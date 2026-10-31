//! A DXF file as its sequence of (group code, value) pairs: ASCII DXF (any line ending, UTF-8 or Windows-1252 text)
//! or binary DXF (R12 one-byte or R13+ two-byte group codes).

use crate::KernelError;

#[derive(Clone, Debug)]
pub enum Value {
    Text(String),
    Real(f64),
    Int(i64),
}

#[derive(Clone, Debug)]
pub struct Pair {
    pub code: i32,
    pub value: Value,
    /// Where the pair starts: the code's line in an ASCII file, its byte offset in a binary one.
    pub at: usize,
    pub binary: bool,
}

const BINARY_SENTINEL: &[u8] = b"AutoCAD Binary DXF\r\n\x1a\0";

fn bad(message: String) -> KernelError {
    KernelError::InvalidInput(message)
}

impl Pair {
    pub fn text(&self) -> String {
        match &self.value {
            Value::Text(s) => s.trim().to_string(),
            Value::Real(r) => r.to_string(),
            Value::Int(i) => i.to_string(),
        }
    }

    pub fn real(&self) -> Result<f64, KernelError> {
        let r = match &self.value {
            Value::Real(r) => *r,
            Value::Int(i) => *i as f64,
            Value::Text(s) => s.trim().parse::<f64>()
                .map_err(|_| bad(format!("{}: group {} should be a number but is '{}'", self.place(), self.code, s.trim())))?,
        };
        if !r.is_finite() {
            return Err(bad(format!("{}: group {} is not a finite number", self.place(), self.code)));
        }
        Ok(r)
    }

    pub fn int(&self) -> Result<i64, KernelError> {
        match &self.value {
            Value::Int(i) => Ok(*i),
            Value::Real(r) if r.fract() == 0.0 && r.is_finite() => Ok(*r as i64),
            Value::Text(s) => {
                let t = s.trim();
                t.parse::<i64>().or_else(|_| match t.parse::<f64>() {
                    Ok(r) if r.fract() == 0.0 && r.is_finite() => Ok(r as i64),
                    _ => Err(()),
                }).map_err(|_| bad(format!("{}: group {} should be a whole number but is '{t}'", self.place(), self.code)))
            }
            _ => Err(bad(format!("{}: group {} should be a whole number", self.place(), self.code))),
        }
    }

    pub fn place(&self) -> String {
        if self.binary { format!("byte {}", self.at) } else { format!("line {}", self.at) }
    }
}

/// Splits a whole DXF file into its pairs.
pub fn pairs(bytes: &[u8]) -> Result<Vec<Pair>, KernelError> {
    if bytes.starts_with(BINARY_SENTINEL) {
        return binary(bytes);
    }
    let bytes = bytes.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(bytes);
    let mut lines = Vec::new();
    for raw in bytes.split(|&b| b == b'\n') {
        let raw = raw.strip_suffix(b"\r").unwrap_or(raw);
        lines.push(raw);
    }
    // A file written with lone CR line ends arrives as one line.
    if lines.len() < 4 && bytes.contains(&b'\r') {
        lines = bytes.split(|&b| b == b'\r').collect();
    }
    let mut out = Vec::with_capacity(lines.len() / 2);
    let mut i = 0;
    while i + 1 < lines.len() {
        let code_text = decode(lines[i]);
        let code_text = code_text.trim();
        if code_text.is_empty() && out.is_empty() {
            i += 1;
            continue;
        }
        let code = code_text.parse::<i32>()
            .map_err(|_| bad(format!("line {}: '{code_text}' is not a DXF group code", i + 1)))?;
        out.push(Pair { code, value: Value::Text(unescape(&decode(lines[i + 1]))), at: i + 1, binary: false });
        i += 2;
    }
    Ok(out)
}

/// UTF-8 when the bytes are valid UTF-8, otherwise Windows-1252 (the code page of almost every pre-2007 DXF).
fn decode(bytes: &[u8]) -> String {
    if let Ok(s) = std::str::from_utf8(bytes) {
        return s.to_string();
    }
    const HIGH: [char; 32] = [
        '\u{20AC}', '\u{81}', '\u{201A}', '\u{192}', '\u{201E}', '\u{2026}', '\u{2020}', '\u{2021}',
        '\u{2C6}', '\u{2030}', '\u{160}', '\u{2039}', '\u{152}', '\u{8D}', '\u{17D}', '\u{8F}',
        '\u{90}', '\u{2018}', '\u{2019}', '\u{201C}', '\u{201D}', '\u{2022}', '\u{2013}', '\u{2014}',
        '\u{2DC}', '\u{2122}', '\u{161}', '\u{203A}', '\u{153}', '\u{9D}', '\u{17E}', '\u{178}',
    ];
    bytes.iter().map(|&b| match b {
        0x80..=0x9F => HIGH[(b - 0x80) as usize],
        _ => b as char,
    }).collect()
}

/// DXF writes characters outside the code page as \U+XXXX.
fn unescape(s: &str) -> String {
    if !s.contains("\\U+") && !s.contains("\\u+") {
        return s.to_string();
    }
    let chars: Vec<char> = s.chars().collect();
    let mut out = String::with_capacity(s.len());
    let mut i = 0;
    while i < chars.len() {
        if chars[i] == '\\' && i + 7 <= chars.len() && (chars[i + 1] == 'U' || chars[i + 1] == 'u') && chars[i + 2] == '+' {
            let hex: String = chars[i + 3..i + 7].iter().collect();
            if let Some(c) = u32::from_str_radix(&hex, 16).ok().and_then(char::from_u32) {
                out.push(c);
                i += 7;
                continue;
            }
        }
        out.push(chars[i]);
        i += 1;
    }
    out
}

#[derive(Clone, Copy)]
enum Kind {
    Text,
    Real,
    I16,
    I32,
    I64,
    Bool,
    Chunk,
}

/// The value type of a group code (DXF Reference, "Group Code Value Types").
fn kind(code: i32) -> Kind {
    match code {
        10..=59 | 110..=149 | 210..=239 | 460..=469 | 1010..=1059 => Kind::Real,
        60..=79 | 170..=179 | 270..=289 | 370..=389 | 400..=409 | 1060..=1070 => Kind::I16,
        90..=99 | 420..=429 | 440..=459 | 1071 => Kind::I32,
        160..=169 => Kind::I64,
        290..=299 => Kind::Bool,
        310..=319 | 1004 => Kind::Chunk,
        _ => Kind::Text,
    }
}

fn binary(bytes: &[u8]) -> Result<Vec<Pair>, KernelError> {
    let mut at = BINARY_SENTINEL.len();
    // R13 and later write two-byte group codes: the first pair is then 00 00 'S' (SECTION), R12's is 00 'S'.
    let wide = bytes.get(at + 1) == Some(&0);
    let take = |at: &mut usize, n: usize| -> Result<&[u8], KernelError> {
        let s = bytes.get(*at..*at + n).ok_or_else(|| bad(format!("binary DXF ends early at byte {}", *at)))?;
        *at += n;
        Ok(s)
    };
    let mut out = Vec::new();
    while at < bytes.len() {
        let start = at;
        let code = if wide {
            i16::from_le_bytes(take(&mut at, 2)?.try_into().unwrap()) as i32
        } else {
            match take(&mut at, 1)?[0] {
                255 => i16::from_le_bytes(take(&mut at, 2)?.try_into().unwrap()) as i32,
                b => b as i32,
            }
        };
        let value = match kind(code) {
            Kind::Text => {
                let end = bytes[at..].iter().position(|&b| b == 0)
                    .ok_or_else(|| bad(format!("binary DXF: unterminated text at byte {at}")))?;
                let s = unescape(&decode(&bytes[at..at + end]));
                at += end + 1;
                Value::Text(s)
            }
            Kind::Real => Value::Real(f64::from_le_bytes(take(&mut at, 8)?.try_into().unwrap())),
            Kind::I16 => Value::Int(i16::from_le_bytes(take(&mut at, 2)?.try_into().unwrap()) as i64),
            Kind::I32 => Value::Int(i32::from_le_bytes(take(&mut at, 4)?.try_into().unwrap()) as i64),
            Kind::I64 => Value::Int(i64::from_le_bytes(take(&mut at, 8)?.try_into().unwrap())),
            Kind::Bool => Value::Int(take(&mut at, 1)?[0] as i64),
            Kind::Chunk => {
                let n = take(&mut at, 1)?[0] as usize;
                take(&mut at, n)?;
                Value::Text(String::new())
            }
        };
        let eof = code == 0 && matches!(&value, Value::Text(s) if s == "EOF");
        out.push(Pair { code, value, at: start, binary: true });
        if eof {
            break;
        }
    }
    Ok(out)
}

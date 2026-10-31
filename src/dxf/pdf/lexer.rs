//! PDF tokens and objects (ISO 32000-1 7.2 lexical conventions, 7.3 objects): white space and comments, numbers,
//! literal strings with escapes and balanced parentheses, hex strings, names with #xx, arrays, dictionaries and
//! indirect references `N G R`.

use std::rc::Rc;

use crate::KernelError;

#[derive(Clone, Debug, PartialEq)]
pub enum Obj {
    Null,
    Bool(bool),
    Int(i64),
    Real(f64),
    Str(Vec<u8>),
    Name(Vec<u8>),
    Array(Vec<Obj>),
    Dict(Dict),
    Stream(Rc<Stream>),
    Ref(u32, u16),
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Dict(pub Vec<(Vec<u8>, Obj)>);

impl Dict {
    pub fn get(&self, key: &[u8]) -> Option<&Obj> {
        self.0.iter().find(|(k, _)| k.as_slice() == key).map(|(_, v)| v)
    }
}

#[derive(Debug, PartialEq)]
pub struct Stream {
    pub dict: Dict,
    /// The stream's bytes as stored (decrypted when the file is encrypted), before any filter.
    pub raw: Vec<u8>,
}

impl Obj {
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Obj::Int(i) => Some(*i as f64),
            Obj::Real(r) if r.is_finite() => Some(*r),
            _ => None,
        }
    }

    pub fn as_i64(&self) -> Option<i64> {
        match self {
            Obj::Int(i) => Some(*i),
            Obj::Real(r) if r.fract() == 0.0 && r.abs() < 9e15 => Some(*r as i64),
            _ => None,
        }
    }

    pub fn as_name(&self) -> Option<&[u8]> {
        match self {
            Obj::Name(n) => Some(n),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&[u8]> {
        match self {
            Obj::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_dict(&self) -> Option<&Dict> {
        match self {
            Obj::Dict(d) => Some(d),
            Obj::Stream(s) => Some(&s.dict),
            _ => None,
        }
    }

    pub fn as_array(&self) -> Option<&[Obj]> {
        match self {
            Obj::Array(a) => Some(a),
            _ => None,
        }
    }
}

pub fn is_white(c: u8) -> bool {
    matches!(c, 0 | 9 | 10 | 12 | 13 | 32)
}

pub fn is_delim(c: u8) -> bool {
    matches!(c, b'(' | b')' | b'<' | b'>' | b'[' | b']' | b'{' | b'}' | b'/' | b'%')
}

pub fn is_regular(c: u8) -> bool {
    !is_white(c) && !is_delim(c)
}

#[derive(Debug, PartialEq)]
pub enum Tok<'a> {
    /// Value, written without a fraction or exponent, and the rounding it carries: half the last written digit's
    /// place, or single-precision rounding of the value when that is larger.
    Num(f64, bool, f64),
    Name(Vec<u8>),
    Str(Vec<u8>),
    ArrOpen,
    ArrClose,
    DictOpen,
    DictClose,
    /// A keyword or operator, or a stray delimiter.
    Kw(&'a [u8]),
    Eof,
}

pub struct Lexer<'a> {
    pub b: &'a [u8],
    pub pos: usize,
}

/// Deeper nesting than this is refused (arrays and dictionaries in one object).
const MAX_DEPTH: usize = 64;

impl<'a> Lexer<'a> {
    pub fn new(b: &'a [u8], pos: usize) -> Self {
        Lexer { b, pos }
    }

    pub fn skip_white(&mut self) {
        while let Some(&c) = self.b.get(self.pos) {
            if is_white(c) {
                self.pos += 1;
            } else if c == b'%' {
                while let Some(&c) = self.b.get(self.pos) {
                    if c == b'\n' || c == b'\r' {
                        break;
                    }
                    self.pos += 1;
                }
            } else {
                break;
            }
        }
    }

    pub fn next(&mut self) -> Tok<'a> {
        self.skip_white();
        let b = self.b;
        let Some(&c) = b.get(self.pos) else { return Tok::Eof };
        match c {
            b'[' => { self.pos += 1; Tok::ArrOpen }
            b']' => { self.pos += 1; Tok::ArrClose }
            b'<' if b.get(self.pos + 1) == Some(&b'<') => { self.pos += 2; Tok::DictOpen }
            b'>' if b.get(self.pos + 1) == Some(&b'>') => { self.pos += 2; Tok::DictClose }
            b'<' => { self.pos += 1; Tok::Str(self.hex_string()) }
            b'(' => { self.pos += 1; Tok::Str(self.literal_string()) }
            b'/' => { self.pos += 1; Tok::Name(self.name()) }
            b')' | b'>' | b'{' | b'}' => { self.pos += 1; Tok::Kw(&b[self.pos - 1..self.pos]) }
            _ => {
                let start = self.pos;
                while self.pos < b.len() && is_regular(b[self.pos]) {
                    self.pos += 1;
                }
                let word = &b[start..self.pos];
                match number(word) {
                    Some((v, int, q)) => Tok::Num(v, int, q),
                    None => Tok::Kw(word),
                }
            }
        }
    }

    fn hex_string(&mut self) -> Vec<u8> {
        let mut out = Vec::new();
        let mut high: Option<u8> = None;
        while let Some(&c) = self.b.get(self.pos) {
            self.pos += 1;
            let v = match c {
                b'>' => break,
                b'0'..=b'9' => c - b'0',
                b'a'..=b'f' => c - b'a' + 10,
                b'A'..=b'F' => c - b'A' + 10,
                _ => continue,
            };
            match high.take() {
                Some(h) => out.push(h << 4 | v),
                None => high = Some(v),
            }
        }
        if let Some(h) = high {
            out.push(h << 4);
        }
        out
    }

    fn literal_string(&mut self) -> Vec<u8> {
        let b = self.b;
        let mut out = Vec::new();
        let mut depth = 1usize;
        while let Some(&c) = b.get(self.pos) {
            self.pos += 1;
            match c {
                b'(' => { depth += 1; out.push(c); }
                b')' => {
                    depth -= 1;
                    if depth == 0 {
                        break;
                    }
                    out.push(c);
                }
                b'\r' => {
                    if b.get(self.pos) == Some(&b'\n') {
                        self.pos += 1;
                    }
                    out.push(b'\n');
                }
                b'\\' => {
                    let Some(&e) = b.get(self.pos) else { break };
                    self.pos += 1;
                    match e {
                        b'n' => out.push(b'\n'),
                        b'r' => out.push(b'\r'),
                        b't' => out.push(b'\t'),
                        b'b' => out.push(8),
                        b'f' => out.push(12),
                        b'0'..=b'7' => {
                            let mut v = (e - b'0') as u32;
                            for _ in 0..2 {
                                match b.get(self.pos) {
                                    Some(&d @ b'0'..=b'7') => { v = v * 8 + (d - b'0') as u32; self.pos += 1; }
                                    _ => break,
                                }
                            }
                            out.push(v as u8);
                        }
                        b'\r' => {
                            if b.get(self.pos) == Some(&b'\n') {
                                self.pos += 1;
                            }
                        }
                        b'\n' => {}
                        other => out.push(other),
                    }
                }
                _ => out.push(c),
            }
        }
        out
    }

    fn name(&mut self) -> Vec<u8> {
        let b = self.b;
        let mut out = Vec::new();
        while self.pos < b.len() && is_regular(b[self.pos]) {
            let c = b[self.pos];
            if c == b'#' {
                let h = |x: u8| (x as char).to_digit(16);
                if let (Some(a), Some(d)) = (b.get(self.pos + 1).and_then(|&x| h(x)), b.get(self.pos + 2).and_then(|&x| h(x))) {
                    out.push((a * 16 + d) as u8);
                    self.pos += 3;
                    continue;
                }
            }
            out.push(c);
            self.pos += 1;
        }
        out
    }
}

/// A PDF number ([+-] digits [. digits], tolerating an exponent and doubled signs that some writers emit).
fn number(w: &[u8]) -> Option<(f64, bool, f64)> {
    let first = *w.first()?;
    if !(first.is_ascii_digit() || matches!(first, b'+' | b'-' | b'.')) {
        return None;
    }
    let mut i = 0;
    let mut negative = false;
    while i < w.len() && matches!(w[i], b'+' | b'-') {
        negative ^= w[i] == b'-';
        i += 1;
    }
    let mut mantissa = 0f64;
    let mut digits = 0;
    let mut decimals: i32 = 0;
    let mut dot = false;
    while i < w.len() {
        match w[i] {
            d @ b'0'..=b'9' => {
                mantissa = mantissa * 10.0 + (d - b'0') as f64;
                digits += 1;
                if dot {
                    decimals += 1;
                }
            }
            b'.' if !dot => dot = true,
            _ => break,
        }
        i += 1;
    }
    let mut exponent: i32 = 0;
    let mut has_exponent = false;
    if digits > 0 && i < w.len() && matches!(w[i], b'e' | b'E') {
        let mut j = i + 1;
        let mut eneg = false;
        if j < w.len() && matches!(w[j], b'+' | b'-') {
            eneg = w[j] == b'-';
            j += 1;
        }
        let start = j;
        let mut e: i32 = 0;
        while j < w.len() && w[j].is_ascii_digit() {
            e = (e * 10 + (w[j] - b'0') as i32).min(400);
            j += 1;
        }
        if j > start {
            exponent = if eneg { -e } else { e };
            has_exponent = true;
            i = j;
        }
    }
    if digits == 0 {
        // "-", "." and the like: some writers emit them for zero.
        return if i == w.len() { Some((0.0, true, 0.5)) } else { None };
    }
    if i != w.len() {
        // Trailing junk such as "1.2.3": keep the number written so far.
        if w[i..].iter().any(|c| c.is_ascii_alphabetic() && !matches!(c, b'e' | b'E')) {
            return None;
        }
    }
    let scale = exponent - decimals;
    let mut v = if scale >= 0 { mantissa * 10f64.powi(scale) } else { mantissa / 10f64.powi(-scale) };
    if negative {
        v = -v;
    }
    if !v.is_finite() {
        return None;
    }
    let int = !dot && !has_exponent && v.abs() < 9e15;
    // Many writers hold coordinates in single precision and print more digits than that carries.
    let rounding = (0.5 * 10f64.powi(scale.clamp(-300, 300))).max(v.abs() * f64::from(f32::EPSILON) / 2.0);
    Some((v, int, rounding))
}

fn bad(m: impl Into<String>) -> KernelError {
    KernelError::InvalidInput(m.into())
}

/// One object starting at the lexer's position. `refs` allows `N G R` (not in content streams).
pub fn parse(lx: &mut Lexer, refs: bool) -> Result<Obj, KernelError> {
    let tok = lx.next();
    build(lx, tok, refs, 0)?.ok_or_else(|| bad("a PDF object was expected"))
}

/// The object a token begins, or None for a keyword other than true/false/null.
pub fn build(lx: &mut Lexer, tok: Tok, refs: bool, depth: usize) -> Result<Option<Obj>, KernelError> {
    if depth > MAX_DEPTH {
        return Err(bad("PDF objects are nested too deeply"));
    }
    Ok(Some(match tok {
        Tok::Num(v, int, _) => {
            if int && refs && v >= 0.0 && v <= u32::MAX as f64 {
                let save = lx.pos;
                if let Tok::Num(g, true, _) = lx.next() {
                    if (0.0..=65535.0).contains(&g) && lx.next() == Tok::Kw(b"R") {
                        return Ok(Some(Obj::Ref(v as u32, g as u16)));
                    }
                }
                lx.pos = save;
            }
            if int { Obj::Int(v as i64) } else { Obj::Real(v) }
        }
        Tok::Name(n) => Obj::Name(n),
        Tok::Str(s) => Obj::Str(s),
        Tok::ArrOpen => {
            let mut items = Vec::new();
            loop {
                let t = lx.next();
                match t {
                    Tok::ArrClose => break,
                    Tok::Eof => return Err(bad("a PDF array is not closed")),
                    Tok::DictClose => return Err(bad("a PDF array is closed by '>>'")),
                    t => {
                        if let Some(o) = build(lx, t, refs, depth + 1)? {
                            items.push(o);
                        }
                    }
                }
            }
            Obj::Array(items)
        }
        Tok::DictOpen => {
            let mut entries: Vec<(Vec<u8>, Obj)> = Vec::new();
            loop {
                match lx.next() {
                    Tok::DictClose => break,
                    Tok::Eof => return Err(bad("a PDF dictionary is not closed")),
                    Tok::Name(k) => {
                        let save = lx.pos;
                        let t = lx.next();
                        if t == Tok::DictClose {
                            entries.push((k, Obj::Null));
                            break;
                        }
                        if matches!(t, Tok::Name(_)) || matches!(t, Tok::Kw(_)) {
                            // A value may be a name; a keyword other than true/false/null means a missing value.
                            match build(lx, t, refs, depth + 1)? {
                                Some(o) => entries.push((k, o)),
                                None => {
                                    lx.pos = save;
                                    let _ = lx.next();
                                    entries.push((k, Obj::Null));
                                }
                            }
                        } else if let Some(o) = build(lx, t, refs, depth + 1)? {
                            entries.push((k, o));
                        }
                    }
                    Tok::ArrClose => return Err(bad("a PDF dictionary is closed by ']'")),
                    t => {
                        // A non-name key: skip the stray object.
                        let _ = build(lx, t, refs, depth + 1)?;
                    }
                }
            }
            Obj::Dict(Dict(entries))
        }
        Tok::Kw(b"true") => Obj::Bool(true),
        Tok::Kw(b"false") => Obj::Bool(false),
        Tok::Kw(b"null") => Obj::Null,
        Tok::Kw(_) => return Ok(None),
        Tok::ArrClose | Tok::DictClose | Tok::Eof => return Ok(None),
    }))
}

/// A text string (7.9.2.2): UTF-16BE with its byte order mark, UTF-8 with its mark, else PDFDocEncoding (read as
/// Latin-1, which it matches for the printable characters).
pub fn text(s: &[u8]) -> String {
    if let Some(rest) = s.strip_prefix(&[0xFE, 0xFF]) {
        let units: Vec<u16> = rest.chunks_exact(2).map(|c| u16::from_be_bytes([c[0], c[1]])).collect();
        return String::from_utf16_lossy(&units);
    }
    if let Some(rest) = s.strip_prefix(&[0xEF, 0xBB, 0xBF]) {
        return String::from_utf8_lossy(rest).into_owned();
    }
    s.iter().map(|&b| b as char).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one(s: &str) -> Obj {
        parse(&mut Lexer::new(s.as_bytes(), 0), true).unwrap()
    }

    #[test]
    fn objects() {
        assert_eq!(one("  % comment\n 12 0 R"), Obj::Ref(12, 0));
        assert_eq!(one("-.5"), Obj::Real(-0.5));
        assert_eq!(one("(a(b)c\\)\\101\\\nd)"), Obj::Str(b"a(b)c)Ad".to_vec()));
        assert_eq!(one("<48 6>"), Obj::Str(b"H`".to_vec()));
        assert_eq!(one("/A#20B"), Obj::Name(b"A B".to_vec()));
        let d = one("<</Type/Page/Kids[1 0 R 2 0 R]/N 3/Sub<</X true>>/Empty>>");
        let d = d.as_dict().unwrap();
        assert_eq!(d.get(b"Kids").unwrap().as_array().unwrap(), &[Obj::Ref(1, 0), Obj::Ref(2, 0)]);
        assert_eq!(d.get(b"N"), Some(&Obj::Int(3)));
        assert_eq!(d.get(b"Empty"), Some(&Obj::Null));
    }
}

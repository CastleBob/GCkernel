//! IGES 5.3 fixed records and Hollerith fields, independently written from
//! sections 2.2.3 and 2.2.4. No interpretation of unknown geometry as display data.
use super::{bad, Result};
use std::collections::BTreeMap;
pub const MAX_BYTES: usize = 64 * 1024 * 1024;
#[derive(Clone, Debug)]
pub enum Value {
    Number(f64),
    Text(String),
    Empty,
}
impl Value {
    pub fn number(&self) -> Result<f64> {
        if let Self::Number(v) = self {
            Ok(*v)
        } else {
            bad("expected numeric parameter")
        }
    }
}
pub struct Entity {
    pub kind: usize,
    pub form: usize,
    pub transform: usize,
    pub label: String,
    pub values: Vec<Value>,
}
pub struct File {
    pub entities: BTreeMap<usize, Entity>,
    pub scale: f64,
}
fn integer(bytes: &[u8]) -> Result<usize> {
    let s = std::str::from_utf8(bytes)
        .map_err(|_| super::error("non-ASCII field"))?
        .trim();
    if s.is_empty() {
        return Ok(0);
    }
    s.parse::<usize>()
        .map_err(|_| super::error("invalid integer field"))
}
fn hollerith(data: &[u8], at: &mut usize) -> Result<Option<Vec<u8>>> {
    while data.get(*at) == Some(&b' ') {
        *at += 1;
    }
    let start = *at;
    let mut end = start;
    while data.get(end).is_some_and(u8::is_ascii_digit) {
        end += 1;
    }
    if end == start || !matches!(data.get(end), Some(b'H' | b'h')) {
        return Ok(None);
    }
    let count = integer(&data[start..end])?;
    if count > MAX_BYTES || end + 1 + count > data.len() {
        return bad("truncated Hollerith string");
    }
    *at = end + 1 + count;
    Ok(Some(data[end + 1..*at].to_vec()))
}
fn values(data: &[u8], delimiter: u8, terminator: u8) -> Result<Vec<Value>> {
    let mut out = Vec::new();
    let mut at = 0;
    loop {
        if out.len() >= 4_000_000 {
            return bad("parameter budget exceeded");
        }
        while data.get(at) == Some(&b' ') {
            at += 1;
        }
        if at == data.len() {
            return bad("missing parameter terminator");
        }
        let value = if let Some(text) = hollerith(data, &mut at)? {
            Value::Text(String::from_utf8(text).map_err(|_| super::error("invalid string"))?)
        } else {
            let start = at;
            while at < data.len() && data[at] != delimiter && data[at] != terminator {
                at += 1;
            }
            let text = std::str::from_utf8(&data[start..at]).unwrap().trim();
            if text.is_empty() {
                Value::Empty
            } else {
                let n = text
                    .replace(['D', 'd'], "E")
                    .parse::<f64>()
                    .map_err(|_| super::error("invalid real parameter"))?;
                if !n.is_finite() {
                    return bad("nonfinite real parameter");
                }
                Value::Number(n)
            }
        };
        out.push(value);
        while data.get(at) == Some(&b' ') {
            at += 1;
        }
        match data.get(at) {
            Some(v) if *v == delimiter => at += 1,
            Some(v) if *v == terminator => {
                if data[at + 1..].iter().any(|&v| v != b' ') {
                    return bad("trailing parameter data");
                }
                return Ok(out);
            }
            _ => return bad("invalid parameter delimiter"),
        }
    }
}
pub fn parse(bytes: &[u8]) -> Result<File> {
    if bytes.len() > MAX_BYTES || bytes.is_empty() {
        return bad("empty or oversized file");
    }
    let mut sections: [Vec<Vec<u8>>; 5] = std::array::from_fn(|_| Vec::new());
    let mut previous = 0;
    let lines = bytes.split(|&v| v == b'\n').collect::<Vec<_>>();
    for (i, line) in lines.iter().enumerate() {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if i + 1 == lines.len() && line.is_empty() {
            continue;
        }
        if line.len() != 80 || line.iter().any(|b| !matches!(b, 32..=126)) {
            return bad("IGES ASCII records must contain exactly 80 printable columns");
        }
        let section = match line[72] {
            b'S' => 0,
            b'G' => 1,
            b'D' => 2,
            b'P' => 3,
            b'T' => 4,
            _ => return bad("invalid section code"),
        };
        if section < previous || section > previous + 1 {
            return bad("invalid section order");
        }
        previous = section;
        if integer(&line[73..80])? != sections[section].len() + 1 {
            return bad("nonconsecutive section sequence");
        }
        sections[section].push(line.to_vec());
    }
    if sections.iter().any(Vec::is_empty)
        || sections[2].len() % 2 != 0
        || sections[4].len() != 1
        || sections[2].len() > 1_000_000
    {
        return bad("missing or invalid sections");
    }
    for i in 0..4 {
        let t = &sections[4][0];
        if t[8 * i] != b"SGDP"[i] || integer(&t[8 * i + 1..8 * i + 8])? != sections[i].len() {
            return bad("terminate section counts disagree");
        }
    }
    let global = sections[1]
        .iter()
        .flat_map(|l| l[..72].iter().copied())
        .collect::<Vec<_>>();
    let mut at = 0;
    let delimiter = if let Some(s) = hollerith(&global, &mut at)? {
        if s.len() != 1 {
            return bad("invalid parameter separator");
        }
        s[0]
    } else {
        b','
    };
    if global.get(at) != Some(&delimiter) {
        return bad("invalid global parameter separator");
    }
    at += 1;
    let terminator = if let Some(s) = hollerith(&global, &mut at)? {
        if s.len() != 1 {
            return bad("invalid record separator");
        }
        s[0]
    } else {
        b';'
    };
    if delimiter == terminator
        || delimiter.is_ascii_alphanumeric()
        || terminator.is_ascii_alphanumeric()
        || delimiter == b' '
        || terminator == b' '
    {
        return bad("invalid delimiters");
    }
    let g = values(&global, delimiter, terminator)?;
    if g.len() < 24 || g.len() > 26 {
        return bad("invalid IGES global parameter count");
    }
    let model = if matches!(g[12], Value::Empty) {
        1.
    } else {
        g[12].number()?
    };
    let unit = if matches!(g[13], Value::Empty) {
        1.
    } else {
        g[13].number()?
    };
    let unit = if unit == 3. {
        let Value::Text(name) = &g[14] else {
            return bad("named units require a units string");
        };
        match name.as_str() {
            "IN" | "INCH" => 1.,
            "MM" => 2.,
            "FT" => 4.,
            "MI" => 5.,
            "M" => 6.,
            "KM" => 7.,
            "MIL" => 8.,
            "UM" => 9.,
            "CM" => 10.,
            "UIN" => 11.,
            _ => return super::unsupported("unknown named IGES units"),
        }
    } else {
        unit
    };
    let inches = match unit {
        1. => 1.,
        2. => 1. / 25.4,
        4. => 12.,
        5. => 63360.,
        6. => 1000. / 25.4,
        7. => 1_000_000. / 25.4,
        8. => 0.001,
        9. => 0.001 / 25.4,
        10. => 10. / 25.4,
        11. => 1e-6,
        _ => return super::unsupported("custom/unknown IGES units"),
    };
    let scale = inches / model;
    if model <= 0. || !scale.is_finite() || scale <= 0. {
        return bad("invalid model space scale");
    }
    let mut entities = BTreeMap::new();
    let mut consumed = vec![false; sections[3].len()];
    for (index, pair) in sections[2].chunks(2).enumerate() {
        let a = &pair[0];
        let b = &pair[1];
        let id = 2 * index + 1;
        let field = |line: &[u8], i: usize| integer(&line[i * 8..(i + 1) * 8]);
        let kind = field(a, 0)?;
        let start = field(a, 1)?;
        let count = field(b, 3)?;
        if kind != field(b, 0)? || start == 0 || count == 0 || start - 1 + count > sections[3].len()
        {
            return bad("invalid directory parameter range");
        }
        if field(a, 2)? != 0 {
            return super::unsupported("IGES external structure reference");
        }
        let mut data = Vec::new();
        for i in start - 1..start - 1 + count {
            let line = &sections[3][i];
            if consumed[i] || integer(&line[64..72])? != id {
                return bad("parameter records have ambiguous directory ownership");
            }
            consumed[i] = true;
            data.extend_from_slice(&line[..64]);
        }
        let params = values(&data, delimiter, terminator)?;
        if params[0].number()? != kind as f64 {
            return bad("entity type disagrees with directory");
        }
        entities.insert(
            id,
            Entity {
                kind,
                form: field(b, 4)?,
                transform: field(a, 6)?,
                label: std::str::from_utf8(&b[56..64]).unwrap().trim().to_owned(),
                values: params[1..].to_vec(),
            },
        );
    }
    if consumed.iter().any(|v| !*v) {
        return bad("orphan parameter record");
    }
    Ok(File { entities, scale })
}

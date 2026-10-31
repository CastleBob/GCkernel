//! Independently written Part 21 syntax reader. References: ISO 10303-21 public
//! grammar, not another translator's implementation. No geometry decisions here.
use super::{invalid, unsupported, Result};
use crate::KernelError;
use std::collections::BTreeMap;

pub const MAX_BYTES: usize = 64 * 1024 * 1024;
const MAX_VALUES: usize = 4_000_000;
const MAX_ENTITIES: usize = 500_000;

#[derive(Clone, Debug)]
pub enum Value {
    Number(f64),
    Text(String),
    Enum(String),
    Ref(u64),
    List(Vec<Value>),
    Typed(String, Vec<Value>),
    Omitted,
    Derived,
    Binary,
}
impl Value {
    pub fn number(&self) -> Result<f64> {
        if let Self::Number(n) = self {
            Ok(*n)
        } else {
            invalid("expected a number")
        }
    }
    pub fn text(&self) -> Result<&str> {
        if let Self::Text(s) = self {
            Ok(s)
        } else {
            invalid("expected a string")
        }
    }
    pub fn reference(&self) -> Result<u64> {
        if let Self::Ref(n) = self {
            Ok(*n)
        } else {
            invalid("expected an entity reference")
        }
    }
    pub fn list(&self) -> Result<&[Value]> {
        if let Self::List(a) = self {
            Ok(a)
        } else {
            invalid("expected a list")
        }
    }
    pub fn enumeration(&self) -> Result<&str> {
        if let Self::Enum(s) = self {
            Ok(s)
        } else {
            invalid("expected an enumeration")
        }
    }
    pub fn flag(&self) -> Result<bool> {
        match self.enumeration()? {
            "T" => Ok(true),
            "F" => Ok(false),
            _ => invalid("expected .T. or .F."),
        }
    }
    fn refs(&self, out: &mut Vec<u64>) {
        match self {
            Self::Ref(n) => out.push(*n),
            Self::List(a) | Self::Typed(_, a) => {
                for v in a {
                    v.refs(out);
                }
            }
            _ => {}
        }
    }
}
#[derive(Clone, Debug)]
pub struct Entity {
    pub components: BTreeMap<String, Vec<Value>>,
}
impl Entity {
    pub fn has(&self, name: &str) -> bool {
        self.components.contains_key(name)
    }
    pub fn args(&self, name: &str, count: usize) -> Result<&[Value]> {
        let Some(a) = self.components.get(name) else {
            return unsupported(format!("expected {name}, found {}", self.kind()));
        };
        if a.len() != count {
            return invalid(format!(
                "{name} requires {count} attributes, found {}",
                a.len()
            ));
        }
        Ok(a)
    }
    pub fn kind(&self) -> String {
        self.components
            .keys()
            .cloned()
            .collect::<Vec<_>>()
            .join("+")
    }
}
pub struct Exchange {
    pub entities: BTreeMap<u64, Entity>,
    pub schemas: Vec<String>,
}

struct Parser<'a> {
    b: &'a [u8],
    at: usize,
    values: usize,
}
impl<'a> Parser<'a> {
    fn error<T>(&self, message: &str) -> Result<T> {
        invalid(format!("STEP byte {}: {message}", self.at))
    }
    fn space(&mut self) -> Result<()> {
        loop {
            while self.b.get(self.at).is_some_and(u8::is_ascii_whitespace) {
                self.at += 1;
            }
            if self.b.get(self.at..self.at + 2) != Some(b"/*") {
                return Ok(());
            }
            self.at += 2;
            while self.b.get(self.at..self.at + 2) != Some(b"*/") {
                if self.at >= self.b.len() {
                    return self.error("unterminated comment");
                }
                self.at += 1;
            }
            self.at += 2;
        }
    }
    fn peek(&mut self) -> Result<u8> {
        self.space()?;
        self.b
            .get(self.at)
            .copied()
            .ok_or_else(|| KernelError::InvalidInput(format!("STEP ended at byte {}", self.at)))
    }
    fn punctuation(&mut self, c: u8) -> Result<()> {
        if self.peek()? != c {
            return self.error(&format!("expected '{}'", c as char));
        }
        self.at += 1;
        Ok(())
    }
    fn word(&mut self) -> Result<String> {
        self.space()?;
        let start = self.at;
        if !self.b.get(start).is_some_and(u8::is_ascii_alphabetic) {
            return self.error("expected a keyword");
        }
        while self
            .b
            .get(self.at)
            .is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_' || *c == b'-')
        {
            self.at += 1;
        }
        Ok(std::str::from_utf8(&self.b[start..self.at])
            .unwrap()
            .to_ascii_uppercase())
    }
    fn keyword(&mut self, wanted: &str) -> Result<()> {
        let actual = self.word()?;
        if actual != wanted {
            return self.error(&format!("expected {wanted}, found {actual}"));
        }
        Ok(())
    }
    fn id(&mut self) -> Result<u64> {
        self.punctuation(b'#')?;
        let start = self.at;
        while self.b.get(self.at).is_some_and(u8::is_ascii_digit) {
            self.at += 1;
        }
        let id = std::str::from_utf8(&self.b[start..self.at])
            .unwrap()
            .parse::<u64>()
            .map_err(|_| {
                KernelError::InvalidInput(format!("invalid STEP entity identifier at byte {start}"))
            })?;
        if id == 0 || id > 9_007_199_254_740_991 {
            return self.error("entity identifier outside supported positive exact-integer range");
        }
        Ok(id)
    }
    fn string(&mut self) -> Result<String> {
        self.punctuation(b'\'')?;
        let start = self.at;
        let mut out = String::new();
        loop {
            if self.at - start > 32767 {
                return self.error("STEP string length limit exceeded");
            }
            let Some(&c) = self.b.get(self.at) else {
                return self.error("unterminated string");
            };
            self.at += 1;
            match c {
                b'\'' => {
                    if self.b.get(self.at) == Some(&b'\'') {
                        self.at += 1;
                        out.push('\'');
                    } else {
                        return Ok(out);
                    }
                }
                b'\\' => {
                    if self.b.get(self.at) == Some(&b'\\') {
                        self.at += 1;
                        out.push('\\');
                        continue;
                    }
                    let remaining = &self.b[self.at..];
                    if remaining.starts_with(b"X2\\") || remaining.starts_with(b"X4\\") {
                        let width = if remaining[1] == b'2' { 4 } else { 8 };
                        self.at += 3;
                        let mut codes = Vec::new();
                        while self.b.get(self.at..self.at + 4) != Some(b"\\X0\\") {
                            let Some(s) = self.b.get(self.at..self.at + width) else {
                                return self.error("unterminated Unicode escape");
                            };
                            let h = std::str::from_utf8(s)
                                .ok()
                                .and_then(|s| u32::from_str_radix(s, 16).ok())
                                .ok_or_else(|| {
                                    KernelError::InvalidInput(
                                        "invalid STEP Unicode hex digits".into(),
                                    )
                                })?;
                            codes.push(h);
                            self.at += width;
                        }
                        self.at += 4;
                        if codes.is_empty() {
                            return self.error("empty Unicode escape");
                        }
                        if width == 4 {
                            for c in char::decode_utf16(codes.iter().map(|c| *c as u16)) {
                                out.push(c.map_err(|_| {
                                    KernelError::InvalidInput(
                                        "invalid STEP UTF-16 surrogate".into(),
                                    )
                                })?);
                            }
                        } else {
                            for c in codes {
                                out.push(char::from_u32(c).ok_or_else(|| {
                                    KernelError::InvalidInput("invalid STEP Unicode scalar".into())
                                })?);
                            }
                        }
                    } else if remaining.starts_with(b"X\\") {
                        self.at += 2;
                        let Some(hex) = self.b.get(self.at..self.at + 2) else {
                            return self.error("short hexadecimal string escape");
                        };
                        let n = std::str::from_utf8(hex)
                            .ok()
                            .and_then(|s| u8::from_str_radix(s, 16).ok())
                            .ok_or_else(|| {
                                KernelError::InvalidInput("invalid STEP string escape".into())
                            })?;
                        out.push(n as char);
                        self.at += 2;
                    } else {
                        return unsupported(format!(
                            "STEP string control directive at byte {}",
                            self.at - 1
                        ));
                    }
                }
                32..=126 => out.push(c as char),
                _ => {
                    let width = if c < 224 {
                        2
                    } else if c < 240 {
                        3
                    } else {
                        4
                    };
                    let bytes = self
                        .b
                        .get(self.at - 1..self.at - 1 + width)
                        .ok_or_else(|| {
                            KernelError::InvalidInput("short STEP UTF-8 sequence".into())
                        })?;
                    let rest = std::str::from_utf8(bytes)
                        .map_err(|_| KernelError::InvalidInput("invalid STEP UTF-8".into()))?;
                    let ch = rest.chars().next().unwrap();
                    out.push(ch);
                    self.at += ch.len_utf8() - 1;
                }
            }
        }
    }
    fn list(&mut self, depth: usize) -> Result<Vec<Value>> {
        self.punctuation(b'(')?;
        let mut values = Vec::new();
        if self.peek()? == b')' {
            self.at += 1;
            return Ok(values);
        }
        loop {
            values.push(self.value(depth + 1)?);
            match self.peek()? {
                b')' => {
                    self.at += 1;
                    return Ok(values);
                }
                b',' => self.at += 1,
                _ => return self.error("expected comma or closing parenthesis"),
            }
        }
    }
    fn value(&mut self, depth: usize) -> Result<Value> {
        self.values += 1;
        if depth > 64 || self.values > MAX_VALUES {
            return self.error("STEP nesting/value limit exceeded");
        }
        match self.peek()? {
            b'#' => Ok(Value::Ref(self.id()?)),
            b'\'' => Ok(Value::Text(self.string()?)),
            b'(' => Ok(Value::List(self.list(depth)?)),
            b'$' => {
                self.at += 1;
                Ok(Value::Omitted)
            }
            b'*' => {
                self.at += 1;
                Ok(Value::Derived)
            }
            b'.' => {
                self.at += 1;
                let s = self.word()?;
                self.punctuation(b'.')?;
                Ok(Value::Enum(s))
            }
            b'"' => {
                self.at += 1;
                let start = self.at;
                while self.b.get(self.at).is_some_and(u8::is_ascii_hexdigit) {
                    self.at += 1;
                }
                if self.at == start || !matches!(self.b[start], b'0'..=b'3') {
                    return self.error("invalid binary encoding");
                }
                self.punctuation(b'"')?;
                Ok(Value::Binary)
            }
            b'+' | b'-' | b'0'..=b'9' => {
                let start = self.at;
                if matches!(self.b[self.at], b'+' | b'-') {
                    self.at += 1;
                }
                let digits = self.at;
                while self.b.get(self.at).is_some_and(u8::is_ascii_digit) {
                    self.at += 1;
                }
                if self.at == digits {
                    return self.error("number has no digits");
                }
                if self.b.get(self.at) == Some(&b'.') {
                    self.at += 1;
                    while self.b.get(self.at).is_some_and(u8::is_ascii_digit) {
                        self.at += 1;
                    }
                    if matches!(self.b.get(self.at), Some(b'E' | b'e')) {
                        self.at += 1;
                        if matches!(self.b.get(self.at), Some(b'+' | b'-')) {
                            self.at += 1;
                        }
                        let exponent = self.at;
                        while self.b.get(self.at).is_some_and(u8::is_ascii_digit) {
                            self.at += 1;
                        }
                        if self.at == exponent {
                            return self.error("empty exponent");
                        }
                    }
                }
                let n = std::str::from_utf8(&self.b[start..self.at])
                    .unwrap()
                    .parse::<f64>()
                    .map_err(|_| {
                        KernelError::InvalidInput(format!("invalid STEP number at byte {start}"))
                    })?;
                if !n.is_finite() {
                    return self.error("non-finite number");
                }
                Ok(Value::Number(n))
            }
            _ => {
                let name = self.word()?;
                Ok(Value::Typed(name, self.list(depth)?))
            }
        }
    }
    fn entity(&mut self) -> Result<Entity> {
        let complex = self.peek()? == b'(';
        if complex {
            self.at += 1;
        }
        let mut components = BTreeMap::new();
        loop {
            let name = self.word()?;
            let args = self.list(0)?;
            if components.insert(name, args).is_some() {
                return self.error("duplicate complex entity component");
            }
            if !complex {
                break;
            }
            if self.peek()? == b')' {
                self.at += 1;
                break;
            }
        }
        self.punctuation(b';')?;
        Ok(Entity { components })
    }
}

pub fn parse(bytes: &[u8]) -> Result<Exchange> {
    if bytes.len() > MAX_BYTES {
        return invalid("STEP file exceeds the 64 MiB limit");
    }
    // Part 21 section 5.2: control octets are ignored even inside tokens/strings.
    // In particular CAD exporters wrap long keywords and strings at line breaks.
    let filtered: Vec<u8> = bytes
        .iter()
        .copied()
        .filter(|b| *b >= 32 && *b != 127)
        .collect();
    std::str::from_utf8(&filtered)
        .map_err(|_| KernelError::InvalidInput("invalid STEP UTF-8 encoding".into()))?;
    let bytes = filtered.as_slice();
    let mut p = Parser {
        b: bytes,
        at: 0,
        values: 0,
    };
    p.keyword("ISO-10303-21")?;
    p.punctuation(b';')?;
    p.keyword("HEADER")?;
    p.punctuation(b';')?;
    let mut header = BTreeMap::new();
    loop {
        let name = p.word()?;
        if name == "ENDSEC" {
            p.punctuation(b';')?;
            break;
        }
        let args = p.list(0)?;
        p.punctuation(b';')?;
        if header.insert(name, args).is_some() {
            return invalid("duplicate STEP header record");
        }
    }
    for required in ["FILE_DESCRIPTION", "FILE_NAME", "FILE_SCHEMA"] {
        if !header.contains_key(required) {
            return invalid(format!("STEP header lacks {required}"));
        }
    }
    let schema = &header["FILE_SCHEMA"];
    if schema.len() != 1 {
        return invalid("FILE_SCHEMA requires one list");
    }
    let schemas = schema[0]
        .list()?
        .iter()
        .map(|v| Ok(v.text()?.to_owned()))
        .collect::<Result<Vec<_>>>()?;
    if schemas.is_empty() {
        return invalid("FILE_SCHEMA is empty");
    }
    p.keyword("DATA")?;
    if p.peek()? != b';' {
        return unsupported("named/multiple STEP DATA sections");
    }
    p.punctuation(b';')?;
    let mut entities = BTreeMap::new();
    loop {
        if p.peek()? != b'#' {
            p.keyword("ENDSEC")?;
            p.punctuation(b';')?;
            break;
        }
        if entities.len() >= MAX_ENTITIES {
            return invalid("STEP entity limit exceeded");
        }
        let id = p.id()?;
        p.punctuation(b'=')?;
        if entities.insert(id, p.entity()?).is_some() {
            return invalid(format!("duplicate STEP entity #{id}"));
        }
    }
    p.keyword("END-ISO-10303-21")?;
    p.punctuation(b';')?;
    p.space()?;
    if p.at != bytes.len() {
        return p.error("trailing content after STEP terminator");
    }
    for (id, e) in &entities {
        let mut refs = Vec::new();
        for a in e.components.values() {
            for v in a {
                v.refs(&mut refs);
            }
        }
        for r in refs {
            if !entities.contains_key(&r) {
                return invalid(format!("STEP #{id} references missing #{r}"));
            }
        }
    }
    Ok(Exchange { entities, schemas })
}

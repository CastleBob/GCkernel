//! The document's objects (ISO 32000-1 7.5 file structure): the header, `startxref`, classic cross-reference tables
//! (7.5.4) with their trailers and /Prev chains, cross-reference streams (7.5.8) including hybrid files' /XRefStm,
//! object streams (7.5.7), and, when the table is missing or wrong, a rebuilt table from scanning the file for
//! `N G obj`. Streams are decrypted as they are read (7.6) and decoded on request within a size budget.

use std::collections::{HashMap, HashSet};
use std::rc::Rc;

use super::crypt::{Crypt, Method};
use super::filters;
use super::lexer::{self, Dict, Lexer, Obj, Stream, Tok};
use crate::KernelError;

/// All stream decoding together may produce at most this much data.
pub const DECODE_BUDGET: usize = 512 << 20;
const MAX_RESOLVE_DEPTH: usize = 32;
const MAX_XREF_SECTIONS: usize = 1000;

#[derive(Clone, Copy, Debug)]
enum Entry {
    At(usize),
    InStream(u32, u32),
}

pub struct Doc<'a> {
    data: &'a [u8],
    entries: HashMap<u32, Entry>,
    /// Object positions found by scanning, used when the table's entry is wrong (built on first need).
    scanned: Option<HashMap<u32, usize>>,
    /// Object streams found by scanning (rebuilt files), searched for objects the table lacks.
    scanned_streams: Vec<u32>,
    pub trailer: Dict,
    cache: HashMap<u32, Rc<Obj>>,
    object_streams: HashMap<u32, Rc<(Vec<u8>, Vec<(u32, usize)>)>>,
    crypt: Option<Crypt>,
    encrypt_num: Option<u32>,
    budget: usize,
    depth: usize,
    pub notes: Vec<String>,
    pub version: String,
}

fn bad(m: impl Into<String>) -> KernelError {
    KernelError::InvalidInput(m.into())
}

fn find(hay: &[u8], needle: &[u8], from: usize) -> Option<usize> {
    if from >= hay.len() || needle.is_empty() {
        return None;
    }
    hay[from..].windows(needle.len()).position(|w| w == needle).map(|p| p + from)
}

fn rfind(hay: &[u8], needle: &[u8]) -> Option<usize> {
    hay.windows(needle.len()).rposition(|w| w == needle)
}

impl<'a> Doc<'a> {
    pub fn open(data: &'a [u8]) -> Result<Doc<'a>, KernelError> {
        let header = find(&data[..data.len().min(1024)], b"%PDF-", 0).ok_or_else(|| bad("the file is not a PDF"))?;
        let version: String = data[header + 5..].iter().take(8).take_while(|c| c.is_ascii_digit() || **c == b'.')
            .map(|&c| c as char).collect();
        let mut doc = Doc {
            data, entries: HashMap::new(), scanned: None, scanned_streams: Vec::new(), trailer: Dict::default(),
            cache: HashMap::new(), object_streams: HashMap::new(), crypt: None, encrypt_num: None,
            budget: DECODE_BUDGET, depth: 0, notes: Vec::new(), version,
        };
        // Decryption is set up before the catalog is checked: it may sit in an encrypted object stream. A password
        // refusal stands; a damaged encryption dictionary may be the table's fault, so the file is then scanned.
        let mut from_table = doc.load_table(header).is_ok();
        if from_table {
            match doc.setup_crypt() {
                Ok(()) => from_table = doc.catalog_ok(),
                Err(e @ KernelError::Unsupported(_)) => return Err(e),
                Err(_) => from_table = false,
            }
        }
        if !from_table {
            doc.entries.clear();
            doc.trailer = Dict::default();
            doc.cache.clear();
            doc.object_streams.clear();
            doc.crypt = None;
            doc.encrypt_num = None;
            doc.rebuild()?;
            doc.notes.push("the PDF's cross-reference table was missing or damaged; its objects were found by scanning the file".into());
            doc.setup_crypt()?;
            if !doc.catalog_ok() {
                return Err(bad("the PDF has no readable document catalog"));
            }
        }
        if let Some(c) = &doc.crypt {
            let d = c.description.clone();
            doc.notes.push(d);
        }
        Ok(doc)
    }

    fn catalog_ok(&mut self) -> bool {
        let Some(root) = self.trailer.get(b"Root").cloned() else { return false };
        match self.resolve(&root) {
            Ok(o) => o.as_dict().is_some_and(|d| d.get(b"Pages").is_some()),
            Err(_) => false,
        }
    }

    fn setup_crypt(&mut self) -> Result<(), KernelError> {
        let Some(enc) = self.trailer.get(b"Encrypt").cloned() else { return Ok(()) };
        if let Obj::Ref(n, _) = enc {
            self.encrypt_num = Some(n);
        }
        let enc = self.resolve(&enc)?;
        let Some(enc_dict) = enc.as_dict() else { return Err(bad("the PDF's /Encrypt entry is not a dictionary")) };
        // The encryption dictionary's own entries may be references (CF is normally direct).
        let mut direct = Dict::default();
        for (k, v) in &enc_dict.0 {
            direct.0.push((k.clone(), (*self.resolve(v)?).clone()));
        }
        let id0 = match self.trailer.get(b"ID").cloned() {
            Some(id) => {
                let id = self.resolve(&id)?;
                id.as_array().and_then(|a| a.first()).and_then(Obj::as_str).map(<[u8]>::to_vec).unwrap_or_default()
            }
            None => Vec::new(),
        };
        self.crypt = Some(Crypt::new(&direct, &id0)?);
        // Objects read so far were not decrypted.
        self.cache.clear();
        self.object_streams.clear();
        Ok(())
    }

    // ---- The cross-reference table ---------------------------------------------------------------------------------

    fn load_table(&mut self, header: usize) -> Result<(), KernelError> {
        let data = self.data;
        let tail_from = data.len().saturating_sub(4096);
        let sx = rfind(&data[tail_from..], b"startxref").map(|p| p + tail_from).ok_or_else(|| bad("no startxref"))?;
        let mut lx = Lexer::new(data, sx + 9);
        let Tok::Num(off, true, _) = lx.next() else { return Err(bad("startxref has no offset")) };
        let mut at = off as usize;
        let mut seen = HashSet::new();
        let mut first = true;
        while seen.insert(at) && seen.len() <= MAX_XREF_SECTIONS {
            let trailer = match self.section(at).or_else(|e| if header > 0 { self.section(at + header) } else { Err(e) }) {
                Ok(t) => t,
                Err(e) if first => return Err(e),
                Err(_) => {
                    self.notes.push("an older cross-reference section of the PDF could not be read; it was skipped".into());
                    break;
                }
            };
            if first {
                self.trailer = trailer.clone();
                first = false;
            } else {
                for (k, v) in trailer.0.iter() {
                    if self.trailer.get(k).is_none() && k.as_slice() != b"Prev" && k.as_slice() != b"XRefStm" {
                        self.trailer.0.push((k.clone(), v.clone()));
                    }
                }
            }
            match trailer.get(b"Prev").and_then(Obj::as_i64) {
                Some(p) if p >= 0 => at = p as usize,
                _ => break,
            }
        }
        if self.entries.is_empty() {
            return Err(bad("empty cross-reference table"));
        }
        Ok(())
    }

    /// One cross-reference section (table or stream) at `at`; returns its trailer dictionary.
    fn section(&mut self, at: usize) -> Result<Dict, KernelError> {
        let data = self.data;
        let mut lx = Lexer::new(data, at);
        lx.skip_white();
        if data[lx.pos.min(data.len())..].starts_with(b"xref") {
            lx.pos += 4;
            loop {
                match lx.next() {
                    Tok::Num(start, true, _) => {
                        let Tok::Num(count, true, _) = lx.next() else { return Err(bad("bad xref subsection")) };
                        if start < 0.0 || count < 0.0 || count > data.len() as f64 / 10.0 {
                            return Err(bad("bad xref subsection"));
                        }
                        for i in 0..count as u64 {
                            let (Tok::Num(off, true, _), Tok::Num(_, true, _)) = (lx.next(), lx.next()) else {
                                return Err(bad("bad xref entry"));
                            };
                            let kind = lx.next();
                            let num = start as u64 + i;
                            if num > u32::MAX as u64 {
                                return Err(bad("bad xref entry"));
                            }
                            if kind == Tok::Kw(b"n") && off > 0.0 {
                                self.entries.entry(num as u32).or_insert(Entry::At(off as usize));
                            } else if kind != Tok::Kw(b"f") && kind != Tok::Kw(b"n") {
                                return Err(bad("bad xref entry"));
                            }
                        }
                    }
                    Tok::Kw(b"trailer") => {
                        let t = lexer::parse(&mut lx, true)?;
                        let Obj::Dict(t) = t else { return Err(bad("the trailer is not a dictionary")) };
                        // A hybrid file's cross-reference stream holds the compressed objects.
                        if let Some(x) = t.get(b"XRefStm").and_then(Obj::as_i64) {
                            if x >= 0 && self.stream_section(x as usize).is_err() {
                                self.notes.push("a hybrid PDF's cross-reference stream could not be read".into());
                            }
                        }
                        return Ok(t);
                    }
                    _ => return Err(bad("the xref table has no trailer")),
                }
            }
        }
        self.stream_section(at)
    }

    fn stream_section(&mut self, at: usize) -> Result<Dict, KernelError> {
        let (_, obj) = self.parse_at(at, None)?;
        let Obj::Stream(s) = obj else { return Err(bad("the cross-reference stream is not a stream")) };
        if s.dict.get(b"Type").and_then(Obj::as_name) != Some(b"XRef") {
            return Err(bad("not a cross-reference stream"));
        }
        let decoded = self.decode(&s)?;
        let w: Vec<usize> = s.dict.get(b"W").and_then(Obj::as_array).unwrap_or(&[]).iter()
            .map(|o| o.as_i64().unwrap_or(0).clamp(0, 8) as usize).collect();
        if w.len() < 3 {
            return Err(bad("the cross-reference stream has no /W"));
        }
        let size = s.dict.get(b"Size").and_then(Obj::as_i64).unwrap_or(0).max(0);
        let index: Vec<i64> = match s.dict.get(b"Index").and_then(Obj::as_array) {
            Some(a) => a.iter().filter_map(Obj::as_i64).collect(),
            None => vec![0, size],
        };
        let row = w[0] + w[1] + w[2];
        if row == 0 {
            return Err(bad("the cross-reference stream has an empty /W"));
        }
        let field = |r: &[u8], from: usize, n: usize| -> u64 { r[from..from + n].iter().fold(0u64, |a, &b| a << 8 | b as u64) };
        let mut rows = decoded.chunks_exact(row);
        for pair in index.chunks_exact(2) {
            let (start, count) = (pair[0].max(0) as u64, pair[1].max(0) as u64);
            for i in 0..count {
                let Some(r) = rows.next() else { break };
                let kind = if w[0] == 0 { 1 } else { field(r, 0, w[0]) };
                let a = field(r, w[0], w[1]);
                let b = field(r, w[0] + w[1], w[2]);
                let num = start + i;
                if num > u32::MAX as u64 {
                    break;
                }
                let entry = match kind {
                    1 if a > 0 && a < self.data.len() as u64 => Entry::At(a as usize),
                    2 if a <= u32::MAX as u64 && b <= u32::MAX as u64 => Entry::InStream(a as u32, b as u32),
                    _ => continue,
                };
                self.entries.entry(num as u32).or_insert(entry);
            }
        }
        Ok(s.dict.clone())
    }

    /// Scans the whole file for `N G obj` and trailers (later definitions replace earlier ones).
    fn rebuild(&mut self) -> Result<(), KernelError> {
        let data = self.data;
        let positions = scan_objects(data);
        let mut trailer: Option<Dict> = None;
        let mut from_xref_stream: Option<Dict> = None;
        let mut catalog: Option<u32> = None;
        for (&num, &at) in positions.iter() {
            self.entries.insert(num, Entry::At(at));
        }
        let mut order: Vec<(usize, u32)> = positions.iter().map(|(n, a)| (*a, *n)).collect();
        order.sort();
        for &(at, num) in &order {
            // Peek at the object's dictionary only (the stream data is not needed here).
            let mut lx = Lexer::new(data, at);
            let (_, _, _) = (lx.next(), lx.next(), lx.next());
            if lx.next() != Tok::DictOpen {
                continue;
            }
            lx.pos -= 2;
            let Ok(Obj::Dict(d)) = lexer::parse(&mut lx, true) else { continue };
            match d.get(b"Type").and_then(Obj::as_name) {
                Some(b"XRef") => {
                    if d.get(b"Root").is_some() {
                        from_xref_stream = Some(d);
                    }
                }
                Some(b"ObjStm") => self.scanned_streams.push(num),
                Some(b"Catalog") => catalog = Some(num),
                _ => {}
            }
        }
        let mut from = 0;
        while let Some(p) = find(data, b"trailer", from) {
            from = p + 7;
            let mut lx = Lexer::new(data, p + 7);
            if let Ok(Obj::Dict(d)) = lexer::parse(&mut lx, true) {
                if d.get(b"Root").is_some() {
                    trailer = Some(d);
                }
            }
        }
        let mut t = trailer.or(from_xref_stream).unwrap_or_default();
        let root_ok = match t.get(b"Root") {
            Some(Obj::Ref(n, _)) => positions.contains_key(n) || !self.scanned_streams.is_empty(),
            _ => false,
        };
        if !root_ok {
            if let Some(n) = catalog {
                t.0.retain(|(k, _)| k.as_slice() != b"Root");
                t.0.push((b"Root".to_vec(), Obj::Ref(n, 0)));
            }
        }
        self.scanned = Some(positions);
        self.trailer = t;
        Ok(())
    }

    // ---- Objects ---------------------------------------------------------------------------------------------------

    /// Follows references (a direct object is returned as it is).
    pub fn resolve(&mut self, o: &Obj) -> Result<Rc<Obj>, KernelError> {
        match o {
            Obj::Ref(n, _) => self.object(*n),
            other => Ok(Rc::new(other.clone())),
        }
    }

    /// A dictionary entry, resolved.
    pub fn get(&mut self, d: &Dict, key: &[u8]) -> Result<Option<Rc<Obj>>, KernelError> {
        match d.get(key) {
            Some(o) => {
                let r = self.resolve(o)?;
                Ok((!matches!(*r, Obj::Null)).then_some(r))
            }
            None => Ok(None),
        }
    }

    pub fn object(&mut self, num: u32) -> Result<Rc<Obj>, KernelError> {
        if let Some(o) = self.cache.get(&num) {
            return Ok(o.clone());
        }
        if self.depth >= MAX_RESOLVE_DEPTH {
            return Err(bad("the PDF's objects refer to each other too deeply"));
        }
        self.depth += 1;
        let r = self.load(num);
        self.depth -= 1;
        let o = Rc::new(r?);
        self.cache.insert(num, o.clone());
        Ok(o)
    }

    fn load(&mut self, num: u32) -> Result<Obj, KernelError> {
        let entry = match self.entries.get(&num) {
            Some(e) => Some(*e),
            None => self.search_streams(num)?,
        };
        let obj = match entry {
            None => return Ok(Obj::Null),
            Some(Entry::InStream(stm, index)) => return self.from_object_stream(num, stm, index),
            Some(Entry::At(at)) => match self.parse_at(at, Some(num)) {
                Ok((generation, o)) => (generation, o),
                Err(e) => {
                    // The table's offset is wrong: look where scanning finds the object.
                    let scanned = self.scanned.get_or_insert_with(|| scan_objects(self.data));
                    match scanned.get(&num).copied() {
                        Some(p) if p != at => self.parse_at(p, Some(num))?,
                        _ => return Err(e),
                    }
                }
            },
        };
        let (generation, mut o) = obj;
        if self.crypt.is_some() && self.encrypt_num != Some(num) {
            o = self.decrypt(o, num, generation);
        }
        Ok(o)
    }

    fn search_streams(&mut self, num: u32) -> Result<Option<Entry>, KernelError> {
        if self.scanned_streams.is_empty() {
            return Ok(None);
        }
        for stm in std::mem::take(&mut self.scanned_streams) {
            if let Ok(s) = self.object_stream(stm) {
                for (i, (n, _)) in s.1.iter().enumerate() {
                    self.entries.entry(*n).or_insert(Entry::InStream(stm, i as u32));
                }
            }
        }
        Ok(self.entries.get(&num).copied())
    }

    fn object_stream(&mut self, stm: u32) -> Result<Rc<(Vec<u8>, Vec<(u32, usize)>)>, KernelError> {
        if let Some(s) = self.object_streams.get(&stm) {
            return Ok(s.clone());
        }
        let o = self.object(stm)?;
        let Obj::Stream(s) = &*o else { return Err(bad(format!("object stream {stm} is not a stream"))) };
        let data = self.decode(s)?;
        let n = s.dict.get(b"N").and_then(Obj::as_i64).unwrap_or(0).clamp(0, 1 << 22) as usize;
        let first = s.dict.get(b"First").and_then(Obj::as_i64).unwrap_or(0).max(0) as usize;
        let mut lx = Lexer::new(&data, 0);
        let mut index = Vec::with_capacity(n.min(1 << 16));
        for _ in 0..n {
            let (Tok::Num(a, true, _), Tok::Num(b, true, _)) = (lx.next(), lx.next()) else { break };
            if a < 0.0 || b < 0.0 || a > u32::MAX as f64 {
                break;
            }
            index.push((a as u32, first + b as usize));
        }
        let parsed = Rc::new(((*data).clone(), index));
        self.object_streams.insert(stm, parsed.clone());
        Ok(parsed)
    }

    /// Objects inside an object stream are not encrypted on their own (the stream was).
    fn from_object_stream(&mut self, num: u32, stm: u32, index: u32) -> Result<Obj, KernelError> {
        if stm == num {
            return Ok(Obj::Null);
        }
        let s = self.object_stream(stm)?;
        let at = match s.1.get(index as usize) {
            Some((n, at)) if *n == num => *at,
            _ => match s.1.iter().find(|(n, _)| *n == num) {
                Some((_, at)) => *at,
                None => return Ok(Obj::Null),
            },
        };
        if at >= s.0.len() {
            return Ok(Obj::Null);
        }
        let mut lx = Lexer::new(&s.0, at);
        lexer::parse(&mut lx, true)
    }

    /// The indirect object `N G obj ... endobj` at `at` (checking the number when given), with its stream.
    fn parse_at(&mut self, at: usize, num: Option<u32>) -> Result<(u16, Obj), KernelError> {
        let data = self.data;
        if at >= data.len() {
            return Err(bad("an object offset is past the end of the file"));
        }
        let mut lx = Lexer::new(data, at);
        let (Tok::Num(n, true, _), Tok::Num(g, true, _), Tok::Kw(b"obj")) = (lx.next(), lx.next(), lx.next()) else {
            return Err(bad("no object at the offset the table gives"));
        };
        if num.is_some_and(|want| want as f64 != n) {
            return Err(bad("the object at the table's offset has another number"));
        }
        let this = n as u32;
        let generation = g.clamp(0.0, 65535.0) as u16;
        let obj = lexer::parse(&mut lx, true)?;
        if lx.next() != Tok::Kw(b"stream") {
            return Ok((generation, obj));
        }
        let Obj::Dict(dict) = obj else { return Err(bad("a stream without a dictionary")) };
        // The data starts after the end-of-line that follows "stream" (CR LF or LF; a lone CR is tolerated).
        let mut start = lx.pos;
        if data.get(start) == Some(&b'\r') {
            start += 1;
        }
        if data.get(start) == Some(&b'\n') {
            start += 1;
        }
        let length = match dict.get(b"Length") {
            Some(Obj::Ref(r, _)) if *r != this => self.object(*r).ok().and_then(|o| o.as_i64()),
            Some(o) => o.as_i64(),
            None => None,
        };
        let ends_well = |end: usize| -> bool {
            let mut l = Lexer::new(data, end);
            l.skip_white();
            data.get(l.pos..).is_some_and(|r| r.starts_with(b"endstream"))
        };
        let end = match length {
            Some(l) if l >= 0 && (l as u64) <= (data.len() - start) as u64 && ends_well(start + l as usize) => start + l as usize,
            _ => {
                let found = find(data, b"endstream", start);
                let end = match found {
                    Some(e) => e,
                    None => find(data, b"endobj", start).unwrap_or(data.len()),
                };
                let mut end = end;
                if end > start && data[end - 1] == b'\n' {
                    end -= 1;
                }
                if end > start && data[end - 1] == b'\r' {
                    end -= 1;
                }
                self.note("a stream's /Length in the PDF was missing or wrong; its end was found by searching for 'endstream'");
                end
            }
        };
        Ok((generation, Obj::Stream(Rc::new(Stream { dict, raw: data[start..end].to_vec() }))))
    }

    fn note(&mut self, s: &str) {
        if !self.notes.iter().any(|n| n == s) {
            self.notes.push(s.to_string());
        }
    }

    // ---- Encryption ------------------------------------------------------------------------------------------------

    fn decrypt(&self, o: Obj, num: u32, generation: u16) -> Obj {
        let Some(c) = &self.crypt else { return o };
        match o {
            Obj::Stream(s) => {
                let d = &s.dict;
                if d.get(b"Type").and_then(Obj::as_name) == Some(b"XRef") {
                    return Obj::Stream(s);
                }
                let dict = decrypt_strings(c, Obj::Dict(d.clone()), num, generation, 0);
                let Obj::Dict(dict) = dict else { return Obj::Stream(s) };
                let is_metadata = dict.get(b"Type").and_then(Obj::as_name) == Some(b"Metadata");
                // A Crypt filter first in the chain names the stream's own crypt filter.
                let first_filter = match dict.get(b"Filter") {
                    Some(Obj::Name(n)) => Some(n.as_slice()),
                    Some(Obj::Array(a)) => a.first().and_then(Obj::as_name),
                    _ => None,
                };
                let method = if first_filter == Some(b"Crypt") {
                    let parms = match dict.get(b"DecodeParms") {
                        Some(Obj::Dict(p)) => Some(p),
                        Some(Obj::Array(a)) => a.first().and_then(Obj::as_dict),
                        _ => None,
                    };
                    c.filter(parms.and_then(|p| p.get(b"Name")).and_then(Obj::as_name).unwrap_or(b"Identity"))
                } else if is_metadata && !c.encrypt_metadata {
                    Method::Identity
                } else {
                    c.stream
                };
                let raw = c.decrypt(method, num, generation, &s.raw);
                Obj::Stream(Rc::new(Stream { dict, raw }))
            }
            other => decrypt_strings(c, other, num, generation, 0),
        }
    }

    // ---- Decoding --------------------------------------------------------------------------------------------------

    /// The stream's data with all its filters applied.
    pub fn decode(&mut self, s: &Stream) -> Result<Rc<Vec<u8>>, KernelError> {
        let filters: Vec<Obj> = match self.get(&s.dict, b"Filter")? {
            None => Vec::new(),
            Some(f) => match &*f {
                Obj::Array(a) => a.clone(),
                other => vec![other.clone()],
            },
        };
        let parms: Vec<Obj> = match self.get(&s.dict, b"DecodeParms")? {
            None => Vec::new(),
            Some(p) => match &*p {
                Obj::Array(a) => a.clone(),
                other => vec![other.clone()],
            },
        };
        if filters.is_empty() {
            return Ok(Rc::new(s.raw.clone()));
        }
        let mut data = s.raw.clone();
        for (i, f) in filters.iter().enumerate() {
            let name = self.resolve(f)?;
            let Some(name) = name.as_name() else { return Err(bad("a stream filter is not a name")) };
            let p = match parms.get(i) {
                Some(p) => {
                    let p = self.resolve(p)?;
                    let mut direct = Dict::default();
                    if let Some(d) = p.as_dict() {
                        for (k, v) in &d.0 {
                            direct.0.push((k.clone(), (*self.resolve(v)?).clone()));
                        }
                    }
                    Some(direct)
                }
                None => None,
            };
            let mut notes = Vec::new();
            data = filters::apply(name, p.as_ref(), &data, self.budget, &mut notes)?;
            for n in notes {
                self.note(&n);
            }
        }
        self.budget = self.budget.saturating_sub(data.len());
        if self.budget == 0 {
            return Err(bad(format!("the PDF's streams expand past {} MB together; the file was refused", DECODE_BUDGET >> 20)));
        }
        Ok(Rc::new(data))
    }
}

fn decrypt_strings(c: &Crypt, o: Obj, num: u32, generation: u16, depth: usize) -> Obj {
    if depth > 64 {
        return o;
    }
    match o {
        Obj::Str(s) => Obj::Str(c.decrypt(c.string, num, generation, &s)),
        Obj::Array(a) => Obj::Array(a.into_iter().map(|x| decrypt_strings(c, x, num, generation, depth + 1)).collect()),
        Obj::Dict(d) => Obj::Dict(Dict(d.0.into_iter().map(|(k, v)| (k, decrypt_strings(c, v, num, generation, depth + 1))).collect())),
        other => other,
    }
}

/// Positions of `N G obj` throughout the file, the last definition of each number winning.
fn scan_objects(data: &[u8]) -> HashMap<u32, usize> {
    let mut out = HashMap::new();
    let mut from = 0;
    while let Some(p) = find(data, b"obj", from) {
        from = p + 3;
        if data.get(p + 3).is_some_and(|&c| lexer::is_regular(c)) {
            continue;
        }
        // Back over white space, the generation, white space, the number.
        let mut i = p;
        let digits_back = |i: &mut usize| -> usize {
            let end = *i;
            while *i > 0 && data[*i - 1].is_ascii_digit() {
                *i -= 1;
            }
            end - *i
        };
        let white_back = |i: &mut usize| -> usize {
            let end = *i;
            while *i > 0 && lexer::is_white(data[*i - 1]) {
                *i -= 1;
            }
            end - *i
        };
        if white_back(&mut i) == 0 || digits_back(&mut i) == 0 || white_back(&mut i) == 0 {
            continue;
        }
        let num_end = i;
        if digits_back(&mut i) == 0 || num_end - i > 10 {
            continue;
        }
        if i > 0 && lexer::is_regular(data[i - 1]) {
            continue;
        }
        let Ok(text) = std::str::from_utf8(&data[i..num_end]) else { continue };
        let Ok(n) = text.parse::<u32>() else { continue };
        out.insert(n, i);
    }
    out
}

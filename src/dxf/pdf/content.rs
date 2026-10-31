//! Content streams (ISO 32000-1 8.2-8.10, 14.6 marked content, 8.11 optional content): the graphics state stack and
//! current transformation matrix, path construction and painting, Form XObjects, inline images (whose data is
//! skipped exactly), text objects, shadings and images counted but not drawn, and optional-content groups as layers.

use std::collections::HashMap;
use std::rc::Rc;

use super::lexer::{self, Dict, Lexer, Obj, Tok};
use super::shapes::{self, Seg};
use super::xref::Doc;
use super::Affine;
use super::super::{fit, Curve, Drawing, P2};
use crate::KernelError;

/// Operators interpreted across the page and all its forms.
const MAX_OPERATORS: usize = 20_000_000;
/// Curves one page may produce.
const MAX_CURVES: usize = 5_000_000;
const MAX_FORM_DEPTH: usize = 16;
/// Form XObject uses one page may make (forms drawing forms many times over).
const MAX_FORM_USES: usize = 1_000_000;
const MAX_STACK: usize = 4096;
const MAX_OPERANDS: usize = 8192;

/// Optional content: each group's layer and whether it is on in the default configuration.
pub struct Oc {
    pub groups: HashMap<u32, (usize, bool)>,
}

/// The state a marked-content sequence (or a form) puts its content in.
#[derive(Clone, Copy)]
pub struct Mark {
    pub layer: usize,
    /// The layer that hides this content, when one does.
    pub hidden_by: Option<usize>,
}

pub struct Painter<'p, 'a> {
    pub doc: &'p mut Doc<'a>,
    pub d: &'p mut Drawing,
    pub oc: &'p Oc,
    group: usize,
    operators: usize,
    form_uses: usize,
    forms: Vec<u32>,
    form_data: HashMap<u32, Rc<Vec<u8>>>,
    hidden: Vec<usize>,
}

struct Sub {
    start: P2,
    segs: Vec<Seg>,
}

#[derive(Default)]
struct Path {
    subs: Vec<Sub>,
    /// The current point and the rounding of the numbers that placed it (a segment takes the coarser of its start's
    /// and its own operator's).
    cur: Option<(P2, f64)>,
    /// Whether the last subpath takes more segments (false after h or re).
    open: bool,
}

impl Path {
    fn move_to(&mut self, p: P2, q: f64) {
        self.subs.push(Sub { start: p, segs: Vec::new() });
        self.cur = Some((p, q));
        self.open = true;
    }

    /// The point a new segment starts from (starting a subpath after h, or at p for a segment with no current point).
    fn from(&mut self, p: P2, q: f64) -> (P2, f64) {
        let (c, cq) = self.cur.unwrap_or((p, q));
        if !self.open || self.subs.is_empty() {
            self.move_to(c, cq);
        }
        (c, cq)
    }

    fn line(&mut self, p: P2, q: f64) {
        let (c, cq) = self.from(p, q);
        if let Some(s) = self.subs.last_mut() {
            s.segs.push(Seg::Line(c, p, q.max(cq)));
        }
        self.cur = Some((p, q));
    }

    fn cubic(&mut self, p1: P2, p2: P2, p3: P2, q: f64) {
        let (c, cq) = self.from(p3, q);
        if let Some(s) = self.subs.last_mut() {
            s.segs.push(Seg::Cubic([c, p1, p2, p3], q.max(cq)));
        }
        self.cur = Some((p3, q));
    }

    fn close(&mut self) {
        if !self.open {
            return;
        }
        if let (Some(s), Some((c, q))) = (self.subs.last_mut(), self.cur) {
            if fit::dist(c, s.start) > 1e-12 * (1.0 + fit::len(c)) && !s.segs.is_empty() {
                s.segs.push(Seg::Line(c, s.start, q));
            }
            self.cur = Some((s.start, q));
        }
        self.open = false;
    }

    /// Filling closes every subpath of two or more segments.
    fn close_all(&mut self) {
        for s in &mut self.subs {
            if s.segs.len() < 2 {
                continue;
            }
            let (end, q) = match s.segs[s.segs.len() - 1] {
                Seg::Line(_, b, q) => (b, q),
                Seg::Cubic(p, q) => (p[3], q),
            };
            if fit::dist(end, s.start) > 1e-12 * (1.0 + fit::len(end)) {
                s.segs.push(Seg::Line(end, s.start, q));
            }
        }
    }
}

fn bad(m: impl Into<String>) -> KernelError {
    KernelError::InvalidInput(m.into())
}

impl<'p, 'a> Painter<'p, 'a> {
    pub fn new(doc: &'p mut Doc<'a>, d: &'p mut Drawing, oc: &'p Oc) -> Self {
        Painter { doc, d, oc, group: 0, operators: 0, form_uses: 0, forms: Vec::new(), form_data: HashMap::new(), hidden: Vec::new() }
    }

    /// Layers whose content was left out because they are off.
    pub fn hidden_layers(&self) -> &[usize] {
        &self.hidden
    }

    fn mark_hidden(&mut self, m: Mark) {
        if let Some(l) = m.hidden_by {
            if !self.hidden.contains(&l) {
                self.hidden.push(l);
            }
        }
    }

    fn skip(&mut self, kind: &str, m: Mark) {
        if m.hidden_by.is_none() {
            self.d.skip(kind);
        } else {
            self.mark_hidden(m);
        }
    }

    pub fn run(&mut self, content: &[u8], resources: &Dict, ctm: Affine, outer: Mark, depth: usize) -> Result<(), KernelError> {
        let mut lx = Lexer::new(content, 0);
        let mut operands: Vec<Obj> = Vec::new();
        let mut quanta: Vec<f64> = Vec::new();
        let mut nums: Vec<f64> = Vec::new();
        let mut ctm = ctm;
        let mut stack: Vec<Affine> = Vec::new();
        let mut overflow = 0usize;
        let mut path = Path::default();
        let mut marks: Vec<Mark> = Vec::new();
        loop {
            let tok = lx.next();
            let op = match tok {
                Tok::Eof => break,
                Tok::Kw(op) => op,
                Tok::Num(v, int, q) => {
                    operands.push(if int { Obj::Int(v as i64) } else { Obj::Real(v) });
                    quanta.push(q);
                    if operands.len() > MAX_OPERANDS {
                        operands.clear();
                        quanta.clear();
                    }
                    continue;
                }
                Tok::ArrClose | Tok::DictClose => continue,
                other => {
                    match lexer::build(&mut lx, other, false, 0) {
                        Ok(Some(o)) => {
                            operands.push(o);
                            quanta.push(0.0);
                        }
                        Ok(None) => {}
                        Err(_) => {
                            self.doc_note("a content stream in the PDF ended inside an array or dictionary; the rest of it was skipped");
                            break;
                        }
                    }
                    if operands.len() > MAX_OPERANDS {
                        operands.clear();
                        quanta.clear();
                    }
                    continue;
                }
            };
            self.operators += 1;
            if self.operators > MAX_OPERATORS {
                return Err(bad(format!("the PDF page has more than {MAX_OPERATORS} drawing operators; it was refused")));
            }
            let mark = marks.last().copied().unwrap_or(outer);
            nums.clear();
            nums.extend(operands.iter().filter_map(Obj::as_f64));
            // Writers drop trailing zeros, so the operand with the most decimals shows the precision they used.
            // Zeros ("0.00000000000000000") say nothing about it.
            let q_user = operands.iter().zip(&quanta).filter(|(o, q)| **q > 0.0 && o.as_f64().is_some_and(|v| v != 0.0))
                .map(|(_, q)| *q).fold(f64::INFINITY, f64::min);
            let q_user = if q_user.is_finite() { q_user } else { 0.0 };
            let q = q_user * ctm.scale();
            let pt = |i: usize| -> P2 { ctm.point([nums[i], nums[i + 1]]) };
            match op {
                b"q" => {
                    if stack.len() < MAX_STACK { stack.push(ctm) } else { overflow += 1 }
                }
                b"Q" => {
                    if overflow > 0 {
                        overflow -= 1;
                    } else if let Some(m) = stack.pop() {
                        ctm = m;
                    }
                }
                b"cm" if nums.len() >= 6 => {
                    let m = Affine([nums[0], nums[1], nums[2], nums[3], nums[4], nums[5]]);
                    if m.0.iter().all(|v| v.is_finite()) {
                        ctm = ctm.then(&m);
                    }
                }
                b"m" if nums.len() >= 2 => path.move_to(pt(0), q),
                b"l" if nums.len() >= 2 => path.line(pt(0), q),
                b"c" if nums.len() >= 6 => path.cubic(pt(0), pt(2), pt(4), q),
                b"v" if nums.len() >= 4 => {
                    let p0 = path.cur.map(|c| c.0).unwrap_or_else(|| pt(0));
                    path.cubic(p0, pt(0), pt(2), q)
                }
                b"y" if nums.len() >= 4 => {
                    let p3 = pt(2);
                    path.cubic(pt(0), p3, p3, q)
                }
                b"h" => path.close(),
                b"re" if nums.len() >= 4 => {
                    let (x, y, w, h) = (nums[0], nums[1], nums[2], nums[3]);
                    path.move_to(ctm.point([x, y]), q);
                    path.line(ctm.point([x + w, y]), q);
                    path.line(ctm.point([x + w, y + h]), q);
                    path.line(ctm.point([x, y + h]), q);
                    path.close();
                }
                b"S" | b"s" | b"f" | b"F" | b"f*" | b"B" | b"B*" | b"b" | b"b*" => {
                    if matches!(op, b"s" | b"b" | b"b*") {
                        path.close();
                    }
                    if op != b"S" && op != b"s" {
                        path.close_all();
                    }
                    let p = std::mem::take(&mut path);
                    self.paint(p, mark)?;
                }
                b"n" => path = Path::default(),
                b"BT" => self.skip("TEXT", mark),
                b"sh" => self.skip("SHADING", mark),
                b"BI" => {
                    skip_inline_image(&mut lx);
                    self.skip("IMAGE", mark);
                }
                b"Do" => {
                    if let Some(Obj::Name(name)) = operands.first() {
                        let name = name.clone();
                        self.xobject(&name, resources, ctm, mark, depth)?;
                    }
                }
                b"BMC" => marks.push(mark),
                b"BDC" => {
                    let m = match (operands.first(), operands.get(1)) {
                        (Some(Obj::Name(tag)), Some(props)) if tag.as_slice() == b"OC" => {
                            let props = props.clone();
                            self.oc_mark(&props, resources, mark)?
                        }
                        _ => mark,
                    };
                    marks.push(m);
                }
                b"EMC" => {
                    marks.pop();
                }
                _ => {}
            }
            operands.clear();
            quanta.clear();
        }
        Ok(())
    }

    fn doc_note(&mut self, s: &str) {
        if !self.doc.notes.iter().any(|n| n == s) {
            self.doc.notes.push(s.to_string());
        }
    }

    fn paint(&mut self, path: Path, mark: Mark) -> Result<(), KernelError> {
        if mark.hidden_by.is_some() {
            if path.subs.iter().any(|s| !s.segs.is_empty()) {
                self.mark_hidden(mark);
            }
            return Ok(());
        }
        let mut any = false;
        for s in &path.subs {
            for geometry in shapes::convert(&s.segs) {
                self.d.curves.push(Curve { geometry, layer: mark.layer, color: 256, group: self.group });
                any = true;
            }
        }
        if any {
            self.group += 1;
        }
        if self.d.curves.len() > MAX_CURVES {
            return Err(bad(format!("the PDF page draws more than {MAX_CURVES} pieces; it was refused")));
        }
        Ok(())
    }

    /// The layer and visibility of an /OC marked-content sequence or form.
    fn oc_mark(&mut self, props: &Obj, resources: &Dict, outer: Mark) -> Result<Mark, KernelError> {
        let (obj, num) = match props {
            Obj::Name(n) => {
                let Some(table) = self.doc.get(resources, b"Properties")? else { return Ok(outer) };
                let Some(table) = table.as_dict() else { return Ok(outer) };
                match table.get(n) {
                    Some(Obj::Ref(r, _)) => (self.doc.object(*r)?, Some(*r)),
                    Some(o) => (Rc::new(o.clone()), None),
                    None => return Ok(outer),
                }
            }
            Obj::Ref(r, _) => (self.doc.object(*r)?, Some(*r)),
            other => (Rc::new(other.clone()), None),
        };
        self.oc_of(&obj, num, outer)
    }

    fn oc_of(&mut self, obj: &Obj, num: Option<u32>, outer: Mark) -> Result<Mark, KernelError> {
        let Some(d) = obj.as_dict() else { return Ok(outer) };
        let group_state = |painter: &mut Self, n: Option<u32>, d: &Dict| -> (usize, bool) {
            if let Some(&(layer, on)) = n.and_then(|n| painter.oc.groups.get(&n)) {
                return (layer, on);
            }
            let name = d.get(b"Name").and_then(Obj::as_str).map(lexer::text).unwrap_or_default();
            (if name.trim().is_empty() { outer.layer } else { painter.d.layer_index(&name) }, true)
        };
        let (layer, on, off_layer) = if d.get(b"Type").and_then(Obj::as_name) == Some(b"OCMD") {
            let members: Vec<Obj> = match d.get(b"OCGs") {
                Some(Obj::Array(a)) => a.clone(),
                Some(o) => vec![o.clone()],
                None => Vec::new(),
            };
            let mut states = Vec::new();
            for m in members.iter().take(256) {
                let n = if let Obj::Ref(r, _) = m { Some(*r) } else { None };
                let g = self.doc.resolve(m)?;
                if let Some(gd) = g.as_dict() {
                    let gd = gd.clone();
                    states.push(group_state(self, n, &gd));
                }
            }
            if states.is_empty() {
                return Ok(outer);
            }
            let policy = d.get(b"P").and_then(Obj::as_name).unwrap_or(b"AnyOn");
            let on = match policy {
                b"AllOn" => states.iter().all(|s| s.1),
                b"AnyOff" => states.iter().any(|s| !s.1),
                b"AllOff" => states.iter().all(|s| !s.1),
                _ => states.iter().any(|s| s.1),
            };
            let off = states.iter().find(|s| !s.1).map_or(states[0].0, |s| s.0);
            (states[0].0, on, off)
        } else {
            let (layer, on) = group_state(self, num, d);
            (layer, on, layer)
        };
        Ok(Mark {
            layer,
            hidden_by: outer.hidden_by.or((!on).then_some(off_layer)),
        })
    }

    fn xobject(&mut self, name: &[u8], resources: &Dict, ctm: Affine, mark: Mark, depth: usize) -> Result<(), KernelError> {
        let Some(table) = self.doc.get(resources, b"XObject")? else { return Ok(()) };
        let Some(table) = table.as_dict() else { return Ok(()) };
        let (obj, num) = match table.get(name) {
            Some(Obj::Ref(r, _)) => (self.doc.object(*r)?, Some(*r)),
            Some(o) => (Rc::new(o.clone()), None),
            None => {
                self.doc_note("a drawing refers to a picture or form the PDF does not have; it was left out");
                return Ok(());
            }
        };
        let Obj::Stream(s) = &*obj else { return Ok(()) };
        match s.dict.get(b"Subtype").and_then(Obj::as_name) {
            Some(b"Image") => {
                self.skip("IMAGE", mark);
                Ok(())
            }
            Some(b"Form") => {
                let mut mark = mark;
                if let Some(oc) = s.dict.get(b"OC") {
                    let oc = oc.clone();
                    let n = if let Obj::Ref(r, _) = oc { Some(r) } else { None };
                    let o = self.doc.resolve(&oc)?;
                    mark = self.oc_of(&o, n, mark)?;
                }
                if mark.hidden_by.is_some() {
                    self.mark_hidden(mark);
                    return Ok(());
                }
                self.form_uses += 1;
                if self.form_uses > MAX_FORM_USES {
                    return Err(bad(format!("the PDF page uses forms more than {MAX_FORM_USES} times; it was refused")));
                }
                if depth >= MAX_FORM_DEPTH || num.is_some_and(|n| self.forms.contains(&n)) {
                    self.doc_note("a form in the PDF uses itself or is nested too deeply; the repeated use was left out");
                    return Ok(());
                }
                let matrix = match self.doc.get(&s.dict, b"Matrix")? {
                    Some(m) => {
                        let v: Vec<f64> = m.as_array().unwrap_or(&[]).iter().filter_map(Obj::as_f64).collect();
                        if v.len() == 6 { Affine([v[0], v[1], v[2], v[3], v[4], v[5]]) } else { Affine::IDENTITY }
                    }
                    None => Affine::IDENTITY,
                };
                let own = self.doc.get(&s.dict, b"Resources")?;
                let res: Dict = match own.as_deref().and_then(Obj::as_dict) {
                    Some(r) => r.clone(),
                    None => resources.clone(),
                };
                let data = match num.and_then(|n| self.form_data.get(&n).cloned()) {
                    Some(d) => d,
                    None => {
                        let d = self.doc.decode(s)?;
                        if let Some(n) = num {
                            self.form_data.insert(n, d.clone());
                        }
                        d
                    }
                };
                if let Some(n) = num {
                    self.forms.push(n);
                }
                let r = self.run(&data, &res, ctm.then(&matrix), mark, depth + 1);
                if num.is_some() {
                    self.forms.pop();
                }
                r
            }
            Some(b"PS") => {
                self.skip("POSTSCRIPT", mark);
                Ok(())
            }
            _ => Ok(()),
        }
    }
}

/// Skips an inline image (8.9.7) from after BI to after EI: the dictionary, then the data, whose length is known
/// when it is unfiltered (or given by /L), else found as white space, "EI", white space or the end, followed by
/// text that looks like operators.
fn skip_inline_image(lx: &mut Lexer) {
    let mut dict = Dict::default();
    loop {
        match lx.next() {
            Tok::Kw(b"ID") => break,
            Tok::Eof => return,
            Tok::Name(k) => {
                let t = lx.next();
                match lexer::build(lx, t, false, 0) {
                    Ok(Some(v)) => dict.0.push((k, v)),
                    Ok(None) => {}
                    Err(_) => return,
                }
            }
            _ => {}
        }
    }
    let b = lx.b;
    let mut start = lx.pos;
    if b.get(start).is_some_and(|&c| lexer::is_white(c)) {
        start += 1;
    }
    let get = |a: &[u8], l: &[u8]| dict.get(a).or_else(|| dict.get(l));
    let filtered = get(b"F", b"Filter").is_some_and(|f| !matches!(f, Obj::Array(a) if a.is_empty()));
    let ends_at = |end: usize| -> Option<usize> {
        let mut i = end;
        while i < b.len() && lexer::is_white(b[i]) {
            i += 1;
        }
        (b.get(i..i + 2) == Some(b"EI") && b.get(i + 2).is_none_or(|&c| !lexer::is_regular(c))).then_some(i + 2)
    };
    let known = get(b"L", b"Length").and_then(Obj::as_i64).or_else(|| {
        if filtered {
            return None;
        }
        let w = get(b"W", b"Width")?.as_i64()?;
        let h = get(b"H", b"Height")?.as_i64()?;
        let mask = matches!(get(b"IM", b"ImageMask"), Some(Obj::Bool(true)));
        let bpc = if mask { 1 } else { get(b"BPC", b"BitsPerComponent").and_then(Obj::as_i64).unwrap_or(8) };
        let comps = if mask { 1 } else {
            match get(b"CS", b"ColorSpace") {
                Some(Obj::Name(n)) => match n.as_slice() {
                    b"G" | b"DeviceGray" | b"CalGray" | b"I" | b"Indexed" => 1,
                    b"RGB" | b"DeviceRGB" | b"CalRGB" => 3,
                    b"CMYK" | b"DeviceCMYK" => 4,
                    _ => return None,
                },
                Some(Obj::Array(a)) if matches!(a.first().and_then(Obj::as_name), Some(b"I" | b"Indexed")) => 1,
                _ => return None,
            }
        };
        if !(0..=1 << 20).contains(&w) || !(0..=1 << 20).contains(&h) || !(1..=16).contains(&bpc) {
            return None;
        }
        Some(h * ((w * comps * bpc + 7) / 8))
    });
    if let Some(n) = known {
        if n >= 0 {
            if let Some(after) = start.checked_add(n as usize).and_then(ends_at) {
                lx.pos = after;
                return;
            }
        }
    }
    // Search: "EI" between white space (or the end), with plain text after it.
    let mut i = start;
    while i + 1 < b.len() {
        if &b[i..i + 2] == b"EI" && (i == 0 || lexer::is_white(b[i - 1]))
            && b.get(i + 2).is_none_or(|&c| !lexer::is_regular(c))
        {
            let after = &b[i + 2..b.len().min(i + 66)];
            if after.iter().all(|&c| lexer::is_white(c) || (0x20..0x7F).contains(&c)) {
                lx.pos = i + 2;
                return;
            }
        }
        i += 1;
    }
    lx.pos = b.len();
}

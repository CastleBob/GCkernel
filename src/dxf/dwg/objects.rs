//! DWG objects (ODA specification chapter 20): the common object and entity data, the classes, layers and block
//! headers, text styles, and every curve and text entity decoded into the DXF reader's entity records (DXF group
//! codes).

use std::collections::HashMap;

use super::bits::{bad, Bits, Ver, R};
use super::Source;
use crate::dxf::codes::{unescape, Pair, Value};
use crate::dxf::read::{Block, Entity};
use crate::dxf::{Drawing, Layer, TextStyle};

/// Fixed object type numbers (spec 20.3).
mod t {
    pub const TEXT: u16 = 1;
    pub const ATTRIB: u16 = 2;
    pub const BLOCK: u16 = 4;
    pub const ENDBLK: u16 = 5;
    pub const SEQEND: u16 = 6;
    pub const INSERT: u16 = 7;
    pub const MINSERT: u16 = 8;
    pub const VERTEX_2D: u16 = 0x0a;
    pub const VERTEX_PFACE_FACE: u16 = 0x0e;
    pub const POLYLINE_2D: u16 = 0x0f;
    pub const POLYLINE_3D: u16 = 0x10;
    pub const VERTEX_3D: u16 = 0x0b;
    pub const ARC: u16 = 0x11;
    pub const CIRCLE: u16 = 0x12;
    pub const LINE: u16 = 0x13;
    pub const POLYLINE_PFACE: u16 = 0x1d;
    pub const POLYLINE_MESH: u16 = 0x1e;
    pub const ELLIPSE: u16 = 0x23;
    pub const SPLINE: u16 = 0x24;
    pub const BLOCK_HEADER: u16 = 0x31;
    pub const LAYER: u16 = 0x33;
    pub const SHAPEFILE: u16 = 0x35;
    pub const MTEXT: u16 = 0x2c;
    pub const LWPOLYLINE: u16 = 0x4d;
}

/// Fixed type numbers that are entities, with their DXF names.
fn fixed_entity(kind: u16) -> Option<&'static str> {
    Some(match kind {
        1 => "TEXT",
        2 => "ATTRIB",
        3 => "ATTDEF",
        4 => "BLOCK",
        5 => "ENDBLK",
        6 => "SEQEND",
        7 | 8 => "INSERT",
        0x0a..=0x0e => "VERTEX",
        0x0f | 0x10 | 0x1d | 0x1e => "POLYLINE",
        0x11 => "ARC",
        0x12 => "CIRCLE",
        0x13 => "LINE",
        0x14..=0x1a => "DIMENSION",
        0x1b => "POINT",
        0x1c => "3DFACE",
        0x1f => "SOLID",
        0x20 => "TRACE",
        0x21 => "SHAPE",
        0x22 => "VIEWPORT",
        0x23 => "ELLIPSE",
        0x24 => "SPLINE",
        0x25 => "REGION",
        0x26 => "3DSOLID",
        0x27 => "BODY",
        0x28 => "RAY",
        0x29 => "XLINE",
        0x2b => "OLEFRAME",
        0x2c => "MTEXT",
        0x2d => "LEADER",
        0x2e => "TOLERANCE",
        0x2f => "MLINE",
        0x4a => "OLE2FRAME",
        0x4d => "LWPOLYLINE",
        0x4e => "HATCH",
        0x1f2 => "ACAD_PROXY_ENTITY",
        _ => return None,
    })
}

struct Class {
    name: String,
    entity: bool,
}

/// One decoded object: its type, handle and the three streams positioned after the common data.
pub(super) struct Obj<'a> {
    pub kind: u16,
    pub handle: u64,
    pub data: Bits<'a>,
    pub refs: Bits<'a>,
    pub strings: Option<Bits<'a>>,
    pub entity: bool,
    // Entity common data.
    pub entmode: u8,
    pub reactors: usize,
    pub xdic: bool,
    pub nolinks: bool,
    pub bylayer_lt: bool,
    pub color: i32,
    pub color_book: bool,
    pub lt_flags: u8,
    pub ps_flags: u8,
    pub mat_flags: u8,
    pub styles: usize,
    pub invisible: bool,
    /// Extended data: (first bit, byte count) of each application's items.
    pub eed: Vec<(usize, usize)>,
}

struct LayerInfo {
    name: String,
}

struct BlockInfo {
    name: String,
    /// The upper-cased name, made unique with the handle when names repeat (anonymous blocks such as *U, *I).
    key: String,
    base: [f64; 3],
    xref: bool,
    /// R2004+: the entities in order; earlier versions find them by owner.
    owned: Option<Vec<u64>>,
}

pub(super) struct Reader<'a> {
    src: &'a Source,
    classes: HashMap<u16, Class>,
    layers: HashMap<u64, LayerInfo>,
    blocks: HashMap<u64, BlockInfo>,
    /// Text style handle -> name.
    styles: HashMap<u64, String>,
}

/// Counts read from a damaged file are refused past this.
const MAX_ITEMS: usize = 10_000_000;

impl<'a> Reader<'a> {
    pub fn new(src: &'a Source) -> R<Reader<'a>> {
        let mut r = Reader {
            src, classes: HashMap::new(), layers: HashMap::new(), blocks: HashMap::new(), styles: HashMap::new(),
        };
        r.read_classes();
        Ok(r)
    }

    fn ver(&self) -> Ver {
        self.src.ver
    }

    /// Class definitions (spec chapter 10 and 5.8). A damaged class list only loses names for custom types.
    fn read_classes(&mut self) {
        let data = &self.src.classes;
        let ver = self.ver();
        let mut b = Bits::new(data, ver);
        let mut strings: Option<Bits> = None;
        if ver >= Ver::R2007 {
            // RL size in bits of the class data, from the start of the data; the string stream ends there.
            let Ok(bits) = b.rl() else { return };
            let end = bits as usize;
            strings = string_stream(&b, end);
            b.end = strings.as_ref().map_or(end, |s| s.pos).min(b.end);
        }
        if ver >= Ver::R2004 {
            if b.bs().is_err() || b.rc().is_err() || b.rc().is_err() || b.b().is_err() {
                return;
            }
        }
        while b.left() >= 32 {
            let mut one = || -> R<(u16, Class)> {
                let number = b.bs()? as u16;
                b.bs()?;
                b.tv(&mut strings)?;
                b.tv(&mut strings)?;
                let name = b.tv(&mut strings)?;
                b.b()?;
                let item = b.bs()? as u16;
                if ver >= Ver::R2004 {
                    // Object count, DWG version, maintenance version and two unknowns, all bit longs.
                    for _ in 0..5 {
                        b.bl()?;
                    }
                }
                Ok((number, Class { name: name.trim().to_ascii_uppercase(), entity: item == 0x1f2 }))
            };
            match one() {
                Ok((number, class)) => {
                    self.classes.insert(number, class);
                }
                Err(_) => break,
            }
        }
    }

    fn is_entity(&self, kind: u16) -> bool {
        if kind >= 500 {
            return self.classes.get(&kind).is_some_and(|c| c.entity);
        }
        fixed_entity(kind).is_some()
    }

    fn type_name(&self, kind: u16) -> String {
        if let Some(c) = self.classes.get(&kind) {
            if kind >= 500 {
                return c.name.clone();
            }
        }
        fixed_entity(kind).map(str::to_string).unwrap_or_else(|| format!("DWG type {kind}"))
    }

    /// The object with this handle, common data read (spec 20.1, 20.2, 20.4.1).
    pub fn object(&self, handle: u64) -> R<Obj<'a>> {
        let Some(&offset) = self.src.map.get(&handle) else {
            return bad(format!("object {handle:X} is not in the object map"));
        };
        let buf: &'a [u8] = &self.src.objects;
        let ver = self.ver();
        let mut at = offset;
        let size = super::bits::ms(buf, &mut at)?;
        let handle_bits = if ver >= Ver::R2010 { super::bits::mc(buf, &mut at, false)? as usize } else { 0 };
        let Some(slice) = buf.get(at..at.saturating_add(size)).filter(|s| s.len() == size) else {
            return bad(format!("object {handle:X} runs past the end of its section"));
        };
        let mut b = Bits::new(slice, ver);
        let kind = b.ot()?;
        let mut data_end = if ver >= Ver::R2010 {
            if handle_bits > size * 8 {
                return bad(format!("object {handle:X}: handle stream larger than the object"));
            }
            Some(size * 8 - handle_bits)
        } else if ver >= Ver::R2000 {
            Some(b.rl()? as usize)
        } else {
            None
        };
        let own = b.h(0)?.value;
        // Extended entity data: (size, application handle, bytes) until a zero size.
        let mut eed = Vec::new();
        for _ in 0..100_000 {
            let n = b.bs()? as u16 as usize;
            if n == 0 {
                break;
            }
            b.h(0)?;
            eed.push((b.pos, n));
            b.skip(n * 8)?;
        }
        let entity = self.is_entity(kind);
        if entity && b.b()? {
            let n = if ver >= Ver::R2010 { b.bll()? as usize } else { b.rl()? as usize };
            b.skip(n.saturating_mul(8))?;
        }
        if ver <= Ver::R14 {
            data_end = Some(b.rl()? as usize);
        }
        let data_end = data_end.unwrap_or(0);
        if data_end > size * 8 || data_end < b.pos {
            return bad(format!("object {handle:X}: its data size does not fit"));
        }
        let strings = if ver >= Ver::R2007 { string_stream(&b, data_end) } else { None };
        let refs = b.window(data_end, size * 8);
        // R2007+: the last data bit flags the string stream; the data ends where the strings (or that bit) start.
        b.end = match &strings {
            Some(s) => s.pos.min(data_end),
            None if ver >= Ver::R2007 => data_end.saturating_sub(1).max(b.pos),
            None => data_end,
        };
        let mut o = Obj {
            kind, handle: own, data: b, refs, strings, entity, entmode: 0, reactors: 0, xdic: true, nolinks: true,
            bylayer_lt: true, color: 256, color_book: false, lt_flags: 0, ps_flags: 0, mat_flags: 0, styles: 0,
            invisible: false, eed,
        };
        let b = &mut o.data;
        if entity {
            o.entmode = b.bb()?;
            o.reactors = b.bl()? as i64 as usize;
            if ver >= Ver::R2004 {
                o.xdic = !b.b()?;
            }
            if ver >= Ver::R2013 {
                b.b()?;
            }
            if ver <= Ver::R14 {
                o.bylayer_lt = b.b()?;
            }
            if ver < Ver::R2004 {
                o.nolinks = b.b()?;
            }
            if ver < Ver::R2004 {
                o.color = b.bs()? as i32;
            } else {
                let n = b.bs()? as u16;
                o.color = (n & 0x1ff) as i32;
                if n & 0x8000 != 0 {
                    let rgb = b.bl()? as u32;
                    let index = super::bits::method_index(rgb);
                    o.color = if (1..=255).contains(&index) { index } else { 256 };
                }
                if n & 0x4000 != 0 {
                    o.color_book = true;
                }
                if n & 0x2000 != 0 {
                    b.bl()?;
                }
            }
            b.bd()?;
            if ver >= Ver::R2000 {
                o.lt_flags = b.bb()?;
                o.ps_flags = b.bb()?;
            }
            if ver >= Ver::R2007 {
                o.mat_flags = b.bb()?;
                b.rc()?;
            }
            if ver >= Ver::R2010 {
                o.styles = usize::from(b.b()?) + usize::from(b.b()?) + usize::from(b.b()?);
            }
            o.invisible = b.bs()? & 1 != 0;
            if ver >= Ver::R2000 {
                b.rc()?;
            }
        } else {
            o.reactors = b.bl()? as i64 as usize;
            if ver >= Ver::R2004 {
                o.xdic = !b.b()?;
            }
            if ver >= Ver::R2013 {
                b.b()?;
            }
        }
        Ok(o)
    }

    /// The common entity handles (spec 20.4.2): gives (owner when stored, layer); leaves `refs` at the entity's own.
    fn entity_refs(&self, o: &mut Obj) -> R<(Option<u64>, u64)> {
        let ver = self.ver();
        let h = o.handle;
        let r = &mut o.refs;
        let owner = if o.entmode == 0 { Some(r.h(h)?.value) } else { None };
        r.skip_h(o.reactors)?;
        if o.xdic {
            r.h(h)?;
        }
        let mut layer = 0;
        if ver <= Ver::R14 {
            layer = r.h(h)?.value;
            if !o.bylayer_lt {
                r.h(h)?;
            }
        }
        if ver <= Ver::R2000 && !o.nolinks {
            r.h(h)?;
            r.h(h)?;
        }
        if ver >= Ver::R2004 && o.color_book {
            r.h(h)?;
        }
        if ver >= Ver::R2000 {
            layer = r.h(h)?.value;
            if o.lt_flags == 3 {
                r.h(h)?;
            }
            if ver >= Ver::R2007 && o.mat_flags == 3 {
                r.h(h)?;
            }
            if o.ps_flags == 3 {
                r.h(h)?;
            }
            r.skip_h(o.styles)?;
        }
        Ok((owner, layer))
    }

    /// Non-entity common handles: owner, reactors, extension dictionary.
    fn object_refs(&self, o: &mut Obj) -> R<u64> {
        let h = o.handle;
        let owner = o.refs.h(h)?.value;
        o.refs.skip_h(o.reactors)?;
        if o.xdic {
            o.refs.h(h)?;
        }
        Ok(owner)
    }

    /// Every LAYER object (spec 20.4.54) into the drawing's layer table.
    pub fn layers(&mut self, d: &mut Drawing) -> R<()> {
        let handles: Vec<u64> = self.src.order.clone();
        for handle in handles {
            let Ok(mut o) = self.object(handle) else { continue };
            if o.kind == t::SHAPEFILE {
                if let Err(e) = self.style(&mut o, handle, d) {
                    d.warn_text(format!("a text style could not be decoded ({e}); its texts name no font"));
                }
                continue;
            }
            if o.kind != t::LAYER {
                continue;
            }
            let ver = self.ver();
            let b = &mut o.data;
            let name = b.table_entry(&mut o.strings)?;
            // The spec names the second flag (bit 2 of the R2000 values) "on", but visible layers have it clear, so it
            // is read as "off" (the twin DXFs agree).
            let (frozen, off, locked);
            if ver <= Ver::R14 {
                frozen = b.b()?;
                off = b.b()?;
                b.b()?;
                locked = b.b()?;
            } else {
                let values = b.bs()?;
                frozen = values & 1 != 0;
                off = values & 2 != 0;
                locked = values & 8 != 0;
            }
            let color = b.cmc(&mut o.strings)?;
            let name = name.trim().to_string();
            if name.is_empty() {
                continue;
            }
            let index = d.layer_index(&name);
            d.layers[index] = Layer {
                name: d.layers[index].name.clone(),
                color: if color.abs() >= 1 && color.abs() <= 255 { color.abs() } else { 7 },
                frozen,
                off: off || color < 0,
                locked,
            };
            self.layers.insert(handle, LayerInfo { name });
        }
        Ok(())
    }

    /// A text style (spec 20.4.56, SHAPEFILE): name, flags, font file; the TrueType family from the extended data.
    fn style(&mut self, o: &mut Obj<'a>, handle: u64, d: &mut Drawing) -> R<()> {
        let ver = self.ver();
        let b = &mut o.data;
        let name = b.table_entry(&mut o.strings)?;
        b.b()?; // vertical
        let shape = b.b()?;
        b.bd()?; // fixed height
        b.bd()?; // width factor
        b.bd()?; // oblique angle
        b.rc()?; // generation
        b.bd()?; // last height
        let file = b.tv(&mut o.strings)?;
        let name = name.trim().to_string();
        self.styles.insert(handle, name.clone());
        if name.is_empty() || shape {
            return Ok(());
        }
        let family = o.eed.iter().find_map(|&(start, n)| eed_string(&o.data.window(start, start + n * 8), ver))
            .unwrap_or_default();
        d.styles.push(TextStyle { name, file: file.trim().to_string(), family: family.trim().to_string() });
        Ok(())
    }

    fn style_name(&self, handle: u64) -> String {
        self.styles.get(&handle).cloned().unwrap_or_else(|| "Standard".to_string())
    }

    /// An insert's ATTRIB (spec 20.4.4) as a DXF-style child record; refused ones carry their reason as the kind.
    fn attrib(&self, handle: u64) -> Entity {
        match self.attrib_record(handle) {
            Ok(Some(e)) => e,
            Ok(None) => Entity { kind: "SEQEND".into(), pairs: Vec::new(), children: Vec::new() },
            Err(e) => refused_attrib(format!("an attribute could not be decoded ({e})")),
        }
    }

    fn attrib_record(&self, handle: u64) -> R<Option<Entity>> {
        let mut o = self.object(handle)?;
        if o.kind != t::ATTRIB {
            return Ok(None);
        }
        let (_, layer) = self.entity_refs(&mut o)?;
        let ver = self.ver();
        let mut e = Rec { kind: "ATTRIB".into(), pairs: Vec::new(), children: Vec::new() };
        e.text(8, &self.layer_name(layer));
        e.int(62, o.color as i64);
        if o.invisible {
            e.int(60, 1);
        }
        let h = o.handle;
        let b = &mut o.data;
        text_data(b, ver, &mut o.strings, &mut e, 74)?;
        if ver >= Ver::R2010 {
            b.rc()?; // version
        }
        if ver >= Ver::R2018 {
            let kind = b.rc()?;
            if kind == 2 || kind == 4 {
                return Ok(Some(refused_attrib("a multi-line attribute (DWG 2018), which is not read".into())));
            }
        }
        e.text(2, &b.tv(&mut o.strings)?);
        b.bs()?; // field length
        e.int(70, b.rc()? as i64);
        let style = o.refs.h(h)?.value;
        e.text(7, &self.style_name(style));
        Ok(Some(e.done()))
    }

    /// Block headers (spec 20.4.52).
    fn block_headers(&mut self) -> R<()> {
        let ver = self.ver();
        for &handle in &self.src.order {
            let Ok(mut o) = self.object(handle) else { continue };
            if o.kind != t::BLOCK_HEADER {
                continue;
            }
            let b = &mut o.data;
            let name = b.table_entry(&mut o.strings)?;
            b.b()?;
            b.b()?;
            let is_xref = b.b()?;
            let overlaid = b.b()?;
            if ver >= Ver::R2000 {
                b.b()?;
            }
            let owned_count = if ver >= Ver::R2004 { b.bl()?.max(0) as usize } else { 0 };
            let base = b.p3()?;
            let xref = is_xref || overlaid;
            let mut owned = None;
            if ver >= Ver::R2004 && !xref {
                self.object_refs(&mut o)?;
                let r = &mut o.refs;
                r.h(handle)?; // NULL
                r.h(handle)?; // BLOCK entity
                if owned_count > r.left() / 8 || owned_count > MAX_ITEMS {
                    return bad(format!("block '{name}' owns {owned_count} entities, more than its data holds"));
                }
                let mut list = Vec::with_capacity(owned_count);
                for _ in 0..owned_count {
                    list.push(r.h(handle)?.value);
                }
                owned = Some(list);
            }
            self.blocks.insert(handle, BlockInfo { name: name.trim().to_string(), key: String::new(), base, xref, owned });
        }
        let mut seen: HashMap<String, usize> = HashMap::new();
        for b in self.blocks.values() {
            *seen.entry(b.name.to_ascii_uppercase()).or_default() += 1;
        }
        for (handle, b) in self.blocks.iter_mut() {
            let upper = b.name.to_ascii_uppercase();
            b.key = if seen[&upper] > 1 { format!("{upper}|{handle:X}") } else { upper };
        }
        Ok(())
    }

    fn layer_name(&self, handle: u64) -> String {
        self.layers.get(&handle).map_or_else(|| "0".to_string(), |l| l.name.clone())
    }

    /// Model-space entities (paper space ones flagged 67 = 1) and every block definition.
    pub fn entities(&mut self, d: &mut Drawing) -> R<(HashMap<String, Block>, Vec<Entity>)> {
        self.block_headers()?;
        let ver = self.ver();
        // Before R2004, entities are found by their owner: model and paper space by entity mode, block contents by
        // the owner handle, in object map order.
        let mut by_owner: HashMap<u64, Vec<u64>> = HashMap::new();
        let mut model_list = Vec::new();
        let mut paper_list = Vec::new();
        if ver < Ver::R2004 {
            for &handle in &self.src.order {
                let Ok(mut o) = self.object(handle) else { continue };
                if !o.entity {
                    continue;
                }
                let Ok((owner, _)) = self.entity_refs(&mut o) else { continue };
                match (o.entmode, owner) {
                    (2, _) => model_list.push(handle),
                    (1, _) => paper_list.push(handle),
                    (0, Some(owner)) => by_owner.entry(owner).or_default().push(handle),
                    _ => {}
                }
            }
        }
        let names: Vec<(u64, String)> = self.blocks.iter().map(|(h, b)| (*h, b.name.clone())).collect();
        let mut model = Vec::new();
        let mut blocks = HashMap::new();
        let mut failures = 0usize;
        let mut first_failure = String::new();
        for (handle, name) in names {
            let upper = name.to_ascii_uppercase();
            let info = &self.blocks[&handle];
            let list: Vec<u64> = match &info.owned {
                Some(list) => list.clone(),
                None if upper == "*MODEL_SPACE" && ver < Ver::R2004 => model_list.clone(),
                None if upper == "*PAPER_SPACE" && ver < Ver::R2004 => paper_list.clone(),
                None => by_owner.get(&handle).cloned().unwrap_or_default(),
            };
            let (base, xref) = (info.base, info.xref);
            let paper = upper == "*PAPER_SPACE";
            let mut records = Vec::with_capacity(list.len());
            for h in list {
                match self.record(h, paper, &by_owner) {
                    Ok(Some(e)) => records.push(e),
                    Ok(None) => {}
                    Err(e) => {
                        failures += 1;
                        if first_failure.is_empty() {
                            first_failure = e.to_string();
                        }
                    }
                }
            }
            if upper == "*MODEL_SPACE" {
                // Model space first, so the groups follow the drawing order.
                records.append(&mut model);
                model = records;
            } else if paper {
                model.extend(records);
            } else {
                blocks.insert(self.blocks[&handle].key.clone(), Block { base, entities: records, xref });
            }
        }
        if failures > 0 {
            d.warn(format!("{failures} entities could not be decoded and were left out (the first: {first_failure})"));
        }
        Ok((blocks, model))
    }

    /// One entity as a DXF-style record; None for objects that are not drawn on their own.
    fn record(&self, handle: u64, paper: bool, by_owner: &HashMap<u64, Vec<u64>>) -> R<Option<Entity>> {
        let mut o = self.object(handle)?;
        if !o.entity {
            return Ok(None);
        }
        if matches!(o.kind, t::ATTRIB | t::BLOCK | t::ENDBLK | t::SEQEND)
            || (t::VERTEX_2D..=t::VERTEX_PFACE_FACE).contains(&o.kind)
        {
            return Ok(None);
        }
        let (_, layer) = self.entity_refs(&mut o)?;
        let kind = self.type_name(o.kind);
        let mut e = Rec { kind, pairs: Vec::new(), children: Vec::new() };
        e.text(8, &self.layer_name(layer));
        e.int(62, o.color as i64);
        if o.invisible {
            e.int(60, 1);
        }
        if paper {
            e.int(67, 1);
        }
        let ver = self.ver();
        let custom = if o.kind >= 500 { self.type_name(o.kind) } else { String::new() };
        // Custom entities that are block references (e.g. ACIDBLOCKREFERENCE) begin with the INSERT data.
        let shape = match custom.as_str() {
            "LWPOLYLINE" => t::LWPOLYLINE,
            c if c.ends_with("BLOCKREFERENCE") => t::INSERT,
            _ => o.kind,
        };
        let h = o.handle;
        let b = &mut o.data;
        match shape {
            t::LINE => {
                let (p, q);
                if ver <= Ver::R14 {
                    p = b.p3()?;
                    q = b.p3()?;
                } else {
                    let flat = b.b()?;
                    let x0 = b.rd()?;
                    let x1 = b.dd(x0)?;
                    let y0 = b.rd()?;
                    let y1 = b.dd(y0)?;
                    let (z0, z1) = if flat { (0.0, 0.0) } else { let z0 = b.rd()?; (z0, b.dd(z0)?) };
                    p = [x0, y0, z0];
                    q = [x1, y1, z1];
                }
                b.bt()?;
                let n = b.be()?;
                e.point(10, p);
                e.point(11, q);
                e.point(210, n);
            }
            t::CIRCLE | t::ARC => {
                let c = b.p3()?;
                let r = b.bd()?;
                b.bt()?;
                let n = b.be()?;
                e.point(10, c);
                e.real(40, r);
                e.point(210, n);
                if o.kind == t::ARC {
                    e.real(50, b.bd()?.to_degrees());
                    e.real(51, b.bd()?.to_degrees());
                }
            }
            t::ELLIPSE => {
                e.point(10, b.p3()?);
                e.point(11, b.p3()?);
                e.point(210, b.p3()?);
                e.real(40, b.bd()?);
                e.real(41, b.bd()?);
                e.real(42, b.bd()?);
            }
            t::SPLINE => spline(b, ver, &mut e)?,
            t::INSERT | t::MINSERT => {
                e.kind = "INSERT".into();
                let p = b.p3()?;
                let s = if ver <= Ver::R14 {
                    [b.bd()?, b.bd()?, b.bd()?]
                } else {
                    match b.bb()? {
                        3 => [1.0; 3],
                        1 => [1.0, b.dd(1.0)?, b.dd(1.0)?],
                        2 => { let x = b.rd()?; [x; 3] }
                        _ => { let x = b.rd()?; [x, b.dd(x)?, b.dd(x)?] }
                    }
                };
                let rotation = b.bd()?;
                let n = b.p3()?;
                // The owned attribute count is stored only when attributes follow.
                let attributes = b.b()?;
                let owned = if attributes && ver >= Ver::R2004 { b.bl()?.max(0) as usize } else { 0 };
                e.point(10, p);
                e.real(41, s[0]);
                e.real(42, s[1]);
                e.real(43, s[2]);
                e.real(50, rotation.to_degrees());
                e.point(210, n);
                if o.kind == t::MINSERT {
                    e.int(70, b.bs()? as i64);
                    e.int(71, b.bs()? as i64);
                    e.real(44, b.bd()?);
                    e.real(45, b.bd()?);
                }
                let block = o.refs.h(h)?.value;
                let Some(info) = self.blocks.get(&block) else {
                    return bad(format!("insert {h:X} names block {block:X}, which is not a block header"));
                };
                e.text(2, &info.key);
                if attributes {
                    // R2004 on lists them; before, they are the entities the insert owns. A damaged list loses the
                    // attributes, not the insert.
                    let list: R<Vec<u64>> = if ver >= Ver::R2004 {
                        if owned > o.refs.left() / 8 || owned > MAX_ITEMS {
                            bad(format!("insert {h:X} owns {owned} attributes, more than its data holds"))
                        } else {
                            (0..owned).map(|_| o.refs.h(h).map(|x| x.value)).collect()
                        }
                    } else {
                        Ok(by_owner.get(&h).cloned().unwrap_or_default())
                    };
                    match list {
                        Ok(list) => for a in list {
                            let child = self.attrib(a);
                            if child.kind != "SEQEND" {
                                e.children.push(child);
                            }
                        },
                        Err(err) => e.children.push(refused_attrib(format!("attributes could not be found ({err})"))),
                    }
                }
            }
            t::POLYLINE_2D | t::POLYLINE_3D | t::POLYLINE_PFACE | t::POLYLINE_MESH => {
                let flags: i64;
                let mut elevation = 0.0;
                let mut normal = [0.0, 0.0, 1.0];
                match o.kind {
                    t::POLYLINE_2D => {
                        flags = b.bs()? as u16 as i64 & !(8 | 16 | 64);
                        b.bs()?;
                        b.bd()?;
                        b.bd()?;
                        b.bt()?;
                        elevation = b.bd()?;
                        normal = b.be()?;
                    }
                    t::POLYLINE_3D => {
                        let spline = b.rc()?;
                        let closed = b.rc()?;
                        flags = 8 | i64::from(closed & 1) | if spline & 3 != 0 { 4 } else { 0 };
                    }
                    t::POLYLINE_PFACE => flags = 64,
                    _ => flags = 16,
                }
                e.int(70, flags);
                e.point(10, [0.0, 0.0, elevation]);
                e.point(210, normal);
                if flags & (16 | 64) == 0 {
                    let vertices: Vec<u64> = if ver >= Ver::R2004 {
                        let count = b.bl()?.max(0) as usize;
                        if count > o.refs.left() / 8 {
                            return bad(format!("polyline {h:X} owns {count} vertices, more than its data holds"));
                        }
                        (0..count).map(|_| o.refs.h(h).map(|x| x.value)).collect::<R<_>>()?
                    } else {
                        by_owner.get(&h).cloned().unwrap_or_default()
                    };
                    for v in vertices {
                        e.children.push(self.vertex(v)?);
                    }
                }
            }
            t::LWPOLYLINE => {
                e.kind = "LWPOLYLINE".into();
                lwpolyline(b, ver, &mut e)?;
            }
            t::TEXT | t::MTEXT => {
                // A text that cannot be decoded stays a skipped TEXT or MTEXT (as for a curves-only read), with the
                // reason in 999 for the text reader.
                let mut decode = || -> R<()> {
                    if o.kind == t::TEXT {
                        text_data(b, ver, &mut o.strings, &mut e, 73)?;
                    } else {
                        e.point(10, b.p3()?);
                        e.point(210, b.p3()?);
                        e.point(11, b.p3()?);
                        e.real(41, b.bd()?);
                        if ver >= Ver::R2007 {
                            b.bd()?; // reference rectangle height
                        }
                        e.real(40, b.bd()?);
                        e.int(71, b.bs()? as i64);
                        e.int(72, b.bs()? as i64);
                        b.bd()?; // extents height
                        b.bd()?; // extents width
                        e.text(1, &unescape(&b.tv(&mut o.strings)?));
                    }
                    let style = o.refs.h(h)?.value;
                    e.text(7, &self.style_name(style));
                    Ok(())
                };
                if let Err(err) = decode() {
                    e.text(999, &format!("{} {h:X}: {err}", e.kind));
                }
            }
            _ => {}
        }
        Ok(Some(e.done()))
    }

    fn vertex(&self, handle: u64) -> R<Entity> {
        let mut o = self.object(handle)?;
        let mut e = Rec { kind: "VERTEX".into(), pairs: Vec::new(), children: Vec::new() };
        let b = &mut o.data;
        match o.kind {
            t::VERTEX_2D => {
                let flags = b.rc()?;
                let p = b.p3()?;
                let start = b.bd()?;
                if start >= 0.0 {
                    b.bd()?;
                }
                let bulge = b.bd()?;
                e.int(70, flags as i64);
                e.point(10, p);
                e.real(42, bulge);
            }
            t::VERTEX_3D => {
                let flags = b.rc()?;
                e.int(70, flags as i64 | 32);
                e.point(10, b.p3()?);
            }
            k => return bad(format!("polyline vertex {handle:X} has type {k}")),
        }
        Ok(e.done())
    }
}

/// An attribute the reader does not give as text, its reason as the record kind (the text reader reports it).
fn refused_attrib(reason: String) -> Entity {
    Entity { kind: reason, pairs: Vec::new(), children: Vec::new() }
}

/// The data TEXT, ATTRIB and ATTDEF share (spec 20.4.3) as DXF codes; `vertical` is 73 for TEXT, 74 for ATTRIB.
fn text_data<'a>(b: &mut Bits<'a>, ver: Ver, strings: &mut Option<Bits<'a>>, e: &mut Rec, vertical: i32) -> R<()> {
    let (elevation, insertion, alignment, normal, rotation, height, width, value, generation, horizontal, upright);
    if ver <= Ver::R14 {
        elevation = b.bd()?;
        insertion = [b.rd()?, b.rd()?];
        alignment = Some([b.rd()?, b.rd()?]);
        normal = b.p3()?;
        b.bd()?; // thickness
        b.bd()?; // oblique angle
        rotation = b.bd()?;
        height = b.bd()?;
        width = b.bd()?;
        value = b.tv(strings)?;
        generation = b.bs()?;
        horizontal = b.bs()?;
        upright = b.bs()?;
    } else {
        // Data flags: a set bit means the value is not stored (and takes its default).
        let flags = b.rc()?;
        elevation = if flags & 0x01 == 0 { b.rd()? } else { 0.0 };
        insertion = [b.rd()?, b.rd()?];
        alignment = if flags & 0x02 == 0 {
            let x = b.dd(insertion[0])?;
            Some([x, b.dd(insertion[1])?])
        } else {
            None
        };
        normal = b.be()?;
        b.bt()?; // thickness
        if flags & 0x04 == 0 {
            b.rd()?; // oblique angle
        }
        rotation = if flags & 0x08 == 0 { b.rd()? } else { 0.0 };
        height = b.rd()?;
        width = if flags & 0x10 == 0 { b.rd()? } else { 1.0 };
        value = b.tv(strings)?;
        generation = if flags & 0x20 == 0 { b.bs()? } else { 0 };
        horizontal = if flags & 0x40 == 0 { b.bs()? } else { 0 };
        upright = if flags & 0x80 == 0 { b.bs()? } else { 0 };
    }
    e.point(10, [insertion[0], insertion[1], elevation]);
    if let Some(a) = alignment {
        e.point(11, [a[0], a[1], elevation]);
    }
    e.point(210, normal);
    e.real(50, rotation.to_degrees());
    e.real(40, height);
    e.real(41, width);
    e.text(1, &unescape(&value));
    e.int(71, generation as i64);
    e.int(72, horizontal as i64);
    e.int(vertical, upright as i64);
    Ok(())
}

/// The first string item (code 0, DXF 1000) of one application's extended data (spec chapter 28).
fn eed_string(b: &Bits, ver: Ver) -> Option<String> {
    let mut b = b.clone();
    for _ in 0..1000 {
        let code = b.rc().ok()?;
        match code {
            0 => {
                if ver >= Ver::R2007 {
                    let n = b.rs().ok()? as usize;
                    let units: Vec<u16> = (0..n).map(|_| b.rs()).collect::<R<_>>().ok()?;
                    return Some(String::from_utf16_lossy(&units).trim_end_matches('\0').to_string());
                }
                let n = b.rc().ok()? as usize;
                b.rs().ok()?; // code page
                let bytes: Vec<u8> = (0..n).map(|_| b.rc()).collect::<R<_>>().ok()?;
                return Some(crate::dxf::codes::decode(&bytes));
            }
            2 => b.skip(8).ok()?,
            3 | 5 | 40..=42 => b.skip(64).ok()?,
            4 => {
                let n = b.rc().ok()? as usize;
                b.skip(n * 8).ok()?;
            }
            10..=13 => b.skip(192).ok()?,
            70 => b.skip(16).ok()?,
            71 => b.skip(32).ok()?,
            _ => return None,
        }
    }
    None
}

fn spline(b: &mut Bits, ver: Ver, e: &mut Rec) -> R<()> {
    let scenario = b.bl()?;
    if ver >= Ver::R2013 {
        b.bl()?;
        b.bl()?;
    }
    let degree = b.bl()?;
    let (mut flags, mut knots, mut control, mut fit, mut weighted) = (8i64, 0usize, 0usize, 0usize, false);
    let mut tangents = None;
    match scenario {
        2 => {
            b.bd()?;
            let t0 = b.p3()?;
            let t1 = b.p3()?;
            tangents = Some((t0, t1));
            fit = b.count(3 * 2, "fit point")?;
        }
        1 => {
            let rational = b.b()?;
            let closed = b.b()?;
            let periodic = b.b()?;
            flags |= i64::from(closed) | i64::from(periodic) << 1 | i64::from(rational) << 2;
            b.bd()?;
            b.bd()?;
            knots = b.count(2, "knot")?;
            control = b.count(6, "control point")?;
            weighted = b.b()?;
        }
        s => return bad(format!("spline kind {s}")),
    }
    e.int(70, flags);
    e.int(71, degree as i64);
    e.int(72, knots as i64);
    e.int(73, control as i64);
    e.int(74, fit as i64);
    for _ in 0..knots {
        e.real(40, b.bd()?);
    }
    for _ in 0..control {
        e.point(10, b.p3()?);
        if weighted {
            e.real(41, b.bd()?);
        }
    }
    for _ in 0..fit {
        e.point(11, b.p3()?);
    }
    if let Some((t0, t1)) = tangents {
        if t0.iter().any(|v| *v != 0.0) && t1.iter().any(|v| *v != 0.0) {
            e.point(12, t0);
            e.point(13, t1);
        }
    }
    Ok(())
}

fn lwpolyline(b: &mut Bits, ver: Ver, e: &mut Rec) -> R<()> {
    let flag = b.bs()? as u16;
    if flag & 4 != 0 {
        b.bd()?;
    }
    let elevation = if flag & 8 != 0 { b.bd()? } else { 0.0 };
    if flag & 2 != 0 {
        b.bd()?;
    }
    let normal = if flag & 1 != 0 { b.p3()? } else { [0.0, 0.0, 1.0] };
    let n = b.count(2, "vertex")?;
    let bulges = if flag & 16 != 0 { b.count(2, "bulge")? } else { 0 };
    let ids = if ver >= Ver::R2010 && flag & 1024 != 0 { b.count(2, "vertex id")? } else { 0 };
    let widths = if flag & 32 != 0 { b.count(4, "width")? } else { 0 };
    let mut points = Vec::with_capacity(n);
    for i in 0..n {
        if ver <= Ver::R14 || i == 0 {
            points.push([b.rd()?, b.rd()?]);
        } else {
            let prev = points[i - 1];
            let x = b.dd(prev[0])?;
            let y = b.dd(prev[1])?;
            points.push([x, y]);
        }
    }
    let mut bulge = Vec::with_capacity(bulges);
    for _ in 0..bulges {
        bulge.push(b.bd()?);
    }
    for _ in 0..ids {
        b.bl()?;
    }
    for _ in 0..widths {
        b.bd()?;
        b.bd()?;
    }
    e.int(90, n as i64);
    e.int(70, i64::from(flag & 512 != 0));
    e.real(38, elevation);
    e.point(210, normal);
    for (i, p) in points.iter().enumerate() {
        e.real(10, p[0]);
        e.real(20, p[1]);
        if let Some(&g) = bulge.get(i) {
            if g != 0.0 {
                e.real(42, g);
            }
        }
    }
    Ok(())
}

/// The R2007+ string stream of an object whose data ends at bit `end` (spec 20.1): a flag bit just before `end`,
/// then going backwards a 15-bit size (two shorts for 30 bits) and the strings.
fn string_stream<'a>(b: &Bits<'a>, end: usize) -> Option<Bits<'a>> {
    if end < 1 {
        return None;
    }
    let mut p = end - 1;
    let mut flag = b.window(p, end);
    if !flag.b().ok()? {
        return None;
    }
    p = p.checked_sub(16)?;
    let mut size = b.window(p, p + 16).rs().ok()? as usize;
    if size & 0x8000 != 0 {
        p = p.checked_sub(16)?;
        let hi = b.window(p, p + 16).rs().ok()? as usize;
        size = (size & 0x7fff) | (hi << 15);
    }
    let start = p.checked_sub(size)?;
    Some(b.window(start, p))
}

/// An entity record under construction.
struct Rec {
    kind: String,
    pairs: Vec<Pair>,
    children: Vec<Entity>,
}

impl Rec {
    fn push(&mut self, code: i32, value: Value) {
        self.pairs.push(Pair { code, value, at: 0, binary: true });
    }
    fn real(&mut self, code: i32, v: f64) {
        self.push(code, Value::Real(v));
    }
    fn int(&mut self, code: i32, v: i64) {
        self.push(code, Value::Int(v));
    }
    fn text(&mut self, code: i32, v: &str) {
        self.push(code, Value::Text(v.to_string()));
    }
    fn point(&mut self, code: i32, p: [f64; 3]) {
        self.real(code, p[0]);
        self.real(code + 10, p[1]);
        self.real(code + 20, p[2]);
    }
    fn done(self) -> Entity {
        Entity { kind: self.kind, pairs: self.pairs, children: self.children }
    }
}

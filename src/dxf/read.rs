//! DXF to `Drawing`: sections, the layer table, block definitions and model-space entities, every curve mapped to
//! world XY through its object coordinate system and any (nested, arrayed, mirrored) block inserts.

use std::collections::HashMap;

use super::codes::{pairs, Pair};
use super::text::{mtext_plain, text_plain};
use super::{fit, Curve, Drawing, Geometry, Layer, Text, TextKind, TextStyle, P2};
use crate::KernelError;

type P3 = [f64; 3];

/// One entity record: its type and pairs; a POLYLINE's VERTEX records (or an INSERT's ATTRIBs) as children. The
/// DWG reader builds the same records (DXF group codes) so both formats share everything below.
pub(crate) struct Entity {
    pub(crate) kind: String,
    pub(crate) pairs: Vec<Pair>,
    pub(crate) children: Vec<Entity>,
}

/// A block definition; `blocks` maps are keyed by the upper-cased block name.
pub(crate) struct Block {
    pub(crate) base: P3,
    pub(crate) entities: Vec<Entity>,
    pub(crate) xref: bool,
}

/// x -> m x + t.
#[derive(Clone, Copy)]
struct Affine {
    m: [[f64; 3]; 3],
    t: P3,
}

impl Affine {
    const IDENTITY: Affine = Affine { m: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], t: [0.0; 3] };

    fn point(&self, p: P3) -> P3 {
        let v = self.vector(p);
        [v[0] + self.t[0], v[1] + self.t[1], v[2] + self.t[2]]
    }

    fn vector(&self, v: P3) -> P3 {
        let m = &self.m;
        [
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
        ]
    }

    /// self after other: x -> self(other(x)).
    fn then(&self, other: &Affine) -> Affine {
        let mut m = [[0.0; 3]; 3];
        for (i, row) in m.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = (0..3).map(|k| self.m[i][k] * other.m[k][j]).sum();
            }
        }
        Affine { m, t: self.point(other.t) }
    }

    fn translate(t: P3) -> Affine {
        Affine { t, ..Affine::IDENTITY }
    }

    fn scale(s: P3) -> Affine {
        Affine { m: [[s[0], 0.0, 0.0], [0.0, s[1], 0.0], [0.0, 0.0, s[2]]], t: [0.0; 3] }
    }

    fn rotate_z(angle: f64) -> Affine {
        let (s, c) = angle.sin_cos();
        Affine { m: [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], t: [0.0; 3] }
    }

    /// Object coordinates to world for extrusion direction n (DXF Reference, "Arbitrary Axis Algorithm").
    fn ocs(n: P3) -> Affine {
        let n = unit(n).unwrap_or([0.0, 0.0, 1.0]);
        if n == [0.0, 0.0, 1.0] {
            return Affine::IDENTITY;
        }
        let ax = if n[0].abs() < 1.0 / 64.0 && n[1].abs() < 1.0 / 64.0 {
            cross([0.0, 1.0, 0.0], n)
        } else {
            cross([0.0, 0.0, 1.0], n)
        };
        let ax = unit(ax).unwrap_or([1.0, 0.0, 0.0]);
        let ay = unit(cross(n, ax)).unwrap_or([0.0, 1.0, 0.0]);
        Affine { m: [[ax[0], ay[0], n[0]], [ax[1], ay[1], n[1]], [ax[2], ay[2], n[2]]], t: [0.0; 3] }
    }
}

/// The counter-clockwise sweep from angle a to angle b (both in [0, 2 pi)); the whole turn when they meet.
fn sweep_of(a: f64, b: f64) -> f64 {
    let sweep = (b - a).rem_euclid(std::f64::consts::TAU);
    if sweep < 1e-12 { std::f64::consts::TAU } else { sweep }
}

fn cross(a: P3, b: P3) -> P3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn norm(a: P3) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn unit(a: P3) -> Option<P3> {
    let n = norm(a);
    (n > 0.0 && n.is_finite()).then(|| [a[0] / n, a[1] / n, a[2] / n])
}

/// Reads a DXF file's bytes.
pub fn read(bytes: &[u8]) -> Result<Drawing, KernelError> {
    read_parts(bytes, false)
}

/// Only the header and layer table (for asking about units before an import).
pub fn read_header(bytes: &[u8]) -> Result<Drawing, KernelError> {
    read_parts(bytes, true)
}

fn read_parts(bytes: &[u8], header_only: bool) -> Result<Drawing, KernelError> {
    let all = pairs(bytes)?;
    let mut drawing = Drawing::default();
    let mut blocks: HashMap<String, Block> = HashMap::new();
    let mut model: Vec<Entity> = Vec::new();
    let mut sections = 0;
    let mut i = 0;
    while i < all.len() {
        let p = &all[i];
        if p.code == 0 && p.text().eq_ignore_ascii_case("SECTION") {
            sections += 1;
            let name = all.get(i + 1).filter(|n| n.code == 2).map(|n| n.text().to_ascii_uppercase()).unwrap_or_default();
            let start = i + 2;
            let end = (start..all.len())
                .find(|&k| all[k].code == 0 && all[k].text().eq_ignore_ascii_case("ENDSEC"))
                .unwrap_or(all.len());
            let body = &all[start.min(end)..end];
            match name.as_str() {
                "HEADER" => header(body, &mut drawing),
                "TABLES" => tables(body, &mut drawing)?,
                "BLOCKS" if !header_only => definitions(body, &mut blocks)?,
                "ENTITIES" if !header_only => model.extend(records(body)),
                _ => {}
            }
            i = end + 1;
        } else if p.code == 0 && p.text().eq_ignore_ascii_case("EOF") {
            break;
        } else {
            i += 1;
        }
    }
    if sections == 0 {
        return Err(KernelError::InvalidInput("this is not a DXF file (it has no SECTION)".into()));
    }
    assemble(drawing, &blocks, &model)
}

/// Model-space entities (and the blocks they insert) as the drawing seen from above. `drawing` carries the header and
/// layer table already read; entities flagged 67 = 1 are paper space.
pub(crate) fn assemble(mut drawing: Drawing, blocks: &HashMap<String, Block>, model: &[Entity])
    -> Result<Drawing, KernelError>
{
    if drawing.layers.is_empty() {
        drawing.layer_index("0");
    }

    let mut x = Explode { blocks, stack: Vec::new(), notes: Notes::default() };
    let mut paper = Vec::new();
    for (group, entity) in model.iter().enumerate() {
        if flag(entity, 67) == 1 {
            // A layout's VIEWPORT is how paper space shows model space, not drawing content.
            if entity.kind != "VIEWPORT" {
                paper.push((group, entity));
            }
            continue;
        }
        x.entity(&mut drawing, entity, &Affine::IDENTITY, None, 256, group)?;
        too_many(&drawing)?;
    }
    if drawing.curves.is_empty() && !paper.is_empty() {
        // A sheet-style export: the drawing is on the layout and model space is empty.
        for (group, entity) in &paper {
            x.entity(&mut drawing, entity, &Affine::IDENTITY, None, 256, *group)?;
            too_many(&drawing)?;
        }
        if !drawing.curves.is_empty() {
            drawing.warn("model space is empty; the curves on the paper-space layout were read".into());
            paper.clear();
        }
    }
    x.notes.paper = paper.len();
    x.notes.finish(&mut drawing);
    Ok(drawing)
}

fn header(body: &[Pair], d: &mut Drawing) {
    let mut name = String::new();
    for p in body {
        if p.code == 9 {
            name = p.text().to_ascii_uppercase();
            continue;
        }
        match (name.as_str(), p.code) {
            ("$ACADVER", 1) => d.version = p.text(),
            ("$INSUNITS", 70) => d.insunits = p.int().ok().map(|v| v as i32),
            ("$MEASUREMENT", 70) => d.measurement = p.int().ok().map(|v| v as i32),
            _ => {}
        }
    }
}

fn tables(body: &[Pair], d: &mut Drawing) -> Result<(), KernelError> {
    let all = records(body);
    for record in all.iter().filter(|r| r.kind == "STYLE") {
        let name = text(record, 2).unwrap_or_default();
        // Flag 1: a shape file (symbols), not a text font.
        if name.is_empty() || flag(record, 70) & 1 != 0 {
            continue;
        }
        // The TrueType family is the first 1000 string of the ACAD extended data.
        let mut family = String::new();
        let mut acad = false;
        for p in &record.pairs {
            match p.code {
                1001 => acad = p.text().eq_ignore_ascii_case("ACAD"),
                1000 if acad && family.is_empty() => family = p.text(),
                _ => {}
            }
        }
        d.styles.push(TextStyle { name, file: text(record, 3).unwrap_or_default(), family });
    }
    for record in all.iter().filter(|r| r.kind == "LAYER") {
        let name = text(record, 2).unwrap_or_default();
        if name.is_empty() {
            continue;
        }
        let color = int(record, 62)?.unwrap_or(7) as i32;
        let flags = int(record, 70)?.unwrap_or(0);
        let index = d.layer_index(&name);
        d.layers[index] = Layer {
            name: d.layers[index].name.clone(),
            color: color.abs(),
            frozen: flags & 1 != 0,
            off: color < 0,
            locked: flags & 4 != 0,
        };
    }
    Ok(())
}

fn definitions(body: &[Pair], blocks: &mut HashMap<String, Block>) -> Result<(), KernelError> {
    let mut current: Option<(String, Block)> = None;
    for record in records(body) {
        match record.kind.as_str() {
            "BLOCK" => {
                let name = text(&record, 2).or_else(|| text(&record, 3)).unwrap_or_default();
                let flags = int(&record, 70)?.unwrap_or(0);
                let base = point(&record, 10)?.unwrap_or([0.0; 3]);
                current = Some((name.to_ascii_uppercase(), Block { base, entities: Vec::new(), xref: flags & 4 != 0 }));
            }
            "ENDBLK" => {
                if let Some((name, block)) = current.take() {
                    blocks.insert(name, block);
                }
            }
            _ => {
                if let Some((_, block)) = current.as_mut() {
                    block.entities.push(record);
                }
            }
        }
    }
    if let Some((name, block)) = current.take() {
        blocks.insert(name, block);
    }
    Ok(())
}

/// Splits a section body into entity records, attaching VERTEX/ATTRIB runs (to SEQEND) to their POLYLINE/INSERT.
fn records(body: &[Pair]) -> Vec<Entity> {
    let mut flat: Vec<Entity> = Vec::new();
    for p in body {
        if p.code == 0 {
            flat.push(Entity { kind: p.text().to_ascii_uppercase(), pairs: Vec::new(), children: Vec::new() });
        } else if let Some(last) = flat.last_mut() {
            last.pairs.push(p.clone());
        }
    }
    let mut out: Vec<Entity> = Vec::with_capacity(flat.len());
    let mut iter = flat.into_iter().peekable();
    while let Some(mut e) = iter.next() {
        if e.kind == "POLYLINE" || e.kind == "INSERT" {
            let child = if e.kind == "POLYLINE" { "VERTEX" } else { "ATTRIB" };
            while iter.peek().is_some_and(|n| n.kind == child) {
                e.children.push(iter.next().unwrap());
            }
            if iter.peek().is_some_and(|n| n.kind == "SEQEND") {
                iter.next();
            }
        }
        out.push(e);
    }
    out
}

fn find<'a>(e: &'a Entity, code: i32) -> Option<&'a Pair> {
    e.pairs.iter().find(|p| p.code == code)
}

fn text(e: &Entity, code: i32) -> Option<String> {
    find(e, code).map(|p| p.text())
}

fn real(e: &Entity, code: i32) -> Result<Option<f64>, KernelError> {
    find(e, code).map(|p| p.real()).transpose()
}

fn int(e: &Entity, code: i32) -> Result<Option<i64>, KernelError> {
    find(e, code).map(|p| p.int()).transpose()
}

fn flag(e: &Entity, code: i32) -> i64 {
    int(e, code).ok().flatten().unwrap_or(0)
}

fn point(e: &Entity, code: i32) -> Result<Option<P3>, KernelError> {
    let Some(x) = real(e, code)? else { return Ok(None) };
    Ok(Some([x, real(e, code + 10)?.unwrap_or(0.0), real(e, code + 20)?.unwrap_or(0.0)]))
}

fn extrusion(e: &Entity) -> Result<P3, KernelError> {
    Ok(point(e, 210)?.unwrap_or([0.0, 0.0, 1.0]))
}

/// Every value of a repeated point code, in order (LWPOLYLINE vertices, SPLINE control and fit points).
fn points(e: &Entity, code: i32) -> Result<Vec<P3>, KernelError> {
    let mut out: Vec<P3> = Vec::new();
    for p in &e.pairs {
        if p.code == code {
            out.push([p.real()?, 0.0, 0.0]);
        } else if p.code == code + 10 {
            if let Some(last) = out.last_mut() { last[1] = p.real()?; }
        } else if p.code == code + 20 {
            if let Some(last) = out.last_mut() { last[2] = p.real()?; }
        }
    }
    Ok(out)
}

fn reals(e: &Entity, code: i32) -> Result<Vec<f64>, KernelError> {
    e.pairs.iter().filter(|p| p.code == code).map(|p| p.real()).collect()
}

#[derive(Default)]
struct Notes {
    paper: usize,
    off_plane: usize,
    edge_on: usize,
    zero: usize,
    text_edge_on: usize,
    text_no_height: usize,
    text_mirrored: usize,
    text_skewed: usize,
    text_vertical: usize,
    text_bad: usize,
    text_bad_first: String,
}

impl Notes {
    fn bad_text(&mut self, reason: String) {
        self.text_bad += 1;
        if self.text_bad_first.is_empty() {
            self.text_bad_first = reason;
        }
    }

    fn finish(&self, d: &mut Drawing) {
        if self.paper > 0 {
            d.warn(format!("{} entities on a paper-space layout were left out (only model space is read)", self.paper));
        }
        if self.off_plane > 0 {
            d.warn(format!("{} curves are not in the XY plane; they are read as seen from above", self.off_plane));
        }
        if self.edge_on > 0 {
            d.warn(format!("{} arcs or circles stand on edge when seen from above and were left out", self.edge_on));
        }
        if self.zero > 0 {
            d.warn(format!("{} zero-length curves were left out", self.zero));
        }
        if self.text_bad > 0 {
            d.warn_text(format!("{} texts could not be read and were left out (the first: {})", self.text_bad,
                self.text_bad_first));
        }
        if self.text_edge_on > 0 {
            d.warn_text(format!("{} texts stand on edge when seen from above and were left out", self.text_edge_on));
        }
        if self.text_no_height > 0 {
            d.warn_text(format!("{} texts have no height and were left out", self.text_no_height));
        }
        if self.text_mirrored > 0 {
            d.warn_text(format!("{} texts are mirrored in the drawing; they are read readable, over the same place",
                self.text_mirrored));
        }
        if self.text_skewed > 0 {
            d.warn_text(format!("{} texts are slanted by a non-uniformly scaled insert; they are read upright", self.text_skewed));
        }
        if self.text_vertical > 0 {
            d.warn_text(format!("{} MTEXTs are written top to bottom; they are read left to right", self.text_vertical));
        }
    }
}

struct Explode<'a> {
    blocks: &'a HashMap<String, Block>,
    stack: Vec<String>,
    notes: Notes,
}

const MAX_DEPTH: usize = 64;
const MAX_ARRAY: i64 = 100_000;
const MAX_CURVES: usize = 5_000_000;

fn too_many(d: &Drawing) -> Result<(), KernelError> {
    if d.curves.len() > MAX_CURVES {
        return Err(KernelError::Unsupported(format!("the drawing expands to more than {MAX_CURVES} curves")));
    }
    if d.texts.len() > MAX_CURVES {
        return Err(KernelError::Unsupported(format!("the drawing expands to more than {MAX_CURVES} texts")));
    }
    Ok(())
}

/// An entity's layer and colour: layer 0 and colour 0 (BYBLOCK) take the enclosing insert's.
fn resolve(d: &mut Drawing, e: &Entity, parent: Option<usize>, parent_color: i32) -> Result<(usize, i32), KernelError> {
    let layer_name = text(e, 8).unwrap_or_else(|| "0".into());
    let layer = match parent {
        Some(p) if layer_name.trim() == "0" || layer_name.trim().is_empty() => p,
        _ => d.layer_index(&layer_name),
    };
    let color = match int(e, 62)?.unwrap_or(256) as i32 {
        0 => parent_color,
        c if c < 0 => c.abs(),
        c => c,
    };
    Ok((layer, color))
}

impl Explode<'_> {
    /// One entity under world transform `m`; `parent` is the enclosing insert's (layer, colour).
    fn entity(&mut self, d: &mut Drawing, e: &Entity, m: &Affine, parent: Option<usize>, parent_color: i32, group: usize)
        -> Result<(), KernelError>
    {
        if flag(e, 60) == 1 {
            return Ok(());
        }
        let (layer, color) = resolve(d, e, parent, parent_color)?;
        if e.kind == "INSERT" {
            return self.insert(d, e, m, layer, color, group);
        }
        if e.kind == "TEXT" || e.kind == "MTEXT" {
            // Still counted as not a curve; the text itself goes to `texts`.
            d.skip(&e.kind);
            if let Err(err) = self.text(d, e, m, layer, color, group) {
                self.notes.bad_text(err.to_string());
            }
            return Ok(());
        }
        let mut out = Emit { d, layer, color, group, notes: &mut self.notes };
        match e.kind.as_str() {
            "LINE" => {
                let a = point(e, 10)?.unwrap_or([0.0; 3]);
                let b = point(e, 11)?.unwrap_or([0.0; 3]);
                out.line(m.point(a), m.point(b));
            }
            "CIRCLE" | "ARC" => {
                let ocs = m.then(&Affine::ocs(extrusion(e)?));
                let c = point(e, 10)?.unwrap_or([0.0; 3]);
                let r = real(e, 40)?.unwrap_or(0.0).abs();
                if r <= 0.0 {
                    out.notes.zero += 1;
                    return Ok(());
                }
                let (t0, t1) = if e.kind == "CIRCLE" {
                    (0.0, std::f64::consts::TAU)
                } else {
                    // Counter-clockwise from 50 to 51; equal angles draw the whole circle.
                    let a0 = real(e, 50)?.unwrap_or(0.0).to_radians().rem_euclid(std::f64::consts::TAU);
                    let a1 = real(e, 51)?.unwrap_or(0.0).to_radians().rem_euclid(std::f64::consts::TAU);
                    (a0, a0 + sweep_of(a0, a1))
                };
                out.conic(ocs.point(c), ocs.vector([r, 0.0, 0.0]), ocs.vector([0.0, r, 0.0]), t0, t1);
            }
            "ELLIPSE" => {
                let c = point(e, 10)?.unwrap_or([0.0; 3]);
                let major = point(e, 11)?.unwrap_or([0.0; 3]);
                let ratio = real(e, 40)?.unwrap_or(1.0);
                let n = unit(extrusion(e)?).unwrap_or([0.0, 0.0, 1.0]);
                let minor = cross(n, major);
                let minor = [minor[0] * ratio, minor[1] * ratio, minor[2] * ratio];
                let t0 = real(e, 41)?.unwrap_or(0.0).rem_euclid(std::f64::consts::TAU);
                let t1 = t0 + sweep_of(t0, real(e, 42)?.unwrap_or(std::f64::consts::TAU).rem_euclid(std::f64::consts::TAU));
                if norm(major) <= 0.0 || ratio <= 0.0 {
                    out.notes.zero += 1;
                    return Ok(());
                }
                out.conic(m.point(c), m.vector(major), m.vector(minor), t0, t1);
            }
            "LWPOLYLINE" => {
                let ocs = m.then(&Affine::ocs(extrusion(e)?));
                let elevation = real(e, 38)?.unwrap_or(0.0);
                let mut vertices: Vec<(P3, f64)> = Vec::new();
                for p in &e.pairs {
                    match p.code {
                        10 => vertices.push(([p.real()?, 0.0, elevation], 0.0)),
                        20 => if let Some(v) = vertices.last_mut() { v.0[1] = p.real()?; },
                        42 => if let Some(v) = vertices.last_mut() { v.1 = p.real()?; },
                        _ => {}
                    }
                }
                out.polyline(&ocs, &vertices, flag(e, 70) & 1 != 0);
            }
            "POLYLINE" => {
                let flags = flag(e, 70);
                if flags & (16 | 64) != 0 {
                    out.d.skip("POLYLINE mesh");
                    return Ok(());
                }
                let three_d = flags & 8 != 0;
                let ocs = if three_d { *m } else { m.then(&Affine::ocs(extrusion(e)?)) };
                let elevation = point(e, 10)?.map(|p| p[2]).unwrap_or(0.0);
                let mut vertices = Vec::new();
                for v in &e.children {
                    if flag(v, 70) & 16 != 0 {
                        continue; // spline frame control point: not drawn
                    }
                    let mut p = point(v, 10)?.unwrap_or([0.0; 3]);
                    if !three_d { p[2] = elevation; }
                    vertices.push((p, if three_d { 0.0 } else { real(v, 42)?.unwrap_or(0.0) }));
                }
                out.polyline(&ocs, &vertices, flags & 1 != 0);
            }
            "SPLINE" => out.spline(e, m)?,
            kind => out.d.skip(kind),
        }
        Ok(())
    }

    fn insert(&mut self, d: &mut Drawing, e: &Entity, m: &Affine, layer: usize, color: i32, group: usize)
        -> Result<(), KernelError>
    {
        // Attribute values sit in the insert's own coordinates, not the block's.
        for attrib in &e.children {
            if attrib.kind != "ATTRIB" {
                // The DWG reader's attributes it could not give as text, the reason as their kind.
                self.notes.bad_text(attrib.kind.clone());
                continue;
            }
            if flag(attrib, 60) == 1 || flag(attrib, 70) & 1 != 0 {
                continue; // invisible
            }
            let placed = resolve(d, attrib, Some(layer), color).and_then(|(l, c)| self.text(d, attrib, m, l, c, group));
            if let Err(err) = placed {
                self.notes.bad_text(err.to_string());
            }
        }
        let name = text(e, 2).unwrap_or_default();
        let key = name.to_ascii_uppercase();
        let Some(block) = self.blocks.get(&key) else {
            d.warn(format!("block '{name}' is inserted but not defined; it was left out"));
            return Ok(());
        };
        if block.xref {
            d.warn(format!("block '{name}' is an external reference; its drawing is not in this file and was left out"));
            return Ok(());
        }
        if self.stack.contains(&key) || self.stack.len() >= MAX_DEPTH {
            d.warn(format!("block '{name}' contains itself; the repeated insert was left out"));
            return Ok(());
        }
        let ins = point(e, 10)?.unwrap_or([0.0; 3]);
        let s = [real(e, 41)?.unwrap_or(1.0), real(e, 42)?.unwrap_or(1.0), real(e, 43)?.unwrap_or(1.0)];
        if s[0] == 0.0 || s[1] == 0.0 {
            d.warn(format!("an insert of block '{name}' has zero scale and was left out"));
            return Ok(());
        }
        let rotation = real(e, 50)?.unwrap_or(0.0).to_radians();
        let columns = int(e, 70)?.unwrap_or(1).max(1);
        let rows = int(e, 71)?.unwrap_or(1).max(1);
        if columns.saturating_mul(rows) > MAX_ARRAY {
            return Err(KernelError::Unsupported(format!(
                "an insert of block '{name}' repeats it {columns} x {rows} times (more than {MAX_ARRAY})")));
        }
        let (dx, dy) = (real(e, 44)?.unwrap_or(0.0), real(e, 45)?.unwrap_or(0.0));
        let ocs = Affine::ocs(extrusion(e)?);
        let place = m.then(&ocs).then(&Affine::translate(ins)).then(&Affine::rotate_z(rotation));
        let local = Affine::scale(s).then(&Affine::translate([-block.base[0], -block.base[1], -block.base[2]]));
        self.stack.push(key);
        for row in 0..rows {
            for column in 0..columns {
                let offset = Affine::translate([column as f64 * dx, row as f64 * dy, 0.0]);
                let total = place.then(&offset).then(&local);
                for child in &block.entities {
                    self.entity(d, child, &total, Some(layer), color, group)?;
                    too_many(d)?;
                }
            }
        }
        self.stack.pop();
        Ok(())
    }
}

impl Explode<'_> {
    /// A TEXT, MTEXT or ATTRIB under world transform `m`, as seen from above.
    fn text(&mut self, d: &mut Drawing, e: &Entity, m: &Affine, layer: usize, color: i32, group: usize)
        -> Result<(), KernelError>
    {
        if let Some(reason) = text(e, 999) {
            // The DWG reader could not decode this text's data.
            self.notes.bad_text(reason);
            return Ok(());
        }
        let kind = match e.kind.as_str() {
            "MTEXT" => TextKind::MText,
            "ATTRIB" => TextKind::Attrib,
            _ => TextKind::Text,
        };
        let value = if kind == TextKind::MText {
            // Text over 250 characters comes in 3 chunks, the last piece in 1.
            let mut s: String = e.pairs.iter().filter(|p| p.code == 3).map(|p| p.raw()).collect();
            s.push_str(&find(e, 1).map(|p| p.raw()).unwrap_or_default());
            mtext_plain(&s).trim_end().to_string()
        } else {
            text_plain(&find(e, 1).map(|p| p.raw()).unwrap_or_default()).trim_end().to_string()
        };
        if value.trim().is_empty() {
            return Ok(()); // draws nothing
        }
        let height = real(e, 40)?.unwrap_or(0.0);
        if !(height > 0.0) {
            self.notes.text_no_height += 1;
            return Ok(());
        }
        let style = text(e, 7).filter(|s| !s.is_empty()).unwrap_or_else(|| "Standard".into());
        let p10 = point(e, 10)?.unwrap_or([0.0; 3]);
        let (at, end, x, up, mut halign, valign, mut width, mut box_width);
        if kind == TextKind::MText {
            // MTEXT's point and its 11 direction are world (here: block) coordinates. Without 11, 50 turns the
            // object coordinate system's X axis; it is read in degrees like every other entity's angle (the DXF
            // Reference lists radians for this one code), and 11 wins when both are given.
            let n = unit(extrusion(e)?).unwrap_or([0.0, 0.0, 1.0]);
            let direction = match point(e, 11)? {
                Some(v) if norm(v) > 0.0 => v,
                _ => {
                    let r = real(e, 50)?.unwrap_or(0.0).to_radians();
                    Affine::ocs(n).vector([r.cos(), r.sin(), 0.0])
                }
            };
            let direction = unit(direction).unwrap_or([1.0, 0.0, 0.0]);
            let attachment = flag(e, 71).clamp(1, 9) as i32;
            halign = (attachment - 1) % 3;
            valign = [3, 2, 1][((attachment - 1) / 3) as usize];
            if flag(e, 72) == 3 {
                self.notes.text_vertical += 1;
            }
            at = m.point(p10);
            end = at;
            x = m.vector(direction);
            up = m.vector(cross(n, direction));
            width = 1.0;
            box_width = real(e, 41)?.unwrap_or(0.0).max(0.0);
        } else {
            let h72 = flag(e, 72).clamp(0, 5) as i32;
            let v = flag(e, if kind == TextKind::Attrib { 74 } else { 73 }).clamp(0, 3) as i32;
            let p11 = point(e, 11)?;
            let mut r = real(e, 50)?.unwrap_or(0.0).to_radians();
            let (anchor, second) = match (h72, v) {
                (3 | 5, _) => (p10, p11.unwrap_or(p10)),
                (0, 0) => (p10, p10),
                _ => (p11.unwrap_or(p10), p11.unwrap_or(p10)),
            };
            if (h72 == 3 || h72 == 5) && (second[0] - anchor[0]).hypot(second[1] - anchor[1]) > 0.0 {
                // Aligned and fit text run from the first point to the second.
                r = (second[1] - anchor[1]).atan2(second[0] - anchor[0]);
            }
            (halign, valign) = if h72 == 4 { (1, 2) } else { (h72, v) };
            let generation = flag(e, 71);
            let (sn, cs) = r.sin_cos();
            let mut lx = [cs, sn, 0.0];
            let mut ly = [-sn, cs, 0.0];
            if generation & 2 != 0 {
                lx = [-lx[0], -lx[1], 0.0]; // backward
            }
            if generation & 4 != 0 {
                ly = [-ly[0], -ly[1], 0.0]; // upside down
            }
            let frame = m.then(&Affine::ocs(extrusion(e)?));
            at = frame.point(anchor);
            end = frame.point(second);
            x = frame.vector(lx);
            up = frame.vector(ly);
            width = real(e, 41)?.filter(|w| *w > 0.0).unwrap_or(1.0);
            box_width = 0.0;
        }
        // The frame as seen from above: the baseline direction and the part of "up" square to it.
        let x2: P2 = [x[0], x[1]];
        let u2: P2 = [up[0], up[1]];
        let lx = fit::len(x2);
        let scale = norm(x).max(norm(up));
        if scale <= 0.0 || lx <= 1e-9 * scale {
            self.notes.text_edge_on += 1;
            return Ok(());
        }
        let xh = [x2[0] / lx, x2[1] / lx];
        let square = fit::cross(xh, u2);
        if square.abs() <= 1e-9 * scale {
            self.notes.text_edge_on += 1;
            return Ok(());
        }
        if fit::dot(xh, u2).abs() > 1e-6 * square.abs() {
            self.notes.text_skewed += 1;
        }
        let mirrored = square < 0.0;
        let mut rotation = xh[1].atan2(xh[0]);
        let (mut at2, mut end2): (P2, P2) = ([at[0], at[1]], [end[0], end[1]]);
        if mirrored {
            // Readable over the same place: the baseline turned round, so the text's far end is where its start was.
            self.notes.text_mirrored += 1;
            rotation += std::f64::consts::PI;
            halign = match halign { 0 => 2, 2 => 0, h => h };
            if halign == 3 || halign == 5 {
                std::mem::swap(&mut at2, &mut end2);
            }
        }
        width *= lx / square.abs();
        box_width *= lx;
        let rotation = rotation.rem_euclid(std::f64::consts::TAU);
        let rotation = if rotation >= std::f64::consts::TAU - 1e-12 { 0.0 } else { rotation };
        let font = d.style_font(&style);
        d.texts.push(Text {
            kind, text: value, at: at2, end: end2, height: height * square.abs(), rotation, width, halign, valign,
            box_width, style, font, mirrored, layer, color, group,
        });
        Ok(())
    }
}

struct Emit<'a> {
    d: &'a mut Drawing,
    layer: usize,
    color: i32,
    group: usize,
    notes: &'a mut Notes,
}

impl Emit<'_> {
    fn push(&mut self, geometry: Geometry) {
        self.d.curves.push(Curve { geometry, layer: self.layer, color: self.color, group: self.group });
    }

    fn line(&mut self, a: P3, b: P3) {
        if (a[2] - b[2]).abs() > 1e-9 * (1.0 + norm(a).max(norm(b))) {
            self.notes.off_plane += 1;
        }
        let (a2, b2) = ([a[0], a[1]], [b[0], b[1]]);
        if fit::dist(a2, b2) <= 1e-12 * (1.0 + fit::len(a2).max(fit::len(b2))) {
            self.notes.zero += 1;
            return;
        }
        self.push(Geometry::Line { a: a2, b: b2 });
    }

    /// The curve c + u cos t + v sin t, t0 <= t <= t1, as seen from above: a circle, an arc or an elliptic arc.
    fn conic(&mut self, c: P3, u: P3, v: P3, t0: f64, t1: f64) {
        let (c2, u2, v2): (P2, P2, P2) = ([c[0], c[1]], [u[0], u[1]], [v[0], v[1]]);
        let size = norm(u).max(norm(v));
        if size <= 0.0 {
            self.notes.zero += 1;
            return;
        }
        let area = fit::cross(u2, v2);
        if area.abs() <= 1e-9 * norm(u) * norm(v) {
            self.notes.edge_on += 1;
            return;
        }
        if u[2].abs() + v[2].abs() > 1e-9 * size {
            self.notes.off_plane += 1;
        }
        let full = t1 - t0 >= std::f64::consts::TAU - 1e-12;
        let (lu, lv) = (fit::len(u2), fit::len(v2));
        let circular = (lu - lv).abs() <= 1e-9 * size && fit::dot(u2, v2).abs() <= 1e-9 * size * size;
        if circular && full {
            self.push(Geometry::Circle { c: c2, r: lu });
        } else if circular {
            let p0 = [u2[0] * t0.cos() + v2[0] * t0.sin(), u2[1] * t0.cos() + v2[1] * t0.sin()];
            let start = p0[1].atan2(p0[0]);
            self.push(Geometry::Arc { c: c2, r: lu, start, sweep: (t1 - t0) * area.signum() });
        } else {
            self.push(Geometry::Ellipse { c: c2, u: u2, v: v2, t0, t1 });
        }
    }

    /// Vertices in object coordinates, each with the bulge of the segment that leaves it.
    fn polyline(&mut self, ocs: &Affine, vertices: &[(P3, f64)], closed: bool) {
        let n = vertices.len();
        let segments = if closed { n } else { n.saturating_sub(1) };
        for i in 0..segments {
            let (p0, bulge) = vertices[i];
            let p1 = vertices[(i + 1) % n].0;
            let chord = ((p1[0] - p0[0]).powi(2) + (p1[1] - p0[1]).powi(2)).sqrt();
            if bulge.abs() < 1e-12 || chord <= 0.0 {
                if chord <= 0.0 && (p1[2] - p0[2]).abs() <= 0.0 {
                    continue; // repeated vertex
                }
                self.line(ocs.point(p0), ocs.point(p1));
                continue;
            }
            // Bulge b = tan(sweep / 4); the centre sits (chord / 2)(1 - b^2) / (2 b) left of the chord's midpoint.
            let sweep = 4.0 * bulge.atan();
            let (dx, dy) = ((p1[0] - p0[0]) / chord, (p1[1] - p0[1]) / chord);
            let offset = chord / 2.0 * (1.0 - bulge * bulge) / (2.0 * bulge);
            let centre = [(p0[0] + p1[0]) / 2.0 - dy * offset, (p0[1] + p1[1]) / 2.0 + dx * offset, p0[2]];
            let r = chord / 2.0 * (1.0 + bulge * bulge) / (2.0 * bulge.abs());
            let a0 = (p0[1] - centre[1]).atan2(p0[0] - centre[0]);
            let (v, t0) = if sweep > 0.0 { ([0.0, r, 0.0], a0) } else { ([0.0, -r, 0.0], -a0) };
            self.conic(ocs.point(centre), ocs.vector([r, 0.0, 0.0]), ocs.vector(v), t0, t0 + sweep.abs());
        }
    }

    fn spline(&mut self, e: &Entity, m: &Affine) -> Result<(), KernelError> {
        let degree = flag(e, 71);
        let control = points(e, 10)?;
        let fit_points = points(e, 11)?;
        let knots = reals(e, 40)?;
        let mut weights = reals(e, 41)?;
        let project = |p: P3| { let w = m.point(p); [w[0], w[1]] };
        let lifted = |ps: &[P3]| ps.iter().filter(|p| m.point(**p)[2].abs() > 1e-9 * (1.0 + norm(**p))).count();
        if !(1..=25).contains(&degree) {
            self.d.warn(format!("a SPLINE of degree {degree} was left out"));
            return Ok(());
        }
        let degree = degree as usize;
        let points2: Vec<P2>;
        if control.len() > degree {
            if knots.len() != control.len() + degree + 1 {
                self.d.warn(format!(
                    "a SPLINE has {} knots for {} control points of degree {degree}; it was left out",
                    knots.len(), control.len()));
                return Ok(());
            }
            if lifted(&control) > 0 { self.notes.off_plane += 1; }
            points2 = control.iter().map(|p| project(*p)).collect();
            if weights.len() != control.len() || weights.iter().all(|w| (*w - 1.0).abs() < 1e-15) {
                weights.clear();
            } else if weights.iter().any(|w| !(*w > 0.0)) {
                self.d.warn("a rational SPLINE with a weight that is not positive was left out".into());
                return Ok(());
            }
        } else if fit_points.len() >= 2 {
            if lifted(&fit_points) > 0 { self.notes.off_plane += 1; }
            let through: Vec<P2> = fit_points.iter().map(|p| project(*p)).collect();
            let tangent = |code| -> Result<Option<P2>, KernelError> {
                Ok(point(e, code)?.map(|t| { let w = m.vector(t); [w[0], w[1]] }))
            };
            let ends = match (tangent(12)?, tangent(13)?) {
                (Some(a), Some(b)) => Some((a, b)),
                _ => None,
            };
            let Some(spline) = fit::interpolate(&through, degree.min(3), ends) else {
                self.d.warn("a SPLINE given only by fit points could not be built (repeated points); it was left out".into());
                return Ok(());
            };
            self.push(spline);
            return Ok(());
        } else {
            self.notes.zero += 1;
            return Ok(());
        }
        if knots.windows(2).any(|w| w[1] < w[0]) || knots[degree] >= knots[knots.len() - degree - 1] {
            self.d.warn("a SPLINE with an invalid knot vector was left out".into());
            return Ok(());
        }
        self.push(Geometry::Spline { degree, knots, points: points2, weights });
        Ok(())
    }
}

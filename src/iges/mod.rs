//! Independent IGES 5.3 manifold-solid intake. Explicit topology only: no sewing,
//! healing, implicit units or conversion of loose sheets into solids.
mod geometry;
mod syntax;
pub mod wire;
use crate::{
    geom::{Curve, Surface, Vec3},
    mass,
    neutral::{self, Body, Edge, Face},
    topology, KernelError,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Read,
    path::Path,
};
type Result<T> = std::result::Result<T, KernelError>;
fn error(s: impl Into<String>) -> KernelError {
    KernelError::InvalidInput(format!("IGES: {}", s.into()))
}
fn bad<T>(s: impl Into<String>) -> Result<T> {
    Err(error(s))
}
fn unsupported<T>(s: impl Into<String>) -> Result<T> {
    Err(KernelError::Unsupported(format!("IGES: {}", s.into())))
}
#[derive(Clone, Debug)]
pub struct ImportedBody {
    pub entity_id: usize,
    pub name: String,
    pub source_length_to_inches: f64,
    pub body: Body,
    pub volume: f64,
}
struct Reader<'a> {
    v: &'a [syntax::Value],
    at: usize,
}
impl<'a> Reader<'a> {
    fn next(&mut self) -> Result<f64> {
        let v = self
            .v
            .get(self.at)
            .ok_or_else(|| error("truncated entity parameters"))?
            .number()?;
        self.at += 1;
        Ok(v)
    }
    fn count(&mut self) -> Result<usize> {
        let n = self.next()?;
        if n < 0. || n.fract() != 0. || n > 1_000_000. {
            return bad("invalid integer parameter");
        }
        Ok(n as usize)
    }
    fn flag(&mut self) -> Result<bool> {
        match self.next()? {
            0. => Ok(false),
            1. => Ok(true),
            _ => bad("invalid logical parameter"),
        }
    }
    fn point(&mut self) -> Result<Vec3> {
        Ok(Vec3::new(self.next()?, self.next()?, self.next()?))
    }
    fn finish(&mut self) -> Result<()> {
        if self.at == self.v.len() {
            return Ok(());
        }
        // The optional associativity/property lists are explicitly empty.
        if self.v.len() - self.at == 2 && self.count()? == 0 && self.count()? == 0 {
            return Ok(());
        }
        unsupported("nonempty associativity/property lists or trailing parameters")
    }
}
struct Document {
    file: syntax::File,
    tolerance: f64,
    used: std::cell::RefCell<BTreeSet<usize>>,
}
impl Document {
    fn entity(&self, id: usize) -> Result<&syntax::Entity> {
        self.used.borrow_mut().insert(id);
        self.file
            .entities
            .get(&id)
            .ok_or_else(|| error(format!("missing directory entry {id}")))
    }
    fn reader(&self, id: usize, kind: usize) -> Result<Reader<'_>> {
        let e = self.entity(id)?;
        if e.kind != kind {
            return bad(format!(
                "directory {id} has type {}, expected {kind}",
                e.kind
            ));
        }
        Ok(Reader {
            v: &e.values,
            at: 0,
        })
    }
    fn topology_entity(&self, id: usize, kind: usize) -> Result<Reader<'_>> {
        let e = self.entity(id)?;
        if e.form != 1 || e.transform != 0 {
            return unsupported(format!(
                "type {kind} requires form 1 and no topology-level transform"
            ));
        }
        self.reader(id, kind)
    }
    fn body(&self, id: usize) -> Result<ImportedBody> {
        let entity = self.entity(id)?;
        if entity.form != 0 {
            return unsupported("manifold solid form must be zero");
        }
        let mut solid = self.reader(id, 186)?;
        let shell = solid.count()?;
        let shell_sense = solid.flag()?;
        if solid.count()? != 0 {
            return unsupported("void shells");
        }
        solid.finish()?;
        let mut r = self.topology_entity(shell, 514)?;
        let count = r.count()?;
        if count == 0 || count > 100000 {
            return bad("invalid shell face count");
        }
        let mut body = Body::default();
        let mut vertices = BTreeMap::new();
        let mut edges = BTreeMap::new();
        let mut faces = BTreeSet::new();
        for _ in 0..count {
            let id = r.count()?;
            let same = r.flag()? == shell_sense;
            if !faces.insert(id) {
                return bad("repeated face in shell");
            }
            let mut f = self.topology_entity(id, 510)?;
            let support = f.count()?;
            let n = f.count()?;
            if n == 0 || n > 2048 || !f.flag()? {
                return unsupported("face requires an explicit outer loop and at most 2048 loops");
            }
            let (surface, uv_domain) = self.surface(support)?;
            let mut face = Face {
                surface,
                sense: same,
                loops: Vec::new(),
            };
            let mut pcurves = Vec::new();
            let mut any_uv = false;
            let mut all_uv = true;
            for _ in 0..n {
                let loop_id = f.count()?;
                let mut lp = self.topology_entity(loop_id, 508)?;
                let n = lp.count()?;
                if n == 0 || n > 2048 {
                    return bad("invalid loop length");
                }
                let mut loop_edges = Vec::new();
                let mut loop_uv = Vec::new();
                for _ in 0..n {
                    if lp.count()? != 0 {
                        return unsupported("vertex-only loop use/pole");
                    }
                    let list = lp.count()?;
                    let index = lp.count()?;
                    let forward = lp.flag()?;
                    let key = (list, index);
                    let edge = if let Some(&e) = edges.get(&key) {
                        e
                    } else {
                        let mut er = self.topology_entity(list, 504)?;
                        let count = er.count()?;
                        if index == 0 || index > count || count > 100000 {
                            return bad("invalid edge list index");
                        }
                        let mut selected = None;
                        for i in 1..=count {
                            let curve = er.count()?;
                            let sv = er.count()?;
                            let si = er.count()?;
                            let tv = er.count()?;
                            let ti = er.count()?;
                            if i == index {
                                selected = Some((curve, sv, si, tv, ti));
                            }
                        }
                        er.finish()?;
                        let (curve, sv, si, tv, ti) = selected.unwrap();
                        let mut ends = [0; 2];
                        for (i, key) in [(sv, si), (tv, ti)].into_iter().enumerate() {
                            ends[i] = if let Some(&v) = vertices.get(&key) {
                                v
                            } else {
                                let mut vr = self.topology_entity(key.0, 502)?;
                                let count = vr.count()?;
                                if key.1 == 0 || key.1 > count || count > 100000 {
                                    return bad("invalid vertex list index");
                                }
                                let mut point = None;
                                for j in 1..=count {
                                    let p = vr.point()?;
                                    if j == key.1 {
                                        point = Some(p);
                                    }
                                }
                                vr.finish()?;
                                let index = body.vertices.len();
                                body.vertices.push(point.unwrap().scale(self.file.scale));
                                vertices.insert(key, index);
                                index
                            };
                        }
                        let curve = self.curve(curve, false)?;
                        if curve.start().distance(body.vertices[ends[0]]) > self.tolerance
                            || curve.end().distance(body.vertices[ends[1]]) > self.tolerance
                        {
                            return unsupported(
                                "edge vertices do not match the bounded curve endpoints",
                            );
                        }
                        let e = body.edges.len();
                        body.edges.push(Edge {
                            curve,
                            start: ends[0],
                            end: ends[1],
                        });
                        edges.insert(key, e);
                        e
                    };
                    loop_edges.push((edge, forward));
                    match lp.count()? {
                        0 => {
                            all_uv = false;
                        }
                        1 => {
                            lp.flag()?;
                            let curve = self.curve(lp.count()?, true)?;
                            let Some([u0, u1, v0, v1]) = uv_domain else {
                                return unsupported("parameter curves on analytic IGES supports");
                            };
                            let (p, w) = crate::spline_curve::controls(&curve)?;
                            if p.iter().any(|p| p.z != 0.) {
                                return bad("parameter curve is not two dimensional");
                            }
                            loop_uv.push(Curve::RationalBezier {
                                p: p.map(|p| {
                                    Vec3::new((p.x - u0) / (u1 - u0), (p.y - v0) / (v1 - v0), 0.)
                                }),
                                w,
                            });
                            any_uv = true;
                        }
                        _ => return unsupported("multiple parameter curves per edge use"),
                    }
                }
                lp.finish()?;
                // IGES loop direction is relative to the underlying surface;
                // face-use and shell-use orientations act on the whole loop.
                if !same {
                    loop_edges.reverse();
                    for (_, forward) in &mut loop_edges {
                        *forward = !*forward;
                    }
                    loop_uv.reverse();
                    loop_uv = loop_uv.into_iter().map(|c| c.reversed()).collect();
                }
                face.loops.push(loop_edges);
                pcurves.push(loop_uv);
            }
            f.finish()?;
            if any_uv {
                if !all_uv {
                    return unsupported("partially specified spline pcurve boundary");
                }
                if let Surface::Spline { data } = &mut face.surface {
                    data.pcurves = Some(pcurves);
                }
            }
            body.faces.push(face);
        }
        r.finish()?;
        body = self.motion(entity.transform)?.body(&body);
        let data = neutral::write_body(&body);
        let mut reader = neutral::Reader::with_tolerance(&data, self.tolerance);
        let body = reader.body()?;
        reader.finish()?;
        topology::check(&body, self.tolerance)?;
        let mut reached = BTreeSet::from([0]);
        let mut todo = vec![0];
        let mut owners = vec![Vec::new(); body.edges.len()];
        for (i, f) in body.faces.iter().enumerate() {
            for &(e, _) in f.loops.iter().flatten() {
                owners[e].push(i);
            }
        }
        while let Some(i) = todo.pop() {
            for &(e, _) in body.faces[i].loops.iter().flatten() {
                for &other in &owners[e] {
                    if reached.insert(other) {
                        todo.push(other);
                    }
                }
            }
        }
        if reached.len() != body.faces.len() {
            return bad("disconnected closed shell");
        }
        for f in &body.faces {
            if let Surface::Plane { n, .. } = f.surface {
                for (i, lp) in f.loops.iter().enumerate() {
                    let curves = lp
                        .iter()
                        .map(|&(e, forward)| {
                            if forward {
                                body.edges[e].curve.clone()
                            } else {
                                body.edges[e].curve.reversed()
                            }
                        })
                        .collect::<Vec<_>>();
                    if crate::extrude::signed_area(&curves, n)
                        * if f.sense { 1. } else { -1. }
                        * if i == 0 { 1. } else { -1. }
                        <= 0.
                    {
                        return bad("planar loop winding disagrees with face orientation");
                    }
                }
            }
        }
        let volume = mass::measure(&body, self.tolerance)?.volume;
        if volume <= 0. {
            return bad("nonpositive solid volume");
        }
        Ok(ImportedBody {
            entity_id: id,
            name: entity.label.clone(),
            source_length_to_inches: self.file.scale,
            body,
            volume,
        })
    }
}
pub fn read(bytes: &[u8], tolerance: f64) -> Result<Vec<ImportedBody>> {
    if !tolerance.is_finite() || tolerance <= 0. {
        return bad("invalid tolerance");
    }
    let file = syntax::parse(bytes)?;
    for (&id, e) in &file.entities {
        if !matches!(
            e.kind,
            100 | 110 | 116 | 123 | 124 | 126 | 128 | 186 | 190 | 192 | 502 | 504 | 508 | 510 | 514
        ) {
            return unsupported(format!("directory {id}, entity type {}", e.kind));
        }
    }
    let doc = Document {
        file,
        tolerance,
        used: Default::default(),
    };
    let ids = doc
        .file
        .entities
        .iter()
        .filter(|(_, e)| e.kind == 186)
        .map(|(&id, _)| id)
        .collect::<Vec<_>>();
    if ids.is_empty() || ids.len() > 1024 {
        return unsupported(
            "requires 1..1024 explicit manifold solids; loose surfaces are not sewn",
        );
    }
    let bodies = ids
        .into_iter()
        .map(|id| doc.body(id))
        .collect::<Result<Vec<_>>>()?;
    for id in doc.file.entities.keys() {
        if !doc.used.borrow().contains(id) {
            return unsupported(format!(
                "unreferenced directory entity {id}; partial geometry import is not allowed"
            ));
        }
    }
    Ok(bodies)
}
pub fn read_file(path: &Path, tolerance: f64) -> Result<Vec<ImportedBody>> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)
        .map_err(|e| error(e.to_string()))?
        .take(syntax::MAX_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| error(e.to_string()))?;
    read(&bytes, tolerance)
}

//! Wire protocol v1 between GCmaker and gckernel-worker (agreed with Astra, battle plan 2026-10-27).
//! All little-endian.
//!
//! Request: magic "GCKQ", u32 version (1), u32 op, u32 units (1 = inch), f64 linear tolerance,
//!          u32 count, then count f64.
//! Reply:   magic "GCKR", u32 version (1), u32 status (0 ok, 1 unsupported, 2 invalid input, 3 internal),
//!          u32 count, then count f64 when ok, else count UTF-8 message bytes.
//!
//! Ops: 1 check (a body: validated and echoed back, followed by its V, E, F, L counts),
//!      2 extrude (a profile -> the prism body, followed by its exact volume),
//!      3 cut (two bodies back to back -> body less tool),
//!      4 measure (a body -> its volume, surface area, and 1 when both are closed form, 0 when by quadrature),
//!      5 rebuild prism (a body -> the kernel's own extrusion of it, when it is a straight prism; rebuild.rs),
//!      6 tessellate (a chord tolerance then a body -> its display mesh and edge polylines; tessellate.rs, the
//!        contract with Astra's viewport),
//!      7 transform (one rigid motion and a batch of items -> the items moved; frame.rs, `transform`),
//!      8 revolve (axis point, unit axis, unit seam direction, n, then n (radius, height) points -> the solid of
//!        revolution, followed by its exact volume; revolve.rs),
//!      9 general revolve (axis point, unit axis, unit seam direction, signed angle, loop count, then per loop a segment
//!        count and segments [kind, r, z, ...]: 0 line, 1 arc (centre r, z, direction +1/-1), 2 ellipse arc (centre, u, v,
//!        direction), 3 cubic Bezier (two control points) -> the solid followed by its exact volume; revolve.rs),
//!     10 named extrude (op 2's request -> op 2's reply, then the provenance envelope; provenance.rs),
//!     11 named revolve (op 9's request -> op 9's reply, then the provenance envelope; provenance.rs).

use crate::frame::{Frame, Point2, Point3, Rigid};
use crate::geom::Vec3;
use crate::neutral::{put_curve, put_surface, write_body, Reader};
use crate::{boolean, extrude, mass, provenance, rebuild, revolve, tessellate, topology, KernelError};

pub const VERSION: u32 = 1;
pub const UNITS_INCH: u32 = 1;

pub const OP_CHECK: u32 = 1;
pub const OP_EXTRUDE: u32 = 2;
pub const OP_CUT: u32 = 3;
pub const OP_MEASURE: u32 = 4;
pub const OP_REBUILD_PRISM: u32 = 5;
pub const OP_TESSELLATE: u32 = 6;
pub const OP_TRANSFORM: u32 = 7;
pub const OP_REVOLVE: u32 = 8;
pub const OP_REVOLVE_GENERAL: u32 = 9;
pub const OP_EXTRUDE_NAMED: u32 = 10;
pub const OP_REVOLVE_NAMED: u32 = 11;

pub struct Request {
    pub op: u32,
    pub units: u32,
    pub tolerance: f64,
    pub data: Vec<f64>,
}

fn u32_at(bytes: &[u8], at: usize) -> Result<u32, KernelError> {
    bytes.get(at..at + 4).map(|b| u32::from_le_bytes(b.try_into().unwrap()))
        .ok_or_else(|| KernelError::InvalidInput("request ended early".into()))
}

pub fn parse_request(bytes: &[u8]) -> Result<Request, KernelError> {
    if bytes.get(0..4) != Some(b"GCKQ") { return Err(KernelError::InvalidInput("not a GC Kernel request".into())); }
    let version = u32_at(bytes, 4)?;
    if version != VERSION { return Err(KernelError::Unsupported(format!("request version {version}"))); }
    let op = u32_at(bytes, 8)?;
    let units = u32_at(bytes, 12)?;
    let tolerance = bytes.get(16..24).map(|b| f64::from_le_bytes(b.try_into().unwrap()))
        .ok_or_else(|| KernelError::InvalidInput("request ended early".into()))?;
    let count = u32_at(bytes, 24)? as usize;
    let payload = &bytes[28..];
    if payload.len() != count * 8 {
        return Err(KernelError::InvalidInput(format!("request says {count} values but carries {} bytes", payload.len())));
    }
    let data = payload.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
    Ok(Request { op, units, tolerance, data })
}

pub fn encode_request(op: u32, tolerance: f64, data: &[f64]) -> Vec<u8> {
    let mut out = b"GCKQ".to_vec();
    for v in [VERSION, op, UNITS_INCH] { out.extend_from_slice(&v.to_le_bytes()); }
    out.extend_from_slice(&tolerance.to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    for v in data { out.extend_from_slice(&v.to_le_bytes()); }
    out
}

/// Runs one request.
pub fn handle(request: &Request) -> Result<Vec<f64>, KernelError> {
    if request.units != UNITS_INCH {
        return Err(KernelError::Unsupported(format!("units code {} (v1 is inches; convert at capture)", request.units)));
    }
    let tol = request.tolerance;
    if !tol.is_finite() || tol <= 0.0 { return Err(KernelError::InvalidInput("the tolerance must be positive".into())); }
    let mut reader = Reader::with_tolerance(&request.data, tol);
    match request.op {
        OP_CHECK => {
            let body = reader.body()?;
            reader.finish()?;
            let report = topology::check(&body, tol)?;
            let mut out = write_body(&body);
            out.extend_from_slice(&[report.vertices as f64, report.edges as f64, report.faces as f64, report.loops as f64]);
            Ok(out)
        }
        OP_EXTRUDE => {
            let profile = reader.profile()?;
            reader.finish()?;
            let body = extrude::extrude(&profile, tol)?;
            topology::check(&body, tol).map_err(|e| KernelError::Internal(format!("extrude built an invalid solid: {e}")))?;
            let mut out = write_body(&body);
            out.push(extrude::volume(&profile));
            Ok(out)
        }
        OP_EXTRUDE_NAMED => {
            let profile = reader.profile()?;
            reader.finish()?;
            let (body, prov) = extrude::extrude_named(&profile, tol)?;
            topology::check(&body, tol).map_err(|e| KernelError::Internal(format!("extrude built an invalid solid: {e}")))?;
            let mut out = write_body(&body);
            out.push(extrude::volume(&profile));
            provenance::write(&mut out, &prov);
            Ok(out)
        }
        OP_CUT => {
            let body = reader.body()?;
            let tool = reader.body()?;
            reader.finish()?;
            topology::check(&body, tol)?;
            topology::check(&tool, tol)?;
            let result = boolean::cut(&body, &tool, tol)?;
            // Complete removal is an explicit, successful empty result (Astra's acceptance cases).
            if !result.faces.is_empty() {
                topology::check(&result, tol).map_err(|e| KernelError::Internal(format!("cut built an invalid solid: {e}")))?;
            }
            Ok(write_body(&result))
        }
        OP_MEASURE => {
            let body = reader.body()?;
            reader.finish()?;
            topology::check(&body, tol)?;
            let m = mass::measure(&body, tol)?;
            Ok(vec![m.volume, m.area, if m.exact { 1.0 } else { 0.0 }])
        }
        OP_REBUILD_PRISM => {
            let body = reader.body()?;
            reader.finish()?;
            topology::check(&body, tol)?;
            let rebuilt = rebuild::rebuild_prism(&body, tol)?;
            topology::check(&rebuilt, tol).map_err(|e| KernelError::Internal(format!("rebuild built an invalid solid: {e}")))?;
            Ok(write_body(&rebuilt))
        }
        OP_TESSELLATE => {
            let chord = reader.next()?;
            let body = reader.body()?;
            reader.finish()?;
            topology::check(&body, tol)?;
            Ok(tessellate::write_mesh(&tessellate::tessellate(&body, chord)?))
        }
        OP_TRANSFORM => transform(&request.data, tol),
        OP_REVOLVE => {
            let (p, a, x) = (reader.point()?, reader.point()?, reader.point()?);
            let n = reader.count()?;
            let mut points = Vec::with_capacity(n);
            for _ in 0..n { points.push((reader.next()?, reader.next()?)); }
            reader.finish()?;
            let section = revolve::Section { p, a, x, points };
            let body = revolve::revolve(&section, tol)?;
            topology::check(&body, tol).map_err(|e| KernelError::Internal(format!("revolve built an invalid solid: {e}")))?;
            let mut out = write_body(&body);
            out.push(revolve::volume(&section));
            Ok(out)
        }
        OP_REVOLVE_GENERAL | OP_REVOLVE_NAMED => {
            let (p, a, x) = (reader.point()?, reader.point()?, reader.point()?);
            let angle = reader.next()?;
            let loop_count = reader.count()?;
            let mut loops = Vec::with_capacity(loop_count);
            for _ in 0..loop_count {
                let n = reader.count()?;
                let mut segments = Vec::with_capacity(n);
                for _ in 0..n {
                    let kind = reader.next()?;
                    let (r, z) = (reader.next()?, reader.next()?);
                    let direction = |v: f64| if v == 1.0 { Ok(true) } else if v == -1.0 { Ok(false) }
                        else { Err(KernelError::InvalidInput(format!("a segment direction must be 1 or -1, not {v}"))) };
                    let kind = match kind {
                        k if k == 0.0 => revolve::Kind::Line,
                        k if k == 1.0 => {
                            let (cr, cz, d) = (reader.next()?, reader.next()?, reader.next()?);
                            revolve::Kind::Arc { cr, cz, ccw: direction(d)? }
                        }
                        k if k == 2.0 => {
                            let v: Vec<f64> = (0..7).map(|_| reader.next()).collect::<Result<_, _>>()?;
                            revolve::Kind::Ellipse { cr: v[0], cz: v[1], ur: v[2], uz: v[3], vr: v[4], vz: v[5], forward: direction(v[6])? }
                        }
                        k if k == 3.0 => {
                            let v: Vec<f64> = (0..4).map(|_| reader.next()).collect::<Result<_, _>>()?;
                            revolve::Kind::Bezier { c1: (v[0], v[1]), c2: (v[2], v[3]) }
                        }
                        k => return Err(KernelError::InvalidInput(format!("unknown section segment kind {k}"))),
                    };
                    segments.push(revolve::Segment { r, z, kind });
                }
                loops.push(segments);
            }
            reader.finish()?;
            let rev = revolve::Revolution { p, a, x, angle, loops };
            let (body, prov) = revolve::revolve_general_named(&rev, tol)?;
            topology::check(&body, tol).map_err(|e| KernelError::Internal(format!("revolve built an invalid solid: {e}")))?;
            let mut out = write_body(&body);
            out.push(revolve::volume_of(&rev, tol)?);
            if request.op == OP_REVOLVE_NAMED { provenance::write(&mut out, &prov); }
            Ok(out)
        }
        op => Err(KernelError::Unsupported(format!("op {op}"))),
    }
}

pub fn encode_reply(result: &Result<Vec<f64>, KernelError>) -> Vec<u8> {
    let mut out = b"GCKR".to_vec();
    out.extend_from_slice(&VERSION.to_le_bytes());
    match result {
        Ok(values) => {
            out.extend_from_slice(&0u32.to_le_bytes());
            out.extend_from_slice(&(values.len() as u32).to_le_bytes());
            for v in values { out.extend_from_slice(&v.to_le_bytes()); }
        }
        Err(e) => {
            let (status, message) = match e {
                KernelError::Unsupported(m) => (1u32, m),
                KernelError::InvalidInput(m) => (2u32, m),
                KernelError::Internal(m) => (3u32, m),
            };
            out.extend_from_slice(&status.to_le_bytes());
            out.extend_from_slice(&(message.len() as u32).to_le_bytes());
            out.extend_from_slice(message.as_bytes());
        }
    }
    out
}

/// Op 7 TRANSFORM (contract with Astra, GCm Handoff/Claude-GcKernel-FilletCuts-Request.md, 2026-10-28).
/// Request: [1, mode (0 apply T, 1 apply T inverse), T (m row-major 9, t 3), N, N items of [kind, len, data[len]]].
/// Kinds: 0 point (3), 1 vector (3), 2 frame (o x y, 9), 3 rigid (12; reply: T applied after it), 4 curve (neutral),
/// 5 body (neutral v1), 6 frame-local point (frame 9, x, y; reply kind 0), 7 surface (neutral).
/// Reply: [1, N, N items of [kind, len, data]]. Any bad item fails the whole batch with its index.
pub fn transform(data: &[f64], tol: f64) -> Result<Vec<f64>, KernelError> {
    let bad = |m: String| KernelError::InvalidInput(m);
    let mut at = 0usize;
    fn take(data: &[f64], at: &mut usize, what: &str) -> Result<f64, KernelError> {
        let v = *data.get(*at).ok_or_else(|| KernelError::InvalidInput(format!("transform data ended early ({what})")))?;
        *at += 1;
        if v.is_finite() { Ok(v) } else { Err(KernelError::InvalidInput(format!("transform data is not finite ({what})"))) }
    }
    let next = |at: &mut usize, what: &str| take(data, at, what);
    if next(&mut at, "version")? != 1.0 { return Err(bad("unsupported transform version".into())); }
    let mode = next(&mut at, "mode")?;
    if mode != 0.0 && mode != 1.0 { return Err(bad(format!("transform mode must be 0 or 1, not {mode}"))); }
    let mut t = [0.0; 12];
    for v in t.iter_mut() { *v = next(&mut at, "transform")?; }
    let mut motion = rigid(&t).map_err(|e| bad(format!("the transform: {e}")))?;
    if mode == 1.0 { motion = motion.inverse(); }
    let count = next(&mut at, "item count")?;
    if count < 0.0 || count != count.floor() || count > 1e7 { return Err(bad(format!("bad item count {count}"))); }
    let mut out = vec![1.0, count];
    for index in 0..count as usize {
        let kind = next(&mut at, "item kind")?;
        let len = next(&mut at, "item length")?;
        if len < 0.0 || len != len.floor() { return Err(bad(format!("item {index}: bad length {len}"))); }
        let len = len as usize;
        let end = at.checked_add(len).filter(|&e| e <= data.len())
            .ok_or_else(|| bad(format!("item {index}: its length {len} runs past the data")))?;
        let item = &data[at..end];
        at = end;
        let (kind, reply) = transform_item(&motion, kind, item, tol).map_err(|e| bad(format!("item {index}: {e}")))?;
        out.push(kind);
        out.push(reply.len() as f64);
        out.extend_from_slice(&reply);
    }
    if at != data.len() { return Err(bad(format!("transform data has {} values left over", data.len() - at))); }
    Ok(out)
}

fn rigid(v: &[f64]) -> Result<Rigid, KernelError> {
    Rigid::from_rows([[v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]], Vec3::new(v[9], v[10], v[11]))
}

fn frame(v: &[f64]) -> Result<Frame, KernelError> {
    Frame::new(Point3::new(v[0], v[1], v[2]), Vec3::new(v[3], v[4], v[5]), Vec3::new(v[6], v[7], v[8]))
}

fn transform_item(motion: &Rigid, kind: f64, item: &[f64], tol: f64) -> Result<(f64, Vec<f64>), KernelError> {
    let exact = |n: usize| -> Result<(), KernelError> {
        if item.len() == n { Ok(()) } else { Err(KernelError::InvalidInput(format!("kind {kind} needs {n} values, not {}", item.len()))) }
    };
    if item.iter().any(|v| !v.is_finite()) { return Err(KernelError::InvalidInput("not finite".into())); }
    let p3 = |p: Point3| vec![p.x, p.y, p.z];
    let v3 = |v: Vec3| vec![v.x, v.y, v.z];
    match kind {
        k if k == 0.0 => { exact(3)?; Ok((0.0, p3(motion.point(Point3::new(item[0], item[1], item[2]))))) }
        k if k == 1.0 => { exact(3)?; Ok((1.0, v3(motion.vector(Vec3::new(item[0], item[1], item[2]))))) }
        k if k == 2.0 => {
            exact(9)?;
            let f = motion.frame(&frame(item)?);
            let mut o = p3(f.origin());
            o.extend(v3(f.x()));
            o.extend(v3(f.y()));
            Ok((2.0, o))
        }
        k if k == 3.0 => {
            exact(12)?;
            let (m, t) = rigid(item)?.then(motion).rows();
            let mut o: Vec<f64> = m.iter().flatten().copied().collect();
            o.extend(v3(t));
            Ok((3.0, o))
        }
        k if k == 4.0 => {
            let mut r = Reader::with_tolerance(item, tol);
            let c = r.curve()?;
            r.finish()?;
            let mut o = Vec::new();
            put_curve(&mut o, &motion.curve(&c));
            Ok((4.0, o))
        }
        k if k == 5.0 => {
            let mut r = Reader::with_tolerance(item, tol);
            let body = r.body()?;
            r.finish()?;
            topology::check(&body, tol)?;
            let moved = motion.body(&body);
            topology::check(&moved, tol).map_err(|e| KernelError::Internal(format!("the moved body is not valid: {e}")))?;
            Ok((5.0, write_body(&moved)))
        }
        k if k == 6.0 => {
            exact(11)?;
            let p = frame(&item[..9])?.point2(Point2::new(item[9], item[10]));
            Ok((0.0, p3(motion.point(p))))
        }
        k if k == 7.0 => {
            let mut r = Reader::with_tolerance(item, tol);
            let s = r.surface()?;
            r.finish()?;
            let mut o = Vec::new();
            put_surface(&mut o, &motion.surface(&s));
            Ok((7.0, o))
        }
        k => Err(KernelError::InvalidInput(format!("unknown item kind {k}"))),
    }
}

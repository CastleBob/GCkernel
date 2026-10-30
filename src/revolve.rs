//! Solids of revolution (Claude 2026-10-29; Bob: tapered Pole stock, then the whole revolve subsystem). A closed
//! section in the half-plane through an axis, turned about it by a signed angle (a full turn or part of one).
//!
//! Section: one outer loop and any number of holes in (radius, height) coordinates (radius along the unit `x`, height
//! along the unit axis `a`, from the axis point `p`). Segments are lines or circular arcs (ellipse arcs and cubic Bezier
//! curves are part of the op 9 contract and refused by name until their meridian surfaces exist). Exact faces: a line
//! along the axis direction sweeps a cylinder, a sloped line a cone (an apex where it reaches the axis), a line across the
//! axis a plane annulus or disc; an arc centred on the axis sweeps a sphere, an arc centred off it a torus. A part turn
//! adds two planar caps in the meridian planes at 0 and at the angle. Segments lying on the axis bound nothing (in a full
//! turn) or join the two caps (in a part turn).
//!
//! Topology: circle (or arc) edges where faces meet, one vertex per point on the axis. Curved faces touching the axis in a
//! full turn, or in a part turn of more than half a turn, are split at the middle angle, so every face has distinct
//! boundary angles around its pole (a chart for display and mass). Each face uses its first point's circle forward and its
//! second point's reversed: with outer loops counter-clockwise and holes clockwise, material lies left of travel and the
//! outward normal is (dz, -dr), which makes the loops run counter-clockwise about it.

use crate::geom::{Curve, Surface, Vec3};
use crate::meridian::Meridian;
use crate::neutral::{Body, Edge, Face};
use crate::provenance::*;
use crate::KernelError;
use std::f64::consts::{PI, TAU};

/// The op 8 input: axis point `p`, unit axis `a`, unit `x` perpendicular to it (the seam's half-plane), and one loop of
/// (radius, height) points joined by straight lines, revolved a full turn.
#[derive(Clone, Debug)]
pub struct Section {
    pub p: Vec3,
    pub a: Vec3,
    pub x: Vec3,
    pub points: Vec<(f64, f64)>,
}

/// One section segment, from its start point to the next segment's start.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Kind {
    Line,
    /// A circular arc about (cr, cz), counter-clockwise in (radius, height) when `ccw`.
    Arc { cr: f64, cz: f64, ccw: bool },
    /// c + u cos t + v sin t in (radius, height); t increasing when `forward` (R4: refused until built).
    Ellipse { cr: f64, cz: f64, ur: f64, uz: f64, vr: f64, vz: f64, forward: bool },
    /// Cubic Bezier with control points c1, c2 (R4: refused until built).
    Bezier { c1: (f64, f64), c2: (f64, f64) },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Segment {
    pub r: f64,
    pub z: f64,
    pub kind: Kind,
}

/// The op 9 input: axis frame, signed angle (0 < |angle| <= 2 pi; negative turns about -a) and loops of segments.
#[derive(Clone, Debug)]
pub struct Revolution {
    pub p: Vec3,
    pub a: Vec3,
    pub x: Vec3,
    pub angle: f64,
    pub loops: Vec<Vec<Segment>>,
}

impl Revolution {
    pub fn from_section(s: &Section) -> Revolution {
        Revolution { p: s.p, a: s.a, x: s.x, angle: TAU,
            loops: vec![s.points.iter().map(|&(r, z)| Segment { r, z, kind: Kind::Line }).collect()] }
    }
}

fn invalid<T>(m: impl Into<String>) -> Result<T, KernelError> { Err(KernelError::InvalidInput(m.into())) }

/// A segment resolved in (radius, height): start, end and, for an arc, its centre, radius, start angle and signed sweep.
#[derive(Clone, Copy, Debug)]
struct Piece {
    s: (f64, f64),
    e: (f64, f64),
    arc: Option<(f64, f64, f64, f64, f64)>, // cr, cz, rho, alpha0, signed sweep
    /// An ellipse or Bezier meridian with its parameter interval, increasing along travel.
    mer: Option<(Meridian, f64, f64)>,
}

impl Piece {
    fn reversed(&self) -> Piece {
        Piece { s: self.e, e: self.s, arc: self.arc.map(|(cr, cz, rho, a0, sw)| (cr, cz, rho, a0 + sw, -sw)),
            mer: self.mer.map(|(m, t0, t1)| match m {
                Meridian::Bezier(_) => (m.reversed(), 1.0 - t1, 1.0 - t0),
                Meridian::Ellipse { .. } => (m.reversed(), -t1, -t0),
            }) }
    }
    /// The point at fraction f in [0, 1] of the piece.
    fn at(&self, f: f64) -> (f64, f64) {
        if let Some((m, t0, t1)) = self.mer { return m.eval(t0 + (t1 - t0) * f); }
        match self.arc {
            None => (self.s.0 + (self.e.0 - self.s.0) * f, self.s.1 + (self.e.1 - self.s.1) * f),
            Some((cr, cz, rho, a0, sw)) => (cr + rho * (a0 + sw * f).cos(), cz + rho * (a0 + sw * f).sin()),
        }
    }
    /// A box holding the piece over fractions [f0, f1] (exact for lines and arcs, conservative for meridians).
    fn bounds(&self, f0: f64, f1: f64) -> (f64, f64, f64, f64) {
        if let Some((m, t0, t1)) = self.mer { return m.bounds(t0 + (t1 - t0) * f0, t0 + (t1 - t0) * f1); }
        let mut pts = vec![self.at(f0), self.at(f1)];
        if let Some((cr, cz, rho, a0, sw)) = self.arc {
            let (lo, hi) = { let (x, y) = (a0 + sw * f0, a0 + sw * f1); (x.min(y), x.max(y)) };
            for k in -8..=8 {
                let t = k as f64 * PI / 2.0;
                if t > lo && t < hi { pts.push((cr + rho * t.cos(), cz + rho * t.sin())); }
            }
        }
        pts.iter().fold((f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
            |b, p| (b.0.min(p.0), b.1.min(p.1), b.2.max(p.0), b.3.max(p.1)))
    }
    /// Twice the signed area contribution: the loop integral of (r dz - z dr).
    fn area2(&self) -> f64 {
        if let Some((m, t0, t1)) = self.mer { return m.area2(t0, t1); }
        match self.arc {
            None => self.s.0 * self.e.1 - self.e.0 * self.s.1,
            Some((cr, cz, rho, a0, sw)) => {
                let f = |t: f64| rho * (cr * t.sin() - cz * t.cos()) + rho * rho * t;
                f(a0 + sw) - f(a0)
            }
        }
    }
    /// The loop integral of r^2 dz (Pappus: pi times its loop sum is the full-turn volume).
    fn r2dz(&self) -> f64 {
        if let Some((m, t0, t1)) = self.mer { return m.r2dz(t0, t1); }
        match self.arc {
            None => {
                let ((r0, z0), (r1, z1)) = (self.s, self.e);
                (z1 - z0) * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0
            }
            // r = cr + rho cos t, dz = rho cos t dt.
            Some((cr, _, rho, a0, sw)) => {
                let f = |t: f64| rho * (cr * cr * t.sin() + 2.0 * cr * rho * (t / 2.0 + (2.0 * t).sin() / 4.0)
                    + rho * rho * (t.sin() - t.sin().powi(3) / 3.0));
                f(a0 + sw) - f(a0)
            }
        }
    }
    /// The smallest radius strictly inside the segment (an arc passing the angle pi), if any.
    fn inner_min_radius(&self) -> Option<f64> {
        if let Some((m, t0, t1)) = self.mer { return m.inner_min_radius(t0, t1); }
        let (cr, _, rho, a0, sw) = self.arc?;
        let (lo, hi) = if sw > 0.0 { (a0, a0 + sw) } else { (a0 + sw, a0) };
        let t = PI + ((lo - PI) / TAU).ceil() * TAU;
        (t > lo && t < hi).then_some(cr - rho)
    }
    fn mid(&self) -> (f64, f64) {
        if self.mer.is_some() { return self.at(0.5); }
        match self.arc {
            None => ((self.s.0 + self.e.0) / 2.0, (self.s.1 + self.e.1) / 2.0),
            Some((cr, cz, rho, a0, sw)) => (cr + rho * (a0 + sw / 2.0).cos(), cz + rho * (a0 + sw / 2.0).sin()),
        }
    }
}

fn resolve(seg: &Segment, next: &Segment, tol: f64, at: &str) -> Result<Piece, KernelError> {
    let (s, e) = ((seg.r, seg.z), (next.r, next.z));
    if (e.0 - s.0).hypot(e.1 - s.1) <= tol { return invalid(format!("{at} has no length")); }
    match seg.kind {
        Kind::Line => Ok(Piece { s, e, arc: None, mer: None }),
        Kind::Arc { cr, cz, ccw } => {
            if !(cr.is_finite() && cz.is_finite()) { return invalid(format!("{at}: the arc centre is not finite")); }
            let (rho, rho_end) = ((s.0 - cr).hypot(s.1 - cz), (e.0 - cr).hypot(e.1 - cz));
            if rho <= tol || (rho - rho_end).abs() > tol {
                return invalid(format!("{at}: the arc's ends are not on one radius ({rho}, {rho_end})"));
            }
            let (a0, a1) = ((s.1 - cz).atan2(s.0 - cr), (e.1 - cz).atan2(e.0 - cr));
            let sw = if ccw { (a1 - a0).rem_euclid(TAU) } else { -(a0 - a1).rem_euclid(TAU) };
            Ok(Piece { s, e, arc: Some((cr, cz, rho, a0, sw)), mer: None })
        }
        Kind::Ellipse { cr, cz, ur, uz, vr, vz, forward } => {
            let vals = [cr, cz, ur, uz, vr, vz];
            if vals.iter().any(|v| !v.is_finite()) { return invalid(format!("{at}: the ellipse is not finite")); }
            let det = ur * vz - uz * vr;
            if det.abs() <= 1e-12 * ur.hypot(uz) * vr.hypot(vz) || ur.hypot(uz) <= tol || vr.hypot(vz) <= tol {
                return invalid(format!("{at}: the ellipse's semi-diameters must be independent"));
            }
            let m = Meridian::Ellipse { c: (cr, cz), u: (ur, uz), v: (vr, vz) };
            // Travel always runs with increasing parameter: a backwards arc is the reversed ellipse.
            let m = if forward { m } else { m.reversed() };
            let t0 = m.param_of(s, 0.0);
            let mut t1 = m.param_of(e, t0);
            if t1 <= t0 { t1 += TAU; }
            for (q, t) in [(s, t0), (e, t1)] {
                let on = m.eval(t);
                if (on.0 - q.0).hypot(on.1 - q.1) > tol { return invalid(format!("{at}: an end is not on the ellipse")); }
            }
            Ok(Piece { s, e, arc: None, mer: Some((m, t0, t1)) })
        }
        Kind::Bezier { c1, c2 } => {
            if !(c1.0.is_finite() && c1.1.is_finite() && c2.0.is_finite() && c2.1.is_finite()) {
                return invalid(format!("{at}: the Bezier controls are not finite"));
            }
            Ok(Piece { s, e, arc: None, mer: Some((Meridian::Bezier([s, c1, c2, e]), 0.0, 1.0)) })
        }
    }
}

/// Proves two pieces stay farther apart than tol (except at the given shared points: pairs of fractions (on p, on q)
/// where they meet as neighbours), by subdividing until boxes separate (sound: boxes hold the pieces). False when the
/// pieces come within tol (boxes smaller than tol still touching) or the proof does not finish.
fn apart(p: &Piece, q: &Piece, shared: &[(f64, f64)], tol: f64) -> bool {
    let size = { let b = p.bounds(0.0, 1.0); (b.2 - b.0).hypot(b.3 - b.1) }.max({ let b = q.bounds(0.0, 1.0); (b.2 - b.0).hypot(b.3 - b.1) });
    let exclusion = (10.0 * tol).max(1e-9 * size);
    let mut stack = vec![((0.0, 1.0), (0.0, 1.0), 0usize)];
    let diag = |b: (f64, f64, f64, f64)| (b.2 - b.0).hypot(b.3 - b.1);
    while let Some(((a0, a1), (b0, b1), depth)) = stack.pop() {
        let (bp, bq) = (p.bounds(a0, a1), q.bounds(b0, b1));
        let gap = (bq.0 - bp.2).max(bp.0 - bq.2).max(0.0).hypot((bq.1 - bp.3).max(bp.1 - bq.3).max(0.0));
        if gap > tol { continue; }
        let at_shared = shared.iter().any(|&(fp, fq)| fp >= a0 && fp <= a1 && fq >= b0 && fq <= b1);
        if at_shared && diag(bp) <= exclusion && diag(bq) <= exclusion { continue; }
        if depth > 60 || (diag(bp) <= tol && diag(bq) <= tol) { return false; }
        if diag(bp) >= diag(bq) {
            let m = (a0 + a1) / 2.0;
            stack.push(((a0, m), (b0, b1), depth + 1));
            stack.push(((m, a1), (b0, b1), depth + 1));
        } else {
            let m = (b0 + b1) / 2.0;
            stack.push(((a0, a1), (b0, m), depth + 1));
            stack.push(((a0, a1), (m, b1), depth + 1));
        }
        if stack.len() > 200_000 { return false; }
    }
    true
}

/// Is angle t on the arc from a0 through the signed sweep sw?
fn within(t: f64, a0: f64, sw: f64) -> bool {
    let d = if sw >= 0.0 { (t - a0).rem_euclid(TAU) } else { (a0 - t).rem_euclid(TAU) };
    d <= sw.abs()
}

fn on_piece(piece: &Piece, pt: (f64, f64)) -> bool {
    match piece.arc {
        None => {
            let (dx, dy) = (piece.e.0 - piece.s.0, piece.e.1 - piece.s.1);
            let t = ((pt.0 - piece.s.0) * dx + (pt.1 - piece.s.1) * dy) / (dx * dx + dy * dy);
            (0.0..=1.0).contains(&t)
        }
        Some((cr, cz, _, a0, sw)) => within((pt.1 - cz).atan2(pt.0 - cr), a0, sw),
    }
}

/// The point of a piece nearest q (exact: a line's clamped foot; an arc's radial point when its angle is on the arc,
/// otherwise the nearer end).
fn closest_on(piece: &Piece, q: (f64, f64)) -> (f64, f64) {
    match piece.arc {
        None => {
            let (dx, dy) = (piece.e.0 - piece.s.0, piece.e.1 - piece.s.1);
            let t = (((q.0 - piece.s.0) * dx + (q.1 - piece.s.1) * dy) / (dx * dx + dy * dy)).clamp(0.0, 1.0);
            (piece.s.0 + t * dx, piece.s.1 + t * dy)
        }
        Some((cr, cz, rho, a0, sw)) => {
            let t = (q.1 - cz).atan2(q.0 - cr);
            if within(t, a0, sw) && (q.0 - cr).hypot(q.1 - cz) > 0.0 { (cr + rho * t.cos(), cz + rho * t.sin()) }
            else if (q.0 - piece.s.0).hypot(q.1 - piece.s.1) <= (q.0 - piece.e.0).hypot(q.1 - piece.e.1) { piece.s }
            else { piece.e }
        }
    }
}

/// Points where the two pieces' supports cross and both pieces pass (exact line/circle algebra).
fn crossings(p: &Piece, q: &Piece) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    match (p.arc, q.arc) {
        (None, None) => {
            let (d1, d2) = ((p.e.0 - p.s.0, p.e.1 - p.s.1), (q.e.0 - q.s.0, q.e.1 - q.s.1));
            let den = d1.0 * d2.1 - d1.1 * d2.0;
            if den != 0.0 {
                let t = ((q.s.0 - p.s.0) * d2.1 - (q.s.1 - p.s.1) * d2.0) / den;
                out.push((p.s.0 + t * d1.0, p.s.1 + t * d1.1));
            }
        }
        (None, Some((cr, cz, rho, _, _))) | (Some((cr, cz, rho, _, _)), None) => {
            let l = if p.arc.is_none() { p } else { q };
            let (dx, dy) = (l.e.0 - l.s.0, l.e.1 - l.s.1);
            let (fx, fy) = (l.s.0 - cr, l.s.1 - cz);
            let (aa, bb, cc) = (dx * dx + dy * dy, 2.0 * (fx * dx + fy * dy), fx * fx + fy * fy - rho * rho);
            let disc = bb * bb - 4.0 * aa * cc;
            if disc >= 0.0 {
                for sg in [-1.0, 1.0] {
                    let t = (-bb + sg * disc.sqrt()) / (2.0 * aa);
                    out.push((l.s.0 + t * dx, l.s.1 + t * dy));
                }
            }
        }
        (Some((c1r, c1z, r1, _, _)), Some((c2r, c2z, r2, _, _))) => {
            let (dx, dy) = (c2r - c1r, c2z - c1z);
            let d = dx.hypot(dy);
            if d > 0.0 && d <= r1 + r2 && d >= (r1 - r2).abs() {
                let along = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
                let h = (r1 * r1 - along * along).max(0.0).sqrt();
                let (mx, my) = (c1r + along * dx / d, c1z + along * dy / d);
                out.push((mx + h * dy / d, my - h * dx / d));
                out.push((mx - h * dy / d, my + h * dx / d));
            } else if d == 0.0 && r1 == r2 {
                // One circle: shared ends or overlap of the two arcs.
                out.extend([p.s, p.e, q.s, q.e, p.mid(), q.mid()]);
            }
        }
    }
    out.into_iter().filter(|&pt| on_piece(p, pt) && on_piece(q, pt)).collect()
}

/// The least distance between two pieces, ignoring contact at `shared` (a neighbour's common end) when given.
fn distance(p: &Piece, q: &Piece, shared: Option<(f64, f64)>, tol: f64) -> f64 {
    let near = |pt: (f64, f64)| shared.is_some_and(|sh| (pt.0 - sh.0).hypot(pt.1 - sh.1) <= tol);
    let mut best = f64::INFINITY;
    let mut consider = |u: (f64, f64), v: (f64, f64)| {
        if near(u) && near(v) { return; }
        best = best.min((u.0 - v.0).hypot(u.1 - v.1));
    };
    for (s, o) in [(p, q), (q, p)] {
        for end in [s.s, s.e] { consider(end, closest_on(o, end)); }
        if let Some((cr, cz, _, _, _)) = s.arc {
            let f = closest_on(o, (cr, cz));
            consider(closest_on(s, f), f);
            // The far side too: the point of o farthest along the line of centres (arc-arc and arc-line extremes).
            if let Some((c2r, c2z, rho2, _, _)) = o.arc {
                let (dx, dy) = (cr - c2r, cz - c2z);
                let d = dx.hypot(dy);
                if d > 0.0 {
                    for sg in [-1.0, 1.0] {
                        let pt = (c2r + sg * rho2 * dx / d, c2z + sg * rho2 * dy / d);
                        if on_piece(o, pt) { consider(closest_on(s, pt), pt); }
                    }
                }
            }
        }
    }
    for x in crossings(p, q) { consider(x, x); }
    best
}

/// Crossings of the ray from q towards +r with a piece (half-open in height, so a shared end counts once).
fn ray_crossings(pc: &Piece, q: (f64, f64)) -> usize {
    let counts = |e0: f64, e1: f64| (e0 > q.1) != (e1 > q.1);
    if pc.mer.is_some() {
        let n = 256;
        let mut found = 0;
        for k in 0..n {
            let (f0, f1) = (k as f64 / n as f64, (k + 1) as f64 / n as f64);
            let (z0, z1) = (pc.at(f0).1, pc.at(f1).1);
            if !counts(z0, z1) { continue; }
            let (mut lo, mut hi) = (f0, f1);
            for _ in 0..60 {
                let m = (lo + hi) / 2.0;
                if counts(pc.at(lo).1, pc.at(m).1) { hi = m } else { lo = m }
            }
            found += usize::from(pc.at((lo + hi) / 2.0).0 > q.0);
        }
        return found;
    }
    match pc.arc {
        None => {
            let ((r0, z0), (r1, z1)) = (pc.s, pc.e);
            if !counts(z0, z1) { return 0; }
            usize::from(r0 + (q.1 - z0) * (r1 - r0) / (z1 - z0) > q.0)
        }
        Some((cr, cz, rho, a0, sw)) => {
            // Split the arc at its top and bottom so each part is monotone in height, then count each part like a line.
            let mut cuts = vec![a0, a0 + sw];
            for k in -2..=2 {
                for base in [PI / 2.0, -PI / 2.0] {
                    let t = base + k as f64 * TAU;
                    if (t - a0) * (t - a0 - sw) < 0.0 { cuts.push(t); }
                }
            }
            cuts.sort_by(|u, v| if sw > 0.0 { u.total_cmp(v) } else { v.total_cmp(u) });
            let mut n = 0;
            for w in cuts.windows(2) {
                let (z0, z1) = (cz + rho * w[0].sin(), cz + rho * w[1].sin());
                if !counts(z0, z1) { continue; }
                let s = ((q.1 - cz) / rho).clamp(-1.0, 1.0).asin();
                // The crossing on this monotone part: the angle with sin = s between w0 and w1.
                let (lo, hi) = (w[0].min(w[1]), w[0].max(w[1]));
                let cands = [s, PI - s];
                let t = (-3..=3).flat_map(|k| cands.map(|c| c + k as f64 * TAU)).find(|t| *t >= lo - 1e-15 && *t <= hi + 1e-15);
                if let Some(t) = t { n += usize::from(cr + rho * t.cos() > q.0); }
            }
            n
        }
    }
}

/// Exact volume: |angle| / 2 times the loop integral of r^2 dz over every loop (holes subtract by orientation).
pub fn volume_of(rev: &Revolution, tol: f64) -> Result<f64, KernelError> {
    let (loops, _) = pieces(rev, tol)?;
    let sum: f64 = loops.iter().flatten().map(Piece::r2dz).sum();
    Ok((rev.angle.abs() / 2.0 * sum).abs())
}

/// Exact volume of an op 8 section (a full turn of straight segments).
pub fn volume(section: &Section) -> f64 {
    let pts = &section.points;
    let sum: f64 = (0..pts.len()).map(|i| {
        let ((r0, z0), (r1, z1)) = (pts[i], pts[(i + 1) % pts.len()]);
        (z1 - z0) * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0
    }).sum();
    (PI * sum).abs()
}

/// The loops resolved, validated and oriented (outer counter-clockwise first, holes clockwise), and each built loop's
/// request loop and whether it was reversed.
fn pieces(rev: &Revolution, tol: f64) -> Result<(Vec<Vec<Piece>>, Vec<(usize, bool)>), KernelError> {
    if !(tol.is_finite() && tol > 0.0) { return invalid("the tolerance must be positive and finite"); }
    if rev.loops.is_empty() { return invalid("a revolved section needs a loop"); }
    let mut loops = Vec::new();
    for (l, lp) in rev.loops.iter().enumerate() {
        if lp.len() < 2 { return invalid(format!("section loop {l} needs at least two segments")); }
        for (i, sg) in lp.iter().enumerate() {
            if !(sg.r.is_finite() && sg.z.is_finite()) { return invalid(format!("loop {l} point {i} is not finite")); }
            if sg.r < 0.0 || (sg.r > 0.0 && sg.r <= tol) {
                return invalid(format!("loop {l} point {i} must lie on the axis (radius 0) or clear of it (radius {})", sg.r));
            }
        }
        let mut ps = Vec::with_capacity(lp.len());
        for i in 0..lp.len() { ps.push(resolve(&lp[i], &lp[(i + 1) % lp.len()], tol, &format!("loop {l} segment {i}"))?); }
        let twice: f64 = ps.iter().map(Piece::area2).sum();
        let r2dz: f64 = ps.iter().map(Piece::r2dz).sum();
        if !twice.is_finite() || !r2dz.is_finite() || !(rev.angle * r2dz).is_finite() {
            return invalid("the revolved section is too large to measure (its area or volume overflows)");
        }
        if twice.abs() <= tol * tol { return invalid(format!("section loop {l} encloses no area")); }
        loops.push((twice.abs(), ps, twice > 0.0));
    }
    let outer = (0..loops.len()).max_by(|&i, &j| loops[i].0.total_cmp(&loops[j].0)).unwrap();
    loops.swap(0, outer);
    let mut origin: Vec<(usize, bool)> = (0..loops.len()).map(|i| (i, false)).collect();
    origin.swap(0, outer);
    let mut out: Vec<Vec<Piece>> = Vec::new();
    for (k, (_, ps, ccw)) in loops.into_iter().enumerate() {
        origin[k].1 = ccw != (k == 0);
        out.push(if ccw == (k == 0) { ps } else { ps.iter().rev().map(Piece::reversed).collect() });
    }
    // Every segment stays in r >= 0, touching the axis only at its ends and only on the outer loop; an arc centred off the
    // axis may not reach it (that would be a horn torus apex).
    for (l, lp) in out.iter().enumerate() {
        for (i, pc) in lp.iter().enumerate() {
            if l > 0 && (pc.s.0 == 0.0 || pc.e.0 == 0.0) { return invalid(format!("hole loop {l} touches the axis")); }
            if let Some(m) = pc.inner_min_radius() {
                if m <= tol { return invalid(format!("loop {l} segment {i} touches or crosses the axis")); }
            }
            if let Some((cr, _, _, _, _)) = pc.arc {
                if cr < -tol { return invalid(format!("loop {l} segment {i}: the arc's centre is on the far side of the axis")); }
            }
        }
        // A point on the axis must be where the boundary runs along or across it (a neighbour on the axis, or a flat disc):
        // two curved walls meeting only there would pinch the solid to a point (two shells at one vertex).
        let n = lp.len();
        for i in 0..n {
            if lp[i].s.0 != 0.0 { continue; }
            let along = |pc: &Piece| pc.arc.is_none() && pc.mer.is_none() && (pc.s.0 == 0.0 && pc.e.0 == 0.0 || (pc.e.1 - pc.s.1).abs() <= tol);
            if !along(&lp[(i + n - 1) % n]) && !along(&lp[i]) {
                return invalid(format!("loop {l} point {i} touches the axis only at a point (the solid would pinch there)"));
            }
        }
    }
    // Simple: no two segments meet except neighbours at their shared point; loops apart.
    let all: Vec<(usize, usize, Piece)> = out.iter().enumerate()
        .flat_map(|(l, lp)| lp.iter().enumerate().map(move |(i, pc)| (l, i, *pc))).collect();
    for u in 0..all.len() {
        for v in u + 1..all.len() {
            let ((l1, i1, p1), (l2, i2, p2)) = (all[u], all[v]);
            let n = out[l1].len();
            let after = l1 == l2 && (i1 + 1) % n == i2;
            let before = l1 == l2 && (i2 + 1) % n == i1;
            if p1.mer.is_some() || p2.mer.is_some() {
                let mut shared = Vec::new();
                if after { shared.push((1.0, 0.0)); }
                if before { shared.push((0.0, 1.0)); }
                if !apart(&p1, &p2, &shared, tol) {
                    return invalid(format!("section segments (loop {l1} segment {i1}) and (loop {l2} segment {i2}) cross, touch or could not be proved apart"));
                }
                continue;
            }
            if after && before {
                // A two-segment loop: the pieces share both ends; they must not overlap or meet elsewhere.
                if (p1.mid().0 - p2.mid().0).hypot(p1.mid().1 - p2.mid().1) <= tol
                    || crossings(&p1, &p2).iter().any(|&x| (x.0 - p1.s.0).hypot(x.1 - p1.s.1) > tol && (x.0 - p1.e.0).hypot(x.1 - p1.e.1) > tol) {
                    return invalid(format!("loop {l1} segments {i1} and {i2} overlap or cross"));
                }
                continue;
            }
            let shared = if after { Some(p1.e) } else if before { Some(p1.s) } else { None };
            if distance(&p1, &p2, shared, tol) <= tol {
                return invalid(format!("section segments (loop {l1} segment {i1}) and (loop {l2} segment {i2}) cross or touch"));
            }
            // Neighbours folding back onto each other share more than the end: their tangents there are opposite.
            if shared.is_some() && (p1.mid().0 - closest_on(&p2, p1.mid()).0).hypot(p1.mid().1 - closest_on(&p2, p1.mid()).1) <= tol {
                return invalid(format!("loop {l1} segments {i1} and {i2} overlap"));
            }
        }
    }
    for (l, lp) in out.iter().enumerate().skip(1) {
        let q = lp[0].mid();
        if out[0].iter().map(|pc| ray_crossings(pc, q)).sum::<usize>() % 2 == 0 {
            return invalid(format!("hole loop {l} is not inside the outer loop"));
        }
    }
    Ok((out, origin))
}

pub fn revolve(section: &Section, tol: f64) -> Result<Body, KernelError> {
    if section.points.len() < 3 { return invalid("a revolved section needs at least three points"); }
    revolve_general(&Revolution::from_section(section), tol)
}

/// The op 9 construction.
pub fn revolve_general(rev: &Revolution, tol: f64) -> Result<Body, KernelError> {
    revolve_general_named(rev, tol).map(|(body, _)| body)
}

/// The op 9 construction and the provenance of its every face, edge and vertex (provenance.rs roles 10-21, request
/// loop/segment indices; angle 0 is the start, the requested signed angle the end).
pub fn revolve_general_named(rev: &Revolution, tol: f64) -> Result<(Body, Provenance), KernelError> {
    let (p, mut a, x) = (rev.p, rev.a, rev.x);
    // Stated as what must hold, so a NaN fails every test (Codex acceptance 2026-10-29).
    if !(p.is_finite() && a.is_finite() && x.is_finite()
        && (a.length() - 1.0).abs() <= 1e-9 && (x.length() - 1.0).abs() <= 1e-9 && a.dot(x).abs() <= 1e-9) {
        return invalid("a revolve needs an axis point, a unit axis and a unit seam direction at right angles");
    }
    if !(rev.angle.is_finite() && rev.angle != 0.0 && rev.angle.abs() <= TAU + 1e-12) {
        return invalid("the revolve angle must satisfy 0 < |angle| <= 2 pi");
    }
    // A negative turn is a positive turn about -a with heights mirrored (reversing each loop's orientation, which the
    // pieces re-orient).
    let mut rev = rev.clone();
    if rev.angle < 0.0 {
        a = a.scale(-1.0);
        for lp in &mut rev.loops {
            for sg in lp.iter_mut() {
                sg.z = -sg.z;
                sg.kind = match sg.kind {
                    Kind::Arc { cr, cz, ccw } => Kind::Arc { cr, cz: -cz, ccw: !ccw },
                    Kind::Ellipse { cr, cz, ur, uz, vr, vz, forward } => Kind::Ellipse { cr, cz: -cz, ur, uz: -uz, vr, vz: -vz, forward },
                    Kind::Bezier { c1, c2 } => Kind::Bezier { c1: (c1.0, -c1.1), c2: (c2.0, -c2.1) },
                    k => k,
                };
            }
        }
        rev.angle = -rev.angle;
    }
    let (loops, origin) = pieces(&rev, tol)?;
    let full = (rev.angle - TAU).abs() <= 1e-12;
    let phi = if full { TAU } else { rev.angle };
    let y = a.cross(x);
    let rho_at = |t: f64| x.scale(t.cos()).add(y.scale(t.sin()));
    let at = |r: f64, z: f64, t: f64| p.add(a.scale(z)).add(rho_at(t).scale(r));
    let split_faces = full || phi > PI;
    let mid = if full { PI } else { phi / 2.0 };
    let curved = |pc: &Piece| pc.arc.is_some() || pc.mer.is_some() || (pc.e.1 - pc.s.1).abs() > tol;
    let touches = |pc: &Piece| pc.s.0 == 0.0 || pc.e.0 == 0.0;

    let mut body = Body::default();
    let mut prov = Provenance::default();
    // Built piece i of built loop l is request segment seg(l, i); its start point is request segment point(l, i)'s start.
    let seg = |l: usize, i: usize| { let (rl, rv) = origin[l]; let m = loops[l].len(); (rl, if rv { m - 1 - i } else { i }) };
    let point = |l: usize, i: usize| { let (rl, rv) = origin[l]; let m = loops[l].len(); (rl, if rv { (m - i) % m } else { i }) };
    let angle_role = |t: f64, at_start: u32, at_end: u32, at_mid: u32| if t == 0.0 { at_start } else if !full && t == phi { at_end } else { at_mid };
    // Each point clear of the axis: vertices at its boundary angles and the arcs between them (a full turn closes back to
    // the vertex at 0). Each axis point: one vertex.
    struct Rim { angles: Vec<f64>, verts: Vec<usize>, arcs: Vec<usize> }
    let mut rims: Vec<Vec<Option<Rim>>> = Vec::new();
    let mut axis_vertex: Vec<Vec<usize>> = Vec::new();
    for (l, lp) in loops.iter().enumerate() {
        let n = lp.len();
        let (mut lr, mut la) = (Vec::with_capacity(n), vec![usize::MAX; n]);
        for i in 0..n {
            let (r, z) = lp[i].s;
            let (rl, rs) = point(l, i);
            if r == 0.0 {
                la[i] = body.vertices.len();
                body.vertices.push(p.add(a.scale(z)));
                prov.vertices.push(Record::new(AXIS_POINT, &[rl, rs]));
                lr.push(None);
                continue;
            }
            let prev = &lp[(i + n - 1) % n];
            let split = split_faces && ((curved(prev) && touches(prev)) || (curved(&lp[i]) && touches(&lp[i])));
            let angles: Vec<f64> = match (full, split) {
                (true, false) => vec![0.0],
                (true, true) => vec![0.0, PI],
                (false, false) => vec![0.0, phi],
                (false, true) => vec![0.0, mid, phi],
            };
            let verts: Vec<usize> = angles.iter().map(|&t| {
                body.vertices.push(at(r, z, t));
                prov.vertices.push(Record::new(angle_role(t, START_POINT, END_POINT, MID_POINT), &[rl, rs]));
                body.vertices.len() - 1
            }).collect();
            let ends: Vec<f64> = if full { angles.iter().cloned().chain([TAU]).collect() } else { angles.clone() };
            let centre = p.add(a.scale(z));
            let mut arcs = Vec::new();
            for k in 0..ends.len() - 1 {
                let end = if full && k + 1 == angles.len() { verts[0] } else { verts[k + 1] };
                arcs.push(body.edges.len());
                body.edges.push(Edge { curve: Curve::Arc { c: centre, n: a, x, r, t0: ends[k], t1: ends[k + 1] }, start: verts[k], end });
                prov.edges.push(Record::new(PARALLEL, &[rl, rs, if ends.len() == 2 { 0 } else { k + 1 }]));
            }
            lr.push(Some(Rim { angles, verts, arcs }));
        }
        rims.push(lr);
        axis_vertex.push(la);
    }
    let vertex_at = |l: usize, i: usize, t: f64| -> usize {
        match &rims[l][i] {
            None => axis_vertex[l][i],
            Some(rim) => {
                let t = if full && t == TAU { 0.0 } else { t };
                rim.verts[rim.angles.iter().position(|&u| u == t).expect("a boundary angle of this point")]
            }
        }
    };
    // Arcs of a point between two boundary angles, in increasing angle.
    let arcs_between = |l: usize, i: usize, t0: f64, t1: f64| -> Vec<usize> {
        match &rims[l][i] {
            None => Vec::new(),
            Some(rim) => {
                let ends: Vec<f64> = if full { rim.angles.iter().cloned().chain([TAU]).collect() } else { rim.angles.clone() };
                (0..ends.len() - 1).filter(|&k| ends[k] >= t0 && ends[k + 1] <= t1).map(|k| rim.arcs[k]).collect()
            }
        }
    };
    // A piece's meridian curve at angle t, from its start to its end.
    let meridian = |pc: &Piece, t: f64| -> Curve {
        if let Some((m, t0, t1)) = pc.mer {
            // The map (r, z) -> 3D at angle t is affine, so control points and conjugate semi-diameters map exactly.
            let lift = |v: (f64, f64)| rho_at(t).scale(v.0).add(a.scale(v.1));
            return match m {
                Meridian::Bezier(_) => Curve::Bezier { p: m.sub_bezier(t0, t1).map(|q| at(q.0, q.1, t)) },
                Meridian::Ellipse { c, u, v } => Curve::Conic { c: at(c.0, c.1, t), u: lift(u), v: lift(v), t0, t1 },
            };
        }
        match pc.arc {
            None => Curve::Line { a: at(pc.s.0, pc.s.1, t), b: at(pc.e.0, pc.e.1, t) },
            Some((cr, cz, rho, _, sw)) => {
                let c = at(cr, cz, t);
                let rr = rho_at(t);
                // Counter-clockwise in (radius, height) turns rho towards a: the arc runs about rho x a.
                let n = if sw > 0.0 { rr.cross(a) } else { a.cross(rr) };
                Curve::Arc { c, n, x: at(pc.s.0, pc.s.1, t).sub(c).scale(1.0 / rho), r: rho, t0: 0.0, t1: sw.abs() }
            }
        }
    };
    let (mut cap0, mut cap1): (Vec<Vec<(usize, bool)>>, Vec<Vec<(usize, bool)>>) = (Vec::new(), Vec::new());
    for (l, lp) in loops.iter().enumerate() {
        let n = lp.len();
        let (mut c0, mut c1) = (Vec::new(), Vec::new());
        for i in 0..n {
            let j = (i + 1) % n;
            let pc = lp[i];
            let (rl, rs) = seg(l, i);
            if pc.arc.is_none() && pc.mer.is_none() && pc.s.0 == 0.0 && pc.e.0 == 0.0 {
                // On the axis: nothing in a full turn; the caps' shared edge in a part turn.
                if !full {
                    let e = body.edges.len();
                    body.edges.push(Edge { curve: meridian(&pc, 0.0), start: axis_vertex[l][i], end: axis_vertex[l][j] });
                    prov.edges.push(Record::new(AXIS_SEGMENT, &[rl, rs]));
                    c0.push((e, true));
                    c1.push((e, false));
                }
                continue;
            }
            let planar = !curved(&pc);
            let split = split_faces && !planar && touches(&pc);
            let bounds: Vec<f64> = match (full, split) {
                (true, false) => vec![0.0, TAU],
                (true, true) => vec![0.0, PI, TAU],
                (false, false) => vec![0.0, phi],
                (false, true) => vec![0.0, mid, phi],
            };
            let surface = match pc.arc {
                // Its natural normal points right of the meridian's travel, which is outward here.
                None if pc.mer.is_some() => Surface::Revolution { p, a, x, m: pc.mer.unwrap().0 },
                None if planar => Surface::Plane { p: p.add(a.scale(pc.s.1)), n: a, x },
                None => {
                    let ((r0, z0), (r1, z1)) = (pc.s, pc.e);
                    if r0 == r1 { Surface::Cylinder { p: p.add(a.scale(z0)), a, x, r: r0 } } else {
                        // A cone's reference point is at a positive-radius end (an apex is only a vertex).
                        let (rb, zb) = if r0 > 0.0 { (r0, z0) } else { (r1, z1) };
                        Surface::Cone { p: p.add(a.scale(zb)), a, x, r: rb, k: (r1 - r0) / (z1 - z0) }
                    }
                }
                Some((cr, cz, rho, _, _)) => {
                    if cr.abs() <= tol { Surface::Sphere { c: p.add(a.scale(cz)), a, x, r: rho } }
                    else { Surface::Torus { c: p.add(a.scale(cz)), a, x, big: cr, r: rho } }
                }
            };
            let sense = match pc.arc {
                None if pc.mer.is_some() => true,
                None if planar => pc.e.0 < pc.s.0,
                None => pc.e.1 > pc.s.1,
                Some((_, _, _, _, sw)) => sw > 0.0,
            };
            if full && planar {
                // An annulus or disc bounded by whole circles: the first point's forward, the second's reversed.
                let first: Vec<(usize, bool)> = arcs_between(l, i, 0.0, TAU).into_iter().map(|e| (e, true)).collect();
                let second: Vec<(usize, bool)> = arcs_between(l, j, 0.0, TAU).into_iter().rev().map(|e| (e, false)).collect();
                let ordered = if pc.s.0 >= pc.e.0 { [first, second] } else { [second, first] };
                body.faces.push(Face { surface, sense, loops: ordered.into_iter().filter(|lp| !lp.is_empty()).collect() });
                prov.faces.push(Record::new(SWEPT, &[rl, rs, 0]));
                continue;
            }
            // Meridian edges at the boundary angles (a full turn's last is its first).
            let mut me = Vec::new();
            for (k, &t) in bounds.iter().enumerate() {
                if full && k + 1 == bounds.len() { me.push(me[0]); continue; }
                me.push(body.edges.len());
                body.edges.push(Edge { curve: meridian(&pc, t), start: vertex_at(l, i, t), end: vertex_at(l, j, t) });
                prov.edges.push(Record::new(angle_role(t, START_MERIDIAN, END_MERIDIAN, MID_MERIDIAN), &[rl, rs]));
            }
            if !full {
                c0.push((me[0], true));
                c1.push((*me.last().unwrap(), false));
            }
            for k in 0..bounds.len() - 1 {
                let (t0, t1) = (bounds[k], bounds[k + 1]);
                let mut lpf: Vec<(usize, bool)> = arcs_between(l, i, t0, t1).into_iter().map(|e| (e, true)).collect();
                lpf.push((me[k + 1], true));
                lpf.extend(arcs_between(l, j, t0, t1).into_iter().rev().map(|e| (e, false)));
                lpf.push((me[k], false));
                body.faces.push(Face { surface: surface.clone(), sense, loops: vec![lpf] });
                prov.faces.push(Record::new(SWEPT, &[rl, rs, if bounds.len() == 2 { 0 } else { k + 1 }]));
            }
        }
        c1.reverse();
        cap0.push(c0);
        cap1.push(c1);
    }
    if !full {
        // Caps in the meridian planes: at 0 the outward normal is -y (the section loops run counter-clockwise about it);
        // at the angle it is the direction of turn there (loops reversed).
        body.faces.push(Face { surface: Surface::Plane { p, n: y, x: a }, sense: false, loops: cap0 });
        body.faces.push(Face { surface: Surface::Plane { p, n: a.cross(rho_at(phi)), x: a }, sense: true, loops: cap1 });
        prov.faces.push(Record::new(START_CAP, &[]));
        prov.faces.push(Record::new(END_CAP, &[]));
    }
    Ok((body, prov))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{mass, tessellate, topology};
    use std::f64::consts::PI;

    const TOL: f64 = 1e-6;

    fn axis(points: Vec<(f64, f64)>) -> Section {
        Section { p: Vec3::new(0.0, 0.0, 0.0), a: Vec3::new(0.0, 1.0, 0.0), x: Vec3::new(1.0, 0.0, 0.0), points }
    }

    /// GCmaker's default Pole (10 OD base 30 long, taper to 4 OD, top 30 long, 90 overall, wall .25), as
    /// SketchFactory_Pole lays out its outer then inner points.
    fn pole(wall: f64) -> Section {
        let (rb, rt) = (5.0, 2.0);
        let mut pts = vec![(rb, 0.0), (rb, 30.0), (rt, 60.0), (rt, 90.0)];
        if wall > 0.0 {
            pts.extend([(rt - wall, 90.0), (rt - wall, 60.0), (rb - wall, 30.0), (rb - wall, 0.0)]);
        } else {
            pts.extend([(0.0, 90.0), (0.0, 0.0)]);
        }
        axis(pts)
    }

    /// Volume of a hollow frustum stack by hand: pi h (R0^2 + R0 R1 + R1^2) / 3 per piece, outer less inner.
    fn frusta(pieces: &[(f64, f64, f64)]) -> f64 { pieces.iter().map(|&(r0, r1, h)| PI * h * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0).sum() }

    #[test]
    fn hollow_pole_is_a_closed_exact_solid_with_the_frustum_volume() {
        let s = pole(0.25);
        let body = revolve(&s, TOL).unwrap();
        let report = topology::check(&body, TOL).unwrap();
        assert_eq!((report.faces, report.edges, report.vertices), (8, 14, 8));
        let expected = frusta(&[(5.0, 5.0, 30.0), (5.0, 2.0, 30.0), (2.0, 2.0, 30.0)])
            - frusta(&[(4.75, 4.75, 30.0), (4.75, 1.75, 30.0), (1.75, 1.75, 30.0)]);
        assert!((volume(&s) - expected).abs() <= 1e-10 * expected);
        let m = mass::measure(&body, TOL).unwrap();
        assert!(m.exact);
        assert!((m.volume - expected).abs() <= 1e-9 * expected, "{} vs {expected}", m.volume);
        // Lateral area: cylinders 2 pi r h, cones pi (r0 + r1) slant; caps the annuli.
        let slant = (9.0f64 + 900.0).sqrt();
        let area = 2.0 * PI * 30.0 * (5.0 + 2.0 + 4.75 + 1.75) + PI * slant * (7.0 + 6.5)
            + PI * (25.0 - 4.75 * 4.75) + PI * (4.0 - 1.75 * 1.75);
        assert!((m.area - area).abs() <= 1e-9 * area, "{} vs {area}", m.area);
    }

    #[test]
    fn solid_pole_and_the_reversed_section_give_the_same_solid() {
        let mut s = pole(0.0);
        let body = revolve(&s, TOL).unwrap();
        topology::check(&body, TOL).unwrap();
        let expected = frusta(&[(5.0, 5.0, 30.0), (5.0, 2.0, 30.0), (2.0, 2.0, 30.0)]);
        assert!((mass::measure(&body, TOL).unwrap().volume - expected).abs() <= 1e-9 * expected);
        s.points.reverse();
        let reversed = revolve(&s, TOL).unwrap();
        assert!((mass::measure(&reversed, TOL).unwrap().volume - expected).abs() <= 1e-9 * expected);
    }

    #[test]
    fn a_tilted_axis_far_away_measures_the_same() {
        let mut s = pole(0.25);
        s.p = Vec3::new(1000.0, -250.0, 40.0);
        s.a = Vec3::new(1.0, 2.0, 2.0).unit().unwrap();
        s.x = Vec3::new(2.0, -1.0, 0.0).unit().unwrap();
        let body = revolve(&s, TOL).unwrap();
        topology::check(&body, TOL).unwrap();
        let v = volume(&s);
        assert!((mass::measure(&body, TOL).unwrap().volume - v).abs() <= 1e-8 * v);
    }

    #[test]
    fn display_mesh_lies_on_the_surfaces_within_the_chord() {
        let body = revolve(&pole(0.25), TOL).unwrap();
        let mesh = tessellate::tessellate(&body, 0.01).unwrap();
        assert!(!mesh.triangles.is_empty());
        // Every display point on the outer taper is on the cone (radius 5 - (h - 30) / 10 between heights 30 and 60).
        for q in &mesh.positions {
            let r = q.x.hypot(q.z);
            if q.y > 30.0 + 1e-9 && q.y < 60.0 - 1e-9 && r > 3.0 && r < 5.0 - 1e-9 {
                let cone = 5.0 - (q.y - 30.0) * 0.1;
                if (r - cone).abs() < 0.2 { assert!((r - cone).abs() <= 1e-9, "{r} vs {cone} at {q:?}"); }
            }
        }
        // Seams are one support each, so not display edges; the eight circles are.
        assert_eq!(mesh.display_edges.len(), 8);
    }

    #[test]
    fn cones_cross_the_wire_and_cuts_on_them_are_refused_by_name() {
        let body = revolve(&pole(0.25), TOL).unwrap();
        let wire = crate::neutral::write_body(&body);
        let mut reader = crate::neutral::Reader::with_tolerance(&wire, TOL);
        let back = reader.body().unwrap();
        reader.finish().unwrap();
        assert_eq!(crate::neutral::write_body(&back), wire);
        assert!(back.faces.iter().filter(|f| matches!(f.surface, Surface::Cone { .. })).count() == 2);
        match crate::boolean::cut(&body, &body, TOL) {
            Err(KernelError::Unsupported(m)) => assert!(m.contains("conical faces"), "{m}"),
            other => panic!("expected a refusal, got {other:?}"),
        }
    }

    /// Builds, checks, measures (exact) and meshes a general revolve against a volume and area found by hand.
    fn solid(rev: &Revolution, volume: f64, area: Option<f64>) -> Body {
        let body = revolve_general(rev, TOL).unwrap_or_else(|e| panic!("{e:?}"));
        topology::check(&body, TOL).unwrap_or_else(|e| panic!("{e:?}"));
        let exact = volume_of(rev, TOL).unwrap();
        assert!((exact - volume).abs() <= 1e-9 * volume, "construction {exact} vs {volume}");
        let m = mass::measure(&body, TOL).unwrap_or_else(|e| panic!("{e:?}"));
        assert!(m.exact);
        assert!((m.volume - volume).abs() <= 1e-9 * volume, "measured {} vs {volume}", m.volume);
        if let Some(area) = area { assert!((m.area - area).abs() <= 1e-9 * area, "area {} vs {area}", m.area); }
        let mesh = tessellate::tessellate(&body, 0.005).unwrap_or_else(|e| panic!("{e:?}"));
        assert!(!mesh.triangles.is_empty());
        body
    }

    fn general(angle: f64, loops: Vec<Vec<Segment>>) -> Revolution {
        Revolution { p: Vec3::new(0.0, 0.0, 0.0), a: Vec3::new(0.0, 1.0, 0.0), x: Vec3::new(1.0, 0.0, 0.0), angle, loops }
    }
    fn line(r: f64, z: f64) -> Segment { Segment { r, z, kind: Kind::Line } }
    fn arc(r: f64, z: f64, cr: f64, cz: f64, ccw: bool) -> Segment { Segment { r, z, kind: Kind::Arc { cr, cz, ccw } } }

    #[test]
    fn cone_to_its_apex_full_and_part_turns() {
        // A solid cone, radius 2, height 5: pi r^2 h / 3; lateral pi r slant, base pi r^2.
        let (r, h) = (2.0, 5.0);
        let section = vec![line(0.0, 0.0), line(r, 0.0), line(0.0, h)];
        let full = PI * r * r * h / 3.0;
        let area = PI * r * (r * r + h * h).sqrt() + PI * r * r;
        solid(&general(TAU, vec![section.clone()]), full, Some(area));
        // Part turns: the volume scales; the caps (two triangles) add 2 (r h / 2).
        for angle in [PI / 3.0, 1.5 * PI, -0.75 * PI] {
            let k = angle.abs() / TAU;
            solid(&general(angle, vec![section.clone()]), full * k, Some(area * k + r * h));
        }
        // Op 8 admits the apex too.
        let body = revolve(&axis(vec![(0.0, 0.0), (2.0, 0.0), (0.0, 5.0)]), TOL).unwrap();
        assert!((mass::measure(&body, TOL).unwrap().volume - full).abs() <= 1e-9 * full);
    }

    #[test]
    fn sphere_ball_and_hemispherical_part_turns() {
        // A ball of radius 3 from a half disc: the semicircle centred on the axis, closed along the axis.
        let rr = 3.0;
        let section = vec![arc(0.0, -rr, 0.0, 0.0, true), line(0.0, rr)];
        solid(&general(TAU, vec![section.clone()]), 4.0 / 3.0 * PI * rr.powi(3), Some(4.0 * PI * rr * rr));
        solid(&general(PI / 2.0, vec![section.clone()]), PI / 3.0 * rr.powi(3), Some(PI * rr * rr + PI * rr * rr));
        // A hollow sphere (shell 3 to 2): outer arc up, line in to the inner pole, inner arc down, line out.
        let shell = vec![arc(0.0, -3.0, 0.0, 0.0, true), line(0.0, 3.0), line(0.0, 2.0), arc(0.0, 2.0, 0.0, 0.0, false), line(0.0, -2.0)];
        let _ = shell;
        let shell = vec![arc(0.0, -3.0, 0.0, 0.0, true), line(0.0, 3.0), arc(0.0, 2.0, 0.0, 0.0, false), line(0.0, -2.0)];
        solid(&general(TAU, vec![shell]), 4.0 / 3.0 * PI * (27.0 - 8.0), Some(4.0 * PI * (9.0 + 4.0)));
    }

    #[test]
    fn ring_torus_from_a_circular_section_and_a_hole() {
        // A full circle of radius 1 at radius 4 (two arcs): volume 2 pi^2 R r^2, area 4 pi^2 R r.
        let (big, r) = (4.0, 1.0);
        let ring = vec![arc(big + r, 0.0, big, 0.0, true), arc(big - r, 0.0, big, 0.0, true)];
        solid(&general(TAU, vec![ring.clone()]), 2.0 * PI * PI * big * r * r, Some(4.0 * PI * PI * big * r));
        solid(&general(PI / 2.0, vec![ring.clone()]), PI * PI * big * r * r / 2.0, Some(PI * PI * big * r + 2.0 * PI * r * r));
        // A square tube 3..6 x -2..2 with that circle as a hole: annular prism less the torus.
        let outer = vec![line(3.0 - 0.5, -2.0), line(6.0, -2.0), line(6.0, 2.0), line(2.5, 2.0)];
        let volume = PI * (36.0 - 6.25) * 4.0 - 2.0 * PI * PI * big * r * r;
        solid(&general(TAU, vec![outer, ring]), volume, None);
    }

    #[test]
    fn a_trimmed_spindle_arc_and_a_negative_turn_keep_their_world_placement() {
        // An outer bulge: arc centred at radius 1 with radius 2 (r >= big, a spindle torus), only its far side used.
        let bulge = vec![line(2.0, -1.0), arc(2.0, -1.0, 1.0, 0.0, true)];
        let s = (2.0f64).hypot(0.0);
        let _ = s;
        let pts = vec![line(1.0, -(3.0f64).sqrt()), arc(1.0 + 1.0, -(3.0f64).sqrt(), 1.0, 0.0, true), line(1.0, (3.0f64).sqrt())];
        let _ = (bulge, pts);
        // Lens: the chord r = 2 between heights -sqrt(3) and sqrt(3) on the circle (r - 1)^2 + z^2 = 4, bulging outwards.
        let h = (3.0f64).sqrt();
        let lens = vec![line(2.0, h), arc(2.0, -h, 1.0, 0.0, true)];
        let body = revolve_general(&general(TAU, vec![lens.clone()]), TOL).unwrap();
        topology::check(&body, TOL).unwrap();
        let v = volume_of(&general(TAU, vec![lens.clone()]), TOL).unwrap();
        assert!((mass::measure(&body, TOL).unwrap().volume - v).abs() <= 1e-9 * v);
        tessellate::tessellate(&body, 0.005).unwrap();
        // A quarter turn backwards lies at negative angles: every vertex has z <= 0 (y = a x x is +z here... x along +X,
        // a along +Y, so y = a x x = -Z; a positive turn goes towards -Z and a negative one towards +Z).
        let back = revolve_general(&general(-PI / 2.0, vec![vec![line(1.0, 0.0), line(2.0, 0.0), line(2.0, 1.0), line(1.0, 1.0)]]), TOL).unwrap();
        assert!(back.vertices.iter().all(|q| q.z >= -1e-12 && q.y >= -1e-12 && q.y <= 1.0 + 1e-12));
        let ahead = revolve_general(&general(PI / 2.0, vec![vec![line(1.0, 0.0), line(2.0, 0.0), line(2.0, 1.0), line(1.0, 1.0)]]), TOL).unwrap();
        assert!(ahead.vertices.iter().all(|q| q.z <= 1e-12));
    }

    #[test]
    fn refuses_ellipse_ends_off_the_curve_and_zero_angles() {
        // The ellipse c (1, 1), u (0, -1), v (1, 0) passes (1, 0) and (2, 1); (2, 1.5) is off it.
        let ellipse = vec![Segment { r: 1.0, z: 0.0, kind: Kind::Ellipse { cr: 1.0, cz: 1.0, ur: 0.0, uz: -1.0, vr: 1.0, vz: 0.0, forward: true } },
            line(2.0, 1.5), line(1.0, 1.5)];
        assert!(matches!(revolve_general(&general(TAU, vec![ellipse]), TOL), Err(KernelError::InvalidInput(_))));
        let horn = vec![line(0.0, 0.0), arc(2.0, 0.0, 1.0, 0.0, true), line(0.0, 0.0)];
        assert!(revolve_general(&general(TAU, vec![horn]), TOL).is_err());
        assert!(matches!(revolve_general(&general(0.0, vec![vec![line(1.0, 0.0), line(2.0, 0.0), line(2.0, 1.0)]]), TOL), Err(KernelError::InvalidInput(_))));
    }

    #[test]
    fn refuses_apexes_crossings_and_points_just_off_the_axis() {
        let bow = axis(vec![(1.0, 0.0), (3.0, 2.0), (3.0, 0.0), (1.0, 2.0)]);
        assert!(matches!(revolve(&bow, TOL), Err(KernelError::InvalidInput(_))));
        let near = axis(vec![(1e-9, 0.0), (2.0, 0.0), (2.0, 1.0), (1e-9, 1.0)]);
        assert!(matches!(revolve(&near, TOL), Err(KernelError::InvalidInput(_))));
        let negative = axis(vec![(-1.0, 0.0), (2.0, 0.0), (2.0, 1.0)]);
        assert!(matches!(revolve(&negative, TOL), Err(KernelError::InvalidInput(_))));
        let huge = axis(vec![(1.0, 0.0), (1e200, 0.0), (1e200, 1e200), (1.0, 1e200)]);
        assert!(matches!(revolve(&huge, TOL), Err(KernelError::InvalidInput(_))));
        for bad in [Vec3::new(f64::NAN, 0.0, 0.0), Vec3::new(0.0, f64::NAN, 0.0)] {
            let mut s = pole(0.25);
            s.a = Vec3::new(0.0, 1.0, 0.0).add(bad);
            assert!(matches!(revolve(&s, TOL), Err(KernelError::InvalidInput(_))));
            let mut s = pole(0.25);
            s.x = Vec3::new(1.0, 0.0, 0.0).add(bad);
            assert!(matches!(revolve(&s, TOL), Err(KernelError::InvalidInput(_))));
        }
    }
}



#[cfg(test)]
mod diag {
    use super::*;
    #[test]
    fn diag_spindle() {
        let h = (3.0f64).sqrt();
        let rev = Revolution { p: Vec3::new(0.0, 0.0, 0.0), a: Vec3::new(0.0, 1.0, 0.0), x: Vec3::new(1.0, 0.0, 0.0), angle: TAU,
            loops: vec![vec![Segment { r: 0.0, z: -h, kind: Kind::Arc { cr: 1.0, cz: 0.0, ccw: true } }, Segment { r: 0.0, z: h, kind: Kind::Line }]] };
        let body = revolve_general(&rev, 1e-7).unwrap();
        println!("result {:?}", crate::tessellate::tessellate(&body, 0.04).err());
    }
}

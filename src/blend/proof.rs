//! Proofs that a blend leaves a valid boundary (Claude 2026-11-01): clearance between curves, winding numbers and
//! clearance of a curve from a region. Each works on pieces that are proved to lie inside a ball:
//! - lines, circles and ellipses c + u cos t + v sin t by their speed bound (|C'| <= |b - a|, r, |u| + |v|): the piece
//!   over [t0, t1] lies within speed (t1 - t0) / 2 of its midpoint;
//! - Bezier and positive-weight rational Bezier pieces by their control hull (they lie inside it), split by de Casteljau.
//! A proof that runs out of work is reported as not proved, never assumed.

use crate::geom::{Curve, Vec3};

const BUDGET: usize = 400_000;

#[derive(Clone, Copy)]
enum Piece {
    Span { t0: f64, t1: f64 },
    Hull { p: [Vec3; 4], w: [f64; 4] },
}

struct Sampler<'a> {
    curve: &'a Curve,
    speed: f64,
}

impl<'a> Sampler<'a> {
    fn of(curve: &'a Curve) -> Result<(Self, Piece), String> {
        let (speed, piece) = match *curve {
            Curve::Line { a, b } => (b.sub(a).length(), Piece::Span { t0: 0.0, t1: 1.0 }),
            Curve::Arc { r, t0, t1, .. } => (r, Piece::Span { t0, t1 }),
            Curve::Conic { u, v, t0, t1, .. } => (u.length() + v.length(), Piece::Span { t0, t1 }),
            Curve::Bezier { p } => (0.0, Piece::Hull { p, w: [1.0; 4] }),
            Curve::RationalBezier { p, w } => (0.0, Piece::Hull { p, w }),
            Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } => return Err("an intersection-curve edge has no clearance bound here".into()),
        };
        Ok((Self { curve, speed }, piece))
    }
    fn ball(&self, piece: &Piece) -> (Vec3, f64) {
        match *piece {
            Piece::Span { t0, t1 } => (self.curve.eval((t0 + t1) / 2.0).0, self.speed * (t1 - t0) / 2.0),
            Piece::Hull { p, .. } => {
                let c = p.iter().fold(Vec3::new(0.0, 0.0, 0.0), |s, q| s.add(*q)).scale(0.25);
                (c, p.iter().map(|q| q.distance(c)).fold(0.0, f64::max))
            }
        }
    }
    fn ends(&self, piece: &Piece) -> (Vec3, Vec3) {
        match *piece {
            Piece::Span { t0, t1 } => (self.curve.eval(t0).0, self.curve.eval(t1).0),
            Piece::Hull { p, .. } => (p[0], p[3]),
        }
    }
}

fn halves(piece: &Piece) -> (Piece, Piece) {
    match *piece {
        Piece::Span { t0, t1 } => {
            let m = (t0 + t1) / 2.0;
            (Piece::Span { t0, t1: m }, Piece::Span { t0: m, t1 })
        }
        Piece::Hull { p, w } => {
            let [l, r] = crate::rational::split(p, w, 0.5);
            (Piece::Hull { p: l.0, w: l.1 }, Piece::Hull { p: r.0, w: r.1 })
        }
    }
}

/// Whether two curves stay more than `gap` apart (Ok(true)), come within a few `gap` of each other (Ok(false)), or
/// could not be decided (Err).
pub(crate) fn apart(a: &Curve, b: &Curve, gap: f64) -> Result<bool, String> {
    let (sa, pa) = Sampler::of(a)?;
    let (sb, pb) = Sampler::of(b)?;
    let mut pending = vec![(pa, pb)];
    let mut work = 0;
    while let Some((x, y)) = pending.pop() {
        work += 1;
        if work > BUDGET { return Err("the clearance proof ran out of work".into()); }
        let ((cx, rx), (cy, ry)) = (sa.ball(&x), sb.ball(&y));
        if cx.distance(cy) > rx + ry + gap { continue; }
        if rx <= gap && ry <= gap { return Ok(false); }
        if rx >= ry {
            let (l, r) = halves(&x);
            pending.push((l, y));
            pending.push((r, y));
        } else {
            let (l, r) = halves(&y);
            pending.push((x, l));
            pending.push((x, r));
        }
    }
    Ok(true)
}

/// The winding number of point q (in the loop's plane, normal n) about a closed loop of curves. Every piece is split
/// until its ball excludes q, so the angle it sweeps about q is less than a half turn and equals the principal angle
/// between its ends.
pub(crate) fn winding(q: Vec3, n: Vec3, lp: &[Curve]) -> Result<f64, String> {
    let mut total = 0.0;
    let mut work = 0;
    for c in lp {
        let (s, p) = Sampler::of(c)?;
        let mut pending = vec![p];
        while let Some(piece) = pending.pop() {
            work += 1;
            if work > BUDGET { return Err("the winding proof ran out of work".into()); }
            let (centre, radius) = s.ball(&piece);
            if centre.distance(q) <= radius {
                if radius < 1e-12 * centre.length().max(1.0) { return Err("a point lies on the loop".into()); }
                let (l, r) = halves(&piece);
                pending.push(r);
                pending.push(l);
                continue;
            }
            let (from, to) = s.ends(&piece);
            let (f, t) = (from.sub(q), to.sub(q));
            total += f.cross(t).dot(n).atan2(f.dot(t));
        }
    }
    Ok(total / std::f64::consts::TAU)
}

/// Whether the curve keeps clear of a region: `test(centre, radius)` answers Some(true) when the ball is proved
/// outside the region, Some(false) when a small ball is inside it, None to split further.
pub(crate) fn clear_of(c: &Curve, test: &dyn Fn(Vec3, f64) -> Option<bool>) -> Result<bool, String> {
    let (s, p) = Sampler::of(c)?;
    let mut pending = vec![p];
    let mut work = 0;
    while let Some(piece) = pending.pop() {
        work += 1;
        if work > BUDGET { return Err("the clearance proof ran out of work".into()); }
        let (centre, radius) = s.ball(&piece);
        match test(centre, radius) {
            Some(true) => continue,
            Some(false) => return Ok(false),
            None => {
                let (l, r) = halves(&piece);
                pending.push(l);
                pending.push(r);
            }
        }
    }
    Ok(true)
}

/// The parameters in [t0, t1] where c + u cos t + v sin t crosses the line through l0 along ld (in the plane of normal
/// n): f(t) = (c - l0).m + (u.m) cos t + (v.m) sin t = 0, m = n x ld, solved in closed form.
pub(crate) fn trig_line_roots(c: Vec3, u: Vec3, v: Vec3, t0: f64, t1: f64, l0: Vec3, ld: Vec3, n: Vec3) -> Vec<f64> {
    let m = n.cross(ld);
    let (k, a, b) = (c.sub(l0).dot(m), u.dot(m), v.dot(m));
    let r = a.hypot(b);
    if r <= 1e-300 || k.abs() > r { return Vec::new(); }
    let phi = b.atan2(a);
    let delta = (-k / r).clamp(-1.0, 1.0).acos();
    let mut out = Vec::new();
    for base in [phi + delta, phi - delta] {
        let mut t = base + ((t0 - base) / std::f64::consts::TAU).ceil() * std::f64::consts::TAU;
        while t <= t1 {
            out.push(t);
            t += std::f64::consts::TAU;
        }
    }
    out
}

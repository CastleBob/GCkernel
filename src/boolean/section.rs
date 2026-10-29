//! Through cuts on right-prism stock as booleans of their sections (Claude, 2026-10-29).
//!
//! When the body is a right prism along a cardinal axis (two caps normal to it, every other face parallel to it) and
//! the tool is a right prism along the same axis spanning the body's whole extent, the result is the body's section
//! less the tool's section, swept over the body's extent. The 2D difference follows the public subdivide-and-
//! classify scheme (Martinez, Rueda and Feito, Computers & Geosciences 2009), written here for line and circular-arc
//! loops: every edge is split at every meeting with the other section (closed-form line/line, line/circle and
//! circle/circle points, plus endpoints lying on the other boundary), meeting points within the tolerance are one
//! vertex (Bob 2026-10-28: what agrees within the tolerance is the same), coincident pieces are decided once per pair,
//! and every other piece is kept or dropped by the winding of its midpoint. The prism property of each operand is
//! proved by its section area times its extent equalling its measured volume. Disconnected material is refused, as
//! in the other domains. Independent implementation from published mathematics; no third-party source.
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Profile};
use crate::KernelError;
use std::f64::consts::PI;

type Result<T> = std::result::Result<T, KernelError>;

fn unsupported(message: &str) -> KernelError {
    KernelError::Unsupported(message.into())
}

const LIMIT: usize = 2000;

#[derive(Clone, Copy, Debug, PartialEq)]
struct P {
    x: f64,
    y: f64,
}

impl P {
    fn new(x: f64, y: f64) -> P { P { x, y } }
    fn add(self, o: P) -> P { P::new(self.x + o.x, self.y + o.y) }
    fn sub(self, o: P) -> P { P::new(self.x - o.x, self.y - o.y) }
    fn scale(self, k: f64) -> P { P::new(self.x * k, self.y * k) }
    fn dot(self, o: P) -> f64 { self.x * o.x + self.y * o.y }
    fn cross(self, o: P) -> f64 { self.x * o.y - self.y * o.x }
    fn len(self) -> f64 { self.x.hypot(self.y) }
    fn dist(self, o: P) -> f64 { self.sub(o).len() }
    fn polar(c: P, r: f64, t: f64) -> P { P::new(c.x + r * t.cos(), c.y + r * t.sin()) }
}

/// A section edge: a line a -> b, a circular arc about c from angle a0 through the signed sweep (positive is
/// counter-clockwise), or an elliptic arc c + u cos t + v sin t (u, v independent semi-diameters) from t = a0 through
/// the signed sweep (Claude 2026-10-29, exact ellipses).
#[derive(Clone, Copy, Debug)]
enum Seg {
    Line { a: P, b: P },
    Arc { c: P, r: f64, a0: f64, sweep: f64 },
    Conic { c: P, u: P, v: P, a0: f64, sweep: f64 },
}

/// A point's coordinates (c, n) with w = c u + n v (Cramer in the plane): the ellipse's unit-circle coordinates.
fn unit_coords(u: P, v: P, w: P) -> (f64, f64) {
    let det = u.cross(v);
    (w.cross(v) / det, u.cross(w) / det)
}

/// A curve c + u cos t + v sin t over t in [a0, a0 + sweep] (either sign): its box, from the ends and every interior
/// extreme (x' = 0 where tan t = vx / ux, likewise y).
fn harmonic_box(c: P, u: P, v: P, a0: f64, sweep: f64) -> (P, P) {
    let at = |t: f64| P::new(c.x + u.x * t.cos() + v.x * t.sin(), c.y + u.y * t.cos() + v.y * t.sin());
    let (lo_t, hi_t) = if sweep >= 0. { (a0, a0 + sweep) } else { (a0 + sweep, a0) };
    let mut points = vec![at(lo_t), at(hi_t)];
    for base in [v.x.atan2(u.x), v.y.atan2(u.y)] {
        for k in -4..=4 {
            let t = base + k as f64 * PI;
            if t > lo_t && t < hi_t {
                points.push(at(t));
            }
        }
    }
    let lo = P::new(points.iter().map(|q| q.x).fold(f64::INFINITY, f64::min), points.iter().map(|q| q.y).fold(f64::INFINITY, f64::min));
    let hi = P::new(points.iter().map(|q| q.x).fold(f64::NEG_INFINITY, f64::max), points.iter().map(|q| q.y).fold(f64::NEG_INFINITY, f64::max));
    (lo, hi)
}

impl Seg {
    fn at(&self, s: f64) -> P {
        match *self {
            Seg::Line { a, b } => a.add(b.sub(a).scale(s)),
            Seg::Arc { c, r, a0, sweep } => P::polar(c, r, a0 + s * sweep),
            Seg::Conic { c, u, v, a0, sweep } => {
                let t = a0 + s * sweep;
                c.add(u.scale(t.cos())).add(v.scale(t.sin()))
            }
        }
    }
    fn start(&self) -> P { self.at(0.) }
    fn end(&self) -> P { self.at(1.) }
    /// Unit direction of travel at s.
    fn tangent(&self, s: f64) -> P {
        match *self {
            Seg::Line { a, b } => b.sub(a).scale(1. / b.dist(a)),
            Seg::Arc { a0, sweep, .. } => {
                let t = a0 + s * sweep;
                P::new(-t.sin(), t.cos()).scale(sweep.signum())
            }
            Seg::Conic { u, v, a0, sweep, .. } => {
                let t = a0 + s * sweep;
                let d = u.scale(-t.sin()).add(v.scale(t.cos())).scale(sweep.signum());
                d.scale(1. / d.len())
            }
        }
    }
    fn length(&self) -> f64 {
        match *self {
            Seg::Line { a, b } => a.dist(b),
            Seg::Arc { r, sweep, .. } => r * sweep.abs(),
            // Used only for closeness of parameters along the edge: a fine chord sum.
            Seg::Conic { .. } => (0..64).map(|i| self.at(i as f64 / 64.).dist(self.at((i + 1) as f64 / 64.))).sum(),
        }
    }
    fn reversed(&self) -> Seg {
        match *self {
            Seg::Line { a, b } => Seg::Line { a: b, b: a },
            Seg::Arc { c, r, a0, sweep } => Seg::Arc { c, r, a0: a0 + sweep, sweep: -sweep },
            Seg::Conic { c, u, v, a0, sweep } => Seg::Conic { c, u, v, a0: a0 + sweep, sweep: -sweep },
        }
    }
    /// The parameter of the point of the supporting curve nearest p (outside [0, 1] beyond the ends).
    fn param(&self, p: P) -> f64 {
        match *self {
            Seg::Line { a, b } => p.sub(a).dot(b.sub(a)) / b.sub(a).dot(b.sub(a)),
            Seg::Arc { c, a0, sweep, .. } => {
                let phi = (p.y - c.y).atan2(p.x - c.x);
                let d = ((phi - a0) * sweep.signum()).rem_euclid(2. * PI);
                // Past the end: nearer the end, or (beyond the gap's middle) just before the start.
                if d <= sweep.abs() + (2. * PI - sweep.abs()) * 0.5 { d / sweep.abs() } else { (d - 2. * PI) / sweep.abs() }
            }
            // The ellipse's own parameter of p (exact on the ellipse), wrapped as for the arc.
            Seg::Conic { c, u, v, a0, sweep } => {
                let (cc, nn) = unit_coords(u, v, p.sub(c));
                let phi = nn.atan2(cc);
                let d = ((phi - a0) * sweep.signum()).rem_euclid(2. * PI);
                if d <= sweep.abs() + (2. * PI - sweep.abs()) * 0.5 { d / sweep.abs() } else { (d - 2. * PI) / sweep.abs() }
            }
        }
    }
    /// Distance from p to this edge (its nearest point within the edge).
    fn distance(&self, p: P) -> f64 {
        let s = self.param(p);
        if (0. ..=1.).contains(&s) {
            match *self {
                Seg::Line { a, b } => (b.sub(a).cross(p.sub(a)) / b.dist(a)).abs(),
                Seg::Arc { c, r, .. } => (p.dist(c) - r).abs(),
                // The point of the ellipse at p's own parameter: exact on the ellipse, a close upper bound near it.
                Seg::Conic { .. } => p.dist(self.at(s)),
            }
        } else {
            p.dist(self.start()).min(p.dist(self.end()))
        }
    }
    /// Its part of a closed loop's signed area: the integral of (x dy - y dx) / 2.
    fn area(&self) -> f64 {
        match *self {
            Seg::Line { a, b } => a.cross(b) * 0.5,
            Seg::Arc { c, r, a0, sweep } => {
                let a1 = a0 + sweep;
                0.5 * (r * c.x * (a1.sin() - a0.sin()) - r * c.y * (a1.cos() - a0.cos()) + r * r * sweep)
            }
            Seg::Conic { c, u, v, a0, sweep } => {
                let a1 = a0 + sweep;
                let (dc, ds) = (a1.cos() - a0.cos(), a1.sin() - a0.sin());
                0.5 * (c.x * (u.y * dc + v.y * ds) - c.y * (u.x * dc + v.x * ds) + u.cross(v) * sweep)
            }
        }
    }
    /// The edge's box (rigorous: ends and interior extremes).
    fn bounds(&self) -> (P, P) {
        match *self {
            Seg::Line { a, b } => (P::new(a.x.min(b.x), a.y.min(b.y)), P::new(a.x.max(b.x), a.y.max(b.y))),
            Seg::Arc { c, r, a0, sweep } => harmonic_box(c, P::new(r, 0.), P::new(0., r), a0, sweep),
            Seg::Conic { c, u, v, a0, sweep } => harmonic_box(c, u, v, a0, sweep),
        }
    }
}

fn loop_area(lp: &[Seg]) -> f64 { lp.iter().map(Seg::area).sum() }

/// Points where two edges meet (candidates; each is checked to lie on both within the tolerance by the caller).
/// An ellipse meeting an arc or another ellipse (a quartic) is refused unless their boxes keep them apart, or it is
/// the same ellipse (coincident: its ends are the meetings).
fn meetings(s: &Seg, t: &Seg, tol: f64) -> Result<Vec<P>> {
    let mut points = vec![];
    for (e, other) in [(s, t), (t, s)] {
        for q in [e.start(), e.end()] {
            if other.distance(q) <= tol {
                points.push(q);
            }
        }
    }
    match (*s, *t) {
        (Seg::Line { a, b }, Seg::Line { a: c, b: d }) => {
            let (u, v) = (b.sub(a), d.sub(c));
            let den = u.cross(v);
            if den.abs() > 1e-12 * u.len() * v.len() {
                points.push(a.add(u.scale(c.sub(a).cross(v) / den)));
            }
        }
        (Seg::Line { a, b }, Seg::Arc { c, r, .. }) | (Seg::Arc { c, r, .. }, Seg::Line { a, b }) => {
            let (d, f) = (b.sub(a), a.sub(c));
            let (qa, qb, qc) = (d.dot(d), 2. * d.dot(f), f.dot(f) - r * r);
            let disc = qb * qb - 4. * qa * qc;
            if disc >= 0. {
                let root = disc.sqrt();
                for k in [(-qb - root) / (2. * qa), (-qb + root) / (2. * qa)] {
                    points.push(a.add(d.scale(k)));
                }
            } else {
                // Within the tolerance of tangency: the nearest point of the line to the centre.
                let k = -qb / (2. * qa);
                let q = a.add(d.scale(k));
                if (q.dist(c) - r).abs() <= tol {
                    points.push(q);
                }
            }
        }
        (Seg::Arc { c: c1, r: r1, .. }, Seg::Arc { c: c2, r: r2, .. }) => {
            let d = c1.dist(c2);
            if d > tol && d <= r1 + r2 + tol && d >= (r1 - r2).abs() - tol {
                let along = (r1 * r1 - r2 * r2 + d * d) / (2. * d);
                let h = (r1 * r1 - along * along).max(0.).sqrt();
                let u = c2.sub(c1).scale(1. / d);
                let base = c1.add(u.scale(along));
                points.push(base.add(P::new(-u.y, u.x).scale(h)));
                points.push(base.add(P::new(u.y, -u.x).scale(h)));
            }
        }
        // In the ellipse's unit-circle coordinates the line stays a line with the same parameter: a circle-line
        // quadratic, exact.
        (Seg::Line { a, b }, Seg::Conic { c, u, v, .. }) | (Seg::Conic { c, u, v, .. }, Seg::Line { a, b }) => {
            let (ax, ay) = unit_coords(u, v, a.sub(c));
            let (bx, by) = unit_coords(u, v, b.sub(c));
            let (d, f) = (P::new(bx - ax, by - ay), P::new(ax, ay));
            let (qa, qb, qc) = (d.dot(d), 2. * d.dot(f), f.dot(f) - 1.);
            let disc = qb * qb - 4. * qa * qc;
            if disc >= 0. {
                let root = disc.sqrt();
                for k in [(-qb - root) / (2. * qa), (-qb + root) / (2. * qa)] {
                    points.push(a.add(b.sub(a).scale(k)));
                }
            } else {
                points.push(a.add(b.sub(a).scale(-qb / (2. * qa))));
            }
        }
        (Seg::Conic { c, u, v, .. }, other @ (Seg::Arc { .. } | Seg::Conic { .. }))
        | (other @ Seg::Arc { .. }, Seg::Conic { c, u, v, .. }) => {
            // Same support, proven analytically (Codex review: sample points can lie on two different ellipses): the
            // other curve is A + B cos t + C sin t in this ellipse's unit-circle coordinates, so c^2 + n^2 - 1 =
            // k0 + k1 cos t + k2 sin t + k3 cos 2t + k4 sin 2t vanishes identically only if every k does.
            let (oc, ou, ov) = match other {
                Seg::Arc { c: oc, r, .. } => (oc, P::new(r, 0.), P::new(0., r)),
                Seg::Conic { c: oc, u: ou, v: ov, .. } => (oc, ou, ov),
                Seg::Line { .. } => unreachable!("lines meet ellipses above"),
            };
            let (a2, b2, c2) = (unit_coords(u, v, oc.sub(c)), unit_coords(u, v, ou), unit_coords(u, v, ov));
            let dot = |x: (f64, f64), y: (f64, f64)| x.0 * y.0 + x.1 * y.1;
            let k = [dot(a2, a2) + (dot(b2, b2) + dot(c2, c2)) / 2. - 1., 2. * dot(a2, b2), 2. * dot(a2, c2),
                (dot(b2, b2) - dot(c2, c2)) / 2., dot(b2, c2)];
            let same = k.iter().map(|x| x.abs()).sum::<f64>() * u.len().max(v.len()) / 2. <= tol;
            if !same {
                let (l1, h1) = s.bounds();
                let (l2, h2) = t.bounds();
                let apart = h1.x + tol < l2.x || h2.x + tol < l1.x || h1.y + tol < l2.y || h2.y + tol < l1.y;
                if !apart {
                    return Err(unsupported("an ellipse meeting an arc or another ellipse is not supported yet"));
                }
            }
        }
    }
    points.retain(|q| q.x.is_finite() && q.y.is_finite() && s.distance(*q) <= tol && t.distance(*q) <= tol);
    Ok(points)
}

/// Winding number of the loops about p, by a ray from p in a fixed direction off the axes. `None` when the ray
/// passes within roundoff of an edge end (the caller turns the ray and asks again).
fn winding_along(loops: &[Vec<Seg>], p: P, angle: f64) -> Option<i32> {
    let (cs, sn) = (angle.cos(), angle.sin());
    let local = |q: P| P::new((q.x - p.x) * cs + (q.y - p.y) * sn, -(q.x - p.x) * sn + (q.y - p.y) * cs);
    let mut w = 0;
    for seg in loops.iter().flatten() {
        match *seg {
            Seg::Line { a, b } => {
                let (a, b) = (local(a), local(b));
                if (a.y <= 0.) != (b.y <= 0.) {
                    let x = a.x + (b.x - a.x) * (-a.y) / (b.y - a.y);
                    if x > 0. {
                        w += if b.y > a.y { 1 } else { -1 };
                    }
                }
            }
            Seg::Arc { c, r, a0, sweep } => {
                let c = local(c);
                let s = -c.y / r;
                if s.abs() >= 1. {
                    continue;
                }
                for t in [s.asin(), PI - s.asin()] {
                    if c.x + r * t.cos() <= 0. {
                        continue;
                    }
                    // Back to the arc's own angle (the local frame is turned by `angle`).
                    let k = ((t + angle - a0) * sweep.signum()).rem_euclid(2. * PI) / sweep.abs();
                    if k < 1e-9 || (k - 1.).abs() < 1e-9 {
                        return None;
                    }
                    if k < 1. {
                        // Upward when the local y grows along the travel: d/dt sin = cos.
                        w += if t.cos() * sweep.signum() > 0. { 1 } else { -1 };
                    }
                }
            }
            Seg::Conic { c, u, v, a0, sweep } => {
                // Turned into the ray's frame: y(t) = cy + uy cos t + vy sin t = R cos(t - phi) = 0 crosses the ray.
                let c = local(c);
                let turn = |q: P| P::new(q.x * cs + q.y * sn, -q.x * sn + q.y * cs);
                let (u, v) = (turn(u), turn(v));
                let r = u.y.hypot(v.y);
                if r == 0. {
                    continue;
                }
                let k = -c.y / r;
                if k.abs() >= 1. {
                    continue;
                }
                let phi = v.y.atan2(u.y);
                for t in [phi + k.acos(), phi - k.acos()] {
                    if c.x + u.x * t.cos() + v.x * t.sin() <= 0. {
                        continue;
                    }
                    let along = ((t - a0) * sweep.signum()).rem_euclid(2. * PI) / sweep.abs();
                    if along < 1e-9 || (along - 1.).abs() < 1e-9 {
                        return None;
                    }
                    if along < 1. {
                        let dy = -u.y * t.sin() + v.y * t.cos();
                        w += if dy * sweep.signum() > 0. { 1 } else { -1 };
                    }
                }
            }
        }
    }
    Some(w)
}

fn winding(loops: &[Vec<Seg>], p: P) -> Result<i32> {
    for angle in [0.5172, 1.3291, 2.4127, 3.8813] {
        if let Some(w) = winding_along(loops, p, angle) {
            return Ok(w);
        }
    }
    Err(unsupported("section winding is ambiguous at a classification point"))
}

/// Orients loops so the region's winding is 1 inside and 0 outside: outer loops counter-clockwise, holes clockwise.
fn orient(loops: Vec<Vec<Seg>>) -> Result<Vec<Vec<Seg>>> {
    let mut out = Vec::with_capacity(loops.len());
    for (i, lp) in loops.iter().enumerate() {
        let probe = lp[0].at(0.5);
        let mut depth = 0;
        for (j, other) in loops.iter().enumerate() {
            if i != j && winding(std::slice::from_ref(other), probe)? != 0 {
                depth += 1;
            }
        }
        let ccw = loop_area(lp) > 0.;
        out.push(if ccw == (depth % 2 == 0) { lp.clone() } else { lp.iter().rev().map(Seg::reversed).collect() });
    }
    Ok(out)
}

struct Piece {
    v0: usize,
    v1: usize,
    seg: Seg,
}

/// Every edge split where it meets the other section, with shared vertex numbers.
fn split(own: &[Vec<Seg>], other: &[Vec<Seg>], pool: &mut Vec<P>, tol: f64) -> Result<Vec<Piece>> {
    let vertex = |q: P, pool: &mut Vec<P>| -> usize {
        if let Some(i) = pool.iter().position(|v| v.dist(q) <= tol) {
            return i;
        }
        pool.push(q);
        pool.len() - 1
    };
    let mut pieces = vec![];
    for seg in own.iter().flatten() {
        let len = seg.length();
        let mut cuts = vec![(0., vertex(seg.start(), pool)), (1., vertex(seg.end(), pool))];
        for o in other.iter().flatten() {
            for q in meetings(seg, o, tol)? {
                cuts.push((seg.param(q).clamp(0., 1.), vertex(q, pool)));
            }
        }
        cuts.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut kept: Vec<(f64, usize)> = vec![];
        for cut in cuts {
            match kept.last() {
                Some(last) if last.1 == cut.1 && (cut.0 - last.0) * len <= 2. * tol => {}
                _ => kept.push(cut),
            }
        }
        for w in kept.windows(2) {
            let ((s0, v0), (s1, v1)) = (w[0], w[1]);
            if v0 == v1 && (s1 - s0) * len <= 2. * tol {
                continue;
            }
            let piece = match *seg {
                Seg::Line { .. } => Seg::Line { a: pool[v0], b: pool[v1] },
                Seg::Arc { c, r, a0, sweep } => Seg::Arc { c, r, a0: a0 + s0 * sweep, sweep: (s1 - s0) * sweep },
                Seg::Conic { c, u, v, a0, sweep } => Seg::Conic { c, u, v, a0: a0 + s0 * sweep, sweep: (s1 - s0) * sweep },
            };
            pieces.push(Piece { v0, v1, seg: piece });
        }
    }
    Ok(pieces)
}

/// The regularized difference a - b of two oriented sections.
fn difference(a: &[Vec<Seg>], b: &[Vec<Seg>], tol: f64) -> Result<Vec<Vec<Seg>>> {
    let mut pool = vec![];
    let pa = split(a, b, &mut pool, tol)?;
    let pb = split(b, a, &mut pool, tol)?;
    let mut kept: Vec<Piece> = vec![];
    let mut b_overlap = vec![false; pb.len()];
    for p in &pa {
        let mid = p.seg.at(0.5);
        let twin = pb.iter().position(|q| {
            ((q.v0 == p.v0 && q.v1 == p.v1) || (q.v0 == p.v1 && q.v1 == p.v0)) && q.seg.at(0.5).dist(mid) <= tol
        });
        if let Some(j) = twin {
            b_overlap[j] = true;
            // Shared boundary: the same way round, b's material lies on a's side and removes it; the other way
            // round the regions only touch and a's boundary stays.
            // Direction from the tangents at the shared midpoint (a closed piece starts and ends at one vertex).
            if pb[j].seg.tangent(0.5).dot(p.seg.tangent(0.5)) < 0. {
                kept.push(Piece { v0: p.v0, v1: p.v1, seg: p.seg });
            }
            continue;
        }
        if b.iter().flatten().any(|s| s.distance(mid) <= tol) {
            return Err(unsupported("section boundaries meet along an unmatched piece"));
        }
        if winding(b, mid)? == 0 {
            kept.push(Piece { v0: p.v0, v1: p.v1, seg: p.seg });
        }
    }
    for (j, q) in pb.iter().enumerate() {
        if b_overlap[j] {
            continue;
        }
        let mid = q.seg.at(0.5);
        if a.iter().flatten().any(|s| s.distance(mid) <= tol) {
            return Err(unsupported("section boundaries meet along an unmatched piece"));
        }
        if winding(a, mid)? != 0 {
            kept.push(Piece { v0: q.v1, v1: q.v0, seg: q.seg.reversed() });
        }
    }
    if kept.len() > 4 * LIMIT {
        return Err(unsupported("section difference exceeds its work limit"));
    }
    // Chain the kept pieces into loops. Where several leave one vertex, take the first one clockwise from the way
    // back along the arriving piece: each loop then bounds one material sector.
    let mut used = vec![false; kept.len()];
    let mut loops = vec![];
    let mut visited = vec![false; pool.len()];
    for first in 0..kept.len() {
        if used[first] {
            continue;
        }
        let mut lp = vec![];
        let mut at = first;
        loop {
            used[at] = true;
            // A vertex on the boundary twice is a pinch: material meeting itself at a point of the section (an edge
            // shared by four faces in the solid), e.g. a hole tangent to the outline. Not a manifold result.
            if std::mem::replace(&mut visited[kept[at].v0], true) {
                return Err(unsupported("section cut touches a boundary at a point (non-manifold result)"));
            }
            lp.push(kept[at].seg);
            let end = kept[at].v1;
            let back = kept[at].seg.tangent(1.).scale(-1.);
            let back_angle = back.y.atan2(back.x);
            let mut best: Option<(f64, usize)> = None;
            for (k, q) in kept.iter().enumerate() {
                // The loop's own first piece competes too: choosing it closes the loop.
                if q.v0 != end || (used[k] && k != first) {
                    continue;
                }
                let out = q.seg.tangent(0.);
                let mut turn = (back_angle - out.y.atan2(out.x)).rem_euclid(2. * PI);
                if turn <= 1e-12 {
                    turn = 2. * PI;
                }
                if let Some((b, _)) = best {
                    if (turn - b).abs() <= 1e-9 {
                        return Err(unsupported("section pieces leave a vertex along one direction"));
                    }
                }
                if best.is_none_or(|(b, _)| turn < b) {
                    best = Some((turn, k));
                }
            }
            match best {
                Some((_, k)) if k == first => break,
                Some((_, k)) => at = k,
                None => return Err(unsupported("section difference does not close")),
            }
        }
        loops.push(merge(lp, tol));
    }
    Ok(loops)
}

/// Joins consecutive collinear lines and consecutive arcs of one circle, so no seam vertex remains on a straight or
/// round boundary.
fn merge(lp: Vec<Seg>, tol: f64) -> Vec<Seg> {
    let join = |s: &Seg, t: &Seg| -> Option<Seg> {
        match (*s, *t) {
            (Seg::Line { a, b }, Seg::Line { a: _, b: d }) => {
                let (u, v) = (b.sub(a), d.sub(b));
                (u.cross(v).abs() <= tol * u.len().max(v.len()) && u.dot(v) > 0.).then_some(Seg::Line { a, b: d })
            }
            (Seg::Arc { c, r, a0, sweep }, Seg::Arc { c: c2, r: r2, sweep: s2, .. }) => {
                (c.dist(c2) <= tol && (r - r2).abs() <= tol && sweep.signum() == s2.signum() && (sweep + s2).abs() < 2. * PI)
                    .then_some(Seg::Arc { c, r, a0, sweep: sweep + s2 })
            }
            (Seg::Conic { c, u, v, a0, sweep }, Seg::Conic { c: c2, u: u2, v: v2, sweep: s2, .. }) => {
                (c.dist(c2) <= tol && u.dist(u2) <= tol && v.dist(v2) <= tol && sweep.signum() == s2.signum()
                    && (sweep + s2).abs() < 2. * PI)
                    .then_some(Seg::Conic { c, u, v, a0, sweep: sweep + s2 })
            }
            _ => None,
        }
    };
    let mut out: Vec<Seg> = vec![];
    for s in lp {
        if let Some(j) = out.last().and_then(|l| join(l, &s)) {
            *out.last_mut().unwrap() = j;
        } else {
            out.push(s);
        }
    }
    while out.len() > 1 {
        match join(out.last().unwrap(), &out[0]) {
            Some(j) => {
                out[0] = j;
                out.pop();
            }
            None => break,
        }
    }
    out
}

fn axes(k: usize) -> (usize, usize) { ((k + 1) % 3, (k + 2) % 3) }
fn unit(i: usize) -> Vec3 { Vec3::new((i == 0) as u8 as f64, (i == 1) as u8 as f64, (i == 2) as u8 as f64) }
fn coord(v: Vec3, i: usize) -> f64 { [v.x, v.y, v.z][i] }

/// A right prism along axis k: its section loops (in the (u, v) plane of `axes(k)`) and its extent along k.
fn prism(body: &Body, k: usize, tol: f64) -> Result<Option<(Vec<Vec<Seg>>, f64, f64)>> {
    let (u, v) = axes(k);
    let mut caps: Vec<(usize, f64, f64)> = vec![]; // face, level, outward sign
    for (i, f) in body.faces.iter().enumerate() {
        match f.surface {
            Surface::Plane { p, n, .. } => {
                let along = coord(n, k);
                if along.abs() >= 1. - 1e-12 && coord(n, u).abs() <= 1e-12 && coord(n, v).abs() <= 1e-12 {
                    caps.push((i, coord(p, k), along.signum() * if f.sense { 1. } else { -1. }));
                } else if along.abs() > 1e-12 {
                    return Ok(None);
                }
            }
            Surface::Cylinder { a, .. } | Surface::EllipticCylinder { a, .. } => {
                if coord(a, k).abs() < 1. - 1e-12 {
                    return Ok(None);
                }
            }
        }
    }
    let lo = caps.iter().map(|c| c.1).fold(f64::INFINITY, f64::min);
    let hi = caps.iter().map(|c| c.1).fold(f64::NEG_INFINITY, f64::max);
    if caps.iter().any(|c| (c.1 - lo).abs() > tol && (c.1 - hi).abs() > tol) || hi - lo <= tol {
        return Ok(None);
    }
    // The bottom section: every bottom cap's boundary, with seams shared by two bottom faces cancelled.
    let mut uses = std::collections::BTreeMap::new();
    let bottom: Vec<_> = caps.iter().filter(|c| (c.1 - lo).abs() <= tol).collect();
    if bottom.iter().any(|c| c.2 > 0.) {
        return Ok(None);
    }
    for c in &bottom {
        for (e, forward) in body.faces[c.0].loops.iter().flatten() {
            uses.entry(*e).or_insert_with(Vec::new).push(*forward);
        }
    }
    let to2 = |q: Vec3| P::new(coord(q, u), coord(q, v));
    let mut segs = vec![];
    for (e, flags) in &uses {
        if flags.len() == 2 && flags[0] != flags[1] {
            continue;
        }
        if flags.len() != 1 {
            return Err(unsupported("prism section edge is used more than once"));
        }
        let seg = match body.edges[*e].curve {
            Curve::Line { a, b } => Seg::Line { a: to2(a), b: to2(b) },
            Curve::Arc { c, n, x, r, t0, t1 } => {
                let side = coord(n, k);
                if side.abs() < 1. - 1e-12 {
                    return Err(unsupported("prism section arc is not in the cap plane"));
                }
                let start = Curve::point_at_angle(c, n, x, r, t0);
                let center = to2(c);
                let s = to2(start).sub(center);
                Seg::Arc { c: center, r, a0: s.y.atan2(s.x), sweep: (t1 - t0) * side.signum() }
            }
            Curve::Conic { c, u: cu, v: cv, t0, t1 } => {
                if coord(cu, k).abs() > 1e-9 * cu.length() || coord(cv, k).abs() > 1e-9 * cv.length() {
                    return Err(unsupported("prism section ellipse is not in the cap plane"));
                }
                Seg::Conic { c: to2(c), u: P::new(coord(cu, u), coord(cu, v)), v: P::new(coord(cv, u), coord(cv, v)), a0: t0, sweep: t1 - t0 }
            }
            _ => return Err(unsupported("prism section edge is not a line, circular arc or ellipse")),
        };
        segs.push(if flags[0] { seg } else { seg.reversed() });
    }
    if segs.len() > LIMIT {
        return Err(unsupported("prism section exceeds its work limit"));
    }
    // Chain the boundary pieces into loops by their ends.
    let mut used = vec![false; segs.len()];
    let mut loops = vec![];
    for first in 0..segs.len() {
        if used[first] {
            continue;
        }
        used[first] = true;
        let mut lp = vec![segs[first]];
        while lp.last().unwrap().end().dist(lp[0].start()) > tol {
            let end = lp.last().unwrap().end();
            let next = (0..segs.len()).find(|j| !used[*j] && segs[*j].start().dist(end) <= tol)
                .ok_or_else(|| unsupported("prism section loop does not close"))?;
            used[next] = true;
            lp.push(segs[next]);
        }
        loops.push(merge(lp, tol));
    }
    let loops = orient(loops)?;
    // The prism proof: section area times extent is the body's measured volume.
    let area: f64 = loops.iter().map(|l| loop_area(l)).sum();
    let volume = crate::mass::measure(body, tol)?.volume;
    if (area * (hi - lo) - volume).abs() > 1e-9 * volume.abs().max(1.) {
        return Ok(None);
    }
    Ok(Some((loops, lo, hi)))
}

/// Whether the tool is made of caps normal to axis k and side faces parallel to it.
fn prism_along(tool: &Body, k: usize) -> bool {
    let (u, v) = axes(k);
    let mut caps = 0;
    tool.faces.iter().all(|f| match f.surface {
        Surface::Plane { n, .. } => {
            if coord(n, k).abs() >= 1. - 1e-12 && coord(n, u).abs() <= 1e-12 && coord(n, v).abs() <= 1e-12 {
                caps += 1;
                true
            } else {
                coord(n, k).abs() <= 1e-12
            }
        }
        Surface::Cylinder { a, .. } | Surface::EllipticCylinder { a, .. } => coord(a, k).abs() >= 1. - 1e-12,
    }) && caps >= 2
}

/// Domain selection, before anything is built (rule proposed to Codex 2026-10-29): both operands are right prisms
/// along one axis with the tool spanning the body, and either the tool's section is not strictly inside the body's
/// (it meets or crosses the outline or a hole) or the body is no longer a plain box. Clear interior cuts in a plain
/// box keep their existing domains.
pub(super) fn selects(body: &Body, tool: &Body, tol: f64) -> bool {
    for k in (0..3).filter(|&k| prism_along(tool, k) && prism_along(body, k)) {
        let (Ok(Some((a, lo, hi))), Ok(Some((b, tlo, thi)))) = (prism(body, k, tol), prism(tool, k, tol)) else {
            continue;
        };
        if tlo > lo + tol || thi < hi - tol {
            continue;
        }
        let meets = a.iter().flatten().any(|s| b.iter().flatten().any(|t| meetings(s, t, tol).map_or(true, |m| !m.is_empty())));
        // No other domain cuts with ellipses.
        let elliptic = b.iter().flatten().chain(a.iter().flatten()).any(|s| matches!(s, Seg::Conic { .. }));
        // A plain box keeps its passage domains for interior cuts; a prism already cut to another outline does not
        // fit them (Bob's p9999: a hole after its corner cut).
        let plain_box = body.faces.len() == 6
            && body.faces.iter().all(|f| matches!(f.surface, Surface::Plane { n, .. } if [n.x, n.y, n.z].iter().filter(|c| **c != 0.).count() == 1));
        return meets || elliptic || !plain_box || !matches!(winding(&a, b[0][0].start()), Ok(1));
    }
    false
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    if !tol.is_finite() || tol <= 0. {
        return Err(KernelError::InvalidInput("invalid section cut tolerance".into()));
    }
    crate::topology::check(body, tol)?;
    crate::topology::check(tool, tol)?;
    // The axis: the first along which both are right prisms and the tool spans the body (a box tool is a prism
    // along all three; the body decides). Chosen from the geometry before anything is built.
    let mut chart = None;
    for k in (0..3).filter(|&k| prism_along(tool, k)) {
        if let (Some((a, lo, hi)), Some((b, tlo, thi))) = (prism(body, k, tol)?, prism(tool, k, tol)?) {
            if tlo <= lo + tol && thi >= hi - tol {
                chart = Some((k, a, lo, hi, b));
                break;
            }
        }
    }
    let (k, a, lo, hi, b) = chart.ok_or_else(|| unsupported("section cut needs a right-prism body and a tool spanning it along one axis"))?;
    let result = difference(&a, &b, tol)?;
    if result.is_empty() {
        return Ok(Body::default());
    }
    let outer = result.iter().filter(|l| loop_area(l) > 0.).count();
    if outer != 1 {
        return Err(unsupported(if outer == 0 {
            "section difference has holes without material"
        } else {
            "section cut would leave disconnected material components"
        }));
    }
    let (u, v) = axes(k);
    let to3 = |q: P| {
        let mut w = [0.; 3];
        w[u] = q.x;
        w[v] = q.y;
        w[k] = lo;
        Vec3::new(w[0], w[1], w[2])
    };
    let curves = result
        .iter()
        .map(|lp| {
            lp.iter()
                .map(|s| match *s {
                    Seg::Line { a, b } => Curve::Line { a: to3(a), b: to3(b) },
                    Seg::Arc { c, r, a0, sweep } if sweep > 0. => {
                        Curve::Arc { c: to3(c), n: unit(k), x: unit(u), r, t0: a0, t1: a0 + sweep }
                    }
                    // Clockwise in (u, v): counter-clockwise about -k, where the angle from u runs the other way.
                    Seg::Arc { c, r, a0, sweep } => {
                        Curve::Arc { c: to3(c), n: unit(k).scale(-1.), x: unit(u), r, t0: -a0, t1: -a0 - sweep }
                    }
                    // c + cu cos t + cv sin t; run backwards it is c + cu cos s - cv sin s for s = -t increasing.
                    Seg::Conic { c, u: cu, v: cv, a0, sweep } => {
                        let along = |q: P| { let mut w = [0.; 3]; w[u] = q.x; w[v] = q.y; Vec3::new(w[0], w[1], w[2]) };
                        if sweep > 0. {
                            Curve::Conic { c: to3(c), u: along(cu), v: along(cv), t0: a0, t1: a0 + sweep }
                        } else {
                            Curve::Conic { c: to3(c), u: along(cu), v: along(cv).scale(-1.), t0: -a0, t1: -a0 - sweep }
                        }
                    }
                })
                .collect()
        })
        .collect();
    crate::extrude::extrude(&Profile { loops: curves, n: unit(k), low: 0., high: hi - lo }, tol)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mass::measure;

    const TOL: f64 = 1e-6;

    fn prism_of(loops: Vec<Vec<Curve>>, low: f64, high: f64) -> Body {
        crate::extrude::extrude(&Profile { loops, n: Vec3::new(0., 0., 1.), low, high }, TOL).unwrap()
    }
    fn polygon(points: &[(f64, f64)]) -> Vec<Curve> {
        (0..points.len())
            .map(|i| {
                let (a, b) = (points[i], points[(i + 1) % points.len()]);
                Curve::Line { a: Vec3::new(a.0, a.1, 0.), b: Vec3::new(b.0, b.1, 0.) }
            })
            .collect()
    }
    fn circle(cx: f64, cy: f64, r: f64) -> Vec<Curve> {
        vec![Curve::Arc { c: Vec3::new(cx, cy, 0.), n: Vec3::new(0., 0., 1.), x: Vec3::new(1., 0., 0.), r, t0: 0., t1: 2. * PI }]
    }
    fn plate() -> Body { prism_of(vec![polygon(&[(0., 0.), (10., 0.), (10., 8., ), (0., 8.)])], 0., 0.5) }
    fn removed(body: &Body, tool: &Body) -> f64 {
        let result = cut(body, tool, TOL).unwrap();
        crate::topology::check(&result, TOL).unwrap();
        measure(body, TOL).unwrap().volume - measure(&result, TOL).unwrap().volume
    }

    #[test]
    fn corner_triangle_with_sides_on_the_plate_faces() {
        // Bob's p9999 shape: the triangle's two legs lie on x = 10 and y = 8, exactly and one ulp off.
        for (xe, ye) in [(10., 8.), (10f64.next_down(), 8f64.next_down()), (10f64.next_up(), 8f64.next_up())] {
            let tool = prism_of(vec![polygon(&[(6., ye), (xe, 5.), (xe, ye)])], -0.001, 0.501);
            assert!((removed(&plate(), &tool) - 0.5 * 4. * 3. * 0.5).abs() < 1e-12);
        }
    }

    #[test]
    fn notch_with_an_arc_across_an_edge_both_ways_round() {
        // A D notch: x 3..5 from below the plate to y 1, then a half circle about (4, 1).
        let d = |reverse: bool| {
            let mut lp = vec![
                Curve::Line { a: Vec3::new(3., -1., 0.), b: Vec3::new(5., -1., 0.) },
                Curve::Line { a: Vec3::new(5., -1., 0.), b: Vec3::new(5., 1., 0.) },
                Curve::Arc { c: Vec3::new(4., 1., 0.), n: Vec3::new(0., 0., 1.), x: Vec3::new(1., 0., 0.), r: 1., t0: 0., t1: PI },
                Curve::Line { a: Vec3::new(3., 1., 0.), b: Vec3::new(3., -1., 0.) },
            ];
            if reverse {
                lp = lp.iter().rev().map(Curve::reversed).collect();
            }
            prism_of(vec![lp], -1., 2.)
        };
        for reverse in [false, true] {
            assert!((removed(&plate(), &d(reverse)) - (2. + PI / 2.) * 0.5).abs() < 1e-12);
        }
    }

    #[test]
    fn existing_hole_cut_across_and_cleared_holes_stay() {
        // Plate with a round hole of radius 1 at (5, 4); a slab x 4..12 cuts across half of it.
        let holed = prism_of(vec![polygon(&[(0., 0.), (10., 0.), (10., 8.), (0., 8.)]), circle(5., 4., 1.).iter().rev().map(Curve::reversed).collect()], 0., 0.5);
        let slab = prism_of(vec![polygon(&[(5., -1.), (12., -1.), (12., 9.), (5., 9.)])], -1., 1.);
        let expected = (5. * 8. - PI / 2.) * 0.5;
        assert!((removed(&holed, &slab) - expected).abs() < 1e-12);
        // A hole elsewhere is untouched: cutting the corner leaves it as it was.
        let corner = prism_of(vec![polygon(&[(8., -1.), (11., -1.), (11., 2.), (8., 2.)])], -1., 1.);
        assert!((removed(&holed, &corner) - 2. * 2. * 0.5).abs() < 1e-12);
    }

    fn fixture(name: &str) -> Body {
        let bytes = std::fs::read(format!("{}/tests/fixtures/section/{name}.f64", env!("CARGO_MANIFEST_DIR"))).unwrap();
        let values: Vec<f64> = bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
        crate::neutral::Reader::with_tolerance(&values, TOL).body().unwrap()
    }

    /// Bob's NC1 plate p9999 (Steel Structures; worker inputs from the GC Kernel view, 2026-10-29): its corner-waste
    /// triangle, legs on the plate's x and y faces one ulp inside them, removes half the legs' product times the
    /// plate thickness (6.772047244094487 x 6.853149606299212 x 0.3748031496062992 plate).
    #[test]
    fn bobs_p9999_corner_cut() {
        let (body, tool) = (fixture("p9999/body"), fixture("p9999/tool"));
        let (w, l, h) = (6.772047244094487, 6.853149606299212, 0.3748031496062992);
        let expected = 0.5 * (w - 3.5015748031496057) * (l - 3.1503937007874008) * h;
        assert!((removed(&body, &tool) - expected).abs() < 1e-12, "{}", removed(&body, &tool) - expected);
    }

    #[test]
    fn whole_removal_is_empty_and_severing_is_refused() {
        let all = prism_of(vec![polygon(&[(-1., -1.), (11., -1.), (11., 9.), (-1., 9.)])], -1., 1.);
        assert!(cut(&plate(), &all, TOL).unwrap().faces.is_empty());
        let band = prism_of(vec![polygon(&[(4., -1.), (6., -1.), (6., 9.), (4., 9.)])], -1., 1.);
        assert!(matches!(cut(&plate(), &band, TOL), Err(KernelError::Unsupported(_))));
        let shallow = prism_of(vec![polygon(&[(8., -1.), (11., -1.), (11., 2.), (8., 2.)])], 0.25, 1.);
        assert!(matches!(cut(&plate(), &shallow, TOL), Err(KernelError::Unsupported(_))));
    }
}

//! Checks that a body is a closed, consistently oriented solid. Occurrence based (Astra's review): a seam
//! edge appears in one face's loop twice, once each way, and every occurrence is counted.

use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::Body;
use crate::KernelError;

/// The counts behind the checks; Euler-Poincare is a diagnostic, not a requirement (holes, genus).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TopologyReport {
    pub vertices: usize,
    pub edges: usize,
    pub faces: usize,
    pub loops: usize,
}

impl TopologyReport {
    /// V - E + F - (L - F) = 2 (S - G) for one shell: the genus it implies.
    pub fn genus(&self) -> f64 {
        let chi = self.vertices as f64 - self.edges as f64 + self.faces as f64 - (self.loops as f64 - self.faces as f64);
        1.0 - chi / 2.0
    }
}

fn bad<T>(message: String) -> Result<T, KernelError> { Err(KernelError::InvalidInput(message)) }

/// Every edge used exactly once forward and once reversed; every loop connects head to tail and closes;
/// every edge's curve runs from its start vertex to its end vertex, and lies on the surface of every face it bounds
/// (all within `tolerance`; Astra's review: a plane moved off its boundary must be reported, never repaired).
pub fn check(body: &Body, tolerance: f64) -> Result<TopologyReport, KernelError> {
    if body.faces.is_empty() { return bad("the body has no faces".into()); }
    for (i, e) in body.edges.iter().enumerate() {
        let (s, t) = (body.vertices[e.start], body.vertices[e.end]);
        if e.curve.start().distance(s) > tolerance || e.curve.end().distance(t) > tolerance {
            return bad(format!("edge {i}'s curve does not run from its start vertex to its end vertex"));
        }
        if let Curve::Line { .. } = e.curve {
            if e.start == e.end { return bad(format!("edge {i} is a line that starts and ends at one vertex")); }
        }
    }
    let mut forward = vec![0usize; body.edges.len()];
    let mut backward = vec![0usize; body.edges.len()];
    let mut loops = 0;
    for (f, face) in body.faces.iter().enumerate() {
        if face.loops.is_empty() { return bad(format!("face {f} has no loop")); }
        for (l, lp) in face.loops.iter().enumerate() {
            loops += 1;
            if lp.is_empty() { return bad(format!("face {f} loop {l} is empty")); }
            let ends = |&(edge, fwd): &(usize, bool)| {
                let e = &body.edges[edge];
                if fwd { (e.start, e.end) } else { (e.end, e.start) }
            };
            for k in 0..lp.len() {
                let (_, tail) = ends(&lp[k]);
                let (head, _) = ends(&lp[(k + 1) % lp.len()]);
                if tail != head { return bad(format!("face {f} loop {l} breaks after segment {k}")); }
                let (edge, fwd) = lp[k];
                if fwd { forward[edge] += 1 } else { backward[edge] += 1 }
            }
        }
    }
    for i in 0..body.edges.len() {
        if forward[i] != 1 || backward[i] != 1 {
            return bad(format!("edge {i} is used {} times forward and {} times reversed (a closed solid uses each once each way)", forward[i], backward[i]));
        }
    }
    for (f, face) in body.faces.iter().enumerate() {
        for &(e, _) in face.loops.iter().flatten() {
            if let Some(miss) = off_surface(&body.edges[e].curve, &face.surface, tolerance) {
                return bad(format!("edge {e} leaves face {f}'s surface by {miss}"));
            }
        }
    }
    Ok(TopologyReport { vertices: body.vertices.len(), edges: body.edges.len(), faces: body.faces.len(), loops })
}

/// Distance from a point to a surface.
fn distance(q: Vec3, surface: &Surface) -> f64 {
    match *surface {
        Surface::Plane { p, n, .. } => q.sub(p).dot(n).abs(),
        Surface::Cylinder { p, a, r, .. } => {
            let w = q.sub(p);
            (w.sub(a.scale(w.dot(a))).length() - r).abs()
        }
    }
}

/// A bound on how far the whole curve leaves the surface, when it may exceed `tolerance` (Astra's analytic checks,
/// not samples):
/// - line or arc/ellipse c + u cos t + v sin t on a plane: |(c - p).n| + sqrt((u.n)^2 + (v.n)^2) (a line's ends);
/// - line on a cylinder: the squared radius is quadratic along it, so its ends and the clamped stationary point;
/// - arc/ellipse on a cylinder: the squared-radius residual is a trigonometric polynomial of degree two; the sum of
///   its coefficient magnitudes bounds it, against tolerance (2 r - tolerance), inward and outward;
/// - a CylinderPair lies on its two cylinders by construction: it may bound a face only on one of them.
fn off_surface(curve: &Curve, surface: &Surface, tolerance: f64) -> Option<f64> {
    let ellipse = |c: Vec3, u: Vec3, v: Vec3| (c, u, v);
    let (c, u, v) = match *curve {
        Curve::Line { a, b } => {
            let worst = match *surface {
                Surface::Plane { .. } => distance(a, surface).max(distance(b, surface)),
                Surface::Cylinder { p, a: axis, .. } => {
                    let perp = |w: Vec3| w.sub(axis.scale(w.dot(axis)));
                    let (w0, d) = (perp(a.sub(p)), perp(b.sub(a)));
                    let t = if d.dot(d) > 0.0 { (-w0.dot(d) / d.dot(d)).clamp(0.0, 1.0) } else { 0.0 };
                    [0.0, 1.0, t].iter().map(|&s| distance(a.add(b.sub(a).scale(s)), surface)).fold(0.0, f64::max)
                }
            };
            return (worst > tolerance).then_some(worst);
        }
        Curve::Arc { c, n, x, r, .. } => ellipse(c, x.scale(r), n.cross(x).scale(r)),
        Curve::Conic { c, u, v, .. } => ellipse(c, u, v),
        Curve::CylinderPair { c, a, rho, p, b, r, .. } => {
            let same = |q: Vec3, axis: Vec3, radius: f64| match *surface {
                Surface::Cylinder { p: fp, a: fa, r: fr, .. } => {
                    let w = q.sub(fp);
                    fa.cross(axis).length() <= 1e-9 && w.sub(fa.scale(w.dot(fa))).length() <= tolerance && (fr - radius).abs() <= tolerance
                }
                Surface::Plane { .. } => false,
            };
            return if same(c, a, rho) || same(p, b, r) { None } else { Some(f64::INFINITY) };
        }
    };
    let worst = match *surface {
        Surface::Plane { p, n, .. } => c.sub(p).dot(n).abs() + (u.dot(n).powi(2) + v.dot(n).powi(2)).sqrt(),
        Surface::Cylinder { p, a, r, .. } => {
            // |perp(c - p + u cos t + v sin t)|^2 - r^2 = k0 + k1 cos t + k2 sin t + k3 cos 2t + k4 sin 2t.
            let perp = |w: Vec3| w.sub(a.scale(w.dot(a)));
            let (w, pu, pv) = (perp(c.sub(p)), perp(u), perp(v));
            let k = [w.dot(w) + (pu.dot(pu) + pv.dot(pv)) / 2.0 - r * r, 2.0 * w.dot(pu), 2.0 * w.dot(pv),
                     (pu.dot(pu) - pv.dot(pv)) / 2.0, pu.dot(pv)];
            let residual: f64 = k.iter().map(|x| x.abs()).sum();
            // |rho^2 - r^2| <= tol (2 r - tol) keeps rho within tol of r on both sides.
            if r > tolerance && residual <= tolerance * (2.0 * r - tolerance) { 0.0 } else { (residual / r).max(tolerance * 2.0) }
        }
    };
    (worst > tolerance).then_some(worst)
}

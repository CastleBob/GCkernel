//! Checks that a body is a closed, consistently oriented solid. Occurrence based (Astra's review): a seam
//! edge appears in one face's loop twice, once each way, and every occurrence is counted.

use crate::geom::{Curve, Surface, Vec3};
use crate::meridian::Meridian;
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
        if let Surface::Spline { ref data } = face.surface {
            data.check_face(body, face, tolerance)?;
            continue;
        }
        for &(e, _) in face.loops.iter().flatten() {
            if let Some(miss) = off_surface(&body.edges[e].curve, &face.surface, tolerance) {
                if let Curve::Bezier { .. } | Curve::RationalBezier { .. } = body.edges[e].curve {
                    let measured = (0..=32).map(|i| distance(body.edges[e].curve.eval(i as f64 / 32.0).0, &face.surface)).fold(0.0, f64::max);
                    return bad(format!("edge {e} Bezier span cannot be certified on face {f} within {tolerance} inches; sampled departure {measured} inches (samples are not an acceptance proof)"));
                }
                return bad(format!("edge {e} leaves face {f}'s surface by {miss}"));
            }
        }
    }
    Ok(TopologyReport { vertices: body.vertices.len(), edges: body.edges.len(), faces: body.faces.len(), loops })
}

/// Distance from a point to a surface.
fn distance(q: Vec3, surface: &Surface) -> f64 {
    match *surface {
        Surface::Spline { .. } => f64::INFINITY,
        Surface::Plane { p, n, .. } => q.sub(p).dot(n).abs(),
        Surface::Cylinder { p, a, r, .. } => {
            let w = q.sub(p);
            (w.sub(a.scale(w.dot(a))).length() - r).abs()
        }
        Surface::EllipticCylinder { p, a, u, v } => {
            // Its unit-circle coordinates (c, n): scaled back onto c^2 + n^2 = 1 they give a point of the surface, so
            // this bounds the distance from above (equal on the surface).
            let (c, n) = unit_circle(p, a, u, v, q);
            let rho = c.hypot(n);
            if rho == 0.0 { return f64::INFINITY; }
            let s = q.sub(p).dot(a);
            q.distance(p.add(u.scale(c / rho)).add(v.scale(n / rho)).add(a.scale(s)))
        }
        // In the half-plane through the axis the cone is the line rho = r + k s: the distance to it.
        Surface::Cone { p, a, r, k, .. } => {
            let w = q.sub(p);
            let s = w.dot(a);
            (w.sub(a.scale(s)).length() - (r + k * s)).abs() / (1.0 + k * k).sqrt()
        }
        Surface::Sphere { c, r, .. } => (q.distance(c) - r).abs(),
        // In the half-plane through the axis the torus is the circle of radius r about (big, 0).
        Surface::Torus { c, a, big, r, .. } => {
            let w = q.sub(c);
            let h = w.dot(a);
            ((w.sub(a.scale(h)).length() - big).hypot(h) - r).abs()
        }
        // Along the axis from the surface point at the same angle and radius: an upper bound of the distance (equal
        // in the axis direction), so a check against it is never looser than the true distance.
        Surface::Conoid { c, a, x, rho, p, b, r, branch } => {
            let w = q.sub(c);
            let y = a.cross(x);
            let t = w.dot(y).atan2(w.dot(x));
            let h = crate::geom::pair_height(c, a, x, rho, p, b, r, branch, t).0;
            (w.dot(a) - h).abs()
        }
        // In the half-plane through the axis: the distance to the meridian over its whole domain.
        Surface::Revolution { p, a, m, .. } => {
            let w = q.sub(p);
            let h = w.dot(a);
            let domain = match m { Meridian::Bezier(_) => (0.0, 1.0), Meridian::Ellipse { .. } => (0.0, std::f64::consts::TAU) };
            m.distance((w.sub(a.scale(h)).length(), h), domain.0, domain.1)
        }
    }
}

/// Whether the pair (c, a, rho; p, b, r, branch) defines the conoid `surface`: the same parameter axis line (either way:
/// reversing a flips the branch, the height times a staying), the same radii and the same root cylinder axis line. The
/// parameter origin (x) and the anchors' places along their axes do not change the points.
#[allow(clippy::too_many_arguments)]
fn same_conoid(c: Vec3, a: Vec3, rho: f64, p: Vec3, b: Vec3, r: f64, branch: f64, surface: &Surface, tolerance: f64) -> bool {
    let Surface::Conoid { c: sc, a: sa, rho: srho, p: sp, b: sb, r: sr, branch: sbranch, .. } = *surface else { return false };
    let perp = |w: Vec3, d: Vec3| w.sub(d.scale(w.dot(d)));
    let turn = a.dot(sa).signum();
    a.cross(sa).length() <= 1e-9 && perp(c.sub(sc), sa).length() <= tolerance && (rho - srho).abs() <= tolerance
        && b.cross(sb).length() <= 1e-9 && perp(p.sub(sp), sb).length() <= tolerance && (r - sr).abs() <= tolerance
        && branch * turn == sbranch
}

/// A curve's (radius, height) form in the half-plane through the axis containing it, if it lies in one: the unit radial
/// direction of that half-plane and the curve's points' coordinates (affine in the curve's own data), with the largest
/// departure from the plane.
fn meridian_plane(points: &[Vec3], p: Vec3, a: Vec3) -> Option<(Vec3, f64)> {
    let perp = |w: Vec3| w.sub(a.scale(w.dot(a)));
    let far = points.iter().map(|&q| perp(q.sub(p))).max_by(|x, y| x.length().total_cmp(&y.length()))?;
    let l = far.length();
    if l == 0.0 { return None; }
    let rho = far.scale(1.0 / l);
    let n = a.cross(rho);
    let off = points.iter().map(|&q| q.sub(p).dot(n).abs()).fold(0.0, f64::max);
    Some((rho, off))
}

/// A point's coordinates (c, n) in an elliptic cylinder's frame across its axis: q - p = c u + n v + s a.
fn unit_circle(p: Vec3, a: Vec3, u: Vec3, v: Vec3, q: Vec3) -> (f64, f64) {
    Surface::elliptic_coordinates(a, u, v, q.sub(p))
}

/// A bound on how far the whole curve leaves the surface, when it may exceed `tolerance` (Astra's analytic checks,
/// not samples):
/// - line or arc/ellipse c + u cos t + v sin t on a plane: |(c - p).n| + sqrt((u.n)^2 + (v.n)^2) (a line's ends);
/// - line on a cylinder: the squared radius is quadratic along it, so its ends and the clamped stationary point;
/// - arc/ellipse on a cylinder: the squared-radius residual is a trigonometric polynomial of degree two; the sum of
///   its coefficient magnitudes bounds it, against tolerance (2 r - tolerance), inward and outward;
/// - a CylinderPair lies on its two cylinders by construction: it may bound a face only on one of them.
pub(crate) fn off_surface(curve: &Curve, surface: &Surface, tolerance: f64) -> Option<f64> {
    let ellipse = |c: Vec3, u: Vec3, v: Vec3| (c, u, v);
    let (c, u, v) = match *curve {
        Curve::Line { a, b } => {
            let worst = match *surface {
                Surface::Spline { .. } => f64::INFINITY,
                Surface::Plane { .. } => distance(a, surface).max(distance(b, surface)),
                Surface::Cylinder { p, a: axis, .. } => {
                    let perp = |w: Vec3| w.sub(axis.scale(w.dot(axis)));
                    let (w0, d) = (perp(a.sub(p)), perp(b.sub(a)));
                    let t = if d.dot(d) > 0.0 { (-w0.dot(d) / d.dot(d)).clamp(0.0, 1.0) } else { 0.0 };
                    [0.0, 1.0, t].iter().map(|&s| distance(a.add(b.sub(a).scale(s)), surface)).fold(0.0, f64::max)
                }
                // Only a generator (a line along the axis) lies on an elliptic cylinder.
                Surface::EllipticCylinder { a: axis, .. } => {
                    let d = b.sub(a);
                    let across = d.sub(axis.scale(d.dot(axis))).length();
                    across.max(distance(a, surface)).max(distance(b, surface))
                }
                // Only a generator lies on a cone: both ends on it in one half-plane through the axis. Across that half-
                // plane the line departs linearly (its far end's offset bounds it); within it the radius and the cone's
                // are both linear, so the ends bound the rest.
                // The half-plane is taken through the end farther from the axis, so a generator may end at the apex.
                Surface::Cone { p, a: axis, .. } => {
                    let perp = |w: Vec3| w.sub(axis.scale(w.dot(axis)));
                    let (mut e0, mut e1) = (perp(a.sub(p)), perp(b.sub(p)));
                    if e1.length() > e0.length() { std::mem::swap(&mut e0, &mut e1); }
                    let l0 = e0.length();
                    if l0 <= tolerance || e1.dot(e0) < -tolerance * l0 { f64::INFINITY } else {
                        let off = e1.cross(e0.scale(1.0 / l0)).length();
                        distance(a, surface).max(distance(b, surface)) + off
                    }
                }
                // No straight line lies on a sphere or (in our admitted faces) a torus or general meridian surface.
                Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } => f64::INFINITY,
                // A ruling of a conoid: square to the axis, both ends at one angle about it (the same half-plane) and on
                // the surface; the surface is straight along it, so the ends bound the rest.
                Surface::Conoid { c, a: axis, .. } => {
                    let perp = |w: Vec3| w.sub(axis.scale(w.dot(axis)));
                    let (mut e0, mut e1) = (perp(a.sub(c)), perp(b.sub(c)));
                    if e1.length() > e0.length() { std::mem::swap(&mut e0, &mut e1); }
                    let l0 = e0.length();
                    let along = b.sub(a).dot(axis).abs();
                    if l0 <= tolerance || e1.dot(e0) < -tolerance * l0 { f64::INFINITY } else {
                        let off = e1.cross(e0.scale(1.0 / l0)).length();
                        distance(a, surface).max(distance(b, surface)) + off + along
                    }
                }
            };
            return (worst > tolerance).then_some(worst);
        }
        Curve::Arc { c, n, x, r, .. } => ellipse(c, x.scale(r), n.cross(x).scale(r)),
        Curve::Conic { c, u, v, .. } => ellipse(c, u, v),
        Curve::CylinderPair { c, a, rho, p, b, r, branch, .. } => {
            let same = |q: Vec3, axis: Vec3, radius: f64| match *surface {
                Surface::Cylinder { p: fp, a: fa, r: fr, .. } => {
                    let w = q.sub(fp);
                    fa.cross(axis).length() <= 1e-9 && w.sub(fa.scale(w.dot(fa))).length() <= tolerance && (fr - radius).abs() <= tolerance
                }
                Surface::Plane { .. } | Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. }
                | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => false,
            };
            // On a conoid: the conoid of this very pair (its rim at s = rho).
            if same_conoid(c, a, rho, p, b, r, branch, surface, tolerance) { return None; }
            return if same(c, a, rho) || same(p, b, r) { None } else { Some(f64::INFINITY) };
        }
        // On the cylinder of radius s about the pair's own axis, and on the conoid of the same pair; nothing else.
        Curve::RadialPair { c, a, rho, p, b, r, branch, s, .. } => {
            let on_cylinder = match *surface {
                Surface::Cylinder { p: fp, a: fa, r: fr, .. } => {
                    let w = c.sub(fp);
                    fa.cross(a).length() <= 1e-9 && w.sub(fa.scale(w.dot(fa))).length() <= tolerance && (fr - s).abs() <= tolerance
                }
                _ => false,
            };
            return if on_cylinder || same_conoid(c, a, rho, p, b, r, branch, surface, tolerance) { None } else { Some(f64::INFINITY) };
        }
        // A Bezier curve lies within its control points' hull: on a plane, their distances bound it; on a general
        // surface of revolution with a Bezier meridian, it must be that meridian (or its reverse) in one half-plane.
        Curve::RationalBezier { p: cp, w } => {
            if crate::rational::validate(cp, w).is_err() { return Some(f64::INFINITY); }
            let worst = match *surface {
                Surface::Spline { .. } => f64::INFINITY,
                Surface::Plane { p, n, .. } => cp.iter().map(|q| { let d = q.sub(p).dot(n).abs(); if d.is_finite() { d } else { f64::INFINITY } }).fold(0.0, f64::max),
                Surface::Cylinder { p, a, r, .. } => {
                    let q = cp.map(|q| { let d = q.sub(p); d.sub(a.scale(d.dot(a))) });
                    if r > tolerance && rational_radius_bound(q, w, r*r, tolerance*(2.0*r-tolerance)) { 0.0 }
                    else { f64::INFINITY }
                },
                Surface::EllipticCylinder { p, a, u, v } => {
                    let q = cp.map(|q| { let (c, s) = unit_circle(p, a, u, v, q); Vec3::new(c, s, 0.0) });
                    let delta = tolerance / (u.dot(u)+v.dot(v)).sqrt();
                    if delta < 1.0 && rational_radius_bound(q, w, 1.0, delta*(2.0-delta)) { 0.0 }
                    else { f64::INFINITY }
                },
                _ => f64::INFINITY,
            };
            return (!worst.is_finite() || worst > tolerance).then_some(worst);
        }
        Curve::Bezier { p: cp } => {
            let worst = match *surface {
                Surface::Spline { .. } => f64::INFINITY,
                Surface::Plane { p, n, .. } => cp.iter().map(|q| q.sub(p).dot(n).abs()).fold(0.0, f64::max),
                Surface::Cylinder { p, a, r, .. } => {
                    let q = cp.map(|q| { let w = q.sub(p); w.sub(a.scale(w.dot(a))) });
                    if r > tolerance && bezier_squared_radius_bound(q, r * r, tolerance * (2.0 * r - tolerance)) { 0.0 }
                    else { f64::INFINITY }
                },
                Surface::EllipticCylinder { p, a, u, v } => {
                    let q = cp.map(|q| { let (c, s) = unit_circle(p, a, u, v, q); Vec3::new(c, s, 0.0) });
                    let sigma = (u.dot(u) + v.dot(v)).sqrt();
                    let delta = tolerance / sigma;
                    if delta < 1.0 && bezier_squared_radius_bound(q, 1.0, delta * (2.0 - delta)) { 0.0 }
                    else { f64::INFINITY }
                },
                Surface::Revolution { p, a, m: Meridian::Bezier(b), .. } => match meridian_plane(&cp, p, a) {
                    None => f64::INFINITY,
                    Some((rho, off)) => {
                        let rz: Vec<(f64, f64)> = cp.iter().map(|q| { let w = q.sub(p); (w.dot(rho), w.dot(a)) }).collect();
                        let gap = |order: [usize; 4]| (0..4).map(|k| (rz[k].0 - b[order[k]].0).hypot(rz[k].1 - b[order[k]].1)).fold(0.0, f64::max);
                        off + gap([0, 1, 2, 3]).min(gap([3, 2, 1, 0]))
                    }
                },
                _ => f64::INFINITY,
            };
            return (worst > tolerance).then_some(worst);
        }
    };
    let worst = match *surface {
        Surface::Spline { .. } => f64::INFINITY,
        // No circle or ellipse lies on a conoid whose height varies (the admitted faces are bounded by rulings and rims).
        Surface::Conoid { .. } => f64::INFINITY,
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
        Surface::EllipticCylinder { p, a, u: su, v: sv } => {
            // In the surface's unit-circle coordinates the curve is A + B cos t + C sin t (linear in the point), so
            // c^2 + n^2 - 1 = k0 + k1 cos t + k2 sin t + k3 cos 2t + k4 sin 2t; the sum of |k| bounds it. A residual e
            // moves a point at most max(|u|, |v|) |sqrt(1 + e) - 1| <= max(|u|, |v|) |e| / 2 off the surface. The curve
            // must also lie across the axis (a section at one height): its semi-diameters' parts along the axis count.
            let at = |q: Vec3| unit_circle(p, a, su, sv, q);
            let lin = |w: Vec3| unit_circle(Vec3::new(0.0, 0.0, 0.0), a, su, sv, w);
            let (big_a, big_b, big_c) = (at(c), lin(u), lin(v));
            let dot = |x: (f64, f64), y: (f64, f64)| x.0 * y.0 + x.1 * y.1;
            let k = [dot(big_a, big_a) + (dot(big_b, big_b) + dot(big_c, big_c)) / 2.0 - 1.0, 2.0 * dot(big_a, big_b),
                     2.0 * dot(big_a, big_c), (dot(big_b, big_b) - dot(big_c, big_c)) / 2.0, dot(big_b, big_c)];
            let residual: f64 = k.iter().map(|x| x.abs()).sum();
            let tilt = u.dot(a).abs().max(v.dot(a).abs());
            residual * su.length().max(sv.length()) / 2.0 + tilt
        }
        // A circle across the axis at height s = (c - p).a, radius r + k s there: the cylinder's residual at that radius,
        // plus the curve's tilt along the axis (which moves both height and required radius).
        Surface::Cone { p, a, r: r0, k, .. } => {
            let r = r0 + k * c.sub(p).dot(a);
            if r <= tolerance { f64::INFINITY } else {
                let perp = |w: Vec3| w.sub(a.scale(w.dot(a)));
                let (w, pu, pv) = (perp(c.sub(p)), perp(u), perp(v));
                let k2 = [w.dot(w) + (pu.dot(pu) + pv.dot(pv)) / 2.0 - r * r, 2.0 * w.dot(pu), 2.0 * w.dot(pv),
                          (pu.dot(pu) - pv.dot(pv)) / 2.0, pu.dot(pv)];
                let residual: f64 = k2.iter().map(|x| x.abs()).sum();
                let tilt = u.dot(a).abs().max(v.dot(a).abs());
                let radial = if r > tolerance && residual <= tolerance * (2.0 * r - tolerance) { 0.0 } else { (residual / r).max(tolerance * 2.0) };
                radial + tilt * (1.0 + k.abs())
            }
        }
        // Any circle on a sphere: |X - c|^2 - r^2 is a degree-two trigonometric polynomial; the sum of its coefficients'
        // magnitudes bounds it, and a residual e moves a point at most |e| / (2 r) (for small e) off the sphere.
        Surface::Sphere { c: sc, r, .. } => {
            let w = c.sub(sc);
            let k = [w.dot(w) + (u.dot(u) + v.dot(v)) / 2.0 - r * r, 2.0 * w.dot(u), 2.0 * w.dot(v),
                     (u.dot(u) - v.dot(v)) / 2.0, u.dot(v)];
            let residual: f64 = k.iter().map(|x| x.abs()).sum();
            if residual <= tolerance * (2.0 * r - tolerance) { 0.0 } else { (residual / r).max(tolerance * 2.0) }
        }
        // A torus holds two families of circles our faces use: across the axis (a parallel: centre on the axis, its
        // radius and height on the meridian circle) and in a plane through the axis (a meridian: centre on the centre
        // circle, radius r). Anything else is refused.
        Surface::Revolution { p, a, m, .. } => {
            let perp = |w: Vec3| w.sub(a.scale(w.dot(a)));
            let w = c.sub(p);
            let tilt = u.dot(a).abs().max(v.dot(a).abs());
            if tilt <= tolerance && perp(w).length() <= tolerance {
                // A parallel: a true circle about the axis whose (radius, height) is on the meridian.
                let (pu, pv) = (perp(u), perp(v));
                let rho = (pu.length() + pv.length()) / 2.0;
                let roundness = (pu.length() - pv.length()).abs() + pu.dot(pv).abs() / rho.max(1e-300);
                let domain = match m { Meridian::Bezier(_) => (0.0, 1.0), Meridian::Ellipse { .. } => (0.0, std::f64::consts::TAU) };
                m.distance((rho, w.dot(a)), domain.0, domain.1) + roundness + perp(w).length()
            } else if let Meridian::Ellipse { c: mc, u: mu, v: mv } = m {
                // An ellipse meridian: the conic lies in a half-plane through the axis and on the meridian ellipse
                // (its points in the meridian's unit-circle coordinates satisfy c^2 + n^2 = 1 identically).
                match meridian_plane(&[c, c.add(u), c.sub(u), c.add(v), c.sub(v)], p, a) {
                    None => f64::INFINITY,
                    Some((rho, off)) => {
                        let det = mu.0 * mv.1 - mu.1 * mv.0;
                        let unit = |x: f64, y: f64| ((x * mv.1 - y * mv.0) / det, (mu.0 * y - mu.1 * x) / det);
                        let (wr, wz) = (w.dot(rho) - mc.0, w.dot(a) - mc.1);
                        let (big_a, big_b, big_c) = (unit(wr, wz), unit(u.dot(rho), u.dot(a)), unit(v.dot(rho), v.dot(a)));
                        let dot = |x: (f64, f64), y: (f64, f64)| x.0 * y.0 + x.1 * y.1;
                        let k = [dot(big_a, big_a) + (dot(big_b, big_b) + dot(big_c, big_c)) / 2.0 - 1.0, 2.0 * dot(big_a, big_b),
                                 2.0 * dot(big_a, big_c), (dot(big_b, big_b) - dot(big_c, big_c)) / 2.0, dot(big_b, big_c)];
                        let residual: f64 = k.iter().map(|x| x.abs()).sum();
                        off + residual * mu.0.hypot(mu.1).max(mv.0.hypot(mv.1)) / 2.0
                    }
                }
            } else { f64::INFINITY }
        }
        Surface::Torus { c: tc, a, x: _, big, r } => {
            let perp = |w: Vec3| w.sub(a.scale(w.dot(a)));
            let w = c.sub(tc);
            let tilt = u.dot(a).abs().max(v.dot(a).abs());
            if tilt <= tolerance && perp(w).length() <= tolerance {
                // A parallel at height h with radius rho: (rho - big)^2 + h^2 = r^2, and a true circle about the axis.
                let (pu, pv) = (perp(u), perp(v));
                let rho = (pu.length() + pv.length()) / 2.0;
                let roundness = (pu.length() - pv.length()).abs() + pu.dot(pv).abs() / rho.max(1e-300);
                let h = w.dot(a);
                ((rho - big).hypot(h) - r).abs() + roundness + perp(w).length()
            } else {
                // A meridian: the circle's plane contains the axis, its centre on the centre circle, radius r.
                let n = u.cross(v);
                let ln = n.length();
                if ln == 0.0 { f64::INFINITY } else {
                    let n = n.scale(1.0 / ln);
                    let plane_off = n.dot(a).abs() * w.length().max(r) + n.dot(w).abs();
                    let centre_off = (perp(w).length() - big).abs() + w.dot(a).abs();
                    let round = (u.length() - r).abs() + (v.length() - r).abs() + u.dot(v).abs() / r;
                    plane_off + centre_off + round
                }
            }
        }
    };
    (worst > tolerance).then_some(worst)
}

/// Bound |q(t).q(t) - radius_squared| on the ENTIRE cubic span. Products of
/// Bernstein polynomials give a degree-six Bernstein polynomial whose range is
/// inside its coefficients' hull. Subdivision tightens that bound without
/// changing the stored curve. This is admission, never fitting or projection.
/// A failed bound (including work/roundoff limits) refuses the curve.
fn bezier_squared_radius_bound(q: [Vec3; 4], radius_squared: f64, allowed: f64) -> bool {
    let c3 = [1.0, 3.0, 3.0, 1.0];
    let c6 = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    let mut pending = vec![(q, 0)];
    let mut visited = 0;
    while let Some((q, depth)) = pending.pop() {
        visited += 1;
        if visited > 8192 {
            return false;
        }
        let mut worst = 0.0f64;
        for k in 0..=6usize {
            let mut sum = 0.0;
            let mut magnitude = radius_squared.abs();
            for i in k.saturating_sub(3)..=k.min(3) {
                let j = k - i;
                let term = c3[i] * c3[j] / c6[k] * q[i].dot(q[j]);
                sum += term;
                magnitude += term.abs();
            }
            let bound = (sum - radius_squared).abs() + 128.0 * f64::EPSILON * magnitude;
            if !bound.is_finite() {
                return false;
            }
            worst = worst.max(bound);
        }
        if worst <= allowed {
            continue;
        }
        if depth >= 24 {
            return false;
        }
        // Actual endpoints/midpoint are rejection witnesses, never acceptance samples.
        let left = crate::geom::bezier_piece(q, 0.0, 0.5);
        let right = crate::geom::bezier_piece(q, 0.5, 1.0);
        if [q[0], q[3], left[3]]
            .iter()
            .any(|p| (p.dot(*p) - radius_squared).abs() > allowed)
        {
            return false;
        }
        pending.push((right, depth + 1));
        pending.push((left, depth + 1));
    }
    true
}

/// For positive weights, |Q.Q - r^2 W^2| / min(W)^2 bounds the radial
/// residual over a whole span. Products use degree-six Bernstein coefficients.
/// Subdivision tightens the bound; samples may reject, never admit.
fn rational_radius_bound(p: [Vec3; 4], w: [f64; 4], r2: f64, allowed: f64) -> bool {
    let max = w.iter().copied().fold(0.0, f64::max);
    let mut pending = vec![(p, w.map(|v| v/max), 0)];
    let (c3, c6) = ([1., 3., 3., 1.], [1., 6., 15., 20., 15., 6., 1.]);
    let mut visits = 0;
    while let Some((p, w, depth)) = pending.pop() {
        visits += 1;
        if visits > 8192 { return false; }
        let q: [Vec3; 4] = std::array::from_fn(|i| p[i].scale(w[i]));
        let min_w = w.iter().copied().fold(f64::INFINITY, f64::min);
        let mut worst = 0.0f64;
        for k in 0usize..=6 {
            let (mut sum, mut magnitude) = (0.0, 0.0);
            for i in k.saturating_sub(3)..=k.min(3) {
                let j = k-i;
                let factor = c3[i]*c3[j]/c6[k];
                let (a, b) = (q[i].dot(q[j]), r2*w[i]*w[j]);
                sum += factor*(a-b);
                magnitude += factor*((q[i].x*q[j].x).abs()+(q[i].y*q[j].y).abs()+(q[i].z*q[j].z).abs()+b.abs());
            }
            let bound = (sum.abs()+256.0*f64::EPSILON*magnitude)/(min_w*min_w);
            if !bound.is_finite() { return false; }
            worst = worst.max(bound);
        }
        if worst <= allowed { continue; }
        if depth >= 24 { return false; }
        let [left, right] = crate::rational::split(p, w, 0.5);
        if [p[0], p[3], left.0[3]].iter().any(|q| (q.dot(*q)-r2).abs() > allowed) { return false; }
        pending.push((right.0, right.1, depth+1));
        pending.push((left.0, left.1, depth+1));
    }
    true
}

#[cfg(test)]
mod spline_boundary_tests {
    use super::*;
    #[test]
    fn whole_span_check_rejects_bulge_hidden_by_endpoints_and_midpoint() {
        let p = [
            Vec3::new(2., 0., 0.),
            Vec3::new(2.1, 0., 1.),
            Vec3::new(1.9, 0., 2.),
            Vec3::new(2., 0., 3.),
        ];
        let curve = Curve::Bezier { p };
        let surface = Surface::Cylinder {
            p: Vec3::new(0., 0., 0.),
            a: Vec3::new(0., 0., 1.),
            x: Vec3::new(1., 0., 0.),
            r: 2.,
        };
        for t in [0., 0.5, 1.] {
            assert!(distance(curve.eval(t).0, &surface) < 1e-12);
        }
        assert!(off_surface(&curve, &surface, 1e-6).is_some());
        let straight = Curve::Bezier {
            p: p.map(|q| Vec3::new(2., 0., q.z)),
        };
        assert!(off_surface(&straight, &surface, 1e-6).is_none());
    }
    #[test]
    fn rational_hull_and_denominator_bounds_reject_hidden_departures() {
        let origin = Vec3::new(0.,0.,0.);
        let axis = Vec3::new(0.,0.,1.);
        let surface = Surface::Cylinder { p: origin, a: axis, x: Vec3::new(1.,0.,0.), r: 2. };
        let p = [Vec3::new(2.,0.,0.), Vec3::new(2.1,0.,1.), Vec3::new(1.9,0.,2.), Vec3::new(2.,0.,3.)];
        let curve = Curve::RationalBezier { p, w: [1.;4] };
        for t in [0.,0.5,1.] { assert!(distance(curve.eval(t).0,&surface) < 1e-12); }
        assert!(off_surface(&curve,&surface,1e-6).is_some());
        let fixture = crate::step::read(include_bytes!("../tests/fixtures/rational_cylinder.step"),1e-6).unwrap();
        let (p,w) = fixture.bodies[0].body.edges.iter().find_map(|e| match e.curve {
            Curve::RationalBezier { p,w } => Some((p,w)), _ => None }).unwrap();
        // Affinely stretch an exact circular rational quadrant into an exact ellipse.
        let p = p.map(|p| Vec3::new(1.5*p.x,0.5*p.y,p.z));
        let ellipse = Surface::EllipticCylinder { p: origin, a: axis, u: Vec3::new(3.,0.,0.), v: Vec3::new(0.,1.,0.) };
        assert!(off_surface(&Curve::RationalBezier { p,w },&ellipse,1e-6).is_none());
        let mut bad = p; bad[1].x += 1e-4;
        assert!(off_surface(&Curve::RationalBezier { p: bad,w },&ellipse,1e-6).is_some());
        let plane = Surface::Plane { p: origin, n: axis, x: Vec3::new(1.,0.,0.) };
        assert!(off_surface(&Curve::RationalBezier { p,w },&plane,1e-6).is_none());
        bad = p; bad[1].z += 1e-4;
        assert!(off_surface(&Curve::RationalBezier { p: bad,w },&plane,1e-6).is_some());
    }
    #[test]
    fn elliptic_cylinder_admission_uses_actual_axes_and_rejects_radial_error() {
        let p = [
            Vec3::new(3., 0., 0.),
            Vec3::new(3., 0., 1.),
            Vec3::new(3., 0., 2.),
            Vec3::new(3., 0., 3.),
        ];
        let surface = Surface::EllipticCylinder {
            p: Vec3::new(0., 0., 0.),
            a: Vec3::new(0., 0., 1.),
            u: Vec3::new(3., 0., 0.),
            v: Vec3::new(0., 1., 0.),
        };
        assert!(off_surface(&Curve::Bezier { p }, &surface, 1e-6).is_none());
        let mut bad = p;
        bad[1].x += 0.0001;
        assert!(off_surface(&Curve::Bezier { p: bad }, &surface, 1e-6).is_some());
    }
}

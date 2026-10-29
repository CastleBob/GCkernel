//! Exact volume and surface area of a neutral body from its boundary (Claude 2026-10-27; an independent check for
//! every operation, Astra's acceptance request). Divergence theorem: volume = 1/3 of the flux of x through the faces.
//!
//! Plane face: x . n is the constant p . n on it, so its share is (p . n) area / 3, the area by Green's theorem over
//! its loops (extrude::signed_area, exact for lines, arcs and conics).
//!
//! Cylinder face (axis point P, axis a, frame X, Y = a x X, radius r): x = P + s a + r rho(theta), outward normal
//! sense * rho. With f(theta) = r (x . rho) = r (P.X cos theta + P.Y sin theta + r), Green's theorem in (theta, s)
//! gives the flux as -(loop integral of s f dtheta) and the area as -sense (loop integral of s r dtheta), whatever
//! the face's sense, because the loops run counter-clockwise about the outward normal. dtheta is a closed form on
//! the whole cylinder, so seams need no bookkeeping. Rulings (lines along the axis) contribute nothing; circles and
//! ellipse arcs on the cylinder are c + u cos t + v sin t with theta = theta_u + e t, so every integral is closed form.
//! A CylinderPair edge has no closed form: its integrals use adaptive Gauss-Legendre quadrature (smooth integrand,
//! converged to about 1e-14 relative), and the result says so (`exact` false).

use crate::extrude::signed_area;
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::Body;
use crate::KernelError;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassProperties {
    pub volume: f64,
    pub area: f64,
    /// False when any integral was by quadrature (CylinderPair edges) rather than closed form.
    pub exact: bool,
}

fn gauss_nodes() -> &'static [(f64, f64)] {
    static NODES: std::sync::OnceLock<Vec<(f64, f64)>> = std::sync::OnceLock::new();
    NODES.get_or_init(|| {
        // Roots and weights of the degree-16 Legendre polynomial, by Newton's method.
        let n = 16;
        (0..n).map(|i| {
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let mut dp = 0.0;
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, x);
                for k in 2..=n { let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64; p0 = p1; p1 = p2; }
                dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
                let dx = p1 / dp;
                x -= dx;
                if dx.abs() < 1e-16 { break; }
            }
            (x, 2.0 / ((1.0 - x * x) * dp * dp))
        }).collect()
    })
}

fn gauss(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let (m, h) = ((a + b) / 2.0, (b - a) / 2.0);
    gauss_nodes().iter().map(|&(x, w)| w * f(m + h * x)).sum::<f64>() * h
}

/// Integral of a smooth f over [a, b]: 16-point Gauss-Legendre, a panel halved until its halves agree with it within
/// `absolute` (shared out by length) or 1e-13 relative. At most 4096 panels: past that the integral is refused, not
/// estimated (Astra).
pub fn quadrature(f: &dyn Fn(f64) -> f64, a: f64, b: f64, absolute: f64) -> Result<f64, String> {
    const PANELS: usize = 4096;
    let mut pending = vec![(a, b, gauss(f, a, b))];
    let (mut total, mut panels) = (0.0, 1usize);
    while let Some((lo, hi, whole)) = pending.pop() {
        let m = (lo + hi) / 2.0;
        let (left, right) = (gauss(f, lo, m), gauss(f, m, hi));
        let allowed = (absolute * (hi - lo) / (b - a)).max(1e-13 * (left.abs() + right.abs()));
        if (left + right - whole).abs() <= allowed {
            total += left + right;
        } else {
            panels += 1;
            if panels > PANELS { return Err(format!("quadrature did not converge within {PANELS} panels")); }
            pending.push((lo, m, left));
            pending.push((m, hi, right));
        }
    }
    Ok(total)
}

/// A loop segment moved by `shift` (to the measuring frame), run the way the loop runs.
fn segment(body: &Body, edge: usize, forward: bool, shift: Vec3) -> Curve {
    let c = body.edges[edge].curve.translated(shift);
    if forward { c } else { c.reversed() }
}

/// Volume and area; a body whose faces point inwards has negative volume. Everything is measured in a frame at the
/// vertices' centroid: each curve and surface is moved there once, so a small part far from the world origin loses
/// nothing to cancellation (Astra). The body itself is not changed.
/// Curves that do not lie on their cylinder (beyond `tolerance`) are refused.
pub fn measure(body: &Body, tolerance: f64) -> Result<MassProperties, KernelError> {
    let (mut volume, mut area) = (0.0, 0.0);
    let exact = !body.edges.iter().any(|e| matches!(e.curve, Curve::CylinderPair { .. }))
        && !body.faces.iter().any(|f| matches!(f.surface, Surface::EllipticCylinder { .. }));
    let n = body.vertices.len().max(1) as f64;
    let o = body.vertices.iter().fold(Vec3::new(0.0, 0.0, 0.0), |s, v| s.add(*v)).scale(1.0 / n);
    // The body's size, for the quadrature's absolute tolerance.
    let size = body.vertices.iter().map(|v| v.distance(o)).fold(0.0, f64::max).max(1e-9);
    let shift = o.scale(-1.0);
    let o = Vec3::new(0.0, 0.0, 0.0);
    for (f, face) in body.faces.iter().enumerate() {
        let sense = if face.sense { 1.0 } else { -1.0 };
        match &face.surface {
            Surface::Plane { p, n, .. } => {
                let outward = n.scale(sense);
                let a: f64 = face.loops.iter()
                    .map(|lp| signed_area(&lp.iter().map(|&(e, fwd)| segment(body, e, fwd, shift)).collect::<Vec<_>>(), outward))
                    .sum();
                area += a;
                volume += p.add(shift).dot(outward) * a / 3.0;
            }
            Surface::Cylinder { p, a, x, r } => {
                let p = &p.add(shift);
                let y = a.cross(*x);
                let (px, py) = (p.sub(o).dot(*x), p.sub(o).dot(y));
                let (mut flux, mut sweep) = (0.0, 0.0);
                for &(e, fwd) in face.loops.iter().flatten() {
                    let (j, k) = cylinder_terms(&segment(body, e, fwd, shift), *p, o, *a, *x, y, *r, px, py, tolerance, size)
                        .map_err(|m| KernelError::InvalidInput(format!("face {f}, edge {e}: {m}")))?;
                    flux += j;
                    sweep += k;
                }
                volume -= flux / 3.0;
                area -= sense * sweep;
            }
            Surface::EllipticCylinder { p, a, u, v } => {
                let p = p.add(shift);
                let (mut flux, mut sweep) = (0.0, 0.0);
                for &(e, fwd) in face.loops.iter().flatten() {
                    let (j, k) = elliptic_terms(&segment(body, e, fwd, shift), p, *a, *u, *v, size)
                        .map_err(|m| KernelError::InvalidInput(format!("face {f}, edge {e}: {m}")))?;
                    flux += j;
                    sweep += k;
                }
                volume -= flux / 3.0;
                area -= sense * sweep;
            }
        }
    }
    if !volume.is_finite() || !area.is_finite() {
        return Err(KernelError::Unsupported("a boundary integral could not be evaluated".into()));
    }
    Ok(MassProperties { volume, area, exact })
}

/// For one boundary curve on the cylinder: (integral of s f dtheta, integral of s r dtheta). s is measured along the
/// axis from the body's centroid `o`, not from the axis point `p`, which may sit far away (Astra): f dtheta and
/// r dtheta are closed forms, so a constant added to s changes no face's loop integrals.
#[allow(clippy::too_many_arguments)]
fn cylinder_terms(curve: &Curve, p: Vec3, o: Vec3, a: Vec3, x: Vec3, y: Vec3, r: f64, px: f64, py: f64, tol: f64, size: f64)
    -> Result<(f64, f64), String> {
    let (c, u, v, t0, t1) = match *curve {
        Curve::Line { a: s, b: e } => {
            let d = e.sub(s);
            if d.cross(a).length() > tol * d.length().max(1.0) { return Err("a line on a cylinder must run along its axis".into()); }
            return Ok((0.0, 0.0));
        }
        Curve::Arc { c, n, x: ax, r: ar, t0, t1 } => (c, ax.scale(ar), n.cross(ax).scale(ar), t0, t1),
        Curve::Conic { c, u, v, t0, t1 } => (c, u, v, t0, t1),
        Curve::CylinderPair { .. } => return curve_on_cylinder_terms(curve, p, o, a, x, y, r, px, py, tol, size),
    };
    let perp = |w: Vec3| w.sub(a.scale(w.dot(a)));
    let (dc, up, vp) = (c.sub(p), perp(u), perp(v));
    if perp(dc).length() > tol || (up.length() - r).abs() > tol || (vp.length() - r).abs() > tol || up.dot(vp).abs() > tol * r {
        return Err("a circle or ellipse that does not lie on its cylinder".into());
    }
    let e = up.cross(vp).dot(a).signum();
    let theta_u = up.dot(y).atan2(up.dot(x));
    let (ca, sa) = (theta_u.cos(), theta_u.sin());
    let (big_a, big_b) = (px * ca + py * sa, e * (py * ca - px * sa));
    let (sc, cu, cv) = (c.sub(o).dot(a), u.dot(a), v.dot(a));
    let span = t1 - t0;
    let ic = t1.sin() - t0.sin();
    let is = t0.cos() - t1.cos();
    let s2 = (2.0 * t1).sin() - (2.0 * t0).sin();
    let (icc, iss, ics) = (span / 2.0 + s2 / 4.0, span / 2.0 - s2 / 4.0, (t1.sin().powi(2) - t0.sin().powi(2)) / 2.0);
    let j = e * r * (sc * big_a * ic + sc * big_b * is + sc * r * span + cu * big_a * icc + (cu * big_b + cv * big_a) * ics
        + cu * r * ic + cv * big_b * iss + cv * r * is);
    let k = e * r * (sc * span + cu * ic + cv * is);
    Ok((j, k))
}

/// For one boundary curve on an elliptic cylinder X = C(t) + s a, C(t) = p + u cos t + v sin t (Claude 2026-10-29): the
/// flux integrand X.(C' x a) = C(t).(C'(t) x a) =: g(t) does not depend on s, so (as for the cylinder) the face's flux is
/// the loop integral of s g(t) dt and its area that of s |C'(t)| dt. Along the curve t is recovered exactly from the
/// point's coordinates (c, n) across the axis, dt = (c dn - n dc) / (c^2 + n^2); both integrals by quadrature.
fn elliptic_terms(curve: &Curve, p: Vec3, a: Vec3, u: Vec3, v: Vec3, size: f64) -> Result<(f64, f64), String> {
    let dual = |w: Vec3| Surface::elliptic_coordinates(a, u, v, w);
    let (t0, t1) = curve.domain();
    let parts = |tau: f64| {
        let (q, dq) = curve.eval(tau);
        let s = q.dot(a);
        let (c, n) = dual(q.sub(p));
        let (dc, dn) = dual(dq);
        let rho2 = c * c + n * n;
        let t = n.atan2(c);
        let dt = if rho2 > 0.0 { (c * dn - n * dc) / rho2 } else { 0.0 };
        let cp = u.scale(-t.sin()).add(v.scale(t.cos()));
        let ct = p.add(u.scale(t.cos())).add(v.scale(t.sin()));
        (s, ct.dot(cp.cross(a)), cp.length(), dt)
    };
    let j = quadrature(&|tau| { let (s, g, _, dt) = parts(tau); s * g * dt }, t0, t1, 1e-14 * size.powi(3))?;
    let k = quadrature(&|tau| { let (s, _, len, dt) = parts(tau); s * len * dt }, t0, t1, 1e-14 * size.powi(2))?;
    Ok((j, k))
}

/// The same two integrals for any curve on the cylinder, by quadrature: with w = point - P, s = (point - o) . a,
/// cos theta = w . x / r, sin theta = w . y / r and dtheta/dt = (w.x (w'.y) - w.y (w'.x)) / r^2.
#[allow(clippy::too_many_arguments)]
fn curve_on_cylinder_terms(curve: &Curve, p: Vec3, o: Vec3, a: Vec3, x: Vec3, y: Vec3, r: f64, px: f64, py: f64, tol: f64,
    size: f64) -> Result<(f64, f64), String> {
    let (t0, t1) = curve.domain();
    for i in 0..=8 {
        let w = curve.eval(t0 + (t1 - t0) * i as f64 / 8.0).0.sub(p);
        if (w.sub(a.scale(w.dot(a))).length() - r).abs() > tol { return Err("a curve that does not lie on its cylinder".into()); }
    }
    // w . x as (q - o) . x - (p - o) . x: the far axis point enters once, as px and py, not at every t.
    let parts = |t: f64| {
        let (near, d) = curve.eval_from(o, t);
        let (wx, wy) = (near.dot(x) - px, near.dot(y) - py);
        (near.dot(a), wx, wy, (wx * d.dot(y) - wy * d.dot(x)) / (r * r))
    };
    let j = quadrature(&|t| { let (s, wx, wy, dtheta) = parts(t); s * r * (px * wx / r + py * wy / r + r) * dtheta }, t0, t1,
        1e-14 * size.powi(3))?;
    let k = quadrature(&|t| { let (s, _, _, dtheta) = parts(t); s * r * dtheta }, t0, t1, 1e-14 * size.powi(2))?;
    Ok((j, k))
}

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
//!
//! Radial conoid face (Claude 2026-10-31): X = C + s rho(t) + h(t) a, natural normal N = X_t x X_s = -s a + h'(t) v(t),
//! v = a x rho. X.N = -s (C.a + h) + h' (C.v), with antiderivative in s G = -(s^2 / 2)(C.a + h) + s h' (C.v); as on the
//! cylinder, Green's theorem in (t, s) makes the outward flux -(loop integral of G dt) whatever the face's sense, and
//! |N| = sqrt(s^2 + h'^2) has antiderivative A = (s / 2) sqrt(s^2 + h'^2) + (h'^2 / 2) asinh(s / |h'|), the area
//! -sense (loop integral of A dt). Rulings have dt = 0; every other boundary curve by quadrature in its own parameter.

use crate::extrude::signed_area;
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::Body;
use crate::KernelError;
use std::f64::consts::{PI, TAU};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassProperties {
    pub volume: f64,
    pub area: f64,
    /// False when any integral was by quadrature (CylinderPair edges) rather than closed form.
    pub exact: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FaceArea {
    pub area: f64,
    pub exact: bool,
}

/// Areas in body face order, using the same trimmed-surface integrals as whole-body mass.
/// As with `measure`, the caller must validate the body topology first.
pub fn face_areas(body: &Body, tolerance: f64) -> Result<Vec<FaceArea>, KernelError> {
    let mut faces = Vec::with_capacity(body.faces.len());
    measure_impl(body, tolerance, Some(&mut faces))?;
    Ok(faces)
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
    measure_impl(body, tolerance, None)
}

fn measure_impl(body: &Body, tolerance: f64, mut faces: Option<&mut Vec<FaceArea>>) -> Result<MassProperties, KernelError> {
    let (mut total_volume, mut total_area) = (0.0, 0.0);
    let exact = !body.edges.iter().any(|e| matches!(e.curve, Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } | Curve::RationalBezier { .. }))
        && !body.faces.iter().any(|f| matches!(f.surface, Surface::EllipticCylinder { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. })
            || matches!(f.surface, Surface::Cylinder { .. })
                && f.loops.iter().flatten().any(|&(e, _)| matches!(body.edges[e].curve, Curve::Bezier { .. })));
    let n = body.vertices.len().max(1) as f64;
    let o = body.vertices.iter().fold(Vec3::new(0.0, 0.0, 0.0), |s, v| s.add(*v)).scale(1.0 / n);
    // The body's size, for the quadrature's absolute tolerance.
    let size = body.vertices.iter().map(|v| v.distance(o)).fold(0.0, f64::max).max(1e-9);
    let shift = o.scale(-1.0);
    let o = Vec3::new(0.0, 0.0, 0.0);
    for (f, face) in body.faces.iter().enumerate() {
        let (mut volume, mut area) = (0.0, 0.0);
        let sense = if face.sense { 1.0 } else { -1.0 };
        match &face.surface {
            Surface::Spline { data } => {
                data.check_face(body, face, tolerance)?;
                let (v, a) = data.measure_face(body, face, tolerance, shift, sense, size)?;
                volume += v;
                area += a;
            }
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
            Surface::Cone { p, a, x, r, k } => {
                let p = p.add(shift);
                let y = a.cross(*x);
                let (mut flux, mut sweep) = (0.0, 0.0);
                for &(e, fwd) in face.loops.iter().flatten() {
                    let (j, s) = cone_terms(&segment(body, e, fwd, shift), p, *a, *x, y, *r, *k, tolerance)
                        .map_err(|m| KernelError::InvalidInput(format!("face {f}, edge {e}: {m}")))?;
                    flux += j;
                    sweep += s;
                }
                // An apex is a boundary segment of the (theta, s) chart that no edge carries: its share is added here.
                for (q, ta, tb) in pole_jumps(body, face, shift, p, *a, *x, y, sense, tolerance) {
                    let (j, s) = cone_parallel(q.sub(p).dot(*a), ta, tb, p, *a, *x, y, *r, *k);
                    flux += j;
                    sweep += s;
                }
                volume -= flux / 3.0;
                area -= sense * sweep;
            }
            Surface::Revolution { p, a, x, m } => {
                // X = P + z(u) a + r(u) rho(theta): X.N = r z' (P.rho) + r (r z' - r' (z + P.a)) (meridian.rs flux_hats),
                // |N| = r |M'|. Parallels have u constant (closed-form flux, area by quadrature); meridians theta.
                let p = p.add(shift);
                let y = a.cross(*x);
                let pa = p.dot(*a);
                let edges: Vec<(usize, bool)> = face.loops.iter().flatten().copied().collect();
                let rz = |q: Vec3| { let w = q.sub(p); let h = w.dot(*a); (w.sub(a.scale(h)).length(), h) };
                // The meridian parameter near which every parallel's lies (an ellipse's is periodic): a meridian edge's
                // midpoint.
                let reference = edges.iter().find_map(|&(e, fwd)| {
                    let seg = segment(body, e, fwd, shift);
                    if parallel_u(&seg, p, *a, 0.0, tolerance).is_some() { return None; }
                    let (t0, t1) = seg.domain();
                    Some(m.param_of(rz(seg.eval((t0 + t1) / 2.0).0), 0.0))
                }).unwrap_or(0.0);
                let u_at = |q: Vec3| m.param_of(rz(q), reference);
                let terms = |u: f64, ta: f64, tb: f64| -> Result<(f64, f64), KernelError> {
                    let (a_hat, b_hat) = m.flux_hats(u, pa);
                    let g_hat = m.area_hat(u, size).map_err(|e| KernelError::Unsupported(format!("face {f}: {e}")))?;
                    let rho_int = p.dot(*x) * (tb.sin() - ta.sin()) - p.dot(y) * (tb.cos() - ta.cos());
                    Ok((a_hat * rho_int + b_hat * (tb - ta), g_hat * (tb - ta)))
                };
                let (mut flux, mut sweep) = (0.0, 0.0);
                for &(e, fwd) in &edges {
                    let seg = segment(body, e, fwd, shift);
                    let (cc, cu, cv, t0, t1) = match seg {
                        Curve::Arc { c, n, x: ax, r: ar, t0, t1 } => (c, ax.scale(ar), n.cross(ax).scale(ar), t0, t1),
                        Curve::Conic { c, u, v, t0, t1 } => (c, u, v, t0, t1),
                        Curve::Bezier { .. } => continue, // a meridian: theta constant
                        _ => return Err(KernelError::InvalidInput(format!("face {f}, edge {e}: not a parallel or meridian of its surface of revolution"))),
                    };
                    if parallel_u(&seg, p, *a, 0.0, tolerance).is_none() { continue; } // a meridian ellipse
                    let perp = |q: Vec3| q.sub(a.scale(q.dot(*a)));
                    let (up, vp) = (perp(cu), perp(cv));
                    let dir = up.cross(vp).dot(*a).signum();
                    let theta_u = up.dot(y).atan2(up.dot(*x));
                    let (j, k) = terms(u_at(cc.add(up)), theta_u + dir * t0, theta_u + dir * t1)?;
                    flux += j;
                    sweep += k;
                }
                for (q, ta, tb) in pole_jumps(body, face, shift, p, *a, *x, y, sense, tolerance) {
                    let (j, k) = terms(u_at(q), ta, tb)?;
                    flux += j;
                    sweep += k;
                }
                volume -= flux / 3.0;
                area -= sense * sweep;
            }
            Surface::Sphere { c, a, x, r } | Surface::Torus { c, a, x, r, .. } => {
                let big = if let Surface::Torus { big, .. } = face.surface { big } else { 0.0 };
                let c = c.add(shift);
                let y = a.cross(*x);
                // Each parallel's meridian angle u; on a torus a face's parallels lie within one turn of u, so the branch
                // is cut at the widest gap between them (a sphere's u is always in [-pi/2, pi/2]).
                let edges: Vec<(usize, bool)> = face.loops.iter().flatten().copied().collect();
                let u_of = |q: Vec3| { let w = q.sub(c); let h = w.dot(*a); h.atan2(w.sub(a.scale(h)).length() - big) };
                let mut us: Vec<Option<f64>> = edges.iter().map(|&(e, fwd)| parallel_u(&segment(body, e, fwd, shift), c, *a, big, tolerance)).collect();
                // A face spans less than a turn of u: every parallel is taken within half a turn of the u at a meridian's
                // midpoint (the meridian runs between the face's ends), which fixes the branch on a torus.
                let reference = edges.iter().zip(&us).find(|(_, u)| u.is_none()).map(|(&(e, fwd), _)| {
                    let seg = segment(body, e, fwd, shift);
                    let (t0, t1) = seg.domain();
                    u_of(seg.eval((t0 + t1) / 2.0).0)
                });
                if let Some(reference) = reference {
                    for u in us.iter_mut().flatten() { *u = reference + (*u - reference + PI).rem_euclid(TAU) - PI; }
                }
                let (mut flux, mut sweep) = (0.0, 0.0);
                for (&(e, fwd), u) in edges.iter().zip(&us) {
                    let (j, k) = circle_family_terms(&segment(body, e, fwd, shift), *u, c, *a, *x, y, big, *r, tolerance)
                        .map_err(|m| KernelError::InvalidInput(format!("face {f}, edge {e}: {m}")))?;
                    flux += j;
                    sweep += k;
                }
                // Poles (a sphere's points on the axis) close the chart like an apex.
                for (q, ta, tb) in pole_jumps(body, face, shift, c, *a, *x, y, sense, tolerance) {
                    let (j, k) = circle_parallel(u_of(q), ta, tb, c, *a, *x, y, big, *r);
                    flux += j;
                    sweep += k;
                }
                volume -= flux / 3.0;
                area -= sense * sweep;
            }
            Surface::Conoid { c, a, x, rho, p, b, r, branch } => {
                let (c, p) = (c.add(shift), p.add(shift));
                let (mut flux, mut sweep) = (0.0, 0.0);
                for &(e, fwd) in face.loops.iter().flatten() {
                    let (j, k) = conoid_terms(&segment(body, e, fwd, shift), c, *a, *x, *rho, p, *b, *r, *branch, tolerance, size)
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
        if !volume.is_finite() || !area.is_finite() {
            return Err(KernelError::Unsupported(format!("face {f}: a boundary integral could not be evaluated")));
        }
        if let Some(output) = faces.as_mut() {
            let face_exact = !matches!(face.surface, Surface::EllipticCylinder { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. })
                && !face.loops.iter().flatten().any(|&(edge, _)| {
                    matches!(body.edges[edge].curve, Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } | Curve::RationalBezier { .. })
                        || matches!(face.surface, Surface::Cylinder { .. }) && matches!(body.edges[edge].curve, Curve::Bezier { .. })
                });
            output.push(FaceArea { area, exact: face_exact });
        }
        total_volume += volume;
        total_area += area;
    }
    if !total_volume.is_finite() || !total_area.is_finite() {
        return Err(KernelError::Unsupported("a boundary integral could not be evaluated".into()));
    }
    Ok(MassProperties { volume: total_volume, area: total_area, exact })
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
        Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } => return curve_on_cylinder_terms(curve, p, o, a, x, y, r, px, py, tol, size),
        Curve::Bezier { .. } | Curve::RationalBezier { .. } => {
            if crate::topology::off_surface(curve, &Surface::Cylinder { p, a, x, r }, tol).is_some() {
                return Err("Bezier boundary is not proved within cylinder tolerance".into());
            }
            // Integrate the exact cylindrical chart using the curve's angular
            // coordinate. The source Bezier stays unchanged, within the proved
            // existing geometric tolerance. This result uses quadrature.
            return elliptic_terms(curve, p, a, x.scale(r), y.scale(r), size);
        },
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

/// For one boundary curve on a radial conoid (module notes): (loop integral of G dt, loop integral of A dt), the measuring
/// origin at 0 (c and p already moved there). Along a curve t and s come from its point (angle and distance about the
/// axis), dt = (w.x w'.y - w.y w'.x) / s^2; a ruling (a line square to the axis at one angle) gives nothing.
#[allow(clippy::too_many_arguments)]
fn conoid_terms(curve: &Curve, c: Vec3, a: Vec3, x: Vec3, rho: f64, p: Vec3, b: Vec3, r: f64, branch: f64, tol: f64, size: f64)
    -> Result<(f64, f64), String> {
    let y = a.cross(x);
    if let Curve::Line { a: s, b: e } = *curve {
        let (ws, we) = (s.sub(c), e.sub(c));
        let (ts, te) = (ws.dot(y).atan2(ws.dot(x)), we.dot(y).atan2(we.dot(x)));
        let turn = (te - ts + PI).rem_euclid(TAU) - PI;
        if turn.abs() * ws.length().max(we.length()) > tol || e.sub(s).dot(a).abs() > tol {
            return Err("a line on a conoid must be one of its rulings".into());
        }
        return Ok((0.0, 0.0));
    }
    let (t0, t1) = curve.domain();
    let ca = c.dot(a);
    let parts = |tau: f64| {
        let (q, d) = curve.eval(tau);
        let w = q.sub(c);
        let (wx, wy) = (w.dot(x), w.dot(y));
        let s2 = wx * wx + wy * wy;
        let t = wy.atan2(wx);
        let dt = (wx * d.dot(y) - wy * d.dot(x)) / s2;
        let (h, dh) = crate::geom::pair_height(c, a, x, rho, p, b, r, branch, t);
        (s2.sqrt(), t, dt, h, dh, w.dot(a))
    };
    for i in 0..=8 {
        let (s, _, _, h, _, along) = parts(t0 + (t1 - t0) * i as f64 / 8.0);
        if !(s > tol) || (along - h).abs() > tol { return Err("a curve that does not lie on its conoid".into()); }
    }
    let j = quadrature(&|tau| {
        let (s, t, dt, h, dh, _) = parts(tau);
        let v = x.scale(-t.sin()).add(y.scale(t.cos()));
        (-(s * s / 2.0) * (ca + h) + s * dh * c.dot(v)) * dt
    }, t0, t1, 1e-14 * size.powi(3))?;
    let k = quadrature(&|tau| {
        let (s, _, dt, _, dh, _) = parts(tau);
        let g = dh.abs();
        let a_hat = s / 2.0 * (s * s + dh * dh).sqrt() + if g > 0.0 { dh * dh / 2.0 * (s / g).asinh() } else { 0.0 };
        a_hat * dt
    }, t0, t1, 1e-14 * size.powi(2))?;
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

/// For one boundary curve on a cone X = P + s a + R(s) rho(theta), R = r + k s (Claude 2026-10-29): the natural normal
/// X_theta x X_s = R (rho - k a), so X.N = R(s) g(theta) with g = P.rho - k (P.a) + r (the s terms cancel), and
/// |N| = R sqrt(1 + k^2). As for the cylinder, Green's theorem in (theta, s) turns the face integrals into loop integrals
/// of G(s) g dtheta and G(s) sqrt(1 + k^2) dtheta, G = r s + k s^2 / 2 (s from P). Generators have dtheta = 0; circles
/// across the axis have s constant, so both are closed form. Other curves on a cone are refused.
#[allow(clippy::too_many_arguments)]
fn cone_terms(curve: &Curve, p: Vec3, a: Vec3, x: Vec3, y: Vec3, r: f64, k: f64, tol: f64) -> Result<(f64, f64), String> {
    let perp = |w: Vec3| w.sub(a.scale(w.dot(a)));
    let (c, u, v, t0, t1) = match *curve {
        Curve::Line { a: s, b: e } => {
            // A generator lies in one half-plane through the axis: measured from its end farther from the axis (an apex
            // end has no direction of its own).
            let (ws, we) = (perp(s.sub(p)), perp(e.sub(p)));
            let w = if we.length() > ws.length() { we } else { ws };
            let d = perp(e.sub(s));
            if w.length() <= tol || d.cross(w.scale(1.0 / w.length())).length() > tol * e.sub(s).length().max(1.0) {
                return Err("a line on a cone must be a generator".into());
            }
            return Ok((0.0, 0.0));
        }
        Curve::Arc { c, n, x: ax, r: ar, t0, t1 } => (c, ax.scale(ar), n.cross(ax).scale(ar), t0, t1),
        Curve::Conic { c, u, v, t0, t1 } => (c, u, v, t0, t1),
        Curve::ConePair { .. } => return curve_on_cone_terms(curve, p, a, x, y, r, k, tol),
        Curve::CylinderPair { .. } | Curve::RadialPair { .. } => return Err("cylinder intersection curves on a cone are not supported".into()),
        Curve::Bezier { .. } | Curve::RationalBezier { .. } => return Err("a Bezier curve on a cone is not supported".into()),
    };
    let sc = c.sub(p).dot(a);
    let radius = r + k * sc;
    let (up, vp) = (perp(u), perp(v));
    if u.dot(a).abs() > tol || v.dot(a).abs() > tol || perp(c.sub(p)).length() > tol || radius <= tol
        || (up.length() - radius).abs() > tol || (vp.length() - radius).abs() > tol || up.dot(vp).abs() > tol * radius {
        return Err("a circle that does not lie across its cone".into());
    }
    let e = up.cross(vp).dot(a).signum();
    let theta_u = up.dot(y).atan2(up.dot(x));
    Ok(cone_parallel(sc, theta_u + e * t0, theta_u + e * t1, p, a, x, y, r, k))
}

/// The same two loop integrals for any curve on the cone (Claude 2026-10-31, a hole's ConePair rim), by quadrature in
/// the curve's own parameter: with w = point - P, s = w.a and theta its angle about the axis, dtheta/dt = (w.x w'.y -
/// w.y w'.x) / |w'|^2 across the axis, the integrands G(s) g(theta) and G(s) sqrt(1 + k^2) times dtheta/dt.
#[allow(clippy::too_many_arguments)]
fn curve_on_cone_terms(curve: &Curve, p: Vec3, a: Vec3, x: Vec3, y: Vec3, r: f64, k: f64, tol: f64) -> Result<(f64, f64), String> {
    let (t0, t1) = curve.domain();
    let (px, py, pa) = (p.dot(x), p.dot(y), p.dot(a));
    let parts = |t: f64| {
        let (q, d) = curve.eval(t);
        let w = q.sub(p);
        let (wx, wy, s) = (w.dot(x), w.dot(y), w.dot(a));
        let rr = wx * wx + wy * wy;
        let theta = wy.atan2(wx);
        (s, theta, (wx * d.dot(y) - wy * d.dot(x)) / rr, rr.sqrt())
    };
    for i in 0..=8 {
        let (s, _, _, rad) = parts(t0 + (t1 - t0) * i as f64 / 8.0);
        if (rad - (r + k * s)).abs() > tol * (1.0 + k * k).sqrt() { return Err("a curve that does not lie on its cone".into()); }
    }
    let size = curve.start().sub(p).length().max(1e-9);
    let big_g = |s: f64| r * s + k * s * s / 2.0;
    let j = quadrature(&|t| { let (s, th, dth, _) = parts(t); big_g(s) * (px * th.cos() + py * th.sin() + r - k * pa) * dth }, t0, t1,
        1e-14 * size.powi(3))?;
    let area = quadrature(&|t| { let (s, _, dth, _) = parts(t); big_g(s) * (1.0 + k * k).sqrt() * dth }, t0, t1, 1e-14 * size.powi(2))?;
    Ok((j, area))
}

/// The cone's loop integrals along a parallel at height s (from P) between angles ta and tb.
#[allow(clippy::too_many_arguments)]
fn cone_parallel(sc: f64, ta: f64, tb: f64, p: Vec3, a: Vec3, x: Vec3, y: Vec3, r: f64, k: f64) -> (f64, f64) {
    let (px, py, pa) = (p.dot(x), p.dot(y), p.dot(a));
    let g = px * (tb.sin() - ta.sin()) - py * (tb.cos() - ta.cos()) + (r - k * pa) * (tb - ta);
    let big_g = r * sc + k * sc * sc / 2.0;
    (big_g * g, big_g * (1.0 + k * k).sqrt() * (tb - ta))
}

/// Where a face's loop passes through the axis (a cone's apex, a sphere's pole), the (theta, u) chart's boundary runs
/// along the pole from the arriving meridian's angle to the departing one's, which no edge carries. The jump's size is
/// the angle between the two meridians (revolve splits faces so it is at most half a turn); its direction follows the
/// loop's counter-clockwise turn about the outward normal: backwards at the pole of greatest u, forwards at the least,
/// reversed when the face's sense is. Returns (pole point, angle from, angle to) in the measuring frame.
#[allow(clippy::too_many_arguments)]
fn pole_jumps(body: &Body, face: &crate::neutral::Face, shift: Vec3, o: Vec3, a: Vec3, x: Vec3, y: Vec3, sense: f64, tol: f64)
    -> Vec<(Vec3, f64, f64)> {
    let mut out = Vec::new();
    let theta = |q: Vec3| { let d = q.sub(o); d.dot(y).atan2(d.dot(x)) };
    let radial = |q: Vec3| { let d = q.sub(o); d.sub(a.scale(d.dot(a))).length() };
    let mid = |e: usize, fwd: bool| { let c = segment(body, e, fwd, shift); let (t0, t1) = c.domain(); c.eval((t0 + t1) / 2.0).0 };
    for lp in &face.loops {
        for k in 0..lp.len() {
            let (e_in, f_in) = lp[k];
            let (e_out, f_out) = lp[(k + 1) % lp.len()];
            let v = if f_in { body.edges[e_in].end } else { body.edges[e_in].start };
            let q = body.vertices[v].add(shift);
            if radial(q) > tol { continue; }
            let (m_in, m_out) = (mid(e_in, f_in), mid(e_out, f_out));
            let (ta, tb) = (theta(m_in), theta(m_out));
            let turn = (tb - ta + PI).rem_euclid(TAU) - PI;
            let size = if turn <= -PI { PI } else { turn.abs() };
            // The pole is the face's extreme in u: above the arriving meridian's midpoint means the top.
            let top = q.sub(m_in).dot(a) > 0.0;
            let dir = if top { -1.0 } else { 1.0 } * sense;
            out.push((q, ta, ta + dir * size));
        }
    }
    out
}

/// A parallel's meridian angle on a sphere (big = 0) or torus: atan2(h, rho - big) for a circle across the axis with
/// centre on it at height h and radius rho; None for any other curve.
fn parallel_u(curve: &Curve, c: Vec3, a: Vec3, big: f64, tol: f64) -> Option<f64> {
    let (cc, u, v) = match *curve {
        Curve::Arc { c, n, x, r, .. } => (c, x.scale(r), n.cross(x).scale(r)),
        Curve::Conic { c, u, v, .. } => (c, u, v),
        _ => return None,
    };
    let w = cc.sub(c);
    let perp = |q: Vec3| q.sub(a.scale(q.dot(a)));
    if u.dot(a).abs() > tol || v.dot(a).abs() > tol || perp(w).length() > tol { return None; }
    Some(w.dot(a).atan2(perp(u).length() - big))
}

/// For one boundary curve on a sphere (big = 0) or torus X = C + R(u) rho(theta) + r sin(u) a, R = big + r cos u (Claude
/// 2026-10-29): X_theta x X_u = R r (cos u rho + sin u a), so X.N = A(u) (C.rho) + B(u) with A = R r cos u and
/// B = R r (big cos u + r + (C.a) sin u), and |N| = R r. Green's theorem in (theta, u) turns the face integrals into loop
/// integrals of (A^(u) (C.rho) + B^(u)) dtheta and r (big u + r sin u) dtheta (antiderivatives in u, closed form). Along a
/// parallel u is constant; along a meridian theta is; anything else is refused.
#[allow(clippy::too_many_arguments)]
fn circle_family_terms(curve: &Curve, u: Option<f64>, c: Vec3, a: Vec3, x: Vec3, y: Vec3, big: f64, r: f64, tol: f64)
    -> Result<(f64, f64), String> {
    let (cc, cu, cv, t0, t1) = match *curve {
        Curve::Arc { c, n, x: ax, r: ar, t0, t1 } => (c, ax.scale(ar), n.cross(ax).scale(ar), t0, t1),
        Curve::Conic { c, u, v, t0, t1 } => (c, u, v, t0, t1),
        _ => return Err("only circles lie on a sphere or torus face's boundary here".into()),
    };
    let Some(u) = u else {
        // A meridian: its plane contains the axis (theta constant along it), so it contributes nothing.
        let n = cu.cross(cv);
        let ln = n.length();
        if ln == 0.0 || n.scale(1.0 / ln).dot(a).abs() > tol || n.scale(1.0 / ln).dot(cc.sub(c)).abs() > tol {
            return Err("a circle on a sphere or torus that is neither a parallel nor a meridian".into());
        }
        return Ok((0.0, 0.0));
    };
    let perp = |q: Vec3| q.sub(a.scale(q.dot(a)));
    let (up, vp) = (perp(cu), perp(cv));
    let e = up.cross(vp).dot(a).signum();
    let theta_u = up.dot(y).atan2(up.dot(x));
    Ok(circle_parallel(u, theta_u + e * t0, theta_u + e * t1, c, a, x, y, big, r))
}

/// The sphere's or torus's loop integrals along a parallel at meridian angle u between angles ta and tb.
#[allow(clippy::too_many_arguments)]
fn circle_parallel(u: f64, ta: f64, tb: f64, c: Vec3, a: Vec3, x: Vec3, y: Vec3, big: f64, r: f64) -> (f64, f64) {
    let (cx, cy, ca) = (c.dot(x), c.dot(y), c.dot(a));
    let (s, co) = (u.sin(), u.cos());
    let a_hat = r * big * s + r * r * (u / 2.0 + (2.0 * u).sin() / 4.0);
    let b_hat = r * big * big * s + r * r * big * u - r * ca * big * co + r * r * big * (u / 2.0 + (2.0 * u).sin() / 4.0)
        + r * r * r * s + r * r * ca * s * s / 2.0;
    let rho_int = cx * (tb.sin() - ta.sin()) - cy * (tb.cos() - ta.cos());
    (a_hat * rho_int + b_hat * (tb - ta), r * (big * u + r * s) * (tb - ta))
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

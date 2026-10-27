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

use crate::extrude::signed_area;
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::Body;
use crate::KernelError;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassProperties {
    pub volume: f64,
    pub area: f64,
}

fn segment(body: &Body, edge: usize, forward: bool) -> Curve {
    let c = &body.edges[edge].curve;
    if forward { c.clone() } else { c.reversed() }
}

/// Volume and area; a body whose faces point inwards has negative volume. The flux is taken about the vertices'
/// centroid, not the world origin, so a small part far from the origin loses nothing to cancellation (Astra).
/// Curves that do not lie on their cylinder (beyond `tolerance`) are refused.
pub fn measure(body: &Body, tolerance: f64) -> Result<MassProperties, KernelError> {
    let (mut volume, mut area) = (0.0, 0.0);
    let n = body.vertices.len().max(1) as f64;
    let o = body.vertices.iter().fold(Vec3::new(0.0, 0.0, 0.0), |s, v| s.add(*v)).scale(1.0 / n);
    for (f, face) in body.faces.iter().enumerate() {
        let sense = if face.sense { 1.0 } else { -1.0 };
        match &face.surface {
            Surface::Plane { p, n, .. } => {
                let outward = n.scale(sense);
                let a: f64 = face.loops.iter()
                    .map(|lp| signed_area(&lp.iter().map(|&(e, fwd)| segment(body, e, fwd)).collect::<Vec<_>>(), outward))
                    .sum();
                area += a;
                volume += p.sub(o).dot(outward) * a / 3.0;
            }
            Surface::Cylinder { p, a, x, r } => {
                let y = a.cross(*x);
                let (px, py) = (p.sub(o).dot(*x), p.sub(o).dot(y));
                let (mut flux, mut sweep) = (0.0, 0.0);
                for &(e, fwd) in face.loops.iter().flatten() {
                    let (j, k) = cylinder_terms(&segment(body, e, fwd), *p, *a, *x, y, *r, px, py, tolerance)
                        .map_err(|m| KernelError::InvalidInput(format!("face {f}, edge {e}: {m}")))?;
                    flux += j;
                    sweep += k;
                }
                volume -= flux / 3.0;
                area -= sense * sweep;
            }
        }
    }
    Ok(MassProperties { volume, area })
}

/// For one boundary curve on the cylinder: (integral of s f dtheta, integral of s r dtheta).
#[allow(clippy::too_many_arguments)]
fn cylinder_terms(curve: &Curve, p: Vec3, a: Vec3, x: Vec3, y: Vec3, r: f64, px: f64, py: f64, tol: f64)
    -> Result<(f64, f64), String> {
    let (c, u, v, t0, t1) = match *curve {
        Curve::Line { a: s, b: e } => {
            let d = e.sub(s);
            if d.cross(a).length() > tol * d.length().max(1.0) { return Err("a line on a cylinder must run along its axis".into()); }
            return Ok((0.0, 0.0));
        }
        Curve::Arc { c, n, x: ax, r: ar, t0, t1 } => (c, ax.scale(ar), n.cross(ax).scale(ar), t0, t1),
        Curve::Conic { c, u, v, t0, t1 } => (c, u, v, t0, t1),
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
    let (sc, cu, cv) = (dc.dot(a), u.dot(a), v.dot(a));
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

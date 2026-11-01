//! One blended edge's exact geometry (Claude 2026-11-01): its two rails (where the blend meets each face), the blend
//! surface and its sense.
//!
//! Straight edge between two planes. With the edge's direction d as it runs in face A, outward normals nA, nB and the
//! in-face directions wA = nA x d, wB = nB x (-d) (each points from the edge into its face), the faces' opening angle
//! phi has cos phi = wA.wB. A chamfer meets A at distance sA and B at sB (the plane through both rails); a fillet of
//! radius r touches both faces at sA = sB = r / tan(phi / 2), its axis r inside A from the touching line (inside the
//! material for a convex edge, outside it for a concave one) and parallel to d: a circular cylinder. A distance-angle
//! chamfer (distance s on face F at angle theta to F) meets the other face at s sin theta / sin(theta + phi) (the sine
//! rule in the triangle the chamfer cuts from the cross-section).
//!
//! Circle where a plane meets a cylinder of radius R at right angles. In the half-plane through the axis the faces are
//! the lines h = 0 and rho = R; a chamfer is the line from (R + alpha sP, 0) to (R, beta sC) (alpha, beta = +-1: which
//! way each face runs from the edge), revolved: a cone. A fillet of radius q is the quarter circle about
//! (R + alpha q, beta q) revolved: a torus with centre-circle radius R + alpha q and tube radius q.

use super::{Item, Kind};
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::Body;
use crate::KernelError;
use std::f64::consts::PI;

pub(crate) fn unsupported<T>(m: impl Into<String>) -> Result<T, KernelError> { Err(KernelError::Unsupported(m.into())) }

#[derive(Clone, Debug)]
pub(crate) enum Shape {
    /// Rails o + w s + t d on each face; a fillet's axis point (at o's section) when it is one.
    Planar { o: Vec3, d: Vec3, wa: Vec3, wb: Vec3, sa: f64, sb: f64, axis: Option<Vec3> },
    /// The edge circle's centre c0 (on the plane), the plane's outward normal ah, the arc's own frame (n, x, t0, t1),
    /// the cylinder radius, the ways each face runs from the edge (alpha radially on the plane, beta along ah on the
    /// cylinder) and the setbacks on the plane (sp) and on the cylinder (sc).
    Rim { c0: Vec3, ah: Vec3, n: Vec3, x: Vec3, t0: f64, t1: f64, radius: f64, alpha: f64, beta: f64, sp: f64, sc: f64, plane: usize, cylinder: usize },
}

#[derive(Clone, Debug)]
pub(crate) struct EdgeBlend {
    pub edge: usize,
    /// The face in which the edge runs forward, and the other.
    pub a: usize,
    pub b: usize,
    pub radius: Option<f64>,
    pub convex: bool,
    pub shape: Shape,
    pub surface: Surface,
    pub sense: bool,
}

/// A plane face's point and outward unit normal.
pub(crate) fn plane_of(body: &Body, f: usize) -> Option<(Vec3, Vec3)> {
    match body.faces[f].surface {
        Surface::Plane { p, n, .. } => Some((p, if body.faces[f].sense { n } else { n.scale(-1.0) })),
        _ => None,
    }
}

/// The setbacks on faces a and b for an opening angle phi between the faces' in-face directions.
fn setbacks(item: &Item, a: usize, b: usize, phi: f64) -> Result<(f64, f64), KernelError> {
    let e = item.edge;
    let on = |face: usize, here: f64, there: f64| -> Result<(f64, f64), KernelError> {
        if face == a { Ok((here, there)) } else if face == b { Ok((there, here)) } else {
            Err(KernelError::InvalidInput(format!("edge {e}: face {face} is not one of its two faces ({a}, {b})")))
        }
    };
    match item.kind {
        Kind::Fillet { radius } => { let s = radius / (phi / 2.0).tan(); Ok((s, s)) }
        Kind::Chamfer { distance } => Ok((distance, distance)),
        Kind::ChamferDistances { first, second, face } => on(face, first, second),
        Kind::ChamferAngle { distance, angle, face } => {
            if angle + phi >= PI - 1e-9 {
                return unsupported(format!("edge {e}: a chamfer at {angle} radians to its face never meets the other face (the faces open at {phi} radians)"));
            }
            on(face, distance, distance * angle.sin() / (angle + phi).sin())
        }
    }
}

pub(crate) fn plan_edge(body: &Body, item: &Item, a: usize, b: usize, tol: f64) -> Result<EdgeBlend, KernelError> {
    let e = item.edge;
    if a == b { return unsupported(format!("edge {e} is a seam within face {a}; only an edge between two faces can be blended")); }
    let radius = match item.kind { Kind::Fillet { radius } => Some(radius), _ => None };
    let edge = &body.edges[e];
    match (&edge.curve, plane_of(body, a), plane_of(body, b)) {
        (Curve::Line { a: p0, b: p1 }, Some((_, na)), Some((_, nb))) => {
            let d = p1.sub(*p0).unit()?;
            let (wa, wb) = (na.cross(d), nb.cross(d.scale(-1.0)));
            let cross = na.cross(nb);
            if cross.length() <= 1e-9 {
                return unsupported(format!("edge {e} lies between faces {a} and {b}, which meet flat (no crease to blend)"));
            }
            let convex = cross.dot(d) > 0.0;
            let phi = wa.dot(wb).clamp(-1.0, 1.0).acos();
            let (sa, sb) = setbacks(item, a, b, phi)?;
            let (pa, pb) = (p0.add(wa.scale(sa)), p0.add(wb.scale(sb)));
            let sigma = if convex { 1.0 } else { -1.0 };
            let (surface, sense, axis) = match radius {
                Some(r) => {
                    let c = pa.sub(na.scale(sigma * r));
                    if c.distance(pb.sub(nb.scale(sigma * r))) > tol {
                        return Err(KernelError::Internal(format!("edge {e}: the fillet axis is not at its radius from both faces")));
                    }
                    (Surface::Cylinder { p: c, a: d, x: na.scale(sigma), r }, convex, Some(c))
                }
                None => {
                    let mut n = d.cross(pb.sub(pa)).unit()?;
                    if n.dot(na.add(nb)) < 0.0 { n = n.scale(-1.0); }
                    (Surface::Plane { p: pa, n, x: d }, true, None)
                }
            };
            Ok(EdgeBlend { edge: e, a, b, radius, convex, shape: Shape::Planar { o: *p0, d, wa, wb, sa, sb, axis }, surface, sense })
        }
        (Curve::Arc { c, n, x, r, t0, t1 }, ..) => {
            let (plane, cylinder) = match (plane_of(body, a), plane_of(body, b)) {
                (Some(_), None) => (a, b),
                (None, Some(_)) => (b, a),
                _ => return unsupported(format!("edge {e}: a circular edge can be blended only where a flat face meets a cylinder")),
            };
            let (pp, ah) = plane_of(body, plane).expect("a plane face");
            let Surface::Cylinder { p: cp, a: ca, r: cr, .. } = body.faces[cylinder].surface else {
                return unsupported(format!("edge {e}: a circular edge can be blended only where a flat face meets a cylinder (face {cylinder} is not one)"));
            };
            let perp = |w: Vec3| w.sub(ca.scale(w.dot(ca)));
            if ah.cross(ca).length() > 1e-9 || n.cross(ca).length() > 1e-9 {
                return unsupported(format!("edge {e}: the flat face {plane} is not at right angles to the cylinder {cylinder} (the blend would not be a torus or cone)"));
            }
            if perp(c.sub(cp)).length() > tol || (r - cr).abs() > tol || c.sub(pp).dot(ah).abs() > tol {
                return Err(KernelError::InvalidInput(format!("edge {e}: its circle is not where face {plane} meets the cylinder {cylinder}")));
            }
            let (x0, tangent) = edge.curve.eval(*t0);
            let te = tangent.unit()?;
            let tp = if plane == a { te } else { te.scale(-1.0) };
            let rho = perp(x0.sub(*c)).unit()?;
            let nc = if body.faces[cylinder].sense { rho } else { rho.scale(-1.0) };
            let (wp, wc) = (ah.cross(tp), nc.cross(tp.scale(-1.0)));
            let (alpha, beta) = (wp.dot(rho).signum(), wc.dot(ah).signum());
            let (na, nb) = if plane == a { (ah, nc) } else { (nc, ah) };
            let convex = na.cross(nb).dot(te) > 0.0;
            let (sa, sb) = setbacks(item, a, b, PI / 2.0)?;
            let (sp, sc) = if plane == a { (sa, sb) } else { (sb, sa) };
            let rp = r + alpha * sp;
            if rp <= tol {
                return unsupported(format!("edge {e}: the blend ({sp} across face {plane}) reaches the cylinder's axis and consumes the flat face"));
            }
            let (surface, sense) = match radius {
                Some(q) => (Surface::Torus { c: c.add(ah.scale(beta * q)), a: ah, x: *x, big: rp, r: q }, convex),
                None => {
                    let k = -alpha * sp / (beta * sc);
                    let normal = rho.sub(ah.scale(k));
                    (Surface::Cone { p: *c, a: ah, x: *x, r: rp, k }, normal.dot(na.add(nb)) > 0.0)
                }
            };
            Ok(EdgeBlend { edge: e, a, b, radius, convex, surface, sense,
                shape: Shape::Rim { c0: *c, ah, n: *n, x: *x, t0: *t0, t1: *t1, radius: *r, alpha, beta, sp, sc, plane, cylinder } })
        }
        _ => unsupported(format!("edge {e}: blends are supported on straight edges between two flat faces and on circles where a flat face meets a cylinder at right angles")),
    }
}

impl EdgeBlend {
    /// The rail on a planar blend's face as (point, direction).
    pub(crate) fn rail_line(&self, face: usize) -> (Vec3, Vec3) {
        match self.shape {
            Shape::Planar { o, d, wa, wb, sa, sb, .. } => {
                if face == self.a { (o.add(wa.scale(sa)), d) } else { (o.add(wb.scale(sb)), d) }
            }
            Shape::Rim { .. } => unreachable!("a rim blend has circular rails"),
        }
    }
    /// The rail on `face` from `from` to `to`, running the edge's way.
    pub(crate) fn rail(&self, face: usize, from: Vec3, to: Vec3) -> Curve {
        match self.shape {
            Shape::Planar { .. } => Curve::Line { a: from, b: to },
            Shape::Rim { c0, ah, n, x, t0, t1, radius, alpha, beta, sp, sc, plane, .. } => {
                if face == plane { Curve::Arc { c: c0, n, x, r: radius + alpha * sp, t0, t1 } }
                else { Curve::Arc { c: c0.add(ah.scale(beta * sc)), n, x, r: radius, t0, t1 } }
            }
        }
    }
}

/// The short arc about c from `from` to `to` (less than a half turn).
pub(crate) fn arc_between(c: Vec3, from: Vec3, to: Vec3) -> Result<Curve, KernelError> {
    let (f, t) = (from.sub(c), to.sub(c));
    let r = f.length();
    let n = f.cross(t);
    if n.length() <= 1e-12 * r * r {
        return unsupported("a blend boundary arc would be a half turn or none (ambiguous)");
    }
    let angle = n.length().atan2(f.dot(t));
    Ok(Curve::Arc { c, n: n.unit()?, x: f.unit()?, r, t0: 0.0, t1: angle })
}

/// Where the cylinder (axis point c, unit axis d, radius r) crosses the plane (point p, unit normal m), from `from` to
/// `to` (both on both, less than a half turn apart round the axis): a circle when the plane is square to the axis,
/// else the ellipse centre + u cos t + v sin t with u = r e1, v = r e2 + lambda d (e1 = unit(d x m), e2 = d x e1,
/// lambda = -r (e2.m) / (d.m)): its part across the axis is the circle r (e1 cos t + e2 sin t), and it lies in the plane.
pub(crate) fn cylinder_section(c: Vec3, d: Vec3, r: f64, p: Vec3, m: Vec3, from: Vec3, to: Vec3) -> Result<Curve, KernelError> {
    let dm = d.dot(m);
    if dm.abs() <= 1e-9 { return unsupported("a blend's cylinder runs along the face it should end on"); }
    let centre = c.add(d.scale(p.sub(c).dot(m) / dm));
    if d.cross(m).length() <= 1e-12 { return arc_between(centre, from, to); }
    let e1 = d.cross(m).unit()?;
    let e2 = d.cross(e1);
    let lambda = -r * e2.dot(m) / dm;
    let (u, v) = (e1.scale(r), e2.scale(r).add(d.scale(lambda)));
    let angle = |q: Vec3| { let w = q.sub(centre); w.dot(e2).atan2(w.dot(e1)) };
    let tf = angle(from);
    let delta = (angle(to) - tf + PI).rem_euclid(2.0 * PI) - PI;
    if delta.abs() <= 1e-12 || delta.abs() >= PI - 1e-12 {
        return unsupported("a blend's end section would be a half turn or none (ambiguous)");
    }
    Ok(if delta > 0.0 { Curve::Conic { c: centre, u, v, t0: tf, t1: tf + delta } }
       else { Curve::Conic { c: centre, u, v, t0: tf + delta, t1: tf }.reversed() })
}

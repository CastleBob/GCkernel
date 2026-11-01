//! What happens at a vertex where blends end (Claude 2026-11-01). Only vertices where three faces meet are handled;
//! k of the three edges there are blended:
//! - k = 1: the blend runs out on the third (end) face. Each rail stops where it meets the end face's edge on its face;
//!   those edges are trimmed there and the end face gains the blend's section by its plane (a line for a chamfer, a
//!   circle or ellipse arc for a fillet).
//! - k = 2: the two blends meet each other. On the face they share their rails meet at T; on the other two faces they
//!   must reach the unblended edge at one point Q (refused otherwise), which trims it. The joint T-Q is the line where
//!   two chamfer planes meet, the ellipse where two equal fillet cylinders meet (they are mirror images in the plane
//!   through their axes' crossing, T and Q), or, along one plane-cylinder rim, the torus meridian / cone generator.
//! - k = 3: a corner patch: the plane through the three face corners (chamfers), or for equal fillets the sphere about
//!   the point at the radius from all three faces, bounded by its three great arcs (admitted when one face is at right
//!   angles to the other two, so the arcs are its frame's equator and meridians).

use super::edge::{arc_between, cylinder_section, plane_of, unsupported, EdgeBlend, Shape};
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge};
use crate::KernelError;
use std::collections::HashMap;

/// One end of an edge: the edge and whether it is the edge's end (else its start).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) struct End { pub edge: usize, pub at_end: bool }

/// A face's corner at a vertex: loop `lp` arrives by `inn` and leaves by `out`, segment `k` being the one leaving.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Corner { pub face: usize, pub lp: usize, pub k: usize, pub inn: End, pub out: End }

#[derive(Default)]
pub(crate) struct Built { pub vertices: Vec<Vec3>, pub edges: Vec<Edge> }

impl Built {
    pub(crate) fn vertex(&mut self, p: Vec3) -> usize { self.vertices.push(p); self.vertices.len() - 1 }
    pub(crate) fn edge(&mut self, curve: Curve, start: usize, end: usize) -> usize {
        self.edges.push(Edge { curve, start, end });
        self.edges.len() - 1
    }
}

/// The new topology at one vertex.
#[derive(Default)]
pub(crate) struct Plan {
    /// Unblended edge ends: their new vertex.
    pub point: HashMap<End, usize>,
    /// Blended edge ends on each of their faces: the rail's vertex there.
    pub rail: HashMap<(End, usize), usize>,
    /// Blended edge ends: the edge closing the blend face there.
    pub chain: HashMap<End, usize>,
    /// Edges inserted at a face corner (face, loop, segment), with their direction there.
    pub insert: HashMap<(usize, usize, usize), (usize, bool)>,
    /// A corner patch: its surface, sense and edges.
    pub corner: Option<(Surface, bool, Vec<usize>)>,
}

fn meet(l1: (Vec3, Vec3), l2: (Vec3, Vec3), tol: f64, what: &str) -> Result<Vec3, KernelError> {
    let ((p1, d1), (p2, d2)) = (l1, l2);
    let n = d1.cross(d2);
    let nn = n.dot(n);
    if nn <= 1e-18 { return unsupported(format!("{what}: the lines run parallel")); }
    let w = p2.sub(p1);
    let s = w.cross(d2).dot(n) / nn;
    let t = w.cross(d1).dot(n) / nn;
    let (a, b) = (p1.add(d1.scale(s)), p2.add(d2.scale(t)));
    if a.distance(b) > tol { return unsupported(format!("{what}: the lines do not meet")); }
    Ok(a.add(b).scale(0.5))
}

fn line_of(body: &Body, e: usize) -> Result<(Vec3, Vec3), KernelError> {
    match body.edges[e].curve {
        Curve::Line { a, b } => Ok((a, b.sub(a).unit()?)),
        _ => unsupported(format!("edge {e} is not straight")),
    }
}

pub(crate) fn plan(v: usize, corners: &[Corner], blends: &HashMap<usize, EdgeBlend>, body: &Body, built: &mut Built, tol: f64)
    -> Result<Plan, KernelError> {
    if corners.len() != 3 {
        return unsupported(format!("vertex {v} joins {} face corners; a blend may end only at a vertex where three faces meet", corners.len()));
    }
    let ends: Vec<End> = corners.iter().map(|c| c.inn).collect();
    let selected: Vec<End> = ends.iter().copied().filter(|x| blends.contains_key(&x.edge)).collect();
    let rims = selected.iter().filter(|s| matches!(blends[&s.edge].shape, Shape::Rim { .. })).count();
    if rims > 0 {
        if rims != selected.len() {
            return unsupported(format!("vertex {v}: a plane-cylinder blend meets a blend of another kind"));
        }
        return rim(v, corners, &ends, &selected, blends, body, built, tol);
    }
    planar(v, corners, &ends, &selected, blends, body, built, tol)
}

fn find<'a>(corners: &'a [Corner], f: impl Fn(&Corner) -> bool) -> &'a Corner {
    corners.iter().find(|c| f(c)).expect("the corner exists at a three-face vertex")
}

#[allow(clippy::too_many_arguments)]
fn planar(v: usize, corners: &[Corner], ends: &[End], selected: &[End], blends: &HashMap<usize, EdgeBlend>, body: &Body,
    built: &mut Built, tol: f64) -> Result<Plan, KernelError> {
    let faces: Vec<usize> = corners.iter().map(|c| c.face).collect();
    if faces[0] == faces[1] || faces[1] == faces[2] || faces[0] == faces[2] {
        return unsupported(format!("vertex {v}: a face meets itself there"));
    }
    if faces.iter().any(|&f| plane_of(body, f).is_none()) || ends.iter().any(|x| !matches!(body.edges[x.edge].curve, Curve::Line { .. })) {
        return unsupported(format!("vertex {v}: a blend between flat faces may end only where three flat faces meet along straight edges"));
    }
    let mut plan = Plan::default();
    match selected.len() {
        1 => {
            let s = selected[0];
            let eb = &blends[&s.edge];
            for c in corners.iter().filter(|c| c.inn == s || c.out == s) {
                let other = if c.inn == s { c.out } else { c.inn };
                let p = meet(eb.rail_line(c.face), line_of(body, other.edge)?, tol,
                    &format!("vertex {v}: the blend of edge {} meeting edge {}", s.edge, other.edge))?;
                let id = built.vertex(p);
                plan.point.insert(other, id);
                plan.rail.insert((s, c.face), id);
            }
            let h = find(corners, |c| c.inn != s && c.out != s);
            let (from, to) = (plan.point[&h.inn], plan.point[&h.out]);
            let (hp, hn) = plane_of(body, h.face).expect("a plane");
            let (pf, pt) = (built.vertices[from], built.vertices[to]);
            let curve = match (eb.radius, &eb.shape) {
                (Some(r), Shape::Planar { axis: Some(c), d, .. }) => cylinder_section(*c, *d, r, hp, hn, pf, pt)?,
                _ => Curve::Line { a: pf, b: pt },
            };
            let id = built.edge(curve, from, to);
            plan.chain.insert(s, id);
            plan.insert.insert((h.face, h.lp, h.k), (id, true));
        }
        2 => {
            let common = *find(corners, |c| selected.contains(&c.inn) && selected.contains(&c.out));
            let (s1, s2) = (common.inn, common.out);
            let u = *ends.iter().find(|x| !selected.contains(x)).expect("one unblended end");
            let (b1, b2) = (&blends[&s1.edge], &blends[&s2.edge]);
            if b1.radius.is_some() != b2.radius.is_some() {
                return unsupported(format!("vertex {v}: a fillet (edge {}) and a chamfer (edge {}) meet there", if b1.radius.is_some() { s1.edge } else { s2.edge }, if b1.radius.is_some() { s2.edge } else { s1.edge }));
            }
            if b1.convex != b2.convex {
                return unsupported(format!("vertex {v}: a convex and a concave blend (edges {} and {}) meet there", s1.edge, s2.edge));
            }
            if let (Some(r1), Some(r2)) = (b1.radius, b2.radius) {
                if (r1 - r2).abs() > tol {
                    return unsupported(format!("vertex {v}: fillets of different radii ({r1}, {r2}) meet there; their cylinders do not meet in a plane curve"));
                }
            }
            let t = meet(b1.rail_line(common.face), b2.rail_line(common.face), tol, &format!("vertex {v}: the blends of edges {} and {}", s1.edge, s2.edge))?;
            let ca = *find(corners, |c| (c.inn == s1 || c.out == s1) && (c.inn == u || c.out == u));
            let cb = *find(corners, |c| (c.inn == s2 || c.out == s2) && (c.inn == u || c.out == u));
            let what = format!("vertex {v}: the blends meeting edge {}", u.edge);
            let qa = meet(b1.rail_line(ca.face), line_of(body, u.edge)?, tol, &what)?;
            let qb = meet(b2.rail_line(cb.face), line_of(body, u.edge)?, tol, &what)?;
            if qa.distance(qb) > tol {
                return unsupported(format!("vertex {v}: the blends of edges {} and {} reach edge {} at different points ({} apart); their faces meet the shared face at different angles or the blends differ in size there", s1.edge, s2.edge, u.edge, qa.distance(qb)));
            }
            let curve = match (&b1.shape, &b2.shape, b1.radius) {
                (Shape::Planar { axis: Some(c1), d: d1, .. }, Shape::Planar { axis: Some(c2), d: d2, .. }, Some(r)) => mitre(*c1, *d1, *c2, *d2, r, t, qa, tol, v)?,
                _ => Curve::Line { a: t, b: qa },
            };
            let (ti, qi) = (built.vertex(t), built.vertex(qa));
            plan.rail.insert((s1, common.face), ti);
            plan.rail.insert((s2, common.face), ti);
            plan.rail.insert((s1, ca.face), qi);
            plan.rail.insert((s2, cb.face), qi);
            plan.point.insert(u, qi);
            let id = built.edge(curve, ti, qi);
            plan.chain.insert(s1, id);
            plan.chain.insert(s2, id);
        }
        3 => {
            let all: Vec<&EdgeBlend> = selected.iter().map(|s| &blends[&s.edge]).collect();
            if all.iter().any(|b| b.radius.is_some() != all[0].radius.is_some()) {
                return unsupported(format!("vertex {v}: fillets and chamfers meet at one corner"));
            }
            if all.iter().any(|b| b.convex != all[0].convex) {
                return unsupported(format!("vertex {v}: convex and concave blends meet at one corner"));
            }
            let mut point_of = HashMap::new();
            for c in corners {
                let (bi, bo) = (&blends[&c.inn.edge], &blends[&c.out.edge]);
                let t = meet(bi.rail_line(c.face), bo.rail_line(c.face), tol, &format!("vertex {v}: the blends on face {}", c.face))?;
                let id = built.vertex(t);
                point_of.insert(c.face, id);
                plan.rail.insert((c.inn, c.face), id);
                plan.rail.insert((c.out, c.face), id);
            }
            let normals: Vec<Vec3> = corners.iter().map(|c| plane_of(body, c.face).expect("a plane").1).collect();
            let sum = normals.iter().fold(Vec3::new(0.0, 0.0, 0.0), |s, n| s.add(*n));
            let tp: Vec<Vec3> = corners.iter().map(|c| built.vertices[point_of[&c.face]]).collect();
            let (surface, sense, centre) = match all[0].radius {
                None => {
                    let mut n = tp[1].sub(tp[0]).cross(tp[2].sub(tp[0])).unit()
                        .map_err(|_| KernelError::Unsupported(format!("vertex {v}: the chamfers' corners lie in a line")))?;
                    if n.dot(sum) < 0.0 { n = n.scale(-1.0); }
                    let w = tp[1].sub(tp[0]);
                    (Surface::Plane { p: tp[0], n, x: w.sub(n.scale(w.dot(n))).unit()? }, true, None)
                }
                Some(r) => {
                    if all.iter().any(|b| (b.radius.unwrap() - r).abs() > tol) {
                        return unsupported(format!("vertex {v}: fillets of different radii meet at one corner"));
                    }
                    let sigma = if all[0].convex { 1.0 } else { -1.0 };
                    let pole = (0..3).find(|&i| (0..3).all(|j| j == i || normals[i].dot(normals[j]).abs() <= 1e-12))
                        .ok_or_else(|| KernelError::Unsupported(format!("vertex {v}: a three-fillet corner is supported only where one face is at right angles to the other two (its sphere patch is then bounded by its equator and meridians)")))?;
                    let other = (pole + 1) % 3;
                    let (n1, n2, n3) = (normals[0], normals[1], normals[2]);
                    let level = |i: usize| plane_of(body, corners[i].face).unwrap().0.dot(normals[i]) - sigma * r;
                    let det = n1.dot(n2.cross(n3));
                    let c = n2.cross(n3).scale(level(0)).add(n3.cross(n1).scale(level(1))).add(n1.cross(n2).scale(level(2))).scale(1.0 / det);
                    for i in 0..3 {
                        if tp[i].distance(c.add(normals[i].scale(sigma * r))) > tol {
                            return Err(KernelError::Internal(format!("vertex {v}: the fillet corner point on face {} is not where the sphere touches it", corners[i].face)));
                        }
                    }
                    (Surface::Sphere { c, a: normals[pole].scale(sigma), x: normals[other].scale(sigma), r }, all[0].convex, Some(c))
                }
            };
            let mut edges = Vec::new();
            for &s in selected {
                let ca = find(corners, |c| c.out == s);
                let cb = find(corners, |c| c.inn == s);
                let (from, to) = (point_of[&ca.face], point_of[&cb.face]);
                let (pf, pt) = (built.vertices[from], built.vertices[to]);
                let curve = match centre { Some(c) => arc_between(c, pf, pt)?, None => Curve::Line { a: pf, b: pt } };
                let id = built.edge(curve, from, to);
                plan.chain.insert(s, id);
                edges.push(id);
            }
            plan.corner = Some((surface, sense, edges));
        }
        _ => unreachable!("a planned vertex has a blended edge"),
    }
    Ok(plan)
}

/// The ellipse where two equal fillet cylinders (axes c1 + s d1, c2 + s d2) meet, from t to q: they must cross, and be
/// mirror images in the plane through the crossing, t and q; the section of the first by that plane is then the
/// section of the second.
#[allow(clippy::too_many_arguments)]
fn mitre(c1: Vec3, d1: Vec3, c2: Vec3, d2: Vec3, r: f64, t: Vec3, q: Vec3, tol: f64, v: usize) -> Result<Curve, KernelError> {
    let p = meet((c1, d1), (c2, d2), tol, &format!("vertex {v}: the two fillets' axes"))?;
    let m = t.sub(p).cross(q.sub(p)).unit()
        .map_err(|_| KernelError::Unsupported(format!("vertex {v}: the fillets' joint is degenerate")))?;
    let mirrored = d1.sub(m.scale(2.0 * m.dot(d1)));
    if mirrored.cross(d2).length() > 1e-9 {
        return unsupported(format!("vertex {v}: the two fillets are not mirror images across their joint; their meeting curve is not a plane ellipse"));
    }
    cylinder_section(c1, d1, r, p, m, t, q)
}

/// A vertex on a plane-cylinder rim: the rim continues through it (both rim edges blended alike) and the third edge is
/// a straight seam along the cylinder.
#[allow(clippy::too_many_arguments)]
fn rim(v: usize, corners: &[Corner], ends: &[End], selected: &[End], blends: &HashMap<usize, EdgeBlend>, body: &Body,
    built: &mut Built, tol: f64) -> Result<Plan, KernelError> {
    if selected.len() != 2 {
        return unsupported(format!("vertex {v}: a plane-cylinder blend must run on round its rim (both rim edges at the vertex blended alike)"));
    }
    let common = *corners.iter().find(|c| selected.contains(&c.inn) && selected.contains(&c.out))
        .ok_or_else(|| KernelError::Unsupported(format!("vertex {v}: the rim edges there do not share a face")))?;
    let (s1, s2) = (common.inn, common.out);
    let u = *ends.iter().find(|x| !selected.contains(x)).expect("one unblended end");
    let (b1, b2) = (&blends[&s1.edge], &blends[&s2.edge]);
    let (Shape::Rim { c0, ah, radius, alpha, beta, sp, sc, plane, cylinder, .. }, Shape::Rim { c0: k0, ah: ak, radius: rk, alpha: al, beta: be, sp: spk, sc: sck, plane: pk, cylinder: cyk, .. }) = (&b1.shape, &b2.shape) else { unreachable!() };
    let same_cylinder = || {
        let (Surface::Cylinder { p, a, r, .. }, Surface::Cylinder { p: q, a: b, r: s, .. }) = (&body.faces[*cylinder].surface, &body.faces[*cyk].surface) else { return false };
        let w = q.sub(*p);
        a.cross(*b).length() <= 1e-9 && w.sub(a.scale(w.dot(*a))).length() <= tol && (r - s).abs() <= tol
    };
    let alike = plane == pk && c0.distance(*k0) <= tol && ah.sub(*ak).length() <= 1e-9 && (radius - rk).abs() <= tol && alpha == al && beta == be
        && (sp - spk).abs() <= tol && (sc - sck).abs() <= tol && b1.radius.is_some() == b2.radius.is_some() && same_cylinder();
    if !alike || common.face != *plane {
        return unsupported(format!("vertex {v}: the blends of rim edges {} and {} differ (kind, size or faces); a rim blend must be the same all round", s1.edge, s2.edge));
    }
    let (sl, sd) = line_of(body, u.edge).map_err(|_| KernelError::Unsupported(format!("vertex {v}: the rim's third edge {} is not a straight seam along the cylinder", u.edge)))?;
    if sd.cross(*ah).length() > 1e-9 || (sl.sub(*c0).sub(ah.scale(sl.sub(*c0).dot(*ah))).length() - radius).abs() > tol {
        return unsupported(format!("vertex {v}: the rim's third edge {} is not a straight seam along the cylinder", u.edge));
    }
    let at = body.vertices[v];
    let w = at.sub(*c0);
    let rho = w.sub(ah.scale(w.dot(*ah))).unit()?;
    let t = c0.add(rho.scale(radius + alpha * sp));
    let q = at.add(ah.scale(beta * sc));
    let curve = match b1.radius {
        Some(r) => arc_between(c0.add(ah.scale(beta * r)).add(rho.scale(radius + alpha * r)), t, q)?,
        None => Curve::Line { a: t, b: q },
    };
    let ca = *find(corners, |c| (c.inn == s1 || c.out == s1) && (c.inn == u || c.out == u));
    let cb = *find(corners, |c| (c.inn == s2 || c.out == s2) && (c.inn == u || c.out == u));
    let (ti, qi) = (built.vertex(t), built.vertex(q));
    let mut plan = Plan::default();
    plan.rail.insert((s1, common.face), ti);
    plan.rail.insert((s2, common.face), ti);
    plan.rail.insert((s1, ca.face), qi);
    plan.rail.insert((s2, cb.face), qi);
    plan.point.insert(u, qi);
    let id = built.edge(curve, ti, qi);
    plan.chain.insert(s1, id);
    plan.chain.insert(s2, id);
    Ok(plan)
}

//! A blend that runs out past a short end face onto the curved face beyond it (Claude 2026-11-03, Bob's L3x3x1/4
//! angle: a fillet along the leg's bottom end edge, whose end face at the toe is 0.0625 flat and then the toe round).
//!
//! The configuration, classified before anything is built (`classify`): a straight blend between two flat faces A and
//! B ends at vertex V on a flat end face H (one blended edge there, three faces). On one side, X (A or B), the rail
//! meets the line of H's edge g on X beyond g's far vertex W. At W three faces meet: X, H and a circular cylinder H2
//! whose axis b is square to X; H2 is a strip bounded by the seam with H (a line along b), its edge on X (an arc), a
//! second line along b and a second arc, and lies on one side of the plane through its axis along the blend.
//!
//! The blend then keeps its radius and runs on to its end faces, as a rolling-ball blend ending on its end faces does
//! (the idea follows solidcore's generic-blend notes, docs/generic-blend.md, Apache-2.0; no code was copied): g and W
//! are consumed; H gains the blend's section by H's plane from the other side's point Q to S, where that section
//! crosses the seam; H2 gains the blend surface's intersection with the cylinder from S to P, where X's rail meets
//! H2's arc on X. For a fillet that intersection is two circular cylinders with square axes (a `CylinderPair`); for a
//! chamfer a plane section of the cylinder (an ellipse arc).
//!
//! Proofs before return: P is the only point of X's rail on H2's arc; S is the only crossing of the section with the
//! seam, inside it; the `CylinderPair`'s discriminant stays positive (closed form); on H2's strip the new curve runs
//! monotonically across (from the seam to P) and down (from S's height to X), so it stays inside the strip and meets
//! its boundary only at S and P. The changed flat faces get the usual proofs in `mod.rs`.

use super::edge::{cylinder_section, plane_of, unsupported, EdgeBlend, Shape};
use super::proof;
use super::vertex::{Corner, End};
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::Body;
use crate::KernelError;
use std::collections::HashMap;
use std::f64::consts::{PI, TAU};

/// What lies beyond the short end edge on one side.
pub(crate) struct Beyond {
    /// The end face's edge on X (consumed) and its far vertex (consumed).
    pub g: usize,
    pub w: usize,
    /// The curved face beyond, its corner at W and the ends at W of its seam with H and its arc on X.
    pub h2: Corner,
    pub seam: End,
    pub arc: End,
    /// H2's cylinder.
    cp: Vec3,
    b: Vec3,
    radius: f64,
    /// H2's second line along the axis (a point on it) and its second arc's height.
    other: Vec3,
    far: f64,
}

/// Where the blend's side on X, running past H's edge `g_at_v` (its end at V), goes beyond: Some when the
/// configuration above holds, None when it does not (the caller then keeps the plain end, which refuses by name).
pub(crate) fn classify(eb: &EdgeBlend, x: usize, h: usize, g_at_v: End, body: &Body, all: &HashMap<usize, Vec<Corner>>,
    blends: &HashMap<usize, EdgeBlend>) -> Option<Beyond> {
    let Shape::Planar { d, .. } = eb.shape else { return None };
    let g = g_at_v.edge;
    let w = if g_at_v.at_end { body.edges[g].start } else { body.edges[g].end };
    let corners = all.get(&w)?;
    if corners.len() != 3 { return None; }
    let on = |f: usize| corners.iter().find(|c| c.face == f).copied();
    let (xc, hc) = (on(x)?, on(h)?);
    let h2 = *corners.iter().find(|c| c.face != x && c.face != h)?;
    let other_end = |c: &Corner| if c.inn.edge == g { c.out } else { c.inn };
    if xc.inn.edge != g && xc.out.edge != g { return None; }
    if hc.inn.edge != g && hc.out.edge != g { return None; }
    let (arc, seam) = (other_end(&xc), other_end(&hc));
    if !((h2.inn == arc && h2.out == seam) || (h2.inn == seam && h2.out == arc)) { return None; }
    if blends.contains_key(&arc.edge) || blends.contains_key(&seam.edge) { return None; }
    let Surface::Cylinder { p: cp, a: b, r: radius, .. } = body.faces[h2.face].surface else { return None };
    let (_, xn) = plane_of(body, x)?;
    if xn.cross(b).length() > 1e-9 || d.dot(b).abs() > 1e-12 { return None; }
    let along_axis = |e: usize| match body.edges[e].curve {
        Curve::Line { a, b: q } => q.sub(a).unit().ok().filter(|u| u.cross(b).length() <= 1e-9).map(|_| a),
        _ => None,
    };
    let square_arc = |e: usize| match body.edges[e].curve {
        Curve::Arc { c, n, r, .. } => {
            let w = c.sub(cp);
            (n.cross(b).length() <= 1e-9 && w.sub(b.scale(w.dot(b))).length() <= 1e-9 * (1.0 + radius) && (r - radius).abs() <= 1e-9 * (1.0 + radius))
                .then(|| w.dot(b))
        }
        _ => None,
    };
    along_axis(seam.edge)?;
    square_arc(arc.edge)?;
    // H2 is a strip: one loop of the seam, the arc, a second line along the axis and a second arc.
    let lp = &body.faces[h2.face].loops;
    if lp.len() != 1 || lp[0].len() != 4 { return None; }
    let rest: Vec<usize> = lp[0].iter().map(|&(e, _)| e).filter(|&e| e != seam.edge && e != arc.edge).collect();
    if rest.len() != 2 { return None; }
    let (other, far) = match (along_axis(rest[0]), square_arc(rest[1]), along_axis(rest[1]), square_arc(rest[0])) {
        (Some(o), Some(f), ..) | (_, _, Some(o), Some(f)) => (o, f),
        _ => return None,
    };
    Some(Beyond { g, w, h2, seam, arc, cp, b, radius, other, far })
}

/// The new points and curves of a run-out.
pub(crate) struct Runout {
    pub s: Vec3,
    pub p: Vec3,
    /// On H, from Q to S.
    pub section: Curve,
    /// On H2, from S to P.
    pub beyond: Curve,
}

/// f(t) = k + a cos t + b sin t over [t0, t1]: its least and greatest values, and whether it is monotone (no stationary
/// point strictly inside).
fn sinus(k: f64, a: f64, b: f64, t0: f64, t1: f64) -> (f64, f64, bool) {
    let f = |t: f64| k + a * t.cos() + b * t.sin();
    let (mut lo, mut hi) = (f(t0).min(f(t1)), f(t0).max(f(t1)));
    let mut monotone = true;
    if a.hypot(b) > 0.0 {
        let base = b.atan2(a);
        let mut t = base + ((t0 - base) / PI).ceil() * PI;
        while t <= t1 {
            if t > t0 + 1e-9 && t < t1 - 1e-9 { monotone = false; }
            lo = lo.min(f(t));
            hi = hi.max(f(t));
            t += PI;
        }
    }
    (lo, hi, monotone)
}

/// An arc or ellipse arc as c + u cos t + v sin t over [t0, t1].
fn trig(curve: &Curve) -> Option<(Vec3, Vec3, Vec3, f64, f64)> {
    match *curve {
        Curve::Arc { c, n, x, r, t0, t1 } => Some((c, x.scale(r), n.cross(x).scale(r), t0, t1)),
        Curve::Conic { c, u, v, t0, t1 } => Some((c, u, v, t0, t1)),
        _ => None,
    }
}

/// The run-out on side X: `q` is the other side's point on H, `past` where X's rail meets the line of g (beyond W),
/// `(hp, hn)` H's plane.
#[allow(clippy::too_many_arguments)]
pub(crate) fn build(v: usize, eb: &EdgeBlend, x: usize, by: &Beyond, body: &Body, q: Vec3, past: Vec3, (hp, hn): (Vec3, Vec3),
    tol: f64) -> Result<Runout, KernelError> {
    let Shape::Planar { o, d, axis, wa, wb, sa, sb } = eb.shape else { unreachable!("classified as a straight blend") };
    let e = eb.edge;
    let (cp, b, big) = (by.cp, by.b, by.radius);
    let wp = body.vertices[by.w];
    let refuse = |m: String| unsupported::<Runout>(format!("vertex {v}: the blend of edge {e} running past its end face onto face {}: {m}", by.h2.face));
    // Room across: the side on X must stop short of the curved face's far line (with the reach and the width the
    // message says what size fits).
    let (wx, reach) = if x == eb.a { (wa, sa) } else { (wb, sb) };
    let width = by.other.sub(o).dot(wx);
    if reach >= width - tol {
        return unsupported(format!("edge {e} is too wide for its end: its side on face {x} reaches {reach} across, past the curved face beyond the end ({width} across)"));
    }
    // P: X's rail meets H2's arc on X (a line and a circle in X's plane).
    let Curve::Arc { c: ac, n: an, x: ax, r: ar, t0: a0, t1: a1 } = body.edges[by.arc.edge].curve else { unreachable!("classified as an arc") };
    let (ox, _) = eb.rail_line(x);
    let wo = ox.sub(ac);
    let (wd, disc) = (wo.dot(d), wo.dot(d).powi(2) - wo.dot(wo) + ar * ar);
    if disc <= 0.0 { return refuse(format!("its side on face {x} misses the curved face")); }
    let ay = an.cross(ax);
    let hits: Vec<Vec3> = [-wd - disc.sqrt(), -wd + disc.sqrt()].iter().map(|&s| ox.add(d.scale(s))).filter(|pt| {
        let r = pt.sub(ac);
        let t = a0 + (r.dot(ay).atan2(r.dot(ax)) - a0).rem_euclid(TAU);
        (t - a0) * ar > tol && (a1 - t) * ar > tol
    }).collect();
    let [p] = hits[..] else { return refuse(format!("its side on face {x} meets the curved face's edge {} {} times", by.arc.edge, hits.len())) };
    // S: the section by H's plane crosses the seam (the line through W along b, in H).
    let (s, section) = match (eb.radius, axis) {
        (Some(r), Some(c)) => {
            let full = cylinder_section(c, d, r, hp, hn, q, past)?;
            let (cc, u, vv, t0, t1) = trig(&full).expect("a section is an arc or an ellipse arc");
            let size = u.length().max(vv.length());
            let roots: Vec<f64> = proof::trig_line_roots(cc, u, vv, t0, t1, wp, b, hn).into_iter()
                .filter(|&t| (t - t0) * size > tol && (t1 - t) * size > tol).collect();
            let [t] = roots[..] else { return refuse(format!("its end section crosses the seam (edge {}) {} times", by.seam.edge, roots.len())) };
            let s = Curve::conic_point(cc, u, vv, t);
            (s, cylinder_section(c, d, r, hp, hn, q, s)?)
        }
        _ => {
            let dir = past.sub(q);
            let n = dir.cross(b);
            let nn = n.dot(n);
            if nn <= 1e-18 { return refuse("its section runs along the seam".into()); }
            let lambda = wp.sub(q).cross(b).dot(n) / nn;
            if lambda * dir.length() <= tol || (1.0 - lambda) * dir.length() <= tol { return refuse("its section does not cross the seam".into()); }
            let s = q.add(dir.scale(lambda));
            (s, Curve::Line { a: q, b: s })
        }
    };
    let seam_far = { let se = &body.edges[by.seam.edge]; body.vertices[if by.seam.at_end { se.start } else { se.end }] };
    let seam_len = seam_far.distance(wp);
    let along = s.sub(wp).dot(seam_far.sub(wp)) / seam_len;
    if along <= tol || along >= seam_len - tol { return refuse(format!("it crosses the seam outside it ({along} along a seam {seam_len} long)")); }
    // Across (e) and along the axis (b) on H2, and the side of the axis plane through d.
    let ec = b.cross(d).unit()?;
    let side = wp.sub(cp).dot(d).signum();
    if wp.sub(cp).dot(d).abs() <= tol { return refuse("the seam lies on the plane through the cylinder's axis along the blend".into()); }
    for &arc_edge in body.faces[by.h2.face].loops[0].iter().map(|(x, _)| x).filter(|&&x| matches!(body.edges[x].curve, Curve::Arc { .. })) {
        let (c, u, vv, t0, t1) = trig(&body.edges[arc_edge].curve).expect("an arc");
        let (lo, hi, _) = sinus(c.sub(cp).dot(d), u.dot(d), vv.dot(d), t0, t1);
        if (if side > 0.0 { lo } else { -hi }) < -tol { return refuse("the curved face passes round its axis".into()); }
    }
    // The curve on H2 from S to P, with its across and along-axis coordinates as k + a cos t + b sin t over [t0, t1].
    let (beyond, coord_e, coord_z, span) = match (eb.radius, axis) {
        (Some(r), Some(c)) => {
            let xs = s.sub(c).sub(d.scale(s.sub(c).dot(d))).unit()?;
            let mut a = d;
            let mut tp = { let y = a.cross(xs); let w = p.sub(c); w.dot(y).atan2(w.dot(xs)) };
            if tp < 0.0 { a = d.scale(-1.0); tp = -tp; }
            if tp * r <= tol || tp >= PI { return refuse("the fillet's end on the curved face is degenerate".into()); }
            let y = a.cross(xs);
            let coord = |m: Vec3| (c.sub(cp).dot(m), xs.dot(m) * r, y.dot(m) * r);
            let (ce, cz) = (coord(ec), coord(b));
            // The pair's discriminant is r_H2^2 - (across)^2 (the axes are square): positive while |across| < r_H2.
            let (lo, hi, _) = sinus(ce.0, ce.1, ce.2, 0.0, tp);
            if lo.abs().max(hi.abs()) >= big - tol { return refuse("the fillet reaches round the curved face's axis".into()); }
            let pair = |branch: f64| Curve::CylinderPair { c, a, x: xs, rho: r, p: cp, b, r: big, branch, t0: 0.0, t1: tp };
            let fits = |cv: &Curve| cv.eval(0.0).0.distance(s) <= tol && cv.eval(tp).0.distance(p) <= tol;
            let found: Vec<Curve> = [1.0, -1.0].into_iter().map(pair).filter(fits).collect();
            let [cv] = &found[..] else { return Err(KernelError::Internal(format!("vertex {v}: the fillet's curve on face {} does not join S and P on one branch", by.h2.face))) };
            (cv.clone(), ce, cz, (0.0, tp))
        }
        _ => {
            let Surface::Plane { p: pp, n: pn, .. } = eb.surface else { unreachable!("a chamfer is a plane") };
            let cv = cylinder_section(cp, b, big, pp, pn, s, p)?;
            let (c, u, vv, t0, t1) = trig(&cv).expect("a section is an arc or an ellipse arc");
            let coord = |m: Vec3| (c.sub(cp).dot(m), u.dot(m), vv.dot(m));
            let cd = coord(d);
            let (lo, hi, _) = sinus(cd.0, cd.1, cd.2, t0, t1);
            if (if side > 0.0 { lo } else { -hi }) < -tol { return refuse("the chamfer reaches round the curved face's axis".into()); }
            (cv, coord(ec), coord(b), (t0, t1))
        }
    };
    // Inside the strip: across from the seam (at S) toward the second line, along the axis from S's height down to X
    // (at P), both monotone; so the curve meets the strip's boundary only at S and P.
    let (_, _, mono_e) = sinus(coord_e.0, coord_e.1, coord_e.2, span.0, span.1);
    let (_, _, mono_z) = sinus(coord_z.0, coord_z.1, coord_z.2, span.0, span.1);
    let at = |k: (f64, f64, f64), t: f64| k.0 + k.1 * t.cos() + k.2 * t.sin();
    let (ys, yp, yo) = (at(coord_e, span.0), at(coord_e, span.1), by.other.sub(cp).dot(ec));
    let (zs, zx, zf) = (at(coord_z, span.0), wp.sub(cp).dot(b), by.far);
    let across = (yo - ys).signum();
    let down = (zf - zx).signum();
    if !mono_e || !mono_z { return refuse("its curve on the curved face turns back".into()); }
    if across * (yp - ys) <= tol { return refuse("its curve on the curved face does not run across it".into()); }
    if across * (yo - yp) <= tol {
        return refuse(format!("it runs across the whole curved face (it reaches {} across; the face is {} wide)", (yp - ys).abs(), (yo - ys).abs()));
    }
    if down * (zs - zx) <= tol || down * (zf - zs) <= tol { return refuse("its curve on the curved face leaves the face's height".into()); }
    Ok(Runout { s, p, section, beyond })
}

//! Edge fillets and chamfers (Claude 2026-11-01, Bob: "continue with ... the fillets and the chamfers"). Exact
//! analytic blends; anything outside the supported domain is refused by name, and nothing partial is returned.
//!
//! Supported exactly:
//! - straight edges between two flat faces, convex or concave: a chamfer (equal distances, two distances, or a
//!   distance and an angle) is a plane, a fillet a circular cylinder;
//! - circles where a flat face meets a cylinder at right angles (hole rims, bar and pipe ends, boss bases): a chamfer
//!   is a cone, a fillet a torus; the rim must be blended all round, alike;
//! - their ends at vertices where three faces meet (see `vertex`): running out on a flat end face (or past a short one
//!   onto the cylinder strip beyond it: `runout`), two blends meeting
//!   in a mitre (chamfer planes in a line; equal fillets in an ellipse), three meeting at a corner (a plane triangle
//!   for chamfers; a sphere patch for equal fillets where one face is square to the other two).
//! - convex outer chamfers along a tangent straight/circular rim on one end plane, with matching setbacks on
//!   the end and skin. The plane and cone join along a common generator. Mixed fillets and non-tangent joins refuse.
//!
//! Every result is proved before it is returned: the rails and trims stay inside their faces, each changed flat face's
//! curves are proved apart (except where they join) and its loops nest as before, no other edge of the body enters a
//! blend's region, and `topology::check` passes. A blend that would consume a face, cross another loop or cut into
//! another feature is refused by name.
//!
//! Sources (public mathematics only): the rolling-ball construction of a constant-radius blend between two surfaces
//! (the ball's centre on both offset surfaces; for planes the offset planes' line, for a plane and a coaxial cylinder
//! the offset circle), plane sections of a circular cylinder (circles and ellipses), the intersection of two equal
//! circular cylinders with crossing axes (two plane ellipses, in the planes of mirror symmetry), the sine rule for the
//! distance-angle chamfer, and Lipschitz/convex-hull bounds for the clearance proofs (README.md).

mod edge;
mod proof;
mod runout;
mod vertex;

use crate::extrude::signed_area;
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Face};
use crate::{topology, KernelError};
use edge::{plan_edge, unsupported, EdgeBlend, Shape};
use std::collections::{BTreeSet, HashMap};
use vertex::{Built, Corner, End, Plan};

/// What to do to one edge. Distances are in the body's unit; angles in radians.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Kind {
    /// A constant-radius fillet.
    Fillet { radius: f64 },
    /// A chamfer meeting both faces at `distance` from the edge.
    Chamfer { distance: f64 },
    /// A chamfer meeting face `face` (input body face index; one of the edge's two faces) at `first` and the other face
    /// at `second`.
    ChamferDistances { first: f64, second: f64, face: usize },
    /// A chamfer meeting face `face` at `distance` and making `angle` with it.
    ChamferAngle { distance: f64, angle: f64, face: usize },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Item {
    /// The edge's index in the input body.
    pub edge: usize,
    pub kind: Kind,
}

/// Where a result face came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Origin {
    /// An input face, kept (trimmed where blends meet it). Wire role 0.
    Face(usize),
    /// The blend of an input edge. Wire role 1.
    Edge(usize),
    /// The corner patch where three blends meet at an input vertex. Wire role 2.
    Vertex(usize),
}

pub struct Blended {
    pub body: Body,
    /// One per result face.
    pub origins: Vec<Origin>,
}

fn invalid<T>(m: impl Into<String>) -> Result<T, KernelError> { Err(KernelError::InvalidInput(m.into())) }

fn validate(body: &Body, items: &[Item]) -> Result<(), KernelError> {
    if items.is_empty() { return invalid("no edges to blend"); }
    let mut seen = BTreeSet::new();
    for item in items {
        let e = item.edge;
        if e >= body.edges.len() { return invalid(format!("edge {e} is not in the body ({} edges)", body.edges.len())); }
        if !seen.insert(e) { return invalid(format!("edge {e} is listed twice")); }
        let positive = |x: f64, what: &str| if x.is_finite() && x > 0.0 { Ok(()) } else { invalid(format!("edge {e}: the {what} must be positive, not {x}")) };
        match item.kind {
            Kind::Fillet { radius } => positive(radius, "radius")?,
            Kind::Chamfer { distance } => positive(distance, "distance")?,
            Kind::ChamferDistances { first, second, face } => {
                positive(first, "first distance")?;
                positive(second, "second distance")?;
                if face >= body.faces.len() { return invalid(format!("edge {e}: face {face} is not in the body")); }
            }
            Kind::ChamferAngle { distance, angle, face } => {
                positive(distance, "distance")?;
                if !(angle.is_finite() && angle > 0.0 && angle < std::f64::consts::PI) {
                    return invalid(format!("edge {e}: the chamfer angle must lie strictly between 0 and pi radians, not {angle}"));
                }
                if face >= body.faces.len() { return invalid(format!("edge {e}: face {face} is not in the body")); }
            }
        }
    }
    Ok(())
}

/// Fillets and chamfers the listed edges of a closed body. The body must pass `topology::check`.
pub fn blend(body: &Body, items: &[Item], tol: f64) -> Result<Blended, KernelError> {
    if !(tol.is_finite() && tol > 0.0) { return invalid("the tolerance must be positive"); }
    validate(body, items)?;
    topology::check(body, tol)?;
    let mut used: Vec<[Option<usize>; 2]> = vec![[None, None]; body.edges.len()];
    for (f, face) in body.faces.iter().enumerate() {
        for &(e, fwd) in face.loops.iter().flatten() { used[e][if fwd { 0 } else { 1 }] = Some(f); }
    }
    let mut blends: HashMap<usize, EdgeBlend> = HashMap::new();
    for item in items {
        let [a, b] = used[item.edge];
        blends.insert(item.edge, plan_edge(body, item, a.expect("checked"), b.expect("checked"), tol)?);
    }
    let selected: BTreeSet<usize> = blends.keys().copied().collect();
    // Every face corner, by vertex.
    let mut at_vertex: HashMap<usize, Vec<Corner>> = HashMap::new();
    for (f, face) in body.faces.iter().enumerate() {
        for (l, lp) in face.loops.iter().enumerate() {
            let n = lp.len();
            for k in 0..n {
                let (pe, pf) = lp[(k + n - 1) % n];
                let (ne, nf) = lp[k];
                let inn = End { edge: pe, at_end: pf };
                let out = End { edge: ne, at_end: !nf };
                let v = if pf { body.edges[pe].end } else { body.edges[pe].start };
                at_vertex.entry(v).or_default().push(Corner { face: f, lp: l, k, inn, out });
            }
        }
    }
    let touched: BTreeSet<usize> = selected.iter().flat_map(|&e| [body.edges[e].start, body.edges[e].end]).collect();
    let mut built = Built::default();
    let vmap: Vec<Option<usize>> = (0..body.vertices.len())
        .map(|v| (!touched.contains(&v)).then(|| built.vertex(body.vertices[v]))).collect();
    let mut plans: HashMap<usize, Plan> = HashMap::new();
    for &v in &touched {
        plans.insert(v, vertex::plan(v, &at_vertex[&v], &at_vertex, &blends, body, &mut built, tol)?);
    }
    // What the plans move, insert and consume, by edge end and face corner (a run-out reaches past its own vertex).
    let moved: HashMap<End, usize> = plans.values().flat_map(|p| p.point.iter().map(|(&x, &id)| (x, id))).collect();
    let inserts: HashMap<(usize, usize, usize), (usize, bool)> = plans.values().flat_map(|p| p.insert.iter().map(|(&k, &x)| (k, x))).collect();
    let consumed: BTreeSet<usize> = plans.values().flat_map(|p| p.consumed_edges.iter().copied()).collect();
    let gone: Vec<(usize, usize)> = plans.iter().flat_map(|(&v, p)| p.consumed_vertices.iter().map(move |&w| (w, v))).collect();
    let end_vertex = |x: End| if x.at_end { body.edges[x.edge].end } else { body.edges[x.edge].start };
    // Unblended edges, trimmed where a blend meets them.
    let mut is_new: Vec<bool> = vec![true; built.edges.len()];
    let mut edge_map: Vec<Option<usize>> = vec![None; body.edges.len()];
    for (e, edge) in body.edges.iter().enumerate() {
        if selected.contains(&e) || consumed.contains(&e) { continue; }
        let id_of = |x: End| match moved.get(&x) {
            Some(&id) => Ok((id, true)),
            None => vmap[end_vertex(x)].map(|id| (id, false))
                .ok_or_else(|| KernelError::Internal(format!("edge {e}: its end at a blended vertex was not planned"))),
        };
        let ((s, ts), (t, tt)) = (id_of(End { edge: e, at_end: false })?, id_of(End { edge: e, at_end: true })?);
        let curve = if !ts && !tt { edge.curve.clone() } else if let Curve::Arc { c, n, x, r, t0, t1 } = edge.curve {
            // A run-out trims the curved face's arc on its flat face: the new ends' angles, inside the old span.
            let y = n.cross(x);
            let angle = |p: Vec3| { let w = p.sub(c); t0 + (w.dot(y).atan2(w.dot(x)) - t0 + 1e-12).rem_euclid(std::f64::consts::TAU) - 1e-12 };
            let (u0, u1) = (if ts { angle(built.vertices[s]) } else { t0 }, if tt { angle(built.vertices[t]) } else { t1 });
            if u0 < t0 - 1e-12 || u1 > t1 + 1e-12 || (u1 - u0) * r <= tol {
                return unsupported(format!("edge {e} would be consumed: the blend reaching it passes its end"));
            }
            Curve::Arc { c, n, x, r, t0: u0, t1: u1 }
        } else {
            let Curve::Line { a, b } = edge.curve else { unreachable!("vertex plans trim only straight edges and run-out arcs") };
            let (ps, pt) = (built.vertices[s], built.vertices[t]);
            let (d, len) = (b.sub(a), b.distance(a));
            let along = |p: Vec3| p.sub(a).dot(d) / (len * len);
            let (u0, u1) = (along(ps), along(pt));
            if u0 < -1e-12 || u1 > 1.0 + 1e-12 || (u1 - u0) * len <= tol {
                // How far along this edge the blends reach in all; with the length it says what size would fit.
                let reach = len * (1.0 - (u1 - u0));
                return unsupported(format!("edge {e} would be consumed: the blends reaching it pass its end (it is {len} long; they reach {reach} along it)"));
            }
            Curve::Line { a: ps, b: pt }
        };
        let id = built.edge(curve, s, t);
        is_new.push(ts || tt);
        edge_map[e] = Some(id);
    }
    // Rails.
    let mut rails: HashMap<(usize, usize), usize> = HashMap::new();
    for &e in &selected {
        let eb = &blends[&e];
        let (sv, ev) = (body.edges[e].start, body.edges[e].end);
        for face in [eb.a, eb.b] {
            let s = plans[&sv].rail[&(End { edge: e, at_end: false }, face)];
            let t = plans[&ev].rail[&(End { edge: e, at_end: true }, face)];
            let (ps, pt) = (built.vertices[s], built.vertices[t]);
            if let Shape::Planar { d, .. } = eb.shape {
                if pt.sub(ps).dot(d) <= tol {
                    return unsupported(format!("edge {e}: the blend's sides cross on face {face}; the edge is too short for it"));
                }
            }
            let id = built.edge(eb.rail(face, ps, pt), s, t);
            rails.insert((e, face), id);
        }
    }
    is_new.resize(built.edges.len(), true);
    // Faces: the input faces with their new loops, then one blend face per edge, then corner patches.
    let mut faces = Vec::new();
    let mut origins = Vec::new();
    let mut changed = Vec::new();
    for (f, face) in body.faces.iter().enumerate() {
        let mut loops = Vec::new();
        for (l, lp) in face.loops.iter().enumerate() {
            let mut nl = Vec::new();
            for (k, &(e, fwd)) in lp.iter().enumerate() {
                if let Some(&(id, dir)) = inserts.get(&(f, l, k)) { nl.push((id, dir)); }
                if consumed.contains(&e) { continue; }
                nl.push((if selected.contains(&e) { rails[&(e, f)] } else { edge_map[e].expect("mapped") }, fwd));
            }
            loops.push(nl);
        }
        changed.push(loops.iter().flatten().any(|&(e, _)| is_new[e]));
        faces.push(Face { surface: face.surface.clone(), sense: face.sense, loops });
        origins.push(Origin::Face(f));
    }
    let mut chain_use: HashMap<usize, bool> = HashMap::new();
    for &e in &selected {
        let eb = &blends[&e];
        let (sv, ev) = (body.edges[e].start, body.edges[e].end);
        let (s_end, e_end) = (End { edge: e, at_end: false }, End { edge: e, at_end: true });
        let ra = rails[&(e, eb.a)];
        let rb = rails[&(e, eb.b)];
        let cs = walk(&built, &plans[&sv].chain[&s_end], plans[&sv].rail[&(s_end, eb.a)], e)?;
        let ce = walk(&built, &plans[&ev].chain[&e_end], plans[&ev].rail[&(e_end, eb.b)], e)?;
        for &(c, d) in cs.iter().chain(&ce) { chain_use.insert(c, d); }
        let mut lp = vec![(ra, false)];
        lp.extend(cs);
        lp.push((rb, true));
        lp.extend(ce);
        faces.push(Face { surface: eb.surface.clone(), sense: eb.sense, loops: vec![lp] });
        origins.push(Origin::Edge(e));
        changed.push(false);
    }
    for &v in &touched {
        let Some((surface, sense, edges)) = &plans[&v].corner else { continue };
        // Each patch edge runs against its use in the blend face beside it; chained head to tail.
        let directed: Vec<(usize, bool)> = edges.iter().map(|&e| (e, !chain_use[&e])).collect();
        let head = |(e, d): (usize, bool)| if d { built.edges[e].start } else { built.edges[e].end };
        let tail = |(e, d): (usize, bool)| if d { built.edges[e].end } else { built.edges[e].start };
        let mut lp = vec![directed[0]];
        while lp.len() < directed.len() {
            let last = tail(*lp.last().unwrap());
            let next = directed.iter().copied().find(|&x| head(x) == last && !lp.contains(&x))
                .ok_or_else(|| KernelError::Internal(format!("vertex {v}: the corner patch's edges do not close")))?;
            lp.push(next);
        }
        faces.push(Face { surface: surface.clone(), sense: *sense, loops: vec![lp] });
        origins.push(Origin::Vertex(v));
        changed.push(false);
    }
    let mut result = Body { vertices: built.vertices, edges: built.edges, faces };
    if !gone.is_empty() { drop_unused_vertices(&mut result); }
    // Proofs.
    for (f, face) in result.faces.iter().enumerate() {
        if !changed[f] { continue; }
        match face.surface {
            Surface::Plane { n, .. } => check_flat_face(&result, f, n.scale(if face.sense { 1.0 } else { -1.0 }), &is_new, tol)?,
            Surface::Cylinder { .. } => {}
            _ => return Err(KernelError::Internal(format!("face {f} changed but is neither flat nor a rim's cylinder"))),
        }
    }
    for eb in blends.values() {
        if let Shape::Rim { .. } = eb.shape { check_rim_cylinder(&result, eb, &is_new, tol)?; }
    }
    // Each straight blend's reach along its edge: its blend face's corners (every end section, mitre and corner arc
    // lies in the cross-section triangle's plane cut through them, so between their values).
    let mut reach: HashMap<usize, (f64, f64)> = HashMap::new();
    for (f, origin) in origins.iter().enumerate() {
        let Origin::Edge(e) = *origin else { continue };
        let Shape::Planar { o, d, .. } = blends[&e].shape else { continue };
        let along = |p: Vec3| p.sub(o).dot(d);
        let points = result.faces[f].loops.iter().flatten().flat_map(|&(x, _)| [result.edges[x].start, result.edges[x].end])
            .map(|v| along(result.vertices[v])).chain([along(body.vertices[body.edges[e].start]), along(body.vertices[body.edges[e].end])]);
        reach.insert(e, points.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), a| (lo.min(a), hi.max(a))));
    }
    check_interference(body, &blends, &selected, &gone, &reach, tol)?;
    topology::check(&result, tol).map_err(|e| KernelError::Internal(format!("the blend built an invalid solid: {e}")))?;
    Ok(Blended { body: result, origins })
}

fn directed(body: &Body, e: usize, fwd: bool) -> Curve {
    if fwd { body.edges[e].curve.clone() } else { body.edges[e].curve.reversed() }
}

/// A changed flat face: every new or trimmed curve is proved apart from every other curve it does not join, joints do
/// not fold back or cross again, and the loops nest as one outer loop (counter-clockwise) holding its holes.
fn check_flat_face(body: &Body, f: usize, normal: Vec3, is_new: &[bool], tol: f64) -> Result<(), KernelError> {
    let face = &body.faces[f];
    let segs: Vec<(usize, usize, usize, Curve)> = face.loops.iter().enumerate()
        .flat_map(|(l, lp)| lp.iter().enumerate().map(move |(k, &(e, fwd))| (l, k, e, fwd)))
        .map(|(l, k, e, fwd)| (l, k, e, directed(body, e, fwd))).collect();
    let refuse = |m: String| KernelError::Unsupported(format!("face {f}: {m}"));
    for i in 0..segs.len() {
        for j in i + 1..segs.len() {
            let ((li, ki, ei, ci), (lj, kj, ej, cj)) = (&segs[i], &segs[j]);
            if !(is_new[*ei] || is_new[*ej]) { continue; }
            let len = face.loops[*li].len();
            let consecutive = li == lj && (kj - ki == 1 || (*ki == 0 && *kj == len - 1));
            if consecutive {
                let (first, second) = if kj - ki == 1 { (ci, cj) } else { (cj, ci) };
                joint_ok(first, second, normal, tol).map_err(|m| refuse(format!("edges {ei} and {ej}: {m}")))?;
            } else {
                match proof::apart(ci, cj, tol) {
                    Ok(true) => {}
                    Ok(false) => return Err(refuse(format!("the blend makes edges {ei} and {ej} meet or cross (the blend runs into another boundary of the face)"))),
                    Err(m) => return Err(refuse(format!("could not prove edges {ei} and {ej} apart: {m}"))),
                }
            }
        }
    }
    let curves: Vec<Vec<Curve>> = face.loops.iter().map(|lp| lp.iter().map(|&(e, fwd)| directed(body, e, fwd)).collect()).collect();
    let areas: Vec<f64> = curves.iter().map(|c| signed_area(c, normal)).collect();
    let outer: Vec<usize> = (0..areas.len()).filter(|&i| areas[i] > 0.0).collect();
    if outer.len() != 1 { return Err(refuse("the blend leaves its loops without one outer boundary".into())); }
    for (h, lp) in curves.iter().enumerate() {
        if h == outer[0] { continue; }
        let q = lp[0].start();
        for (l, other) in curves.iter().enumerate() {
            if l == h { continue; }
            let w = proof::winding(q, normal, other).map_err(|m| refuse(format!("loop nesting: {m}")))?;
            let want = if l == outer[0] { 1.0 } else { 0.0 };
            if (w - want).abs() > 1e-6 { return Err(refuse(format!("the blend leaves loop {h} outside the face or inside another hole"))); }
        }
    }
    Ok(())
}

/// Two curves joined head to tail at a new vertex: they must not fold back on each other, nor cross again.
fn joint_ok(first: &Curve, second: &Curve, normal: Vec3, tol: f64) -> Result<(), String> {
    let din = first.eval(first.domain().1).1.unit().map_err(|e| e.to_string())?;
    let dout = second.eval(second.domain().0).1.unit().map_err(|e| e.to_string())?;
    if din.dot(dout) < -1.0 + 1e-9 { return Err("they fold back on each other".into()); }
    let trig = |c: &Curve| match *c {
        Curve::Arc { c, n, x, r, t0, t1 } => Some((c, x.scale(r), n.cross(x).scale(r), t0, t1, r)),
        Curve::Conic { c, u, v, t0, t1 } => Some((c, u, v, t0, t1, u.length().max(v.length()))),
        _ => None,
    };
    match (first, second) {
        (Curve::Line { .. }, Curve::Line { .. }) => Ok(()),
        (Curve::Line { a, b }, other) | (other, Curve::Line { a, b }) => {
            let (c, u, v, s0, s1, size) = trig(other).ok_or("no proof for this joint")?;
            let ld = b.sub(*a).unit().map_err(|e| e.to_string())?;
            let len = b.distance(*a);
            let joint = if std::ptr::eq(other, second) { s0 } else { s1 };
            for t in proof::trig_line_roots(c, u, v, s0, s1, *a, ld, normal) {
                if (t - joint).abs() * size <= tol { continue; }
                let p = Curve::conic_point(c, u, v, t);
                let along = p.sub(*a).dot(ld);
                if along > -tol && along < len + tol { return Err("they cross again".into()); }
            }
            Ok(())
        }
        _ => match (trig(first), trig(second)) {
            (Some((c1, u1, v1, ..)), Some((c2, u2, v2, ..))) if c1.distance(c2) <= tol
                && (u1.length() - u2.length()).abs() <= tol && (v1.length() - v2.length()).abs() <= tol
                && (u1.length() - v1.length()).abs() <= tol => Ok(()), // one circle continuing
            _ => Err("no proof for this joint".into()),
        },
    }
}

/// A rim blend's cylinder face: every curve of it other than the new rails and trimmed seams must lie beyond the rail
/// (farther from the flat face along the cylinder).
fn check_rim_cylinder(body: &Body, eb: &EdgeBlend, is_new: &[bool], tol: f64) -> Result<(), KernelError> {
    let Shape::Rim { c0, ah, beta, sc, cylinder, .. } = eb.shape else { return Ok(()) };
    for &(e, _) in body.faces[cylinder].loops.iter().flatten() {
        if is_new[e] { continue; }
        let h = |p: Vec3| beta * p.sub(c0).dot(ah);
        let least = match body.edges[e].curve {
            Curve::Line { a, b } => h(a).min(h(b)),
            Curve::Arc { c, n, x, r, .. } => h(c) - r * (x.dot(ah).abs() + n.cross(x).dot(ah).abs()),
            Curve::Conic { c, u, v, .. } => h(c) - u.dot(ah).abs() - v.dot(ah).abs(),
            Curve::Bezier { p } | Curve::RationalBezier { p, .. } => p.iter().map(|&q| h(q)).fold(f64::INFINITY, f64::min),
            Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } => return unsupported(format!("edge {}: the rim blend's cylinder face {cylinder} also holds an intersection curve (edge {e}); a clearance proof for it is not available", eb.edge)),
        };
        if least <= sc + tol {
            return unsupported(format!("edge {}: the blend reaches {sc} along the cylinder face {cylinder}, into its other boundary (input edge near {least})", eb.edge));
        }
    }
    Ok(())
}

/// No edge of the input body other than those at the blended chain's own vertices may enter a blend's region: the
/// triangle the blend cuts from (or adds to) the cross-section, swept along the edge (for a rim, the annulus it sweeps).
fn check_interference(body: &Body, blends: &HashMap<usize, EdgeBlend>, selected: &BTreeSet<usize>, gone: &[(usize, usize)],
    reach: &HashMap<usize, (f64, f64)>, tol: f64) -> Result<(), KernelError> {
    // The chains of blended edges joined at vertices.
    let mut parent: Vec<usize> = (0..body.vertices.len()).collect();
    fn root(p: &mut Vec<usize>, mut v: usize) -> usize { while p[v] != v { p[v] = p[p[v]]; v = p[v]; } v }
    for &e in selected {
        let (a, b) = (root(&mut parent, body.edges[e].start), root(&mut parent, body.edges[e].end));
        parent[a] = b;
    }
    // A run-out's consumed vertex joins its vertex's chain: the edges there are trimmed on faces whose new curves are
    // proved (the flat ones in check_flat_face, the curved one in runout::build).
    for &(w, v) in gone {
        let (a, b) = (root(&mut parent, w), root(&mut parent, v));
        parent[a] = b;
    }
    let touched: BTreeSet<usize> = selected.iter().flat_map(|&e| [body.edges[e].start, body.edges[e].end])
        .chain(gone.iter().map(|&(w, _)| w)).collect();
    for &e in selected {
        let eb = &blends[&e];
        let chain = root(&mut parent, body.edges[e].start);
        let near: Vec<bool> = (0..body.vertices.len()).map(|v| touched.contains(&v) && root(&mut parent, v) == chain).collect();
        let test: Box<dyn Fn(Vec3, f64) -> Option<bool>> = match eb.shape {
            Shape::Planar { o, d, wa, wb, sa, sb, .. } => {
                let (pa, pb) = (o.add(wa.scale(sa)), o.add(wb.scale(sb)));
                let (lo, hi) = reach[&e];
                // Side planes of the cross-section triangle (o, pa, pb), outward normals across d.
                let inside = o.add(pa).add(pb).scale(1.0 / 3.0);
                let sides: Vec<(Vec3, Vec3)> = [(o, pa), (pa, pb), (pb, o)].iter().map(|&(p, q)| {
                    let mut n = d.cross(q.sub(p)).unit().unwrap_or(d);
                    if n.dot(inside.sub(p)) > 0.0 { n = n.scale(-1.0); }
                    (p, n)
                }).collect();
                Box::new(move |c: Vec3, r: f64| {
                    let along = c.sub(o).dot(d);
                    let worst = sides.iter().map(|&(p, n)| c.sub(p).dot(n)).fold(f64::NEG_INFINITY, f64::max)
                        .max(lo - tol - along).max(along - hi - tol);
                    if worst > r { Some(true) } else if r <= tol / 4.0 { Some(false) } else { None }
                })
            }
            Shape::Rim { c0, ah, n, x, t0, t1, radius, alpha, beta, sp, sc, .. } => {
                let (r0, r1) = (radius.min(radius + alpha * sp), radius.max(radius + alpha * sp));
                let (h0, h1) = (0f64.min(beta * sc), 0f64.max(beta * sc));
                // A <= half-turn rim occupies a convex angular sector, not the
                // complete ring. Both inward radial planes contain that sector.
                // Account for any admitted arc-axis/skin-axis residual over the
                // axial range; no curve is exempted from the clearance proof.
                let sector = if t1 > t0 && t1 - t0 <= std::f64::consts::PI {
                    let y = n.cross(x);
                    let start = x.scale(t0.cos()).add(y.scale(t0.sin()));
                    let end = x.scale(t1.cos()).add(y.scale(t1.sin()));
                    let left = n.cross(start).unit()?;
                    let right = end.cross(n).unit()?;
                    Some((left, right, h0.abs().max(h1.abs())))
                } else { None };
                Box::new(move |c: Vec3, r: f64| {
                    let w = c.sub(c0);
                    let h = w.dot(ah);
                    let rho = w.sub(ah.scale(h)).length();
                    let mut worst = (h0 - h).max(h - h1).max(r0 - rho).max(rho - r1);
                    if let Some((left, right, axial)) = sector {
                        worst = worst.max(-w.dot(left) - axial * ah.dot(left).abs() - tol)
                            .max(-w.dot(right) - axial * ah.dot(right).abs() - tol);
                    }
                    if worst > r { Some(true) } else if r <= tol / 4.0 { Some(false) } else { None }
                })
            }
        };
        for (x, other) in body.edges.iter().enumerate() {
            if selected.contains(&x) || near[other.start] || near[other.end] { continue; }
            match proof::clear_of(&other.curve, &*test) {
                Ok(true) => {}
                Ok(false) => return unsupported(format!("the blend of edge {e} would run into edge {x} of the body")),
                Err(m) => return unsupported(format!("could not prove the blend of edge {e} clear of edge {x}: {m}")),
            }
        }
    }
    Ok(())
}

/// A blend face's closing edges at one end, from the rail vertex `from`, each with the direction it runs.
fn walk(built: &Built, chain: &[usize], from: usize, e: usize) -> Result<Vec<(usize, bool)>, KernelError> {
    let mut at = from;
    let mut out: Vec<(usize, bool)> = Vec::new();
    for _ in 0..chain.len() {
        let next = chain.iter().copied()
            .find(|&c| !out.iter().any(|&(o, _)| o == c) && (built.edges[c].start == at || built.edges[c].end == at))
            .ok_or_else(|| KernelError::Internal(format!("edge {e}: the blend face's end edges do not chain")))?;
        let fwd = built.edges[next].start == at;
        at = if fwd { built.edges[next].end } else { built.edges[next].start };
        out.push((next, fwd));
    }
    Ok(out)
}

/// Removes vertices no edge uses (a run-out's consumed vertex), keeping the others in order.
fn drop_unused_vertices(body: &mut Body) {
    let mut used = vec![false; body.vertices.len()];
    for e in &body.edges {
        used[e.start] = true;
        used[e.end] = true;
    }
    let mut map = vec![usize::MAX; body.vertices.len()];
    let mut kept = Vec::new();
    for (i, &p) in body.vertices.iter().enumerate() {
        if used[i] {
            map[i] = kept.len();
            kept.push(p);
        }
    }
    body.vertices = kept;
    for e in &mut body.edges {
        e.start = map[e.start];
        e.end = map[e.end];
    }
}

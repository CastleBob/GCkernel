//! Union domain: a boss standing on one planar face of the body (Claude 2026-11-02, op 24): a SolidWorks-style boss
//! extrude or revolve sketched on a plate, a web, a flange or a leg, anywhere on any body.
//!
//! Admission, decided from the geometry before anything is built:
//! 1. the base: every tool face that is a plane on the plane of a body face F (within the tolerance) facing the opposite
//!    way. Their union (edges shared by two base faces cancelled, e.g. the halves of a revolved disc) is chained into
//!    one outer loop and any holes (a ring boss). The base faces must lie on one plane;
//! 2. the base lies strictly inside F: in F's plane (its own frame, x and n x x) no edge of the base meets an edge of F
//!    (the section domain's closed-form meetings, the tolerance included), the base's outer loop lies inside F's region
//!    (winding 1) and no loop of F lies inside the base's outer loop;
//! 3. the boss lies above F's plane and meets it only in its base: every tool vertex at a height of at least -tol, those
//!    within the tolerance of the plane being base vertices; every other tool edge (a line, arc or ellipse arc) has no
//!    interior extreme of height at or below the tolerance; every other tool face is a plane not on F's plane, or a
//!    cylinder, elliptic cylinder, cone, sphere or torus about an axis square to F (their heights reach their extremes
//!    on their edges);
//! 4. nothing else of the body reaches into the boss's room: every body face other than F has a conservative analytic
//!    box (in F's frame) at or below the plane, above the boss's highest point, or beside the boss's box.
//! Then the union is exact and local: the body is kept entity for entity; F gains the base's outer loop as a new inner
//! loop (the base's own edges, run as the base ran them); each base hole becomes a new face on F's surface; the tool's
//! other faces, edges and vertices are added unchanged; base seams are dropped. Nothing is recomputed. Volume is
//! checked: body plus tool.
//! Refused by name (op 24 then tries the prism stack domain when both are right prisms along one axis): a base touching
//! or crossing F's boundary or one of its holes, a boss reaching below F or touching F's plane outside its base, body
//! geometry within the boss's room, surfaces whose bounds are not proved (general revolution, spline, conoid faces).
use super::*;
use crate::neutral::Face;
use crate::provenance::{Provenance, Record, FROM_BODY_EDGE, FROM_BODY_FACE, FROM_BODY_VERTEX, FROM_TOOL_EDGE, FROM_TOOL_FACE, FROM_TOOL_VERTEX};

/// The admitted union: F, the base faces, the base's outer loop and holes (tool edge occurrences as the base runs them).
pub(in crate::boolean) struct Plan {
    face: usize,
    base: Vec<usize>,
    outer: Vec<(usize, bool)>,
    holes: Vec<Vec<(usize, bool)>>,
}

/// A planar face's anchor, outward normal and in-plane direction.
pub(super) fn plane(face: &Face) -> Option<(Vec3, Vec3, Vec3)> {
    match face.surface {
        Surface::Plane { p, n, x } => Some((p, if face.sense { n } else { n.scale(-1.) }, x)),
        _ => None,
    }
}

/// The range of c0 + a cos t + b sin t over [t0, t1], and its values at interior extremes (t = atan2(b, a) + j pi).
fn harmonic(c0: f64, a: f64, b: f64, t0: f64, t1: f64) -> (f64, f64, Vec<f64>) {
    let at = |t: f64| c0 + a * t.cos() + b * t.sin();
    let (v0, v1) = (at(t0), at(t1));
    let (mut lo, mut hi) = (v0.min(v1), v0.max(v1));
    let mut inner = vec![];
    if a != 0. || b != 0. {
        let base = b.atan2(a);
        let margin = 1e-12 * (t1 - t0).abs().max(1.);
        for j in -6..=6 {
            let t = base + j as f64 * PI;
            if t > t0.min(t1) + margin && t < t0.max(t1) - margin {
                let v = at(t);
                lo = lo.min(v);
                hi = hi.max(v);
                inner.push(v);
            }
        }
    }
    (lo, hi, inner)
}

/// The range of (q - o).d over a curve, and the values at interior extremes (exact for lines, arcs and ellipse arcs;
/// None for other curves, which are bounded by `world_range`).
pub(super) fn exact_range(curve: &Curve, o: Vec3, d: Vec3) -> Option<(f64, f64, Vec<f64>)> {
    match *curve {
        Curve::Line { a, b } => {
            let (x, y) = (a.sub(o).dot(d), b.sub(o).dot(d));
            Some((x.min(y), x.max(y), vec![]))
        }
        Curve::Arc { c, n, x, r, t0, t1 } => {
            let y = n.cross(x);
            Some(harmonic(c.sub(o).dot(d), r * x.dot(d), r * y.dot(d), t0, t1))
        }
        Curve::Conic { c, u, v, t0, t1 } => Some(harmonic(c.sub(o).dot(d), u.dot(d), v.dot(d), t0, t1)),
        _ => None,
    }
}

/// A conservative range of (q - o).d over any curve: exact where `exact_range` is, else its world box's corners.
pub(super) fn world_range(curve: &Curve, o: Vec3, d: Vec3) -> (f64, f64) {
    if let Some((lo, hi, _)) = exact_range(curve, o, d) {
        return (lo, hi);
    }
    let (lo, hi) = super::super::passage::curve_bounds(curve);
    let mut range = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..8 {
        let q = Vec3::new(if i & 1 == 0 { lo[0] } else { hi[0] }, if i & 2 == 0 { lo[1] } else { hi[1] }, if i & 4 == 0 { lo[2] } else { hi[2] });
        let v = q.sub(o).dot(d);
        range = (range.0.min(v), range.1.max(v));
    }
    range
}

pub(super) fn kind(s: &Surface) -> &'static str {
    match s {
        Surface::Plane { .. } => "planar",
        Surface::Cylinder { .. } => "cylindrical",
        Surface::EllipticCylinder { .. } => "elliptic cylindrical",
        Surface::Cone { .. } => "conical",
        Surface::Sphere { .. } => "spherical",
        Surface::Torus { .. } => "toroidal",
        Surface::Revolution { .. } => "general revolved",
        Surface::Spline { .. } => "spline",
        Surface::Conoid { .. } => "conoid",
    }
}

/// A conservative box of a face in the frame (o; x, y, n): per direction, the range of (q - o).d over the face.
pub(super) fn face_box(body: &Body, face: &Face, o: Vec3, frame: [Vec3; 3]) -> Result<[(f64, f64); 3]> {
    let edges: Vec<&Curve> = face.loops.iter().flatten().map(|&(e, _)| &body.edges[e].curve).collect();
    let along = |p: Vec3, a: Vec3| edges.iter().map(|c| world_range(c, p, a)).fold((f64::INFINITY, f64::NEG_INFINITY), |r, x| (r.0.min(x.0), r.1.max(x.1)));
    // An axis segment p + s a, s in the face's axial range, widened by `reach(d)` across the axis.
    let band = |p: Vec3, a: Vec3, reach: &dyn Fn(Vec3) -> f64| -> [(f64, f64); 3] {
        let (s0, s1) = along(p, a);
        frame.map(|d| {
            let (x, y) = (p.sub(o).dot(d) + s0 * a.dot(d), p.sub(o).dot(d) + s1 * a.dot(d));
            (x.min(y) - reach(d), x.max(y) + reach(d))
        })
    };
    let across = |a: Vec3, d: Vec3| (1. - a.dot(d).powi(2)).max(0.).sqrt();
    Ok(match face.surface {
        // A plane region reaches its extremes on its boundary.
        Surface::Plane { .. } => frame.map(|d| along(o, d)),
        Surface::Cylinder { p, a, r, .. } => band(p, a, &|d| r * across(a, d)),
        Surface::EllipticCylinder { p, a, u, v } => band(p, a, &|d| u.dot(d).hypot(v.dot(d))),
        Surface::Cone { p, a, r, k, .. } => {
            let (s0, s1) = along(p, a);
            let reach = (r + k * s0).abs().max((r + k * s1).abs());
            band(p, a, &|d| reach * across(a, d))
        }
        Surface::Sphere { c, r, .. } => frame.map(|d| (c.sub(o).dot(d) - r, c.sub(o).dot(d) + r)),
        Surface::Torus { c, big, r, .. } => frame.map(|d| (c.sub(o).dot(d) - big - r, c.sub(o).dot(d) + big + r)),
        ref s => {
            return Err(unsupported(&format!(
                "union: the body has a {} face; its clearance from the boss cannot be proved yet",
                kind(s)
            )))
        }
    })
}

/// A curve in F's plane frame as a section edge, run as the occurrence runs it.
pub(super) fn seg2(curve: &Curve, forward: bool, o: Vec3, x: Vec3, y: Vec3, n: Vec3) -> Result<Seg> {
    let to2 = |q: Vec3| P::new(q.sub(o).dot(x), q.sub(o).dot(y));
    let seg = match *curve {
        Curve::Line { a, b } => Seg::Line { a: to2(a), b: to2(b) },
        Curve::Arc { c, n: an, r, t0, t1, .. } => {
            if an.cross(n).length() > 1e-9 {
                return Err(unsupported("union: an arc of the face or base does not lie in its plane"));
            }
            let s = curve.start().sub(c);
            Seg::Arc { c: to2(c), r, a0: s.dot(y).atan2(s.dot(x)), sweep: (t1 - t0) * an.dot(n).signum() }
        }
        Curve::Conic { c, u, v, t0, t1 } => Seg::Conic { c: to2(c), u: P::new(u.dot(x), u.dot(y)), v: P::new(v.dot(x), v.dot(y)), a0: t0, sweep: t1 - t0 },
        _ => return Err(unsupported("union: the face or the base has an edge that is not a line, arc or ellipse arc")),
    };
    Ok(if forward { seg } else { seg.reversed() })
}

/// The base: the tool faces on a body face's plane facing it, the body faces they rest on (candidates for F), and the
/// base's boundary loops (tool edge occurrences as the base faces run them). Ok(None): no such tool face.
pub(super) struct Base {
    pub(super) faces: Vec<usize>,
    pub(super) candidates: Vec<usize>,
    pub(super) loops: Vec<Vec<(usize, bool)>>,
}

pub(super) fn base_of(body: &Body, tool: &Body, tol: f64) -> Result<Option<Base>> {
    let mut pairs = vec![];
    for (i, f) in body.faces.iter().enumerate() {
        let Some((p, n, _)) = plane(f) else { continue };
        for (j, t) in tool.faces.iter().enumerate() {
            let Some((q, m, _)) = plane(t) else { continue };
            if n.add(m).length() <= 1e-9 && q.sub(p).dot(n).abs() <= tol {
                pairs.push((i, j));
            }
        }
    }
    if pairs.is_empty() {
        return Ok(None);
    }
    let mut base: Vec<usize> = pairs.iter().map(|p| p.1).collect();
    base.sort();
    base.dedup();
    let (p0, n0, _) = plane(&body.faces[pairs[0].0]).unwrap();
    if base.iter().any(|&j| { let (q, m, _) = plane(&tool.faces[j]).unwrap(); m.add(n0).length() > 1e-9 || q.sub(p0).dot(n0).abs() > tol }) {
        return Err(unsupported("union: the boss rests on more than one plane of the body"));
    }
    // The base's boundary: occurrences not cancelled by another base face.
    let mut count = vec![0usize; tool.edges.len()];
    for &j in &base {
        for &(e, _) in tool.faces[j].loops.iter().flatten() {
            count[e] += 1;
        }
    }
    let boundary: Vec<(usize, bool)> = base.iter().flat_map(|&j| tool.faces[j].loops.iter().flatten().copied()).filter(|&(e, _)| count[e] == 1).collect();
    let ends = |&(e, w): &(usize, bool)| { let edge = &tool.edges[e]; if w { (edge.start, edge.end) } else { (edge.end, edge.start) } };
    let mut starts = std::collections::HashMap::new();
    for (i, o) in boundary.iter().enumerate() {
        if starts.insert(ends(o).0, i).is_some() {
            return Err(unsupported("union: the boss's base meets itself at a point"));
        }
    }
    let mut seen = vec![false; boundary.len()];
    let mut loops: Vec<Vec<(usize, bool)>> = vec![];
    for first in 0..boundary.len() {
        if seen[first] {
            continue;
        }
        let mut lp = vec![];
        let mut i = first;
        loop {
            if std::mem::replace(&mut seen[i], true) {
                return Err(KernelError::InvalidInput("union: the boss's base outline does not close".into()));
            }
            lp.push(boundary[i]);
            i = *starts.get(&ends(&boundary[i]).1).ok_or_else(|| KernelError::InvalidInput("union: the boss's base outline is open".into()))?;
            if i == first {
                break;
            }
        }
        loops.push(lp);
    }
    let mut candidates: Vec<usize> = pairs.iter().map(|p| p.0).collect();
    candidates.sort();
    candidates.dedup();
    Ok(Some(Base { faces: base, candidates, loops }))
}

/// Step 3: the boss lies above the plane (o, n) and meets it only in its base.
pub(super) fn above(tool: &Body, base: &[usize], o: Vec3, n: Vec3, tol: f64) -> Result<()> {
    let mut base_vertices = std::collections::HashSet::new();
    let mut base_edges = std::collections::HashSet::new();
    for &j in base {
        for &(e, _) in tool.faces[j].loops.iter().flatten() {
            base_edges.insert(e);
            base_vertices.insert(tool.edges[e].start);
            base_vertices.insert(tool.edges[e].end);
        }
    }
    let mut in_use = vec![false; tool.vertices.len()];
    for e in &tool.edges {
        in_use[e.start] = true;
        in_use[e.end] = true;
    }
    for (vi, q) in tool.vertices.iter().enumerate() {
        if !in_use[vi] {
            continue; // e.g. a revolve's axis point under a disc
        }
        let h = q.sub(o).dot(n);
        if h < -tol {
            return Err(unsupported("union: the boss reaches below the face it stands on (a boss into the body is not this domain)"));
        }
        if h <= tol && !base_vertices.contains(&vi) {
            return Err(unsupported("union: the boss touches the face's plane outside its base"));
        }
    }
    for (ei, e) in tool.edges.iter().enumerate() {
        if base_edges.contains(&ei) {
            continue;
        }
        let (_, _, inner) = exact_range(&e.curve, o, n).ok_or_else(|| unsupported("union: the boss has an edge whose height above the face is not bounded yet"))?;
        if inner.iter().any(|&h| h <= tol) {
            return Err(unsupported("union: the boss dips to or below the face's plane along an edge"));
        }
    }
    for (j, t) in tool.faces.iter().enumerate() {
        if base.contains(&j) {
            continue;
        }
        let square = |a: Vec3| a.cross(n).length() <= 1e-9;
        let ok = match t.surface {
            Surface::Plane { p, n: m, .. } => !(m.cross(n).length() <= 1e-9 && p.sub(o).dot(n).abs() <= tol),
            Surface::Cylinder { a, .. } | Surface::EllipticCylinder { a, .. } | Surface::Cone { a, .. } | Surface::Sphere { a, .. } | Surface::Torus { a, .. } => square(a),
            _ => false,
        };
        if !ok {
            return Err(unsupported(&format!("union: the boss has a {} face whose height above the face it stands on is not bounded by its edges", kind(&t.surface))));
        }
    }
    Ok(())
}

/// Step 4: no body face other than F reaches into the boss's room. Returns the boss's box in the frame (x, y, n).
pub(super) fn room(body: &Body, tool: &Body, face: usize, o: Vec3, frame: [Vec3; 3], tol: f64) -> Result<[(f64, f64); 3]> {
    let mut reach = [(f64::INFINITY, f64::NEG_INFINITY); 3];
    for t in &tool.faces {
        let b = face_box(tool, t, o, frame)?;
        for d in 0..3 {
            reach[d] = (reach[d].0.min(b[d].0), reach[d].1.max(b[d].1));
        }
    }
    let top = reach[2].1;
    for (gi, g) in body.faces.iter().enumerate() {
        if gi == face {
            continue;
        }
        let b = face_box(body, g, o, frame)?;
        let clear = b[2].1 <= tol || b[2].0 >= top + tol || (0..2).any(|d| b[d].1 < reach[d].0 - tol || b[d].0 > reach[d].1 + tol);
        if !clear {
            return Err(unsupported("union: another part of the body comes within the boss's room above the face (the boss would run into it)"));
        }
    }
    Ok(reach)
}

/// Admission (see the module notes). Ok(None): no tool face lies on a body face's plane facing it (not this domain).
pub(in crate::boolean) fn admit(body: &Body, tool: &Body, tol: f64) -> Result<Option<Plan>> {
    let Some(Base { faces: base, candidates, loops }) = base_of(body, tool, tol)? else { return Ok(None) };
    // 2. Strictly inside one face F on that plane. The frame: F's outward normal n, its x, y = n x x.
    let mut last_reason = None;
    let mut chosen = None;
    for &fi in &candidates {
        let (o, n, x) = plane(&body.faces[fi]).unwrap();
        let y = n.cross(x);
        let segs = |lp: &[(usize, bool)], b: &Body| lp.iter().map(|&(e, w)| seg2(&b.edges[e].curve, w, o, x, y, n)).collect::<Result<Vec<Seg>>>();
        let face_loops: Vec<Vec<Seg>> = body.faces[fi].loops.iter().map(|lp| segs(lp, body)).collect::<Result<_>>()?;
        let base_loops: Vec<Vec<Seg>> = loops.iter().map(|lp| segs(lp, tool)).collect::<Result<_>>()?;
        // The base outline runs clockwise about n (its faces look the other way): the outer has negative area here.
        let areas: Vec<f64> = base_loops.iter().map(|l| loop_area(l)).collect();
        let outers: Vec<usize> = (0..areas.len()).filter(|&i| areas[i] < 0.).collect();
        if outers.len() != 1 {
            return Err(unsupported("union: the boss's base is not one region (one outline with any holes inside it)"));
        }
        let mut touches = false;
        for s in face_loops.iter().flatten() {
            for t in base_loops.iter().flatten() {
                if !meetings(s, t, tol)?.is_empty() {
                    touches = true;
                }
            }
        }
        if touches {
            last_reason = Some("union: the boss's base touches or crosses the edge of the face it stands on, or a hole in it");
            continue;
        }
        let outer = &base_loops[outers[0]];
        if winding(&face_loops, outer[0].start())? != 1 {
            last_reason = Some("union: the boss's base is not on the face it is sketched on");
            continue;
        }
        if face_loops.iter().any(|fl| winding(std::slice::from_ref(outer), fl[0].start()).map_or(true, |w| w != 0)) {
            last_reason = Some("union: a hole or edge of the face lies under the boss's base");
            continue;
        }
        let outer_loop = loops[outers[0]].clone();
        let holes = loops.iter().enumerate().filter(|(i, _)| *i != outers[0]).map(|(_, l)| l.clone()).collect();
        chosen = Some(Plan { face: fi, base: base.clone(), outer: outer_loop, holes });
        break;
    }
    let plan = chosen.ok_or_else(|| unsupported(last_reason.unwrap_or("union: the boss's base is not inside one face of the body")))?;
    let (o, n, x) = plane(&body.faces[plan.face]).unwrap();
    // 3. The boss above the plane, meeting it only in its base; 4. body clearance in F's frame.
    above(tool, &plan.base, o, n, tol)?;
    room(body, tool, plan.face, o, [x, n.cross(x), n], tol)?;
    Ok(Some(plan))
}

/// The admitted union and its history (op 25 roles 30-39: the body kept entity for entity, F modified or split by base
/// holes, the tool's faces other than the base and every tool edge and vertex they use kept unchanged).
pub(in crate::boolean) fn unite(body: &Body, tool: &Body, plan: &Plan, tol: f64) -> Result<(Body, Provenance)> {
    let mut result = body.clone();
    let mut prov = Provenance::default();
    let face_trim = plan.face;
    for i in 0..body.faces.len() {
        prov.faces.push(Record::new(FROM_BODY_FACE, &[i, (i == face_trim) as usize]));
    }
    for i in 0..body.edges.len() {
        prov.edges.push(Record::new(FROM_BODY_EDGE, &[i, 0]));
    }
    for i in 0..body.vertices.len() {
        prov.vertices.push(Record::new(FROM_BODY_VERTEX, &[i]));
    }
    // Tool edges and vertices the kept faces use (the base boundary among them; base seams are not).
    let mut used = vec![false; tool.edges.len()];
    for (j, t) in tool.faces.iter().enumerate() {
        if !plan.base.contains(&j) {
            for &(e, _) in t.loops.iter().flatten() {
                used[e] = true;
            }
        }
    }
    let mut vmap = vec![usize::MAX; tool.vertices.len()];
    for (e, edge) in tool.edges.iter().enumerate() {
        if used[e] {
            for v in [edge.start, edge.end] {
                if vmap[v] == usize::MAX {
                    vmap[v] = usize::MAX - 1;
                }
            }
        }
    }
    for (v, q) in tool.vertices.iter().enumerate() {
        if vmap[v] != usize::MAX {
            vmap[v] = result.vertices.len();
            result.vertices.push(*q);
            prov.vertices.push(Record::new(FROM_TOOL_VERTEX, &[v]));
        }
    }
    let mut emap = vec![usize::MAX; tool.edges.len()];
    for (e, edge) in tool.edges.iter().enumerate() {
        if used[e] {
            emap[e] = result.edges.len();
            result.edges.push(crate::neutral::Edge { curve: edge.curve.clone(), start: vmap[edge.start], end: vmap[edge.end] });
            prov.edges.push(Record::new(FROM_TOOL_EDGE, &[e, 0]));
        }
    }
    let moved = |lp: &[(usize, bool)]| -> Vec<(usize, bool)> { lp.iter().map(|&(e, w)| (emap[e], w)).collect() };
    for (j, t) in tool.faces.iter().enumerate() {
        if plan.base.contains(&j) {
            continue;
        }
        result.faces.push(Face { surface: t.surface.clone(), sense: t.sense, loops: t.loops.iter().map(|l| moved(l)).collect() });
        prov.faces.push(Record::new(FROM_TOOL_FACE, &[j, 0]));
    }
    if plan.outer.iter().chain(plan.holes.iter().flatten()).any(|&(e, _)| !used[e]) {
        return Err(KernelError::Internal("union: a base outline edge bounds no other boss face".into()));
    }
    result.faces[plan.face].loops.push(moved(&plan.outer));
    let (surface, sense) = (body.faces[plan.face].surface.clone(), body.faces[plan.face].sense);
    for hole in &plan.holes {
        result.faces.push(Face { surface: surface.clone(), sense, loops: vec![moved(hole)] });
        prov.faces.push(Record::new(FROM_BODY_FACE, &[plan.face, 1]));
    }
    crate::topology::check(&result, tol).map_err(|e| KernelError::Internal(format!("union built an invalid solid: {e}")))?;
    let (mb, mt, mr) = (crate::mass::measure(body, tol)?, crate::mass::measure(tool, tol)?, crate::mass::measure(&result, tol)?);
    let expected = mb.volume + mt.volume;
    if (mr.volume - expected).abs() > tol * mr.area + 1e-9 * expected.abs().max(1.) {
        return Err(KernelError::Internal(format!("union built {}, body and boss hold {expected}", mr.volume)));
    }
    Ok((result, prov))
}

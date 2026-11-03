//! Boundary tidying for unions (Claude 2026-11-02). Exact topology edits only: no point is moved and no surface is
//! fitted.
//! - `join_faces`: faces of one group (the caller's proof that they lie on one support with one outward side) are joined
//!   across every edge they share. Two faces become one by splicing their loops at the shared edge; an edge met twice in
//!   one face is a spike (removed), a split of one loop into two (removed: a hole appears), or, on a periodic support,
//!   the seam of a loop pair winding round the axis (kept, so every cylinder band keeps one seam as the extruder makes).
//! - `join_edges`: two edges on one support meeting at a vertex that only they use, and bounding the same two faces, are
//!   one edge (collinear lines; arcs of one circle; arcs of one ellipse).
//! - `compact`: drops unused edges and vertices and renumbers.
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;
use std::f64::consts::TAU;

type Result<T> = std::result::Result<T, KernelError>;

fn locate(face: &Face, edge: usize) -> Vec<(usize, usize)> {
    let mut out = vec![];
    for (l, lp) in face.loops.iter().enumerate() {
        for (i, &(e, _)) in lp.iter().enumerate() {
            if e == edge {
                out.push((l, i));
            }
        }
    }
    out
}

/// Joins the faces of each group across their shared edges (see the module notes). `wraps(edges, loop)` says whether
/// a loop winds round the periodic support of the face's group (only asked for a loop that would split in two).
pub(super) fn join_faces(body: &mut Body, group: &[Option<usize>], wraps: &dyn Fn(&[Edge], &[(usize, bool)]) -> bool) -> Result<()> {
    let n = body.faces.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut f: usize) -> usize {
        while p[f] != f {
            p[f] = p[p[f]];
            f = p[f];
        }
        f
    }
    let mut uses = vec![vec![]; body.edges.len()];
    for (f, face) in body.faces.iter().enumerate() {
        for &(e, w) in face.loops.iter().flatten() {
            uses[e].push((f, w));
        }
    }
    let mut faces: Vec<Option<Face>> = std::mem::take(&mut body.faces).into_iter().map(Some).collect();
    for (e, used) in uses.iter().enumerate() {
        let [(f1, w1), (f2, w2)] = used[..] else { continue };
        if w1 == w2 || group[f1].is_none() || group[f1] != group[f2] {
            continue;
        }
        let (r1, r2) = (find(&mut parent, f1), find(&mut parent, f2));
        if r1 != r2 {
            let b = faces[r2].take().ok_or_else(|| KernelError::Internal("union join: a face was lost".into()))?;
            let a = faces[r1].as_mut().ok_or_else(|| KernelError::Internal("union join: a face was lost".into()))?;
            let (Some(&(la, ia)), Some(&(lb, ib))) = (locate(a, e).first(), locate(&b, e).first()) else {
                return Err(KernelError::Internal("union join: a shared edge is missing from its face".into()));
            };
            let (lpa, lpb) = (a.loops[la].clone(), &b.loops[lb]);
            let mut spliced = lpa[..ia].to_vec();
            spliced.extend_from_slice(&lpb[ib + 1..]);
            spliced.extend_from_slice(&lpb[..ib]);
            spliced.extend_from_slice(&lpa[ia + 1..]);
            a.loops[la] = spliced;
            for (l, lp) in b.loops.iter().enumerate() {
                if l != lb {
                    a.loops.push(lp.clone());
                }
            }
            parent[r2] = r1;
            continue;
        }
        let a = faces[r1].as_mut().ok_or_else(|| KernelError::Internal("union join: a face was lost".into()))?;
        let found = locate(a, e);
        let [(l1, i1), (l2, i2)] = found[..] else {
            return Err(KernelError::Internal("union join: an edge shared within one face is not met twice".into()));
        };
        if l1 == l2 {
            let lp = a.loops[l1].clone();
            let m = lp.len();
            let (i, j) = (i1.min(i2), i1.max(i2));
            if j == i + 1 || (i == 0 && j == m - 1) {
                // A spike: the edge out and straight back.
                let kept: Vec<_> = lp.iter().enumerate().filter(|(x, _)| *x != i && *x != j).map(|(_, o)| *o).collect();
                if kept.is_empty() {
                    a.loops.remove(l1);
                } else {
                    a.loops[l1] = kept;
                }
                continue;
            }
            let inner = lp[i + 1..j].to_vec();
            let mut outer = lp[j + 1..].to_vec();
            outer.extend_from_slice(&lp[..i]);
            if wraps(&body.edges, &inner) && wraps(&body.edges, &outer) {
                continue; // the seam of a band round a periodic support
            }
            a.loops[l1] = outer;
            a.loops.push(inner);
        } else {
            // Two loops of the face joined through the edge.
            let (lpa, lpb) = (a.loops[l1].clone(), a.loops[l2].clone());
            let mut joined = lpa[..i1].to_vec();
            joined.extend_from_slice(&lpb[i2 + 1..]);
            joined.extend_from_slice(&lpb[..i2]);
            joined.extend_from_slice(&lpa[i1 + 1..]);
            a.loops[l1] = joined;
            a.loops.remove(l2);
        }
    }
    body.faces = faces.into_iter().flatten().collect();
    Ok(())
}

/// The curve run from `at` (one of its ends' vertices): itself or reversed.
fn from_vertex(edge: &Edge, at: usize) -> Curve {
    if edge.start == at { edge.curve.clone() } else { edge.curve.reversed() }
}

/// c1 then c2 (c1 ends where c2 starts) as one curve of their common support, when they have one.
fn join_curves(c1: &Curve, c2: &Curve, tol: f64) -> Option<Curve> {
    match (c1, c2) {
        (Curve::Line { a, b }, Curve::Line { a: _, b: d }) => {
            let (u, w) = (b.sub(*a), d.sub(*a));
            let l = w.length();
            if !(l > 0.) || u.dot(w) <= 0. || d.sub(*b).dot(w) <= 0. {
                return None;
            }
            // The joint lies on the joined line within the tolerance.
            (u.cross(w).length() / l <= tol).then_some(Curve::Line { a: *a, b: *d })
        }
        (Curve::Arc { c, n, x, r, t0, t1 }, Curve::Arc { c: c2, n: n2, r: r2, t0: s0, t1: s1, .. }) => {
            if c.distance(*c2) > tol || (r - r2).abs() > tol || n.dot(*n2) < 1. - 1e-12 {
                return None;
            }
            let total = (t1 - t0) + (s1 - s0);
            (total <= TAU + 1e-12).then(|| Curve::Arc { c: *c, n: *n, x: *x, r: *r, t0: *t0, t1: t0 + total.min(TAU) })
        }
        (Curve::Conic { c, u, v, t0, t1 }, Curve::Conic { c: c2, u: u2, v: v2, t0: s0, t1: s1 }) => {
            if c.distance(*c2) > tol || u.distance(*u2) > tol || v.distance(*v2) > tol {
                return None;
            }
            // The second continues the first's parameter (mod a turn).
            let gap = (s0 - t1).rem_euclid(TAU);
            if gap.min(TAU - gap) > 1e-9 {
                return None;
            }
            let total = (t1 - t0) + (s1 - s0);
            (total <= TAU + 1e-12).then(|| Curve::Conic { c: *c, u: *u, v: *v, t0: *t0, t1: t0 + total.min(TAU) })
        }
        _ => None,
    }
}

/// Joins edges through vertices only they use (see the module notes), until none is left.
pub(super) fn join_edges(body: &mut Body, tol: f64) {
    loop {
        let mut incident = vec![vec![]; body.vertices.len()];
        for (i, e) in body.edges.iter().enumerate() {
            incident[e.start].push(i);
            incident[e.end].push(i);
        }
        let mut faces_of = vec![vec![]; body.edges.len()];
        for (f, face) in body.faces.iter().enumerate() {
            for &(e, _) in face.loops.iter().flatten() {
                faces_of[e].push(f);
            }
        }
        for list in &mut faces_of {
            list.sort();
        }
        let mut changed = false;
        for v in 0..body.vertices.len() {
            let [e1, e2] = incident[v][..] else { continue };
            if e1 == e2 || faces_of[e1] != faces_of[e2] || faces_of[e1].len() != 2 || faces_of[e1][0] == faces_of[e1][1] {
                continue;
            }
            let (a, b) = (&body.edges[e1], &body.edges[e2]);
            let (far1, far2) = (if a.start == v { a.end } else { a.start }, if b.start == v { b.end } else { b.start });
            let into = from_vertex(a, v).reversed();
            let out = from_vertex(b, v);
            let Some(joined) = join_curves(&into, &out, tol) else { continue };
            // A closed line is no edge; a closed arc keeps its one vertex.
            if far1 == far2 && matches!(joined, Curve::Line { .. }) {
                continue;
            }
            // The joined edge runs far1 -> far2, as the first edge was run into v.
            let (old_start, old_end) = (body.edges[e1].start, body.edges[e1].end);
            body.edges[e1] = Edge { curve: joined, start: far1, end: far2 };
            for face in &mut body.faces {
                for lp in &mut face.loops {
                    let Some(i) = lp.iter().position(|o| o.0 == e1) else { continue };
                    let m = lp.len();
                    let (next, prev) = ((i + 1) % m, (i + m - 1) % m);
                    // This face runs e1 into v then e2 (the joined edge's way), or e2 into v then e1 out (against it).
                    let into_v = (if lp[i].1 { old_end } else { old_start }) == v;
                    if into_v && lp[next].0 == e2 {
                        lp[i].1 = true;
                        lp.remove(next);
                    } else if !into_v && lp[prev].0 == e2 {
                        lp[i].1 = false;
                        lp.remove(prev);
                    }
                }
            }
            // e2 is now unused: a closed placeholder at the far vertex, dropped by `compact`.
            body.edges[e2] = Edge { curve: Curve::Line { a: body.vertices[far2], b: body.vertices[far2] }, start: usize::MAX, end: usize::MAX };
            changed = true;
            break;
        }
        if !changed {
            break;
        }
        // Drop the placeholder before the next pass reads the incidences.
        let _ = compact(body);
    }
}

/// Drops edges no face uses and vertices no edge uses, and renumbers. Returns the old -> new edge and vertex maps
/// (usize::MAX: dropped).
pub(super) fn compact(body: &mut Body) -> (Vec<usize>, Vec<usize>) {
    let mut used_e = vec![false; body.edges.len()];
    for &(e, _) in body.faces.iter().flat_map(|f| f.loops.iter().flatten()) {
        used_e[e] = true;
    }
    let mut used_v = vec![false; body.vertices.len()];
    for (i, e) in body.edges.iter().enumerate() {
        if used_e[i] {
            used_v[e.start] = true;
            used_v[e.end] = true;
        }
    }
    let mut vm = vec![usize::MAX; body.vertices.len()];
    let mut vertices = vec![];
    for (i, q) in body.vertices.iter().enumerate() {
        if used_v[i] {
            vm[i] = vertices.len();
            vertices.push(*q);
        }
    }
    let mut em = vec![usize::MAX; body.edges.len()];
    let mut edges = vec![];
    for (i, e) in body.edges.iter().enumerate() {
        if used_e[i] {
            em[i] = edges.len();
            edges.push(Edge { curve: e.curve.clone(), start: vm[e.start], end: vm[e.end] });
        }
    }
    for face in &mut body.faces {
        for o in face.loops.iter_mut().flatten() {
            o.0 = em[o.0];
        }
    }
    body.vertices = vertices;
    body.edges = edges;
    (em, vm)
}

/// Twice the vector area of a closed planar loop of lines, arcs and ellipse arcs (exact: the integral of q x dq).
pub(super) fn vector_area(edges: &[Edge], lp: &[(usize, bool)]) -> Option<Vec3> {
    let mut sum = Vec3::new(0., 0., 0.);
    for &(e, forward) in lp {
        let curve = if forward { edges[e].curve.clone() } else { edges[e].curve.reversed() };
        let term = match curve {
            Curve::Line { a, b } => a.cross(b),
            Curve::Arc { c, n, r, t0, t1, .. } => c.cross(curve.end().sub(curve.start())).add(n.scale(r * r * (t1 - t0))),
            Curve::Conic { c, u, v, t0, t1 } => c.cross(curve.end().sub(curve.start())).add(u.cross(v).scale(t1 - t0)),
            _ => return None,
        };
        sum = sum.add(term);
    }
    Some(sum)
}

/// Puts each planar face's outer loop (the largest area counter-clockwise about its outward normal) first.
pub(super) fn outer_first(body: &mut Body) {
    for face in &mut body.faces {
        let Surface::Plane { n, .. } = face.surface else { continue };
        if face.loops.len() < 2 {
            continue;
        }
        let out = if face.sense { n } else { n.scale(-1.) };
        let areas: Vec<f64> = face.loops.iter().map(|lp| vector_area(&body.edges, lp).map_or(f64::NEG_INFINITY, |a| a.dot(out))).collect();
        if let Some(best) = (0..areas.len()).max_by(|&i, &j| areas[i].total_cmp(&areas[j])) {
            face.loops.swap(0, best);
        }
    }
}

/// Whether the faces form one connected shell through their shared edges.
pub(super) fn connected(body: &Body) -> bool {
    let n = body.faces.len();
    if n == 0 {
        return false;
    }
    let mut owners = vec![vec![]; body.edges.len()];
    for (f, face) in body.faces.iter().enumerate() {
        for &(e, _) in face.loops.iter().flatten() {
            owners[e].push(f);
        }
    }
    let mut seen = vec![false; n];
    let mut stack = vec![0];
    seen[0] = true;
    let mut adjacency = vec![vec![]; n];
    for list in &owners {
        for &a in list {
            for &b in list {
                if a != b {
                    adjacency[a].push(b);
                }
            }
        }
    }
    while let Some(f) = stack.pop() {
        for &g in &adjacency[f] {
            if !seen[g] {
                seen[g] = true;
                stack.push(g);
            }
        }
    }
    seen.iter().all(|s| *s)
}

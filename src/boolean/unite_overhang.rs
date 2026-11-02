//! Union domain: a boss on one planar face F of a body that is not a layered prism stack, its base crossing F's edge
//! (Claude 2026-11-02, op 24): a round boss over an I-beam flange's tip, a pad over an angle leg's toe.
//!
//! Admission, decided from the geometry before anything is built, after the face and stack domains declined:
//! 1. the base (as the face domain: tool faces on the plane of a body face, facing it) touches exactly one body face F
//!    on that plane, and its boundary crosses F's boundary transversally: in F's plane (frame x, n x x) both outlines
//!    are split at their closed-form meetings (`section::split`, points within the tolerance one vertex); a piece of
//!    the base running along a piece of F (flush) is refused; F's pieces are classified inside or outside the base and
//!    the base's inside or outside F by the winding of their midpoints. The base must rest on F (some of it inside F)
//!    and hang over (some of it outside);
//! 2. the boss lies above F's plane and meets it only in its base, and nothing of the body reaches into its room (the
//!    face domain's checks);
//! 3. under the overhang there is nothing of the body at the plane: no other body face on the plane near the base;
//!    every other body face reaching the plane near the base a plane or a cylinder whose edges reach the plane only at
//!    F's own vertices (lines' ends; arcs' and ellipse arcs' ends and interior extremes, in closed form);
//!    no base corner on F's outline (a vertex of both bodies).
//! Then the union is local and exact: F becomes F less the base (one face per region); the overhang (the base less F)
//! becomes the boss's underside, on the base's own plane; F's edges and the base's edges that are crossed are split at
//! the crossings (pieces lifted from the arrangement onto F's plane: lines between their vertices, arcs and ellipse
//! arcs on their own supports), the faces beside them (F's neighbours, the boss's walls) take the pieces; every other
//! entity of both bodies is kept unchanged. Volume checked: body plus tool.
//! History (op 25, roles 30-42): every kept entity its source [i, 0]; F's regions [F, 1]; faces beside a split edge
//! [i, 1]; the underside [base face, 1] (one base face only, else op 25 refuses); split edges' pieces [e, 1]; each
//! crossing EDGE_CROSSING [body edge, tool edge, n] (n along the tool edge).
use super::unite_face::{above, base_of, face_box, plane, room, seg2, Base};
use super::unite_stack::chain;
use super::*;
use crate::boolean::unite_tidy;
use crate::neutral::{Edge, Face};
use crate::provenance::{Provenance, Record, EDGE_CROSSING, FROM_BODY_EDGE, FROM_BODY_FACE, FROM_BODY_VERTEX, FROM_TOOL_EDGE, FROM_TOOL_FACE, FROM_TOOL_VERTEX};

/// The admitted overhang: F, the base faces, F's frame, the arrangement and its classification.
pub(in crate::boolean) struct Plan {
    face: usize,
    base: Vec<usize>,
    o: Vec3,
    x: Vec3,
    y: Vec3,
    n: Vec3,
    pool: Vec<P>,
    /// F's occurrences (edge, way) in F's loop order, flattened; the base's as the base runs them reversed (B2: the
    /// base region on the left in F's frame), flattened.
    f_occ: Vec<(usize, bool)>,
    b_occ: Vec<(usize, bool)>,
    /// Pieces: (occurrence, piece, inside the other region).
    f_pieces: Vec<(usize, Piece, bool)>,
    b_pieces: Vec<(usize, Piece, bool)>,
    /// Pool point -> (tool?, vertex) of the bodies, when it is one of their vertices.
    known: Vec<Option<(bool, usize)>>,
}

/// A piece of the arrangement as a 3D curve on F's plane, run the piece's way.
fn lift(seg: &Seg, a: Vec3, b: Vec3, o: Vec3, x: Vec3, y: Vec3, n: Vec3) -> Curve {
    let at = |q: P| o.add(x.scale(q.x)).add(y.scale(q.y));
    let dir = |q: P| x.scale(q.x).add(y.scale(q.y));
    match *seg {
        Seg::Line { .. } => Curve::Line { a, b },
        Seg::Arc { c, r, a0, sweep } => {
            let start = dir(P::new(a0.cos(), a0.sin()));
            if sweep > 0. {
                Curve::Arc { c: at(c), n, x: start, r, t0: 0., t1: sweep }
            } else {
                Curve::Arc { c: at(c), n: n.scale(-1.), x: start, r, t0: 0., t1: -sweep }
            }
        }
        Seg::Conic { c, u, v, a0, sweep } => {
            if sweep > 0. {
                Curve::Conic { c: at(c), u: dir(u), v: dir(v), t0: a0, t1: a0 + sweep }
            } else {
                Curve::Conic { c: at(c), u: dir(u), v: dir(v).scale(-1.), t0: -a0, t1: -a0 - sweep }
            }
        }
    }
}

/// Groups chained loops (region on the left) into faces: each outer loop (positive area) with the holes it is the
/// innermost outline of.
fn regions(bits: &[(usize, usize, Seg)], vertices: usize) -> Result<Vec<Vec<Vec<usize>>>> {
    let loops = chain(bits, vertices)?;
    let area = |lp: &Vec<usize>| lp.iter().map(|&x| bits[x].2.area()).sum::<f64>();
    let (outers, holes): (Vec<&Vec<usize>>, Vec<&Vec<usize>>) = loops.iter().partition(|lp| area(lp) > 0.);
    let mut owned: Vec<Vec<Vec<usize>>> = outers.iter().map(|o| vec![(*o).clone()]).collect();
    for h in holes {
        let probe = bits[h[0]].2.at(0.5);
        let mut best: Option<(f64, usize)> = None;
        for (o, outer) in outers.iter().enumerate() {
            let segs: Vec<Seg> = outer.iter().map(|&x| bits[x].2).collect();
            if winding(&[segs], probe)? != 0 && best.is_none_or(|(a, _)| area(outer) < a) {
                best = Some((area(outer), o));
            }
        }
        let (_, o) = best.ok_or_else(|| KernelError::Internal("overhang union: a hole lies in no outline".into()))?;
        owned[o].push(h.clone());
    }
    Ok(owned)
}

/// Admission (see the module notes). Ok(None): not an overhang of a base on a face (no base, or no crossing).
pub(in crate::boolean) fn admit(body: &Body, tool: &Body, tol: f64) -> Result<Option<Plan>> {
    let Some(Base { faces: base, candidates, loops }) = base_of(body, tool, tol)? else { return Ok(None) };
    let mut chosen = None;
    for &fi in &candidates {
        let (o, n, x) = plane(&body.faces[fi]).unwrap();
        let y = n.cross(x);
        let f_occ: Vec<(usize, bool)> = body.faces[fi].loops.iter().flatten().copied().collect();
        let b_occ: Vec<(usize, bool)> = loops.iter().flat_map(|lp| lp.iter().rev().map(|&(e, w)| (e, !w))).collect();
        let f_loops: Vec<Vec<Seg>> = body.faces[fi].loops.iter().map(|lp| lp.iter().map(|&(e, w)| seg2(&body.edges[e].curve, w, o, x, y, n)).collect::<Result<_>>()).collect::<Result<_>>()?;
        let b_loops: Vec<Vec<Seg>> = loops.iter().map(|lp| lp.iter().rev().map(|&(e, w)| seg2(&tool.edges[e].curve, !w, o, x, y, n)).collect::<Result<_>>()).collect::<Result<_>>()?;
        let flat_tags = |ls: &[Vec<Seg>], tool: bool| {
            let mut k = 0;
            ls.iter().map(|l| l.iter().map(|_| { k += 1; Tag { tool, seg: k - 1, s0: 0., s1: 1. } }).collect::<Vec<_>>()).collect::<Vec<_>>()
        };
        let (f_tags, b_tags) = (flat_tags(&f_loops, false), flat_tags(&b_loops, true));
        let mut pool = vec![];
        let pf = split(&f_loops, &f_tags, &b_loops, &mut pool, tol)?;
        let pb = split(&b_loops, &b_tags, &f_loops, &mut pool, tol)?;
        if pf.len() + pb.len() > 4 * LIMIT {
            return Err(unsupported("overhang union exceeds its work limit"));
        }
        for p in &pf {
            let mid = p.seg.at(0.5);
            if pb.iter().any(|q| ((q.v0 == p.v0 && q.v1 == p.v1) || (q.v0 == p.v1 && q.v1 == p.v0)) && q.seg.at(0.5).dist(mid) <= tol) {
                return Err(unsupported("union: the boss's edge runs along the face's edge (flush) on a body that is not a layered prism stack; not supported yet"));
            }
        }
        let mut f_pieces = vec![];
        for p in pf {
            let mid = p.seg.at(0.5);
            if b_loops.iter().flatten().any(|s| s.distance(mid) <= tol) {
                return Err(unsupported("overhang union: outlines meet along an unmatched piece"));
            }
            let inside = winding(&b_loops, mid)? != 0;
            f_pieces.push((p.tag.seg, p, inside));
        }
        let mut b_pieces = vec![];
        for p in pb {
            let mid = p.seg.at(0.5);
            if f_loops.iter().flatten().any(|s| s.distance(mid) <= tol) {
                return Err(unsupported("overhang union: outlines meet along an unmatched piece"));
            }
            let inside = winding(&f_loops, mid)? != 0;
            b_pieces.push((p.tag.seg, p, inside));
        }
        let rests = b_pieces.iter().any(|b| b.2) || f_pieces.iter().any(|f| f.2);
        let hangs = b_pieces.iter().any(|b| !b.2);
        let crosses = f_pieces.iter().any(|f| f.2);
        if !rests {
            continue;
        }
        if chosen.is_some() {
            return Err(unsupported("union: the boss rests on more than one face of the body"));
        }
        if !(hangs && crosses) {
            return Ok(None); // inside F (the face domain's) or not crossing: not an overhang
        }
        // The bodies' vertices among the pool points.
        let mut known: Vec<Option<(bool, usize)>> = vec![None; pool.len()];
        let mut mark = |q: Vec3, owner: (bool, usize)| -> Result<()> {
            let p2 = P::new(q.sub(o).dot(x), q.sub(o).dot(y));
            if let Some(i) = pool.iter().position(|v| v.dist(p2) <= tol) {
                match known[i] {
                    Some(k) if k != owner => return Err(unsupported("union: a corner of the boss's base lies on the face's outline")),
                    _ => known[i] = Some(owner),
                }
            }
            Ok(())
        };
        for &(e, _) in &f_occ {
            for v in [body.edges[e].start, body.edges[e].end] {
                mark(body.vertices[v], (false, v))?;
            }
        }
        for &(e, _) in &b_occ {
            for v in [tool.edges[e].start, tool.edges[e].end] {
                mark(tool.vertices[v], (true, v))?;
            }
        }
        chosen = Some(Plan { face: fi, base: base.clone(), o, x, y, n, pool, f_occ, b_occ, f_pieces, b_pieces, known });
    }
    let Some(plan) = chosen else { return Ok(None) };
    let (o, x, y, n) = (plan.o, plan.x, plan.y, plan.n);
    above(tool, &plan.base, o, n, tol)?;
    let reach = room(body, tool, plan.face, o, [x, y, n], tol)?;
    // 3. Nothing of the body at the plane under the overhang: near the boss (its box in the plane, tolerance widened),
    // no other face on the plane, and faces reaching it only at F's vertices.
    let near = |b: &[(f64, f64); 3]| (0..2).all(|d| b[d].1 >= reach[d].0 - tol && b[d].0 <= reach[d].1 + tol);
    let f_edges: std::collections::HashSet<usize> = plan.f_occ.iter().map(|o| o.0).collect();
    let f_vertices: std::collections::HashSet<usize> = plan.f_occ.iter().flat_map(|&(e, _)| [body.edges[e].start, body.edges[e].end]).collect();
    for (gi, g) in body.faces.iter().enumerate() {
        if gi == plan.face {
            continue;
        }
        let b = face_box(body, g, o, [x, y, n])?;
        if b[2].1 < -tol || !near(&b) {
            continue;
        }
        if !matches!(g.surface, Surface::Plane { .. } | Surface::Cylinder { .. }) {
            return Err(unsupported("union: a curved part of the body reaches the plane under the boss's overhang"));
        }
        if b[2].0 >= -tol {
            return Err(unsupported("union: the boss rests on more than one face of the body"));
        }
        for &(e, _) in g.loops.iter().flatten() {
            if f_edges.contains(&e) {
                continue;
            }
            let edge = &body.edges[e];
            let (_, _, inner) = super::unite_face::exact_range(&edge.curve, o, n)
                .ok_or_else(|| unsupported("union: an edge near the boss's overhang is not bounded yet"))?;
            let touches_inside = inner.iter().any(|&h| h >= -tol);
            let touches_end = [edge.start, edge.end].iter().any(|&v| body.vertices[v].sub(o).dot(n) >= -tol && !f_vertices.contains(&v));
            if touches_inside || touches_end {
                return Err(unsupported("union: another part of the body reaches the plane under the boss's overhang"));
            }
        }
    }
    Ok(Some(plan))
}

/// The admitted overhang union and its history (see the module notes).
pub(in crate::boolean) fn unite(body: &Body, tool: &Body, plan: &Plan, tol: f64) -> Result<(Body, Provenance)> {
    let (o, x, y, n) = (plan.o, plan.x, plan.y, plan.n);
    let mut r = body.clone();
    let (tv, te) = (r.vertices.len(), r.edges.len());
    r.vertices.extend(tool.vertices.iter().copied());
    r.edges.extend(tool.edges.iter().map(|e| Edge { curve: e.curve.clone(), start: e.start + tv, end: e.end + tv }));
    // Pool point -> result vertex.
    let mut vrec: Vec<Record> = (0..body.vertices.len()).map(|i| Record::new(FROM_BODY_VERTEX, &[i]))
        .chain((0..tool.vertices.len()).map(|j| Record::new(FROM_TOOL_VERTEX, &[j]))).collect();
    let mut erec: Vec<Record> = (0..body.edges.len()).map(|i| Record::new(FROM_BODY_EDGE, &[i, 0]))
        .chain((0..tool.edges.len()).map(|j| Record::new(FROM_TOOL_EDGE, &[j, 0]))).collect();
    // A crossing's body edge and tool edge, and its parameter along the tool edge (B2 runs it backwards).
    let crossing = |i: usize| -> Option<(usize, usize, f64)> {
        let f = plan.f_pieces.iter().find(|(_, p, _)| p.v0 == i || p.v1 == i)?;
        let (k, p, _) = plan.b_pieces.iter().find(|(_, p, _)| p.v0 == i || p.v1 == i)?;
        let s = if p.v0 == i { p.tag.s0 } else { p.tag.s1 };
        let (t, w) = plan.b_occ[*k];
        // The parameter along the tool edge's own direction: B2's s runs the occurrence (t, w).
        let along = if w { s } else { 1. - s };
        Some((plan.f_occ[f.0].0, t, along))
    };
    let crossings: Vec<Option<(usize, usize, f64)>> = (0..plan.pool.len()).map(|i| if plan.known[i].is_none() { crossing(i) } else { None }).collect();
    let mut vmap = vec![usize::MAX; plan.pool.len()];
    for (i, k) in plan.known.iter().enumerate() {
        vmap[i] = match *k {
            Some((false, v)) => v,
            Some((true, v)) => tv + v,
            None => {
                let (e, t, s) = crossings[i].ok_or_else(|| KernelError::Internal("overhang union: a new vertex is no crossing".into()))?;
                let n = crossings.iter().flatten().filter(|c| c.0 == e && c.1 == t && c.2 < s).count();
                vrec.push(Record::new(EDGE_CROSSING, &[e, t, n]));
                r.vertices.push(o.add(x.scale(plan.pool[i].x)).add(y.scale(plan.pool[i].y)));
                r.vertices.len() - 1
            }
        };
    }
    // Each occurrence's pieces as result occurrences run the occurrence's way: kept whole, or new edges per piece.
    let pieces_of = |list: &[(usize, Piece, bool)], k: usize| -> Vec<usize> {
        let mut ids: Vec<usize> = (0..list.len()).filter(|&i| list[i].0 == k).collect();
        ids.sort_by(|&a, &b| list[a].1.tag.s0.min(list[a].1.tag.s1).total_cmp(&list[b].1.tag.s0.min(list[b].1.tag.s1)));
        ids
    };
    let mut f_piece_occ = vec![(0, true); plan.f_pieces.len()];
    let mut b_piece_occ = vec![(0, true); plan.b_pieces.len()];
    // Replacements in the faces beside: (edge, way) -> the occurrences that take its place, in its order.
    let mut replace: std::collections::HashMap<(usize, bool), Vec<(usize, bool)>> = std::collections::HashMap::new();
    for (side, occs, list, out, offset) in [
        (false, &plan.f_occ, &plan.f_pieces, &mut f_piece_occ, 0usize),
        (true, &plan.b_occ, &plan.b_pieces, &mut b_piece_occ, te),
    ] {
        for (k, &(e, w)) in occs.iter().enumerate() {
            let ids = pieces_of(list, k);
            let whole = ids.len() == 1 && list[ids[0]].1.tag.whole();
            if whole {
                out[ids[0]] = (e + offset, w);
                continue;
            }
            let mut along = vec![];
            for &i in &ids {
                let p = &list[i].1;
                let (s, t) = (vmap[p.v0], vmap[p.v1]);
                let curve = lift(&p.seg, r.vertices[s], r.vertices[t], o, x, y, n);
                r.edges.push(Edge { curve, start: s, end: t });
                erec.push(if side { Record::new(FROM_TOOL_EDGE, &[e, 1]) } else { Record::new(FROM_BODY_EDGE, &[e, 1]) });
                out[i] = (r.edges.len() - 1, true);
                along.push((r.edges.len() - 1, true));
            }
            // F's neighbour runs F's edge the other way; the boss's wall runs a base edge as B2 does (the base's
            // opposite).
            if side {
                replace.insert((e + offset, w), along);
            } else {
                replace.insert((e, !w), along.iter().rev().map(|&(id, f)| (id, !f)).collect());
            }
        }
    }
    let swap = |lp: &[(usize, bool)]| -> Vec<(usize, bool)> {
        lp.iter().flat_map(|o| replace.get(o).cloned().unwrap_or_else(|| vec![*o])).collect()
    };
    let mut faces: Vec<Face> = vec![];
    // Regions: F less the base (F's pieces outside, base pieces inside F run backwards), and the overhang (base pieces
    // outside F, F's pieces inside the base run backwards), both chained with the region on the left.
    let regions_of = |keep_f: bool| -> Result<Vec<Vec<Vec<(usize, bool)>>>> {
        let mut bits = vec![];
        let mut occ = vec![];
        for (i, (_, p, inside)) in plan.f_pieces.iter().enumerate() {
            if *inside != keep_f {
                let o = f_piece_occ[i];
                if keep_f {
                    bits.push((p.v0, p.v1, p.seg));
                    occ.push(o);
                } else {
                    bits.push((p.v1, p.v0, p.seg.reversed()));
                    occ.push((o.0, !o.1));
                }
            }
        }
        for (i, (_, p, inside)) in plan.b_pieces.iter().enumerate() {
            if *inside == keep_f {
                let o = b_piece_occ[i];
                if keep_f {
                    bits.push((p.v1, p.v0, p.seg.reversed()));
                    occ.push((o.0, !o.1));
                } else {
                    bits.push((p.v0, p.v1, p.seg));
                    occ.push(o);
                }
            }
        }
        Ok(regions(&bits, plan.pool.len())?.into_iter().map(|face| face.into_iter().map(|lp| lp.into_iter().map(|i| occ[i]).collect()).collect()).collect())
    };
    let kept = regions_of(true)?;
    let hanging = regions_of(false)?;
    if kept.is_empty() || hanging.is_empty() {
        return Err(KernelError::Internal("overhang union: a region vanished".into()));
    }
    let mut frec = vec![];
    let changed = |lp: &[(usize, bool)]| lp.iter().any(|o| replace.contains_key(o));
    let (fs, fsense) = (body.faces[plan.face].surface.clone(), body.faces[plan.face].sense);
    for (i, g) in body.faces.iter().enumerate() {
        if i == plan.face {
            faces.push(Face { surface: fs.clone(), sense: fsense, loops: kept[0].clone() });
            frec.push(Record::new(FROM_BODY_FACE, &[i, 1]));
        } else {
            faces.push(Face { surface: g.surface.clone(), sense: g.sense, loops: g.loops.iter().map(|l| swap(l)).collect() });
            frec.push(Record::new(FROM_BODY_FACE, &[i, g.loops.iter().any(|l| changed(l)) as usize]));
        }
    }
    for extra in &kept[1..] {
        faces.push(Face { surface: fs.clone(), sense: fsense, loops: extra.clone() });
        frec.push(Record::new(FROM_BODY_FACE, &[plan.face, 1]));
    }
    for (j, t) in tool.faces.iter().enumerate() {
        if plan.base.contains(&j) {
            continue;
        }
        let moved: Vec<Vec<(usize, bool)>> = t.loops.iter().map(|l| l.iter().map(|&(e, w)| (e + te, w)).collect::<Vec<_>>()).collect();
        frec.push(Record::new(FROM_TOOL_FACE, &[j, moved.iter().any(|l| changed(l)) as usize]));
        faces.push(Face { surface: t.surface.clone(), sense: t.sense, loops: moved.iter().map(|l| swap(l)).collect() });
    }
    // The underside looks along -n: its loops run the other way round.
    let (bs, bsense) = (tool.faces[plan.base[0]].surface.clone(), tool.faces[plan.base[0]].sense);
    for face in hanging {
        let loops = face.into_iter().map(|lp| lp.into_iter().rev().map(|(e, w)| (e, !w)).collect()).collect();
        faces.push(Face { surface: bs.clone(), sense: bsense, loops });
        frec.push(Record::new(FROM_TOOL_FACE, &[plan.base[0], 1]));
    }
    r.faces = faces;
    let (emap, vmap2) = unite_tidy::compact(&mut r);
    let mut history = Provenance { faces: frec, edges: vec![Record::new(0, &[]); r.edges.len()], vertices: vec![Record::new(0, &[]); r.vertices.len()] };
    for (old, &new) in emap.iter().enumerate() {
        if new != usize::MAX {
            history.edges[new] = erec[old].clone();
        }
    }
    for (old, &new) in vmap2.iter().enumerate() {
        if new != usize::MAX {
            history.vertices[new] = vrec[old].clone();
        }
    }
    // An underside made from several base faces (a revolved disc's halves) has more than one parent.
    if plan.base.len() > 1 {
        history = Provenance::default();
    }
    unite_tidy::outer_first(&mut r);
    crate::topology::check(&r, tol).map_err(|e| KernelError::Internal(format!("overhang union built an invalid solid: {e}")))?;
    if !unite_tidy::connected(&r) {
        return Err(KernelError::Internal("overhang union built separate shells".into()));
    }
    let (mb, mt, mr) = (crate::mass::measure(body, tol)?, crate::mass::measure(tool, tol)?, crate::mass::measure(&r, tol)?);
    let expected = mb.volume + mt.volume;
    if (mr.volume - expected).abs() > tol * mr.area + 1e-9 * expected.abs().max(1.) {
        return Err(KernelError::Internal(format!("overhang union built {}, body and boss hold {expected}", mr.volume)));
    }
    Ok((r, history))
}

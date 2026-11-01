//! Blind pockets in right-prism stock (Claude 2026-10-31, Bob's plates: a sketch profile cut part way into a Flat).
//!
//! Domain, decided from the geometry before anything is built: the body is a right prism along a cardinal axis k
//! (section A over [lo, hi], `section::prism`), the tool a right prism along k (section B over [tlo, thi]) that leaves
//! the body through its top cap (thi >= hi - tol) and ends inside it (lo + tol < tlo < hi - tol); a pocket from the
//! bottom is the same turned half a turn about a section axis. The result is exactly the stack of two prisms, A over
//! [lo, d] and A - B over [d, hi] (d = tlo), as one B-rep:
//! - the bottom cap A at lo, the top cap A - B at hi and the floor A n B at d (facing up);
//! - one wall per segment of A, from lo up to hi where the segment lies outside B and to d where it lies inside B
//!   (a stepped outline, its vertical steps where B's boundary meets it), so no seam crosses a wall;
//! - one pocket wall per piece of B's boundary inside A, from d to hi, facing into the pocket.
//! The 2D arrangement is the section domain's own: every edge split at its closed-form meetings with the other
//! section, meeting points within the tolerance one vertex, coincident pieces decided once by their directions
//! (B along A the same way: A's piece borders the floor; the other way: it rises to the top), every other piece
//! classified by its midpoint's winding. Loops chain by the least clockwise turn; a pinch (a vertex on one level's
//! boundary twice) is refused as non-manifold. Surfaces are the segments' own planes, cylinders and elliptic
//! cylinders along k. The result's volume is checked against area(A) (d - lo) + area(A - B) (hi - d).
use super::*;
use crate::neutral::{Edge, Face};

/// One piece of the arrangement: its 2D vertices and curve.
#[derive(Clone, Copy)]
struct Bit {
    v0: usize,
    v1: usize,
    seg: Seg,
}

/// A section loop's material lies on its left; the outward side face of a piece: (surface, sense).
fn wall_surface(seg: &Seg, k: usize, at: f64) -> Result<(Surface, bool)> {
    let (u, v) = axes(k);
    let lift = |q: P| { let mut w = [0.; 3]; w[u] = q.x; w[v] = q.y; w[k] = at; Vec3::new(w[0], w[1], w[2]) };
    let dir = |q: P| { let mut w = [0.; 3]; w[u] = q.x; w[v] = q.y; Vec3::new(w[0], w[1], w[2]) };
    Ok(match *seg {
        Seg::Line { a, b } => {
            let t = b.sub(a);
            let l = t.len();
            if !(l > 0.) {
                return Err(unsupported("blind pocket: a section line has no length"));
            }
            let t = t.scale(1. / l);
            (Surface::Plane { p: lift(a), n: dir(P::new(t.y, -t.x)), x: dir(t) }, true)
        }
        Seg::Arc { c, r, sweep, .. } => (Surface::Cylinder { p: lift(c), a: unit(k), x: unit(u), r }, sweep > 0.),
        Seg::Conic { c, u: cu, v: cv, sweep, .. } => {
            let turn = cu.cross(cv);
            let (eu, ev) = if turn > 0. { (cu, cv) } else { (cu, cv.scale(-1.)) };
            // Counter-clockwise travel about k keeps the material inside: the outward normal points away.
            let ccw = sweep.signum() * turn.signum() > 0.;
            (Surface::EllipticCylinder { p: lift(c), a: unit(k), u: dir(eu), v: dir(ev) }, ccw)
        }
    })
}

/// The piece as a 3D curve at level `at`, run the piece's way.
fn curve_at(seg: &Seg, k: usize, at: f64) -> Curve {
    let (u, v) = axes(k);
    let lift = |q: P| { let mut w = [0.; 3]; w[u] = q.x; w[v] = q.y; w[k] = at; Vec3::new(w[0], w[1], w[2]) };
    let along = |q: P| { let mut w = [0.; 3]; w[u] = q.x; w[v] = q.y; Vec3::new(w[0], w[1], w[2]) };
    match *seg {
        Seg::Line { a, b } => Curve::Line { a: lift(a), b: lift(b) },
        Seg::Arc { c, r, a0, sweep } if sweep > 0. => Curve::Arc { c: lift(c), n: unit(k), x: unit(u), r, t0: a0, t1: a0 + sweep },
        Seg::Arc { c, r, a0, sweep } => Curve::Arc { c: lift(c), n: unit(k).scale(-1.), x: unit(u), r, t0: -a0, t1: -a0 - sweep },
        Seg::Conic { c, u: cu, v: cv, a0, sweep } if sweep > 0. => Curve::Conic { c: lift(c), u: along(cu), v: along(cv), t0: a0, t1: a0 + sweep },
        Seg::Conic { c, u: cu, v: cv, a0, sweep } => Curve::Conic { c: lift(c), u: along(cu), v: along(cv).scale(-1.), t0: -a0, t1: -a0 - sweep },
    }
}

/// Chains pieces (run their own way, material on the left) into loops: at each vertex the least clockwise turn from
/// the way back. Returns indices into `bits`. A vertex met twice on one level is a pinch: refused.
fn chain(bits: &[Bit], vertices: usize) -> Result<Vec<Vec<usize>>> {
    let mut used = vec![false; bits.len()];
    let mut visited = vec![false; vertices];
    let mut loops = vec![];
    for first in 0..bits.len() {
        if used[first] {
            continue;
        }
        let mut lp = vec![];
        let mut at = first;
        loop {
            used[at] = true;
            if std::mem::replace(&mut visited[bits[at].v0], true) {
                return Err(unsupported("blind pocket touches a boundary at a point (non-manifold result)"));
            }
            lp.push(at);
            let end = bits[at].v1;
            let back = bits[at].seg.tangent(1.).scale(-1.);
            let back_angle = back.y.atan2(back.x);
            let mut best: Option<(f64, usize)> = None;
            for (j, q) in bits.iter().enumerate() {
                if q.v0 != end || (used[j] && j != first) {
                    continue;
                }
                let out = q.seg.tangent(0.);
                let mut turn = (back_angle - out.y.atan2(out.x)).rem_euclid(2. * PI);
                if turn <= 1e-12 {
                    turn = 2. * PI;
                }
                if let Some((b, _)) = best {
                    if (turn - b).abs() <= 1e-9 {
                        return Err(unsupported("blind pocket pieces leave a vertex along one direction"));
                    }
                }
                if best.is_none_or(|(b, _)| turn < b) {
                    best = Some((turn, j));
                }
            }
            match best {
                Some((_, j)) if j == first => break,
                Some((_, j)) => at = j,
                None => return Err(unsupported("blind pocket section does not close")),
            }
        }
        loops.push(lp);
    }
    Ok(loops)
}

/// The pocket for a tool leaving through the top cap (see the module notes); `a` and `b` the body's and tool's prisms.
pub(super) fn cut(body: &Body, a: &Prism, b: &Prism, k: usize, tol: f64) -> Result<Body> {
    let (lo, hi, d) = (a.lo, a.hi, b.lo);
    if !(b.hi >= hi - tol && d > lo + tol && d < hi - tol) {
        return Err(unsupported("blind pocket: the tool must leave through the top and end inside the stock"));
    }
    // The arrangement, one A segment at a time (its pieces in order along it).
    let mut pool = vec![];
    let mut walls: Vec<Vec<Bit>> = vec![];
    let mut a_loops: Vec<Vec<usize>> = vec![];
    for (lp, tags) in a.loops.iter().zip(&a.tags) {
        let mut ids = vec![];
        for (seg, tag) in lp.iter().zip(tags) {
            let pieces = split(&[vec![*seg]], &[vec![*tag]], &b.loops, &mut pool, tol)?;
            ids.push(walls.len());
            walls.push(pieces.iter().map(|p| Bit { v0: p.v0, v1: p.v1, seg: p.seg }).collect());
        }
        a_loops.push(ids);
    }
    let pb: Vec<Bit> = split(&b.loops, &b.tags, &a.loops, &mut pool, tol)?.iter().map(|p| Bit { v0: p.v0, v1: p.v1, seg: p.seg }).collect();
    if walls.iter().map(Vec::len).sum::<usize>() + pb.len() > 4 * LIMIT {
        return Err(unsupported("blind pocket exceeds its work limit"));
    }
    // Classify: each A piece inside B (its wall stops at d) or outside (it rises to hi); B pieces inside A bound the
    // pocket. A B piece on an A piece is decided with it.
    let mut b_twin = vec![false; pb.len()];
    let mut inside: Vec<Vec<bool>> = vec![];
    for pieces in &walls {
        let mut row = vec![];
        for p in pieces {
            let mid = p.seg.at(0.5);
            let twin = pb.iter().position(|q| {
                ((q.v0 == p.v0 && q.v1 == p.v1) || (q.v0 == p.v1 && q.v1 == p.v0)) && q.seg.at(0.5).dist(mid) <= tol
            });
            if let Some(j) = twin {
                b_twin[j] = true;
                row.push(pb[j].seg.tangent(0.5).dot(p.seg.tangent(0.5)) > 0.);
                continue;
            }
            if b.loops.iter().flatten().any(|s| s.distance(mid) <= tol) {
                return Err(unsupported("blind pocket: section boundaries meet along an unmatched piece"));
            }
            row.push(winding(&b.loops, mid)? != 0);
        }
        inside.push(row);
    }
    let mut pocket = vec![];
    for (j, q) in pb.iter().enumerate() {
        if b_twin[j] {
            continue;
        }
        let mid = q.seg.at(0.5);
        if a.loops.iter().flatten().any(|s| s.distance(mid) <= tol) {
            return Err(unsupported("blind pocket: section boundaries meet along an unmatched piece"));
        }
        if winding(&a.loops, mid)? != 0 {
            pocket.push(*q);
        }
    }
    if pocket.is_empty() && inside.iter().flatten().all(|i| !i) {
        return Err(unsupported("blind pocket: the tool does not meet the stock's section"));
    }
    // 3D vertices: a 2D vertex at a level (0 lo, 1 d, 2 hi).
    let levels = [lo, d, hi];
    let (u, v) = axes(k);
    let mut vertices: Vec<Vec3> = vec![];
    let mut vmap = std::collections::HashMap::new();
    let mut vertex = |w: usize, l: usize, vertices: &mut Vec<Vec3>| -> usize {
        *vmap.entry((w, l)).or_insert_with(|| {
            let mut c = [0.; 3];
            c[u] = pool[w].x;
            c[v] = pool[w].y;
            c[k] = levels[l];
            vertices.push(Vec3::new(c[0], c[1], c[2]));
            vertices.len() - 1
        })
    };
    let mut edges: Vec<Edge> = vec![];
    let add = |curve: Curve, s: usize, e: usize, edges: &mut Vec<Edge>| { edges.push(Edge { curve, start: s, end: e }); edges.len() - 1 };
    // Level edges.
    let mut bottom = vec![];
    for (s, pieces) in walls.iter().enumerate() {
        let seg = a.loops.iter().flatten().nth(s).copied().unwrap();
        let (w0, w1) = (pieces[0].v0, pieces[pieces.len() - 1].v1);
        let (s0, s1) = (vertex(w0, 0, &mut vertices), vertex(w1, 0, &mut vertices));
        bottom.push(add(curve_at(&seg, k, lo), s0, s1, &mut edges));
    }
    let mut piece_edge: Vec<Vec<(usize, usize)>> = vec![]; // per A piece: (edge, level)
    for (pieces, row) in walls.iter().zip(&inside) {
        let mut out = vec![];
        for (p, &ins) in pieces.iter().zip(row) {
            let l = if ins { 1 } else { 2 };
            let (s0, s1) = (vertex(p.v0, l, &mut vertices), vertex(p.v1, l, &mut vertices));
            out.push((add(curve_at(&p.seg, k, levels[l]), s0, s1, &mut edges), l));
        }
        piece_edge.push(out);
    }
    let mut pocket_edges = vec![]; // (bottom at d, top at hi), both run B's way
    for q in &pocket {
        let (b0, b1) = (vertex(q.v0, 1, &mut vertices), vertex(q.v1, 1, &mut vertices));
        let (t0, t1) = (vertex(q.v0, 2, &mut vertices), vertex(q.v1, 2, &mut vertices));
        pocket_edges.push((add(curve_at(&q.seg, k, d), b0, b1, &mut edges), add(curve_at(&q.seg, k, hi), t0, t1, &mut edges)));
    }
    // Wall outlines as requests: level edges and vertical runs (2D vertex, from level, to level).
    #[derive(Clone, Copy)]
    enum Step {
        Edge(usize, bool),
        Rise(usize, usize, usize),
    }
    let mut outlines: Vec<(Surface, bool, Vec<Step>)> = vec![];
    for (s, pieces) in walls.iter().enumerate() {
        let seg = a.loops.iter().flatten().nth(s).copied().unwrap();
        let (surface, sense) = wall_surface(&seg, k, lo)?;
        let pe = &piece_edge[s];
        let m = pieces.len();
        let mut steps = vec![Step::Edge(bottom[s], true), Step::Rise(pieces[m - 1].v1, 0, pe[m - 1].1)];
        for i in (0..m).rev() {
            steps.push(Step::Edge(pe[i].0, false));
            if i > 0 && pe[i].1 != pe[i - 1].1 {
                steps.push(Step::Rise(pieces[i].v0, pe[i].1, pe[i - 1].1));
            }
        }
        steps.push(Step::Rise(pieces[0].v0, pe[0].1, 0));
        outlines.push((surface, sense, steps));
    }
    for (q, &(eb, et)) in pocket.iter().zip(&pocket_edges) {
        // The pocket wall bounds A - B, whose boundary runs this piece backwards: its outward side faces into B.
        let (surface, sense) = wall_surface(&q.seg.reversed(), k, lo)?;
        outlines.push((surface, sense, vec![Step::Edge(eb, false), Step::Rise(q.v0, 1, 2), Step::Edge(et, true), Step::Rise(q.v1, 2, 1)]));
    }
    // Vertical edges: unit spans (w, 0-1) and (w, 1-2), joined into one lo..hi edge where the same two faces use both.
    let mut users: std::collections::HashMap<(usize, usize), Vec<usize>> = std::collections::HashMap::new();
    for (f, (_, _, steps)) in outlines.iter().enumerate() {
        for st in steps {
            if let Step::Rise(w, a0, a1) = *st {
                for l in a0.min(a1)..a0.max(a1) {
                    users.entry((w, l)).or_default().push(f);
                }
            }
        }
    }
    let mut vertical: std::collections::HashMap<(usize, usize, usize), usize> = std::collections::HashMap::new();
    let mut keys: Vec<_> = users.keys().copied().collect();
    keys.sort();
    for (w, l) in keys {
        if vertical.contains_key(&(w, l, l + 1)) || (l == 1 && vertical.contains_key(&(w, 0, 2))) {
            continue;
        }
        let mut faces = users[&(w, l)].clone();
        faces.sort();
        if faces.len() != 2 {
            return Err(unsupported("blind pocket: a vertical edge is not shared by exactly two faces"));
        }
        let joined = l == 0 && users.get(&(w, 1)).is_some_and(|g| { let mut g = g.clone(); g.sort(); g == faces });
        let (l0, l1) = if joined { (0, 2) } else { (l, l + 1) };
        let (s, e) = (vertex(w, l0, &mut vertices), vertex(w, l1, &mut vertices));
        let id = add(Curve::Line { a: vertices[s], b: vertices[e] }, s, e, &mut edges);
        vertical.insert((w, l0, l1), id);
    }
    let rise = |w: usize, a0: usize, a1: usize| -> Result<Vec<(usize, bool)>> {
        let up = a1 > a0;
        let (lo_l, hi_l) = (a0.min(a1), a0.max(a1));
        let mut out = vec![];
        let mut l = lo_l;
        while l < hi_l {
            if let Some(&e) = vertical.get(&(w, 0, 2)).filter(|_| l == 0 && hi_l == 2) {
                out.push((e, up));
                l = 2;
            } else {
                let e = *vertical.get(&(w, l, l + 1)).ok_or_else(|| KernelError::Internal("blind pocket: a vertical edge is missing".into()))?;
                out.push((e, up));
                l += 1;
            }
        }
        if !up {
            out.reverse();
        }
        Ok(out)
    };
    let mut faces = vec![];
    for (surface, sense, steps) in outlines {
        let mut lp = vec![];
        for st in steps {
            match st {
                Step::Edge(e, f) => lp.push((e, f)),
                Step::Rise(w, a0, a1) => lp.extend(rise(w, a0, a1)?),
            }
        }
        faces.push(Face { surface, sense, loops: vec![lp] });
    }
    // The caps: the bottom (A, facing down: A's loops backwards), the top (A - B) and the floor (A n B), chained.
    let cap = |level: f64, n: Vec3| Surface::Plane { p: { let mut c = [0.; 3]; c[k] = level; Vec3::new(c[0], c[1], c[2]) }, n, x: unit(u) };
    faces.push(Face {
        surface: cap(lo, unit(k).scale(-1.)),
        sense: true,
        loops: a_loops.iter().map(|ids| ids.iter().rev().map(|s| (bottom[*s], false)).collect()).collect(),
    });
    let mut top_bits = vec![];
    let mut floor_bits = vec![];
    for (s, pieces) in walls.iter().enumerate() {
        for (i, p) in pieces.iter().enumerate() {
            let (e, l) = piece_edge[s][i];
            if l == 2 { top_bits.push((*p, e, true)) } else { floor_bits.push((*p, e, true)) }
        }
    }
    for (q, &(eb, et)) in pocket.iter().zip(&pocket_edges) {
        floor_bits.push((*q, eb, true));
        top_bits.push((Bit { v0: q.v1, v1: q.v0, seg: q.seg.reversed() }, et, false));
    }
    let mut area_top = 0.;
    for (bits, level) in [(&top_bits, hi), (&floor_bits, d)] {
        if bits.is_empty() {
            continue;
        }
        let plain: Vec<Bit> = bits.iter().map(|b| b.0).collect();
        let loops = chain(&plain, pool.len())?;
        if level == hi {
            area_top = loops.iter().map(|l| l.iter().map(|i| plain[*i].seg.area()).sum::<f64>()).sum();
        }
        faces.push(Face {
            surface: cap(level, unit(k)),
            sense: true,
            loops: loops.iter().map(|l| l.iter().map(|i| (bits[*i].1, bits[*i].2)).collect()).collect(),
        });
    }
    let out = Body { vertices, edges, faces };
    crate::topology::check(&out, tol)?;
    // The stack's volume: A over [lo, d] and A - B over [d, hi].
    let area_a: f64 = a.loops.iter().map(|l| loop_area(l)).sum();
    let expected = area_a * (d - lo) + area_top * (hi - d);
    // Points within the tolerance are one vertex, so a face may sit up to the tolerance off its edges: allowed that
    // times the boundary's area, and no more.
    let mass = crate::mass::measure(&out, tol)?;
    let got = mass.volume;
    if (got - expected).abs() > tol * mass.area + 1e-9 * expected.abs().max(1.) {
        return Err(KernelError::Internal(format!("blind pocket built {got}, its layers hold {expected}")));
    }
    let _ = body;
    Ok(out)
}

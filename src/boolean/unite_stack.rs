//! Union of two right prisms along one cardinal axis (Claude 2026-11-02, op 24): a boss overhanging a plate's edge,
//! flush with it, a bar through a plate, two overlapping blocks.
//!
//! Domain, decided from the geometry before anything is built: both operands are right prisms along the same cardinal
//! axis k (`section::prism`: caps on two levels normal to k, every other face a plane, cylinder or elliptic cylinder
//! parallel to k, area times extent equal to the measured volume), sections A over [alo, ahi] and B over [blo, bhi].
//! Levels within the tolerance are one level (the body's value is kept). The union is the stack of layers between
//! consecutive levels, each layer's section the union of the sections of the operands spanning it:
//! - the arrangement is the section domain's own (`section::split`): every edge split at its closed-form meetings with
//!   the other section, meeting points within the tolerance one vertex, coincident pieces decided once by their
//!   directions, every other piece classified by the winding of its midpoint. Each arrangement span knows on each of
//!   its sides whether it is inside A and inside B;
//! - a wall stands on a span in each layer where material lies on one side only; its outward side is the empty one;
//! - at each level the face looking up is the material below and not above, the face looking down the material above
//!   and not below (chained by the least clockwise turn; a vertex met twice on one level is refused as non-manifold);
//! - walls on one support facing one way are joined into one face across the edges they share, and edges on one
//!   support through a vertex only they use are joined (`unite_tidy`), so no seam is left on a plane or a cylinder
//!   (one seam per full cylinder band, as the extruder makes).
//! Surfaces are the operands' own (a wall the side face of its span's operand, a level face the cap of the operand
//! ending there); curves are the spans lifted to their levels. Refused by name: a layer with no material (the boss does
//! not reach the body), an edge or vertical line met by four faces (bodies touching only along a line), disconnected
//! material, ellipses meeting arcs or ellipses (as the section domain). The result's volume is checked against the
//! layers' areas times their thicknesses, computed from the arrangement.
use super::pocket::{curve_at, wall_surface};
use super::*;
use crate::boolean::unite_tidy;
use crate::neutral::{Edge, Face};
use std::collections::HashMap;

/// One span of the arrangement: its vertices and curve, whether A and B lie on its left and right, and the section
/// segment of each operand whose boundary it lies on.
struct Span {
    v0: usize,
    v1: usize,
    seg: Seg,
    a: [bool; 2],
    b: [bool; 2],
    from_a: Option<usize>,
    from_b: Option<usize>,
}

fn separate() -> KernelError {
    unsupported("the boss does not reach the body: their union would be two separate solids (one solid only)")
}

/// The first cardinal axis along which both are right prisms.
fn chart(body: &Body, tool: &Body, tol: f64) -> Result<Option<(usize, Prism, Prism)>> {
    for k in 0..3 {
        if !(prism_along(body, k) && prism_along(tool, k)) {
            continue;
        }
        if let (Some(a), Some(b)) = (prism(body, k, tol)?, prism(tool, k, tol)?) {
            return Ok(Some((k, a, b)));
        }
    }
    Ok(None)
}

/// Whether this domain takes the pair (both right prisms along one cardinal axis).
pub(in crate::boolean) fn selects(body: &Body, tool: &Body, tol: f64) -> bool {
    matches!(chart(body, tool, tol), Ok(Some(_)))
}

/// One support within the tolerance: collinear lines, arcs of one circle, arcs of one ellipse.
fn same_support(s: &Seg, t: &Seg, tol: f64) -> bool {
    match (*s, *t) {
        (Seg::Line { a, b }, Seg::Line { a: c, b: d }) => {
            let u = b.sub(a);
            let l = u.len();
            l > 0. && [c, d].iter().all(|q| (u.cross(q.sub(a)) / l).abs() <= tol)
        }
        (Seg::Arc { c, r, .. }, Seg::Arc { c: c2, r: r2, .. }) => c.dist(c2) <= tol && (r - r2).abs() <= tol,
        (Seg::Conic { c, u, v, .. }, Seg::Conic { c: c2, u: u2, v: v2, .. }) => c.dist(c2) <= tol && u.dist(u2) <= tol && v.dist(v2) <= tol,
        _ => false,
    }
}

/// Which way an oriented span faces on its support (+1 or -1): along the support's reference line, or round it.
fn facing(oriented: &Seg, reference: &Seg) -> i8 {
    match (*oriented, *reference) {
        (Seg::Line { a, b }, Seg::Line { a: c, b: d }) => if b.sub(a).dot(d.sub(c)) > 0. { 1 } else { -1 },
        (Seg::Arc { sweep, .. }, _) => if sweep > 0. { 1 } else { -1 },
        (Seg::Conic { u, v, sweep, .. }, _) => if sweep * u.cross(v) > 0. { 1 } else { -1 },
        _ => 1,
    }
}

/// Chains oriented pieces (material on the left) into loops by the least clockwise turn at each vertex; a vertex met
/// twice on one level is a pinch (refused).
fn chain(bits: &[(usize, usize, Seg)], vertices: usize) -> Result<Vec<Vec<usize>>> {
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
            if std::mem::replace(&mut visited[bits[at].0], true) {
                return Err(unsupported("the union would touch itself at a point of its outline (non-manifold result)"));
            }
            lp.push(at);
            let end = bits[at].1;
            let back = bits[at].2.tangent(1.).scale(-1.);
            let back_angle = back.y.atan2(back.x);
            let mut best: Option<(f64, usize)> = None;
            for (j, q) in bits.iter().enumerate() {
                if q.0 != end || (used[j] && j != first) {
                    continue;
                }
                let out = q.2.tangent(0.);
                let mut turn = (back_angle - out.y.atan2(out.x)).rem_euclid(2. * PI);
                if turn <= 1e-12 {
                    turn = 2. * PI;
                }
                if let Some((b, _)) = best {
                    if (turn - b).abs() <= 1e-9 {
                        return Err(unsupported("union outline pieces leave a vertex along one direction"));
                    }
                }
                if best.is_none_or(|(b, _)| turn < b) {
                    best = Some((turn, j));
                }
            }
            match best {
                Some((_, j)) if j == first => break,
                Some((_, j)) => at = j,
                None => return Err(unsupported("union outline does not close")),
            }
        }
        loops.push(lp);
    }
    Ok(loops)
}

/// The 3D entities of the stack, made once each: vertices (2D vertex, level), level edges (span, level) run the span's
/// way, rises (2D vertex, layer) from the layer's lower level to its upper.
struct Build<'a> {
    pool: &'a [P],
    spans: &'a [Span],
    levels: &'a [f64],
    k: usize,
    body: Body,
    vertex: HashMap<(usize, usize), usize>,
    level_edge: HashMap<(usize, usize), usize>,
    rise: HashMap<(usize, usize), usize>,
}

impl Build<'_> {
    fn vertex(&mut self, w: usize, j: usize) -> usize {
        if let Some(&v) = self.vertex.get(&(w, j)) {
            return v;
        }
        let (u, v) = axes(self.k);
        let mut c = [0.; 3];
        c[u] = self.pool[w].x;
        c[v] = self.pool[w].y;
        c[self.k] = self.levels[j];
        self.body.vertices.push(Vec3::new(c[0], c[1], c[2]));
        let id = self.body.vertices.len() - 1;
        self.vertex.insert((w, j), id);
        id
    }
    fn level(&mut self, s: usize, j: usize) -> usize {
        if let Some(&e) = self.level_edge.get(&(s, j)) {
            return e;
        }
        let (v0, v1) = (self.spans[s].v0, self.spans[s].v1);
        let (a, b) = (self.vertex(v0, j), self.vertex(v1, j));
        self.body.edges.push(Edge { curve: curve_at(&self.spans[s].seg, self.k, self.levels[j]), start: a, end: b });
        let id = self.body.edges.len() - 1;
        self.level_edge.insert((s, j), id);
        id
    }
    fn rise(&mut self, w: usize, l: usize) -> usize {
        if let Some(&e) = self.rise.get(&(w, l)) {
            return e;
        }
        let (a, b) = (self.vertex(w, l), self.vertex(w, l + 1));
        self.body.edges.push(Edge { curve: Curve::Line { a: self.body.vertices[a], b: self.body.vertices[b] }, start: a, end: b });
        let id = self.body.edges.len() - 1;
        self.rise.insert((w, l), id);
        id
    }
}

/// The union (see the module notes).
pub(in crate::boolean) fn union(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    let (k, a, b) = chart(body, tool, tol)?.ok_or_else(|| unsupported("the prism union needs both bodies to be right prisms along one axis"))?;
    let mut levels = vec![a.lo, a.hi];
    for z in [b.lo, b.hi] {
        if levels.iter().all(|l| (l - z).abs() > tol) {
            levels.push(z);
        }
    }
    levels.sort_by(f64::total_cmp);
    let at = |z: f64| levels.iter().position(|l| (l - z).abs() <= tol).unwrap();
    let (a0, a1, b0, b1) = (at(a.lo), at(a.hi), at(b.lo), at(b.hi));
    let layers = levels.len() - 1;
    let active = |l: usize| (a0 <= l && l < a1, b0 <= l && l < b1);
    if (0..layers).any(|l| active(l) == (false, false)) {
        return Err(separate());
    }
    // The arrangement of both sections.
    let mut b_tags = b.tags.clone();
    b_tags.iter_mut().flatten().for_each(|t| t.tool = true);
    let mut pool = vec![];
    let pa = split(&a.loops, &a.tags, &b.loops, &mut pool, tol)?;
    let pb = split(&b.loops, &b_tags, &a.loops, &mut pool, tol)?;
    if pa.len() + pb.len() > 4 * LIMIT {
        return Err(unsupported("prism union exceeds its work limit"));
    }
    let mut spans = vec![];
    let mut twinned = vec![false; pb.len()];
    for p in &pa {
        let mid = p.seg.at(0.5);
        let twin = pb.iter().position(|q| ((q.v0 == p.v0 && q.v1 == p.v1) || (q.v0 == p.v1 && q.v1 == p.v0)) && q.seg.at(0.5).dist(mid) <= tol);
        let (inside_b, from_b) = match twin {
            Some(j) => {
                twinned[j] = true;
                let same = pb[j].seg.tangent(0.5).dot(p.seg.tangent(0.5)) > 0.;
                ([same, !same], Some(pb[j].tag.seg))
            }
            None => {
                if b.loops.iter().flatten().any(|s| s.distance(mid) <= tol) {
                    return Err(unsupported("prism union: section boundaries meet along an unmatched piece"));
                }
                let w = winding(&b.loops, mid)? != 0;
                ([w, w], None)
            }
        };
        spans.push(Span { v0: p.v0, v1: p.v1, seg: p.seg, a: [true, false], b: inside_b, from_a: Some(p.tag.seg), from_b });
    }
    for (j, q) in pb.iter().enumerate() {
        if twinned[j] {
            continue;
        }
        let mid = q.seg.at(0.5);
        if a.loops.iter().flatten().any(|s| s.distance(mid) <= tol) {
            return Err(unsupported("prism union: section boundaries meet along an unmatched piece"));
        }
        let w = winding(&a.loops, mid)? != 0;
        spans.push(Span { v0: q.v0, v1: q.v1, seg: q.seg, a: [w, w], b: [true, false], from_a: None, from_b: Some(q.tag.seg) });
    }
    let solid = |l: usize, s: &Span, side: usize| {
        let (ia, ib) = active(l);
        (ia && s.a[side]) || (ib && s.b[side])
    };
    // Supports, and wall groups (support, facing).
    let mut support = vec![0; spans.len()];
    let mut reps: Vec<usize> = vec![];
    for i in 0..spans.len() {
        support[i] = match reps.iter().position(|&r| same_support(&spans[r].seg, &spans[i].seg, tol)) {
            Some(s) => s,
            None => {
                reps.push(i);
                reps.len() - 1
            }
        };
    }
    // A wall's surface: the side face of the operand whose boundary carries it (A first), else the span's own.
    let side_face = |s: &Span, left: bool, l: usize| -> Option<(Surface, bool)> {
        let (ia, ib) = active(l);
        let side = if left { 0 } else { 1 };
        let from = |seg: Option<usize>, prism: &Prism, owner: &Body| {
            seg.filter(|&g| g != MIXED).map(|g| prism.origins[g].face).filter(|&f| f != usize::MAX).map(|f| (owner.faces[f].surface.clone(), owner.faces[f].sense))
        };
        if ia && s.a[side] && s.a[0] != s.a[1] {
            if let Some(found) = from(s.from_a, &a, body) {
                return Some(found);
            }
        }
        if ib && s.b[side] && s.b[0] != s.b[1] {
            return from(s.from_b, &b, tool);
        }
        None
    };
    let mut build = Build { pool: &pool, spans: &spans, levels: &levels, k, body: Body::default(), vertex: HashMap::new(), level_edge: HashMap::new(), rise: HashMap::new() };
    let mut faces = vec![];
    let mut group: Vec<Option<usize>> = vec![];
    let mut groups: HashMap<(usize, i8), (usize, Option<(Surface, bool)>)> = HashMap::new();
    let mut cells = vec![];
    for (i, s) in spans.iter().enumerate() {
        for l in 0..layers {
            let (left, right) = (solid(l, s, 0), solid(l, s, 1));
            if left == right {
                continue;
            }
            let oriented = if left { s.seg } else { s.seg.reversed() };
            let key = (support[i], facing(&oriented, &spans[reps[support[i]]].seg));
            let count = groups.len();
            let entry = groups.entry(key).or_insert((count, None));
            if entry.1.is_none() {
                entry.1 = side_face(s, left, l);
            }
            cells.push((i, l, left, entry.0, key));
        }
    }
    let mut surfaces: HashMap<usize, (Surface, bool)> = HashMap::new();
    for &(i, l, left, g, key) in &cells {
        if !surfaces.contains_key(&g) {
            let found = groups[&key].1.clone();
            let oriented = if left { spans[i].seg } else { spans[i].seg.reversed() };
            surfaces.insert(g, match found { Some(f) => f, None => wall_surface(&oriented, k, levels[0])? });
        }
        let s = &spans[i];
        let lp = if left {
            vec![(build.level(i, l), true), (build.rise(s.v1, l), true), (build.level(i, l + 1), false), (build.rise(s.v0, l), false)]
        } else {
            vec![(build.level(i, l), false), (build.rise(s.v0, l), true), (build.level(i, l + 1), true), (build.rise(s.v1, l), false)]
        };
        let (surface, sense) = surfaces[&g].clone();
        faces.push(Face { surface, sense, loops: vec![lp] });
        group.push(Some(g));
    }
    // The level faces, and each layer's area for the volume check.
    let cap_face = |prism: &Prism, owner: &Body, top: bool| -> Option<(Surface, bool)> {
        let list = if top { &prism.top } else { &prism.bottom };
        list.first().map(|&f| (owner.faces[f].surface.clone(), owner.faces[f].sense))
    };
    for j in 0..=layers {
        for up in [true, false] {
            let region = |s: &Span, side: usize| {
                let below = j > 0 && solid(j - 1, s, side);
                let above = j < layers && solid(j, s, side);
                if up { below && !above } else { above && !below }
            };
            let mut bits = vec![];
            let mut which = vec![];
            for (i, s) in spans.iter().enumerate() {
                let (left, right) = (region(s, 0), region(s, 1));
                if left != right {
                    bits.push(if left { (s.v0, s.v1, s.seg) } else { (s.v1, s.v0, s.seg.reversed()) });
                    which.push((i, left));
                }
            }
            if bits.is_empty() {
                continue;
            }
            let loops = chain(&bits, pool.len())?;
            let area = |lp: &Vec<usize>| lp.iter().map(|&x| bits[x].2.area()).sum::<f64>();
            let (outers, holes): (Vec<&Vec<usize>>, Vec<&Vec<usize>>) = loops.iter().partition(|lp| area(lp) > 0.);
            let mut owned: Vec<Vec<&Vec<usize>>> = outers.iter().map(|o| vec![*o]).collect();
            for h in holes {
                let probe = bits[h[0]].2.at(0.5);
                let mut best: Option<(f64, usize)> = None;
                for (o, outer) in outers.iter().enumerate() {
                    let segs: Vec<Seg> = outer.iter().map(|&x| bits[x].2).collect();
                    if winding(&[segs], probe)? != 0 && best.is_none_or(|(a, _)| area(outer) < a) {
                        best = Some((area(outer), o));
                    }
                }
                let (_, o) = best.ok_or_else(|| KernelError::Internal("prism union: a hole lies in no outline".into()))?;
                owned[o].push(h);
            }
            // The cap of an operand ending (up) or starting (down) at this level; A first.
            let (surface, sense) = if up {
                if j == a1 { cap_face(&a, body, true) } else { cap_face(&b, tool, true) }
            } else if j == a0 {
                cap_face(&a, body, false)
            } else {
                cap_face(&b, tool, false)
            }
            .ok_or_else(|| KernelError::Internal("prism union: a level face has no operand cap".into()))?;
            for face_loops in owned {
                let mut out = vec![];
                for lp in face_loops {
                    let mut occ: Vec<(usize, bool)> = lp.iter().map(|&x| (build.level(which[x].0, j), which[x].1)).collect();
                    if !up {
                        occ = occ.iter().rev().map(|&(e, w)| (e, !w)).collect();
                    }
                    out.push(occ);
                }
                faces.push(Face { surface: surface.clone(), sense, loops: out });
                group.push(None);
            }
        }
    }
    let mut result = build.body;
    result.faces = faces;
    // Every edge must bound exactly two faces, once each way: otherwise the bodies only touch along a line.
    let mut uses = vec![(0usize, 0usize); result.edges.len()];
    for &(e, w) in result.faces.iter().flat_map(|f| f.loops.iter().flatten()) {
        if w { uses[e].0 += 1 } else { uses[e].1 += 1 }
    }
    if uses.iter().any(|&u| u != (1, 1)) {
        return Err(unsupported("the bodies would touch only along an edge (non-manifold result)"));
    }
    let wraps = |edges: &[Edge], lp: &[(usize, bool)]| {
        let net: f64 = lp
            .iter()
            .map(|&(e, w)| {
                let turn = match edges[e].curve {
                    Curve::Arc { n, t0, t1, .. } => (t1 - t0) * coord(n, k).signum(),
                    Curve::Conic { u: cu, v: cv, t0, t1, .. } => (t1 - t0) * coord(cu.cross(cv), k).signum(),
                    _ => 0.,
                };
                if w { turn } else { -turn }
            })
            .sum();
        net.abs() > PI
    };
    unite_tidy::join_faces(&mut result, &group, &wraps)?;
    unite_tidy::compact(&mut result);
    unite_tidy::join_edges(&mut result, tol);
    unite_tidy::outer_first(&mut result);
    crate::topology::check(&result, tol).map_err(|e| KernelError::Internal(format!("prism union built an invalid solid: {e}")))?;
    if !unite_tidy::connected(&result) {
        return Err(separate());
    }
    // The stack's volume: each layer's union area times its thickness.
    let mut expected = 0.;
    for l in 0..layers {
        let mut bits = vec![];
        for s in &spans {
            let (left, right) = (solid(l, s, 0), solid(l, s, 1));
            if left != right {
                bits.push(if left { s.seg } else { s.seg.reversed() });
            }
        }
        expected += bits.iter().map(Seg::area).sum::<f64>() * (levels[l + 1] - levels[l]);
    }
    let mass = crate::mass::measure(&result, tol)?;
    if (mass.volume - expected).abs() > tol * mass.area + 1e-9 * expected.abs().max(1.) {
        return Err(KernelError::Internal(format!("prism union built {}, its layers hold {expected}", mass.volume)));
    }
    Ok(result)
}

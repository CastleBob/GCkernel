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

/// One prism of a layered operand: its section (material left of every loop), its extent along k, each segment's tag
/// (for `split`) and the face it was read from.
struct Slab {
    loops: Vec<Vec<Seg>>,
    tags: Vec<Vec<Tag>>,
    lo: f64,
    hi: f64,
    faces: Vec<usize>,
}

/// A layered operand: its slabs and its caps (level, looking up, face).
struct Layered {
    slabs: Vec<Slab>,
    caps: Vec<(f64, bool, usize)>,
}

/// One span of the arrangement: its vertices and curve, and per slab whether the slab lies on its left and right and,
/// when the span lies on the slab's boundary, the slab's piece it came from.
struct Span {
    v0: usize,
    v1: usize,
    seg: Seg,
    member: Vec<[bool; 2]>,
    from: Vec<Option<usize>>,
}

fn separate() -> KernelError {
    unsupported("the boss does not reach the body: their union would be two separate solids (one solid only)")
}

/// Pairs each downward crossing with the next upward one along a periodic support (gap in (0, 2 pi]).
fn pair_round(downs: &[P], ups: &[P], param: &dyn Fn(P) -> f64, forward: bool) -> Result<Vec<(f64, f64)>> {
    let mut used = vec![false; ups.len()];
    let mut out = vec![];
    for d in downs {
        let a0 = param(*d);
        let gap = |q: &P| {
            let g = ((param(*q) - a0) * if forward { 1. } else { -1. }).rem_euclid(2. * PI);
            if g <= 1e-12 { 2. * PI } else { g }
        };
        let j = (0..ups.len())
            .filter(|&j| !used[j])
            .min_by(|&i, &j| gap(&ups[i]).total_cmp(&gap(&ups[j])))
            .ok_or_else(|| unsupported("layered union: a wall's crossings do not pair"))?;
        used[j] = true;
        let g = gap(&ups[j]);
        out.push((a0, if forward { g } else { -g }));
    }
    Ok(out)
}

/// The body as a stack of right prisms along k (see the module notes), or None when it is not one.
fn layered(body: &Body, k: usize, tool: bool, tol: f64) -> Result<Option<Layered>> {
    let (u, v) = axes(k);
    let kk = unit(k);
    let square = |w: Vec3| coord(w, k).abs() >= 1. - 1e-12 && coord(w, u).abs() <= 1e-12 && coord(w, v).abs() <= 1e-12;
    let mut caps = vec![];
    let mut walls = vec![];
    for (i, f) in body.faces.iter().enumerate() {
        match f.surface {
            Surface::Plane { p, n, .. } => {
                if square(n) {
                    let out = coord(n, k) * if f.sense { 1. } else { -1. };
                    caps.push((coord(p, k), out > 0., i));
                } else if coord(n, k).abs() <= 1e-12 {
                    walls.push(i);
                } else {
                    return Ok(None);
                }
            }
            Surface::Cylinder { a, .. } | Surface::EllipticCylinder { a, .. } if square(a) => walls.push(i),
            _ => return Ok(None),
        }
    }
    // Wall edges run square to k (at one level) or along it.
    let level = |c: &Curve| -> bool {
        match *c {
            Curve::Line { a, b } => (coord(a, k) - coord(b, k)).abs() <= tol,
            Curve::Arc { n, .. } => square(n),
            Curve::Conic { u: cu, v: cv, .. } => coord(cu, k).abs() <= 1e-9 * cu.length() && coord(cv, k).abs() <= 1e-9 * cv.length(),
            _ => false,
        }
    };
    let vertical = |c: &Curve| matches!(*c, Curve::Line { a, b } if b.sub(a).cross(kk).length() <= 1e-12 * b.sub(a).length());
    for &w in &walls {
        for &(e, _) in body.faces[w].loops.iter().flatten() {
            let c = &body.edges[e].curve;
            if !level(c) && !vertical(c) {
                return Ok(None);
            }
        }
    }
    let mut levels: Vec<f64> = vec![];
    for &(z, _, _) in &caps {
        if levels.iter().all(|l| (l - z).abs() > tol) {
            levels.push(z);
        }
    }
    levels.sort_by(f64::total_cmp);
    if levels.len() < 2 || levels.len() > 64 {
        return Ok(None);
    }
    let to2 = |q: Vec3| P::new(coord(q, u), coord(q, v));
    let mut slabs = vec![];
    let mut area_volume = 0.;
    for l in 0..levels.len() - 1 {
        let (lo, hi) = (levels[l], levels[l + 1]);
        let mid = 0.5 * (lo + hi);
        // Every wall's pieces at the middle: from each vertical edge crossing it downward to the next crossing it upward
        // along the wall's travel (k x outward: material on the left).
        let mut pieces: Vec<(Seg, usize)> = vec![];
        for &w in &walls {
            let face = &body.faces[w];
            let mut crossings: Vec<(P, bool)> = vec![];
            for &(e, fwd) in face.loops.iter().flatten() {
                let edge = &body.edges[e];
                if !vertical(&edge.curve) {
                    continue;
                }
                let (s, t) = if fwd { (edge.start, edge.end) } else { (edge.end, edge.start) };
                let (zs, zt) = (coord(body.vertices[s], k), coord(body.vertices[t], k));
                if zs.min(zt) < mid && zs.max(zt) > mid {
                    crossings.push((to2(body.vertices[s]), zt > zs));
                }
            }
            if crossings.is_empty() {
                continue;
            }
            let downs: Vec<P> = crossings.iter().filter(|c| !c.1).map(|c| c.0).collect();
            let ups: Vec<P> = crossings.iter().filter(|c| c.1).map(|c| c.0).collect();
            if downs.len() != ups.len() {
                return Ok(None);
            }
            match face.surface {
                Surface::Plane { p, n, .. } => {
                    let out = if face.sense { n } else { n.scale(-1.) };
                    let t = kk.cross(out);
                    let (o, t2) = (to2(p), P::new(coord(t, u), coord(t, v)));
                    let along = |q: P| q.sub(o).dot(t2);
                    let mut all: Vec<(f64, bool, P)> = crossings.iter().map(|c| (along(c.0), c.1, c.0)).collect();
                    all.sort_by(|a, b| a.0.total_cmp(&b.0));
                    for pair in all.chunks(2) {
                        let [(_, false, a), (_, true, b)] = pair[..] else { return Ok(None) };
                        pieces.push((Seg::Line { a, b }, w));
                    }
                }
                Surface::Cylinder { p, r, .. } => {
                    let c = to2(p);
                    let angle = |q: P| (q.y - c.y).atan2(q.x - c.x);
                    for (a0, sweep) in pair_round(&downs, &ups, &angle, face.sense)? {
                        pieces.push((Seg::Arc { c, r, a0, sweep }, w));
                    }
                }
                Surface::EllipticCylinder { p, u: eu, v: ev, .. } => {
                    let c = to2(p);
                    let (u2, v2) = (P::new(coord(eu, u), coord(eu, v)), P::new(coord(ev, u), coord(ev, v)));
                    let forward = face.sense == (u2.cross(v2) > 0.);
                    let param = |q: P| {
                        let (cc, nn) = unit_coords(u2, v2, q.sub(c));
                        nn.atan2(cc)
                    };
                    for (a0, sweep) in pair_round(&downs, &ups, &param, forward)? {
                        pieces.push((Seg::Conic { c, u: u2, v: v2, a0, sweep }, w));
                    }
                }
                _ => return Ok(None),
            }
        }
        if pieces.len() > LIMIT {
            return Err(unsupported("layered union: a section exceeds its work limit"));
        }
        // Chain the pieces into loops by their ends, then join pieces of one support.
        let mut used = vec![false; pieces.len()];
        let (mut loops, mut tags) = (vec![], vec![]);
        let whole = |i: usize| Tag { tool, seg: i, s0: 0., s1: 1. };
        for first in 0..pieces.len() {
            if used[first] {
                continue;
            }
            used[first] = true;
            let (mut lp, mut tg) = (vec![pieces[first].0], vec![whole(first)]);
            while lp.last().unwrap().end().dist(lp[0].start()) > tol {
                let end = lp.last().unwrap().end();
                let Some(next) = (0..pieces.len()).find(|&j| !used[j] && pieces[j].0.start().dist(end) <= tol) else { return Ok(None) };
                used[next] = true;
                lp.push(pieces[next].0);
                tg.push(whole(next));
            }
            let (lp, tg) = merge(lp, tg, tol);
            loops.push(lp);
            tags.push(tg);
        }
        let area: f64 = loops.iter().map(|l| loop_area(l)).sum();
        if loops.is_empty() || area <= 0. {
            return Ok(None);
        }
        area_volume += area * (hi - lo);
        slabs.push(Slab { loops, tags, lo, hi, faces: pieces.iter().map(|p| p.1).collect() });
    }
    // The stack proof: the layers' areas times their thicknesses are the body's measured volume.
    let volume = crate::mass::measure(body, tol)?.volume;
    if (area_volume - volume).abs() > 1e-9 * volume.abs().max(1.) {
        return Ok(None);
    }
    Ok(Some(Layered { slabs, caps }))
}

/// The first cardinal axis along which both are layered.
fn chart(body: &Body, tool: &Body, tol: f64) -> Result<Option<(usize, Layered, Layered)>> {
    for k in 0..3 {
        if let (Ok(Some(a)), Ok(Some(b))) = (layered(body, k, false, tol), layered(tool, k, true, tol)) {
            return Ok(Some((k, a, b)));
        }
    }
    Ok(None)
}

/// Whether this domain takes the pair (both layered along one cardinal axis).
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
    let (k, a, b) = chart(body, tool, tol)?.ok_or_else(|| unsupported("the layered prism union needs both bodies to be stacks of right prisms along one axis"))?;
    // Every slab, the body's first; levels within the tolerance are one (the body's value kept).
    let slabs: Vec<(&Slab, &Body)> = a.slabs.iter().map(|s| (s, body)).chain(b.slabs.iter().map(|s| (s, tool))).collect();
    let mut levels: Vec<f64> = vec![];
    for (s, _) in &slabs {
        for z in [s.lo, s.hi] {
            if levels.iter().all(|l| (l - z).abs() > tol) {
                levels.push(z);
            }
        }
    }
    levels.sort_by(f64::total_cmp);
    let at = |z: f64| levels.iter().position(|l| (l - z).abs() <= tol).unwrap();
    let spans_of: Vec<(usize, usize)> = slabs.iter().map(|(s, _)| (at(s.lo), at(s.hi))).collect();
    let layers = levels.len() - 1;
    let active = |slab: usize, l: usize| spans_of[slab].0 <= l && l < spans_of[slab].1;
    if (0..layers).any(|l| (0..slabs.len()).all(|s| !active(s, l))) {
        return Err(separate());
    }
    // The arrangement of every slab's section against all the others'.
    let mut pool = vec![];
    let mut pieces = vec![];
    for (i, (s, _)) in slabs.iter().enumerate() {
        let others: Vec<Vec<Seg>> = slabs.iter().enumerate().filter(|(j, _)| *j != i).flat_map(|(_, (o, _))| o.loops.iter().cloned()).collect();
        for p in split(&s.loops, &s.tags, &others, &mut pool, tol)? {
            pieces.push((i, p));
        }
        if pieces.len() > 4 * LIMIT {
            return Err(unsupported("layered prism union exceeds its work limit"));
        }
    }
    let mut spans: Vec<Span> = vec![];
    for (i, p) in &pieces {
        let mid = p.seg.at(0.5);
        let found = spans.iter().position(|q| ((q.v0 == p.v0 && q.v1 == p.v1) || (q.v0 == p.v1 && q.v1 == p.v0)) && q.seg.at(0.5).dist(mid) <= tol);
        let j = match found {
            Some(j) => j,
            None => {
                spans.push(Span { v0: p.v0, v1: p.v1, seg: p.seg, member: vec![[false; 2]; slabs.len()], from: vec![None; slabs.len()] });
                spans.len() - 1
            }
        };
        // The slab lies left of its own boundary pieces: on the span's left when they run the same way.
        let same = spans[j].seg.tangent(0.5).dot(p.seg.tangent(0.5)) > 0.;
        if spans[j].from[*i].is_some() {
            return Err(unsupported("layered prism union: a section runs over itself"));
        }
        spans[j].member[*i] = [same, !same];
        spans[j].from[*i] = Some(p.tag.seg);
    }
    for span in &mut spans {
        let mid = span.seg.at(0.5);
        for (i, (s, _)) in slabs.iter().enumerate() {
            if span.from[i].is_some() {
                continue;
            }
            if s.loops.iter().flatten().any(|g| g.distance(mid) <= tol) {
                return Err(unsupported("layered prism union: section boundaries meet along an unmatched piece"));
            }
            let w = winding(&s.loops, mid)? != 0;
            span.member[i] = [w, w];
        }
    }
    let solid = |l: usize, s: &Span, side: usize| (0..slabs.len()).any(|i| active(i, l) && s.member[i][side]);
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
    // A wall's surface: the face its span was read from, of a slab active there with material on the wall's side (the
    // body's first), else the span's own.
    let side_face = |s: &Span, left: bool, l: usize| -> Option<(Surface, bool)> {
        let side = if left { 0 } else { 1 };
        (0..slabs.len()).find_map(|i| {
            let (slab, owner) = slabs[i];
            let g = s.from[i].filter(|&g| g != MIXED && active(i, l) && s.member[i][side])?;
            let f = &owner.faces[slab.faces[g]];
            Some((f.surface.clone(), f.sense))
        })
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
    // A cap of an operand at a level, looking the given way (the body's first).
    let cap_face = |level: f64, up: bool| -> Option<(Surface, bool)> {
        [(&a, body), (&b, tool)].iter().find_map(|(op, owner)| {
            op.caps.iter().find(|c| (c.0 - level).abs() <= tol && c.1 == up).map(|c| (owner.faces[c.2].surface.clone(), owner.faces[c.2].sense))
        })
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
            let (surface, sense) = cap_face(levels[j], up).ok_or_else(|| KernelError::Internal("layered prism union: a level face has no operand cap".into()))?;
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
    // One shell: separate pieces (or an enclosed void, whose inner shell this domain does not build as one) are refused.
    if !unite_tidy::connected(&result) {
        return Err(unsupported("the union would leave separate solids or enclose a void (one solid without voids only)"));
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

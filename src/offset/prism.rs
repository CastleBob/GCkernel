//! Offsets and shells of straight prisms (Claude 2026-11-01). A prism is its profile A (the top cap's loops) swept
//! through its thickness L, so its offsets have closed forms:
//!
//! - inward by e (rounded or sharp): every point at least e inside A x [-L, 0] is (A eroded by e) x [-L + e, -e]; a
//!   prism of the planar offset. Sharp outward is the same with the sharp planar offset, L + 2d thick.
//! - outward by d, rounded (the true offset, the Minkowski sum with a ball): the distance to A x I is
//!   sqrt(dist_A^2 + dist_I^2), so the side wall is the planar offset A + d swept through [-L, 0] (rounding arcs at
//!   convex corners sweep cylinders about the old rulings) and each cap band is a tube of radius d about the cap's
//!   boundary: a quarter cylinder about each line, a quarter torus about each arc, a sphere about each convex corner,
//!   and at a concave line-line corner the two quarter cylinders meet along an ellipse in the corner's bisecting plane
//!   (two equal cylinders whose axes cross: s^2 sin^2(beta) + h^2 = d^2). The caps are A moved by d.
//!
//! Every layer is built face by face from shared vertices and edges; shells join an outer and an inner (void) layer,
//! with flat rims where a removed cap leaves the shell open.

use super::planar::{self, area2, check_nesting, check_simple, dir, offset_loop, right, Corner, Elem, OffsetLoop, Seg2, Tag, P2};
use super::Source;
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;
use std::f64::consts::FRAC_PI_2;

fn unsupported<T>(m: String) -> Result<T, KernelError> { Err(KernelError::Unsupported(m)) }

/// A prism found in a body, with every profile element's source face, edges and vertices.
pub struct Prism {
    pub o: Vec3,
    pub n: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub length: f64,
    pub top: usize,
    pub bottom: usize,
    /// The top cap's loops in the (u, v) frame, material on the left; element j of loop l is tagged Seg(l, j).
    pub loops: Vec<Vec<Elem>>,
    pub side_face: Vec<Vec<usize>>,
    pub top_edge: Vec<Vec<usize>>,
    pub bottom_edge: Vec<Vec<usize>>,
    /// The ruling, top vertex and bottom vertex at the start of element j.
    pub ruling: Vec<Vec<usize>>,
    pub top_vertex: Vec<Vec<usize>>,
    pub bottom_vertex: Vec<Vec<usize>>,
}

impl Prism {
    fn point(&self, q: P2, z: f64) -> Vec3 { self.o.add(self.u.scale(q.0)).add(self.v.scale(q.1)).add(self.n.scale(z)) }
    fn vector(&self, q: P2) -> Vec3 { self.u.scale(q.0).add(self.v.scale(q.1)) }
    fn curve(&self, s: &Seg2, z: f64) -> Curve {
        match *s {
            Seg2::Line { a, b } => Curve::Line { a: self.point(a, z), b: self.point(b, z) },
            Seg2::Arc { c, r, a0, sw } => Curve::Arc {
                c: self.point(c, z), n: if sw > 0.0 { self.n } else { self.n.scale(-1.0) }, x: self.vector(dir(a0)), r, t0: 0.0, t1: sw.abs(),
            },
        }
    }
    pub fn source(&self, role: &Role) -> Source {
        match *role {
            Role::Side(Tag::Seg(l, j)) => Source::Face(self.side_face[l][j]),
            Role::Side(Tag::Corner(l, j)) => Source::Edge(self.ruling[l][j]),
            Role::Side(Tag::Axis) => unreachable!("prism profiles have no axis"),
            Role::Cap { top } | Role::Rim { top } => Source::Face(if top { self.top } else { self.bottom }),
            Role::Pipe { top, l, j } => Source::Edge(if top { self.top_edge[l][j] } else { self.bottom_edge[l][j] }),
            Role::Ball { top, l, j } => Source::Vertex(if top { self.top_vertex[l][j] } else { self.bottom_vertex[l][j] }),
        }
    }
}

fn same_curve(a: &Curve, b: &Curve, tol: f64) -> bool {
    let mid = |c: &Curve| { let (t0, t1) = c.domain(); c.eval((t0 + t1) / 2.0).0 };
    let ends = (a.start().distance(b.start()) <= tol && a.end().distance(b.end()) <= tol)
        || (a.start().distance(b.end()) <= tol && a.end().distance(b.start()) <= tol);
    ends && mid(a).distance(mid(b)) <= tol
}

fn outward(body: &Body, f: usize) -> Option<Vec3> {
    match body.faces[f].surface {
        Surface::Plane { n, .. } => Some(if body.faces[f].sense { n } else { n.scale(-1.0) }),
        _ => None,
    }
}

/// The body as a straight prism whose caps include every face in `caps` (empty: any caps), or why it is not one.
pub fn recognize(body: &Body, caps: &[usize], tol: f64) -> Result<Prism, KernelError> {
    let mut why = String::from("it has no plane face to be a cap");
    let mut proved = false;
    for top in 0..body.faces.len() {
        let Some(n) = outward(body, top) else { continue };
        match prove(body, top, n, tol) {
            Ok(p) => {
                if caps.iter().all(|&f| f == p.top || f == p.bottom) { return Ok(p); }
                proved = true;
            }
            Err(KernelError::Unsupported(m)) => if !proved { why = m },
            Err(e) => return Err(e),
        }
    }
    if proved { return unsupported("a straight prism whose removed faces are not the two caps of one sweep (only caps can be removed)".into()); }
    unsupported(format!("not a straight prism: {why}"))
}

fn prove(body: &Body, top: usize, n: Vec3, tol: f64) -> Result<Prism, KernelError> {
    let bottoms: Vec<usize> = (0..body.faces.len()).filter(|&f| outward(body, f).is_some_and(|m| m.add(n).length() < 1e-9)).collect();
    let [bottom] = bottoms[..] else { return unsupported(format!("{} faces face away from face {top}; a prism has one bottom", bottoms.len())) };
    let (Surface::Plane { p: pt, .. }, Surface::Plane { p: pb, .. }) = (&body.faces[top].surface, &body.faces[bottom].surface) else { unreachable!() };
    let length = pt.sub(*pb).dot(n);
    if !(length > tol) { return unsupported(format!("faces {top} and {bottom} are not apart")); }
    let swept = |s: &Surface| match *s {
        Surface::Plane { n: m, .. } => m.dot(n).abs() < 1e-9,
        Surface::Cylinder { a, .. } | Surface::EllipticCylinder { a, .. } => a.cross(n).length() < 1e-9,
        _ => false,
    };
    for (f, face) in body.faces.iter().enumerate() {
        if f != top && f != bottom && !swept(&face.surface) { return unsupported(format!("face {f} does not run straight through")); }
    }
    let o = *pt;
    let u = n.perpendicular();
    let v = n.cross(u);
    let flat = |q: Vec3| { let d = q.sub(o); (d.dot(u), d.dot(v)) };
    let mut prism = Prism { o, n, u, v, length, top, bottom, loops: Vec::new(), side_face: Vec::new(), top_edge: Vec::new(),
        bottom_edge: Vec::new(), ruling: Vec::new(), top_vertex: Vec::new(), bottom_vertex: Vec::new() };
    let mut used = vec![false; body.faces.len()];
    let mut bottom_edges = Vec::new();
    for (l, lp) in body.faces[top].loops.iter().enumerate() {
        let (mut elems, mut sides, mut tops, mut bots) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let (mut rules, mut tvs, mut bvs) = (Vec::new(), Vec::new(), Vec::new());
        for (j, &(e, fwd)) in lp.iter().enumerate() {
            let edge = &body.edges[e];
            let curve = if fwd { edge.curve.clone() } else { edge.curve.reversed() };
            let start = if fwd { edge.start } else { edge.end };
            let owners: Vec<usize> = (0..body.faces.len())
                .filter(|&f| f != top && body.faces[f].loops.iter().flatten().any(|&(g, _)| g == e)).collect();
            let [side] = owners[..] else { return unsupported(format!("top edge {e} bounds {} side faces", owners.len())) };
            if side == bottom || used[side] { return unsupported(format!("top edge {e} has no side face of its own")); }
            used[side] = true;
            let face = &body.faces[side];
            let [slp] = &face.loops[..] else { return unsupported(format!("side face {side} has {} loops", face.loops.len())) };
            if slp.len() != 4 { return unsupported(format!("side face {side} has {} edges, not 4", slp.len())); }
            let is_ruling = |g: usize| matches!(body.edges[g].curve, Curve::Line { a, b }
                if (b.sub(a).length() - length).abs() <= tol && b.sub(a).cross(n).length() <= tol);
            let Some(r) = slp.iter().map(|&(g, _)| g).find(|&g| g != e && is_ruling(g) && (body.edges[g].start == start || body.edges[g].end == start)) else {
                return unsupported(format!("side face {side} has no ruling through the thickness at the start of top edge {e}"));
            };
            let bv = if body.edges[r].start == start { body.edges[r].end } else { body.edges[r].start };
            let others: Vec<usize> = slp.iter().map(|&(g, _)| g).filter(|&g| g != e && !is_ruling(g)).collect();
            let [b] = others[..] else { return unsupported(format!("side face {side} is not bounded by two rulings and two caps' edges")) };
            if !same_curve(&body.edges[b].curve, &curve.translated(n.scale(-length)), tol) {
                return unsupported(format!("side face {side}'s bottom edge is not its top edge moved through the thickness"));
            }
            let seg = match curve {
                Curve::Line { a, b } => Seg2::Line { a: flat(a), b: flat(b) },
                Curve::Arc { c, n: an, x, r, t0, t1 } => {
                    if an.cross(n).length() > 1e-9 { return unsupported(format!("top edge {e} is an arc out of the cap's plane")); }
                    let q = Curve::point_at_angle(c, an, x, r, t0);
                    let (cc, s) = (flat(c), flat(q));
                    Seg2::Arc { c: cc, r, a0: (s.1 - cc.1).atan2(s.0 - cc.0), sw: (t1 - t0) * an.dot(n).signum() }
                }
                Curve::Conic { .. } => return unsupported(format!(
                    "side face {side} is an elliptic cylinder: an ellipse's offset is not an ellipse (no exact offset surface)")),
                _ => return unsupported(format!("top edge {e} is not a line or circular arc")),
            };
            elems.push(Elem { seg, tag: Tag::Seg(l, j), corner: Tag::Corner(l, j), mirrored: false, removed: false });
            sides.push(side);
            tops.push(e);
            bots.push(b);
            rules.push(r);
            tvs.push(start);
            bvs.push(bv);
            bottom_edges.push(b);
        }
        prism.loops.push(elems);
        prism.side_face.push(sides);
        prism.top_edge.push(tops);
        prism.bottom_edge.push(bots);
        prism.ruling.push(rules);
        prism.top_vertex.push(tvs);
        prism.bottom_vertex.push(bvs);
    }
    if let Some(f) = (0..body.faces.len()).find(|&f| f != top && f != bottom && !used[f]) {
        return unsupported(format!("face {f} is not the sweep of a top edge"));
    }
    let mut on_bottom: Vec<usize> = body.faces[bottom].loops.iter().flatten().map(|&(g, _)| g).collect();
    on_bottom.sort_unstable();
    bottom_edges.sort_unstable();
    if on_bottom != bottom_edges { return unsupported("the bottom cap is not the top cap moved through the thickness".into()); }
    // A regression check on the structural proof: the swept profile has the body's volume.
    let expected = length * profile_area(&prism.loops);
    let measured = crate::mass::measure(body, tol)?.volume;
    if (expected - measured).abs() > 1e-9 * expected.abs().max(measured.abs()) {
        return unsupported(format!("the body is not the sweep of its top (volume {measured} vs {expected})"));
    }
    Ok(prism)
}

/// The material area of a profile (outer loops count, holes subtract by orientation).
pub fn profile_area(loops: &[Vec<Elem>]) -> f64 { loops.iter().map(|l| area2(l)).sum::<f64>() / 2.0 }

/// What a built face is, in the source profile's terms.
#[derive(Clone, Copy, Debug)]
pub enum Role {
    /// Swept from a profile element (its tag).
    Side(Tag),
    Cap { top: bool },
    /// The rounding tube about a cap edge (source element j of loop l).
    Pipe { top: bool, l: usize, j: usize },
    /// The rounding sphere about a cap corner.
    Ball { top: bool, l: usize, j: usize },
    /// A flat rim where a removed cap left the shell open.
    Rim { top: bool },
}

/// One layer: `base` swept from z_lo to z_hi, rounded outward by `d` (0: the base itself); an open end leaves its
/// boundary edges unused (for a rim), a closed end is a cap (with a rounding band when d > 0).
pub struct Layer {
    pub body: Body,
    pub roles: Vec<Role>,
    pub open_hi: Vec<Vec<usize>>,
    pub open_lo: Vec<Vec<usize>>,
    pub volume: f64,
}

/// The planar offset of every loop (checked simple), or the loops themselves when d = 0.
pub fn offset_profile(loops: &[Vec<Elem>], d: f64, rounded: bool, tol: f64, what: &str) -> Result<Vec<OffsetLoop>, KernelError> {
    let out: Vec<OffsetLoop> = if d == 0.0 {
        loops.iter().map(|lp| OffsetLoop {
            elems: lp.clone(), seg_at: (0..lp.len()).collect(), corner_at: vec![None; lp.len()], corners: vec![Corner::Still; lp.len()],
            turns: (0..lp.len()).map(|j| { let (a, b) = (lp[(j + lp.len() - 1) % lp.len()].seg.tangent_end(), lp[j].seg.tangent_start());
                planar::cross(a, b).atan2(planar::dot(a, b)) }).collect(),
        }).collect()
    } else {
        loops.iter().enumerate().map(|(l, lp)| offset_loop(lp, &|_| d, rounded, tol, &format!("{what} loop {l}")))
            .collect::<Result<_, _>>()?
    };
    if d != 0.0 {
        let elems: Vec<Vec<Elem>> = out.iter().map(|o| o.elems.clone()).collect();
        check_simple(&elems, tol, what)?;
        check_nesting(loops, &elems, what)?;
        // Each loop keeps its sense (an outer loop stays counter-clockwise, a hole clockwise).
        for (l, (o, s)) in elems.iter().zip(loops).enumerate() {
            if area2(o).signum() != area2(s).signum() { return unsupported(format!("{what}: loop {l} turns inside out")); }
        }
    }
    Ok(out)
}

struct Builder<'a> {
    prism: &'a Prism,
    body: Body,
    roles: Vec<Role>,
}

impl Builder<'_> {
    fn vertex(&mut self, q: P2, z: f64) -> usize { self.body.vertices.push(self.prism.point(q, z)); self.body.vertices.len() - 1 }
    fn edge(&mut self, curve: Curve, start: usize, end: usize) -> usize {
        self.body.edges.push(Edge { curve, start, end });
        self.body.edges.len() - 1
    }
    fn face(&mut self, surface: Surface, sense: bool, loops: Vec<Vec<(usize, bool)>>, role: Role) {
        self.body.faces.push(Face { surface, sense, loops });
        self.roles.push(role);
    }
}

/// Builds one layer of the prism's frame (`prism` supplies the frame and nothing else here).
pub fn layer(prism: &Prism, base: &[Vec<Elem>], d: f64, z_lo: f64, z_hi: f64, lo_open: bool, hi_open: bool, tol: f64, what: &str)
    -> Result<Layer, KernelError> {
    let rings = offset_profile(base, d, true, tol, what)?;
    let mut b = Builder { prism, body: Body::default(), roles: Vec::new() };
    let (n, height) = (prism.n, z_hi - z_lo);
    let mut out_hi = Vec::new();
    let mut out_lo = Vec::new();
    let (mut cap_hi, mut cap_lo) = (Vec::new(), Vec::new());
    let (mut area_base, mut area_off, mut perimeter, mut k2) = (0.0, 0.0, 0.0, 0.0);
    for (l, ring) in rings.iter().enumerate() {
        let (src, o) = (&base[l], &ring.elems);
        let (m, kk) = (src.len(), o.len());
        area_base += area2(src) / 2.0;
        area_off += area2(o) / 2.0;
        let v_lo: Vec<usize> = o.iter().map(|e| b.vertex(e.seg.start(), z_lo)).collect();
        let v_hi: Vec<usize> = o.iter().map(|e| b.vertex(e.seg.start(), z_hi)).collect();
        let mut lo_e = Vec::with_capacity(kk);
        let mut hi_e = Vec::with_capacity(kk);
        let mut rule = Vec::with_capacity(kk);
        for k in 0..kk {
            let k1 = (k + 1) % kk;
            lo_e.push(b.edge(prism.curve(&o[k].seg, z_lo), v_lo[k], v_lo[k1]));
            hi_e.push(b.edge(prism.curve(&o[k].seg, z_hi), v_hi[k], v_hi[k1]));
            rule.push(b.edge(Curve::Line { a: prism.point(o[k].seg.start(), z_lo), b: prism.point(o[k].seg.start(), z_hi) }, v_lo[k], v_hi[k]));
        }
        for k in 0..kk {
            let (surface, sense) = match o[k].seg {
                Seg2::Line { a, b: e } => {
                    let t = planar::sub(e, a);
                    let t = planar::scale(t, 1.0 / planar::len(t));
                    (Surface::Plane { p: prism.point(a, z_lo), n: prism.vector(right(t)), x: prism.vector(t) }, true)
                }
                Seg2::Arc { c, r, sw, .. } => (Surface::Cylinder { p: prism.point(c, z_lo), a: n, x: prism.u, r }, sw > 0.0),
            };
            let lp = vec![(lo_e[k], true), (rule[(k + 1) % kk], true), (hi_e[k], false), (rule[k], false)];
            b.face(surface, sense, vec![lp], Role::Side(o[k].tag));
        }
        if d > 0.0 {
            perimeter += src.iter().map(|e| e.seg.length()).sum::<f64>();
            for j in 0..m {
                if let Seg2::Arc { sw, .. } = src[j].seg { k2 += sw / 2.0; }
                k2 += match ring.corners[j] {
                    Corner::Round => ring.turns[j] / 2.0,
                    Corner::Sharp => -(ring.turns[j].abs() / 2.0).tan(),
                    Corner::Tangent | Corner::Still => 0.0,
                };
            }
        }
        for (top, z, v_mid, mid_e) in [(true, z_hi, &v_hi, &hi_e), (false, z_lo, &v_lo, &lo_e)] {
            let open = if top { hi_open } else { lo_open };
            if open {
                if top { out_hi.push(mid_e.clone()) } else { out_lo.push(mid_e.clone()) }
                continue;
            }
            if d == 0.0 {
                let lp: Vec<(usize, bool)> = if top { mid_e.iter().map(|&e| (e, true)).collect() } else { mid_e.iter().rev().map(|&e| (e, false)).collect() };
                if top { cap_hi.push(lp) } else { cap_lo.push(lp) }
                continue;
            }
            // The rounding band: cap corners, cap edges, the edges up (or down) the band at each corner.
            let zc = if top { z + d } else { z - d };
            let up = if top { n } else { n.scale(-1.0) };
            let corner_v: Vec<usize> = src.iter().map(|e| b.vertex(e.seg.start(), zc)).collect();
            let cap_e: Vec<usize> = (0..m).map(|j| b.edge(prism.curve(&src[j].seg, zc), corner_v[j], corner_v[(j + 1) % m])).collect();
            // Per corner: (edge from the end of the previous element's offset, edge from the start of this one's).
            let mut climbs: Vec<(usize, usize)> = Vec::with_capacity(m);
            for j in 0..m {
                let prev = (j + m - 1) % m;
                let q = src[j].seg.start();
                let centre = prism.point(q, z);
                let quarter = |h: P2| {
                    let h3 = prism.vector(h);
                    Curve::Arc { c: centre, n: h3.cross(up), x: h3, r: d, t0: 0.0, t1: FRAC_PI_2 }
                };
                let k_out = ring.seg_at[j];
                match (ring.corner_at[j], ring.corners[j]) {
                    (Some(ka), _) => {
                        let e_in = b.edge(quarter(right(src[prev].seg.tangent_end())), v_mid[ka], corner_v[j]);
                        let e_out = b.edge(quarter(right(src[j].seg.tangent_start())), v_mid[k_out], corner_v[j]);
                        climbs.push((e_in, e_out));
                    }
                    (None, Corner::Tangent) => {
                        let e = b.edge(quarter(right(src[j].seg.tangent_start())), v_mid[k_out], corner_v[j]);
                        climbs.push((e, e));
                    }
                    (None, Corner::Sharp) => {
                        if !matches!((src[prev].seg, src[j].seg), (Seg2::Line { .. }, Seg2::Line { .. })) {
                            return unsupported(format!(
                                "{what}: loop {l} corner {j} is a sharp concave corner beside an arc; the rounding tubes there meet in a \
                                 curve that is not a conic (only line-line concave corners are supported)"));
                        }
                        let w = planar::sub(o[k_out].seg.start(), q);
                        let curve = Curve::Conic { c: centre, u: prism.vector(w), v: up.scale(d), t0: 0.0, t1: FRAC_PI_2 };
                        let e = b.edge(curve, v_mid[k_out], corner_v[j]);
                        climbs.push((e, e));
                    }
                    (None, c) => return Err(KernelError::Internal(format!("{what}: corner {j} offset as {c:?} without a rounding arc"))),
                }
            }
            for j in 0..m {
                let next = (j + 1) % m;
                let k = ring.seg_at[j];
                let (climb_start, climb_end) = (climbs[j].1, climbs[next].0);
                let centre = prism.point(src[j].seg.start(), z);
                let surface = match src[j].seg {
                    Seg2::Line { a, b: e } => {
                        let t = planar::sub(e, a);
                        Surface::Cylinder { p: centre, a: prism.vector(planar::scale(t, 1.0 / planar::len(t))), x: n, r: d }
                    }
                    Seg2::Arc { c, r, .. } => Surface::Torus { c: prism.point(c, z), a: n, x: prism.u, big: r, r: d },
                };
                let lp = if top {
                    vec![(mid_e[k], true), (climb_end, true), (cap_e[j], false), (climb_start, false)]
                } else {
                    vec![(mid_e[k], false), (climb_start, true), (cap_e[j], true), (climb_end, false)]
                };
                b.face(surface, true, vec![lp], Role::Pipe { top, l, j });
                if let Some(ka) = ring.corner_at[j] {
                    let sphere = Surface::Sphere { c: centre, a: n, x: prism.u, r: d };
                    let (e_in, e_out) = climbs[j];
                    let lp = if top {
                        vec![(mid_e[ka], true), (e_out, true), (e_in, false)]
                    } else {
                        vec![(mid_e[ka], false), (e_in, true), (e_out, false)]
                    };
                    b.face(sphere, true, vec![lp], Role::Ball { top, l, j });
                }
            }
            let lp: Vec<(usize, bool)> = if top { cap_e.iter().map(|&e| (e, true)).collect() } else { cap_e.iter().rev().map(|&e| (e, false)).collect() };
            if top { cap_hi.push(lp) } else { cap_lo.push(lp) }
        }
    }
    let zc_hi = if d > 0.0 { z_hi + d } else { z_hi };
    let zc_lo = if d > 0.0 { z_lo - d } else { z_lo };
    let anchor = base[0][0].seg.start();
    if !hi_open {
        b.face(Surface::Plane { p: prism.point(anchor, zc_hi), n, x: prism.u }, true, cap_hi, Role::Cap { top: true });
    }
    if !lo_open {
        b.face(Surface::Plane { p: prism.point(anchor, zc_lo), n, x: prism.u }, false, cap_lo, Role::Cap { top: false });
    }
    // The volume in closed form: the wall swept through the height, and per closed end the band
    // integral over h of area(A + sqrt(d^2 - h^2)) = |A| d + P pi d^2 / 4 + 2 K d^3 / 3, area(A + rho) = |A| + P rho + K rho^2.
    let mut volume = height * area_off;
    if d > 0.0 {
        let predicted = area_base + perimeter * d + k2 * d * d;
        if (predicted - area_off).abs() > 1e-9 * area_off.abs().max(1.0) {
            return Err(KernelError::Internal(format!("{what}: the offset profile's area {area_off} is not the predicted {predicted}")));
        }
        let closed = usize::from(!hi_open) + usize::from(!lo_open);
        volume += closed as f64 * (area_base * d + perimeter * std::f64::consts::PI * d * d / 4.0 + 2.0 * k2 * d * d * d / 3.0);
    }
    Ok(Layer { body: b.body, roles: b.roles, open_hi: out_hi, open_lo: out_lo, volume })
}

/// `inner` turned inside out (a void) and joined to `outer`, with a rim per loop at each open end.
pub fn join(prism: &Prism, outer: Layer, inner: Layer, z_top: f64, z_bottom: f64) -> (Body, Vec<Role>, f64) {
    let volume = outer.volume - inner.volume;
    let (mut body, mut roles) = (outer.body, outer.roles);
    let eo = body.edges.len();
    super::append(&mut body, &super::flipped(&inner.body));
    roles.extend(inner.roles);
    let anchor = prism.loops[0][0].seg.start();
    for (top, outer_loops, inner_loops) in [(true, &outer.open_hi, &inner.open_hi), (false, &outer.open_lo, &inner.open_lo)] {
        for (l, (ol, il)) in outer_loops.iter().zip(inner_loops.iter()).enumerate() {
            let (o_lp, i_lp): (Vec<(usize, bool)>, Vec<(usize, bool)>) = if top {
                (ol.iter().map(|&e| (e, true)).collect(), il.iter().rev().map(|&e| (e + eo, false)).collect())
            } else {
                (ol.iter().rev().map(|&e| (e, false)).collect(), il.iter().map(|&e| (e + eo, true)).collect())
            };
            let loops = if l == 0 { vec![o_lp, i_lp] } else { vec![i_lp, o_lp] };
            let z = if top { z_top } else { z_bottom };
            body.faces.push(Face { surface: Surface::Plane { p: prism.point(anchor, z), n: prism.n, x: prism.u }, sense: top, loops });
            roles.push(Role::Rim { top });
        }
    }
    (body, roles, volume)
}

/// The prism offset by `distance` (positive outward), rounded or sharp.
pub fn offset(prism: &Prism, distance: f64, sharp: bool, tol: f64) -> Result<(Body, Vec<Role>, f64), KernelError> {
    let len = prism.length;
    let layer = if distance > 0.0 && !sharp {
        layer(prism, &prism.loops, distance, -len, 0.0, false, false, tol, "the outward offset")?
    } else {
        let e = -distance;
        if distance < 0.0 && !(len - 2.0 * e > tol) {
            return unsupported(format!("an inward offset of {e} consumes the prism's thickness {len}"));
        }
        let rings = offset_profile(&prism.loops, distance, !sharp, tol, "the offset profile")?;
        let profile: Vec<Vec<Elem>> = rings.into_iter().map(|r| r.elems).collect();
        let (lo, hi) = (-len - distance, distance);
        layer(prism, &profile, 0.0, lo, hi, false, false, tol, "the offset")?
    };
    Ok((layer.body, layer.roles, layer.volume))
}

/// The prism shelled open at the removed caps: wall thickness t inward (t > 0) or outward (t < 0).
pub fn shell_open(prism: &Prism, t: f64, top_open: bool, bottom_open: bool, tol: f64) -> Result<(Body, Vec<Role>, f64), KernelError> {
    let len = prism.length;
    let (outer, inner) = if t > 0.0 {
        let outer = layer(prism, &prism.loops, 0.0, -len, 0.0, bottom_open, top_open, tol, "the shell")?;
        let z_lo = if bottom_open { -len } else { -len + t };
        let z_hi = if top_open { 0.0 } else { -t };
        if !(z_hi - z_lo > tol) { return unsupported(format!("a wall of {t} consumes the prism's thickness {len}")); }
        let rings = offset_profile(&prism.loops, -t, true, tol, "the shell's inner wall")?;
        let profile: Vec<Vec<Elem>> = rings.into_iter().map(|r| r.elems).collect();
        (outer, layer(prism, &profile, 0.0, z_lo, z_hi, bottom_open, top_open, tol, "the shell's cavity")?)
    } else {
        let outer = layer(prism, &prism.loops, -t, -len, 0.0, bottom_open, top_open, tol, "the shell's outer wall")?;
        (outer, layer(prism, &prism.loops, 0.0, -len, 0.0, bottom_open, top_open, tol, "the shell's cavity")?)
    };
    Ok(join(prism, outer, inner, 0.0, -len))
}

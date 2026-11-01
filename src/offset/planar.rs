//! The planar engine behind offsets and shells (Claude 2026-11-01): closed loops of lines and circular arcs in a plane,
//! material on the left of travel (outer loops counter-clockwise, holes clockwise), offset by a signed distance per
//! element (positive away from the material).
//!
//! Exact construction from plane geometry: a line moves parallel to itself, an arc keeps its centre and changes its
//! radius by the distance (growing on a convex arc, shrinking on a concave one). At a corner the moved neighbours either
//! meet where their supports cross (the corner opens or closes sharp), or, when the corner opens and the offset is
//! rounded, are joined by an arc of the distance's radius about the old corner point: the set of points at that distance
//! from the corner (the Minkowski sum with a disc). Corners whose turn is within tolerance of straight are tangent: the
//! moved ends agree within tolerance and simply join.
//!
//! Nothing is repaired. An element that would vanish or turn back (a concave arc no larger than the distance, a short
//! side eaten by its neighbours), neighbours whose moved supports do not meet, and loops that touch or cross after the
//! move (the offset changed the topology) are refused by name.

use crate::KernelError;
use std::f64::consts::{PI, TAU};

pub type P2 = (f64, f64);

pub fn add(a: P2, b: P2) -> P2 { (a.0 + b.0, a.1 + b.1) }
pub fn sub(a: P2, b: P2) -> P2 { (a.0 - b.0, a.1 - b.1) }
pub fn scale(a: P2, k: f64) -> P2 { (a.0 * k, a.1 * k) }
pub fn dot(a: P2, b: P2) -> f64 { a.0 * b.0 + a.1 * b.1 }
pub fn cross(a: P2, b: P2) -> f64 { a.0 * b.1 - a.1 * b.0 }
pub fn len(a: P2) -> f64 { a.0.hypot(a.1) }
pub fn dist(a: P2, b: P2) -> f64 { len(sub(a, b)) }
pub fn dir(t: f64) -> P2 { (t.cos(), t.sin()) }
/// The unit normal to the right of travel along unit `t` (away from the material).
pub fn right(t: P2) -> P2 { (t.1, -t.0) }
fn unit(a: P2) -> P2 { let l = len(a); (a.0 / l, a.1 / l) }
/// An angle difference wrapped into (-pi, pi].
pub fn wrap(t: f64) -> f64 { let w = (t + PI).rem_euclid(TAU) - PI; if w <= -PI { w + TAU } else { w } }

/// A line from a to b, or an arc about c of radius r from angle a0 through the signed sweep sw (counter-clockwise
/// positive, 0 < |sw| <= 2 pi).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Seg2 {
    Line { a: P2, b: P2 },
    Arc { c: P2, r: f64, a0: f64, sw: f64 },
}

impl Seg2 {
    pub fn start(&self) -> P2 {
        match *self { Seg2::Line { a, .. } => a, Seg2::Arc { c, r, a0, .. } => add(c, scale(dir(a0), r)) }
    }
    pub fn end(&self) -> P2 {
        match *self { Seg2::Line { b, .. } => b, Seg2::Arc { c, r, a0, sw } => add(c, scale(dir(a0 + sw), r)) }
    }
    /// The point at fraction f of the element.
    pub fn at(&self, f: f64) -> P2 {
        match *self {
            Seg2::Line { a, b } => add(a, scale(sub(b, a), f)),
            Seg2::Arc { c, r, a0, sw } => add(c, scale(dir(a0 + sw * f), r)),
        }
    }
    pub fn mid(&self) -> P2 { self.at(0.5) }
    /// Unit tangent along travel at the start and at the end.
    pub fn tangent_start(&self) -> P2 {
        match *self {
            Seg2::Line { a, b } => unit(sub(b, a)),
            Seg2::Arc { a0, sw, .. } => scale((-a0.sin(), a0.cos()), sw.signum()),
        }
    }
    pub fn tangent_end(&self) -> P2 {
        match *self {
            Seg2::Line { a, b } => unit(sub(b, a)),
            Seg2::Arc { a0, sw, .. } => { let t = a0 + sw; scale((-t.sin(), t.cos()), sw.signum()) }
        }
    }
    pub fn reversed(&self) -> Seg2 {
        match *self {
            Seg2::Line { a, b } => Seg2::Line { a: b, b: a },
            Seg2::Arc { c, r, a0, sw } => Seg2::Arc { c, r, a0: a0 + sw, sw: -sw },
        }
    }
    pub fn length(&self) -> f64 {
        match *self { Seg2::Line { a, b } => dist(a, b), Seg2::Arc { r, sw, .. } => r * sw.abs() }
    }
    /// The loop integral of (x dy - y dx) along the element: twice its share of the enclosed signed area.
    pub fn area2(&self) -> f64 {
        match *self {
            Seg2::Line { a, b } => cross(a, b),
            Seg2::Arc { c, r, a0, sw } => {
                let a1 = a0 + sw;
                r * c.0 * (a1.sin() - a0.sin()) - r * c.1 * (a1.cos() - a0.cos()) + r * r * sw
            }
        }
    }
    /// The same element mirrored in the line x = 0 (a meridian's reflection in the axis), run the same way.
    pub fn mirrored(&self) -> Seg2 {
        match *self {
            Seg2::Line { a, b } => Seg2::Line { a: (-a.0, a.1), b: (-b.0, b.1) },
            Seg2::Arc { c, r, a0, sw } => Seg2::Arc { c: (-c.0, c.1), r, a0: PI - a0, sw: -sw },
        }
    }
}

/// Where an element came from, in the source's own terms: segment j of source loop l, the corner at the start of that
/// segment, or a piece of the axis (a meridian's closing segment).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Tag { Seg(usize, usize), Corner(usize, usize), Axis }

/// A loop element with its origin. `corner` names the source corner at the element's start (the tag a rounding arc
/// made there takes); `mirrored` marks the reflected half of a meridian; `removed` an element on a removed face.
#[derive(Clone, Copy, Debug)]
pub struct Elem {
    pub seg: Seg2,
    pub tag: Tag,
    pub corner: Tag,
    pub mirrored: bool,
    pub removed: bool,
}

/// How a corner was offset.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Corner {
    /// Both neighbours unmoved (distance 0).
    Still,
    /// The turn is within tolerance of straight: the moved ends join.
    Tangent,
    /// A rounding arc about the corner point.
    Round,
    /// The moved supports meet at a point.
    Sharp,
}

/// An offset loop: its elements in order, and for each source element and corner where it went.
#[derive(Clone, Debug)]
pub struct OffsetLoop {
    pub elems: Vec<Elem>,
    /// Index in `elems` of each source element.
    pub seg_at: Vec<usize>,
    /// Index in `elems` of the rounding arc at each source corner (the start of source element j), if any.
    pub corner_at: Vec<Option<usize>>,
    pub corners: Vec<Corner>,
    /// The signed turn at each source corner (counter-clockwise positive: a convex corner for material on the left).
    pub turns: Vec<f64>,
}

fn unsupported<T>(m: String) -> Result<T, KernelError> { Err(KernelError::Unsupported(m)) }

/// A moved element's support: a line through a point along a unit direction, or a circle.
#[derive(Clone, Copy, Debug)]
enum Support { Line { p: P2, d: P2 }, Circle { c: P2, r: f64 } }

fn moved(seg: &Seg2, delta: f64, tol: f64, name: &str) -> Result<(Seg2, Support), KernelError> {
    match *seg {
        Seg2::Line { a, b } => {
            let t = unit(sub(b, a));
            let k = scale(right(t), delta);
            let (a, b) = (add(a, k), add(b, k));
            Ok((Seg2::Line { a, b }, Support::Line { p: a, d: t }))
        }
        Seg2::Arc { c, r, a0, sw } => {
            let r2 = r + delta * sw.signum();
            if !(r2 > tol) {
                return unsupported(format!("{name}: an arc of radius {r} would shrink to radius {r2} (the offset removes it)"));
            }
            Ok((Seg2::Arc { c, r: r2, a0, sw }, Support::Circle { c, r: r2 }))
        }
    }
}

/// Points where two supports cross (tangent contact within tolerance counts once).
fn crossings(s: Support, t: Support, tol: f64) -> Vec<P2> {
    match (s, t) {
        (Support::Line { p, d }, Support::Line { p: q, d: e }) => {
            let den = cross(d, e);
            if den.abs() <= 1e-12 { return Vec::new(); }
            let s = cross(sub(q, p), e) / den;
            vec![add(p, scale(d, s))]
        }
        (Support::Line { p, d }, Support::Circle { c, r }) | (Support::Circle { c, r }, Support::Line { p, d }) => {
            let w = sub(p, c);
            let b = dot(d, w);
            let disc = b * b - (dot(w, w) - r * r);
            if disc < -2.0 * r * tol { return Vec::new(); }
            let h = disc.max(0.0).sqrt();
            if h == 0.0 { vec![add(p, scale(d, -b))] } else { vec![add(p, scale(d, -b - h)), add(p, scale(d, -b + h))] }
        }
        (Support::Circle { c, r }, Support::Circle { c: c2, r: r2 }) => {
            let dd = dist(c, c2);
            if dd <= 1e-12 * (1.0 + r) { return Vec::new(); }
            let u = scale(sub(c2, c), 1.0 / dd);
            let a = (r * r - r2 * r2 + dd * dd) / (2.0 * dd);
            let h2 = r * r - a * a;
            if h2 < -2.0 * r.min(r2) * tol { return Vec::new(); }
            let h = h2.max(0.0).sqrt();
            let base = add(c, scale(u, a));
            let n = (-u.1, u.0);
            if h == 0.0 { vec![base] } else { vec![add(base, scale(n, h)), add(base, scale(n, -h))] }
        }
    }
}

/// The element re-cut to run from `s` to `e` along its support, refused when it would vanish or turn back.
fn trimmed(seg: &Seg2, s: P2, e: P2, tol: f64, name: &str) -> Result<Seg2, KernelError> {
    match *seg {
        Seg2::Line { a, b } => {
            if dot(sub(e, s), sub(b, a)) <= 0.0 || dist(s, e) <= tol {
                return unsupported(format!("{name}: the offset makes this side vanish (its neighbours' offsets meet across it)"));
            }
            Ok(Seg2::Line { a: s, b: e })
        }
        Seg2::Arc { c, r, a0, sw } => {
            let ds = wrap((s.1 - c.1).atan2(s.0 - c.0) - a0);
            let de = wrap((e.1 - c.1).atan2(e.0 - c.0) - (a0 + sw));
            let sw2 = sw + de - ds;
            if sw2.signum() != sw.signum() || sw2.abs() * r <= tol || sw2.abs() > TAU + 1e-9 {
                return unsupported(format!("{name}: the offset makes this arc vanish (its neighbours' offsets meet across it)"));
            }
            Ok(Seg2::Arc { c, r, a0: a0 + ds, sw: sw2 })
        }
    }
}

/// Offsets one loop. `delta(elem)` is each element's distance (positive away from the material); `rounded` rounds the
/// corners that open (exact offset), otherwise every corner is the meeting of the moved supports (sharp). `name` and
/// the corner tags identify the loop in messages and in the result.
pub fn offset_loop(lp: &[Elem], delta: &dyn Fn(&Elem) -> f64, rounded: bool, tol: f64, name: &str)
    -> Result<OffsetLoop, KernelError> {
    let m = lp.len();
    let mut moved_segs = Vec::with_capacity(m);
    for (i, el) in lp.iter().enumerate() {
        moved_segs.push(moved(&el.seg, delta(el), tol, &format!("{name} element {i}"))?);
    }
    // Per corner j (the start of element j): the end of element j - 1 and the start of element j after the move.
    let mut ends: Vec<(P2, P2)> = Vec::with_capacity(m);
    let mut arcs: Vec<Option<Seg2>> = Vec::with_capacity(m);
    let mut corners = Vec::with_capacity(m);
    let mut turns = Vec::with_capacity(m);
    for j in 0..m {
        let prev = (j + m - 1) % m;
        let (t_in, t_out) = (lp[prev].seg.tangent_end(), lp[j].seg.tangent_start());
        let turn = cross(t_in, t_out).atan2(dot(t_in, t_out));
        turns.push(turn);
        let (d_in, d_out) = (delta(&lp[prev]), delta(&lp[j]));
        let (e_in, s_out) = (moved_segs[prev].0.end(), moved_segs[j].0.start());
        let v = lp[j].seg.start();
        let meet = |e_in: P2, s_out: P2| scale(add(e_in, s_out), 0.5);
        if d_in == d_out && d_in == 0.0 {
            ends.push((v, v));
            arcs.push(None);
            corners.push(Corner::Still);
            continue;
        }
        if d_in == d_out && (d_in * turn).abs() <= tol && dist(e_in, s_out) <= tol {
            let q = meet(e_in, s_out);
            ends.push((q, q));
            arcs.push(None);
            corners.push(Corner::Tangent);
            continue;
        }
        if d_in == d_out && rounded && d_in * turn > 0.0 {
            let d = d_in;
            let start = scale(right(t_in), d.signum());
            arcs.push(Some(Seg2::Arc { c: v, r: d.abs(), a0: start.1.atan2(start.0), sw: turn }));
            ends.push((e_in, s_out));
            corners.push(Corner::Round);
            continue;
        }
        let near = meet(e_in, s_out);
        let q = crossings(moved_segs[prev].1, moved_segs[j].1, tol).into_iter()
            .min_by(|a, b| dist(*a, near).total_cmp(&dist(*b, near)));
        let Some(q) = q else {
            return unsupported(format!("{name} corner {j}: the offset faces either side do not meet"));
        };
        ends.push((q, q));
        arcs.push(None);
        corners.push(Corner::Sharp);
    }
    let mut out = OffsetLoop { elems: Vec::new(), seg_at: Vec::new(), corner_at: Vec::new(), corners, turns };
    for j in 0..m {
        if let Some(arc) = arcs[j] {
            out.corner_at.push(Some(out.elems.len()));
            out.elems.push(Elem { seg: arc, tag: lp[j].corner, corner: lp[j].corner, mirrored: lp[j].mirrored, removed: false });
        } else {
            out.corner_at.push(None);
        }
        let (s, e) = (ends[j].1, ends[(j + 1) % m].0);
        let seg = trimmed(&moved_segs[j].0, s, e, tol, &format!("{name} element {j}"))?;
        out.seg_at.push(out.elems.len());
        out.elems.push(Elem { seg, ..lp[j] });
    }
    Ok(out)
}

/// Twice a loop's signed area.
pub fn area2(lp: &[Elem]) -> f64 { lp.iter().map(|e| e.seg.area2()).sum() }

/// Is angle t on the arc from a0 through sw (with a small angular slack)?
fn on_arc(a0: f64, sw: f64, t: f64, slack: f64) -> bool {
    let rel = if sw > 0.0 { (t - a0).rem_euclid(TAU) } else { (a0 - t).rem_euclid(TAU) };
    rel <= sw.abs() + slack || rel >= TAU - slack
}

fn support(seg: &Seg2) -> Support {
    match *seg {
        Seg2::Line { a, b } => Support::Line { p: a, d: unit(sub(b, a)) },
        Seg2::Arc { c, r, .. } => Support::Circle { c, r },
    }
}

fn on_piece(seg: &Seg2, q: P2, tol: f64) -> bool {
    match *seg {
        Seg2::Line { a, b } => {
            let d = sub(b, a);
            let s = dot(sub(q, a), d) / dot(d, d);
            s >= -tol / len(d) && s <= 1.0 + tol / len(d)
        }
        Seg2::Arc { c, r, a0, sw } => on_arc(a0, sw, (q.1 - c.1).atan2(q.0 - c.0), tol / r),
    }
}

/// Distance from a point to an element (exact).
pub fn point_distance(q: P2, seg: &Seg2) -> f64 {
    match *seg {
        Seg2::Line { a, b } => {
            let d = sub(b, a);
            let s = (dot(sub(q, a), d) / dot(d, d)).clamp(0.0, 1.0);
            dist(q, add(a, scale(d, s)))
        }
        Seg2::Arc { c, r, a0, sw } => {
            let w = sub(q, c);
            if len(w) > 0.0 && on_arc(a0, sw, w.1.atan2(w.0), 0.0) { (len(w) - r).abs() }
            else { dist(q, seg.start()).min(dist(q, seg.end())) }
        }
    }
}

/// The least distance between two elements that do not cross: at an end of one, or where both are square to the line
/// joining them (a line's foot from a circle's centre; two arcs along their line of centres).
fn piece_distance(p: &Seg2, q: &Seg2) -> f64 {
    let mut best = [point_distance(p.start(), q), point_distance(p.end(), q), point_distance(q.start(), p), point_distance(q.end(), p)]
        .into_iter().fold(f64::INFINITY, f64::min);
    let mut pair = |x: P2, y: P2| best = best.min(dist(x, y));
    match (*p, *q) {
        (Seg2::Line { a, b }, Seg2::Arc { c, r, a0, sw }) | (Seg2::Arc { c, r, a0, sw }, Seg2::Line { a, b }) => {
            let d = sub(b, a);
            let s = dot(sub(c, a), d) / dot(d, d);
            if s > 0.0 && s < 1.0 {
                let f = add(a, scale(d, s));
                let n = unit((-d.1, d.0));
                for k in [1.0, -1.0] {
                    let t = (n.1 * k).atan2(n.0 * k);
                    if on_arc(a0, sw, t, 0.0) { pair(f, add(c, scale(n, k * r))); }
                }
            }
        }
        (Seg2::Arc { c, r, a0, sw }, Seg2::Arc { c: c2, r: r2, a0: b0, sw: sw2 }) => {
            let dd = dist(c, c2);
            if dd <= 1e-12 * (r + r2) {
                // Concentric: |r - r2| wherever their angles overlap.
                if [a0, a0 + sw].iter().any(|&t| on_arc(b0, sw2, t, 0.0)) || [b0, b0 + sw2].iter().any(|&t| on_arc(a0, sw, t, 0.0)) {
                    best = best.min((r - r2).abs());
                }
            } else {
                let u = scale(sub(c2, c), 1.0 / dd);
                for k1 in [1.0, -1.0] {
                    for k2 in [1.0, -1.0] {
                        let (t1, t2) = ((u.1 * k1).atan2(u.0 * k1), (u.1 * k2).atan2(u.0 * k2));
                        if on_arc(a0, sw, t1, 0.0) && on_arc(b0, sw2, t2, 0.0) {
                            pair(add(c, scale(u, k1 * r)), add(c2, scale(u, k2 * r2)));
                        }
                    }
                }
            }
        }
        _ => {}
    }
    best
}

/// Every element keeps clear of every other except its neighbours, which meet only at their shared end (without folding
/// back): the loops are simple and apart. `what` names the loops in the refusal.
pub fn check_simple(loops: &[Vec<Elem>], tol: f64, what: &str) -> Result<(), KernelError> {
    let all: Vec<(usize, usize, Seg2)> = loops.iter().enumerate()
        .flat_map(|(l, lp)| lp.iter().enumerate().map(move |(i, e)| (l, i, e.seg))).collect();
    for u in 0..all.len() {
        for w in u + 1..all.len() {
            let ((l1, i1, p), (l2, i2, q)) = (all[u], all[w]);
            let n = loops[l1].len();
            let after = l1 == l2 && (i1 + 1) % n == i2;
            let before = l1 == l2 && (i2 + 1) % n == i1;
            let mut shared = Vec::new();
            if after { shared.push(p.end()); }
            if before { shared.push(p.start()); }
            let fail = || unsupported(format!(
                "{what}: elements {i1} of loop {l1} and {i2} of loop {l2} touch or cross (the offset changes the shape's topology)"));
            let cross_pts = crossings(support(&p), support(&q), tol);
            for x in cross_pts {
                if on_piece(&p, x, tol) && on_piece(&q, x, tol) && shared.iter().all(|s| dist(*s, x) > 4.0 * tol) {
                    return fail();
                }
            }
            if after || before {
                // Folding back onto a neighbour: one's middle lies on the other.
                if point_distance(p.mid(), &q) <= tol || point_distance(q.mid(), &p) <= tol { return fail(); }
                continue;
            }
            if piece_distance(&p, &q) <= tol { return fail(); }
        }
    }
    Ok(())
}

/// Is q inside the loop (odd crossings of the ray from q towards +x)? Exact: lines by the half-open rule in y, arcs
/// split at their highest and lowest points into pieces monotone in y, each crossed at most once.
pub fn inside(q: P2, lp: &[Elem]) -> bool {
    let mut odd = false;
    let mut cross_monotone = |s: P2, e: P2, x_at: &dyn Fn(f64) -> f64| {
        if (s.1 > q.1) != (e.1 > q.1) && x_at(q.1) > q.0 { odd = !odd; }
    };
    for el in lp {
        match el.seg {
            Seg2::Line { a, b } => cross_monotone(a, b, &|y| a.0 + (b.0 - a.0) * (y - a.1) / (b.1 - a.1)),
            Seg2::Arc { c, r, a0, sw } => {
                // Angles where the arc turns in y: pi/2 + k pi inside (a0, a0 + sw).
                let (lo, hi) = if sw > 0.0 { (a0, a0 + sw) } else { (a0 + sw, a0) };
                let mut cuts = vec![lo];
                let mut k = ((lo - PI / 2.0) / PI).floor() + 1.0;
                while PI / 2.0 + k * PI < hi { cuts.push(PI / 2.0 + k * PI); k += 1.0; }
                cuts.push(hi);
                for w in cuts.windows(2) {
                    let (t0, t1) = (w[0], w[1]);
                    let side = ((t0 + t1) / 2.0).cos().signum();
                    let at = |t: f64| add(c, scale(dir(t), r));
                    cross_monotone(at(t0), at(t1), &|y| c.0 + side * (r * r - (y - c.1).powi(2)).max(0.0).sqrt());
                }
            }
        }
    }
    odd
}

/// The offset loops nest as the source loops did (a hole that grew past its outer loop without touching it, or two
/// outer loops that swallowed one another, changed the topology).
pub fn check_nesting(before: &[Vec<Elem>], after: &[Vec<Elem>], what: &str) -> Result<(), KernelError> {
    for i in 0..before.len() {
        for j in 0..before.len() {
            if i == j { continue; }
            if inside(before[i][0].seg.mid(), &before[j]) != inside(after[i][0].seg.mid(), &after[j]) {
                return unsupported(format!("{what}: loops {i} and {j} no longer nest as they did (the offset changes the shape's topology)"));
            }
        }
    }
    Ok(())
}

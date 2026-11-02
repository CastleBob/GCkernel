//! Radial bores through a cylindrical wall at any clock angle (Claude, 2026-10-29; Bob's Pipe Drill Test: holes turned
//! 340 and 320 degrees about the pipe were refused, the transverse domains taking only coordinate-axis cylinders).
//!
//! A full round cutter whose axis is perpendicular to a pipe's axis (any angle about it, any offset that keeps the
//! hole clear of tangency) and which starts outside the pipe and ends either in its bore (one wall) or beyond the far
//! wall (both walls). Each pierced wall gets an inner loop on its outer and inner cylinder faces, the exact
//! intersection of the cutter's cylinder with that wall (a CylinderPair edge, wire type 4), and the cutter's own
//! cylinder between the two loops, joined by one straight seam. Nothing else in the body changes: no partition, no
//! seams elsewhere. Admission is decided from the geometry before anything is built, and refuses explicitly when the
//! hole would meet a wall's other boundaries (ends, seam, other holes) or any other face, or when a wall is split.
//! Clearance tests use rigorous boxes: arcs and conics by their extrema, pair edges by the perpendicular-cylinder
//! bound below; an inconclusive box refuses. Independent implementation from elementary geometry.
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;
use std::f64::consts::TAU;

type Result<T> = std::result::Result<T, KernelError>;

fn unsupported(message: &str) -> KernelError {
    KernelError::Unsupported(message.into())
}

fn xyz(v: Vec3) -> [f64; 3] { [v.x, v.y, v.z] }

/// The cutter: start cap centre `s` (outside the pipe), unit direction `d` into the pipe, radius, length.
struct Cutter {
    s: Vec3,
    d: Vec3,
    rho: f64,
    len: f64,
}

/// The pipe: axis point `p`, unit axis `b`, the outer wall face (radius `ro`) and inner wall face (radius `ri`).
struct Pipe {
    p: Vec3,
    b: Vec3,
    outer: usize,
    ro: f64,
    inner: usize,
    ri: f64,
}

fn perp(v: Vec3, b: Vec3) -> Vec3 { v.sub(b.scale(v.dot(b))) }

/// A full round cutter: one cylinder face, two plane caps normal to its axis, two full-circle cap edges on the axis.
pub(super) fn cutter(tool: &Body, tol: f64) -> Option<(Vec3, Vec3, f64)> {
    if tool.faces.len() != 3 || tool.edges.len() != 3 || tool.vertices.len() != 2 {
        return None;
    }
    let (a, r) = tool.faces.iter().find_map(|f| match f.surface {
        Surface::Cylinder { a, r, .. } => Some((a, r)),
        _ => None,
    })?;
    let caps: Vec<Vec3> = tool
        .edges
        .iter()
        .filter_map(|e| match e.curve {
            Curve::Arc { c, n, r: ar, t0, t1, .. }
                if (ar - r).abs() <= tol && (t1 - t0 - TAU).abs() <= 1e-9 && n.cross(a).length() <= 1e-9 =>
            {
                Some(c)
            }
            _ => None,
        })
        .collect();
    if caps.len() != 2 || caps[1].sub(caps[0]).cross(a).length() > tol {
        return None;
    }
    Some((caps[0], caps[1], r))
}

/// The pipe the cutter crosses radially, with the cutter oriented from outside in; None unless every radial-bore
/// condition holds (selection is decided here, before anything is built).
fn admit(body: &Body, tool: &Body, tol: f64) -> Option<(Cutter, Pipe)> {
    let (c0, c1, rho) = cutter(tool, tol)?;
    let axis = c1.sub(c0).unit().ok()?;
    let mut found: Option<(Cutter, Pipe)> = None;
    for (i, f) in body.faces.iter().enumerate() {
        let Surface::Cylinder { p, a: b, r: ro, .. } = f.surface else { continue };
        if !f.sense || b.dot(axis).abs() > 1e-9 {
            continue;
        }
        // Its inner wall: a coaxial cylinder of smaller radius facing the axis.
        let inner: Vec<(usize, f64)> = body
            .faces
            .iter()
            .enumerate()
            .filter_map(|(j, g)| match g.surface {
                Surface::Cylinder { p: q, a: c, r, .. }
                    if !g.sense && r < ro - tol && c.cross(b).length() <= 1e-9 && perp(q.sub(p), b).length() <= tol =>
                {
                    Some((j, r))
                }
                _ => None,
            })
            .collect();
        let [(j, ri)] = inner[..] else { continue };
        // Outside in: the cap farther from the pipe's axis is the start.
        let dist = |q: Vec3| perp(q.sub(p), b).length();
        let (s, e) = if dist(c0) >= dist(c1) { (c0, c1) } else { (c1, c0) };
        let d = e.sub(s).unit().ok()?;
        // A cap disk spans b and w = d x b; a point of it is (q - p).d along d and offset + rho sin(phi) along w, so
        // its distance from the pipe axis is at least |(q - p).d| and at most sqrt(((q - p).d)^2 + (offset + rho)^2).
        let offset = s.sub(p).dot(d.cross(b)).abs();
        let (sd, ed) = (s.sub(p).dot(d), e.sub(p).dot(d));
        let starts_outside = sd < 0. && -sd > ro + tol;
        let in_bore = (ed * ed + (offset + rho).powi(2)).sqrt() < ri - tol;
        let past_far_wall = ed > ro + tol;
        if rho + offset >= ri - tol || !starts_outside || !(in_bore || past_far_wall) {
            continue;
        }
        if found.is_some() {
            return None; // two pipes would qualify: ambiguous, not chosen
        }
        found = Some((
            Cutter { s, d, rho, len: e.sub(s).length() },
            Pipe { p, b, outer: i, ro, inner: j, ri },
        ));
    }
    found
}

/// Selected only when every radial-bore condition holds, clearances included: anything else keeps the other domains.
pub(super) fn selects(body: &Body, tool: &Body, tol: f64) -> bool {
    admit(body, tool, tol).is_some_and(|(k, pipe)| {
        let body = reseamed(body, &k, &pipe, tol);
        clear(&body, tool, &k, &pipe, tol).is_ok()
    })
}

/// The angular interval a perpendicular round hole occupies on a wall of radius r about the pipe axis (frame e1,
/// e2 = b x e1): the hole's centre direction is branch * a (the smaller root h lies toward -a), and every point of
/// the loop is within |offset| + rho of that generator sideways, so within asin((|offset| + rho) / r) of its angle.
/// None when that does not bound it (the hole not perpendicular, or as wide as the wall).
#[allow(clippy::too_many_arguments)]
fn hole_interval(c: Vec3, a: Vec3, rho: f64, branch: f64, pipe: &Pipe, r: f64, e1: Vec3, e2: Vec3) -> Option<(f64, f64)> {
    if a.dot(pipe.b).abs() > 1e-9 {
        return None;
    }
    let offset = c.sub(pipe.p).dot(pipe.b.cross(a)).abs();
    let reach = (offset + rho) / r;
    if !(reach < 1.) {
        return None;
    }
    let dir = a.scale(branch);
    let centre = dir.dot(e2).atan2(dir.dot(e1));
    let half = reach.asin();
    Some((centre - half, centre + half))
}

/// Bob's Pipe 1.25 12 Cross Holes: a side hole at 3:00 lands on the wall's seam, the generator where the
/// periodic face's parameter starts (the profile circles start at +X). The seam is not geometry: the same face
/// with the same rings and the seam on another generator is the same solid. So when the new hole's loops reach a
/// plain wall's seam (one full ring at each end, a straight seam used both ways, nothing else at its two vertices),
/// the seam moves to the middle of the widest angular gap between that wall's holes, old and new; the rings are
/// the same circles started at the new generator, and the wall surface's parameter starts there too. The surfaces' geometry, radii
/// and every other edge are unchanged. Decided from the geometry before the cut; when no clear generator exists the
/// wall keeps its seam and the cut is refused by name as before.
fn reseamed(body: &Body, k: &Cutter, pipe: &Pipe, tol: f64) -> Body {
    let mut out = body.clone();
    let b = pipe.b;
    let end = k.s.add(k.d.scale(k.len));
    let both = perp(end.sub(pipe.p), b).length() > pipe.ri;
    for (face, r) in [(pipe.outer, pipe.ro), (pipe.inner, pipe.ri)] {
        let Surface::Cylinder { x: e1, .. } = body.faces[face].surface else { continue };
        let Ok(e1) = perp(e1, b).unit() else { continue };
        let e2 = b.cross(e1);
        // The new hole's intervals on this wall: entry side (branch -1) and, through both walls, the exit side.
        let mut busy = vec![];
        for branch in if both { vec![-1., 1.] } else { vec![-1.] } {
            match hole_interval(k.s, k.d, k.rho, branch, pipe, r, e1, e2) {
                Some(i) => busy.push(i),
                None => return out,
            }
        }
        // Its plain outer loop: two full rings and one seam used both ways.
        let f = &body.faces[face];
        let Some(outer) = f.loops.iter().position(|l| l.len() == 4) else { continue };
        let lp = &f.loops[outer];
        let lines: Vec<usize> =
            lp.iter().filter(|(e, _)| matches!(body.edges[*e].curve, Curve::Line { .. })).map(|(e, _)| *e).collect();
        let rings: Vec<usize> =
            lp.iter().filter(|(e, _)| matches!(body.edges[*e].curve, Curve::Arc { .. })).map(|(e, _)| *e).collect();
        if lines.len() != 2 || lines[0] != lines[1] || rings.len() != 2 || rings[0] == rings[1] {
            continue;
        }
        let seam = lines[0];
        let ring_ok = |e: usize| match body.edges[e].curve {
            Curve::Arc { c, n, r: ar, t0, t1, .. } => {
                (t1 - t0 - TAU).abs() <= 1e-9
                    && n.cross(b).length() <= 1e-9
                    && (ar - r).abs() <= tol
                    && perp(c.sub(pipe.p), b).length() <= tol
                    && body.edges[e].start == body.edges[e].end
            }
            _ => false,
        };
        if !rings.iter().all(|e| ring_ok(*e)) {
            continue;
        }
        let ends = [body.edges[seam].start, body.edges[seam].end];
        let ring_vertices = [body.edges[rings[0]].start, body.edges[rings[1]].start];
        if !(ends.contains(&ring_vertices[0]) && ends.contains(&ring_vertices[1]) && ends[0] != ends[1]) {
            continue;
        }
        // Nothing else may end at the seam's vertices.
        if body.edges.iter().enumerate().any(|(i, e)| {
            i != seam && !rings.contains(&i) && (ends.contains(&e.start) || ends.contains(&e.end))
        }) {
            continue;
        }
        // Does the new hole reach the seam? (Its angle within the hole's interval widened by tol / r.)
        let Curve::Line { a: sa, .. } = body.edges[seam].curve else { continue };
        let q = perp(sa.sub(pipe.p), b);
        let seam_angle = q.dot(e2).atan2(q.dot(e1));
        let margin = tol / r;
        let inside = |(lo, hi): (f64, f64), t: f64| {
            let d = (t - (lo + hi) * 0.5 + std::f64::consts::PI).rem_euclid(TAU) - std::f64::consts::PI;
            d.abs() <= (hi - lo) * 0.5 + margin
        };
        if !busy.iter().any(|i| inside(*i, seam_angle)) {
            continue;
        }
        // The wall's other holes.
        let mut known = true;
        for (i, l) in f.loops.iter().enumerate() {
            if i == outer {
                continue;
            }
            match (l.as_slice(), &body.edges[l[0].0].curve) {
                ([_], Curve::CylinderPair { c, a, rho, branch, t0, t1, .. }) if (t1 - t0 - TAU).abs() <= 1e-9 => {
                    match hole_interval(*c, *a, *rho, *branch, pipe, r, e1, e2) {
                        Some(i) => busy.push(i),
                        None => known = false,
                    }
                }
                _ => known = false,
            }
        }
        if !known {
            continue;
        }
        // The widest gap between the busy intervals, round the circle.
        let mut spans: Vec<(f64, f64)> = busy
            .iter()
            .map(|(lo, hi)| {
                let l = lo.rem_euclid(TAU);
                (l, l + (hi - lo))
            })
            .collect();
        spans.sort_by(|x, y| x.0.total_cmp(&y.0));
        // Each interval once more a turn later, so a span crossing angle 2 pi covers the start of the next turn.
        let first = spans[0].0;
        let turned: Vec<(f64, f64)> = spans.iter().map(|(l, h)| (l + TAU, h + TAU)).collect();
        spans.extend(turned);
        let mut best: Option<(f64, f64)> = None;
        let mut reach = spans[0].1;
        for &(next, top) in &spans[1..] {
            if next > first + TAU {
                break;
            }
            if next > reach && best.is_none_or(|(w, _)| next - reach > w) {
                best = Some((next - reach, (next + reach) * 0.5));
            }
            reach = reach.max(top);
        }
        let Some((width, angle)) = best else { continue };
        if r * (width * 0.5).sin() <= 4. * tol {
            continue;
        }
        let e = e1.scale(angle.cos()).add(e2.scale(angle.sin()));
        for ring in &rings {
            let v = out.edges[*ring].start;
            if let Curve::Arc { c, n, x, t0, .. } = &mut out.edges[*ring].curve {
                let start = c.add(e.scale(r));
                *x = e.scale(t0.cos()).sub(n.cross(e).scale(t0.sin()));
                out.vertices[v] = start;
            }
        }
        let (va, vb) = (out.edges[seam].start, out.edges[seam].end);
        out.edges[seam].curve = Curve::Line { a: out.vertices[va], b: out.vertices[vb] };
        // The wall's own parameter starts on the seam too (readers that place a periodic face's seam at its
        // parameter origin, as GCmaker's display import does, then see the same face).
        if let Surface::Cylinder { x, .. } = &mut out.faces[face].surface {
            *x = e;
        }
    }
    out
}

/// A rigorous box of a curve: arcs and conics by their extrema, a pair edge whose cylinders are perpendicular by the
/// bound below, anything else none (the caller refuses).
fn bounds(curve: &Curve) -> Option<([f64; 3], [f64; 3])> {
    match *curve {
        Curve::CylinderPair { c, a, rho, p, b, r, branch, .. } => {
            if a.dot(b).abs() > 1e-9 {
                return None;
            }
            // Q - c is perpendicular to a and its parts along b and w = b x a are at most rho; the point is on the
            // root cylinder, so its part along a is branch sqrt(r^2 - w^2) with w its (bounded) part along w.
            let w = b.cross(a);
            let (cb, cw) = (c.sub(p).dot(b), c.sub(p).dot(w));
            let (w0, w1) = (cw - rho, cw + rho);
            let m = w0.abs().max(w1.abs());
            if m >= r {
                return None;
            }
            let near = (r * r - m * m).sqrt();
            let along = if branch > 0. { [near, r] } else { [-r, -near] };
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for sb in [cb - rho, cb + rho] {
                for sw in [w0, w1] {
                    for sa in along {
                        let q = xyz(p.add(b.scale(sb)).add(w.scale(sw)).add(a.scale(sa)));
                        for k in 0..3 {
                            lo[k] = lo[k].min(q[k]);
                            hi[k] = hi[k].max(q[k]);
                        }
                    }
                }
            }
            Some((lo, hi))
        }
        _ => Some(super::passage::curve_bounds(curve)),
    }
}

/// A RadialPair edge (a normalized pipe end's rim) in 128 pieces, each with a rigorous box (the piece of its circle of
/// radius s by its extrema, plus the pair's proved height range along a: geom::pair_height_range) and its angular
/// interval in the frame (e1, e2) across the axis. None for any other curve, or when a height range is not proved.
fn radial_pieces(curve: &Curve, e1: Vec3, e2: Vec3) -> Option<Vec<(([f64; 3], [f64; 3]), (f64, f64))>> {
    let Curve::RadialPair { c, a, x, rho, p, b, r, branch, t0, t1, s } = *curve else { return None };
    let n = 128;
    let y = a.cross(x);
    let angle = |t: f64| { let u = x.scale(t.cos()).add(y.scale(t.sin())); u.dot(e2).atan2(u.dot(e1)) };
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let (ta, tb) = (t0 + (t1 - t0) * i as f64 / n as f64, t0 + (t1 - t0) * (i + 1) as f64 / n as f64);
        let (mut lo, mut hi) = super::passage::curve_bounds(&Curve::Arc { c, n: a, x, r: s, t0: ta, t1: tb });
        let (hl, hh) = crate::geom::pair_height_range(c, a, x, rho, p, b, r, branch, ta, tb, 16)?;
        for (k, ak) in [a.x, a.y, a.z].into_iter().enumerate() {
            lo[k] += (ak * hl).min(ak * hh);
            hi[k] += (ak * hl).max(ak * hh);
        }
        let (aa, ab) = (angle(ta), angle(tb));
        let ab = aa + (ab - aa + std::f64::consts::PI).rem_euclid(TAU) - std::f64::consts::PI;
        out.push(((lo, hi), (aa.min(ab), aa.max(ab))));
    }
    Some(out)
}

/// Whether two angular intervals (lo, hi) on the circle meet, within `margin`.
fn arcs_meet(i: (f64, f64), j: (f64, f64), margin: f64) -> bool {
    let d = ((i.0 + i.1) * 0.5 - (j.0 + j.1) * 0.5 + std::f64::consts::PI).rem_euclid(TAU) - std::f64::consts::PI;
    d.abs() <= (i.1 - i.0) * 0.5 + (j.1 - j.0) * 0.5 + margin
}

fn apart(a: ([f64; 3], [f64; 3]), b: ([f64; 3], [f64; 3]), gap: f64) -> bool {
    (0..3).any(|k| a.1[k] + gap < b.0[k] || b.1[k] + gap < a.0[k])
}

fn union(boxes: impl Iterator<Item = ([f64; 3], [f64; 3])>) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for (l, h) in boxes {
        for k in 0..3 {
            lo[k] = lo[k].min(l[k]);
            hi[k] = hi[k].max(h[k]);
        }
    }
    (lo, hi)
}

/// The walls crossed (outside in: face, radius, branch) and their loop curves, after every clearance check: each wall
/// one face per radius with its boundary clear of the loops (ends, seam, earlier holes) and the loops strictly between
/// its two end rings along the axis; every other face clear of the whole cutter.
fn clear(body: &Body, tool: &Body, k: &Cutter, pipe: &Pipe, tol: f64) -> Result<(Vec<(usize, f64, f64)>, Vec<Curve>)> {
    let b = pipe.b;
    // The walls crossed, outside in: (face, radius, branch). The cutter's axis meets a wall at h = -beta +- root.
    let mut walls = vec![(pipe.outer, pipe.ro, -1.), (pipe.inner, pipe.ri, -1.)];
    let end = k.s.add(k.d.scale(k.len));
    if perp(end.sub(pipe.p), b).length() > pipe.ri {
        walls.extend([(pipe.inner, pipe.ri, 1.), (pipe.outer, pipe.ro, 1.)]);
    }
    let loop_curve = |r: f64, branch: f64| Curve::CylinderPair {
        c: k.s, a: k.d, x: b, rho: k.rho, p: pipe.p, b, r, branch, t0: 0., t1: TAU,
    };
    let loops: Vec<Curve> = walls.iter().map(|&(_, r, branch)| loop_curve(r, branch)).collect();
    let loop_box = union(loops.iter().map(|c| bounds(c).expect("perpendicular by admission")));
    let along_box = |bx: ([f64; 3], [f64; 3])| {
        let along = |q: [f64; 3]| (q[0] - pipe.p.x) * b.x + (q[1] - pipe.p.y) * b.y + (q[2] - pipe.p.z) * b.z;
        let corners = (0..8).map(|m| [if m & 1 == 0 { bx.0[0] } else { bx.1[0] }, if m & 2 == 0 { bx.0[1] } else { bx.1[1] }, if m & 4 == 0 { bx.0[2] } else { bx.1[2] }]);
        let values: Vec<f64> = corners.map(along).collect();
        (values.iter().cloned().fold(f64::INFINITY, f64::min), values.iter().cloned().fold(f64::NEG_INFINITY, f64::max))
    };
    for &(face, r) in [(pipe.outer, pipe.ro), (pipe.inner, pipe.ri)].iter() {
        let f = &body.faces[face];
        let edges: Vec<usize> = f.loops.iter().flatten().map(|(e, _)| *e).collect();
        let mut rings: Vec<(f64, f64)> = vec![];
        // A normalized end's rim (Claude 2026-10-31) is bounded piece by piece: every piece clear of the loops, and the
        // pieces at the hole's own angles (its intervals on this wall) taken together as that end's place along the axis.
        let Surface::Cylinder { x: wx, .. } = f.surface else { unreachable!() };
        let e1 = perp(wx, b).unit()?;
        let e2 = b.cross(e1);
        let busy: Vec<(f64, f64)> = walls.iter().filter(|w| w.0 == face)
            .filter_map(|&(_, _, branch)| hole_interval(k.s, k.d, k.rho, branch, pipe, r, e1, e2)).collect();
        for &e in &edges {
            if let Some(pieces) = radial_pieces(&body.edges[e].curve, e1, e2) {
                let mut span: Option<(f64, f64)> = None;
                for (bx, angles) in pieces {
                    if !apart(bx, loop_box, tol) {
                        return Err(unsupported("radial bore meets the normalized end of the pipe wall"));
                    }
                    if busy.iter().any(|i| arcs_meet(*i, angles, tol / r)) {
                        let (l, h) = along_box(bx);
                        span = Some(span.map_or((l, h), |(sl, sh)| (sl.min(l), sh.max(h))));
                    }
                }
                rings.extend(span);
                continue;
            }
            let bx = bounds(&body.edges[e].curve).ok_or_else(|| unsupported("radial bore: a wall edge has no bounding box"))?;
            if !apart(bx, loop_box, tol) {
                return Err(unsupported("radial bore meets a boundary of the pipe wall (an end, its seam or another hole)"));
            }
            if !matches!(body.edges[e].curve, Curve::Line { .. }) {
                rings.push(along_box(bx));
            }
        }
        let corners = (0..8).map(|m| [if m & 1 == 0 { loop_box.0[0] } else { loop_box.1[0] }, if m & 2 == 0 { loop_box.0[1] } else { loop_box.1[1] }, if m & 4 == 0 { loop_box.0[2] } else { loop_box.1[2] }]);
        let along: Vec<f64> = corners.map(|q| (q[0] - pipe.p.x) * b.x + (q[1] - pipe.p.y) * b.y + (q[2] - pipe.p.z) * b.z).collect();
        let (lo, hi) = (along.iter().cloned().fold(f64::INFINITY, f64::min), along.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
        if !rings.iter().any(|r| r.1 < lo - tol) || !rings.iter().any(|r| r.0 > hi + tol) {
            return Err(unsupported("radial bore is not between the pipe wall's two ends"));
        }
    }
    // Every other face stays clear of the whole cutter.
    let cutter_box = union(tool.edges.iter().map(|e| super::passage::curve_bounds(&e.curve)));
    for (i, f) in body.faces.iter().enumerate() {
        if i == pipe.outer || i == pipe.inner {
            continue;
        }
        for (e, _) in f.loops.iter().flatten() {
            let boxes = match radial_pieces(&body.edges[*e].curve, b.perpendicular(), b.cross(b.perpendicular())) {
                Some(pieces) => pieces.into_iter().map(|(bx, _)| bx).collect(),
                None => vec![bounds(&body.edges[*e].curve).ok_or_else(|| unsupported("radial bore: an edge has no bounding box"))?],
            };
            if boxes.into_iter().any(|bx| !apart(bx, cutter_box, tol)) {
                return Err(unsupported("radial bore meets another face of the body"));
            }
        }
    }
    Ok((walls, loops))
}

/// Every planar face's loops enclose positive area about its outward normal (a face turned inside out, or a loop run
/// the wrong way, is malformed input, not a cut to make). The topology check alone does not see this.
pub(super) fn planar_faces_outward(body: &Body) -> Result<()> {
    for (i, f) in body.faces.iter().enumerate() {
        if let Surface::Plane { n, .. } = f.surface {
            let outward = n.scale(if f.sense { 1. } else { -1. });
            let area: f64 = f
                .loops
                .iter()
                .map(|lp| {
                    let curves: Vec<Curve> = lp
                        .iter()
                        .map(|(e, fwd)| if *fwd { body.edges[*e].curve.clone() } else { body.edges[*e].curve.reversed() })
                        .collect();
                    crate::extrude::signed_area(&curves, outward)
                })
                .sum();
            if !(area > 0.) {
                return Err(KernelError::InvalidInput(format!("planar face {i} is inside out")));
            }
        }
    }
    Ok(())
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    cut_traced(body, tool, tol).map(|(result, _)| result)
}

/// The radial bore and its history (trace::bored): the body kept, each wall crossed gaining its meeting loop, each wall
/// passage a seam and the face from the cutter's cylinder. No history when the bore moved a wall's seam (the seam,
/// rings and their vertices are then the same names on moved geometry; refused rather than called unchanged).
pub(super) fn cut_traced(body: &Body, tool: &Body, tol: f64) -> Result<(Body, super::History)> {
    crate::topology::check(body, tol)?;
    crate::topology::check(tool, tol)?;
    planar_faces_outward(body)?;
    planar_faces_outward(tool)?;
    let (k, pipe) = admit(body, tool, tol).ok_or_else(|| unsupported("not a radial bore through a pipe wall"))?;
    let original = body;
    let body = &reseamed(body, &k, &pipe, tol);
    let moved = crate::neutral::write_body(body) != crate::neutral::write_body(original);
    crate::topology::check(body, tol)?;
    let (walls, loops) = clear(body, tool, &k, &pipe, tol)?;
    let b = pipe.b;
    // Build: per wall pair (entry and exit of one wall), two loop edges, one seam, one cutter-cylinder face.
    let mut out = body.clone();
    let mut hole_flags = vec![];
    for (w, &(face, _, _)) in walls.iter().enumerate() {
        let curve = loops[w].clone();
        let (q, tangent) = curve.eval(0.);
        let v = out.vertices.len();
        out.vertices.push(q);
        let e = out.edges.len();
        out.edges.push(Edge { curve, start: v, end: v });
        // A hole runs clockwise about the wall's outward normal: the hole (toward the cutter's axis) on the right.
        let radial = perp(q.sub(pipe.p), b).unit()?;
        let normal = radial.scale(if out.faces[face].sense { 1. } else { -1. });
        let foot = k.s.add(k.d.scale(q.sub(k.s).dot(k.d)));
        let inward = perp(foot.sub(q), normal);
        let forward = tangent.cross(normal).dot(inward) > 0.;
        out.faces[face].loops.push(vec![(e, forward)]);
        hole_flags.push((e, forward, v));
    }
    for pair in hole_flags.chunks_exact(2) {
        let ((e_in, f_in, v_in), (e_out, f_out, v_out)) = (pair[0], pair[1]);
        let seam = out.edges.len();
        out.edges.push(Edge {
            curve: Curve::Line { a: out.vertices[v_in], b: out.vertices[v_out] },
            start: v_in,
            end: v_out,
        });
        // The cutter's wall faces into the hole; it runs each loop the other way from its wall.
        out.faces.push(Face {
            surface: Surface::Cylinder { p: k.s, a: k.d, x: b, r: k.rho },
            sense: false,
            loops: vec![vec![(e_in, !f_in), (seam, true), (e_out, !f_out), (seam, false)]],
        });
    }
    crate::topology::check(&out, tol)?;
    let history = if moved {
        Err("cut history is not recorded for a radial bore that moves the pipe wall's seam yet".to_string())
    } else {
        match tool.faces.iter().position(|f| matches!(f.surface, Surface::Cylinder { .. })) {
            Some(j) => super::trace::bored(body, &walls.iter().map(|w| w.0).collect::<Vec<_>>(), j),
            None => Err("cut history is not recorded for this cut: the cutter has no cylinder face".to_string()),
        }
    };
    Ok((out, history))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mass::measure;
    use crate::neutral::Profile;

    const TOL: f64 = 1e-6;

    fn circle(r: f64) -> Curve {
        Curve::Arc { c: Vec3::new(0., 0., 0.), n: Vec3::new(0., 0., 1.), x: Vec3::new(1., 0., 0.), r, t0: 0., t1: TAU }
    }
    /// A pipe along Y (outer radius ro, inner ri, length 10), as the shape factory places it (centred on the Y axis).
    fn pipe(ro: f64, ri: f64) -> Body {
        let lp = |r: f64| vec![Curve::Arc { c: Vec3::new(0., 0., 0.), n: Vec3::new(0., -1., 0.), x: Vec3::new(1., 0., 0.), r, t0: 0., t1: TAU }];
        crate::extrude::extrude(&Profile { loops: vec![lp(ro), lp(ri)], n: Vec3::new(0., 1., 0.), low: 0., high: 10. }, TOL).unwrap()
    }
    /// A round cutter of radius a at axial position y, turned `deg` about the pipe axis from +Z, from height `top`
    /// above the axis down to `bottom` (negative = past the axis).
    fn bore(a: f64, y: f64, deg: f64, top: f64, bottom: f64) -> Body {
        let t = deg.to_radians();
        let mut body = crate::extrude::extrude(&Profile { loops: vec![vec![circle(a)]], n: Vec3::new(0., 0., 1.), low: bottom, high: top }, TOL).unwrap();
        // Turn the Z-axis cutter about Y by deg and move it to y.
        let turn = |v: Vec3| Vec3::new(v.x * t.cos() + v.z * t.sin(), v.y, -v.x * t.sin() + v.z * t.cos());
        let place = |v: Vec3| turn(v).add(Vec3::new(0., y, 0.));
        for v in &mut body.vertices {
            *v = place(*v);
        }
        for e in &mut body.edges {
            e.curve = match e.curve {
                Curve::Line { a, b } => Curve::Line { a: place(a), b: place(b) },
                Curve::Arc { c, n, x, r, t0, t1 } => Curve::Arc { c: place(c), n: turn(n), x: turn(x), r, t0, t1 },
                _ => unreachable!(),
            };
        }
        for f in &mut body.faces {
            f.surface = match f.surface {
                Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => unreachable!("no elliptic-cylinder or revolved faces in these test bodies"),
                Surface::Plane { p, n, x } => Surface::Plane { p: place(p), n: turn(n), x: turn(x) },
                Surface::Cylinder { p, a, x, r } => Surface::Cylinder { p: place(p), a: turn(a), x: turn(x), r },
            };
        }
        body
    }
    /// Removed through one wall: integral of 2 sqrt(a^2 - x^2) (sqrt(ro^2 - x^2) - sqrt(ri^2 - x^2)) dx over (-a, a).
    fn one_wall(a: f64, ro: f64, ri: f64) -> f64 {
        let n = 400_000;
        (0..n)
            .map(|i| {
                let x = -a + (i as f64 + 0.5) * 2. * a / n as f64;
                2. * (a * a - x * x).sqrt() * ((ro * ro - x * x).sqrt() - (ri * ri - x * x).sqrt()) * 2. * a / n as f64
            })
            .sum()
    }

    #[test]
    fn bobs_turned_holes_cut_one_wall_each_at_any_angle() {
        // Bob's Pipe Drill Test: OD 3.5, wall .216, 3/8 holes at 0, 340 and 320 degrees, one wall deep.
        let (ro, ri, a) = (1.75, 1.75 - 0.216, 0.1875);
        let mut body = pipe(ro, ri);
        let start = measure(&body, TOL).unwrap().volume;
        // As GCmaker places them: from the surface plus the 0.001 boolean overshoot, 1.0 deep.
        for (k, (y, deg)) in [(4., 0.), (6., 340.), (8., 320.)].into_iter().enumerate() {
            let tool = bore(a, y, deg, ro + 0.001, ro - 1.0);
            assert!(selects(&body, &tool, TOL));
            body = crate::boolean::cut(&body, &tool, TOL).unwrap();
            let removed = start - measure(&body, TOL).unwrap().volume;
            assert!((removed - (k + 1) as f64 * one_wall(a, ro, ri)).abs() < 1e-8, "after hole {}: {removed}", k + 1);
        }
    }

    #[test]
    fn through_both_walls_removes_two_walls() {
        let (ro, ri, a) = (1.75, 1.5, 0.25);
        let body = pipe(ro, ri);
        let start = measure(&body, TOL).unwrap().volume;
        let result = crate::boolean::cut(&body, &bore(a, 5., 37., 2.5, -2.5), TOL).unwrap();
        let removed = start - measure(&result, TOL).unwrap().volume;
        assert!((removed - 2. * one_wall(a, ro, ri)).abs() < 1e-8, "{removed}");
    }

    #[test]
    fn a_hole_on_the_seam_moves_the_seam_to_the_widest_clear_gap() {
        // The profile circles start at +X (90 degrees from +Z): a hole there, one wall and through both walls, and
        // with holes already on both sides of the seam; each removes its exact wall integrals, the seam moved.
        let (ro, ri, a) = (1.75, 1.5, 0.25);
        let body = pipe(ro, ri);
        let start = measure(&body, TOL).unwrap().volume;
        let one = crate::boolean::cut(&body, &bore(a, 5., 90., 2.5, 0.5), TOL).unwrap();
        assert!((start - measure(&one, TOL).unwrap().volume - one_wall(a, ro, ri)).abs() < 1e-8);
        let both = crate::boolean::cut(&body, &bore(a, 5., 90., 2.5, -2.5), TOL).unwrap();
        assert!((start - measure(&both, TOL).unwrap().volume - 2. * one_wall(a, ro, ri)).abs() < 1e-8);
        let mut holed = body.clone();
        for deg in [60., 120., 0.] {
            holed = crate::boolean::cut(&holed, &bore(a, 5., deg, 2.5, 0.5), TOL).unwrap();
        }
        let next = crate::boolean::cut(&holed, &bore(a, 5., 90., 2.5, 0.5), TOL).unwrap();
        crate::topology::check(&next, TOL).unwrap();
        assert!((start - measure(&next, TOL).unwrap().volume - 4. * one_wall(a, ro, ri)).abs() < 1e-8);
        // Overlapping an existing hole is still refused, seam or not.
        assert!(matches!(crate::boolean::cut(&one, &bore(a, 5.2, 90., 2.5, 0.5), TOL), Err(KernelError::Unsupported(_))));
    }

    #[test]
    fn holes_meeting_ends_seams_or_each_other_and_wall_pockets_are_refused() {
        let (ro, ri, a) = (1.75, 1.5, 0.25);
        let body = pipe(ro, ri);
        // Ending inside the wall is a pocket, not this domain.
        assert!(!selects(&body, &bore(a, 5., 10., 2.5, 1.6), TOL));
        // Across the pipe's end.
        assert!(matches!(crate::boolean::cut(&body, &bore(a, 0.1, 10., 2.5, 0.5), TOL), Err(KernelError::Unsupported(_))));
        // A second hole overlapping the first.
        let once = crate::boolean::cut(&body, &bore(a, 5., 10., 2.5, 0.5), TOL).unwrap();
        assert!(matches!(crate::boolean::cut(&once, &bore(a, 5.2, 10., 2.5, 0.5), TOL), Err(KernelError::Unsupported(_))));
    }
}

//! Offsets and shells of full solids of revolution (Claude 2026-11-01). The distance from a point to a solid turned
//! about an axis is the distance, in the point's own half-plane through the axis, to the solid's meridian section (the
//! nearest point of a rotationally symmetric set lies at the same angle). So every offset or shell of such a solid is
//! the solid of revolution of the planar offset of its section, taken in the whole plane through the axis (the section
//! and its mirror image, so corners on the axis offset correctly) and then cut back to the half-plane r >= 0.
//!
//! The section's lines and arcs move to lines and arcs: a plane across the axis stays a plane, a cylinder stays a
//! coaxial cylinder, a cone a coaxial cone, a sphere a concentric sphere, a torus a torus with its tube radius changed;
//! a rounded corner (a convex edge circle under an outward offset) becomes a torus about that circle, or a sphere about
//! an apex on the axis. The kernel's revolve (revolve.rs) builds the result exactly from the offset section.
//!
//! A shell opened through removed faces is the section less its cavity: the section offset by the wall everywhere
//! except along the removed faces' meridian segments, the cavity's walls reversed and joined to the rims the removed
//! faces leave.

use super::planar::{self, area2, check_simple, dist, offset_loop, Elem, Seg2, Tag, P2};
use super::Source;
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::Body;
use crate::revolve::{self, Kind, Revolution, Segment};
use crate::KernelError;
use std::f64::consts::TAU;

fn unsupported<T>(m: String) -> Result<T, KernelError> { Err(KernelError::Unsupported(m)) }

/// Where a section corner came from: the body's circle edges there, or its vertex on the axis.
#[derive(Clone, Debug)]
pub enum CornerSource { Parallels(Vec<usize>), AxisVertex(Option<usize>) }

/// A full solid of revolution found in a body.
pub struct Revolved {
    pub p: Vec3,
    pub a: Vec3,
    pub x: Vec3,
    /// Section loops in (r, z) (r from the axis, z along it from p), material on the left; element j of loop l is
    /// Seg(l, j) or, closing the loop along the axis, Axis.
    pub loops: Vec<Vec<Elem>>,
    /// Each element's source faces (two when the body splits the surface in halves; none on the axis).
    pub faces: Vec<Vec<Vec<usize>>>,
    /// The source of the corner at each element's start.
    pub corners: Vec<Vec<CornerSource>>,
}

impl Revolved {
    fn rz(&self, q: Vec3) -> P2 {
        let w = q.sub(self.p);
        let h = w.dot(self.a);
        (w.sub(self.a.scale(h)).length(), h)
    }
    fn radial(&self, q: Vec3) -> Vec3 { let w = q.sub(self.p); w.sub(self.a.scale(w.dot(self.a))) }
}

enum RawKind { Line, Arc { c: P2, r: f64, hint: P2 } }

struct Raw {
    ends: [P2; 2],
    kind: RawKind,
    faces: Vec<usize>,
}

/// The unit radial direction (in the plane across the axis) at which a face mostly lies: the mean of its edges' middles
/// and vertices' directions, for telling the halves of a split surface apart.
fn heading(body: &Body, face: usize, rev: &Revolved) -> Option<Vec3> {
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    for &(e, _) in body.faces[face].loops.iter().flatten() {
        let c = &body.edges[e].curve;
        let (t0, t1) = c.domain();
        for q in [c.eval((t0 + t1) / 2.0).0, c.start()] {
            let w = rev.radial(q);
            if w.length() > 1e-9 * (1.0 + q.length()) { sum = sum.add(w.scale(1.0 / w.length())); }
        }
    }
    (sum.length() > 1e-9).then(|| sum.scale(1.0 / sum.length()))
}

fn pick(candidates: &[usize], want: Option<Vec3>, of: &dyn Fn(usize) -> Option<Vec3>) -> usize {
    if candidates.len() == 1 { return candidates[0]; }
    let Some(w) = want else { return candidates[0] };
    *candidates.iter().max_by(|&&i, &&j| {
        let s = |k: usize| of(k).map_or(-2.0, |h| h.dot(w));
        s(i).total_cmp(&s(j))
    }).unwrap()
}

/// The body as a full solid of revolution of lines and arcs, or why it is not one.
pub fn recognize(body: &Body, tol: f64) -> Result<Revolved, KernelError> {
    let axis = body.faces.iter().find_map(|f| match f.surface {
        Surface::Cylinder { p, a, x, .. } | Surface::Cone { p, a, x, .. } | Surface::Torus { c: p, a, x, .. } => Some((p, a, x)),
        _ => None,
    }).or_else(|| body.faces.iter().find_map(|f| match f.surface { Surface::Sphere { c, a, x, .. } => Some((c, a, x)), _ => None }));
    let Some((p, a, x)) = axis else { return unsupported("not a solid of revolution: no face turns about an axis".into()) };
    let mut rev = Revolved { p, a, x, loops: Vec::new(), faces: Vec::new(), corners: Vec::new() };
    let on_axis = |q: Vec3| { let w = q.sub(p); w.sub(a.scale(w.dot(a))).length() <= tol };
    let parallel = |c: &Curve| match *c {
        Curve::Arc { c, n, .. } => n.cross(a).length() < 1e-9 && on_axis(c),
        _ => false,
    };
    let mut raws: Vec<Raw> = Vec::new();
    for (f, face) in body.faces.iter().enumerate() {
        let coaxial = |b: Vec3, q: Vec3| b.cross(a).length() < 1e-9 && on_axis(q);
        let ok = match face.surface {
            Surface::Plane { n, .. } => n.cross(a).length() < 1e-9,
            Surface::Cylinder { p: q, a: b, .. } | Surface::Cone { p: q, a: b, .. } | Surface::Torus { c: q, a: b, .. } => coaxial(b, q),
            Surface::Sphere { c, .. } => on_axis(c),
            Surface::Revolution { .. } => return unsupported(format!(
                "not offset here: face {f} is a general surface of revolution (an ellipse or Bezier meridian), whose offset is not one")),
            Surface::EllipticCylinder { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => return unsupported(format!("not a solid of revolution: face {f} is not turned about an axis")),
        };
        if !ok { return unsupported(format!("not a solid of revolution: face {f} does not turn about the body's axis")); }
        let mut levels: Vec<P2> = Vec::new();
        let mut hint = None;
        let add = |q: P2, levels: &mut Vec<P2>| if !levels.iter().any(|l| dist(*l, q) <= tol) { levels.push(q) };
        for &(e, _) in face.loops.iter().flatten() {
            let c = &body.edges[e].curve;
            if parallel(c) {
                let Curve::Arc { c: centre, r, .. } = *c else { unreachable!() };
                add((r, centre.sub(p).dot(a)), &mut levels);
            } else {
                for q in [c.start(), c.end()] { add(rev.rz(q), &mut levels); }
                let (t0, t1) = c.domain();
                hint = Some(rev.rz(c.eval((t0 + t1) / 2.0).0));
            }
        }
        if let Surface::Plane { .. } = face.surface {
            if levels.len() == 1 && levels[0].0 > tol { let z = levels[0].1; levels.push((0.0, z)); }
        }
        let [l0, l1] = levels[..] else {
            return unsupported(format!("not a solid of revolution: face {f} spans {} meridian levels, not two (a part turn or a cut face)", levels.len()));
        };
        let kind = match face.surface {
            Surface::Plane { .. } | Surface::Cylinder { .. } | Surface::Cone { .. } => RawKind::Line,
            Surface::Sphere { c, r, .. } => RawKind::Arc { c: (0.0, c.sub(p).dot(a)), r, hint: hint.unwrap_or((r, c.sub(p).dot(a))) },
            Surface::Torus { c, big, r, .. } => {
                let Some(h) = hint else { return unsupported(format!("not a solid of revolution: torus face {f} has no meridian edge")) };
                RawKind::Arc { c: (big, c.sub(p).dot(a)), r, hint: h }
            }
            _ => unreachable!(),
        };
        if let RawKind::Arc { c, r, .. } = kind {
            if (dist(l0, c) - r).abs() > tol || (dist(l1, c) - r).abs() > tol {
                return unsupported(format!("not a solid of revolution: face {f}'s boundary circles are not on its meridian"));
            }
        }
        if let Some(raw) = raws.iter_mut().find(|w| {
            let same_ends = (dist(w.ends[0], l0) <= tol && dist(w.ends[1], l1) <= tol) || (dist(w.ends[0], l1) <= tol && dist(w.ends[1], l0) <= tol);
            same_ends && match (&w.kind, &kind) {
                (RawKind::Line, RawKind::Line) => true,
                (RawKind::Arc { c, r, hint }, RawKind::Arc { c: c2, r: r2, hint: h2 }) => dist(*c, *c2) <= tol && (r - r2).abs() <= tol && {
                    // The same arc between these ends: the hints lie on the same side of the chord.
                    let chord = planar::sub(l1, l0);
                    planar::cross(chord, planar::sub(*hint, l0)).signum() == planar::cross(chord, planar::sub(*h2, l0)).signum()
                },
                _ => false,
            }
        }) {
            raw.faces.push(f);
            continue;
        }
        raws.push(Raw { ends: [l0, l1], kind, faces: vec![f] });
    }
    // Chain the meridian segments into loops; ends on the axis join along it.
    let mut points: Vec<P2> = Vec::new();
    let id = |q: P2, points: &mut Vec<P2>| -> usize {
        if let Some(i) = points.iter().position(|w| dist(*w, q) <= tol) { return i; }
        points.push(if q.0 <= tol { (0.0, q.1) } else { q });
        points.len() - 1
    };
    let ends: Vec<[usize; 2]> = raws.iter().map(|w| [id(w.ends[0], &mut points), id(w.ends[1], &mut points)]).collect();
    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); points.len()];
    for (i, e) in ends.iter().enumerate() {
        if e[0] == e[1] { return unsupported(format!("not a solid of revolution: face {} is a whole meridian circle", raws[i].faces[0])); }
        for &q in e { incident[q].push(i); }
    }
    // Axis points with one segment: paired from the top down, each pair a stretch of axis inside the section.
    let mut lonely: Vec<usize> = (0..points.len()).filter(|&q| points[q].0 == 0.0 && incident[q].len() == 1).collect();
    lonely.sort_by(|&i, &j| points[j].1.total_cmp(&points[i].1));
    if lonely.len() % 2 == 1 { return unsupported("not a solid of revolution: the section meets the axis an odd number of times".into()); }
    // (Axis segments are numbered after the raw segments.)
    let mut segs: Vec<[usize; 2]> = ends.clone();
    for pair in lonely.chunks(2) {
        let k = segs.len();
        incident[pair[0]].push(k);
        incident[pair[1]].push(k);
        segs.push([pair[0], pair[1]]);
    }
    if let Some(q) = (0..points.len()).find(|&q| incident[q].len() != 2) {
        return unsupported(format!("not a solid of revolution: {} meridian segments meet at ({}, {})", incident[q].len(), points[q].0, points[q].1));
    }
    let mut used = vec![false; segs.len()];
    let mut loops: Vec<Vec<(usize, bool)>> = Vec::new();
    for first in 0..segs.len() {
        if used[first] { continue; }
        let mut lp = Vec::new();
        let (mut k, mut at) = (first, segs[first][0]);
        loop {
            used[k] = true;
            let fwd = segs[k][0] == at;
            lp.push((k, fwd));
            at = if fwd { segs[k][1] } else { segs[k][0] };
            let Some(&next) = incident[at].iter().find(|&&j| j != k) else { break };
            if next == first { break; }
            if used[next] { return unsupported("not a solid of revolution: the meridian does not close into loops".into()); }
            k = next;
        }
        loops.push(lp);
    }
    // Oriented elements; the largest loop counter-clockwise, every other loop clockwise.
    let seg_of = |k: usize, fwd: bool| -> Seg2 {
        let [i, j] = if fwd { segs[k] } else { [segs[k][1], segs[k][0]] };
        let (s, e) = (points[i], points[j]);
        if k >= raws.len() { return Seg2::Line { a: s, b: e }; }
        match raws[k].kind {
            RawKind::Line => Seg2::Line { a: s, b: e },
            RawKind::Arc { c, r, hint } => {
                let ang = |q: P2| (q.1 - c.1).atan2(q.0 - c.0);
                let (a0, a1, am) = (ang(s), ang(e), ang(hint));
                let ccw = (a1 - a0).rem_euclid(TAU);
                let sw = if (am - a0).rem_euclid(TAU) < ccw { ccw } else { ccw - TAU };
                Seg2::Arc { c, r, a0, sw }
            }
        }
    };
    let mut oriented: Vec<Vec<(usize, Seg2)>> = loops.iter().map(|lp| lp.iter().map(|&(k, fwd)| (k, seg_of(k, fwd))).collect()).collect();
    let area = |lp: &[(usize, Seg2)]| lp.iter().map(|(_, s)| s.area2()).sum::<f64>();
    let outer = (0..oriented.len()).max_by(|&i, &j| area(&oriented[i]).abs().total_cmp(&area(&oriented[j]).abs())).unwrap();
    oriented.swap(0, outer);
    for (l, lp) in oriented.iter_mut().enumerate() {
        if (area(lp) > 0.0) != (l == 0) { *lp = lp.iter().rev().map(|&(k, s)| (k, s.reversed())).collect(); }
    }
    for (l, lp) in oriented.iter().enumerate() {
        let m = lp.len();
        let mut elems = Vec::with_capacity(m);
        let mut faces = Vec::with_capacity(m);
        let mut corners = Vec::with_capacity(m);
        for (j, &(k, seg)) in lp.iter().enumerate() {
            let axis_seg = k >= raws.len();
            elems.push(Elem { seg, tag: if axis_seg { Tag::Axis } else { Tag::Seg(l, j) }, corner: Tag::Corner(l, j), mirrored: false, removed: false });
            faces.push(if axis_seg { Vec::new() } else { raws[k].faces.clone() });
            let q = seg.start();
            corners.push(if q.0 == 0.0 {
                CornerSource::AxisVertex((0..body.vertices.len()).find(|&v| { let w = rev.rz(body.vertices[v]); w.0 <= tol && (w.1 - q.1).abs() <= tol }))
            } else {
                CornerSource::Parallels((0..body.edges.len()).filter(|&e| {
                    let c = &body.edges[e].curve;
                    parallel(c) && matches!(*c, Curve::Arc { c: centre, r, .. } if (r - q.0).abs() <= tol && (centre.sub(p).dot(a) - q.1).abs() <= tol)
                }).collect())
            });
        }
        rev.loops.push(elems);
        rev.faces.push(faces);
        rev.corners.push(corners);
    }
    // A regression check on the structure: the section turned a full turn has the body's volume and area.
    let revolution = revolution(&rev, &rev.loops, tol)?;
    let rebuilt = revolve::revolve_general(&revolution, tol)
        .map_err(|e| KernelError::Unsupported(format!("not a solid of revolution the kernel can rebuild: {e}")))?;
    let (before, after) = (crate::mass::measure(body, tol)?, crate::mass::measure(&rebuilt, tol)?);
    let agree = |x: f64, y: f64| (x - y).abs() <= 1e-9 * x.abs().max(y.abs());
    if !agree(before.volume, after.volume) || !agree(before.area, after.area) {
        return unsupported(format!("not a full solid of revolution of its section (volume {} vs {}, area {} vs {})",
            before.volume, after.volume, before.area, after.area));
    }
    Ok(rev)
}

/// The op 9 request turning half-plane loops a full turn about the solid's axis.
fn revolution(rev: &Revolved, half: &[Vec<Elem>], tol: f64) -> Result<Revolution, KernelError> {
    let snap = |r: f64| if r.abs() <= tol { 0.0 } else { r };
    let loops = half.iter().map(|lp| lp.iter().map(|e| {
        let (r, z) = e.seg.start();
        let kind = match e.seg {
            Seg2::Line { .. } => Kind::Line,
            Seg2::Arc { c, sw, .. } => Kind::Arc { cr: snap(c.0), cz: c.1, ccw: sw > 0.0 },
        };
        Segment { r: snap(r), z, kind }
    }).collect()).collect();
    Ok(Revolution { p: rev.p, a: rev.a, x: rev.x, angle: TAU, loops })
}

/// The section in the whole plane through the axis: each run of the half-plane loop between axis stretches, followed by
/// its mirror image run backwards (a loop clear of the axis keeps its own mirror image as a second loop).
pub fn whole_plane(half: &[Vec<Elem>]) -> Vec<Vec<Elem>> {
    let mirror = |l: usize, j: usize, e: &Elem, m: usize| Elem {
        seg: e.seg.mirrored().reversed(), tag: e.tag, corner: Tag::Corner(l, (j + 1) % m), mirrored: true, removed: e.removed,
    };
    let mut out = Vec::new();
    for (l, lp) in half.iter().enumerate() {
        let m = lp.len();
        let Some(axis) = lp.iter().position(|e| e.tag == Tag::Axis) else {
            out.push(lp.clone());
            out.push(lp.iter().enumerate().rev().map(|(j, e)| mirror(l, j, e, m)).collect());
            continue;
        };
        // Runs start just after an axis element.
        let mut run: Vec<(usize, Elem)> = Vec::new();
        for step in 1..=m {
            let j = (axis + step) % m;
            if lp[j].tag == Tag::Axis {
                if !run.is_empty() {
                    let mut whole: Vec<Elem> = run.iter().map(|&(_, e)| e).collect();
                    whole.extend(run.iter().rev().map(|&(i, e)| mirror(l, i, &e, m)));
                    out.push(whole);
                    run.clear();
                }
            } else {
                run.push((j, lp[j]));
            }
        }
    }
    out
}

/// Whole-plane loops cut back to r >= 0 and closed along the axis.
pub fn half_plane(whole: &[Vec<Elem>], tol: f64) -> Result<Vec<Vec<Elem>>, KernelError> {
    let snap = |q: P2| if q.0.abs() <= tol { (0.0, q.1) } else { q };
    let mut pieces: Vec<Elem> = Vec::new();
    for lp in whole {
        for e in lp {
            let mut cuts: Vec<f64> = Vec::new();
            match e.seg {
                Seg2::Line { a, b } => {
                    if (a.0 > tol && b.0 < -tol) || (a.0 < -tol && b.0 > tol) { cuts.push(a.0 / (a.0 - b.0)); }
                }
                Seg2::Arc { c, r, a0, sw } => {
                    if c.0.abs() < r {
                        let base = (-c.0 / r).acos();
                        for t in [base, -base] {
                            let f = (if sw > 0.0 { (t - a0).rem_euclid(TAU) } else { (a0 - t).rem_euclid(TAU) }) / sw.abs();
                            if f * e.seg.length() > tol && (1.0 - f) * e.seg.length() > tol && f < 1.0 { cuts.push(f); }
                        }
                    }
                }
            }
            cuts.sort_by(f64::total_cmp);
            let mut from = 0.0;
            for f in cuts.into_iter().chain([1.0]) {
                let seg = match e.seg {
                    Seg2::Line { a, b } => {
                        let s = if from == 0.0 { a } else { (0.0, e.seg.at(from).1) };
                        let t = if f == 1.0 { b } else { (0.0, e.seg.at(f).1) };
                        Seg2::Line { a: snap(s), b: snap(t) }
                    }
                    Seg2::Arc { c, r, a0, sw } => Seg2::Arc { c, r, a0: a0 + sw * from, sw: sw * (f - from) },
                };
                from = f;
                if seg.mid().0 > tol { pieces.push(Elem { seg, ..*e }); }
            }
        }
    }
    // Axis stretches: from each arrival (a piece ending on the axis) down to the next departure below it.
    let at_axis = |q: P2| q.0.abs() <= tol;
    let mut stops: Vec<(f64, bool)> = Vec::new(); // (z, arrival)
    for e in &pieces {
        if at_axis(e.seg.end()) { stops.push((e.seg.end().1, true)); }
        if at_axis(e.seg.start()) { stops.push((e.seg.start().1, false)); }
    }
    stops.sort_by(|x, y| y.0.total_cmp(&x.0));
    if stops.len() % 2 == 1 || stops.chunks(2).any(|w| !(w[0].1 && !w[1].1) || w[0].0 - w[1].0 <= tol) {
        return unsupported("the offset section meets the axis in a way that does not close into a solid of revolution".into());
    }
    for w in stops.chunks(2) {
        pieces.push(Elem { seg: Seg2::Line { a: (0.0, w[0].0), b: (0.0, w[1].0) }, tag: Tag::Axis, corner: Tag::Axis, mirrored: false, removed: false });
    }
    chain(pieces, tol, "the offset section")
}

/// Pieces joined end to start into loops (each end meets exactly one start within tolerance).
fn chain(pieces: Vec<Elem>, tol: f64, what: &str) -> Result<Vec<Vec<Elem>>, KernelError> {
    let mut used = vec![false; pieces.len()];
    let mut loops = Vec::new();
    for first in 0..pieces.len() {
        if used[first] { continue; }
        let mut lp = Vec::new();
        let mut k = first;
        loop {
            used[k] = true;
            lp.push(pieces[k]);
            let end = pieces[k].seg.end();
            let next: Vec<usize> = (0..pieces.len()).filter(|&j| dist(pieces[j].seg.start(), end) <= tol && (j == first || !used[j])).collect();
            match next[..] {
                [j] if j == first => break,
                [j] => k = j,
                _ => return unsupported(format!("{what}: its pieces do not join into closed loops ({} continue one end)", next.len())),
            }
        }
        // Each piece starts exactly where the one before it ends.
        for i in 0..lp.len() {
            let e = lp[i].seg.end();
            let n = lp.len();
            let next = &mut lp[(i + 1) % n];
            if let Seg2::Line { a, .. } = &mut next.seg { *a = e; }
        }
        loops.push(lp);
    }
    Ok(loops)
}

/// Turns the half-plane loops into the result body and the source of each of its faces.
fn build(body: &Body, rev: &Revolved, half: &[Vec<Elem>], tol: f64) -> Result<(Body, Vec<Source>, f64), KernelError> {
    let revolution = revolution(rev, half, tol)?;
    let (out, prov) = revolve::revolve_general_named(&revolution, tol)?;
    let volume = revolve::volume_of(&revolution, tol)?;
    let mut sources = Vec::with_capacity(out.faces.len());
    for (f, rec) in prov.faces.iter().enumerate() {
        let (l, j) = (rec.refs[0], rec.refs[1]);
        let want = heading(&out, f, rev);
        let source = match half[l][j].tag {
            Tag::Seg(sl, sj) => Source::Face(pick(&rev.faces[sl][sj], want, &|k| heading(body, k, rev))),
            Tag::Corner(sl, sj) => match &rev.corners[sl][sj] {
                CornerSource::Parallels(edges) if !edges.is_empty() => {
                    let of = |e: usize| { let c = &body.edges[e].curve; let (t0, t1) = c.domain(); let w = rev.radial(c.eval((t0 + t1) / 2.0).0); (w.length() > 0.0).then(|| w.scale(1.0 / w.length())) };
                    Source::Edge(pick(edges, want, &of))
                }
                CornerSource::AxisVertex(Some(v)) => Source::Vertex(*v),
                _ => return Err(KernelError::Internal(format!("section corner ({sl}, {sj}) has no source edge or vertex"))),
            },
            Tag::Axis => return Err(KernelError::Internal("a face was swept from the axis".into())),
        };
        sources.push(source);
    }
    Ok((out, sources, volume))
}

/// Every whole-plane loop offset (per-element distances), checked simple and keeping its sense.
fn offset_all(whole: &[Vec<Elem>], delta: &dyn Fn(&Elem) -> f64, rounded: bool, tol: f64, what: &str) -> Result<Vec<Vec<Elem>>, KernelError> {
    let rings: Vec<Vec<Elem>> = whole.iter().enumerate()
        .map(|(l, lp)| offset_loop(lp, delta, rounded, tol, &format!("{what} loop {l}")).map(|o| o.elems))
        .collect::<Result<_, _>>()?;
    check_simple(&rings, tol, what)?;
    planar::check_nesting(whole, &rings, what)?;
    for (l, (o, s)) in rings.iter().zip(whole).enumerate() {
        if area2(o).signum() != area2(s).signum() { return unsupported(format!("{what}: loop {l} turns inside out")); }
    }
    Ok(rings)
}

pub fn offset(body: &Body, rev: &Revolved, distance: f64, sharp: bool, tol: f64) -> Result<(Body, Vec<Source>, f64), KernelError> {
    let whole = whole_plane(&rev.loops);
    let rings = offset_all(&whole, &|_| distance, !sharp, tol, "the offset section")?;
    let half = half_plane(&rings, tol)?;
    build(body, rev, &half, tol)
}

/// The part of `outer` (a removed element's piece in the outer region) not covered by `inner` (the same element's piece
/// in the cavity, along the same support): the rims the removed face leaves.
fn rims(outer: &Elem, inner: &Elem, tol: f64) -> Result<Vec<Elem>, KernelError> {
    let refuse = || unsupported(
        "a removed face meets a kept face at a concave or tangent edge, where the shell's wall would not end on the removed face".into());
    let frac = |q: P2| -> f64 {
        match outer.seg {
            Seg2::Line { a, b } => { let d = planar::sub(b, a); planar::dot(planar::sub(q, a), d) / planar::dot(d, d) }
            Seg2::Arc { c, a0, sw, .. } => {
                let t = (q.1 - c.1).atan2(q.0 - c.0);
                let rel = if sw > 0.0 { (t - a0).rem_euclid(TAU) } else { (a0 - t).rem_euclid(TAU) };
                let rel = if rel > TAU - 1e-12 { 0.0 } else { rel };
                rel / sw.abs()
            }
        }
    };
    let (f0, f1) = (frac(inner.seg.start()), frac(inner.seg.end()));
    let slack = tol / outer.seg.length();
    if !(f0 >= -slack && f1 <= 1.0 + slack && f1 > f0) { return refuse(); }
    let piece = |from: f64, to: f64, s: P2, e: P2| -> Option<Elem> {
        if (to - from) * outer.seg.length() <= tol { return None; }
        let seg = match outer.seg {
            Seg2::Line { .. } => Seg2::Line { a: s, b: e },
            Seg2::Arc { c, r, a0, sw } => Seg2::Arc { c, r, a0: a0 + sw * from, sw: sw * (to - from) },
        };
        Some(Elem { seg, removed: false, ..*outer })
    };
    Ok([piece(0.0, f0, outer.seg.start(), inner.seg.start()), piece(f1, 1.0, inner.seg.end(), outer.seg.end())].into_iter().flatten().collect())
}

/// The solid shelled to wall t (positive inward), open through the removed faces.
pub fn shell_open(body: &Body, rev: &Revolved, t: f64, removed: &[usize], tol: f64) -> Result<(Body, Vec<Source>, f64), KernelError> {
    let mut half = rev.loops.clone();
    for (l, lp) in half.iter_mut().enumerate() {
        for (j, e) in lp.iter_mut().enumerate() {
            let faces = &rev.faces[l][j];
            let n = faces.iter().filter(|f| removed.contains(f)).count();
            if n > 0 && n < faces.len() {
                return unsupported(format!("faces {faces:?} are pieces of one revolved surface: remove all of them or none"));
            }
            e.removed = n > 0;
        }
    }
    let whole = whole_plane(&half);
    let wall = |e: &Elem| if e.removed { 0.0 } else { -t };
    let away = |e: &Elem| if e.removed { 0.0 } else { t.abs() };
    let (outer, inner) = if t > 0.0 {
        (whole.clone(), offset_all(&whole, &wall, true, tol, "the shell's cavity")?)
    } else {
        (offset_all(&whole, &away, true, tol, "the shell's outer wall")?, whole.clone())
    };
    let mut pieces: Vec<Elem> = Vec::new();
    for e in outer.iter().flatten() {
        if !e.removed { pieces.push(*e); continue; }
        let Some(i) = inner.iter().flatten().find(|i| i.removed && i.tag == e.tag && i.mirrored == e.mirrored) else {
            return Err(KernelError::Internal("a removed face's piece is missing from the cavity".into()));
        };
        pieces.extend(rims(e, i, tol)?);
    }
    pieces.extend(inner.iter().flatten().filter(|e| !e.removed).map(|e| Elem { seg: e.seg.reversed(), ..*e }));
    let region = chain(pieces, tol, "the shell's section")?;
    check_simple(&region, tol, "the shell's section")?;
    let half = half_plane(&region, tol)?;
    build(body, rev, &half, tol)
}

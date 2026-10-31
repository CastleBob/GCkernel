//! Display tessellation of a body (Claude 2026-10-28, for Astra's GC Kernel viewport; contract in GCm
//! GCmaker/Handoff/Claude-GcKernel-Status.md Log 30). Display only: the exact B-rep stays authoritative, and
//! nothing here feeds a later boolean, save or machining.
//!
//! Every edge is cut once into a polyline within the chord tolerance, and both faces of the edge use that same
//! polyline, so the mesh is closed. A planar face is triangulated in its plane (ear clipping; holes joined to the
//! outer loop by bridges). A cylindrical face is triangulated in (r theta, s), its loops unwrapped round the axis
//! (a full face needs its seam edge), then inner edges are split until no triangle spans more angle than the chord
//! tolerance allows. Triangles run counter-clockwise about the outward normal; vertices are repeated per face.
//! A face that cannot be triangulated is an error, never patched.

use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::Body;
use crate::KernelError;
use std::collections::HashMap;
use std::f64::consts::{PI, TAU};

pub struct Mesh {
    pub chord: f64,
    pub positions: Vec<Vec3>,
    pub normals: Vec<Vec3>,
    /// Counter-clockwise about the outward normal.
    pub triangles: Vec<[usize; 3]>,
    /// The body face of each triangle.
    pub faces: Vec<usize>,
    /// Each body edge's polyline, start vertex to end vertex.
    pub edges: Vec<Vec<Vec3>>,
    /// Body edge IDs to draw. Same-support face divisions and periodic seams remain in topology/triangulation.
    pub display_edges: Vec<usize>,
}

/// Only suppress divisions of the SAME analytic support with the same outward sense. Tangent fillet
/// boundaries (plane/cylinder or different cylinders) and even shallow real creases remain visible.
/// This roundoff allowance is independent of display chord and modelling tolerance; no body is changed.
fn same_support(a: &crate::neutral::Face, b: &crate::neutral::Face) -> bool {
    const ROUND: f64 = 64.0 * f64::EPSILON;
    let sign = |sense: bool| if sense { 1.0 } else { -1.0 };
    match (&a.surface, &b.surface) {
        (Surface::Plane { p, n, .. }, Surface::Plane { p: q, n: m, .. }) => {
            let d = q.sub(*p);
            let magnitude = (d.x * n.x).abs() + (d.y * n.y).abs() + (d.z * n.z).abs();
            n.scale(sign(a.sense)).sub(m.scale(sign(b.sense))).length() <= ROUND
                && d.dot(*n).abs() <= ROUND * magnitude.max(1.0)
        }
        (Surface::Cylinder { p, a: n, r, .. }, Surface::Cylinder { p: q, a: m, r: s, .. }) => {
            let d = q.sub(*p);
            a.sense == b.sense && n.cross(*m).length() <= ROUND
                && d.cross(*n).length() <= ROUND * d.length().max(1.0)
                && (r - s).abs() <= ROUND * r.abs().max(s.abs()).max(1.0)
        }
        // One ellipse swept along one axis: same semi-diameters (a seam, or pieces of one extruded wall).
        (Surface::EllipticCylinder { p, a: n, u, v }, Surface::EllipticCylinder { p: q, a: m, u: u2, v: v2 }) => {
            let d = q.sub(*p);
            let scale = u.length().max(v.length()).max(1.0);
            a.sense == b.sense && n.cross(*m).length() <= ROUND
                && d.cross(*n).length() <= ROUND * d.length().max(1.0)
                && u.sub(*u2).length() <= ROUND * scale && v.sub(*v2).length() <= ROUND * scale
        }
        // One cone: the same axis line and orientation-independent radius function (radius at p, slope along the axis).
        (Surface::Cone { p, a: n, r, k, .. }, Surface::Cone { p: q, a: m, r: r2, k: k2, .. }) => {
            let d = p.sub(*q);
            let along = n.dot(*m).signum();
            let there = r2 + k2 * d.dot(*m);
            a.sense == b.sense && n.cross(*m).length() <= ROUND
                && d.cross(*n).length() <= ROUND * d.length().max(1.0)
                && (r - there).abs() <= ROUND * r.abs().max(there.abs()).max(1.0)
                && (k - k2 * along).abs() <= ROUND * k.abs().max(k2.abs()).max(1.0)
        }
        (Surface::Sphere { c, r, .. }, Surface::Sphere { c: c2, r: r2, .. }) => {
            a.sense == b.sense && c.distance(*c2) <= ROUND * c.length().max(1.0) && (r - r2).abs() <= ROUND * r.max(*r2).max(1.0)
        }
        (Surface::Torus { c, a: n, big, r, .. }, Surface::Torus { c: c2, a: m, big: big2, r: r2, .. }) => {
            a.sense == b.sense && n.cross(*m).length() <= ROUND && c.distance(*c2) <= ROUND * c.length().max(1.0)
                && (big - big2).abs() <= ROUND * big.max(*big2).max(1.0) && (r - r2).abs() <= ROUND * r.max(*r2).max(1.0)
        }
        // One meridian swept about one axis: the same meridian data, axis line and orientation (a seam, or pieces of one
        // swept curve). Distinct meridian pieces (different data) stay visible.
        (Surface::Revolution { p, a: n, m, .. }, Surface::Revolution { p: q, a: k, m: m2, .. }) => {
            let d = q.sub(*p);
            let close = |x: f64, y: f64| (x - y).abs() <= ROUND * x.abs().max(y.abs()).max(1.0);
            let same_meridian = match (m, m2) {
                (crate::meridian::Meridian::Ellipse { c, u, v }, crate::meridian::Meridian::Ellipse { c: c2, u: u2, v: v2 }) =>
                    [c.0, c.1, u.0, u.1, v.0, v.1].iter().zip([c2.0, c2.1, u2.0, u2.1, v2.0, v2.1]).all(|(x, y)| close(*x, y)),
                (crate::meridian::Meridian::Bezier(b1), crate::meridian::Meridian::Bezier(b2)) =>
                    b1.iter().zip(b2.iter()).all(|(x, y)| close(x.0, y.0) && close(x.1, y.1)),
                _ => false,
            };
            a.sense == b.sense && n.sub(*k).length() <= ROUND && d.length() <= ROUND * p.length().max(1.0) && same_meridian
        }
        _ => false,
    }
}

fn display_edges(body: &Body) -> Vec<usize> {
    let mut uses = vec![Vec::new(); body.edges.len()];
    for (f, face) in body.faces.iter().enumerate() {
        for lp in &face.loops {
            for &(e, _) in lp { uses[e].push(f); }
        }
    }
    uses.iter().enumerate().filter_map(|(e, faces)| {
        let internal = faces.len() == 2 && same_support(&body.faces[faces[0]], &body.faces[faces[1]]);
        if internal { None } else { Some(e) }
    }).collect()
}

fn unsupported<T>(m: impl Into<String>) -> Result<T, KernelError> { Err(KernelError::Unsupported(m.into())) }

/// The angle step whose chord stays within `chord` of a circle of radius r (at most a quarter turn).
/// Stable for any radius (Codex review: acos(1 - chord / r) rounds to zero when r is huge): the sagitta r (1 - cos(h / 2))
/// = 2 r sin^2(h / 4) inverted as h = 4 asin(sqrt(chord / (2 r))).
fn angle_step(r: f64, chord: f64) -> f64 {
    if !(r > 0.0) || chord >= r { PI / 2.0 } else { (4.0 * (chord / (2.0 * r)).sqrt().asin()).min(PI / 2.0) }
}

/// An edge as points from its start vertex to its end vertex.
fn polyline(body: &Body, edge: usize, chord: f64) -> Result<Vec<Vec3>, KernelError> {
    let e = &body.edges[edge];
    let (t0, t1) = e.curve.domain();
    let at = |t: f64| e.curve.eval(t).0;
    let params: Vec<f64> = match e.curve {
        Curve::Line { .. } => vec![0.0, 1.0],
        Curve::Arc { r, .. } => {
            let steps = ((t1 - t0) / angle_step(r, chord)).ceil().max(1.0);
            if !(steps <= 1_000_000.0) { return unsupported("an arc edge needs over 1000000 display points"); }
            let n = steps as usize;
            (0..=n).map(|k| t0 + (t1 - t0) * k as f64 / n as f64).collect()
        }
        Curve::Conic { u, v, .. } => {
            // |P''| <= sqrt(|u|^2 + |v|^2); a chord over step h departs by at most h^2 |P''| / 8.
            let bend = (u.dot(u) + v.dot(v)).sqrt();
            let step = (8.0 * chord / bend).sqrt().min(PI / 2.0);
            let steps = ((t1 - t0) / step).ceil().max(1.0);
            if !(steps <= 1_000_000.0) { return unsupported("an ellipse edge needs over 1000000 display points"); }
            let n = steps as usize;
            (0..=n).map(|k| t0 + (t1 - t0) * k as f64 / n as f64).collect()
        }
        Curve::RationalBezier { p, w } => {
            crate::rational::validate(p, w)?;
            let mut params = vec![0.0];
            let mut pending = vec![(p, w, 0.0, 1.0, 0)];
            let mut visits = 0;
            while let Some((p, w, a, b, depth)) = pending.pop() {
                visits += 1;
                if visits > 200_000 { return unsupported("rational edge exceeds display subdivision budget"); }
                let d = p[3].sub(p[0]);
                let l = d.dot(d);
                // Positive rational weights keep the ENTIRE curve in its control hull.
                let off = p.iter().map(|q| {
                    let s = if l > 0.0 { (q.sub(p[0]).dot(d)/l).clamp(0.0, 1.0) } else { 0.0 };
                    q.distance(p[0].add(d.scale(s)))
                }).fold(0.0, f64::max);
                if !off.is_finite() { return unsupported("nonfinite rational display bound"); }
                if off <= chord/2.0 { params.push(b); continue; }
                if depth >= 24 { return unsupported("rational edge display bound did not converge"); }
                let [left, right] = crate::rational::split(p, w, 0.5);
                let m = (a+b)/2.0;
                pending.push((right.0, right.1, m, b, depth+1));
                pending.push((left.0, left.1, a, m, depth+1));
            }
            params
        }
        Curve::CylinderPair { .. } | Curve::Bezier { .. } => {
            // No closed-form bend: halve each span until its quarter points sit within chord / 2 of the chord.
            let mut params = vec![t0];
            let mut pending = vec![(t0, t1, 0u32)];
            while let Some((a, b, depth)) = pending.pop() {
                let (pa, pb) = (at(a), at(b));
                let off = [0.25, 0.5, 0.75].iter().map(|f| {
                    let q = at(a + (b - a) * f);
                    let d = pb.sub(pa);
                    let l = d.dot(d);
                    let s = if l > 0.0 { (q.sub(pa).dot(d) / l).clamp(0.0, 1.0) } else { 0.0 };
                    q.distance(pa.add(d.scale(s)))
                }).fold(0.0, f64::max);
                if depth < 3 || (off > chord / 2.0 && depth < 24) {
                    let m = (a + b) / 2.0;
                    pending.push((m, b, depth + 1));
                    pending.push((a, m, depth + 1));
                } else {
                    params.push(b);
                }
                if params.len() > 100_000 { return unsupported("an intersection edge needs over 100000 display points"); }
            }
            params
        }
    };
    let mut points: Vec<Vec3> = params.iter().map(|&t| at(t)).collect();
    // The ends are the edge's vertices exactly, so both faces and the neighbouring edges share them.
    let n = points.len();
    points[0] = body.vertices[e.start];
    points[n - 1] = body.vertices[e.end];
    Ok(points)
}

/// A loop as points, each joint once.
fn loop_points(lines: &[Vec<Vec3>], lp: &[(usize, bool)]) -> Vec<Vec3> {
    let mut out = Vec::new();
    for &(edge, forward) in lp {
        let line = &lines[edge];
        let n = line.len();
        for k in 0..n - 1 {
            out.push(if forward { line[k] } else { line[n - 1 - k] });
        }
    }
    out
}

fn cross2(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 { (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) }
fn area2(p: &[[f64; 2]], idx: &[usize]) -> f64 {
    (0..idx.len()).map(|i| { let (a, b) = (p[idx[i]], p[idx[(i + 1) % idx.len()]]); a[0] * b[1] - a[1] * b[0] }).sum::<f64>() / 2.0
}
fn contains(p: &[[f64; 2]], idx: &[usize], q: [f64; 2]) -> bool {
    let mut inside = false;
    for i in 0..idx.len() {
        let (a, b) = (p[idx[i]], p[idx[(i + 1) % idx.len()]]);
        if (a[1] > q[1]) != (b[1] > q[1]) && q[0] < a[0] + (q[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]) { inside = !inside; }
    }
    inside
}

/// Joins a clockwise hole to a counter-clockwise polygon through a bridge from the hole's rightmost point to a
/// polygon vertex it can see (the standard ray-cast construction).
fn bridge(p: &[[f64; 2]], poly: &mut Vec<usize>, hole: &[usize]) -> Result<(), KernelError> {
    let (hm, &m) = hole.iter().enumerate().max_by(|a, b| p[*a.1][0].total_cmp(&p[*b.1][0])).unwrap();
    let pm = p[m];
    // Nearest crossing of the ray +x from M with the polygon.
    let mut best: Option<(f64, usize)> = None;
    for i in 0..poly.len() {
        let (a, b) = (p[poly[i]], p[poly[(i + 1) % poly.len()]]);
        if (a[1] <= pm[1] && b[1] > pm[1]) || (b[1] <= pm[1] && a[1] > pm[1]) || (a[1] == pm[1] && b[1] == pm[1]) {
            let x = if a[1] == b[1] { a[0].min(b[0]) } else { a[0] + (pm[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]) };
            if x >= pm[0] && best.is_none_or(|(bx, _)| x < bx) { best = Some((x, i)); }
        }
    }
    let Some((ix, i)) = best else { return unsupported("a face's hole lies outside its boundary") };
    let (a, b) = (i, (i + 1) % poly.len());
    let mut target = if p[poly[a]][0] > p[poly[b]][0] { a } else { b };
    let pi = [ix, pm[1]];
    // A polygon vertex inside triangle M, I, P would block the bridge: take the one nearest the ray's direction.
    let tp = p[poly[target]];
    let mut best_angle = f64::INFINITY;
    for (k, &v) in poly.iter().enumerate() {
        let q = p[v];
        if q == pm { continue; }
        // A bridge endpoint may occur more than once after earlier holes were joined. Each occurrence
        // owns a different angular sector: the new bridge must enter the polygon's LOCAL interior.
        let prev = p[poly[(k + poly.len() - 1) % poly.len()]];
        let next = p[poly[(k + 1) % poly.len()]];
        let incoming = cross2(prev, q, pm) >= 0.0;
        let outgoing = cross2(q, next, pm) >= 0.0;
        let local = if cross2(prev, q, next) >= 0.0 { incoming && outgoing } else { incoming || outgoing };
        if !local { continue; }
        let (s1, s2, s3) = (cross2(pm, pi, q), cross2(pi, tp, q), cross2(tp, pm, q));
        let inside = (s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0) || (s1 <= 0.0 && s2 <= 0.0 && s3 <= 0.0);
        if inside && q[0] >= pm[0] {
            let angle = (q[1] - pm[1]).atan2(q[0] - pm[0]).abs();
            if angle < best_angle { best_angle = angle; target = k; }
        }
    }
    if !best_angle.is_finite() { return unsupported("a face hole has no locally visible bridge endpoint"); }
    let mut joined = Vec::with_capacity(poly.len() + hole.len() + 2);
    joined.extend_from_slice(&poly[..=target]);
    for k in 0..=hole.len() { joined.push(hole[(hm + k) % hole.len()]); }
    joined.extend_from_slice(&poly[target..]);
    *poly = joined;
    Ok(())
}

/// Ear clipping of a counter-clockwise polygon (indices into p, possibly repeated at bridges). `short` takes the
/// ear with the shortest new side each time (cylinders: a strip of short triangles, not long fans across the
/// angle), at n^3 cost; otherwise the next ear found.
/// Which straight runs of a chart are straight on the surface too (where a zero-area ear may keep a vertex): any on a
/// plane, those along the axis on a cylinder or cone, none on a sphere or torus (both directions curve).
#[derive(Clone, Copy, PartialEq)]
enum Runs { Any, Axis, None }

fn clip(p: &[[f64; 2]], mut poly: Vec<usize>, short: bool, runs: Runs, out: &mut Vec<[usize; 3]>) -> Result<(), KernelError> {
    let scale = poly.iter().map(|&i| p[i][0].abs().max(p[i][1].abs())).fold(1.0, f64::max);
    let eps = 1e-14 * scale * scale;
    let mut i = 0;
    let mut stall = 0;
    while poly.len() > 3 {
        let n = poly.len();
        // A straight-run vertex (a zero-area ear) is clipped only when no proper ear is left (Codex acceptance: display
        // triangles have area); otherwise it stays and joins proper triangles.
        let proper = |k: usize| {
            let n = poly.len();
            cross2(p[poly[(k + n - 1) % n]], p[poly[k]], p[poly[(k + 1) % n]]).abs() > eps && is_ear(p, &poly, k, eps, runs)
        };
        if short {
            let side = |k: usize| { let (a, c) = (p[poly[(k + n - 1) % n]], p[poly[(k + 1) % n]]); (a[0] - c[0]).hypot(a[1] - c[1]) };
            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by(|&x, &y| side(x).total_cmp(&side(y)).then(x.cmp(&y)));
            i = order.iter().copied().find(|&k| proper(k)).or_else(|| order.iter().copied().find(|&k| is_ear(p, &poly, k, eps, runs)))
                .ok_or_else(|| KernelError::Unsupported("a face could not be triangulated for display".into()))?;
        } else if !proper(i % n) {
            if let Some(k) = (0..n).map(|d| (i + d) % n).find(|&k| proper(k)) { i = k; }
        }
        let (ia, ib, ic) = ((i + n - 1) % n, i % n, (i + 1) % n);
        if is_ear(p, &poly, ib, eps, runs) {
            out.push([poly[ia], poly[ib], poly[ic]]);
            poly.remove(ib);
            stall = 0;
            i = if ib == 0 { 0 } else { ib - 1 };
        } else {
            i += 1;
            stall += 1;
            if stall > 2 * n { return unsupported("a face could not be triangulated for display"); }
        }
    }
    if cross2(p[poly[0]], p[poly[1]], p[poly[2]]) > -eps { out.push([poly[0], poly[1], poly[2]]); }
    else { return unsupported("a face could not be triangulated for display"); }
    Ok(())
}

/// Whether vertex ib of the polygon is an ear: a convex corner with no other vertex in or on its triangle, or a
/// vertex on a straight run (a zero-area triangle keeps it, so the faces either side share every edge point). On a
/// cylinder (`curved`) only a run along the axis is straight: one along the angle is an arc.
fn is_ear(p: &[[f64; 2]], poly: &[usize], ib: usize, eps: f64, runs: Runs) -> bool {
    let n = poly.len();
    let (ia, ic) = ((ib + n - 1) % n, (ib + 1) % n);
    let (a, b, c) = (p[poly[ia]], p[poly[ib]], p[poly[ic]]);
    let turn = cross2(a, b, c);
    if turn.abs() <= eps {
        if runs == Runs::None || (runs == Runs::Axis && (a[0] != b[0] || b[0] != c[0])) { return false; }
        return a != c && (b[0] - a[0]) * (c[0] - b[0]) + (b[1] - a[1]) * (c[1] - b[1]) > 0.0;
    }
    if turn < 0.0 { return false; }
    poly.iter().enumerate().all(|(k, &v)| {
        let q = p[v];
        k == ia || k == ib || k == ic || q == a || q == b || q == c
            || !(cross2(a, b, q) >= 0.0 && cross2(b, c, q) >= 0.0 && cross2(c, a, q) >= 0.0)
    })
}

/// Triangles over the loops (outer counter-clockwise, holes clockwise in p's frame).
fn triangulate(p: &[[f64; 2]], loops: &[Vec<usize>], short: bool) -> Result<Vec<[usize; 3]>, KernelError> {
    triangulate_runs(p, loops, short, if short { Runs::Axis } else { Runs::Any })
}

fn triangulate_runs(p: &[[f64; 2]], loops: &[Vec<usize>], short: bool, runs: Runs) -> Result<Vec<[usize; 3]>, KernelError> {
    let (mut outers, mut holes): (Vec<&Vec<usize>>, Vec<&Vec<usize>>) = (Vec::new(), Vec::new());
    for l in loops {
        if l.len() < 3 { return unsupported("a face loop has fewer than three display points"); }
        if area2(p, l) > 0.0 { outers.push(l) } else { holes.push(l) }
    }
    let polys: Vec<Vec<usize>> = outers.iter().map(|l| l.to_vec()).collect();
    let mut owned: Vec<Vec<&Vec<usize>>> = vec![Vec::new(); polys.len()];
    for h in holes {
        let q = p[h[0]];
        let owner = (0..polys.len()).filter(|&k| contains(p, &polys[k], q))
            .min_by(|&x, &y| area2(p, &polys[x]).total_cmp(&area2(p, &polys[y])));
        let Some(k) = owner else { return unsupported("a face's hole lies outside its boundary") };
        owned[k].push(h);
    }
    let mut out = Vec::new();
    for (k, mut poly) in polys.into_iter().enumerate() {
        let mut hs = owned[k].clone();
        hs.sort_by(|x, y| {
            let mx = |l: &Vec<usize>| l.iter().map(|&i| p[i][0]).fold(f64::NEG_INFINITY, f64::max);
            mx(y).total_cmp(&mx(x))
        });
        for h in hs { bridge(p, &mut poly, h)?; }
        clip(p, poly, short, runs, &mut out)?;
    }
    Ok(out)
}

pub fn tessellate(body: &Body, chord: f64) -> Result<Mesh, KernelError> {
    if !chord.is_finite() || chord <= 0.0 { return Err(KernelError::InvalidInput("the chord tolerance must be positive".into())); }
    let lines: Vec<Vec<Vec3>> = (0..body.edges.len()).map(|e| polyline(body, e, chord)).collect::<Result<_, _>>()?;
    let mut mesh = Mesh { chord, positions: Vec::new(), normals: Vec::new(), triangles: Vec::new(), faces: Vec::new(), edges: lines.clone(), display_edges: display_edges(body) };
    for (f, face) in body.faces.iter().enumerate() {
        let loops3: Vec<Vec<Vec3>> = face.loops.iter().map(|l| loop_points(&lines, l)).collect();
        let mut points3: Vec<Vec3> = loops3.concat();
        let mut index = Vec::new();
        let mut at = 0;
        for l in &loops3 { index.push((at..at + l.len()).collect::<Vec<usize>>()); at += l.len(); }
        let side = if face.sense { 1.0 } else { -1.0 };
        let (tris, normals) = match face.surface {
            Surface::Spline { ref data } => {
                data.check_face(body, face, 1e-6)?;
                let (points, tris, normals) = data.mesh(points3, side, chord)?;
                points3 = points;
                (tris, normals)
            }
            Surface::Plane { p, n, x } => {
                let w = n.scale(side);
                let v = w.cross(x);
                let p2: Vec<[f64; 2]> = points3.iter().map(|q| { let d = q.sub(p); [d.dot(x), d.dot(v)] }).collect();
                let tris = triangulate(&p2, &index, false).map_err(|e| KernelError::Unsupported(format!("planar face {f}: {e}")))?;
                (tris, vec![w; points3.len()])
            }
            Surface::Cylinder { p, a, x, r } => {
                let y = a.cross(x);
                // (r theta, s) with outward radial (sense) giving counter-clockwise loops; inward flips s.
                let flip = side;
                let mut p2 = Vec::with_capacity(points3.len());
                let mut spans = Vec::new();
                for l in &index {
                    let mut prev: Option<f64> = None;
                    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                    for &i in l {
                        let d = points3[i].sub(p);
                        let raw = d.dot(y).atan2(d.dot(x));
                        let theta = match prev { None => raw, Some(t) => t + (raw - t + PI).rem_euclid(TAU) - PI };
                        prev = Some(theta);
                        lo = lo.min(theta);
                        hi = hi.max(theta);
                        p2.push([r * theta, flip * d.dot(a)]);
                    }
                    let (first, last) = (p2[l[0]][0] / r, prev.unwrap());
                    let close = last + (first - last + PI).rem_euclid(TAU) - PI;
                    if (close - first).abs() > PI {
                        return unsupported("a cylindrical face loop wraps the axis without a seam edge");
                    }
                    spans.push((lo, hi));
                }
                // Holes unwrapped a turn away from their outer loop are moved into its range.
                if let Some(k) = (0..index.len()).max_by(|&i, &j| (spans[i].1 - spans[i].0).total_cmp(&(spans[j].1 - spans[j].0))) {
                    let (lo, hi) = spans[k];
                    for (j, l) in index.iter().enumerate() {
                        if j == k { continue; }
                        let mid = (spans[j].0 + spans[j].1) / 2.0;
                        let turns = if mid < lo { ((lo - mid) / TAU).ceil() } else if mid > hi { -((mid - hi) / TAU).ceil() } else { 0.0 };
                        for &i in l { p2[i][0] += turns * TAU * r; }
                    }
                }
                let mut tris = triangulate(&p2, &index, true)?;
                // Split inner edges spanning more angle than the chord allows, lifting new points to the cylinder.
                let limit = r * angle_step(r, chord);
                let boundary: std::collections::HashSet<(usize, usize)> = index.iter()
                    .flat_map(|l| (0..l.len()).map(move |k| { let (u, w) = (l[k], l[(k + 1) % l.len()]); (u.min(w), u.max(w)) }))
                    .collect();
                for _pass in 0..64 {
                    let mut by_edge: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
                    for (t, tri) in tris.iter().enumerate() {
                        for k in 0..3 { let (u, w) = (tri[k], tri[(k + 1) % 3]); by_edge.entry((u.min(w), u.max(w))).or_default().push(t); }
                    }
                    let mut long: Vec<(usize, usize)> = by_edge.keys().copied()
                        .filter(|e| !boundary.contains(e) && (p2[e.0][0] - p2[e.1][0]).abs() > limit).collect();
                    if long.is_empty() { break; }
                    long.sort_by(|e1, e2| (p2[e2.0][0] - p2[e2.1][0]).abs().total_cmp(&(p2[e1.0][0] - p2[e1.1][0]).abs()).then(e1.cmp(e2)));
                    let mut touched = vec![false; tris.len()];
                    for (u, w) in long {
                        let around = &by_edge[&(u, w)];
                        if around.iter().any(|&t| touched[t]) { continue; }
                        let m = p2.len();
                        let mid = [(p2[u][0] + p2[w][0]) / 2.0, (p2[u][1] + p2[w][1]) / 2.0];
                        p2.push(mid);
                        let theta = mid[0] / r;
                        points3.push(p.add(a.scale(flip * mid[1])).add(x.scale(r * theta.cos())).add(y.scale(r * theta.sin())));
                        for &t in around {
                            touched[t] = true;
                            let tri = tris[t];
                            let k = (0..3).find(|&k| { let (s, e) = (tri[k], tri[(k + 1) % 3]); (s == u && e == w) || (s == w && e == u) }).unwrap();
                            let (s, e, o) = (tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]);
                            tris[t] = [s, m, o];
                            tris.push([m, e, o]);
                            touched.push(true);
                        }
                    }
                    if _pass == 63 { return unsupported("a cylindrical face needs too many display triangles"); }
                }
                let normals = points3.iter().map(|q| {
                    let d = q.sub(p);
                    d.sub(a.scale(d.dot(a))).scale(side / r)
                }).collect();
                (tris, normals)
            }
            Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } => {
                let rev = Revolved::of(&face.surface).expect("a surface of revolution");
                let (pts, tris, normals) = revolved_mesh(&rev, &loops3, side, chord, f)?;
                points3 = pts;
                (tris, normals)
            }
            Surface::EllipticCylinder { p, a, u, v } => {
                // As the cylinder, in (R t, s) with t the elliptic parameter (exact on the surface) and R the mean
                // semi-diameter (Claude 2026-10-29). Inner edges longer than the chord allows are split and the new
                // points lifted onto the surface; the chord over a step h is at most h^2 |C''| / 8, |C''| <= sqrt(|u|^2 + |v|^2).
                let (uu, vv) = (u.dot(u), v.dot(v));
                let param = |q: Vec3| {
                    let d = q.sub(p);
                    let (c, n) = Surface::elliptic_coordinates(a, u, v, d);
                    (n.atan2(c), d.dot(a))
                };
                let big_r = (uu.sqrt() + vv.sqrt()) / 2.0;
                let flip = side;
                let mut p2 = Vec::with_capacity(points3.len());
                let mut spans = Vec::new();
                for l in &index {
                    let mut prev: Option<f64> = None;
                    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                    for &i in l {
                        let (raw, along) = param(points3[i]);
                        let theta = match prev { None => raw, Some(t) => t + (raw - t + PI).rem_euclid(TAU) - PI };
                        prev = Some(theta);
                        lo = lo.min(theta);
                        hi = hi.max(theta);
                        p2.push([big_r * theta, flip * along]);
                    }
                    let (first, last) = (p2[l[0]][0] / big_r, prev.unwrap());
                    let close = last + (first - last + PI).rem_euclid(TAU) - PI;
                    if (close - first).abs() > PI {
                        return unsupported("an elliptic-cylinder face loop wraps the axis without a seam edge");
                    }
                    spans.push((lo, hi));
                }
                if let Some(k) = (0..index.len()).max_by(|&i, &j| (spans[i].1 - spans[i].0).total_cmp(&(spans[j].1 - spans[j].0))) {
                    let (lo, hi) = spans[k];
                    for (j, l) in index.iter().enumerate() {
                        if j == k { continue; }
                        let mid = (spans[j].0 + spans[j].1) / 2.0;
                        let turns = if mid < lo { ((lo - mid) / TAU).ceil() } else if mid > hi { -((mid - hi) / TAU).ceil() } else { 0.0 };
                        for &i in l { p2[i][0] += turns * TAU * big_r; }
                    }
                }
                let mut tris = triangulate(&p2, &index, true)?;
                let bend = (uu + vv).sqrt();
                let limit = big_r * (8.0 * chord / bend).sqrt().min(PI / 2.0);
                let boundary: std::collections::HashSet<(usize, usize)> = index.iter()
                    .flat_map(|l| (0..l.len()).map(move |k| { let (s0, w) = (l[k], l[(k + 1) % l.len()]); (s0.min(w), s0.max(w)) }))
                    .collect();
                for _pass in 0..64 {
                    let mut by_edge: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
                    for (t, tri) in tris.iter().enumerate() {
                        for k in 0..3 { let (s0, w) = (tri[k], tri[(k + 1) % 3]); by_edge.entry((s0.min(w), s0.max(w))).or_default().push(t); }
                    }
                    let mut long: Vec<(usize, usize)> = by_edge.keys().copied()
                        .filter(|e| !boundary.contains(e) && (p2[e.0][0] - p2[e.1][0]).abs() > limit).collect();
                    if long.is_empty() { break; }
                    long.sort_by(|e1, e2| (p2[e2.0][0] - p2[e2.1][0]).abs().total_cmp(&(p2[e1.0][0] - p2[e1.1][0]).abs()).then(e1.cmp(e2)));
                    let mut touched = vec![false; tris.len()];
                    for (s0, w) in long {
                        let around = &by_edge[&(s0, w)];
                        if around.iter().any(|&t| touched[t]) { continue; }
                        let m = p2.len();
                        let mid = [(p2[s0][0] + p2[w][0]) / 2.0, (p2[s0][1] + p2[w][1]) / 2.0];
                        p2.push(mid);
                        let theta = mid[0] / big_r;
                        points3.push(p.add(a.scale(flip * mid[1])).add(u.scale(theta.cos())).add(v.scale(theta.sin())));
                        for &t in around {
                            touched[t] = true;
                            let tri = tris[t];
                            let k = (0..3).find(|&k| { let (b0, e0) = (tri[k], tri[(k + 1) % 3]); (b0 == s0 && e0 == w) || (b0 == w && e0 == s0) }).unwrap();
                            let (b0, e0, o) = (tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]);
                            tris[t] = [b0, m, o];
                            tris.push([m, e0, o]);
                            touched.push(true);
                        }
                    }
                    if _pass == 63 { return unsupported("an elliptic-cylinder face needs too many display triangles"); }
                }
                let normals = points3.iter().map(|q| {
                    let (t, _) = param(*q);
                    let tangent = u.scale(-t.sin()).add(v.scale(t.cos()));
                    let n = tangent.cross(a);
                    n.scale(side / n.length())
                }).collect();
                (tris, normals)
            }
        };
        let base = mesh.positions.len();
        mesh.positions.extend_from_slice(&points3);
        mesh.normals.extend(normals);
        for t in tris {
            mesh.triangles.push([base + t[0], base + t[1], base + t[2]]);
            mesh.faces.push(f);
        }
    }
    Ok(mesh)
}

/// A surface of revolution's chart (Claude 2026-10-29): theta round the axis, u along the meridian (a cone's height from
/// its point; a sphere's or torus's meridian angle).
#[derive(Clone, Copy)]
struct Revolved { o: Vec3, a: Vec3, x: Vec3, y: Vec3, cone: Option<(f64, f64)>, circle: Option<(f64, f64)>,
    mer: Option<crate::meridian::Meridian> }

impl Revolved {
    fn of(s: &Surface) -> Option<Revolved> {
        match *s {
            Surface::Cone { p, a, x, r, k } => Some(Revolved { o: p, a, x, y: a.cross(x), cone: Some((r, k)), circle: None, mer: None }),
            Surface::Sphere { c, a, x, r } => Some(Revolved { o: c, a, x, y: a.cross(x), cone: None, circle: Some((0.0, r)), mer: None }),
            Surface::Torus { c, a, x, big, r } => Some(Revolved { o: c, a, x, y: a.cross(x), cone: None, circle: Some((big, r)), mer: None }),
            Surface::Revolution { p, a, x, m } => Some(Revolved { o: p, a, x, y: a.cross(x), cone: None, circle: None, mer: Some(m) }),
            _ => None,
        }
    }
    fn radial(&self, q: Vec3) -> Vec3 { let d = q.sub(self.o); d.sub(self.a.scale(d.dot(self.a))) }
    fn theta(&self, q: Vec3) -> f64 { let d = q.sub(self.o); d.dot(self.y).atan2(d.dot(self.x)) }
    fn u(&self, q: Vec3) -> f64 {
        let h = q.sub(self.o).dot(self.a);
        if let Some(m) = self.mer { return m.param_of((self.radial(q).length(), h), 0.0); }
        match (self.cone, self.circle) {
            (Some(_), _) => h,
            (_, Some((big, _))) => h.atan2(self.radial(q).length() - big),
            _ => unreachable!(),
        }
    }
    fn radius(&self, u: f64) -> f64 {
        if let Some(m) = self.mer { return m.eval(u).0; }
        match (self.cone, self.circle) {
            (Some((r, k)), _) => r + k * u,
            (_, Some((big, rho))) => big + rho * u.cos(),
            _ => unreachable!(),
        }
    }
    fn rho(&self, t: f64) -> Vec3 { self.x.scale(t.cos()).add(self.y.scale(t.sin())) }
    fn point(&self, t: f64, u: f64) -> Vec3 {
        let along = match (self.cone, self.circle, self.mer) {
            (_, _, Some(m)) => m.eval(u).1,
            (Some(_), _, _) => u,
            (_, Some((_, rho)), _) => rho * u.sin(),
            _ => unreachable!(),
        };
        self.o.add(self.a.scale(along)).add(self.rho(t).scale(self.radius(u).max(0.0)))
    }
    /// The unit natural normal (away from the axis, or from the torus's centre circle).
    fn normal(&self, t: f64, u: f64) -> Vec3 {
        if let Some(m) = self.mer {
            let (dr, dz) = m.deriv(u);
            let l = dr.hypot(dz).max(1e-300);
            return self.rho(t).scale(dz / l).sub(self.a.scale(dr / l));
        }
        match (self.cone, self.circle) {
            (Some((_, k)), _) => self.rho(t).sub(self.a.scale(k)).scale(1.0 / (1.0 + k * k).sqrt()),
            (_, Some(_)) => self.rho(t).scale(u.cos()).add(self.a.scale(u.sin())),
            _ => unreachable!(),
        }
    }
    /// Distance from a point to the surface (in the half-plane through the axis: to the cone's line or the torus's
    /// meridian circle; exact).
    fn off(&self, q: Vec3) -> f64 {
        let (h, rad) = (q.sub(self.o).dot(self.a), self.radial(q).length());
        if let Some(m) = self.mer {
            return match m {
                crate::meridian::Meridian::Bezier(_) => m.distance((rad, h), 0.0, 1.0),
                crate::meridian::Meridian::Ellipse { .. } => m.distance((rad, h), 0.0, TAU),
            };
        }
        match (self.cone, self.circle) {
            (Some((r, k)), _) => (rad - (r + k * h)).abs() / (1.0 + k * k).sqrt(),
            (_, Some((big, rho))) => ((rad - big).hypot(h) - rho).abs(),
            _ => unreachable!(),
        }
    }
    /// Length along the surface per unit of u.
    fn u_scale(&self) -> f64 {
        if let Some(m) = self.mer {
            return match m {
                crate::meridian::Meridian::Bezier(_) => m.max_speed(0.0, 1.0),
                crate::meridian::Meridian::Ellipse { .. } => m.max_speed(0.0, TAU),
            }.max(1e-300);
        }
        match (self.cone, self.circle) { (Some((_, k)), _) => (1.0 + k * k).sqrt(), (_, Some((_, rho))) => rho, _ => unreachable!() }
    }
}

/// A face on a cone, sphere or torus as triangles in its (theta, u) chart. Loops are unwrapped by continuity; a point on
/// the axis (an apex or pole) has no angle of its own, so it is entered twice, at the angles of its neighbours (revolve
/// splits faces touching the axis so those neighbours' angles differ by at most half a turn). Inner edges longer than the
/// chord allows (round the axis at the smaller end's radius; along a curved meridian) are split and lifted onto the
/// surface. Returns the points, triangles and outward normals.
fn revolved_mesh(rev: &Revolved, loops3: &[Vec<Vec3>], side: f64, chord: f64, f: usize)
    -> Result<(Vec<Vec3>, Vec<[usize; 3]>, Vec<Vec3>), KernelError> {
    let scale = loops3.iter().flatten().map(|q| q.distance(rev.o)).fold(1.0, f64::max);
    let on_axis = |q: Vec3| rev.radial(q).length() <= 1e-10 * scale;
    let wraps_u = matches!(rev.circle, Some((big, _)) if big > 0.0) || matches!(rev.mer, Some(crate::meridian::Meridian::Ellipse { .. }));
    let mut points3 = Vec::new();
    let mut tu: Vec<(f64, f64)> = Vec::new();
    let mut index = Vec::new();
    let mut spans = Vec::new();
    let near = |prev: f64, raw: f64| prev + (raw - prev + PI).rem_euclid(TAU) - PI;
    for l in loops3 {
        let n = l.len();
        let Some(first) = (0..n).find(|&i| !on_axis(l[i])) else {
            return unsupported(format!("face {f}: a loop lies entirely on the axis"));
        };
        let mut idx = Vec::new();
        let (mut theta, mut u): (Option<f64>, Option<f64>) = (None, None);
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for step in 0..n {
            let q = l[(first + step) % n];
            let raw_u = rev.u(q);
            let uu = match u { Some(p) if wraps_u => near(p, raw_u), _ => raw_u };
            u = Some(uu);
            if on_axis(q) {
                let next = (1..n).map(|k| l[(first + step + k) % n]).find(|&w| !on_axis(w))
                    .ok_or_else(|| KernelError::Unsupported(format!("face {f}: an axis point has no neighbour off the axis")))?;
                let before = theta.expect("the loop starts off the axis");
                // The jump along the pole: its size is the angle to the next meridian (at most half a turn; exactly half
                // is ambiguous by angle alone), its direction the loop's turn: backwards at the top of the chart (greatest
                // u), forwards at the bottom, reversed for an inward face.
                let turn = near(before, rev.theta(next)) - before;
                let size = turn.abs();
                let prev_u = tu.last().expect("the loop starts off the axis").1;
                let dir = if uu > prev_u { -1.0 } else { 1.0 } * side;
                let after = before + dir * size;
                // Entered along the pole every 1/32 turn (one point in space), so each triangle meets the pole at an
                // angle near its own: a triangle reaching a far entry would fold over the pole.
                let steps = ((size / (PI / 32.0)).ceil() as usize).max(1);
                for k in 0..=steps {
                    idx.push(points3.len());
                    points3.push(q);
                    tu.push((before + (after - before) * k as f64 / steps as f64, uu));
                }
                theta = Some(after);
                continue;
            }
            let t = match theta { None => rev.theta(q), Some(p) => near(p, rev.theta(q)) };
            theta = Some(t);
            lo = lo.min(t);
            hi = hi.max(t);
            idx.push(points3.len());
            points3.push(q);
            tu.push((t, uu));
        }
        // A loop must come back to its start angle (a seam closes a whole turn); one that wraps without a seam is refused.
        let (start, last) = (tu[idx[0]].0, theta.unwrap());
        if (last - start).abs() > PI {
            return unsupported("a revolved face loop wraps the axis without a seam edge");
        }
        spans.push((lo, hi));
        index.push(idx);
    }
    // Holes unwrapped a turn away from the outer loop are moved into its range.
    if let Some(k) = (0..index.len()).max_by(|&i, &j| (spans[i].1 - spans[i].0).total_cmp(&(spans[j].1 - spans[j].0))) {
        let (lo, hi) = spans[k];
        for (j, l) in index.iter().enumerate() {
            if j == k { continue; }
            let mid = (spans[j].0 + spans[j].1) / 2.0;
            let turns = if mid < lo { ((lo - mid) / TAU).ceil() } else if mid > hi { -((mid - hi) / TAU).ceil() } else { 0.0 };
            for &i in l { tu[i].0 += turns * TAU; }
        }
    }
    let big_r = tu.iter().map(|&(_, u)| rev.radius(u)).fold(0.0, f64::max);
    if !(big_r > 0.0) { return unsupported(format!("revolved face {f} has no positive radius")); }
    // A cone's generators are straight, so its chart is squashed along them (an affine change: triangulation stays
    // valid) until it is no taller than a tenth of its width: triangles then follow the generators instead of crossing
    // the angle on long diagonals that would need endless splitting.
    let (u_lo, u_hi) = tu.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &(_, u)| (lo.min(u), hi.max(u)));
    let ys = match rev.cone {
        Some(_) => rev.u_scale() * (0.1 * big_r * TAU / (rev.u_scale() * (u_hi - u_lo)).max(1e-300)).min(1.0),
        None => rev.u_scale(),
    };
    let flip = side;
    let mut p2: Vec<[f64; 2]> = tu.iter().map(|&(t, u)| [big_r * t, flip * ys * u]).collect();
    let mut tris = triangulate_runs(&p2, &index, true, if rev.cone.is_some() { Runs::Axis } else { Runs::None })?;
    // Refinement by true distance from the surface (exact for these surfaces): an inner edge whose midpoint departs by
    // more than 0.6 of the chord is split, and so is the longest inner edge of a triangle whose centre departs by more than
    // 0.95 of it. A triangle's centre lies at most 4/9 of its edges' midpoint departures below them, so with boundary
    // edges within the chord (their polylines) every face ends within it. Edges from an apex along a cone are generators
    // (no departure) and are never split.
    let limit_edge = 0.6 * chord;
    let limit_centre = 0.95 * chord;
    let boundary: std::collections::HashSet<(usize, usize)> = index.iter()
        .flat_map(|l| (0..l.len()).map(move |j| { let (s0, w) = (l[j], l[(j + 1) % l.len()]); (s0.min(w), s0.max(w)) }))
        .collect();
    for pass in 0..100 {
        delaunay_flips(&p2, &mut tris, &boundary);
        let mut by_edge: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (t, tri) in tris.iter().enumerate() {
            for j in 0..3 { let (s0, w) = (tri[j], tri[(j + 1) % 3]); by_edge.entry((s0.min(w), s0.max(w))).or_default().push(t); }
        }
        let inner_edge = |e: &(usize, usize)| !boundary.contains(e) && points3[e.0] != points3[e.1];
        let mut split: Vec<(usize, usize)> = by_edge.keys().copied()
            .filter(|e| inner_edge(e) && rev.off(points3[e.0].add(points3[e.1]).scale(0.5)) > limit_edge).collect();
        for tri in &tris {
            let [pa, pb, pc] = tri.map(|i| points3[i]);
            // A triangle facing against the surface folds the chart near a pole: split it in the chart until it does not.
            if pa != pb && pb != pc && pc != pa {
                let normal = tri.iter().fold(Vec3::new(0.0, 0.0, 0.0), |acc, &i| acc.add(rev.normal(tu[i].0, tu[i].1))).scale(side);
                if pb.sub(pa).cross(pc.sub(pa)).dot(normal) <= 0.0 {
                    let chart = |e: &(usize, usize)| (p2[e.0][0] - p2[e.1][0]).hypot(p2[e.0][1] - p2[e.1][1]);
                    let longest = (0..3).map(|j| { let (u0, u1) = (tri[j], tri[(j + 1) % 3]); (u0.min(u1), u0.max(u1)) }).filter(|e| inner_edge(e))
                        .max_by(|e1, e2| chart(e1).total_cmp(&chart(e2)).then(e1.cmp(e2)));
                    if let Some(e) = longest { split.push(e); continue; }
                }
            }
            // Triangles on two entries of one axis point are dropped at the end; their edges are checked above.
            if pa == pb || pb == pc || pc == pa { continue; }
            let centre = pa.add(pb).add(pc).scale(1.0 / 3.0);
            if rev.off(centre) <= limit_centre { continue; }
            let longest = (0..3).map(|j| { let (a, b) = (tri[j], tri[(j + 1) % 3]); (a.min(b), a.max(b)) }).filter(|e| inner_edge(e))
                .max_by(|e1, e2| points3[e1.0].distance(points3[e1.1]).total_cmp(&points3[e2.0].distance(points3[e2.1])).then(e1.cmp(e2)));
            if let Some(e) = longest { split.push(e); }
        }
        split.sort();
        split.dedup();
        if split.is_empty() { break; }
        if pass == 99 || tris.len() > 2_000_000 { return unsupported(format!("revolved face {f} needs too many display triangles")); }
        let len = |e: &(usize, usize)| points3[e.0].distance(points3[e.1]);
        split.sort_by(|e1, e2| len(e2).total_cmp(&len(e1)).then(e1.cmp(e2)));
        let mut touched = vec![false; tris.len()];
        for (s0, w) in split {
            let around = &by_edge[&(s0, w)];
            if around.iter().any(|&t| touched[t]) { continue; }
            let m = p2.len();
            let mid = [(p2[s0][0] + p2[w][0]) / 2.0, (p2[s0][1] + p2[w][1]) / 2.0];
            p2.push(mid);
            let (t, u) = (mid[0] / big_r, mid[1] / (flip * ys));
            points3.push(rev.point(t, u));
            tu.push((t, u));
            for &tr in around {
                touched[tr] = true;
                let tri = tris[tr];
                let j = (0..3).find(|&j| { let (b0, e0) = (tri[j], tri[(j + 1) % 3]); (b0 == s0 && e0 == w) || (b0 == w && e0 == s0) }).unwrap();
                let (b0, e0, o) = (tri[j], tri[(j + 1) % 3], tri[(j + 2) % 3]);
                tris[tr] = [b0, m, o];
                tris.push([m, e0, o]);
                touched.push(true);
            }
        }
    }
    // The two entries of an axis point are one point in space: a triangle using both has no area and is dropped (its other
    // two edges are one segment, shared by the triangles either side).
    tris.retain(|t| { let [a, b, c] = t.map(|i| points3[i]); a != b && b != c && c != a });
    let normals = tu.iter().map(|&(t, u)| rev.normal(t, u).scale(side)).collect();
    Ok((points3, tris, normals))
}

/// Lawson flips in the chart: an inner edge whose opposite vertex lies inside the circumcircle of the triangle across it
/// (and whose quadrilateral is convex) is replaced by the other diagonal, until none is left. Removes slivers the
/// splitting cannot (a chain of points along one meridian). Boundary edges never move.
fn delaunay_flips(p2: &[[f64; 2]], tris: &mut [[usize; 3]], boundary: &std::collections::HashSet<(usize, usize)>) {
    let in_circle = |a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]| {
        let (adx, ady, bdx, bdy, cdx, cdy) = (a[0] - d[0], a[1] - d[1], b[0] - d[0], b[1] - d[1], c[0] - d[0], c[1] - d[1]);
        (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) - (bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady)
            + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady)
    };
    for _ in 0..64 {
        let mut by_edge: HashMap<(usize, usize), Vec<(usize, usize)>> = HashMap::new();
        for (t, tri) in tris.iter().enumerate() {
            for j in 0..3 { let (a, b) = (tri[j], tri[(j + 1) % 3]); by_edge.entry((a.min(b), a.max(b))).or_default().push((t, j)); }
        }
        let mut flipped = vec![false; tris.len()];
        let mut any = false;
        let mut keys: Vec<(usize, usize)> = by_edge.keys().copied().collect();
        keys.sort();
        for e in keys {
            let around = &by_edge[&e];
            if boundary.contains(&e) || around.len() != 2 { continue; }
            let ((t1, j1), (t2, j2)) = (around[0], around[1]);
            if flipped[t1] || flipped[t2] { continue; }
            let (a, b, c) = (tris[t1][j1], tris[t1][(j1 + 1) % 3], tris[t1][(j1 + 2) % 3]);
            let d = tris[t2][(j2 + 2) % 3];
            // t1 = (a, b, c) counter-clockwise; t2 runs b -> a with d opposite.
            if tris[t2][j2] != b || tris[t2][(j2 + 1) % 3] != a { continue; }
            if in_circle(p2[a], p2[b], p2[c], p2[d]) <= 1e-12 * (cross2(p2[a], p2[b], p2[c]).abs() + 1e-300).sqrt().max(1e-300) { continue; }
            // The new triangles (c, a, d) and (d, b, c) must both turn counter-clockwise (a convex quadrilateral).
            if cross2(p2[c], p2[a], p2[d]) <= 0.0 || cross2(p2[d], p2[b], p2[c]) <= 0.0 { continue; }
            tris[t1] = [c, a, d];
            tris[t2] = [d, b, c];
            flipped[t1] = true;
            flipped[t2] = true;
            any = true;
        }
        if !any { break; }
    }
}

/// The worker reply: [1, chord, V, T, E, positions 3V, normals 3V, triangles 4T (i, j, k, face), then per edge:
/// original body edge index (possibly sparse), point count n, 3n points]. Only display edges are exported;
/// triangles and the exact body's edges remain unchanged.
pub fn write_mesh(m: &Mesh) -> Vec<f64> {
    let mut out = vec![1.0, m.chord, m.positions.len() as f64, m.triangles.len() as f64, m.display_edges.len() as f64];
    for v in m.positions.iter().chain(&m.normals) { out.extend_from_slice(&[v.x, v.y, v.z]); }
    for (t, f) in m.triangles.iter().zip(&m.faces) { out.extend_from_slice(&[t[0] as f64, t[1] as f64, t[2] as f64, *f as f64]); }
    for &e in &m.display_edges {
        let line = &m.edges[e];
        out.push(e as f64);
        out.push(line.len() as f64);
        for v in line { out.extend_from_slice(&[v.x, v.y, v.z]); }
    }
    out
}

#[cfg(test)]
mod edge_classification_tests {
    use super::*;
    use crate::neutral::Face;

    #[test]
    fn support_identity_preserves_shallow_creases_offsets_and_fillet_boundaries() {
        let z = Vec3::new(0.0, 0.0, 1.0);
        let x = Vec3::new(1.0, 0.0, 0.0);
        let o = Vec3::new(0.0, 0.0, 0.0);
        let face = |surface| Face { surface, sense: true, loops: Vec::new() };
        let p = face(Surface::Plane { p: o, n: z, x });
        let shifted = face(Surface::Plane { p: x.scale(1e8), n: z, x });
        assert!(same_support(&p, &shifted));
        let reversed = Face { surface: Surface::Plane { p: o, n: z.scale(-1.0), x }, sense: false, loops: Vec::new() };
        assert!(same_support(&p, &reversed));
        assert!(!same_support(&p, &Face { sense: false, ..p.clone() }));
        assert!(!same_support(&p, &face(Surface::Plane { p: o, n: Vec3::new(1e-8, 0.0, 1.0).unit().unwrap(), x })));
        assert!(!same_support(&p, &face(Surface::Plane { p: z.scale(1e-8), n: z, x })));
        let c = face(Surface::Cylinder { p: o, a: z, x, r: 1.0 });
        assert!(same_support(&c, &face(Surface::Cylinder { p: z.scale(1000.0), a: z.scale(-1.0), x, r: 1.0 })));
        assert!(!same_support(&c, &Face { sense: false, ..c.clone() }));
        assert!(!same_support(&c, &face(Surface::Cylinder { p: x.scale(1e-7), a: z, x, r: 1.0 })));
        assert!(!same_support(&c, &face(Surface::Cylinder { p: o, a: z, x, r: 1.0 + 1e-7 })));
        assert!(!same_support(&p, &c), "fillet transitions are intentional edges");
    }
}

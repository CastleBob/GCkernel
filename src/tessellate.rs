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
fn angle_step(r: f64, chord: f64) -> f64 {
    if chord >= r { PI / 2.0 } else { (2.0 * (1.0 - chord / r).acos()).min(PI / 2.0) }
}

/// An edge as points from its start vertex to its end vertex.
fn polyline(body: &Body, edge: usize, chord: f64) -> Result<Vec<Vec3>, KernelError> {
    let e = &body.edges[edge];
    let (t0, t1) = e.curve.domain();
    let at = |t: f64| e.curve.eval(t).0;
    let params: Vec<f64> = match e.curve {
        Curve::Line { .. } => vec![0.0, 1.0],
        Curve::Arc { r, .. } => {
            let n = ((t1 - t0) / angle_step(r, chord)).ceil().max(1.0) as usize;
            (0..=n).map(|k| t0 + (t1 - t0) * k as f64 / n as f64).collect()
        }
        Curve::Conic { u, v, .. } => {
            // |P''| <= sqrt(|u|^2 + |v|^2); a chord over step h departs by at most h^2 |P''| / 8.
            let bend = (u.dot(u) + v.dot(v)).sqrt();
            let step = (8.0 * chord / bend).sqrt().min(PI / 2.0);
            let n = ((t1 - t0) / step).ceil().max(1.0) as usize;
            (0..=n).map(|k| t0 + (t1 - t0) * k as f64 / n as f64).collect()
        }
        Curve::CylinderPair { .. } => {
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
fn clip(p: &[[f64; 2]], mut poly: Vec<usize>, short: bool, out: &mut Vec<[usize; 3]>) -> Result<(), KernelError> {
    let scale = poly.iter().map(|&i| p[i][0].abs().max(p[i][1].abs())).fold(1.0, f64::max);
    let eps = 1e-14 * scale * scale;
    let mut i = 0;
    let mut stall = 0;
    while poly.len() > 3 {
        let n = poly.len();
        if short {
            let side = |k: usize| { let (a, c) = (p[poly[(k + n - 1) % n]], p[poly[(k + 1) % n]]); (a[0] - c[0]).hypot(a[1] - c[1]) };
            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by(|&x, &y| side(x).total_cmp(&side(y)).then(x.cmp(&y)));
            i = order.into_iter().find(|&k| is_ear(p, &poly, k, eps, short)).ok_or_else(|| KernelError::Unsupported(
                "a face could not be triangulated for display".into()))?;
        }
        let (ia, ib, ic) = ((i + n - 1) % n, i % n, (i + 1) % n);
        if is_ear(p, &poly, ib, eps, short) {
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
fn is_ear(p: &[[f64; 2]], poly: &[usize], ib: usize, eps: f64, curved: bool) -> bool {
    let n = poly.len();
    let (ia, ic) = ((ib + n - 1) % n, (ib + 1) % n);
    let (a, b, c) = (p[poly[ia]], p[poly[ib]], p[poly[ic]]);
    let turn = cross2(a, b, c);
    if turn.abs() <= eps {
        if curved && (a[0] != b[0] || b[0] != c[0]) { return false; }
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
        clip(p, poly, short, &mut out)?;
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

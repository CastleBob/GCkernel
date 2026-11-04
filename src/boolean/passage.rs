//! Local admission for a round through-passage in an otherwise unchanged solid.
//! Proves planar entry/exit patch pairs contain the disk with clearance, and
//! conservatively excludes EVERY other surface patch from its finite corridor.
//! Bounds are analytic; rejecting an inconclusive bound never changes geometry.
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Face};
use crate::KernelError;
use std::f64::consts::{PI, TAU};

#[cfg(test)]
#[path = "passage_tests.rs"]
mod tests;

type Result<T> = std::result::Result<T, KernelError>;
fn bad(s: &str) -> KernelError {
    KernelError::InvalidInput(s.into())
}
fn unsupported(s: &str) -> KernelError {
    KernelError::Unsupported(s.into())
}
fn xyz(p: Vec3) -> [f64; 3] {
    [p.x, p.y, p.z]
}
fn vector(p: [f64; 3]) -> Vec3 {
    Vec3::new(p[0], p[1], p[2])
}
fn unit(p: Vec3) -> bool {
    p.is_finite() && (p.length() - 1.).abs() < 1e-10
}
fn sign(forward: bool) -> f64 {
    if forward {
        1.
    } else {
        -1.
    }
}

fn harmonic(curve: &Curve) -> Option<(Vec3, Vec3, Vec3, f64, f64)> {
    match *curve {
        Curve::Line { .. } | Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } | Curve::Bezier { .. } | Curve::RationalBezier { .. } => None,
        Curve::Arc { c, n, x, r, t0, t1 } => Some((c, x.scale(r), n.cross(x).scale(r), t0, t1)),
        Curve::Conic { c, u, v, t0, t1 } => Some((c, u, v, t0, t1)),
    }
}

fn parameters(base: f64, lo: f64, hi: f64) -> Vec<f64> {
    let first = base + ((lo - base) / TAU).ceil() * TAU;
    [first, first + TAU]
        .into_iter()
        .filter(|t| *t >= lo && *t <= hi)
        .collect()
}

pub(super) fn curve_bounds(curve: &Curve) -> ([f64; 3], [f64; 3]) {
    let mut points = vec![curve.start(), curve.end()];
    if let Some((c, u, v, t0, t1)) = harmonic(curve) {
        let a = xyz(u);
        let b = xyz(v);
        for k in 0..3 {
            for root in [b[k].atan2(a[k]), b[k].atan2(a[k]) + PI] {
                for t in parameters(root, t0, t1) {
                    points.push(Curve::conic_point(c, u, v, t));
                }
            }
        }
    }
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in points {
        let p = xyz(p);
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    (lo, hi)
}

fn distance(curve: &Curve, q: Vec3) -> f64 {
    match *curve {
        Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } | Curve::Bezier { .. } | Curve::RationalBezier { .. } => unreachable!("pair edges are refused by passage admission"),
        Curve::Conic { .. } => {
            // A lower distance bound is sufficient for clearance admission.
            // An inconclusive box rejects the cut; it never substitutes a curve.
            let (lo, hi) = curve_bounds(curve);
            let point = xyz(q);
            (0..3)
                .map(|k| (lo[k] - point[k]).max(point[k] - hi[k]).max(0.).powi(2))
                .sum::<f64>()
                .sqrt()
        }
        Curve::Line { a, b } => {
            let d = b.sub(a);
            let t = (q.sub(a).dot(d) / d.dot(d)).clamp(0., 1.);
            q.distance(a.add(d.scale(t)))
        }
        Curve::Arc { c, n, x, r, t0, t1 } => {
            let d = q.sub(c);
            let t = d.dot(n.cross(x)).atan2(d.dot(x));
            let mut closest = q.distance(curve.start()).min(q.distance(curve.end()));
            if !parameters(t, t0, t1).is_empty() {
                let radial = d.sub(n.scale(d.dot(n))).length();
                closest = closest.min(((radial - r).powi(2) + d.dot(n).powi(2)).sqrt());
            }
            closest
        }
    }
}

fn loop_area(body: &Body, lp: &[(usize, bool)], n: Vec3) -> f64 {
    lp.iter()
        .map(|(i, fwd)| {
            let c = &body.edges[*i].curve;
            let integral = match *c {
                Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } | Curve::Bezier { .. } | Curve::RationalBezier { .. } => {
                    unreachable!("pair edges are refused by passage admission")
                }
                Curve::Conic {
                    c: center,
                    u,
                    v,
                    t0,
                    t1,
                } => center
                    .cross(c.end().sub(c.start()))
                    .add(u.cross(v).scale(t1 - t0)),
                Curve::Line { a, b } => a.cross(b),
                Curve::Arc {
                    c: center,
                    n,
                    r,
                    t0,
                    t1,
                    ..
                } => center
                    .cross(c.end().sub(c.start()))
                    .add(n.scale(r * r * (t1 - t0))),
            };
            0.5 * integral.dot(n) * sign(*fwd)
        })
        .sum()
}

/// Exact line/circle ray crossings. Choose a ray that avoids authored vertices
/// so no endpoint ownership is guessed from a near-equality comparison.
fn inside(body: &Body, face: &Face, q: Vec3, n: Vec3, tol: f64) -> Result<bool> {
    let x = n.perpendicular();
    let y = n.cross(x);
    let mut frame = None;
    for k in 0..16 {
        let a = k as f64 * 0.3713906763541037;
        let u = x.scale(a.cos()).add(y.scale(a.sin()));
        let v = n.cross(u);
        if face.loops.iter().flatten().all(|(i, _)| {
            let c = &body.edges[*i].curve;
            c.start().sub(q).dot(v).abs() > tol && c.end().sub(q).dot(v).abs() > tol
        }) {
            frame = Some((u, v));
            break;
        }
    }
    let (u, v) = frame.ok_or_else(|| {
        unsupported("no unambiguous passage classification ray within the work limit")
    })?;
    let mut winding = 0i32;
    for (i, fwd) in face.loops.iter().flatten() {
        let source = &body.edges[*i].curve;
        let c = if *fwd {
            source.clone()
        } else {
            source.reversed()
        };
        match c {
            Curve::Line { a, b } => {
                let ay = a.sub(q).dot(v);
                let by = b.sub(q).dot(v);
                if (ay < 0.) != (by < 0.) {
                    let at = -ay / (by - ay);
                    let hit = a.add(b.sub(a).scale(at));
                    if hit.sub(q).dot(u) > 0. {
                        winding += if by > ay { 1 } else { -1 };
                    }
                }
            }
            _ => {
                let (c, cu, cv, t0, t1) = harmonic(&c).unwrap();
                let a = cu.dot(v);
                let b = cv.dot(v);
                let height = c.sub(q).dot(v);
                let amplitude = a.hypot(b);
                let ratio = -height / amplitude;
                if ratio.abs() >= 1. {
                    continue;
                }
                let phase = b.atan2(a);
                let angle = ratio.acos();
                for root in [phase - angle, phase + angle] {
                    for t in parameters(root, t0, t1) {
                        let hit = Curve::conic_point(c, cu, cv, t);
                        if hit.sub(q).dot(u) > 0. {
                            winding += if -a * t.sin() + b * t.cos() > 0. {
                                1
                            } else {
                                -1
                            };
                        }
                    }
                }
            }
        }
    }
    if !(0..=1).contains(&winding) {
        return Err(bad("passage face loops have inconsistent winding"));
    }
    Ok(winding == 1)
}

fn validate(body: &Body, tol: f64) -> Result<Vec<([f64; 3], [f64; 3])>> {
    if body.faces.len() > 20000 || body.edges.len() > 80000 || body.vertices.len() > 80000 {
        return Err(unsupported("passage input exceeds the topology work limit"));
    }
    if body.vertices.iter().any(|p| !p.is_finite())
        || body
            .edges
            .iter()
            .any(|e| e.start >= body.vertices.len() || e.end >= body.vertices.len())
        || body
            .faces
            .iter()
            .flat_map(|f| &f.loops)
            .flatten()
            .any(|(i, _)| *i >= body.edges.len())
    {
        return Err(bad(
            "passage input has invalid coordinates or topology references",
        ));
    }
    for edge in &body.edges {
        match edge.curve {
            Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } | Curve::Bezier { .. } | Curve::RationalBezier { .. } => {
                return Err(unsupported(
                    "cylinder-pair trims require transverse admission",
                ))
            }
            Curve::Conic { c, u, v, t0, t1 } => {
                if !c.is_finite()
                    || !u.is_finite()
                    || !v.is_finite()
                    || u.length() <= tol
                    || v.length() <= tol
                    || u.cross(v).length() <= 1e-10 * u.length() * v.length()
                    || !t0.is_finite()
                    || !t1.is_finite()
                    || t1 <= t0
                    || t1 - t0 > TAU + 1e-12
                    || t0.abs() > 1e5
                    || t1.abs() > 1e5
                {
                    return Err(bad("invalid passage conic frame"));
                }
            }
            Curve::Line { a, b } => {
                if !a.is_finite() || !b.is_finite() || a.distance(b) <= tol {
                    return Err(bad("invalid passage boundary line"));
                }
            }
            Curve::Arc { c, n, x, r, t0, t1 } => {
                if !c.is_finite()
                    || !unit(n)
                    || !unit(x)
                    || n.dot(x).abs() > 1e-10
                    || !r.is_finite()
                    || r <= tol
                    || !t0.is_finite()
                    || !t1.is_finite()
                    || t1 <= t0
                    || t1 - t0 > TAU + 1e-12
                {
                    return Err(bad("invalid passage boundary arc"));
                }
            }
        }
    }
    crate::topology::check(body, tol)?;
    let mut bounds = Vec::new();
    for face in &body.faces {
        let (p, n, x) = match face.surface {
            Surface::Plane { p, n, x } => (p, n, x),
            Surface::Cylinder { p, a, x, .. } => (p, a, x),
            Surface::EllipticCylinder { .. } => return Err(unsupported("passages through bodies with elliptic-cylinder faces are not supported")),
            Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => return Err(unsupported("passages through bodies with revolved or spline faces are not supported")),
        };
        if !p.is_finite() || !unit(n) || !unit(x) || n.dot(x).abs() > 1e-10 {
            return Err(bad("invalid passage surface frame"));
        }
        let mut low = [f64::INFINITY; 3];
        let mut high = [f64::NEG_INFINITY; 3];
        for (i, _) in face.loops.iter().flatten() {
            let (lo, hi) = curve_bounds(&body.edges[*i].curve);
            for k in 0..3 {
                low[k] = low[k].min(lo[k]);
                high[k] = high[k].max(hi[k]);
            }
        }
        match face.surface {
            Surface::Plane { p, n, .. } => {
                for (i, _) in face.loops.iter().flatten() {
                    let curve = &body.edges[*i].curve;
                    if curve.start().sub(p).dot(n).abs() > tol
                        || curve.end().sub(p).dot(n).abs() > tol
                    {
                        return Err(bad("passage edge is off its plane"));
                    }
                    if let Some((c, u, v, _, _)) = harmonic(curve) {
                        if c.sub(p).dot(n).abs() + u.dot(n).abs() + v.dot(n).abs() > tol {
                            return Err(bad("passage conic is not in its plane"));
                        }
                    }
                }
                let out = n.scale(sign(face.sense));
                for (k, lp) in face.loops.iter().enumerate() {
                    let area = loop_area(body, lp, out);
                    if !area.is_finite()
                        || if k == 0 {
                            area <= tol * tol
                        } else {
                            area >= -tol * tol
                        }
                    {
                        return Err(bad("passage planar loop orientation is invalid"));
                    }
                }
            }
            Surface::EllipticCylinder { .. } => return Err(unsupported("passages through bodies with elliptic-cylinder faces are not supported")),
            Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => return Err(unsupported("passages through bodies with revolved or spline faces are not supported")),
            Surface::Cylinder { p, a, r, .. } => {
                if !r.is_finite() || r <= tol || face.loops.len() != 1 || face.loops[0].len() < 4 {
                    return Err(unsupported(
                        "passage requires simple axial cylindrical patches",
                    ));
                }
                let mut arcs = Vec::new();
                let mut lines = 0;
                for (i, fwd) in &face.loops[0] {
                    match body.edges[*i].curve {
                        Curve::Line { a: from, b: to } => {
                            if to.sub(from).cross(a).length() > tol
                                || (from.sub(p).cross(a).length() - r).abs() > tol
                            {
                                return Err(bad("cylindrical generator is off its surface"));
                            }
                            lines += 1;
                        }
                        _ => {
                            let (c, u, v, t0, t1) = harmonic(&body.edges[*i].curve).unwrap();
                            let radial = |q: Vec3| q.sub(a.scale(q.dot(a)));
                            let q = radial(c.sub(p));
                            let ur = radial(u);
                            let vr = radial(v);
                            // Exact harmonic coefficients bound squared radial error
                            // for the entire trim, not merely at selected samples.
                            let error = (q.dot(q) + (ur.dot(ur) + vr.dot(vr)) * 0.5 - r * r).abs()
                                + 2. * q.dot(ur).abs()
                                + 2. * q.dot(vr).abs()
                                + ((ur.dot(ur) - vr.dot(vr)) * 0.5).abs()
                                + ur.dot(vr).abs();
                            if error > tol * (2. * r + tol)
                                || q.length() > tol
                                || (ur.length() - r).abs() > tol
                                || (vr.length() - r).abs() > tol
                                || ur.dot(vr).abs() > tol * r
                            {
                                return Err(bad("cylindrical trim is off its surface"));
                            }
                            let axial = |t: f64| Curve::conic_point(c, u, v, t).dot(a);
                            let mut low = axial(t0).min(axial(t1));
                            let mut high = axial(t0).max(axial(t1));
                            let phase = v.dot(a).atan2(u.dot(a));
                            for angle in [phase, phase + PI] {
                                for t in parameters(angle, t0, t1) {
                                    low = low.min(axial(t));
                                    high = high.max(axial(t));
                                }
                            }
                            let direction = ur.cross(vr).dot(a) / (r * r) * sign(*fwd);
                            let midpoint =
                                radial(Curve::conic_point(c, u, v, (t0 + t1) * 0.5).sub(p));
                            arcs.push((low, high, direction, midpoint, t1 - t0));
                        }
                    }
                }
                if arcs.len() != 2 || lines < 2 {
                    return Err(unsupported(
                        "cylinder patch must have two conic trims and axial generators",
                    ));
                }
                arcs.sort_by(|a, b| a.0.total_cmp(&b.0));
                if arcs[1].0 - arcs[0].1 <= tol
                    || arcs[0].2 * sign(face.sense) < 0.999999
                    || arcs[1].2 * sign(face.sense) > -0.999999
                {
                    return Err(bad(
                        "cylinder patch axial ordering or orientation is invalid",
                    ));
                }
                let first = &arcs[0];
                let second = &arcs[1];
                if (first.4 - second.4).abs() > 1e-10 || second.3.distance(first.3) > tol {
                    return Err(unsupported("cylinder trims do not bound one axial strip"));
                }
            }
        }
        bounds.push((low, high));
    }
    Ok(bounds)
}

/// Join connected patches on exactly the same axis-normal plane. This removes
/// artificial partition seams, not geometric boundaries: only a shared edge
/// with opposite occurrences on two equal outward planes can be cancelled.
pub(super) fn coalesce_caps(body: &Body, axis: usize, tol: f64) -> Result<Body> {
    validate(body, tol)?;
    let planes: Vec<_> = body
        .faces
        .iter()
        .map(|f| match f.surface {
            Surface::Plane { p, n, .. }
                if xyz(n)[axis].abs() == 1. && (0..3).all(|k| k == axis || xyz(n)[k] == 0.) =>
            {
                Some((xyz(p)[axis], xyz(n)[axis] * sign(f.sense)))
            }
            _ => None,
        })
        .collect();
    let mut owners = vec![Vec::new(); body.edges.len()];
    for (i, face) in body.faces.iter().enumerate() {
        for (edge, forward) in face.loops.iter().flatten() {
            owners[*edge].push((i, *forward));
        }
    }
    let mut adjacency = vec![Vec::new(); body.faces.len()];
    for pair in &owners {
        let [(a, af), (b, bf)] = pair[..] else {
            return Err(bad("cap union requires two edge occurrences"));
        };
        if a != b && af != bf && planes[a].is_some() && planes[a] == planes[b] {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
    }
    let mut visited = vec![false; body.faces.len()];
    let mut membership = vec![usize::MAX; body.faces.len()];
    let mut faces = Vec::new();
    let mut changed = false;
    for seed in 0..body.faces.len() {
        if visited[seed] {
            continue;
        }
        let mut component = vec![seed];
        visited[seed] = true;
        let mut at = 0;
        while at < component.len() {
            for &next in &adjacency[component[at]] {
                if !visited[next] {
                    visited[next] = true;
                    component.push(next);
                }
            }
            at += 1;
        }
        if component.len() == 1 {
            faces.push(body.faces[seed].clone());
            continue;
        }
        changed = true;
        for &i in &component {
            membership[i] = seed;
        }
        let mut boundary = Vec::new();
        for &i in &component {
            for &(edge, forward) in body.faces[i].loops.iter().flatten() {
                if owners[edge]
                    .iter()
                    .all(|(face, _)| membership[*face] == seed)
                {
                    continue;
                }
                boundary.push((edge, forward));
            }
        }
        let mut starts = std::collections::BTreeMap::new();
        let mut ends = std::collections::BTreeSet::new();
        for (i, &(edge, forward)) in boundary.iter().enumerate() {
            let e = &body.edges[edge];
            let (start, end) = if forward {
                (e.start, e.end)
            } else {
                (e.end, e.start)
            };
            if starts.insert(start, i).is_some() || !ends.insert(end) {
                return Err(unsupported("coplanar union boundary branches at a vertex"));
            }
        }
        let mut seen = vec![false; boundary.len()];
        let mut loops = Vec::new();
        for first in 0..boundary.len() {
            if seen[first] {
                continue;
            }
            let mut lp = Vec::new();
            let mut i = first;
            loop {
                if seen[i] {
                    return Err(bad("coplanar union boundary closes on another loop"));
                }
                seen[i] = true;
                let (edge, forward) = boundary[i];
                lp.push((edge, forward));
                let e = &body.edges[edge];
                let end = if forward { e.end } else { e.start };
                i = *starts
                    .get(&end)
                    .ok_or_else(|| bad("coplanar union boundary is open"))?;
                if i == first {
                    break;
                }
            }
            loops.push(lp);
        }
        let mut n = [0.; 3];
        n[axis] = planes[seed].unwrap().1;
        let n = vector(n);
        let outer: Vec<_> = loops
            .iter()
            .enumerate()
            .filter_map(|(i, lp)| (loop_area(body, lp, n) > tol * tol).then_some(i))
            .collect();
        if outer.len() != 1
            || loops
                .iter()
                .any(|lp| loop_area(body, lp, n).abs() <= tol * tol)
        {
            return Err(unsupported(
                "coplanar union needs one outer boundary and resolved holes",
            ));
        }
        loops.swap(0, outer[0]);
        let mut face = body.faces[seed].clone();
        face.loops = loops;
        faces.push(face);
    }
    if !changed {
        return Ok(body.clone());
    }
    // Remove only unused seam edges/vertices and remap surviving references.
    let mut edge_used = vec![false; body.edges.len()];
    for face in &faces {
        for (edge, _) in face.loops.iter().flatten() {
            edge_used[*edge] = true;
        }
    }
    let mut vertex_used = vec![false; body.vertices.len()];
    for (i, edge) in body.edges.iter().enumerate() {
        if edge_used[i] {
            vertex_used[edge.start] = true;
            vertex_used[edge.end] = true;
        }
    }
    let mut result = Body::default();
    let mut vm = vec![0; body.vertices.len()];
    for (i, vertex) in body.vertices.iter().enumerate() {
        if vertex_used[i] {
            vm[i] = result.vertices.len();
            result.vertices.push(*vertex);
        }
    }
    let mut em = vec![0; body.edges.len()];
    for (i, edge) in body.edges.iter().enumerate() {
        if edge_used[i] {
            em[i] = result.edges.len();
            let mut edge = edge.clone();
            edge.start = vm[edge.start];
            edge.end = vm[edge.end];
            result.edges.push(edge);
        }
    }
    for mut face in faces {
        for (edge, _) in face.loops.iter_mut().flatten() {
            *edge = em[*edge];
        }
        result.faces.push(face);
    }
    validate(&result, tol)?;
    Ok(result)
}

/// Seams on a cap plane (Bob's angle corner cut, 2026-11-03): earlier transverse cuts leave a cap split into faces on
/// one plane with one outward side, meeting across edges inside the material. Those edges are not cap boundaries. The
/// axis-normal faces the corridor reaches are grouped by plane and side and joined across the edges they share (exact
/// topology edits, unite_tidy: no point moved, no surface fitted). `None` when the corridor meets no such seam. Decided
/// from the operands before admission; the joined body is then admitted as any other.
pub(super) fn join_cap_seams(
    body: &Body,
    axis: usize,
    center: Vec3,
    radius: f64,
    tool_low: f64,
    tool_high: f64,
    tol: f64,
) -> Result<Option<Body>> {
    let bounds = validate(body, tol)?;
    let c = xyz(center);
    let mut corridor_low = c.map(|v| v - radius);
    let mut corridor_high = c.map(|v| v + radius);
    corridor_low[axis] = tool_low;
    corridor_high[axis] = tool_high;
    let mut planes: Vec<(f64, f64)> = vec![];
    let group: Vec<Option<usize>> = body
        .faces
        .iter()
        .enumerate()
        .map(|(index, face)| {
            let (lo, hi) = bounds[index];
            if (0..3).any(|k| hi[k] < corridor_low[k] - tol || lo[k] > corridor_high[k] + tol) {
                return None;
            }
            let Surface::Plane { p, n, .. } = face.surface else { return None };
            let nn = xyz(n);
            if !cap_normal(nn, axis) {
                return None;
            }
            let (height, side) = (xyz(p)[axis], nn[axis].signum() * sign(face.sense));
            Some(planes.iter().position(|&(h, s)| s == side && (h - height).abs() <= tol).unwrap_or_else(|| {
                planes.push((height, side));
                planes.len() - 1
            }))
        })
        .collect();
    let mut owner: Vec<Option<usize>> = vec![None; body.edges.len()];
    let mut seam = false;
    for (index, face) in body.faces.iter().enumerate() {
        let Some(g) = group[index] else { continue };
        for &(e, _) in face.loops.iter().flatten() {
            match owner[e] {
                Some(other) if other != index && group[other] == Some(g) => seam = true,
                _ => owner[e] = Some(index),
            }
        }
    }
    if !seam {
        return Ok(None);
    }
    let mut joined = body.clone();
    super::unite_tidy::join_faces(&mut joined, &group, &|_, _| false)?;
    let _ = super::unite_tidy::compact(&mut joined);
    super::unite_tidy::join_edges(&mut joined, tol);
    super::unite_tidy::outer_first(&mut joined);
    crate::topology::check(&joined, tol)
        .map_err(|e| KernelError::Internal(format!("joining cap seams built an invalid solid: {e}")))?;
    Ok(Some(joined))
}

/// A planar face's unit normal along the corridor axis, within 1e-12 per component: a cap. Stock built by rotated or
/// filleted profile math carries normals off the axis in the last bits (a web face's 2.6e-17); over any admitted part
/// that is far inside the kernel tolerance, so such a face is the same cap (Bob 10-28: what agrees within the kernel
/// tolerance is the same).
fn cap_normal(nn: [f64; 3], axis: usize) -> bool {
    (nn[axis].abs() - 1.).abs() <= 1e-12 && (0..3).all(|k| k == axis || nn[k].abs() <= 1e-12)
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The distance between segments p0-p1 and q0-q1 (closest points by the clamped normal equations; a degenerate
/// segment is its point).
fn segment_distance(p0: [f64; 3], p1: [f64; 3], q0: [f64; 3], q1: [f64; 3]) -> f64 {
    let (d1, d2, r) = (sub3(p1, p0), sub3(q1, q0), sub3(p0, q0));
    let (a, e, f) = (dot3(d1, d1), dot3(d2, d2), dot3(d2, r));
    let (s, t);
    if a <= f64::EPSILON && e <= f64::EPSILON {
        (s, t) = (0., 0.);
    } else if a <= f64::EPSILON {
        (s, t) = (0., (f / e).clamp(0., 1.));
    } else {
        let c = dot3(d1, r);
        if e <= f64::EPSILON {
            (s, t) = ((-c / a).clamp(0., 1.), 0.);
        } else {
            let b = dot3(d1, d2);
            let denom = a * e - b * b;
            let mut s1 = if denom > 0. { ((b * f - c * e) / denom).clamp(0., 1.) } else { 0. };
            let mut t1 = (b * s1 + f) / e;
            if t1 < 0. {
                t1 = 0.;
                s1 = (-c / a).clamp(0., 1.);
            } else if t1 > 1. {
                t1 = 1.;
                s1 = ((b - c) / a).clamp(0., 1.);
            }
            (s, t) = (s1, t1);
        }
    }
    let cp = [p0[0] + d1[0] * s, p0[1] + d1[1] * s, p0[2] + d1[2] * s];
    let cq = [q0[0] + d2[0] * t, q0[1] + d2[1] * t, q0[2] + d2[2] * t];
    let d = sub3(cp, cq);
    dot3(d, d).sqrt()
}

/// A planar face bounded by straight edges only that lies farther than radius + tol from the corridor's axis segment
/// (from tool_low to tool_high along the axis through the centre): every point of the corridor's cylinder is within
/// the radius of that segment, so the face cannot meet it. The least distance from a segment to a planar region is
/// reached on the region's boundary, where the segment crosses the region, or at a segment end above the region;
/// each is tested. Any curved boundary edge, or a point on a boundary in the inside test, answers false.
#[allow(clippy::too_many_arguments)]
fn clear_of_corridor(body: &Body, face: &Face, n: Vec3, axis: usize, c: [f64; 3], radius: f64, tool_low: f64, tool_high: f64, tol: f64) -> bool {
    let (mut a, mut b) = (c, c);
    a[axis] = tool_low;
    b[axis] = tool_high;
    let limit = radius + tol;
    let mut first = None;
    for &(e, _) in face.loops.iter().flatten() {
        let edge = &body.edges[e];
        let Curve::Line { .. } = edge.curve else { return false };
        let (p, q) = (xyz(body.vertices[edge.start]), xyz(body.vertices[edge.end]));
        first.get_or_insert(p);
        if segment_distance(a, b, p, q) <= limit {
            return false;
        }
    }
    let Some(origin) = first else { return false };
    let nn = xyz(n);
    let height = |x: [f64; 3]| dot3(sub3(x, origin), nn);
    let (ha, hb) = (height(a), height(b));
    let on_region = |x: [f64; 3]| inside(body, face, vector(x), n, tol).unwrap_or(true);
    if ha * hb <= 0. && ha != hb {
        let t = ha / (ha - hb);
        let x = [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
        if on_region(x) {
            return false;
        }
    }
    for (x, h) in [(a, ha), (b, hb)] {
        if h.abs() <= limit && on_region([x[0] - nn[0] * h, x[1] - nn[1] * h, x[2] - nn[2] * h]) {
            return false;
        }
    }
    true
}

/// Ordered (lower height, upper height, lower face, upper face) intervals.
/// A disk wholly inside each cap, alternating signs and exclusion of all other
/// patches prove complete material intervals throughout the finite corridor.
pub(super) fn round_passages(
    body: &Body,
    axis: usize,
    center: Vec3,
    radius: f64,
    tool_low: f64,
    tool_high: f64,
    tol: f64,
) -> Result<Vec<(f64, f64, usize, usize)>> {
    let bounds = validate(body, tol)?;
    let c = xyz(center);
    let mut corridor_low = c.map(|v| v - radius);
    let mut corridor_high = c.map(|v| v + radius);
    corridor_low[axis] = tool_low;
    corridor_high[axis] = tool_high;
    let mut hits = Vec::new();
    for (index, face) in body.faces.iter().enumerate() {
        let (lo, hi) = bounds[index];
        if (0..3).any(|k| hi[k] < corridor_low[k] - tol || lo[k] > corridor_high[k] + tol) {
            continue;
        }
        let Surface::Plane { p, n, .. } = face.surface else {
            return Err(unsupported("bore corridor may intersect a curved face; fillet-crossing cuts require intersection support"));
        };
        let nn = xyz(n);
        if !cap_normal(nn, axis) {
            // A side or sloped face inside the corridor's box need not meet the corridor itself: excluded when it is
            // provably farther than the radius from the corridor's axis segment.
            if clear_of_corridor(body, face, n, axis, c, radius, tool_low, tool_high, tol) {
                continue;
            }
            return Err(unsupported("bore corridor may cross a non-cap planar face"));
        }
        let height = xyz(p)[axis];
        let mut q = c;
        q[axis] = height;
        let q = vector(q);
        let clearance = face
            .loops
            .iter()
            .flatten()
            .map(|(i, _)| distance(&body.edges[*i].curve, q))
            .fold(f64::INFINITY, f64::min);
        if clearance <= radius + tol {
            return Err(unsupported(
                "bore must clear every cap boundary and existing hole",
            ));
        }
        if inside(body, face, q, n.scale(sign(face.sense)), tol)? {
            hits.push((height, nn[axis].signum() * sign(face.sense), index));
        }
    }
    hits.sort_by(|a, b| a.0.total_cmp(&b.0));
    if hits.is_empty() || hits.len() % 2 != 0 || hits.len() > 256 {
        return Err(unsupported(
            "passage requires one to 128 complete planar material intervals",
        ));
    }
    if hits.windows(2).any(|w| w[1].0 - w[0].0 <= tol) {
        return Err(unsupported("planar passage boundaries are not separated"));
    }
    let mut intervals = Vec::new();
    for pair in hits.chunks_exact(2) {
        if pair[0].1 != -1. || pair[1].1 != 1. {
            return Err(unsupported("planar passage must enter and exit each wall"));
        }
        if tool_low > pair[0].0 || tool_high < pair[1].0 {
            return Err(unsupported("passage must span the full material thickness"));
        }
        intervals.push((pair[0].0, pair[1].0, pair[0].2, pair[1].2));
    }
    Ok(intervals)
}

/// Attach the inverted side boundary of a caller-validated cavity, reusing its
/// actual cap loops as holes in the two proven entry/exit faces. The cavity is
/// produced by our own extrusion after passage admission, not guessed from bounds.
pub(super) fn splice_cavity(
    body: &Body,
    cavity: Body,
    axis: usize,
    bottom_face: usize,
    top_face: usize,
    tol: f64,
) -> Result<Body> {
    let mut result = body.clone();
    let vertex_offset = result.vertices.len();
    let edge_offset = result.edges.len();
    let mut cap_count = [0, 0];
    result.vertices.extend(cavity.vertices);
    result.edges.extend(cavity.edges.into_iter().map(|mut e| {
        e.start += vertex_offset;
        e.end += vertex_offset;
        e
    }));
    for mut face in cavity.faces {
        let loops: Vec<Vec<(usize, bool)>> = face
            .loops
            .into_iter()
            .map(|lp| {
                lp.into_iter()
                    .rev()
                    .map(|(i, fwd)| (i + edge_offset, !fwd))
                    .collect()
            })
            .collect();
        if let Surface::Plane { n, .. } = face.surface {
            let a = xyz(n);
            if a[axis].abs() == 1. && (0..3).all(|i| i == axis || a[i] == 0.) {
                let which = if a[axis] * sign(face.sense) < 0. {
                    0
                } else {
                    1
                };
                cap_count[which] += 1;
                result.faces[if which == 0 { bottom_face } else { top_face }]
                    .loops
                    .extend(loops);
                continue;
            }
        }
        face.sense = !face.sense;
        face.loops = loops;
        result.faces.push(face);
    }
    if cap_count != [1, 1] {
        return Err(KernelError::Internal(
            "passage cavity does not have exactly two axial caps".into(),
        ));
    }
    crate::topology::check(&result, tol).map_err(|e| {
        KernelError::Internal(format!("passage splice failed closed topology: {e}"))
    })?;
    validate(&result, tol).map_err(|e| {
        KernelError::Internal(format!("passage splice failed analytic validation: {e}"))
    })?;
    Ok(result)
}

#[cfg(test)]
mod codex_clearance_review {
    use super::*;
    use crate::neutral::Edge;

    #[test]
    fn long_nearly_parallel_crossing_segments_have_zero_distance() {
        let distance = segment_distance([0.,0.,0.], [0.,100000.,0.],
            [-0.0001,0.,0.], [0.0001,100000.,0.]);
        assert!(distance <= 1e-12, "segments intersect at (0,50000,0), distance={distance}");
    }

    #[test]
    fn a_planar_strip_crossing_the_bore_axis_is_not_clear() {
        let vertices = vec![Vec3::new(-0.0001,0.,0.), Vec3::new(-0.00009,0.,0.),
            Vec3::new(0.00011,100000.,0.), Vec3::new(0.0001,100000.,0.)];
        let edges = (0..4).map(|i| Edge { start:i, end:(i+1)%4,
            curve:Curve::Line { a:vertices[i], b:vertices[(i+1)%4] } }).collect();
        let face = Face { surface:Surface::Plane { p:vertices[0], n:Vec3::new(0.,0.,1.), x:Vec3::new(1.,0.,0.) },
            sense:true, loops:vec![(0..4).map(|i|(i,true)).collect()] };
        let body = Body { vertices, edges, faces:vec![face] };
        assert!(inside(&body,&body.faces[0],Vec3::new(0.,47500.,0.),Vec3::new(0.,0.,1.),1e-6).unwrap());
        assert!(!clear_of_corridor(&body,&body.faces[0],Vec3::new(0.,0.,1.),1,[0.,0.,0.],1e-5,0.,100000.,1e-6),
            "the bore axis crosses the interior of this face; it cannot be excluded");
    }
}

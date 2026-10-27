//! Local admission for a round through-passage in an otherwise unchanged solid.
//! Proves two planar entry/exit patches contain the disk with clearance, and
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
        Curve::Line { .. } => None,
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

fn curve_bounds(curve: &Curve) -> ([f64; 3], [f64; 3]) {
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

/// Returns lower height, upper height, lower face, upper face. A circular disk
/// wholly inside both planar patches plus exclusion of all other boundary
/// patches proves one constant material interval over the entire disk.
pub(super) fn round_passage(
    body: &Body,
    axis: usize,
    center: Vec3,
    radius: f64,
    tool_low: f64,
    tool_high: f64,
    tol: f64,
) -> Result<(f64, f64, usize, usize)> {
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
        if nn[axis].abs() != 1. || (0..3).any(|k| k != axis && nn[k] != 0.) {
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
            hits.push((height, nn[axis] * sign(face.sense), index));
        }
    }
    hits.sort_by(|a, b| a.0.total_cmp(&b.0));
    if hits.len() != 2 || hits[0].1 != -1. || hits[1].1 != 1. || hits[1].0 - hits[0].0 <= tol {
        return Err(unsupported(
            "bore must cross exactly one material interval bounded by two planar faces",
        ));
    }
    if tool_low > hits[0].0 || tool_high < hits[1].0 {
        return Err(unsupported("bore must span the full material thickness"));
    }
    Ok((hits[0].0, hits[1].0, hits[0].2, hits[1].2))
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

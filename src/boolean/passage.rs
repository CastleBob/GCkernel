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

fn parameters(base: f64, lo: f64, hi: f64) -> Vec<f64> {
    let first = base + ((lo - base) / TAU).ceil() * TAU;
    [first, first + TAU]
        .into_iter()
        .filter(|t| *t >= lo && *t <= hi)
        .collect()
}

fn curve_bounds(curve: &Curve) -> ([f64; 3], [f64; 3]) {
    let mut points = vec![curve.start(), curve.end()];
    if let Curve::Arc { c, n, x, r, t0, t1 } = *curve {
        let a = xyz(x);
        let b = xyz(n.cross(x));
        for k in 0..3 {
            for root in [b[k].atan2(a[k]), b[k].atan2(a[k]) + PI] {
                for t in parameters(root, t0, t1) {
                    points.push(Curve::point_at_angle(c, n, x, r, t));
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
        Curve::Conic { .. } => unreachable!("conics are refused by passage admission"),
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
            Curve::Conic { .. } => {
                return Err(unsupported("conic-trimmed passages are not admitted yet"))
            }
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
            Curve::Arc { c, n, x, r, t0, t1 } => {
                let a = r * x.dot(v);
                let b = r * n.cross(x).dot(v);
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
                        let hit = Curve::point_at_angle(c, n, x, r, t);
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
            Curve::Conic { .. } => {
                return Err(unsupported("conic-trimmed passages are not admitted yet"))
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
                    if let Curve::Arc { c, n: an, .. } = *curve {
                        if c.sub(p).dot(n).abs() > tol || an.cross(n).length() > 1e-10 {
                            return Err(bad("passage arc is not in its plane"));
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
                if !r.is_finite() || r <= tol || face.loops.len() != 1 || face.loops[0].len() != 4 {
                    return Err(unsupported(
                        "passage requires simple axial cylindrical patches",
                    ));
                }
                let mut arcs = Vec::new();
                let mut lines = 0;
                for (i, fwd) in &face.loops[0] {
                    match body.edges[*i].curve {
                        Curve::Conic { .. } => {
                            return Err(unsupported("conic-trimmed passages are not admitted yet"))
                        }
                        Curve::Line { a: from, b: to } => {
                            if to.sub(from).cross(a).length() > tol
                                || (from.sub(p).cross(a).length() - r).abs() > tol
                            {
                                return Err(bad("cylindrical generator is off its surface"));
                            }
                            lines += 1;
                        }
                        Curve::Arc {
                            c,
                            n,
                            x,
                            r: cr,
                            t0,
                            t1,
                        } => {
                            if c.sub(p).cross(a).length() > tol
                                || n.cross(a).length() > 1e-10
                                || (cr - r).abs() > tol
                            {
                                return Err(bad("cylindrical trim arc is off its surface"));
                            }
                            arcs.push((c.dot(a), n.dot(a) * sign(*fwd), c, n, x, cr, t0, t1));
                        }
                    }
                }
                if arcs.len() != 2 || lines != 2 {
                    return Err(unsupported(
                        "cylinder patch must have two circular trims and two generators",
                    ));
                }
                arcs.sort_by(|a, b| a.0.total_cmp(&b.0));
                if arcs[1].0 - arcs[0].0 <= tol
                    || arcs[0].1 * sign(face.sense) < 0.999999
                    || arcs[1].1 * sign(face.sense) > -0.999999
                {
                    return Err(bad(
                        "cylinder patch axial ordering or orientation is invalid",
                    ));
                }
                let first = &arcs[0];
                let second = &arcs[1];
                let m1 = Curve::point_at_angle(
                    first.2,
                    first.3,
                    first.4,
                    first.5,
                    (first.6 + first.7) * 0.5,
                );
                let m2 = Curve::point_at_angle(
                    second.2,
                    second.3,
                    second.4,
                    second.5,
                    (second.6 + second.7) * 0.5,
                );
                if ((first.7 - first.6) - (second.7 - second.6)).abs() > 1e-10
                    || m2.sub(m1).cross(a).length() > tol
                {
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

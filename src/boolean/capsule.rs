//! Contained capsule-shaped through-slot in a plain box. Admission verifies a
//! complete analytic capsule prism (including oriented topology), then preserves
//! its original cap curves as the hole boundary. No boundary-crossing solver,
//! approximate polygon, coordinate offset or failed-operation retry is used.
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Face, Profile};
use crate::KernelError;
use std::f64::consts::PI;
type Result<T> = std::result::Result<T, KernelError>;
fn invalid(s: &str) -> KernelError {
    KernelError::InvalidInput(s.into())
}
fn unsupported(s: &str) -> KernelError {
    KernelError::Unsupported(s.into())
}
fn xyz(p: Vec3) -> [f64; 3] {
    [p.x, p.y, p.z]
}
fn vec(p: [f64; 3]) -> Vec3 {
    Vec3::new(p[0], p[1], p[2])
}
fn unit(v: Vec3) -> bool {
    v.is_finite() && (v.length() - 1.).abs() < 1e-10
}
fn cardinal(v: Vec3) -> Result<usize> {
    if !unit(v) {
        return Err(invalid("slot axis is not a finite unit direction"));
    }
    let p = xyz(v);
    let axes: Vec<usize> = (0..3).filter(|i| p[*i] != 0.).collect();
    if axes.len() != 1 {
        return Err(unsupported("slot sweep must follow a coordinate axis"));
    }
    Ok(axes[0])
}

fn near(a: Vec3, b: Vec3, tol: f64) -> bool {
    a.distance(b) <= tol
}
fn middle(c: &Curve) -> Vec3 {
    let (lo, hi) = c.domain();
    c.eval((lo + hi) * 0.5).0
}
fn same_curve(a: &Curve, b: &Curve, tol: f64) -> bool {
    if !near(a.start(), b.start(), tol)
        || !near(a.end(), b.end(), tol)
        || !near(middle(a), middle(b), tol)
    {
        return false;
    }
    match (a, b) {
        (Curve::Line { .. }, Curve::Line { .. }) => true,
        (
            Curve::Arc {
                c: ac,
                n: an,
                r: ar,
                t0: a0,
                t1: a1,
                ..
            },
            Curve::Arc {
                c: bc,
                n: bn,
                r: br,
                t0: b0,
                t1: b1,
                ..
            },
        ) => {
            near(*ac, *bc, tol)
                && (ar - br).abs() <= tol
                && an.dot(*bn) > 1. - 1e-10
                && ((a1 - a0) - (b1 - b0)).abs() < 1e-10
        }
        _ => false,
    }
}
fn same_surface(a: &Face, b: &Face, tol: f64) -> bool {
    match (&a.surface, &b.surface) {
        (Surface::Plane { p: ap, n: an, .. }, Surface::Plane { p: bp, n: bn, .. }) => {
            let ao = an.scale(if a.sense { 1. } else { -1. });
            let bo = bn.scale(if b.sense { 1. } else { -1. });
            ao.dot(bo) > 1. - 1e-10 && ap.sub(*bp).dot(*bn).abs() <= tol
        }
        (
            Surface::Cylinder {
                p: ap,
                a: aa,
                r: ar,
                ..
            },
            Surface::Cylinder {
                p: bp,
                a: ba,
                r: br,
                ..
            },
        ) => {
            a.sense == b.sense
                && aa.dot(*ba).abs() > 1. - 1e-10
                && (ar - br).abs() <= tol
                && ap.sub(*bp).cross(*ba).length() <= tol
        }
        _ => false,
    }
}
fn cyclic_equal(a: &[(usize, bool)], b: &[(usize, bool)]) -> bool {
    a.len() == b.len()
        && (0..b.len()).any(|offset| (0..a.len()).all(|i| a[i] == b[(i + offset) % b.len()]))
}

/// Match complete finite B-reps, up to indexing, curve direction, loop start,
/// and equivalent surface frame conventions. This is admission verification,
/// not a conversion: output later uses the original input cap curves.
fn match_prism(actual: &Body, expected: &Body, tol: f64) -> Result<Vec<usize>> {
    let mut vertices = Vec::new();
    let mut seen = vec![false; expected.vertices.len()];
    for p in &actual.vertices {
        let candidates: Vec<usize> = expected
            .vertices
            .iter()
            .enumerate()
            .filter_map(|(i, q)| near(*p, *q, tol).then_some(i))
            .collect();
        if candidates.len() != 1 || seen[candidates[0]] {
            return Err(unsupported(
                "slot vertices do not identify one capsule prism within tolerance",
            ));
        }
        seen[candidates[0]] = true;
        vertices.push(candidates[0]);
    }
    let mut edge_map = Vec::new();
    let mut seen = vec![false; expected.edges.len()];
    for edge in &actual.edges {
        let start = vertices[edge.start];
        let end = vertices[edge.end];
        let mut candidates = Vec::new();
        for (i, e) in expected.edges.iter().enumerate() {
            if start == e.start && end == e.end && same_curve(&edge.curve, &e.curve, tol) {
                candidates.push((i, false));
            }
            if start == e.end && end == e.start && same_curve(&edge.curve, &e.curve.reversed(), tol)
            {
                candidates.push((i, true));
            }
        }
        if candidates.len() != 1 || seen[candidates[0].0] {
            return Err(unsupported(
                "slot edges do not form the supported capsule prism",
            ));
        }
        seen[candidates[0].0] = true;
        edge_map.push(candidates[0]);
    }
    let mut face_map = Vec::new();
    let mut seen = vec![false; expected.faces.len()];
    for face in &actual.faces {
        let boundary: Vec<(usize, bool)> = face.loops[0]
            .iter()
            .map(|(edge, forward)| {
                let (mapped, reverse) = edge_map[*edge];
                (mapped, *forward ^ reverse)
            })
            .collect();
        let candidates: Vec<usize> = expected
            .faces
            .iter()
            .enumerate()
            .filter_map(|(i, f)| {
                (same_surface(face, f, tol) && cyclic_equal(&boundary, &f.loops[0])).then_some(i)
            })
            .collect();
        if candidates.len() != 1 || seen[candidates[0]] {
            return Err(unsupported(
                "slot face geometry or orientation is not the supported prism",
            ));
        }
        seen[candidates[0]] = true;
        face_map.push(candidates[0]);
    }
    Ok(face_map)
}

fn capsule_profile(a: Vec3, b: Vec3, n: Vec3, r: f64) -> Result<Vec<Curve>> {
    let d = b.sub(a).unit()?;
    let side = n.cross(d);
    let am = a.sub(side.scale(r));
    let ap = a.add(side.scale(r));
    let bm = b.sub(side.scale(r));
    let bp = b.add(side.scale(r));
    Ok(vec![
        Curve::Line { a: am, b: bm },
        Curve::Arc {
            c: b,
            n,
            x: side.scale(-1.),
            r,
            t0: 0.,
            t1: PI,
        },
        Curve::Line { a: bp, b: ap },
        Curve::Arc {
            c: a,
            n,
            x: side,
            r,
            t0: 0.,
            t1: PI,
        },
    ])
}

struct Admitted {
    axis: usize,
    n: Vec3,
    low: f64,
    high: f64,
    lower_face: usize,
    enclosing_center: Vec3,
    enclosing_radius: f64,
}
fn admit(tool: &Body, tol: f64) -> Result<Admitted> {
    if tool.vertices.len() != 8
        || tool.edges.len() != 12
        || tool.faces.len() != 6
        || tool
            .faces
            .iter()
            .any(|f| f.loops.len() != 1 || f.loops[0].len() != 4)
    {
        return Err(unsupported("slot cutter requires two semicircles and two straight sides swept into a six-face prism"));
    }
    if tool.vertices.iter().any(|v| !v.is_finite())
        || tool.edges.iter().any(|e| e.start >= 8 || e.end >= 8)
        || tool
            .faces
            .iter()
            .flat_map(|f| &f.loops)
            .flatten()
            .any(|(i, _)| *i >= 12)
    {
        return Err(invalid("invalid slot topology references or vertices"));
    }
    for edge in &tool.edges {
        if !edge.curve.start().is_finite()
            || !edge.curve.end().is_finite()
            || !middle(&edge.curve).is_finite()
        {
            return Err(invalid("nonfinite slot curve"));
        }
        if let Curve::Arc { c, n, x, r, t0, t1 } = edge.curve {
            if !c.is_finite()
                || !unit(n)
                || !unit(x)
                || n.dot(x).abs() > 1e-10
                || !r.is_finite()
                || r <= 0.
                || !t0.is_finite()
                || !t1.is_finite()
                || t1 <= t0
            {
                return Err(invalid("invalid slot arc frame"));
            }
        }
    }
    let mut cylinders = Vec::new();
    for face in &tool.faces {
        let (p, n, x, r) = match face.surface {
            Surface::Plane { p, n, x } => (p, n, x, None),
            Surface::Cylinder { p, a, x, r } => (p, a, x, Some(r)),
        };
        if !p.is_finite()
            || !unit(n)
            || !unit(x)
            || n.dot(x).abs() > 1e-10
            || r.is_some_and(|r| !r.is_finite() || r <= 0.)
        {
            return Err(invalid("invalid slot surface frame"));
        }
        if let Some(r) = r {
            cylinders.push((p, n, r));
        }
    }
    if cylinders.len() != 2 {
        return Err(unsupported("slot requires two semicylindrical walls"));
    }
    let axis = cardinal(cylinders[0].1)?;
    if cardinal(cylinders[1].1)? != axis || cylinders[0].2 != cylinders[1].2 {
        return Err(unsupported("slot end radii and axes must agree"));
    }
    let radius = cylinders[0].2;
    if radius <= tol {
        return Err(unsupported("slot radius is below supported resolution"));
    }
    let low = tool
        .vertices
        .iter()
        .map(|v| xyz(*v)[axis])
        .fold(f64::INFINITY, f64::min);
    let high = tool
        .vertices
        .iter()
        .map(|v| xyz(*v)[axis])
        .fold(f64::NEG_INFINITY, f64::max);
    if !high.is_finite() || !low.is_finite() || high - low <= tol {
        return Err(invalid("slot has no supported axial extent"));
    }
    let mut a = xyz(cylinders[0].0);
    let mut b = xyz(cylinders[1].0);
    a[axis] = low;
    b[axis] = low;
    let a = vec(a);
    let b = vec(b);
    if a.distance(b) <= tol {
        return Err(unsupported("slot end centers must be distinct"));
    }
    let mut n = [0.; 3];
    n[axis] = 1.;
    let n = vec(n);
    let expected = crate::extrude::extrude(
        &Profile {
            loops: vec![capsule_profile(a, b, n, radius)?],
            n,
            low: 0.,
            high: high - low,
        },
        tol,
    )?;
    crate::topology::check(tool, tol)?;
    let mapping = match_prism(tool, &expected, tol)?;
    let expected_lower = expected
        .faces
        .iter()
        .position(|f| {
            matches!(f.surface,Surface::Plane {n:fnorm,..}
        if fnorm.dot(n)*(if f.sense {1.}else{-1.})< -0.999999)
        })
        .unwrap();
    let lower_face = mapping.iter().position(|i| *i == expected_lower).unwrap();
    Ok(Admitted {
        axis,
        n,
        low,
        high,
        lower_face,
        enclosing_center: a.add(b).scale(0.5),
        enclosing_radius: radius + a.distance(b) * 0.5,
    })
}

fn boundary(body: &Body, face: &Face) -> Vec<Curve> {
    face.loops[0]
        .iter()
        .map(|(i, forward)| {
            let c = &body.edges[*i].curve;
            if *forward {
                c.clone()
            } else {
                c.reversed()
            }
        })
        .collect()
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    if !tol.is_finite() || tol <= 0. {
        return Err(invalid("slot tolerance must be finite and positive"));
    }
    let admitted = admit(tool, tol)?;
    let axis = admitted.axis;
    let n = admitted.n;
    if super::endcut::crosses_cylinder(
        body,
        axis,
        admitted.enclosing_center,
        admitted.enclosing_radius,
        admitted.low,
        admitted.high,
        tol,
    ) {
        return super::endcut::cross_cut(body, tool, tol);
    }
    if !(body.vertices.len() == 8 && body.edges.len() == 12 && body.faces.len() == 6) {
        let (low, high, bottom_face, top_face) = super::passage::round_passage(
            body,
            axis,
            admitted.enclosing_center,
            admitted.enclosing_radius,
            admitted.low,
            admitted.high,
            tol,
        )?;
        let shift = n.scale(low - admitted.low);
        let outline = boundary(tool, &tool.faces[admitted.lower_face])
            .into_iter()
            .map(|c| c.translated(shift))
            .collect();
        let cavity = crate::extrude::extrude(
            &Profile {
                loops: vec![outline],
                n,
                low: 0.,
                high: high - low,
            },
            tol,
        )?;
        return super::passage::splice_cavity(body, cavity, axis, bottom_face, top_face, tol);
    }
    let (low, high) = super::rectilinear::plain_box_bounds(body, tol)?;
    if admitted.low > low[axis] || admitted.high < high[axis] {
        return Err(unsupported(
            "slot cutter must span the entire box thickness",
        ));
    }
    let hole = boundary(tool, &tool.faces[admitted.lower_face]);
    // Conservative analytic bounds of the ORIGINAL curves, not just their
    // vertices. A semicircle can leave the box between its endpoints.
    for curve in &hole {
        let (minimum, maximum) = match curve {
            Curve::Conic { .. } | Curve::CylinderPair { .. } => {
                return Err(unsupported("capsule cutters require circular end arcs"))
            }
            Curve::Line { a, b } => {
                let a = xyz(*a);
                let b = xyz(*b);
                (
                    std::array::from_fn(|i| a[i].min(b[i])),
                    std::array::from_fn(|i| a[i].max(b[i])),
                )
            }
            Curve::Arc { c, r, .. } => {
                let c = xyz(*c);
                (c.map(|v| v - r), c.map(|v| v + r))
            }
        };
        if (0..3)
            .filter(|i| *i != axis)
            .any(|i| minimum[i] - low[i] <= tol || high[i] - maximum[i] <= tol)
        {
            return Err(unsupported(
                "slot must remain strictly inside the box side walls",
            ));
        }
    }
    let lower = body
        .faces
        .iter()
        .find(|f| {
            matches!(f.surface,Surface::Plane {p,n:fnorm,..}
        if xyz(p)[axis]==low[axis]&&fnorm.dot(n)*(if f.sense {1.}else{-1.})< -0.999999)
        })
        .ok_or_else(|| invalid("box lower cap is missing"))?;
    let outer: Vec<Curve> = boundary(body, lower)
        .iter()
        .rev()
        .map(Curve::reversed)
        .collect();
    let shift = n.scale(low[axis] - admitted.low);
    let hole = hole.iter().map(|c| c.translated(shift)).collect();
    let result = crate::extrude::extrude(
        &Profile {
            loops: vec![outer, hole],
            n,
            low: 0.,
            high: high[axis] - low[axis],
        },
        tol,
    )?;
    crate::topology::check(&result, tol).map_err(|e| {
        KernelError::Internal(format!("slot result failed topology validation: {e}"))
    })?;
    Ok(result)
}

#[cfg(test)]
#[path = "capsule_tests.rs"]
mod tests;

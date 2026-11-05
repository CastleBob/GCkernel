//! Transverse end cuts on analytic solids with rectangular material sections.
//! In an orthonormal frame, each admitted X band is
//! Ylo(X)..Yhi(X) times Zlo(X)..Zhi(X). Boundaries are planes/cylinders,
//! reconstructed from the actual current B-rep, including prior cuts.
//! This is a bounded geometric domain, not an I-beam recipe or a fallback.
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;
use std::f64::consts::{PI, TAU};

#[path = "endcut_cross.rs"]
mod cross;
#[path = "endcut_split.rs"]
mod split;

type Result<T> = std::result::Result<T, KernelError>;
fn unsupported(s: &str) -> KernelError {
    KernelError::Unsupported(s.into())
}
fn invalid(s: &str) -> KernelError {
    KernelError::InvalidInput(s.into())
}
fn p(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn xyz(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn sign(b: bool) -> f64 {
    if b {
        1.
    } else {
        -1.
    }
}
fn eval(c: &Curve, t: f64) -> Vec3 {
    c.eval(t).0
}
fn domain(c: &Curve) -> (f64, f64) {
    c.domain()
}
fn conic(c: &Curve) -> Option<(Vec3, Vec3, Vec3, f64, f64)> {
    match *c {
        Curve::Line { .. } | Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } | Curve::Bezier { .. } | Curve::RationalBezier { .. } => None,
        Curve::Arc { c, n, x, r, t0, t1 } => Some((c, x.scale(r), n.cross(x).scale(r), t0, t1)),
        Curve::Conic { c, u, v, t0, t1 } => Some((c, u, v, t0, t1)),
    }
}
// Only the X coordinate is used for partition intersections. For admitted
// cylinder pairs a.x=0, so the axial root cannot change X.
fn section_projection(curve: &Curve) -> Option<(Vec3, Vec3, Vec3, f64, f64)> {
    if let Curve::CylinderPair {
        c,
        a,
        x,
        rho,
        t0,
        t1,
        ..
    } = *curve
    {
        Some((c, x.scale(rho), a.cross(x).scale(rho), t0, t1))
    } else {
        conic(curve)
    }
}
fn parameters(t: f64, lo: f64, hi: f64) -> Vec<f64> {
    let first = t + ((lo - t) / TAU).floor() * TAU;
    [first, first + TAU, first + 2. * TAU]
        .into_iter()
        .filter(|s| *s >= lo - 1e-13 && *s <= hi + 1e-13)
        .map(|s| s.clamp(lo, hi))
        .collect()
}
fn roots(c: &Curve, x: f64, tol: f64) -> Vec<Vec3> {
    match *c {
        Curve::Line { a, b } => {
            if a.x == b.x || x < a.x.min(b.x) - tol || x > a.x.max(b.x) + tol {
                return vec![];
            }
            vec![eval(c, ((x - a.x) / (b.x - a.x)).clamp(0., 1.))]
        }
        _ => {
            let (center, u, v, t0, t1) = section_projection(c).unwrap();
            let amplitude = u.x.hypot(v.x);
            if amplitude <= 1e-14 {
                return vec![];
            }
            let ratio = (x - center.x) / amplitude;
            if ratio.abs() > 1. + tol / amplitude {
                return vec![];
            }
            let angle = ratio.clamp(-1., 1.).acos();
            let phase = v.x.atan2(u.x);
            let mut points = Vec::new();
            for theta in [phase - angle, phase + angle] {
                for t in parameters(theta, t0, t1) {
                    let hit = eval(c, t);
                    if points.iter().all(|q: &Vec3| q.distance(hit) > tol) {
                        points.push(hit);
                    }
                }
            }
            points
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Graph {
    Line {
        m: f64,
        b: f64,
    },
    Circle {
        cx: f64,
        cv: f64,
        r: f64,
        branch: f64,
    },
}
impl Graph {
    fn value(self, x: f64) -> f64 {
        match self {
            Self::Line { m, b } => m * x + b,
            Self::Circle { cx, cv, r, branch } => {
                // Preserve the analytic extremal station exactly: subtracting
                // cx again can turn (cx+r)-cx into r-1ulp, inventing a root.
                if x == cx - r || x == cx + r {
                    cv
                } else {
                    cv + branch * (r * r - (x - cx) * (x - cx)).max(0.).sqrt()
                }
            }
        }
    }
    #[cfg(test)]
    fn range(self, a: f64, b: f64) -> (f64, f64) {
        let va = self.value(a);
        let vb = self.value(b);
        let mut lo = va.min(vb);
        let mut hi = va.max(vb);
        if let Self::Circle { cx, .. } = self {
            if cx > a && cx < b {
                let v = self.value(cx);
                lo = lo.min(v);
                hi = hi.max(v);
            }
        }
        (lo, hi)
    }
}
#[derive(Clone, Debug)]
struct Band {
    lo: f64,
    hi: f64,
    y0: Graph,
    y1: Graph,
    z0: Graph,
    z1: Graph,
}

/// Partition stations identify coincident topology within the supplied input
/// tolerance. Supporting curves/surfaces are never moved or approximated.
fn stations(mut values: Vec<f64>, tol: f64) -> Result<Vec<f64>> {
    if values.iter().any(|v| !v.is_finite()) {
        return Err(invalid("nonfinite end-cut partition"));
    }
    values.sort_by(f64::total_cmp);
    let mut result: Vec<f64> = Vec::new();
    for value in values {
        if result.last().is_none_or(|last| (value - *last) > tol) {
            result.push(value);
        }
    }
    if result.len() > 256 {
        return Err(unsupported(
            "end-cut partition exceeds 256 analytic stations",
        ));
    }
    Ok(result)
}
fn events(body: &Body) -> Vec<f64> {
    let mut values: Vec<f64> = body.vertices.iter().map(|p| p.x).collect();
    for edge in &body.edges {
        if let Some((_, u, v, t0, t1)) = section_projection(&edge.curve) {
            let phase = v.x.atan2(u.x);
            for angle in [phase, phase + PI] {
                for t in parameters(angle, t0, t1) {
                    values.push(eval(&edge.curve, t).x);
                }
            }
        }
    }
    values
}

fn validate(body: &Body, tol: f64) -> Result<()> {
    if body.faces.len() > 5000 || body.edges.len() > 20000 || body.vertices.len() > 20000 {
        return Err(unsupported("end-cut input exceeds its topology work limit"));
    }
    if body.vertices.iter().any(|v| !v.is_finite())
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
        return Err(invalid("invalid end-cut topology references"));
    }
    for e in &body.edges {
        if !e.curve.start().is_finite() || !e.curve.end().is_finite() {
            return Err(invalid("nonfinite end-cut curve"));
        }
        if let Curve::CylinderPair {
            c,
            a,
            x,
            rho,
            p: host,
            b,
            r,
            branch,
            t0,
            t1,
        } = e.curve
        {
            let unit = |v: Vec3| v.is_finite() && (v.length() - 1.).abs() < 1e-10;
            if !c.is_finite()
                || !host.is_finite()
                || !unit(a)
                || !unit(x)
                || !unit(b)
                || a.dot(x).abs() > 1e-10
                || !rho.is_finite()
                || rho <= tol
                || !r.is_finite()
                || r <= tol
                || branch.abs() != 1.
                || !t0.is_finite()
                || !t1.is_finite()
                || t1 <= t0
                || t1 - t0 > TAU + 1e-12
                || t0.abs() > 1e5
                || t1.abs() > 1e5
            {
                return Err(invalid("invalid cylinder-pair trim frame"));
            }
            let yz = a.x == 0.
                && b.x == 0.
                && ((a.y == 0. && a.z.abs() == 1. && b.z == 0. && b.y.abs() == 1.)
                    || (a.z == 0. && a.y.abs() == 1. && b.y == 0. && b.z.abs() == 1.));
            if !yz {
                return Err(unsupported(
                    "transverse cylinder-pair trim axes must match the section frame",
                ));
            }
            // Bob 2026-10-28: a trim may end where its cylinders touch within tolerance.
            if crate::geom::pair_trim_touches(c, a, x, rho, host, b, r, t0, t1, tol).is_none() {
                return Err(unsupported(
                    "cylinder-pair trim approaches a tangency or branch singularity",
                ));
            }
        }
        if let Some((c, u, v, t0, t1)) = conic(&e.curve) {
            if !c.is_finite()
                || !u.is_finite()
                || !v.is_finite()
                || u.cross(v).length() < 1e-20
                || !t0.is_finite()
                || !t1.is_finite()
                || t1 <= t0
                || t1 - t0 > TAU + 1e-10
                || t0.abs() > 1e5
                || t1.abs() > 1e5
            {
                return Err(invalid("invalid end-cut conic frame"));
            }
        }
    }
    crate::topology::check(body, tol)?;
    for face in &body.faces {
        let (anchor, n, x) = match face.surface {
            Surface::Plane { p, n, x } => (p, n, x),
            Surface::Cylinder { p, a, x, .. } => (p, a, x),
            Surface::EllipticCylinder { .. } => return Err(unsupported("end cuts across elliptic-cylinder faces are not supported")),
            Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => return Err(unsupported("end cuts across revolved or spline faces are not supported")),
        };
        if !anchor.is_finite()
            || !n.is_finite()
            || !x.is_finite()
            || (n.length() - 1.).abs() > 1e-10
            || (x.length() - 1.).abs() > 1e-10
            || n.dot(x).abs() > 1e-10
        {
            return Err(invalid("invalid end-cut surface frame"));
        }
        for (i, _) in face.loops.iter().flatten() {
            let c = &body.edges[*i].curve;
            if let Curve::CylinderPair {
                c: cp,
                a: ca,
                rho,
                p: hp,
                b: ha,
                r: hr,
                ..
            } = *c
            {
                let Surface::Cylinder { p, a, r, .. } = face.surface else {
                    return Err(unsupported(
                        "cylinder-pair trim must belong to its two supporting cylinders",
                    ));
                };
                let same = |q: Vec3, d: Vec3, radius: f64| {
                    a.cross(d).length() < 1e-10
                        && p.sub(q).cross(a).length() <= tol
                        && (r - radius).abs() <= tol
                };
                if !same(cp, ca, rho) && !same(hp, ha, hr) {
                    return Err(invalid("cylinder-pair trim is off its supporting surface"));
                }
                continue;
            }
            match face.surface {
                Surface::Plane { p, n, .. } => {
                    let miss = if let Some((center, u, v, _, _)) = conic(c) {
                        center.sub(p).dot(n).abs() + u.dot(n).abs() + v.dot(n).abs()
                    } else {
                        c.start()
                            .sub(p)
                            .dot(n)
                            .abs()
                            .max(c.end().sub(p).dot(n).abs())
                    };
                    if miss > tol {
                        return Err(invalid("end-cut curve is off its plane"));
                    }
                }
                Surface::Cylinder { p, a, r, .. } => {
                    if !r.is_finite() || r <= tol {
                        return Err(invalid("invalid end-cut cylinder radius"));
                    }
                    let perpendicular = |v: Vec3| v.sub(a.scale(v.dot(a)));
                    if let Some((center, u, v, _, _)) = conic(c) {
                        let q = perpendicular(center.sub(p));
                        let u = perpendicular(u);
                        let v = perpendicular(v);
                        let bound = (q.dot(q) + (u.dot(u) + v.dot(v)) * 0.5 - r * r).abs()
                            + 2. * q.dot(u).abs()
                            + 2. * q.dot(v).abs()
                            + ((u.dot(u) - v.dot(v)) * 0.5).abs()
                            + u.dot(v).abs();
                        if bound > tol * (2. * r + tol) {
                            return Err(invalid("end-cut conic is off its cylinder"));
                        }
                    } else if perpendicular(c.end().sub(c.start())).length() > tol
                        || (perpendicular(c.start().sub(p)).length() - r).abs() > tol
                    {
                        return Err(invalid("end-cut generator is off its cylinder"));
                    }
                }
                Surface::EllipticCylinder { .. } => {
                    return Err(unsupported("end cuts across elliptic-cylinder faces are not supported"));
                }
                Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => return Err(unsupported("end cuts across revolved or spline faces are not supported")),
            }
        }
    }
    Ok(())
}

struct Hit {
    graph: Graph,
    role: usize,
    lo: f64,
    hi: f64,
}
fn hit(face: &Face, points: &[Vec3], tol: f64) -> Result<Hit> {
    let center = points[0].add(points[1]).scale(0.5);
    let (graph, is_z, out) = match face.surface {
        Surface::Plane { p, n, .. } if n.y == 0. && n.z != 0. => (
            Graph::Line {
                m: -n.x / n.z,
                b: p.z + n.x * p.x / n.z,
            },
            true,
            n.z * sign(face.sense),
        ),
        Surface::Plane { p, n, .. } if n.z == 0. && n.y != 0. => (
            Graph::Line {
                m: -n.x / n.y,
                b: p.y + n.x * p.x / n.y,
            },
            false,
            n.y * sign(face.sense),
        ),
        Surface::Cylinder { p, a, r, .. } if a.y.abs() == 1. && a.x == 0. && a.z == 0. => {
            let branch = if center.z > p.z { 1. } else { -1. };
            (
                Graph::Circle {
                    cx: p.x,
                    cv: p.z,
                    r,
                    branch,
                },
                true,
                branch * sign(face.sense),
            )
        }
        Surface::Cylinder { p, a, r, .. } if a.z.abs() == 1. && a.x == 0. && a.y == 0. => {
            let branch = if center.y > p.y { 1. } else { -1. };
            (
                Graph::Circle {
                    cx: p.x,
                    cv: p.y,
                    r,
                    branch,
                },
                false,
                branch * sign(face.sense),
            )
        }
        _ => {
            return Err(unsupported(
                "end-cut surface does not bound a rectangular transverse section",
            ))
        }
    };
    let coord = |p: Vec3| if is_z { p.z } else { p.y };
    let span = |p: Vec3| if is_z { p.y } else { p.z };
    if points
        .iter()
        .any(|p| (coord(*p) - graph.value(p.x)).abs() > tol)
        || (span(points[0]) - span(points[1])).abs() <= tol
    {
        return Err(unsupported(
            "end-cut section does not have one straight segment per surface",
        ));
    }
    let role = if is_z {
        if out > 0. {
            3
        } else {
            2
        }
    } else if out > 0. {
        1
    } else {
        0
    };
    Ok(Hit {
        graph,
        role,
        lo: span(points[0]).min(span(points[1])),
        hi: span(points[0]).max(span(points[1])),
    })
}

// Independent single-rectangle admission reference for the original rat-hole tests.
#[cfg(test)]
fn section(body: &Body, lo: f64, hi: f64, tol: f64) -> Result<Option<Band>> {
    let mid = lo + (hi - lo) * 0.5;
    let mut hits: [Option<Hit>; 4] = [None, None, None, None];
    let mut count = 0;
    for face in &body.faces {
        if matches!(face.surface,Surface::Plane {n,..} if n.y==0.&&n.z==0.) {
            continue;
        }
        let mut points = Vec::new();
        for (i, _) in face.loops.iter().flatten() {
            for q in roots(&body.edges[*i].curve, mid, tol) {
                if points.iter().all(|p: &Vec3| p.distance(q) > tol) {
                    points.push(q);
                }
            }
        }
        if points.is_empty() {
            continue;
        }
        if points.len() != 2 {
            return Err(unsupported(
                "end-cut needs one rectangular material component in every transverse band",
            ));
        }
        let h = hit(face, &points, tol)?;
        let role = h.role;
        if hits[role].is_some() {
            return Err(unsupported(
                "end-cut section has multiple material intervals or split surface patches",
            ));
        }
        hits[role] = Some(h);
        count += 1;
    }
    if count == 0 {
        return Ok(None);
    }
    if count != 4 {
        return Err(invalid("end-cut section does not close as a rectangle"));
    }
    let [Some(y0), Some(y1), Some(z0), Some(z1)] = hits else {
        unreachable!()
    };
    let ys = (y0.graph.value(mid), y1.graph.value(mid));
    let zs = (z0.graph.value(mid), z1.graph.value(mid));
    if ys.1 - ys.0 <= tol
        || zs.1 - zs.0 <= tol
        || [&y0, &y1]
            .iter()
            .any(|h| (h.lo - zs.0).abs() > tol || (h.hi - zs.1).abs() > tol)
        || [&z0, &z1]
            .iter()
            .any(|h| (h.lo - ys.0).abs() > tol || (h.hi - ys.1).abs() > tol)
    {
        return Err(invalid(
            "end-cut surface sections do not meet at rectangle corners",
        ));
    }
    let b = Band {
        lo,
        hi,
        y0: y0.graph,
        y1: y1.graph,
        z0: z0.graph,
        z1: z1.graph,
    };
    if b.y1.range(lo, hi).0 - b.y0.range(lo, hi).1 <= tol
        || b.z1.range(lo, hi).0 - b.z0.range(lo, hi).1 <= tol
    {
        return Err(unsupported(
            "end-cut cannot prove positive rectangular thickness throughout a band",
        ));
    }
    Ok(Some(b))
}

// Candidate extrema of the difference of two admitted line/circle graphs.
// Circle derivatives are -branch*(x-cx)/sqrt(r^2-(x-cx)^2). Equating
// two and squaring leaves two linear equations; extra candidates are harmless.
fn difference_stations(a: Graph, b: Graph, lo: f64, hi: f64) -> Vec<f64> {
    let mut xs = vec![lo, hi];
    match (a, b) {
        (Graph::Circle { cx: ac, r: ar, .. }, Graph::Circle { cx: bc, r: br, .. }) => {
            for sign in [-1., 1.] {
                let divisor = br + sign * ar;
                if divisor != 0. {
                    xs.push((br * ac + sign * ar * bc) / divisor);
                }
            }
        }
        (Graph::Circle { cx, r, branch, .. }, Graph::Line { m, .. })
        | (Graph::Line { m, .. }, Graph::Circle { cx, r, branch, .. }) => {
            xs.push(cx - branch * r * m / m.hypot(1.));
        }
        _ => {}
    }
    xs.retain(|x| x.is_finite() && *x >= lo && *x <= hi);
    xs
}
#[cfg(test)]
fn below(a: Graph, b: Graph, lo: f64, hi: f64, tol: f64) -> bool {
    difference_stations(a, b, lo, hi)
        .into_iter()
        .all(|x| a.value(x) - b.value(x) <= tol)
}
#[cfg(test)]
fn subtract(mut body: Band, tool: Option<Band>, tol: f64) -> Result<Option<Band>> {
    let Some(tool) = tool else {
        return Ok(Some(body));
    };
    let (lo, hi) = (body.lo, body.hi);
    // Disjointness is proved over the WHOLE interval, never by one sample.
    if below(tool.y1, body.y0, lo, hi, 0.)
        || below(body.y1, tool.y0, lo, hi, 0.)
        || below(tool.z1, body.z0, lo, hi, 0.)
        || below(body.z1, tool.z0, lo, hi, 0.)
    {
        return Ok(Some(body));
    }
    if !below(tool.z0, body.z0, lo, hi, 0.) || !below(body.z1, tool.z1, lo, hi, 0.) {
        return Err(unsupported(
            "transverse end-cut tool must span the entire material depth",
        ));
    }
    if below(tool.y0, body.y0, lo, hi, 0.) && below(body.y1, tool.y1, lo, hi, 0.) {
        return Ok(None);
    }
    if below(tool.y0, body.y0, lo, hi, 0.)
        && below(body.y0, tool.y1, lo, hi, -tol)
        && below(tool.y1, body.y1, lo, hi, -tol)
    {
        body.y0 = tool.y1;
    } else if below(body.y1, tool.y1, lo, hi, 0.)
        && below(body.y0, tool.y0, lo, hi, -tol)
        && below(tool.y0, body.y1, lo, hi, -tol)
    {
        body.y1 = tool.y0;
    } else {
        return Err(unsupported(
            "cut does not remove one provable end interval throughout its band",
        ));
    }
    for y in [body.y0, body.y1] {
        for z in [body.z0, body.z1] {
            if matches!((y, z), (Graph::Circle { .. }, Graph::Circle { .. })) {
                return Err(unsupported("cylinder/cylinder trim curves are not implemented in the transverse end-cut domain"));
            }
        }
    }
    Ok(Some(body))
}

struct Patch {
    surface: Surface,
    sense: bool,
    boundary: Vec<Curve>,
}
fn segment(a: Vec3, b: Vec3) -> Curve {
    Curve::Line { a, b }
}
fn reversed(boundary: Vec<Curve>) -> Vec<Curve> {
    boundary.into_iter().rev().map(|c| c.reversed()).collect()
}

fn circle_angle(cx: f64, r: f64, x: f64) -> f64 {
    if x == cx + r {
        0.
    } else if x == cx - r {
        PI
    } else {
        ((x - cx) / r).clamp(-1., 1.).acos()
    }
}

fn lifted(lo: f64, hi: f64, y: Graph, z: Graph, tol: f64) -> Result<Curve> {
    let (center, u, v, cx, r, slope) = match (y, z) {
        (Graph::Line { .. }, Graph::Line { .. }) => {
            return Ok(segment(
                p(lo, y.value(lo), z.value(lo)),
                p(hi, y.value(hi), z.value(hi)),
            ))
        }
        (Graph::Circle { cx, cv, r, branch }, Graph::Line { m, b }) => (
            p(cx, cv, m * cx + b),
            p(r, 0., m * r),
            p(0., branch * r, 0.),
            cx,
            r,
            m,
        ),
        (Graph::Line { m, b }, Graph::Circle { cx, cv, r, branch }) => (
            p(cx, m * cx + b, cv),
            p(r, m * r, 0.),
            p(0., 0., branch * r),
            cx,
            r,
            m,
        ),
        (
            Graph::Circle {
                cx,
                cv,
                r: rho,
                branch: yb,
            },
            Graph::Circle {
                cx: zx,
                cv: zc,
                r,
                branch: zb,
            },
        ) => {
            let y_chart = Curve::CylinderPair {
                // Choose the cutter axis anchor near this actual intersection,
                // rather than an arbitrary world-zero height.
                c: p(cx, cv, zc),
                a: p(0., 0., yb),
                x: p(1., 0., 0.),
                rho,
                p: p(zx, cv, zc),
                b: p(0., 1., 0.),
                r,
                branch: yb * zb,
                t0: circle_angle(cx, rho, hi),
                t1: circle_angle(cx, rho, lo),
            };
            let z_chart = Curve::CylinderPair {
                c: p(zx, cv, zc),
                a: p(0., -zb, 0.),
                x: p(1., 0., 0.),
                rho: r,
                p: p(cx, cv, zc),
                b: p(0., 0., 1.),
                r: rho,
                branch: -zb * yb,
                t0: circle_angle(zx, r, hi),
                t1: circle_angle(zx, r, lo),
            };
            // Choose a smooth analytic parameter chart before constructing an
            // intersection edge. Both describe the same two supporting surfaces.
            let smooth = |curve: &Curve| {
                let Curve::CylinderPair {
                    c,
                    a,
                    x,
                    rho,
                    p,
                    b,
                    r,
                    t0,
                    t1,
                    ..
                } = *curve
                else {
                    unreachable!()
                };
                // Smooth over the trim; ending where the cylinders touch within tolerance is allowed (Bob 2026-10-28).
                crate::geom::pair_trim_touches(c, a, x, rho, p, b, r, t0, t1, tol).is_some()
            };
            let curve = if smooth(&y_chart) {
                y_chart
            } else if smooth(&z_chart) {
                z_chart
            } else {
                return Err(unsupported(
                    "neither analytic cylinder-intersection chart excludes tangency",
                ));
            };
            return Ok(curve.reversed());
        }
    };
    let t0 = circle_angle(cx, r, hi);
    let t1 = circle_angle(cx, r, lo);
    let curve = if slope == 0. {
        Curve::Arc {
            c: center,
            n: u.cross(v).unit()?,
            x: u.scale(1. / r),
            r,
            t0,
            t1,
        }
    } else {
        Curve::Conic {
            c: center,
            u,
            v,
            t0,
            t1,
        }
    };
    Ok(curve.reversed())
}
fn surface(graph: Graph, is_z: bool, positive: bool) -> Result<(Surface, bool)> {
    match graph {
        Graph::Line { m, b } => {
            let point = if is_z { p(0., 0., b) } else { p(0., b, 0.) };
            let normal = if is_z { p(-m, 0., 1.) } else { p(-m, 1., 0.) };
            let normal = normal.unit()?.scale(sign(positive));
            Ok((
                Surface::Plane {
                    p: point,
                    n: normal,
                    x: normal.perpendicular(),
                },
                true,
            ))
        }
        Graph::Circle { cx, cv, r, branch } => {
            let point = if is_z { p(cx, 0., cv) } else { p(cx, cv, 0.) };
            let axis = if is_z { p(0., 1., 0.) } else { p(0., 0., 1.) };
            Ok((
                Surface::Cylinder {
                    p: point,
                    a: axis,
                    x: p(1., 0., 0.),
                    r,
                },
                positive == (branch > 0.),
            ))
        }
    }
}
// Order: lower/upper Z, then lower/upper Y. Internal cell faces are omitted
// before constructing any curves, so no artificial internal edges are made.
fn selected_patches(b: &Band, exposed: [bool; 4], tol: f64) -> Result<Vec<Patch>> {
    let mut result = Vec::new();
    let (lo, hi) = (b.lo, b.hi);
    for (k, (zg, positive)) in [(b.z0, false), (b.z1, true)].into_iter().enumerate() {
        if !exposed[k] {
            continue;
        }
        let boundary = vec![
            lifted(lo, hi, b.y0, zg, tol)?,
            segment(
                p(hi, b.y0.value(hi), zg.value(hi)),
                p(hi, b.y1.value(hi), zg.value(hi)),
            ),
            lifted(lo, hi, b.y1, zg, tol)?.reversed(),
            segment(
                p(lo, b.y1.value(lo), zg.value(lo)),
                p(lo, b.y0.value(lo), zg.value(lo)),
            ),
        ];
        let (surface, sense) = surface(zg, true, positive)?;
        result.push(Patch {
            surface,
            sense,
            boundary: if positive {
                boundary
            } else {
                reversed(boundary)
            },
        });
    }
    for (k, (yg, positive)) in [(b.y0, false), (b.y1, true)].into_iter().enumerate() {
        if !exposed[k + 2] {
            continue;
        }
        let boundary = vec![
            lifted(lo, hi, yg, b.z0, tol)?,
            segment(
                p(hi, yg.value(hi), b.z0.value(hi)),
                p(hi, yg.value(hi), b.z1.value(hi)),
            ),
            lifted(lo, hi, yg, b.z1, tol)?.reversed(),
            segment(
                p(lo, yg.value(lo), b.z1.value(lo)),
                p(lo, yg.value(lo), b.z0.value(lo)),
            ),
        ];
        let (surface, sense) = surface(yg, false, positive)?;
        result.push(Patch {
            surface,
            sense,
            boundary: if positive {
                reversed(boundary)
            } else {
                boundary
            },
        });
    }
    Ok(result)
}
#[derive(Clone, Copy)]
struct Rect {
    y0: f64,
    y1: f64,
    z0: f64,
    z1: f64,
}
impl Rect {
    fn at(b: &Band, x: f64) -> Self {
        Self {
            y0: b.y0.value(x),
            y1: b.y1.value(x),
            z0: b.z0.value(x),
            z1: b.z1.value(x),
        }
    }
    fn contains(self, y: f64, z: f64) -> bool {
        y > self.y0 && y < self.y1 && z > self.z0 && z < self.z1
    }
}
fn end_patches(left: &[Band], right: &[Band], x: f64, tol: f64) -> Result<Vec<Patch>> {
    let left: Vec<_> = left.iter().map(|b| Rect::at(b, x)).collect();
    let right: Vec<_> = right.iter().map(|b| Rect::at(b, x)).collect();
    let mut ys = Vec::new();
    let mut zs = Vec::new();
    for r in left.iter().chain(&right) {
        ys.extend([r.y0, r.y1]);
        zs.extend([r.z0, r.z1]);
    }
    let ys = stations(ys, tol)?;
    let zs = stations(zs, tol)?;
    let mut result = Vec::new();
    for y in ys.windows(2) {
        for z in zs.windows(2) {
            let midy = (y[0] + y[1]) * 0.5;
            let midz = (z[0] + z[1]) * 0.5;
            let l = left.iter().any(|r| r.contains(midy, midz));
            let r = right.iter().any(|r| r.contains(midy, midz));
            if l == r {
                continue;
            }
            let q = [
                p(x, y[0], z[0]),
                p(x, y[1], z[0]),
                p(x, y[1], z[1]),
                p(x, y[0], z[1]),
            ];
            let boundary = (0..4).map(|i| segment(q[i], q[(i + 1) % 4])).collect();
            result.push(Patch {
                surface: Surface::Plane {
                    p: q[0],
                    n: p(if l { 1. } else { -1. }, 0., 0.),
                    x: p(0., 1., 0.),
                },
                sense: true,
                boundary: if l { boundary } else { reversed(boundary) },
            });
        }
    }
    Ok(result)
}

/// Verify the faces parallel to the partition planes too: their oriented
/// coverage must equal the difference of the adjoining material rectangles.
/// Cross-section sampling alone would otherwise miss duplicate/internal caps.
fn verify_interfaces(body: &Body, xs: &[f64], bands: &[Vec<Band>], tol: f64) -> Result<()> {
    for face in &body.faces {
        if matches!(face.surface,Surface::Plane {n,..} if n.y==0.&&n.z==0.) {
            continue;
        }
        let seen = xs.windows(2).any(|w| {
            let mid = w[0] + (w[1] - w[0]) * 0.5;
            face.loops
                .iter()
                .flatten()
                .any(|(i, _)| !roots(&body.edges[*i].curve, mid, tol).is_empty())
        });
        if !seen {
            return Err(unsupported(
                "an end-cut surface patch is below the admitted partition resolution",
            ));
        }
    }
    let mut faces = Vec::new();
    for face in &body.faces {
        let Surface::Plane { p: anchor, n, .. } = face.surface else {
            continue;
        };
        if n.y != 0. || n.z != 0. {
            continue;
        }
        if face.loops.len() != 1 {
            return Err(unsupported(
                "end-cut transverse faces must be rectangular patches",
            ));
        }
        let mut q = Vec::new();
        let mut area = 0.;
        for (i, forward) in &face.loops[0] {
            let Curve::Line { a, b } = body.edges[*i].curve else {
                return Err(unsupported("end-cut transverse face has a curved boundary"));
            };
            if (a.y - b.y).abs() > tol && (a.z - b.z).abs() > tol {
                return Err(unsupported(
                    "end-cut transverse face has an oblique boundary",
                ));
            }
            area += 0.5 * (a.y * b.z - a.z * b.y) * sign(*forward) * n.x * sign(face.sense);
            q.extend([a, b]);
        }
        let r = Rect {
            y0: q.iter().map(|p| p.y).fold(f64::INFINITY, f64::min),
            y1: q.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max),
            z0: q.iter().map(|p| p.z).fold(f64::INFINITY, f64::min),
            z1: q.iter().map(|p| p.z).fold(f64::NEG_INFINITY, f64::max),
        };
        if (area - (r.y1 - r.y0) * (r.z1 - r.z0)).abs() > tol * 2. * ((r.y1 - r.y0) + (r.z1 - r.z0))
            || area <= 0.
        {
            return Err(invalid(
                "end-cut transverse patch orientation or coverage is invalid",
            ));
        }
        for (i, _) in &face.loops[0] {
            let c = &body.edges[*i].curve;
            let mid = c.start().add(c.end()).scale(0.5);
            if (mid.y - r.y0).abs() > tol
                && (mid.y - r.y1).abs() > tol
                && (mid.z - r.z0).abs() > tol
                && (mid.z - r.z1).abs() > tol
            {
                return Err(unsupported("end-cut transverse patch is not rectangular"));
            }
        }
        let station = xs
            .iter()
            .position(|x| (*x - anchor.x).abs() <= tol)
            .ok_or_else(|| invalid("transverse face has no matching partition station"))?;
        faces.push((station, r, if n.x * sign(face.sense) > 0. { 1 } else { -1 }));
    }
    let mut work = 0usize;
    for (i, x) in xs.iter().enumerate() {
        let left: Vec<_> = if i > 0 {
            bands[i - 1].iter().map(|b| Rect::at(b, *x)).collect()
        } else {
            Vec::new()
        };
        let right: Vec<_> = bands
            .get(i)
            .into_iter()
            .flatten()
            .map(|b| Rect::at(b, *x))
            .collect();
        let patches: Vec<(Rect, i32)> = faces
            .iter()
            .filter(|(station, _, _)| *station == i)
            .map(|(_, r, s)| (*r, *s))
            .collect();
        let mut ys = Vec::new();
        let mut zs = Vec::new();
        for r in left
            .iter()
            .chain(&right)
            .copied()
            .chain(patches.iter().map(|(r, _)| *r))
        {
            ys.extend([r.y0, r.y1]);
            zs.extend([r.z0, r.z1]);
        }
        let ys = stations(ys, tol)?;
        let zs = stations(zs, tol)?;
        for y in ys.windows(2) {
            for z in zs.windows(2) {
                work += patches.len() + 1;
                if work > 2_000_000 {
                    return Err(unsupported(
                        "end-cut transverse coverage exceeds its work limit",
                    ));
                }
                let my = (y[0] + y[1]) * 0.5;
                let mz = (z[0] + z[1]) * 0.5;
                let expected = left.iter().any(|r| r.contains(my, mz)) as i32
                    - right.iter().any(|r| r.contains(my, mz)) as i32;
                let actual: Vec<i32> = patches
                    .iter()
                    .filter_map(|(r, s)| r.contains(my, mz).then_some(*s))
                    .collect();
                if actual.len() != expected.unsigned_abs() as usize
                    || actual.iter().sum::<i32>() != expected
                {
                    return Err(invalid(
                        "end-cut transverse boundary differs from its material sections",
                    ));
                }
            }
        }
    }
    Ok(())
}

fn parameter(curve: &Curve, q: Vec3, tol: f64) -> Option<f64> {
    // A shared B-rep vertex is allowed to differ from an analytic endpoint by
    // the declared input tolerance; classify endpoints before interior roots.
    let (lo, hi) = domain(curve);
    if q.distance(curve.start()) <= tol {
        return Some(lo);
    }
    if q.distance(curve.end()) <= tol {
        return Some(hi);
    }
    match *curve {
        Curve::Line { a, b } => {
            let d = b.sub(a);
            let t = q.sub(a).dot(d) / d.dot(d);
            (t >= -1e-12 && t <= 1. + 1e-12 && q.distance(eval(curve, t.clamp(0., 1.))) <= tol)
                .then_some(t.clamp(0., 1.))
        }
        Curve::CylinderPair {
            c, a, x, t0, t1, ..
        } => {
            let d = q.sub(c);
            parameters(d.dot(a.cross(x)).atan2(d.dot(x)), t0, t1)
                .into_iter()
                .find(|t| eval(curve, *t).distance(q) <= tol)
        }
        _ => {
            let (c, u, v, t0, t1) = conic(curve).unwrap();
            let d = q.sub(c);
            let uu = u.dot(u);
            let uv = u.dot(v);
            let vv = v.dot(v);
            let det = uu * vv - uv * uv;
            let a = (d.dot(u) * vv - d.dot(v) * uv) / det;
            let b = (d.dot(v) * uu - d.dot(u) * uv) / det;
            parameters(b.atan2(a), t0, t1)
                .into_iter()
                .find(|t| eval(curve, *t).distance(q) <= tol)
        }
    }
}
fn slice(curve: &Curve, t0: f64, t1: f64) -> Curve {
    match *curve {
        Curve::RationalBezier { p, w } => {
            let (p, w) = crate::rational::piece(p, w, t0, t1);
            Curve::RationalBezier { p, w }
        },
        Curve::Bezier { p } => Curve::Bezier { p: crate::geom::bezier_piece(p, t0, t1) },
        Curve::Line { .. } => segment(eval(curve, t0), eval(curve, t1)),
        Curve::Arc { c, n, x, r, .. } => Curve::Arc { c, n, x, r, t0, t1 },
        Curve::Conic { c, u, v, .. } => Curve::Conic { c, u, v, t0, t1 },
        Curve::CylinderPair {
            c,
            a,
            x,
            rho,
            p,
            b,
            r,
            branch,
            ..
        } => Curve::CylinderPair {
            c,
            a,
            x,
            rho,
            p,
            b,
            r,
            branch,
            t0,
            t1,
        },
        Curve::RadialPair { c, a, x, rho, p, b, r, branch, s, .. } => Curve::RadialPair { c, a, x, rho, p, b, r, branch, t0, t1, s },
        Curve::ConePair { c, a, x, rho, p, b, r, k, branch, .. } => Curve::ConePair { c, a, x, rho, p, b, r, k, branch, t0, t1 },
    }
}
fn same_curve(a: &Curve, b: &Curve, reverse: bool, tol: f64) -> bool {
    if std::mem::discriminant(a) != std::mem::discriminant(b) {
        return false;
    }
    let (a0, a1) = domain(a);
    let (b0, b1) = domain(b);
    [0., 0.25, 0.5, 0.75, 1.].iter().all(|s| {
        let t = if reverse { 1. - s } else { *s };
        eval(a, a0 + (a1 - a0) * s).distance(eval(b, b0 + (b1 - b0) * t)) <= tol
    })
}
fn assemble(patches: Vec<Patch>, tol: f64) -> Result<Body> {
    let mut body = Body::default();
    for patch in &patches {
        for curve in &patch.boundary {
            for q in [curve.start(), curve.end()] {
                let candidates: Vec<usize> = body
                    .vertices
                    .iter()
                    .enumerate()
                    .filter_map(|(i, v)| (v.distance(q) <= tol).then_some(i))
                    .collect();
                if candidates.len() > 1 {
                    return Err(unsupported(
                        "ambiguous end-cut vertex coincidence within tolerance",
                    ));
                }
                if candidates.is_empty() {
                    body.vertices.push(q);
                }
            }
        }
    }
    let curve_count: usize = patches.iter().map(|p| p.boundary.len()).sum();
    if body.vertices.len().saturating_mul(curve_count) > 20_000_000 {
        return Err(unsupported("end-cut edge assembly exceeds its work limit"));
    }
    let mut edges: std::collections::HashMap<(usize, usize), Vec<usize>> =
        std::collections::HashMap::new();
    for patch in patches {
        let mut boundary = Vec::new();
        for curve in patch.boundary {
            let mut cuts: Vec<(f64, usize)> = body
                .vertices
                .iter()
                .enumerate()
                .filter_map(|(i, q)| parameter(&curve, *q, tol).map(|t| (t, i)))
                .collect();
            cuts.sort_by(|a, b| a.0.total_cmp(&b.0));
            // At a circle extremum a ruled patch can have a collapsed generator.
            // Its endpoint is already represented by the two adjacent curves.
            if cuts.len() == 1
                && matches!(curve, Curve::Line { .. })
                && curve.start().distance(curve.end()) <= tol
            {
                continue;
            }
            if cuts.len() < 2 {
                return Err(KernelError::Internal(format!(
                    "end-cut edge lost an endpoint: {curve:?}; matches {cuts:?}"
                )));
            }
            for pair in cuts.windows(2) {
                let (t0, a) = pair[0];
                let (t1, b) = pair[1];
                if a == b || t1 <= t0 {
                    continue;
                }
                let piece = slice(&curve, t0, t1);
                let key = (a.min(b), a.max(b));
                let mut found = None;
                if let Some(candidates) = edges.get(&key) {
                    for index in candidates {
                        let edge = &body.edges[*index];
                        let reverse = edge.start != a;
                        if same_curve(&piece, &edge.curve, reverse, tol) {
                            if found.is_some() {
                                return Err(unsupported("ambiguous coincident end-cut edges"));
                            }
                            found = Some((*index, !reverse));
                        }
                    }
                }
                let entry = if let Some(found) = found {
                    found
                } else {
                    let index = body.edges.len();
                    body.edges.push(Edge {
                        curve: piece,
                        start: a,
                        end: b,
                    });
                    edges.entry(key).or_default().push(index);
                    (index, true)
                };
                boundary.push(entry);
            }
        }
        body.faces.push(Face {
            surface: patch.surface,
            sense: patch.sense,
            loops: vec![boundary],
        });
    }
    crate::topology::check(&body, tol)
        .map_err(|e| KernelError::Internal(format!("end-cut assembly is not closed: {e}")))?;
    validate(&body, tol).map_err(|e| {
        KernelError::Internal(format!("end-cut result failed analytic checks: {e}"))
    })?;
    Ok(body)
}

#[derive(Clone, Copy)]
struct Frame {
    u: Vec3,
    v: Vec3,
    w: Vec3,
}
impl Frame {
    fn local(self, p: Vec3) -> Vec3 {
        crate::geom::Vec3::new(p.dot(self.u), p.dot(self.v), p.dot(self.w))
    }
    fn world(self, p: Vec3) -> Vec3 {
        self.u
            .scale(p.x)
            .add(self.v.scale(p.y))
            .add(self.w.scale(p.z))
    }
}
fn axis(v: Vec3) -> Result<usize> {
    let a = xyz(v);
    let nonzero: Vec<usize> = (0..3).filter(|i| a[*i] != 0.).collect();
    if nonzero.len() != 1 || a[nonzero[0]].abs() != 1. {
        return Err(unsupported(
            "transverse end-cut cylinder axes must be cardinal",
        ));
    }
    Ok(nonzero[0])
}
fn frame(body: &Body, tool: &Body) -> Result<Frame> {
    let mut tool_axis = None;
    for f in &tool.faces {
        if let Surface::Cylinder { a, .. } = f.surface {
            let a = axis(a)?;
            if tool_axis.is_some_and(|old| old != a) {
                return Err(unsupported("end-cut tool must have one sweep axis"));
            }
            tool_axis = Some(a);
        }
    }
    let w = tool_axis.ok_or_else(|| {
        unsupported("end-cut sweep axis is not specified by a cylindrical tool wall")
    })?;
    frame_for_axis_with(body, Some(tool), w)
}
#[cfg(test)]
fn frame_for_axis(body: &Body, w: usize) -> Result<Frame> {
    frame_for_axis_with(body, None, w)
}
/// The beam axis v is the one longitudinal cylinder axis other than the sweep w. With two (a pipe already coped
/// across another axis, Bob's Pipe 1.25 18 Copes Rotated 90), v is the one for which the other axis's faces lie
/// wholly beyond the tool along v, so the cut runs on a split piece (`split`); chosen from the geometry before
/// the cut, never by retrying a failed one.
fn frame_for_axis_with(body: &Body, tool: Option<&Body>, w: usize) -> Result<Frame> {
    let mut axes = Vec::new();
    for f in &body.faces {
        if let Surface::Cylinder { a, .. } = f.surface {
            let a = axis(a)?;
            if a != w && !axes.contains(&a) {
                axes.push(a);
            }
        }
    }
    let unit = |i: usize| match i {
        0 => p(1., 0., 0.),
        1 => p(0., 1., 0.),
        _ => p(0., 0., 1.),
    };
    let make = |v: usize| Frame { u: unit(v).cross(unit(w)), v: unit(v), w: unit(w) };
    let v = match (axes.as_slice(), tool) {
        (&[v], _) => v,
        (&[], _) => return Err(unsupported("end-cut base has no perpendicular longitudinal cylinder axis")),
        (&[a, b], Some(tool)) => {
            let fits = |v: usize| {
                let f = make(v);
                split::separable(&transform(body, |q| f.local(q)), &transform(tool, |q| f.local(q)))
            };
            match (fits(a), fits(b)) {
                (true, false) => a,
                (false, true) => b,
                _ => return Err(unsupported("end-cut base has more than one longitudinal cylinder axis")),
            }
        }
        _ => return Err(unsupported("end-cut base has more than one longitudinal cylinder axis")),
    };
    Ok(make(v))
}
fn transform(body: &Body, map: impl Fn(Vec3) -> Vec3) -> Body {
    let mut result = body.clone();
    for p in &mut result.vertices {
        *p = map(*p);
    }
    for e in &mut result.edges {
        match &mut e.curve {
            Curve::Line { a, b } => {
                *a = map(*a);
                *b = map(*b);
            }
            Curve::Arc { c, n, x, .. } => {
                *c = map(*c);
                *n = map(*n);
                *x = map(*x);
            }
            Curve::Conic { c, u, v, .. } => {
                *c = map(*c);
                *u = map(*u);
                *v = map(*v);
            }
            Curve::CylinderPair { c, a, x, p, b, .. } | Curve::RadialPair { c, a, x, p, b, .. } | Curve::ConePair { c, a, x, p, b, .. } => {
                *c = map(*c);
                *a = map(*a);
                *x = map(*x);
                *p = map(*p);
                *b = map(*b);
            }
            Curve::Bezier { p } | Curve::RationalBezier { p, .. } => { for q in p.iter_mut() { *q = map(*q); } }
        }
    }
    for f in &mut result.faces {
        match &mut f.surface {
            Surface::Spline { data } => **data = data.mapped(&map),
            Surface::Plane { p, n, x } => {
                *p = map(*p);
                *n = map(*n);
                *x = map(*x);
            }
            Surface::Cylinder { p, a, x, .. } => {
                *p = map(*p);
                *a = map(*a);
                *x = map(*x);
            }
            Surface::EllipticCylinder { p, a, u, v } => {
                *p = map(*p);
                *a = map(*a);
                *u = map(*u);
                *v = map(*v);
            }
            Surface::Cone { p, a, x, .. } | Surface::Sphere { c: p, a, x, .. } | Surface::Torus { c: p, a, x, .. } | Surface::Revolution { p, a, x, .. } => {
                *p = map(*p);
                *a = map(*a);
                *x = map(*x);
            }
            Surface::Conoid { c, a, x, p, b, .. } => {
                *c = map(*c);
                *a = map(*a);
                *x = map(*x);
                *p = map(*p);
                *b = map(*b);
            }
        }
    }
    result
}

/// Geometric domain selection, before either cut executes. A cutter corridor
/// overlapping a perpendicular cylindrical patch needs transverse intersections.
/// No operation error is caught or retried through another algorithm.
pub(super) fn crosses_cylinder(
    body: &Body,
    axis: usize,
    center: Vec3,
    radius: f64,
    low: f64,
    high: f64,
    tol: f64,
) -> bool {
    if body
        .edges
        .iter()
        .any(|e| matches!(e.curve, Curve::CylinderPair { .. }))
    {
        return true;
    }
    let center = xyz(center);
    let mut clo = center.map(|c| c - radius);
    let mut chi = center.map(|c| c + radius);
    clo[axis] = low;
    chi[axis] = high;
    for face in &body.faces {
        let Surface::Cylinder { a, .. } = face.surface else {
            continue;
        };
        if xyz(a)[axis] != 0. {
            continue;
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for (i, _) in face.loops.iter().flatten() {
            let Some(e) = body.edges.get(*i) else {
                return true;
            };
            let (l, h) = super::passage::curve_bounds(&e.curve);
            for k in 0..3 {
                lo[k] = lo[k].min(l[k]);
                hi[k] = hi[k].max(h[k]);
            }
        }
        if (0..3).all(|k| hi[k] >= clo[k] - tol && lo[k] <= chi[k] + tol) {
            return true;
        }
    }
    false
}
/// One direction with every component within 1e-12 of zero made exactly zero, renormalized. Stock built by rotated,
/// tapered or filleted profile math carries plane normals and axes off the coordinate axes in the last bits (a channel
/// flange face's 1.3e-16, a C6x13 inner web's 2.6e-17); the domains read directions exactly, and over any admitted
/// part such a direction is the same within the kernel tolerance (Bob 10-28: what agrees within the kernel tolerance is
/// the same). Points never move.
fn snapped(v: Vec3) -> Vec3 {
    let k = |c: f64| if c.abs() <= 1e-12 { 0. } else { c };
    let w = Vec3::new(k(v.x), k(v.y), k(v.z));
    let l = w.length();
    if l > 0. && w != v { w.scale(1. / l) } else { v }
}

/// A transverse-cut operand with its planes' and cylinders' directions snapped (see `snapped`): the domain's charts and
/// partitions read directions exactly. Applied inside the domain, after selection, so routing is unchanged.
pub(super) fn snap_directions(body: &Body) -> Body {
    let mut out = body.clone();
    for f in &mut out.faces {
        match &mut f.surface {
            Surface::Plane { n, x, .. } => (*n, *x) = (snapped(*n), snapped(*x)),
            Surface::Cylinder { a, x, .. } => (*a, *x) = (snapped(*a), snapped(*x)),
            _ => {}
        }
    }
    out
}

pub(super) fn cross_cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    cross::cut(body, tool, tol)
}
pub(super) fn box_cross_cut(body: &Body, tool: &Body, axis: usize, tol: f64) -> Result<Body> {
    cross::cut_along(body, tool, axis, tol)
}

/// Select a transverse coordinate chart from the planar cutter's normals.
/// This is geometric admission before construction, with no retry if it fails.
pub(super) fn planar_cross_cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    let (body, tool) = (&snap_directions(body), &snap_directions(tool));
    let mut normals = Vec::new();
    for face in &tool.faces {
        let Surface::Plane { p, n, x } = face.surface else {
            return Err(unsupported("planar section cutter has a curved face"));
        };
        if !p.is_finite() || !n.is_finite() || !x.is_finite() || (n.length() - 1.).abs() > 1e-10 {
            return Err(invalid("invalid planar section cutter frame"));
        }
        normals.push(xyz(n));
    }
    let mut cylinder_axes = Vec::new();
    for face in &body.faces {
        if let Surface::Cylinder { a, .. } = face.surface {
            let k = axis(a)?;
            if !cylinder_axes.contains(&k) {
                cylinder_axes.push(k);
            }
        }
    }
    // The stock's own planes must be section graphs in the chart too: a plane turned in both the beam direction v
    // and the sweep w (Bob's earlier miter, normal in X and Y, under a later cardinal cope box) has no axial span
    // there. Decided from the operands before the cut, never by retrying another chart after a failure.
    let body_normals: Vec<[f64; 3]> = body
        .faces
        .iter()
        .filter_map(|f| match f.surface {
            Surface::Plane { n, .. } => Some(xyz(n)),
            _ => None,
        })
        .collect();
    let chosen = (0..3)
        .find(|w| {
            let others: Vec<usize> = cylinder_axes.iter().copied().filter(|k| *k != *w).collect();
            others.len() == 1
                && normals.iter().all(|n| {
                    n[*w] == 0. || ((0..3).all(|k| k == *w || n[k] == 0.) && n[*w].abs() == 1.)
                })
                && body_normals.iter().all(|n| n[others[0]] == 0. || n[*w] == 0.)
        })
        .ok_or_else(|| {
            unsupported("no cardinal section chart for the planar cutter and stock cylinders")
        })?;
    cross::cut_along(body, tool, chosen, tol)
}

#[cfg(test)]
mod tests {
    use super::cross::cut;
    use super::*;
    fn derivative(curve: &Curve, t: f64) -> Vec3 {
        curve.eval(t).1
    }
    // Divergence theorem on the RETURNED boundary. Planar Green integrals are
    // exact; cylinder patches use a boundary integral in angle/axial coordinates.
    // This does not use bands or the body's original cutting recipe.
    fn volume(body: &Body) -> f64 {
        let mut total = 0.;
        // Use a nearby divergence-field origin to avoid cancellation of large
        // coordinate terms when the complete part is translated far from zero.
        let origin = body.vertices[0];
        for face in &body.faces {
            match face.surface {
                Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => unreachable!("no elliptic-cylinder or revolved faces in these test bodies"),
                Surface::Plane { p: anchor, n, .. } => {
                    let outward = n.scale(sign(face.sense));
                    let mut area = 0.;
                    for (i, fwd) in face.loops.iter().flatten() {
                        let curve = &body.edges[*i].curve;
                        let integral = if let Some((c, u, v, t0, t1)) = conic(curve) {
                            c.sub(origin)
                                .cross(curve.end().sub(curve.start()))
                                .add(u.cross(v).scale(t1 - t0))
                        } else {
                            curve.start().sub(origin).cross(curve.end().sub(origin))
                        };
                        area += 0.5 * integral.dot(outward) * sign(*fwd);
                    }
                    assert!(area > 0.);
                    total += anchor.sub(origin).dot(outward) * area / 3.;
                }
                Surface::Cylinder {
                    p: center, a, r, ..
                } => {
                    for (i, fwd) in face.loops.iter().flatten() {
                        let curve = &body.edges[*i].curve;
                        let (t0, t1) = domain(curve);
                        let steps = 1024usize;
                        let f = |t| {
                            let point = eval(curve, t);
                            let d = point.sub(center);
                            let along = point.sub(origin).dot(a);
                            let radial = d.sub(a.scale(d.dot(a)));
                            let delta = derivative(curve, t);
                            let delta = delta.sub(a.scale(delta.dot(a)));
                            -along
                                * (center.sub(origin).dot(radial) + r * r)
                                * a.dot(radial.cross(delta))
                                / (r * r)
                        };
                        let mut integral = f(t0) + f(t1);
                        for j in 1..steps {
                            integral += if j % 2 == 0 { 2. } else { 4. }
                                * f(t0 + (t1 - t0) * j as f64 / steps as f64);
                        }
                        integral *= (t1 - t0) / (3. * steps as f64);
                        total += sign(*fwd) * integral / 3.;
                    }
                }
            }
        }
        total
    }
    fn fixture(name: &str) -> Body {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/rat_hole")
            .join(name);
        let bytes = std::fs::read(path).unwrap();
        let data: Vec<f64> = bytes
            .chunks_exact(8)
            .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
            .collect();
        crate::neutral::Reader::new(&data).body().unwrap()
    }
    #[test]
    fn actual_rat_hole_is_admitted_as_transverse_analytic_bands() {
        let body = fixture("base.f64");
        let tool = fixture("step01.f64");
        let tol = 1e-6;
        validate(&body, tol).unwrap();
        validate(&tool, tol).unwrap();
        let mut values = events(&body);
        values.extend(events(&tool));
        let xs = stations(values, tol).unwrap();
        let mut count = 0;
        for w in xs.windows(2) {
            if let Some(b) = section(&body, w[0], w[1], tol).unwrap() {
                let t = section(&tool, w[0], w[1], tol).unwrap();
                let kept = subtract(b, t, tol).unwrap();
                assert!(kept.is_some());
                count += 1;
            }
        }
        assert!(count >= 7);
    }

    #[test]
    fn both_actual_rat_holes_replay_as_closed_analytic_solids() {
        let base = fixture("base.f64");
        let first = cut(&base, &fixture("step01.f64"), 1e-6).unwrap();
        crate::topology::check(&first, 1e-6).unwrap();
        let second = cut(&first, &fixture("step02.f64"), 1e-6).unwrap();
        crate::topology::check(&second, 1e-6).unwrap();
        // The rat holes' plane-fillet crossings: circles (arcs) since the domain takes the captured near-cardinal
        // planes as exact (snap_directions), conics before.
        let curved = |b: &Body| b.edges.iter().filter(|e| matches!(e.curve, Curve::Conic { .. } | Curve::Arc { .. })).count();
        assert!(curved(&second) > curved(&base));
        assert_eq!(crate::topology::check(&second, 1e-6).unwrap().genus(), 0.);
        // Independent dimensions/area formula for this captured fixture.
        let stock =
            (2. * 0.24 * 4. + (9.99 - 2. * 0.24) * 0.23 + 4. * 0.5_f64.powi(2) * (1. - PI / 4.))
                * 36.;
        let front = 2.580302226686857;
        let removed = front * (1.23 * 0.5 - 2. * PI * 0.5_f64.powi(2) / 4. + 0.25 * 0.23)
            + 0.23 * (0.5 * (front - 0.5) + PI * 0.5_f64.powi(2) / 4.);
        assert!(
            (volume(&base) - stock).abs() < 1e-8,
            "base {} expected {stock}",
            volume(&base)
        );
        assert!(
            (volume(&first) - (stock - removed)).abs() < 1e-8,
            "first {} expected {}",
            volume(&first),
            stock - removed
        );
        assert!(
            (volume(&second) - (stock - 2. * removed)).abs() < 1e-8,
            "second {} expected {}",
            volume(&second),
            stock - 2. * removed
        );
    }

    #[test]
    fn returned_rat_holes_match_production_material_probes() {
        let first = cut(&fixture("base.f64"), &fixture("step01.f64"), 1e-6).unwrap();
        let result = cut(&first, &fixture("step02.f64"), 1e-6).unwrap();
        let manifest = std::fs::read_to_string(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures/rat_hole/manifest.txt"),
        )
        .unwrap();
        let mut checked = 0;
        let mut near_boundary = 0;
        for line in manifest.lines() {
            let Some((coordinates, label)) = line.trim().split_once(": ") else {
                continue;
            };
            if label != "Inside" && label != "Outside" {
                continue;
            }
            let numbers: Vec<f64> = coordinates
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if numbers.len() != 3 {
                continue;
            }
            let q = p(numbers[0], numbers[1], numbers[2]);
            let xs = stations(events(&result), 1e-6).unwrap();
            if xs.iter().any(|x| (*x - q.x).abs() < 2e-6) {
                near_boundary += 1;
                continue;
            }
            let interval = xs.windows(2).find(|w| q.x > w[0] && q.x < w[1]);
            let band = interval.and_then(|w| section(&result, w[0], w[1], 1e-6).unwrap());
            let inside = if let Some(b) = band {
                let r = Rect::at(&b, q.x);
                if [
                    (q.y - r.y0).abs(),
                    (q.y - r.y1).abs(),
                    (q.z - r.z0).abs(),
                    (q.z - r.z1).abs(),
                ]
                .into_iter()
                .any(|v| v < 2e-6)
                {
                    near_boundary += 1;
                    continue;
                }
                r.contains(q.y, q.z)
            } else {
                false
            };
            assert_eq!(inside, label == "Inside", "probe {line}");
            checked += 1;
        }
        assert!(
            checked > 400,
            "only {checked} checked; {near_boundary} near boundaries"
        );
        eprintln!("rat-hole material probes: {checked} checked, {near_boundary} excluded within coordinate/boundary tolerance");
    }

    fn replay(base: Body, first: Body, second: Body, tol: f64) -> Body {
        let operation = |body: &Body, tool: &Body| {
            let mut data = crate::neutral::write_body(body);
            data.extend(crate::neutral::write_body(tool));
            let request = crate::protocol::parse_request(&crate::protocol::encode_request(
                crate::protocol::OP_CUT,
                tol,
                &data,
            ))
            .unwrap();
            let output = crate::protocol::handle(&request).unwrap();
            let mut reader = crate::neutral::Reader::new(&output);
            let body = reader.body().unwrap();
            reader.finish().unwrap();
            body
        };
        operation(&operation(&base, &first), &second)
    }
    #[test]
    fn both_rat_holes_replay_through_public_wire_protocol_without_mutating_inputs() {
        let base = fixture("base.f64");
        let first = fixture("step01.f64");
        let second = fixture("step02.f64");
        let input = crate::neutral::write_body(&base);
        let tool = crate::neutral::write_body(&first);
        let result = replay(base.clone(), first.clone(), second, 1e-6);
        assert_eq!(crate::neutral::write_body(&base), input);
        assert_eq!(crate::neutral::write_body(&first), tool);
        assert_eq!(crate::topology::check(&result, 1e-6).unwrap().genus(), 0.);
        // The public cut, as the wire runs it (operand directions normalized at the boolean entry).
        let direct = crate::boolean::cut(
            &crate::boolean::cut(&base, &first, 1e-6).unwrap(),
            &fixture("step02.f64"),
            1e-6,
        )
        .unwrap();
        let actual = crate::neutral::write_body(&result);
        let expected = crate::neutral::write_body(&direct);
        assert_eq!(actual.len(), expected.len());
        let delta = actual
            .iter()
            .zip(&expected)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        assert!(
            delta < 1e-6,
            "wire/direct maximum scalar difference {delta}"
        );
        assert!((volume(&result) - volume(&direct)).abs() < 1e-8);
    }
    #[test]
    fn transverse_cut_is_invariant_under_all_six_cardinal_axis_frames() {
        let axes = [p(1., 0., 0.), p(0., 1., 0.), p(0., 0., 1.)];
        let expected = volume(&replay(
            fixture("base.f64"),
            fixture("step01.f64"),
            fixture("step02.f64"),
            1e-6,
        ));
        for v in 0..3 {
            for w in 0..3 {
                if v != w {
                    let f = Frame {
                        u: axes[v].cross(axes[w]),
                        v: axes[v],
                        w: axes[w],
                    };
                    let map = |name| transform(&fixture(name), |p| f.world(p));
                    let result =
                        replay(map("base.f64"), map("step01.f64"), map("step02.f64"), 1e-6);
                    assert!((volume(&result) - expected).abs() < 1e-8);
                }
            }
        }
    }
    fn scaled(body: &Body, scale: f64, offset: Vec3) -> Body {
        let mut out = body.clone();
        for v in &mut out.vertices {
            *v = v.scale(scale).add(offset);
        }
        for e in &mut out.edges {
            match &mut e.curve {
                Curve::Line { a, b } => {
                    *a = a.scale(scale).add(offset);
                    *b = b.scale(scale).add(offset);
                }
                Curve::Arc { c, r, .. } => {
                    *c = c.scale(scale).add(offset);
                    *r *= scale;
                }
                Curve::Conic { c, u, v, .. } => {
                    *c = c.scale(scale).add(offset);
                    *u = u.scale(scale);
                    *v = v.scale(scale);
                }
                Curve::CylinderPair { c, p, rho, r, .. } => {
                    *c = c.scale(scale).add(offset);
                    *p = p.scale(scale).add(offset);
                    *rho *= scale;
                    *r *= scale;
                }
                Curve::Bezier { .. } | Curve::RationalBezier { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } => unreachable!("no Bezier or radial-pair edges in these test bodies"),
            }
        }
        for f in &mut out.faces {
            match &mut f.surface {
                Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => unreachable!("no elliptic-cylinder or revolved faces in these test bodies"),
                Surface::Plane { p, .. } => *p = p.scale(scale).add(offset),
                Surface::Cylinder { p, r, .. } => {
                    *p = p.scale(scale).add(offset);
                    *r *= scale;
                }
            }
        }
        out
    }
    #[test]
    fn transverse_cut_preserves_units_and_translation() {
        let expected = volume(&replay(
            fixture("base.f64"),
            fixture("step01.f64"),
            fixture("step02.f64"),
            1e-6,
        ));
        for scale in [0.1, 25.4] {
            let map = |name| scaled(&fixture(name), scale, p(100., -50., 15.));
            let result = replay(
                map("base.f64"),
                map("step01.f64"),
                map("step02.f64"),
                scale * 1e-6,
            );
            let measured = volume(&result) / scale.powi(3);
            assert!(
                (measured - expected).abs() < 1e-6,
                "scale {scale}: measured {measured}, expected {expected}, difference {}",
                measured - expected
            );
        }
    }
    #[test]
    fn partial_depth_and_interior_copes_preserve_thin_skin_refusal() {
        let base = fixture("base.f64");
        let tool = fixture("step01.f64");
        let full_volume = 155.58846611769187;
        let removed = full_volume - 154.58210003470572;
        for (delta, fraction) in [(p(0., 0., 4.), 0.5), (p(0., 5., 0.), 1.)] {
            let moved = scaled(&tool, 1., delta);
            let result = cut(&base, &moved, 1e-6).unwrap();
            crate::topology::check(&result, 1e-6).unwrap();
            let measured = crate::mass::measure(&result, 1e-6).unwrap().volume;
            assert!((measured - (full_volume - removed * fraction)).abs() < 1e-9);
        }
        // A sub-tolerance retained end skin is still material. Admission must
        // reject that interior cut, not erase the skin by treating it as flush.
        let moved = scaled(&tool, 1., p(0., 0.25e-6, 0.));
        assert!(matches!(
            cut(&base, &moved, 1e-6),
            Err(KernelError::Unsupported(_))
        ));
        let mut malformed = base.clone();
        let face = malformed
            .faces
            .iter_mut()
            .find(|f| matches!(f.surface,Surface::Plane {n,..} if n.x.abs()==1.))
            .unwrap();
        face.sense = !face.sense;
        assert!(matches!(
            cut(&malformed, &tool, 1e-6),
            Err(KernelError::InvalidInput(_))
        ));
    }
}

#[cfg(test)]
mod codex_direction_snap_review {
    use super::*;
    use crate::frame::{Point3, Rigid};
    use crate::neutral::Profile;

    fn box_at_angle(length: f64, tol: f64) -> Body {
        let points = [p(0.,0.,0.),p(length,0.,0.),p(length,1.,0.),p(0.,1.,0.)];
        let profile = Profile { loops:vec![(0..4).map(|i| Curve::Line { a:points[i],b:points[(i+1)%4] }).collect()],
            n:p(0.,0.,1.),low:0.,high:1. };
        let body = crate::extrude::extrude(&profile,tol).unwrap();
        let rotated = Rigid::rotation(p(0.,0.,1.),5e-13,Point3::new(0.,0.,0.)).unwrap().body(&body);
        validate(&rotated,tol).unwrap();
        rotated
    }

    #[test]
    fn direction_normalization_must_respect_the_request_distance_tolerance() {
        let tol=1e-10;
        let original=box_at_angle(1000.,tol);
        let normalized=snap_directions(&original);
        let checked=crate::topology::check(&normalized,tol);
        assert!(checked.is_ok(),"valid input must not be changed beyond request tolerance by angular threshold: {checked:?}");
    }

    #[test]
    fn direction_normalization_must_bound_displacement_for_large_admitted_parts() {
        let tol=1e-6;
        let original=box_at_angle(1e8,tol);
        let normalized=snap_directions(&original);
        let checked=crate::topology::check(&normalized,tol);
        assert!(checked.is_ok(),"valid input must not be changed beyond default tolerance: {checked:?}");
    }

    #[test]
    fn a_small_face_can_have_a_distant_valid_plane_origin() {
        let tol=1e-6;
        let mut original=box_at_angle(1.,tol);
        for face in &mut original.faces {
            if let Surface::Plane { p:origin,x,.. }=&mut face.surface {
                *origin=origin.add(x.scale(1e8));
            }
        }
        // Moving a plane's origin within its own plane changes only its parameterization.
        validate(&original,tol).unwrap();
        let normalized=snap_directions(&original);
        let checked=crate::topology::check(&normalized,tol);
        assert!(checked.is_ok(),"a one-inch part also needs a spatial bound relative to each support origin: {checked:?}");
    }
}

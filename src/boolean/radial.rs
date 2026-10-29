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
fn cutter(tool: &Body, tol: f64) -> Option<(Vec3, Vec3, f64)> {
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
    admit(body, tool, tol).is_some_and(|(k, pipe)| clear(body, tool, &k, &pipe, tol).is_ok())
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
    for &face in [pipe.outer, pipe.inner].iter() {
        let f = &body.faces[face];
        let edges: Vec<usize> = f.loops.iter().flatten().map(|(e, _)| *e).collect();
        let mut rings: Vec<(f64, f64)> = vec![];
        for &e in &edges {
            let bx = bounds(&body.edges[e].curve).ok_or_else(|| unsupported("radial bore: a wall edge has no bounding box"))?;
            if !apart(bx, loop_box, tol) {
                return Err(unsupported("radial bore meets a boundary of the pipe wall (an end, its seam or another hole)"));
            }
            if !matches!(body.edges[e].curve, Curve::Line { .. }) {
                let along = |q: [f64; 3]| (q[0] - pipe.p.x) * b.x + (q[1] - pipe.p.y) * b.y + (q[2] - pipe.p.z) * b.z;
                let corners = (0..8).map(|m| [if m & 1 == 0 { bx.0[0] } else { bx.1[0] }, if m & 2 == 0 { bx.0[1] } else { bx.1[1] }, if m & 4 == 0 { bx.0[2] } else { bx.1[2] }]);
                let values: Vec<f64> = corners.map(along).collect();
                rings.push((values.iter().cloned().fold(f64::INFINITY, f64::min), values.iter().cloned().fold(f64::NEG_INFINITY, f64::max)));
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
            let bx = bounds(&body.edges[*e].curve).ok_or_else(|| unsupported("radial bore: an edge has no bounding box"))?;
            if !apart(bx, cutter_box, tol) {
                return Err(unsupported("radial bore meets another face of the body"));
            }
        }
    }
    Ok((walls, loops))
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    crate::topology::check(body, tol)?;
    crate::topology::check(tool, tol)?;
    let (k, pipe) = admit(body, tool, tol).ok_or_else(|| unsupported("not a radial bore through a pipe wall"))?;
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
    Ok(out)
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
                Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } => unreachable!("no elliptic-cylinder or revolved faces in these test bodies"),
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
    fn holes_meeting_ends_seams_or_each_other_and_wall_pockets_are_refused() {
        let (ro, ri, a) = (1.75, 1.5, 0.25);
        let body = pipe(ro, ri);
        // Ending inside the wall is a pocket, not this domain.
        assert!(!selects(&body, &bore(a, 5., 10., 2.5, 1.6), TOL));
        // Across the pipe's end, and on the seam (the profile circles start at +X: 90 degrees from +Z).
        assert!(matches!(crate::boolean::cut(&body, &bore(a, 0.1, 10., 2.5, 0.5), TOL), Err(KernelError::Unsupported(_))));
        assert!(matches!(crate::boolean::cut(&body, &bore(a, 5., 90., 2.5, 0.5), TOL), Err(KernelError::Unsupported(_))));
        // A second hole overlapping the first.
        let once = crate::boolean::cut(&body, &bore(a, 5., 10., 2.5, 0.5), TOL).unwrap();
        assert!(matches!(crate::boolean::cut(&once, &bore(a, 5.2, 10., 2.5, 0.5), TOL), Err(KernelError::Unsupported(_))));
    }
}

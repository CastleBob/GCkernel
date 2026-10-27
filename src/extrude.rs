//! A flat profile swept into a prism: the kernel's first operation and a self-test of the contract.
//! Exact construction: two plane caps, a plane side face per line and a cylinder side face per arc.

use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face, Profile};
use crate::KernelError;

/// Signed area of a closed loop about `n`, exact for lines and arcs (Green's theorem in the plane frame
/// u, v = n x u; an arc contributes 1/2 [r cx (sin b - sin a) - r cy (cos b - cos a) + r^2 (b - a)]).
pub fn signed_area(lp: &[Curve], n: Vec3) -> f64 {
    let u = n.perpendicular();
    let v = n.cross(u);
    let uv = |p: Vec3| (p.dot(u), p.dot(v));
    let mut area = 0.0;
    for c in lp {
        match c {
            Curve::Line { a, b } => {
                let ((x1, y1), (x2, y2)) = (uv(*a), uv(*b));
                area += 0.5 * (x1 * y2 - x2 * y1);
            }
            Curve::Arc { c: center, n: an, x, r, t0, t1 } => {
                let (cx, cy) = uv(*center);
                let sweep = (t1 - t0) * an.dot(n).signum();
                let start = Curve::point_at_angle(*center, *an, *x, *r, *t0);
                let (sx, sy) = uv(start);
                let a = (sy - cy).atan2(sx - cx);
                let b = a + sweep;
                area += 0.5 * (r * cx * (b.sin() - a.sin()) - r * cy * (b.cos() - a.cos()) + r * r * (b - a));
            }
        }
    }
    area
}

fn closes(lp: &[Curve], tolerance: f64) -> Result<(), KernelError> {
    for k in 0..lp.len() {
        if lp[k].end().distance(lp[(k + 1) % lp.len()].start()) > tolerance {
            return Err(KernelError::InvalidInput(format!("profile loop does not close after curve {k}")));
        }
    }
    Ok(())
}

/// The material removed or added by sweeping the profile: outer area less the holes, times the sweep.
pub fn volume(profile: &Profile) -> f64 {
    let areas: Vec<f64> = profile.loops.iter().map(|l| signed_area(l, profile.n).abs()).collect();
    let outer = areas.iter().cloned().fold(0.0, f64::max);
    (2.0 * outer - areas.iter().sum::<f64>()) * (profile.high - profile.low)
}

pub fn extrude(profile: &Profile, tolerance: f64) -> Result<Body, KernelError> {
    let n = profile.n;
    for lp in &profile.loops {
        closes(lp, tolerance)?;
        for c in lp {
            if let Curve::Arc { n: an, .. } = c {
                if an.cross(n).length() > 1e-9 {
                    return Err(KernelError::InvalidInput("a profile arc does not lie in the profile plane".into()));
                }
            }
        }
    }
    // Outer loop first, counter-clockwise about n; holes clockwise.
    let mut loops: Vec<Vec<Curve>> = profile.loops.clone();
    let outer = (0..loops.len())
        .max_by(|&a, &b| signed_area(&loops[a], n).abs().total_cmp(&signed_area(&loops[b], n).abs()))
        .unwrap();
    loops.swap(0, outer);
    for (i, lp) in loops.iter_mut().enumerate() {
        let ccw = signed_area(lp, n) > 0.0;
        if (i == 0) != ccw {
            *lp = lp.iter().rev().map(Curve::reversed).collect();
        }
    }

    let (down, up) = (n.scale(profile.low), n.scale(profile.high));
    let u = n.perpendicular();
    let mut body = Body::default();
    let mut top_loops = Vec::new();
    let mut bottom_loops = Vec::new();
    let mut sides = Vec::new();
    for lp in &loops {
        let m = lp.len();
        let first_vertex = body.vertices.len();
        for c in lp {
            body.vertices.push(c.start().add(down));
            body.vertices.push(c.start().add(up));
        }
        let bottom_v = |j: usize| first_vertex + 2 * (j % m);
        let top_v = |j: usize| first_vertex + 2 * (j % m) + 1;
        let first_edge = body.edges.len();
        for (j, c) in lp.iter().enumerate() {
            body.edges.push(Edge { curve: c.translated(down), start: bottom_v(j), end: bottom_v(j + 1) });
            body.edges.push(Edge { curve: c.translated(up), start: top_v(j), end: top_v(j + 1) });
            body.edges.push(Edge {
                curve: Curve::Line { a: c.start().add(down), b: c.start().add(up) },
                start: bottom_v(j),
                end: top_v(j),
            });
        }
        let bottom_e = |j: usize| first_edge + 3 * (j % m);
        let top_e = |j: usize| first_edge + 3 * (j % m) + 1;
        let up_e = |j: usize| first_edge + 3 * (j % m) + 2;
        top_loops.push((0..m).map(|j| (top_e(j), true)).collect::<Vec<_>>());
        bottom_loops.push((0..m).rev().map(|j| (bottom_e(j), false)).collect::<Vec<_>>());
        for (j, c) in lp.iter().enumerate() {
            // Outward normal of a side face = tangent x n (the solid lies left of the loop's travel).
            let surface_and_sense = match c {
                Curve::Line { a, b } => {
                    let along = b.sub(*a).unit()?;
                    (Surface::Plane { p: a.add(down), n: along.cross(n).unit()?, x: along }, true)
                }
                Curve::Arc { c: center, n: an, r, .. } => {
                    (Surface::Cylinder { p: center.add(down), a: n, x: u, r: *r }, an.dot(n) > 0.0)
                }
            };
            let lp = vec![(bottom_e(j), true), (up_e(j + 1), true), (top_e(j), false), (up_e(j), false)];
            sides.push(Face { surface: surface_and_sense.0, sense: surface_and_sense.1, loops: vec![lp] });
        }
    }
    let anchor = loops[0][0].start();
    body.faces.push(Face { surface: Surface::Plane { p: anchor.add(up), n, x: u }, sense: true, loops: top_loops });
    body.faces.push(Face { surface: Surface::Plane { p: anchor.add(down), n, x: u }, sense: false, loops: bottom_loops });
    body.faces.extend(sides);
    Ok(body)
}

//! A closed half-section revolved a full turn about an axis (Claude 2026-10-29, Bob: tapered Pole stock): the kernel's
//! second construction. Exact faces: a cylinder for each segment along the axis direction at one radius, a cone for each
//! sloped segment, a plane annulus (or disc) for each segment across the axis; circle edges where they meet, one seam
//! generator per curved face. Segments lying on the axis (a solid section) bound nothing.
//!
//! The section is a polygon of (radius, height) points in the half-plane through the axis containing `x`: radius >= 0,
//! simple (no two segments cross or touch except neighbours at their shared point), and only planar segments may reach
//! the axis (a cone running to its apex is refused, as is a radius within the tolerance of zero but not zero).

use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;

/// The revolve input: axis point `p`, unit axis `a`, unit `x` perpendicular to it (the seam's half-plane), and the
/// section's (radius, height) points in order, closed back to the first.
#[derive(Clone, Debug)]
pub struct Section {
    pub p: Vec3,
    pub a: Vec3,
    pub x: Vec3,
    pub points: Vec<(f64, f64)>,
}

fn invalid<T>(m: impl Into<String>) -> Result<T, KernelError> { Err(KernelError::InvalidInput(m.into())) }

/// Twice the signed area in (radius, height), counter-clockwise positive.
fn area2(pts: &[(f64, f64)]) -> f64 {
    (0..pts.len()).map(|i| { let (a, b) = (pts[i], pts[(i + 1) % pts.len()]); a.0 * b.1 - b.0 * a.1 }).sum()
}

/// Exact volume by Pappus through Green's theorem: pi times the loop integral of r^2 dz, each straight segment
/// contributing dz (r0^2 + r0 r1 + r1^2) / 3 (counter-clockwise in (radius, height) for positive volume).
pub fn volume(section: &Section) -> f64 {
    let pts = &section.points;
    let sum: f64 = (0..pts.len()).map(|i| {
        let ((r0, z0), (r1, z1)) = (pts[i], pts[(i + 1) % pts.len()]);
        (z1 - z0) * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0
    }).sum();
    (std::f64::consts::PI * sum).abs()
}

/// Do segments a0-a1 and b0-b1 cross or touch (closed segments)?
fn meet(a0: (f64, f64), a1: (f64, f64), b0: (f64, f64), b1: (f64, f64), tol: f64) -> bool {
    let cross = |o: (f64, f64), p: (f64, f64), q: (f64, f64)| (p.0 - o.0) * (q.1 - o.1) - (p.1 - o.1) * (q.0 - o.0);
    let dist = |q: (f64, f64), s: (f64, f64), e: (f64, f64)| {
        let (dx, dy) = (e.0 - s.0, e.1 - s.1);
        let l = dx * dx + dy * dy;
        let t = if l > 0.0 { (((q.0 - s.0) * dx + (q.1 - s.1) * dy) / l).clamp(0.0, 1.0) } else { 0.0 };
        (q.0 - s.0 - t * dx).hypot(q.1 - s.1 - t * dy)
    };
    let (d1, d2, d3, d4) = (cross(b0, b1, a0), cross(b0, b1, a1), cross(a0, a1, b0), cross(a0, a1, b1));
    ((d1 > 0.0) != (d2 > 0.0) && (d3 > 0.0) != (d4 > 0.0))
        || dist(a0, b0, b1) <= tol || dist(a1, b0, b1) <= tol || dist(b0, a0, a1) <= tol || dist(b1, a0, a1) <= tol
}

pub fn revolve(section: &Section, tol: f64) -> Result<Body, KernelError> {
    let (p, a, x) = (section.p, section.a, section.x);
    if !p.is_finite() || (a.length() - 1.0).abs() > 1e-9 || (x.length() - 1.0).abs() > 1e-9 || a.dot(x).abs() > 1e-9 {
        return invalid("a revolve needs an axis point, a unit axis and a unit seam direction at right angles");
    }
    let mut pts = section.points.clone();
    let n = pts.len();
    if n < 3 { return invalid("a revolved section needs at least three points"); }
    for (i, &(r, z)) in pts.iter().enumerate() {
        if !(r.is_finite() && z.is_finite()) { return invalid(format!("section point {i} is not finite")); }
        if r < 0.0 || (r > 0.0 && r <= tol) {
            return invalid(format!("section point {i} must lie on the axis (radius 0) or clear of it (radius {r})"));
        }
    }
    for i in 0..n {
        let (s, e) = (pts[i], pts[(i + 1) % n]);
        if (e.0 - s.0).hypot(e.1 - s.1) <= tol { return invalid(format!("section segment {i} has no length")); }
        if (s.0 == 0.0) != (e.0 == 0.0) && (e.1 - s.1).abs() > tol {
            return Err(KernelError::Unsupported(format!("section segment {i} runs to the axis on a slope (a cone apex)")));
        }
        for j in i + 1..n {
            if j == i + 1 || (i == 0 && j == n - 1) { continue; }
            if meet(s, e, pts[j], pts[(j + 1) % n], tol) { return invalid(format!("section segments {i} and {j} cross or touch")); }
        }
    }
    let twice = area2(&pts);
    if twice.abs() <= tol * tol { return invalid("the revolved section encloses no area"); }
    if twice < 0.0 { pts.reverse(); }

    let mut body = Body::default();
    // A vertex and a full circle edge (counter-clockwise about a, starting on the seam) at each point clear of the axis.
    let mut circle = vec![usize::MAX; n];
    let mut vertex = vec![usize::MAX; n];
    for (i, &(r, z)) in pts.iter().enumerate() {
        if r == 0.0 { continue; }
        let centre = p.add(a.scale(z));
        vertex[i] = body.vertices.len();
        body.vertices.push(centre.add(x.scale(r)));
        circle[i] = body.edges.len();
        body.edges.push(Edge {
            curve: Curve::Arc { c: centre, n: a, x, r, t0: 0.0, t1: 2.0 * std::f64::consts::PI },
            start: vertex[i],
            end: vertex[i],
        });
    }
    // Each face uses its first point's circle forward and its second point's reversed: with the section counter-clockwise
    // its outward normal is (dz, -dr) in (radius, height), and the loops then run counter-clockwise about it.
    for i in 0..n {
        let j = (i + 1) % n;
        let ((r0, z0), (r1, z1)) = (pts[i], pts[j]);
        if r0 == 0.0 && r1 == 0.0 { continue; }
        let base = p.add(a.scale(z0));
        if (z1 - z0).abs() <= tol {
            let mut loops = Vec::new();
            let (outer, inner) = if r0 >= r1 { (i, j) } else { (j, i) };
            for k in [outer, inner] {
                if circle[k] != usize::MAX { loops.push(vec![(circle[k], k == i)]); }
            }
            body.faces.push(Face { surface: Surface::Plane { p: base, n: a, x }, sense: r1 < r0, loops });
            continue;
        }
        let seam = body.edges.len();
        body.edges.push(Edge { curve: Curve::Line { a: body.vertices[vertex[i]], b: body.vertices[vertex[j]] }, start: vertex[i], end: vertex[j] });
        // s runs from `base` along a; the face lies between heights z0 and z1.
        let surface = if r1 == r0 {
            Surface::Cylinder { p: base, a, x, r: r0 }
        } else {
            Surface::Cone { p: base, a, x, r: r0, k: (r1 - r0) / (z1 - z0) }
        };
        body.faces.push(Face { surface, sense: z1 > z0, loops: vec![vec![(circle[i], true), (seam, true), (circle[j], false), (seam, false)]] });
    }
    Ok(body)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{mass, tessellate, topology};
    use std::f64::consts::PI;

    const TOL: f64 = 1e-6;

    fn axis(points: Vec<(f64, f64)>) -> Section {
        Section { p: Vec3::new(0.0, 0.0, 0.0), a: Vec3::new(0.0, 1.0, 0.0), x: Vec3::new(1.0, 0.0, 0.0), points }
    }

    /// GCmaker's default Pole (10 OD base 30 long, taper to 4 OD, top 30 long, 90 overall, wall .25), as
    /// SketchFactory_Pole lays out its outer then inner points.
    fn pole(wall: f64) -> Section {
        let (rb, rt) = (5.0, 2.0);
        let mut pts = vec![(rb, 0.0), (rb, 30.0), (rt, 60.0), (rt, 90.0)];
        if wall > 0.0 {
            pts.extend([(rt - wall, 90.0), (rt - wall, 60.0), (rb - wall, 30.0), (rb - wall, 0.0)]);
        } else {
            pts.extend([(0.0, 90.0), (0.0, 0.0)]);
        }
        axis(pts)
    }

    /// Volume of a hollow frustum stack by hand: pi h (R0^2 + R0 R1 + R1^2) / 3 per piece, outer less inner.
    fn frusta(pieces: &[(f64, f64, f64)]) -> f64 { pieces.iter().map(|&(r0, r1, h)| PI * h * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0).sum() }

    #[test]
    fn hollow_pole_is_a_closed_exact_solid_with_the_frustum_volume() {
        let s = pole(0.25);
        let body = revolve(&s, TOL).unwrap();
        let report = topology::check(&body, TOL).unwrap();
        assert_eq!((report.faces, report.edges, report.vertices), (8, 14, 8));
        let expected = frusta(&[(5.0, 5.0, 30.0), (5.0, 2.0, 30.0), (2.0, 2.0, 30.0)])
            - frusta(&[(4.75, 4.75, 30.0), (4.75, 1.75, 30.0), (1.75, 1.75, 30.0)]);
        assert!((volume(&s) - expected).abs() <= 1e-10 * expected);
        let m = mass::measure(&body, TOL).unwrap();
        assert!(m.exact);
        assert!((m.volume - expected).abs() <= 1e-9 * expected, "{} vs {expected}", m.volume);
        // Lateral area: cylinders 2 pi r h, cones pi (r0 + r1) slant; caps the annuli.
        let slant = (9.0f64 + 900.0).sqrt();
        let area = 2.0 * PI * 30.0 * (5.0 + 2.0 + 4.75 + 1.75) + PI * slant * (7.0 + 6.5)
            + PI * (25.0 - 4.75 * 4.75) + PI * (4.0 - 1.75 * 1.75);
        assert!((m.area - area).abs() <= 1e-9 * area, "{} vs {area}", m.area);
    }

    #[test]
    fn solid_pole_and_the_reversed_section_give_the_same_solid() {
        let mut s = pole(0.0);
        let body = revolve(&s, TOL).unwrap();
        topology::check(&body, TOL).unwrap();
        let expected = frusta(&[(5.0, 5.0, 30.0), (5.0, 2.0, 30.0), (2.0, 2.0, 30.0)]);
        assert!((mass::measure(&body, TOL).unwrap().volume - expected).abs() <= 1e-9 * expected);
        s.points.reverse();
        let reversed = revolve(&s, TOL).unwrap();
        assert!((mass::measure(&reversed, TOL).unwrap().volume - expected).abs() <= 1e-9 * expected);
    }

    #[test]
    fn a_tilted_axis_far_away_measures_the_same() {
        let mut s = pole(0.25);
        s.p = Vec3::new(1000.0, -250.0, 40.0);
        s.a = Vec3::new(1.0, 2.0, 2.0).unit().unwrap();
        s.x = Vec3::new(2.0, -1.0, 0.0).unit().unwrap();
        let body = revolve(&s, TOL).unwrap();
        topology::check(&body, TOL).unwrap();
        let v = volume(&s);
        assert!((mass::measure(&body, TOL).unwrap().volume - v).abs() <= 1e-8 * v);
    }

    #[test]
    fn display_mesh_lies_on_the_surfaces_within_the_chord() {
        let body = revolve(&pole(0.25), TOL).unwrap();
        let mesh = tessellate::tessellate(&body, 0.01).unwrap();
        assert!(!mesh.triangles.is_empty());
        // Every display point on the outer taper is on the cone (radius 5 - (h - 30) / 10 between heights 30 and 60).
        for q in &mesh.positions {
            let r = q.x.hypot(q.z);
            if q.y > 30.0 + 1e-9 && q.y < 60.0 - 1e-9 && r > 3.0 && r < 5.0 - 1e-9 {
                let cone = 5.0 - (q.y - 30.0) * 0.1;
                if (r - cone).abs() < 0.2 { assert!((r - cone).abs() <= 1e-9, "{r} vs {cone} at {q:?}"); }
            }
        }
        // Seams are one support each, so not display edges; the eight circles are.
        assert_eq!(mesh.display_edges.len(), 8);
    }

    #[test]
    fn cones_cross_the_wire_and_cuts_on_them_are_refused_by_name() {
        let body = revolve(&pole(0.25), TOL).unwrap();
        let wire = crate::neutral::write_body(&body);
        let mut reader = crate::neutral::Reader::with_tolerance(&wire, TOL);
        let back = reader.body().unwrap();
        reader.finish().unwrap();
        assert_eq!(crate::neutral::write_body(&back), wire);
        assert!(back.faces.iter().filter(|f| matches!(f.surface, Surface::Cone { .. })).count() == 2);
        match crate::boolean::cut(&body, &body, TOL) {
            Err(KernelError::Unsupported(m)) => assert!(m.contains("conical faces"), "{m}"),
            other => panic!("expected a refusal, got {other:?}"),
        }
    }

    #[test]
    fn refuses_apexes_crossings_and_points_just_off_the_axis() {
        let apex = axis(vec![(0.0, 0.0), (2.0, 0.0), (0.0, 5.0)]);
        assert!(matches!(revolve(&apex, TOL), Err(KernelError::Unsupported(_))));
        let bow = axis(vec![(1.0, 0.0), (3.0, 2.0), (3.0, 0.0), (1.0, 2.0)]);
        assert!(matches!(revolve(&bow, TOL), Err(KernelError::InvalidInput(_))));
        let near = axis(vec![(1e-9, 0.0), (2.0, 0.0), (2.0, 1.0), (1e-9, 1.0)]);
        assert!(matches!(revolve(&near, TOL), Err(KernelError::InvalidInput(_))));
        let negative = axis(vec![(-1.0, 0.0), (2.0, 0.0), (2.0, 1.0)]);
        assert!(matches!(revolve(&negative, TOL), Err(KernelError::InvalidInput(_))));
    }
}

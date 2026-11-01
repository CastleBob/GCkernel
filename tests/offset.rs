//! Offsets and shells (Claude 2026-11-01): exact volumes from closed forms worked independently of the construction
//! (Steiner's formula for boxes, Pappus for solids of revolution, strip-and-corner sums for plates), every result
//! checked closed and consistent, measured, tessellated, and refusals by name.

use std::f64::consts::PI;

use gc_kernel::extrude::extrude;
use gc_kernel::geom::{Curve, Vec3};
use gc_kernel::neutral::{write_body, Body, Profile, Reader};
use gc_kernel::offset::{offset, shell, Mode, Offset, Source};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_OFFSET, OP_SHELL};
use gc_kernel::revolve::{revolve_general, Kind, Revolution, Segment};
use gc_kernel::{mass, tessellate, topology, KernelError};

const TOL: f64 = 1e-7;
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn p3(x: f64, y: f64) -> Vec3 { Vec3::new(x, y, 0.0) }
fn line(a: (f64, f64), b: (f64, f64)) -> Curve { Curve::Line { a: p3(a.0, a.1), b: p3(b.0, b.1) } }
/// Counter-clockwise arc about (cx, cy) from angle t0 to t1 (degrees).
fn arc(cx: f64, cy: f64, r: f64, t0: f64, t1: f64) -> Curve {
    Curve::Arc { c: p3(cx, cy), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r, t0: t0.to_radians(), t1: t1.to_radians() }
}
fn circle(cx: f64, cy: f64, r: f64) -> Curve { Curve::Arc { c: p3(cx, cy), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r, t0: 0.0, t1: 2.0 * PI } }
fn rectangle(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Curve> {
    vec![line((x0, y0), (x1, y0)), line((x1, y0), (x1, y1)), line((x1, y1), (x0, y1)), line((x0, y1), (x0, y0))]
}
/// A rectangle with corner radius r, counter-clockwise.
fn rounded_rectangle(x0: f64, y0: f64, x1: f64, y1: f64, r: f64) -> Vec<Curve> {
    vec![
        line((x0 + r, y0), (x1 - r, y0)), arc(x1 - r, y0 + r, r, -90.0, 0.0),
        line((x1, y0 + r), (x1, y1 - r)), arc(x1 - r, y1 - r, r, 0.0, 90.0),
        line((x1 - r, y1), (x0 + r, y1)), arc(x0 + r, y1 - r, r, 90.0, 180.0),
        line((x0, y1 - r), (x0, y0 + r)), arc(x0 + r, y0 + r, r, 180.0, 270.0),
    ]
}
fn prism(loops: Vec<Vec<Curve>>, thickness: f64) -> Body {
    let body = extrude(&Profile { loops, n: Z, low: 0.0, high: thickness }, TOL).unwrap();
    topology::check(&body, TOL).unwrap();
    body
}
fn close(a: f64, b: f64, rel: f64) -> bool { (a - b).abs() <= rel * a.abs().max(b.abs()) }

/// Every result: closed, its measured volume and its reported volume the expected one, displayable.
fn accept(r: &Offset, expected: f64, what: &str) {
    topology::check(&r.body, TOL).unwrap_or_else(|e| panic!("{what}: {e}"));
    let m = mass::measure(&r.body, TOL).unwrap();
    assert!(m.exact, "{what}: measured by quadrature");
    assert!(close(r.volume, expected, 1e-10), "{what}: reported {} vs expected {expected}", r.volume);
    assert!(close(m.volume, expected, 1e-10), "{what}: measured {} vs expected {expected}", m.volume);
    assert_eq!(r.sources.len(), r.body.faces.len());
    let mesh = tessellate::tessellate(&r.body, 0.01).unwrap_or_else(|e| panic!("{what}: tessellate {e}"));
    assert!(!mesh.triangles.is_empty());
}

fn count(r: &Offset) -> (usize, usize, usize) {
    r.sources.iter().fold((0, 0, 0), |(f, e, v), s| match s { Source::Face(_) => (f + 1, e, v), Source::Edge(_) => (f, e + 1, v), Source::Vertex(_) => (f, e, v + 1) })
}

fn unsupported(r: Result<Offset, KernelError>, needle: &str) {
    match r {
        Err(KernelError::Unsupported(m)) => assert!(m.contains(needle), "message {m:?} does not name {needle:?}"),
        other => panic!("expected an Unsupported refusal naming {needle:?}, got {:?}", other.map(|o| o.volume)),
    }
}

// ---- prisms ----

#[test]
fn box_rounded_outward_offset_is_steiners_parallel_body() {
    let (a, b, c, d) = (1.0, 2.0, 3.0, 0.25);
    let body = prism(vec![rectangle(0.0, 0.0, a, b)], c);
    let r = offset(&body, d, Mode::Rounded, TOL).unwrap();
    // Steiner: V + S d + M d^2 + 4/3 pi d^3 (faces, quarter cylinders on the 12 edges, octants at the 8 corners).
    let expected = a * b * c + 2.0 * (a * b + b * c + c * a) * d + PI * (a + b + c) * d * d + 4.0 / 3.0 * PI * d * d * d;
    accept(&r, expected, "box outward");
    assert_eq!(count(&r), (6, 12, 8));
    // Exact area too: faces + quarter cylinders (pi d / 2 per unit edge) + one sphere.
    let area = 2.0 * (a * b + b * c + c * a) + 2.0 * PI * d * (a + b + c) + 4.0 * PI * d * d;
    assert!(close(mass::measure(&r.body, TOL).unwrap().area, area, 1e-10));
}

#[test]
fn box_inward_and_sharp_offsets_are_boxes() {
    let (a, b, c) = (1.0, 2.0, 3.0);
    let body = prism(vec![rectangle(0.0, 0.0, a, b)], c);
    let e = 0.2;
    accept(&offset(&body, -e, Mode::Rounded, TOL).unwrap(), (a - 2.0 * e) * (b - 2.0 * e) * (c - 2.0 * e), "box inward");
    accept(&offset(&body, -e, Mode::Sharp, TOL).unwrap(), (a - 2.0 * e) * (b - 2.0 * e) * (c - 2.0 * e), "box inward sharp");
    let r = offset(&body, e, Mode::Sharp, TOL).unwrap();
    accept(&r, (a + 2.0 * e) * (b + 2.0 * e) * (c + 2.0 * e), "box sharp outward");
    assert_eq!(count(&r), (6, 0, 0));
    unsupported(offset(&body, -0.5, Mode::Rounded, TOL), "");
}

#[test]
fn plate_with_holes_offsets_exactly() {
    let (w, h, t, r0, d) = (6.0, 4.0, 0.5, 0.5, 0.125);
    let body = prism(vec![rectangle(0.0, 0.0, w, h), vec![circle(1.5, 2.0, r0)], vec![circle(4.5, 2.0, r0)]], t);
    let area = w * h - 2.0 * PI * r0 * r0;
    // Outward: the wall (rectangle grown with round corners, holes shrunk), the caps moved, quarter cylinders along the
    // outline, a sphere at the corners, and inner quarter tori in the holes (Pappus: quarter disc pi d^2 / 4 about the
    // axis at r0 - 4 d / (3 pi)).
    let wall = (w + 2.0 * d) * (h + 2.0 * d) - (4.0 - PI) * d * d - 2.0 * PI * (r0 - d).powi(2);
    let expected = t * wall + 2.0 * d * area + 2.0 * (w + h) * PI * d * d / 2.0 + 4.0 / 3.0 * PI * d.powi(3)
        + 2.0 * 2.0 * (PI * d * d / 4.0) * 2.0 * PI * (r0 - 4.0 * d / (3.0 * PI));
    let r = offset(&body, d, Mode::Rounded, TOL).unwrap();
    accept(&r, expected, "plate outward");
    // Inward: a plate of the eroded profile (holes grown, corners sharp), thinner by 2 e.
    let e = 0.1;
    let expected = ((w - 2.0 * e) * (h - 2.0 * e) - 2.0 * PI * (r0 + e).powi(2)) * (t - 2.0 * e);
    accept(&offset(&body, -e, Mode::Rounded, TOL).unwrap(), expected, "plate inward");
    // A hole no larger than the distance would vanish; holes 1.0 apart edge to edge merge under a 0.6 inward offset.
    unsupported(offset(&body, 0.5, Mode::Rounded, TOL), "shrink");
    let near = prism(vec![rectangle(0.0, 0.0, w, h), vec![circle(2.0, 2.0, r0)], vec![circle(3.2, 2.0, r0)]], 2.0);
    unsupported(offset(&near, -0.15, Mode::Rounded, TOL), "touch or cross");
}

/// Exact Steiner-type volume for a tangent-continuous (G1) profile with sharp convex corners of turn `convex[i]` and
/// line-line concave corners of turn `concave[i]` (radians, magnitudes): area(A + rho) = A + P rho + K rho^2 with
/// K = sum(arc sweeps) / 2 + sum(convex) / 2 - sum tan(concave / 2), and per cap band integral A d + P pi d^2 / 4 + 2 K d^3 / 3.
fn rounded_prism_volume(area: f64, perimeter: f64, k: f64, t: f64, d: f64) -> f64 {
    t * (area + perimeter * d + k * d * d) + 2.0 * (area * d + perimeter * PI * d * d / 4.0 + 2.0 * k * d.powi(3) / 3.0)
}

#[test]
fn rectangular_tube_offsets_in_and_out() {
    // 4 x 2 tube, wall 0.25, outer corner radius 0.5, inner 0.25 (tangent arcs throughout), 10 long.
    let (w, h, wall, ro, len) = (4.0, 2.0, 0.25, 0.5, 10.0);
    let ri = ro - wall;
    let outer = rounded_rectangle(0.0, 0.0, w, h, ro);
    let inner = rounded_rectangle(wall, wall, w - wall, h - wall, ri);
    let body = prism(vec![outer, inner], len);
    let rect_area = |a: f64, b: f64, r: f64| a * b - (4.0 - PI) * r * r;
    let rect_perimeter = |a: f64, b: f64, r: f64| 2.0 * (a + b) - 8.0 * r + 2.0 * PI * r;
    let area = rect_area(w, h, ro) - rect_area(w - 2.0 * wall, h - 2.0 * wall, ri);
    let perimeter = rect_perimeter(w, h, ro) + rect_perimeter(w - 2.0 * wall, h - 2.0 * wall, ri);
    // The arcs sweep a full turn each way: K = 0.
    let d = 0.1;
    accept(&offset(&body, d, Mode::Rounded, TOL).unwrap(), rounded_prism_volume(area, perimeter, 0.0, len, d), "tube outward");
    // Inward by 0.05: the outer corners shrink to radius 0.45, the inner grow to 0.3; the wall is 0.15.
    let e = 0.05;
    let eroded = rect_area(w - 2.0 * e, h - 2.0 * e, ro - e) - rect_area(w - 2.0 * wall + 2.0 * e, h - 2.0 * wall + 2.0 * e, ri + e);
    accept(&offset(&body, -e, Mode::Rounded, TOL).unwrap(), eroded * (len - 2.0 * e), "tube inward");
    // Inward past the wall's half thickness closes the wall.
    unsupported(offset(&body, -0.13, Mode::Rounded, TOL), "");
}

#[test]
fn angle_channel_and_beam_profiles_offset_outward_exactly() {
    // Angle 4 x 3 x 0.5 with a 0.375 fillet in the inside corner (tangent) and sharp heel and toes.
    let (a, b, t, rf, len) = (4.0, 3.0, 0.5, 0.375, 8.0);
    let angle = vec![
        line((0.0, 0.0), (a, 0.0)), line((a, 0.0), (a, t)), line((a, t), (t + rf, t)),
        Curve::Arc { c: p3(t + rf, t + rf), n: Z.scale(-1.0), x: Vec3::new(0.0, -1.0, 0.0), r: rf, t0: 0.0, t1: PI / 2.0 },
        line((t, t + rf), (t, b)), line((t, b), (0.0, b)), line((0.0, b), (0.0, 0.0)),
    ];
    let body = prism(vec![angle], len);
    let area = a * t + (b - t) * t + rf * rf - PI * rf * rf / 4.0;
    let perimeter = a + t + (a - t - rf) + PI * rf / 2.0 + (b - t - rf) + t + b;
    // Five convex right-angle corners (heel, two per toe... the heel and four toe corners), the concave fillet sweeps -pi/2.
    let k = 5.0 * (PI / 2.0) / 2.0 - (PI / 2.0) / 2.0;
    let d = 0.125;
    let r = offset(&body, d, Mode::Rounded, TOL).unwrap();
    accept(&r, rounded_prism_volume(area, perimeter, k, len, d), "angle outward");
    // The same angle with a sharp inside corner: the tubes about the two inner edges meet along ellipses.
    let sharp_angle = vec![
        line((0.0, 0.0), (a, 0.0)), line((a, 0.0), (a, t)), line((a, t), (t, t)), line((t, t), (t, b)),
        line((t, b), (0.0, b)), line((0.0, b), (0.0, 0.0)),
    ];
    let body = prism(vec![sharp_angle], len);
    let area = a * t + (b - t) * t;
    let perimeter = 2.0 * (a + b);
    let k = 5.0 * (PI / 2.0) / 2.0 - (PI / 4.0).tan();
    accept(&offset(&body, d, Mode::Rounded, TOL).unwrap(), rounded_prism_volume(area, perimeter, k, len, d), "sharp angle outward");
    accept(&offset(&body, -0.1, Mode::Rounded, TOL).unwrap(), {
        // Inward: legs 0.3 thick, the inside corner rounded with radius 0.1 about the old corner.
        let e = 0.1;
        let (a2, b2, t2) = (a - 2.0 * e, b - 2.0 * e, t - 2.0 * e);
        (a2 * t2 + (b2 - t2) * t2 + e * e - PI * e * e / 4.0) * (len - 2.0 * e)
    }, "sharp angle inward");

    // Channel 6 x 2 x 0.25, sharp corners: a C with two concave corners.
    let (cw, ch, ct) = (6.0, 2.0, 0.25);
    let channel = vec![
        line((0.0, 0.0), (ch, 0.0)), line((ch, 0.0), (ch, ct)), line((ch, ct), (ct, ct)), line((ct, ct), (ct, cw - ct)),
        line((ct, cw - ct), (ch, cw - ct)), line((ch, cw - ct), (ch, cw)), line((ch, cw), (0.0, cw)), line((0.0, cw), (0.0, 0.0)),
    ];
    let body = prism(vec![channel], len);
    let area = cw * ch - (cw - 2.0 * ct) * (ch - ct);
    let perimeter = 2.0 * (cw + ch) + 2.0 * (ch - ct);
    let k = 6.0 * (PI / 2.0) / 2.0 - 2.0 * (PI / 4.0).tan();
    accept(&offset(&body, d, Mode::Rounded, TOL).unwrap(), rounded_prism_volume(area, perimeter, k, len, d), "channel outward");

    // W-shape beam 8 x 6, flanges 0.5, web 0.3, root fillets 0.4 (tangent).
    let (bh, bf, tf, tw, rr) = (8.0, 6.0, 0.5, 0.3, 0.4);
    let (xl, xr) = (bf / 2.0 - tw / 2.0, bf / 2.0 + tw / 2.0);
    let fillet = |cx: f64, cy: f64, from: f64| Curve::Arc { c: p3(cx, cy), n: Z.scale(-1.0), x: Vec3::new(from.to_radians().cos(), from.to_radians().sin(), 0.0), r: rr, t0: 0.0, t1: PI / 2.0 };
    let beam = vec![
        line((0.0, 0.0), (bf, 0.0)), line((bf, 0.0), (bf, tf)), line((bf, tf), (xr + rr, tf)), fillet(xr + rr, tf + rr, -90.0),
        line((xr, tf + rr), (xr, bh - tf - rr)), fillet(xr + rr, bh - tf - rr, 180.0), line((xr + rr, bh - tf), (bf, bh - tf)),
        line((bf, bh - tf), (bf, bh)), line((bf, bh), (0.0, bh)), line((0.0, bh), (0.0, bh - tf)), line((0.0, bh - tf), (xl - rr, bh - tf)),
        fillet(xl - rr, bh - tf - rr, 90.0), line((xl, bh - tf - rr), (xl, tf + rr)), fillet(xl - rr, tf + rr, 0.0),
        line((xl - rr, tf), (0.0, tf)), line((0.0, tf), (0.0, 0.0)),
    ];
    let body = prism(vec![beam], len);
    let area = 2.0 * bf * tf + (bh - 2.0 * tf) * tw + 4.0 * (rr * rr - PI * rr * rr / 4.0);
    let perimeter = 2.0 * bf + 4.0 * tf + 2.0 * (bf - tw - 2.0 * rr) + 2.0 * (bh - 2.0 * tf - 2.0 * rr) + 4.0 * PI * rr / 2.0;
    let k = 8.0 * (PI / 2.0) / 2.0 - 4.0 * (PI / 2.0) / 2.0;
    accept(&offset(&body, d, Mode::Rounded, TOL).unwrap(), rounded_prism_volume(area, perimeter, k, len, d), "beam outward");
    // The fillets (radius 0.4) cannot take a 0.5 outward offset: they would vanish.
    unsupported(offset(&body, 0.5, Mode::Rounded, TOL), "shrink");
}

#[test]
fn elliptic_and_arc_concave_corners_are_refused_by_name() {
    let ellipse = Curve::Conic { c: p3(0.0, 0.0), u: Vec3::new(2.0, 0.0, 0.0), v: Vec3::new(0.0, 1.0, 0.0), t0: 0.0, t1: 2.0 * PI };
    let body = prism(vec![vec![ellipse]], 1.0);
    unsupported(offset(&body, 0.1, Mode::Rounded, TOL), "ellipse");
    // A sharp concave corner between an arc and a line (a half-round boss on a plate's edge): the rounding tubes there
    // meet in a curve that is not a conic.
    let boss = vec![
        line((0.0, 0.0), (4.0, 0.0)), line((4.0, 0.0), (4.0, 2.0)), line((4.0, 2.0), (3.0, 2.0)),
        arc(2.0, 2.0, 1.0, 0.0, 180.0), line((1.0, 2.0), (0.0, 2.0)), line((0.0, 2.0), (0.0, 0.0)),
    ];
    unsupported(offset(&prism(vec![boss], 1.0), 0.1, Mode::Rounded, TOL), "not a conic");
    // A half-round notch cut into the edge has convex corners where the arc meets the edge: rounded exactly (spheres
    // there), area(A + rho) = A + P rho + K rho^2 with six right-angle corners and the notch's sweep of -pi.
    let notch = vec![
        line((0.0, 0.0), (4.0, 0.0)), line((4.0, 0.0), (4.0, 2.0)), line((4.0, 2.0), (3.0, 2.0)),
        Curve::Arc { c: p3(2.0, 2.0), n: Z.scale(-1.0), x: Vec3::new(1.0, 0.0, 0.0), r: 1.0, t0: 0.0, t1: PI },
        line((1.0, 2.0), (0.0, 2.0)), line((0.0, 2.0), (0.0, 0.0)),
    ];
    let body = prism(vec![notch], 1.0);
    let (area, perimeter, k) = (8.0 - PI / 2.0, 10.0 + PI, 6.0 * PI / 4.0 - PI / 2.0);
    accept(&offset(&body, 0.1, Mode::Rounded, TOL).unwrap(), rounded_prism_volume(area, perimeter, k, 1.0, 0.1), "notch outward");
    // Inward (a prism of the planar offset) it is exact: the eroded rectangle [e, 4 - e] x [e, 2 - e] less the notch
    // grown to radius 1 + e about (2, 2), whose part below the eroded top edge is a circular segment.
    let e = 0.1;
    let (rr, h0): (f64, f64) = (1.0 + e, e);
    let segment = rr * rr * (h0 / rr).acos() - h0 * (rr * rr - h0 * h0).sqrt();
    let expected = ((4.0 - 2.0 * e) * (2.0 - 2.0 * e) - segment) * (1.0 - 2.0 * e);
    accept(&offset(&body, -e, Mode::Rounded, TOL).unwrap(), expected, "notch inward");
}

// ---- solids of revolution ----

fn revolved(points: &[(f64, f64, Kind)]) -> Body {
    let rev = Revolution { p: Vec3::new(0.0, 0.0, 0.0), a: Z, x: Vec3::new(1.0, 0.0, 0.0), angle: 2.0 * PI,
        loops: vec![points.iter().map(|&(r, z, kind)| Segment { r, z, kind }).collect()] };
    let body = revolve_general(&rev, TOL).unwrap();
    topology::check(&body, TOL).unwrap();
    body
}
fn cylinder(r: f64, h: f64) -> Body { revolved(&[(0.0, 0.0, Kind::Line), (r, 0.0, Kind::Line), (r, h, Kind::Line), (0.0, h, Kind::Line)]) }

#[test]
fn cylinder_offsets_exactly() {
    let (r, h) = (2.0, 5.0);
    let body = cylinder(r, h);
    let d = 0.25;
    // Pappus: the grown cylinder, the two discs moved out, and two quarter tori about the rim circles.
    let expected = PI * (r + d).powi(2) * h + 2.0 * PI * r * r * d + 2.0 * (PI * d * d / 4.0) * 2.0 * PI * (r + 4.0 * d / (3.0 * PI));
    let out = offset(&body, d, Mode::Rounded, TOL).unwrap();
    accept(&out, expected, "cylinder outward");
    assert!(out.sources.iter().any(|s| matches!(s, Source::Edge(_))), "the rim tori come from the rim edges");
    accept(&offset(&body, -d, Mode::Rounded, TOL).unwrap(), PI * (r - d).powi(2) * (h - 2.0 * d), "cylinder inward");
    accept(&offset(&body, d, Mode::Sharp, TOL).unwrap(), PI * (r + d).powi(2) * (h + 2.0 * d), "cylinder sharp");
    unsupported(offset(&body, -2.0, Mode::Rounded, TOL), "");
}

#[test]
fn sphere_torus_and_cone_offsets_exactly() {
    let sphere = revolved(&[(0.0, -2.0, Kind::Arc { cr: 0.0, cz: 0.0, ccw: true }), (0.0, 2.0, Kind::Line)]);
    accept(&offset(&sphere, 0.5, Mode::Rounded, TOL).unwrap(), 4.0 / 3.0 * PI * 2.5f64.powi(3), "sphere outward");
    accept(&offset(&sphere, -0.5, Mode::Rounded, TOL).unwrap(), 4.0 / 3.0 * PI * 1.5f64.powi(3), "sphere inward");
    let (big, tube) = (3.0, 1.0);
    let torus = revolved(&[(big + tube, 0.0, Kind::Arc { cr: big, cz: 0.0, ccw: true }), (big - tube, 0.0, Kind::Arc { cr: big, cz: 0.0, ccw: true })]);
    accept(&offset(&torus, 0.25, Mode::Rounded, TOL).unwrap(), 2.0 * PI * PI * big * 1.25f64.powi(2), "torus outward");
    accept(&offset(&torus, -0.25, Mode::Rounded, TOL).unwrap(), 2.0 * PI * PI * big * 0.75f64.powi(2), "torus inward");
    unsupported(offset(&torus, -1.0, Mode::Rounded, TOL), "shrink");
    // Cone, base radius 3 at z = 0, apex at height 4 (slant 5): offset by d, by Pappus over the offset section.
    let cone = revolved(&[(0.0, 0.0, Kind::Line), (3.0, 0.0, Kind::Line), (0.0, 4.0, Kind::Line)]);
    let d = 0.5;
    // The slant's outward normal is (4, 3) / 5; the section: base line z = -d, an arc about (3, 0) from -90 degrees to
    // the normal's angle, the moved slant, an arc about the apex from the normal's angle to the axis.
    let normal = (0.8f64, 0.6f64);
    let phi = normal.1.atan2(normal.0);
    let pappus_line = |a: (f64, f64), b: (f64, f64)| PI * (b.1 - a.1) * (a.0 * a.0 + a.0 * b.0 + b.0 * b.0) / 3.0;
    // pi times the integral of r^2 dz along an arc about (cr, cz), radius rho, from angle t0 to t1.
    let pappus_arc = |cr: f64, cz: f64, rho: f64, t0: f64, t1: f64| {
        let _ = cz;
        let f = |t: f64| rho * (cr * cr * t.sin() + 2.0 * cr * rho * (t / 2.0 + (2.0 * t).sin() / 4.0) + rho * rho * (t.sin() - t.sin().powi(3) / 3.0));
        PI * (f(t1) - f(t0))
    };
    let (s0, s1) = ((3.0 + d * normal.0, d * normal.1), (d * normal.0, 4.0 + d * normal.1));
    // The loop runs counter-clockwise: up the axis is not part of it (r = 0 contributes nothing).
    let expected = pappus_line((0.0, -d), (3.0, -d)) + pappus_arc(3.0, 0.0, d, -PI / 2.0, phi) + pappus_line(s0, s1)
        + pappus_arc(0.0, 4.0, d, phi, PI / 2.0);
    accept(&offset(&cone, d, Mode::Rounded, TOL).unwrap(), expected, "cone outward");
    // Inward: every face moved in by e; the base up to z = e, the slant in to meet it and the axis (a smaller cone).
    let e = 0.25;
    // The moved slant: 4 r + 3 z = 12 - 5 e; at z = e, r = (12 - 8 e) / 4; at r = 0, z = (12 - 5 e) / 3.
    let (rb, za) = ((12.0 - 8.0 * e) / 4.0, (12.0 - 5.0 * e) / 3.0);
    accept(&offset(&cone, -e, Mode::Rounded, TOL).unwrap(), PI * rb * rb * (za - e) / 3.0, "cone inward");
}

#[test]
fn revolved_pole_with_a_bore_offsets_both_walls() {
    // A tube (bore 1, outside 2, 6 long): its section is clear of the axis.
    let tube = revolved(&[(1.0, 0.0, Kind::Line), (2.0, 0.0, Kind::Line), (2.0, 6.0, Kind::Line), (1.0, 6.0, Kind::Line)]);
    let d = 0.2;
    // Outward: the bore shrinks, the outside grows, the ends move; outer rims get outer quarter tori, bore rims inner.
    let quarter = PI * d * d / 4.0;
    let c = 4.0 * d / (3.0 * PI);
    let expected = PI * ((2.0 + d).powi(2) - (1.0 - d).powi(2)) * 6.0 + 2.0 * PI * (4.0 - 1.0) * d
        + 2.0 * quarter * 2.0 * PI * (2.0 + c) + 2.0 * quarter * 2.0 * PI * (1.0 - c);
    accept(&offset(&tube, d, Mode::Rounded, TOL).unwrap(), expected, "tube outward");
    // Past the bore's radius the bore would close: refused.
    unsupported(offset(&tube, 1.0, Mode::Rounded, TOL), "");
}

/// pi times the loop integral of r^2 dz along a line (Pappus).
fn pappus_line(a: (f64, f64), b: (f64, f64)) -> f64 { PI * (b.1 - a.1) * (a.0 * a.0 + a.0 * b.0 + b.0 * b.0) / 3.0 }
/// The same along an arc about (cr, .) of radius rho from angle t0 to t1 (either way round).
fn pappus_arc(cr: f64, rho: f64, t0: f64, t1: f64) -> f64 {
    let f = |t: f64| rho * (cr * cr * t.sin() + 2.0 * cr * rho * (t / 2.0 + (2.0 * t).sin() / 4.0) + rho * rho * (t.sin() - t.sin().powi(3) / 3.0));
    PI * (f(t1) - f(t0))
}

#[test]
fn filleted_and_grooved_cylinders_offset_their_tori_exactly() {
    let h = PI / 2.0;
    // Radius 2, height 3, the top rim filleted with radius 0.5 (a torus face).
    let filleted = revolved(&[(0.0, 0.0, Kind::Line), (2.0, 0.0, Kind::Line), (2.0, 2.5, Kind::Arc { cr: 1.5, cz: 2.5, ccw: true }),
        (1.5, 3.0, Kind::Line), (0.0, 3.0, Kind::Line)]);
    let d = 0.25;
    // Out: a torus about the bottom rim, the side grown, the fillet's tube grown to 0.75, the top moved up.
    let expected = pappus_arc(2.0, d, -h, 0.0) + pappus_line((2.0 + d, 0.0), (2.0 + d, 2.5)) + pappus_arc(1.5, 0.5 + d, 0.0, h);
    accept(&offset(&filleted, d, Mode::Rounded, TOL).unwrap(), expected, "filleted cylinder outward");
    // In: the bottom corner sharp, the fillet's tube shrunk to 0.3.
    let e = 0.2;
    let expected = pappus_line((2.0 - e, e), (2.0 - e, 2.5)) + pappus_arc(1.5, 0.5 - e, 0.0, h);
    accept(&offset(&filleted, -e, Mode::Rounded, TOL).unwrap(), expected, "filleted cylinder inward");
    unsupported(offset(&filleted, -0.5, Mode::Rounded, TOL), "shrink");

    // Radius 2, height 4, a half-round groove of radius 0.5 about (2, 2): its edges are convex right-angle corners.
    let grooved = revolved(&[(0.0, 0.0, Kind::Line), (2.0, 0.0, Kind::Line), (2.0, 1.5, Kind::Arc { cr: 2.0, cz: 2.0, ccw: false }),
        (2.0, 2.5, Kind::Line), (2.0, 4.0, Kind::Line), (0.0, 4.0, Kind::Line)]);
    let d = 0.2;
    let expected = pappus_arc(2.0, d, -h, 0.0) + pappus_line((2.0 + d, 0.0), (2.0 + d, 1.5)) + pappus_arc(2.0, d, 0.0, h)
        + pappus_arc(2.0, 0.5 - d, -h, -3.0 * h) + pappus_arc(2.0, d, -h, 0.0) + pappus_line((2.0 + d, 2.5), (2.0 + d, 4.0))
        + pappus_arc(2.0, d, 0.0, h);
    let r = offset(&grooved, d, Mode::Rounded, TOL).unwrap();
    accept(&r, expected, "grooved cylinder outward");
    // Rounding tori about the four rim and groove-edge circles, from those edges.
    assert_eq!(count(&r).1, 4);
    unsupported(offset(&grooved, 0.5, Mode::Rounded, TOL), "shrink");
}

#[test]
fn box_offset_sources_name_the_entities_each_face_rolls_about() {
    use gc_kernel::geom::Surface;
    let body = prism(vec![rectangle(0.0, 0.0, 1.0, 2.0)], 3.0);
    let d = 0.25;
    let r = offset(&body, d, Mode::Rounded, TOL).unwrap();
    for (f, s) in r.sources.iter().enumerate() {
        match (*s, &r.body.faces[f].surface) {
            (Source::Face(i), Surface::Plane { p, n, .. }) => {
                let Surface::Plane { p: q, n: m, .. } = body.faces[i].surface else { panic!("face {i} is not a plane") };
                let out = if body.faces[i].sense { m } else { m.scale(-1.0) };
                assert!(n.cross(m).length() < 1e-12 && (p.sub(q).dot(out) - d).abs() < 1e-12, "face {f} is not face {i} moved out by d");
            }
            (Source::Edge(e), Surface::Cylinder { p, a, r: radius, .. }) => {
                let gc_kernel::geom::Curve::Line { a: s, b: t } = body.edges[e].curve else { panic!("edge {e} is not a line") };
                let w = s.sub(*p);
                assert!((radius - d).abs() < 1e-12 && t.sub(s).cross(*a).length() < 1e-12 && w.sub(a.scale(w.dot(*a))).length() < 1e-12,
                    "face {f} is not the tube about edge {e}");
            }
            (Source::Vertex(v), Surface::Sphere { c, r: radius, .. }) => {
                assert!((radius - d).abs() < 1e-12 && c.distance(body.vertices[v]) < 1e-12, "face {f} is not the ball about vertex {v}");
            }
            (s, surface) => panic!("face {f} from {s:?} is a {surface:?}"),
        }
    }
}

// ---- shells ----

#[test]
fn box_shells_open_and_closed() {
    let (a, b, c, t) = (3.0, 2.0, 1.5, 0.2);
    let body = prism(vec![rectangle(0.0, 0.0, a, b)], c);
    let top = (0..body.faces.len()).find(|&f| matches!(body.faces[f].surface, gc_kernel::geom::Surface::Plane { n, .. }
        if (if body.faces[f].sense { n } else { n.scale(-1.0) }).dot(Z) > 0.5)).unwrap();
    let open = shell(&body, t, &[top], TOL).unwrap();
    accept(&open, a * b * c - (a - 2.0 * t) * (b - 2.0 * t) * (c - t), "box open at the top");
    // Closed: the cavity is a void inside the same body.
    let closed = shell(&body, t, &[], TOL).unwrap();
    accept(&closed, a * b * c - (a - 2.0 * t) * (b - 2.0 * t) * (c - 2.0 * t), "box closed shell");
    assert_eq!(closed.body.faces.len(), 12);
    // Removing a side face: the box is a prism along that side's normal too, the cavity running through that wall.
    let side = (0..body.faces.len()).find(|&f| f != top && matches!(body.faces[f].surface, gc_kernel::geom::Surface::Plane { n, .. } if n.dot(Z).abs() < 0.5)).unwrap();
    let gc_kernel::geom::Surface::Plane { n, .. } = body.faces[side].surface else { unreachable!() };
    let (along, across) = if n.x.abs() > 0.5 { (a, b) } else { (b, a) };
    accept(&shell(&body, t, &[side], TOL).unwrap(), a * b * c - (across - 2.0 * t) * (c - 2.0 * t) * (along - t), "box open at a side");
    // Outward: the wall grows outside the box, with rounded bottom edges and corners; the top stays at the box's top.
    let w = 0.25;
    let outer = c * (a * b + 2.0 * (a + b) * w + PI * w * w) + (a * b * w + 2.0 * (a + b) * PI * w * w / 4.0 + 2.0 * PI * w.powi(3) / 3.0);
    accept(&shell(&body, -w, &[top], TOL).unwrap(), outer - a * b * c, "box open, wall outward");
    // Two adjacent faces are not both caps of one prism.
    unsupported(shell(&body, t, &[top, side], TOL), "caps");
}

#[test]
fn cylinder_cup_and_sphere_shells() {
    let (r, h, t) = (2.0, 5.0, 0.25);
    let body = cylinder(r, h);
    let top = (0..body.faces.len()).find(|&f| matches!(body.faces[f].surface, gc_kernel::geom::Surface::Plane { p, .. } if p.z > h - 1e-9)).unwrap();
    let cup = shell(&body, t, &[top], TOL).unwrap();
    accept(&cup, PI * r * r * h - PI * (r - t).powi(2) * (h - t), "cylinder cup");
    // Both ends removed: a tube.
    let bottom = (0..body.faces.len()).find(|&f| matches!(body.faces[f].surface, gc_kernel::geom::Surface::Plane { p, .. } if p.z < 1e-9)).unwrap();
    accept(&shell(&body, t, &[top, bottom], TOL).unwrap(), PI * (r * r - (r - t).powi(2)) * h, "cylinder tube");
    // Outward cup: the wall outside, the bottom rim rounded (a quarter torus), open top flush.
    let w = 0.25;
    let outer = PI * (r + w).powi(2) * h + PI * r * r * w + (PI * w * w / 4.0) * 2.0 * PI * (r + 4.0 * w / (3.0 * PI));
    accept(&shell(&body, -w, &[top], TOL).unwrap(), outer - PI * r * r * h, "cylinder cup outward");
    // A closed spherical shell: outer sphere, inner void sphere.
    let sphere = revolved(&[(0.0, -2.0, Kind::Arc { cr: 0.0, cz: 0.0, ccw: true }), (0.0, 2.0, Kind::Line)]);
    accept(&shell(&sphere, 0.3, &[], TOL).unwrap(), 4.0 / 3.0 * PI * (8.0 - 1.7f64.powi(3)), "hollow sphere");
    // The sphere's surface is split in two halves: removing only one is refused.
    let half = (0..sphere.faces.len()).next().unwrap();
    unsupported(shell(&sphere, 0.3, &[half], TOL), "remove all of them or none");
}

#[test]
fn plate_with_holes_shell_open_at_the_top() {
    let (w, h, th, r0, t) = (6.0, 4.0, 1.0, 0.5, 0.2);
    let body = prism(vec![rectangle(0.0, 0.0, w, h), vec![circle(1.5, 2.0, r0)], vec![circle(4.5, 2.0, r0)]], th);
    let area = w * h - 2.0 * PI * r0 * r0;
    let eroded = (w - 2.0 * t) * (h - 2.0 * t) - 2.0 * PI * (r0 + t).powi(2);
    let top = 0; // extrude's high cap
    let r = shell(&body, t, &[top], TOL).unwrap();
    accept(&r, area * th - eroded * (th - t), "plate shell");
    // Removing a hole's wall is not a cap of the plate.
    let hole = (0..body.faces.len()).find(|&f| matches!(body.faces[f].surface, gc_kernel::geom::Surface::Cylinder { .. })).unwrap();
    unsupported(shell(&body, t, &[hole], TOL), "caps");
}

// ---- the wire ----

fn reply_values(op: u32, data: &[f64]) -> Vec<f64> {
    let req = parse_request(&encode_request(op, TOL, data)).unwrap();
    handle(&req).unwrap()
}

#[test]
fn ops_20_and_21_carry_body_volume_and_sources() {
    let (a, b, c) = (1.0, 2.0, 3.0);
    let body = prism(vec![rectangle(0.0, 0.0, a, b)], c);
    let mut data = vec![0.25, 0.0];
    data.extend(write_body(&body));
    let out = reply_values(OP_OFFSET, &data);
    let mut reader = Reader::new(&out);
    let result = reader.body().unwrap();
    let volume = reader.next().unwrap();
    let faces = reader.count().unwrap();
    assert_eq!(faces, result.faces.len());
    let mut kinds = [0; 3];
    for _ in 0..faces {
        let (kind, index) = (reader.next().unwrap(), reader.next().unwrap());
        kinds[kind as usize] += 1;
        let limit = [body.faces.len(), body.edges.len(), body.vertices.len()][kind as usize];
        assert!(index >= 0.0 && (index as usize) < limit);
    }
    reader.finish().unwrap();
    assert_eq!(kinds, [6, 12, 8]);
    let d = 0.25;
    assert!(close(volume, a * b * c + 2.0 * (a * b + b * c + c * a) * d + PI * (a + b + c) * d * d + 4.0 / 3.0 * PI * d.powi(3), 1e-10));

    let mut data = vec![0.1, 1.0, 0.0];
    data.extend(write_body(&body));
    let out = reply_values(OP_SHELL, &data);
    let mut reader = Reader::new(&out);
    let result = reader.body().unwrap();
    let volume = reader.next().unwrap();
    assert_eq!(reader.count().unwrap(), result.faces.len());
    assert!(close(volume, a * b * c - 0.8 * 1.8 * 2.9, 1e-10));

    // Refusals travel as status 1 with the reason.
    let mut data = vec![5.0, 0.0];
    data.extend(write_body(&prism(vec![rectangle(0.0, 0.0, 1.0, 1.0), vec![circle(0.5, 0.5, 0.25)]], 1.0)));
    let req = parse_request(&encode_request(OP_OFFSET, TOL, &data)).unwrap();
    assert!(matches!(handle(&req), Err(KernelError::Unsupported(_))));
    let mut data = vec![0.1, 2.0];
    data.extend(write_body(&body));
    let req = parse_request(&encode_request(OP_OFFSET, TOL, &data)).unwrap();
    assert!(matches!(handle(&req), Err(KernelError::InvalidInput(_))));
}

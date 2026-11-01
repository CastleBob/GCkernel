//! Fillet and chamfer acceptance (Claude 2026-11-01). Every expected volume is closed form, from the decomposition
//! named at each test; every result must pass topology::check, measure exactly, round-trip through the neutral format
//! and tessellate into a closed, outward mesh whose volume agrees within the chord tolerance.

use std::collections::HashMap;
use std::f64::consts::PI;

use gc_kernel::blend::{blend, Item, Kind, Origin};
use gc_kernel::extrude::extrude;
use gc_kernel::geom::{Curve, Surface, Vec3};
use gc_kernel::neutral::{write_body, Body, Profile, Reader};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_BLEND};
use gc_kernel::revolve::{revolve, Section};
use gc_kernel::{mass, tessellate, topology, KernelError};

const TOL: f64 = 1e-7;
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn p2(x: f64, y: f64) -> Vec3 { Vec3::new(x, y, 0.0) }
fn line(a: (f64, f64), b: (f64, f64)) -> Curve { Curve::Line { a: p2(a.0, a.1), b: p2(b.0, b.1) } }
fn polygon(points: &[(f64, f64)]) -> Vec<Curve> {
    (0..points.len()).map(|i| line(points[i], points[(i + 1) % points.len()])).collect()
}
fn circle(cx: f64, cy: f64, r: f64) -> Curve {
    Curve::Arc { c: p2(cx, cy), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r, t0: 0.0, t1: 2.0 * PI }
}
fn prism(loops: Vec<Vec<Curve>>, high: f64) -> Body {
    extrude(&Profile { loops, n: Z, low: 0.0, high }, TOL).unwrap()
}
fn block(a: f64, b: f64, c: f64) -> Body { prism(vec![polygon(&[(0., 0.), (a, 0.), (a, b), (0., b)])], c) }

/// The edge whose curve passes through `p` at its parameter midpoint.
fn edge_at(body: &Body, p: Vec3) -> usize {
    let found: Vec<usize> = (0..body.edges.len()).filter(|&e| {
        let (t0, t1) = body.edges[e].curve.domain();
        body.edges[e].curve.eval((t0 + t1) / 2.0).0.distance(p) < 1e-9
    }).collect();
    assert_eq!(found.len(), 1, "edges at {p:?}: {found:?}");
    found[0]
}
/// The circle edge with centre c (a full circle's midpoint is opposite its start).
fn circle_at(body: &Body, c: Vec3, r: f64) -> usize {
    let found: Vec<usize> = (0..body.edges.len()).filter(|&e| match body.edges[e].curve {
        Curve::Arc { c: ac, r: ar, .. } => ac.distance(c) < 1e-9 && (ar - r).abs() < 1e-9,
        _ => false,
    }).collect();
    assert_eq!(found.len(), 1, "circles at {c:?}: {found:?}");
    found[0]
}
fn plane_face(body: &Body, outward: Vec3, through: Vec3) -> usize {
    (0..body.faces.len()).find(|&f| match body.faces[f].surface {
        Surface::Plane { p, n, .. } => {
            let n = if body.faces[f].sense { n } else { n.scale(-1.0) };
            n.sub(outward).length() < 1e-9 && through.sub(p).dot(n).abs() < 1e-9
        }
        _ => false,
    }).expect("the face")
}

fn close(got: f64, want: f64) { assert!((got - want).abs() <= 1e-9 * want.abs().max(1.0), "volume {got} vs closed form {want} (off {})", got - want); }

/// Blends, then proves the result: topology, exact measure against `want`, neutral round trip, closed outward mesh.
fn blended(body: &Body, items: &[Item], want: f64) -> Body {
    let result = blend(body, items, TOL).unwrap_or_else(|e| panic!("blend refused: {e}"));
    let out = result.body;
    assert_eq!(result.origins.len(), out.faces.len());
    topology::check(&out, TOL).unwrap();
    let m = mass::measure(&out, TOL).unwrap();
    assert!(m.exact, "the measure should be closed form");
    close(m.volume, want);
    let wire = write_body(&out);
    let mut reader = Reader::with_tolerance(&wire, TOL);
    let back = reader.body().unwrap();
    reader.finish().unwrap();
    topology::check(&back, TOL).unwrap();
    mesh_check(&back, 0.002, m.volume, m.area);
    out
}

fn mesh_check(body: &Body, chord: f64, volume: f64, area: f64) {
    let mesh = tessellate::tessellate(body, chord).unwrap();
    let bits = |x: f64| if x == 0.0 { 0 } else { x.to_bits() };
    let key = |p: Vec3| (bits(p.x), bits(p.y), bits(p.z));
    let mut boundary = HashMap::<_, i32>::new();
    let mut v = 0.0;
    let o = mesh.positions[0];
    for (tri, &face) in mesh.triangles.iter().zip(&mesh.faces) {
        let [a, b, c] = tri.map(|i| mesh.positions[i]);
        for (p, q) in [(a, b), (b, c), (c, a)] {
            if key(p) == key(q) { continue; }
            *boundary.entry((key(p), key(q))).or_default() += 1;
            *boundary.entry((key(q), key(p))).or_default() -= 1;
        }
        let cross = b.sub(a).cross(c.sub(a));
        if cross.length() > 1e-14 {
            let normal = tri.iter().fold(Vec3::new(0.0, 0.0, 0.0), |s, &i| s.add(mesh.normals[i]));
            assert!(cross.dot(normal) > 0.0, "face {face}: inward triangle");
        }
        v += a.sub(o).dot(b.sub(o).cross(c.sub(o))) / 6.0;
    }
    assert!(boundary.values().all(|&n| n == 0), "the mesh is not closed");
    assert!((v - volume).abs() <= chord * area, "mesh volume {v} vs {volume}");
}

fn refused(body: &Body, items: &[Item]) -> String {
    match blend(body, items, TOL) {
        Err(KernelError::Unsupported(m)) => m,
        Err(e) => panic!("expected a named refusal, got {e}"),
        Ok(_) => panic!("expected a refusal"),
    }
}

// The spandrel a fillet of radius r removes between two faces at right angles: area (1 - pi/4) r^2, centroid at
// r (5/6 - pi/4) / (1 - pi/4) from each face.
fn spandrel(r: f64) -> f64 { (1.0 - PI / 4.0) * r * r }
fn spandrel_moment(r: f64) -> f64 { r * r * r * (5.0 / 6.0 - PI / 4.0) }

// ---- 1. chamfers between flat faces ----

#[test]
fn chamfer_one_box_edge_equal_distances() {
    let (a, b, c, d) = (4.0, 3.0, 2.0, 0.5);
    let body = block(a, b, c);
    let e = edge_at(&body, Vec3::new(a / 2.0, 0.0, c));
    // Box less a triangular prism d^2 / 2 by a.
    let out = blended(&body, &[Item { edge: e, kind: Kind::Chamfer { distance: d } }], a * b * c - d * d / 2.0 * a);
    assert_eq!(out.faces.len(), 7);
}

#[test]
fn chamfer_two_distances_runs_the_first_on_the_named_face() {
    let (a, b, c) = (4.0, 3.0, 2.0);
    let body = block(a, b, c);
    let e = edge_at(&body, Vec3::new(a / 2.0, 0.0, c));
    let top = plane_face(&body, Z, Vec3::new(0.0, 0.0, c));
    let out = blended(&body, &[Item { edge: e, kind: Kind::ChamferDistances { first: 0.8, second: 0.3, face: top } }],
        a * b * c - 0.8 * 0.3 / 2.0 * a);
    assert!(out.vertices.iter().any(|v| v.distance(Vec3::new(0.0, 0.8, c)) < 1e-12), "0.8 runs on the top face");
    assert!(out.vertices.iter().any(|v| v.distance(Vec3::new(0.0, 0.0, c - 0.3)) < 1e-12), "0.3 runs down the side");
}

#[test]
fn chamfer_distance_and_angle() {
    let (a, b, c) = (4.0, 3.0, 2.0);
    let body = block(a, b, c);
    let e = edge_at(&body, Vec3::new(a / 2.0, 0.0, c));
    let top = plane_face(&body, Z, Vec3::new(0.0, 0.0, c));
    let (d, angle) = (0.5, PI / 6.0);
    // At right angles the other leg is d tan(angle).
    blended(&body, &[Item { edge: e, kind: Kind::ChamferAngle { distance: d, angle, face: top } }],
        a * b * c - d * d * angle.tan() / 2.0 * a);
}

fn l_block() -> Body { prism(vec![polygon(&[(0., 0.), (4., 0.), (4., 1.), (1., 1.), (1., 3.), (0., 3.)])], 2.0) }
const L_VOLUME: f64 = (4.0 + 2.0) * 2.0;

#[test]
fn chamfer_concave_edge_adds_material() {
    let body = l_block();
    let e = edge_at(&body, Vec3::new(1.0, 1.0, 1.0));
    blended(&body, &[Item { edge: e, kind: Kind::Chamfer { distance: 0.25 } }], L_VOLUME + 0.25 * 0.25 / 2.0 * 2.0);
}

fn top_loop(body: &Body, a: f64, b: f64, c: f64) -> Vec<usize> {
    [(a / 2.0, 0.0), (a, b / 2.0), (a / 2.0, b), (0.0, b / 2.0)].iter().map(|&(x, y)| edge_at(body, Vec3::new(x, y, c))).collect()
}

#[test]
fn chamfer_top_loop_mitres_at_the_corners() {
    let (a, b, c, d) = (4.0, 3.0, 2.0, 0.4);
    let body = block(a, b, c);
    let items: Vec<Item> = top_loop(&body, a, b, c).into_iter().map(|edge| Item { edge, kind: Kind::Chamfer { distance: d } }).collect();
    // Prismatoid: the top slab a b d keeps h/6 (ab + 4 (a-d)(b-d) + (a-2d)(b-2d)).
    let kept = d / 6.0 * (a * b + 4.0 * (a - d) * (b - d) + (a - 2.0 * d) * (b - 2.0 * d));
    blended(&body, &items, a * b * (c - d) + kept);
}

#[test]
fn chamfer_every_box_edge_gives_triangular_corners() {
    let (a, b, c, d) = (4.0, 3.0, 2.0, 0.3);
    let body = block(a, b, c);
    let items: Vec<Item> = (0..body.edges.len()).map(|edge| Item { edge, kind: Kind::Chamfer { distance: d } }).collect();
    let out = blended(&body, &items, {
        // Core, six face slabs d thick, twelve edge prisms keeping d^2 / 2, eight corner cubes keeping 1/6 (the corner
        // plane through the three face corners is x + y + z = 2d inside the d-cube; the edge planes are implied there).
        let (x, y, z) = (a - 2.0 * d, b - 2.0 * d, c - 2.0 * d);
        x * y * z + 2.0 * d * (x * y + x * z + y * z) + 4.0 * (x + y + z) * d * d / 2.0 + 8.0 * d * d * d / 6.0
    });
    assert_eq!(out.faces.len(), 6 + 12 + 8);
    let origins = blend(&body, &items, TOL).unwrap().origins;
    assert_eq!(origins.iter().filter(|o| matches!(o, Origin::Vertex(_))).count(), 8);
    assert_eq!(origins.iter().filter(|o| matches!(o, Origin::Edge(_))).count(), 12);
    assert_eq!(&origins[..6], &(0..6).map(Origin::Face).collect::<Vec<_>>()[..]);
}

// ---- 2. fillets between flat faces ----

#[test]
fn fillet_one_box_edge() {
    let (a, b, c, r) = (4.0, 3.0, 2.0, 0.5);
    let body = block(a, b, c);
    let e = edge_at(&body, Vec3::new(a / 2.0, 0.0, c));
    let out = blended(&body, &[Item { edge: e, kind: Kind::Fillet { radius: r } }], a * b * c - spandrel(r) * a);
    assert!(out.faces.iter().any(|f| matches!(f.surface, Surface::Cylinder { r: q, .. } if q == r)));
}

#[test]
fn fillet_concave_edge_adds_material() {
    let body = l_block();
    let e = edge_at(&body, Vec3::new(1.0, 1.0, 1.0));
    blended(&body, &[Item { edge: e, kind: Kind::Fillet { radius: 0.3 } }], L_VOLUME + spandrel(0.3) * 2.0);
}

#[test]
fn fillet_top_loop_mitres_in_ellipses() {
    let (a, b, c, r) = (4.0, 3.0, 2.0, 0.4);
    let body = block(a, b, c);
    let items: Vec<Item> = top_loop(&body, a, b, c).into_iter().map(|edge| Item { edge, kind: Kind::Fillet { radius: r } }).collect();
    // Each side's spandrel prism is cut by 45 degree mitres: its volume is the area times (side - 2 x centroid offset).
    let out = blended(&body, &items, a * b * c - (spandrel(r) * 2.0 * (a + b) - 8.0 * spandrel_moment(r)));
    assert_eq!(out.edges.iter().filter(|e| matches!(e.curve, Curve::Conic { .. })).count(), 4, "four mitre ellipses");
}

#[test]
fn fillet_every_box_edge_rounds_corners_with_sphere_octants() {
    let (a, b, c, r) = (4.0, 3.0, 2.0, 0.35);
    let body = block(a, b, c);
    let items: Vec<Item> = (0..body.edges.len()).map(|edge| Item { edge, kind: Kind::Fillet { radius: r } }).collect();
    let (x, y, z) = (a - 2.0 * r, b - 2.0 * r, c - 2.0 * r);
    let want = x * y * z + 2.0 * r * (x * y + x * z + y * z) + PI * r * r * (x + y + z) + 4.0 / 3.0 * PI * r * r * r;
    let out = blended(&body, &items, want);
    assert_eq!(out.faces.iter().filter(|f| matches!(f.surface, Surface::Sphere { .. })).count(), 8);
}

#[test]
fn fillet_ending_on_slanted_faces_ends_in_ellipses() {
    // A parallelogram plate: the slanted ends cut the fillet in ellipses; every section across the edge is 4 long.
    let body = prism(vec![polygon(&[(0., 0.), (4., 0.), (5., 2.), (1., 2.)])], 1.5);
    let e = edge_at(&body, Vec3::new(2.0, 0.0, 1.5));
    let r = 0.3;
    let out = blended(&body, &[Item { edge: e, kind: Kind::Fillet { radius: r } }], 8.0 * 1.5 - spandrel(r) * 4.0);
    assert_eq!(out.edges.iter().filter(|e| matches!(e.curve, Curve::Conic { .. })).count(), 2);
}

#[test]
fn fillet_plate_corner_with_slanted_walls() {
    // A triangular plate: at each top corner the top face is square to both walls (which meet at 60 degrees), so the
    // three equal fillets there meet in a sphere patch. Every edge filleted.
    let s = 4.0;
    let h = s * 3f64.sqrt() / 2.0;
    let t = 1.0;
    let r = 0.2;
    let body = prism(vec![polygon(&[(0., 0.), (s, 0.), (s / 2.0, h)])], t);
    let items: Vec<Item> = (0..body.edges.len()).map(|edge| Item { edge, kind: Kind::Fillet { radius: r } }).collect();
    // The rounded prism is the inner prism (the triangle offset in by r, sides s - 2 sqrt(3) r; height t - 2r) grown by
    // a ball of radius r. Steiner: V = V' + S' r + (1/2 sum of edge length x exterior dihedral) r^2 + 4/3 pi r^3; the
    // six horizontal edges turn pi/2, the three vertical ones 2 pi/3.
    let si = s - 2.0 * r * 3f64.sqrt();
    let (ai, per) = (si * si * 3f64.sqrt() / 4.0, 3.0 * si);
    let ti = t - 2.0 * r;
    let want = ai * ti + (2.0 * ai + per * ti) * r + (PI / 2.0 * per + PI * ti) * r * r + 4.0 / 3.0 * PI * r * r * r;
    blended(&body, &items, want);
}

// ---- 3. plane-cylinder edges ----

fn holed_plate(r: f64, t: f64) -> Body { prism(vec![polygon(&[(0., 0.), (4., 0.), (4., 4.), (0., 4.)]), vec![circle(2.0, 2.0, r)]], t) }

#[test]
fn fillet_hole_rim_is_a_torus() {
    let (big, t, q) = (0.5, 1.0, 0.2);
    let body = holed_plate(big, t);
    let e = circle_at(&body, Vec3::new(2.0, 2.0, t), big);
    // Pappus: the spandrel revolved at radius big + its centroid offset.
    let removed = 2.0 * PI * (spandrel(q) * big + spandrel_moment(q));
    let out = blended(&body, &[Item { edge: e, kind: Kind::Fillet { radius: q } }], 16.0 * t - PI * big * big * t - removed);
    assert!(out.faces.iter().any(|f| matches!(f.surface, Surface::Torus { .. })));
}

#[test]
fn chamfer_hole_rim_is_a_cone() {
    let (big, t, d) = (0.5, 1.0, 0.15);
    let body = holed_plate(big, t);
    let e = circle_at(&body, Vec3::new(2.0, 2.0, t), big);
    let removed = 2.0 * PI * (d * d / 2.0) * (big + d / 3.0);
    let out = blended(&body, &[Item { edge: e, kind: Kind::Chamfer { distance: d } }], 16.0 * t - PI * big * big * t - removed);
    assert!(out.faces.iter().any(|f| matches!(f.surface, Surface::Cone { .. })));
}

#[test]
fn fillet_both_rims_of_a_hole() {
    let (big, t, q) = (0.5, 1.0, 0.2);
    let body = holed_plate(big, t);
    let items = [circle_at(&body, Vec3::new(2.0, 2.0, t), big), circle_at(&body, Vec3::new(2.0, 2.0, 0.0), big)]
        .map(|edge| Item { edge, kind: Kind::Fillet { radius: q } });
    let removed = 2.0 * 2.0 * PI * (spandrel(q) * big + spandrel_moment(q));
    blended(&body, &items, 16.0 * t - PI * big * big * t - removed);
}

#[test]
fn round_bar_end_fillet_and_two_distance_chamfer() {
    let (big, len) = (1.0, 3.0);
    let body = prism(vec![vec![circle(0.0, 0.0, big)]], len);
    let top = circle_at(&body, Vec3::new(0.0, 0.0, len), big);
    let q = 0.25;
    let removed = 2.0 * PI * (spandrel(q) * big - spandrel_moment(q));
    blended(&body, &[Item { edge: top, kind: Kind::Fillet { radius: q } }], PI * big * big * len - removed);
    // Chamfer the bottom end: 0.3 across the end face, 0.1 along the bar.
    let bottom = circle_at(&body, Vec3::new(0.0, 0.0, 0.0), big);
    let end = plane_face(&body, Z.scale(-1.0), Vec3::new(0.0, 0.0, 0.0));
    let (sp, sc) = (0.3, 0.1);
    blended(&body, &[Item { edge: bottom, kind: Kind::ChamferDistances { first: sp, second: sc, face: end } }],
        PI * big * big * len - 2.0 * PI * (sp * sc / 2.0) * (big - sp / 3.0));
}

#[test]
fn pipe_end_fillet_inside_and_out() {
    let (ro, ri, len, q) = (1.0, 0.7, 2.0, 0.1);
    let body = prism(vec![vec![circle(0.0, 0.0, ro)], vec![circle(0.0, 0.0, ri)]], len);
    let items = [circle_at(&body, Vec3::new(0.0, 0.0, len), ro), circle_at(&body, Vec3::new(0.0, 0.0, len), ri)]
        .map(|edge| Item { edge, kind: Kind::Fillet { radius: q } });
    let removed = 2.0 * PI * (spandrel(q) * ro - spandrel_moment(q)) + 2.0 * PI * (spandrel(q) * ri + spandrel_moment(q));
    blended(&body, &items, PI * (ro * ro - ri * ri) * len - removed);
}

#[test]
fn hole_rim_in_two_half_arcs_blends_through_both_seams() {
    let (big, t, q) = (0.5, 1.0, 0.2);
    let half = |t0: f64| Curve::Arc { c: p2(2.0, 2.0), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r: big, t0, t1: t0 + PI };
    let body = prism(vec![polygon(&[(0., 0.), (4., 0.), (4., 4.), (0., 4.)]), vec![half(0.0), half(PI)]], t);
    let rims: Vec<usize> = (0..body.edges.len()).filter(|&e| matches!(body.edges[e].curve, Curve::Arc { c, .. } if (c.z - t).abs() < 1e-12)).collect();
    assert_eq!(rims.len(), 2);
    let items: Vec<Item> = rims.iter().map(|&edge| Item { edge, kind: Kind::Fillet { radius: q } }).collect();
    let removed = 2.0 * PI * (spandrel(q) * big + spandrel_moment(q));
    blended(&body, &items, 16.0 * t - PI * big * big * t - removed);
    // Blending only half the rim is refused by name.
    assert!(refused(&body, &items[..1]).contains("round its rim"));
}

#[test]
fn stepped_shaft_shoulder_fillets_concave_and_convex() {
    // A revolved shaft: radius 2 for 1, then radius 1 for 2. The shoulder's inner rim (annulus to the radius-1 wall) is
    // concave: its fillet adds the spandrel revolved at 1 + its centroid offset; the outer rim (radius-2 wall to the
    // annulus) is convex.
    let section = Section { p: Vec3::new(0.0, 0.0, 0.0), a: Z, x: Vec3::new(1.0, 0.0, 0.0),
        points: vec![(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 3.0), (0.0, 3.0)] };
    let body = revolve(&section, TOL).unwrap();
    let base = PI * 4.0 + PI * 2.0;
    let inner = circle_at(&body, Vec3::new(0.0, 0.0, 1.0), 1.0);
    let q = 0.2;
    blended(&body, &[Item { edge: inner, kind: Kind::Fillet { radius: q } }], base + 2.0 * PI * (spandrel(q) * 1.0 + spandrel_moment(q)));
    let outer = circle_at(&body, Vec3::new(0.0, 0.0, 1.0), 2.0);
    blended(&body, &[Item { edge: inner, kind: Kind::Chamfer { distance: q } }, Item { edge: outer, kind: Kind::Fillet { radius: q } }],
        base + 2.0 * PI * (q * q / 2.0) * (1.0 + q / 3.0) - 2.0 * PI * (spandrel(q) * 2.0 - spandrel_moment(q)));
}

// ---- refusals by name ----

#[test]
fn a_fillet_that_consumes_a_face_is_refused() {
    let body = block(4.0, 3.0, 0.5);
    let e = edge_at(&body, Vec3::new(2.0, 0.0, 0.5));
    let m = refused(&body, &[Item { edge: e, kind: Kind::Fillet { radius: 0.6 } }]);
    assert!(m.contains("consumed"), "{m}");
}

#[test]
fn unequal_fillets_and_mixed_blends_at_a_corner_are_refused() {
    let (a, b, c) = (4.0, 3.0, 2.0);
    let body = block(a, b, c);
    let edges = top_loop(&body, a, b, c);
    let m = refused(&body, &[Item { edge: edges[0], kind: Kind::Fillet { radius: 0.3 } }, Item { edge: edges[1], kind: Kind::Fillet { radius: 0.4 } }]);
    assert!(m.contains("different radii"), "{m}");
    let m = refused(&body, &[Item { edge: edges[0], kind: Kind::Fillet { radius: 0.3 } }, Item { edge: edges[1], kind: Kind::Chamfer { distance: 0.3 } }]);
    assert!(m.contains("fillet") && m.contains("chamfer"), "{m}");
}

#[test]
fn a_blend_running_into_a_hole_is_refused() {
    let body = prism(vec![polygon(&[(0., 0.), (4., 0.), (4., 4.), (0., 4.)]), vec![circle(2.0, 0.6, 0.3)]], 1.0);
    let e = edge_at(&body, Vec3::new(2.0, 0.0, 1.0));
    let m = refused(&body, &[Item { edge: e, kind: Kind::Chamfer { distance: 0.4 } }]);
    assert!(m.contains("meet or cross") || m.contains("run into"), "{m}");
    // Clear of the hole it is fine: 0.2 leaves 0.1.
    blended(&body, &[Item { edge: e, kind: Kind::Chamfer { distance: 0.2 } }], 16.0 - PI * 0.09 - 0.02 * 4.0);
}

#[test]
fn a_hole_rim_fillet_deeper_than_the_plate_is_refused() {
    let body = holed_plate(0.5, 0.3);
    let e = circle_at(&body, Vec3::new(2.0, 2.0, 0.3), 0.5);
    let m = refused(&body, &[Item { edge: e, kind: Kind::Fillet { radius: 0.35 } }]);
    assert!(m.contains("consumed"), "{m}");
}

#[test]
fn a_bar_end_fillet_reaching_the_axis_is_refused() {
    let body = prism(vec![vec![circle(0.0, 0.0, 1.0)]], 3.0);
    let e = circle_at(&body, Vec3::new(0.0, 0.0, 3.0), 1.0);
    assert!(refused(&body, &[Item { edge: e, kind: Kind::Fillet { radius: 1.0 } }]).contains("axis"));
}

#[test]
fn edges_on_unsupported_surfaces_are_refused() {
    let ellipse = Curve::Conic { c: p2(0.0, 0.0), u: Vec3::new(2.0, 0.0, 0.0), v: Vec3::new(0.0, 1.0, 0.0), t0: 0.0, t1: 2.0 * PI };
    let body = prism(vec![vec![ellipse]], 1.0);
    let e = (0..body.edges.len()).find(|&e| matches!(body.edges[e].curve, Curve::Conic { c, .. } if c.z == 1.0)).unwrap();
    assert!(refused(&body, &[Item { edge: e, kind: Kind::Fillet { radius: 0.1 } }]).contains("blends are supported"));
}

#[test]
fn bad_requests_are_invalid_input() {
    let body = block(4.0, 3.0, 2.0);
    let e = edge_at(&body, Vec3::new(2.0, 0.0, 2.0));
    for items in [vec![Item { edge: e, kind: Kind::Fillet { radius: -1.0 } }],
                  vec![Item { edge: e, kind: Kind::Chamfer { distance: 0.1 } }, Item { edge: e, kind: Kind::Chamfer { distance: 0.1 } }],
                  vec![Item { edge: 999, kind: Kind::Chamfer { distance: 0.1 } }],
                  vec![Item { edge: e, kind: Kind::ChamferDistances { first: 0.1, second: 0.2, face: plane_face(&body, Z.scale(-1.0), Vec3::new(0., 0., 0.)) } }]] {
        assert!(matches!(blend(&body, &items, TOL), Err(KernelError::InvalidInput(_))), "{items:?}");
    }
}

// ---- the worker op ----

#[test]
fn op_19_blends_and_reports_volume_and_origins() {
    let (a, b, c, r) = (4.0, 3.0, 2.0, 0.5);
    let body = block(a, b, c);
    let e = edge_at(&body, Vec3::new(a / 2.0, 0.0, c));
    let top = plane_face(&body, Z, Vec3::new(0.0, 0.0, c));
    let other = edge_at(&body, Vec3::new(a / 2.0, b, c));
    let mut data = write_body(&body);
    data.extend_from_slice(&[2.0, 0.0, e as f64, r, 2.0, other as f64, 0.3, 0.2, top as f64]);
    let request = parse_request(&encode_request(OP_BLEND, TOL, &data)).unwrap();
    let reply = handle(&request).unwrap();
    let mut reader = Reader::with_tolerance(&reply, TOL);
    let out = reader.body().unwrap();
    let (volume, exact, n) = (reader.next().unwrap(), reader.next().unwrap(), reader.count().unwrap());
    close(volume, a * b * c - spandrel(r) * a - 0.3 * 0.2 / 2.0 * a);
    assert_eq!(exact, 1.0);
    assert_eq!(n, out.faces.len());
    let origins: Vec<(f64, f64)> = (0..n).map(|_| (reader.next().unwrap(), reader.next().unwrap())).collect();
    reader.finish().unwrap();
    assert!(origins.contains(&(1.0, e as f64)) && origins.contains(&(1.0, other as f64)));
    assert_eq!(origins.iter().filter(|o| o.0 == 0.0).count(), body.faces.len());
}

// ---- reflex corners, other angles, more refusals ----

fn l_top_loop(body: &Body) -> Vec<usize> {
    let pts = [(0., 0.), (4., 0.), (4., 1.), (1., 1.), (1., 3.), (0., 3.)];
    (0..6).map(|i| { let (a, b) = (pts[i], pts[(i + 1) % 6]); edge_at(body, Vec3::new((a.0 + b.0) / 2.0, (a.1 + b.1) / 2.0, 2.0)) }).collect()
}

#[test]
fn chamfer_l_top_loop_mitres_at_a_reflex_corner() {
    // 45 degree chamfers: the slab below the top keeps the L offset in by delta at each depth, area A - P delta +
    // (convex - reflex corners) delta^2 = 6 - 14 delta + 4 delta^2; removed = P d^2 / 2 - 4 d^3 / 3.
    let body = l_block();
    let d = 0.3;
    let items: Vec<Item> = l_top_loop(&body).into_iter().map(|edge| Item { edge, kind: Kind::Chamfer { distance: d } }).collect();
    blended(&body, &items, L_VOLUME - (14.0 * d * d / 2.0 - 4.0 * d * d * d / 3.0));
}

#[test]
fn fillet_l_top_loop_mitres_at_a_reflex_corner() {
    // Each side's spandrel prism is shortened by the centroid offset at a convex end, lengthened at a reflex end:
    // 10 convex ends, 2 reflex.
    let body = l_block();
    let r = 0.3;
    let items: Vec<Item> = l_top_loop(&body).into_iter().map(|edge| Item { edge, kind: Kind::Fillet { radius: r } }).collect();
    blended(&body, &items, L_VOLUME - (spandrel(r) * 14.0 - 8.0 * spandrel_moment(r)));
}

#[test]
fn blends_on_a_sixty_degree_edge() {
    // An equilateral prism's vertical edge: the faces open at phi = 60 degrees.
    let (s, len) = (4.0, 2.0);
    let body = prism(vec![polygon(&[(0., 0.), (s, 0.), (s / 2.0, s * 3f64.sqrt() / 2.0)])], len);
    let e = edge_at(&body, Vec3::new(0.0, 0.0, len / 2.0));
    let base = s * s * 3f64.sqrt() / 4.0 * len;
    let r = 0.3;
    // Setback r / tan 30 = r sqrt 3; the kite r t less the sector of (pi - phi).
    let t = r * 3f64.sqrt();
    blended(&body, &[Item { edge: e, kind: Kind::Fillet { radius: r } }], base - (r * t - r * r * (PI - PI / 3.0) / 2.0) * len);
    // Distance 0.4 on the x-axis face at 50 degrees: the other leg by the sine rule; triangle area d s sin(phi) / 2.
    let face = plane_face(&body, Vec3::new(0.0, -1.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
    let (d, angle, phi) = (0.4, 50f64.to_radians(), PI / 3.0);
    let other = d * angle.sin() / (angle + phi).sin();
    blended(&body, &[Item { edge: e, kind: Kind::ChamferAngle { distance: d, angle, face } }], base - d * other * phi.sin() / 2.0 * len);
}

#[test]
fn convex_and_concave_blends_meeting_are_refused() {
    let body = l_block();
    let vertical = edge_at(&body, Vec3::new(1.0, 1.0, 1.0));
    let top = edge_at(&body, Vec3::new(2.5, 1.0, 2.0));
    let m = refused(&body, &[Item { edge: vertical, kind: Kind::Fillet { radius: 0.2 } }, Item { edge: top, kind: Kind::Fillet { radius: 0.2 } }]);
    assert!(m.contains("convex and a concave"), "{m}");
}

#[test]
fn a_rim_fillet_reaching_another_hole_is_refused() {
    let body = prism(vec![polygon(&[(0., 0.), (4., 0.), (4., 4.), (0., 4.)]), vec![circle(1.5, 2.0, 0.4)], vec![circle(2.5, 2.0, 0.4)]], 1.0);
    let e = circle_at(&body, Vec3::new(1.5, 2.0, 1.0), 0.4);
    let m = refused(&body, &[Item { edge: e, kind: Kind::Chamfer { distance: 0.25 } }]);
    assert!(m.contains("meet or cross"), "{m}");
    blended(&body, &[Item { edge: e, kind: Kind::Chamfer { distance: 0.15 } }],
        16.0 - 2.0 * PI * 0.16 - 2.0 * PI * (0.15 * 0.15 / 2.0) * (0.4 + 0.05));
}

#[test]
fn rim_of_a_bored_hole_from_the_boolean_cut() {
    // A hole made by the kernel's own cut (bore domain), as Bob's parts arrive: fillet one rim, chamfer the other.
    let plate = block(4.0, 3.0, 1.0);
    let tool = extrude(&Profile { loops: vec![vec![circle(2.0, 1.5, 0.5)]], n: Z, low: -1.0, high: 2.0 }, TOL).unwrap();
    let body = gc_kernel::boolean::cut(&plate, &tool, TOL).unwrap();
    topology::check(&body, TOL).unwrap();
    let rims: Vec<usize> = (0..body.edges.len()).filter(|&e| matches!(body.edges[e].curve, Curve::Arc { .. })).collect();
    let top: Vec<usize> = rims.iter().copied().filter(|&e| body.edges[e].curve.start().z == 1.0).collect();
    let bottom: Vec<usize> = rims.iter().copied().filter(|&e| body.edges[e].curve.start().z == 0.0).collect();
    let (q, d) = (0.15, 0.1);
    let mut items: Vec<Item> = top.iter().map(|&edge| Item { edge, kind: Kind::Fillet { radius: q } }).collect();
    items.extend(bottom.iter().map(|&edge| Item { edge, kind: Kind::Chamfer { distance: d } }));
    let base = 12.0 - PI * 0.25;
    blended(&body, &items, base - 2.0 * PI * (spandrel(q) * 0.5 + spandrel_moment(q)) - 2.0 * PI * (d * d / 2.0) * (0.5 + d / 3.0));
}

#[test]
fn hole_rims_in_a_rectangular_tube_wall() {
    // Bob's survey case: a round hole through one wall of a 4 x 3 x 1/4 tube (made by the kernel's cut); both rims
    // (outside and inside the wall) filleted. Each is hole-like: the spandrel revolved at 0.5 + its centroid offset.
    let outer = polygon(&[(0., 0.), (4., 0.), (4., 3.), (0., 3.)]);
    let inner = polygon(&[(0.25, 0.25), (3.75, 0.25), (3.75, 2.75), (0.25, 2.75)]);
    let tube = prism(vec![outer, inner], 10.0);
    let y = Vec3::new(0.0, 1.0, 0.0);
    let bore = Curve::Arc { c: Vec3::new(2.0, 0.0, 5.0), n: y, x: Vec3::new(1.0, 0.0, 0.0), r: 0.5, t0: 0.0, t1: 2.0 * PI };
    let tool = extrude(&Profile { loops: vec![vec![bore]], n: y, low: -1.0, high: 0.5 }, TOL).unwrap();
    let body = gc_kernel::boolean::cut(&tube, &tool, TOL).unwrap();
    let rims: Vec<usize> = (0..body.edges.len()).filter(|&e| matches!(body.edges[e].curve, Curve::Arc { .. })).collect();
    assert_eq!(rims.len(), 2);
    let q = 0.05;
    let items: Vec<Item> = rims.iter().map(|&edge| Item { edge, kind: Kind::Fillet { radius: q } }).collect();
    let base = (12.0 - 3.5 * 2.5) * 10.0 - PI * 0.25 * 0.25;
    blended(&body, &items, base - 2.0 * 2.0 * PI * (spandrel(q) * 0.5 + spandrel_moment(q)));
}

#[test]
fn a_blend_cutting_into_an_internal_void_is_refused() {
    // A sealed void just under the top front edge: no face the chamfer touches shows it, only the void's own edges
    // inside the chamfer's region. The void is a second shell: a small box turned inside out.
    let mut body = block(4.0, 3.0, 2.0);
    let void = extrude(&Profile { loops: vec![polygon(&[(1.0, 0.1), (3.0, 0.1), (3.0, 1.0), (1.0, 1.0)])], n: Z, low: 1.85, high: 1.95 }, TOL).unwrap();
    let (nv, ne) = (body.vertices.len(), body.edges.len());
    body.vertices.extend(void.vertices.iter().copied());
    body.edges.extend(void.edges.iter().map(|e| gc_kernel::neutral::Edge { curve: e.curve.clone(), start: e.start + nv, end: e.end + nv }));
    body.faces.extend(void.faces.iter().map(|f| gc_kernel::neutral::Face { surface: f.surface.clone(), sense: !f.sense,
        loops: f.loops.iter().map(|lp| lp.iter().rev().map(|&(e, fwd)| (e + ne, !fwd)).collect()).collect() }));
    topology::check(&body, TOL).unwrap();
    let base = 24.0 - 2.0 * 0.9 * 0.1;
    close(mass::measure(&body, TOL).unwrap().volume, base);
    let e = edge_at(&body, Vec3::new(2.0, 0.0, 2.0));
    let m = refused(&body, &[Item { edge: e, kind: Kind::Chamfer { distance: 0.3 } }]);
    assert!(m.contains("run into edge"), "{m}");
    // A smaller chamfer clears it (0.1 + 0.05 > 0.12).
    blended(&body, &[Item { edge: e, kind: Kind::Chamfer { distance: 0.12 } }], base - 0.12 * 0.12 / 2.0 * 4.0);
}

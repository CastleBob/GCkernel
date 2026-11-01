//! Blind pockets in plates (Claude 2026-10-31): a sketch-profile prism that ends inside a Flat. Each removed volume is
//! the profile's area inside the plate (closed form from its dimensions) times the depth; every result passes the
//! topology check and tessellates; the flat nester's prism profile (op 13) still refuses a pocketed plate.

use gc_kernel::boolean::cut;
use gc_kernel::geom::{Curve, Vec3};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{Body, Profile};
use gc_kernel::topology::check;
use gc_kernel::KernelError;
use std::f64::consts::{PI, TAU};

const TOL: f64 = 1e-6;

fn line(a: (f64, f64), b: (f64, f64), z: f64) -> Curve { Curve::Line { a: Vec3::new(a.0, a.1, z), b: Vec3::new(b.0, b.1, z) } }

/// A polygon of lines at height z.
fn polygon(points: &[(f64, f64)], z: f64) -> Vec<Curve> {
    (0..points.len()).map(|i| line(points[i], points[(i + 1) % points.len()], z)).collect()
}

/// A full circle about (x, y) at height z, counter-clockwise about +z.
fn circle(x: f64, y: f64, r: f64, z: f64) -> Vec<Curve> {
    vec![Curve::Arc { c: Vec3::new(x, y, z), n: Vec3::new(0., 0., 1.), x: Vec3::new(1., 0., 0.), r, t0: 0., t1: TAU }]
}

/// A rounded rectangle centred (x, y), w x h, corner radius r, at height z, counter-clockwise.
fn rounded(x: f64, y: f64, w: f64, h: f64, r: f64, z: f64) -> Vec<Curve> {
    let (a, b) = (w / 2. - r, h / 2. - r);
    let arc = |cx: f64, cy: f64, t0: f64| Curve::Arc { c: Vec3::new(cx, cy, z), n: Vec3::new(0., 0., 1.), x: Vec3::new(1., 0., 0.), r, t0, t1: t0 + PI / 2. };
    vec![
        line((x - a, y - h / 2.), (x + a, y - h / 2.), z),
        arc(x + a, y - b, -PI / 2.),
        line((x + w / 2., y - b), (x + w / 2., y + b), z),
        arc(x + a, y + b, 0.),
        line((x + a, y + h / 2.), (x - a, y + h / 2.), z),
        arc(x - a, y + b, PI / 2.),
        line((x - w / 2., y + b), (x - w / 2., y - b), z),
        arc(x - a, y - b, PI),
    ]
}

fn prism(loops: Vec<Vec<Curve>>, low: f64, high: f64) -> Body {
    gc_kernel::extrude::extrude(&Profile { loops, n: Vec3::new(0., 0., 1.), low, high }, TOL).unwrap()
}

/// Bob's plates: 10 x 8, 0.5 thick, top at z 0.
fn plate() -> Body { prism(vec![polygon(&[(0., 0.), (10., 0.), (10., 8.), (0., 8.)], 0.)], -0.5, 0.) }

/// A pocket tool from just above the top down to `depth` (GCmaker's sketch cut: padding above, the depth below).
fn pocket(loops: Vec<Vec<Curve>>, depth: f64) -> Body { prism(loops, -depth, 0.01) }

fn removed(body: &Body, tool: &Body) -> (f64, Body) {
    let result = cut(body, tool, TOL).unwrap_or_else(|e| panic!("{e}"));
    check(&result, TOL).unwrap();
    gc_kernel::tessellate::tessellate(&result, 0.01).unwrap();
    (measure(body, TOL).unwrap().volume - measure(&result, TOL).unwrap().volume, result)
}

#[test]
fn rounded_rectangle_and_round_pockets_remove_their_area_times_depth() {
    let (got, result) = removed(&plate(), &pocket(vec![rounded(4., 3., 3., 2., 0.4, 0.)], 0.2));
    let exact = (3. * 2. - (4. - PI) * 0.16) * 0.2;
    assert!((got - exact).abs() < 1e-12, "removed {got}, exact {exact}");
    // A floor face at -0.2 and the walls: bottom, top, floor and 4 + 8 side faces.
    assert_eq!(result.faces.len(), 3 + 4 + 8);
    let (got, _) = removed(&plate(), &pocket(vec![circle(7., 5., 1., 0.)], 0.25));
    assert!((got - PI * 0.25).abs() < 1e-12, "removed {got}");
    // Two pockets in turn: the second in a body that is no longer a prism is refused by name, not cut.
    let (_, once) = removed(&plate(), &pocket(vec![circle(7., 5., 1., 0.)], 0.25));
    match cut(&once, &pocket(vec![circle(2., 5., 0.5, 0.)], 0.1), TOL) {
        Err(KernelError::Unsupported(_)) => {}
        other => panic!("expected a named refusal, got {:?}", other.map(|b| b.faces.len())),
    }
}

#[test]
fn a_pocket_over_the_edge_is_a_rebate_and_one_across_the_plate_a_groove() {
    // A rebate: 2 x 3 over the plate's right edge (x 9..11), 0.3 deep: the part inside is 1 x 3.
    let (got, _) = removed(&plate(), &pocket(vec![polygon(&[(9., 2.), (11., 2.), (11., 5.), (9., 5.)], 0.)], 0.3));
    assert!((got - 3. * 0.3).abs() < 1e-12, "removed {got}");
    // Flush with the edge (x 8..10): the same.
    let (got, _) = removed(&plate(), &pocket(vec![polygon(&[(8., 2.), (10., 2.), (10., 5.), (8., 5.)], 0.)], 0.3));
    assert!((got - 2. * 3. * 0.3).abs() < 1e-12, "removed {got}");
    // A groove right across (y 3..4 over x -1..11): the top splits in two, held together by the floor.
    let (got, result) = removed(&plate(), &pocket(vec![polygon(&[(-1., 3.), (11., 3.), (11., 4.), (-1., 4.)], 0.)], 0.1));
    assert!((got - 10. * 0.1).abs() < 1e-12, "removed {got}");
    assert_eq!(check(&result, TOL).unwrap().genus(), 0.);
    // A tool covering the whole top: the plate is thinner.
    let (got, _) = removed(&plate(), &pocket(vec![polygon(&[(-1., -1.), (11., -1.), (11., 9.), (-1., 9.)], 0.)], 0.2));
    assert!((got - 80. * 0.2).abs() < 1e-10, "removed {got}");
}

#[test]
fn counterbore_about_a_through_hole_and_a_pocket_from_the_bottom() {
    // A plate with a 0.5 through hole, counterbored 1.0 diameter 0.2 deep from the top.
    let holed = prism(vec![polygon(&[(0., 0.), (10., 0.), (10., 8.), (0., 8.)], 0.), circle(5., 4., 0.25, 0.)], -0.5, 0.);
    let (got, result) = removed(&holed, &pocket(vec![circle(5., 4., 0.5, 0.)], 0.2));
    assert!((got - PI * (0.25 - 0.0625) * 0.2).abs() < 1e-12, "removed {got}");
    assert_eq!(check(&result, TOL).unwrap().genus(), 1.);
    // From the bottom: a 2 x 1 pocket 0.15 up from z -0.5.
    let tool = prism(vec![polygon(&[(1., 1.), (3., 1.), (3., 2.), (1., 2.)], 0.)], -0.6, -0.35);
    let (got, result) = removed(&plate(), &tool);
    assert!((got - 2. * 0.15).abs() < 1e-12, "removed {got}");
    assert!(result.vertices.iter().any(|v| (v.z + 0.35).abs() < 1e-12));
}

#[test]
fn an_oblique_pocket_and_one_touching_the_outline_at_a_point() {
    // A triangle pocket with oblique sides, 0.25 deep: area 6.
    let (got, _) = removed(&plate(), &pocket(vec![polygon(&[(2., 2.), (6., 2.), (4., 5.)], 0.)], 0.25));
    assert!((got - 6. * 0.25).abs() < 1e-12, "removed {got}");
    // A square standing on its corner at the plate's edge (touching it at one point): non-manifold, refused by name.
    match cut(&plate(), &pocket(vec![polygon(&[(5., 0.), (6., 1.), (5., 2.), (4., 1.)], 0.)], 0.25), TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("point"), "{m}"),
        other => panic!("expected a named refusal, got {:?}", other.map(|b| b.faces.len())),
    }
}

#[test]
fn the_nesters_prism_profile_still_refuses_a_pocketed_plate() {
    let (_, result) = removed(&plate(), &pocket(vec![circle(7., 5., 1., 0.)], 0.25));
    match gc_kernel::rebuild::prism_profile(&result, Vec3::new(0., 0., 1.), TOL) {
        Err(KernelError::Unsupported(m)) => assert!(!m.is_empty()),
        other => panic!("op 13 must refuse a pocketed plate, got {:?}", other.map(|p| p.loops.len())),
    }
}

#[test]
fn round_pocket_from_the_bottom_and_holes_after_a_pocket() {
    // From the bottom: a 1.5 diameter pocket 0.2 up from z -0.5 (cut turned half a turn and back).
    let tool = prism(vec![circle(3., 3., 0.75, 0.)], -0.6, -0.3);
    let (got, result) = removed(&plate(), &tool);
    assert!((got - PI * 0.5625 * 0.2).abs() < 1e-12, "removed {got}");
    assert!(result.vertices.iter().all(|v| v.z >= -0.5 - 1e-12 && v.z <= 1e-12));
    // A through hole after a pocket: in the pocket's floor and elsewhere in the plate.
    let (_, pocketed) = removed(&plate(), &pocket(vec![rounded(4., 3., 3., 2., 0.4, 0.)], 0.2));
    for (x, y) in [(4., 3.), (8., 6.)] {
        let hole = prism(vec![circle(x, y, 0.25, 0.)], -1., 0.5);
        match cut(&pocketed, &hole, TOL) {
            Ok(r) => {
                check(&r, TOL).unwrap();
                let depth = if (x, y) == (4., 3.) { 0.3 } else { 0.5 };
                let got = measure(&pocketed, TOL).unwrap().volume - measure(&r, TOL).unwrap().volume;
                assert!((got - PI * 0.0625 * depth).abs() < 1e-12, "hole at ({x}, {y}) removed {got}");
            }
            Err(KernelError::Unsupported(m)) => eprintln!("hole at ({x}, {y}) after a pocket refused: {m}"),
            Err(e) => panic!("{e}"),
        }
    }
}

#[test]
fn a_box_notch_in_a_tube_end_keeps_its_exact_domain() {
    // GCmaker's Rect Tube sketch cut (captured): an 8 x 6 x 0.25 tube (corner radii 1 and 0.75), notched by a box
    // x 0..3.5 through its height from the end y 0 to y 1. A box tool is not a pocket (the transverse box domain cuts
    // it, exactly as before): removed the section inside x <= 3.5, (19 + pi/2) - (16.75 + 0.28125 pi), times 1.
    let read = |name: &str| {
        let path = format!("{}/tests/fixtures/survey_gaps/rect_tube_sketch_pocket/{name}", env!("CARGO_MANIFEST_DIR"));
        let data: Vec<f64> = std::fs::read(path).unwrap().chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
        gc_kernel::neutral::Reader::with_tolerance(&data, TOL).body().unwrap()
    };
    let (got, _) = removed(&read("body.f64"), &read("tool.f64"));
    let exact = 2.25 + 0.21875 * PI;
    assert!((got - exact).abs() < 1e-10, "removed {got}, section {exact}");
}

//! The kernel rebuilding a straight prism from its cap (Claude): GCmaker's captured I-beam, a holed block, refusals.

use std::f64::consts::PI;

use gc_kernel::extrude::extrude;
use gc_kernel::geom::{Curve, Vec3};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{write_body, Body, Profile, Reader};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_EXTRUDE, OP_PRISM_PROFILE, OP_REBUILD_PRISM};
use gc_kernel::rebuild::{prism_profile, rebuild_prism};
use gc_kernel::topology::check;
use gc_kernel::KernelError;

const TOL: f64 = 1e-7;

fn fixture(name: &str) -> Body {
    let bytes = std::fs::read(format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap();
    let data: Vec<f64> = bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
    let mut r = Reader::new(&data);
    let body = r.body().unwrap();
    r.finish().unwrap();
    body
}

#[test]
fn captured_ibeam_is_rebuilt_by_the_kernel() {
    let base = fixture("rat_hole/base.f64");
    let rebuilt = rebuild_prism(&base, 1e-6).unwrap();
    check(&rebuilt, 1e-6).unwrap();
    let (a, b) = (measure(&base, 1e-6).unwrap(), measure(&rebuilt, 1e-6).unwrap());
    assert!((a.volume - b.volume).abs() < 1e-9 * a.volume && (a.area - b.area).abs() < 1e-9 * a.area);
    // The same through the worker protocol.
    let reply = handle(&parse_request(&encode_request(OP_REBUILD_PRISM, 1e-6, &write_body(&base))).unwrap()).unwrap();
    assert_eq!(reply, write_body(&rebuilt));
}

fn line(a: (f64, f64), b: (f64, f64)) -> Curve { Curve::Line { a: Vec3::new(a.0, a.1, 0.0), b: Vec3::new(b.0, b.1, 0.0) } }

#[test]
fn holed_block_along_any_axis_is_rebuilt() {
    let z = Vec3::new(0.0, 0.0, 1.0);
    let rect = vec![line((0.0, 0.0), (10.0, 0.0)), line((10.0, 0.0), (10.0, 8.0)), line((10.0, 8.0), (0.0, 8.0)), line((0.0, 8.0), (0.0, 0.0))];
    let hole = vec![Curve::Arc { c: Vec3::new(5.0, 4.0, 0.0), n: z, x: Vec3::new(1.0, 0.0, 0.0), r: 1.0, t0: 0.0, t1: 2.0 * PI }];
    let block = extrude(&Profile { loops: vec![rect, hole], n: z, low: -1.0, high: 3.0 }, TOL).unwrap();
    let rebuilt = rebuild_prism(&block, TOL).unwrap();
    assert!((measure(&rebuilt, TOL).unwrap().volume - 4.0 * (80.0 - PI)).abs() < 1e-12);
}

/// A block with a blind pocket is not a straight prism of either cap: refused, not approximated.
#[test]
fn a_non_prism_is_refused() {
    let z = Vec3::new(0.0, 0.0, 1.0);
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64| vec![line((x0, y0), (x1, y0)), line((x1, y0), (x1, y1)), line((x1, y1), (x0, y1)), line((x0, y1), (x0, y0))];
    let block = extrude(&Profile { loops: vec![rect(0.0, 0.0, 10.0, 8.0)], n: z, low: 0.0, high: 4.0 }, TOL).unwrap();
    let pocket = extrude(&Profile { loops: vec![rect(2.0, 2.0, 4.0, 5.0)], n: z, low: 1.0, high: 5.0 }, TOL).unwrap();
    let cut = gc_kernel::boolean::cut(&block, &pocket, TOL).unwrap();
    assert!(matches!(rebuild_prism(&cut, TOL), Err(KernelError::Unsupported(_))));
}

/// A bottom edge bulged out of its side plane is now refused by the input check itself (edge-on-surface).
#[test]
fn a_bottom_edge_bulged_off_its_side_face_is_refused_by_the_input_check() {
    let z = Vec3::new(0.0, 0.0, 1.0);
    let rect = vec![line((0.0, 0.0), (10.0, 0.0)), line((10.0, 0.0), (10.0, 8.0)), line((10.0, 8.0), (0.0, 8.0)), line((0.0, 8.0), (0.0, 0.0))];
    let mut block = extrude(&Profile { loops: vec![rect], n: z, low: 0.0, high: 4.0 }, TOL).unwrap();
    let bottom = block.edges.iter().position(|e| matches!(e.curve, Curve::Line { a, b } if a.z == 0.0 && b.z == 0.0)).unwrap();
    let Curve::Line { a, b } = block.edges[bottom].curve else { unreachable!() };
    let (mid, r) = (a.add(b).scale(0.5), a.distance(b) / 2.0);
    block.edges[bottom].curve = Curve::Arc { c: mid, n: z, x: a.sub(mid).scale(1.0 / r), r, t0: 0.0, t1: PI };
    assert!(matches!(rebuild_prism(&block, TOL), Err(KernelError::InvalidInput(m)) if m.contains("leaves face")));
}

/// A sound body that only the structural proof refuses: a 10 x 8 x 4 box whose y = 0 side is split at z = 2 into
/// two stacked faces. Topology, surfaces and measures are all right; the top edge's side face is bounded by a
/// ruling of 2, not 4 (and its neighbours have five edges), so it is not "one side face per cap edge" and is
/// refused (conservatively), not rebuilt.
#[test]
fn a_split_side_face_is_refused_by_the_structural_proof() {
    use gc_kernel::geom::Surface;
    use gc_kernel::neutral::{Edge, Face};
    let v = |x: f64, y: f64, z: f64| Vec3::new(x, y, z);
    let vertices = vec![v(0.0, 0.0, 0.0), v(10.0, 0.0, 0.0), v(10.0, 8.0, 0.0), v(0.0, 8.0, 0.0), v(0.0, 0.0, 4.0), v(10.0, 0.0, 4.0),
                        v(10.0, 8.0, 4.0), v(0.0, 8.0, 4.0), v(0.0, 0.0, 2.0), v(10.0, 0.0, 2.0)];
    let pairs = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 8), (8, 4), (1, 9), (9, 5), (2, 6), (3, 7), (8, 9)];
    let edges = pairs.iter().map(|&(a, b)| Edge { curve: Curve::Line { a: vertices[a], b: vertices[b] }, start: a, end: b }).collect();
    let plane = |p: Vec3, n: Vec3, x: Vec3| Surface::Plane { p, n, x };
    let (x, y, z) = (v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0));
    let face = |surface, sense, lp: Vec<(usize, bool)>| Face { surface, sense, loops: vec![lp] };
    let faces = vec![
        face(plane(v(0.0, 0.0, 0.0), z, x), false, vec![(3, false), (2, false), (1, false), (0, false)]),
        face(plane(v(0.0, 0.0, 4.0), z, x), true, vec![(4, true), (5, true), (6, true), (7, true)]),
        face(plane(v(0.0, 0.0, 0.0), y.scale(-1.0), x), true, vec![(0, true), (10, true), (14, false), (8, false)]),
        face(plane(v(0.0, 0.0, 0.0), y.scale(-1.0), x), true, vec![(14, true), (11, true), (4, false), (9, false)]),
        face(plane(v(10.0, 0.0, 0.0), x, y), true, vec![(1, true), (12, true), (5, false), (11, false), (10, false)]),
        face(plane(v(0.0, 8.0, 0.0), y, x), true, vec![(2, true), (13, true), (6, false), (12, false)]),
        face(plane(v(0.0, 0.0, 0.0), x.scale(-1.0), y), true, vec![(3, true), (8, true), (9, true), (7, false), (13, false)]),
    ];
    let body = Body { vertices, edges, faces };
    check(&body, TOL).unwrap();
    let m = measure(&body, TOL).unwrap();
    assert!((m.volume - 320.0).abs() < 1e-12 && (m.area - 304.0).abs() < 1e-12);
    match rebuild_prism(&body, TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("side face"), "{m}"),   // refused by the proof
        other => panic!("{other:?}"),
    }
}

/// Astra's probe: both X-side planes of a box moved one inch off their edges. Topology closes and the measures
/// cancel (volume 320, area 304), but the planes no longer hold their boundaries: reported, never repaired.
#[test]
fn planes_moved_off_their_edges_are_reported_not_rebuilt() {
    let z = Vec3::new(0.0, 0.0, 1.0);
    let rect = vec![line((0.0, 0.0), (10.0, 0.0)), line((10.0, 0.0), (10.0, 8.0)), line((10.0, 8.0), (0.0, 8.0)), line((0.0, 8.0), (0.0, 0.0))];
    let mut block = extrude(&Profile { loops: vec![rect], n: z, low: 0.0, high: 4.0 }, TOL).unwrap();
    let mut moved = 0;
    for f in &mut block.faces {
        if let gc_kernel::geom::Surface::Plane { p, n, .. } = &mut f.surface {
            if n.x.abs() > 0.5 { p.x += 1.0; moved += 1; }
        }
    }
    assert_eq!(moved, 2);
    let m = measure(&block, TOL).unwrap();
    assert!((m.volume - 320.0).abs() < 1e-9 && (m.area - 304.0).abs() < 1e-9);   // why measures cannot admit
    assert!(matches!(check(&block, TOL), Err(KernelError::InvalidInput(m)) if m.contains("leaves face")));
    let reply = handle(&parse_request(&encode_request(OP_REBUILD_PRISM, TOL, &write_body(&block))).unwrap());
    assert!(matches!(reply, Err(KernelError::InvalidInput(_))), "{reply:?}");
    // The Rust function itself, not only the protocol, refuses it.
    assert!(matches!(rebuild_prism(&block, TOL), Err(KernelError::InvalidInput(m)) if m.contains("leaves face")));
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Curve> {
    vec![line((x0, y0), (x1, y0)), line((x1, y0), (x1, y1)), line((x1, y1), (x0, y1)), line((x0, y1), (x0, y0))]
}

/// A 6 x 4 x 0.5 plate (top at z = 0.5) drilled through by kernel cuts: a hole, a second hole and a corner notch.
fn cut_plate() -> Body {
    let z = Vec3::new(0.0, 0.0, 1.0);
    let mut plate = extrude(&Profile { loops: vec![rect(0.0, 0.0, 6.0, 4.0)], n: z, low: 0.0, high: 0.5 }, TOL).unwrap();
    let drill = |cx: f64, cy: f64, r: f64| extrude(&Profile {
        loops: vec![vec![Curve::Arc { c: Vec3::new(cx, cy, 0.0), n: z, x: Vec3::new(1.0, 0.0, 0.0), r, t0: 0.0, t1: 2.0 * PI }]],
        n: z, low: -1.0, high: 2.0 }, TOL).unwrap();
    let notch = extrude(&Profile { loops: vec![rect(5.0, 3.0, 7.0, 5.0)], n: z, low: -1.0, high: 2.0 }, TOL).unwrap();
    for tool in [drill(1.5, 2.0, 0.5), drill(3.5, 2.0, 0.25), notch] {
        plate = gc_kernel::boolean::cut(&plate, &tool, TOL).unwrap();
    }
    plate
}

/// Op 13: the plate's own outline and holes, on its top, and sweeping that profile back gives the same solid.
#[test]
fn a_cut_plate_gives_its_profile() {
    let plate = cut_plate();
    let z = Vec3::new(0.0, 0.0, 1.0);
    let profile = prism_profile(&plate, z, TOL).unwrap();
    assert_eq!(profile.loops.len(), 3);
    assert!(profile.loops.iter().flatten().all(|c| (c.start().z - 0.5).abs() < 1e-12 && (c.end().z - 0.5).abs() < 1e-12));
    assert!((profile.high - profile.low - 0.5).abs() < 1e-12 && profile.n.add(z).length() < 1e-12);
    let expected = 0.5 * (24.0 - 1.0 - PI * 0.25 - PI * 0.0625);
    let rebuilt = extrude(&profile, TOL).unwrap();
    assert!((measure(&rebuilt, TOL).unwrap().volume - expected).abs() < 1e-9);
    // Through the worker: the reply is op 2's own request.
    let mut request = vec![0.0, 0.0, 1.0];
    request.extend(write_body(&plate));
    let reply = handle(&parse_request(&encode_request(OP_PRISM_PROFILE, TOL, &request)).unwrap()).unwrap();
    assert_eq!(reply, gc_kernel::neutral::write_profile(&profile));
    let swept = handle(&parse_request(&encode_request(OP_EXTRUDE, TOL, &reply)).unwrap()).unwrap();
    assert!((swept[swept.len() - 1] - expected).abs() < 1e-9);
}

/// The direction is the caller's: a plain box is a prism three ways and gives the cap it is asked for; a direction
/// no cap faces is refused by name.
#[test]
fn the_profile_follows_the_requested_direction() {
    let z = Vec3::new(0.0, 0.0, 1.0);
    let block = extrude(&Profile { loops: vec![rect(0.0, 0.0, 10.0, 8.0)], n: z, low: 0.0, high: 4.0 }, TOL).unwrap();
    let x = prism_profile(&block, Vec3::new(1.0, 0.0, 0.0), TOL).unwrap();
    assert!((x.high - 10.0).abs() < 1e-12 && x.loops[0].iter().all(|c| (c.start().x - 10.0).abs() < 1e-12));
    let tilted = Vec3::new(1.0, 0.0, 1.0).unit().unwrap();
    assert!(matches!(prism_profile(&block, tilted, TOL), Err(KernelError::Unsupported(m)) if m.contains("0 faces face the top")));
}

/// A blind pocket's floor also faces the top: refused by name, never nested as a through cut.
#[test]
fn a_pocketed_plate_is_refused_by_name() {
    let z = Vec3::new(0.0, 0.0, 1.0);
    let block = extrude(&Profile { loops: vec![rect(0.0, 0.0, 10.0, 8.0)], n: z, low: 0.0, high: 4.0 }, TOL).unwrap();
    let pocket = extrude(&Profile { loops: vec![rect(2.0, 2.0, 4.0, 5.0)], n: z, low: 1.0, high: 5.0 }, TOL).unwrap();
    let cut = gc_kernel::boolean::cut(&block, &pocket, TOL).unwrap();
    match prism_profile(&cut, z, TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("a pocket"), "{m}"),
        other => panic!("{other:?}"),
    }
}
/// Rectangle cut-outs through a plate, inside it and across its edge: still one top and one bottom.
#[test]
fn rectangle_cutouts_keep_one_top() {
    let z = Vec3::new(0.0, 0.0, 1.0);
    let mut plate = extrude(&Profile { loops: vec![rect(0.0, 0.0, 6.0, 4.0)], n: z, low: 0.0, high: 0.5 }, TOL).unwrap();
    for (x0, y0, x1, y1) in [(1.0, 1.0, 2.0, 3.0), (3.0, 1.5, 4.5, 2.5), (5.5, -1.0, 7.0, 1.0)] {
        let tool = extrude(&Profile { loops: vec![rect(x0, y0, x1, y1)], n: z, low: -1.0, high: 2.0 }, TOL).unwrap();
        plate = gc_kernel::boolean::cut(&plate, &tool, TOL).unwrap();
        let profile = prism_profile(&plate, z, TOL).unwrap_or_else(|e| panic!("after ({x0},{y0})-({x1},{y1}): {e}"));
        assert!(profile.loops.len() >= 2);
        assert!((measure(&extrude(&profile, TOL).unwrap(), TOL).unwrap().volume - measure(&plate, TOL).unwrap().volume).abs() < 1e-9);
    }
}
/// The seams between a cap's cells are dropped and the cut outline's pieces join: a plate with one rectangle
/// cut-out gives exactly two four-line loops.
#[test]
fn a_cutout_plate_gives_clean_loops() {
    let z = Vec3::new(0.0, 0.0, 1.0);
    let plate = extrude(&Profile { loops: vec![rect(0.0, 0.0, 6.0, 4.0)], n: z, low: 0.0, high: 0.5 }, TOL).unwrap();
    let tool = extrude(&Profile { loops: vec![rect(1.0, 1.0, 2.0, 3.0)], n: z, low: -1.0, high: 2.0 }, TOL).unwrap();
    let cut = gc_kernel::boolean::cut(&plate, &tool, TOL).unwrap();
    let profile = prism_profile(&cut, z, TOL).unwrap();
    let mut sizes: Vec<usize> = profile.loops.iter().map(|l| l.len()).collect();
    sizes.sort_unstable();
    assert_eq!(sizes, vec![4, 4], "{:?}", profile.loops);
}
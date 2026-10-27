//! The kernel rebuilding a straight prism from its cap (Claude): GCmaker's captured I-beam, a holed block, refusals.

use std::f64::consts::PI;

use gc_kernel::extrude::extrude;
use gc_kernel::geom::{Curve, Vec3};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{write_body, Body, Profile, Reader};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_REBUILD_PRISM};
use gc_kernel::rebuild::rebuild_prism;
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

/// The structural proof, not the measures, catches a bottom edge that is not the top edge moved: one bottom line of a
/// block replaced by a half circle between the same two points.
#[test]
fn a_bottom_edge_that_is_not_the_moved_top_edge_is_refused_by_the_proof() {
    let z = Vec3::new(0.0, 0.0, 1.0);
    let rect = vec![line((0.0, 0.0), (10.0, 0.0)), line((10.0, 0.0), (10.0, 8.0)), line((10.0, 8.0), (0.0, 8.0)), line((0.0, 8.0), (0.0, 0.0))];
    let mut block = extrude(&Profile { loops: vec![rect], n: z, low: 0.0, high: 4.0 }, TOL).unwrap();
    let bottom = block.edges.iter().position(|e| matches!(e.curve, Curve::Line { a, b } if a.z == 0.0 && b.z == 0.0)).unwrap();
    let Curve::Line { a, b } = block.edges[bottom].curve else { unreachable!() };
    let (mid, r) = (a.add(b).scale(0.5), a.distance(b) / 2.0);
    block.edges[bottom].curve = Curve::Arc { c: mid, n: z, x: a.sub(mid).scale(1.0 / r), r, t0: 0.0, t1: PI };
    assert!(block.edges[bottom].curve.end().distance(b) < 1e-12);
    match rebuild_prism(&block, TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("bottom edge is not its top edge moved"), "{m}"),
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
}

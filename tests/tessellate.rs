//! Display tessellation (Claude, for Astra's viewport): every mesh must be closed (each directed triangle edge met
//! by its reverse), wound outward (triangle normals agree with the vertex normals), and enclose the exact volume
//! within the chord tolerance times the curved area. Real captured parts and the shop survey replay first.

use std::f64::consts::PI;

use gc_kernel::geom::{Curve, Surface, Vec3};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{write_body, Body, Profile, Reader};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_CUT, OP_TESSELLATE};
use gc_kernel::tessellate::{tessellate, write_mesh, Mesh};
use std::collections::HashMap;

const TOL: f64 = 1e-6;
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn values(path: &str) -> Vec<f64> {
    let bytes = std::fs::read(format!("{}/tests/fixtures/{path}", env!("CARGO_MANIFEST_DIR"))).unwrap();
    bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect()
}
fn read(data: &[f64]) -> Body { Reader::with_tolerance(data, TOL).body().unwrap() }
fn cut(body: &[f64], tool: &[f64]) -> Result<Vec<f64>, String> {
    let data: Vec<f64> = body.iter().chain(tool).copied().collect();
    handle(&parse_request(&encode_request(OP_CUT, TOL, &data)).unwrap()).map_err(|e| e.to_string())
}

/// Checks a mesh; returns its enclosed volume.
fn sound(body: &Body, m: &Mesh) -> f64 {
    let key = |v: Vec3| (v.x.to_bits(), v.y.to_bits(), v.z.to_bits());
    let mut directed: HashMap<_, i64> = HashMap::new();
    let mut volume = 0.0;
    let o = m.positions[0];
    for t in &m.triangles {
        let [a, b, c] = t.map(|i| m.positions[i]);
        for (u, w) in [(a, b), (b, c), (c, a)] {
            if key(u) == key(w) { continue; }
            *directed.entry((key(u), key(w))).or_default() += 1;
            *directed.entry((key(w), key(u))).or_default() -= 1;
        }
        let (a, b, c) = (a.sub(o), b.sub(o), c.sub(o));
        volume += a.dot(b.cross(c)) / 6.0;
        let n = b.sub(a).cross(c.sub(a));
        if n.length() > 1e-12 {
            let avg = t.iter().fold(Vec3::new(0.0, 0.0, 0.0), |s, &i| s.add(m.normals[i]));
            assert!(n.dot(avg) > 0.0, "a triangle winds inward: {a:?} {b:?} {c:?}");
        }
    }
    let open = directed.values().filter(|&&k| k != 0).count();
    assert_eq!(open, 0, "{open} unmatched directed edges");
    assert_eq!(m.edges.len(), body.edges.len());
    assert!(m.display_edges.windows(2).all(|e| e[0] < e[1]));
    assert!(m.display_edges.iter().all(|&e| e < body.edges.len()));
    assert!(m.positions.iter().chain(&m.normals).all(|v| v.is_finite()));
    volume
}

/// Exact volume, mesh volume and their allowed gap (chord times the curved area).
fn check(body: &Body, chord: f64) -> (f64, f64) {
    let original = write_body(body);
    let m = tessellate(body, chord).unwrap();
    assert_eq!(write_body(body), original, "display must not edit exact topology");
    let mesh_volume = sound(body, &m);
    let exact = measure(body, TOL).unwrap().volume;
    let curved: f64 = m.triangles.iter().zip(&m.faces)
        .filter(|(_, f)| matches!(body.faces[**f].surface, Surface::Cylinder { .. }))
        .map(|(t, _)| { let [a, b, c] = t.map(|i| m.positions[i]); b.sub(a).cross(c.sub(a)).length() / 2.0 })
        .sum();
    assert!((mesh_volume - exact).abs() <= chord * curved + 1e-9 * exact.abs(), "mesh {mesh_volume} exact {exact} curved {curved}");
    (exact, mesh_volume)
}

fn line(a: (f64, f64), b: (f64, f64)) -> Curve { Curve::Line { a: Vec3::new(a.0, a.1, 0.0), b: Vec3::new(b.0, b.1, 0.0) } }
fn plate() -> Body {
    let loops = vec![vec![line((0.0, 0.0), (10.0, 0.0)), line((10.0, 0.0), (10.0, 8.0)), line((10.0, 8.0), (0.0, 8.0)), line((0.0, 8.0), (0.0, 0.0))]];
    gc_kernel::extrude::extrude(&Profile { loops, n: Z, low: 0.0, high: 0.5 }, TOL).unwrap()
}
fn rod(cx: f64, cy: f64, r: f64) -> Body {
    let c = Curve::Arc { c: Vec3::new(cx, cy, 0.0), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r, t0: 0.0, t1: 2.0 * PI };
    gc_kernel::extrude::extrude(&Profile { loops: vec![vec![c]], n: Z, low: -1.0, high: 1.5 }, TOL).unwrap()
}

#[test]
fn native_plate_with_separate_hole_and_rotated_slot_displays() {
    for degrees in [0.0_f64, 37.0, -37.0, 90.0, 180.0] {
    let angle = degrees.to_radians();
    let x = Vec3::new(angle.cos(), angle.sin(), 0.0);
    let y = Z.cross(x);
    let p = |u: f64, v: f64| Vec3::new(7.0, 5.0, 0.0).add(x.scale(u)).add(y.scale(v));
    let slot = gc_kernel::extrude::extrude(&Profile { n: Z, low: -1.0, high: 1.5, loops: vec![vec![
        Curve::Line { a: p(-1.0, -0.5), b: p(1.0, -0.5) },
        Curve::Arc { c: p(1.0, 0.0), n: Z, x, r: 0.5, t0: -PI / 2.0, t1: PI / 2.0 },
        Curve::Line { a: p(1.0, 0.5), b: p(-1.0, 0.5) },
        Curve::Arc { c: p(-1.0, 0.0), n: Z, x, r: 0.5, t0: PI / 2.0, t1: 3.0 * PI / 2.0 },
    ]] }, TOL).unwrap();
    let body = gc_kernel::boolean::cut(&plate(), &rod(2.0, 2.0, 0.5), TOL).unwrap();
    let body = gc_kernel::boolean::cut(&body, &slot, TOL).unwrap();
    check(&body, 1e-3);
    }
}

#[test]
fn plate_and_plate_with_hole_are_closed_and_exact() {
    let p = plate();
    let (exact, mesh) = check(&p, 1e-3);
    assert!((exact - 40.0).abs() < 1e-12 && (mesh - 40.0).abs() < 1e-12);
    let holed = gc_kernel::boolean::cut(&p, &rod(5.0, 4.0, 1.0), TOL).unwrap();
    let (exact, mesh) = check(&holed, 1e-3);
    assert!((exact - (40.0 - 0.5 * PI)).abs() < 1e-12);
    // The hole's polygon is inscribed: the mesh keeps a little more material, never less.
    assert!(mesh > exact);
    // A finer chord comes closer, about in proportion.
    let (_, finer) = check(&holed, 1e-5);
    assert!((finer - exact) < (mesh - exact) / 50.0);
}

#[test]
fn plate_with_hole_display_fixture_matches_the_worker() {
    // The first joint fixture for Astra's viewport (GCm Claude-GcKernel-Status.md Log 30).
    let holed = gc_kernel::boolean::cut(&plate(), &rod(5.0, 4.0, 1.0), TOL).unwrap();
    let body = write_body(&holed);
    let mut request = vec![1e-3];
    request.extend_from_slice(&body);
    let reply = handle(&parse_request(&encode_request(OP_TESSELLATE, TOL, &request)).unwrap()).unwrap();
    assert_eq!(reply, write_mesh(&tessellate(&holed, 1e-3).unwrap()));
    let dir = format!("{}/tests/fixtures/display/plate_hole", env!("CARGO_MANIFEST_DIR"));
    let bytes = |v: &[f64]| v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>();
    // Filtering edge strokes must not change the triangulated surface, normals or per-triangle face IDs.
    let previous = values("display/plate_hole/mesh.f64");
    let surface_end = 5 + 6 * reply[2] as usize + 4 * reply[3] as usize;
    assert_eq!(&previous[5..surface_end], &reply[5..surface_end]);
    if std::env::var("GCKERNEL_WRITE_DISPLAY").is_ok() {
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(format!("{dir}/body.f64"), bytes(&body)).unwrap();
        std::fs::write(format!("{dir}/mesh.f64"), bytes(&reply)).unwrap();
    }
    assert_eq!(values("display/plate_hole/body.f64"), body);
    assert_eq!(values("display/plate_hole/mesh.f64"), reply);
}

#[test]
fn display_hides_periodic_and_partition_seams_but_keeps_hole_rims_and_stock_corners() {
    let p = plate();
    let pm = tessellate(&p, 1e-3).unwrap();
    assert_eq!(pm.display_edges.len(), 12, "all sharp box edges are physical");
    let r = rod(5.0, 4.0, 1.0);
    let rm = tessellate(&r, 1e-3).unwrap();
    assert_eq!(rm.display_edges.len(), 2, "draw the two rims, not the cylinder's periodic seam");
    assert!(rm.display_edges.iter().all(|&e| matches!(r.edges[e].curve, Curve::Arc { .. })));
    let holed = gc_kernel::boolean::cut(&p, &r, TOL).unwrap();
    let hm = tessellate(&holed, 1e-3).unwrap();
    assert_eq!(hm.display_edges.iter().filter(|&&e| matches!(holed.edges[e].curve, Curve::Arc { .. })).count(), 2);

    for (fixture, steps) in [("rat_hole", 2), ("tangent_hole", 3), ("parts/rect_tube", 5)] {
        let body = replay(fixture, steps);
        let m = tessellate(&body, 1e-3).unwrap();
        assert!(m.display_edges.len() < body.edges.len(), "{fixture}: expected artificial subdivisions");
        assert!(!m.display_edges.is_empty());
        let reply = write_mesh(&m);
        let mut at = 5 + 6 * m.positions.len() + 4 * m.triangles.len();
        for &id in &m.display_edges {
            assert_eq!(reply[at] as usize, id, "must retain the original topology ID");
            at += 2 + 3 * reply[at + 1] as usize;
        }
        assert_eq!(at, reply.len());
        eprintln!("{fixture}: {} display edges / {} body edges", m.display_edges.len(), body.edges.len());
    }
}

fn replay(dir: &str, steps: usize) -> Body {
    let mut body = values(&format!("{dir}/base.f64"));
    for k in 1..=steps { body = cut(&body, &values(&format!("{dir}/step{k:02}.f64"))).unwrap(); }
    read(&body)
}

#[test]
fn captured_parts_tessellate_closed() {
    check(&replay("rat_hole", 0), 1e-3);
    check(&replay("rat_hole", 2), 1e-3);
    check(&replay("tangent_hole", 3), 1e-3);   // Bob's hole touching the flange face (ends at a double root)
    check(&replay("circle_over_fillet", 1), 1e-3);
    check(&replay("parts/rect_tube", 5), 1e-3);
    check(&replay("parts/str1", 1), 1e-3);
}

#[test]
fn shop_survey_tessellates_closed() {
    let dir = format!("{}/tests/fixtures/survey", env!("CARGO_MANIFEST_DIR"));
    let mut done = 0;
    for entry in std::fs::read_dir(&dir).unwrap() {
        let path = entry.unwrap().path();
        if !path.is_dir() { continue; }
        let name = format!("survey/{}", path.file_name().unwrap().to_string_lossy());
        let mut body = values(&format!("{name}/base.f64"));
        check(&read(&body), 1e-3);
        let mut k = 1;
        while std::path::Path::new(&format!("{dir}/../{name}/step{k:02}.f64")).exists() {
            match cut(&body, &values(&format!("{name}/step{k:02}.f64"))) {
                Ok(next) => body = next,
                Err(_) => break,   // outside the cutting domain today; the survey records which
            }
            k += 1;
        }
        check(&read(&body), 1e-3);
        done += 1;
    }
    assert!(done >= 10, "{done}");
}

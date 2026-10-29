//! Independent axis/apex acceptance for solid revolve. No shared worker/build is changed.
//! Expectations use elementary cone formulae and combinatorial manifold checks.
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;
use gc_kernel::geom::Vec3;
use gc_kernel::neutral::{write_body, Body, Reader};
use gc_kernel::revolve::{revolve, Section};
use gc_kernel::{mass, protocol, tessellate, topology};

const TOL: f64 = 1e-7;
const CHORD: f64 = 0.02;

fn section(points: Vec<(f64, f64)>) -> Section {
    Section { p: Vec3::new(0., 0., 0.), a: Vec3::new(0., 1., 0.),
        x: Vec3::new(1., 0., 0.), points }
}

fn near(actual: f64, expected: f64) {
    assert!(actual.is_finite() && (actual - expected).abs() <= 1e-8 * expected.abs().max(1.),
        "{actual} != analytic {expected}");
}

fn payload(s: &Section) -> Vec<f64> {
    let mut values = vec![s.p.x, s.p.y, s.p.z, s.a.x, s.a.y, s.a.z,
        s.x.x, s.x.y, s.x.z, s.points.len() as f64];
    for &(radius, height) in &s.points { values.extend([radius, height]); }
    values
}

/// Weld only exactly shared positions, including +0/-0. Shared edge tessellation must
/// make adjacent face boundaries identical; a positional tolerance could hide a crack.
fn closed_manifold_mesh(body: &Body, expected_volume: f64, expected_area: f64) {
    let original = write_body(body);
    let mesh = tessellate::tessellate(body, CHORD).expect("apex body must tessellate");
    assert_eq!(original, write_body(body), "display changed authoritative geometry");
    assert!(!mesh.triangles.is_empty());
    assert_eq!(mesh.triangles.len(), mesh.faces.len());
    assert_eq!(mesh.positions.len(), mesh.normals.len());
    assert!(mesh.positions.iter().chain(mesh.normals.iter()).all(|p| p.is_finite()));
    let bits = |x: f64| if x == 0. { 0 } else { x.to_bits() };
    let mut ids = HashMap::new();
    let welded: Vec<usize> = mesh.positions.iter().map(|p| {
        let next = ids.len();
        *ids.entry((bits(p.x), bits(p.y), bits(p.z))).or_insert(next)
    }).collect();
    let mut edge_uses: HashMap<(usize, usize), (usize, i32)> = HashMap::new();
    let mut links: HashMap<usize, HashMap<usize, Vec<usize>>> = HashMap::new();
    let origin = mesh.positions[0];
    let mut volume = 0.;
    for (tri, &face) in mesh.triangles.iter().zip(&mesh.faces) {
        assert!(face < body.faces.len());
        assert!(tri.iter().all(|&i| i < mesh.positions.len()));
        let [a, b, c] = tri.map(|i| mesh.positions[i]);
        let [u, v, w] = tri.map(|i| welded[i]);
        assert!(u != v && v != w && w != u, "collapsed apex triangle on face {face}");
        let cross = b.sub(a).cross(c.sub(a));
        assert!(cross.length() > 1e-13, "zero-area triangle on face {face}");
        let normal = tri.iter().fold(Vec3::new(0., 0., 0.), |n, &i| n.add(mesh.normals[i]));
        assert!(cross.dot(normal) > 0., "triangle orientation opposes outward normal on face {face}");
        volume += a.sub(origin).dot(b.sub(origin).cross(c.sub(origin))) / 6.;
        for (s, e) in [(u, v), (v, w), (w, u)] {
            let entry = edge_uses.entry((s.min(e), s.max(e))).or_default();
            entry.0 += 1;
            entry.1 += if s < e { 1 } else { -1 };
        }
        for (at, left, right) in [(u, v, w), (v, w, u), (w, u, v)] {
            let link = links.entry(at).or_default();
            link.entry(left).or_default().push(right);
            link.entry(right).or_default().push(left);
        }
    }
    assert!(edge_uses.values().all(|&(count, winding)| count == 2 && winding == 0),
        "mesh has an open, multiply-used, or inconsistently oriented edge");
    // A closed manifold vertex, including a cone apex, has one connected cycle as
    // its link. Edge incidence alone misses two closed shells pinched at a vertex.
    for (vertex, link) in &links {
        assert!(link.values().all(|adjacent| adjacent.len() == 2),
            "vertex {vertex} link is not a cycle");
        let mut pending = vec![*link.keys().next().unwrap()];
        let mut visited = HashSet::new();
        while let Some(v) = pending.pop() {
            if visited.insert(v) { pending.extend(link[&v].iter().copied()); }
        }
        assert_eq!(visited.len(), link.len(), "vertex {vertex} joins disconnected fans");
    }
    assert_eq!(links.len() as isize - edge_uses.len() as isize + mesh.triangles.len() as isize,
        2, "a cone or bicone boundary must have sphere topology");
    assert!(volume > 0., "mesh is globally inside out");
    assert!((volume - expected_volume).abs() <= CHORD * expected_area + 1e-8);
}

fn verify(s: &Section, expected_volume: f64, expected_area: f64) {
    let body = revolve(s, TOL).expect("valid cone/apex section must build");
    topology::check(&body, TOL).unwrap();
    let measured = mass::measure(&body, TOL).unwrap();
    assert!(measured.exact, "polygonal cones have closed-form mass");
    near(measured.volume, expected_volume);
    near(measured.area, expected_area);
    let wire = write_body(&body);
    let mut reader = Reader::with_tolerance(&wire, TOL);
    let restored = reader.body().unwrap();
    reader.finish().unwrap();
    topology::check(&restored, TOL).unwrap();
    let roundtrip = write_body(&restored);
    assert_eq!(wire.len(), roundtrip.len());
    for (&a, &b) in wire.iter().zip(&roundtrip) {
        assert!((a - b).abs() <= 1e-12 * a.abs().max(1.));
    }
    near(mass::measure(&restored, TOL).unwrap().volume, expected_volume);
    closed_manifold_mesh(&restored, expected_volume, expected_area);
    let reply = protocol::handle(&protocol::Request {
        op: protocol::OP_REVOLVE, units: 1, tolerance: TOL, data: payload(s),
    }).expect("op 8 must admit the same apex section as the direct constructor");
    near(*reply.last().unwrap(), expected_volume);
    let mut reader = Reader::with_tolerance(&reply[..reply.len() - 1], TOL);
    let from_worker = reader.body().unwrap();
    reader.finish().unwrap();
    near(mass::measure(&from_worker, TOL).unwrap().area, expected_area);
}

#[test]
fn cone_and_inverted_cone_have_exact_mass_and_manifold_apices() {
    for points in [vec![(0., 0.), (2., 0.), (0., 3.)],
                   vec![(0., 0.), (2., 3.), (0., 3.)]] {
        for frame in 0..3 {
            for reverse in [false, true] {
                let mut s = section(points.clone());
                if reverse { s.points.reverse(); }
                if frame > 0 {
                    s.p = Vec3::new(17., -23., 91.);
                    s.a = Vec3::new(1., 2., 2.).unit().unwrap();
                    s.x = Vec3::new(2., -1., 0.).unit().unwrap();
                }
                if frame == 2 { s.a = s.a.scale(-1.); }
                verify(&s, 4. * PI, PI * (4. + 2. * 13_f64.sqrt()));
            }
        }
    }
}

#[test]
fn two_opposing_tips_with_positive_middle_radius_form_one_manifold_solid() {
    for reverse in [false, true] {
        let mut s = section(vec![(0., 0.), (2., 3.), (0., 8.)]);
        s.p = Vec3::new(-43., 19., 7.);
        s.a = Vec3::new(2., -2., 1.).unit().unwrap();
        s.x = Vec3::new(1., 1., 0.).unit().unwrap();
        if reverse { s.points.reverse(); }
        verify(&s, 32. * PI / 3., 2. * PI * (13_f64.sqrt() + 29_f64.sqrt()));
    }
}

#[test]
fn isolated_axis_contact_and_pinched_sections_are_refused() {
    for points in [
        // A simple triangle touching the axis only at an isolated point: its two
        // swept conical walls meet at that point with two disconnected vertex links.
        vec![(0., 0.), (2., -1.), (2., 1.)],
        // Two material lobes touch on the axis inside the closing axis segment.
        vec![(0., 0.), (2., 0.), (2., 2.), (0., 3.), (2., 4.), (2., 6.), (0., 6.)],
    ] {
        for reverse in [false, true] {
            let mut s = section(points.clone());
            if reverse { s.points.reverse(); }
            assert!(revolve(&s, TOL).is_err(), "nonmanifold axis contact was admitted");
            assert!(protocol::handle(&protocol::Request {
                op: protocol::OP_REVOLVE, units: 1, tolerance: TOL, data: payload(&s),
            }).is_err());
        }
    }
}

#[test]
fn direct_revolve_refuses_bad_tolerances_coordinates_and_overflow_before_building() {
    let valid = section(vec![(0., 0.), (2., 0.), (2., 3.), (0., 3.)]);
    for tolerance in [0., -1e-7, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(revolve(&valid, tolerance).is_err(), "bad direct tolerance {tolerance} admitted");
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for item in 0..5 {
            let mut s = valid.clone();
            match item {
                0 => s.p.x = bad,
                1 => s.a.y = bad,
                2 => s.x.z = bad,
                3 => s.points[1].0 = bad,
                _ => s.points[1].1 = bad,
            }
            assert!(revolve(&s, TOL).is_err());
        }
    }
    let enormous = section(vec![(0., 0.), (1e200, 0.), (0., 1e200)]);
    assert!(revolve(&enormous, TOL).is_err(), "overflowing analytic cone admitted");
    // Deliberately do not tessellate enormous radii/tiny chords here: the current
    // subdivision code needs a resource guard before such a regression is safe.
}

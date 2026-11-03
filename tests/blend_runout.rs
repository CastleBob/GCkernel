//! Blends running out past a short end face onto the curved face beyond (Claude 2026-11-03, Bob's L3x3x1/4 angle:
//! the leg's bottom end edge, whose end face at the toe is 0.0625 flat and then the 0.1875 toe round). Expected
//! volumes are the removed cross-section integrated against the toe's outline (Simpson on each smooth piece).

use std::collections::HashMap;
use std::f64::consts::PI;

use gc_kernel::blend::{blend, Item, Kind};
use gc_kernel::extrude::extrude;
use gc_kernel::geom::{Curve, Vec3};
use gc_kernel::neutral::{write_body, Body, Profile, Reader};
use gc_kernel::{mass, tessellate, topology, KernelError};

const TOL: f64 = 1e-7;
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);
const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);

fn p2(x: f64, y: f64) -> Vec3 { Vec3::new(x, y, 0.0) }
fn line(a: (f64, f64), b: (f64, f64)) -> Curve { Curve::Line { a: p2(a.0, a.1), b: p2(b.0, b.1) } }
fn quarter(cx: f64, cy: f64, r: f64, from: f64) -> Curve { Curve::Arc { c: p2(cx, cy), n: Z, x: X, r, t0: from, t1: from + PI / 2.0 } }

/// Bob's angle: legs 3, thickness 0.25, root fillet 0.375, toe rounds 0.1875, 3.5 long, all times `s`.
fn angle(s: f64) -> Body {
    let (w, t, ri, rt, len) = (3.0 * s, 0.25 * s, 0.375 * s, 0.1875 * s, 3.5 * s);
    let profile = vec![
        line((0.0, 0.0), (w, 0.0)),
        line((w, 0.0), (w, t - rt)),
        quarter(w - rt, t - rt, rt, 0.0),
        line((w - rt, t), (t + ri, t)),
        quarter(t + ri, t + ri, ri, PI).reversed(),
        line((t, t + ri), (t, w - rt)),
        quarter(t - rt, w - rt, rt, 0.0),
        line((t - rt, w), (0.0, w)),
        line((0.0, w), (0.0, 0.0)),
    ];
    extrude(&Profile { loops: vec![profile], n: Z, low: 0.0, high: len }, TOL).unwrap()
}

/// The edge whose curve passes through `p` at its parameter midpoint.
fn edge_at(body: &Body, p: Vec3) -> usize {
    let found: Vec<usize> = (0..body.edges.len()).filter(|&e| {
        let (t0, t1) = body.edges[e].curve.domain();
        body.edges[e].curve.eval((t0 + t1) / 2.0).0.distance(p) < 1e-9
    }).collect();
    assert_eq!(found.len(), 1, "edges at {p:?}: {found:?}");
    found[0]
}

/// How far the angle reaches along x at height y above its bottom face (the toe's outline), times `s`.
fn toe(y: f64, s: f64) -> f64 {
    let (w, flat, rt) = (3.0 * s, 0.0625 * s, 0.1875 * s);
    if y <= flat { w } else { w - rt + (rt * rt - (y - flat).powi(2)).sqrt() }
}

/// Simpson's rule on n (even) intervals.
fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;
    let sum: f64 = (0..=n).map(|i| {
        let w = if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
        w * f(a + h * i as f64)
    }).sum();
    sum * h / 3.0
}

/// What a fillet of radius r removes: at height y its spandrel is r - sqrt(r^2 - (r - y)^2) deep, as long as the toe
/// outline; with y = r - r sin t that is r^2 (1 - cos t) cos t toe(r - r sin t) dt, split where y passes the flat.
fn fillet_removed(r: f64, s: f64) -> f64 {
    let f = |t: f64| r * r * (1.0 - t.cos()) * t.cos() * toe(r - r * t.sin(), s);
    let flat = 0.0625 * s;
    if r <= flat { return simpson(f, 0.0, PI / 2.0, 4000); }
    let split = (1.0 - flat / r).asin();
    simpson(f, 0.0, split, 4000) + simpson(f, split, PI / 2.0, 4000)
}

/// What an equal chamfer d removes: at height y a strip d - y deep.
fn chamfer_removed(d: f64, s: f64) -> f64 {
    let f = |y: f64| (d - y) * toe(y, s);
    let flat = 0.0625 * s;
    if d <= flat { return simpson(f, 0.0, d, 4000); }
    simpson(f, 0.0, flat, 4000) + simpson(f, flat, d, 4000)
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

/// Blends the angle's bottom end edge and proves the result: topology, measure against the removed volume, neutral
/// round trip, closed outward mesh.
fn check(s: f64, kind: Kind, removed: f64) {
    let body = angle(s);
    let v0 = mass::measure(&body, TOL).unwrap().volume;
    let e = edge_at(&body, Vec3::new(1.5 * s, 0.0, 0.0));
    let result = blend(&body, &[Item { edge: e, kind }], TOL).unwrap_or_else(|err| panic!("{kind:?} at scale {s}: refused: {err}"));
    let out = result.body;
    assert_eq!(result.origins.len(), out.faces.len());
    topology::check(&out, TOL).unwrap();
    let m = mass::measure(&out, TOL).unwrap();
    let want = v0 - removed;
    assert!((m.volume - want).abs() <= 1e-9 * want, "{kind:?} at scale {s}: volume {} vs {want} (off {})", m.volume, m.volume - want);
    let wire = write_body(&out);
    let mut reader = Reader::with_tolerance(&wire, TOL);
    let back = reader.body().unwrap();
    reader.finish().unwrap();
    topology::check(&back, TOL).unwrap();
    mesh_check(&back, 0.0005 * s, m.volume, m.area);
}

#[test]
fn fillet_runs_out_over_the_toe_round() {
    for s in [1.0, 25.4] {
        for r in [0.05, 0.0625 * 1.5, 0.125, 0.2, 0.24] {
            check(s, Kind::Fillet { radius: r * s }, fillet_removed(r * s, s));
        }
    }
}

#[test]
fn chamfer_runs_out_over_the_toe_round() {
    for s in [1.0, 25.4] {
        for d in [0.05, 0.125, 0.2] {
            check(s, Kind::Chamfer { distance: d * s }, chamfer_removed(d * s, s));
        }
    }
}

#[test]
fn a_blend_wider_than_the_leg_is_refused_by_name() {
    let body = angle(1.0);
    let e = edge_at(&body, Vec3::new(1.5, 0.0, 0.0));
    match blend(&body, &[Item { edge: e, kind: Kind::Fillet { radius: 0.3 } }], TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("too wide for its end") && m.contains("reaches 0.3") && m.contains("(0.25 across)"), "{m}"),
        Err(other) => panic!("expected a named refusal, got {other}"),
        Ok(_) => panic!("expected a refusal"),
    }
}

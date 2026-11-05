//! Geometry primitives and op 7 TRANSFORM (Claude 2026-10-28, contract with Astra): points move and vectors only turn,
//! rigid motions compose in order and invert, frames are admitted only as given (never repaired), every curve kind and
//! a real captured body move exactly, and the worker batch reports a bad item by index without partial results.

use std::f64::consts::PI;

use gc_kernel::frame::{Frame, Point2, Point3, Rigid, Vec2};
use gc_kernel::geom::{pair_disc_coefficients, pair_disc_lower_bound, Curve, Surface, Vec3};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{put_curve, write_body, Body, Reader};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_CUT, OP_TRANSFORM};
use gc_kernel::topology;

const TOL: f64 = 1e-6;

fn values(name: &str) -> Vec<f64> {
    let bytes = std::fs::read(format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap();
    bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect()
}
fn run(op: u32, data: &[f64]) -> Result<Vec<f64>, String> {
    handle(&parse_request(&encode_request(op, TOL, data)).unwrap()).map_err(|e| e.to_string())
}
fn read(data: &[f64]) -> Body { Reader::with_tolerance(data, TOL).body().unwrap() }
fn close(a: Vec3, b: Vec3, eps: f64) -> bool { a.distance(b) <= eps }
fn p(x: f64, y: f64, z: f64) -> Point3 { Point3::new(x, y, z) }
fn v(x: f64, y: f64, z: f64) -> Vec3 { Vec3::new(x, y, z) }

/// A tilted motion used throughout: 37 degrees about a skew unit axis through (2, -1, 5), then a shift.
fn tilted() -> Rigid {
    let axis = v(1.0, 2.0, 2.0).scale(1.0 / 3.0);
    Rigid::rotation(axis, 37f64.to_radians(), p(2.0, -1.0, 5.0)).unwrap().then(&Rigid::translation(v(-3.0, 4.5, 0.25)).unwrap())
}

#[test]
fn points_move_and_vectors_only_turn() {
    let shift = Rigid::translation(v(1.0, 2.0, 3.0)).unwrap();
    assert_eq!(shift.point(p(1.0, 1.0, 1.0)), p(2.0, 3.0, 4.0));
    assert_eq!(shift.vector(v(1.0, 1.0, 1.0)), v(1.0, 1.0, 1.0));
    // A quarter turn about z through (1, 1, 0): (2, 1, 0) goes to (1, 2, 0); the x direction goes to y.
    let quarter = Rigid::rotation(v(0.0, 0.0, 1.0), PI / 2.0, p(1.0, 1.0, 0.0)).unwrap();
    assert!(close(quarter.point(p(2.0, 1.0, 0.0)).to_vec3(), v(1.0, 2.0, 0.0), 1e-15));
    assert!(close(quarter.vector(v(1.0, 0.0, 0.0)), v(0.0, 1.0, 0.0), 1e-15));
    // Point3 and Point2 differences are displacements.
    assert_eq!(p(3.0, 2.0, 1.0).sub(p(1.0, 1.0, 1.0)), v(2.0, 1.0, 0.0));
    assert_eq!(Point2::new(3.0, 2.0).sub(Point2::new(1.0, 1.0)), Vec2::new(2.0, 1.0));
    assert_eq!(Vec2::new(1.0, 0.0).cross(Vec2::new(0.0, 1.0)), 1.0);
    assert!(Vec2::new(0.0, 0.0).unit().is_err());
}

#[test]
fn motions_compose_in_order_and_invert() {
    let r = Rigid::rotation(v(0.0, 0.0, 1.0), PI / 2.0, p(0.0, 0.0, 0.0)).unwrap();
    let t = Rigid::translation(v(1.0, 0.0, 0.0)).unwrap();
    // Rotate then shift, against shift then rotate.
    assert!(close(r.then(&t).point(p(1.0, 0.0, 0.0)).to_vec3(), v(1.0, 1.0, 0.0), 1e-15));
    assert!(close(t.then(&r).point(p(1.0, 0.0, 0.0)).to_vec3(), v(0.0, 2.0, 0.0), 1e-15));
    let a = tilted();
    let b = Rigid::rotation(v(0.0, 1.0, 0.0), -1.1, p(4.0, 0.0, -2.0)).unwrap();
    let c = Rigid::translation(v(0.5, -7.0, 2.0)).unwrap();
    let q = p(3.25, -8.5, 11.0);
    assert!(close(a.then(&b).then(&c).point(q).to_vec3(), a.then(&b.then(&c)).point(q).to_vec3(), 1e-13));
    for k in 0..50 {
        let q = p(k as f64 * 0.37 - 9.0, 13.0 - k as f64 * 0.21, (k * k) as f64 * 0.01);
        let back = a.inverse().point(a.point(q));
        assert!(close(back.to_vec3(), q.to_vec3(), 1e-13 * (1.0 + q.to_vec3().length())), "{back:?} vs {q:?}");
        assert!(close(a.then(&a.inverse()).point(q).to_vec3(), q.to_vec3(), 1e-13 * (1.0 + q.to_vec3().length())));
    }
}

#[test]
fn frames_are_admitted_only_as_given() {
    let (c, s) = (0.6, 0.8);
    let a = Frame::new(p(1.0, 2.0, 3.0), v(c, s, 0.0), v(-s, c, 0.0)).unwrap();
    assert!(close(a.z(), v(0.0, 0.0, 1.0), 1e-15));
    assert!(close(a.point2(Point2::new(1.0, 0.0)).to_vec3(), v(1.6, 2.8, 3.0), 1e-15));
    assert!(close(a.local_point(a.point(p(0.3, -0.2, 0.7))).to_vec3(), v(0.3, -0.2, 0.7), 1e-15));
    for (x, y, why) in [
        (v(2.0, 0.0, 0.0), v(0.0, 1.0, 0.0), "not unit"),
        (v(1.0, 0.0, 0.0), v(0.0, 1.0 + 1e-9, 0.0), "not unit"),
        (v(1.0, 0.0, 0.0), v(1e-9, 1.0, 0.0), "right angles"),
        (v(f64::NAN, 0.0, 0.0), v(0.0, 1.0, 0.0), "not finite"),
    ] {
        let e = Frame::new(p(0.0, 0.0, 0.0), x, y).unwrap_err().to_string();
        assert!(e.contains(why), "{e}");
    }
    // from_frames takes origin to origin and axis to axis.
    let b = Frame::new(p(-4.0, 0.5, 9.0), v(0.0, 0.0, 1.0), v(1.0, 0.0, 0.0)).unwrap();
    let m = Rigid::from_frames(&a, &b);
    assert!(close(m.point(a.origin()).to_vec3(), b.origin().to_vec3(), 1e-14));
    assert!(close(m.vector(a.x()), b.x(), 1e-15) && close(m.vector(a.y()), b.y(), 1e-15) && close(m.vector(a.z()), b.z(), 1e-15));
    assert_eq!(m.frame(&a).origin(), m.point(a.origin()));
}

#[test]
fn only_proper_rotations_are_rigid() {
    let t = v(0.0, 0.0, 0.0);
    for (m, why) in [
        ([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], "reflection"),
        ([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]], "orthonormal"),
        ([[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], "orthonormal"),
        ([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0 + 1e-9]], "orthonormal"),
    ] {
        let e = Rigid::from_rows(m, t).unwrap_err().to_string();
        assert!(e.contains(why), "{e}");
    }
    assert!(Rigid::rotation(v(0.0, 0.0, 2.0), 1.0, p(0.0, 0.0, 0.0)).is_err());
    assert!(Rigid::translation(v(f64::INFINITY, 0.0, 0.0)).is_err());
    let (m, t) = tilted().rows();
    assert!(Rigid::from_rows(m, t).is_ok());
}

/// A real captured part after its cuts: rat_hole (two cuts: lines, arcs, conics) or tangent_hole (three cuts: Bob's
/// hole touching the flange face, with CylinderPair edges, four of them ending where their cylinders touch).
fn replay(part: &str, steps: usize) -> Body {
    let mut body = values(&format!("{part}/base.f64"));
    for step in 1..=steps {
        let tool = values(&format!("{part}/step{step:02}.f64"));
        body = run(OP_CUT, &body.iter().chain(&tool).copied().collect::<Vec<_>>()).unwrap();
    }
    read(&body)
}
fn rat_hole() -> Body { replay("rat_hole", 2) }

#[test]
fn every_curve_kind_and_a_real_body_move_exactly() {
    let m = tilted();
    let body = rat_hole();
    let touching = replay("tangent_hole", 3);
    let mut kinds = [0usize; 5];
    // A genuinely oblique ellipse arc: the rat holes' plane-cylinder crossings are circles (arcs) since cut operands'
    // near-cardinal directions are taken as exact (boolean::snap_directions).
    let ellipse = Curve::Conic { c: Vec3::new(1.0, 2.0, 3.0), u: Vec3::new(2.0, 0.0, 0.5), v: Vec3::new(0.0, 1.0, 0.0), t0: 0.3, t1: 2.1 };
    for c in body.edges.iter().chain(&touching.edges).map(|e| &e.curve).chain(std::iter::once(&ellipse)) {
        let moved = m.curve(c);
        let (t0, t1) = c.domain();
        for k in 0..=8 {
            let t = t0 + (t1 - t0) * k as f64 / 8.0;
            let ((q, d), (mq, md)) = (c.eval(t), moved.eval(t));
            let scale = 1.0 + q.length();
            // A CylinderPair point solves h from sqrt(S): where the cylinders touch (S at a double root, Bob's tangent
            // hole) rounding of 1e-16 in the moved data becomes about 1e-8 in the point, inside the kernel tolerance
            // ("agrees within the kernel tolerance is the same", Bob 2026-10-28). Every other kind moves to 1e-12.
            let (point_eps, tangent_eps) = if matches!(c, Curve::CylinderPair { .. }) { (TOL, 1e-4) } else { (1e-12, 1e-12) };
            assert!(close(mq, m.point(Point3::from_vec3(q)).to_vec3(), point_eps * scale), "{c:?} at {t}");
            assert!(close(md, m.vector(d), tangent_eps * (1.0 + d.length())), "{c:?} tangent at {t}");
        }
        kinds[match c { Curve::Line { .. } => 0, Curve::Arc { .. } => 1, Curve::Conic { .. } => 3, Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } => 4, Curve::Bezier { .. } => 5, Curve::RationalBezier { .. } => panic!("rational curves have dedicated tests") }] += 1;
        if let (Curve::CylinderPair { c, a, x, rho, p, b, r, t0, t1, .. },
            Curve::CylinderPair { c: c0, a: a0, x: x0, p: p0, b: b0, .. }) = (moved, c.clone()) {
            // The discriminant (whether the crossing is real along its trim) is invariant under the motion.
            let before = pair_disc_lower_bound(&pair_disc_coefficients(c0, a0, x0, rho, p0, b0, r), t0, t1, 64);
            let after = pair_disc_lower_bound(&pair_disc_coefficients(c, a, x, rho, p, b, r), t0, t1, 64);
            assert!((before - after).abs() <= 1e-9 * (1.0 + before.abs()), "{before} vs {after}");
        }
    }
    assert!(kinds[0] > 0 && kinds[1] > 0 && kinds[3] > 0 && kinds[4] > 0, "{kinds:?}");

    for part in [&body, &touching] {
        let moved = m.body(part);
        topology::check(&moved, TOL).unwrap();
        let (a, b) = (measure(part, TOL).unwrap(), measure(&moved, TOL).unwrap());
        assert!((a.volume - b.volume).abs() <= 1e-9 * a.volume, "{} vs {}", a.volume, b.volume);
        assert!((a.area - b.area).abs() <= 1e-9 * a.area, "{} vs {}", a.area, b.area);
    }
    let moved = m.body(&body);
    // And back again.
    let back = m.inverse().body(&moved);
    for (u, w) in body.vertices.iter().zip(&back.vertices) { assert!(close(*u, *w, 1e-12)); }

    let s = Surface::Cylinder { p: v(1.0, 0.0, 0.0), a: v(0.0, 0.0, 1.0), x: v(1.0, 0.0, 0.0), r: 2.0 };
    if let Surface::Cylinder { a: ma, r: mr, .. } = m.surface(&s) {
        assert!(close(ma, m.vector(v(0.0, 0.0, 1.0)), 1e-15) && mr == 2.0);
    } else { panic!("the cylinder changed kind"); }
}

/// Op 7: [1, mode, T 12, N, items [kind, len, data]] -> [1, N, items].
fn request(mode: f64, t: &Rigid, items: &[(f64, Vec<f64>)]) -> Vec<f64> {
    let (m, tr) = t.rows();
    let mut d = vec![1.0, mode];
    d.extend(m.iter().flatten());
    d.extend([tr.x, tr.y, tr.z, items.len() as f64]);
    for (k, data) in items { d.push(*k); d.push(data.len() as f64); d.extend(data); }
    d
}

fn items(reply: &[f64]) -> Vec<(f64, Vec<f64>)> {
    assert_eq!(reply[0], 1.0);
    let mut at = 2;
    let mut out = Vec::new();
    for _ in 0..reply[1] as usize {
        let (k, n) = (reply[at], reply[at + 1] as usize);
        out.push((k, reply[at + 2..at + 2 + n].to_vec()));
        at += 2 + n;
    }
    assert_eq!(at, reply.len());
    out
}

#[test]
fn worker_batch_moves_every_kind_and_fails_whole_on_a_bad_item() {
    let m = tilted();
    let f = [1.0, 2.0, 3.0, 0.6, 0.8, 0.0, -0.8, 0.6, 0.0];
    let mut line = Vec::new();
    put_curve(&mut line, &Curve::Line { a: v(0.0, 0.0, 0.0), b: v(1.0, 2.0, 3.0) });
    let body = write_body(&rat_hole());
    let (rm, rt) = Rigid::translation(v(0.0, 0.0, 1.0)).unwrap().rows();
    let rigid: Vec<f64> = rm.iter().flatten().copied().chain([rt.x, rt.y, rt.z]).collect();
    let plane = vec![0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0];
    let all = vec![(0.0, vec![1.0, 2.0, 3.0]), (1.0, vec![1.0, 2.0, 3.0]), (2.0, f.to_vec()), (3.0, rigid), (4.0, line),
        (5.0, body), (6.0, f.iter().copied().chain([2.0, -1.0]).collect()), (7.0, plane)];
    let reply = items(&run(OP_TRANSFORM, &request(0.0, &m, &all)).unwrap());
    let kinds: Vec<f64> = reply.iter().map(|(k, _)| *k).collect();
    assert_eq!(kinds, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 7.0]);
    let vec3 = |d: &[f64]| v(d[0], d[1], d[2]);
    assert!(close(vec3(&reply[0].1), m.point(p(1.0, 2.0, 3.0)).to_vec3(), 1e-14));
    assert!(close(vec3(&reply[1].1), m.vector(v(1.0, 2.0, 3.0)), 1e-14));
    let frame = Frame::new(p(1.0, 2.0, 3.0), v(0.6, 0.8, 0.0), v(-0.8, 0.6, 0.0)).unwrap();
    assert!(close(vec3(&reply[6].1), m.point(frame.point2(Point2::new(2.0, -1.0))).to_vec3(), 1e-14));
    let moved = read(&reply[5].1);
    assert!((measure(&moved, TOL).unwrap().volume - measure(&rat_hole(), TOL).unwrap().volume).abs() < 1e-9);

    // Mode 1 undoes mode 0.
    let back = items(&run(OP_TRANSFORM, &request(1.0, &m, &[(0.0, reply[0].1.clone())])).unwrap());
    assert!(close(vec3(&back[0].1), v(1.0, 2.0, 3.0), 1e-13));

    // A bad item anywhere fails the batch with its index; a wrong length and a mirrored T are refused.
    let bad = vec![(0.0, vec![1.0, 2.0, 3.0]), (2.0, vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0])];
    let e = run(OP_TRANSFORM, &request(0.0, &m, &bad)).unwrap_err();
    assert!(e.contains("item 1") && e.contains("not unit"), "{e}");
    let e = run(OP_TRANSFORM, &request(0.0, &m, &[(0.0, vec![1.0, 2.0])])).unwrap_err();
    assert!(e.contains("item 0") && e.contains("needs 3"), "{e}");
    let mut mirrored = request(0.0, &m, &[(0.0, vec![1.0, 2.0, 3.0])]);
    for k in 2..11 { mirrored[k] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0][k - 2]; }
    assert!(run(OP_TRANSFORM, &mirrored).unwrap_err().contains("reflection"));
    let mut long = request(0.0, &m, &[(0.0, vec![1.0, 2.0, 3.0])]);
    long[16] = 4.0;
    assert!(run(OP_TRANSFORM, &long).unwrap_err().contains("runs past"));
}

//! Foundation tests (Claude): strict neutral reading, exact extrusions, topology, the worker protocol.

use std::f64::consts::PI;

use gc_kernel::extrude::{extrude, signed_area, volume};
use gc_kernel::geom::{Curve, Vec3};
use gc_kernel::neutral::{put_curve, write_body, write_profile, Profile, Reader};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_CHECK, OP_CUT, OP_EXTRUDE};
use gc_kernel::topology::check;
use gc_kernel::KernelError;

const TOL: f64 = 1e-7;
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn line(a: (f64, f64), b: (f64, f64)) -> Curve { Curve::Line { a: Vec3::new(a.0, a.1, 0.0), b: Vec3::new(b.0, b.1, 0.0) } }

fn circle(cx: f64, cy: f64, r: f64) -> Curve {
    Curve::Arc { c: Vec3::new(cx, cy, 0.0), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r, t0: 0.0, t1: 2.0 * PI }
}

fn rectangle(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Curve> {
    vec![line((x0, y0), (x1, y0)), line((x1, y0), (x1, y1)), line((x1, y1), (x0, y1)), line((x0, y1), (x0, y0))]
}

fn arc_curve(values: &[f64]) -> Result<Curve, KernelError> {
    let mut data = vec![1.0];
    data.extend_from_slice(values);
    Reader::new(&data).curve()
}

// ---- strict reading (Astra's review, 2026-10-27) ----

#[test]
fn valid_arc_is_read() {
    assert!(arc_curve(&[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0, 0.0, PI]).is_ok());
}

#[test]
fn nan_direction_is_refused() {
    assert!(Vec3::new(f64::NAN, 0.0, 1.0).unit().is_err());
    assert!(arc_curve(&[0.0, 0.0, 0.0, f64::NAN, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0, 0.0, PI]).is_err());
}

#[test]
fn huge_count_is_refused() {
    let data = [1.0, 1e100, 0.0, 0.0];
    assert!(matches!(Reader::new(&data).body(), Err(KernelError::InvalidInput(_))));
}

#[test]
fn parallel_arc_frame_is_refused() {
    // normal (0,0,1) and reference direction (0,0,1)
    assert!(matches!(arc_curve(&[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.0, 0.0, PI]), Err(KernelError::InvalidInput(_))));
}

#[test]
fn trailing_and_truncated_data_are_refused() {
    let profile = Profile { loops: vec![rectangle(0.0, 0.0, 1.0, 1.0)], n: Z, low: 0.0, high: 1.0 };
    let body = write_body(&extrude(&profile, TOL).unwrap());
    let mut longer = body.clone();
    longer.push(0.0);
    let mut reader = Reader::new(&longer);
    reader.body().unwrap();
    assert!(reader.finish().is_err());
    assert!(Reader::new(&body[..body.len() - 1]).body().is_err());
}

#[test]
fn bad_edge_reference_is_refused() {
    let profile = Profile { loops: vec![rectangle(0.0, 0.0, 1.0, 1.0)], n: Z, low: 0.0, high: 1.0 };
    let mut body = write_body(&extrude(&profile, TOL).unwrap());
    // The first face's first segment edge index sits after the faces' surface/sense/loop counts; corrupt
    // the last value instead (the last segment's flag) to an out-of-range flag, then an edge index.
    let last = body.len() - 2;
    body[last] = 999.0;
    assert!(matches!(Reader::new(&body).body(), Err(KernelError::InvalidInput(_))));
}

#[test]
fn unsupported_curves_say_so() {
    let data = [2.0, 3.0, 4.0, 8.0];
    assert!(matches!(Reader::new(&data).curve(), Err(KernelError::Unsupported(_))));
}

// ---- exact extrusions ----

#[test]
fn box_from_astras_acceptance_baseline() {
    let profile = Profile { loops: vec![rectangle(0.0, 0.0, 10.0, 8.0)], n: Z, low: 0.0, high: 4.0 };
    let body = extrude(&profile, TOL).unwrap();
    let report = check(&body, TOL).unwrap();
    assert_eq!((report.vertices, report.edges, report.faces), (8, 12, 6));
    assert!((report.genus() - 0.0).abs() < 1e-12);
    assert!((volume(&profile) - 320.0).abs() < 1e-9);
}

#[test]
fn clockwise_input_is_normalized() {
    let mut lp = rectangle(0.0, 0.0, 10.0, 8.0);
    lp = lp.iter().rev().map(Curve::reversed).collect();
    assert!(signed_area(&lp, Z) < 0.0);
    let profile = Profile { loops: vec![lp], n: Z, low: 0.0, high: 4.0 };
    check(&extrude(&profile, TOL).unwrap(), TOL).unwrap();
    assert!((volume(&profile) - 320.0).abs() < 1e-9);
}

#[test]
fn box_with_round_hole_is_genus_one() {
    let profile = Profile { loops: vec![vec![circle(5.0, 4.0, 1.0)], rectangle(0.0, 0.0, 10.0, 8.0)], n: Z, low: 0.0, high: 4.0 };
    let body = extrude(&profile, TOL).unwrap();
    let report = check(&body, TOL).unwrap();
    assert!((report.genus() - 1.0).abs() < 1e-12);
    assert!((volume(&profile) - (320.0 - 4.0 * PI)).abs() < 1e-9);
}

#[test]
fn ring_of_two_circles_uses_each_seam_twice() {
    let profile = Profile { loops: vec![vec![circle(0.0, 0.0, 2.0)], vec![circle(0.0, 0.0, 1.5)]], n: Z, low: -1.0, high: 3.0 };
    let body = extrude(&profile, TOL).unwrap();
    check(&body, TOL).unwrap();
    assert!((volume(&profile) - PI * (4.0 - 2.25) * 4.0).abs() < 1e-9);
}

/// An I-beam section (web horizontal) with concave 0.5 fillets where the web meets the flanges.
fn ibeam_section() -> Vec<Curve> {
    let (w, h, tf, tw, r) = (9.99, 4.0, 0.24, 0.23, 0.5);
    let (web_lo, web_hi) = (-h / 2.0 - tw / 2.0 + h / 2.0 - h / 2.0, 0.0);
    let _ = (web_lo, web_hi);
    let (zb, zt) = (-tw / 2.0, tw / 2.0); // web bottom/top relative to the section centre
    let (yb, yt) = (-h / 2.0, h / 2.0);
    let fillet = |cx: f64, cy: f64, from: f64, to: f64| -> Curve {
        // Clockwise about +Z (concave): counter-clockwise about -Z.
        let n = Vec3::new(0.0, 0.0, -1.0);
        let c = Vec3::new(cx, cy, 0.0);
        let x = Vec3::new(from.cos(), -from.sin(), 0.0);
        Curve::Arc { c, n, x, r, t0: 0.0, t1: (from - to).abs() }
    };
    let _ = fillet;
    // Build by points and arcs, counter-clockwise about +Z.
    let mut lp = Vec::new();
    let p = |x: f64, y: f64| (x, y);
    let pts = [p(0.0, yb), p(tf, yb), p(tf, zb - r)];
    lp.push(line(pts[0], pts[1]));
    lp.push(line(pts[1], pts[2]));
    // fillet from (tf, zb - r) to (tf + r, zb): centre (tf + r, zb - r), concave (clockwise about +Z)
    lp.push(concave(tf + r, zb - r, (tf, zb - r), (tf + r, zb), r));
    lp.push(line((tf + r, zb), (w - tf - r, zb)));
    lp.push(concave(w - tf - r, zb - r, (w - tf - r, zb), (w - tf, zb - r), r));
    lp.push(line((w - tf, zb - r), (w - tf, yb)));
    lp.push(line((w - tf, yb), (w, yb)));
    lp.push(line((w, yb), (w, yt)));
    lp.push(line((w, yt), (w - tf, yt)));
    lp.push(line((w - tf, yt), (w - tf, zt + r)));
    lp.push(concave(w - tf - r, zt + r, (w - tf, zt + r), (w - tf - r, zt), r));
    lp.push(line((w - tf - r, zt), (tf + r, zt)));
    lp.push(concave(tf + r, zt + r, (tf + r, zt), (tf, zt + r), r));
    lp.push(line((tf, zt + r), (tf, yt)));
    lp.push(line((tf, yt), (0.0, yt)));
    lp.push(line((0.0, yt), (0.0, yb)));
    let _ = h;
    lp
}

/// A quarter arc from `a` to `b` about `c`, running clockwise about +Z (counter-clockwise about -Z).
fn concave(cx: f64, cy: f64, a: (f64, f64), b: (f64, f64), r: f64) -> Curve {
    let c = Vec3::new(cx, cy, 0.0);
    let n = Vec3::new(0.0, 0.0, -1.0);
    let x = Vec3::new(a.0 - cx, a.1 - cy, 0.0).scale(1.0 / r);
    let to = Vec3::new(b.0 - cx, b.1 - cy, 0.0).scale(1.0 / r);
    // Angle from x to `to` measured counter-clockwise about n.
    let y = n.cross(x);
    let t1 = to.dot(y).atan2(to.dot(x)).rem_euclid(2.0 * PI);
    Curve::Arc { c, n, x, r, t0: 0.0, t1 }
}

fn sampled_area(lp: &[Curve]) -> f64 {
    let mut pts = Vec::new();
    for c in lp {
        match c {
            Curve::Line { a, .. } => pts.push(*a),
            Curve::Arc { c: center, n, x, r, t0, t1 } => {
                for k in 0..2000 {
                    let t = t0 + (t1 - t0) * k as f64 / 2000.0;
                    pts.push(Curve::point_at_angle(*center, *n, *x, *r, t));
                }
            }
        }
    }
    let mut a = 0.0;
    for i in 0..pts.len() {
        let (p, q) = (pts[i], pts[(i + 1) % pts.len()]);
        a += 0.5 * (p.x * q.y - q.x * p.y);
    }
    a
}

#[test]
fn ibeam_section_with_concave_fillets() {
    let lp = ibeam_section();
    let exact = signed_area(&lp, Z);
    assert!((exact - sampled_area(&lp)).abs() < 1e-5, "exact {exact} sampled {}", sampled_area(&lp));
    let profile = Profile { loops: vec![lp], n: Z, low: 0.0, high: 36.0 };
    let body = extrude(&profile, TOL).unwrap();
    let report = check(&body, TOL).unwrap();
    assert_eq!(report.faces, 2 + 16);
    assert!(body.faces.iter().filter(|f| matches!(f.surface, gc_kernel::geom::Surface::Cylinder { .. }) && !f.sense).count() == 4);
}

#[test]
fn open_profile_is_refused() {
    let lp = vec![line((0.0, 0.0), (1.0, 0.0)), line((1.0, 0.0), (1.0, 1.0)), line((1.0, 1.0), (0.0, 0.9))];
    let profile = Profile { loops: vec![lp], n: Z, low: 0.0, high: 1.0 };
    assert!(matches!(extrude(&profile, TOL), Err(KernelError::InvalidInput(_))));
}

// ---- worker protocol ----

fn reply_status(reply: &[u8]) -> u32 { u32::from_le_bytes(reply[8..12].try_into().unwrap()) }

#[test]
fn protocol_extrude_check_and_cut() {
    let profile = Profile { loops: vec![rectangle(0.0, 0.0, 10.0, 8.0)], n: Z, low: 0.0, high: 4.0 };
    let request = parse_request(&encode_request(OP_EXTRUDE, TOL, &write_profile(&profile))).unwrap();
    let out = handle(&request).unwrap();
    assert!((out.last().unwrap() - 320.0).abs() < 1e-9);
    let body = &out[..out.len() - 1];

    let checked = handle(&parse_request(&encode_request(OP_CHECK, TOL, body)).unwrap()).unwrap();
    assert_eq!(&checked[checked.len() - 4..], &[8.0, 12.0, 6.0, 6.0]);

    // A rectangular through hole (Astra's acceptance table): body and tool back to back.
    let tool = Profile { loops: vec![rectangle(2.0, 2.0, 4.0, 5.0)], n: Z, low: -1.0, high: 5.0 };
    let mut two = body.to_vec();
    two.extend_from_slice(&write_body(&extrude(&tool, TOL).unwrap()));
    let cut = handle(&parse_request(&encode_request(OP_CUT, TOL, &two)).unwrap());
    let reply = gc_kernel::protocol::encode_reply(&cut);
    assert_eq!(reply_status(&reply), 0, "{cut:?}");
    let result = Reader::new(cut.as_ref().unwrap()).body().unwrap();
    check(&result, TOL).unwrap();

    // A round tool is outside the boolean's supported domain and says so.
    let round = Profile { loops: vec![vec![circle(5.0, 4.0, 1.0)]], n: Z, low: -1.0, high: 5.0 };
    let mut two = body.to_vec();
    two.extend_from_slice(&write_body(&extrude(&round, TOL).unwrap()));
    let cut = handle(&parse_request(&encode_request(OP_CUT, TOL, &two)).unwrap());
    assert_eq!(reply_status(&gc_kernel::protocol::encode_reply(&cut)), 1, "{cut:?}");
}

#[test]
fn protocol_refuses_other_units_and_bad_headers() {
    let mut request = encode_request(OP_CHECK, TOL, &[]);
    request[12] = 2; // millimetres
    let parsed = parse_request(&request).unwrap();
    assert!(matches!(handle(&parsed), Err(KernelError::Unsupported(_))));
    assert!(parse_request(b"NOPE").is_err());
    let mut short = encode_request(OP_CHECK, TOL, &[1.0, 2.0]);
    short.pop();
    assert!(parse_request(&short).is_err());
}

#[test]
fn put_curve_round_trips_through_the_reader() {
    let mut data = Vec::new();
    put_curve(&mut data, &circle(1.0, 2.0, 3.0));
    let back = Reader::new(&data).curve().unwrap();
    assert!(back.start().distance(Vec3::new(4.0, 2.0, 0.0)) < 1e-12);
}

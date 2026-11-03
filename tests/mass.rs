//! Exact mass properties (Claude): src/mass.rs against closed-form answers, the extruder's own exact volume,
//! the boolean results and GCmaker's captured I-beam.

use std::f64::consts::PI;

use gc_kernel::boolean;
use gc_kernel::extrude::{extrude, volume};
use gc_kernel::geom::{Curve, Vec3};
use gc_kernel::mass::{measure, face_areas};
use gc_kernel::neutral::{write_body, Body, Profile, Reader};
use gc_kernel::protocol::{handle, parse_request, encode_request, OP_MEASURE};

const TOL: f64 = 1e-7;
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn line(a: (f64, f64), b: (f64, f64)) -> Curve { Curve::Line { a: Vec3::new(a.0, a.1, 0.0), b: Vec3::new(b.0, b.1, 0.0) } }
fn circle(cx: f64, cy: f64, r: f64) -> Curve {
    Curve::Arc { c: Vec3::new(cx, cy, 0.0), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r, t0: 0.0, t1: 2.0 * PI }
}
fn rectangle(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Curve> {
    vec![line((x0, y0), (x1, y0)), line((x1, y0), (x1, y1)), line((x1, y1), (x0, y1)), line((x0, y1), (x0, y0))]
}
fn prism(loops: Vec<Vec<Curve>>, low: f64, high: f64) -> (Profile, Body) {
    let profile = Profile { loops, n: Z, low, high };
    let body = extrude(&profile, TOL).unwrap();
    (profile, body)
}
fn close(a: f64, b: f64, tol: f64) { assert!((a - b).abs() <= tol * b.abs().max(1.0), "{a} vs {b}"); }

/// A capsule: centreline (x0, y) to (x1, y), radius r, counter-clockwise.
fn capsule(x0: f64, x1: f64, y: f64, r: f64) -> Vec<Curve> {
    let arc = |cx: f64, from: f64| Curve::Arc { c: Vec3::new(cx, y, 0.0), n: Z, x: Vec3::new(from.cos(), from.sin(), 0.0), r, t0: 0.0, t1: PI };
    vec![line((x0, y - r), (x1, y - r)), arc(x1, -PI / 2.0), line((x1, y + r), (x0, y + r)), arc(x0, PI / 2.0)]
}

#[test]
fn box_and_holed_prisms_are_exact() {
    let (_, b) = prism(vec![rectangle(0.0, 0.0, 10.0, 8.0)], 0.0, 4.0);
    let m = measure(&b, TOL).unwrap();
    close(m.volume, 320.0, 1e-14);
    close(m.area, 304.0, 1e-14);

    let (p, b) = prism(vec![rectangle(0.0, 0.0, 10.0, 8.0), vec![circle(5.0, 4.0, 1.0)]], 0.0, 4.0);
    let m = measure(&b, TOL).unwrap();
    close(m.volume, volume(&p), 1e-14);
    close(m.volume, 4.0 * (80.0 - PI), 1e-14);
    close(m.area, 2.0 * (80.0 - PI) + 4.0 * 36.0 + 2.0 * PI * 4.0, 1e-14);

    let (p, b) = prism(vec![capsule(3.0, 7.0, 4.0, 1.0)], -1.0, 5.0);
    close(measure(&b, TOL).unwrap().volume, volume(&p), 1e-13);
    close(volume(&p), 6.0 * (8.0 + PI), 1e-13);
}

/// Far from the origin, small: the centroid origin keeps full precision.
#[test]
fn small_part_far_away_keeps_its_precision() {
    let (x, y) = (1.0e5, -2.0e5);
    let (p, b) = prism(vec![rectangle(x, y, x + 0.01, y + 0.02), vec![circle(x + 0.005, y + 0.01, 0.002)]], 3.0e4, 3.0e4 + 0.005);
    close(measure(&b, TOL).unwrap().volume, volume(&p), 1e-8);
}

#[test]
fn inside_out_body_has_negative_volume() {
    let (_, mut b) = prism(vec![rectangle(0.0, 0.0, 10.0, 8.0), vec![circle(5.0, 4.0, 1.0)]], 0.0, 4.0);
    for f in &mut b.faces {
        f.sense = !f.sense;
        for lp in &mut f.loops { lp.reverse(); for s in lp.iter_mut() { s.1 = !s.1; } }
    }
    close(measure(&b, TOL).unwrap().volume, -4.0 * (80.0 - PI), 1e-14);
}

/// The oblique-cut cylinder with a conic top (foundation test body): volume 2 pi; area = bottom pi +
/// ellipse pi sqrt(1 + k^2) + wall (integral of 2 + k cos theta) 4 pi.
#[test]
fn oblique_cut_cylinder_with_a_conic_edge() {
    let k = 0.3f64.tan();
    let m = (1.0 + k * k).sqrt();
    let mut data = vec![1.0, 2.0, 3.0, 3.0, 1.0, 0.0, 0.0, 1.0, 0.0, 2.0 + k];
    data.extend_from_slice(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 2.0 * PI, 0.0, 0.0]);
    data.extend_from_slice(&[3.0, 0.0, 0.0, 2.0, 1.0, 0.0, k, 0.0, 1.0, 0.0, 0.0, 2.0 * PI, 1.0, 1.0]);
    data.extend_from_slice(&[0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 2.0 + k, 0.0, 1.0]);
    data.extend_from_slice(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]);
    data.extend_from_slice(&[0.0, 0.0, 0.0, 2.0, -k / m, 0.0, 1.0 / m, 1.0 / m, 0.0, k / m, 1.0, 1.0, 1.0, 1.0, 1.0]);
    data.extend_from_slice(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4.0, 0.0, 1.0, 2.0, 1.0, 1.0, 0.0, 2.0, 0.0]);
    let body = Reader::new(&data).body().unwrap();
    let props = measure(&body, TOL).unwrap();
    close(props.volume, 2.0 * PI, 1e-14);
    close(props.area, PI + PI * m + 4.0 * PI, 1e-14);
}

/// Astra's cuts, measured independently of the boolean code's own checks.
#[test]
fn boolean_results_measure_as_expected() {
    let (_, block) = prism(vec![rectangle(0.0, 0.0, 10.0, 8.0)], 0.0, 4.0);
    let (_, bore) = prism(vec![vec![circle(5.0, 4.0, 1.0)]], -1.0, 5.0);
    let (_, slot) = prism(vec![capsule(3.0, 7.0, 4.0, 1.0)], -1.0, 5.0);
    let (_, hole) = prism(vec![rectangle(2.0, 2.0, 4.0, 5.0)], -1.0, 5.0);
    close(measure(&boolean::cut(&block, &bore, TOL).unwrap(), TOL).unwrap().volume, 320.0 - 4.0 * PI, 1e-13);
    close(measure(&boolean::cut(&block, &slot, TOL).unwrap(), TOL).unwrap().volume, 320.0 - 4.0 * (8.0 + PI), 1e-13);
    close(measure(&boolean::cut(&block, &hole, TOL).unwrap(), TOL).unwrap().volume, 296.0, 1e-13);
}

fn fixture(name: &str) -> Body {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/");
    let bytes = std::fs::read(format!("{path}{name}")).unwrap();
    let data: Vec<f64> = bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
    let mut r = Reader::new(&data);
    let body = r.body().unwrap();
    r.finish().unwrap();
    body
}

/// GCmaker's captured I-beam (planes, concave fillet cylinders, loops as Eyeshot had them): the exact volume is the
/// section area times 36; Eyeshot's mesh (manifest) reads 155.589939, a little high: on concave fillets the chords
/// add material.
#[test]
fn captured_ibeam_base_measures_exactly() {
    let base = fixture("rat_hole/base.f64");
    let m = measure(&base, 1e-6).unwrap();
    let (tf, tw, r, w, h) = (0.24, 0.23, 0.5, 9.99, 4.0);
    let section = 2.0 * tf * h + (w - 2.0 * tf) * tw + 4.0 * r * r * (1.0 - PI / 4.0);
    close(m.volume, section * 36.0, 1e-12);
    assert!(m.volume < 155.589939 && m.volume > 155.589939 * (1.0 - 1e-5), "{}", m.volume);
}

#[test]
fn measure_op_answers_volume_and_area() {
    let (_, b) = prism(vec![rectangle(0.0, 0.0, 10.0, 8.0)], 0.0, 4.0);
    let reply = handle(&parse_request(&encode_request(OP_MEASURE, TOL, &write_body(&b))).unwrap()).unwrap();
    assert_eq!(reply.len(), 3);
    close(reply[0], 320.0, 1e-14);
    close(reply[1], 304.0, 1e-14);
    assert_eq!(reply[2], 1.0);   // closed form
}

#[test]
fn selected_face_areas_include_holes_and_curved_walls_in_body_order() {
    let (_, body) = prism(vec![rectangle(0., 0., 10., 8.), vec![circle(5., 4., 1.)]], 0., 4.);
    let faces = face_areas(&body, TOL).unwrap();
    assert_eq!(faces.len(), body.faces.len());
    for (face, result) in body.faces.iter().zip(&faces) {
        assert!(result.exact);
        match face.surface {
            gc_kernel::geom::Surface::Plane { n, .. } if n.dot(Z).abs() > 0.9 => close(result.area, 80. - PI, 1e-13),
            gc_kernel::geom::Surface::Cylinder { .. } => close(result.area, 8. * PI, 1e-13),
            _ => assert!(result.area == 32. || result.area == 40.),
        }
    }
    close(faces.iter().map(|f| f.area).sum(), measure(&body, TOL).unwrap().area, 1e-14);
    let request = encode_request(gc_kernel::protocol::OP_FACE_AREAS, TOL, &write_body(&body));
    let reply = handle(&parse_request(&request).unwrap()).unwrap();
    assert_eq!(reply.len(), 2 + 2 * faces.len());
    assert_eq!(&reply[..2], &[1., faces.len() as f64]);
    for (i, f) in faces.iter().enumerate() { assert_eq!(&reply[2+2*i..4+2*i], &[f.area, 1.]); }
    let mut bad = body.clone(); bad.faces.pop();
    assert!(handle(&parse_request(&encode_request(gc_kernel::protocol::OP_FACE_AREAS, TOL, &write_body(&bad))).unwrap()).is_err());
}

#[test]
fn elliptical_wall_reports_quadrature_separately_from_its_closed_form_caps() {
    let ellipse = Curve::Conic { c: Vec3::new(0.,0.,0.), u: Vec3::new(2.,0.,0.), v: Vec3::new(0.,1.,0.), t0: 0., t1: 2.*PI };
    let (_, body) = prism(vec![vec![ellipse]], 0., 3.);
    let faces = face_areas(&body, TOL).unwrap();
    let mut walls = 0;
    for (face, area) in body.faces.iter().zip(&faces) {
        match face.surface {
            gc_kernel::geom::Surface::Plane { .. } => { assert!(area.exact); close(area.area, 2.*PI, 1e-12); },
            gc_kernel::geom::Surface::EllipticCylinder { .. } => {
                walls += 1; assert!(!area.exact); close(area.area, 3.*9.688448220547675, 1e-12);
            },
            _ => panic!("unexpected surface"),
        }
    }
    assert_eq!(walls, 1);
    close(faces.iter().map(|a|a.area).sum(), measure(&body, TOL).unwrap().area, 1e-14);
}

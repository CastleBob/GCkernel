use super::*;
use crate::neutral::Profile;
const TOL: f64 = 1e-7;
fn p(x: f64, y: f64) -> Vec3 {
    Vec3::new(x, y, 0.)
}
fn line(a: (f64, f64), b: (f64, f64)) -> Curve {
    Curve::Line {
        a: p(a.0, a.1),
        b: p(b.0, b.1),
    }
}
fn quarter(c: (f64, f64), from: (f64, f64)) -> Curve {
    Curve::Arc {
        c: p(c.0, c.1),
        n: Vec3::new(0., 0., -1.),
        x: p(from.0 - c.0, from.1 - c.1),
        r: 1.,
        t0: 0.,
        t1: PI / 2.,
    }
}
fn ibeam() -> Body {
    // Two 1x6 flanges, 8x1 web, four concave unit-radius root fillets.
    let outline = vec![
        line((0., -3.), (1., -3.)),
        line((1., -3.), (1., -1.5)),
        quarter((2., -1.5), (1., -1.5)),
        line((2., -0.5), (8., -0.5)),
        quarter((8., -1.5), (8., -0.5)),
        line((9., -1.5), (9., -3.)),
        line((9., -3.), (10., -3.)),
        line((10., -3.), (10., 3.)),
        line((10., 3.), (9., 3.)),
        line((9., 3.), (9., 1.5)),
        quarter((8., 1.5), (9., 1.5)),
        line((8., 0.5), (2., 0.5)),
        quarter((2., 1.5), (2., 0.5)),
        line((1., 1.5), (1., 3.)),
        line((1., 3.), (0., 3.)),
        line((0., 3.), (0., -3.)),
    ];
    crate::extrude::extrude(
        &Profile {
            loops: vec![outline],
            n: Vec3::new(0., 0., 1.),
            low: 0.,
            high: 20.,
        },
        TOL,
    )
    .unwrap()
}
fn bore(x: f64, z: f64, r: f64, low: f64, high: f64) -> Body {
    crate::extrude::extrude(
        &Profile {
            loops: vec![vec![Curve::Arc {
                c: Vec3::new(x, 0., z),
                n: Vec3::new(0., 1., 0.),
                x: Vec3::new(1., 0., 0.),
                r,
                t0: 0.,
                t1: TAU,
            }]],
            n: Vec3::new(0., 1., 0.),
            low,
            high,
        },
        TOL,
    )
    .unwrap()
}

fn volume(body: &Body) -> f64 {
    let mut sum = 0.;
    for face in &body.faces {
        match face.surface {
            Surface::EllipticCylinder { .. } | Surface::Cone { .. } => unreachable!("no elliptic-cylinder or conical faces in these test bodies"),
            Surface::Plane { p, n, .. } => {
                let outward = n.scale(sign(face.sense));
                let area: f64 = face
                    .loops
                    .iter()
                    .map(|lp| loop_area(body, lp, outward))
                    .sum();
                sum += p.dot(outward) * area / 3.;
            }
            Surface::Cylinder { a, r, .. } => {
                let vertices: Vec<Vec3> = face
                    .loops
                    .iter()
                    .flatten()
                    .flat_map(|(i, _)| {
                        let e = &body.edges[*i];
                        [body.vertices[e.start], body.vertices[e.end]]
                    })
                    .collect();
                let lo = vertices
                    .iter()
                    .map(|p| p.dot(a))
                    .fold(f64::INFINITY, f64::min);
                let hi = vertices
                    .iter()
                    .map(|p| p.dot(a))
                    .fold(f64::NEG_INFINITY, f64::max);
                let arc = face
                    .loops
                    .iter()
                    .flatten()
                    .find_map(|(i, _)| match body.edges[*i].curve {
                        Curve::Arc {
                            c, n, x, t0, t1, ..
                        } => Some((c, n, x, t0, t1)),
                        _ => None,
                    })
                    .unwrap();
                let (c, n, x, t0, t1) = arc;
                sum += sign(face.sense)
                    * r
                    * (hi - lo)
                    * (c.dot(x) * (t1.sin() - t0.sin())
                        + c.dot(n.cross(x)) * (t0.cos() - t1.cos())
                        + r * (t1 - t0))
                    / 3.;
            }
        }
    }
    sum
}

#[test]
fn web_bore_preserves_all_four_root_fillets_and_removes_exact_material() {
    let body = ibeam();
    let tool = bore(5., 10., 1., -4., 4.);
    let original = crate::neutral::write_body(&body);
    let result = crate::boolean::cut(&body, &tool, TOL).unwrap();
    assert_eq!(original, crate::neutral::write_body(&body));
    assert_eq!(crate::topology::check(&result, TOL).unwrap().genus(), 1.);
    assert!((volume(&body) - (480. - 20. * PI)).abs() < 1e-9);
    assert!((volume(&result) - (480. - 21. * PI)).abs() < 1e-9);
    for (i, face) in body.faces.iter().enumerate() {
        if matches!(face.surface, Surface::Cylinder { .. }) {
            assert_eq!(format!("{face:?}"), format!("{:?}", result.faces[i]));
        }
    }
    assert_eq!(result.vertices[..body.vertices.len()], body.vertices);
}

#[test]
fn separated_web_bores_accumulate_without_changing_previous_geometry() {
    let mut body = ibeam();
    for (i, z) in [2., 5., 8., 11., 14., 17.].into_iter().enumerate() {
        body = crate::boolean::cut(&body, &bore(5., z, 0.4, -4., 4.), TOL).unwrap();
        assert_eq!(
            crate::topology::check(&body, TOL).unwrap().genus(),
            (i + 1) as f64
        );
        assert!((volume(&body) - (480. - 20. * PI - (i + 1) as f64 * PI * 0.16)).abs() < 1e-8);
    }
}

#[test]
fn larger_fillet_crossing_succeeds_but_blind_and_interacting_bores_remain_explicit() {
    let crossing = crate::boolean::cut(&ibeam(), &bore(1.5, 10., 1., -4., 4.), TOL).unwrap();
    assert_eq!(crate::topology::check(&crossing, TOL).unwrap().genus(), 1.);
    for tool in [bore(5., 0., 1., -4., 4.), bore(5., 10., 1., 0., 4.)] {
        let result = crate::boolean::cut(&ibeam(), &tool, TOL);
        assert!(
            matches!(result, Err(KernelError::Unsupported(_))),
            "{:?}",
            result.as_ref().map(|b| b.faces.len())
        );
    }
    let body = crate::boolean::cut(&ibeam(), &bore(5., 10., 1., -4., 4.), TOL).unwrap();
    for z in [10., 11., 12.] {
        assert!(matches!(
            crate::boolean::cut(&body, &bore(5., z, 1., -4., 4.), TOL),
            Err(KernelError::Unsupported(_))
        ));
    }
}

#[test]
fn malformed_unaffected_face_is_not_carried_into_a_new_solid() {
    let mut body = ibeam();
    let f = body
        .faces
        .iter_mut()
        .find(|f| matches!(f.surface, Surface::Cylinder { .. }))
        .unwrap();
    f.sense = !f.sense;
    assert!(matches!(
        crate::boolean::cut(&body, &bore(5., 10., 1., -4., 4.), TOL),
        Err(KernelError::InvalidInput(_))
    ));
    let mut body = ibeam();
    body.edges[0].start = 9999;
    assert!(matches!(
        crate::boolean::cut(&body, &bore(5., 10., 1., -4., 4.), TOL),
        Err(KernelError::InvalidInput(_))
    ));
}

#[test]
fn web_bore_survives_the_wire_contract() {
    let mut data = crate::neutral::write_body(&ibeam());
    data.extend(crate::neutral::write_body(&bore(5., 10., 1., -4., 4.)));
    let request = crate::protocol::parse_request(&crate::protocol::encode_request(
        crate::protocol::OP_CUT,
        TOL,
        &data,
    ))
    .unwrap();
    let reply = crate::protocol::handle(&request).unwrap();
    let result = crate::neutral::Reader::new(&reply).body().unwrap();
    assert!((volume(&result) - (480. - 21. * PI)).abs() < 1e-9);
}

#[test]
fn captured_gcmaker_ibeam_admits_a_clear_web_hole() {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/rat_hole/base.f64");
    let bytes = std::fs::read(path).expect("Claude's captured GCmaker I-beam fixture");
    let data: Vec<f64> = bytes
        .chunks_exact(8)
        .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
        .collect();
    let body = crate::neutral::Reader::new(&data).body().unwrap();
    let n = Vec3::new(0., 0., 1.);
    let tool = crate::extrude::extrude(
        &Profile {
            loops: vec![vec![Curve::Arc {
                c: Vec3::new(5., 18., 0.),
                n,
                x: Vec3::new(1., 0., 0.),
                r: 0.5,
                t0: 0.,
                t1: TAU,
            }]],
            n,
            low: -6.,
            high: 1.,
        },
        TOL,
    )
    .unwrap();
    let result = crate::boolean::cut(&body, &tool, TOL).unwrap();
    assert_eq!(crate::topology::check(&result, TOL).unwrap().genus(), 1.);
    assert!((volume(&body) - volume(&result) - PI * 0.25 * 0.23).abs() < 1e-9);
    assert_eq!(
        result
            .faces
            .iter()
            .filter(|f| matches!(f.surface, Surface::Cylinder { .. }))
            .count(),
        5
    );
}

fn permute(body: &Body, shift: usize) -> Body {
    let map = |p: Vec3| {
        let a = xyz(p);
        vector(std::array::from_fn(|k| a[(k + shift) % 3]))
    };
    let mut result = body.clone();
    for p in &mut result.vertices {
        *p = map(*p);
    }
    for e in &mut result.edges {
        match &mut e.curve {
            Curve::Conic { c, u, v, .. } => {
                *c = map(*c);
                *u = map(*u);
                *v = map(*v);
            }
            Curve::CylinderPair { .. } => panic!("pair edges are outside these passage fixtures"),
            Curve::Line { a, b } => {
                *a = map(*a);
                *b = map(*b);
            }
            Curve::Arc { c, n, x, .. } => {
                *c = map(*c);
                *n = map(*n);
                *x = map(*x);
            }
        }
    }
    for f in &mut result.faces {
        match &mut f.surface {
            Surface::EllipticCylinder { .. } | Surface::Cone { .. } => unreachable!("no elliptic-cylinder or conical faces in these test bodies"),
            Surface::Plane { p, n, x } => {
                *p = map(*p);
                *n = map(*n);
                *x = map(*x);
            }
            Surface::Cylinder { p, a, x, .. } => {
                *p = map(*p);
                *a = map(*a);
                *x = map(*x);
            }
        }
    }
    result
}

#[test]
fn web_passage_works_on_all_coordinate_axes() {
    for shift in 0..3 {
        let body = permute(&ibeam(), shift);
        let tool = permute(&bore(5., 10., 1., -4., 4.), shift);
        let result = crate::boolean::cut(&body, &tool, TOL).unwrap();
        assert!((volume(&result) - (480. - 21. * PI)).abs() < 1e-9);
    }
}

#[test]
fn cap_union_cancels_exact_seams_but_never_merges_near_planes() {
    let points = [(0., 0.), (2., 0.), (2., 2.), (0., 2.), (0., 1.)];
    let profile = Profile {
        loops: vec![(0..points.len())
            .map(|i| line(points[i], points[(i + 1) % points.len()]))
            .collect()],
        n: Vec3::new(0., 0., 1.),
        low: 0.,
        high: 3.,
    };
    let body = crate::extrude::extrude(&profile, TOL).unwrap();
    let merged = coalesce_caps(&body, 0, TOL).unwrap();
    assert_eq!(merged.faces.len(), body.faces.len() - 1);
    assert!((crate::mass::measure(&merged, TOL).unwrap().volume - 12.).abs() < 1e-12);
    let mut near_body = body.clone();
    let face = near_body
        .faces
        .iter_mut()
        .find(|f| matches!(f.surface,Surface::Plane {p,n,..} if p.x==0. && n.x.abs()==1.))
        .unwrap();
    if let Surface::Plane { p, .. } = &mut face.surface {
        p.x += 0.25 * TOL;
    }
    let unchanged = coalesce_caps(&near_body, 0, TOL).unwrap();
    assert_eq!(
        crate::neutral::write_body(&near_body),
        crate::neutral::write_body(&unchanged)
    );
    let mut malformed = body.clone();
    malformed.faces[0].sense = !malformed.faces[0].sense;
    assert!(coalesce_caps(&malformed, 0, TOL).is_err());
}

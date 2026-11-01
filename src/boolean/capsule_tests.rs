use super::*;

const TOL: f64 = 1e-7;

// Integrate the returned analytic boundary, independent of the cut recipe.
fn boundary_volume(body: &Body) -> f64 {
    let mut volume = 0.;
    for face in &body.faces {
        let sign = if face.sense { 1. } else { -1. };
        match face.surface {
            Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => unreachable!("no elliptic-cylinder or revolved faces in these test bodies"),
            Surface::Plane { p, n, .. } => {
                let outward = n.scale(sign);
                let mut area = 0.;
                for (i, forward) in face.loops.iter().flatten() {
                    let curve = &body.edges[*i].curve;
                    let integral = match *curve {
                        Curve::CylinderPair { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } | Curve::Bezier { .. } | Curve::RationalBezier { .. } => {
                            panic!("this independent prism integral excludes cylinder pairs")
                        }
                        Curve::Conic { c, u, v, t0, t1 } => c
                            .cross(curve.end().sub(curve.start()))
                            .add(u.cross(v).scale(t1 - t0)),
                        Curve::Line { a, b } => a.cross(b),
                        Curve::Arc {
                            c, n, r, t0, t1, ..
                        } => c
                            .cross(curve.end().sub(curve.start()))
                            .add(n.scale(r * r * (t1 - t0))),
                    };
                    area += 0.5 * integral.dot(outward) * if *forward { 1. } else { -1. };
                }
                assert!(area > 0.);
                volume += p.dot(outward) * area / 3.;
            }
            Surface::Cylinder { a, r, .. } => {
                let heights: Vec<f64> = face
                    .loops
                    .iter()
                    .flatten()
                    .flat_map(|(i, _)| {
                        let e = &body.edges[*i];
                        [body.vertices[e.start].dot(a), body.vertices[e.end].dot(a)]
                    })
                    .collect();
                let height = heights.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                    - heights.iter().copied().fold(f64::INFINITY, f64::min);
                let arc = face
                    .loops
                    .iter()
                    .flatten()
                    .find_map(|(i, _)| {
                        matches!(body.edges[*i].curve, Curve::Arc { .. })
                            .then_some(&body.edges[*i].curve)
                    })
                    .unwrap();
                let Curve::Arc {
                    c, n, x, t0, t1, ..
                } = *arc
                else {
                    unreachable!()
                };
                let y = n.cross(x);
                volume += sign
                    * r
                    * height
                    * (c.dot(x) * (t1.sin() - t0.sin())
                        + c.dot(y) * (t0.cos() - t1.cos())
                        + r * (t1 - t0))
                    / 3.;
            }
        }
    }
    volume
}

fn fixture(axis: usize, angle: f64, offset: f64, depth: (f64, f64)) -> (Body, Body) {
    let u = (axis + 1) % 3;
    let v = (axis + 2) % 3;
    let point = |s: f64, t: f64| {
        let mut p = [0.; 3];
        p[u] = s;
        p[v] = t;
        vec(p)
    };
    let mut n = [0.; 3];
    n[axis] = 1.;
    let n = vec(n);
    let corners = [point(0., 0.), point(10., 0.), point(10., 8.), point(0., 8.)];
    let outline = (0..4)
        .map(|i| Curve::Line {
            a: corners[i],
            b: corners[(i + 1) % 4],
        })
        .collect();
    let body = crate::extrude::extrude(
        &Profile {
            loops: vec![outline],
            n,
            low: 0.,
            high: 4.,
        },
        TOL,
    )
    .unwrap();
    // Independent 2D fixture with two diameter-length straight sides and arcs.
    let transform = |s: f64, t: f64| {
        point(
            5. + offset + s * angle.cos() - t * angle.sin(),
            4. + s * angle.sin() + t * angle.cos(),
        )
    };
    let radial = |s: f64, t: f64| {
        point(
            s * angle.cos() - t * angle.sin(),
            s * angle.sin() + t * angle.cos(),
        )
    };
    let slot = vec![
        Curve::Line {
            a: transform(-2., -1.),
            b: transform(2., -1.),
        },
        Curve::Arc {
            c: transform(2., 0.),
            n,
            x: radial(0., -1.),
            r: 1.,
            t0: 0.,
            t1: PI,
        },
        Curve::Line {
            a: transform(2., 1.),
            b: transform(-2., 1.),
        },
        Curve::Arc {
            c: transform(-2., 0.),
            n,
            x: radial(0., 1.),
            r: 1.,
            t0: 0.,
            t1: PI,
        },
    ];
    let tool = crate::extrude::extrude(
        &Profile {
            loops: vec![slot],
            n,
            low: depth.0,
            high: depth.1,
        },
        TOL,
    )
    .unwrap();
    (body, tool)
}

#[test]
fn slot_retains_analytic_geometry_and_correct_volume_on_every_axis() {
    for axis in 0..3 {
        for angle in [0., 0.4, 1.2] {
            let (body, tool) = fixture(axis, angle, 0., (-1., 5.));
            let original = crate::neutral::write_body(&body);
            let tool_original = crate::neutral::write_body(&tool);
            let result = crate::boolean::cut(&body, &tool, TOL).unwrap();
            assert_eq!(crate::neutral::write_body(&body), original);
            assert_eq!(crate::neutral::write_body(&tool), tool_original);
            assert_eq!(crate::topology::check(&result, TOL).unwrap().genus(), 1.);
            assert_eq!(result.faces.len(), 10);
            assert_eq!(
                result.faces.iter().filter(|f| f.loops.len() == 2).count(),
                2
            );
            assert_eq!(
                result
                    .faces
                    .iter()
                    .filter(|f| matches!(f.surface,Surface::Cylinder {r,..} if r==1.) && !f.sense)
                    .count(),
                2
            );
            assert!((boundary_volume(&result) - (320. - 4. * (8. + PI))).abs() < 1e-9);
            let data = crate::neutral::write_body(&result);
            let mut reader = crate::neutral::Reader::new(&data);
            let back = reader.body().unwrap();
            reader.finish().unwrap();
            assert!((boundary_volume(&back) - (320. - 4. * (8. + PI))).abs() < 1e-9);
        }
    }
}

#[test]
fn indexing_curve_direction_and_loop_start_do_not_change_the_slot() {
    let (body, mut tool) = fixture(2, 0.4, 0., (-1., 5.));
    for (i, e) in tool.edges.iter_mut().enumerate() {
        if i % 2 == 0 {
            e.curve = e.curve.reversed();
            std::mem::swap(&mut e.start, &mut e.end);
            for face in &mut tool.faces {
                for (index, forward) in face.loops.iter_mut().flatten() {
                    if *index == i {
                        *forward = !*forward;
                    }
                }
            }
        }
    }
    tool.vertices.reverse();
    for e in &mut tool.edges {
        e.start = 7 - e.start;
        e.end = 7 - e.end;
    }
    tool.edges.reverse();
    for face in &mut tool.faces {
        for (i, _) in face.loops.iter_mut().flatten() {
            *i = 11 - *i;
        }
        face.loops[0].rotate_left(1);
    }
    tool.faces.reverse();
    let result = crate::boolean::cut(&body, &tool, TOL).unwrap();
    assert!((boundary_volume(&result) - (320. - 4. * (8. + PI))).abs() < 1e-9);
}

#[test]
fn slot_crossing_cuts_in_the_section_domain_tangent_and_blind_stay_unsupported() {
    // Since 2026-10-29 a slot reaching the outline is a prism section cut. At offset 3 its right half circle is
    // past x = 10: inside remain the straight part x 6..10 (8) and the left half disc (pi / 2). At offset 2 its end
    // touches x = 10 at a point (non-manifold): refused. A blind slot is a pocket since 2026-10-31 (pocket.rs): 2 of
    // its depth inside, (8 + pi) 2 removed.
    let (body, tool) = fixture(2, 0., 3., (-1., 5.));
    let result = crate::boolean::cut(&body, &tool, TOL).unwrap();
    crate::topology::check(&result, TOL).unwrap();
    assert!((boundary_volume(&result) - (320. - 4. * (8. + PI / 2.))).abs() < 1e-9);
    let (body, tool) = fixture(2, 0., 2., (-1., 5.));
    assert!(matches!(crate::boolean::cut(&body, &tool, TOL), Err(KernelError::Unsupported(_))));
    let (body, tool) = fixture(2, 0., 0., (2., 5.));
    let result = crate::boolean::cut(&body, &tool, TOL).unwrap();
    crate::topology::check(&result, TOL).unwrap();
    assert!((boundary_volume(&result) - (320. - 2. * (8. + PI))).abs() < 1e-9);
    let (body, tool) = fixture(2, 0., 0., (0., 4.));
    let result = crate::boolean::cut(&body, &tool, TOL).unwrap();
    assert!((boundary_volume(&result) - (320. - 4. * (8. + PI))).abs() < 1e-9);
    // Recutting the same slot removes nothing more.
    let again = crate::boolean::cut(&result, &tool, TOL).unwrap();
    assert!((boundary_volume(&again) - (320. - 4. * (8. + PI))).abs() < 1e-9);
}

#[test]
fn malformed_caps_and_unequal_end_radii_are_refused() {
    let (body, mut tool) = fixture(2, 0., 0., (-1., 5.));
    let arc = tool
        .edges
        .iter_mut()
        .find(|e| matches!(e.curve, Curve::Arc { .. }))
        .unwrap();
    if let Curve::Arc { r, .. } = &mut arc.curve {
        *r = 1.2;
    }
    assert!(crate::boolean::cut(&body, &tool, TOL).is_err());
    let (_, mut tool) = fixture(2, 0., 0., (-1., 5.));
    let face = tool
        .faces
        .iter_mut()
        .find(|f| matches!(f.surface, Surface::Cylinder { .. }))
        .unwrap();
    if let Surface::Cylinder { r, .. } = &mut face.surface {
        *r = 1.2;
    }
    assert!(matches!(
        crate::boolean::cut(&body, &tool, TOL),
        Err(KernelError::Unsupported(_))
    ));
    let (_, mut tool) = fixture(2, 0., 0., (-1., 5.));
    tool.faces[0].sense = !tool.faces[0].sense;
    assert!(crate::boolean::cut(&body, &tool, TOL).is_err());
    let (_, mut tool) = fixture(2, 0., 0., (-1., 5.));
    tool.edges[0].start = 100;
    assert!(matches!(
        crate::boolean::cut(&body, &tool, TOL),
        Err(KernelError::InvalidInput(_))
    ));
}

#[test]
fn capsule_cut_succeeds_through_protocol() {
    let (body, tool) = fixture(2, 0., 0., (-1., 5.));
    let mut data = crate::neutral::write_body(&body);
    data.extend(crate::neutral::write_body(&tool));
    let request = crate::protocol::parse_request(&crate::protocol::encode_request(
        crate::protocol::OP_CUT,
        TOL,
        &data,
    ))
    .unwrap();
    let response = crate::protocol::handle(&request).unwrap();
    let result = crate::neutral::Reader::new(&response).body().unwrap();
    assert_eq!(crate::topology::check(&result, TOL).unwrap().genus(), 1.);
    assert!((boundary_volume(&result) - (320. - 4. * (8. + PI))).abs() < 1e-9);
}

fn captured_ibeam() -> Body {
    captured_rat_hole("base")
}
fn captured_rat_hole(name: &str) -> Body {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(format!("tests/fixtures/rat_hole/{name}.f64"));
    let bytes = std::fs::read(path).unwrap();
    let data: Vec<f64> = bytes
        .chunks_exact(8)
        .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
        .collect();
    crate::neutral::Reader::new(&data).body().unwrap()
}

#[test]
fn web_slot_after_both_rat_holes_retains_previous_conic_boundaries() {
    let mut base = captured_ibeam();
    for step in ["step01", "step02"] {
        base = crate::boolean::cut(&base, &captured_rat_hole(step), 1e-6).unwrap();
        base = crate::neutral::Reader::new(&crate::neutral::write_body(&base))
            .body()
            .unwrap();
    }
    assert!(base
        .edges
        .iter()
        .any(|e| matches!(e.curve, Curve::Conic { .. })));
    let saved = crate::neutral::write_body(&base);
    let (_, tool) = fixture(2, 0., 0., (-6., 1.));
    let result = crate::boolean::cut(&base, &tool, 1e-6).unwrap();
    assert_eq!(crate::topology::check(&result, 1e-6).unwrap().genus(), 1.);
    let expected = 153.57573395171957 - 0.23 * (8. + PI);
    assert!((crate::mass::measure(&result, 1e-6).unwrap().volume - expected).abs() < 1e-8);
    let n = Vec3::new(0., 0., 1.);
    let bore = crate::extrude::extrude(
        &Profile {
            loops: vec![vec![Curve::Arc {
                c: Vec3::new(5., 14., 0.),
                n,
                x: Vec3::new(1., 0., 0.),
                r: 0.5,
                t0: 0.,
                t1: 2. * PI,
            }]],
            n,
            low: -6.,
            high: 1.,
        },
        1e-6,
    )
    .unwrap();
    let two = crate::boolean::cut(&result, &bore, 1e-6).unwrap();
    assert_eq!(crate::topology::check(&two, 1e-6).unwrap().genus(), 2.);
    assert!(
        (crate::mass::measure(&two, 1e-6).unwrap().volume - (expected - 0.23 * PI * 0.25)).abs()
            < 1e-8
    );
    for (i, edge) in base.edges.iter().enumerate() {
        assert_eq!(format!("{edge:?}"), format!("{:?}", result.edges[i]));
    }
    assert_eq!(saved, crate::neutral::write_body(&base));
    let mut broken = base.clone();
    let edge = broken
        .edges
        .iter_mut()
        .find(|e| matches!(e.curve, Curve::Conic { .. }))
        .unwrap();
    if let Curve::Conic { u, .. } = &mut edge.curve {
        *u = u.scale(1.1);
    }
    assert!(crate::boolean::cut(&broken, &tool, 1e-6).is_err());
}
fn translate(body: &Body, offset: Vec3) -> Body {
    let mut out = body.clone();
    for v in &mut out.vertices {
        *v = v.add(offset);
    }
    for e in &mut out.edges {
        e.curve = e.curve.translated(offset);
    }
    for f in &mut out.faces {
        match &mut f.surface {
            Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => unreachable!("no elliptic-cylinder or revolved faces in these test bodies"),
            Surface::Plane { p, .. } | Surface::Cylinder { p, .. } => *p = p.add(offset),
        }
    }
    out
}
#[test]
fn real_ibeam_supports_separated_slot_bore_slot_sequence() {
    let base = captured_ibeam();
    let original = crate::neutral::write_body(&base);
    let (_, tool) = fixture(2, 0., 0., (-6., 1.));
    let one = crate::boolean::cut(&base, &tool, TOL).unwrap();
    let wire = crate::neutral::write_body(&one);
    let one = crate::neutral::Reader::new(&wire).body().unwrap();
    assert_eq!(crate::topology::check(&one, TOL).unwrap().genus(), 1.);
    assert!((boundary_volume(&base) - boundary_volume(&one) - 0.23 * (8. + PI)).abs() < 1e-9);
    let n = Vec3::new(0., 0., 1.);
    let circle = crate::extrude::extrude(
        &Profile {
            loops: vec![vec![Curve::Arc {
                c: Vec3::new(5., 14., 0.),
                n,
                x: Vec3::new(1., 0., 0.),
                r: 0.5,
                t0: 0.,
                t1: std::f64::consts::TAU,
            }]],
            n,
            low: -6.,
            high: 1.,
        },
        TOL,
    )
    .unwrap();
    let two = crate::boolean::cut(&one, &circle, TOL).unwrap();
    let three = crate::boolean::cut(&two, &translate(&tool, Vec3::new(0., 20., 0.)), TOL).unwrap();
    assert_eq!(crate::topology::check(&three, TOL).unwrap().genus(), 3.);
    assert!(
        (boundary_volume(&base) - boundary_volume(&three) - 0.23 * (2. * (8. + PI) + PI * 0.25))
            .abs()
            < 1e-8
    );
    assert_eq!(crate::neutral::write_body(&base), original);
    for (i, face) in base.faces.iter().enumerate() {
        if matches!(face.surface, Surface::Cylinder { .. }) {
            assert_eq!(format!("{face:?}"), format!("{:?}", three.faces[i]));
        }
    }
}

#[test]
fn separated_web_slots_commute_in_either_order() {
    let base = captured_ibeam();
    let (_, tool) = fixture(2, 0., 0., (-6., 1.));
    let other = translate(&tool, Vec3::new(0., 20., 0.));
    let ab = crate::boolean::cut(
        &crate::boolean::cut(&base, &tool, TOL).unwrap(),
        &other,
        TOL,
    )
    .unwrap();
    let ba = crate::boolean::cut(
        &crate::boolean::cut(&base, &other, TOL).unwrap(),
        &tool,
        TOL,
    )
    .unwrap();
    assert!((boundary_volume(&ab) - boundary_volume(&ba)).abs() < 1e-9);
    assert_eq!(crate::topology::check(&ab, TOL).unwrap().genus(), 2.);
    assert_eq!(crate::topology::check(&ba, TOL).unwrap().genus(), 2.);
}

#[test]
fn web_slot_clearance_and_existing_hole_boundaries_are_enforced() {
    let base = captured_ibeam();
    let (_, tool) = fixture(2, 0., 0., (-6., 1.));
    let one = crate::boolean::cut(&base, &tool, TOL).unwrap();
    for offset in [
        Vec3::new(-2., 0., 0.),
        Vec3::new(0., 0., 0.),
        Vec3::new(0., 1., 0.),
    ] {
        assert!(matches!(
            crate::boolean::cut(&one, &translate(&tool, offset), TOL),
            Err(KernelError::Unsupported(_))
        ));
    }
}

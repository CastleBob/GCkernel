use super::*;

const TOL: f64 = 1e-7;

// Integrate the returned analytic boundary, independent of the cut recipe.
fn boundary_volume(body: &Body) -> f64 {
    let mut volume = 0.;
    for face in &body.faces {
        let sign = if face.sense { 1. } else { -1. };
        match face.surface {
            Surface::Plane { p, n, .. } => {
                let outward = n.scale(sign);
                let mut area = 0.;
                for (i, forward) in face.loops.iter().flatten() {
                    let curve = &body.edges[*i].curve;
                    let integral = match *curve {
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
fn slot_crossing_tangent_blind_and_repeated_cases_are_explicitly_unsupported() {
    for (offset, depth) in [(3., (-1., 5.)), (2., (-1., 5.)), (0., (2., 5.))] {
        let (body, tool) = fixture(2, 0., offset, depth);
        assert!(matches!(
            crate::boolean::cut(&body, &tool, TOL),
            Err(KernelError::Unsupported(_))
        ));
    }
    let (body, tool) = fixture(2, 0., 0., (0., 4.));
    let result = crate::boolean::cut(&body, &tool, TOL).unwrap();
    assert!((boundary_volume(&result) - (320. - 4. * (8. + PI))).abs() < 1e-9);
    assert!(matches!(
        crate::boolean::cut(&result, &tool, TOL),
        Err(KernelError::Unsupported(_))
    ));
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

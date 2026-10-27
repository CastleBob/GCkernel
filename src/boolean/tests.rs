use super::*;

fn box_body(low: [f64; 3], high: [f64; 3]) -> Body {
    let vertices = vec![
        vec(low),
        vec([high[0], low[1], low[2]]),
        vec([high[0], high[1], low[2]]),
        vec([low[0], high[1], low[2]]),
        vec([low[0], low[1], high[2]]),
        vec([high[0], low[1], high[2]]),
        vec(high),
        vec([low[0], high[1], high[2]]),
    ];
    let polygons = [
        ([0, 3, 2, 1], [0., 0., -1.]),
        ([4, 5, 6, 7], [0., 0., 1.]),
        ([0, 4, 7, 3], [-1., 0., 0.]),
        ([1, 2, 6, 5], [1., 0., 0.]),
        ([0, 1, 5, 4], [0., -1., 0.]),
        ([3, 7, 6, 2], [0., 1., 0.]),
    ];
    let mut body = Body {
        vertices,
        edges: Vec::new(),
        faces: Vec::new(),
    };
    let mut edges = BTreeMap::new();
    for (indices, n) in polygons {
        let mut lp = Vec::new();
        for i in 0..4 {
            let a = indices[i];
            let b = indices[(i + 1) % 4];
            let key = (a.min(b), a.max(b));
            let index = *edges.entry(key).or_insert_with(|| {
                let index = body.edges.len();
                body.edges.push(Edge {
                    curve: Curve::Line {
                        a: body.vertices[key.0],
                        b: body.vertices[key.1],
                    },
                    start: key.0,
                    end: key.1,
                });
                index
            });
            lp.push((index, a == key.0));
        }
        let reference = if n[0] != 0. {
            [0., 1., 0.]
        } else {
            [1., 0., 0.]
        };
        body.faces.push(Face {
            surface: Surface::Plane {
                p: body.vertices[indices[0]],
                n: vec(n),
                x: vec(reference),
            },
            sense: true,
            loops: vec![lp],
        });
    }
    body
}

fn base() -> Body {
    box_body([0., 0., 0.], [10., 8., 4.])
}

/// Independent divergence-theorem volume using oriented face triangles, not
/// the boolean operation's cells or membership classifier.
fn volume(body: &Body) -> f64 {
    let origin = body
        .vertices
        .first()
        .copied()
        .unwrap_or(Vec3::new(0., 0., 0.));
    let mut result = 0.;
    for face in &body.faces {
        for lp in &face.loops {
            let points: Vec<Vec3> = lp
                .iter()
                .map(|(i, forward)| {
                    let e = &body.edges[*i];
                    body.vertices[if *forward { e.start } else { e.end }].sub(origin)
                })
                .collect();
            for i in 1..points.len() - 1 {
                result += points[0].dot(points[i].cross(points[i + 1])) / 6.;
            }
        }
    }
    result
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-8 * b.abs().max(1.), "{a} != {b}");
}

/// Independent triangulated +X ray intersections. Queries avoid all authored
/// boundaries; used only as an independent point oracle in tests.
fn inside(body: &Body, p: Vec3) -> bool {
    let direction = Vec3::new(1., 0.137, 0.071);
    let mut distances: Vec<f64> = Vec::new();
    for f in &body.faces {
        for lp in &f.loops {
            let ps: Vec<Vec3> = lp
                .iter()
                .map(|(i, forward)| {
                    let e = &body.edges[*i];
                    body.vertices[if *forward { e.start } else { e.end }]
                })
                .collect();
            // Result faces are convex quads. This independent check is not used to
            // classify arbitrary input polygons or drive the production algorithm.
            for i in 1..ps.len() - 1 {
                let a = ps[0];
                let e1 = ps[i].sub(a);
                let e2 = ps[i + 1].sub(a);
                let h = direction.cross(e2);
                let det = e1.dot(h);
                if det.abs() < 1e-12 {
                    continue;
                }
                let s = p.sub(a);
                let u = s.dot(h) / det;
                if !(0.0..=1.0).contains(&u) {
                    continue;
                }
                let q = s.cross(e1);
                let v = direction.dot(q) / det;
                if v < 0. || u + v > 1. {
                    continue;
                }
                let t = e2.dot(q) / det;
                if t > 1e-10 && distances.iter().all(|d| (*d - t).abs() > 1e-9) {
                    distances.push(t);
                }
            }
        }
    }
    distances.len() % 2 == 1
}

#[test]
fn rectangular_cut_acceptance_table() {
    let cases = [
        ([2., 2., -1.], [4., 5., 5.], 296.),
        ([0., 2., -1.], [2., 5., 5.], 296.),
        ([11., 2., -1.], [12., 5., 5.], 320.),
        ([10., 2., -1.], [12., 5., 5.], 320.),
        ([-1., -1., -1.], [11., 9., 5.], 0.),
        ([2., 2., 2.], [4., 5., 5.], 308.),
    ];
    for (low, high, expected) in cases {
        let a = base();
        let b = box_body(low, high);
        let original = crate::neutral::write_body(&a);
        let result = cut_rectilinear(&a, &b).unwrap();
        near(volume(&result), expected);
        assert_eq!(crate::neutral::write_body(&a), original, "input mutated");
        assert!(result
            .faces
            .iter()
            .all(|f| matches!(f.surface, Surface::Plane { .. })));
        validate_faces(&result).unwrap();
    }
}

#[test]
fn through_hole_and_pocket_have_expected_material() {
    for (low, removed) in [(-1., true), (2., false)] {
        let result = cut_rectilinear(&base(), &box_body([2., 2., low], [4., 5., 5.])).unwrap();
        assert!(inside(&result, Vec3::new(1., 1., 1.)));
        assert_eq!(inside(&result, Vec3::new(3., 3., 1.)), !removed);
        assert!(!inside(&result, Vec3::new(3., 3., 3.)));
        assert!(!inside(&result, Vec3::new(20., 3., 1.)));
    }
}

#[test]
fn repeated_and_overlapping_cuts_preserve_material() {
    let first_tool = box_body([2., 2., -1.], [4., 5., 5.]);
    let first = cut_rectilinear(&base(), &first_tool).unwrap();
    let again = cut_rectilinear(&first, &first_tool).unwrap();
    near(volume(&again), 296.);
    let second = cut_rectilinear(&first, &box_body([3., 3., -1.], [5., 6., 5.])).unwrap();
    near(volume(&second), 320. - 24. - 24. + 8.);
    assert!(!inside(&second, Vec3::new(4.5, 4., 2.)));
    assert!(inside(&second, Vec3::new(1., 1., 2.)));
}

#[test]
fn empty_operands_are_explicit_success() {
    assert!(cut_rectilinear(&Body::default(), &base())
        .unwrap()
        .faces
        .is_empty());
    near(
        volume(&cut_rectilinear(&base(), &Body::default()).unwrap()),
        320.,
    );
}

#[test]
fn internal_void_is_closed_and_oriented() {
    let result = cut_rectilinear(&base(), &box_body([2., 2., 1.], [4., 5., 3.])).unwrap();
    near(volume(&result), 308.);
    assert!(!inside(&result, Vec3::new(3., 3., 2.)));
    assert!(inside(&result, Vec3::new(3., 3., 0.5)));
}

#[test]
fn finite_near_contact_is_distinct_from_exact_contact() {
    let gap = 1e-5;
    let no_contact =
        cut_rectilinear(&base(), &box_body([10. + gap, 2., -1.], [12., 5., 5.])).unwrap();
    let overlap = cut_rectilinear(&base(), &box_body([10. - gap, 2., -1.], [12., 5., 5.])).unwrap();
    near(volume(&no_contact), 320.);
    near(volume(&overlap), 320. - gap * 3. * 4.);
    let unsupported = cut_rectilinear(&base(), &box_body([10. - 1e-9, 2., -1.], [12., 5., 5.]));
    assert!(matches!(unsupported, Err(CutError::Unsupported(_))));
}

#[test]
fn coordinate_translation_and_scale_preserve_geometry() {
    for (scale, offset) in [(1., -1000.), (25.4, 0.), (0.1, 40.)] {
        let transform = |p: [f64; 3]| p.map(|v| v * scale + offset);
        let body = box_body(transform([0., 0., 0.]), transform([10., 8., 4.]));
        let tool = box_body(transform([2., 2., -1.]), transform([4., 5., 5.]));
        near(
            volume(&cut_rectilinear(&body, &tool).unwrap()),
            296. * scale.powi(3),
        );
    }
}

#[test]
fn rejects_bad_topology_and_off_plane_geometry() {
    let mut body = base();
    body.faces[0].loops[0][0].1 ^= true;
    assert!(matches!(
        cut_rectilinear(&body, &base()),
        Err(CutError::InvalidInput(_))
    ));
    let mut body = base();
    body.vertices[0].x = f64::NAN;
    assert!(matches!(
        cut_rectilinear(&body, &base()),
        Err(CutError::InvalidInput(_))
    ));
    let mut body = base();
    body.edges[0].end = 999;
    assert!(matches!(
        cut_rectilinear(&body, &base()),
        Err(CutError::InvalidInput(_))
    ));
    let mut body = base();
    if let Surface::Plane { p, .. } = &mut body.faces[0].surface {
        p.z += 1.;
    }
    assert!(matches!(
        cut_rectilinear(&body, &base()),
        Err(CutError::InvalidInput(_))
    ));
}

#[test]
fn rejects_curves_and_noncardinal_faces() {
    let mut body = base();
    body.faces[0].surface = Surface::Cylinder {
        p: Vec3::new(0., 0., 0.),
        a: Vec3::new(0., 0., 1.),
        x: Vec3::new(1., 0., 0.),
        r: 1.,
    };
    assert!(matches!(
        cut_rectilinear(&body, &base()),
        Err(CutError::Unsupported(_))
    ));
    let mut body = base();
    if let Surface::Plane { n, .. } = &mut body.faces[0].surface {
        n.x = 0.1;
    }
    assert!(matches!(
        cut_rectilinear(&body, &base()),
        Err(CutError::Unsupported(_))
    ));
}

#[test]
fn seeded_rectangular_cuts_match_independent_overlap_volume() {
    let mut state = 0x5a17_u64;
    for _ in 0..60 {
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((state >> 32) % 40) as f64 * 0.25
        };
        let low = [next() - 2., next() - 2., next() - 2.];
        let high = [
            low[0] + next() + 0.25,
            low[1] + next() + 0.25,
            low[2] + next() + 0.25,
        ];
        let overlap = (0..3)
            .map(|i| (high[i].min([10., 8., 4.][i]) - low[i].max(0.)).max(0.))
            .product::<f64>();
        let extent = [10., 8., 4.];
        let severs = (0..3).any(|axis| {
            low[axis] > 0.
                && high[axis] < extent[axis]
                && (0..3)
                    .filter(|i| *i != axis)
                    .all(|i| low[i] <= 0. && high[i] >= extent[i])
        });
        let result = cut_rectilinear(&base(), &box_body(low, high));
        if severs {
            assert!(matches!(result, Err(CutError::Unsupported(_))));
        } else {
            near(volume(&result.unwrap()), 320. - overlap);
        }
    }
}

#[test]
fn severing_cut_reports_unsupported_without_dropping_either_piece() {
    let result = cut_rectilinear(&base(), &box_body([4., -1., -1.], [6., 9., 5.]));
    assert!(matches!(result,Err(CutError::Unsupported(ref m)) if m.contains("disconnected")));
}

#[test]
fn tolerance_is_explicit_and_no_silent_plane_merging_occurs() {
    let tool = box_body([10. - 1e-5, 2., -1.], [12., 5., 5.]);
    near(volume(&cut(&base(), &tool, 1e-7).unwrap()), 320. - 12e-5);
    assert!(matches!(
        cut(&base(), &tool, 1e-4),
        Err(CutError::Unsupported(_))
    ));
    for t in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(matches!(
            cut(&base(), &tool, t),
            Err(CutError::InvalidInput(_))
        ));
    }
}

#[test]
fn wire_protocol_returns_real_cut_and_successful_empty_result() {
    for (tool, expected) in [
        (box_body([2., 2., -1.], [4., 5., 5.]), 296.),
        (box_body([-1., -1., -1.], [11., 9., 5.]), 0.),
    ] {
        let mut data = crate::neutral::write_body(&base());
        data.extend(crate::neutral::write_body(&tool));
        let bytes = crate::protocol::encode_request(crate::protocol::OP_CUT, 1e-7, &data);
        let request = crate::protocol::parse_request(&bytes).unwrap();
        let result = crate::protocol::handle(&request).unwrap();
        let mut reader = crate::neutral::Reader::new(&result);
        let body = reader.body().unwrap();
        reader.finish().unwrap();
        near(volume(&body), expected);
    }
}

/// Build a tube with one planar top/bottom face containing inner loops, as a
/// native CAD solid may represent a through hole (not the cut's tiled output).
fn rectangular_tube() -> Body {
    let mut vertices = Vec::new();
    for z in [0., 4.] {
        for (x, y) in [
            (0., 0.),
            (10., 0.),
            (10., 8.),
            (0., 8.),
            (2., 2.),
            (4., 2.),
            (4., 5.),
            (2., 5.),
        ] {
            vertices.push(Vec3::new(x, y, z));
        }
    }
    let descriptions = vec![
        ([0., 0., -1.], vec![vec![0, 3, 2, 1], vec![4, 5, 6, 7]]),
        ([0., 0., 1.], vec![vec![8, 9, 10, 11], vec![12, 15, 14, 13]]),
        ([0., -1., 0.], vec![vec![0, 1, 9, 8]]),
        ([1., 0., 0.], vec![vec![1, 2, 10, 9]]),
        ([0., 1., 0.], vec![vec![2, 3, 11, 10]]),
        ([-1., 0., 0.], vec![vec![3, 0, 8, 11]]),
        ([0., 1., 0.], vec![vec![4, 12, 13, 5]]),
        ([-1., 0., 0.], vec![vec![5, 13, 14, 6]]),
        ([0., -1., 0.], vec![vec![6, 14, 15, 7]]),
        ([1., 0., 0.], vec![vec![7, 15, 12, 4]]),
    ];
    let mut body = Body {
        vertices,
        edges: Vec::new(),
        faces: Vec::new(),
    };
    let mut edges = BTreeMap::new();
    for (normal, polygons) in descriptions {
        let point = body.vertices[polygons[0][0]];
        let mut loops = Vec::new();
        for polygon in polygons {
            let mut lp = Vec::new();
            for i in 0..polygon.len() {
                let a = polygon[i];
                let b = polygon[(i + 1) % polygon.len()];
                let key = (a.min(b), a.max(b));
                let edge = *edges.entry(key).or_insert_with(|| {
                    let index = body.edges.len();
                    body.edges.push(Edge {
                        curve: Curve::Line {
                            a: body.vertices[key.0],
                            b: body.vertices[key.1],
                        },
                        start: key.0,
                        end: key.1,
                    });
                    index
                });
                lp.push((edge, a == key.0));
            }
            loops.push(lp);
        }
        body.faces.push(Face {
            surface: Surface::Plane {
                p: point,
                n: vec(normal),
                x: if normal[0] != 0. {
                    Vec3::new(0., 1., 0.)
                } else {
                    Vec3::new(1., 0., 0.)
                },
            },
            sense: true,
            loops,
        });
    }
    body
}

#[test]
fn imported_face_holes_survive_followup_cut() {
    let tube = rectangular_tube();
    near(volume(&tube), 296.);
    let result = cut_rectilinear(&tube, &box_body([6., 2., -1.], [8., 5., 5.])).unwrap();
    near(volume(&result), 272.);
    assert!(!inside(&result, Vec3::new(3., 3., 2.)));
    assert!(!inside(&result, Vec3::new(7., 3., 2.)));
    assert!(inside(&result, Vec3::new(5., 3., 2.)));
}

use super::*;
use crate::geom::{Surface, Vec3};
use crate::neutral::Profile;
fn box_body() -> Body {
    let p = [
        Vec3::new(0., 0., 0.),
        Vec3::new(4., 0., 0.),
        Vec3::new(4., 3., 0.),
        Vec3::new(0., 3., 0.),
    ];
    crate::extrude::extrude(
        &Profile {
            loops: vec![(0..4)
                .map(|i| Curve::Line {
                    a: p[i],
                    b: p[(i + 1) % 4],
                })
                .collect()],
            n: Vec3::new(0., 0., 1.),
            low: 0.,
            high: 1.,
        },
        1e-8,
    )
    .unwrap()
}
#[test]
fn unsplit_box_preserves_faces_edges_and_body() {
    let body = box_body();
    let before = crate::neutral::write_body(&body);
    let regions = build(&body, 1e-8).unwrap();
    assert_eq!(regions.faces.len(), 6);
    assert_eq!(regions.edges.len(), 12);
    assert!(regions.faces.iter().all(|g| g.faces.len() == 1));
    assert!(regions.edges.iter().all(|g| g.edges.len() == 1));
    assert_eq!(before, crate::neutral::write_body(&body));
}
#[test]
fn collinear_extrusion_vertices_are_not_measurement_edges_or_corners() {
    let points = [
        Vec3::new(0., 0., 0.),
        Vec3::new(1., 0., 0.),
        Vec3::new(4., 0., 0.),
        Vec3::new(4., 3., 0.),
        Vec3::new(0., 3., 0.),
    ];
    let body = crate::extrude::extrude(
        &Profile {
            loops: vec![(0..5)
                .map(|i| Curve::Line {
                    a: points[i],
                    b: points[(i + 1) % 5],
                })
                .collect()],
            n: Vec3::new(0., 0., 1.),
            low: 0.,
            high: 1.,
        },
        1e-8,
    )
    .unwrap();
    let original = crate::neutral::write_body(&body);
    let regions = build(&body, 1e-8).unwrap();
    assert_eq!(regions.faces.len(), 6);
    assert_eq!(regions.edges.len(), 12);
    assert_eq!(
        regions.edges.iter().filter(|g| g.edges.len() == 2).count(),
        2
    );
    for g in regions.edges.iter().filter(|g| g.edges.len() == 2) {
        assert!((g.curve.end().distance(g.curve.start()) - 4.).abs() < 1e-12);
    }
    assert_eq!(original, crate::neutral::write_body(&body));
}
#[test]
fn circular_fragments_form_one_circle_without_erasing_holes() {
    let z = Vec3::new(0., 0., 1.);
    let x = Vec3::new(1., 0., 0.);
    let center = Vec3::new(0., 0., 0.);
    let arcs: Vec<_> = (0..8)
        .map(|i| Curve::Arc {
            c: center,
            n: z,
            x,
            r: 2.,
            t0: i as f64 * TAU / 8.,
            t1: (i + 1) as f64 * TAU / 8.,
        })
        .collect();
    let body = crate::extrude::extrude(
        &Profile {
            loops: vec![arcs],
            n: z,
            low: 0.,
            high: 1.,
        },
        1e-8,
    )
    .unwrap();
    let regions = build(&body, 1e-8).unwrap();
    assert_eq!(regions.faces.len(), 3);
    assert_eq!(regions.edges.len(), 2);
    for g in &regions.edges {
        assert!(g.closed);
        assert_eq!(g.edges.len(), 8);
        let Curve::Arc { t0, t1, r, .. } = g.curve else {
            panic!()
        };
        assert!((t1 - t0 - TAU).abs() < 1e-12);
        assert_eq!(r, 2.);
    }
}
#[test]
fn analytic_support_nearness_is_not_permission_to_merge() {
    let p = Vec3::new(0., 0., 0.);
    let x = Vec3::new(1., 0., 0.);
    let z = Vec3::new(0., 0., 1.);
    let line = Curve::Line { a: p, b: x };
    assert!(!same_curve(
        &line,
        &Curve::Line {
            a: p,
            b: Vec3::new(1., 1e-8, 0.)
        }
    ));
    let circle = Curve::Arc {
        c: p,
        n: z,
        x,
        r: 1.,
        t0: 0.,
        t1: 1.,
    };
    assert!(!same_curve(
        &circle,
        &Curve::Arc {
            c: Vec3::new(1e-8, 0., 0.),
            n: z,
            x,
            r: 1.,
            t0: 1.,
            t1: 2.
        }
    ));
    assert!(!same_curve(
        &circle,
        &Curve::Arc {
            c: p,
            n: z,
            x,
            r: 1. + 1e-8,
            t0: 1.,
            t1: 2.
        }
    ));
    let far = Vec3::new(1e12, 0., 0.);
    assert!(!same_curve(
        &circle.translated(far),
        &circle.translated(far.add(x))
    ));
}
#[test]
fn malformed_input_is_refused_before_indexing() {
    let mut body = box_body();
    body.edges[0].start = usize::MAX;
    assert!(build(&body, 1e-8).is_err());
    assert!(build(&box_body(), f64::NAN).is_err());
}

#[test]
fn disconnected_coplanar_faces_stay_separate() {
    let mut body = box_body();
    let mut other = box_body();
    let shift = Vec3::new(10., 0., 0.);
    let (nv, ne) = (body.vertices.len(), body.edges.len());
    for p in &mut other.vertices {
        *p = p.add(shift);
    }
    for edge in &mut other.edges {
        edge.start += nv;
        edge.end += nv;
        edge.curve = edge.curve.translated(shift);
    }
    for face in &mut other.faces {
        let Surface::Plane { p, .. } = &mut face.surface else {
            panic!()
        };
        *p = p.add(shift);
        for (edge, _) in face.loops.iter_mut().flatten() {
            *edge += ne;
        }
    }
    body.vertices.extend(other.vertices);
    body.edges.extend(other.edges);
    body.faces.extend(other.faces);
    let regions = build(&body, 1e-8).unwrap();
    assert_eq!(regions.faces.len(), 12);
    assert_eq!(regions.edges.len(), 24);
    assert!(regions.faces.iter().all(|face| face.faces.len() == 1));
}

#[test]
fn real_angle_top_and_slot_are_whole_regions() {
    let bytes = include_bytes!("../../tests/fixtures/selection/angle.f64");
    let values: Vec<f64> = bytes
        .chunks_exact(8)
        .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
        .collect();
    let body = crate::neutral::Reader::new(&values).body().unwrap();
    let before = crate::neutral::write_body(&body);
    let regions = build(&body, 1e-7).unwrap();
    println!(
        "Angle: {} raw faces -> {} regions; {} raw edges -> {} boundaries",
        body.faces.len(),
        regions.faces.len(),
        body.edges.len(),
        regions.edges.len()
    );
    let top: Vec<_> = body
        .faces
        .iter()
        .enumerate()
        .filter_map(|(i, f)| match f.surface {
            Surface::Plane { p, n, .. } if p.z == 0. && n.z == 1. => Some(i),
            _ => None,
        })
        .collect();
    assert!(top.len() > 10);
    let top_id = regions.face_region[top[0]];
    assert!(top.iter().all(|&f| regions.face_region[f] == top_id));
    assert_eq!(regions.faces[top_id].faces.len(), top.len());
    let areas = crate::mass::face_areas(&body, 1e-7).unwrap();
    let top_area: f64 = top.iter().map(|&f| areas[f].area).sum();
    assert!(
        (top_area - (108. - std::f64::consts::PI * 0.5 - 1.)).abs() < 1e-9,
        "{top_area}"
    );
    for (center, count) in [
        (Vec3::new(1.1464466094067263, 3.646446609406726, 0.), 4),
        (Vec3::new(1.8535533905932737, 4.353553390593274, 0.), 4),
    ] {
        let raw: Vec<_> = body
            .faces
            .iter()
            .enumerate()
            .filter_map(|(i, f)| match f.surface {
                Surface::Cylinder { p, r, .. } if p.distance(center) < 1e-12 && r == 0.5 => Some(i),
                _ => None,
            })
            .collect();
        assert_eq!(raw.len(), count);
        assert!(raw
            .iter()
            .all(|&f| regions.face_region[f] == regions.face_region[raw[0]]));
        let arcs:Vec<_>=regions.edges.iter().filter(|e|matches!(e.curve,Curve::Arc{c,r,..} if (c.x-center.x).abs()<1e-12&&(c.y-center.y).abs()<1e-12&&r==0.5)).collect();
        assert_eq!(arcs.len(), 2);
        for edge in arcs {
            let Curve::Arc { t0, t1, .. } = edge.curve else {
                panic!()
            };
            assert!((t1 - t0 - std::f64::consts::PI).abs() < 1e-12);
            assert_eq!(edge.edges.len(), 4);
        }
    }
    assert!(
        before == crate::neutral::write_body(&body),
        "Selection changed the raw body"
    );
}

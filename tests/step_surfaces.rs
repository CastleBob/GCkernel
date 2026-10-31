use gc_kernel::{
    frame::{Point3, Rigid},
    geom::{Surface, Vec3},
    mass, neutral, step, tessellate, topology,
};
const POLY: &str = include_str!("fixtures/SplineSurfaceBox.step");
const RATIONAL: &str = include_str!("fixtures/RationalSplineSurfaceBox.step");
const TOL: f64 = 1e-6;
#[test]
fn curved_tensor_faces_preserve_closed_solid_mass_mesh_and_frames() {
    for source in [
        POLY,
        RATIONAL,
        include_str!("fixtures/CurvedBoundarySplineBox.step"),
        include_str!("fixtures/MultispanSplineSurfaceBox.step"),
    ] {
        let item = step::read(source.as_bytes(), TOL).unwrap().bodies.remove(0);
        let body = item.body;
        let expected = mass::measure(&body, TOL).unwrap().volume;
        if source == POLY {
            assert!((expected - 25.5).abs() < 1e-9, "{expected}");
        } else {
            assert!(expected > 25.5 && expected < 27.);
        }
        let motion = Rigid::rotation(Vec3::new(0., 1., 0.), 0.37, Point3::new(0., 0., 0.))
            .unwrap()
            .then(&Rigid::translation(Vec3::new(200., -100., 30.)).unwrap());
        for body in [body.clone(), motion.body(&body)] {
            assert_eq!(topology::check(&body, TOL).unwrap().genus(), 0.);
            let m = mass::measure(&body, TOL).unwrap();
            assert!(!m.exact);
            assert!((m.volume - expected).abs() < 1e-8);
            let wire = neutral::write_body(&body);
            let mut r = neutral::Reader::with_tolerance(&wire, TOL);
            let restored = r.body().unwrap();
            r.finish().unwrap();
            topology::check(&restored, TOL).unwrap();
            let mesh = tessellate::tessellate(&body, 0.01).unwrap();
            let o = mesh.positions[0];
            let volume: f64 = mesh
                .triangles
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(|i| mesh.positions[i].sub(o));
                    a.dot(b.cross(c)) / 6.
                })
                .sum();
            assert!(
                (volume - expected).abs() < 0.1,
                "mesh {volume}, expected {expected}"
            );
            // Every canonical edge-polyline segment occurs twice geometrically.
            let positions = &mesh.positions;
            for line in &mesh.edges {
                for pair in line.windows(2) {
                    let count = mesh
                        .triangles
                        .iter()
                        .flat_map(|t| {
                            (0..3).map(move |i| (positions[t[i]], positions[t[(i + 1) % 3]]))
                        })
                        .filter(|(a, b)| {
                            (*a == pair[0] && *b == pair[1]) || (*a == pair[1] && *b == pair[0])
                        })
                        .count();
                    assert_eq!(count, 2, "unshared boundary {pair:?}");
                }
            }
        }
    }
}
#[test]
fn surface_refusals_are_atomic_and_keep_supported_source_usable() {
    for (from, to) in [
        ("(4,4),(4,4)", "(3,5),(4,4)"),
        ("3,3,((#800", "4,3,((#800"),
        (".UNSPECIFIED.,.F.,.F.,.F.", ".UNSPECIFIED.,.T.,.F.,.F."),
        (
            "#804=CARTESIAN_POINT('',(0.66666666666666663,0.0,4.0))",
            "#804=CARTESIAN_POINT('',(0.66666666666666663,0.0,4.01))",
        ),
    ] {
        assert!(POLY.contains(from));
        assert!(
            step::read(POLY.replace(from, to).as_bytes(), TOL).is_err(),
            "{to}"
        );
    }
    assert!(step::read(RATIONAL.replace("1.25", "-1.25").as_bytes(), TOL).is_err());
    let mut body = step::read(POLY.as_bytes(), TOL)
        .unwrap()
        .bodies
        .remove(0)
        .body;
    let f = body
        .faces
        .iter_mut()
        .find(|f| matches!(f.surface, Surface::Spline { .. }))
        .unwrap();
    f.sense = !f.sense;
    assert!(topology::check(&body, TOL).is_err());
    step::read(POLY.as_bytes(), TOL).unwrap();
}
#[test]
fn vertex_and_cartesian_selectors_trim_without_changing_source_geometry() {
    let source = include_str!("fixtures/vertex_trimmed_spline_box.step");
    for source in [
        source.to_owned(),
        source.replace("#600,.T.)", "#102,.T.)").replace(
            "(PARAMETER_VALUE(-4.),#1),(PARAMETER_VALUE(16.),#2),.T.,.PARAMETER.",
            "(#1),(#2),.T.,.CARTESIAN.",
        ),
    ] {
        let b = step::read(source.as_bytes(), TOL)
            .unwrap()
            .bodies
            .remove(0)
            .body;
        assert!((mass::measure(&b, TOL).unwrap().volume - 24.).abs() < 1e-10);
    }
}

#[test]
fn vertex_trimmed_rational_seam_preserves_the_analytic_solid() {
    let mut source = include_str!("fixtures/trimmed_rational_cylinder.step").to_owned();
    for i in 0..4 {
        for (trim, basis) in [(900 + i, 11), (910 + i, 12)] {
            let sense = if i % 2 == 0 { "T" } else { "F" };
            source = source.replace(&format!(",#{trim},.{sense}.)"), &format!(",#{basis},.T.)"));
        }
    }
    let b = step::read(source.as_bytes(), TOL)
        .unwrap()
        .bodies
        .remove(0)
        .body;
    assert!((mass::measure(&b, TOL).unwrap().volume - 16. * std::f64::consts::PI).abs() < 1e-9);
    topology::check(&b, TOL).unwrap();
    tessellate::tessellate(&b, 0.005).unwrap();
}

#[test]
fn reversed_surface_chart_with_false_face_sense_preserves_solid() {
    let mut body = step::read(POLY.as_bytes(), TOL)
        .unwrap()
        .bodies
        .remove(0)
        .body;
    for f in &mut body.faces {
        if let Surface::Spline { data } = &mut f.surface {
            for p in &mut data.patches {
                let before = p.clone();
                for i in 0..4 {
                    for j in 0..4 {
                        p.p[i * 4 + j] = before.p[(3 - i) * 4 + j];
                        p.w[i * 4 + j] = before.w[(3 - i) * 4 + j];
                    }
                }
            }
            f.sense = false;
        }
    }
    topology::check(&body, TOL).unwrap();
    assert!((mass::measure(&body, TOL).unwrap().volume - 25.5).abs() < 1e-9);
    tessellate::tessellate(&body, 0.01).unwrap();
}

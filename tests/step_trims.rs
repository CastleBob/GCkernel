use gc_kernel::{
    frame::{Point3, Rigid},
    geom::Vec3,
    mass, neutral, step, tessellate, topology,
};
use std::f64::consts::PI;

const SOURCE: &str = include_str!("fixtures/trimmed_rational_cylinder.step");
const TOL: f64 = 1e-6;

#[test]
fn interior_and_seam_trims_form_an_exact_shared_closed_solid() {
    let item = step::read(SOURCE.as_bytes(), TOL).unwrap().bodies.remove(0);
    let body = item.body;
    assert_eq!(
        (body.vertices.len(), body.edges.len(), body.faces.len()),
        (16, 17, 3)
    );
    let motion = Rigid::rotation(Vec3::new(0., 1., 0.), 0.37, Point3::new(0., 0., 0.))
        .unwrap()
        .then(&Rigid::translation(Vec3::new(200., -100., 30.)).unwrap());
    for body in [body.clone(), motion.body(&body)] {
        assert_eq!(topology::check(&body, TOL).unwrap().genus(), 0.);
        let m = mass::measure(&body, TOL).unwrap();
        assert!((m.volume - 16. * PI).abs() < 1e-9, "{m:?}");
        assert!((m.area - 24. * PI).abs() < 1e-9, "{m:?}");
        let wire = neutral::write_body(&body);
        let mut reader = neutral::Reader::with_tolerance(&wire, TOL);
        topology::check(&reader.body().unwrap(), TOL).unwrap();
        reader.finish().unwrap();
        let mesh = tessellate::tessellate(&body, 0.005).unwrap();
        let origin = mesh.positions[0];
        let volume: f64 = mesh
            .triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| mesh.positions[i].sub(origin));
                a.dot(b.cross(c)) / 6.
            })
            .sum();
        assert!((volume - 16. * PI).abs() < 0.15, "display volume {volume}");
    }
}

#[test]
fn parameter_preference_checks_redundant_cartesian_selector_without_repair() {
    let consistent = SOURCE.replace(
        "#900=TRIMMED_CURVE('',#11,(PARAMETER_VALUE(-5.)),",
        "#900=TRIMMED_CURVE('',#11,(#1,PARAMETER_VALUE(-5.)),",
    );
    step::read(consistent.as_bytes(), TOL).unwrap();
    // Point #5 is the origin, off the curve. Parameter is authoritative but the
    // contradictory Cartesian selector must still fail the strict contract.
    let bad = consistent.replace("(#1,PARAMETER_VALUE(-5.))", "(#5,PARAMETER_VALUE(-5.))");
    assert!(step::read(bad.as_bytes(), TOL)
        .unwrap_err()
        .to_string()
        .contains("parameter/Cartesian trim mismatch"));
    assert!(step::read(
        consistent
            .replace(".PARAMETER.", ".UNSPECIFIED.")
            .as_bytes(),
        TOL
    )
    .is_err());
    step::read(
        SOURCE.replace(".PARAMETER.", ".UNSPECIFIED.").as_bytes(),
        TOL,
    )
    .unwrap();
}

#[test]
fn unsupported_or_contradictory_trims_fail_atomically_and_leave_original_usable() {
    for (from, to, reason) in [
        (
            "PARAMETER_VALUE(-5.)",
            "PARAMETER_VALUE(-8.)",
            "outside the source knot domain",
        ),
        (
            "#903=TRIMMED_CURVE('',#11,(PARAMETER_VALUE(-5.)),(PARAMETER_VALUE(9.5))",
            "#903=TRIMMED_CURVE('',#11,(PARAMETER_VALUE(-5.)),(PARAMETER_VALUE(-5.))",
            "equal STEP trim parameters",
        ),
        (".PARAMETER.", ".CARTESIAN.", "Cartesian-preferred"),
        (
            "(PARAMETER_VALUE(-5.))",
            "(PARAMETER_VALUE(-5.),PARAMETER_VALUE(-5.))",
            "duplicate",
        ),
        ("(PARAMETER_VALUE(-5.))", "()", "trim set"),
        (
            "PARAMETER_VALUE(-5.)",
            "LENGTH_MEASURE(-5.)",
            "invalid STEP trim selector",
        ),
        (
            ".CIRCULAR_ARC.,.T.,.F.",
            ".CIRCULAR_ARC.,.F.,.F.",
            "trim order",
        ),
        (
            "#900=TRIMMED_CURVE('',#11,",
            "#900=TRIMMED_CURVE('',#900,",
            "not a supported spline",
        ),
        (
            "#900=TRIMMED_CURVE('',#11,",
            "#900=TRIMMED_CURVE('',#901,",
            "not a supported spline",
        ),
        (
            "#15=EDGE_CURVE('',#3,#1100,#900,.T.)",
            "#15=EDGE_CURVE('',#3,#1100,#900,.F.)",
            "EDGE_CURVE vertices",
        ),
        (
            "#2108=CARTESIAN_POINT('',(2.,0.,0.))",
            "#2108=CARTESIAN_POINT('',(2.01,0.,0.))",
            "nonclosing seam",
        ),
    ] {
        assert!(SOURCE.contains(from), "{from}");
        let changed = SOURCE.replace(from, to);
        let error = step::read(changed.as_bytes(), TOL).unwrap_err().to_string();
        assert!(error.contains(reason), "expected {reason}, got {error}");
    }
    assert!(step::read(SOURCE.as_bytes(), TOL).is_ok());
}

#[test]
fn trim_sense_edge_sense_and_face_orientation_compose_independently() {
    for trim_forward in [false, true] {
        for edge_forward in [false, true] {
            let mut source = SOURCE.to_owned();
            if !trim_forward {
                source = source.replace(
                    "#900=TRIMMED_CURVE('',#11,(PARAMETER_VALUE(-5.)),(PARAMETER_VALUE(-0.5)),.T.",
                    "#900=TRIMMED_CURVE('',#11,(PARAMETER_VALUE(-0.5)),(PARAMETER_VALUE(-5.)),.F.",
                );
            }
            let backwards = trim_forward != edge_forward;
            let (v0, v1) = if backwards { (1100, 3) } else { (3, 1100) };
            source = source.replace(
                "#15=EDGE_CURVE('',#3,#1100,#900,.T.);",
                &format!(
                    "#15=EDGE_CURVE('',#{v0},#{v1},#900,.{}.);",
                    if edge_forward { "T" } else { "F" }
                ),
            );
            if backwards {
                source = source
                    .replace(
                        "ORIENTED_EDGE('',*,*,#15,.T.)",
                        "ORIENTED_EDGE('',*,*,#15,.SWAP.)",
                    )
                    .replace(
                        "ORIENTED_EDGE('',*,*,#15,.F.)",
                        "ORIENTED_EDGE('',*,*,#15,.T.)",
                    )
                    .replace(
                        "ORIENTED_EDGE('',*,*,#15,.SWAP.)",
                        "ORIENTED_EDGE('',*,*,#15,.F.)",
                    );
            }
            let body = step::read(source.as_bytes(), TOL).unwrap().bodies.remove(0);
            assert!((body.volume - 16. * PI).abs() < 1e-9);
        }
    }
}

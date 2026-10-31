use gc_kernel::{
    boolean,
    frame::{Point3, Rigid},
    geom::{Curve, Vec3},
    mass, neutral, rational, step, tessellate, topology,
};
use std::f64::consts::PI;

const TOL: f64 = 1e-6;
const SOURCE: &str = include_str!("fixtures/rational_cylinder.step");

#[test]
fn complete_rational_solid_roundtrips_moves_measures_and_tessellates() {
    let imported = step::read(SOURCE.as_bytes(), TOL).unwrap();
    let body = &imported.bodies[0].body;
    assert_eq!(
        (body.vertices.len(), body.edges.len(), body.faces.len()),
        (8, 9, 3)
    );
    assert_eq!(
        body.edges
            .iter()
            .filter(|e| matches!(e.curve, Curve::RationalBezier { .. }))
            .count(),
        8
    );
    let motion = Rigid::rotation(
        Vec3::new(1., 2., 2.).scale(1. / 3.),
        0.73,
        Point3::new(0., 0., 0.),
    )
    .unwrap()
    .then(&Rigid::translation(Vec3::new(1000., -400., 600.)).unwrap());
    for b in [body.clone(), motion.body(body)] {
        let wire = neutral::write_body(&b);
        let mut reader = neutral::Reader::with_tolerance(&wire, TOL);
        let b = reader.body().unwrap();
        reader.finish().unwrap();
        // Existing neutral frames normalize unit directions; rational data itself is bit-preserved.
        assert!(neutral::write_body(&b)
            .iter()
            .zip(&wire)
            .all(|(a, b)| (a - b).abs() < 1e-12));
        assert_eq!(topology::check(&b, TOL).unwrap().genus(), 0.0);
        let m = mass::measure(&b, TOL).unwrap();
        assert!(!m.exact);
        assert!((m.volume - 16. * PI).abs() < 1e-9, "{m:?}");
        assert!((m.area - 24. * PI).abs() < 1e-9, "{m:?}");
        let mesh = tessellate::tessellate(&b, 0.005).unwrap();
        assert!(!mesh.triangles.is_empty());
        assert_eq!(mesh.edges.len(), 9);
        // Independent triangle divergence sum checks orientation and complete display solid.
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
        for (e, line) in b.edges.iter().zip(&mesh.edges) {
            for i in 0..=1000 {
                let q = e.curve.eval(i as f64 / 1000.).0;
                let gap = line
                    .windows(2)
                    .map(|ab| {
                        let d = ab[1].sub(ab[0]);
                        let t = (q.sub(ab[0]).dot(d) / d.dot(d)).clamp(0., 1.);
                        q.distance(ab[0].add(d.scale(t)))
                    })
                    .fold(f64::INFINITY, f64::min);
                assert!(gap <= 0.005 + 1e-10, "display gap {gap}");
            }
        }
        assert!(boolean::cut(&b, &b, TOL)
            .unwrap_err()
            .to_string()
            .contains("rational Bezier"));
    }
}

#[test]
fn rational_derivative_subdivision_reversal_and_weight_scaling_preserve_geometry() {
    let body = step::read(SOURCE.as_bytes(), TOL)
        .unwrap()
        .bodies
        .remove(0)
        .body;
    let (p, w) = body
        .edges
        .iter()
        .find_map(|e| match e.curve {
            Curve::RationalBezier { p, w } => Some((p, w)),
            _ => None,
        })
        .unwrap();
    let c = Curve::RationalBezier { p, w };
    let (cp, cw) = rational::piece(p, w, 0.2, 0.8);
    let piece = Curve::RationalBezier { p: cp, w: cw };
    for i in 1..100 {
        let t = i as f64 / 100.;
        let (q, d) = c.eval(t);
        assert!((q.x.hypot(q.y) - 2.).abs() < 2e-15);
        let numerical = c.eval(t + 1e-5).0.sub(c.eval(t - 1e-5).0).scale(0.5e5);
        assert!(d.distance(numerical) < 1e-8);
        assert!(q.distance(c.reversed().eval(1. - t).0) < 2e-15);
        assert!(q.distance(rational::eval(p, w.map(|w| w * 1e200), t).0) < 2e-15);
        assert!(piece.eval(t).0.distance(c.eval(0.2 + 0.6 * t).0) < 2e-15);
    }
    let mut wire = Vec::new();
    neutral::put_curve(&mut wire, &c);
    for value in [0., -1., f64::NAN, f64::INFINITY, 1e-20] {
        let mut bad = wire.clone();
        bad[13] = value;
        assert!(neutral::Reader::with_tolerance(&bad, TOL).curve().is_err());
    }
}

#[test]
fn source_weight_errors_and_strict_surface_mismatch_are_atomic_failures() {
    for (old, new) in [
        ("0.7071067811865475244", "0."),
        ("0.7071067811865475244", "-1."),
        ("0.7071067811865475244", "1E-20"),
        ("0.7071067811865475244", "0.708"),
        ("RATIONAL_B_SPLINE_CURVE((1.,", "RATIONAL_B_SPLINE_CURVE(("),
    ] {
        assert!(
            step::read(SOURCE.replace(old, new).as_bytes(), TOL).is_err(),
            "{new}"
        );
    }
    assert!(step::read(SOURCE.as_bytes(), TOL).is_ok());
}

use gc_kernel::{mass, neutral, step, tessellate, topology};
const TOL: f64 = 1e-6;
#[test]
fn independent_analytic_and_spline_fixtures_survive_export() {
    for source in [
        include_str!("fixtures/vertex_trimmed_spline_box.step"),
        include_str!("fixtures/RationalSplineSurfaceBox.step"),
        include_str!("fixtures/MultispanSplineSurfaceBox.step"),
        include_str!("fixtures/CurvedBoundarySplineBox.step"),
        include_str!("fixtures/trimmed_rational_cylinder.step"),
    ] {
        let body = step::read(source.as_bytes(), TOL)
            .unwrap()
            .bodies
            .remove(0)
            .body;
        let name = "O'Brien \\ pipe Ω 🧰";
        let encoded = step::write(&[(name, &body)], TOL).unwrap();
        assert_eq!(encoded, step::write(&[(name, &body)], TOL).unwrap());
        let restored = step::read(encoded.as_bytes(), TOL)
            .unwrap()
            .bodies
            .remove(0);
        assert_eq!(restored.name, name);
        assert_eq!(restored.source_length_to_inches, 1.);
        let (a, b) = (
            mass::measure(&body, TOL).unwrap(),
            mass::measure(&restored.body, TOL).unwrap(),
        );
        assert!((a.volume - b.volume).abs() < 1e-9);
        assert!((a.area - b.area).abs() < 1e-9);
        assert_eq!(
            topology::check(&body, TOL).unwrap(),
            topology::check(&restored.body, TOL).unwrap()
        );
        tessellate::tessellate(&restored.body, 0.02).unwrap();
    }
}
#[test]
fn batch_export_is_atomic_and_leaves_input_unchanged() {
    let b = step::read(
        include_bytes!("fixtures/vertex_trimmed_spline_box.step"),
        TOL,
    )
    .unwrap()
    .bodies
    .remove(0)
    .body;
    let before = neutral::write_body(&b);
    let text = step::write(&[("first", &b), ("second", &b)], TOL).unwrap();
    let read = step::read(text.as_bytes(), TOL).unwrap();
    assert_eq!(read.bodies.len(), 2);
    assert_eq!(read.bodies[0].name, "first");
    assert_eq!(read.bodies[1].name, "second");
    let mut bad = b.clone();
    bad.faces.pop();
    assert!(step::write(&[("good", &b), ("bad", &bad)], TOL).is_err());
    assert!(step::write(&[], TOL).is_err());
    assert!(step::write(&[("good", &b)], f64::NAN).is_err());
    assert_eq!(before, neutral::write_body(&b));
}

#[test]
fn export_protocol_rejects_truncation_and_preserves_batch() {
    let body = step::read(
        include_bytes!("fixtures/vertex_trimmed_spline_box.step"),
        TOL,
    )
    .unwrap()
    .bodies
    .remove(0)
    .body;
    let values = neutral::write_body(&body);
    let mut data = vec![1., 1., 1., b'A' as f64, values.len() as f64];
    data.extend(values);
    let reply = step::wire::write_op(&data, TOL).unwrap();
    assert_eq!(reply[0], 1.);
    assert_eq!(reply[1] as usize, reply.len() - 2);
    let bytes = reply[2..].iter().map(|&v| v as u8).collect::<Vec<_>>();
    assert_eq!(step::read(&bytes, TOL).unwrap().bodies[0].name, "A");
    for n in 0..data.len() {
        assert!(step::wire::write_op(&data[..n], TOL).is_err());
    }
    for bad in [f64::NAN, -1., 0.5, 256.] {
        let mut packet = data.clone();
        packet[3] = bad;
        assert!(step::wire::write_op(&packet, TOL).is_err());
    }
    data.push(0.);
    assert!(step::wire::write_op(&data, TOL).is_err());
}

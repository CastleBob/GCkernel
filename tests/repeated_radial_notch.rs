use gc_kernel::{
    boolean,
    geom::{Curve, Vec3},
    neutral::{write_body, Body, Profile, Reader},
    revolve::{revolve, Section},
    topology,
};
use std::f64::consts::TAU;
const TOL: f64 = 1e-6;
fn stock(bevel: bool) -> Body {
    let mut points = vec![(1.5, 0.)];
    if bevel {
        points.extend([(1.5, 1.45), (1.4, 1.5)]);
    } else {
        points.push((1.5, 1.5));
    }
    points.extend([(1.375, 1.5), (1.375, 0.)]);
    revolve(
        &Section {
            p: Vec3::new(0., 0., 0.),
            a: Vec3::new(0., 1., 0.),
            x: Vec3::new(1., 0., 0.),
            points,
        },
        TOL,
    )
    .unwrap()
}
fn tool(radius: f64) -> Body {
    tool_with_reach(radius, 1.501)
}
fn tool_with_reach(radius: f64, high: f64) -> Body {
    gc_kernel::extrude::extrude(
        &Profile {
            loops: vec![vec![Curve::Arc {
                c: Vec3::new(0., 0., 0.),
                n: Vec3::new(0., 0., 1.),
                x: Vec3::new(1., 0., 0.),
                r: radius,
                t0: 0.,
                t1: TAU,
            }]],
            n: Vec3::new(0., 0., 1.),
            low: 0.5,
            high,
        },
        TOL,
    )
    .unwrap()
}
fn near_bevel(height: f64, skin: f64) -> Body {
    revolve(
        &Section {
            p: Vec3::new(0., 0., 0.),
            a: Vec3::new(0., 1., 0.),
            x: Vec3::new(1., 0., 0.),
            points: vec![
                (1.5 - height, 0.),
                (1.5, skin),
                (1.5, 3.),
                (1.3, 3.),
                (1.3, 0.),
            ],
        },
        TOL,
    )
    .unwrap()
}
#[test]
fn repeated_near_bevel_notches_preserve_the_supplied_body() {
    for (height, skin, rho) in [
        (0.1, 0.1 * 30_f64.to_radians().tan(), 0.5),
        (0.1, 0.1, 0.38),
        (0.05, 0.2, 0.4),
    ] {
        let t = tool(rho);
        let mut body = boolean::cut(&near_bevel(height, skin), &t, TOL).unwrap();
        assert_eq!(
            (body.faces.len(), body.edges.len(), body.vertices.len()),
            (6, 16, 11)
        );
        let original = write_body(&body);
        for _ in 0..4 {
            body = boolean::cut(&body, &t, TOL).unwrap();
            assert_eq!(write_body(&body), original);
            topology::check(&body, TOL).unwrap();
        }
    }
}
#[test]
fn near_bevel_repeat_proof_survives_wire_roundtrip_and_cardinal_placement() {
    for (rows, translation) in [
        (
            [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            Vec3::new(0., 0., 0.),
        ),
        (
            [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
            Vec3::new(4., -8., 3.),
        ),
        (
            [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
            Vec3::new(0., 0., 0.),
        ),
    ] {
        let rigid = gc_kernel::frame::Rigid::from_rows(rows, translation).unwrap();
        let t = rigid.body(&tool_with_reach(0.5, 2.));
        let source = rigid.body(&near_bevel(0.1, 0.25));
        let body = boolean::cut(&source, &t, TOL).unwrap();
        let data = write_body(&body);
        let mut reader = Reader::with_tolerance(&data, TOL);
        let back = reader.body().unwrap();
        reader.finish().unwrap();
        let repeated = boolean::cut(&back, &t, TOL)
            .unwrap_or_else(|e| panic!("rows={rows:?} translation={translation:?}: {e:?}"));
        assert_eq!(write_body(&repeated), data);
    }
}
#[test]
fn a_wire_normalized_translated_boundary_does_not_get_approximately_matched() {
    let rigid = gc_kernel::frame::Rigid::from_rows(
        [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
        Vec3::new(4., -8., 3.),
    )
    .unwrap();
    let t = rigid.body(&tool(0.5));
    let body = boolean::cut(&rigid.body(&near_bevel(0.1, 0.25)), &t, TOL).unwrap();
    let data = write_body(&body);
    assert_eq!(write_body(&boolean::cut(&body, &t, TOL).unwrap()), data);
    let mut reader = Reader::with_tolerance(&data, TOL);
    let back = reader.body().unwrap();
    reader.finish().unwrap();
    // The existing wire reader normalizes a 0.9999999999999999 direction to
    // 1.0 here. An exact construction certificate must not wave that through.
    let normalized = write_body(&back);
    assert_ne!(normalized, data);
    assert!(boolean::cut(&back, &t, TOL).is_err());
    assert_eq!(write_body(&back), normalized);
}
#[test]
fn near_bevel_changes_to_tool_trim_or_cone_are_not_skipped() {
    let t = tool(0.5);
    let body = boolean::cut(&near_bevel(0.1, 0.1), &t, TOL).unwrap();
    for radius in [0.5 + 1e-10, 0.501, 0.499] {
        let other = tool(radius);
        let before = (write_body(&body), write_body(&other));
        assert!(boolean::cut(&body, &other, TOL).is_err());
        assert_eq!(before, (write_body(&body), write_body(&other)));
    }
    for vertex in 5..body.vertices.len() {
        let mut altered = body.clone();
        altered.vertices[vertex].x += 1e-8;
        topology::check(&altered, TOL).unwrap();
        let before = write_body(&altered);
        assert!(
            boolean::cut(&altered, &t, TOL).is_err(),
            "trim vertex {vertex}"
        );
        assert_eq!(before, write_body(&altered));
    }
    let mut altered = body.clone();
    for face in &mut altered.faces {
        if let gc_kernel::geom::Surface::Cone { k, .. } = &mut face.surface {
            *k += 1e-10;
        }
    }
    assert!(topology::check(&altered, TOL).is_err()); // ConePair no longer lies on the altered support.
    let before = write_body(&altered);
    assert!(boolean::cut(&altered, &t, TOL).is_err());
    assert_eq!(before, write_body(&altered));
}
#[test]
fn repeated_flat_notches_preserve_every_supplied_neutral_value() {
    for bevel in [false, true] {
        let t = tool(0.38);
        let mut body = boolean::cut(&stock(bevel), &t, TOL).unwrap();
        let original = write_body(&body);
        let source = write_body(&t);
        for _ in 0..4 {
            body = boolean::cut(&body, &t, TOL).unwrap();
            assert_eq!(write_body(&body), original);
            assert_eq!(write_body(&t), source);
            topology::check(&body, TOL).unwrap();
        }
    }
}
#[test]
fn the_certificate_uses_geometry_after_neutral_round_trip_and_cardinal_transform() {
    let t = tool(0.38);
    let body = boolean::cut(&stock(true), &t, TOL).unwrap();
    let data = write_body(&body);
    let mut reader = Reader::with_tolerance(&data, TOL);
    let back = reader.body().unwrap();
    reader.finish().unwrap();
    assert_eq!(write_body(&boolean::cut(&back, &t, TOL).unwrap()), data);
    let rigid = gc_kernel::frame::Rigid::from_rows(
        [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
        Vec3::new(4., 8., -3.),
    )
    .unwrap();
    // Construct after transformation so the certificate can retain exact values.
    let tt = rigid.body(&t);
    let moved = boolean::cut(&rigid.body(&stock(false)), &tt, TOL).unwrap();
    assert_eq!(
        write_body(&boolean::cut(&moved, &tt, TOL).unwrap()),
        write_body(&moved)
    );
}
#[test]
fn near_matching_tools_and_altered_notch_boundaries_are_not_skipped() {
    let t = tool(0.38);
    let body = boolean::cut(&stock(true), &t, TOL).unwrap();
    for radius in [0.381, 0.38 + 1e-10, 0.379] {
        let other = tool(radius);
        let saved = (write_body(&body), write_body(&other));
        assert!(matches!(
            boolean::cut(&body, &other, TOL),
            Err(gc_kernel::KernelError::Unsupported(_))
        ));
        assert_eq!(saved, (write_body(&body), write_body(&other)));
    }
    let mut altered = body.clone();
    let last = altered.vertices.len() - 1;
    altered.vertices[last].x += 1e-8;
    topology::check(&altered, TOL).unwrap(); // Within general topology tolerance, but not the exact certificate.
    let saved = write_body(&altered);
    assert!(matches!(
        boolean::cut(&altered, &t, TOL),
        Err(gc_kernel::KernelError::Unsupported(_))
    ));
    assert_eq!(write_body(&altered), saved);
}
#[test]
fn a_no_change_result_still_reports_its_actual_history_limit() {
    for (source, radius) in [(stock(true), 0.38), (near_bevel(0.1, 0.1), 0.5)] {
        let t = tool(radius);
        let body = boolean::cut(&source, &t, TOL).unwrap();
        let (same, history) = boolean::cut_recorded(&body, &t, TOL).unwrap();
        assert_eq!(write_body(&same), write_body(&body));
        assert!(history.unwrap_err().contains("UnchangedRadialNotch"));
        assert!(
            matches!(boolean::cut_named(&body,&t,TOL),Err(gc_kernel::KernelError::Unsupported(e)) if e.contains("UnchangedRadialNotch"))
        );
    }
}

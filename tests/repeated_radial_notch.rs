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
            high: 1.501,
        },
        TOL,
    )
    .unwrap()
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
    let t = tool(0.38);
    let body = boolean::cut(&stock(true), &t, TOL).unwrap();
    let (same, history) = boolean::cut_recorded(&body, &t, TOL).unwrap();
    assert_eq!(write_body(&same), write_body(&body));
    assert!(history.unwrap_err().contains("UnchangedRadialNotch"));
    assert!(
        matches!(boolean::cut_named(&body,&t,TOL),Err(gc_kernel::KernelError::Unsupported(e)) if e.contains("UnchangedRadialNotch"))
    );
}

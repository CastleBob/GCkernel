//! Independent dimension-based oracle, not a mass result reused as an expected value.
use gc_kernel::{
    boolean,
    geom::{Curve, Surface, Vec3},
    mass,
    neutral::{write_body, Body, Profile, Reader},
    revolve::{Kind, Revolution, Segment},
    tessellate, topology,
};
use std::f64::consts::{PI, TAU};
const TOL: f64 = 1e-6;
fn sector(
    ri: f64,
    ro: f64,
    lo: f64,
    hi: f64,
    angle: f64,
    center: f64,
    origin: Vec3,
    axis: Vec3,
    seam: Vec3,
) -> Body {
    let start = center - angle / 2.;
    let x = seam
        .scale(start.cos())
        .add(axis.cross(seam).scale(start.sin()));
    gc_kernel::revolve::revolve_general(
        &Revolution {
            p: origin,
            a: axis,
            x,
            angle,
            loops: vec![[(ri, lo), (ro, lo), (ro, hi), (ri, hi)]
                .into_iter()
                .map(|(r, z)| Segment {
                    r,
                    z,
                    kind: Kind::Line,
                })
                .collect()],
        },
        TOL,
    )
    .unwrap()
}
fn canonical(ri: f64, ro: f64, lo: f64, hi: f64, angle: f64) -> Body {
    sector(
        ri,
        ro,
        lo,
        hi,
        angle,
        0.,
        Vec3::new(0., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    )
}
fn closed_mesh(body: &Body, expected: f64, scale: f64) {
    let mesh = tessellate::tessellate(body, 0.01 * scale).unwrap();
    use std::collections::HashMap;
    let mut ids = HashMap::new();
    let mut remap = Vec::new();
    for p in &mesh.positions {
        let key = [p.x, p.y, p.z].map(|v| if v == 0. { 0 } else { v.to_bits() });
        let count = ids.len();
        remap.push(*ids.entry(key).or_insert(count));
    }
    let mut edges: HashMap<(usize, usize), (usize, i32)> = HashMap::new();
    let mut volume = 0.;
    for [a, b, c] in &mesh.triangles {
        volume += mesh.positions[*a].dot(mesh.positions[*b].cross(mesh.positions[*c])) / 6.;
        for (u, v) in [(*a, *b), (*b, *c), (*c, *a)] {
            let (u, v) = (remap[u], remap[v]);
            assert_ne!(u, v);
            let value = edges.entry((u.min(v), u.max(v))).or_default();
            value.0 += 1;
            value.1 += if u < v { 1 } else { -1 };
        }
    }
    assert!(
        edges.values().all(|v| *v == (2, 0)),
        "mesh has an open or inconsistently oriented edge"
    );
    assert!(
        volume > 0.
            && (volume - expected).abs() < mass::measure(body, TOL).unwrap().area * mesh.chord,
        "mesh {volume} vs {expected}"
    );
}

#[test]
fn seven_sweep_sizes_in_both_directions_have_analytic_volume_and_closed_meshes() {
    let body = canonical(4., 5., 0., 10., TAU);
    let before = write_body(&body);
    for degrees in [0.01_f64, 30., 90., 180., 270., 350., 359.99] {
        for sign in [-1., 1.] {
            let angle = degrees.to_radians() * sign;
            let tool = canonical(3.9999, 5.0001, 3.9999, 6.0001, angle);
            let tool_before = write_body(&tool);
            let result = boolean::cut(&body, &tool, TOL)
                .unwrap_or_else(|e| panic!("{degrees} * {sign}: {e}"));
            let expected = 90. * PI - 0.5 * angle.abs() * 9. * 2.0002;
            let report = topology::check(&result, TOL).unwrap();
            assert_eq!((report.faces, report.edges, report.vertices), (8, 18, 12));
            assert!((mass::measure(&result, TOL).unwrap().volume - expected).abs() < 1e-8);
            let area = TAU * 9. * 10. + TAU * 9. - 9. * angle.abs() * 2.0002
                + angle.abs() * 9.
                + 2. * 2.0002;
            assert!(
                (mass::measure(&result, TOL).unwrap().area - area).abs() < 1e-8,
                "surface area at {degrees}"
            );
            // Check the physical clock angles, which a volume-only oracle cannot distinguish.
            for e in &result.edges {
                if let Curve::Arc { c, n, x, r, t0, t1 } = e.curve {
                    if t1 - t0 < TAU {
                        for t in [t0, t1] {
                            let direction = x.scale(t.cos()).add(n.cross(x).scale(t.sin()));
                            let expected_directions = [-angle.abs() / 2., angle.abs() / 2.]
                                .map(|u| Vec3::new(u.sin(), 0., u.cos()));
                            assert!(
                                expected_directions
                                    .iter()
                                    .any(|d| d.distance(direction) < 1e-12),
                                "window is at the wrong clock angle: {direction:?}"
                            );
                            assert!((c.y - 3.9999).abs() < 1e-12 || (c.y - 6.0001).abs() < 1e-12);
                            assert!(r == 4. || r == 5.);
                        }
                    }
                }
            }
            closed_mesh(&result, expected, 1.);
            let wire = write_body(&result);
            let again = Reader::with_tolerance(&wire, TOL).body().unwrap();
            topology::check(&again, TOL).unwrap();
            assert!((mass::measure(&again, TOL).unwrap().volume - expected).abs() < 1e-8);
            assert_eq!(before, write_body(&body));
            assert_eq!(tool_before, write_body(&tool));
        }
    }
}

#[test]
fn recorded_cut_reports_missing_history_and_strict_named_cut_refuses() {
    let body = canonical(4., 5., 0., 10., TAU);
    let tool = canonical(3.9, 5.1, 4., 6., PI / 2.);
    let ordinary = boolean::cut(&body, &tool, TOL).unwrap();
    let (recorded, history) = boolean::cut_recorded(&body, &tool, TOL).unwrap();
    assert_eq!(write_body(&ordinary), write_body(&recorded));
    assert!(history.unwrap_err().contains("AnnularWindow"));
    assert!(boolean::cut_named(&body, &tool, TOL)
        .unwrap_err()
        .to_string()
        .contains("AnnularWindow"));
}

#[test]
fn end_notches_on_both_ends_keep_analytic_volume_area_and_closed_meshes() {
    let body = canonical(4., 5., 0., 10., TAU);
    let before = write_body(&body);
    for degrees in [0.01_f64, 30., 90., 180., 270., 350., 359.99] {
        for sign in [-1., 1.] {
            for upper in [false, true] {
                let (lo, hi) = if upper {
                    (7.9999, 12.0001)
                } else {
                    (-2.0001, 2.0001)
                };
                let tool = canonical(3.9999, 5.0001, lo, hi, degrees.to_radians() * sign);
                let tool_before = write_body(&tool);
                let result = boolean::cut(&body, &tool, TOL)
                    .unwrap_or_else(|e| panic!("{degrees} * {sign}, upper {upper}: {e}"));
                let theta = degrees.to_radians();
                let h = 2.0001;
                let expected = 90. * PI - 0.5 * theta * 9. * h;
                let area = TAU * 9. * 10. + TAU * 9. - 9. * theta * h + 2. * h;
                let props = mass::measure(&result, TOL).unwrap();
                assert!(
                    (props.volume - expected).abs() < 1e-8,
                    "{degrees}: {} vs {expected}",
                    props.volume
                );
                assert!(
                    (props.area - area).abs() < 1e-8,
                    "{degrees}: area {} vs {area}",
                    props.area
                );
                let report = topology::check(&result, TOL).unwrap();
                assert_eq!((report.faces, report.edges, report.vertices), (7, 18, 12));
                closed_mesh(&result, expected, 1.);
                let wire = write_body(&result);
                let round = Reader::with_tolerance(&wire, TOL).body().unwrap();
                topology::check(&round, TOL).unwrap();
                assert!((mass::measure(&round, TOL).unwrap().volume - expected).abs() < 1e-8);
                assert_eq!(write_body(&body), before);
                assert_eq!(write_body(&tool), tool_before);
            }
        }
    }
}

#[test]
fn end_notches_preserve_rotated_positions_cardinal_frames_and_scale() {
    for scale in [0.1, 1., 25.4] {
        for upper in [false, true] {
            let o = Vec3::new(7., 11., 13.).scale(scale);
            let a = Vec3::new(0., 0., -1.);
            let x = Vec3::new(1., 0., 0.);
            let body = sector(4. * scale, 5. * scale, 0., 10. * scale, TAU, 0., o, a, x);
            let (lo, hi) = if upper { (8., 13.) } else { (-3., 2.) };
            let tool = sector(
                3.9 * scale,
                5.1 * scale,
                lo * scale,
                hi * scale,
                1.5 * PI,
                0.73,
                o,
                a,
                x,
            );
            let result = boolean::cut(&body, &tool, TOL).unwrap();
            let expected = (90. * PI - 0.5 * 1.5 * PI * 9. * 2.) * scale.powi(3);
            assert!(
                (mass::measure(&result, TOL).unwrap().volume - expected).abs()
                    < 1e-8 * scale.powi(3)
            );
            closed_mesh(&result, expected, scale);
            // Every new radial face lies at an authored boundary angle, and its
            // cut reaches the original stock end rather than the tool's far cap.
        for f in &result.faces {
            if let Surface::Plane{n,..}=f.surface {if n.dot(a)==0. {
                assert!([0.73-0.75*PI,0.73+0.75*PI].into_iter().any(|t| {
                    let radial=x.scale(t.cos()).add(a.cross(x).scale(t.sin()));
                    let normal=a.cross(radial);
                    n.distance(normal)<1e-12 || n.distance(normal.scale(-1.))<1e-12
                }));
            }}
        }
            for p in &result.vertices {
                let h = p.sub(o).dot(a) / scale;
                assert!(h >= -1e-12 && h <= 10. + 1e-12);
            }
            let (again, history) = boolean::cut_recorded(&body, &tool, TOL).unwrap();
            assert_eq!(write_body(&result), write_body(&again));
            assert!(history.unwrap_err().contains("AnnularWindow"));
        }
    }
}

#[test]
fn extruded_pipe_representation_is_admitted_and_malformed_supports_refuse() {
    let ring = |r| {
        vec![Curve::Arc {
            c: Vec3::new(0., 0., 0.),
            n: Vec3::new(0., -1., 0.),
            x: Vec3::new(1., 0., 0.),
            r,
            t0: 0.,
            t1: TAU,
        }]
    };
    let body = gc_kernel::extrude::extrude(
        &Profile {
            loops: vec![ring(5.), ring(4.)],
            n: Vec3::new(0., 1., 0.),
            low: 0.,
            high: 10.,
        },
        TOL,
    )
    .unwrap();
    let tool = canonical(3.9, 5.1, 4., 6., PI / 2.);
    let result = boolean::cut(&body, &tool, TOL).unwrap();
    assert!((mass::measure(&result, TOL).unwrap().volume - (90. * PI - 4.5 * PI)).abs() < 1e-8);
    let mut bad = tool.clone();
    let plane = bad
        .faces
        .iter_mut()
        .find(|f| matches!(f.surface,Surface::Plane{n,..} if n.y==1.))
        .unwrap();
    if let Surface::Plane { p, .. } = &mut plane.surface {
        p.y += TOL / 4.;
    }
    // This is within the generic topology tolerance, but not the domain's exact recipe.
    topology::check(&bad, TOL).unwrap();
    assert!(boolean::cut(&body, &bad, TOL)
        .unwrap_err()
        .to_string()
        .contains("exact ring plane"));
    let mut inside = tool.clone();
    let wall = inside
        .faces
        .iter_mut()
        .find(|f| matches!(f.surface, Surface::Cylinder { .. }) && f.sense)
        .unwrap();
    wall.sense = false;
    assert!(boolean::cut(&body, &inside, TOL)
        .unwrap_err()
        .to_string()
        .contains("outward senses"));
}

#[test]
fn translated_cardinal_frames_units_and_rotated_windows() {
    for scale in [0.1, 1., 25.4] {
        for (axis, seam) in [
            (Vec3::new(1., 0., 0.), Vec3::new(0., 1., 0.)),
            (Vec3::new(0., 0., -1.), Vec3::new(1., 0., 0.)),
        ] {
            let o = Vec3::new(7., 11., 13.).scale(scale);
            let body = sector(
                4. * scale,
                5. * scale,
                0.,
                10. * scale,
                TAU,
                0.,
                o,
                axis,
                seam,
            );
            let tool = sector(
                3.9 * scale,
                5.1 * scale,
                4. * scale,
                6. * scale,
                1.5 * PI,
                0.73,
                o,
                axis,
                seam,
            );
            let result = boolean::cut(&body, &tool, TOL).unwrap();
            let expected = (90. * PI - 0.5 * 1.5 * PI * 9. * 2.) * scale.powi(3);
            assert!(
                (mass::measure(&result, TOL).unwrap().volume - expected).abs()
                    < 1e-8 * scale.powi(3)
            );
            closed_mesh(&result, expected, scale);
        }
    }
}

#[test]
fn insufficient_reach_touching_ends_and_previous_windows_refuse_without_mutation() {
    let body = canonical(4., 5., 0., 10., TAU);
    let saved = write_body(&body);
    for (ri, ro, lo, hi) in [
        (4., 5.1, 4., 6.),
        (3.9, 5., 4., 6.),
        (3.9, 5.1, 0., 6.),
        (3.9, 5.1, 4., 10.),
        (3.9, 5.1, -1., 11.),
    ] {
        let tool = canonical(ri, ro, lo, hi, PI / 2.);
        let err = boolean::cut(&body, &tool, TOL).unwrap_err().to_string();
        assert!(err.contains("annular pipe window"), "{err}");
        assert_eq!(saved, write_body(&body));
    }
    let tool = canonical(3.9, 5.1, 4., 6., PI / 2.);
    let result = boolean::cut(&body, &tool, TOL).unwrap();
    let saved = write_body(&result);
    assert!(boolean::cut(&result, &tool, TOL)
        .unwrap_err()
        .to_string()
        .contains("plain pipe"));
    assert_eq!(saved, write_body(&result));
}

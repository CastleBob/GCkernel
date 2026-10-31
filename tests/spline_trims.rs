use gc_kernel::{
    extrude,
    frame::{Point3, Rigid},
    geom::{Curve, Surface, Vec3},
    mass,
    neutral::{self, Body, Profile},
    spline_surface::{Patch, SplineSurface},
    tessellate, topology,
};
const TOL: f64 = 1e-6;
#[test]
fn exported_trimmed_roofs_roundtrip_and_can_supply_native_reference_files() {
    for hole in [false, true] {
        for multi in [false, true] {
            let b = solid(hole, multi);
            let text = gc_kernel::step::write(&[("Trimmed roof", &b)], TOL).unwrap();
            let restored = gc_kernel::step::read(text.as_bytes(), TOL)
                .unwrap()
                .bodies
                .remove(0)
                .body;
            assert!(
                (mass::measure(&restored, TOL).unwrap().volume
                    - mass::measure(&b, TOL).unwrap().volume)
                    .abs()
                    < 1e-9
            );
            topology::check(&restored, TOL).unwrap();
            tessellate::tessellate(&restored, 0.02).unwrap();
            if let Ok(directory) = std::env::var("GC_STEP_EXPORT_FIXTURES") {
                std::fs::write(
                    std::path::Path::new(&directory).join(format!(
                        "TrimmedRoof-{}-{}.step",
                        if hole { "hole" } else { "solid" },
                        if multi { "multi" } else { "single" }
                    )),
                    text,
                )
                .unwrap();
            }
        }
    }
}
fn rectangle(a: f64, b: f64, c: f64, d: f64, reverse: bool) -> Vec<Curve> {
    let mut p = vec![
        Vec3::new(a, c, 0.),
        Vec3::new(b, c, 0.),
        Vec3::new(b, d, 0.),
        Vec3::new(a, d, 0.),
    ];
    if reverse {
        p.reverse();
    }
    (0..4)
        .map(|i| Curve::Line {
            a: p[i],
            b: p[(i + 1) % 4],
        })
        .collect()
}
fn solid(hole: bool, multispan: bool) -> Body {
    let mut loops = vec![rectangle(0.2, 1.8, 0.3, 2.7, false)];
    if hole {
        loops.push(rectangle(0.8, 1.2, 0.9, 2.1, true));
    }
    let mut body = extrude::extrude(
        &Profile {
            loops,
            n: Vec3::new(0., 0., 1.),
            low: 0.,
            high: 4.,
        },
        TOL,
    )
    .unwrap();
    for e in &mut body.edges {
        let a = e.curve.start();
        let b = e.curve.end();
        if a.z == 4. && b.z == 4. {
            let dx = b.x - a.x;
            e.curve = Curve::Bezier {
                p: std::array::from_fn(|i| {
                    let t = i as f64 / 3.;
                    Vec3::new(
                        a.x + dx * t,
                        a.y + (b.y - a.y) * t,
                        4. + (a.x * a.x
                            + 2. * a.x * dx * t
                            + dx * dx * (i * (i.saturating_sub(1))) as f64 / 6.)
                            / 4.,
                    )
                }),
            };
        } else {
            e.curve = Curve::Line {
                a: lift(a),
                b: lift(b),
            };
        }
    }
    for v in &mut body.vertices {
        *v = lift(*v);
    }
    let patch = Patch {
        p: std::array::from_fn(|k| {
            let (i, j) = (k / 4, k % 4);
            Vec3::new(
                2. * i as f64 / 3.,
                j as f64,
                4. + (i * (i.saturating_sub(1))) as f64 / 6.,
            )
        }),
        w: [1.; 16],
    };
    let data = if multispan {
        SplineSurface {
            u: vec![0., 0.5, 1.],
            v: vec![0., 0.4, 1.],
            patches: vec![
                patch.piece(0., 0.5, 0., 0.4),
                patch.piece(0., 0.5, 0.4, 1.),
                patch.piece(0.5, 1., 0., 0.4),
                patch.piece(0.5, 1., 0.4, 1.),
            ],
        }
    } else {
        SplineSurface {
            u: vec![0., 1.],
            v: vec![0., 1.],
            patches: vec![patch],
        }
    };
    data.validate().unwrap();
    let face = body
        .faces
        .iter_mut()
        .find(|f| matches!(f.surface,Surface::Plane{n,..} if n.z>0.99 && f.sense))
        .unwrap();
    face.surface = Surface::Spline {
        data: Box::new(data),
    };
    body
}
fn lift(p: Vec3) -> Vec3 {
    if p.z == 4. {
        Vec3::new(p.x, p.y, 4. + p.x * p.x / 4.)
    } else {
        p
    }
}
fn volume(a: f64, b: f64, c: f64, d: f64) -> f64 {
    (d - c) * (4. * (b - a) + (b.powi(3) - a.powi(3)) / 12.)
}
fn horizontal_area(a: f64, b: f64, c: f64, d: f64) -> f64 {
    let f = |x: f64| x / 2. * (1. + x * x / 4.).sqrt() + (x / 2.).asinh();
    (d - c) * (b - a + f(b) - f(a))
}
fn wall_area(a: f64, b: f64, c: f64, d: f64) -> f64 {
    2. * (4. * (b - a) + (b.powi(3) - a.powi(3)) / 12.) + (d - c) * (8. + (a * a + b * b) / 4.)
}
#[test]
fn trims_and_holes_preserve_mass_shared_mesh_and_rigid_frames() {
    for hole in [false, true] {
        for multi in [false, true] {
            let body = solid(hole, multi);
            let expected =
                volume(0.2, 1.8, 0.3, 2.7) - if hole { volume(0.8, 1.2, 0.9, 2.1) } else { 0. };
            let expected_area = horizontal_area(0.2, 1.8, 0.3, 2.7)
                + wall_area(0.2, 1.8, 0.3, 2.7)
                + if hole {
                    wall_area(0.8, 1.2, 0.9, 2.1) - horizontal_area(0.8, 1.2, 0.9, 2.1)
                } else {
                    0.
                };
            let rigid = Rigid::rotation(Vec3::new(0., 1., 0.), 0.37, Point3::new(0., 0., 0.))
                .unwrap()
                .then(&Rigid::translation(Vec3::new(20., -10., 3.)).unwrap());
            for b in [body.clone(), rigid.body(&body)] {
                topology::check(&b, TOL).unwrap();
                let m = mass::measure(&b, TOL).unwrap();
                assert!(
                    (m.volume - expected).abs() < 1e-9,
                    "{} vs {expected}",
                    m.volume
                );
                assert!(
                    (m.area - expected_area).abs() < 1e-9,
                    "area {} vs {expected_area}",
                    m.area
                );
                let data = neutral::write_body(&b);
                let restored = neutral::Reader::with_tolerance(&data, TOL).body().unwrap();
                topology::check(&restored, TOL).unwrap();
                let mesh = tessellate::tessellate(&b, 0.02).unwrap();
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
                    (volume - expected).abs() < 0.05,
                    "display {volume} vs {expected}"
                );
                for line in &mesh.edges {
                    for pair in line.windows(2) {
                        let count = mesh
                            .triangles
                            .iter()
                            .flat_map(|t| (0..3).map(move |k| (t[k], t[(k + 1) % 3])))
                            .filter(|&(a, b)| {
                                (mesh.positions[a] == pair[0] && mesh.positions[b] == pair[1])
                                    || (mesh.positions[b] == pair[0]
                                        && mesh.positions[a] == pair[1])
                            })
                            .count();
                        assert_eq!(count, 2, "unshared {pair:?}");
                    }
                }
            }
        }
    }
}
#[test]
fn mismatches_winding_outside_and_crossing_loops_are_refused() {
    let original = solid(true, false);
    let face = original
        .faces
        .iter()
        .position(|f| matches!(f.surface, Surface::Spline { .. }))
        .unwrap();
    for kind in 0..4 {
        let mut b = original.clone();
        match kind {
            0 => {
                let e = b.faces[face].loops[0][0].0;
                if let Curve::Bezier { p } = &mut b.edges[e].curve {
                    p[1].z += 1e-3;
                }
            }
            1 => {
                b.faces[face].sense = false;
            }
            2 => {
                b.faces[face].loops[1].reverse();
            }
            _ => {
                let lp = b.faces[face].loops[1].clone();
                b.faces[face].loops.push(lp);
            }
        }
        assert!(topology::check(&b, TOL).is_err(), "kind {kind}");
    }
    topology::check(&original, TOL).unwrap();
}

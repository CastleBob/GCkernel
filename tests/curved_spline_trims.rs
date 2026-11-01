use gc_kernel::{
    extrude,
    geom::{Curve, Surface, Vec3},
    mass,
    neutral::{self, Body, Edge, Face, Profile},
    spline_surface::{Patch, SplineSurface},
    step, tessellate, topology,
};
const TOL: f64 = 1e-6;

#[test]
fn curved_multispan_trims_keep_mass_shared_mesh_and_step_geometry() {
    for (name, mut body) in [
        ("CurvedMultispanRoof.step", curved_roof()),
        ("CircularMultispanHole.step", circular_hole()),
    ] {
        let expected = mass::measure(&body, TOL).unwrap();
        for face in &mut body.faces {
            if let Surface::Spline { data } = &mut face.surface {
                let original = data.patches[0].clone();
                data.u = vec![0., 0.35, 0.5, 0.8, 1.];
                data.v = vec![0., 0.05, 0.4, 0.5, 0.6, 1.];
                data.patches = data
                    .u
                    .windows(2)
                    .flat_map(|u| {
                        data.v
                            .windows(2)
                            .map(|v| original.piece(u[0], u[1], v[0], v[1]))
                    })
                    .collect();
            }
        }
        topology::check(&body, TOL).unwrap();
        let measured = mass::measure(&body, TOL).unwrap();
        assert!((measured.volume - expected.volume).abs() < 1e-9);
        assert!((measured.area - expected.area).abs() < 1e-8);
        let mesh = tessellate::tessellate(&body, 0.01).unwrap();
        assert!(!mesh.triangles.is_empty());
        let bytes = step::write(&[(name, &body)], TOL).unwrap();
        let round = step::read(bytes.as_bytes(), TOL).unwrap();
        assert!((round.bodies[0].volume - expected.volume).abs() < 1e-9);
        if let Ok(dir) = std::env::var("GC_STEP_EXPORT_FIXTURES") {
            std::fs::write(std::path::Path::new(&dir).join(name), bytes).unwrap();
        }
        // A defective interior patch must not escape whole-edge checks just
        // because the trim's endpoints lie on other patches.
        let face = body
            .faces
            .iter_mut()
            .find(|f| matches!(f.surface, Surface::Spline { .. }))
            .unwrap();
        if let Surface::Spline { data } = &mut face.surface {
            data.patches[6].p[5].z += 0.01;
        }
        assert!(topology::check(&body, TOL).is_err());
    }
}
fn graph() -> SplineSurface {
    SplineSurface {
        u: vec![0., 1.],
        v: vec![0., 1.],
        patches: vec![Patch {
            p: std::array::from_fn(|k| {
                let (i, j) = (k / 4, k % 4);
                Vec3::new(
                    2. * i as f64 / 3.,
                    j as f64,
                    4. + (i * i.saturating_sub(1)) as f64 / 6.,
                )
            }),
            w: [1.; 16],
        }],
    }
}
fn curved_roof() -> Body {
    let corners = [
        Vec3::new(0.2, 0.3, 0.),
        Vec3::new(1.8, 0.3, 0.),
        Vec3::new(1.8, 2.7, 0.),
        Vec3::new(0.2, 2.7, 0.),
    ];
    let bottom: [Curve; 4] = std::array::from_fn(|k| {
        let (a, b) = (corners[k], corners[(k + 1) % 4]);
        Curve::Bezier {
            p: std::array::from_fn(|i| {
                let mut p = a.add(b.sub(a).scale(i as f64 / 3.));
                if k == 0 {
                    p.y += 0.6 * (-(i as f64) / 3. + (i * i.saturating_sub(1)) as f64 / 6.);
                }
                p
            }),
        }
    });
    let top: [Curve; 4] = std::array::from_fn(|k| {
        let Curve::Bezier { p } = bottom[k] else {
            unreachable!()
        };
        let (a, d) = (p[0].x, p[3].x - p[0].x);
        Curve::Bezier {
            p: std::array::from_fn(|i| {
                Vec3::new(
                    p[i].x,
                    p[i].y,
                    4. + (a * a
                        + 2. * a * d * i as f64 / 3.
                        + d * d * (i * i.saturating_sub(1)) as f64 / 6.)
                        / 4.,
                )
            }),
        }
    });
    let mut b = Body::default();
    b.vertices.extend(corners);
    b.vertices.extend(top.iter().map(Curve::start));
    for i in 0..4 {
        b.edges.push(Edge {
            curve: bottom[i].clone(),
            start: i,
            end: (i + 1) % 4,
        });
    }
    for i in 0..4 {
        b.edges.push(Edge {
            curve: top[i].clone(),
            start: i + 4,
            end: (i + 1) % 4 + 4,
        });
    }
    for i in 0..4 {
        b.edges.push(Edge {
            curve: Curve::Line {
                a: b.vertices[i],
                b: b.vertices[i + 4],
            },
            start: i,
            end: i + 4,
        });
    }
    b.faces.push(Face {
        surface: Surface::Plane {
            p: Vec3::new(0., 0., 0.),
            n: Vec3::new(0., 0., 1.),
            x: Vec3::new(1., 0., 0.),
        },
        sense: false,
        loops: vec![(0..4).rev().map(|i| (i, false)).collect()],
    });
    b.faces.push(Face {
        surface: Surface::Spline {
            data: Box::new(graph()),
        },
        sense: true,
        loops: vec![(4..8).map(|i| (i, true)).collect()],
    });
    for i in 0..4 {
        let Curve::Bezier { p: low } = bottom[i] else {
            unreachable!()
        };
        let Curve::Bezier { p: high } = top[i] else {
            unreachable!()
        };
        let data = SplineSurface {
            u: vec![0., 1.],
            v: vec![0., 1.],
            patches: vec![Patch {
                p: std::array::from_fn(|k| {
                    low[k / 4].add(high[k / 4].sub(low[k / 4]).scale((k % 4) as f64 / 3.))
                }),
                w: [1.; 16],
            }],
        };
        data.validate().unwrap();
        b.faces.push(Face {
            surface: Surface::Spline {
                data: Box::new(data),
            },
            sense: true,
            loops: vec![vec![
                (i, true),
                (8 + (i + 1) % 4, true),
                (4 + i, false),
                (8 + i, false),
            ]],
        });
    }
    b
}
fn circular_hole() -> Body {
    let p = [
        Vec3::new(0.2, 0.3, 0.),
        Vec3::new(1.8, 0.3, 0.),
        Vec3::new(1.8, 2.7, 0.),
        Vec3::new(0.2, 2.7, 0.),
    ];
    let loops = vec![
        (0..4)
            .map(|i| Curve::Line {
                a: p[i],
                b: p[(i + 1) % 4],
            })
            .collect(),
        vec![Curve::Arc {
            c: Vec3::new(1., 1.5, 0.),
            n: Vec3::new(0., 0., -1.),
            x: Vec3::new(1., 0., 0.),
            r: 0.3,
            t0: 0.,
            t1: std::f64::consts::TAU,
        }],
    ];
    let mut b = extrude::extrude(
        &Profile {
            loops,
            n: Vec3::new(0., 0., 1.),
            low: 0.,
            high: 4.,
        },
        TOL,
    )
    .unwrap();
    let mut data = graph();
    for patch in &mut data.patches {
        for p in &mut patch.p {
            p.z = 4.;
        }
    }
    b.faces
        .iter_mut()
        .find(|f| matches!(f.surface,Surface::Plane{n,..} if n.z>0.99 && f.sense))
        .unwrap()
        .surface = Surface::Spline {
        data: Box::new(data),
    };
    b
}
#[test]
fn curved_uv_roof_and_round_hole_preserve_analytic_solids() {
    let simpson = |f: &dyn Fn(f64) -> f64| {
        let n = 2048;
        (0..=n)
            .map(|i| {
                f(i as f64 / n as f64)
                    * if i == 0 || i == n {
                        1.
                    } else if i % 2 == 0 {
                        2.
                    } else {
                        4.
                    }
            })
            .sum::<f64>()
            / (3. * n as f64)
    };
    let roof_area = 4.
        + 2.4 * (4. + 0.01 + 4. + 0.81)
        + (4. * 1.6 + (1.8f64.powi(3) - 0.2f64.powi(3)) / 12.)
        + simpson(&|t| {
            let x = 0.2 + 1.6 * t;
            1.6 * (1. + x * x / 4.).sqrt() * (2.4 + 0.6 * t * (1. - t))
        })
        + simpson(&|t| {
            let x = 0.2 + 1.6 * t;
            (4. + x * x / 4.) * (1.6f64.powi(2) + 0.36 * (2. * t - 1.).powi(2)).sqrt()
        });
    let roof_volume = 2.4 * (4. * 1.6 + (1.8f64.powi(3) - 0.2f64.powi(3)) / 12.)
        + 1.6 * 0.6 * (4. / 6. + 0.2 * 0.2 / 24. + 0.2 * 1.6 / 24. + 1.6 * 1.6 / 80.);
    for (name, b, expected) in [
        ("CurvedUvRoof.step", curved_roof(), roof_volume),
        (
            "CircularTrimSplineFace.step",
            circular_hole(),
            4. * (1.6 * 2.4 - std::f64::consts::PI * 0.09),
        ),
    ] {
        topology::check(&b, TOL).unwrap();
        let m = mass::measure(&b, TOL).unwrap();
        let expected_area = if name == "CurvedUvRoof.step" {
            roof_area
        } else {
            2. * (1.6 * 2.4 - std::f64::consts::PI * 0.09) + 4. * (8. + std::f64::consts::TAU * 0.3)
        };
        assert!(
            (m.area - expected_area).abs() < 1e-8,
            "area {} vs {expected_area}",
            m.area
        );
        assert!(
            (m.volume - expected).abs() < 1e-9,
            "{} vs {expected}",
            m.volume
        );
        let data = neutral::write_body(&b);
        let restored = neutral::Reader::with_tolerance(&data, TOL).body().unwrap();
        topology::check(&restored, TOL).unwrap();
        let mesh = tessellate::tessellate(&b, 0.02).unwrap();
        for line in &mesh.edges {
            for pair in line.windows(2) {
                let count = mesh
                    .triangles
                    .iter()
                    .flat_map(|t| (0..3).map(move |i| (t[i], t[(i + 1) % 3])))
                    .filter(|&(a, b)| {
                        (mesh.positions[a] == pair[0] && mesh.positions[b] == pair[1])
                            || (mesh.positions[a] == pair[1] && mesh.positions[b] == pair[0])
                    })
                    .count();
                assert_eq!(count, 2);
            }
        }
        let text = step::write(&[(name, &b)], TOL).unwrap();
        let restored = step::read(text.as_bytes(), TOL)
            .unwrap()
            .bodies
            .remove(0)
            .body;
        assert!((mass::measure(&restored, TOL).unwrap().volume - expected).abs() < 1e-9);
        if let Ok(dir) = std::env::var("GC_STEP_EXPORT_FIXTURES") {
            std::fs::write(std::path::Path::new(&dir).join(name), text).unwrap();
        }
    }
}

#[test]
fn curved_loop_crossings_touches_and_support_mismatches_are_refused() {
    let mut data = graph();
    for patch in &mut data.patches {
        for p in &mut patch.p {
            p.z = 4.;
        }
    }
    let curve = Curve::Bezier {
        p: [(0.2, 0.2), (0.95, 0.8), (0.05, 0.8), (0.8, 0.2)]
            .map(|(u, v)| Vec3::new(2. * u, 3. * v, 4.)),
    };
    let a = curve.start();
    let b = curve.end();
    let face = Face {
        surface: Surface::Spline {
            data: Box::new(data.clone()),
        },
        sense: true,
        loops: vec![vec![(0, true), (1, true)]],
    };
    let local = Body {
        vertices: vec![a, b],
        edges: vec![
            Edge {
                curve,
                start: 0,
                end: 1,
            },
            Edge {
                curve: Curve::Line { a: b, b: a },
                start: 1,
                end: 0,
            },
        ],
        faces: vec![face.clone()],
    };
    assert!(data.check_face(&local, &face, TOL).is_err());
    let mut touching = circular_hole();
    for e in &mut touching.edges {
        if let Curve::Arc { r, .. } = &mut e.curve {
            *r = 0.8;
        }
    }
    let face = touching
        .faces
        .iter()
        .find(|f| matches!(f.surface, Surface::Spline { .. }))
        .unwrap();
    let Surface::Spline { data } = &face.surface else {
        unreachable!()
    };
    assert!(data.check_face(&touching, face, TOL).is_err());
    let mut mismatch = curved_roof();
    if let Curve::Bezier { p } = &mut mismatch.edges[4].curve {
        p[1].z += 0.001;
    }
    let face = &mismatch.faces[1];
    let Surface::Spline { data } = &face.surface else {
        unreachable!()
    };
    assert!(data.check_face(&mismatch, face, TOL).is_err());
    topology::check(&curved_roof(), TOL).unwrap();
}

#[test]
fn curved_trim_mass_and_display_survive_rigid_frames() {
    use gc_kernel::frame::{Point3, Rigid};
    let transform = Rigid::rotation(Vec3::new(0., 1., 0.), 0.4, Point3::new(0., 0., 0.))
        .unwrap()
        .then(&Rigid::translation(Vec3::new(20., -30., 40.)).unwrap());
    for body in [curved_roof(), circular_hole()] {
        let original = mass::measure(&body, TOL).unwrap();
        let moved = transform.body(&body);
        topology::check(&moved, TOL).unwrap();
        let m = mass::measure(&moved, TOL).unwrap();
        assert!((m.volume - original.volume).abs() < 1e-8);
        assert!((m.area - original.area).abs() < 1e-8);
        tessellate::tessellate(&moved, 0.02).unwrap();
    }
}

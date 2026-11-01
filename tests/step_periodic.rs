use gc_kernel::{
    geom::{Curve, Surface, Vec3},
    mass,
    neutral::{Body, Edge, Face},
    spline_surface::{Patch, SplineSurface},
    step, tessellate, topology,
};
const TOL: f64 = 1e-6;
fn periodic_body() -> Body {
    let corners = [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(-1., 0., 0.),
        Vec3::new(0., -1., 0.),
    ];
    let curves: Vec<_> = (0..4)
        .map(|i| {
            let a = corners[i].add(corners[(i + 1) % 4]).scale(0.5);
            let b = corners[(i + 1) % 4];
            let c = b.add(corners[(i + 2) % 4]).scale(0.5);
            Curve::Bezier {
                p: [
                    a,
                    a.scale(1. / 3.).add(b.scale(2. / 3.)),
                    b.scale(2. / 3.).add(c.scale(1. / 3.)),
                    c,
                ],
            }
        })
        .collect();
    let mut body = Body::default();
    for z in [0., 3.] {
        for c in &curves {
            body.vertices.push(c.start().add(Vec3::new(0., 0., z)));
        }
    }
    for level in 0..2 {
        for i in 0..4 {
            let Curve::Bezier { p } = curves[i] else {
                unreachable!()
            };
            body.edges.push(Edge {
                curve: Curve::Bezier {
                    p: p.map(|p| p.add(Vec3::new(0., 0., 3. * level as f64))),
                },
                start: 4 * level + i,
                end: 4 * level + (i + 1) % 4,
            });
        }
    }
    body.edges.push(Edge {
        curve: Curve::Line {
            a: body.vertices[0],
            b: body.vertices[4],
        },
        start: 0,
        end: 4,
    });
    for level in 0..2 {
        body.faces.push(Face {
            surface: Surface::Plane {
                p: Vec3::new(0., 0., 3. * level as f64),
                n: Vec3::new(0., 0., 1.),
                x: Vec3::new(1., 0., 0.),
            },
            sense: level == 1,
            loops: vec![if level == 0 {
                (0..4).rev().map(|i| (i, false)).collect()
            } else {
                (4..8).map(|i| (i, true)).collect()
            }],
        });
    }
    let line = |u0, v0, u1, v1| Curve::Line {
        a: Vec3::new(u0, v0, 0.),
        b: Vec3::new(u1, v1, 0.),
    };
    let mut uv: Vec<_> = (0..4)
        .map(|i| line(i as f64 / 4., 0., (i + 1) as f64 / 4., 0.))
        .collect();
    uv.push(line(1., 0., 1., 1.));
    uv.extend(
        (0..4)
            .rev()
            .map(|i| line((i + 1) as f64 / 4., 1., i as f64 / 4., 1.)),
    );
    uv.push(line(0., 1., 0., 0.));
    let data = SplineSurface {
        u: vec![0., 0.25, 0.5, 0.75, 1.],
        v: vec![0., 1.],
        pcurves: Some(vec![uv]),
        patches: curves
            .iter()
            .map(|c| {
                let Curve::Bezier { p } = c else {
                    unreachable!()
                };
                Patch {
                    p: std::array::from_fn(|k| p[k / 4].add(Vec3::new(0., 0., (k % 4) as f64))),
                    w: [1.; 16],
                }
            })
            .collect(),
    };
    let mut lp: Vec<_> = (0..4).map(|i| (i, true)).collect();
    lp.push((8, true));
    lp.extend((4..8).rev().map(|i| (i, false)));
    lp.push((8, false));
    body.faces.push(Face {
        surface: Surface::Spline {
            data: Box::new(data),
        },
        sense: true,
        loops: vec![lp],
    });
    body
}
#[test]
fn closed_seam_and_unclamped_periodic_surface_preserve_the_original_solid() {
    let body = periodic_body();
    topology::check(&body, TOL).unwrap();
    assert!((mass::measure(&body, TOL).unwrap().volume - 5.).abs() < 1e-9);
    assert!(!tessellate::tessellate(&body, 0.01)
        .unwrap()
        .triangles
        .is_empty());
    let clamped = step::write(&[("Periodic spline prism", &body)], TOL).unwrap();
    let surface = clamped
        .lines()
        .find(|l| l.contains("B_SPLINE_SURFACE("))
        .unwrap();
    let id = surface.split('=').next().unwrap();
    let mut controls = String::new();
    let mut rows = Vec::new();
    let xy = [(1., 0.), (0., 1.), (-1., 0.), (0., -1.), (1., 0.), (0., 1.)];
    for (i, (x, y)) in xy.into_iter().enumerate() {
        for j in 0..2 {
            controls += &format!(
                "#{}=CARTESIAN_POINT('',({x},{y},{}));\n",
                900000 + i * 2 + j,
                3 * j
            );
        }
        rows.push(format!("(#{},#{})", 900000 + i * 2, 900001 + i * 2));
    }
    let definition=format!("{id}=B_SPLINE_SURFACE_WITH_KNOTS('',2,1,({}),.UNSPECIFIED.,.T.,.F.,.F.,(1,1,1,1,1,1,1,1,1),(2,2),(-0.5,-0.25,0.,0.25,0.5,0.75,1.,1.25,1.5),(0.,1.),.UNIFORM_KNOTS.);",rows.join(","));
    let periodic = clamped.replace(surface, &format!("{controls}{definition}"));
    for (name, text) in [
        ("ClosedSeamSpline.step", clamped),
        ("PeriodicSpline.step", periodic.clone()),
    ] {
        let imported = step::read(text.as_bytes(), TOL).unwrap();
        assert!((imported.bodies[0].volume - 5.).abs() < 1e-9);
        let mesh = tessellate::tessellate(&imported.bodies[0].body, 0.01).unwrap();
        let o = mesh.positions[0];
        let volume: f64 = mesh
            .triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| mesh.positions[i].sub(o));
                a.dot(b.cross(c)) / 6.
            })
            .sum();
        assert!((volume - 5.).abs() < 0.1);
        if let Ok(dir) = std::env::var("GC_STEP_EXPORT_FIXTURES") {
            std::fs::write(std::path::Path::new(&dir).join(name), text).unwrap();
        }
    }
    let broken = periodic.replace(
        "#900010=CARTESIAN_POINT('',(0,1,0))",
        "#900010=CARTESIAN_POINT('',(0,1.1,0))",
    );
    assert_ne!(broken, periodic);
    assert!(step::read(broken.as_bytes(), TOL)
        .unwrap_err()
        .to_string()
        .contains("does not close"));
}

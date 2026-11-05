use gc_kernel::{
    boolean,
    geom::{Curve, Vec3},
    mass,
    neutral::{write_body, Body, Profile, Reader},
    revolve::{revolve, Section},
    tessellate, topology,
};
use std::f64::consts::{PI, TAU};
const TOL: f64 = 1e-6;
fn stock(bevel: bool) -> Body {
    let h = 0.0625;
    let s = h * 30_f64.to_radians().tan();
    let mut points = vec![(1.5, 0.)];
    if bevel {
        points.extend([(1.5, 1.5 - s), (1.5 - h, 1.5)]);
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
fn tool(y: f64, rho: f64) -> Body {
    gc_kernel::extrude::extrude(
        &Profile {
            loops: vec![vec![Curve::Arc {
                c: Vec3::new(0., y, 0.),
                n: Vec3::new(0., 0., 1.),
                x: Vec3::new(1., 0., 0.),
                r: rho,
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
fn removed() -> f64 {
    let rho = 0.38;
    let n = 2048;
    let h = PI / n as f64;
    let f = |t: f64| {
        let x = rho * t.sin();
        rho * rho * t.cos().powi(2) * ((2.25 - x * x).sqrt() - (1.375_f64.powi(2) - x * x).sqrt())
    };
    let mut sum = f(-PI / 2.) + f(PI / 2.);
    for i in 1..n {
        sum += (if i % 2 == 0 { 2. } else { 4. }) * f(-PI / 2. + i as f64 * h);
    }
    sum * h / 3.
}
fn verify(out: &Body, want: f64) {
    topology::check(out, TOL).unwrap();
    let m = mass::measure(out, TOL).unwrap();
    assert!((m.volume - want).abs() < 1e-9, "{} vs {want}", m.volume);
    let wire = write_body(out);
    let mut r = Reader::with_tolerance(&wire, TOL);
    let back = r.body().unwrap();
    r.finish().unwrap();
    topology::check(&back, TOL).unwrap();
    let mesh = tessellate::tessellate(&back, 0.001).unwrap();
    let bits = |v: f64| if v == 0. { 0 } else { v.to_bits() };
    let key = |p: Vec3| (bits(p.x), bits(p.y), bits(p.z));
    let mut boundary = std::collections::HashMap::<_, i32>::new();
    let mut volume = 0.;
    let origin = mesh.positions[0];
    for tri in &mesh.triangles {
        let [a, b, c] = tri.map(|v| mesh.positions[v]);
        for (p, q) in [(a, b), (b, c), (c, a)] {
            if key(p) == key(q) {
                continue;
            }
            *boundary.entry((key(p), key(q))).or_default() += 1;
            *boundary.entry((key(q), key(p))).or_default() -= 1;
        }
        let normal = tri
            .iter()
            .fold(Vec3::new(0., 0., 0.), |n, &v| n.add(mesh.normals[v]));
        let cross = b.sub(a).cross(c.sub(a));
        if cross.length() > 1e-14 {
            assert!(cross.dot(normal) > 0.);
        }
        volume += a.sub(origin).dot(b.sub(origin).cross(c.sub(origin))) / 6.;
    }
    assert!(boundary.values().all(|v| *v == 0), "mesh is open");
    assert!((volume - want).abs() <= mesh.chord * m.area);
}
#[test]
fn centered_flat_end_notch_preserves_a_distant_bevel() {
    for bevel in [false, true] {
        let b = stock(bevel);
        let t = tool(0., 0.38);
        let original = (write_body(&b), write_body(&t));
        let out = boolean::cut(&b, &t, TOL).unwrap();
        let mut expected = PI * (2.25 - 1.375_f64.powi(2)) * 1.5 - removed();
        if bevel {
            let h = 0.0625;
            let s = h * 30_f64.to_radians().tan();
            expected -= PI * h * s * (1.5 - h / 3.);
        }
        verify(&out, expected);
        for (fi, face) in b.faces.iter().enumerate() {
            if matches!(face.surface, gc_kernel::geom::Surface::Cone { .. })
                || matches!(face.surface,gc_kernel::geom::Surface::Plane{p,..} if p.y==1.5)
            {
                assert_eq!(out.faces[fi].loops, face.loops);
                assert_eq!(
                    format!("{:?}", out.faces[fi].surface),
                    format!("{:?}", face.surface)
                );
                for &(ei, _) in face.loops.iter().flatten() {
                    assert_eq!(format!("{:?}", out.edges[ei]), format!("{:?}", b.edges[ei]));
                }
            }
        }
        assert_eq!(write_body(&b), original.0);
        assert_eq!(write_body(&t), original.1);
    }
}

#[test]
fn flat_notches_work_at_both_ends_and_in_cardinal_rigid_frames() {
    let b = stock(false);
    let wanted = mass::measure(&b, TOL).unwrap().volume - removed();
    let frames = [
        [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
        [[-1., 0., 0.], [0., -1., 0.], [0., 0., 1.]],
    ];
    for y in [0., 1.5] {
        for rows in frames {
            let rigid = gc_kernel::frame::Rigid::from_rows(rows, Vec3::new(8., -16., 32.)).unwrap();
            let moved = rigid.body(&b);
            let cutter = rigid.body(&tool(y, 0.38));
            let out = boolean::cut(&moved, &cutter, TOL).unwrap();
            verify(&out, wanted);
        }
    }
}

#[test]
fn existing_openings_partial_depth_seams_and_bevel_intersections_refuse_without_changes() {
    let b = stock(true);
    let t = tool(0., 0.38);
    let done = boolean::cut(&b, &t, TOL).unwrap();
    let mid = boolean::cut(&b, &tool(0.75, 0.1), TOL).unwrap();
    let near_bevel = revolve(
        &Section {
            p: Vec3::new(0., 0., 0.),
            a: Vec3::new(0., 1., 0.),
            x: Vec3::new(1., 0., 0.),
            points: vec![(1.4, 0.), (1.5, 0.1), (1.5, 1.5), (1.375, 1.5), (1.375, 0.)],
        },
        TOL,
    )
    .unwrap();
    let turn = gc_kernel::frame::Rigid::from_rows(
        [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
        Vec3::new(0., 0., 0.),
    )
    .unwrap();
    let shallow = gc_kernel::extrude::extrude(
        &Profile {
            loops: vec![vec![Curve::Arc {
                c: Vec3::new(0., 0., 0.),
                n: Vec3::new(0., 0., 1.),
                x: Vec3::new(1., 0., 0.),
                r: 0.38,
                t0: 0.,
                t1: TAU,
            }]],
            n: Vec3::new(0., 0., 1.),
            low: 1.45,
            high: 1.501,
        },
        TOL,
    )
    .unwrap();
    for (body, cutter) in [
        (&done, t.clone()),
        (&mid, t.clone()),
        (&near_bevel, t.clone()),
        (&b, turn.body(&t)),
        (&b, shallow),
        (&b, tool(0., 1.38)),
        (&b, tool(0.001, 0.38)),
    ] {
        let original = (write_body(body), write_body(&cutter));
        assert!(matches!(
            boolean::cut(body, &cutter, TOL),
            Err(gc_kernel::KernelError::Unsupported(_))
        ));
        assert_eq!(original, (write_body(body), write_body(&cutter)));
    }
}

#[test]
fn missing_persistent_history_is_explicit_for_the_new_boundary() {
    let b = stock(false);
    let t = tool(0., 0.38);
    assert!(
        matches!(boolean::cut_named(&b,&t,TOL),Err(gc_kernel::KernelError::Unsupported(e)) if e.contains("RadialNotch"))
    );
    let (out, history) = boolean::cut_recorded(&b, &t, TOL).unwrap();
    assert!(history.unwrap_err().contains("RadialNotch"));
    verify(&out, mass::measure(&b, TOL).unwrap().volume - removed());
}

#[test]
fn inverted_cutter_is_not_treated_as_a_positive_cutting_volume() {
    let b = stock(false);
    let mut t = tool(0., 0.38);
    for f in &mut t.faces {
        f.sense = !f.sense;
        for lp in &mut f.loops {
            lp.reverse();
            for (_, fw) in lp {
                *fw = !*fw;
            }
        }
    }
    topology::check(&t, TOL).unwrap();
    let original = (write_body(&b), write_body(&t));
    assert!(
        matches!(boolean::cut(&b,&t,TOL),Err(gc_kernel::KernelError::Unsupported(e)) if e.contains("inverted"))
    );
    assert_eq!(original, (write_body(&b), write_body(&t)));
}

#[test]
fn scaled_stock_and_tools_preserve_the_dimensionless_cut() {
    fn scaled(body: &Body, k: f64) -> Body {
        let mut out = body.clone();
        for v in &mut out.vertices {
            *v = v.scale(k);
        }
        for e in &mut out.edges {
            match &mut e.curve {
                Curve::Arc { c, r, .. } => {
                    *c = c.scale(k);
                    *r *= k;
                }
                Curve::Line { a, b } => {
                    *a = a.scale(k);
                    *b = b.scale(k);
                }
                _ => panic!("fixture must only contain arcs and lines"),
            }
        }
        for f in &mut out.faces {
            match &mut f.surface {
                gc_kernel::geom::Surface::Plane { p, .. } => *p = p.scale(k),
                gc_kernel::geom::Surface::Cylinder { p, r, .. }
                | gc_kernel::geom::Surface::Cone { p, r, .. } => {
                    *p = p.scale(k);
                    *r *= k;
                }
                _ => panic!("fixture must be analytic axisymmetric stock"),
            }
        }
        out
    }
    for bevel in [false, true] {
        let b = stock(bevel);
        let t = tool(0., 0.38);
        let want = mass::measure(&b, TOL).unwrap().volume - removed();
        for k in [0.25_f64, 4.] {
            let b = scaled(&b, k);
            let t = scaled(&t, k);
            verify(&boolean::cut(&b, &t, TOL).unwrap(), want * k.powi(3));
        }
    }
}

use gc_kernel::{
    boolean,
    geom::{Curve, Surface, Vec3},
    mass,
    neutral::{write_body, Body, Profile, Reader},
    revolve::{revolve, Section},
    tessellate, topology,
};
use std::f64::consts::{PI, TAU};
const TOL: f64 = 1e-6;
fn stock(height: f64, skin: f64) -> Body {
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
fn tool(rho: f64) -> Body {
    gc_kernel::extrude::extrude(
        &Profile {
            loops: vec![vec![Curve::Arc {
                c: Vec3::new(0., 0., 0.),
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
// Independent disk-strip integration, split at the cone/cylinder transition.
fn reference(height: f64, skin: f64, rho: f64) -> f64 {
    let primitive = |r: f64, a: f64| a * (r * r - a * a).sqrt() + r * r * (a / r).asin();
    let f = |t: f64| {
        let y = rho * t.sin();
        let width = rho * t.cos();
        let r = if y < skin {
            1.5 - height + height * y / skin
        } else {
            1.5
        };
        (primitive(r, width) - primitive(1.3, width)) * rho * t.cos()
    };
    let integrate = |lo: f64, hi: f64| {
        let n = 2048;
        let step = (hi - lo) / n as f64;
        let mut sum = f(lo) + f(hi);
        for i in 1..n {
            sum += (if i % 2 == 0 { 2. } else { 4. }) * f(lo + i as f64 * step);
        }
        sum * step / 3.
    };
    let beta = (skin / rho).asin();
    PI * (1.5_f64.powi(2) - 1.3_f64.powi(2)) * 3.
        - PI * height * skin * (1.5 - height / 3.)
        - integrate(0., beta)
        - integrate(beta, PI / 2.)
}
#[test]
fn notch_crosses_the_end_cone_and_cylinder_without_regenerating_stock() {
    for (height, skin, rho) in [
        (0.1, 0.1 * 30_f64.to_radians().tan(), 0.5),
        (0.1, 0.1, 0.38),
        (0.05, 0.2, 0.4),
    ] {
        let body = stock(height, skin);
        let cutter = tool(rho);
        let original = (write_body(&body), write_body(&cutter));
        let out = boolean::cut(&body, &cutter, TOL).unwrap();
        let wanted = reference(height, skin, rho);
        verify(&out, wanted);
        assert_eq!(
            out.edges
                .iter()
                .filter(|e| matches!(e.curve, Curve::ConePair { .. }))
                .count(),
            2
        );
        assert_eq!(
            out.edges
                .iter()
                .filter(|e| matches!(e.curve, Curve::CylinderPair { .. }))
                .count(),
            2
        );
        for (fi, face) in body.faces.iter().enumerate() {
            if matches!(face.surface,Surface::Plane{p,..} if p.y==3.) {
                assert_eq!(out.faces[fi].loops, face.loops);
                assert_eq!(
                    format!("{:?}", out.faces[fi].surface),
                    format!("{:?}", face.surface)
                );
                for &(ei, _) in face.loops.iter().flatten() {
                    assert_eq!(
                        format!("{:?}", out.edges[ei]),
                        format!("{:?}", body.edges[ei])
                    );
                }
            }
        }
        assert_eq!(original, (write_body(&body), write_body(&cutter)));
        if rho == 0.5 {
            assert!((wanted - 5.17435183987158620).abs() < 1e-12);
        }
    }
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
fn near_bevel_notch_survives_cardinal_frames_and_the_opposite_end() {
    let skin = 0.1 * 30_f64.to_radians().tan();
    let body = stock(0.1, skin);
    let cutter = tool(0.5);
    let wanted = reference(0.1, skin, 0.5);
    for rows in [
        [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
        [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
    ] {
        let rigid = gc_kernel::frame::Rigid::from_rows(rows, Vec3::new(0., 0., 0.)).unwrap();
        let b = rigid.body(&body);
        let t = rigid.body(&cutter);
        verify(&boolean::cut(&b, &t, TOL).unwrap(), wanted);
    }
    // Integer near-transition coordinates remain exactly representable after translation.
    let b = stock(0.1, 0.25);
    let rigid = gc_kernel::frame::Rigid::from_rows(
        [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
        Vec3::new(4., -8., 3.),
    )
    .unwrap();
    verify(
        &boolean::cut(&rigid.body(&b), &rigid.body(&cutter), TOL).unwrap(),
        reference(0.1, 0.25, 0.5),
    );
}
#[test]
fn near_bevel_refusals_preserve_both_operands() {
    let b = stock(0.1, 0.1);
    let t = tool(0.5);
    let done = boolean::cut(&b, &t, TOL).unwrap();
    let seam_turn = gc_kernel::frame::Rigid::from_rows(
        [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
        Vec3::new(0., 0., 0.),
    )
    .unwrap();
    let offset = gc_kernel::frame::Rigid::from_rows(
        [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        Vec3::new(0., 0.001, 0.),
    )
    .unwrap();
    for (body, cutter) in [
        (done, tool(0.501)), // Exact repeats have a separate whole-boundary proof.
        (stock(0.1, 0.5), t.clone()),
        (stock(0.1, 0.6), t.clone()),
        (b.clone(), tool(1.31)),
        (b.clone(), seam_turn.body(&t)),
        (b.clone(), offset.body(&t)),
    ] {
        let original = (write_body(&body), write_body(&cutter));
        assert!(matches!(
            boolean::cut(&body, &cutter, TOL),
            Err(gc_kernel::KernelError::Unsupported(_))
        ));
        assert_eq!(original, (write_body(&body), write_body(&cutter)));
    }
}
#[test]
fn near_bevel_history_is_explicit_and_never_fabricated() {
    let b = stock(0.1, 0.1);
    let t = tool(0.5);
    assert!(matches!(boolean::cut_named(&b,&t,TOL),
        Err(gc_kernel::KernelError::Unsupported(e)) if e.contains("RadialNotch")));
    let (out, history) = boolean::cut_recorded(&b, &t, TOL).unwrap();
    assert!(history.unwrap_err().contains("RadialNotch"));
    verify(&out, reference(0.1, 0.1, 0.5));
}
#[test]
fn scaled_near_bevel_notches_keep_the_same_dimensionless_boundary() {
    let b = stock(0.1, 0.1);
    let t = tool(0.5);
    let wanted = reference(0.1, 0.1, 0.5);
    for scale in [0.25, 4.] {
        let out = boolean::cut(&scaled(&b, scale), &scaled(&t, scale), TOL * scale).unwrap();
        verify(&out, wanted * scale.powi(3));
    }
}

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

//! The CylinderPair curve (Claude, layout agreed with Astra): a cutter cylinder crossing a host cylinder.

use std::f64::consts::PI;

use gc_kernel::geom::{Curve, Surface, Vec3};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{put_curve, write_body, Body, Edge, Face, Reader};
use gc_kernel::topology::check;
use gc_kernel::KernelError;

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);
const O: Vec3 = Vec3::new(0.0, 0.0, 0.0);
const RHO: f64 = 0.4;
const R: f64 = 0.5;

/// Cutter radius .4 along z through the origin, host radius .5 along y through the origin, upper branch:
/// z = sqrt(R^2 - x^2) over the cutter circle (the radii of GCmaker's circle-over-fillet fixture).
fn plug_edge() -> Curve {
    Curve::CylinderPair { c: O, a: Z, x: X, rho: RHO, p: O, b: Y, r: R, branch: 1.0, t0: 0.0, t1: 2.0 * PI }
}

fn on_both(q: Vec3) -> (f64, f64) { ((q.x * q.x + q.y * q.y).sqrt() - RHO, (q.x * q.x + q.z * q.z).sqrt() - R) }

#[test]
fn points_lie_on_both_cylinders_and_tangents_are_derivatives() {
    let c = plug_edge();
    for i in 0..=64 {
        let t = 2.0 * PI * i as f64 / 64.0;
        let (q, d) = c.eval(t);
        let (e1, e2) = on_both(q);
        assert!(e1.abs() < 1e-15 && e2.abs() < 1e-15, "{e1} {e2}");
        assert!(q.z > 0.0);
        let h = 1e-6;
        let numeric = c.eval(t + h).0.sub(c.eval(t - h).0).scale(0.5 / h);
        assert!(numeric.distance(d) < 1e-8, "tangent at {t}");
    }
}

/// Reversal flips the cutter axis and the branch (Astra): the same points, run the other way, at interior points.
#[test]
fn reversal_is_the_same_curve_backwards_at_interior_points() {
    let oblique = Curve::CylinderPair {
        c: Vec3::new(0.1, -0.2, 0.3), a: Vec3::new(0.0, 0.6, 0.8), x: X, rho: 0.3,
        p: Vec3::new(0.0, 0.0, 0.0), b: Vec3::new(0.0, 0.0, 1.0).sub(Vec3::new(0.0, 0.0, 0.0)), r: 2.0,
        branch: -1.0, t0: 0.3, t1: 2.1,
    };
    for c in [plug_edge(), oblique] {
        let r = c.reversed();
        let (t0, t1) = c.domain();
        assert_eq!(r.domain(), (-t1, -t0));
        for i in 0..=20 {
            let t = t0 + (t1 - t0) * i as f64 / 20.0;
            let (q, d) = c.eval(t);
            let (qr, dr) = r.eval(-t);
            assert!(q.distance(qr) < 1e-14, "point at {t}");
            assert!(d.add(dr).length() < 1e-12, "tangent at {t}");
        }
        assert!(r.start().distance(c.end()) < 1e-14 && r.end().distance(c.start()) < 1e-14);
        let d = Vec3::new(3.0, -1.0, 2.0);
        assert!(c.translated(d).eval(t0 + 0.4).0.distance(c.eval(t0 + 0.4).0.add(d)) < 1e-14);
    }
}

fn pair_values(c: &Curve) -> Vec<f64> { let mut v = Vec::new(); put_curve(&mut v, c); v }

fn read(values: &[f64]) -> Result<Curve, KernelError> {
    let mut r = Reader::new(values);
    let c = r.curve()?;
    r.finish()?;
    Ok(c)
}

#[test]
fn reads_strictly_and_round_trips() {
    let v = pair_values(&plug_edge());
    assert_eq!(v.len(), 21);
    assert_eq!(pair_values(&read(&v).unwrap()), v);
    let refuse = |change: &dyn Fn(&mut Vec<f64>), why: &str| {
        let mut w = v.clone();
        change(&mut w);
        assert!(matches!(read(&w), Err(KernelError::InvalidInput(_))), "{why}");
    };
    refuse(&|w| w[18] = 0.5, "branch not +-1");
    refuse(&|w| w[10] = 0.0, "zero cutter radius");
    refuse(&|w| { w[14] = 0.0; w[15] = 0.0; w[16] = 1.0; }, "parallel axes");
    refuse(&|w| w[17] = 0.39, "host thinner than the cutter: the discriminant goes negative");
    refuse(&|w| { w[19] = 1.0; w[20] = 1.0; }, "empty interval");
    refuse(&|w| w[13] = f64::NAN, "non-finite");
    refuse(&|w| { w.pop(); }, "truncated");
    // Just thicker than the cutter: valid, the discriminant is small but proved positive.
    let mut thin = v.clone();
    thin[17] = 0.4001;
    assert!(read(&thin).is_ok());
}

/// The plug: inside the cutter, above z = 0, below the host. Bottom disc, cutter wall, a patch of the host.
fn plug() -> Body {
    let top = plug_edge();
    let v0 = Vec3::new(RHO, 0.0, 0.0);
    let v1 = top.start();
    let bottom = Curve::Arc { c: O, n: Z, x: X, r: RHO, t0: 0.0, t1: 2.0 * PI };
    Body {
        vertices: vec![v0, v1],
        edges: vec![
            Edge { curve: bottom, start: 0, end: 0 },
            Edge { curve: top, start: 1, end: 1 },
            Edge { curve: Curve::Line { a: v0, b: v1 }, start: 0, end: 1 },
        ],
        faces: vec![
            Face { surface: Surface::Plane { p: O, n: Z, x: X }, sense: false, loops: vec![vec![(0, false)]] },
            Face { surface: Surface::Cylinder { p: O, a: Z, x: X, r: RHO }, sense: true,
                   loops: vec![vec![(0, true), (2, true), (1, false), (2, false)]] },
            Face { surface: Surface::Cylinder { p: O, a: Y, x: Z, r: R }, sense: true, loops: vec![vec![(1, true)]] },
        ],
    }
}

/// Simpson on the substitution x = rho sin(phi), so the square-root end behaviour disappears.
fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;
    (0..=n).map(|i| f(a + h * i as f64) * if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 }).sum::<f64>() * h / 3.0
}

#[test]
fn plug_checks_and_measures_by_quadrature() {
    let body = plug();
    check(&body, 1e-9).unwrap();
    let data = write_body(&body);
    let mut r = Reader::new(&data);
    assert_eq!(write_body(&r.body().unwrap()), data);
    r.finish().unwrap();

    let m = measure(&body, 1e-9).unwrap();
    assert!(!m.exact);
    // Volume: integral over the disc of sqrt(R^2 - x^2) = int 2 rho^2 cos^2(phi) sqrt(R^2 - rho^2 sin^2 phi) dphi.
    let volume = simpson(|f| 2.0 * RHO * RHO * f.cos().powi(2) * (R * R - RHO * RHO * f.sin().powi(2)).sqrt(), -PI / 2.0, PI / 2.0, 20000);
    // Area: disc + wall (rho times the height round the cutter) + host patch (R / sqrt(R^2 - x^2) over the disc).
    let wall = RHO * simpson(|t| (R * R - RHO * RHO * t.cos().powi(2)).sqrt(), 0.0, 2.0 * PI, 20000);
    let patch = simpson(|f| 2.0 * RHO * RHO * f.cos().powi(2) * R / (R * R - RHO * RHO * f.sin().powi(2)).sqrt(), -PI / 2.0, PI / 2.0, 20000);
    assert!((m.volume - volume).abs() < 1e-12, "{} vs {volume}", m.volume);
    assert!((m.area - (PI * RHO * RHO + wall + patch)).abs() < 1e-12, "{} vs {}", m.area, PI * RHO * RHO + wall + patch);
}

/// Astra's case: the same plug moved far away, its cylinders' axis points anchored far along their axes. The measure
/// matches the plug at home and finishes promptly (axial distances from the body's centroid, bounded work).
#[test]
fn translated_plug_with_far_axis_anchors_measures_the_same() {
    let home = measure(&plug(), 1e-9).unwrap();
    let d = Vec3::new(1.0e4, -3.0e4, 2.0e4);
    let mut far = plug();
    for v in &mut far.vertices { *v = v.add(d); }
    for e in &mut far.edges { e.curve = e.curve.translated(d); }
    // The host's axis point may sit anywhere on its axis (exact: only its part across b enters). The cutter's c
    // may not: h then cancels at that size, so writers keep c near the curve.
    if let Curve::CylinderPair { p, b, .. } = &mut far.edges[1].curve {
        *p = p.add(b.scale(7.0e5));
    }
    for f in &mut far.faces {
        f.surface = match f.surface.clone() {
            Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } => unreachable!("no elliptic-cylinder or conical faces in these test bodies"),
            Surface::Plane { p, n, x } => Surface::Plane { p: p.add(d), n, x },
            Surface::Cylinder { p, a, x, r } => Surface::Cylinder { p: p.add(d).add(a.scale(-9.0e5)), a, x, r },
        };
    }
    check(&far, 1e-7).unwrap();
    let clock = std::time::Instant::now();
    let m = measure(&far, 1e-7).unwrap();
    assert!(clock.elapsed().as_secs_f64() < 2.0, "{:?}", clock.elapsed());
    assert!((m.volume - home.volume).abs() < 1e-9 * home.volume, "{} vs {}", m.volume, home.volume);
    assert!((m.area - home.area).abs() < 1e-9 * home.area, "{} vs {}", m.area, home.area);
}

#[test]
fn quadrature_refuses_what_it_cannot_converge() {
    use gc_kernel::mass::quadrature;
    assert!((quadrature(&|t| t.cos(), 0.0, 1.0, 1e-15).unwrap() - 1f64.sin()).abs() < 1e-15);
    assert!(quadrature(&|t| (1.0 / t).sin(), 1e-9, 1.0, 1e-15).is_err());
}

/// Astra's probe: only the parameter cylinder's anchor moved to z = 1e8 along its own axis. The discriminant in
/// Lagrange form (alpha r^2 - |D' x a'|^2) keeps the curve: the points stay on both cylinders and on the plug's top.
#[test]
fn far_parameter_anchor_keeps_the_curve() {
    let near = plug_edge();
    let Curve::CylinderPair { a, x, rho, p, b, r, branch, t0, t1, .. } = near.clone() else { unreachable!() };
    let far = Curve::CylinderPair { c: Vec3::new(0.0, 0.0, 1.0e8), a, x, rho, p, b, r, branch, t0, t1 };
    for i in 0..=16 {
        let t = 2.0 * PI * i as f64 / 16.0;
        let (q, qn) = (far.eval(t).0, near.eval(t).0);
        assert!(q.distance(qn) < 1e-6, "{q:?} vs {qn:?}");
    }
}

/// The analytic membership bounds refuse curves just off their face: an arc 1e-5 too big for its cylinder, an arc
/// tilted out of its plane, a CylinderPair on a plane face.
#[test]
fn curves_just_off_their_surfaces_are_refused() {
    let mut body = plug();
    let Curve::Arc { r, .. } = &mut body.edges[0].curve else { unreachable!() };
    *r += 1e-5;
    body.vertices[0].x += 1e-5;
    body.edges[2].curve = Curve::Line { a: body.vertices[0], b: body.vertices[1] };
    assert!(matches!(check(&body, 1e-7), Err(KernelError::InvalidInput(m)) if m.contains("leaves face")));

    let mut tilted = plug();
    let Curve::Arc { n, .. } = &mut tilted.edges[0].curve else { unreachable!() };
    *n = Vec3::new(0.0, 1e-4, 1.0).unit().unwrap();
    assert!(matches!(check(&tilted, 1e-7), Err(KernelError::InvalidInput(m)) if m.contains("leaves face")));

    let mut on_plane = plug();
    on_plane.faces[2].surface = Surface::Plane { p: Vec3::new(0.0, 0.0, 0.3), n: Z, x: X };
    assert!(matches!(check(&on_plane, 1e-7), Err(KernelError::InvalidInput(m)) if m.contains("leaves face 2")));
}

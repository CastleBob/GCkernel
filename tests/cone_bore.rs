//! Holes through conical walls (Claude 2026-10-31, tapered Pole stock, hollow and solid). A round cutter square to the
//! pole's axis meets each cone wall in a ConePair loop (wire type 8, not a conic). Every removed volume is checked
//! against an independent quadrature from the dimensions: over the hole's disk (x across, y along the pole) the depth
//! of material the hole passes through, sqrt(Ro(y)^2 - x^2) - sqrt(Ri(y)^2 - x^2) per wall (Ro, Ri the walls' radii at
//! that height), by Gauss-Legendre in the disk's radius and the trapezoid rule round it (smooth, periodic). Every
//! result passes the topology check, survives the neutral round trip and tessellates.

use gc_kernel::boolean::cut;
use gc_kernel::geom::{Curve, Vec3};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{Body, Profile};
use gc_kernel::revolve::{revolve, Section};
use gc_kernel::topology::check;
use gc_kernel::KernelError;
use std::f64::consts::TAU;

const TOL: f64 = 1e-6;

/// A pole along +Y: a 5 radius base to height 30, tapering to radius 2 at 60, straight to 90; hollow with `wall`
/// (the inner surface parallel) or solid.
fn pole(wall: f64) -> Body {
    let (rb, rt) = (5.0, 2.0);
    let mut pts = vec![(rb, 0.0), (rb, 30.0), (rt, 60.0), (rt, 90.0)];
    if wall > 0.0 {
        pts.extend([(rt - wall, 90.0), (rt - wall, 60.0), (rb - wall, 30.0), (rb - wall, 0.0)]);
    } else {
        pts.extend([(0.0, 90.0), (0.0, 0.0)]);
    }
    revolve(&Section { p: Vec3::new(0., 0., 0.), a: Vec3::new(0., 1., 0.), x: Vec3::new(1., 0., 0.), points: pts }, TOL).unwrap()
}

/// The outer radius at height y on the taper (30..60).
fn outer(y: f64) -> f64 { 5.0 - (y - 30.0) * 0.1 }

/// A round cutter of radius rho square to the pole: centred at (x0, y0) across it, turned `deg` about the pole from +Z,
/// from `top` down to `bottom` along its axis (heights above the pole's axis, negative past it).
fn bore(rho: f64, x0: f64, y0: f64, deg: f64, top: f64, bottom: f64) -> Body {
    let circle = vec![Curve::Arc { c: Vec3::new(x0, y0, 0.), n: Vec3::new(0., 0., 1.), x: Vec3::new(1., 0., 0.), r: rho, t0: 0., t1: TAU }];
    let tool = gc_kernel::extrude::extrude(&Profile { loops: vec![circle], n: Vec3::new(0., 0., 1.), low: bottom, high: top }, TOL).unwrap();
    let t = deg.to_radians();
    let turn = gc_kernel::frame::Rigid::from_rows([[t.cos(), 0., t.sin()], [0., 1., 0.], [-t.sin(), 0., t.cos()]], Vec3::new(0., 0., 0.)).unwrap();
    turn.body(&tool)
}

/// The integral over the hole's disk (centre x0, y0, radius rho) of f(x, y): Gauss-Legendre (16 points on each of 16
/// panels) in the radius, the trapezoid rule (512 points) round the disk.
fn over_disk(x0: f64, y0: f64, rho: f64, f: &dyn Fn(f64, f64) -> f64) -> f64 {
    let nodes: Vec<(f64, f64)> = {
        // Legendre roots and weights by Newton's method, independently of the kernel's own table.
        let n = 16;
        (0..n).map(|i| {
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let mut dp = 0.;
            for _ in 0..100 {
                let (mut p0, mut p1) = (1., x);
                for k in 2..=n { let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64; p0 = p1; p1 = p2; }
                dp = n as f64 * (x * p1 - p0) / (x * x - 1.);
                let dx = p1 / dp;
                x -= dx;
                if dx.abs() < 1e-16 { break; }
            }
            (x, 2. / ((1. - x * x) * dp * dp))
        }).collect()
    };
    let (panels, m) = (16, 512);
    let mut total = 0.;
    for p in 0..panels {
        let (a, b) = (rho * p as f64 / panels as f64, rho * (p + 1) as f64 / panels as f64);
        for &(u, w) in &nodes {
            let r = (a + b) / 2. + (b - a) / 2. * u;
            let ring: f64 = (0..m).map(|j| { let t = TAU * j as f64 / m as f64; f(x0 + r * t.cos(), y0 + r * t.sin()) }).sum::<f64>() * TAU / m as f64;
            total += w * (b - a) / 2. * r * ring;
        }
    }
    total
}

fn removed(body: &Body, tool: &Body) -> (f64, Body) {
    let result = cut(body, tool, TOL).unwrap_or_else(|e| panic!("{e}"));
    check(&result, TOL).unwrap();
    // The neutral round trip (the reader takes directions to unit length, so a turned cutter's axis may move an ulp).
    let wire = gc_kernel::neutral::write_body(&result);
    let mut reader = gc_kernel::neutral::Reader::with_tolerance(&wire, TOL);
    let back = reader.body().unwrap();
    reader.finish().unwrap();
    check(&back, TOL).unwrap();
    assert!((measure(&back, TOL).unwrap().volume - measure(&result, TOL).unwrap().volume).abs() < 1e-12);
    gc_kernel::tessellate::tessellate(&result, 0.01).unwrap();
    (measure(body, TOL).unwrap().volume - measure(&result, TOL).unwrap().volume, result)
}

fn wall(r: f64, x: f64) -> f64 { (r * r - x * x).sqrt() }

#[test]
fn a_hole_through_one_tapered_wall_removes_its_disk_integral() {
    let body = pole(0.25);
    // A 1 inch hole at height 45 (outer radius 3.5), from 6 above the axis down to 1 above it (in the bore).
    let (got, result) = removed(&body, &bore(0.5, 0., 45., 0., 6., 1.));
    let exact = over_disk(0., 45., 0.5, &|x, y| wall(outer(y), x) - wall(outer(y) - 0.25, x));
    assert!((got - exact).abs() < 1e-9, "removed {got}, disk integral {exact}");
    assert_eq!(result.edges.iter().filter(|e| matches!(e.curve, Curve::ConePair { .. })).count(), 2);
}

#[test]
fn through_both_walls_off_centre_and_turned() {
    let body = pole(0.25);
    // A 0.75 hole offset 0.6 across the pole at height 50, turned 130 degrees about it, right through.
    let (got, _) = removed(&body, &bore(0.375, 0.6, 50., 130., 6., -6.));
    let exact = 2. * over_disk(0.6, 50., 0.375, &|x, y| wall(outer(y), x) - wall(outer(y) - 0.25, x));
    assert!((got - exact).abs() < 1e-9, "removed {got}, disk integral {exact}");
}

#[test]
fn solid_tapered_pole_through_and_blind() {
    let body = pole(0.);
    let (got, _) = removed(&body, &bore(0.5, 0., 40., 0., 6., -6.));
    let exact = over_disk(0., 40., 0.5, &|x, y| 2. * wall(outer(y), x));
    assert!((got - exact).abs() < 1e-9, "through: removed {got}, disk integral {exact}");
    // Blind: from the surface down to 1.5 above the axis (outer radius at 40 is 4).
    let (got, result) = removed(&body, &bore(0.5, 0., 40., 0., 6., 1.5));
    let exact = over_disk(0., 40., 0.5, &|x, y| wall(outer(y), x) - 1.5);
    assert!((got - exact).abs() < 1e-9, "blind: removed {got}, disk integral {exact}");
    assert!(result.faces.iter().any(|f| matches!(f.surface, gc_kernel::geom::Surface::Plane { n, .. } if (n.z - 1.).abs() < 1e-12)));
}

#[test]
fn several_holes_in_turn_and_one_on_the_straight_section() {
    let mut body = pole(0.25);
    let start = measure(&body, TOL).unwrap().volume;
    let mut exact = 0.;
    for (y, deg) in [(35., 0.), (42., 200.), (55., 120.)] {
        body = cut(&body, &bore(0.25, 0., y, deg, 6., 0.5), TOL).unwrap_or_else(|e| panic!("hole at {y}: {e}"));
        check(&body, TOL).unwrap();
        exact += over_disk(0., y, 0.25, &|x, yy| wall(outer(yy), x) - wall(outer(yy) - 0.25, x));
    }
    let got = start - measure(&body, TOL).unwrap().volume;
    assert!((got - exact).abs() < 1e-9, "three holes removed {got}, integrals {exact}");
    // The straight base (radius 5, a cylinder) of the same pole takes the same domain.
    let (got, _) = removed(&body, &bore(0.25, 0., 10., 0., 6., 3.));
    let exact = over_disk(0., 10., 0.25, &|x, _| wall(5., x) - wall(4.75, x));
    assert!((got - exact).abs() < 1e-9, "base hole removed {got}, integral {exact}");
}

#[test]
fn cone_wall_holes_outside_the_domain_are_refused_by_name() {
    let body = pole(0.25);
    let refused = |tool: Body| match cut(&body, &tool, TOL) {
        Err(KernelError::Unsupported(m)) => m,
        other => panic!("expected a named refusal, got {:?}", other.map(|b| b.faces.len())),
    };
    // Across the ring where the taper starts (height 30).
    assert!(refused(bore(0.5, 0., 30.2, 0., 6., 1.)).contains("axial window"));
    // Ending inside the wall.
    assert!(refused(bore(0.5, 0., 45., 0., 6., 3.4)).contains("inside the wall"));
    // Not square to the pole: a cutter along the pole's own axis.
    let along = gc_kernel::extrude::extrude(&Profile {
        loops: vec![vec![Curve::Arc { c: Vec3::new(0., 0., 4.), n: Vec3::new(0., 1., 0.), x: Vec3::new(1., 0., 0.), r: 0.2, t0: 0., t1: TAU }]],
        n: Vec3::new(0., 1., 0.), low: 40., high: 50.,
    }, TOL).unwrap();
    assert!(refused(along).contains("square"));
    // On the wall's seam (the revolve's seam lies at +X, 90 degrees from +Z): refused, the seam is not moved on cones.
    assert!(refused(bore(0.5, 0., 45., 90., 6., 1.)).contains("seam"));
    // Two holes overlapping: the second meets the first.
    let once = cut(&body, &bore(0.5, 0., 45., 0., 6., 1.), TOL).unwrap();
    match cut(&once, &bore(0.5, 0.3, 45.4, 0., 6., 1.), TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("earlier hole"), "{m}"),
        other => panic!("expected a named refusal, got {:?}", other.map(|b| b.faces.len())),
    }
}

//! Bob's linear-nesting parts that the GCmaker native survey of 2026-10-31 found stopping in the kernel (Claude). Each
//! fixture is the worker's own input from the GC Kernel view: the body before the refused cut and the extruded tool.
//! Every removed volume is checked against an independent value from the parts' saved dimensions (closed form where
//! there is one, otherwise a quadrature of the analytic section written here, not the kernel's B-rep or mass code),
//! and every result must pass the topology check.

use gc_kernel::geom::{Curve, Vec3};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{Body, Profile, Reader};
use gc_kernel::topology::check;
use gc_kernel::KernelError;
use std::f64::consts::{PI, TAU};

const TOL: f64 = 1e-6;

fn read(case: &str, name: &str) -> Body {
    let path = format!("{}/tests/fixtures/survey_gaps/{case}/{name}", env!("CARGO_MANIFEST_DIR"));
    let bytes = std::fs::read(path).unwrap();
    let data: Vec<f64> = bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
    Reader::with_tolerance(&data, TOL).body().unwrap()
}

/// Cuts the case and returns (removed volume, result).
fn removed(case: &str) -> (f64, Body) {
    let (body, tool) = (read(case, "body.f64"), read(case, "tool.f64"));
    let result = gc_kernel::boolean::cut(&body, &tool, TOL).unwrap_or_else(|e| panic!("{case}: {e}"));
    let report = check(&result, TOL).unwrap_or_else(|e| panic!("{case}: result topology {e:?}"));
    assert!(report.faces > 0);
    (measure(&body, TOL).unwrap().volume - measure(&result, TOL).unwrap().volume, result)
}

/// Adaptive Simpson, to `eps` absolute.
fn simpson(f: &dyn Fn(f64) -> f64, a: f64, b: f64, eps: f64) -> f64 {
    fn step(f: &dyn Fn(f64) -> f64, a: f64, b: f64, fa: f64, fm: f64, fb: f64, whole: f64, eps: f64, depth: u32) -> f64 {
        let m = (a + b) / 2.;
        let (lm, rm) = ((a + m) / 2., (m + b) / 2.);
        let (flm, frm) = (f(lm), f(rm));
        let left = (m - a) / 6. * (fa + 4. * flm + fm);
        let right = (b - m) / 6. * (fm + 4. * frm + fb);
        if depth == 0 || (left + right - whole).abs() <= 15. * eps {
            return left + right + (left + right - whole) / 15.;
        }
        step(f, a, m, fa, flm, fm, left, eps / 2., depth - 1) + step(f, m, b, fm, frm, fb, right, eps / 2., depth - 1)
    }
    let (fa, fm, fb) = (f(a), f((a + b) / 2.), f(b));
    step(f, a, b, fa, fm, fb, (b - a) / 6. * (fa + 4. * fm + fb), eps, 48)
}

/// Integral over the pipe annulus (ri..ro, all the way round) of g(x, z) dA: Simpson in r, the trapezoid rule in
/// angle (exact to roundoff for a smooth periodic integrand with this many points).
fn annulus(ri: f64, ro: f64, g: &dyn Fn(f64, f64) -> f64) -> f64 {
    let n = 1024;
    simpson(
        &|r| (0..n).map(|j| g(r * (TAU * j as f64 / n as f64).cos(), r * (TAU * j as f64 / n as f64).sin())).sum::<f64>() * TAU / n as f64 * r,
        ri,
        ro,
        1e-13,
    )
}

// L3X3X1/4 (heel at x 0, z 0; horizontal leg z in [-0.25, 0] to its toe at x 3, toe radius 0.125 about
// (2.875, -0.125) on the underside). Bob's cope box: x in [1.25, 2.95], z from below the leg to above it, y in
// [-0.45, 1.75] from the end. It takes the leg's whole thickness over x in [1.25, 2.875] and, over [2.875, 2.95], the
// toe section 0.125 + sqrt(r^2 - u^2), u = x - 2.875; the stock starts at y = 0.
fn angle_cope() -> f64 {
    let r: f64 = 0.125;
    let u: f64 = 0.075;
    let toe = 0.125 * u + (u * (r * r - u * u).sqrt() + r * r * (u / r).asin()) / 2.;
    1.75 * (0.25 * (2.875 - 1.25) + toe)
}

#[test]
fn angle_cope_box_spanning_no_stock_axis_removes_its_section() {
    // A: the box spans none of the stock's axes; it runs in the transverse chart (previously refused).
    let (got, _) = removed("angle_cope_one_end");
    assert!((got - angle_cope()).abs() < 1e-12, "removed {got}, closed form {}", angle_cope());
}

#[test]
fn angle_cope_after_a_miter_chooses_the_chart_the_miter_fits() {
    // A: the earlier miter's plane turns in X and Y, so the chart sweeping X cannot hold it; Z is chosen.
    let (got, _) = removed("angle_cope_after_miter");
    assert!((got - angle_cope()).abs() < 1e-12, "removed {got}, closed form {}", angle_cope());
}

#[test]
fn angle_second_cope_joins_the_first_copes_seams() {
    // A: the first cope's result keeps coplanar cell seams (the heel face split at y 1.75); their section segments
    // are one boundary.
    let (got, _) = removed("angle_second_cope");
    assert!((got - angle_cope()).abs() < 1e-12, "removed {got}, closed form {}", angle_cope());
}

/// One wall of a pipe (ro, ri) pierced square by a round hole of radius a: the integral of
/// 2 sqrt(a^2 - x^2) (sqrt(ro^2 - x^2) - sqrt(ri^2 - x^2)) over (-a, a), with x = a sin t.
fn one_wall(a: f64, ro: f64, ri: f64) -> f64 {
    simpson(
        &|t: f64| {
            let s = a * t.sin();
            2. * a * a * t.cos() * t.cos() * ((ro * ro - s * s).sqrt() - (ri * ri - s * s).sqrt())
        },
        -PI / 2.,
        PI / 2.,
        1e-15,
    )
}

#[test]
fn side_hole_on_the_pipe_seam_moves_the_seam_and_cuts_one_wall() {
    // B: 1.25 pipe (OD 1.66, wall 0.14), half-inch hole at 3:00 where the wall's seam is; one wall deep.
    let (got, result) = removed("pipe_side_hole_on_seam");
    let exact = one_wall(0.25, 0.83, 0.69);
    assert!((got - exact).abs() < 1e-10, "removed {got}, one wall {exact}");
    // Two holes now: genus 3 (the bore and two wall holes).
    assert_eq!(check(&result, TOL).unwrap().genus(), 3.);
}

#[test]
fn square_cope_after_a_cope_turned_across_another_axis_runs_on_a_split_piece() {
    // B: cope OD 1.664 square across the pipe at y = 0 (axis X), after the far end's cope across Z at y = 18. The
    // removed region is y in [0, sqrt(Rc^2 - z^2)] over the pipe's annulus.
    let (got, _) = removed("pipe_cope_after_turned_cope");
    let rc: f64 = 0.832;
    let exact = annulus(0.69, 0.83, &|_x, z| (rc * rc - z * z).sqrt());
    assert!((got - exact).abs() < 1e-10, "removed {got}, section integral {exact}");
}

#[test]
fn pipe_cross_slot_after_a_miter_runs_on_a_split_piece() {
    // D: a 0.5 wide, 2.0 long slot across the top of the pipe, from z 0.33 up, after the miter at the other end. Its
    // removed volume is the integral over x of the slot's width times the top wall's height above z 0.33.
    let (got, _) = removed("pipe_slot_after_miter");
    let (ro, ri, z0): (f64, f64, f64) = (0.83, 0.69, 0.33);
    let column = |x: f64| {
        let top = (ro * ro - x * x).max(0.).sqrt();
        let bottom = if x * x < ri * ri { (ri * ri - x * x).sqrt().max(z0) } else { z0 };
        (top - bottom).max(0.)
    };
    let width = |x: f64| if x <= 0.75 { 0.5 } else { 2. * (0.0625 - (x - 0.75) * (x - 0.75)).max(0.).sqrt() };
    let f = |x: f64| column(x) * width(x);
    let knees = [0., (ri * ri - z0 * z0).sqrt(), ri, 0.75, (ro * ro - z0 * z0).sqrt()];
    let exact = 2. * knees.windows(2).map(|w| simpson(&f, w[0], w[1], 1e-14)).sum::<f64>();
    assert!((got - exact).abs() < 1e-9, "removed {got}, section integral {exact}");
}

/// The larger root y of the line (x, ., z) on the cylinder radius rc about the axis through `q` along `d`.
fn y_on(x: f64, z: f64, q: Vec3, d: Vec3, rc: f64) -> f64 {
    // |w|^2 - (w.d)^2 = rc^2 with w = (x, y, z) - q: quadratic in y.
    let w0 = Vec3::new(x - q.x, -q.y, z - q.z);
    let (a, b, c) = (1. - d.y * d.y, 2. * (w0.y - w0.dot(d) * d.y), w0.dot(w0) - w0.dot(d).powi(2) - rc * rc);
    (-b + (b * b - 4. * a * c).sqrt()) / (2. * a)
}

#[test]
fn pipe_coped_at_sixty_degrees_to_a_larger_pipe() {
    // C: Bob's 1.25 pipe coped to a 2.375 OD run meeting it at 60 degrees: the cutter's axis (0, -sin 30, cos 30)
    // through the pipe end's centre. Removed: the integral over the annulus of the cope's reach y+ (the end at y 0).
    let (got, result) = removed("pipe_cope_60");
    let d = Vec3::new(0., -0.5, 0.75_f64.sqrt());
    let exact = annulus(0.69, 0.83, &|x, z| y_on(x, z, Vec3::new(0., 0., 0.), d, 1.1875));
    assert!((got - exact).abs() < 1e-10, "removed {got}, annulus integral {exact}");
    assert_eq!(result.faces.len(), 4);
    assert!(result.edges.iter().filter(|e| matches!(e.curve, Curve::CylinderPair { .. })).count() == 2);
}

/// A plain pipe along +Y from y 0 to `length` (outer ro, inner ri), as the shape factory places it.
fn pipe(ro: f64, ri: f64, length: f64) -> Body {
    let ring = |r: f64| vec![Curve::Arc { c: Vec3::new(0., 0., 0.), n: Vec3::new(0., -1., 0.), x: Vec3::new(1., 0., 0.), r, t0: 0., t1: TAU }];
    gc_kernel::extrude::extrude(&Profile { loops: vec![ring(ro), ring(ri)], n: Vec3::new(0., 1., 0.), low: 0., high: length }, TOL).unwrap()
}

/// A round cutter of radius rc along unit d, its axis through q, from s0 to s1 along d.
fn cutter(q: Vec3, d: Vec3, rc: f64, s0: f64, s1: f64) -> Body {
    let x = d.perpendicular();
    let start = q.add(d.scale(s0));
    gc_kernel::extrude::extrude(
        &Profile { loops: vec![vec![Curve::Arc { c: start, n: d, x, r: rc, t0: 0., t1: TAU }]], n: d, low: 0., high: s1 - s0 },
        TOL,
    )
    .unwrap()
}

#[test]
fn oblique_copes_at_any_turn_and_offset_remove_their_annulus_integral() {
    // C, general: a 3.5 x 0.216 pipe coped by a 6 inch cutter at 40 degrees, turned 25 degrees about the pipe and
    // offset 0.3 from its axis, at the y = 0 end; and the same at the far end, the cutter reversed.
    let (ro, ri, rc) = (1.75, 1.534, 3.0);
    let (tilt, turn) = (40_f64.to_radians(), 25_f64.to_radians());
    let d = Vec3::new(turn.sin() * tilt.cos(), tilt.sin(), turn.cos() * tilt.cos());
    let body = pipe(ro, ri, 20.);
    let start = measure(&body, TOL).unwrap().volume;
    let q = Vec3::new(0.3 * turn.cos(), 0., -0.3 * turn.sin());
    let result = gc_kernel::boolean::cut(&body, &cutter(q, d, rc, -12., 12.), TOL).unwrap();
    check(&result, TOL).unwrap();
    let got = start - measure(&result, TOL).unwrap().volume;
    let exact = annulus(ri, ro, &|x, z| y_on(x, z, q, d, rc));
    assert!((got - exact).abs() < 1e-9, "removed {got}, annulus integral {exact}");
    // The far end: mirror y about 10.
    let far = Vec3::new(q.x, 20., q.z);
    let dm = Vec3::new(d.x, -d.y, d.z);
    let result = gc_kernel::boolean::cut(&body, &cutter(far, dm.scale(-1.), rc, -12., 12.), TOL).unwrap();
    check(&result, TOL).unwrap();
    let got = start - measure(&result, TOL).unwrap().volume;
    assert!((got - exact).abs() < 1e-9, "far end removed {got}, annulus integral {exact}");
}

#[test]
fn oblique_copes_outside_the_domain_are_refused_by_name() {
    let (ro, ri, rc) = (1.75, 1.534, 3.0);
    let tilt = 40_f64.to_radians();
    let d = Vec3::new(0., tilt.sin(), tilt.cos());
    let body = pipe(ro, ri, 20.);
    let refused = |tool: Body| match gc_kernel::boolean::cut(&body, &tool, TOL) {
        Err(KernelError::Unsupported(m)) => m,
        other => panic!("expected a named refusal, got {other:?}"),
    };
    // Too short: its end cap lies among the shared points.
    assert!(refused(cutter(Vec3::new(0., 0., 0.), d, rc, -1., 12.)).contains("cutter's ends"));
    // Not holding the whole end: a small cutter through the end.
    assert!(refused(cutter(Vec3::new(0., 0., 0.), d, 1.0, -12., 12.)).contains("whole pipe end"));
    // A pipe that already has a hole: not plain.
    let holed = {
        let bore = cutter(Vec3::new(0., 10., 0.), Vec3::new(0., 0., -1.), 0.25, -3., 0.);
        gc_kernel::boolean::cut(&body, &bore, TOL).unwrap()
    };
    match gc_kernel::boolean::cut(&holed, &cutter(Vec3::new(0., 0., 0.), d, rc, -12., 12.), TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("plain pipe"), "{m}"),
        other => panic!("expected a named refusal, got {other:?}"),
    }
}

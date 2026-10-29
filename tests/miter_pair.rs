//! Second miter on stock already mitered at its other end (Bob's linear-nesting parts, GCmaker native survey
//! 2026-10-29: HSS2x2 36 Miter Both Trapezoid, L3 Angle 48 Miter Both Ends). The bodies and tools are the worker's
//! own inputs from the GC Kernel view (body after the first miter, then the second miter's extruded wedge).
//!
//! The wedge removes { 0 <= y <= x - c } from the stock, c = 0.00508 (its oblique face passes through
//! (3.00004, 2.99496) and (2.00004, 1.99496) at 45 degrees). The removed volume is therefore the section integral
//! of (x - c) h(x) dx over x > c, h the stock section height at x, from the shapes' saved dimensions:
//! HSS2X2X3/16 outer corner radius 0.25, inner 0.0625; L3X3X1/4 root fillet 0.25, toe radii 0.125.
//! The circular pieces are integrated exactly: with u = x - a,
//! integral (u + a - c) sqrt(r^2 - u^2) du = -(r^2 - u^2)^1.5 / 3 + (a - c) (u sqrt(r^2 - u^2) + r^2 asin(u / r)) / 2.

use gc_kernel::mass::measure;
use gc_kernel::neutral::{Body, Reader};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_CUT};
use gc_kernel::topology::check;

const TOL: f64 = 1e-6;
const C: f64 = 0.00508;

fn values(name: &str) -> Vec<f64> {
    let bytes = std::fs::read(format!("{}/tests/fixtures/miter_pair/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap();
    bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect()
}

fn read(data: &[f64]) -> Body { Reader::with_tolerance(data, TOL).body().unwrap() }

fn cut(body: &[f64], tool: &[f64]) -> Vec<f64> {
    let data: Vec<f64> = body.iter().chain(tool).copied().collect();
    handle(&parse_request(&encode_request(OP_CUT, TOL, &data)).unwrap()).unwrap()
}

/// integral over [x0, x1] of (x - C) k dx.
fn line(x0: f64, x1: f64, k: f64) -> f64 { k * ((x1 - C).powi(2) - (x0 - C).powi(2)) / 2. }

/// integral over [x0, x1] of (x - C) sqrt(r^2 - (x - a)^2) dx, exactly.
fn circle(x0: f64, x1: f64, a: f64, r: f64) -> f64 {
    let f = |x: f64| {
        let u = ((x - a) / r).clamp(-1., 1.) * r;
        let s = (r * r - u * u).max(0.).sqrt();
        -(s * s * s) / 3. + (a - C) * (u * s + r * r * (u / r).asin()) / 2.
    };
    f(x1) - f(x0)
}

/// Cuts Bob's second miter and checks the result: sound closed topology, and removed volume equal to `removed`.
fn second_miter(part: &str, removed: f64) {
    let body = values(&format!("{part}/body.f64"));
    let after = cut(&body, &values(&format!("{part}/tool.f64")));
    let (before, result) = (read(&body), read(&after));
    check(&result, TOL).unwrap_or_else(|e| panic!("{part}: result topology {e:?}"));
    let got = measure(&before, TOL).unwrap().volume - measure(&result, TOL).unwrap().volume;
    assert!((got - removed).abs() < 1e-9, "{part}: removed {got}, section integral {removed}");
}

#[test]
fn hss_tube_second_miter_removes_its_section_integral() {
    // Outer: 2 x 2, corners R 0.25; inner: x in [0.1875, 1.8125], 1.625 high, corners 0.0625 about the same centres.
    let outer = line(C, 0.25, 1.5) + 2. * circle(C, 0.25, 0.25, 0.25) + line(0.25, 1.75, 2.)
        + line(1.75, 2., 1.5) + 2. * circle(1.75, 2., 1.75, 0.25);
    let inner = line(0.1875, 0.25, 1.5) + 2. * circle(0.1875, 0.25, 0.25, 0.0625) + line(0.25, 1.75, 1.625)
        + line(1.75, 1.8125, 1.5) + 2. * circle(1.75, 1.8125, 1.75, 0.0625);
    second_miter("hss_trapezoid", outer - inner);
}

#[test]
fn angle_second_miter_removes_its_section_integral() {
    // Heel at x 0, z 0: vertical leg x in [0, 0.25] to z -3 (toe radius 0.125 about (0.125, -2.875)), horizontal leg
    // z in [-0.25, 0] to x 3 (toe about (2.875, -0.125)), root fillet radius 0.25 about (0.5, -0.5).
    let removed = line(C, 0.125, 3.)
        + line(0.125, 0.25, 2.875) + circle(0.125, 0.25, 0.125, 0.125)
        + line(0.25, 0.5, 0.5) - circle(0.25, 0.5, 0.5, 0.25)
        + line(0.5, 2.875, 0.25)
        + line(2.875, 3., 0.125) + circle(2.875, 3., 2.875, 0.125);
    second_miter("angle_both", removed);
}

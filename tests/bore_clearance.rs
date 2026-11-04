//! Bores through stock whose caps are off the axis in the last bits (Claude 2026-11-05, Bob 11-04: "fill in any gaps
//! that EyeShot is doing"). Valley Industrial B24 (C6x13 channel, NC1 import): its first 0.875 hole through the web. The
//! channel's inner web face carries the normal (-2.6e-17, 0, -1) from the filleted, tapered stock profile; the round
//! passage now takes it as the cap it is, and a planar side face inside the corridor's box is excluded when it lies
//! provably farther than the radius from the corridor's axis.

use gc_kernel::boolean::{cut, cut_recorded};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{Body, Reader};

const TOL: f64 = 1e-6;

fn read(name: &str) -> Body {
    let bytes = std::fs::read(format!("{}/tests/fixtures/parts/channel_web_hole/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap();
    let values: Vec<f64> = bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
    Reader::with_tolerance(&values, TOL).body().unwrap()
}

#[test]
fn a_web_hole_through_a_near_cardinal_cap_is_cut() {
    let (body, tool) = (read("body.f64"), read("tool.f64"));
    let before = measure(&body, TOL).unwrap().volume;
    // The rebuild's named cut (op 26) and the plain cut give the same body.
    let (recorded, _) = cut_recorded(&body, &tool, TOL).unwrap();
    let after = measure(&recorded, TOL).unwrap().volume;
    assert!((measure(&cut(&body, &tool, TOL).unwrap(), TOL).unwrap().volume - after).abs() < 1e-12);
    // The 0.875 hole (r = 0.4374 in) through the 11.1 mm web.
    let removed = before - after;
    let r = 0.43740157480314956_f64;
    let web = 11.1 / 25.4;
    assert!((removed - std::f64::consts::PI * r * r * web).abs() < 1e-9, "removed {removed}");
}

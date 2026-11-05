//! A transverse end cut on a body whose end face is stepped (Claude 2026-11-05, Bob 11-04: "fill in any gaps that
//! EyeShot is doing"). Bob Cody 12004SL0 (C12x20.7 channel, NC1 import), its second AK end cutter: the first left the
//! end face an L-shaped outline (two material sections meeting at one station). The transverse coverage check took
//! every face parallel to the partition planes as a rectangle and refused it ("end-cut transverse patch orientation or
//! coverage is invalid"); it now takes any outline of edges parallel to the axes, sampled on all its coordinates.

use gc_kernel::boolean::{cut, cut_recorded};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{Body, Reader};

const TOL: f64 = 1e-6;

fn read(name: &str) -> Body {
    let bytes = std::fs::read(format!("{}/tests/fixtures/parts/channel_stepped_end/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap();
    let values: Vec<f64> = bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
    Reader::with_tolerance(&values, TOL).body().unwrap()
}

#[test]
fn an_end_cut_beside_a_stepped_end_face_is_cut() {
    let (body, tool) = (read("body.f64"), read("tool.f64"));
    let (recorded, _) = cut_recorded(&body, &tool, TOL).unwrap();
    // One closed solid. The end cut takes away the material round the first cutter's opening: genus 1 to 0.
    let genus = |b: &Body| gc_kernel::topology::check(b, TOL).unwrap().genus();
    assert!(genus(&recorded) <= genus(&body));
    let (before, after) = (measure(&body, TOL).unwrap().volume, measure(&recorded, TOL).unwrap().volume);
    assert!((measure(&cut(&body, &tool, TOL).unwrap(), TOL).unwrap().volume - after).abs() < 1e-12);
    // Material is removed, none added; the value is this fixture's result. The whole part (all five AK cutters) is
    // GC Maker's production solid within its mesh error.
    assert!(after < before, "{after} of {before}");
    assert!((after - 40.157983135437796).abs() < 1e-8, "{after}");
}

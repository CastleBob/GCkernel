//! A transverse end cut on stock whose faces are off the coordinate axes in the last bits (Claude 2026-11-05, Bob
//! 11-04: "fill in any gaps that EyeShot is doing"). Arc Steel B2025 (C12x20.7 channel, NC1 import), its first AK end
//! cutter: the channel's flange end faces carry normals (-1, 1.3e-16, 0) from the tapered, filleted stock profile, so
//! the transverse partition read them as non-transverse and refused ("an end-cut surface patch is below the admitted
//! partition resolution"). The transverse domain now takes directions within 1e-12 of an axis as that axis
//! (endcut::snap_directions).

use gc_kernel::boolean::{cut, cut_recorded};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{Body, Reader};

const TOL: f64 = 1e-6;

fn read(name: &str) -> Body {
    let bytes = std::fs::read(format!("{}/tests/fixtures/parts/channel_end_partition/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap();
    let values: Vec<f64> = bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
    Reader::with_tolerance(&values, TOL).body().unwrap()
}

#[test]
fn a_channel_end_cut_with_near_cardinal_end_faces_is_cut() {
    let (body, tool) = (read("body.f64"), read("tool.f64"));
    let _ = gc_kernel::boolean::take_equivalences();
    let (recorded, _) = cut_recorded(&body, &tool, TOL).unwrap();
    // Bob 10-28: the status says so. The cut reports the equivalence it used.
    let notes = gc_kernel::boolean::take_equivalences();
    assert!(notes.iter().any(|n| n.contains("taken as exact")), "{notes:?}");
    let plain = cut(&body, &tool, TOL).unwrap();
    let genus = |b: &Body| gc_kernel::topology::check(b, TOL).unwrap().genus();
    assert_eq!(genus(&recorded), genus(&body));
    let (before, after) = (measure(&body, TOL).unwrap().volume, measure(&recorded, TOL).unwrap().volume);
    assert!((measure(&plain, TOL).unwrap().volume - after).abs() < 1e-12);
    // Material is removed at the end, none added; the value is this fixture's result (281.576718641662 in3).
    assert!(after < before, "{after} of {before}");
    assert!((after - 281.57671864166207).abs() < 1e-8, "{after}");
}

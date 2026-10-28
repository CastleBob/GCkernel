//! Touching within tolerance (Bob 2026-10-28: what agrees within the tolerance is the same): a hole whose wall
//! reaches exactly the I-beam flange's inside face touches both root fillets there, so its CylinderPair edges end at
//! a double root of their discriminant (Claude, geom::pair_trim_touches).

use std::f64::consts::PI;

use gc_kernel::geom::{pair_trim_touches, Curve, Vec3};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{Body, Reader};
use gc_kernel::protocol::{handle, parse_request, encode_request, OP_CUT};

const TOL: f64 = 1e-6;

fn values(name: &str) -> Vec<f64> {
    let bytes = std::fs::read(format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap();
    bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect()
}

fn cut(body: &[f64], tool: &[f64]) -> Vec<f64> {
    let data: Vec<f64> = body.iter().chain(tool).copied().collect();
    handle(&parse_request(&encode_request(OP_CUT, TOL, &data)).unwrap()).unwrap()
}

fn read(data: &[f64]) -> Body { Reader::with_tolerance(data, TOL).body().unwrap() }

/// Bob's captured part (GCmaker 2026-10-28): rat holes, then a hole at x 2.332133337554964, radius 2.092133337554964,
/// so its wall reaches x = .24 + 2e-16. Removed by the hole, an independent integral: the web .23 pi R^2 plus
/// 2 x integral over x .24-.74 of 2 sqrt(R^2 - (x - cx)^2) (.5 - sqrt(.25 - (x - .74)^2)) = 3.29255763471708.
#[test]
fn bobs_hole_touching_the_flange_face_cuts_exactly() {
    let mut body = values("tangent_hole/base.f64");
    for step in ["step01", "step02"] { body = cut(&body, &values(&format!("tangent_hole/{step}.f64"))); }
    let before = measure(&read(&body), TOL).unwrap().volume;
    let after = cut(&body, &values("tangent_hole/step03.f64"));
    let result = read(&after);
    let removed = before - measure(&result, TOL).unwrap().volume;
    assert!((removed - 3.29255763471708).abs() < 1e-9, "removed {removed}");

    // Four hole/fillet edges (two fillets, two sides of the hole), each touching at the flange end.
    let mut touching = 0;
    for e in &result.edges {
        if let Curve::CylinderPair { c, a, x, rho, p, b, r, t0, t1, .. } = e.curve {
            let ends = pair_trim_touches(c, a, x, rho, p, b, r, t0, t1, TOL).unwrap();
            if ends.0 || ends.1 {
                touching += 1;
                let at = if ends.0 { t0 } else { t1 };
                let (q, d) = e.curve.eval(at);
                assert!((q.x - 0.24).abs() < 1e-9 && d.is_finite() && d.length() > 0.1, "{q:?} {d:?}");
            }
        }
    }
    assert_eq!(touching, 4);
    // The strict reader (no tolerance) refuses those edges; the request's tolerance admits them.
    assert!(Reader::new(&after).body().is_err());
}

/// The hole/fillet pair in Bob's numbers: rulings x = cx + R cos t of the hole meet the fillet (axis y through
/// x .74, z -2.615, radius .5) where |x - .74| <= .5; at t = pi, x = .24 and S has a double root.
fn hole_on_fillet(cx: f64, rho: f64, t0: f64, t1: f64) -> Option<(bool, bool)> {
    let (c, a, x) = (Vec3::new(cx, 6.442907, -2.615), Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0));
    let (p, b) = (Vec3::new(0.74, 0.0, -2.615), Vec3::new(0.0, 1.0, 0.0));
    pair_trim_touches(c, a, x, rho, p, b, 0.5, t0, t1, TOL)
}

#[test]
fn a_touch_is_admitted_only_at_an_end_and_only_as_a_double_root() {
    let r = 2.092133337554964;
    assert_eq!(hole_on_fillet(0.24 + r, r, PI - 0.7, PI), Some((false, true)));
    assert_eq!(hole_on_fillet(0.24 + r, r, PI - 0.7, PI - 0.1), Some((false, false)));
    // Touching inside the trim (a corner within the edge): refused.
    assert_eq!(hole_on_fillet(0.24 + r, r, PI - 0.7, PI + 0.7), None);
    // Overlapping the tangent line by 1e-7 in (inside tolerance): still a touch; by 1e-4 in: refused. A 1e-7 gap
    // leaves S proved positive: an ordinary, steep trim.
    assert_eq!(hole_on_fillet(0.24 + r - 1e-7, r, PI - 0.7, PI), Some((false, true)));
    assert_eq!(hole_on_fillet(0.24 + r - 1e-4, r, PI - 0.7, PI), None);
    assert_eq!(hole_on_fillet(0.24 + r + 1e-7, r, PI - 0.7, PI), Some((false, false)));
    // A simple root (the hole's ruling grazes the fillet while the surfaces cross, Astra's large-hole case): the
    // speed is unbounded in this chart, refused.
    let graze = (5.0f64 / 7.0).acos();
    assert_eq!(hole_on_fillet(0.74, 0.7, graze, PI / 2.0), None);
}

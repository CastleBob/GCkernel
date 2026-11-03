//! Bores across cap seams (Claude 2026-11-03, Bob's angle corner cutout): a transverse cut leaves the vertical leg's flat
//! faces split on one plane; a later side hole crossing those seams is admitted by joining them (passage::join_cap_seams),
//! with no point moved, and records no history for that cut.

use gc_kernel::mass::measure;
use gc_kernel::neutral::{Body, Reader};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_CUT, OP_CUT_RECORDED};

const TOL: f64 = 1e-6;

fn values(name: &str) -> Vec<f64> {
    let bytes = std::fs::read(format!("{}/tests/fixtures/parts/angle_corner_cutout/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap();
    bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect()
}
fn run(op: u32, data: &[f64]) -> Result<Vec<f64>, String> {
    handle(&parse_request(&encode_request(op, TOL, data)).unwrap()).map_err(|e| e.to_string())
}
fn read(data: &[f64]) -> Body { Reader::with_tolerance(data, TOL).body().unwrap() }

#[test]
fn side_hole_across_the_corner_cut_seams_completes_the_angle() {
    let mut body = values("base.f64");
    let mut volumes = vec![measure(&read(&body), TOL).unwrap().volume];
    for step in 1..=4 {
        let tool = values(&format!("step{step:02}.f64"));
        body = run(OP_CUT, &body.iter().chain(&tool).copied().collect::<Vec<_>>()).unwrap();
        volumes.push(measure(&read(&body), TOL).unwrap().volume);
    }
    // Each cut removes material; the last removes exactly the 1 inch hole through the 1/4 inch leg.
    assert!(volumes.windows(2).all(|w| w[1] < w[0]));
    let hole = std::f64::consts::PI * 0.25 * 0.25;
    assert!((volumes[3] - volumes[4] - hole).abs() < 1e-9, "{volumes:?}");
    assert!((volumes[4] - 3.656972).abs() < 1e-5, "{volumes:?}");
}

#[test]
fn seamed_bore_runs_once_and_records_no_history() {
    let mut body = values("base.f64");
    for step in 1..=3 {
        let tool = values(&format!("step{step:02}.f64"));
        body = run(OP_CUT, &body.iter().chain(&tool).copied().collect::<Vec<_>>()).unwrap();
    }
    let request: Vec<f64> = body.iter().chain(&values("step04.f64")).copied().collect();
    let plain = run(OP_CUT, &request).unwrap();
    let recorded = run(OP_CUT_RECORDED, &request).unwrap();
    let n = plain.len();
    assert_eq!(recorded[..n], plain[..]);
    assert_eq!(recorded[n], 0.0);
    let reason: String = recorded[n + 2..].iter().map(|&b| b as u8 as char).collect();
    assert!(reason.contains("cap seams"), "{reason}");
}

//! Exact saved SL0 final-cut operands, captured by GCcore 5ed55c6 at e0a0105.
//! This records a current refusal, not successful completed channel geometry.
use gc_kernel::{
    neutral::{write_body, Body, Reader},
    KernelError,
};
const TOL: f64 = 1e-6;
fn body(text: &str) -> Body {
    let data: Vec<f64> = text
        .split_whitespace()
        .map(|s| s.parse().unwrap())
        .collect();
    Reader::with_tolerance(&data, TOL).body().unwrap()
}
#[test]
fn saved_channel_final_cut_retains_valid_operands_on_explicit_coincidence_refusal() {
    let before = body(include_str!("fixtures/channel-sl0-before.brep.txt"));
    let tool = body(include_str!("fixtures/channel-sl0-tool.brep.txt"));
    gc_kernel::topology::check(&before, TOL).unwrap();
    gc_kernel::topology::check(&tool, TOL).unwrap();
    let original = (write_body(&before), write_body(&tool));
    let error = gc_kernel::boolean::cut(&before, &tool, TOL).unwrap_err();
    assert!(
        matches!(error,KernelError::Unsupported(ref e) if e=="transverse graph coincidence is unresolved"),
        "{error}"
    );
    assert_eq!(write_body(&before), original.0);
    assert_eq!(write_body(&tool), original.1);
}

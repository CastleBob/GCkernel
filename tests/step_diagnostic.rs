use gc_kernel::step::{self, diagnostic};

#[test]
fn diagnostic_refuses_malformed_input_and_invalid_tolerance() {
    assert!(diagnostic::inspect(b"broken", 1e-6).is_err());
    assert!(diagnostic::inspect(b"", f64::NAN).is_err());
    assert!(diagnostic::inspect(b"", 0.).is_err());
}

#[test]
fn original_seven_failures_are_present_in_the_original_basis() {
    let Ok(root) = std::env::var("GC_STEP_CORPUS") else {
        return;
    };
    for (name, pairs) in [
        ("Ellipse Pipe.step", 12),
        ("JaysPipe.step", 12),
        ("New Pipe.step", 20),
        ("Pipe 7in Big Ellipse.step", 4),
        ("Pipe_10in_Slot_3holes.stp", 32),
        ("test.step", 12),
        ("TubePipe.step", 8),
    ] {
        let path = std::path::Path::new(&root).join(name);
        let before = std::fs::read(&path).unwrap();
        let report = diagnostic::inspect(&before, 1e-6).unwrap();
        assert!(!report.admitted, "{name}");
        assert!(report.admission_error.is_some());
        assert_eq!(report.mismatches.len(), pairs, "{name}");
        assert_eq!(report.body_errors.len(), 1, "{name}");
        for r in report.mismatches {
            assert!(
                r.source_observed_departure.unwrap() > r.tolerance,
                "{name} #{:?}",
                r.edge
            );
            assert!(
                r.conversion_difference.unwrap() < 1e-10,
                "{name} #{:?}",
                r.edge
            );
            assert!((r.observed_departure - r.source_observed_departure.unwrap()).abs() < 1e-10);
            assert!(r.solid > 0 && r.face > 0 && r.edge > 0 && r.curve > 0 && r.surface > 0);
        }
        assert!(step::read(&before, 1e-6).is_err());
        assert_eq!(before, std::fs::read(&path).unwrap());
    }
}

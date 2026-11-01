//! Every failed edge/support pair, with source IDs and independent source-basis
//! observations. Samples diagnose failures; they never authorize admission.
fn main() {
    println!("file,solid,face,edge,curve,surface,tolerance_inches,observed_inches,source_basis_inches,conversion_difference_inches,source_to_inches");
    for path in std::env::args().skip(1) {
        match gc_kernel::step::diagnostic::inspect_file(std::path::Path::new(&path), 1e-6) {
            Ok(report) => {
                eprintln!(
                    "{path}: admitted={}, mismatches={}, {}",
                    report.admitted,
                    report.mismatches.len(),
                    report.admission_error.unwrap_or_default()
                );
                for r in report.mismatches {
                    println!(
                        "\"{}\",{},{},{},{},{},{:.17e},{:.17e},{},{},{:.17e}",
                        path.replace('"', "\"\""),
                        r.solid,
                        r.face,
                        r.edge,
                        r.curve,
                        r.surface,
                        r.tolerance,
                        r.observed_departure,
                        r.source_observed_departure
                            .map(|x| format!("{x:.17e}"))
                            .unwrap_or_default(),
                        r.conversion_difference
                            .map(|x| format!("{x:.17e}"))
                            .unwrap_or_default(),
                        r.source_length_to_inches
                    );
                }
            }
            Err(e) => eprintln!("{path}: diagnostic failed: {e}"),
        }
    }
}

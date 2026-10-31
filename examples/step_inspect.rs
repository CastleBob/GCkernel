//! Explicit development intake: cargo run --example step_inspect -- file.step ...
//! Reports exact imported bodies; never changes the files or installs a worker.
fn main() {
    let mut failed = false;
    let mut args = std::env::args().skip(1).peekable();
    let tolerance = if args.peek().is_some_and(|a| a == "--tolerance") {
        args.next();
        match args
            .next()
            .and_then(|s| s.parse::<f64>().ok())
            .filter(|n| n.is_finite() && *n > 0.0)
        {
            Some(n) => n,
            None => {
                eprintln!("--tolerance requires a positive finite number in inches");
                std::process::exit(2);
            }
        }
    } else {
        1e-6
    };
    println!("Requested tolerance: {tolerance} inches (no automatic retries)");
    for path in args {
        match gc_kernel::step::read_file(std::path::Path::new(&path), tolerance) {
            Ok(import) => {
                println!("OK {path}: {} bodies", import.bodies.len());
                for b in import.bodies {
                    println!(
                        "  #{} {:?}: V/E/F={}/{}/{}, source-to-inch={}, volume={} in^3",
                        b.entity_id,
                        b.name,
                        b.body.vertices.len(),
                        b.body.edges.len(),
                        b.body.faces.len(),
                        b.source_length_to_inches,
                        b.volume
                    );
                }
            }
            Err(e) => {
                failed = true;
                println!("REFUSED {path}: {e}");
            }
        }
    }
    if failed {
        std::process::exit(1);
    }
}

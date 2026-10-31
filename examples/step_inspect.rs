//! Explicit development intake: cargo run --example step_inspect -- file.step ...
//! Reports exact imported bodies; never changes the files or installs a worker.
fn main() {
    let mut failed = false;
    for path in std::env::args().skip(1) {
        match gc_kernel::step::read_file(std::path::Path::new(&path), 1e-6) {
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

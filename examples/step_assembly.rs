fn main() {
    for path in std::env::args().skip(1) {
        match gc_kernel::step::assembly::read_file(std::path::Path::new(&path), 1e-6) {
            Ok(items) => {
                println!("{path}: {} occurrences", items.len());
                for i in items {
                    println!(
                        "  {:?}: solid #{}, volume {:.12}",
                        i.path, i.item.entity_id, i.item.volume
                    );
                }
            }
            Err(error) => eprintln!("{path}: {error}"),
        }
    }
}

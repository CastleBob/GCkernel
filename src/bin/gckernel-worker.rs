//! One GC Kernel request per process: the request on stdin, the reply on stdout. A panic is caught and
//! answered as an internal error; nothing unwinds out of the process into GCmaker.

use std::io::{Read, Write};

use gc_kernel::protocol::{encode_reply, handle, parse_request};
use gc_kernel::KernelError;

fn main() {
    let mut input = Vec::new();
    let result = match std::io::stdin().read_to_end(&mut input) {
        Err(e) => Err(KernelError::InvalidInput(format!("could not read the request: {e}"))),
        Ok(_) => std::panic::catch_unwind(|| parse_request(&input).and_then(|r| handle(&r)))
            .unwrap_or_else(|panic| {
                let what = panic.downcast_ref::<&str>().map(|s| s.to_string())
                    .or_else(|| panic.downcast_ref::<String>().cloned())
                    .unwrap_or_else(|| "unknown".into());
                Err(KernelError::Internal(format!("the kernel panicked: {what}")))
            }),
    };
    let mut out = std::io::stdout().lock();
    let _ = out.write_all(&encode_reply(&result));
    let _ = out.flush();
}

//! GC Kernel: GCmaker's own solid-modelling kernel (Bob 2026-10-27). Clean room: our own code from public
//! mathematics, no third-party geometry source. Plan and ownership: GCm GCmaker/Handoff/Claude-GcKernel-BattlePlan.md
//! (foundation: Claude; src/boolean/: Astra).
//!
//! GCmaker talks to it through `gckernel-worker.exe`, one request per process (see `protocol`), so a
//! failure here cannot affect an open GCmaker.

pub mod boolean;
pub mod extrude;
pub mod geom;
pub mod mass;
pub mod neutral;
pub mod protocol;
pub mod topology;

/// Why an operation produced no result. Nothing partial is ever returned with an error.
#[derive(Clone, Debug, PartialEq)]
pub enum KernelError {
    /// Valid input the kernel cannot handle yet (a surface type, a case outside a supported domain).
    Unsupported(String),
    /// Malformed or inconsistent input.
    InvalidInput(String),
    /// A fault in the kernel itself.
    Internal(String),
}

impl std::fmt::Display for KernelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KernelError::Unsupported(m) => write!(f, "unsupported: {m}"),
            KernelError::InvalidInput(m) => write!(f, "invalid input: {m}"),
            KernelError::Internal(m) => write!(f, "internal error: {m}"),
        }
    }
}

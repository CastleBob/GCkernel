//! Independent boolean operations. Owner: Astra/Codex.
//! The initial supported domain is closed axis-aligned rectilinear solids.
mod rectilinear;

pub fn cut(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
) -> Result<crate::neutral::Body, crate::KernelError> {
    rectilinear::cut(body, tool, tolerance)
}

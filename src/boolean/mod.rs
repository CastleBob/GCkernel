//! Independent boolean operations. Owner: Astra/Codex.
//! Supported domains are rectilinear solids and contained circular through-bores
//! in plain axis-aligned boxes. Admission selects one operation before execution;
//! failure never triggers another algorithm.
mod bore;
mod rectilinear;

pub fn cut(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
) -> Result<crate::neutral::Body, crate::KernelError> {
    if tool
        .faces
        .iter()
        .any(|f| matches!(f.surface, crate::geom::Surface::Cylinder { .. }))
    {
        bore::cut(body, tool, tolerance)
    } else {
        rectilinear::cut(body, tool, tolerance)
    }
}

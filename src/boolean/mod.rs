//! Independent boolean operations. Owner: Astra/Codex.
//! Supported domains are rectilinear solids, contained circular through-bores,
//! and capsule through-slots in plain axis-aligned boxes. Admission selects one operation before execution;
//! failure never triggers another algorithm.
mod bore;
mod capsule;
mod rectilinear;

pub fn cut(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
) -> Result<crate::neutral::Body, crate::KernelError> {
    let cylinders = tool
        .faces
        .iter()
        .filter(|f| matches!(f.surface, crate::geom::Surface::Cylinder { .. }))
        .count();
    match cylinders {
        0 => rectilinear::cut(body, tool, tolerance),
        1 => bore::cut(body, tool, tolerance),
        2 => capsule::cut(body, tool, tolerance),
        _ => Err(crate::KernelError::Unsupported(
            "this curved cutter is outside the supported cut domains".into(),
        )),
    }
}

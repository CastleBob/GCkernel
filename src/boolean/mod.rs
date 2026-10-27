//! Independent boolean operations. Owner: Astra/Codex.
//! Supported domains are rectilinear solids, clear circular passages in planar
//! regions of solids, and capsule through-slots in plain axis-aligned boxes. Admission selects one operation before execution;
//! failure never triggers another algorithm.
mod bore;
mod capsule;
mod endcut;
mod passage;
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
        1 if tool.vertices.len() == 2 && tool.edges.len() == 3 && tool.faces.len() == 3 => {
            bore::cut(body, tool, tolerance)
        }
        1 => endcut::cut(body, tool, tolerance),
        2 => capsule::cut(body, tool, tolerance),
        _ => Err(crate::KernelError::Unsupported(
            "this curved cutter is outside the supported cut domains".into(),
        )),
    }
}

//! Independent boolean operations. Owner: Astra/Codex.
//! Supported domains are rectilinear solids, clear planar passages, and bounded
//! transverse intersections of analytic plane/cylinder solids. Admission selects one operation before execution;
//! failure never triggers another algorithm.
mod bore;
mod box_cut;
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
        0 if body
            .faces
            .iter()
            .any(|f| matches!(f.surface, crate::geom::Surface::Cylinder { .. })) =>
        {
            box_cut::cut(body, tool, tolerance)
        }
        0 => rectilinear::cut(body, tool, tolerance),
        1 if tool.vertices.len() == 2 && tool.edges.len() == 3 && tool.faces.len() == 3 => {
            bore::cut(body, tool, tolerance)
        }
        1 => endcut::cross_cut(body, tool, tolerance),
        2 => capsule::cut(body, tool, tolerance),
        _ => Err(crate::KernelError::Unsupported(
            "this curved cutter is outside the supported cut domains".into(),
        )),
    }
}

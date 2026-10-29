//! Independent boolean operations. Owner: Claude (from 2026-10-29; previously Astra/Codex).
//! Supported domains are rectilinear solids, clear planar passages, and bounded
//! transverse intersections of analytic plane/cylinder solids. Admission selects one operation before execution;
//! failure never triggers another algorithm.
mod bore;
mod box_cut;
mod capsule;
mod endcut;
mod passage;
mod radial;
mod rectilinear;
mod section;

pub fn cut(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
) -> Result<crate::neutral::Body, crate::KernelError> {
    // No boolean domain handles conical faces yet (tapered Pole stock): refused by name before any selector reads them.
    if body.faces.iter().chain(&tool.faces).any(|f| matches!(f.surface, crate::geom::Surface::Cone { .. })) {
        return Err(crate::KernelError::Unsupported("cuts on bodies with conical faces are not supported yet".into()));
    }
    if body.faces.iter().chain(&tool.faces).any(|f| matches!(f.surface, crate::geom::Surface::Sphere { .. } | crate::geom::Surface::Torus { .. } | crate::geom::Surface::Revolution { .. })) {
        return Err(crate::KernelError::Unsupported("cuts on bodies with spherical, toroidal or general revolved faces are not supported yet".into()));
    }
    // Through cuts reaching a prism's outline or holes: the section domain (selected before execution).
    if section::selects(body, tool, tolerance) {
        return section::cut(body, tool, tolerance);
    }
    // Round cutters crossing a pipe wall radially, at any clock angle (selected before execution).
    if radial::selects(body, tool, tolerance) {
        return radial::cut(body, tool, tolerance);
    }
    let cylinders = tool
        .faces
        .iter()
        .filter(|f| matches!(f.surface, crate::geom::Surface::Cylinder { .. }))
        .count();
    match cylinders {
        0 if body
            .faces
            .iter()
            .any(|f| matches!(f.surface, crate::geom::Surface::Cylinder { .. }))
            || tool.faces.iter().any(|f| {
                matches!(f.surface, crate::geom::Surface::Plane { n, .. }
                if [n.x,n.y,n.z].iter().filter(|v| **v != 0.).count() != 1)
            }) =>
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

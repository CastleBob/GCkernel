//! Independent boolean operations. Owner: Claude (from 2026-10-29; previously Astra/Codex).
//! Supported domains are rectilinear solids, clear planar passages, and bounded
//! transverse intersections of analytic plane/cylinder solids. Admission selects one operation before execution;
//! failure never triggers another algorithm.
mod bore;
mod box_cut;
mod capsule;
mod cone_bore;
mod cope;
mod endcut;
mod normal_end;
mod passage;
mod radial;
mod rectilinear;
mod section;
mod unite;
mod unite_tidy;

pub use unite::{union, union_named};

#[derive(Clone, Copy, Debug)]
enum Domain {
    Section,
    Radial,
    BoxCut,
    Rectilinear,
    Bore,
    Cross,
    Capsule,
    Cope,
    NormalEnd,
    ConeBore,
}

/// The exact drill the bore domain takes: one cylinder wall between two caps, one seam.
fn drill(tool: &crate::neutral::Body) -> bool {
    tool.vertices.len() == 2 && tool.edges.len() == 3 && tool.faces.len() == 3
}

/// The one domain a cut runs in, chosen from the operands before anything is built.
fn select(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
) -> Result<Domain, crate::KernelError> {
    if body.faces.iter().chain(&tool.faces).any(|f| matches!(f.surface, crate::geom::Surface::Spline { .. })) {
        return Err(crate::KernelError::Unsupported("cuts on general spline surfaces are not supported yet".into()));
    }
    if body.edges.iter().chain(&tool.edges).any(|e| matches!(e.curve, crate::geom::Curve::RationalBezier { .. })) {
        return Err(crate::KernelError::Unsupported("cuts on rational Bezier boundaries are not supported yet".into()));
    }
    // Conical faces (tapered Pole stock): only radial bores through the walls (cone_bore.rs, Claude 2026-10-31); its
    // admission names why anything else is refused, before any other selector reads them.
    if tool.faces.iter().any(|f| matches!(f.surface, crate::geom::Surface::Cone { .. })) {
        return Err(crate::KernelError::Unsupported("cuts with tools that have conical faces are not supported yet".into()));
    }
    if body.faces.iter().any(|f| matches!(f.surface, crate::geom::Surface::Cone { .. })) {
        cone_bore::admit(body, tool, tolerance)?;
        return Ok(Domain::ConeBore);
    }
    if body.faces.iter().chain(&tool.faces).any(|f| matches!(f.surface, crate::geom::Surface::Sphere { .. } | crate::geom::Surface::Torus { .. } | crate::geom::Surface::Revolution { .. })) {
        return Err(crate::KernelError::Unsupported("cuts on bodies with spherical, toroidal or general revolved faces are not supported yet".into()));
    }
    // A normalized-end waste band (op 23's tool) takes its own domain; a body with a normalized end only the radial bore
    // domain, which bounds the end's rims piece by piece (Claude 2026-10-31).
    if tool.faces.iter().any(|f| matches!(f.surface, crate::geom::Surface::Conoid { .. })) {
        return Ok(Domain::NormalEnd);
    }
    if body.faces.iter().any(|f| matches!(f.surface, crate::geom::Surface::Conoid { .. }))
        || body.edges.iter().chain(&tool.edges).any(|e| matches!(e.curve, crate::geom::Curve::RadialPair { .. })) {
        if radial::selects(body, tool, tolerance) {
            return Ok(Domain::Radial);
        }
        return Err(crate::KernelError::Unsupported(
            "cuts on a normalized pipe end are supported only for radial bores clear of the end".into(),
        ));
    }
    // Through cuts reaching a prism's outline or holes: the section domain (selected before execution).
    if section::selects(body, tool, tolerance) {
        return Ok(Domain::Section);
    }
    // Round cutters crossing a pipe wall radially, at any clock angle (selected before execution).
    if radial::selects(body, tool, tolerance) {
        return Ok(Domain::Radial);
    }
    // A round cutter oblique to a pipe's axis (a cope at an angle): neither the bore nor the transverse chart.
    if cope::selects(body, tool, tolerance) {
        return Ok(Domain::Cope);
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
            Ok(Domain::BoxCut)
        }
        0 => Ok(Domain::Rectilinear),
        1 if drill(tool) => Ok(Domain::Bore),
        1 => Ok(Domain::Cross),
        2 => Ok(Domain::Capsule),
        _ => Err(crate::KernelError::Unsupported(
            "this curved cutter is outside the supported cut domains".into(),
        )),
    }
}

pub fn cut(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
) -> Result<crate::neutral::Body, crate::KernelError> {
    match select(body, tool, tolerance)? {
        Domain::Section => section::cut(body, tool, tolerance),
        Domain::Radial => radial::cut(body, tool, tolerance),
        Domain::BoxCut => box_cut::cut(body, tool, tolerance),
        Domain::Rectilinear => rectilinear::cut(body, tool, tolerance),
        Domain::Bore => bore::cut(body, tool, tolerance),
        Domain::Cross => endcut::cross_cut(body, tool, tolerance),
        Domain::Capsule => capsule::cut(body, tool, tolerance),
        Domain::Cope => cope::cut(body, tool, tolerance),
        Domain::NormalEnd => normal_end::cut(body, tool, tolerance),
        Domain::ConeBore => cone_bore::cut(body, tool, tolerance),
    }
}

/// Op 12: the cut, the recorded history of its every output entity (provenance.rs roles 30-39) and the outcome of every
/// input entity. A domain that does not record its history refuses by name; nothing is matched after the fact.
pub fn cut_named(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
) -> Result<(crate::neutral::Body, crate::provenance::Provenance, Vec<crate::provenance::Outcome>), crate::KernelError> {
    use crate::provenance::*;
    let (result, history) = match select(body, tool, tolerance)? {
        Domain::Bore => bore::cut_named(body, tool, tolerance)?,
        Domain::Section => section::cut_named(body, tool, tolerance)?,
        domain => {
            return Err(crate::KernelError::Unsupported(format!(
                "cut history is not recorded for the {domain:?} cut domain yet"
            )))
        }
    };
    let counts = (history.faces.len(), history.edges.len(), history.vertices.len());
    if counts != (result.faces.len(), result.edges.len(), result.vertices.len()) {
        return Err(crate::KernelError::Internal(format!("cut history names {counts:?} entities of another count")));
    }
    // Generated entities must be told apart by their records alone.
    for list in [&history.faces, &history.edges, &history.vertices] {
        for (i, r) in list.iter().enumerate() {
            if matches!(r.role, MEETING | TOOL_SEAM | SEAM_POINT | MEETING_END | TOOL_EDGE_END) && list[i + 1..].contains(r) {
                return Err(crate::KernelError::Unsupported(format!("cut history would name two entities {r:?}")));
            }
        }
    }
    let inputs = [body.faces.len(), body.edges.len(), body.vertices.len(), tool.faces.len(), tool.edges.len(), tool.vertices.len()];
    let outcomes = outcomes(&history, inputs);
    Ok((result, history, outcomes))
}

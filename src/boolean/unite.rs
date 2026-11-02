//! Union, ops 24 and 25 (Claude 2026-11-02, Bob: boss extrudes and revolves in the GC Kernel view). The same
//! philosophy as the cut: the union is solved by the first exact domain that recognizes the configuration, chosen from
//! the operands before anything is built; a domain's failure is returned, never retried in another; anything else is
//! refused by name. Domains, in order:
//! 1. Face (`section::unite_face`): a boss standing on one planar face of any body, its base strictly inside that face
//!    and nothing of the body in its room. The body is kept entity for entity. History recorded (op 25).
//! 2. Stack (`section::unite_stack`): both operands right prisms along one cardinal axis (a boss overhanging or flush
//!    with a plate's edge, a bar through a plate, overlapping blocks). Rebuilt from the arrangement of both sections,
//!    coplanar walls and caps joined. History not recorded yet: op 25 refuses by name, op 24 succeeds.
//! 3. Overhang (`section::unite_overhang`): a boss on one planar face of a body that is not a layered stack, its base
//!    crossing the face's edge (a round boss over an I-beam flange tip). History not recorded yet.
//! A boss that does not reach the body (two solids) is refused by name.
use super::section::{unite_face, unite_overhang, unite_stack};
use crate::neutral::Body;
use crate::provenance::{outcomes, Outcome, Provenance};
use crate::KernelError;

enum Domain {
    Face(unite_face::Plan),
    Stack,
    Overhang(unite_overhang::Plan),
}

fn unsupported(message: &str) -> KernelError {
    KernelError::Unsupported(message.into())
}

/// Whether the bodies' world boxes are apart by more than the tolerance (then they cannot touch).
fn apart(body: &Body, tool: &Body, tol: f64) -> bool {
    let bounds = |b: &Body| {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for e in &b.edges {
            let (l, h) = super::passage::curve_bounds(&e.curve);
            for d in 0..3 {
                lo[d] = lo[d].min(l[d]);
                hi[d] = hi[d].max(h[d]);
            }
        }
        (lo, hi)
    };
    // Edge boxes bound plane faces only; any curved face makes this test say nothing.
    let planar = |b: &Body| b.faces.iter().all(|f| matches!(f.surface, crate::geom::Surface::Plane { .. }));
    if !(planar(body) && planar(tool)) {
        return false;
    }
    let ((l1, h1), (l2, h2)) = (bounds(body), bounds(tool));
    (0..3).any(|d| h1[d] + tol < l2[d] || h2[d] + tol < l1[d])
}

fn select(body: &Body, tool: &Body, tol: f64) -> Result<Domain, KernelError> {
    use crate::geom::{Curve, Surface};
    if body.faces.is_empty() || tool.faces.is_empty() {
        return Err(KernelError::InvalidInput("a union needs two non-empty bodies".into()));
    }
    if body.faces.iter().chain(&tool.faces).any(|f| matches!(f.surface, Surface::Spline { .. })) {
        return Err(unsupported("unions with general spline surfaces are not supported yet"));
    }
    if body.edges.iter().chain(&tool.edges).any(|e| matches!(e.curve, Curve::RationalBezier { .. } | Curve::Bezier { .. })) {
        return Err(unsupported("unions with Bezier boundaries are not supported yet"));
    }
    let face = unite_face::admit(body, tool, tol);
    if let Ok(Some(plan)) = face {
        return Ok(Domain::Face(plan));
    }
    if unite_stack::selects(body, tool, tol) {
        return Ok(Domain::Stack);
    }
    match unite_overhang::admit(body, tool, tol) {
        Ok(Some(plan)) => return Ok(Domain::Overhang(plan)),
        Err(reason) => return Err(reason),
        Ok(None) => {}
    }
    match face {
        Err(reason) => Err(reason),
        _ if apart(body, tool, tol) => Err(unsupported("the boss does not reach the body: their union would be two separate solids (one solid only)")),
        _ => Err(unsupported(
            "this union is outside the supported union domains (a boss standing on one planar face of the body; two right prisms along one axis)",
        )),
    }
}

/// Op 24: the union of two closed solids (body, then tool), one solid.
pub fn union(body: &Body, tool: &Body, tol: f64) -> Result<Body, KernelError> {
    match select(body, tool, tol)? {
        Domain::Face(plan) => unite_face::unite(body, tool, &plan, tol).map(|(b, _)| b),
        Domain::Stack => unite_stack::union(body, tool, tol),
        Domain::Overhang(plan) => unite_overhang::unite(body, tool, &plan, tol),
    }
}

/// Op 25: the union, the recorded history of every output entity (provenance.rs roles 30-39: body = the first operand,
/// tool = the second) and the outcome of every input entity. A domain that does not record its history refuses by name.
pub fn union_named(body: &Body, tool: &Body, tol: f64) -> Result<(Body, Provenance, Vec<Outcome>), KernelError> {
    let (result, history) = match select(body, tool, tol)? {
        Domain::Face(plan) => unite_face::unite(body, tool, &plan, tol)?,
        Domain::Stack => return Err(unsupported("union history is not recorded for the prism stack union domain yet")),
        Domain::Overhang(_) => return Err(unsupported("union history is not recorded for overhanging bosses yet")),
    };
    let counts = (history.faces.len(), history.edges.len(), history.vertices.len());
    if counts != (result.faces.len(), result.edges.len(), result.vertices.len()) {
        return Err(KernelError::Internal(format!("union history names {counts:?} entities of another count")));
    }
    let inputs = [body.faces.len(), body.edges.len(), body.vertices.len(), tool.faces.len(), tool.edges.len(), tool.vertices.len()];
    let outcomes = outcomes(&history, inputs);
    Ok((result, history, outcomes))
}

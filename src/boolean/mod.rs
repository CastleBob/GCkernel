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
mod trace;
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
        // A clear box passage through a body with no curved faces splices its cavity (box_cut.rs, which records its
        // history) instead of rebuilding cells (Claude 2026-11-02).
        0 if box_cut::clear_box_passage(body, tool, tolerance) => Ok(Domain::BoxCut),
        0 => Ok(Domain::Rectilinear),
        1 if drill(tool) => Ok(Domain::Bore),
        1 => Ok(Domain::Cross),
        2 => Ok(Domain::Capsule),
        _ => Err(crate::KernelError::Unsupported(
            "this curved cutter is outside the supported cut domains".into(),
        )),
    }
}

thread_local! {
    static EQUIVALENCES: std::cell::RefCell<Vec<String>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Records a within-tolerance equivalence a domain used for the current cut (Bob 10-28: "the status says so").
pub(crate) fn note(text: String) {
    EQUIVALENCES.with(|e| e.borrow_mut().push(text));
}

/// The within-tolerance equivalences the cuts run on this thread have used since the last call, in order, and clears
/// them: the caller's status reports them.
pub fn take_equivalences() -> Vec<String> {
    EQUIVALENCES.with(|e| std::mem::take(&mut *e.borrow_mut()))
}

pub fn cut(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
) -> Result<crate::neutral::Body, crate::KernelError> {
    run(select(body, tool, tolerance)?, body, tool, tolerance)
}

/// The cut in its selected domain, with no history.
fn run(
    domain: Domain,
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
) -> Result<crate::neutral::Body, crate::KernelError> {
    match domain {
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

/// A domain's recorded history of one cut, or why it records none. The body never depends on which.
pub(crate) type History = std::result::Result<crate::provenance::Provenance, String>;

/// A cut's records and input outcomes, or why there are none.
pub type Recorded = std::result::Result<(crate::provenance::Provenance, Vec<crate::provenance::Outcome>), String>;

/// The domains that record history (for some or all of their cuts).
fn records(domain: Domain) -> bool {
    matches!(domain, Domain::Bore | Domain::Section | Domain::Radial | Domain::ConeBore | Domain::Capsule | Domain::BoxCut)
}

/// The cut in its one selected domain and that domain's history. `strict` (op 12) refuses a domain that records nothing
/// before cutting. A domain recording history computes the body once, the same body as op 3, whether or not it can
/// record this cut's history.
fn traced(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
    strict: bool,
) -> Result<(crate::neutral::Body, History), crate::KernelError> {
    let domain = select(body, tool, tolerance)?;
    let none = format!("cut history is not recorded for the {domain:?} cut domain yet");
    if strict && !records(domain) {
        return Err(crate::KernelError::Unsupported(none));
    }
    Ok(match domain {
        Domain::Bore => bore::cut_traced(body, tool, tolerance)?,
        Domain::Section => section::cut_traced(body, tool, tolerance)?,
        Domain::Radial => radial::cut_traced(body, tool, tolerance)?,
        Domain::ConeBore => cone_bore::cut_traced(body, tool, tolerance)?,
        Domain::Capsule => capsule::cut_traced(body, tool, tolerance)?,
        Domain::BoxCut => box_cut::cut_traced(body, tool, tolerance)?,
        domain => (run(domain, body, tool, tolerance)?, Err(none)),
    })
}

/// The admitted history's checks and the input outcomes. A count mismatch is a defect; two generated entities with one
/// name make the history unrecordable.
fn finish(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    result: &crate::neutral::Body,
    history: History,
) -> Result<Recorded, crate::KernelError> {
    use crate::provenance::*;
    let history = match history {
        Ok(history) => history,
        Err(reason) => return Ok(Err(reason)),
    };
    let counts = (history.faces.len(), history.edges.len(), history.vertices.len());
    if counts != (result.faces.len(), result.edges.len(), result.vertices.len()) {
        return Err(crate::KernelError::Internal(format!("cut history names {counts:?} entities of another count")));
    }
    // Generated entities must be told apart by their records alone.
    for list in [&history.faces, &history.edges, &history.vertices] {
        for (i, r) in list.iter().enumerate() {
            if matches!(r.role, MEETING | TOOL_SEAM | SEAM_POINT | MEETING_END | TOOL_EDGE_END | MEETING_POINT) && list[i + 1..].contains(r) {
                return Ok(Err(format!("cut history would name two entities {r:?}")));
            }
        }
    }
    let inputs = [body.faces.len(), body.edges.len(), body.vertices.len(), tool.faces.len(), tool.edges.len(), tool.vertices.len()];
    let outcomes = outcomes(&history, inputs);
    Ok(Ok((history, outcomes)))
}

/// Op 12: the cut, the recorded history of its every output entity (provenance.rs roles 30-41) and the outcome of every
/// input entity. A domain that does not record its history refuses by name; nothing is matched after the fact.
pub fn cut_named(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
) -> Result<(crate::neutral::Body, crate::provenance::Provenance, Vec<crate::provenance::Outcome>), crate::KernelError> {
    let (result, history) = traced(body, tool, tolerance, true)?;
    let (history, outcomes) = finish(body, tool, &result, history)?.map_err(crate::KernelError::Unsupported)?;
    Ok((result, history, outcomes))
}

/// Op 26 (Claude 2026-11-02, persistent names on every stock): op 3's cut, run once in its one selected domain, with
/// op 12's history when that domain records this cut's history, otherwise the reason it does not. The body is op 3's
/// either way; a failed cut fails. Nothing is retried, matched or guessed: the host keeps names only while every cut
/// in its history carries records.
pub fn cut_recorded(
    body: &crate::neutral::Body,
    tool: &crate::neutral::Body,
    tolerance: f64,
) -> Result<(crate::neutral::Body, Recorded), crate::KernelError> {
    let (result, history) = traced(body, tool, tolerance, false)?;
    let recorded = finish(body, tool, &result, history)?;
    Ok((result, recorded))
}

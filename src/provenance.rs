//! Operation provenance (Claude 2026-10-30, contract with Codex in GCm GCmaker/Handoff/Claude-GcKernel-PersistentReferences.md):
//! for every output face, edge and vertex, the role it plays and the request-local indices it was made from. The kernel
//! never names anything; the host maps its own persistent identities to request indices and back.
//!
//! Envelope appended to the named ops' replies: start marker, schema version, the count of values that follow up to the end
//! marker, the output counts, then one record [role, k, refs...] per face, edge and vertex in body order, then the end marker.

pub const START: f64 = 9731.0;
pub const END: f64 = -9731.0;
pub const VERSION: f64 = 1.0;

// Op 10 (named extrude). Loop and segment indices are the request's, before the extruder's reordering and reversal.
pub const LOW_CAP: u32 = 1;
pub const HIGH_CAP: u32 = 2;
/// [loop, seg]
pub const SIDE: u32 = 3;
/// [loop, seg]
pub const LOW_SEGMENT: u32 = 4;
/// [loop, seg]
pub const HIGH_SEGMENT: u32 = 5;
/// The ruling through segment seg's start point: [loop, seg].
pub const RULING: u32 = 6;
/// Segment seg's start point at the low cap: [loop, seg].
pub const LOW_POINT: u32 = 7;
/// Segment seg's start point at the high cap: [loop, seg].
pub const HIGH_POINT: u32 = 8;

// Op 11 (named revolve, op 9's request). Request loop/segment indices, before the reordering, reversal and negative-turn
// mirroring. Angles: start = 0 (the seam direction x), end = the requested signed angle (part turns only), middle = the
// split of faces touching the axis (pi in a full turn, half the angle in a part turn of more than pi). part: 0 unsplit,
// 1 from the start to the middle, 2 from the middle to the end (or back to the start in a full turn).
/// The surface swept by segment seg: [loop, seg, part].
pub const SWEPT: u32 = 10;
pub const START_CAP: u32 = 11;
pub const END_CAP: u32 = 12;
/// The circle or arc swept by segment seg's start point: [loop, seg, part].
pub const PARALLEL: u32 = 13;
/// Segment seg at the start angle (in a full turn, the seam): [loop, seg].
pub const START_MERIDIAN: u32 = 14;
/// Segment seg at the end angle: [loop, seg].
pub const END_MERIDIAN: u32 = 15;
/// Segment seg at the middle angle: [loop, seg].
pub const MID_MERIDIAN: u32 = 16;
/// Segment seg's start point, on the axis: [loop, seg].
pub const AXIS_POINT: u32 = 17;
/// Segment seg's start point at the start angle: [loop, seg].
pub const START_POINT: u32 = 18;
/// Segment seg's start point at the end angle: [loop, seg].
pub const END_POINT: u32 = 19;
/// Segment seg's start point at the middle angle: [loop, seg].
pub const MID_POINT: u32 = 20;
/// Segment seg lying on the axis, shared by the two caps of a part turn: [loop, seg].
pub const AXIS_SEGMENT: u32 = 21;

// Op 12 (named cut, op 3's request: body then tool). Refs are the request bodies' own entity indices. trimmed: 0 when the
// output has its source's support and the same bounding edges (by provenance), 1 otherwise.
/// From body face i: [i, trimmed].
pub const FROM_BODY_FACE: u32 = 30;
/// From tool face j: [j, trimmed].
pub const FROM_TOOL_FACE: u32 = 31;
/// From body edge i: [i, trimmed].
pub const FROM_BODY_EDGE: u32 = 32;
/// The n-th edge generated where body face i meets tool face j (n = 0, 1, ... in the domain's recorded order): [i, j, n].
pub const MEETING: u32 = 33;
/// From tool edge j: [j, trimmed].
pub const FROM_TOOL_EDGE: u32 = 34;
/// Body vertex i: [i].
pub const FROM_BODY_VERTEX: u32 = 35;
/// The seam line a domain places on the face it makes from tool face j, starting where body face i meets it: [i, j].
pub const TOOL_SEAM: u32 = 37;
/// The seam point on the meeting edge of body face i and tool face j: [i, j].
pub const SEAM_POINT: u32 = 38;
/// Tool vertex j: [j].
pub const FROM_TOOL_VERTEX: u32 = 39;
/// The end on body face c of meeting edge [i, j, n]: [i, j, n, c].
pub const MEETING_END: u32 = 40;
/// The end on body face c of the edge from tool edge j: [j, c].
pub const TOOL_EDGE_END: u32 = 41;

// Op 12 outcome codes, one per input entity.
pub const UNCHANGED: u32 = 1;
pub const MODIFIED: u32 = 2;
pub const SPLIT: u32 = 3;
pub const DELETED: u32 = 4;
/// A tool entity contributing nothing to the result's boundary.
pub const CONSUMED: u32 = 5;

/// What became of one input entity: its code and the output entities (of its own kind) descending from it.
#[derive(Clone, Debug, PartialEq)]
pub struct Outcome {
    pub code: u32,
    pub outputs: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub role: u32,
    pub refs: Vec<usize>,
}

impl Record {
    pub fn new(role: u32, refs: &[usize]) -> Self { Self { role, refs: refs.to_vec() } }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Provenance {
    pub faces: Vec<Record>,
    pub edges: Vec<Record>,
    pub vertices: Vec<Record>,
}

/// Each input entity's outcome (body faces, edges, vertices, then the tool's; `inputs` = their counts), from the outputs
/// whose records descend from it (roles 30, 31, 32, 34, 35, 39): none = deleted (body) or consumed (tool), one = unchanged
/// or modified by its trimmed flag (vertices: unchanged), more = split.
pub fn outcomes(p: &Provenance, inputs: [usize; 6]) -> Vec<Outcome> {
    let lists = [
        (&p.faces, FROM_BODY_FACE, DELETED),
        (&p.edges, FROM_BODY_EDGE, DELETED),
        (&p.vertices, FROM_BODY_VERTEX, DELETED),
        (&p.faces, FROM_TOOL_FACE, CONSUMED),
        (&p.edges, FROM_TOOL_EDGE, CONSUMED),
        (&p.vertices, FROM_TOOL_VERTEX, CONSUMED),
    ];
    let mut out = Vec::new();
    for ((records, role, none), &count) in lists.into_iter().zip(&inputs) {
        for i in 0..count {
            let outputs: Vec<usize> = (0..records.len()).filter(|&o| records[o].role == role && records[o].refs[0] == i).collect();
            let code = match outputs.len() {
                0 => none,
                1 if records[outputs[0]].refs.get(1).map_or(true, |&t| t == 0) => UNCHANGED,
                1 => MODIFIED,
                _ => SPLIT,
            };
            out.push(Outcome { code, outputs });
        }
    }
    out
}

/// The envelope, appended to `out`.
pub fn write(out: &mut Vec<f64>, p: &Provenance) {
    envelope(out, records(p));
}

/// Op 12's envelope: the records, then the input counts and one outcome per input entity in that order.
pub fn write_cut(out: &mut Vec<f64>, p: &Provenance, inputs: [usize; 6], outcomes: &[Outcome]) {
    let mut body = records(p);
    body.extend(inputs.iter().map(|&n| n as f64));
    for o in outcomes {
        body.push(o.code as f64);
        body.push(o.outputs.len() as f64);
        body.extend(o.outputs.iter().map(|&i| i as f64));
    }
    envelope(out, body);
}

fn records(p: &Provenance) -> Vec<f64> {
    let mut body = vec![p.faces.len() as f64, p.edges.len() as f64, p.vertices.len() as f64];
    for r in p.faces.iter().chain(&p.edges).chain(&p.vertices) {
        body.push(r.role as f64);
        body.push(r.refs.len() as f64);
        body.extend(r.refs.iter().map(|&i| i as f64));
    }
    body
}

fn envelope(out: &mut Vec<f64>, body: Vec<f64>) {
    out.push(START);
    out.push(VERSION);
    out.push(body.len() as f64);
    out.extend(body);
    out.push(END);
}

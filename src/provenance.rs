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

/// The envelope, appended to `out`.
pub fn write(out: &mut Vec<f64>, p: &Provenance) {
    let mut body = vec![p.faces.len() as f64, p.edges.len() as f64, p.vertices.len() as f64];
    for r in p.faces.iter().chain(&p.edges).chain(&p.vertices) {
        body.push(r.role as f64);
        body.push(r.refs.len() as f64);
        body.extend(r.refs.iter().map(|&i| i as f64));
    }
    out.push(START);
    out.push(VERSION);
    out.push(body.len() as f64);
    out.extend(body);
    out.push(END);
}

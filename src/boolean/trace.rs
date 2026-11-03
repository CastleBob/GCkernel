//! Op 12 history shared by the domains that keep the body and attach a cavity (Claude 2026-11-02, persistent names on
//! every stock; contract with Codex in GCm GCmaker/Handoff/Claude-GcKernel-PersistentReferences.md). Everything here is
//! read from the construction itself: which tool segment each cavity segment was built from, which body face each
//! cavity cap joined. Nothing is matched by position after the fact; where a construction fact is missing the history
//! is refused with the reason (the cut itself is unaffected).
use crate::neutral::Body;
use crate::provenance::*;
use std::collections::HashMap;

pub(super) type Refusal = String;

fn refuse<T>(m: &str) -> Result<T, Refusal> {
    Err(format!("cut history is not recorded for this cut: {m}"))
}

/// The face other than `not` holding edge `e`, when exactly one does.
pub(super) fn other_face(body: &Body, e: usize, not: usize) -> Result<usize, Refusal> {
    let faces: Vec<usize> = (0..body.faces.len())
        .filter(|&f| f != not && body.faces[f].loops.iter().flatten().any(|(x, _)| *x == e))
        .collect();
    match faces[..] {
        [f] => Ok(f),
        _ => refuse("a prism edge has no single side face"),
    }
}

/// The edge at vertex `v` that is not in face `cap`: a right prism's ruling there, when exactly one.
pub(super) fn ruling_at(body: &Body, cap: usize, v: usize) -> Result<usize, Refusal> {
    let in_cap: Vec<usize> = body.faces[cap].loops.iter().flatten().map(|(e, _)| *e).collect();
    let r: Vec<usize> = (0..body.edges.len())
        .filter(|e| !in_cap.contains(e) && (body.edges[*e].start == v || body.edges[*e].end == v))
        .collect();
    match r[..] {
        [e] => Ok(e),
        _ => refuse("a prism vertex has no single ruling"),
    }
}

/// The one edge face `side` shares with face `cap`.
pub(super) fn shared_edge(body: &Body, side: usize, cap: usize) -> Result<usize, Refusal> {
    let in_cap: Vec<usize> = body.faces[cap].loops.iter().flatten().map(|(e, _)| *e).collect();
    let r: Vec<usize> = body.faces[side].loops.iter().flatten().map(|(e, _)| *e).filter(|e| in_cap.contains(e)).collect();
    match r[..] {
        [e] => Ok(e),
        _ => refuse("a prism side face has no single cap edge"),
    }
}

/// What a prism operand's cap loop segment was built from: its side face and the ruling through its start.
#[derive(Clone, Copy, Debug)]
pub(super) struct Segment {
    pub side: usize,
    pub ruling: usize,
}

/// The segments of loop 0 of `solid`'s cap face `cap`, in that loop's order, each oriented as the loop runs it.
pub(super) fn cap_segments(solid: &Body, cap: usize) -> Result<Vec<Segment>, Refusal> {
    let Some(lp) = solid.faces[cap].loops.first() else { return refuse("a prism cap has no loop") };
    lp.iter()
        .map(|&(e, forward)| {
            let start = if forward { solid.edges[e].start } else { solid.edges[e].end };
            Ok(Segment { side: other_face(solid, e, cap)?, ruling: ruling_at(solid, cap, start)? })
        })
        .collect()
}

/// Numbers the meetings of one body face and one tool face, n = 0, 1, ... in the order they are recorded.
#[derive(Default)]
pub(super) struct Meetings(HashMap<(usize, usize), usize>);

impl Meetings {
    pub(super) fn next(&mut self, face: usize, tool_face: usize) -> usize {
        let n = self.0.entry((face, tool_face)).or_insert(0);
        *n += 1;
        *n - 1
    }
}

/// The history of a body kept whole: every face, edge and vertex from itself, unchanged.
pub(super) fn kept(body: &Body) -> Provenance {
    Provenance {
        faces: (0..body.faces.len()).map(|i| Record::new(FROM_BODY_FACE, &[i, 0])).collect(),
        edges: (0..body.edges.len()).map(|i| Record::new(FROM_BODY_EDGE, &[i, 0])).collect(),
        vertices: (0..body.vertices.len()).map(|i| Record::new(FROM_BODY_VERTEX, &[i])).collect(),
    }
}

/// passage::splice_cavity's records: the cavity (op 10's extrusion of one loop whose segment s was built from tool
/// segment `tool[s]`) joins its low cap's loop to body face `bottom` and its high cap's to `top`, both then modified, and
/// appends its vertices, edges and side faces in its own order. Side s is the tool's side face; its low and high copies
/// the meetings of the two body faces with that side; the ruling the tool's ruling cut to the wall; its ends that ruling's
/// ends on the two body faces.
pub(super) fn splice(
    history: &mut Provenance,
    cavity: &Provenance,
    tool: &[Segment],
    bottom: usize,
    top: usize,
    meetings: &mut Meetings,
) -> Result<(), Refusal> {
    let segment = |r: &Record| -> Result<Segment, Refusal> {
        match (r.refs.first(), r.refs.get(1)) {
            (Some(0), Some(&s)) if s < tool.len() => Ok(tool[s]),
            _ => refuse("a cavity segment has no tool segment"),
        }
    };
    history.faces[bottom] = Record::new(FROM_BODY_FACE, &[bottom, 1]);
    history.faces[top] = Record::new(FROM_BODY_FACE, &[top, 1]);
    for r in &cavity.faces {
        match r.role {
            LOW_CAP | HIGH_CAP => {}
            SIDE => history.faces.push(Record::new(FROM_TOOL_FACE, &[segment(r)?.side, 1])),
            _ => return refuse("an unexpected cavity face"),
        }
    }
    for r in &cavity.edges {
        let s = segment(r)?;
        history.edges.push(match r.role {
            LOW_SEGMENT => Record::new(MEETING, &[bottom, s.side, meetings.next(bottom, s.side)]),
            HIGH_SEGMENT => Record::new(MEETING, &[top, s.side, meetings.next(top, s.side)]),
            RULING => Record::new(FROM_TOOL_EDGE, &[s.ruling, 1]),
            _ => return refuse("an unexpected cavity edge"),
        });
    }
    for r in &cavity.vertices {
        let s = segment(r)?;
        history.vertices.push(match r.role {
            LOW_POINT => Record::new(TOOL_EDGE_END, &[s.ruling, bottom]),
            HIGH_POINT => Record::new(TOOL_EDGE_END, &[s.ruling, top]),
            _ => return refuse("an unexpected cavity vertex"),
        });
    }
    Ok(())
}

/// The records of a round cutter bored through curved walls (radial.rs, cone_bore.rs): the body is kept, each wall the
/// cutter crosses gains one closed loop (one vertex, one edge: the meeting of that wall and the cutter's cylinder `j`, n
/// counting that wall's crossings along the cutter), and each pair of crossings one seam (keyed by the pair's entry wall)
/// and one cylinder face from the cutter's. `crossings` are (wall face) in construction order, the vertex and edge of
/// each appended in that order, then per pair its seam edge and its face.
pub(super) fn bored(body: &Body, crossings: &[usize], j: usize) -> Result<Provenance, Refusal> {
    if crossings.len() % 2 != 0 {
        return refuse("a bore crossing without its pair");
    }
    let mut history = kept(body);
    let mut meetings = Meetings::default();
    for &face in crossings {
        let n = meetings.next(face, j);
        history.faces[face] = Record::new(FROM_BODY_FACE, &[face, 1]);
        history.vertices.push(Record::new(MEETING_POINT, &[face, j, n]));
        history.edges.push(Record::new(MEETING, &[face, j, n]));
    }
    for pair in crossings.chunks_exact(2) {
        history.edges.push(Record::new(TOOL_SEAM, &[pair[0], j]));
        history.faces.push(Record::new(FROM_TOOL_FACE, &[j, 1]));
    }
    Ok(history)
}

/// The tool segments of a box cutter's passage outline built from its bounds (box_cut.rs): segment i runs from
/// `corners[i]` to `corners[i + 1]` across the cutter's cap `cap` (corners given in the cap's two in-plane coordinates
/// u, v; each corner IS one of the cutter's cap vertices, its bounds being made of them). Side i is the cutter face
/// holding the cap edge joining those two vertices, its ruling the one at corner i.
pub(super) fn box_segments(tool: &Body, cap: usize, u: usize, v: usize, corners: &[[f64; 2]]) -> Result<Vec<Segment>, Refusal> {
    let at = |q: crate::geom::Vec3| { let c = [q.x, q.y, q.z]; [c[u], c[v]] };
    let in_cap: Vec<usize> = tool.faces[cap].loops.iter().flatten().map(|(e, _)| *e).collect();
    let mut cap_vertices: Vec<usize> = in_cap.iter().flat_map(|&e| [tool.edges[e].start, tool.edges[e].end]).collect();
    cap_vertices.sort();
    cap_vertices.dedup();
    let vertex = |c: [f64; 2]| -> Result<usize, Refusal> {
        let found: Vec<usize> = cap_vertices.iter().copied().filter(|&i| at(tool.vertices[i]) == c).collect();
        match found[..] { [i] => Ok(i), _ => refuse("a cutter corner is not one cap vertex") }
    };
    (0..corners.len())
        .map(|i| {
            let (a, b) = (vertex(corners[i])?, vertex(corners[(i + 1) % corners.len()])?);
            let edges: Vec<usize> = in_cap.iter().copied()
                .filter(|&e| { let x = &tool.edges[e]; (x.start, x.end) == (a, b) || (x.start, x.end) == (b, a) })
                .collect();
            let [edge] = edges[..] else { return refuse("two cutter corners are not joined by one cap edge") };
            Ok(Segment { side: other_face(tool, edge, cap)?, ruling: ruling_at(tool, cap, a)? })
        })
        .collect()
}

/// The records of a plain box rebuilt by op 10's extrusion with one hole (capsule.rs on a box): profile loop 0 is the
/// body's low cap loop run backwards (segment s is that loop's segment L - 1 - s reversed), loop 1 the tool's cap loop
/// `tool_cap` in its order. The box's caps (`low`, `high`) are modified; its sides, their cap edges, rulings and corners
/// are whole and unchanged; the hole's walls, meetings, rulings and their ends come from the tool.
pub(super) fn box_with_hole(
    body: &Body,
    low: usize,
    high: usize,
    tool: &Body,
    tool_cap: usize,
    records: &Provenance,
) -> Result<Provenance, Refusal> {
    let Some(lp) = body.faces[low].loops.first() else { return refuse("the box cap has no loop") };
    let n = lp.len();
    let start = |t: usize| { let (e, f) = lp[t % n]; if f { body.edges[e].start } else { body.edges[e].end } };
    let hole = cap_segments(tool, tool_cap)?;
    let mut meetings = Meetings::default();
    // Body segment of profile segment s: (its cap edge, side face, the corner it starts at, the ruling there).
    let outer = |s: usize| -> Result<(usize, usize, usize, usize), Refusal> {
        if s >= n {
            return refuse("a rebuilt box segment has no box segment");
        }
        let t = n - 1 - s;
        let edge = lp[t].0;
        let corner = start(t + 1);
        Ok((edge, other_face(body, edge, low)?, corner, ruling_at(body, low, corner)?))
    };
    let tool_seg = |s: usize| hole.get(s).copied().ok_or_else(|| "cut history is not recorded for this cut: a hole segment has no tool segment".to_string());
    let mut out = Provenance::default();
    for r in &records.faces {
        out.faces.push(match (r.role, r.refs.first()) {
            (LOW_CAP, _) => Record::new(FROM_BODY_FACE, &[low, 1]),
            (HIGH_CAP, _) => Record::new(FROM_BODY_FACE, &[high, 1]),
            (SIDE, Some(0)) => Record::new(FROM_BODY_FACE, &[outer(r.refs[1])?.1, 0]),
            (SIDE, Some(1)) => Record::new(FROM_TOOL_FACE, &[tool_seg(r.refs[1])?.side, 1]),
            _ => return refuse("an unexpected rebuilt face"),
        });
    }
    for r in &records.edges {
        out.edges.push(match (r.role, r.refs[0]) {
            (LOW_SEGMENT, 0) => Record::new(FROM_BODY_EDGE, &[outer(r.refs[1])?.0, 0]),
            (HIGH_SEGMENT, 0) => Record::new(FROM_BODY_EDGE, &[shared_edge(body, outer(r.refs[1])?.1, high)?, 0]),
            (RULING, 0) => Record::new(FROM_BODY_EDGE, &[outer(r.refs[1])?.3, 0]),
            (LOW_SEGMENT, 1) => { let t = tool_seg(r.refs[1])?; Record::new(MEETING, &[low, t.side, meetings.next(low, t.side)]) }
            (HIGH_SEGMENT, 1) => { let t = tool_seg(r.refs[1])?; Record::new(MEETING, &[high, t.side, meetings.next(high, t.side)]) }
            (RULING, 1) => Record::new(FROM_TOOL_EDGE, &[tool_seg(r.refs[1])?.ruling, 1]),
            _ => return refuse("an unexpected rebuilt edge"),
        });
    }
    for r in &records.vertices {
        out.vertices.push(match (r.role, r.refs[0]) {
            (LOW_POINT, 0) => Record::new(FROM_BODY_VERTEX, &[outer(r.refs[1])?.2]),
            (HIGH_POINT, 0) => {
                let (_, _, corner, ruling) = outer(r.refs[1])?;
                let e = &body.edges[ruling];
                Record::new(FROM_BODY_VERTEX, &[if e.start == corner { e.end } else { e.start }])
            }
            (LOW_POINT, 1) => Record::new(TOOL_EDGE_END, &[tool_seg(r.refs[1])?.ruling, low]),
            (HIGH_POINT, 1) => Record::new(TOOL_EDGE_END, &[tool_seg(r.refs[1])?.ruling, high]),
            _ => return refuse("an unexpected rebuilt vertex"),
        });
    }
    Ok(out)
}

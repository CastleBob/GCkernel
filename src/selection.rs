//! Read-only selection regions over a validated B-rep. Raw topology is never rewritten.
//! Shared boundaries between identical oriented analytic supports join face regions.
//! Degree-two boundaries on one line/circle join only with the same two face regions.
use crate::{geom::Curve, neutral::Body, KernelError};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;
const ROUND: f64 = 64. * f64::EPSILON;

#[derive(Clone, Debug)]
pub struct FaceRegion {
    pub faces: Vec<usize>,
}
#[derive(Clone, Debug)]
pub struct EdgeRegion {
    /// Ordered source edges with traversal direction. No approximate reconstructed polyline.
    pub edges: Vec<(usize, bool)>,
    pub curve: Curve,
    pub start: usize,
    pub end: usize,
    pub closed: bool,
}
#[derive(Clone, Debug)]
pub struct SelectionTopology {
    pub faces: Vec<FaceRegion>,
    pub edges: Vec<EdgeRegion>,
    pub face_region: Vec<usize>,
    /// None marks an internal same-support division, not a selectable boundary.
    pub edge_region: Vec<Option<usize>>,
}
fn bad(message: impl Into<String>) -> KernelError {
    KernelError::InvalidInput(message.into())
}
fn root(parent: &[usize], mut i: usize) -> usize {
    while parent[i] != i {
        i = parent[i];
    }
    i
}
fn same_curve(a: &Curve, b: &Curve) -> bool {
    match (a, b) {
        (Curve::Line { a: p, b: q }, Curve::Line { a: r, b: s }) => {
            let u = q.sub(*p);
            let v = s.sub(*r);
            let delta = r.sub(*p);
            let (lu, lv) = (u.length(), v.length());
            lu > 0.
                && lv > 0.
                && u.cross(v).length() <= ROUND * lu * lv
                && delta.cross(u).length() <= ROUND * lu * delta.length().max(1.)
        }
        (
            Curve::Arc { c, n, r, .. },
            Curve::Arc {
                c: d, n: m, r: s, ..
            },
        ) => {
            c.distance(*d) <= ROUND * r.abs().max(s.abs()).max(1.)
                && n.cross(*m).length() <= ROUND
                && (r - s).abs() <= ROUND * r.abs().max(s.abs()).max(1.)
        }
        _ => false,
    }
}

/// Validates first, then constructs a detached selection catalog. This is neither
/// solid healing nor persistent topological naming. IDs belong to this body revision.
pub fn build(body: &Body, tolerance: f64) -> Result<SelectionTopology, KernelError> {
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(bad(
            "Selection topology requires a positive finite tolerance.",
        ));
    }
    // The reader checks indices and frames before topology::check indexes them.
    let data = crate::neutral::write_body(body);
    let mut reader = crate::neutral::Reader::with_tolerance(&data, tolerance);
    reader.body()?;
    reader.finish()?;
    crate::topology::check(body, tolerance)?;
    let mut uses = vec![Vec::new(); body.edges.len()];
    for (f, face) in body.faces.iter().enumerate() {
        for &(e, _) in face.loops.iter().flatten() {
            uses[e].push(f);
        }
    }
    let mut parents: Vec<_> = (0..body.faces.len()).collect();
    for adjacent in &uses {
        let [a, b] = adjacent.as_slice() else {
            return Err(bad(
                "A selection boundary needs exactly two oriented face uses.",
            ));
        };
        let (a, b) = (root(&parents, *a), root(&parents, *b));
        if a == b {
            continue;
        }
        let (first, second) = (a.min(b), a.max(b));
        // Compare to one representative so a sequence of roundoff-near supports cannot drift.
        if crate::tessellate::same_support(&body.faces[first], &body.faces[second])
            && (0..body.faces.len())
                .filter(|&f| root(&parents, f) == second)
                .all(|f| crate::tessellate::same_support(&body.faces[first], &body.faces[f]))
        {
            parents[second] = first;
        }
    }
    let mut face_groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for f in 0..body.faces.len() {
        face_groups.entry(root(&parents, f)).or_default().push(f);
    }
    let faces: Vec<_> = face_groups
        .into_values()
        .map(|faces| FaceRegion { faces })
        .collect();
    let mut face_region = vec![0; body.faces.len()];
    for (id, region) in faces.iter().enumerate() {
        for &f in &region.faces {
            face_region[f] = id;
        }
    }
    let mut adjacency = vec![Vec::new(); body.vertices.len()];
    let mut boundary = Vec::new();
    let mut supports = vec![(0, 0); body.edges.len()];
    for (id, adjacent) in uses.iter().enumerate() {
        let (a, b) = (face_region[adjacent[0]], face_region[adjacent[1]]);
        if a == b {
            continue;
        }
        supports[id] = (a.min(b), a.max(b));
        boundary.push(id);
        adjacency[body.edges[id].start].push(id);
        adjacency[body.edges[id].end].push(id);
    }
    let mut links = vec![[None; 2]; body.edges.len()];
    for (vertex, edges) in adjacency.iter().enumerate() {
        if edges.len() != 2 || edges[0] == edges[1] {
            continue;
        }
        let (a, b) = (edges[0], edges[1]);
        if supports[a] == supports[b] && same_curve(&body.edges[a].curve, &body.edges[b].curve) {
            let end_a = usize::from(body.edges[a].end == vertex);
            let end_b = usize::from(body.edges[b].end == vertex);
            links[a][end_a] = Some(b);
            links[b][end_b] = Some(a);
        }
    }
    let mut unvisited: BTreeSet<_> = boundary.into_iter().collect();
    let mut edges = Vec::new();
    let mut edge_region = vec![None; body.edges.len()];
    while let Some(&seed) = unvisited.first() {
        let mut component = BTreeSet::new();
        let mut queue = vec![seed];
        while let Some(e) = queue.pop() {
            if component.insert(e) {
                queue.extend(links[e].iter().flatten().copied());
            }
        }
        let start_edge = component
            .iter()
            .copied()
            .find(|&e| links[e].iter().any(Option::is_none))
            .unwrap_or(seed);
        let first = &body.edges[start_edge];
        let start = if links[start_edge][0].is_none() {
            first.start
        } else {
            first.end
        };
        let (mut at, mut edge) = (start, start_edge);
        let mut ordered = Vec::new();
        loop {
            if !unvisited.remove(&edge) {
                return Err(bad("Selection boundary traversal repeats before closing."));
            }
            let raw = &body.edges[edge];
            let forward = raw.start == at;
            ordered.push((edge, forward));
            at = if forward { raw.end } else { raw.start };
            let next = links[edge][usize::from(forward)];
            match next {
                None => break,
                Some(next) if next == start_edge => break,
                Some(next) => edge = next,
            }
        }
        if ordered.len() != component.len() {
            return Err(bad("Selection boundary is not one uninterrupted chain."));
        }
        let closed = at == start;
        let curve = chain_curve(body, &ordered, closed)?;
        let id = edges.len();
        for &(e, _) in &ordered {
            edge_region[e] = Some(id);
        }
        edges.push(EdgeRegion {
            edges: ordered,
            curve,
            start,
            end: at,
            closed,
        });
    }
    Ok(SelectionTopology {
        faces,
        edges,
        face_region,
        edge_region,
    })
}

fn oriented(body: &Body, edge: (usize, bool)) -> Curve {
    let curve = &body.edges[edge.0].curve;
    if edge.1 {
        curve.clone()
    } else {
        curve.reversed()
    }
}
fn chain_curve(body: &Body, edges: &[(usize, bool)], closed: bool) -> Result<Curve, KernelError> {
    let first = oriented(body, edges[0]);
    if edges.len() == 1 {
        return Ok(first);
    }
    match &first {
        Curve::Line { a, b } => {
            if closed {
                return Err(bad(
                    "A straight selection boundary cannot form a closed loop.",
                ));
            }
            let direction = b.sub(*a).unit()?;
            let mut end = *b;
            for &edge in &edges[1..] {
                let next = oriented(body, edge);
                let Curve::Line { a: p, b: q } = next else {
                    return Err(bad("Mixed curve kinds in a straight selection boundary."));
                };
                let length = q.distance(p);
                let scale = length.max(end.distance(*a)).max(1.);
                if p.distance(end) > ROUND * scale
                    || q.sub(p).dot(direction) <= 0.
                    || q.sub(*a).cross(direction).length() > ROUND * scale
                {
                    return Err(bad("Straight selection boundary does not have one continuous analytic support."));
                }
                end = q;
            }
            Ok(Curve::Line { a: *a, b: end })
        }
        Curve::Arc { c, n, x, r, t0, t1 } => {
            let mut sweep = t1 - t0;
            for &edge in &edges[1..] {
                let next = oriented(body, edge);
                let Curve::Arc {
                    n: m, t0: s, t1: t, ..
                } = &next
                else {
                    return Err(bad("Mixed curve kinds in a circular selection boundary."));
                };
                let combined = Curve::Arc {
                    c: *c,
                    n: *n,
                    x: *x,
                    r: *r,
                    t0: *t0,
                    t1: *t0 + sweep,
                };
                if !same_curve(&first, &next)
                    || n.sub(*m).length() > ROUND
                    || combined.end().distance(next.start())
                        > ROUND * r.max(1.) * edges.len() as f64
                {
                    return Err(bad("Circular selection boundary does not have one continuous analytic support."));
                }
                sweep += t - s;
            }
            if sweep > TAU + ROUND * edges.len() as f64 {
                return Err(bad("Circular selection boundary overlaps itself."));
            }
            if closed && (sweep - TAU).abs() > ROUND * edges.len() as f64 {
                return Err(bad(
                    "Closed circular selection boundary does not cover one turn.",
                ));
            }
            Ok(Curve::Arc {
                c: *c,
                n: *n,
                x: *x,
                r: *r,
                t0: *t0,
                t1: *t0 + if closed { TAU } else { sweep },
            })
        }
        _ => Err(KernelError::Unsupported(
            "Only line and circle boundary fragments have a unified measurement curve.".into(),
        )),
    }
}

#[cfg(test)]
mod tests;

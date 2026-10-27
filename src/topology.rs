//! Checks that a body is a closed, consistently oriented solid. Occurrence based (Astra's review): a seam
//! edge appears in one face's loop twice, once each way, and every occurrence is counted.

use crate::geom::Curve;
use crate::neutral::Body;
use crate::KernelError;

/// The counts behind the checks; Euler-Poincare is a diagnostic, not a requirement (holes, genus).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TopologyReport {
    pub vertices: usize,
    pub edges: usize,
    pub faces: usize,
    pub loops: usize,
}

impl TopologyReport {
    /// V - E + F - (L - F) = 2 (S - G) for one shell: the genus it implies.
    pub fn genus(&self) -> f64 {
        let chi = self.vertices as f64 - self.edges as f64 + self.faces as f64 - (self.loops as f64 - self.faces as f64);
        1.0 - chi / 2.0
    }
}

fn bad<T>(message: String) -> Result<T, KernelError> { Err(KernelError::InvalidInput(message)) }

/// Every edge used exactly once forward and once reversed; every loop connects head to tail and closes;
/// every edge's curve runs from its start vertex to its end vertex (within `tolerance`).
pub fn check(body: &Body, tolerance: f64) -> Result<TopologyReport, KernelError> {
    if body.faces.is_empty() { return bad("the body has no faces".into()); }
    for (i, e) in body.edges.iter().enumerate() {
        let (s, t) = (body.vertices[e.start], body.vertices[e.end]);
        if e.curve.start().distance(s) > tolerance || e.curve.end().distance(t) > tolerance {
            return bad(format!("edge {i}'s curve does not run from its start vertex to its end vertex"));
        }
        if let Curve::Line { .. } = e.curve {
            if e.start == e.end { return bad(format!("edge {i} is a line that starts and ends at one vertex")); }
        }
    }
    let mut forward = vec![0usize; body.edges.len()];
    let mut backward = vec![0usize; body.edges.len()];
    let mut loops = 0;
    for (f, face) in body.faces.iter().enumerate() {
        if face.loops.is_empty() { return bad(format!("face {f} has no loop")); }
        for (l, lp) in face.loops.iter().enumerate() {
            loops += 1;
            if lp.is_empty() { return bad(format!("face {f} loop {l} is empty")); }
            let ends = |&(edge, fwd): &(usize, bool)| {
                let e = &body.edges[edge];
                if fwd { (e.start, e.end) } else { (e.end, e.start) }
            };
            for k in 0..lp.len() {
                let (_, tail) = ends(&lp[k]);
                let (head, _) = ends(&lp[(k + 1) % lp.len()]);
                if tail != head { return bad(format!("face {f} loop {l} breaks after segment {k}")); }
                let (edge, fwd) = lp[k];
                if fwd { forward[edge] += 1 } else { backward[edge] += 1 }
            }
        }
    }
    for i in 0..body.edges.len() {
        if forward[i] != 1 || backward[i] != 1 {
            return bad(format!("edge {i} is used {} times forward and {} times reversed (a closed solid uses each once each way)", forward[i], backward[i]));
        }
    }
    Ok(TopologyReport { vertices: body.vertices.len(), edges: body.edges.len(), faces: body.faces.len(), loops })
}

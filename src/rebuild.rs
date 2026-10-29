//! The kernel's own base solid (Claude 2026-10-27): a body that is a straight prism is rebuilt by the kernel's
//! extruder from one cap's loops, so later cuts start from kernel geometry rather than from the host CAD's body.
//!
//! A prism here, proved from the topology and geometry, not from aggregate measures (Astra): two plane caps with
//! opposite outward normals d and -d, at distance L; every edge of the top cap bounds exactly one side face, a plane
//! parallel to d or a cylinder along d, whose one loop is that edge, a straight ruling of length L along d, the same
//! curve moved by -L d, and a ruling back; every side face is reached that way exactly once; and the bottom cap's
//! edges are exactly those moved curves. Equal volume and area of the rebuilt body (mass::measure, 1e-9 relative)
//! are then a regression check on top. Anything else is Unsupported, never a guess.

use crate::extrude::extrude;
use crate::geom::{Curve, Surface, Vec3};
use crate::mass::measure;
use crate::neutral::{Body, Profile};
use crate::KernelError;

const PARALLEL: f64 = 1e-9;
const AGREE: f64 = 1e-9;

fn unsupported<T>(m: &str) -> Result<T, KernelError> { Err(KernelError::Unsupported(format!("rebuild: {m}"))) }

/// The rebuilt prism. The input is checked first (topology and every edge on its faces), whoever calls this:
/// a malformed body is reported, never rebuilt into a sound one (Astra).
pub fn rebuild_prism(body: &Body, tolerance: f64) -> Result<Body, KernelError> {
    crate::topology::check(body, tolerance)?;
    let outward = |i: usize| match body.faces[i].surface {
        Surface::Plane { n, .. } => Some(if body.faces[i].sense { n } else { n.scale(-1.0) }),
        _ => None,
    };
    // Every pair of plane faces with opposite normals is a candidate; the sides decide.
    for top in 0..body.faces.len() {
        let Some(d) = outward(top) else { continue };
        let Some(bottom) = (0..body.faces.len()).find(|&j| outward(j).is_some_and(|m| m.add(d).length() < PARALLEL)) else {
            continue;
        };
        let sides_ok = body.faces.iter().enumerate().filter(|(i, _)| *i != top && *i != bottom).all(|(_, f)| match f.surface {
            Surface::Plane { n, .. } => n.dot(d).abs() < PARALLEL,
            Surface::Cylinder { a, .. } | Surface::EllipticCylinder { a, .. } => a.cross(d).length() < PARALLEL,
            Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } => false,
        });
        if !sides_ok { continue; }
        let caps = [top, bottom].map(|i| match body.faces[i].surface { Surface::Plane { p, .. } => p, _ => unreachable!() });
        let length = caps[0].sub(caps[1]).dot(d);
        if length <= tolerance { continue; }
        prove_prism(body, top, bottom, d, length, tolerance)?;
        let loops = body.faces[top].loops.iter()
            .map(|lp| lp.iter().map(|&(e, fwd)| {
                let c = &body.edges[e].curve;
                if fwd { c.clone() } else { c.reversed() }
            }).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let rebuilt = extrude(&Profile { loops, n: d.scale(-1.0), low: 0.0, high: length }, tolerance)?;
        let (before, after) = (measure(body, tolerance)?, measure(&rebuilt, tolerance)?);
        let agree = |a: f64, b: f64| (a - b).abs() <= AGREE * a.abs().max(b.abs()).max(1e-300);
        if !agree(before.volume, after.volume) || !agree(before.area, after.area) {
            return unsupported(&format!(
                "the body is not the straight prism of its cap (volume {} vs {}, area {} vs {})",
                before.volume, after.volume, before.area, after.area));
        }
        return Ok(rebuilt);
    }
    unsupported("the body is not a straight prism")
}

/// Same geometric curve, in either direction: ends and middles within `tol`.
fn same_curve(a: &Curve, b: &Curve, tol: f64) -> bool {
    let mid = |c: &Curve| { let (t0, t1) = c.domain(); c.eval((t0 + t1) / 2.0).0 };
    let ends = (a.start().distance(b.start()) <= tol && a.end().distance(b.end()) <= tol)
        || (a.start().distance(b.end()) <= tol && a.end().distance(b.start()) <= tol);
    ends && mid(a).distance(mid(b)) <= tol
}

fn ruling(curve: &Curve, d: Vec3, length: f64, tol: f64) -> bool {
    matches!(curve, Curve::Line { a, b } if (b.sub(*a).length() - length).abs() <= tol && b.sub(*a).cross(d).length() <= tol)
}

/// The structural proof described at the top of the file.
fn prove_prism(body: &Body, top: usize, bottom: usize, d: Vec3, length: f64, tol: f64) -> Result<(), KernelError> {
    let shift = d.scale(-length);
    let mut side_used = vec![false; body.faces.len()];
    let mut moved_edges = Vec::new();
    for &(e, _) in body.faces[top].loops.iter().flatten() {
        let owners: Vec<usize> = (0..body.faces.len())
            .filter(|&f| f != top && body.faces[f].loops.iter().flatten().any(|&(g, _)| g == e))
            .collect();
        let [side] = owners[..] else { return unsupported(&format!("top edge {e} bounds {} side faces", owners.len())) };
        if side == bottom || side_used[side] { return unsupported(&format!("top edge {e} does not have a side face of its own")); }
        side_used[side] = true;
        let face = &body.faces[side];
        let [lp] = &face.loops[..] else { return unsupported(&format!("side face {side} has {} loops", face.loops.len())) };
        if lp.len() != 4 { return unsupported(&format!("side face {side} has {} edges, not 4", lp.len())); }
        let at = lp.iter().position(|&(g, _)| g == e).unwrap();
        let [r1, moved, r2] = [1, 2, 3].map(|k| lp[(at + k) % 4].0);
        if !ruling(&body.edges[r1].curve, d, length, tol) || !ruling(&body.edges[r2].curve, d, length, tol) {
            return unsupported(&format!("side face {side} is not bounded by rulings of length {length}"));
        }
        if !same_curve(&body.edges[moved].curve, &body.edges[e].curve.translated(shift), tol) {
            return unsupported(&format!("side face {side}'s bottom edge is not its top edge moved along the sweep"));
        }
        moved_edges.push(moved);
    }
    if let Some(f) = (0..body.faces.len()).find(|&f| f != top && f != bottom && !side_used[f]) {
        return unsupported(&format!("face {f} is not the sweep of a top cap edge"));
    }
    let mut bottom_edges: Vec<usize> = body.faces[bottom].loops.iter().flatten().map(|&(g, _)| g).collect();
    bottom_edges.sort_unstable();
    moved_edges.sort_unstable();
    if bottom_edges != moved_edges { return unsupported("the bottom cap is not the top cap moved along the sweep"); }
    Ok(())
}

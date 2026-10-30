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
            Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } => false,
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

/// Op 13 (Claude 2026-10-30, Bob: nest a Flat part as the kernel builds it): the profile whose straight sweep along
/// `-d` is exactly `body`, its loops on the top (outer counter-clockwise about `d`, holes clockwise).
///
/// A cut may leave a cap in several coplanar pieces (the rectilinear cells), so a cap here is every plane face
/// facing `d` (or `-d`); all of them must lie in one plane. Edges between two pieces of a cap are seams, not
/// outline. The proof is `prove_prism`'s on the outline: each top outline edge bounds exactly one side face (a plane
/// parallel to `d` or a cylinder along `d`), whose one loop is that edge, a ruling of the thickness, the edge moved
/// down, and a ruling back; every side face is reached once; the bottom outline is exactly the moved edges. The
/// outline chains into loops (two loops touching at a point are refused), runs of collinear lines join, and the
/// profile's own sweep must agree with the body in volume and area. Everything else is refused by name.
pub fn prism_profile(body: &Body, d: Vec3, tolerance: f64) -> Result<Profile, KernelError> {
    crate::topology::check(body, tolerance)?;
    let facing = |m: Vec3| (0..body.faces.len()).filter(|&i| outward(body, i).is_some_and(|n| n.sub(m).length() < PARALLEL)).collect::<Vec<_>>();
    let (tops, bottoms) = (facing(d), facing(d.scale(-1.0)));
    if tops.is_empty() || bottoms.is_empty() {
        return unsupported(&format!("{} faces face the top and {} the bottom; a plate needs both", tops.len(), bottoms.len()));
    }
    let height = |f: usize| match body.faces[f].surface { Surface::Plane { p, .. } => p.dot(d), _ => unreachable!() };
    for cap in [&tops, &bottoms] {
        let (lo, hi) = cap.iter().map(|&f| height(f)).fold((f64::MAX, f64::MIN), |(lo, hi), h| (lo.min(h), hi.max(h)));
        if hi - lo > tolerance {
            return unsupported(&format!(
                "faces facing the {} lie {} apart; a straight-through plate has one top and one bottom (a pocket or a \
                 step does not)", if std::ptr::eq(cap, &tops) { "top" } else { "bottom" }, hi - lo));
        }
    }
    let length = height(tops[0]) - height(bottoms[0]);
    if length <= tolerance { return unsupported(&format!("the plate is {length} thick")); }
    let is_cap = |f: usize| tops.contains(&f) || bottoms.contains(&f);
    if let Some(f) = (0..body.faces.len()).filter(|&f| !is_cap(f)).find(|&f| !swept(&body.faces[f].surface, d)) {
        return unsupported(&format!("face {f} ({}) does not run straight through the plate", kind(&body.faces[f].surface)));
    }
    let top_outline = outline(body, &tops)?;
    let bottom_outline = outline(body, &bottoms)?;

    let shift = d.scale(-length);
    let mut side_used = vec![false; body.faces.len()];
    let mut moved_edges = Vec::new();
    for &(e, _) in &top_outline {
        let owners: Vec<usize> = (0..body.faces.len())
            .filter(|&f| !is_cap(f) && body.faces[f].loops.iter().flatten().any(|&(g, _)| g == e))
            .collect();
        let [side] = owners[..] else { return unsupported(&format!("top edge {e} bounds {} side faces", owners.len())) };
        if side_used[side] { return unsupported(&format!("top edge {e} does not have a side face of its own")); }
        side_used[side] = true;
        let face = &body.faces[side];
        let [lp] = &face.loops[..] else { return unsupported(&format!("side face {side} has {} loops", face.loops.len())) };
        if lp.len() != 4 { return unsupported(&format!("side face {side} has {} edges, not 4", lp.len())); }
        let at = lp.iter().position(|&(g, _)| g == e).unwrap();
        let [r1, moved, r2] = [1, 2, 3].map(|k| lp[(at + k) % 4].0);
        if !ruling(&body.edges[r1].curve, d, length, tolerance) || !ruling(&body.edges[r2].curve, d, length, tolerance) {
            return unsupported(&format!("side face {side} is not bounded by rulings of the plate's thickness {length}"));
        }
        if !same_curve(&body.edges[moved].curve, &body.edges[e].curve.translated(shift), tolerance) {
            return unsupported(&format!("side face {side}'s bottom edge is not its top edge moved through the plate"));
        }
        moved_edges.push(moved);
    }
    if let Some(f) = (0..body.faces.len()).find(|&f| !is_cap(f) && !side_used[f]) {
        return unsupported(&format!("face {f} is not the sweep of a top outline edge"));
    }
    let mut bottom_edges: Vec<usize> = bottom_outline.iter().map(|&(g, _)| g).collect();
    bottom_edges.sort_unstable();
    moved_edges.sort_unstable();
    if bottom_edges != moved_edges { return unsupported("the bottom outline is not the top outline moved through the plate"); }

    let loops = chain(body, &top_outline, tolerance)?;
    let profile = Profile { loops, n: d.scale(-1.0), low: 0.0, high: length };
    let rebuilt = extrude(&profile, tolerance)?;
    let (before, after) = (measure(body, tolerance)?, measure(&rebuilt, tolerance)?);
    let agree = |a: f64, b: f64| (a - b).abs() <= AGREE * a.abs().max(b.abs()).max(1e-300);
    if !agree(before.volume, after.volume) || !agree(before.area, after.area) {
        return unsupported(&format!(
            "the body is not the straight sweep of its top (volume {} vs {}, area {} vs {})",
            before.volume, after.volume, before.area, after.area));
    }
    Ok(profile)
}

fn swept(s: &Surface, d: Vec3) -> bool {
    match *s {
        Surface::Plane { n, .. } => n.dot(d).abs() < PARALLEL,
        Surface::Cylinder { a, .. } | Surface::EllipticCylinder { a, .. } => a.cross(d).length() < PARALLEL,
        Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } => false,
    }
}

/// A cap's outline: the uses (edge, forward) of edges that bound exactly one of its pieces. An edge used by two
/// pieces is a seam between them; any other count is refused.
fn outline(body: &Body, cap: &[usize]) -> Result<Vec<(usize, bool)>, KernelError> {
    let mut uses: Vec<Vec<bool>> = vec![Vec::new(); body.edges.len()];
    for &f in cap { for &(e, fwd) in body.faces[f].loops.iter().flatten() { uses[e].push(fwd); } }
    let mut out = Vec::new();
    for (e, u) in uses.iter().enumerate() {
        match u[..] {
            [] => {}
            [fwd] => out.push((e, fwd)),
            [a, b] if a != b => {}
            _ => return unsupported(&format!("edge {e} is used {} times by one cap", u.len())),
        }
    }
    Ok(out)
}

/// The outline's oriented curves chained end to start into closed loops, collinear runs of lines joined.
fn chain(body: &Body, outline: &[(usize, bool)], tol: f64) -> Result<Vec<Vec<Curve>>, KernelError> {
    let ends = |&(e, fwd): &(usize, bool)| { let g = &body.edges[e]; if fwd { (g.start, g.end) } else { (g.end, g.start) } };
    let mut from: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for (k, u) in outline.iter().enumerate() {
        if from.insert(ends(u).0, k).is_some() {
            return unsupported(&format!("two outline loops touch at vertex {}", ends(u).0));
        }
    }
    let mut done = vec![false; outline.len()];
    let mut loops = Vec::new();
    for first in 0..outline.len() {
        if done[first] { continue; }
        let mut lp = Vec::new();
        let mut k = first;
        while !done[k] {
            done[k] = true;
            let (e, fwd) = outline[k];
            let c = &body.edges[e].curve;
            lp.push(if fwd { c.clone() } else { c.reversed() });
            let Some(&next) = from.get(&ends(&outline[k]).1) else {
                return unsupported(&format!("the outline is open at vertex {}", ends(&outline[k]).1));
            };
            k = next;
        }
        if k != first { return unsupported("the outline does not close into loops"); }
        loops.push(join_lines(lp, tol));
    }
    Ok(loops)
}

/// Consecutive lines along one straight line become one (a cap cell's seams cut the outline's edges into pieces).
fn join_lines(lp: Vec<Curve>, tol: f64) -> Vec<Curve> {
    let along = |p: &Curve, q: &Curve| match (p, q) {
        (Curve::Line { a, b }, Curve::Line { b: c, .. }) => {
            let u = b.sub(*a);
            let len = u.length();
            len > 0.0 && u.scale(1.0 / len).cross(c.sub(*a)).length() <= tol && c.sub(*a).dot(u) > len
        }
        _ => false,
    };
    let mut out: Vec<Curve> = Vec::with_capacity(lp.len());
    for c in lp {
        if let Some(last) = out.last_mut() {
            if along(last, &c) {
                let (Curve::Line { a, .. }, Curve::Line { b, .. }) = (last.clone(), &c) else { unreachable!() };
                *last = Curve::Line { a, b: *b };
                continue;
            }
        }
        out.push(c);
    }
    // The loop's first and last lines may be one run too.
    if out.len() > 2 && along(&out[out.len() - 1], &out[0]) {
        let (Curve::Line { a, .. }, Curve::Line { b, .. }) = (out.pop().unwrap(), &out[0]) else { unreachable!() };
        out[0] = Curve::Line { a, b: *b };
    }
    out
}

fn outward(body: &Body, i: usize) -> Option<Vec3> {
    match body.faces[i].surface {
        Surface::Plane { n, .. } => Some(if body.faces[i].sense { n } else { n.scale(-1.0) }),
        _ => None,
    }
}

fn kind(s: &Surface) -> &'static str {
    match s {
        Surface::Plane { .. } => "plane",
        Surface::Cylinder { .. } => "cylinder",
        Surface::EllipticCylinder { .. } => "elliptic cylinder",
        Surface::Cone { .. } => "cone",
        Surface::Sphere { .. } => "sphere",
        Surface::Torus { .. } => "torus",
        Surface::Revolution { .. } => "surface of revolution",
    }
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

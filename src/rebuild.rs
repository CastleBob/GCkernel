//! The kernel's own base solid (Claude 2026-10-27): a body that is a straight prism is rebuilt by the kernel's
//! extruder from one cap's loops, so later cuts start from kernel geometry rather than from the host CAD's body.
//!
//! A prism here: two plane caps with opposite outward normals d and -d, each one face, and every other face a plane
//! parallel to d or a cylinder along d. The rebuilt body is accepted only when its volume and area (mass::measure)
//! agree with the input's to 1e-9 relative; anything else is Unsupported, never a guess.

use crate::extrude::extrude;
use crate::geom::Surface;
use crate::mass::measure;
use crate::neutral::{Body, Profile};
use crate::KernelError;

const PARALLEL: f64 = 1e-9;
const AGREE: f64 = 1e-9;

fn unsupported<T>(m: &str) -> Result<T, KernelError> { Err(KernelError::Unsupported(format!("rebuild: {m}"))) }

/// The rebuilt prism.
pub fn rebuild_prism(body: &Body, tolerance: f64) -> Result<Body, KernelError> {
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
            Surface::Cylinder { a, .. } => a.cross(d).length() < PARALLEL,
        });
        if !sides_ok { continue; }
        let caps = [top, bottom].map(|i| match body.faces[i].surface { Surface::Plane { p, .. } => p, _ => unreachable!() });
        let length = caps[0].sub(caps[1]).dot(d);
        if length <= tolerance { continue; }
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

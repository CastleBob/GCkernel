//! Normalized pipe ends (Claude 2026-10-31, Bob's Normalize End; GCmaker's EndConditionSideCompositeBuilder and
//! PipeCopeNormalizer): a radial (vertical) torch cuts one boundary per angle round the pipe, so the normalized end keeps,
//! at each angle, whichever wall's crossing with the mating cylinder lies deeper in the material, and the end face follows
//! that height radially across the wall. Its exact face is the radial conoid of the deeper wall's CylinderPair
//! (`Surface::Conoid`). This module builds the waste tool for one round mating cylinder (op 23); boolean::cut removes it
//! from a plain pipe (boolean/normal_end.rs).
//!
//! Domain, refused by name otherwise: the mating cylinder crosses both walls all the way round (the pair's discriminant
//! proved positive), and one wall is proved at least as deep as the other at every angle, within the tolerance (what
//! agrees within the tolerance is the same: a square cope's two crossings meet at the pipe's sides). A cope whose deeper
//! wall changes round the pipe has a two-piece envelope, and is refused.
//!
//! Proofs (independent, elementary): with D(t) = h_in(t) - h_out(t), each height's second derivative is bounded
//! (geom::pair_bend_bound), so on each of n steps D stays above its chord's lesser end less max |D''| step^2 / 8.
use crate::geom::{pair_bend_bound, pair_disc_coefficients, pair_disc_lower_bound, pair_height, pair_height_range, Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;
use std::f64::consts::TAU;

type Result<T> = std::result::Result<T, KernelError>;

fn unsupported<T>(message: &str) -> Result<T> {
    Err(KernelError::Unsupported(message.into()))
}

/// The pipe and its mating cylinder, in one frame: pipe axis point `o`, unit `d` from the cut end into the material,
/// unit `x` square to it (where the seams go), the inner and outer wall radii; the mating cylinder's axis point `p`,
/// unit axis `b` and radius `rc`.
#[derive(Clone, Copy, Debug)]
pub struct Envelope {
    pub o: Vec3,
    pub d: Vec3,
    pub x: Vec3,
    pub ri: f64,
    pub ro: f64,
    pub p: Vec3,
    pub b: Vec3,
    pub rc: f64,
}

const STEPS: usize = 16384;

/// The wall whose crossing is the normalized envelope (the deeper one, toward d), proved at every angle.
pub fn deeper_wall(e: &Envelope, tol: f64) -> Result<f64> {
    if !(e.ri > 0.0 && e.ro > e.ri + tol && e.rc > 0.0) {
        return Err(KernelError::InvalidInput("a normalized end needs 0 < inner radius < outer radius and a positive mating radius".into()));
    }
    if (e.d.length() - 1.0).abs() > 1e-12 || (e.x.length() - 1.0).abs() > 1e-12 || e.d.dot(e.x).abs() > 1e-12
        || (e.b.length() - 1.0).abs() > 1e-12 {
        return Err(KernelError::InvalidInput("a normalized end needs unit axes and a seam direction square to the pipe".into()));
    }
    if e.b.cross(e.d).length() < 1e-9 {
        return unsupported("Normalize End: the mating cylinder runs along the pipe");
    }
    let mut floors = [0.0; 2];
    let mut bends = [0.0; 2];
    for (i, rho) in [e.ri, e.ro].into_iter().enumerate() {
        let k = pair_disc_coefficients(e.o, e.d, e.x, rho, e.p, e.b, e.rc);
        floors[i] = pair_disc_lower_bound(&k, 0.0, TAU, STEPS);
        if !(floors[i] > tol * tol) {
            return unsupported("Normalize End: the mating cylinder does not cross both pipe walls all the way round (part of the end would stay uncut)");
        }
        bends[i] = pair_bend_bound(e.d, rho, e.b, &k, floors[i]);
    }
    let step = TAU / STEPS as f64;
    let pad = (bends[0] + bends[1]) * step * step / 8.0;
    let gap = |t: f64| pair_height(e.o, e.d, e.x, e.ri, e.p, e.b, e.rc, 1.0, t).0 - pair_height(e.o, e.d, e.x, e.ro, e.p, e.b, e.rc, 1.0, t).0;
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..=STEPS {
        let g = gap(step * i as f64);
        lo = lo.min(g);
        hi = hi.max(g);
    }
    if lo - pad >= -tol {
        Ok(e.ri)
    } else if hi + pad <= tol {
        Ok(e.ro)
    } else {
        unsupported("Normalize End: the deeper wall changes round the pipe (a two-piece envelope); not supported natively yet")
    }
}

/// The waste band GCmaker's normalized cut removes, exact: radii sa..sb about the pipe axis, from level `far` along d
/// (beyond the end, outside the material) up to the conoid of the deeper wall. Returns it and its volume, the integral
/// (sb^2 - sa^2) / 2 times that of (h(t) - far) round the turn (by quadrature).
pub fn waste_tool(e: &Envelope, sa: f64, sb: f64, far: f64, tol: f64) -> Result<(Body, f64)> {
    let rho = deeper_wall(e, tol)?;
    if !(sa > tol && sa < e.ri - tol && sb > e.ro + tol) {
        return Err(KernelError::InvalidInput("a normalized end's waste band must reach past both walls (0 < inner < pipe bore, outer > pipe)".into()));
    }
    let (o, d, x) = (e.o, e.d, e.x);
    let pair = |s: f64| Curve::RadialPair { c: o, a: d, x, rho, p: e.p, b: e.b, r: e.rc, branch: 1.0, t0: 0.0, t1: TAU, s };
    let (hmin, _) = pair_height_range(o, d, x, rho, e.p, e.b, e.rc, 1.0, 0.0, TAU, STEPS)
        .ok_or_else(|| KernelError::Unsupported("Normalize End: the envelope height is not proved defined all the way round".into()))?;
    if !(far < hmin - tol) {
        return Err(KernelError::InvalidInput("a normalized end's waste must start beyond the envelope".into()));
    }
    let h0 = pair_height(o, d, x, rho, e.p, e.b, e.rc, 1.0, 0.0).0;
    let vertices = vec![
        o.add(x.scale(sa)).add(d.scale(h0)),
        o.add(x.scale(sb)).add(d.scale(h0)),
        o.add(x.scale(sa)).add(d.scale(far)),
        o.add(x.scale(sb)).add(d.scale(far)),
    ];
    let ring = |r: f64| Curve::Arc { c: o.add(d.scale(far)), n: d, x, r, t0: 0.0, t1: TAU };
    let edges = vec![
        Edge { curve: pair(sa), start: 0, end: 0 },
        Edge { curve: pair(sb), start: 1, end: 1 },
        Edge { curve: Curve::Line { a: vertices[0], b: vertices[1] }, start: 0, end: 1 },
        Edge { curve: ring(sa), start: 2, end: 2 },
        Edge { curve: ring(sb), start: 3, end: 3 },
        Edge { curve: Curve::Line { a: vertices[2], b: vertices[0] }, start: 2, end: 0 },
        Edge { curve: Curve::Line { a: vertices[3], b: vertices[1] }, start: 3, end: 1 },
    ];
    let conoid = Surface::Conoid { c: o, a: d, x, rho, p: e.p, b: e.b, r: e.rc, branch: 1.0 };
    let faces = vec![
        // The top: outward along d, against the conoid's natural normal (-s d + h' v).
        Face { surface: conoid, sense: false, loops: vec![vec![(2, true), (1, true), (2, false), (0, false)]] },
        // The far cap: outward -d; the outer ring counter-clockwise about -d, the inner one a hole.
        Face { surface: Surface::Plane { p: o.add(d.scale(far)), n: d.scale(-1.0), x }, sense: true, loops: vec![vec![(4, false)], vec![(3, true)]] },
        // The inner cylinder faces the axis; the outer one away from it.
        Face { surface: Surface::Cylinder { p: o, a: d, x, r: sa }, sense: false, loops: vec![vec![(5, true), (0, true), (5, false), (3, false)]] },
        Face { surface: Surface::Cylinder { p: o, a: d, x, r: sb }, sense: true, loops: vec![vec![(4, true), (6, true), (1, false), (6, false)]] },
    ];
    let body = Body { vertices, edges, faces };
    crate::topology::check(&body, tol)?;
    let size = sb.max(far.abs()).max(1e-9);
    let integral = crate::mass::quadrature(&|t| pair_height(o, d, x, rho, e.p, e.b, e.rc, 1.0, t).0 - far, 0.0, TAU, 1e-14 * size)
        .map_err(KernelError::Unsupported)?;
    Ok((body, (sb * sb - sa * sa) / 2.0 * integral))
}

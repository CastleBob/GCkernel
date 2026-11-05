//! A plain pipe end cut by a normalized-end waste band (Claude 2026-10-31, Bob's Normalize End; the tool comes from
//! normalize::waste_tool, worker op 23).
//!
//! Domain, decided from the geometry before anything is built: the tool is a band about an axis, radii sa < sb, bounded
//! by a radial conoid (its top), a plane square to the axis (its far cap) and two coaxial cylinders; the body is a plain
//! pipe (outer and inner wall, an annular cap at each end) on the same axis line, its walls strictly inside the band
//! (sa < ri, ro < sb, by more than the tolerance); the conoid's height round the whole turn is proved strictly between
//! the pipe's two end planes, and the tool's far cap lies beyond one end, outside the pipe. Then pipe less tool is,
//! exactly, the pipe from the conoid to its other end: each wall's new end is the conoid's rim on it (a RadialPair at
//! the wall's radius, one closed loop) and the new end face is the conoid between the rims. The far cap, both wall
//! surfaces and the far rings are kept exactly; the seams lie on the far ring's generator. Anything else is refused by
//! name; nothing is retried.
//! A second opposite normalized end is also admitted when the first complete conoid
//! and the new conoid have strictly separated axial ranges. Its entire far cap,
//! rims and radial seam are retained. This does not recut an existing conoid.
use crate::geom::{pair_height_range, Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;
use std::f64::consts::TAU;

type Result<T> = std::result::Result<T, KernelError>;

fn unsupported(message: &str) -> KernelError {
    KernelError::Unsupported(format!("normalized pipe end: {message}"))
}

fn perp(v: Vec3, b: Vec3) -> Vec3 {
    v.sub(b.scale(v.dot(b)))
}

/// The band tool: its conoid face's surface, inner and outer radius, and its far cap's level along the conoid's axis.
struct Band {
    surface: Surface,
    sa: f64,
    sb: f64,
    far: f64,
}

#[derive(Clone, Copy)]
struct End {
    face: usize,
    low: f64,
    high: f64,
    outer: usize,
    inner: usize,
    normalized: bool,
}

fn band(tool: &Body, tol: f64) -> Result<Band> {
    if tool.faces.len() != 4 {
        return Err(unsupported("the tool is not a normalized-end waste band (four faces)"));
    }
    let conoids: Vec<&Face> = tool.faces.iter().filter(|f| matches!(f.surface, Surface::Conoid { .. })).collect();
    let [top] = conoids[..] else { return Err(unsupported("the tool has more than one conoid face")) };
    let Surface::Conoid { c, a, .. } = top.surface else { unreachable!() };
    let mut radii = vec![];
    let mut far = None;
    for f in &tool.faces {
        match f.surface {
            Surface::Conoid { .. } => {}
            Surface::Cylinder { p, a: ca, r, .. } => {
                if ca.cross(a).length() > 1e-9 || perp(p.sub(c), a).length() > tol {
                    return Err(unsupported("the waste band's walls are not coaxial with its conoid"));
                }
                radii.push(r);
            }
            Surface::Plane { p, n, .. } => {
                if n.cross(a).length() > 1e-9 || far.is_some() {
                    return Err(unsupported("the waste band's far cap is not one plane square to its axis"));
                }
                far = Some(p.sub(c).dot(a));
            }
            _ => return Err(unsupported("the waste band has a face that is not its conoid, a wall or its cap")),
        }
    }
    let (Some(far), &[r0, r1]) = (far, &radii[..]) else { return Err(unsupported("the waste band needs two walls and a far cap")) };
    let (sa, sb) = (r0.min(r1), r0.max(r1));
    if !(sb > sa + tol) {
        return Err(unsupported("the waste band's walls coincide"));
    }
    Ok(Band { surface: top.surface.clone(), sa, sb, far })
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    crate::topology::check(body, tol)?;
    crate::topology::check(tool, tol)?;
    super::radial::planar_faces_outward(body)?;
    super::radial::planar_faces_outward(tool)?;
    let band = band(tool, tol)?;
    let Surface::Conoid { c, a, rho, p: cp, b: cb, r: cr, branch, .. } = band.surface else { unreachable!() };
    if body.faces.len() != 4 {
        return Err(unsupported("needs a plain pipe (two walls and an annular cap at each end)"));
    }
    let mut outer = None;
    let mut inner = None;
    for (i, f) in body.faces.iter().enumerate() {
        if let Surface::Cylinder { p, a: wa, r, .. } = f.surface {
            let slot = if f.sense { &mut outer } else { &mut inner };
            if slot.is_some() {
                return Err(unsupported("needs a plain pipe (two walls and an annular cap at each end)"));
            }
            *slot = Some((i, p, wa, r));
        }
    }
    let (Some((oi, p, b, ro)), Some((ii, pi, bi, ri))) = (outer, inner) else {
        return Err(unsupported("needs a plain pipe (two walls and an annular cap at each end)"));
    };
    if bi.cross(b).length() > 1e-9 || perp(pi.sub(p), b).length() > tol || !(ri < ro - tol) || ri <= tol {
        return Err(unsupported("the pipe's walls are not coaxial with the inner one inside"));
    }
    if b.cross(a).length() > 1e-9 || perp(p.sub(c), a).length() > tol {
        return Err(unsupported("the waste band is not on the pipe's axis"));
    }
    if !(band.sa < ri - tol && band.sb > ro + tol) {
        return Err(unsupported("the waste band does not reach past both walls"));
    }
    // Each wall: one band of two complete rims and a twice-used axial seam.
    let wall = |face: usize, r: f64| -> Result<[usize; 2]> {
        let f = &body.faces[face];
        let [lp] = &f.loops[..] else { return Err(unsupported("a pipe wall already has holes or cuts")) };
        let rings: Vec<usize> = lp.iter().filter(|(e, _)| matches!(body.edges[*e].curve, Curve::Arc { .. } | Curve::RadialPair { .. })).map(|(e, _)| *e).collect();
        let lines: Vec<usize> = lp.iter().filter(|(e, _)| matches!(body.edges[*e].curve, Curve::Line { .. })).map(|(e, _)| *e).collect();
        if lp.len() != 4 || rings.len() != 2 || lines.len() != 2 || lines[0] != lines[1] || rings[0] == rings[1] {
            return Err(unsupported("a pipe wall is not one plain band"));
        }
        for e in &rings {
            let (rc, n, ar, t0, t1) = match body.edges[*e].curve {
                Curve::Arc { c, n, r, t0, t1, .. } => (c, n, r, t0, t1),
                Curve::RadialPair { c, a, s, t0, t1, .. } => (c, a, s, t0, t1),
                _ => unreachable!(),
            };
            if (t1 - t0 - TAU).abs() > 1e-9 || n.cross(a).length() > 1e-9 || (ar - r).abs() > tol || perp(rc.sub(c), a).length() > tol
                || body.edges[*e].start != body.edges[*e].end {
                return Err(unsupported("a pipe wall's end is not a full ring"));
            }
        }
        Ok([rings[0], rings[1]])
    };
    let (wo, wi) = (wall(oi, ro)?, wall(ii, ri)?);
    // Each end supplies a containing axial range. A normalized far end must
    // have exactly its two RadialPair rims and a twice-used radial seam.
    let mut ends = vec![];
    for cf in (0..4).filter(|i| *i != oi && *i != ii) {
        let f = &body.faces[cf];
        let ring_of = |w: &[usize; 2]| w.iter().copied().find(|e| f.loops.iter().flatten().any(|(fe, _)| fe == e));
        let (Some(o_e), Some(i_e)) = (ring_of(&wo), ring_of(&wi)) else {
            return Err(unsupported("a pipe end does not close its walls"));
        };
        let (low, high, normalized) = match f.surface {
            Surface::Plane { p: q, n, .. } => {
                if n.cross(a).length() > 1e-9 || f.loops.len() != 2 || f.loops.iter().any(|l| l.len() != 1)
                    || ![o_e, i_e].iter().all(|e| matches!(body.edges[*e].curve, Curve::Arc { .. })) {
                    return Err(unsupported("a pipe end is not a plain annulus square to the pipe"));
                }
                let h = q.sub(c).dot(a);
                (h, h, false)
            }
            Surface::Conoid { c: fc, a: fa, x: fx, rho: frho, p: fp, b: fb, r: fr, branch: fbranch } => {
                // This new domain carries exact parallel frames; it never snaps
                // a tilted retained conoid onto the new end's axis.
                if fa != a && fa != a.scale(-1.0) {
                    return Err(unsupported("the retained normalized end does not have the same axis direction"));
                }
                if perp(fc.sub(c), a).length() > tol || f.loops.len() != 1 || f.loops[0].len() != 4 {
                    return Err(unsupported("the retained normalized end is not one complete conoid"));
                }
                for (edge, radius) in [(o_e, ro), (i_e, ri)] {
                    let Curve::RadialPair { c: ec, a: ea, x: ex, rho: erho, p: ep, b: eb, r: er, branch: ebranch, s: es, .. } = body.edges[edge].curve else {
                        return Err(unsupported("the retained normalized end has a non-radial rim"));
                    };
                    if ec != fc || ea != fa || ex != fx || erho != frho || ep != fp || eb != fb || er != fr || ebranch != fbranch || (es-radius).abs() > tol {
                        return Err(unsupported("the retained conoid and its rim recipes differ"));
                    }
                }
                let seam: Vec<_> = f.loops[0].iter().filter(|(e, _)| *e != o_e && *e != i_e).collect();
                if seam.len() != 2 || seam[0].0 != seam[1].0 || seam[0].1 == seam[1].1
                    || !matches!(body.edges[seam[0].0].curve, Curve::Line { .. }) {
                    return Err(unsupported("the retained conoid lacks one twice-used radial seam"));
                }
                let (h0, h1) = pair_height_range(fc, fa, fx, frho, fp, fb, fr, fbranch, 0., TAU, 16384)
                    .ok_or_else(|| unsupported("the retained conoid's height is not proved defined all the way round"))?;
                let sign = fa.dot(a);
                let origin = fc.sub(c).dot(a);
                let (q0,q1) = (origin+sign*h0, origin+sign*h1);
                (q0.min(q1), q0.max(q1), true)
            }
            _ => return Err(unsupported("a pipe end is neither planar nor one complete radial conoid")),
        };
        ends.push(End {face:cf, low, high, outer:o_e, inner:i_e, normalized});
    }
    if ends.iter().filter(|e| e.normalized).count() > 1 {
        return Err(unsupported("both pipe ends are already normalized"));
    }
    ends.sort_by(|x, y| x.low.total_cmp(&y.low));
    let (lo_end, hi_end) = (ends[0], ends[1]);
    let (hmin, hmax) = pair_height_range(c, a, band_x(&band.surface), rho, cp, cb, cr, branch, 0.0, TAU, 16384)
        .ok_or_else(|| unsupported("the conoid's height is not proved defined all the way round"))?;
    // Global separation is sufficient for every angular/radial point. If ranges
    // overlap, refuse even when a tighter pointwise proof might admit the cut.
    if !(hmin > lo_end.high + tol && hmax < hi_end.low - tol) {
        return Err(unsupported("the envelope reaches a pipe end (part of the end would stay uncut, or the pipe would be cut through)"));
    }
    let (far_end, toward) = if band.far < lo_end.low - tol && !lo_end.normalized {
        (hi_end, 1.0)
    } else if band.far > hi_end.high + tol && !hi_end.normalized {
        (lo_end, -1.0)
    } else {
        return Err(unsupported("the waste band's far cap is not beyond a planar pipe end"));
    };
    // Rebuild: the far cap and its rings kept; each wall trimmed back to its rim; the conoid between the rims.
    let mut vertices = vec![];
    let mut edges: Vec<Edge> = vec![];
    let mut vmap = std::collections::HashMap::new();
    let mut keep = |v: usize, vertices: &mut Vec<Vec3>| *vmap.entry(v).or_insert_with(|| {
        vertices.push(body.vertices[v]);
        vertices.len() - 1
    });
    let mut emap = std::collections::HashMap::new();
    let mut retained = vec![far_end.outer, far_end.inner];
    for &(edge, _) in body.faces[far_end.face].loops.iter().flatten() {
        if !retained.contains(&edge) { retained.push(edge); }
    }
    for ring in retained {
        let e = &body.edges[ring];
        let (s, t) = (keep(e.start, &mut vertices), keep(e.end, &mut vertices));
        edges.push(Edge { curve: e.curve.clone(), start: s, end: t });
        emap.insert(ring, edges.len() - 1);
    }
    // The seams on the far outer ring's generator: the rims' parameter starts there.
    let v_far = body.edges[far_end.outer].start;
    let x = perp(body.vertices[v_far].sub(c), a).unit()?;
    let mut rims = vec![];
    for (r, far_ring) in [(ro, far_end.outer), (ri, far_end.inner)] {
        let curve = Curve::RadialPair { c, a, x, rho, p: cp, b: cb, r: cr, branch, t0: 0.0, t1: TAU, s: r };
        let q = curve.start();
        let vf = body.edges[far_ring].start;
        if perp(body.vertices[vf].sub(q), a).length() > tol {
            return Err(unsupported("the pipe's far rings do not start on one generator"));
        }
        vertices.push(q);
        let v = vertices.len() - 1;
        edges.push(Edge { curve, start: v, end: v });
        let rim = edges.len() - 1;
        let vf = keep(vf, &mut vertices);
        edges.push(Edge { curve: Curve::Line { a: q, b: vertices[vf] }, start: v, end: vf });
        let seam = edges.len() - 1;
        let n = match body.edges[far_ring].curve {
            Curve::Arc { n, .. } => n,
            Curve::RadialPair { a, .. } => a,
            _ => unreachable!(),
        };
        // Whether the far ring runs the way t does (counter-clockwise about a).
        rims.push((rim, seam, emap[&far_ring], n.dot(a) > 0.0, v));
    }
    let (rim_o, seam_o, far_o, far_o_ccw, v_o) = rims[0];
    let (rim_i, seam_i, far_i, far_i_ccw, v_i) = rims[1];
    edges.push(Edge { curve: Curve::Line { a: vertices[v_i], b: vertices[v_o] }, start: v_i, end: v_o });
    let ruling = edges.len() - 1;
    let far_face = Face {
        surface: body.faces[far_end.face].surface.clone(),
        sense: body.faces[far_end.face].sense,
        loops: body.faces[far_end.face].loops.iter().map(|l| l.iter().map(|(e, f)| (emap[e], *f)).collect()).collect(),
    };
    // In each wall's chart (t round a, height along a; natural normal away from the axis) the outer wall runs
    // counter-clockwise, the inner one clockwise; the conoid (natural normal -s a + h' v) with the material above it
    // (toward = +1) has its natural normal outward.
    let (outer_face, inner_face, end_face) = if toward > 0.0 {
        (
            vec![(rim_o, true), (seam_o, true), (far_o, !far_o_ccw), (seam_o, false)],
            vec![(seam_i, true), (far_i, far_i_ccw), (seam_i, false), (rim_i, false)],
            (true, vec![(rim_i, true), (ruling, true), (rim_o, false), (ruling, false)]),
        )
    } else {
        (
            vec![(far_o, far_o_ccw), (seam_o, false), (rim_o, false), (seam_o, true)],
            vec![(seam_i, false), (rim_i, true), (seam_i, true), (far_i, !far_i_ccw)],
            (false, vec![(ruling, true), (rim_o, true), (ruling, false), (rim_i, false)]),
        )
    };
    let Surface::Conoid { c: sc, a: sa_, rho: sr, p: sp, b: sb_, r: srr, branch: sbr, .. } = band.surface else { unreachable!() };
    let faces = vec![
        far_face,
        Face { surface: body.faces[oi].surface.clone(), sense: true, loops: vec![outer_face] },
        Face { surface: body.faces[ii].surface.clone(), sense: false, loops: vec![inner_face] },
        Face { surface: Surface::Conoid { c: sc, a: sa_, x, rho: sr, p: sp, b: sb_, r: srr, branch: sbr }, sense: end_face.0, loops: vec![end_face.1] },
    ];
    let out = Body { vertices, edges, faces };
    crate::topology::check(&out, tol)?;
    Ok(out)
}

fn band_x(s: &Surface) -> Vec3 {
    match *s {
        Surface::Conoid { x, .. } => x,
        _ => unreachable!(),
    }
}

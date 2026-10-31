//! A pipe end coped by a round cutter at an oblique angle (Claude, 2026-10-31; Bob's Pipe 1.25 12 Coped 60 To
//! 2.375: a 1.25 pipe coped to a 2.375 OD run meeting it at 60 degrees, so the cutter's axis is turned 30 degrees
//! from square across the pipe and neither the bore nor the transverse chart (cardinal axes only) can take it).
//!
//! Domain, decided from the geometry before anything is built: a plain pipe (an outer and a coaxial inner wall, an
//! annular cap at each end) and a full round cutter (one cylinder, two caps) whose axis is oblique to the pipe's
//! (neither parallel nor square to it), whose cylinder holds one whole end annulus strictly inside it and clears
//! the other end's whole disk, and whose finite length covers every point the two solid cylinders share. Then each
//! generator of each wall is inside the cutter at the coped end and outside it at the far end, crossing its
//! cylinder exactly once between them (a line meets a solid cylinder in one interval), so each wall's new end is
//! one closed loop round the pipe: the exact CylinderPair (wire type 4) on that wall with the cutter's cylinder,
//! branch toward the material, its discriminant proved positive all the way round. The result: the far cap
//! unchanged, each wall trimmed back to its loop (its seam on the far ring's generator), and the cutter's cylinder
//! between the two loops facing into the cope. Anything else is refused by name; nothing is retried.
//!
//! Bounds (independent, from elementary geometry): with phi the angle between the axes and the two axes' common
//! normal n, a point within Ro of the pipe axis and within Rc of the cutter's has |w.b'| <= (Rc + Ro |cos phi|) /
//! sin phi along the pipe and |w.d| <= (Ro + Rc |cos phi|) / sin phi along the cutter (w measured from the
//! axes' closest points, b' and d the in-plane parts). The squared distance from the cutter's axis is a
//! trigonometric polynomial of degree two round any circle across the pipe: exact from eight samples, bounded
//! with the second-derivative chord bound (geom::pair_disc_lower_bound).
use crate::geom::{pair_disc_lower_bound, pair_eval, Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;
use std::f64::consts::TAU;

type Result<T> = std::result::Result<T, KernelError>;

fn unsupported(message: &str) -> KernelError {
    KernelError::Unsupported(message.into())
}

fn perp(v: Vec3, b: Vec3) -> Vec3 {
    v.sub(b.scale(v.dot(b)))
}

/// Selected when the tool is a full round cutter oblique to a pipe in the body (one outward cylinder with a coaxial
/// inward one). The cut itself admits the rest or refuses by name.
pub(super) fn selects(body: &Body, tool: &Body, tol: f64) -> bool {
    let Some((c0, c1, _)) = super::radial::cutter(tool, tol) else { return false };
    let Ok(d) = c1.sub(c0).unit() else { return false };
    body.faces.iter().any(|f| match f.surface {
        Surface::Cylinder { p, a, r, .. } if f.sense => {
            let cos = a.dot(d).abs();
            cos > 1e-9
                && cos < 1. - 1e-9
                && body.faces.iter().any(|g| {
                    matches!(g.surface, Surface::Cylinder { p: q, a: b, r: s, .. }
                        if !g.sense && s < r - tol && b.cross(a).length() <= 1e-9 && perp(q.sub(p), a).length() <= tol)
                })
        }
        _ => false,
    })
}

/// The squared distance from the line (q, d) round the circle centre o, radius r, in the plane spanned by e1, e2:
/// its trigonometric coefficients, exact from eight samples (degree two).
fn distance_coefficients(o: Vec3, e1: Vec3, e2: Vec3, r: f64, q: Vec3, d: Vec3) -> [f64; 5] {
    let mut k = [0.0; 5];
    for i in 0..8 {
        let t = std::f64::consts::PI * i as f64 / 4.0;
        let w = perp(o.add(e1.scale(r * t.cos())).add(e2.scale(r * t.sin())).sub(q), d);
        let s = w.dot(w);
        k[0] += s / 8.0;
        k[1] += s * t.cos() / 4.0;
        k[2] += s * t.sin() / 4.0;
        k[3] += s * (2.0 * t).cos() / 4.0;
        k[4] += s * (2.0 * t).sin() / 4.0;
    }
    k
}

fn neg(k: [f64; 5]) -> [f64; 5] {
    k.map(|v| -v)
}

/// The least squared distance from the line (q, d) over the disk centre o, radius r, plane e1, e2 (the plane not
/// parallel to d): the quadratic's minimum if it lies in the disk, else the boundary circle's proved lower bound.
fn disk_min(o: Vec3, e1: Vec3, e2: Vec3, r: f64, q: Vec3, d: Vec3) -> f64 {
    let (g, a, b) = (perp(o.sub(q), d), perp(e1, d), perp(e2, d));
    let (aa, ab, bb) = (a.dot(a), a.dot(b), b.dot(b));
    let det = aa * bb - ab * ab;
    let (ga, gb) = (g.dot(a), g.dot(b));
    let u = (-ga * bb + gb * ab) / det;
    let v = (-gb * aa + ga * ab) / det;
    if det > 0. && u.hypot(v) <= r {
        let w = g.add(a.scale(u)).add(b.scale(v));
        return w.dot(w);
    }
    pair_disc_lower_bound(&distance_coefficients(o, e1, e2, r, q, d), 0., TAU, 4096)
}

/// A plain wall: one loop of two full rings and one seam used both ways. Its (far ring, seam) given the far end's
/// ring set.
struct Wall {
    face: usize,
    r: f64,
    rings: [usize; 2],
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    crate::topology::check(body, tol)?;
    crate::topology::check(tool, tol)?;
    super::radial::planar_faces_outward(body)?;
    super::radial::planar_faces_outward(tool)?;
    let (c0, c1, rc) =
        super::radial::cutter(tool, tol).ok_or_else(|| unsupported("oblique cope needs a full round cutter"))?;
    let length = c1.sub(c0).length();
    let d = c1.sub(c0).unit()?;
    if body.faces.len() != 4 {
        return Err(unsupported("oblique cope needs a plain pipe (two walls and an annular cap at each end)"));
    }
    // The walls.
    let mut outer = None;
    let mut inner = None;
    for (i, f) in body.faces.iter().enumerate() {
        if let Surface::Cylinder { p, a, r, .. } = f.surface {
            let slot = if f.sense { &mut outer } else { &mut inner };
            if slot.is_some() {
                return Err(unsupported("oblique cope needs a plain pipe (two walls and an annular cap at each end)"));
            }
            *slot = Some((i, p, a, r));
        }
    }
    let (Some((oi, p, b, ro)), Some((ii, pi, bi, ri))) = (outer, inner) else {
        return Err(unsupported("oblique cope needs a plain pipe (two walls and an annular cap at each end)"));
    };
    if bi.cross(b).length() > 1e-9 || perp(pi.sub(p), b).length() > tol || !(ri < ro - tol) || ri <= tol {
        return Err(unsupported("oblique cope: the pipe's walls are not coaxial with the inner one inside"));
    }
    let cos = b.dot(d);
    if !(cos.abs() > 1e-9 && cos.abs() < 1. - 1e-9) {
        return Err(unsupported("oblique cope: the cutter is not oblique to the pipe"));
    }
    let sin = (1. - cos * cos).sqrt();
    // Each wall's loop: two full rings and a seam used both ways.
    let wall = |face: usize, r: f64| -> Result<Wall> {
        let f = &body.faces[face];
        let [lp] = &f.loops[..] else { return Err(unsupported("oblique cope: a pipe wall already has holes or cuts")) };
        let rings: Vec<usize> = lp.iter().filter(|(e, _)| matches!(body.edges[*e].curve, Curve::Arc { .. })).map(|(e, _)| *e).collect();
        let lines: Vec<usize> = lp.iter().filter(|(e, _)| matches!(body.edges[*e].curve, Curve::Line { .. })).map(|(e, _)| *e).collect();
        if lp.len() != 4 || rings.len() != 2 || lines.len() != 2 || lines[0] != lines[1] || rings[0] == rings[1] {
            return Err(unsupported("oblique cope: a pipe wall is not one plain band"));
        }
        for e in &rings {
            let Curve::Arc { c, n, r: ar, t0, t1, .. } = body.edges[*e].curve else { unreachable!() };
            if (t1 - t0 - TAU).abs() > 1e-9 || n.cross(b).length() > 1e-9 || (ar - r).abs() > tol || perp(c.sub(p), b).length() > tol {
                return Err(unsupported("oblique cope: a pipe wall's end is not a full ring"));
            }
        }
        Ok(Wall { face, r, rings: [rings[0], rings[1]] })
    };
    let (wo, wi) = (wall(oi, ro)?, wall(ii, ri)?);
    // The two ends along b, each with its cap face and its outer and inner ring.
    let caps: Vec<usize> = (0..4).filter(|i| *i != oi && *i != ii).collect();
    let mut ends = Vec::new();
    for &cf in &caps {
        let f = &body.faces[cf];
        let Surface::Plane { p: q, n, .. } = f.surface else {
            return Err(unsupported("oblique cope: a pipe end is not planar"));
        };
        if n.cross(b).length() > 1e-9 || f.loops.len() != 2 || f.loops.iter().any(|l| l.len() != 1) {
            return Err(unsupported("oblique cope: a pipe end is not a plain annulus square to the pipe"));
        }
        let level = q.sub(p).dot(b);
        let ring_of = |w: &Wall| {
            w.rings.iter().copied().find(|e| f.loops.iter().any(|l| l[0].0 == *e))
        };
        let (Some(ro_e), Some(ri_e)) = (ring_of(&wo), ring_of(&wi)) else {
            return Err(unsupported("oblique cope: a pipe end does not close its walls"));
        };
        ends.push((cf, level, ro_e, ri_e));
    }
    let centre = |level: f64| p.add(b.scale(level));
    let e1 = b.perpendicular();
    let e2 = b.cross(e1);
    let margin = 2. * rc * tol + tol * tol;
    // Coped end: its whole annulus strictly inside the cutter's cylinder (the outer ring bounds the convex distance).
    // Far end: its whole disk strictly outside it.
    let inside = |level: f64| -pair_disc_lower_bound(&neg(distance_coefficients(centre(level), e1, e2, ro, c0, d)), 0., TAU, 4096) < rc * rc - margin;
    let outside = |level: f64| disk_min(centre(level), e1, e2, ro, c0, d) > rc * rc + margin;
    let (coped, far) = if inside(ends[0].1) && outside(ends[1].1) {
        (ends[0], ends[1])
    } else if inside(ends[1].1) && outside(ends[0].1) {
        (ends[1], ends[0])
    } else {
        return Err(unsupported(
            "oblique cope must hold one whole pipe end inside the cutter and clear the other end",
        ));
    };
    // b' from the coped end into the material.
    let toward = if far.1 > coped.1 { b } else { b.scale(-1.) };
    // The axes' closest points: the shared region of the two solid cylinders is bounded about them.
    let w0 = c0.sub(p);
    let (bd, wb, wd) = (b.dot(d), w0.dot(b), w0.dot(d));
    let den = 1. - bd * bd;
    let sp = (wb - bd * wd) / den; // along b from p
    let sc = (bd * wb - wd) / den; // along d from c0
    // A shared point is O_b + alpha b + beta d + gamma n with |alpha| <= Rc / sin, |beta| <= Ro / sin, and lies in
    // the pipe only on the material side of the coped end: sigma (alpha + beta cos) >= sigma (coped - sp). Its place
    // along the cutter, alpha cos + beta from O_d, is extreme at a vertex of that clipped box.
    let (big_a, big_b) = (rc / sin, ro / sin);
    let sigma = toward.dot(b);
    let floor = coped.1 - sp;
    let mut corners = Vec::new();
    for alpha in [-big_a, big_a] {
        for beta in [-big_b, big_b] {
            corners.push((alpha, beta));
        }
        corners.push((alpha, (floor - alpha) / cos));
    }
    for beta in [-big_b, big_b] {
        corners.push((floor - beta * cos, beta));
    }
    let along: Vec<f64> = corners
        .into_iter()
        .filter(|(a, bt)| {
            a.abs() <= big_a * (1. + 1e-12) && bt.abs() <= big_b * (1. + 1e-12) && sigma * (a + bt * cos) >= sigma * floor - 1e-12
        })
        .map(|(a, bt)| sc + a * cos + bt)
        .collect();
    let (lo_d, hi_d) = along.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), v| (l.min(*v), h.max(*v)));
    if along.is_empty() || lo_d < tol || hi_d > length - tol {
        return Err(unsupported("oblique cope: the cutter's ends are inside the pipe"));
    }
    // Each wall's loop: the pair on the wall (parameter cylinder) with the cutter's (root), branch toward b'.
    let mut out = body.clone();
    let start = centre(coped.1);
    let mut loops = Vec::new();
    for (w, far_ring) in [(&wo, far.2), (&wi, far.3)] {
        let v_far = out.edges[far_ring].start;
        let x = perp(out.vertices[v_far].sub(p), b).unit()?;
        let k = crate::geom::pair_disc_coefficients(start, toward, x, w.r, c0, d, rc);
        if pair_disc_lower_bound(&k, 0., TAU, 4096) <= tol * tol {
            return Err(unsupported("oblique cope: a wall's loop is not proved clear of tangency"));
        }
        let curve = Curve::CylinderPair { c: start, a: toward, x, rho: w.r, p: c0, b: d, r: rc, branch: 1., t0: 0., t1: TAU };
        let q = pair_eval(start, toward, x, w.r, c0, d, rc, 1., 0.).0.add(start);
        if (perp(out.vertices[v_far].sub(q), b).length()) > tol {
            return Err(KernelError::Internal("oblique cope: a seam is not a generator".into()));
        }
        loops.push((curve, q, v_far, far_ring));
    }
    // Rebuild: keep the far cap and its two rings; drop the coped cap, its rings and the old seams.
    let mut vertices = Vec::new();
    let mut edges = Vec::new();
    let mut vmap = std::collections::HashMap::new();
    let mut keep_vertex = |v: usize, vertices: &mut Vec<Vec3>| *vmap.entry(v).or_insert_with(|| {
        vertices.push(out.vertices[v]);
        vertices.len() - 1
    });
    let mut emap = std::collections::HashMap::new();
    for ring in [far.2, far.3] {
        let e = &out.edges[ring];
        let (s, t) = (keep_vertex(e.start, &mut vertices), keep_vertex(e.end, &mut vertices));
        edges.push(Edge { curve: e.curve.clone(), start: s, end: t });
        emap.insert(ring, edges.len() - 1);
    }
    let mut new = Vec::new();
    for (curve, q, v_far, far_ring) in loops {
        vertices.push(q);
        let v = vertices.len() - 1;
        edges.push(Edge { curve, start: v, end: v });
        let pair = edges.len() - 1;
        let vf = keep_vertex(v_far, &mut vertices);
        edges.push(Edge { curve: Curve::Line { a: q, b: vertices[vf] }, start: v, end: vf });
        let seam = edges.len() - 1;
        let Curve::Arc { n, .. } = out.edges[far_ring].curve else { unreachable!() };
        new.push((pair, seam, emap[&far_ring], n.dot(toward) > 0.));
    }
    let far_face = Face {
        surface: out.faces[far.0].surface.clone(),
        sense: out.faces[far.0].sense,
        loops: out.faces[far.0].loops.iter().map(|l| l.iter().map(|(e, f)| (emap[e], *f)).collect()).collect(),
    };
    // Outer wall (outward radial): the pair counter-clockwise about b' (the material on its b' side), the far ring
    // the other way. Inner wall (outward toward the axis): both reversed.
    let (po, so, ro_far, ro_ccw) = new[0];
    let (pin, si, ri_far, ri_ccw) = new[1];
    let outer_face = Face {
        surface: out.faces[wo.face].surface.clone(),
        sense: true,
        loops: vec![vec![(po, true), (so, true), (ro_far, !ro_ccw), (so, false)]],
    };
    let inner_face = Face {
        surface: out.faces[wi.face].surface.clone(),
        sense: false,
        loops: vec![vec![(pin, false), (si, true), (ri_far, ri_ccw), (si, false)]],
    };
    // The cutter's cylinder between the loops, facing its axis; its parameter starts where the pipe's axis leaves it.
    let axis_root = pair_eval(start, toward, e1, 0., c0, d, rc, 1., 0.).0.add(start);
    let x = perp(axis_root.sub(c0), d).unit()?;
    let cope_face = Face {
        surface: Surface::Cylinder { p: c0, a: d, x, r: rc },
        sense: false,
        loops: vec![vec![(po, false)], vec![(pin, true)]],
    };
    out = Body { vertices, edges, faces: vec![far_face, outer_face, inner_face, cope_face] };
    crate::topology::check(&out, tol)?;
    Ok(out)
}

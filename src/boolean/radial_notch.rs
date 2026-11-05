//! Centered radial notch through one flat annular pipe end. Complete simple
//! axisymmetric boundaries are admitted before splicing; no stock is rebuilt.
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;
use std::f64::consts::{PI, TAU};
type Result<T> = std::result::Result<T, KernelError>;
fn no(s: &str) -> KernelError {
    KernelError::Unsupported(format!("radial end notch: {s}"))
}
fn parallel(a: Vec3, b: Vec3) -> bool {
    a.cross(b).length() == 0.
}

/// A full annular cap and a cutter whose axis lies exactly in that plane and
/// passes through its axis point. This only selects a domain; admission below
/// proves the whole boundary and finite cutter before any geometry changes.
fn cap(body: &Body, tool: &Body, tol: f64) -> Option<(usize, Vec3, Vec3)> {
    let (s, e, _) = super::radial::cutter(tool, tol)?;
    let mut found = None;
    for (fi, f) in body.faces.iter().enumerate() {
        let Surface::Plane { n, .. } = f.surface else {
            continue;
        };
        if f.loops.len() != 2 || f.loops.iter().any(|lp| lp.len() != 1) {
            continue;
        }
        let Curve::Arc {
            c, n: axis, t0, t1, ..
        } = body.edges.get(f.loops[0][0].0)?.curve
        else {
            continue;
        };
        if t1 - t0 != TAU || !parallel(n, axis) {
            continue;
        }
        let Curve::Arc {
            c: c2,
            t0: a,
            t1: b,
            ..
        } = body.edges.get(f.loops[1][0].0)?.curve
        else {
            continue;
        };
        if c2 != c
            || b - a != TAU
            || s.sub(c).dot(n) != 0.
            || e.sub(c).dot(n) != 0.
            || !parallel(s.sub(c), e.sub(c))
        {
            continue;
        }
        if found.is_some() {
            return None;
        }
        found = Some((fi, c, n.scale(if f.sense { -1. } else { 1. })));
    }
    found
}
pub(super) fn selects(body: &Body, tool: &Body, tol: f64) -> bool {
    // Domain selection proves eligibility before execution. In particular,
    // existing blind-pocket domains keep cutters that do not reach the bore.
    cap(body, tool, tol).is_some() && admission(body, tool, tol).is_ok()
}

struct Band {
    face: usize,
    rings: [(usize, f64); 2],
}
fn band(body: &Body, fi: usize, c: Vec3, u: Vec3) -> Result<Band> {
    let f = &body.faces[fi];
    let [lp] = &f.loops[..] else {
        return Err(no("wall has earlier openings"));
    };
    if lp.len() != 4 {
        return Err(no("wall is not one full band"));
    }
    let mut rings = vec![];
    let mut seams = vec![];
    for &(ei, fw) in lp {
        let e = &body.edges[ei];
        match e.curve {
            Curve::Arc {
                c: q, n, t0, t1, ..
            } if parallel(q.sub(c), u) && parallel(n, u) && t1 - t0 == TAU && e.start == e.end => {
                rings.push((ei, q.sub(c).dot(u), fw, n.dot(u)))
            }
            Curve::Line { .. } => seams.push((ei, fw)),
            _ => return Err(no("wall needs full coaxial rings and a straight seam")),
        }
    }
    if rings.len() != 2 || seams.len() != 2 || seams[0].0 != seams[1].0 || seams[0].1 == seams[1].1
    {
        return Err(no(
            "wall band connectivity is not two rings and one repeated seam",
        ));
    }
    rings.sort_by(|a, b| a.1.total_cmp(&b.1));
    let edge = &body.edges[seams[0].0];
    let ends = [body.edges[rings[0].0].start, body.edges[rings[1].0].start];
    if ends[0] == ends[1] || !ends.contains(&edge.start) || !ends.contains(&edge.end) {
        return Err(no("band seam does not join its rings"));
    }
    for (i, &(_, _, fw, dir)) in rings.iter().enumerate() {
        if (fw == (dir > 0.)) != ((i == 0) == f.sense) {
            return Err(no("wall loop is inside out"));
        }
    }
    Ok(Band {
        face: fi,
        rings: [(rings[0].0, rings[0].1), (rings[1].0, rings[1].1)],
    })
}

struct Plan {
    cap: usize,
    c: Vec3,
    u: Vec3,
    d: Vec3,
    x: Vec3,
    rho: f64,
    walls: [Band; 2],
    radii: [f64; 2],
    near_cone: Option<Band>,
}
fn plan(body: &Body, tool: &Body, tol: f64) -> Result<Plan> {
    let (cap, c, u) =
        cap(body, tool, tol).ok_or_else(|| no("cutter is not centered on one flat annular end"))?;
    if body.faces.len() != 4 && body.faces.len() != 5 {
        return Err(no("requires a simple pipe with at most one outer bevel"));
    }
    let (a, b, rho) =
        super::radial::cutter(tool, tol).ok_or_else(|| no("not a full round cutter"))?;
    let (s, e) = if a.sub(c).length() > b.sub(c).length() {
        (a, b)
    } else {
        (b, a)
    };
    let d = s.sub(c).unit()?;
    let x = u.cross(d);
    let tool_wall = tool
        .faces
        .iter()
        .position(|f| matches!(f.surface, Surface::Cylinder { .. }))
        .ok_or_else(|| no("cutter has no cylindrical wall"))?;
    if !tool.faces[tool_wall].sense {
        return Err(no("cutter cylindrical wall is inverted"));
    }
    let tool_band = band(tool, tool_wall, c, d)?;
    for (ei, _) in tool_band.rings {
        if !matches!(tool.edges[ei].curve,Curve::Arc{r,..} if r==rho) {
            return Err(no("cutter cap radius differs from its wall"));
        }
    }
    for f in &tool.faces {
        if let Surface::Plane { p, n, .. } = f.surface {
            if f.loops.len() != 1 || f.loops[0].len() != 1 || !parallel(n, d) {
                return Err(no("cutter cap is not one full axis-normal disk"));
            }
            let outward = n.scale(if f.sense { 1. } else { -1. }).dot(d);
            if !((p.sub(s).dot(d) == 0. && outward > 0.) || (p.sub(e).dot(d) == 0. && outward < 0.))
            {
                return Err(no("cutter cap points into its material"));
            }
        }
    }
    let mut cyl = vec![];
    let mut cones = vec![];
    let mut caps = vec![];
    for (fi, f) in body.faces.iter().enumerate() {
        match f.surface {
            Surface::Cylinder { p, a, r, .. } if parallel(a, u) && parallel(p.sub(c), u) => {
                cyl.push((fi, r))
            }
            Surface::Cone { p, a, .. } if parallel(a, u) && parallel(p.sub(c), u) => cones.push(fi),
            Surface::Plane { p, n, .. } if parallel(n, u) => caps.push((fi, p.sub(c).dot(u))),
            _ => return Err(no("boundary is not a simple coaxial pipe")),
        }
    }
    if cyl.len() != 2 || caps.len() != 2 || cones.len() > 1 {
        return Err(no(
            "requires two cylindrical walls, two annular ends, and at most one cone",
        ));
    }
    cyl.sort_by(|a, b| a.1.total_cmp(&b.1));
    let (ri, ro) = (cyl[0].1, cyl[1].1);
    if body.faces[cyl[0].0].sense
        || !body.faces[cyl[1].0].sense
        || !(ri > rho + tol && ro > ri + tol)
    {
        return Err(no("wall radii/senses or cutter clearance are invalid"));
    }
    if !(s.sub(c).dot(d) > ro + tol
        && e.sub(c).dot(d) > 0.
        && (e.sub(c).dot(d).powi(2) + rho * rho).sqrt() < ri - tol)
    {
        return Err(no(
            "cutter must start outside and finish strictly inside the near bore",
        ));
    }
    let inner = band(body, cyl[0].0, c, u)?;
    let outer = band(body, cyl[1].0, c, u)?;
    let length = inner.rings[1].1;
    let transition = outer.rings[1].1;
    let near = outer.rings[0].1;
    if inner.rings[0].1 != 0.
        || !(near >= 0. && near < rho - tol)
        || !(transition > rho + tol && length >= transition)
    {
        return Err(no(
            "hole reaches a wall transition or does not cross its flat end",
        ));
    }
    for &(fi, z) in &caps {
        let f = &body.faces[fi];
        if f.loops.len() != 2 || f.loops.iter().any(|lp| lp.len() != 1) || !(z == 0. || z == length)
        {
            return Err(no("end is not one full annulus at the wall limits"));
        }
        for lp in &f.loops {
            let Curve::Arc {
                c: q, n, t0, t1, ..
            } = body.edges[lp[0].0].curve
            else {
                return Err(no("end has a noncircular rim"));
            };
            if q != c.add(u.scale(z)) || !parallel(n, u) || t1 - t0 != TAU {
                return Err(no("end rim is not its exact full coaxial circle"));
            }
        }
    }
    let mut near_cone = None;
    if cones.is_empty() {
        if transition != length || near != 0. {
            return Err(no("straight walls have different end levels"));
        }
    } else {
        let cone = band(body, cones[0], c, u)?;
        let radius = |ei: usize| match body.edges[ei].curve {
            Curve::Arc { r, .. } => r,
            _ => unreachable!(),
        };
        if !body.faces[cone.face].sense {
            return Err(no("outer cone is inverted"));
        }
        if near == 0. {
            if cone.rings[0].1 != transition
                || cone.rings[1].1 != length
                || !(length > transition + tol)
                || radius(cone.rings[0].0) != ro
                || !(radius(cone.rings[1].0) > ri + tol && radius(cone.rings[1].0) < ro)
            {
                return Err(no("distant cone is not an outer end bevel"));
            }
        } else {
            let r0 = radius(cone.rings[0].0);
            if transition != length
                || cone.rings[0].1 != 0.
                || cone.rings[1] != outer.rings[0]
                || !(near > tol && r0 > ri + tol && r0 < ro)
                || radius(cone.rings[1].0) != ro
            {
                return Err(no(
                    "near cone does not join the flat end to the outer cylinder",
                ));
            }
            // The increasing radius is bounded below by r0 > ri > rho;
            // therefore both cone/cutter intersections have positive discriminant.
            let Surface::Cone {
                p: cp,
                a: ca,
                r: cr,
                k,
                ..
            } = body.faces[cone.face].surface
            else {
                unreachable!()
            };
            if !(k * ca.dot(u) > 0.)
                || cr + k * c.sub(cp).dot(ca) != r0
                || cr + k * c.add(u.scale(near)).sub(cp).dot(ca) != ro
            {
                return Err(no(
                    "near cone radius/frame does not agree with the end ring",
                ));
            }
            near_cone = Some(cone);
        }
    }
    // No hidden extra boundaries or unused topology: these counts are exact
    // for the fully checked two-band shell, plus one optional conical band.
    let expected = if cones.is_empty() { (6, 4) } else { (8, 5) };
    if (body.edges.len(), body.vertices.len()) != expected {
        return Err(no(
            "additional edges or vertices are outside the supported shell",
        ));
    }
    Ok(Plan {
        cap,
        c,
        u,
        d,
        x,
        rho,
        walls: [inner, outer],
        radii: [ri, ro],
        near_cone,
    })
}

fn add_edge(out: &mut Body, curve: Curve, start: usize, end: usize) -> usize {
    let ei = out.edges.len();
    out.edges.push(Edge { curve, start, end });
    ei
}
/// Join a degree-two boundary from an oriented first edge. Ambiguity is an
/// internal construction error, never repaired by proximity or point welding.
fn loop_from(out: &Body, edges: &[usize], first: (usize, bool)) -> Result<Vec<(usize, bool)>> {
    let ends = |e: usize, f: bool| {
        if f {
            (out.edges[e].start, out.edges[e].end)
        } else {
            (out.edges[e].end, out.edges[e].start)
        }
    };
    let (start, mut at) = ends(first.0, first.1);
    let mut lp = vec![first];
    while lp.len() < edges.len() {
        let next: Vec<_> = edges
            .iter()
            .copied()
            .filter(|e| !lp.iter().any(|v| v.0 == *e))
            .filter_map(|e| {
                if out.edges[e].start == at {
                    Some((e, true))
                } else if out.edges[e].end == at {
                    Some((e, false))
                } else {
                    None
                }
            })
            .collect();
        let [v] = &next[..] else {
            return Err(KernelError::Internal(
                "radial notch boundary does not have one continuation".into(),
            ));
        };
        at = ends(v.0, v.1).1;
        lp.push(*v);
    }
    if at != start {
        return Err(KernelError::Internal(
            "radial notch boundary is not closed".into(),
        ));
    }
    Ok(lp)
}

type Split = (Vec3, Vec3, Vec3, f64, f64, f64, f64, f64);
fn admission(body: &Body, tool: &Body, tol: f64) -> Result<(Plan, Vec<Split>)> {
    crate::topology::check(body, tol)?;
    crate::topology::check(tool, tol)?;
    super::radial::planar_faces_outward(body)?;
    super::radial::planar_faces_outward(tool)?;
    let p = plan(body, tool, tol)?;
    // Check both periodic seams before constructing anything. A ring is split
    // only when its removed minor arc stays strictly inside its existing chart.
    let mut splits = vec![];
    let mut rings: Vec<_> = p
        .walls
        .iter()
        .zip(p.radii)
        .map(|(band, r)| {
            (
                band.rings[0].0,
                r,
                (p.rho * p.rho - band.rings[0].1.powi(2)).sqrt(),
            )
        })
        .collect();
    if let Some(cone) = &p.near_cone {
        let ring = cone.rings[0].0;
        let Curve::Arc { r, .. } = body.edges[ring].curve else {
            unreachable!()
        };
        rings.push((ring, r, p.rho));
    }
    for (ring, r, width) in rings {
        let Curve::Arc {
            c,
            n,
            x,
            r: ar,
            t0,
            t1,
        } = body.edges[ring].curve
        else {
            unreachable!()
        };
        if ar != r {
            return Err(no("end ring and wall radii differ"));
        }
        let center = (p.d.dot(n.cross(x)).atan2(p.d.dot(x)) - t0).rem_euclid(TAU) + t0;
        let half = (width / r).asin();
        let (lo, hi) = (center - half, center + half);
        if !(r * (lo - t0) > tol && r * (t1 - hi) > tol) {
            return Err(no("notch meets the existing periodic seam"));
        }
        splits.push((c, n, x, r, t0, t1, lo, hi));
    }
    Ok((p, splits))
}
pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    let (p, splits) = admission(body, tool, tol)?;
    if p.near_cone.is_some() {
        return cut_near_bevel(body, &p, &splits, tol);
    }
    let mut out = body.clone();
    let mut pieces = vec![];
    for (i, band) in p.walls.iter().enumerate() {
        let (c, n, x, r, t0, t1, lo, hi) = splits[i];
        let pair = Curve::CylinderPair {
            c: p.c,
            a: p.d,
            x: p.x,
            rho: p.rho,
            p: p.c,
            b: p.u,
            r,
            branch: 1.,
            t0: 0.,
            t1: PI,
        };
        let plus = out.vertices.len();
        out.vertices.push(pair.eval(0.).0);
        let minus = out.vertices.len();
        out.vertices.push(pair.eval(PI).0);
        let forward = n.dot(p.u) < 0.; // ring lo -> hi follows or reverses the pair
        let (vlow, vhigh) = if forward {
            (plus, minus)
        } else {
            (minus, plus)
        };
        let ring = band.rings[0].0;
        let seam = out.edges[ring].start;
        out.edges[ring] = Edge {
            curve: Curve::Arc {
                c,
                n,
                x,
                r,
                t0,
                t1: lo,
            },
            start: seam,
            end: vlow,
        };
        let last = add_edge(
            &mut out,
            Curve::Arc {
                c,
                n,
                x,
                r,
                t0: hi,
                t1,
            },
            vhigh,
            seam,
        );
        let meeting = add_edge(&mut out, pair, plus, minus);
        let old = &body.faces[band.face].loops[0];
        let ring_fw = old.iter().find(|v| v.0 == ring).unwrap().1;
        let path = [(ring, true), (meeting, forward), (last, true)];
        let mut lp = vec![];
        for &(ei, fw) in old {
            if ei != ring {
                lp.push((ei, fw));
            } else if fw {
                lp.extend(path);
            } else {
                lp.extend(path.iter().rev().map(|(e, f)| (*e, !*f)));
            }
        }
        out.faces[band.face].loops = vec![lp];
        pieces.push((ring, last, meeting, plus, minus, forward == ring_fw));
    }
    let mut sides = vec![];
    for side in [3, 4] {
        let v = |i: usize| if side == 3 { pieces[i].3 } else { pieces[i].4 };
        let (a, b) = (v(0), v(1));
        let curve = Curve::Line {
            a: out.vertices[a],
            b: out.vertices[b],
        };
        sides.push(add_edge(&mut out, curve, a, b));
    }
    let outer_ring_fw = body.faces[p.cap]
        .loops
        .iter()
        .flatten()
        .find(|v| v.0 == pieces[1].0)
        .unwrap()
        .1;
    let cap_edges = [
        pieces[0].0,
        pieces[0].1,
        pieces[1].0,
        pieces[1].1,
        sides[0],
        sides[1],
    ];
    out.faces[p.cap].loops = vec![loop_from(&out, &cap_edges, (pieces[1].0, outer_ring_fw))?];
    let hole_edges = [pieces[0].2, pieces[1].2, sides[0], sides[1]];
    let hole_loop = loop_from(&out, &hole_edges, (pieces[1].2, !pieces[1].5))?;
    out.faces.push(Face {
        surface: Surface::Cylinder {
            p: p.c,
            a: p.d,
            x: p.x,
            r: p.rho,
        },
        sense: false,
        loops: vec![hole_loop],
    });
    crate::topology::check(&out, tol)?;
    super::radial::planar_faces_outward(&out)?;
    Ok(out)
}

/// A near bevel has two cone/cutter edges joined by the central cylindrical
/// edge. The supplied cone's seam becomes unnecessary once the notch joins
/// both of its rings: reuse that edge slot, retaining all distant edge indices.
fn cut_near_bevel(body: &Body, p: &Plan, splits: &[Split], tol: f64) -> Result<Body> {
    let cone = p.near_cone.as_ref().unwrap();
    let beta = (cone.rings[1].1 / p.rho).asin();
    let mut out = body.clone();
    // (original ring, second complement arc, intersection, plus, minus,
    // intersection direction in its wall).
    let mut pieces = vec![];
    for (i, band) in p.walls.iter().enumerate() {
        let (c, n, x, r, t0, t1, lo, hi) = splits[i];
        let (a, b) = if i == 0 { (0., PI) } else { (beta, PI - beta) };
        let pair = Curve::CylinderPair {
            c: p.c,
            a: p.d,
            x: p.x,
            rho: p.rho,
            p: p.c,
            b: p.u,
            r,
            branch: 1.,
            t0: a,
            t1: b,
        };
        let plus = out.vertices.len();
        out.vertices.push(pair.eval(a).0);
        let minus = out.vertices.len();
        out.vertices.push(pair.eval(b).0);
        let forward = n.dot(p.u) < 0.;
        let (vlow, vhigh) = if forward {
            (plus, minus)
        } else {
            (minus, plus)
        };
        let ring = band.rings[0].0;
        let seam = out.edges[ring].start;
        out.edges[ring] = Edge {
            curve: Curve::Arc {
                c,
                n,
                x,
                r,
                t0,
                t1: lo,
            },
            start: seam,
            end: vlow,
        };
        let last = add_edge(
            &mut out,
            Curve::Arc {
                c,
                n,
                x,
                r,
                t0: hi,
                t1,
            },
            vhigh,
            seam,
        );
        let meeting = add_edge(&mut out, pair, plus, minus);
        let old = &body.faces[band.face].loops[0];
        let ring_fw = old.iter().find(|v| v.0 == ring).unwrap().1;
        let path = [(ring, true), (meeting, forward), (last, true)];
        let mut lp = vec![];
        for &(ei, fw) in old {
            if ei != ring {
                lp.push((ei, fw));
            } else if fw {
                lp.extend(path);
            } else {
                lp.extend(path.iter().rev().map(|(e, f)| (*e, !*f)));
            }
        }
        out.faces[band.face].loops = vec![lp];
        pieces.push((ring, last, meeting, plus, minus, forward == ring_fw));
    }
    let Surface::Cone {
        p: cp,
        a: ca,
        r: cr,
        k,
        ..
    } = body.faces[cone.face].surface
    else {
        unreachable!()
    };
    let conic = |a, b| Curve::ConePair {
        c: p.c,
        a: p.d,
        x: p.x,
        rho: p.rho,
        p: cp,
        b: ca,
        r: cr,
        k,
        branch: 1.,
        t0: a,
        t1: b,
    };
    let first_curve = conic(0., beta);
    let last_curve = conic(PI - beta, PI);
    let plus = out.vertices.len();
    out.vertices.push(first_curve.eval(0.).0);
    let minus = out.vertices.len();
    out.vertices.push(last_curve.eval(PI).0);
    let seam_edge = body.faces[cone.face].loops[0]
        .iter()
        .find(|(ei, _)| matches!(body.edges[*ei].curve, Curve::Line { .. }))
        .unwrap()
        .0;
    out.edges[seam_edge] = Edge {
        curve: first_curve,
        start: plus,
        end: pieces[1].3,
    };
    let last_cone = add_edge(&mut out, last_curve, pieces[1].4, minus);
    let (c, n, x, r, t0, t1, lo, hi) = splits[2];
    let forward = n.dot(p.u) < 0.;
    let (vlow, vhigh) = if forward {
        (plus, minus)
    } else {
        (minus, plus)
    };
    let ring = cone.rings[0].0;
    let seam = out.edges[ring].start;
    out.edges[ring] = Edge {
        curve: Curve::Arc {
            c,
            n,
            x,
            r,
            t0,
            t1: lo,
        },
        start: seam,
        end: vlow,
    };
    let last = add_edge(
        &mut out,
        Curve::Arc {
            c,
            n,
            x,
            r,
            t0: hi,
            t1,
        },
        vhigh,
        seam,
    );
    let cone_edges = [ring, last, pieces[1].0, pieces[1].1, seam_edge, last_cone];
    let ring_fw = body.faces[cone.face].loops[0]
        .iter()
        .find(|v| v.0 == ring)
        .unwrap()
        .1;
    out.faces[cone.face].loops = vec![loop_from(&out, &cone_edges, (ring, ring_fw))?];

    let mut sides = vec![];
    for (a, b) in [(pieces[0].3, plus), (pieces[0].4, minus)] {
        let line = Curve::Line {
            a: out.vertices[a],
            b: out.vertices[b],
        };
        sides.push(add_edge(&mut out, line, a, b));
    }
    let cap_edges = [pieces[0].0, pieces[0].1, ring, last, sides[0], sides[1]];
    let cap_fw = body.faces[p.cap]
        .loops
        .iter()
        .flatten()
        .find(|v| v.0 == ring)
        .unwrap()
        .1;
    out.faces[p.cap].loops = vec![loop_from(&out, &cap_edges, (ring, cap_fw))?];
    let hole_edges = [
        pieces[0].2,
        pieces[1].2,
        seam_edge,
        last_cone,
        sides[0],
        sides[1],
    ];
    let hole_loop = loop_from(&out, &hole_edges, (pieces[1].2, !pieces[1].5))?;
    out.faces.push(Face {
        surface: Surface::Cylinder {
            p: p.c,
            a: p.d,
            x: p.x,
            r: p.rho,
        },
        sense: false,
        loops: vec![hole_loop],
    });
    crate::topology::check(&out, tol)?;
    super::radial::planar_faces_outward(&out)?;
    Ok(out)
}

/// A bounded whole-boundary certificate, independent of feature history.
/// If B equals the admitted S\T boundary exactly, set difference gives
/// B\T = (S\T)\T = B. Only the original B is returned by this domain.
/// The temporary inverse splice is proof data, never replacement stock.
pub(super) fn already_removed(body: &Body, tool: &Body, tol: f64) -> bool {
    fn prove(body: &Body, tool: &Body, tol: f64) -> Result<bool> {
        let counts = (body.faces.len(), body.edges.len(), body.vertices.len());
        let near_bevel = counts == (6, 16, 11);
        if counts != (5, 12, 8) && counts != (6, 14, 9) && !near_bevel {
            return Ok(false);
        }
        crate::topology::check(body, tol)?;
        crate::topology::check(tool, tol)?;
        let base_edges = body.edges.len() - if near_bevel { 8 } else { 6 };
        let base_vertices = body.vertices.len() - if near_bevel { 6 } else { 4 };
        let base_faces = body.faces.len() - 1;
        let mut proof = body.clone();
        let mut rings = vec![];
        for i in 0..2 {
            let last = base_edges + 2 * i;
            let meeting = last + 1;
            if !matches!(body.edges[meeting].curve, Curve::CylinderPair { .. }) {
                return Ok(false);
            }
            let Curve::Arc {
                c,
                n,
                x,
                r,
                t0: hi,
                t1,
            } = body.edges[last].curve
            else {
                return Ok(false);
            };
            let walls: Vec<_> = body.faces[..base_faces]
                .iter()
                .enumerate()
                .filter(|(_, f)| f.loops.len() == 1 && f.loops[0].iter().any(|v| v.0 == meeting))
                .collect();
            let [(fi, face)] = &walls[..] else {
                return Ok(false);
            };
            let arcs: Vec<_> = face.loops[0]
                .iter()
                .copied()
                .filter(|(ei, _)| {
                    *ei < base_edges
                        && matches!(body.edges[*ei].curve,
                    Curve::Arc { c: q, n: a, x: xx, r: rr, t0, t1: lo }
                    if q==c && a==n && xx==x && rr==r && lo<hi && t1-t0==TAU)
                })
                .collect();
            let [(ring, fw)] = &arcs[..] else {
                return Ok(false);
            };
            let Curve::Arc { t0, .. } = body.edges[*ring].curve else {
                unreachable!()
            };
            let seam = body.edges[*ring].start;
            if seam != body.edges[last].end {
                return Ok(false);
            }
            proof.edges[*ring] = Edge {
                curve: Curve::Arc { c, n, x, r, t0, t1 },
                start: seam,
                end: seam,
            };
            let mut restored = vec![];
            let mut inserted = false;
            for &(ei, dir) in &face.loops[0] {
                if [*ring, last, meeting].contains(&ei) {
                    if !inserted {
                        restored.push((*ring, *fw));
                        inserted = true;
                    }
                } else {
                    restored.push((ei, dir));
                }
            }
            proof.faces[*fi].loops = vec![restored];
            rings.push(*ring);
        }
        if rings[0] == rings[1] {
            return Ok(false);
        }
        if near_bevel {
            // The third split ring bounds the end cap and the cone. The
            // cone's old seam slot now holds its first intersection segment.
            let last = base_edges + 5;
            let Curve::Arc {
                c,
                n,
                x,
                r,
                t0: hi,
                t1,
            } = body.edges[last].curve
            else {
                return Ok(false);
            };
            let cones: Vec<_> = body.faces[..base_faces]
                .iter()
                .enumerate()
                .filter(|(_, f)| {
                    matches!(f.surface, Surface::Cone { .. })
                        && f.loops.len() == 1
                        && f.loops[0].iter().any(|v| v.0 == last)
                })
                .collect();
            let [(fi, face)] = &cones[..] else {
                return Ok(false);
            };
            let arcs: Vec<_> = face.loops[0]
                .iter()
                .copied()
                .filter(|(ei, _)| {
                    *ei < base_edges
                        && matches!(body.edges[*ei].curve,
                    Curve::Arc { c: q, n: a, x: xx, r: rr, t0, t1: lo }
                    if q == c && a == n && xx == x && rr == r && lo < hi && t1 - t0 == TAU)
                })
                .collect();
            let [(ring, fw)] = &arcs[..] else {
                return Ok(false);
            };
            let seams: Vec<_> = face.loops[0]
                .iter()
                .filter(|(ei, _)| {
                    *ei < base_edges && matches!(body.edges[*ei].curve, Curve::ConePair { .. })
                })
                .collect();
            let [seam] = &seams[..] else {
                return Ok(false);
            };
            let upper: Vec<_> = face.loops[0]
                .iter()
                .copied()
                .filter(|v| v.0 == rings[1])
                .collect();
            let [upper] = &upper[..] else {
                return Ok(false);
            };
            if rings.contains(ring) {
                return Ok(false);
            }
            let Curve::Arc { t0, .. } = body.edges[*ring].curve else {
                unreachable!();
            };
            let start = body.edges[*ring].start;
            let end = proof.edges[rings[1]].start;
            if start != body.edges[last].end || start >= base_vertices || end >= base_vertices {
                return Ok(false);
            }
            proof.edges[*ring] = Edge {
                curve: Curve::Arc { c, n, x, r, t0, t1 },
                start,
                end: start,
            };
            // Only a proof seam is restored. It is replaced during replay,
            // and is never returned as geometry from this no-change domain.
            proof.edges[seam.0] = Edge {
                curve: Curve::Line {
                    a: body.vertices[start],
                    b: body.vertices[end],
                },
                start,
                end,
            };
            proof.faces[*fi].loops =
                vec![vec![(*ring, *fw), (seam.0, true), *upper, (seam.0, false)]];
            rings[1] = *ring;
        }
        let offset = if near_bevel { 6 } else { 4 };
        let sides = [base_edges + offset, base_edges + offset + 1];
        let caps: Vec<_> = body.faces[..base_faces]
            .iter()
            .enumerate()
            .filter(|(_, f)| {
                matches!(f.surface, Surface::Plane { .. })
                    && sides
                        .iter()
                        .all(|ei| f.loops.iter().flatten().any(|v| v.0 == *ei))
            })
            .collect();
        let [(cap, face)] = &caps[..] else {
            return Ok(false);
        };
        let mut cap_loops = vec![];
        for &ring in rings.iter().rev() {
            let references: Vec<_> = face
                .loops
                .iter()
                .flatten()
                .copied()
                .filter(|v| v.0 == ring)
                .collect();
            let [item] = &references[..] else {
                return Ok(false);
            };
            cap_loops.push(vec![*item]);
        }
        proof.faces[*cap].loops = cap_loops;
        proof.faces.truncate(base_faces);
        proof.edges.truncate(base_edges);
        proof.vertices.truncate(base_vertices);
        // Invalid inverse references must be rejected before topology::check indexes them.
        if proof
            .edges
            .iter()
            .any(|e| e.start >= base_vertices || e.end >= base_vertices)
            || proof
                .faces
                .iter()
                .flat_map(|f| f.loops.iter().flatten())
                .any(|v| v.0 >= base_edges)
        {
            return Ok(false);
        }
        let (plan, _) = admission(&proof, tool, tol)?;
        if plan.near_cone.is_some() != near_bevel {
            return Ok(false);
        }
        // Direct single-domain construction used as the certificate equation,
        // not a retry of a Boolean against the supplied body.
        let certified = cut(&proof, tool, tol)?;
        Ok(crate::neutral::write_body(&certified) == crate::neutral::write_body(body))
    }
    prove(body, tool, tol).unwrap_or(false)
}

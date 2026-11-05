//! One annular-sector window through both radii of a plain pipe. The complete
//! operand boundaries are admitted before construction. This domain has no
//! history: its periodic wall seams are placed in the complement of the window.
//! No previous feature is discarded; a body with any previous opening refuses.
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::revolve::{Kind, Revolution, Segment};
use crate::KernelError;
use std::f64::consts::TAU;
type Result<T> = std::result::Result<T, KernelError>;
fn no(s: &str) -> KernelError {
    KernelError::Unsupported(format!("annular pipe window: {s}"))
}
fn parallel(a: Vec3, b: Vec3) -> bool {
    a.cross(b).length() == 0.
}

pub(super) fn selects(body: &Body, tool: &Body) -> bool {
    let axes: Vec<_> = tool
        .faces
        .iter()
        .filter_map(|f| match f.surface {
            Surface::Cylinder { p, a, .. } => Some((p, a)),
            _ => None,
        })
        .collect();
    axes.len() == 2
        && parallel(axes[0].1, axes[1].1)
        && parallel(axes[0].0.sub(axes[1].0), axes[0].1)
        && tool
            .edges
            .iter()
            .any(|e| matches!(e.curve,Curve::Arc{t0,t1,..} if t1>t0 && t1-t0<TAU))
        && body
            .faces
            .iter()
            .any(|f| matches!(f.surface, Surface::Cylinder { a, .. } if parallel(a, axes[0].1)))
}

/// Two coaxial walls, outward at the larger radius, inward at the smaller one.
fn walls(b: &Body, c: Vec3, a: Vec3) -> Result<[(usize, f64); 2]> {
    let mut w = Vec::new();
    for (i, f) in b.faces.iter().enumerate() {
        match f.surface {
            Surface::Cylinder { p, a: d, r, .. } => {
                if !parallel(a, d) || !parallel(p.sub(c), a) {
                    return Err(no("wall axes are not exactly coaxial"));
                }
                w.push((i, r));
            }
            Surface::Plane { .. } => {}
            _ => return Err(no("only plane and cylinder supports are admitted")),
        }
    }
    if w.len() != 2 {
        return Err(no("requires exactly two cylindrical walls"));
    }
    w.sort_by(|p, q| p.1.total_cmp(&q.1));
    if !(w[0].1 > 0. && w[1].1 > w[0].1) || b.faces[w[0].0].sense || !b.faces[w[1].0].sense {
        return Err(no("wall radii or outward senses do not bound an annulus"));
    }
    Ok([w[0], w[1]])
}

struct Sector {
    c: Vec3,
    a: Vec3,
    x: Vec3,
    t0: f64,
    t1: f64,
    ri: f64,
    ro: f64,
    low: f64,
    high: f64,
}

fn sector(b: &Body, tol: f64) -> Result<Sector> {
    if (b.faces.len(), b.edges.len(), b.vertices.len()) != (6, 12, 8) {
        return Err(no("cutter is not one six-face annular sector"));
    }
    let Some((c, a, x, t0, t1)) = b.edges.iter().find_map(|e| match e.curve {
        Curve::Arc {
            c, n, x, t0, t1, ..
        } => Some((c, n, x, t0, t1)),
        _ => None,
    }) else {
        return Err(no("sector has no circular boundary"));
    };
    if !(t1 > t0 && t1 - t0 < TAU) {
        return Err(no(
            "sector angle must be strictly between zero and a full turn",
        ));
    }
    let w = walls(b, c, a)?;
    let (ri, ro) = (w[0].1, w[1].1);
    if ri * (t1 - t0).min(TAU - (t1 - t0)) <= 4. * tol {
        return Err(no(
            "angular opening or its complement is below the clearance budget",
        ));
    }
    let mut arcs = Vec::new();
    for (i, e) in b.edges.iter().enumerate() {
        match e.curve {
            Curve::Arc {
                c: q,
                n,
                x: y,
                r,
                t0: u,
                t1: v,
            } => {
                if n != a
                    || y != x
                    || u != t0
                    || v != t1
                    || !parallel(q.sub(c), a)
                    || (r != ri && r != ro)
                {
                    return Err(no(
                        "all four sector arcs must share one exact angular recipe",
                    ));
                }
                arcs.push((i, if r == ri { 0 } else { 1 }, q.sub(c).dot(a)));
            }
            Curve::Line { .. } => {}
            _ => return Err(no("sector boundaries must be lines and circular arcs")),
        }
    }
    if arcs.len() != 4 {
        return Err(no("sector needs four circular edges"));
    }
    let low = arcs.iter().map(|r| r.2).fold(f64::INFINITY, f64::min);
    let high = arcs.iter().map(|r| r.2).fold(f64::NEG_INFINITY, f64::max);
    if high - low <= tol {
        return Err(no("sector height is below tolerance"));
    }
    // Every vertex is one (inner/outer, low/high, start/end) corner, uniquely.
    let mut tags = vec![None; 8];
    for &(e, r, z) in &arcs {
        if z != low && z != high {
            return Err(no("sector arcs do not lie on two axial levels"));
        }
        let z = if z == low { 0 } else { 1 };
        for (v, t) in [(b.edges[e].start, 0), (b.edges[e].end, 1)] {
            let tag = (r, z, t);
            if tags[v].is_some() || tags.contains(&Some(tag)) {
                return Err(no("sector corners are duplicated"));
            }
            tags[v] = Some(tag);
        }
    }
    let mut axial = 0;
    let mut radial = 0;
    for e in &b.edges {
        if !matches!(e.curve, Curve::Line { .. }) {
            continue;
        }
        let (Some((r, z, t)), Some((s, h, u))) = (tags[e.start], tags[e.end]) else {
            return Err(no("unclassified sector corner"));
        };
        if t != u {
            return Err(no("sector straight edge crosses its angular interval"));
        }
        if r == s && z != h {
            axial += 1;
        } else if z == h && r != s {
            radial += 1;
        } else {
            return Err(no("sector has a diagonal or collapsed edge"));
        }
    }
    if (axial, radial) != (4, 4) {
        return Err(no("sector needs four radial and four axial sides"));
    }
    let mut end_caps = 0;
    let mut radial_caps = 0;
    for (fi, f) in b.faces.iter().enumerate() {
        let [lp] = &f.loops[..] else {
            return Err(no("sector face has additional loops"));
        };
        if lp.len() != 4 {
            return Err(no("sector face is not one four-edge patch"));
        }
        let corners: Vec<_> = lp
            .iter()
            .map(|(e, _)| tags[b.edges[*e].start].unwrap())
            .chain(lp.iter().map(|(e, _)| tags[b.edges[*e].end].unwrap()))
            .collect();
        match f.surface {
            Surface::Cylinder { .. } => {
                let r = if fi == w[0].0 { 0 } else { 1 };
                if !corners.iter().all(|v| v.0 == r) {
                    return Err(no("sector wall mixes radii"));
                }
                for &(ei, forward) in lp {
                    if let Some(&(_, _, z)) = arcs.iter().find(|v| v.0 == ei) {
                        if forward != ((z == low) == f.sense) {
                            return Err(no("sector cylindrical loop is inside out"));
                        }
                    }
                }
            }
            Surface::Plane { p, n, .. } => {
                if parallel(n, a) && corners.iter().all(|v| v.1 == corners[0].1) {
                    let z = if corners[0].1 == 0 { low } else { high };
                    if (n != a && n != a.scale(-1.)) || p.sub(c.add(a.scale(z))).dot(a) != 0. {
                        return Err(no("sector axial cap differs from its exact ring plane"));
                    }
                    end_caps += 1;
                } else if corners.iter().all(|v| v.2 == corners[0].2) {
                    let t = if corners[0].2 == 0 { t0 } else { t1 };
                    let normal = a.cross(x.scale(t.cos()).add(a.cross(x).scale(t.sin())));
                    if (n != normal && n != normal.scale(-1.)) || p.sub(c).dot(n) != 0. {
                        return Err(no(
                            "sector radial cap differs from its exact meridian plane",
                        ));
                    }
                    radial_caps += 1;
                } else {
                    return Err(no("sector plane is not an axial or radial cap"));
                }
            }
            _ => unreachable!(),
        }
    }
    if (end_caps, radial_caps) != (2, 2) {
        return Err(no("sector cap arrangement is incomplete"));
    }
    Ok(Sector {
        c,
        a,
        x,
        t0,
        t1,
        ri,
        ro,
        low,
        high,
    })
}

fn plain_pipe(b: &Body, s: &Sector, tol: f64) -> Result<(f64, f64, f64, f64)> {
    if (b.faces.len(), b.edges.len(), b.vertices.len()) != (4, 6, 4) {
        return Err(no("stock must be a plain pipe without earlier openings"));
    }
    let w = walls(b, s.c, s.a)?;
    let (ri, ro) = (w[0].1, w[1].1);
    if ro - ri <= tol {
        return Err(no("pipe wall thickness is below tolerance"));
    }
    let mut levels = Vec::new();
    for &(fi, r) in &w {
        let [lp] = &b.faces[fi].loops[..] else {
            return Err(no("pipe wall has additional loops"));
        };
        if lp.len() != 4 {
            return Err(no("pipe wall is not one complete band"));
        }
        let mut rings = Vec::new();
        let mut lines = Vec::new();
        for &(ei, fw) in lp {
            let e = &b.edges[ei];
            match e.curve {
                Curve::Arc {
                    c,
                    n,
                    r: ar,
                    t0,
                    t1,
                    ..
                } if parallel(n, s.a)
                    && parallel(c.sub(s.c), s.a)
                    && ar == r
                    && t1 - t0 == TAU
                    && e.start == e.end =>
                {
                    rings.push((ei, c.sub(s.c).dot(s.a), fw, n.dot(s.a)))
                }
                Curve::Line { .. } => lines.push((ei, fw)),
                _ => return Err(no("pipe rim is not a complete coaxial circle")),
            }
        }
        if rings.len() != 2
            || lines.len() != 2
            || lines[0].0 != lines[1].0
            || lines[0].1 == lines[1].1
        {
            return Err(no("pipe band lacks two rings and one repeated seam"));
        }
        rings.sort_by(|p, q| p.1.total_cmp(&q.1));
        if rings[1].1 - rings[0].1 <= tol {
            return Err(no("pipe height is below tolerance"));
        }
        for (k, &(_, _, fw, dir)) in rings.iter().enumerate() {
            let positive = fw == (dir > 0.);
            if positive != ((k == 0) == b.faces[fi].sense) {
                return Err(no("pipe cylindrical loop is inside out"));
            }
        }
        levels.extend(rings.iter().map(|p| p.1));
    }
    if levels[0] != levels[2] || levels[1] != levels[3] {
        return Err(no("pipe ends do not share exact axial levels"));
    }
    for (fi, f) in b.faces.iter().enumerate() {
        if w.iter().any(|v| v.0 == fi) {
            continue;
        }
        let Surface::Plane { p, n, .. } = f.surface else {
            return Err(no("pipe end is not planar"));
        };
        if !parallel(n, s.a) || f.loops.len() != 2 || f.loops.iter().any(|lp| lp.len() != 1) {
            return Err(no("pipe end is not one complete annulus"));
        }
        let mut radii = Vec::new();
        for lp in &f.loops {
            let Curve::Arc { c, r, .. } = b.edges[lp[0].0].curve else {
                return Err(no("pipe cap has a noncircular edge"));
            };
            if p.sub(c).dot(n) != 0. {
                return Err(no("pipe cap differs from its exact ring plane"));
            }
            radii.push(r);
        }
        radii.sort_by(f64::total_cmp);
        if radii != vec![ri, ro] {
            return Err(no("pipe cap does not join the two walls"));
        }
    }
    Ok((ri, ro, levels[0], levels[1]))
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    crate::topology::check(body, tol)?;
    crate::topology::check(tool, tol)?;
    super::radial::planar_faces_outward(body)?;
    super::radial::planar_faces_outward(tool)?;
    let s = sector(tool, tol)?;
    let (ri, ro, low, high) = plain_pipe(body, &s, tol)?;
    if !(s.ri < ri - tol && s.ro > ro + tol) {
        return Err(no("cutter must extend beyond both wall radii"));
    }
    if !(s.low > low + tol && s.high < high - tol) {
        return Err(no("window must be strictly clear of both pipe ends"));
    }
    let phi = s.t1 - s.t0;
    let gap = (TAU - phi) / 2.;
    let turn = s.t0 - gap;
    let x = s.x.scale(turn.cos()).add(s.a.cross(s.x).scale(turn.sin()));
    let mut out = crate::revolve::revolve_general(
        &Revolution {
            p: s.c,
            a: s.a,
            x,
            angle: TAU,
            loops: vec![[(ri, low), (ro, low), (ro, high), (ri, high)]
                .into_iter()
                .map(|(r, z)| Segment {
                    r,
                    z,
                    kind: Kind::Line,
                })
                .collect()],
        },
        tol,
    )?;
    // The canonical seam lies halfway round the complement; the hole is wholly
    // inside one unwrapped cylindrical chart, including for sweeps above pi.
    let angles = [gap, gap + phi];
    let radial = |t: f64| x.scale(t.cos()).add(s.a.cross(x).scale(t.sin()));
    let mut v = [[[0usize; 2]; 2]; 2];
    for (r, rr) in [ri, ro].into_iter().enumerate() {
        for (z, zz) in [s.low, s.high].into_iter().enumerate() {
            for (t, tt) in angles.into_iter().enumerate() {
                v[r][z][t] = out.vertices.len();
                out.vertices
                    .push(s.c.add(s.a.scale(zz)).add(radial(tt).scale(rr)));
            }
        }
    }
    let mut arc = [[0; 2]; 2];
    let mut axis = [[0; 2]; 2];
    let mut ray = [[0; 2]; 2];
    let add = |out: &mut Body, curve: Curve, start: usize, end: usize| {
        let i = out.edges.len();
        out.edges.push(Edge { curve, start, end });
        i
    };
    for (r, rr) in [ri, ro].into_iter().enumerate() {
        for (z, zz) in [s.low, s.high].into_iter().enumerate() {
            arc[r][z] = add(
                &mut out,
                Curve::Arc {
                    c: s.c.add(s.a.scale(zz)),
                    n: s.a,
                    x,
                    r: rr,
                    t0: angles[0],
                    t1: angles[1],
                },
                v[r][z][0],
                v[r][z][1],
            );
        }
    }
    for r in 0..2 {
        for t in 0..2 {
            let (start, end) = (v[r][0][t], v[r][1][t]);
            let curve = Curve::Line {
                a: out.vertices[start],
                b: out.vertices[end],
            };
            axis[r][t] = add(&mut out, curve, start, end);
        }
    }
    for z in 0..2 {
        for t in 0..2 {
            let (start, end) = (v[0][z][t], v[1][z][t]);
            let curve = Curve::Line {
                a: out.vertices[start],
                b: out.vertices[end],
            };
            ray[z][t] = add(&mut out, curve, start, end);
        }
    }
    for f in &mut out.faces {
        if let Surface::Cylinder { r, .. } = f.surface {
            f.loops.push(if r == ro {
                vec![
                    (axis[1][0], true),
                    (arc[1][1], true),
                    (axis[1][1], false),
                    (arc[1][0], false),
                ]
            } else {
                vec![
                    (arc[0][0], true),
                    (axis[0][1], true),
                    (arc[0][1], false),
                    (axis[0][0], false),
                ]
            });
        }
    }
    out.faces.push(Face {
        surface: Surface::Plane {
            p: s.c.add(s.a.scale(s.low)),
            n: s.a,
            x,
        },
        sense: true,
        loops: vec![vec![
            (ray[0][0], true),
            (arc[1][0], true),
            (ray[0][1], false),
            (arc[0][0], false),
        ]],
    });
    out.faces.push(Face {
        surface: Surface::Plane {
            p: s.c.add(s.a.scale(s.high)),
            n: s.a,
            x,
        },
        sense: false,
        loops: vec![vec![
            (ray[1][0], false),
            (arc[0][1], true),
            (ray[1][1], true),
            (arc[1][1], false),
        ]],
    });
    out.faces.push(Face {
        surface: Surface::Plane {
            p: s.c,
            n: s.a.cross(radial(angles[0])),
            x: s.a,
        },
        sense: true,
        loops: vec![vec![
            (axis[0][0], true),
            (ray[1][0], true),
            (axis[1][0], false),
            (ray[0][0], false),
        ]],
    });
    out.faces.push(Face {
        surface: Surface::Plane {
            p: s.c,
            n: s.a.cross(radial(angles[1])),
            x: s.a,
        },
        sense: false,
        loops: vec![vec![
            (axis[1][1], true),
            (ray[1][1], false),
            (axis[0][1], false),
            (ray[0][1], true),
        ]],
    });
    crate::topology::check(&out, tol)?;
    super::radial::planar_faces_outward(&out)?;
    Ok(out)
}

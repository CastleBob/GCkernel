//! Transverse cuts on stock that already carries features the section chart cannot represent (Claude, 2026-10-31).
//!
//! Bob's Pipe 1.25 18 Copes Rotated 90 (a cope at y = 0 after one turned 90 degrees at y = 18) and Pipe 1.25 24 Slot
//! And Miter (a cross slot at y = 8 after a miter at y = 0): the transverse chart (x = u stations, y = v the beam
//! axis, z = w the sweep) needs every face to be a graph y(x) or z(x), but the far feature's faces are not (a
//! cylinder along x, a plane turned in both y and z). They are far from the new cut along the beam axis, so the
//! body is split by the plane y = c between them, where it is a plain prism along y: the near piece is cut in the
//! chart and glued back to the far piece, untouched, along y = c.
//!
//! Admission, from the geometry before anything is built: every face outside the chart lies wholly on one side of
//! y = c and the whole tool (its surfaces' full y extent where they bend in y) on the other, with c in the middle
//! of the widest gap between every vertex and edge extremum there; every face that y = c crosses is prismatic
//! along y (a plane with n.y = 0 or a cylinder along y) and is crossed by exactly two straight edge uses. The
//! split adds a cap on each side; the cut runs on the near piece only; gluing removes the caps, splits the open
//! edges along y = c at each other's vertices and joins them where they are the same curve (one vertex within
//! tolerance, Bob's rule), and must close every one of them, or the cut is refused by name. The far piece keeps
//! its own faces, curves and surfaces exactly; the faces y = c crossed keep a seam there. No retry: when this
//! admission fails the cut is refused.
use super::*;

/// Can the transverse chart represent this face: a section graph y(x) or z(x), or an x-normal plane?
fn compatible(face: &Face) -> bool {
    match face.surface {
        Surface::Plane { n, .. } => n.y == 0. || n.z == 0.,
        Surface::Cylinder { a, .. } => {
            (a.y.abs() == 1. && a.x == 0. && a.z == 0.) || (a.z.abs() == 1. && a.x == 0. && a.y == 0.)
        }
        _ => false,
    }
}

/// A rigorous y range of a curve: lines by their ends, conics by their extrema, a cylinder pair by whichever of
/// its two cylinders has its axis across y (every point of it lies within that cylinder's radius of its axis).
fn curve_y(curve: &Curve) -> Option<(f64, f64)> {
    match *curve {
        Curve::Line { a, b } => Some((a.y.min(b.y), a.y.max(b.y))),
        Curve::CylinderPair { c, a, rho, p, b, r, .. } => {
            if b.y == 0. {
                Some((p.y - r, p.y + r))
            } else if a.y == 0. {
                Some((c.y - rho, c.y + rho))
            } else {
                None
            }
        }
        Curve::Bezier { p } => Some(p.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), q| (l.min(q.y), h.max(q.y)))),
        Curve::RationalBezier { p, .. } => {
            Some(p.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), q| (l.min(q.y), h.max(q.y))))
        }
        _ => {
            let (c, u, v, t0, t1) = conic(curve)?;
            let (s, e) = (curve.start().y, curve.end().y);
            let (mut lo, mut hi) = (s.min(e), s.max(e));
            let phase = v.y.atan2(u.y);
            for angle in [phase, phase + PI] {
                for t in parameters(angle, t0, t1) {
                    let y = c.y + u.y * t.cos() + v.y * t.sin();
                    lo = lo.min(y);
                    hi = hi.max(y);
                }
            }
            Some((lo, hi))
        }
    }
}

/// The y values where a curve has a vertex or a y extremum (the split plane keeps clear of them).
fn curve_events(curve: &Curve, out: &mut Vec<f64>) {
    out.push(curve.start().y);
    out.push(curve.end().y);
    if let Some((c, u, v, t0, t1)) = conic(curve) {
        let phase = v.y.atan2(u.y);
        for angle in [phase, phase + PI] {
            for t in parameters(angle, t0, t1) {
                out.push(c.y + u.y * t.cos() + v.y * t.sin());
            }
        }
    } else if let Some((lo, hi)) = curve_y(curve) {
        out.push(lo);
        out.push(hi);
    }
}

/// A face's y range: its boundary, and for a cylinder whose axis is across y its whole surface (its interior can
/// bend past its boundary there).
fn face_y(body: &Body, face: &Face) -> Result<(f64, f64)> {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for (e, _) in face.loops.iter().flatten() {
        let (l, h) = curve_y(&body.edges[*e].curve)
            .ok_or_else(|| unsupported("transverse split: an edge has no y bound"))?;
        lo = lo.min(l);
        hi = hi.max(h);
    }
    match face.surface {
        // y is linear: its extremes on a plane, or on a cylinder along y, are on the boundary.
        Surface::Plane { .. } => {}
        Surface::Cylinder { a, .. } if a.y.abs() == 1. && a.x == 0. && a.z == 0. => {}
        // Across y: the whole surface lies within r of the axis.
        Surface::Cylinder { p, a, r, .. } if a.y == 0. => {
            lo = lo.min(p.y - r);
            hi = hi.max(p.y + r);
        }
        _ => return Err(unsupported("transverse split: a face is neither a plane nor a cylinder along or across the stock")),
    }
    Ok((lo, hi))
}

/// Does the chart hold the body directly, or after a split (the faces outside it wholly beyond the tool along y)?
pub(super) fn separable(body: &Body, tool: &Body) -> bool {
    body.faces.iter().all(compatible) || reaches(body, tool).is_ok_and(|r| r.is_some())
}

/// The tool's and the out-of-chart faces' y ranges when they are apart: (gap low, gap high, tool below).
fn reaches(body: &Body, tool: &Body) -> Result<Option<(f64, f64, bool)>> {
    let (mut tlo, mut thi) = (f64::INFINITY, f64::NEG_INFINITY);
    for face in &tool.faces {
        let (l, h) = face_y(tool, face)?;
        tlo = tlo.min(l);
        thi = thi.max(h);
    }
    let (mut ilo, mut ihi) = (f64::INFINITY, f64::NEG_INFINITY);
    for face in body.faces.iter().filter(|f| !compatible(f)) {
        let (l, h) = face_y(body, face)?;
        ilo = ilo.min(l);
        ihi = ihi.max(h);
    }
    Ok(if thi < ilo {
        Some((thi, ilo, true))
    } else if ihi < tlo {
        Some((ihi, tlo, false))
    } else {
        None
    })
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64, local: fn(&Body, &Body, f64) -> Result<Body>) -> Result<Body> {
    if body.faces.iter().all(compatible) {
        return local(body, tool, tol);
    }
    let Some((gap_lo, gap_hi, near_low)) = reaches(body, tool)? else {
        return Err(unsupported(
            "transverse cut meets faces its section chart cannot represent (another feature's faces are in the cut's reach along the stock)",
        ));
    };
    // The split plane: the middle of the widest gap between every vertex and edge extremum inside (gap_lo, gap_hi).
    let mut events = vec![gap_lo, gap_hi];
    for edge in &body.edges {
        curve_events(&edge.curve, &mut events);
    }
    events.retain(|y| *y >= gap_lo && *y <= gap_hi);
    events.sort_by(f64::total_cmp);
    let (width, c) = events
        .windows(2)
        .map(|w| (w[1] - w[0], (w[0] + w[1]) * 0.5))
        .fold((0., 0.), |best, g| if g.0 > best.0 { g } else { best });
    if width <= 8. * tol {
        return Err(unsupported("transverse split: no clear plane between the cut and the stock's other features"));
    }
    let (lower, upper) = split(body, c, tol)?;
    let (near, far) = if near_low { (lower, upper) } else { (upper, lower) };
    if !near.faces.iter().all(compatible) {
        return Err(unsupported("transverse split: the near piece still has faces outside the section chart"));
    }
    let cut = local(&near, tool, tol)?;
    if cut.faces.is_empty() {
        // The near piece is wholly removed: the far piece, with its cap, is the result.
        return Ok(far);
    }
    let outward = if near_low { 1. } else { -1. };
    glue(&cut, &far, c, outward, tol)
}

/// One oriented use of an edge piece, its side of the plane (true = lower) and its ends (vertex indices).
#[derive(Clone, Copy)]
struct Use {
    edge: usize,
    forward: bool,
    lower: bool,
}

/// Splits `body` by the plane y = c into its lower (y <= c) and upper pieces, each closed by a cap at y = c.
fn split(body: &Body, c: f64, tol: f64) -> Result<(Body, Body)> {
    let mut vertices = body.vertices.clone();
    // Each original edge's pieces: (lower piece, upper piece) as indices into `edges`, None where it has none, and
    // for a crossing edge whether its start is the lower end.
    let mut edges: Vec<Edge> = Vec::new();
    let mut pieces: Vec<(Option<usize>, Option<usize>, bool)> = Vec::new();
    for edge in &body.edges {
        let (lo, hi) = curve_y(&edge.curve).ok_or_else(|| unsupported("transverse split: an edge has no y bound"))?;
        if hi < c {
            edges.push(edge.clone());
            pieces.push((Some(edges.len() - 1), None, true));
        } else if lo > c {
            edges.push(edge.clone());
            pieces.push((None, Some(edges.len() - 1), true));
        } else {
            let Curve::Line { a, b } = edge.curve else {
                return Err(unsupported("transverse split: the split plane crosses a curved edge"));
            };
            let t = (c - a.y) / (b.y - a.y);
            let q = a.add(b.sub(a).scale(t));
            let v = vertices.len();
            vertices.push(Vec3::new(q.x, c, q.z));
            let q = vertices[v];
            let first = Edge { curve: Curve::Line { a, b: q }, start: edge.start, end: v };
            let second = Edge { curve: Curve::Line { a: q, b }, start: v, end: edge.end };
            let start_low = a.y < c;
            edges.push(first);
            edges.push(second);
            let (f, s) = (edges.len() - 2, edges.len() - 1);
            pieces.push(if start_low { (Some(f), Some(s), true) } else { (Some(s), Some(f), false) });
        }
    }
    let crossing = |e: usize| pieces[e].0.is_some() && pieces[e].1.is_some();
    let mut lower_faces: Vec<Face> = Vec::new();
    let mut upper_faces: Vec<Face> = Vec::new();
    // Section edges (index in `edges`), as the lower side faces use them forward.
    let mut sections: Vec<usize> = Vec::new();
    for face in &body.faces {
        let crossings: usize = face.loops.iter().flatten().filter(|(e, _)| crossing(*e)).count();
        if crossings == 0 {
            // Wholly one side: decided by any of its edges.
            let (e, _) = face.loops[0][0];
            let low = pieces[e].1.is_none();
            let map = |(e, f): &(usize, bool)| (if low { pieces[*e].0.unwrap() } else { pieces[*e].1.unwrap() }, *f);
            let moved = Face {
                surface: face.surface.clone(),
                sense: face.sense,
                loops: face.loops.iter().map(|l| l.iter().map(map).collect()).collect(),
            };
            if low { lower_faces.push(moved) } else { upper_faces.push(moved) }
            continue;
        }
        let prismatic = match face.surface {
            Surface::Plane { n, .. } => n.y == 0.,
            Surface::Cylinder { a, .. } => a.y.abs() == 1. && a.x == 0. && a.z == 0.,
            _ => false,
        };
        if !prismatic || crossings != 2 {
            return Err(unsupported(
                "transverse split: a face the split plane crosses is not a plain prism face crossed twice",
            ));
        }
        let mut lower_loops = Vec::new();
        let mut upper_loops = Vec::new();
        for lp in &face.loops {
            if !lp.iter().any(|(e, _)| crossing(*e)) {
                let (e, _) = lp[0];
                let low = pieces[e].1.is_none();
                let mapped: Vec<(usize, bool)> = lp
                    .iter()
                    .map(|(e, f)| (if low { pieces[*e].0.unwrap() } else { pieces[*e].1.unwrap() }, *f))
                    .collect();
                if low { lower_loops.push(mapped) } else { upper_loops.push(mapped) }
                continue;
            }
            // The loop as pieces in traversal order.
            let mut uses: Vec<Use> = Vec::new();
            for &(e, forward) in lp {
                let (l, u, start_low) = pieces[e];
                match (l, u) {
                    (Some(l), None) => uses.push(Use { edge: l, forward, lower: true }),
                    (None, Some(u)) => uses.push(Use { edge: u, forward, lower: false }),
                    (Some(l), Some(u)) => {
                        // Forward, the start's piece comes first.
                        let (a, b) = if start_low { ((l, true), (u, false)) } else { ((u, false), (l, true)) };
                        let order = if forward { [a, b] } else { [b, a] };
                        for (edge, lower) in order {
                            uses.push(Use { edge, forward, lower });
                        }
                    }
                    (None, None) => unreachable!(),
                }
            }
            // Rotate so it starts with the first lower piece after an upper one.
            let n = uses.len();
            let start = (0..n)
                .find(|i| uses[*i].lower && !uses[(i + n - 1) % n].lower)
                .ok_or_else(|| invalid("transverse split: a crossed loop has no lower run"))?;
            uses.rotate_left(start);
            let k = uses.iter().position(|u| !u.lower).unwrap();
            if uses[k..].iter().any(|u| u.lower) {
                return Err(unsupported("transverse split: a crossed loop has more than one lower run"));
            }
            let end_of = |u: &Use| if u.forward { edges[u.edge].end } else { edges[u.edge].start };
            let start_of = |u: &Use| if u.forward { edges[u.edge].start } else { edges[u.edge].end };
            let (p_in, p_out) = (start_of(&uses[0]), end_of(&uses[k - 1]));
            // The section edge from p_out to p_in with the lower material on its left (about the outward normal).
            let (a, b) = (vertices[p_out], vertices[p_in]);
            let curve = match face.surface {
                Surface::Plane { n, .. } => {
                    if p_in == p_out {
                        return Err(invalid("transverse split: a planar face is crossed at one point"));
                    }
                    let outward = n.scale(sign(face.sense));
                    if outward.cross(b.sub(a)).y >= 0. {
                        return Err(invalid("transverse split: a crossed planar face is inconsistently oriented"));
                    }
                    Curve::Line { a, b }
                }
                Surface::Cylinder { p, r, .. } => {
                    let center = Vec3::new(p.x, c, p.z);
                    let s = sign(face.sense);
                    let axis = Vec3::new(0., -s, 0.);
                    let (da, db) = (a.sub(center), b.sub(center));
                    if (da.length() - r).abs() > tol || (db.length() - r).abs() > tol {
                        return Err(invalid("transverse split: a crossing point is off its cylinder"));
                    }
                    let x = da.scale(1. / da.length());
                    let yv = axis.cross(x);
                    let mut sweep = db.dot(yv).atan2(db.dot(x));
                    if p_in == p_out {
                        sweep = TAU;
                    } else if sweep <= 0. {
                        sweep += TAU;
                    }
                    Curve::Arc { c: center, n: axis, x, r, t0: 0., t1: sweep }
                }
                _ => unreachable!(),
            };
            edges.push(Edge { curve, start: p_out, end: p_in });
            let section = edges.len() - 1;
            sections.push(section);
            let mut lower_loop: Vec<(usize, bool)> = uses[..k].iter().map(|u| (u.edge, u.forward)).collect();
            lower_loop.push((section, true));
            let mut upper_loop: Vec<(usize, bool)> = uses[k..].iter().map(|u| (u.edge, u.forward)).collect();
            upper_loop.push((section, false));
            lower_loops.push(lower_loop);
            upper_loops.push(upper_loop);
        }
        // Keep the outer loop first: the crossed loop is the outer one of each piece.
        lower_loops.rotate_right(1);
        upper_loops.rotate_right(1);
        lower_faces.push(Face { surface: face.surface.clone(), sense: face.sense, loops: lower_loops });
        upper_faces.push(Face { surface: face.surface.clone(), sense: face.sense, loops: upper_loops });
    }
    // The caps: the lower one (outward +y) runs every section edge backward, the upper one forward.
    let chain = |forward: bool| -> Result<Vec<Vec<(usize, bool)>>> {
        let ends = |e: usize| if forward { (edges[e].start, edges[e].end) } else { (edges[e].end, edges[e].start) };
        let mut left: Vec<usize> = sections.clone();
        let mut loops = Vec::new();
        while let Some(first) = left.pop() {
            let mut lp = vec![(first, forward)];
            let (head, mut at) = ends(first);
            while at != head {
                let next: Vec<usize> = left.iter().copied().filter(|e| ends(*e).0 == at).collect();
                let [e] = next[..] else {
                    return Err(unsupported("transverse split: the section at the split plane does not chain into loops"));
                };
                left.retain(|x| *x != e);
                lp.push((e, forward));
                at = ends(e).1;
            }
            loops.push(lp);
        }
        Ok(loops)
    };
    let cap = |forward: bool, outward: f64| -> Result<Face> {
        let mut loops = chain(forward)?;
        let normal = Vec3::new(0., outward, 0.);
        let area = |lp: &Vec<(usize, bool)>| {
            let curves: Vec<Curve> = lp
                .iter()
                .map(|(e, f)| if *f { edges[*e].curve.clone() } else { edges[*e].curve.reversed() })
                .collect();
            crate::extrude::signed_area(&curves, normal)
        };
        if loops.iter().filter(|l| area(l) > 0.).count() != 1 {
            return Err(unsupported("transverse split: the section at the split plane is not one region"));
        }
        loops.sort_by(|a, b| area(b).total_cmp(&area(a)));
        Ok(Face {
            surface: Surface::Plane { p: Vec3::new(0., c, 0.), n: Vec3::new(0., 1., 0.), x: Vec3::new(0., 0., 1.) },
            sense: outward > 0.,
            loops,
        })
    };
    lower_faces.push(cap(false, 1.)?);
    upper_faces.push(cap(true, -1.)?);
    let lower = compact(&vertices, &edges, lower_faces);
    let upper = compact(&vertices, &edges, upper_faces);
    crate::topology::check(&lower, tol).map_err(|e| KernelError::Internal(format!("transverse split lower piece: {e}")))?;
    crate::topology::check(&upper, tol).map_err(|e| KernelError::Internal(format!("transverse split upper piece: {e}")))?;
    Ok((lower, upper))
}

/// The body of `faces`, keeping only the edges and vertices they use.
fn compact(vertices: &[Vec3], edges: &[Edge], faces: Vec<Face>) -> Body {
    let mut vmap = vec![usize::MAX; vertices.len()];
    let mut emap = vec![usize::MAX; edges.len()];
    let mut body = Body::default();
    for face in &faces {
        for (e, _) in face.loops.iter().flatten() {
            if emap[*e] != usize::MAX {
                continue;
            }
            let mut edge = edges[*e].clone();
            for v in [&mut edge.start, &mut edge.end] {
                if vmap[*v] == usize::MAX {
                    vmap[*v] = body.vertices.len();
                    body.vertices.push(vertices[*v]);
                }
                *v = vmap[*v];
            }
            emap[*e] = body.edges.len();
            body.edges.push(edge);
        }
    }
    for mut face in faces {
        for lp in &mut face.loops {
            for (e, _) in lp.iter_mut() {
                *e = emap[*e];
            }
        }
        body.faces.push(face);
    }
    body
}

/// Joins the cut near piece and the untouched far piece along y = c: both caps go, the open edges along the plane
/// are split at each other's vertices and joined where they are the same curve.
fn glue(near: &Body, far: &Body, c: f64, outward: f64, tol: f64) -> Result<Body> {
    let is_cap = |f: &Face, out: f64| {
        matches!(f.surface, Surface::Plane { p, n, .. }
            if n.x == 0. && n.z == 0. && (p.y - c).abs() <= tol && n.y * sign(f.sense) * out > 0.)
    };
    let mut vertices = near.vertices.clone();
    let mut edges = near.edges.clone();
    let mut faces: Vec<Face> = near.faces.iter().filter(|f| !is_cap(f, outward)).cloned().collect();
    let (nv, ne) = (vertices.len(), edges.len());
    vertices.extend(far.vertices.iter().copied());
    edges.extend(far.edges.iter().map(|e| Edge { curve: e.curve.clone(), start: e.start + nv, end: e.end + nv }));
    let far_caps = far.faces.iter().filter(|f| is_cap(f, -outward)).count();
    if far_caps != 1 {
        return Err(internal("transverse split: the far piece's cap is not one face"));
    }
    for f in far.faces.iter().filter(|f| !is_cap(f, -outward)) {
        let mut f = f.clone();
        for lp in &mut f.loops {
            for (e, _) in lp.iter_mut() {
                *e += ne;
            }
        }
        faces.push(f);
    }
    // Open edges: used once.
    let mut uses = vec![0usize; edges.len()];
    for (e, _) in faces.iter().flat_map(|f| f.loops.iter().flatten()) {
        uses[*e] += 1;
    }
    let open: Vec<usize> = (0..edges.len()).filter(|e| uses[*e] == 1).collect();
    for e in &open {
        let curve = &edges[*e].curve;
        let (lo, hi) = curve_y(curve).ok_or_else(|| internal("transverse split: an open edge has no y bound"))?;
        if lo < c - tol || hi > c + tol {
            return Err(unsupported("transverse split: the cut reached the split plane"));
        }
    }
    // One vertex where two agree within tolerance (Bob's rule); the near piece's position is kept.
    let mut vmap: Vec<usize> = (0..vertices.len()).collect();
    let open_vertices: Vec<usize> = {
        let mut v: Vec<usize> = open.iter().flat_map(|e| [edges[*e].start, edges[*e].end]).collect();
        v.sort();
        v.dedup();
        v
    };
    for &v in open_vertices.iter().filter(|v| **v >= nv) {
        let matches: Vec<usize> = open_vertices
            .iter()
            .copied()
            .filter(|w| *w < nv && vertices[*w].distance(vertices[v]) <= tol)
            .collect();
        match matches[..] {
            [] => {}
            [w] => vmap[v] = w,
            _ => return Err(unsupported("transverse split: ambiguous vertices along the split plane")),
        }
    }
    for e in &mut edges {
        e.start = vmap[e.start];
        e.end = vmap[e.end];
    }
    let mut points: Vec<usize> = open_vertices.iter().map(|v| vmap[*v]).collect();
    points.sort();
    points.dedup();
    // Split every open edge at the other open vertices lying on it; each piece replaces the edge in its loop.
    let mut replace: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for &e in &open {
        let curve = edges[e].curve.clone();
        let (t0, t1) = domain(&curve);
        let mut cuts: Vec<(f64, usize)> = vec![(t0, edges[e].start), (t1, edges[e].end)];
        for &v in &points {
            if v == edges[e].start || v == edges[e].end {
                continue;
            }
            if let Some(t) = parameter(&curve, vertices[v], tol) {
                if t > t0 && t < t1 {
                    cuts.push((t, v));
                }
            }
        }
        if cuts.len() == 2 {
            continue;
        }
        cuts.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut list = Vec::new();
        for w in cuts.windows(2) {
            edges.push(Edge { curve: slice(&curve, w[0].0, w[1].0), start: w[0].1, end: w[1].1 });
            list.push(edges.len() - 1);
        }
        replace.insert(e, list);
    }
    for f in &mut faces {
        for lp in &mut f.loops {
            let mut out = Vec::with_capacity(lp.len());
            for &(e, forward) in lp.iter() {
                match replace.get(&e) {
                    Some(list) if forward => out.extend(list.iter().map(|p| (*p, true))),
                    Some(list) => out.extend(list.iter().rev().map(|p| (*p, false))),
                    None => out.push((e, forward)),
                }
            }
            *lp = out;
        }
    }
    // Pair the open pieces: a far piece joins the near piece with the same ends that is the same curve run the
    // other way round in its loop.
    let mut uses: std::collections::HashMap<usize, (usize, bool)> = std::collections::HashMap::new();
    let mut count = std::collections::HashMap::new();
    for (e, forward) in faces.iter().flat_map(|f| f.loops.iter().flatten()) {
        *count.entry(*e).or_insert(0usize) += 1;
        uses.insert(*e, (*e, *forward));
    }
    let open: Vec<usize> = count.iter().filter(|(_, n)| **n == 1).map(|(e, _)| *e).collect();
    let (near_open, far_open): (Vec<usize>, Vec<usize>) = open.iter().partition(|e| **e < ne || replaced_from(&replace, **e, ne));
    let mut join: std::collections::HashMap<usize, (usize, bool)> = std::collections::HashMap::new();
    for &f in &far_open {
        let (fs, fe) = (edges[f].start, edges[f].end);
        let candidates: Vec<(usize, bool)> = near_open
            .iter()
            .copied()
            .filter_map(|n| {
                let (ns, nn) = (edges[n].start, edges[n].end);
                if (ns, nn) == (fs, fe) && same_curve(&edges[f].curve, &edges[n].curve, false, tol) {
                    Some((n, true))
                } else if (ns, nn) == (fe, fs) && same_curve(&edges[f].curve, &edges[n].curve, true, tol) {
                    Some((n, false))
                } else {
                    None
                }
            })
            .collect();
        let [(n, same)] = candidates[..] else {
            return Err(unsupported("transverse split: the cut piece and the far piece do not meet edge for edge at the split plane"));
        };
        if uses[&f].1 == (uses[&n].1 == same) {
            return Err(internal("transverse split: joined edges run the same way in both pieces"));
        }
        join.insert(f, (n, same));
    }
    if join.len() != far_open.len() || join.len() != near_open.len() {
        return Err(unsupported("transverse split: the cut piece and the far piece do not meet edge for edge at the split plane"));
    }
    for f in &mut faces {
        for lp in &mut f.loops {
            for (e, forward) in lp.iter_mut() {
                if let Some((n, same)) = join.get(e) {
                    *forward = if *same { *forward } else { !*forward };
                    *e = *n;
                }
            }
        }
    }
    let body = compact(&vertices, &edges, faces);
    crate::topology::check(&body, tol).map_err(|e| KernelError::Internal(format!("transverse split glue: {e}")))?;
    Ok(body)
}

/// Is a split piece one of the near piece's edges' pieces?
fn replaced_from(replace: &std::collections::HashMap<usize, Vec<usize>>, piece: usize, ne: usize) -> bool {
    replace.iter().any(|(e, list)| *e < ne && list.contains(&piece))
}

fn internal(m: &str) -> KernelError {
    KernelError::Internal(m.into())
}

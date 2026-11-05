//! Multiple transverse material intervals for through cuts crossing cylinders.
//! Shares analytic boundary construction with the end-cut domain. Input sections
//! are read from the current B-rep; no prior stock or operation recipe is used.
use super::*;

fn positive(lower: Graph, upper: Graph, lo: f64, hi: f64, tol: f64) -> bool {
    // At a computed intersection station, evaluating two independently stored
    // graphs can leave a few ulps of endpoint residual. Bound only that arithmetic
    // residual, using the terms before cancellation; do not use the modelling
    // tolerance to erase an interior crossing or a narrow material interval.
    let evaluation_size = |g: Graph, x: f64| match g {
        Graph::Line { m, b } => (m * x).abs() + b.abs(),
        Graph::Circle { cv, r, .. } => cv.abs() + r,
    };
    [lo, hi].into_iter().all(|x| {
        let roundoff =
            (16. * f64::EPSILON * (evaluation_size(lower, x) + evaluation_size(upper, x))).min(tol);
        upper.value(x) - lower.value(x) >= -roundoff
    }) && upper.value((lo + hi) * 0.5) - lower.value((lo + hi) * 0.5) > tol
        && difference_stations(lower, upper, lo, hi)
            .into_iter()
            .filter(|x| *x > lo && *x < hi)
            .all(|x| upper.value(x) - lower.value(x) > tol)
}

// Intersections of the actual supporting graphs partition intervals in which
// their ordering is constant. Squared equations supply candidates only: both
// circle branches are checked before accepting a station. The downstream
// whole-band ordering proof remains mandatory.
fn intersections(a: Graph, b: Graph, lo: f64, hi: f64, tol: f64) -> Result<Vec<f64>> {
    let mut candidates = Vec::new();
    match (a, b) {
        (Graph::Line { m: am, b: ab }, Graph::Line { m: bm, b: bb }) => {
            if am != bm {
                candidates.push(((bb - ab) / (am - bm), false));
            }
        }
        (Graph::Circle { .. }, Graph::Line { .. }) => return intersections(b, a, lo, hi, tol),
        (Graph::Line { m, b }, Graph::Circle { cx, cv, r, .. }) => {
            let d = m * cx + b - cv;
            let aa = 1. + m * m;
            let discriminant = aa * r * r - d * d;
            if discriminant >= 0. {
                let root = discriminant.sqrt();
                // Stable quadratic roots in the coordinate relative to cx.
                let q = -m * d - root.copysign(m * d);
                if q == 0. {
                    candidates.push((cx, true));
                } else {
                    candidates.push((cx + q / aa, discriminant == 0.));
                    candidates.push((cx + (d - r) * (d + r) / q, discriminant == 0.));
                }
            }
        }
        (
            Graph::Circle {
                cx: ax,
                cv: ay,
                r: ar,
                ..
            },
            Graph::Circle {
                cx: bx,
                cv: by,
                r: br,
                ..
            },
        ) => {
            let (dx, dy) = (bx - ax, by - ay);
            let d = dx.hypot(dy);
            if d != 0. && d <= ar + br && d >= (ar - br).abs() {
                let along = ((ar - br) * (ar + br) + d * d) / (2. * d);
                let height2 = (ar - along) * (ar + along);
                if height2 < 0. {
                    return Err(unsupported("unresolved circle intersection discriminant"));
                }
                let base_x = ax + along * dx / d;
                let offset = height2.sqrt() * dy / d;
                candidates.push((base_x - offset, height2 == 0.));
                candidates.push((base_x + offset, height2 == 0.));
            }
            // Concentric equal circles either coincide or meet at their X
            // extrema, which are already input events; unequal ones never meet.
        }
    }
    let mut result = Vec::new();
    for (x, tangent) in candidates {
        if !x.is_finite() {
            return Err(unsupported("nonfinite transverse intersection"));
        }
        if x <= lo || x >= hi {
            continue;
        }
        if (a.value(x) - b.value(x)).abs() > tol {
            continue; // Candidate belongs to the other circle branch.
        }
        if tangent {
            return Err(unsupported(
                "tangent transverse graphs inside an active band",
            ));
        }
        result.push(x);
    }
    Ok(result)
}

fn intersection_stations(
    body: &Body,
    tool: &Body,
    xs: &[f64],
    extrema: &[f64],
    tol: f64,
    work: &mut usize,
) -> Result<Vec<f64>> {
    let mut crossings = Vec::new();
    for band in xs.windows(2) {
        let mid = band[0] + (band[1] - band[0]) * 0.5;
        let mut hits = boundary_hits(body, mid, tol)?;
        hits.extend(boundary_hits(tool, mid, tol)?);
        for is_z in [false, true] {
            let mut graphs = Vec::new();
            for hit in hits.iter().filter(|h| (h.role >= 2) == is_z) {
                spend(work, graphs.len() + 1)?;
                if !graphs.contains(&hit.graph) {
                    graphs.push(hit.graph);
                }
            }
            for i in 0..graphs.len() {
                for j in i + 1..graphs.len() {
                    spend(work, 1)?;
                    crossings.extend(intersections(graphs[i], graphs[j], band[0], band[1], tol)?);
                }
            }
        }
    }
    let mut refined = stations(xs.iter().chain(&crossings).copied().collect(), tol)?;
    // A station within tolerance of an analytic crossing IS that crossing (Bob
    // 2026-10-28: what agrees within the tolerance is the same). Merging keeps
    // the smaller value, which can be an authored vertex a few ulps short of the
    // crossing; the band would then start past it and fail its ordering proof.
    // Take the crossing, as extremal stations take the extremum. Extremal
    // stations stay put; crossings that disagree beyond roundoff are left to
    // the ordering proof, which refuses them.
    for x in refined.iter_mut() {
        if extrema.contains(x) {
            continue;
        }
        let near: Vec<f64> = crossings.iter().copied().filter(|c| (c - *x).abs() <= tol).collect();
        let (Some(lo), Some(hi)) = (
            near.iter().copied().reduce(f64::min),
            near.iter().copied().reduce(f64::max),
        ) else {
            continue;
        };
        if hi - lo <= 64. * f64::EPSILON * (1. + x.abs()) {
            *x = lo;
        }
    }
    if refined.windows(2).any(|w| w[1] - w[0] <= tol) {
        return Err(unsupported(
            "ambiguous transverse crossing stations within tolerance",
        ));
    }
    Ok(refined)
}

/// Two stored graphs of ONE supporting surface (Claude 2026-10-31, Bob's aa5: an NC1 angle preparation is built on
/// the adjoining leg's exact miter plane, and its tool face and the earlier cut's face are that plane stored through
/// different points and normals, 1e-13 apart): they are the same graph when they agree within evaluation roundoff,
/// 64 machine epsilons of the terms before cancellation, at both ends of the band (two lines: their difference is
/// linear, so the ends bound it) or have the same circle to that roundoff. The first is kept. Distinct graphs closer
/// than the tolerance but beyond roundoff are still refused.
fn same_graph(a: Graph, b: Graph, lo: f64, hi: f64) -> bool {
    let size = |g: Graph, x: f64| match g {
        Graph::Line { m, b } => (m * x).abs() + b.abs(),
        Graph::Circle { cv, r, .. } => cv.abs() + r,
    };
    match (a, b) {
        (Graph::Line { .. }, Graph::Line { .. }) => [lo, hi]
            .into_iter()
            .all(|x| (a.value(x) - b.value(x)).abs() <= 64. * f64::EPSILON * (size(a, x) + size(b, x))),
        (Graph::Circle { cx, cv, r, branch }, Graph::Circle { cx: dx, cv: dv, r: dr, branch: db }) => {
            let ulps = 64. * f64::EPSILON * (cx.abs() + dx.abs() + cv.abs() + dv.abs() + r + dr);
            branch == db && (cx - dx).abs() <= ulps && (cv - dv).abs() <= ulps && (r - dr).abs() <= ulps
        }
        _ => false,
    }
}

fn ordered_graphs(mut graphs: Vec<Graph>, lo: f64, hi: f64, tol: f64) -> Result<Vec<Graph>> {
    let mid = (lo + hi) * 0.5;
    graphs.sort_by(|a, b| a.value(mid).total_cmp(&b.value(mid)));
    let mut out: Vec<Graph> = Vec::new();
    for g in graphs {
        if let Some(&last) = out.last() {
            if (g.value(mid) - last.value(mid)).abs() <= tol {
                if g == last || same_graph(last, g, lo, hi) {
                    continue;
                }
                return Err(unsupported("transverse graph coincidence is unresolved"));
            }
            if !positive(last, g, lo, hi, tol) {
                return Err(unsupported(
                    "transverse graphs cross or touch inside a band",
                ));
            }
        }
        out.push(g);
    }
    if out.len() > 129 {
        return Err(unsupported(
            "transverse arrangement has too many coordinates",
        ));
    }
    Ok(out)
}

#[cfg(test)]
fn sections(body: &Body, lo: f64, hi: f64, tol: f64) -> Result<Vec<Band>> {
    sections_with_budget(body, lo, hi, tol, &mut 0)
}
fn spend(work: &mut usize, amount: usize) -> Result<()> {
    *work = work.saturating_add(amount);
    if *work > 2_000_000 {
        return Err(unsupported(
            "transverse arrangement exceeds its work budget",
        ));
    }
    Ok(())
}
fn boundary_hits(body: &Body, mid: f64, tol: f64) -> Result<Vec<Hit>> {
    let mut hits = Vec::new();
    for face in &body.faces {
        if matches!(face.surface, Surface::Plane { n, .. } if n.y==0. && n.z==0.) {
            continue;
        }
        let mut points = Vec::new();
        for (i, forward) in face.loops.iter().flatten() {
            let curve = &body.edges[*i].curve;
            for q in roots(curve, mid, tol) {
                let t = parameter(curve, q, tol)
                    .ok_or_else(|| unsupported("cannot locate oriented section crossing"))?;
                let tangent = curve.eval(t).1.scale(sign(*forward));
                let normal = match face.surface {
                    Surface::Plane { n, .. } => n,
                    Surface::EllipticCylinder { .. } => {
                        return Err(unsupported("transverse sections across elliptic-cylinder faces are not supported"));
                    }
                    Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => return Err(unsupported("transverse sections across revolved or spline faces are not supported")),
                    Surface::Cylinder { p, a, .. } => {
                        let d = q.sub(p);
                        d.sub(a.scale(d.dot(a)))
                            .unit()
                            .map_err(|_| invalid("section on cylinder axis"))?
                    }
                }
                .scale(sign(face.sense));
                let inward = normal.cross(tangent);
                if points
                    .iter()
                    .any(|(p, _): &(Vec3, Vec3)| p.distance(q) <= tol)
                {
                    return Err(unsupported("coincident oriented section crossings"));
                }
                points.push((q, inward));
            }
        }
        if points.is_empty() {
            continue;
        }
        // A whole cylindrical wall can contribute both circle branches. Split
        // by its actual radial coordinate; each branch must be one axial strip.
        let coordinate = match face.surface {
            Surface::Cylinder { a, .. } if a.y.abs() == 1. => Some(2),
            Surface::Cylinder { a, .. } if a.z.abs() == 1. => Some(1),
            _ => None,
        };
        let mut groups: Vec<Vec<(Vec3, Vec3)>> = Vec::new();
        for point in points {
            let found = groups.iter().position(|g| {
                coordinate.is_none_or(|k| (xyz(g[0].0)[k] - xyz(point.0)[k]).abs() <= tol)
            });
            if let Some(i) = found {
                groups[i].push(point);
            } else {
                groups.push(vec![point]);
            }
        }
        for mut group in groups {
            let span = match face.surface {
                Surface::Plane { n, .. } if n.y == 0. => 1,
                Surface::Plane { n, .. } if n.z == 0. => 2,
                Surface::Cylinder { a, .. } if a.y.abs() == 1. => 1,
                Surface::Cylinder { a, .. } if a.z.abs() == 1. => 2,
                _ => return Err(unsupported("section surface has no axial span")),
            };
            group.sort_by(|a, b| xyz(a.0)[span].total_cmp(&xyz(b.0)[span]));
            if group.len() % 2 != 0 {
                return Err(unsupported(
                    "transverse section has an odd number of oriented boundaries",
                ));
            }
            for pair in group.chunks_exact(2) {
                if xyz(pair[0].1)[span] <= 0. || xyz(pair[1].1)[span] >= 0. {
                    return Err(invalid(
                        "transverse face interval orientation is inconsistent",
                    ));
                }
                hits.push(hit(face, &[pair[0].0, pair[1].0], tol)?);
            }
        }
    }
    if hits.len() > 1024 {
        return Err(unsupported("too many transverse section segments"));
    }
    Ok(join_seams(hits, tol))
}
/// An earlier transverse cut leaves coplanar partition seams: one plane or cylinder split into faces that meet
/// along an edge (Bob's L3 Angle 8 Coped Both Ends: the heel face split at the first cope's end, y 1.75). Their
/// section segments lie on the IDENTICAL graph, face the same way and share the seam's crossing point, so
/// together they are one boundary segment: they are joined, exactly. A segment pair that differs in graph or
/// orientation, or abuts only within a gap, is left alone (and refused later as before).
fn join_seams(mut hits: Vec<Hit>, tol: f64) -> Vec<Hit> {
    hits.sort_by(|a, b| a.lo.total_cmp(&b.lo));
    let mut joined: Vec<Hit> = Vec::with_capacity(hits.len());
    for h in hits {
        match joined
            .iter_mut()
            .find(|p| p.role == h.role && p.graph == h.graph && (h.lo - p.hi).abs() <= tol)
        {
            Some(p) => p.hi = p.hi.max(h.hi),
            None => joined.push(h),
        }
    }
    joined
}
fn sections_with_budget(
    body: &Body,
    lo: f64,
    hi: f64,
    tol: f64,
    work: &mut usize,
) -> Result<Vec<Band>> {
    let mid = lo + (hi - lo) * 0.5;
    let hits = boundary_hits(body, mid, tol)?;
    if hits.is_empty() {
        return Ok(Vec::new());
    }
    let ys = ordered_graphs(
        hits.iter()
            .filter(|h| h.role < 2)
            .map(|h| h.graph)
            .collect(),
        lo,
        hi,
        tol,
    )?;
    let zs = ordered_graphs(
        hits.iter()
            .filter(|h| h.role >= 2)
            .map(|h| h.graph)
            .collect(),
        lo,
        hi,
        tol,
    )?;
    if ys.len() < 2 || zs.len() < 2 || ys.len() * zs.len() > 4096 {
        return Err(unsupported(
            "transverse arrangement dimensions exceed limits",
        ));
    }
    let (ny, nz) = (ys.len() - 1, zs.len() - 1);
    spend(work, ny * nz * 2 + hits.len() * (ny + nz))?;
    let mut ye = vec![vec![0_i32; nz]; ny + 1];
    let mut ze = vec![vec![0_i32; nz + 1]; ny];
    for h in &hits {
        let (coords, spans) = if h.role < 2 { (&ys, &zs) } else { (&zs, &ys) };
        let index = |v: f64, graphs: &[Graph]| {
            graphs
                .iter()
                .position(|g| (g.value(mid) - v).abs() <= tol)
                .ok_or_else(|| invalid("section boundary misses an arrangement coordinate"))
        };
        let k = index(h.graph.value(mid), coords)?;
        let start = index(h.lo, spans)?;
        let end = index(h.hi, spans)?;
        if start >= end {
            return Err(invalid("section boundary has no positive span"));
        }
        for j in start..end {
            let count = if h.role < 2 {
                &mut ye[k][j]
            } else {
                &mut ze[j][k]
            };
            if *count != 0 {
                return Err(invalid("duplicate or overlapping section boundaries"));
            }
            *count = if h.role % 2 == 0 { 1 } else { -1 };
        }
    }
    let mut material = vec![vec![false; nz]; ny];
    for z in 0..nz {
        let mut winding = 0;
        for y in 0..=ny {
            winding += ye[y][z];
            if winding != 0 && winding != 1 {
                return Err(invalid("inconsistent transverse Y winding"));
            }
            if y < ny {
                material[y][z] = winding == 1;
            }
        }
        if winding != 0 {
            return Err(invalid("unclosed transverse Y boundary"));
        }
    }
    for y in 0..ny {
        let mut winding = 0;
        for z in 0..=nz {
            winding += ze[y][z];
            if winding != 0 && winding != 1 {
                return Err(invalid("inconsistent transverse Z winding"));
            }
            if z < nz && material[y][z] != (winding == 1) {
                return Err(invalid(
                    "transverse boundary coverage disagrees between axes",
                ));
            }
        }
        if winding != 0 {
            return Err(invalid("unclosed transverse Z boundary"));
        }
    }
    let mut result = Vec::new();
    for y in 0..ny {
        for z in 0..nz {
            if material[y][z] {
                result.push(Band {
                    lo,
                    hi,
                    y0: ys[y],
                    y1: ys[y + 1],
                    z0: zs[z],
                    z1: zs[z + 1],
                });
            }
        }
    }
    Ok(result)
}

fn subtract_sections(
    body: &[Band],
    tool: &[Band],
    lo: f64,
    hi: f64,
    tol: f64,
    work: &mut usize,
) -> Result<Vec<Band>> {
    if body.is_empty() {
        return Ok(Vec::new());
    }
    let ys = ordered_graphs(
        body.iter().chain(tool).flat_map(|b| [b.y0, b.y1]).collect(),
        lo,
        hi,
        tol,
    )?;
    let zs = ordered_graphs(
        body.iter().chain(tool).flat_map(|b| [b.z0, b.z1]).collect(),
        lo,
        hi,
        tol,
    )?;
    if ys.len() * zs.len() > 4096 {
        return Err(unsupported(
            "transverse subtraction arrangement exceeds limits",
        ));
    }
    spend(work, ys.len() * zs.len() * (body.len() + tool.len() + 1))?;
    let mid = (lo + hi) * 0.5;
    let mut result = Vec::new();
    for y in ys.windows(2) {
        for z in zs.windows(2) {
            let (my, mz) = (
                (y[0].value(mid) + y[1].value(mid)) * 0.5,
                (z[0].value(mid) + z[1].value(mid)) * 0.5,
            );
            if body.iter().any(|b| Rect::at(b, mid).contains(my, mz))
                && !tool.iter().any(|b| Rect::at(b, mid).contains(my, mz))
            {
                result.push(Band {
                    lo,
                    hi,
                    y0: y[0],
                    y1: y[1],
                    z0: z[0],
                    z1: z[1],
                });
            }
        }
    }
    if result.len() > 128 {
        return Err(unsupported("too many transverse material cells"));
    }
    Ok(result)
}

fn neighbors(a: &Band, b: &Band, tol: f64) -> Option<(usize, usize)> {
    let mid = (a.lo + a.hi) * 0.5;
    let (a, b) = (Rect::at(a, mid), Rect::at(b, mid));
    if (a.y0 - b.y0).abs() <= tol && (a.y1 - b.y1).abs() <= tol {
        if (a.z0 - b.z1).abs() <= tol {
            return Some((0, 1));
        }
        if (a.z1 - b.z0).abs() <= tol {
            return Some((1, 0));
        }
    }
    if (a.z0 - b.z0).abs() <= tol && (a.z1 - b.z1).abs() <= tol {
        if (a.y0 - b.y1).abs() <= tol {
            return Some((2, 3));
        }
        if (a.y1 - b.y0).abs() <= tol {
            return Some((3, 2));
        }
    }
    None
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    let (body, tool) = (&super::snap_directions(body, tol), &super::snap_directions(tool, tol));
    cut_in_frame(body, tool, tol, frame(body, tool)?)
}
pub(super) fn cut_along(body: &Body, tool: &Body, axis: usize, tol: f64) -> Result<Body> {
    let (body, tool) = (&super::snap_directions(body, tol), &super::snap_directions(tool, tol));
    cut_in_frame(body, tool, tol, frame_for_axis_with(body, Some(tool), axis)?)
}
fn cut_in_frame(body: &Body, tool: &Body, tol: f64, f: Frame) -> Result<Body> {
    if !tol.is_finite() || tol <= 0. {
        return Err(invalid("invalid transverse cut tolerance"));
    }
    // Translate geometry, never directions, to a nearby computation origin.
    // This avoids evaluating tiny circular roofs against world coordinates.
    let origin = body.vertices.first().copied().unwrap_or(p(0., 0., 0.));
    let shift = |body: &Body, by: Vec3| {
        let mut out = body.clone();
        for v in &mut out.vertices {
            *v = v.add(by);
        }
        for e in &mut out.edges {
            e.curve = e.curve.translated(by);
        }
        for f in &mut out.faces {
            match &mut f.surface {
                Surface::Spline { data } => **data = data.mapped(|p| p.add(by)),
                Surface::Plane { p, .. } | Surface::Cylinder { p, .. } | Surface::EllipticCylinder { p, .. } | Surface::Cone { p, .. } | Surface::Sphere { c: p, .. } | Surface::Torus { c: p, .. } | Surface::Revolution { p, .. } => *p = p.add(by),
                Surface::Conoid { c, p, .. } => { *c = c.add(by); *p = p.add(by); }
            }
        }
        out
    };
    let body = transform(&shift(body, origin.scale(-1.)), |p| f.local(p));
    let tool = transform(&shift(tool, origin.scale(-1.)), |p| f.local(p));
    let result = super::split::cut(&body, &tool, tol, cut_local)?;
    Ok(shift(&transform(&result, |p| f.world(p)), origin))
}
/// The transverse cut in the local chart (x stations, y the beam axis, z the sweep).
fn cut_local(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    let (body, tool) = (body.clone(), tool.clone());
    validate(&body, tol)?;
    validate(&tool, tol)?;
    let mut xs = events(&body);
    xs.extend(events(&tool));
    let mut xs = stations(xs, tol)?;
    // An authored vertex may differ from its curve endpoint within the input
    // tolerance. At an extremum choose the analytic station, so a round tool
    // closes to zero width there instead of inventing a tiny transverse cap.
    let mut extrema = Vec::new();
    for operand in [&body, &tool] {
        for edge in &operand.edges {
            if let Some((c, u, v, t0, t1)) = section_projection(&edge.curve) {
                if u.x.hypot(v.x) == 0. {
                    continue;
                }
                // The analytic extremum c.x +- r for a circle whose plane holds the station direction (a cylinder
                // across x), not the circle evaluated through its own frame: a wall arc left by an earlier cut
                // starts elsewhere round its circle, and evaluating it at its extremum gave 0.8299999999999998 for
                // cx + r = 0.83 (Bob's Pipe 1.25 18 Copes Rotated 90), whose two ulps acos turns into 2e-8 of angle,
                // 1e-8 off the vertex and 7.6e-8 in^3. Otherwise the evaluated extremum, as before.
                let exact = match edge.curve {
                    Curve::Arc { n, r, .. } if n.x == 0. => Some(r),
                    Curve::CylinderPair { a, rho, .. } if a.x == 0. => Some(rho),
                    _ => None,
                };
                let phase = v.x.atan2(u.x);
                for (k, angle) in [phase, phase + PI].into_iter().enumerate() {
                    for t in parameters(angle, t0, t1) {
                        extrema.push(match exact {
                            Some(r) => c.x + if k == 0 { r } else { -r },
                            None => eval(&edge.curve, t).x,
                        });
                    }
                }
            }
        }
    }
    for x in &mut xs {
        if let Some(value) = extrema.iter().find(|v| (**v - *x).abs() <= tol) {
            *x = *value;
        }
    }
    if xs.windows(2).any(|w| w[1] - w[0] <= tol) {
        return Err(unsupported(
            "ambiguous transverse extremal stations within tolerance",
        ));
    }
    if xs.len().saturating_mul(body.edges.len() + tool.edges.len()) > 2_000_000 {
        return Err(unsupported("transverse recognition exceeds its work limit"));
    }
    let mut arrangement_work = 0;
    xs = intersection_stations(&body, &tool, &xs, &extrema, tol, &mut arrangement_work)?;
    if xs.len().saturating_mul(body.edges.len() + tool.edges.len()) > 2_000_000 {
        return Err(unsupported(
            "refined transverse recognition exceeds its work limit",
        ));
    }
    let mut original = Vec::new();
    let mut cutters = Vec::new();
    let mut bands = Vec::new();
    for x in xs.windows(2) {
        let b = sections_with_budget(&body, x[0], x[1], tol, &mut arrangement_work)?;
        let t = sections_with_budget(&tool, x[0], x[1], tol, &mut arrangement_work)?;
        let remaining = subtract_sections(&b, &t, x[0], x[1], tol, &mut arrangement_work)?;
        original.push(b);
        cutters.push(t);
        bands.push(remaining);
    }
    verify_interfaces(&body, &xs, &original, tol)?;
    verify_interfaces(&tool, &xs, &cutters, tol)?;
    let total: usize = bands.iter().map(Vec::len).sum();
    if total == 0 {
        return Ok(Body::default());
    }
    if total > 4096 {
        return Err(unsupported("transverse result exceeds its patch limit"));
    }
    // Material connectivity follows positive-area overlaps at band interfaces.
    // Separate intervals in one section may reconnect around a through hole.
    let mut adjacency = vec![Vec::new(); total];
    let mut offsets = vec![0];
    for band in &bands {
        offsets.push(offsets.last().unwrap() + band.len());
    }
    for i in 1..bands.len() {
        for (j, left) in bands[i - 1].iter().enumerate() {
            let l = Rect::at(left, xs[i]);
            for (k, right) in bands[i].iter().enumerate() {
                let r = Rect::at(right, xs[i]);
                if l.y1.min(r.y1) - l.y0.max(r.y0) > tol && l.z1.min(r.z1) - l.z0.max(r.z0) > tol {
                    let a = offsets[i - 1] + j;
                    let b = offsets[i] + k;
                    adjacency[a].push(b);
                    adjacency[b].push(a);
                }
            }
        }
    }
    for (i, cells) in bands.iter().enumerate() {
        for j in 0..cells.len() {
            for k in j + 1..cells.len() {
                if neighbors(&cells[j], &cells[k], tol).is_some() {
                    let (a, b) = (offsets[i] + j, offsets[i] + k);
                    adjacency[a].push(b);
                    adjacency[b].push(a);
                }
            }
        }
    }
    // A removed interval collapsing against an exterior X boundary creates
    // four incident sheets along one edge, not a manifold solid. Detect this
    // geometric tangency before assembly; ordinary interior hole tips have
    // material on the opposite side and are unaffected by this check.
    for i in 0..xs.len() {
        let left = if i > 0 { bands[i - 1].as_slice() } else { &[] };
        let right = bands.get(i).map(Vec::as_slice).unwrap_or(&[]);
        let side = if left.is_empty() {
            right
        } else if right.is_empty() {
            left
        } else {
            continue;
        };
        for j in 0..side.len() {
            for k in j + 1..side.len() {
                let (a, b) = (Rect::at(&side[j], xs[i]), Rect::at(&side[k], xs[i]));
                let mid = (side[j].lo + side[j].hi) * 0.5;
                let (am, bm) = (Rect::at(&side[j], mid), Rect::at(&side[k], mid));
                let y_contact = ((a.y1 - b.y0).abs() <= tol && bm.y0 - am.y1 > tol)
                    || ((b.y1 - a.y0).abs() <= tol && am.y0 - bm.y1 > tol);
                let z_contact = ((a.z1 - b.z0).abs() <= tol && bm.z0 - am.z1 > tol)
                    || ((b.z1 - a.z0).abs() <= tol && am.z0 - bm.z1 > tol);
                if (y_contact && a.z1.min(b.z1) - a.z0.max(b.z0) > tol)
                    || (z_contact && a.y1.min(b.y1) - a.y0.max(b.y0) > tol)
                {
                    return Err(unsupported(
                        "cut meets an exterior boundary tangentially along an edge",
                    ));
                }
            }
        }
    }
    let mut reached = vec![false; total];
    let mut stack = vec![0];
    reached[0] = true;
    while let Some(i) = stack.pop() {
        for j in &adjacency[i] {
            if !reached[*j] {
                reached[*j] = true;
                stack.push(*j);
            }
        }
    }
    if reached.iter().any(|r| !r) {
        return Err(unsupported(
            "transverse cut would leave disconnected material components",
        ));
    }
    let mut result = Vec::new();
    for cells in &bands {
        let mut exposed = vec![[true; 4]; cells.len()];
        for i in 0..cells.len() {
            for j in i + 1..cells.len() {
                if let Some((a, b)) = neighbors(&cells[i], &cells[j], tol) {
                    exposed[i][a] = false;
                    exposed[j][b] = false;
                }
            }
        }
        for (i, b) in cells.iter().enumerate() {
            result.extend(selected_patches(b, exposed[i], tol)?);
        }
    }
    for i in 0..xs.len() {
        let left = if i > 0 { bands[i - 1].as_slice() } else { &[] };
        let right = bands.get(i).map(Vec::as_slice).unwrap_or(&[]);
        result.extend(end_patches(left, right, xs[i], tol)?);
    }
    assemble(result, tol)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn analytic_intersection_events_check_branches_and_keep_tangencies_explicit() {
        let horizontal = |b| Graph::Line { m: 0., b };
        let circle = |cx, cv, r, branch| Graph::Circle { cx, cv, r, branch };
        let upper = circle(0., 0., 1., 1.);
        let mut roots = intersections(horizontal(0.5), upper, -1., 1., 1e-9).unwrap();
        roots.sort_by(f64::total_cmp);
        assert_eq!(roots.len(), 2);
        assert!((roots[0] + 0.75_f64.sqrt()).abs() < 1e-14);
        assert!((roots[1] - 0.75_f64.sqrt()).abs() < 1e-14);
        assert!(intersections(horizontal(-0.5), upper, -1., 1., 1e-9)
            .unwrap()
            .is_empty());
        assert!(intersections(horizontal(1.), upper, -1., 1., 1e-9).is_err());
        assert!(intersections(horizontal(1.01), upper, -1., 1., 1e-9)
            .unwrap()
            .is_empty());
        let roots = intersections(upper, circle(1., 0., 1., 1.), 0., 1., 1e-9).unwrap();
        assert!(!roots.is_empty() && roots.iter().all(|x| (*x - 0.5).abs() < 1e-14));
        assert!(intersections(upper, circle(1., 0., 1., -1.), 0., 1., 1e-9)
            .unwrap()
            .is_empty());
        assert!(
            intersections(upper, circle(0., 0., 0.9, 1.), -0.9, 0.9, 1e-9)
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            intersections(horizontal(2.), Graph::Line { m: 2., b: 0. }, 0., 2., 1e-9).unwrap(),
            vec![1.]
        );
        assert!(intersections(horizontal(2.), horizontal(1.), 0., 2., 1e-9)
            .unwrap()
            .is_empty());
        // A distinct sub-tolerance separation is not a coincident graph, and
        // a sub-tolerance interior reversal is still a failed ordering proof.
        assert!(ordered_graphs(vec![horizontal(0.), horizontal(1e-12)], -1., 1., 1e-6).is_err());
        assert!(!positive(
            horizontal(0.),
            Graph::Line { m: 1., b: 0. },
            -1e-10,
            1.,
            1e-6
        ));
        let lower = horizontal(-0.4721303058753803);
        let upper = Graph::Line {
            m: 6.001499981278495,
            b: 38.42878924596892,
        };
        let x = intersections(lower, upper, -7., -6., 1e-6).unwrap()[0];
        assert!(positive(lower, upper, x, x + 0.1, 1e-6));
        assert!(positive(upper, lower, x - 0.1, x, 1e-6));
    }
    fn survey(case: &str, name: &str) -> Body {
        let bytes = std::fs::read(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join(format!("tests/fixtures/survey/{case}/{name}.f64")),
        )
        .unwrap();
        let values: Vec<_> = bytes
            .chunks_exact(8)
            .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
            .collect();
        crate::neutral::Reader::new(&values).body().unwrap()
    }
    fn part(case: &str, name: &str) -> Body {
        let bytes = std::fs::read(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join(format!("tests/fixtures/parts/{case}/{name}.f64")),
        )
        .unwrap();
        let values: Vec<_> = bytes
            .chunks_exact(8)
            .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
            .collect();
        crate::neutral::Reader::new(&values).body().unwrap()
    }
    #[test]
    fn actual_rect_tube_miter_matches_section_moments_and_four_steps_replay() {
        let base = part("rect_tube", "base");
        let tool = part("rect_tube", "step01");
        let (m, intercept) = tool
            .faces
            .iter()
            .find_map(|f| match f.surface {
                Surface::Plane { p, n, .. } if n.x != 0. && n.y != 0. => {
                    Some((-n.x / n.y, p.y + n.x * p.x / n.y))
                }
                _ => None,
            })
            .unwrap();
        // Independent dimensional reference: 8x6 tube, outer corner radius 1,
        // inner .75, length 12. Integrate the clipped linear end height against
        // each outer/inner strip's area and first X moment. No B-rep sections or
        // kernel mass integration are used for this reference.
        let mut expected = 0.;
        for (lo, hi, c, r, constant, sign) in [
            (0., 1., 1., 1., 4., 1.),
            (1., 7., 0., 0., 6., 1.),
            (7., 8., 7., 1., 4., 1.),
            (0.25, 1., 1., 0.75, 4., -1.),
            (1., 7., 0., 0., 5.5, -1.),
            (7., 7.75, 7., 0.75, 4., -1.),
        ] {
            let mut xs = vec![lo, hi];
            for y in [0., 12.] {
                let x = (y - intercept) / m;
                if x > lo && x < hi {
                    xs.push(x);
                }
            }
            xs.sort_by(f64::total_cmp);
            for band in xs.windows(2) {
                let y = m * (band[0] + band[1]) * 0.5 + intercept;
                let (slope, b) = if y <= 0. {
                    (0., 0.)
                } else if y >= 12. {
                    (0., 12.)
                } else {
                    (m, intercept)
                };
                let primitive = |x: f64| {
                    let mut value = constant * (slope * x * x * 0.5 + b * x);
                    if r != 0. {
                        let t = (x - c).clamp(-r, r);
                        let root = (r * r - t * t).max(0.).sqrt();
                        let area = (t * root + r * r * (t / r).asin()) * 0.5;
                        let moment = c * area - root.powi(3) / 3.;
                        value += 2. * (slope * moment + b * area);
                    }
                    value
                };
                expected += sign * (primitive(band[1]) - primitive(band[0]));
            }
        }
        let result = operation(&base, &tool, 1e-6);
        assert!((crate::mass::measure(&result, 1e-6).unwrap().volume - expected).abs() < 1e-9);
        let f = frame_for_axis(&base, 2).unwrap();
        let local = transform(&result, |q| f.local(q));
        let xs = stations(events(&local), 1e-6).unwrap();
        let bands: Vec<_> = xs
            .windows(2)
            .map(|w| sections(&local, w[0], w[1], 1e-6).unwrap())
            .collect();
        for ix in 0..41 {
            for iy in 0..37 {
                for iz in 0..31 {
                    let q = p(
                        -0.213 + 0.217 * ix as f64,
                        -0.191 + 0.353 * iy as f64,
                        -6.119 + 0.211 * iz as f64,
                    );
                    let d2 =
                        (q.x - q.x.clamp(1., 7.)).powi(2) + (q.z - q.z.clamp(-5., -1.)).powi(2);
                    let expected = d2 < 1.
                        && d2 > 0.75_f64.powi(2)
                        && q.y > 0.
                        && q.y < 12.
                        && q.y < m * q.x + intercept;
                    let ql = f.local(q);
                    let actual = xs
                        .windows(2)
                        .position(|w| ql.x > w[0] && ql.x < w[1])
                        .is_some_and(|i| {
                            bands[i]
                                .iter()
                                .any(|b| Rect::at(b, ql.x).contains(ql.y, ql.z))
                        });
                    assert_eq!(actual, expected, "miter probe {q:?}");
                }
            }
        }
        // Independent native mesh references recorded by Claude, 9ebdf7d.
        // Their tessellation error is larger than our analytic integral error.
        let mut body = result;
        for (step, native_volume) in [(2, 58.662702784), (3, 55.725486031), (4, 54.783702062)] {
            let tool = part("rect_tube", &format!("step{step:02}"));
            let before = crate::neutral::write_body(&body);
            let result = operation(&body, &tool, 1e-6);
            assert_eq!(before, crate::neutral::write_body(&body));
            crate::topology::check(&result, 1e-6).unwrap();
            assert!(
                (crate::mass::measure(&result, 1e-6).unwrap().volume - native_volume).abs() < 0.001
            );
            body = result;
        }
        let base = part("str1", "base");
        let result = operation(&base, &part("str1", "step01"), 1e-6);
        crate::topology::check(&result, 1e-6).unwrap();
        assert!(
            (crate::mass::measure(&result, 1e-6).unwrap().volume - 112.845646699).abs() < 0.001
        );
        // Its next wedges were refused as "two different stored lines"; they are one plane stored twice, 1 to 2 ulps
        // apart (Claude 2026-10-31, same_graph), so all six steps now replay, each within the reference's mesh
        // allowance (tests/fixtures/parts/str1/reference.md: Eyeshot's step-by-step mesh volumes).
        let mut body = result;
        for (step, native_volume) in [(2, 110.414764497), (3, 106.701206509), (4, 101.28852158), (5, 101.231382506), (6, 98.800405878)] {
            let result = operation(&body, &part("str1", &format!("step{step:02}")), 1e-6);
            crate::topology::check(&result, 1e-6).unwrap();
            let volume = crate::mass::measure(&result, 1e-6).unwrap().volume;
            assert!((volume - native_volume).abs() < 0.001, "STR-1 step {step}: {volume} vs {native_volume}");
            body = result;
        }
    }
    #[test]
    fn rectangular_tube_multiple_walls_preserve_the_unaffected_stock() {
        let base = survey("rect_tube_hole_both_walls", "base");
        let n = p(1., 0., 0.);
        let slot = vec![
            segment(p(0., -1., 11.625), p(0., 1., 11.625)),
            Curve::Arc {
                c: p(0., 1., 12.),
                n,
                x: p(0., 0., -1.),
                r: 0.375,
                t0: 0.,
                t1: PI,
            },
            segment(p(0., 1., 12.375), p(0., -1., 12.375)),
            Curve::Arc {
                c: p(0., -1., 12.),
                n,
                x: p(0., 0., 1.),
                r: 0.375,
                t0: 0.,
                t1: PI,
            },
        ];
        let corners = [
            p(0., -0.5, 11.),
            p(0., 0.5, 11.),
            p(0., 0.5, 13.),
            p(0., -0.5, 13.),
        ];
        let rectangle = (0..4)
            .map(|i| segment(corners[i], corners[(i + 1) % 4]))
            .collect();
        let make = |outline| {
            crate::extrude::extrude(
                &crate::neutral::Profile {
                    loops: vec![outline],
                    n,
                    low: -3.,
                    high: 3.,
                },
                1e-6,
            )
            .unwrap()
        };
        let original_volume = 72. + 4.5 * PI;
        let original = crate::neutral::write_body(&base);
        for (tool, removed) in [
            (survey("rect_tube_hole_both_walls", "step01"), PI * 0.125),
            (make(slot), 0.5 * (1.5 + PI * 0.375_f64.powi(2))),
            (make(rectangle), 1.),
        ] {
            let out = operation(&base, &tool, 1e-6);
            assert_eq!(crate::topology::check(&out, 1e-6).unwrap().genus(), 3.);
            assert!(
                (crate::mass::measure(&out, 1e-6).unwrap().volume - (original_volume - removed))
                    .abs()
                    < 1e-8
            );
            assert_eq!(&out.vertices[..base.vertices.len()], &base.vertices);
            for (i, f) in base.faces.iter().enumerate() {
                if matches!(f.surface, Surface::Cylinder { .. }) {
                    assert_eq!(format!("{f:?}"), format!("{:?}", out.faces[i]));
                }
            }
            let next = operation(&out, &placed(&tool, 1., p(0., 0., 6.)), 1e-6);
            assert_eq!(crate::topology::check(&next, 1e-6).unwrap().genus(), 5.);
            assert!(
                (crate::mass::measure(&next, 1e-6).unwrap().volume
                    - (original_volume - 2. * removed))
                    .abs()
                    < 1e-8
            );
            assert!(matches!(
                crate::boolean::cut(&out, &tool, 1e-6),
                Err(KernelError::Unsupported(_))
            ));
            for (scale, turns) in [(0.1_f64, 0), (1., 1), (25.4, 2)] {
                let map = |q: Vec3| {
                    let a = xyz(q);
                    p(a[turns], a[(turns + 1) % 3], a[(turns + 2) % 3])
                };
                let offset = p(100., -50., 25.);
                let moved = operation(
                    &placed(&transform(&base, map), scale, offset),
                    &placed(&transform(&tool, map), scale, offset),
                    scale * 1e-6,
                );
                assert!(
                    (crate::mass::measure(&moved, scale * 1e-6).unwrap().volume / scale.powi(3)
                        - (original_volume - removed))
                        .abs()
                        < 1e-8
                );
            }
        }
        // Finite tools ending inside one wall or starting inside the other must
        // not accidentally pair a gap's exit/entry faces as solid material.
        let tool = survey("rect_tube_hole_one_wall", "step01");
        for offset in [-3.4, 0.4] {
            assert!(matches!(
                crate::boolean::cut(&base, &placed(&tool, 1., p(offset, 0., 0.)), 1e-6),
                Err(KernelError::Unsupported(_))
            ));
        }
        assert_eq!(crate::neutral::write_body(&base), original);
    }
    #[test]
    fn square_and_exact_planar_miter_ends_match_centroid_integrals() {
        for (case, area) in [
            ("pipe_square_end", 0.9375 * PI),
            ("rect_tube_square_end", 3. + 0.1875 * PI),
        ] {
            let base = survey(case, "base");
            let square = survey(case, "step01");
            let n = p(1., 0., 1.).unit().unwrap();
            let u = p(n.z, 0., -n.x);
            let v = p(0., 1., 0.);
            let c = p(0., 0., 22.);
            let corners = [
                c.sub(u.scale(10.)).sub(v.scale(10.)),
                c.add(u.scale(10.)).sub(v.scale(10.)),
                c.add(u.scale(10.)).add(v.scale(10.)),
                c.sub(u.scale(10.)).add(v.scale(10.)),
            ];
            let miter = crate::extrude::extrude(
                &crate::neutral::Profile {
                    loops: vec![(0..4)
                        .map(|i| segment(corners[i], corners[(i + 1) % 4]))
                        .collect()],
                    n,
                    low: 0.,
                    high: 10.,
                },
                1e-6,
            )
            .unwrap();
            for (name, tool, height) in [("square", &square, 23.), ("miter", &miter, 22.)] {
                let result = operation(&base, tool, 1e-6);
                let volume = crate::mass::measure(&result, 1e-6).unwrap().volume;
                eprintln!("{case} {name}: {volume:.14}");
                assert!((volume - area * height).abs() < 1e-8);
                assert_eq!(crate::topology::check(&result, 1e-6).unwrap().genus(), 1.);
                let f = frame_for_axis(&base, if name == "square" { 0 } else { 1 }).unwrap();
                let local = transform(&result, |q| f.local(q));
                let xs = stations(events(&local), 1e-6).unwrap();
                let bs: Vec<_> = xs
                    .windows(2)
                    .map(|w| sections(&local, w[0], w[1], 1e-6).unwrap())
                    .collect();
                for ix in 0..25 {
                    for iy in 0..25 {
                        for iz in 0..21 {
                            let q = p(
                                -2.171 + ix as f64 * 0.1813,
                                -2.137 + iy as f64 * 0.1771,
                                19.137 + iz as f64 * 0.293,
                            );
                            let stock = if case.starts_with("pipe") {
                                let d = q.x * q.x + q.y * q.y;
                                d > 1.75_f64.powi(2) && d < 4.
                            } else {
                                let d = (q.x.abs() - 1.5).max(0.).powi(2)
                                    + (q.y.abs() - 1.5).max(0.).powi(2);
                                d > 0.0625 && d < 0.25
                            };
                            let expected =
                                stock && q.z < height - if name == "miter" { q.x } else { 0. };
                            let ql = f.local(q);
                            let actual = xs
                                .windows(2)
                                .position(|w| ql.x > w[0] && ql.x < w[1])
                                .is_some_and(|i| {
                                    bs[i].iter().any(|b| Rect::at(b, ql.x).contains(ql.y, ql.z))
                                });
                            assert_eq!(actual, expected, "{case} {name}: {q:?}");
                        }
                    }
                }
                let deeper = operation(&result, &placed(tool, 1., p(0., 0., -1.)), 1e-6);
                assert!(
                    (crate::mass::measure(&deeper, 1e-6).unwrap().volume - area * (height - 1.))
                        .abs()
                        < 1e-8
                );
                for (scale, turns) in [(0.1_f64, 0), (1., 1), (25.4, 2)] {
                    let map = |q: Vec3| {
                        let a = xyz(q);
                        p(a[turns], a[(turns + 1) % 3], a[(turns + 2) % 3])
                    };
                    let offset = p(1e5, -1e5, 1e5);
                    // Move the miter down before this test: it must not make a
                    // numerically ambiguous near-contact with the old end cap.
                    let shift = if name == "miter" { -1. } else { 0. };
                    let tool = placed(tool, 1., p(0., 0., shift));
                    let moved = operation(
                        &placed(&transform(&base, map), scale, offset),
                        &placed(&transform(&tool, map), scale, offset),
                        scale * 1e-6,
                    );
                    assert!(
                        (crate::mass::measure(&moved, scale * 1e-6).unwrap().volume
                            / scale.powi(3)
                            - area * (height + shift))
                            .abs()
                            < 1e-6
                    );
                }
                if let Ok(folder) = std::env::var("GCKERNEL_ASTRA_END_EXPORT") {
                    let folder = std::path::PathBuf::from(folder);
                    std::fs::create_dir_all(&folder).unwrap();
                    for (suffix, b) in [("base", &base), ("tool", tool), ("result", &result)] {
                        let bytes: Vec<_> = crate::neutral::write_body(b)
                            .into_iter()
                            .flat_map(f64::to_le_bytes)
                            .collect();
                        std::fs::write(folder.join(format!("{case}_{name}_{suffix}.f64")), bytes)
                            .unwrap();
                    }
                }
            }
            let hole_case = if case.starts_with("pipe") {
                "pipe_hole_one_wall"
            } else {
                "rect_tube_hole_one_wall"
            };
            let hole = survey(hole_case, "step01");
            let removed = if case.starts_with("pipe") {
                22.5 * PI - 70.48768294431949
            } else {
                PI * 0.0625
            };
            for out in [
                operation(&operation(&base, &hole, 1e-6), &square, 1e-6),
                operation(&operation(&base, &square, 1e-6), &hole, 1e-6),
            ] {
                assert!(
                    (crate::mass::measure(&out, 1e-6).unwrap().volume - (23. * area - removed))
                        .abs()
                        < 1e-8
                );
            }
            let tangent = operation(&base, &placed(&square, 1., p(0., 0., 1.)), 1e-6);
            assert!(
                (crate::mass::measure(&tangent, 1e-6).unwrap().volume - area * 24.).abs() < 1e-8
            );
            assert!(matches!(
                crate::boolean::cut(&base, &placed(&square, 1., p(0., 0., -14.)), 1e-6),
                Err(KernelError::Unsupported(_))
            ));
        }
        // The captured miter's cap has a tiny nonzero third normal component. The declared normalization contract
        // (boolean::snap_directions, Claude 2026-11-05 under Bob's 10-28 rule): directions within 1e-12 of a
        // coordinate axis are that axis, so the miter cuts to one valid solid.
        for case in ["pipe_miter_45", "rect_tube_miter_45"] {
            let base = survey(case, "base");
            let out = crate::boolean::cut(&base, &survey(case, "step01"), 1e-6).unwrap();
            assert_eq!(crate::topology::check(&out, 1e-6).unwrap().genus(), crate::topology::check(&base, 1e-6).unwrap().genus());
            let (before, after) = (crate::mass::measure(&base, 1e-6).unwrap().volume, crate::mass::measure(&out, 1e-6).unwrap().volume);
            assert!(after > 0. && after < before, "{case}: {after} of {before}");
        }
    }
    #[test]
    fn near_distinct_graphs_and_interior_contacts_are_not_merged() {
        let line = Graph::Line { m: 0., b: 0. };
        assert_eq!(
            ordered_graphs(vec![line, line], 0., 1., 1e-6)
                .unwrap()
                .len(),
            1
        );
        assert!(matches!(
            ordered_graphs(vec![line, Graph::Line { m: 0., b: 5e-7 }], 0., 1., 1e-6),
            Err(KernelError::Unsupported(_))
        ));
        // One plane stored two ways (aa5's miter: the cut's face and the preparation wedge's), apart only by the
        // roundoff of their normals and points: one graph. 1e-9 apart is beyond roundoff and still refused.
        let (miter, wedge) = (Graph::Line { m: 0.99999999999999986, b: -44.003543307086616 },
            Graph::Line { m: 1.0000000000000048, b: -44.00354330708683 });
        assert_eq!(ordered_graphs(vec![miter, wedge], 43.3, 44.01, 1e-6).unwrap().len(), 1);
        assert!(ordered_graphs(vec![miter, Graph::Line { m: 1., b: -44.003543307086616 + 1e-9 }], 43.3, 44.01, 1e-6).is_err());
        let arc = Graph::Circle {
            cx: 0.,
            cv: 0.,
            r: 1.,
            branch: 1.,
        };
        // Both endpoints lie below this line, but the circle rises above it
        // in the interior. Endpoint ordering alone would miss both crossings.
        let tilted = Graph::Line { m: 0.1, b: 0.95 };
        assert!(!below(arc, tilted, -0.5, 0.5, 0.));
        assert!(matches!(
            ordered_graphs(vec![arc, tilted], -0.5, 0.5, 1e-6),
            Err(KernelError::Unsupported(_))
        ));
    }
    #[test]
    fn corner_cut_and_curved_blind_pockets_match_independent_integrals() {
        let integral = |f: &dyn Fn(f64) -> f64, lo: f64, hi: f64| {
            let n = 4096;
            let h = (hi - lo) / n as f64;
            let mut sum = f(lo) + f(hi);
            for i in 1..n {
                sum += if i % 2 == 0 { 2. } else { 4. } * f(lo + h * i as f64);
            }
            sum * h / 3.
        };
        let corner_removed = integral(&|t| 0.125 * t.cos().powi(2), -PI / 2., -PI / 6.)
            + integral(
                &|t| {
                    let y = 1.5 + 0.25 * t.sin();
                    2. * (0.25 - (y - 1.75).powi(2)).sqrt()
                        * ((0.25 - (y - 1.5).powi(2)).sqrt() - 0.25 * t.cos())
                        * 0.25
                        * t.cos()
                },
                0.,
                PI / 2.,
            )
            + integral(
                &|t| {
                    let y = 1.5 + 0.5 * t.sin();
                    0.5 * (0.25 - (y - 1.75).powi(2)).sqrt() * t.cos().powi(2)
                },
                PI / 6.,
                PI / 2.,
            );
        let pipe_removed = integral(
            &|t| 0.5 * t.cos().powi(2) * ((4. - 0.25 * t.sin().powi(2)).sqrt() - 1.9),
            -PI / 2.,
            PI / 2.,
        );
        let transition = ((1.75_f64.powi(2) - 1.7_f64.powi(2)).sqrt() / 0.5).asin();
        let clipped_wall = |t: f64| {
            let y = 0.5 * t.sin();
            0.5 * t.cos().powi(2)
                * ((4. - y * y).sqrt() - (1.75_f64.powi(2) - y * y).sqrt().max(1.7))
        };
        let breakthrough_removed = 2.
            * (integral(&clipped_wall, 0., transition)
                + integral(&clipped_wall, transition, PI / 2.));
        let beam_expected = (155.58846611769187 + reference_volume()) * 0.5;
        for (name, base, tool, expected, genus) in [
            (
                "corner",
                survey("rect_tube_hole_over_corner", "base"),
                survey("rect_tube_hole_over_corner", "step01"),
                72. + 4.5 * PI - corner_removed,
                2.,
            ),
            (
                "pipe pocket",
                survey("pipe_hole_one_wall", "base"),
                placed(&survey("pipe_hole_one_wall", "step01"), 1., p(0.9, 0., 0.)),
                22.5 * PI - pipe_removed,
                1.,
            ),
            (
                "pipe partial breakthrough",
                survey("pipe_hole_one_wall", "base"),
                placed(&survey("pipe_hole_one_wall", "step01"), 1., p(0.7, 0., 0.)),
                22.5 * PI - breakthrough_removed,
                2.,
            ),
            (
                "beam pocket",
                fixture("base"),
                placed(&fixture("step01"), 1., p(0., 0., 4.)),
                beam_expected,
                0.,
            ),
        ] {
            let result = operation(&base, &tool, 1e-6);
            let actual = crate::mass::measure(&result, 1e-6).unwrap().volume;
            eprintln!("{name}: {actual:.14}, independent {expected:.14}");
            assert!((actual - expected).abs() < 1e-8, "{name}");
            assert_eq!(
                crate::topology::check(&result, 1e-6).unwrap().genus(),
                genus
            );
            let frame = frame(&base, &tool).unwrap();
            let local = transform(&result, |q| frame.local(q));
            let xs = stations(events(&local), 1e-6).unwrap();
            let bs: Vec<_> = xs
                .windows(2)
                .map(|w| sections(&local, w[0], w[1], 1e-6).unwrap())
                .collect();
            for ix in 0..23 {
                for iy in 0..25 {
                    for iz in 0..19 {
                        let (q, stock, in_tool) = if name == "corner" {
                            let q = p(
                                1.013 + ix as f64 * 0.047,
                                1.21931 + iy as f64 * 0.04077,
                                11.413 + iz as f64 * 0.061,
                            );
                            let d = (q.x.abs() - 1.5).max(0.).powi(2)
                                + (q.y.abs() - 1.5).max(0.).powi(2);
                            (
                                q,
                                d > 0.0625 && d < 0.25,
                                q.x > 1.5 && (q.y - 1.75).powi(2) + (q.z - 12.).powi(2) < 0.25,
                            )
                        } else if name.starts_with("pipe") {
                            let q = p(
                                1.613 + ix as f64 * 0.023,
                                -0.613 + iy as f64 * 0.051,
                                11.413 + iz as f64 * 0.061,
                            );
                            let r2 = q.x * q.x + q.y * q.y;
                            (
                                q,
                                r2 > 1.75_f64.powi(2) && r2 < 4.,
                                q.x > (if name == "pipe pocket" { 1.9 } else { 1.7 })
                                    && q.y * q.y + (q.z - 12.).powi(2) < 0.25,
                            )
                        } else {
                            let q = p(
                                0.313 + ix as f64 * 0.041,
                                4.413 + iy as f64 * 0.051,
                                -2.713 + iz as f64 * 0.081,
                            );
                            let (z0, z1) = if q.x < 0.74 {
                                let root = (0.25 - (q.x - 0.74).powi(2)).sqrt();
                                (-2.615 + root, -1.385 - root)
                            } else {
                                (-2.115, -1.885)
                            };
                            (
                                q,
                                q.z > z0 && q.z < z1,
                                q.z > -2. && (q.x - 0.74).powi(2) + (q.y - 5.).powi(2) < 0.16,
                            )
                        };
                        let l = frame.local(q);
                        let actual_inside = xs
                            .windows(2)
                            .position(|w| l.x > w[0] && l.x < w[1])
                            .is_some_and(|i| {
                                bs[i].iter().any(|b| Rect::at(b, l.x).contains(l.y, l.z))
                            });
                        assert_eq!(actual_inside, stock && !in_tool, "{name} probe {q:?}");
                    }
                }
            }
            let offset = if name == "beam pocket" {
                p(0., 10., 0.)
            } else {
                p(0., 0., 6.)
            };
            let again = operation(&result, &placed(&tool, 1., offset), 1e-6);
            let original = if name == "corner" {
                72. + 4.5 * PI
            } else if name.starts_with("pipe") {
                22.5 * PI
            } else {
                155.58846611769187
            };
            assert!(
                (crate::mass::measure(&again, 1e-6).unwrap().volume - (2. * expected - original))
                    .abs()
                    < 1e-8
            );
            for (scale, turns) in [(0.1_f64, 0), (1., 1), (25.4, 2)] {
                let map = |q: Vec3| {
                    let a = xyz(q);
                    p(a[turns], a[(turns + 1) % 3], a[(turns + 2) % 3])
                };
                let offset = p(1e5, -1e5, 1e5);
                let placed_result = operation(
                    &placed(&transform(&base, map), scale, offset),
                    &placed(&transform(&tool, map), scale, offset),
                    scale * 1e-6,
                );
                assert!(
                    (crate::mass::measure(&placed_result, scale * 1e-6)
                        .unwrap()
                        .volume
                        / scale.powi(3)
                        - expected)
                        .abs()
                        < 1e-6,
                    "{name} scale {scale}"
                );
            }
        }
    }
    #[test]
    fn twenty_one_pipe_holes_replay_without_losing_prior_openings() {
        let case = "pipe_hole_one_wall";
        let mut body = survey(case, "base");
        let tool = survey(case, "step01");
        let removed = 22.5 * PI - 70.48768294431949;
        for step in 0..21 {
            let z = 1. + step as f64 * 1.1;
            body = operation(&body, &placed(&tool, 1., p(0., 0., z - 12.)), 1e-6);
            assert_eq!(
                crate::topology::check(&body, 1e-6).unwrap().genus(),
                (step + 2) as f64
            );
            assert!(
                (crate::mass::measure(&body, 1e-6).unwrap().volume
                    - (22.5 * PI - (step + 1) as f64 * removed))
                    .abs()
                    < 1e-8
            );
        }
    }
    #[test]
    fn off_center_pipe_holes_match_independent_wall_integrals() {
        let base = survey("pipe_hole_one_wall", "base");
        for radius in [0.1_f64, 0.25, 0.5, 0.8] {
            for cy in [-1.3_f64, -0.7, 0., 0.4, 1.2] {
                if cy.abs() + radius >= 1.75 {
                    continue;
                }
                let n = p(1., 0., 0.);
                let tool = crate::extrude::extrude(
                    &crate::neutral::Profile {
                        loops: vec![vec![Curve::Arc {
                            c: p(0., cy, 12.),
                            n,
                            x: p(0., 1., 0.),
                            r: radius,
                            t0: 0.,
                            t1: TAU,
                        }]],
                        n,
                        low: 0.,
                        high: 3.,
                    },
                    1e-6,
                )
                .unwrap();
                let result = operation(&base, &tool, 1e-6);
                let integrand = |t: f64| {
                    let y = cy + radius * t.sin();
                    2. * radius
                        * radius
                        * t.cos().powi(2)
                        * ((4. - y * y).sqrt() - (1.75_f64.powi(2) - y * y).sqrt())
                };
                let count = 4096;
                let h = PI / count as f64;
                let mut sum = integrand(-PI / 2.) + integrand(PI / 2.);
                for i in 1..count {
                    sum += if i % 2 == 0 { 2. } else { 4. } * integrand(-PI / 2. + h * i as f64);
                }
                let expected = 22.5 * PI - sum * h / 3.;
                assert!(
                    (crate::mass::measure(&result, 1e-6).unwrap().volume - expected).abs() < 1e-8,
                    "r {radius} cy {cy}"
                );
                assert_eq!(crate::topology::check(&result, 1e-6).unwrap().genus(), 2.);
            }
        }
    }
    #[test]
    fn surveyed_pipe_cuts_admit_oriented_annular_sections() {
        let integrate = |f: &dyn Fn(f64) -> f64, lo: f64, hi: f64| {
            let n = 4096;
            let h = (hi - lo) / n as f64;
            let mut sum = f(lo) + f(hi);
            for i in 1..n {
                sum += if i % 2 == 0 { 2. } else { 4. } * f(lo + h * i as f64);
            }
            sum * h / 3.
        };
        // Independent radial wall thickness of the dimensioned 4 OD x .25 pipe.
        let wall = |y: f64| (4. - y * y).sqrt() - (1.75_f64.powi(2) - y * y).sqrt();
        let hole = integrate(
            &|t| 0.5 * t.cos().powi(2) * wall(0.5 * t.sin()),
            -PI / 2.,
            PI / 2.,
        );
        let slot = 0.75 * integrate(&wall, -1., 1.)
            + 4. * 0.375_f64.powi(2)
                * integrate(
                    &|t| t.cos().powi(2) * wall(1. + 0.375 * t.sin()),
                    0.,
                    PI / 2.,
                );
        let original_volume = 22.5 * PI;
        for (case, removed, holes) in [
            ("pipe_hole_one_wall", hole, 1),
            ("pipe_hole_both_walls", 2. * hole, 2),
            ("pipe_slot_one_wall", slot, 1),
        ] {
            let base = survey(case, "base");
            let tool = survey(case, "step01");
            let out = operation(&base, &tool, 1e-6);
            let expected = original_volume - removed;
            assert!(
                (crate::mass::measure(&out, 1e-6).unwrap().volume - expected).abs() < 1e-8,
                "{case}"
            );
            assert_eq!(
                crate::topology::check(&out, 1e-6).unwrap().genus(),
                (1 + holes) as f64
            );
            let f = frame(&base, &tool).unwrap();
            let local = transform(&out, |q| f.local(q));
            let xs = stations(events(&local), 1e-6).unwrap();
            let bs: Vec<_> = xs
                .windows(2)
                .map(|w| sections(&local, w[0], w[1], 1e-6).unwrap())
                .collect();
            for ix in 0..23 {
                for iy in 0..21 {
                    for z in [10.1, 10.8, 11.3, 11.8, 12.13, 12.7, 13.3, 13.9] {
                        let q = p(-2.15 + ix as f64 * 0.19, -2.13 + iy as f64 * 0.21, z);
                        let radial = q.x * q.x + q.y * q.y;
                        let mut expected_inside = radial > 1.75_f64.powi(2) && radial < 4.;
                        let in_profile = if case == "pipe_slot_one_wall" {
                            (q.y.abs() - 1.).max(0.).powi(2) + (q.z - 12.).powi(2)
                                < 0.375_f64.powi(2)
                        } else {
                            q.y * q.y + (q.z - 12.).powi(2) < 0.25
                        };
                        if in_profile && (holes == 2 || q.x > 1.) {
                            expected_inside = false;
                        }
                        let qlocal = f.local(q);
                        let actual_inside = xs
                            .windows(2)
                            .position(|w| qlocal.x > w[0] && qlocal.x < w[1])
                            .is_some_and(|i| {
                                bs[i]
                                    .iter()
                                    .any(|b| Rect::at(b, qlocal.x).contains(qlocal.y, qlocal.z))
                            });
                        assert_eq!(actual_inside, expected_inside, "{case} probe {q:?}");
                    }
                }
            }
            let next = operation(&out, &placed(&tool, 1., p(0., 0., 6.)), 1e-6);
            assert_eq!(
                crate::topology::check(&next, 1e-6).unwrap().genus(),
                (1 + 2 * holes) as f64
            );
            assert!(
                (crate::mass::measure(&next, 1e-6).unwrap().volume
                    - (original_volume - 2. * removed))
                    .abs()
                    < 1e-8
            );
            for scale in [0.1_f64, 1., 25.4] {
                let offset = p(1e5, -1e5, 1e5);
                let moved = operation(
                    &placed(&base, scale, offset),
                    &placed(&tool, scale, offset),
                    scale * 1e-6,
                );
                assert!(
                    (crate::mass::measure(&moved, scale * 1e-6).unwrap().volume / scale.powi(3)
                        - expected)
                        .abs()
                        < 1e-6,
                    "{case} scale {scale}"
                );
            }
            let unit = |i| match i {
                0 => p(1., 0., 0.),
                1 => p(0., 1., 0.),
                _ => p(0., 0., 1.),
            };
            for w in 0..3 {
                for v in 0..3 {
                    if v == w {
                        continue;
                    }
                    let f = Frame {
                        u: unit(v).cross(unit(w)),
                        v: unit(v),
                        w: unit(w),
                    };
                    let rotated = operation(
                        &transform(&base, |q| f.world(q)),
                        &transform(&tool, |q| f.world(q)),
                        1e-6,
                    );
                    assert!(
                        (crate::mass::measure(&rotated, 1e-6).unwrap().volume - expected).abs()
                            < 1e-8
                    );
                }
            }
            let mut bad = base.clone();
            let cap = bad
                .faces
                .iter_mut()
                .find(|f| matches!(f.surface, Surface::Plane { .. }))
                .unwrap();
            cap.sense = !cap.sense;
            assert!(crate::boolean::cut(&bad, &tool, 1e-6).is_err());
        }
    }
    fn fixture(name: &str) -> Body {
        let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join(format!("tests/fixtures/circle_over_fillet/{name}.f64"));
        let data = std::fs::read(file).unwrap();
        let values: Vec<_> = data
            .chunks_exact(8)
            .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
            .collect();
        crate::neutral::Reader::new(&values).body().unwrap()
    }
    #[test]
    fn captured_fillet_circle_has_proven_interior_section_intervals() {
        let body = fixture("base");
        let tool = fixture("step01");
        validate(&body, 1e-6).unwrap();
        validate(&tool, 1e-6).unwrap();
        let mut xs = events(&body);
        xs.extend(events(&tool));
        let xs = stations(xs, 1e-6).unwrap();
        let mut bs = Vec::new();
        let mut ts = Vec::new();
        let mut split = 0;
        for x in xs.windows(2) {
            let b = sections(&body, x[0], x[1], 1e-6).unwrap();
            let t = sections(&tool, x[0], x[1], 1e-6).unwrap();
            let cut = subtract_sections(&b, &t, x[0], x[1], 1e-6, &mut 0).unwrap();
            split += usize::from(cut.len() > b.len());
            bs.push(b);
            ts.push(t);
        }
        verify_interfaces(&body, &xs, &bs, 1e-6).unwrap();
        verify_interfaces(&tool, &xs, &ts, 1e-6).unwrap();
        assert!(split >= 2);
    }

    #[test]
    fn captured_circle_crossing_a_fillet_returns_closed_analytic_boundary() {
        let body = fixture("base");
        let tool = fixture("step01");
        let result = cut(&body, &tool, 1e-6).unwrap();
        assert_eq!(crate::topology::check(&result, 1e-6).unwrap().genus(), 1.);
        assert!(result
            .edges
            .iter()
            .any(|e| matches!(e.curve, Curve::CylinderPair { .. })));
        let m = crate::mass::measure(&result, 1e-6).unwrap();
        assert!(
            (m.volume - reference_volume()).abs() < 1e-9,
            "{} expected {}",
            m.volume,
            reference_volume()
        );
        assert!(!m.exact);
        eprintln!(
            "fillet circle volume: {:.14}; reference {:.14}",
            m.volume,
            reference_volume()
        );
    }
    // Independent dimension integral over the disk: x=.74+.4*cos(theta).
    // This uses stock dimensions, not returned faces or section reconstruction.
    fn reference_volume() -> f64 {
        let n = 4096;
        let step = PI / n as f64;
        let integrand = |t: f64| {
            let width = 0.32 * t.sin().powi(2);
            let thickness = if t <= PI / 2. {
                0.23
            } else {
                1.23 - 2. * (0.25 - 0.16 * t.cos().powi(2)).sqrt()
            };
            width * thickness
        };
        let mut sum = integrand(0.) + integrand(PI);
        for i in 1..n {
            sum += if i % 2 == 0 { 2. } else { 4. } * integrand(i as f64 * step);
        }
        155.58846611769187 - sum * step / 3.
    }

    #[test]
    fn nonconcentric_fillet_crossings_use_geometry_not_fixture_dimensions() {
        for (cx, radius) in [(0.55_f64, 0.15_f64), (0.9, 0.3), (0.7, 0.35)] {
            let n = p(0., 0., 1.);
            let tool = crate::extrude::extrude(
                &crate::neutral::Profile {
                    loops: vec![vec![Curve::Arc {
                        c: p(cx, 12., 0.),
                        n,
                        x: p(1., 0., 0.),
                        r: radius,
                        t0: 0.,
                        t1: TAU,
                    }]],
                    n,
                    low: -6.,
                    high: 1.,
                },
                1e-6,
            )
            .unwrap();
            let result = operation(&fixture("base"), &tool, 1e-6);
            let integrand = |t: f64| {
                let x = cx + radius * t.cos();
                let height = if x >= 0.74 {
                    0.23
                } else {
                    1.23 - 2. * (0.25 - (x - 0.74).powi(2)).sqrt()
                };
                2. * radius * radius * t.sin().powi(2) * height
            };
            let split = ((0.74 - cx) / radius).clamp(-1., 1.).acos();
            let mut removed = 0.;
            for (a, b) in [(0., split), (split, PI)] {
                let n = 4096;
                let step = (b - a) / n as f64;
                let mut sum = integrand(a) + integrand(b);
                for i in 1..n {
                    sum += if i % 2 == 0 { 2. } else { 4. } * integrand(a + step * i as f64);
                }
                removed += sum * step / 3.;
            }
            let expected = 155.58846611769187 - removed;
            assert!(
                (crate::mass::measure(&result, 1e-6).unwrap().volume - expected).abs() < 1e-8,
                "cx {cx}, r {radius}"
            );
        }
    }

    #[test]
    fn large_circle_crosses_the_vertical_fillet_end_with_a_smooth_chart() {
        let n = p(0., 0., 1.);
        let radius = 0.7;
        let tool = crate::extrude::extrude(
            &crate::neutral::Profile {
                loops: vec![vec![Curve::Arc {
                    c: p(0.74, 12., 0.),
                    n,
                    x: p(1., 0., 0.),
                    r: radius,
                    t0: 0.,
                    t1: TAU,
                }]],
                n,
                low: -6.,
                high: 1.,
            },
            1e-6,
        )
        .unwrap();
        let result = operation(&fixture("base"), &tool, 1e-6);
        assert_eq!(crate::topology::check(&result, 1e-6).unwrap().genus(), 1.);
        assert!(result
            .edges
            .iter()
            .any(|e| matches!(e.curve,Curve::CylinderPair {a,..} if a.y.abs()==1.)));
        // Independent stock-section integral. Use the host circle angle in the
        // fillet strip, avoiding the square-root endpoint in x coordinates.
        let flange_area = -0.5 * (radius * radius - 0.25).sqrt()
            + radius * radius * (-0.5 / radius).asin()
            + radius * radius * PI / 2.;
        let web_area = radius * radius * PI / 2.;
        let f = |t: f64| {
            2. * (radius * radius - 0.25 * t.cos().powi(2)).sqrt()
                * 0.5
                * t.sin()
                * (1.23 - t.sin())
        };
        let steps = 4096;
        let h = PI / 2. / steps as f64;
        let mut sum = f(PI / 2.) + f(PI);
        for i in 1..steps {
            sum += if i % 2 == 0 { 2. } else { 4. } * f(PI / 2. + h * i as f64);
        }
        let expected = 155.58846611769187 - 4. * flange_area - 0.23 * web_area - sum * h / 3.;
        let actual = crate::mass::measure(&result, 1e-6).unwrap().volume;
        assert!(
            (actual - expected).abs() < 1e-8,
            "{actual} expected {expected}"
        );
        eprintln!("large fillet crossing volume {actual:.14}");
        // Independent stock dimensions across the entire flange depth, not only
        // the thin web. This detects shape errors that could cancel in volume.
        for (shape, is_cut) in [(&fixture("base"), false), (&result, true)] {
            let xs = stations(events(shape), 1e-6).unwrap();
            let bands: Vec<_> = xs
                .windows(2)
                .map(|w| sections(shape, w[0], w[1], 1e-6).unwrap())
                .collect();
            for ix in 0..31 {
                for iy in 0..37 {
                    for iz in 0..18 {
                        let q = p(
                            0.035 + ix as f64 * 0.047,
                            11.313 + iy as f64 * 0.039,
                            -3.95 + iz as f64 * 0.233,
                        );
                        let (z0, z1) = if q.x < 0.24 {
                            (-4., 0.)
                        } else if q.x < 0.74 {
                            let root = (0.25 - (q.x - 0.74).powi(2)).sqrt();
                            (-2.615 + root, -1.385 - root)
                        } else {
                            (-2.115, -1.885)
                        };
                        let expected = q.z > z0
                            && q.z < z1
                            && (!is_cut
                                || (q.x - 0.74).powi(2) + (q.y - 12.).powi(2) > radius * radius);
                        let inside = xs
                            .windows(2)
                            .position(|w| q.x > w[0] && q.x < w[1])
                            .is_some_and(|i| {
                                bands[i].iter().any(|b| Rect::at(b, q.x).contains(q.y, q.z))
                            });
                        assert_eq!(inside, expected, "large cut={is_cut} probe {q:?}");
                    }
                }
            }
        }
        let unit = |i| match i {
            0 => p(1., 0., 0.),
            1 => p(0., 1., 0.),
            _ => p(0., 0., 1.),
        };
        for w in 0..3 {
            for v in 0..3 {
                if w == v {
                    continue;
                }
                let frame = Frame {
                    u: unit(v).cross(unit(w)),
                    v: unit(v),
                    w: unit(w),
                };
                let rotated = operation(
                    &transform(&fixture("base"), |q| frame.world(q)),
                    &transform(&tool, |q| frame.world(q)),
                    1e-6,
                );
                assert!(
                    (crate::mass::measure(&rotated, 1e-6).unwrap().volume - expected).abs() < 1e-8
                );
            }
        }
        for scale in [0.1_f64, 1., 25.4] {
            let offset = p(1e5, -1e5, 1e5);
            let moved = operation(
                &placed(&fixture("base"), scale, offset),
                &placed(&tool, scale, offset),
                scale * 1e-6,
            );
            assert!(
                (crate::mass::measure(&moved, scale * 1e-6).unwrap().volume / scale.powi(3)
                    - expected)
                    .abs()
                    < 1e-6
            );
        }
        let second = operation(&result, &placed(&tool, 1., p(0., 10., 0.)), 1e-6);
        assert_eq!(crate::topology::check(&second, 1e-6).unwrap().genus(), 2.);
        assert!(
            (crate::mass::measure(&second, 1e-6).unwrap().volume
                - (2. * expected - 155.58846611769187))
                .abs()
                < 1e-8
        );
    }

    fn operation(body: &Body, tool: &Body, tol: f64) -> Body {
        let mut data = crate::neutral::write_body(body);
        data.extend(crate::neutral::write_body(tool));
        let wire = crate::protocol::encode_request(crate::protocol::OP_CUT, tol, &data);
        let request = crate::protocol::parse_request(&wire).unwrap();
        let out = crate::protocol::handle(&request).unwrap();
        let mut reader = crate::neutral::Reader::new(&out);
        let result = reader.body().unwrap();
        reader.finish().unwrap();
        result
    }
    fn placed(body: &Body, scale: f64, offset: Vec3) -> Body {
        let mut out = body.clone();
        for v in &mut out.vertices {
            *v = v.scale(scale).add(offset);
        }
        for e in &mut out.edges {
            match &mut e.curve {
                Curve::Line { a, b } => {
                    *a = a.scale(scale).add(offset);
                    *b = b.scale(scale).add(offset);
                }
                Curve::Arc { c, r, .. } => {
                    *c = c.scale(scale).add(offset);
                    *r *= scale;
                }
                Curve::Conic { c, u, v, .. } => {
                    *c = c.scale(scale).add(offset);
                    *u = u.scale(scale);
                    *v = v.scale(scale);
                }
                Curve::CylinderPair { c, p, rho, r, .. } => {
                    *c = c.scale(scale).add(offset);
                    *p = p.scale(scale).add(offset);
                    *rho *= scale;
                    *r *= scale;
                }
                Curve::Bezier { .. } | Curve::RationalBezier { .. } | Curve::RadialPair { .. } | Curve::ConePair { .. } => unreachable!("no Bezier or radial-pair edges in these test bodies"),
            }
        }
        for f in &mut out.faces {
            match &mut f.surface {
                Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => unreachable!("no elliptic-cylinder or revolved faces in these test bodies"),
                Surface::Plane { p, .. } => *p = p.scale(scale).add(offset),
                Surface::Cylinder { p, r, .. } => {
                    *p = p.scale(scale).add(offset);
                    *r *= scale;
                }
            }
        }
        out
    }
    #[test]
    fn separated_fillet_circles_replay_through_the_worker_contract() {
        let base = fixture("base");
        let tool = fixture("step01");
        let saved = crate::neutral::write_body(&base);
        let first = operation(&base, &tool, 1e-6);
        let second = operation(&first, &placed(&tool, 1., p(0., 10., 0.)), 1e-6);
        assert_eq!(crate::topology::check(&second, 1e-6).unwrap().genus(), 2.);
        let expected = 2. * reference_volume() - 155.58846611769187;
        assert!((crate::mass::measure(&second, 1e-6).unwrap().volume - expected).abs() < 1e-8);
        assert_eq!(crate::neutral::write_body(&base), saved);
    }
    #[test]
    fn fillet_circle_survives_all_six_cardinal_axis_frames() {
        let unit = |i| match i {
            0 => p(1., 0., 0.),
            1 => p(0., 1., 0.),
            _ => p(0., 0., 1.),
        };
        for w in 0..3 {
            for v in 0..3 {
                if w == v {
                    continue;
                }
                let f = Frame {
                    u: unit(v).cross(unit(w)),
                    v: unit(v),
                    w: unit(w),
                };
                let body = transform(&fixture("base"), |p| f.world(p));
                let tool = transform(&fixture("step01"), |p| f.world(p));
                let result = operation(&body, &tool, 1e-6);
                assert!(
                    (crate::mass::measure(&result, 1e-6).unwrap().volume - reference_volume())
                        .abs()
                        < 1e-8
                );
            }
        }
    }
    #[test]
    fn fillet_circle_preserves_scaled_and_translated_geometry() {
        for (scale, offset) in [
            (0.1, p(100., -50., 15.)),
            (25.4, p(100., -50., 15.)),
            (1., p(1e5, -1e5, 1e5)),
        ] {
            let body = placed(&fixture("base"), scale, offset);
            let tool = placed(&fixture("step01"), scale, offset);
            eprintln!("scaled fillet cut: scale {scale}, offset {offset:?}");
            let result = operation(&body, &tool, scale * 1e-6);
            eprintln!("scaled fillet cut returned; measuring boundary");
            let measured =
                crate::mass::measure(&result, scale * 1e-6).unwrap().volume / scale.powi(3);
            assert!(
                (measured - reference_volume()).abs() < 1e-6,
                "{scale} {offset:?}: {measured}"
            );
        }
    }
    #[test]
    fn translated_fillet_geometry_passes_independently_of_measurement() {
        for scale in [0.1, 1., 25.4] {
            let offset = p(1e5, -1e5, 1e5);
            let body = placed(&fixture("base"), scale, offset);
            let tool = placed(&fixture("step01"), scale, offset);
            let result = operation(&body, &tool, scale * 1e-6);
            assert_eq!(
                crate::topology::check(&result, scale * 1e-6)
                    .unwrap()
                    .genus(),
                1.
            );
        }
    }

    #[test]
    fn singular_and_malformed_crossing_cuts_leave_inputs_untouched() {
        let base = fixture("base");
        let original = crate::neutral::write_body(&base);
        // A concentric radius equal to the root fillet reaches the flange face at the fillets' tangent line: hole and
        // fillets touch there. Claude 2026-10-28 (Bob's rule: within tolerance it touches) admits it; removed, by an
        // independent integral, .23 pi / 4 + 2 (pi / 16 - 1 / 6) = .1825 pi - 1 / 3. Within 1e-7, not 1e-9: the
        // fillet's tangent line (x .24000000000000044) and the hole's extreme (.24) are one station, and the hole's
        // angle at it is acos(-1 + 9e-16) = pi - 4e-8, so its edges end 2e-8 short of the exact touch point (the
        // kernel's points stay within 5e-16 of both cylinders). Bob's captured hole, whose extreme is the station,
        // matches its integral within 1e-13 (tests/touch.rs).
        let n = p(0., 0., 1.);
        let tangent = crate::extrude::extrude(
            &crate::neutral::Profile {
                loops: vec![vec![Curve::Arc {
                    c: p(0.74, 5., 0.),
                    n,
                    x: p(1., 0., 0.),
                    r: 0.5,
                    t0: 0.,
                    t1: TAU,
                }]],
                n,
                low: -6.,
                high: 1.,
            },
            1e-6,
        )
        .unwrap();
        let touching = crate::boolean::cut(&base, &tangent, 1e-6).unwrap();
        crate::topology::check(&touching, 1e-6).unwrap();
        let removed = crate::mass::measure(&base, 1e-6).unwrap().volume
            - crate::mass::measure(&touching, 1e-6).unwrap().volume;
        assert!((removed - (0.1825 * PI - 1. / 3.)).abs() < 1e-7, "{removed}");
        let mut malformed = operation(&base, &fixture("step01"), 1e-6);
        let e = malformed
            .edges
            .iter_mut()
            .find(|e| matches!(e.curve, Curve::CylinderPair { .. }))
            .unwrap();
        if let Curve::CylinderPair { branch, .. } = &mut e.curve {
            *branch = 0.;
        }
        assert!(matches!(
            crate::boolean::cut(&malformed, &fixture("step01"), 1e-6),
            Err(KernelError::InvalidInput(_))
        ));
        assert_eq!(crate::neutral::write_body(&base), original);
    }

    #[test]
    fn capsule_crossing_a_root_fillet_retains_analytic_intersections() {
        let base = fixture("base");
        let n = p(0., 0., 1.);
        let radius = 0.4;
        // A Y-directed capsule: lower/upper semicircle centers (.74,5),(.74,7).
        let outline = vec![
            Curve::Arc {
                c: p(0.74, 5., 0.),
                n,
                x: p(-1., 0., 0.),
                r: radius,
                t0: 0.,
                t1: PI,
            },
            segment(p(1.14, 5., 0.), p(1.14, 7., 0.)),
            Curve::Arc {
                c: p(0.74, 7., 0.),
                n,
                x: p(1., 0., 0.),
                r: radius,
                t0: 0.,
                t1: PI,
            },
            segment(p(0.34, 7., 0.), p(0.34, 5., 0.)),
        ];
        let tool = crate::extrude::extrude(
            &crate::neutral::Profile {
                loops: vec![outline],
                n,
                low: -6.,
                high: 1.,
            },
            1e-6,
        )
        .unwrap();
        let result = operation(&base, &tool, 1e-6);
        assert_eq!(crate::topology::check(&result, 1e-6).unwrap().genus(), 1.);
        // Capsule ends together remove one circle. The intervening 2-inch
        // rectangle crosses .4 of fillet and .4 of web; integrate sqrt exactly.
        let circular_integral = 0.5 * (0.4 * 0.3 + 0.25 * (0.4_f64 / 0.5).asin());
        let removed_strip = 2. * (0.4 * 1.23 - 2. * circular_integral + 0.4 * 0.23);
        let expected = reference_volume() - removed_strip;
        assert!((crate::mass::measure(&result, 1e-6).unwrap().volume - expected).abs() < 1e-8);
    }

    #[test]
    fn rat_holes_fillet_circle_web_slot_and_bore_preserve_all_prior_cuts() {
        let mut body = fixture("base");
        for name in ["step01", "step02"] {
            let bytes = std::fs::read(
                std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join(format!("tests/fixtures/rat_hole/{name}.f64")),
            )
            .unwrap();
            let values: Vec<_> = bytes
                .chunks_exact(8)
                .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
                .collect();
            body = operation(
                &body,
                &crate::neutral::Reader::new(&values).body().unwrap(),
                1e-6,
            );
        }
        body = operation(&body, &fixture("step01"), 1e-6);
        let n = p(0., 0., 1.);
        let slot = vec![
            segment(p(3., 3., 0.), p(7., 3., 0.)),
            Curve::Arc {
                c: p(7., 4., 0.),
                n,
                x: p(0., -1., 0.),
                r: 1.,
                t0: 0.,
                t1: PI,
            },
            segment(p(7., 5., 0.), p(3., 5., 0.)),
            Curve::Arc {
                c: p(3., 4., 0.),
                n,
                x: p(0., 1., 0.),
                r: 1.,
                t0: 0.,
                t1: PI,
            },
        ];
        let hole = vec![Curve::Arc {
            c: p(5., 14., 0.),
            n,
            x: p(1., 0., 0.),
            r: 0.5,
            t0: 0.,
            t1: TAU,
        }];
        for outline in [slot, hole] {
            let tool = crate::extrude::extrude(
                &crate::neutral::Profile {
                    loops: vec![outline],
                    n,
                    low: -6.,
                    high: 1.,
                },
                1e-6,
            )
            .unwrap();
            body = operation(&body, &tool, 1e-6);
        }
        let expected = 153.57573395171957
            - (155.58846611769187 - reference_volume())
            - 0.23 * (8. + PI + PI * 0.25);
        assert_eq!(crate::topology::check(&body, 1e-6).unwrap().genus(), 3.);
        assert!((crate::mass::measure(&body, 1e-6).unwrap().volume - expected).abs() < 1e-8);
    }

    #[test]
    fn captured_fillet_circle_matches_every_production_material_probe() {
        let result = cut(&fixture("base"), &fixture("step01"), 1e-6).unwrap();
        let xs = stations(events(&result), 1e-6).unwrap();
        let bands: Vec<_> = xs
            .windows(2)
            .map(|w| sections(&result, w[0], w[1], 1e-6).unwrap())
            .collect();
        let manifest = std::fs::read_to_string(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures/circle_over_fillet/manifest.txt"),
        )
        .unwrap();
        let mut checked = 0;
        for line in manifest.lines() {
            let Some((coordinates, label)) = line.trim().split_once(": ") else {
                continue;
            };
            if label != "Inside" && label != "Outside" {
                continue;
            }
            let numbers: Vec<f64> = coordinates
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if numbers.len() != 3 {
                continue;
            }
            let q = p(numbers[0], numbers[1], numbers[2]);
            let interval = xs.windows(2).position(|w| q.x > w[0] && q.x < w[1]);
            let inside = interval
                .is_some_and(|i| bands[i].iter().any(|b| Rect::at(b, q.x).contains(q.y, q.z)));
            assert_eq!(inside, label == "Inside", "probe {line}");
            checked += 1;
        }
        assert_eq!(checked, 315);
    }
    #[test]
    fn interior_rat_hole_reconnects_material_around_the_opening() {
        let body = fixture("base");
        let bytes = std::fs::read(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures/rat_hole/step01.f64"),
        )
        .unwrap();
        let values: Vec<_> = bytes
            .chunks_exact(8)
            .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
            .collect();
        let mut tool = crate::neutral::Reader::new(&values).body().unwrap();
        let offset = p(0., 5., 0.);
        for v in &mut tool.vertices {
            *v = v.add(offset);
        }
        for e in &mut tool.edges {
            e.curve = e.curve.translated(offset);
        }
        for f in &mut tool.faces {
            match &mut f.surface {
                Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => unreachable!("no elliptic-cylinder or revolved faces in these test bodies"),
                Surface::Plane { p, .. } | Surface::Cylinder { p, .. } => *p = p.add(offset),
            }
        }
        let result = cut(&body, &tool, 1e-6).unwrap();
        assert_eq!(crate::topology::check(&result, 1e-6).unwrap().genus(), 1.);
        assert!(
            (crate::mass::measure(&result, 1e-6).unwrap().volume - 154.58210003470572).abs() < 1e-8
        );
        let offset = p(0., 10., 0.);
        for v in &mut tool.vertices {
            *v = v.add(offset);
        }
        for e in &mut tool.edges {
            e.curve = e.curve.translated(offset);
        }
        for f in &mut tool.faces {
            match &mut f.surface {
                Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } | Surface::Spline { .. } | Surface::Conoid { .. } => unreachable!("no elliptic-cylinder or revolved faces in these test bodies"),
                Surface::Plane { p, .. } | Surface::Cylinder { p, .. } => *p = p.add(offset),
            }
        }
        let wire = crate::neutral::write_body(&result);
        let restored = crate::neutral::Reader::new(&wire).body().unwrap();
        let second = cut(&restored, &tool, 1e-6).unwrap();
        assert_eq!(crate::topology::check(&second, 1e-6).unwrap().genus(), 2.);
        assert!(
            (crate::mass::measure(&second, 1e-6).unwrap().volume - 153.57573395171957).abs() < 1e-8
        );
    }
}

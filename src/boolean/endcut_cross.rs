//! Multiple transverse material intervals for through cuts crossing cylinders.
//! Shares analytic boundary construction with the end-cut domain. Input sections
//! are read from the current B-rep; no prior stock or operation recipe is used.
use super::*;

fn positive(lower: Graph, upper: Graph, lo: f64, hi: f64, tol: f64) -> bool {
    if below(lower, upper, lo, hi, -tol) {
        return true;
    }
    // A circle closes to zero width at its extremal stations but has positive
    // width throughout the open band. This is a proved analytic degeneracy.
    match (lower, upper) {
        (
            Graph::Circle {
                cx,
                cv,
                r,
                branch: -1.,
            },
            Graph::Circle {
                cx: bx,
                cv: bv,
                r: br,
                branch: 1.,
            },
        ) if cx == bx && cv == bv && r == br => {
            lo >= cx - r - tol
                && hi <= cx + r + tol
                && upper.value((lo + hi) * 0.5) - lower.value((lo + hi) * 0.5) > tol
        }
        _ => false,
    }
}

fn sections(body: &Body, lo: f64, hi: f64, tol: f64) -> Result<Vec<Band>> {
    let mid = lo + (hi - lo) * 0.5;
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
    let mut used = vec![false; hits.len()];
    let mut result = Vec::new();
    for (i, low) in hits.iter().enumerate().filter(|(_, h)| h.role == 0) {
        let yl = low.graph.value(mid);
        let matches: Vec<_> = hits
            .iter()
            .enumerate()
            .filter(|(_, h)| {
                h.role == 1
                    && h.graph.value(mid) > yl + tol
                    && (h.lo - low.lo).abs() <= tol
                    && (h.hi - low.hi).abs() <= tol
            })
            .collect();
        let (j, high) = matches
            .into_iter()
            .min_by(|(_, a), (_, b)| a.graph.value(mid).total_cmp(&b.graph.value(mid)))
            .ok_or_else(|| invalid("transverse section has no matching upper wall"))?;
        let yh = high.graph.value(mid);
        let mut roofs = Vec::new();
        for (role, z) in [(2, low.lo), (3, low.hi)] {
            let candidates: Vec<_> = hits
                .iter()
                .enumerate()
                .filter(|(_, h)| {
                    h.role == role
                        && (h.graph.value(mid) - z).abs() <= tol
                        && (h.lo - yl).abs() <= tol
                        && (h.hi - yh).abs() <= tol
                })
                .collect();
            if candidates.len() != 1 {
                return Err(unsupported("transverse section roof is split or ambiguous"));
            }
            roofs.push(candidates[0]);
        }
        for k in [i, j, roofs[0].0, roofs[1].0] {
            if used[k] {
                return Err(invalid("overlapping transverse material sections"));
            }
            used[k] = true;
        }
        let b = Band {
            lo,
            hi,
            y0: low.graph,
            y1: high.graph,
            z0: roofs[0].1.graph,
            z1: roofs[1].1.graph,
        };
        if !positive(b.y0, b.y1, lo, hi, tol) || !positive(b.z0, b.z1, lo, hi, tol) {
            return Err(unsupported(
                "cannot prove positive transverse material thickness",
            ));
        }
        result.push(b);
    }
    if used.iter().any(|u| !u) {
        return Err(invalid("unpaired transverse boundary segments"));
    }
    result.sort_by(|a, b| a.y0.value(mid).total_cmp(&b.y0.value(mid)));
    for i in 0..result.len() {
        for j in i + 1..result.len() {
            let (a, b) = (&result[i], &result[j]);
            if !(below(a.y1, b.y0, lo, hi, 0.)
                || below(b.y1, a.y0, lo, hi, 0.)
                || below(a.z1, b.z0, lo, hi, 0.)
                || below(b.z1, a.z0, lo, hi, 0.))
            {
                return Err(unsupported(
                    "transverse material intervals must remain separated",
                ));
            }
        }
    }
    Ok(result)
}

fn subtract_interval(body: Band, tool: &Band, tol: f64) -> Result<Vec<Band>> {
    let (lo, hi) = (body.lo, body.hi);
    if below(tool.y1, body.y0, lo, hi, 0.)
        || below(body.y1, tool.y0, lo, hi, 0.)
        || below(tool.z1, body.z0, lo, hi, 0.)
        || below(body.z1, tool.z0, lo, hi, 0.)
    {
        return Ok(vec![body]);
    }
    if !below(tool.z0, body.z0, lo, hi, 0.) || !below(body.z1, tool.z1, lo, hi, 0.) {
        return Err(unsupported(
            "transverse cutter must cover the entire material depth",
        ));
    }
    let covers_low = below(tool.y0, body.y0, lo, hi, 0.);
    let covers_high = below(body.y1, tool.y1, lo, hi, 0.);
    let inner_low = below(body.y0, tool.y0, lo, hi, -tol) && below(tool.y0, body.y1, lo, hi, -tol);
    let inner_high = below(body.y0, tool.y1, lo, hi, -tol) && below(tool.y1, body.y1, lo, hi, -tol);
    let mut result = Vec::new();
    if !(covers_low || inner_low) || !(covers_high || inner_high) {
        return Err(unsupported(
            "cutter crosses a material boundary inside a transverse band",
        ));
    }
    if inner_low {
        let mut b = body.clone();
        b.y1 = tool.y0;
        result.push(b);
    }
    if inner_high {
        let mut b = body;
        b.y0 = tool.y1;
        result.push(b);
    }
    Ok(result)
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    cut_in_frame(body, tool, tol, frame(body, tool)?)
}
pub(super) fn cut_along(body: &Body, tool: &Body, axis: usize, tol: f64) -> Result<Body> {
    cut_in_frame(body, tool, tol, frame_for_axis(body, axis)?)
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
                Surface::Plane { p, .. } | Surface::Cylinder { p, .. } => *p = p.add(by),
            }
        }
        out
    };
    let body = transform(&shift(body, origin.scale(-1.)), |p| f.local(p));
    let tool = transform(&shift(tool, origin.scale(-1.)), |p| f.local(p));
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
            if let Some((_, u, v, t0, t1)) = section_projection(&edge.curve) {
                if u.x.hypot(v.x) == 0. {
                    continue;
                }
                let phase = v.x.atan2(u.x);
                for angle in [phase, phase + PI] {
                    for t in parameters(angle, t0, t1) {
                        extrema.push(eval(&edge.curve, t).x);
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
    let mut original = Vec::new();
    let mut cutters = Vec::new();
    let mut bands = Vec::new();
    for x in xs.windows(2) {
        let b = sections(&body, x[0], x[1], tol)?;
        let t = sections(&tool, x[0], x[1], tol)?;
        let mut remaining = b.clone();
        for cutter in &t {
            let mut next = Vec::new();
            for stock in remaining {
                next.extend(subtract_interval(stock, cutter, tol)?);
            }
            if next.len() > 128 {
                return Err(unsupported("too many transverse material intervals"));
            }
            remaining = next;
        }
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
    for b in bands.iter().flatten() {
        result.extend(patches(b, tol)?);
    }
    for i in 0..xs.len() {
        let left = if i > 0 { bands[i - 1].as_slice() } else { &[] };
        let right = bands.get(i).map(Vec::as_slice).unwrap_or(&[]);
        result.extend(end_patches(left, right, xs[i], tol)?);
    }
    let result = assemble(result, tol)?;
    Ok(shift(&transform(&result, |p| f.world(p)), origin))
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let case = "rect_tube_hole_over_corner";
        assert!(matches!(
            crate::boolean::cut(&survey(case, "base"), &survey(case, "step01"), 1e-6),
            Err(KernelError::Unsupported(_))
        ));
        let base = survey("pipe_hole_one_wall", "base");
        let original = crate::neutral::write_body(&base);
        let blind = placed(&survey("pipe_hole_one_wall", "step01"), 1., p(0.9, 0., 0.));
        assert!(matches!(
            crate::boolean::cut(&base, &blind, 1e-6),
            Err(KernelError::Unsupported(_))
        ));
        assert_eq!(crate::neutral::write_body(&base), original);
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
            for b in &b {
                for t in &t {
                    let cut = subtract_interval(b.clone(), t, 1e-6).unwrap();
                    split += usize::from(cut.len() == 2);
                }
            }
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
            }
        }
        for f in &mut out.faces {
            match &mut f.surface {
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
    fn singular_blind_and_malformed_crossing_cuts_leave_inputs_untouched() {
        let base = fixture("base");
        let original = crate::neutral::write_body(&base);
        // A concentric radius equal to the root fillet meets a branch singularity.
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
        assert!(matches!(
            crate::boolean::cut(&base, &tangent, 1e-6),
            Err(KernelError::Unsupported(_))
        ));
        let blind = placed(&fixture("step01"), 1., p(0., 0., 4.));
        assert!(matches!(
            crate::boolean::cut(&base, &blind, 1e-6),
            Err(KernelError::Unsupported(_))
        ));
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

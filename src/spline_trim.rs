//! Certified polygonal trim domains on tensor surfaces with an affine inverse
//! chart. The chart is established from the entire homogeneous control net;
//! edge/support agreement uses Bernstein composition, never point sampling.
use crate::{
    geom::{Curve, Vec3},
    neutral::{Body, Face},
    spline_surface::SplineSurface,
    KernelError,
};
type Result<T> = std::result::Result<T, KernelError>;
fn reject<T>(s: &str) -> Result<T> {
    Err(KernelError::Unsupported(s.into()))
}

#[derive(Clone, Copy)]
pub(crate) struct Chart {
    o: Vec3,
    u: Vec3,
    v: Vec3,
}
impl Chart {
    pub(crate) fn point(&self, p: Vec3) -> [f64; 2] {
        let d = p.sub(self.o);
        [d.dot(self.u), d.dot(self.v)]
    }
}
fn choose(n: usize, k: usize) -> f64 {
    (0..k.min(n - k)).fold(1., |a, i| a * (n - i) as f64 / (i + 1) as f64)
}
fn product(a: &[f64], b: &[f64]) -> Vec<f64> {
    let (n, m) = (a.len() - 1, b.len() - 1);
    let mut out = vec![0.; n + m + 1];
    for i in 0..=n {
        for j in 0..=m {
            out[i + j] += a[i] * b[j] * choose(n, i) * choose(m, j) / choose(n + m, i + j);
        }
    }
    out
}
fn cubic_basis(a: &[f64], h: &[f64]) -> [Vec<f64>; 4] {
    let b: Vec<_> = h.iter().zip(a).map(|(h, a)| h - a).collect();
    std::array::from_fn(|i| {
        let mut q = vec![choose(3, i)];
        for k in 0..3 {
            q = product(&q, if k < i { a } else { &b });
        }
        q
    })
}

/// An affine inverse exists only if both homogeneous identities hold across
/// every patch. Near rank-deficiency and inconclusive roundoff are refused.
pub(crate) fn chart(s: &SplineSurface) -> Option<Chart> {
    let o = s.eval(0., 0.).0;
    let a = s.eval(1., 0.).0.sub(o);
    let b = s.eval(0., 1.).0.sub(o);
    let n = a.cross(b);
    let n2 = n.dot(n);
    if !n2.is_finite() || n2 <= 1e-28 {
        return None;
    }
    let mut c = Chart {
        o,
        u: b.cross(n).scale(1. / n2),
        v: n.cross(a).scale(1. / n2),
    };
    let mut third = (0., Vec3::new(0., 0., 0.), [0., 0.]);
    for (k, p) in s.patches.iter().enumerate() {
        let (i, j) = (k / (s.v.len() - 1), k % (s.v.len() - 1));
        // Only constant-weight nets supply control-point parameter coordinates
        // for selecting the third basis vector. The final identity is general.
        if p.w.iter().any(|w| *w != p.w[0]) {
            continue;
        }
        for x in 0..4 {
            for y in 0..4 {
                let d = p.p[x * 4 + y].sub(o);
                let det = d.dot(n);
                if det.abs() > third.0 {
                    third = (
                        det.abs(),
                        d,
                        [
                            s.u[i] + (s.u[i + 1] - s.u[i]) * x as f64 / 3.,
                            s.v[j] + (s.v[j + 1] - s.v[j]) * y as f64 / 3.,
                        ],
                    );
                }
            }
        }
    }
    if third.0 > 1e-10 * n.length() * a.length().max(b.length()) {
        let d = third.1;
        let det = n.dot(d);
        c.u = b.cross(d).add(n.scale(third.2[0])).scale(1. / det);
        c.v = d.cross(a).add(n.scale(third.2[1])).scale(1. / det);
    }
    if !c.u.is_finite() || !c.v.is_finite() {
        return None;
    }
    let cond = (c.u.length() + c.v.length()) * a.length().max(b.length());
    if cond > 1e6 {
        return None;
    }
    for (k, p) in s.patches.iter().enumerate() {
        let (i, j) = (k / (s.v.len() - 1), k % (s.v.len() - 1));
        let wmax = p.w.iter().copied().fold(0., f64::max);
        let wmin = p.w.iter().copied().fold(f64::INFINITY, f64::min) / wmax;
        for x in 0usize..5 {
            for y in 0usize..5 {
                let mut residual = [0.; 2];
                for a in x.saturating_sub(1)..=x.min(3) {
                    for b in y.saturating_sub(1)..=y.min(3) {
                        let t = c.point(p.p[a * 4 + b]);
                        let factor = choose(3, a) * choose(3, b) / (choose(4, x) * choose(4, y))
                            * p.w[a * 4 + b]
                            / wmax;
                        residual[0] +=
                            factor * (t[0] - (if x - a == 0 { s.u[i] } else { s.u[i + 1] }));
                        residual[1] +=
                            factor * (t[1] - (if y - b == 0 { s.v[j] } else { s.v[j + 1] }));
                    }
                }
                if residual
                    .iter()
                    .any(|r| !r.is_finite() || r.abs() / wmin > 256. * f64::EPSILON * cond.max(1.))
                {
                    return None;
                }
            }
        }
    }
    Some(c)
}

fn cross(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}
fn inside(p: [f64; 2], lp: &[[f64; 2]]) -> bool {
    let mut result = false;
    for i in 0..lp.len() {
        let (a, b) = (lp[i], lp[(i + 1) % lp.len()]);
        if (a[1] > p[1]) != (b[1] > p[1])
            && p[0] < a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
        {
            result = !result;
        }
    }
    result
}

/// Positive-weight rational curve composed with one tensor patch, in relative
/// coordinates. The numerator of the difference has degree 21. Its convex hull
/// bounds the WHOLE edge, including excursions between evaluation points.
fn residual(s: &SplineSurface, c: &Chart, curve: &Curve, i: usize, j: usize) -> Result<f64> {
    let (p, mut w) = crate::spline_curve::controls(curve)?;
    let max = w.iter().copied().fold(0., f64::max);
    w = w.map(|x| x / max);
    let uv = p.map(|p| c.point(p));
    let u: Vec<_> = (0..4)
        .map(|k| w[k] * (uv[k][0] - s.u[i]) / (s.u[i + 1] - s.u[i]))
        .collect();
    let v: Vec<_> = (0..4)
        .map(|k| w[k] * (uv[k][1] - s.v[j]) / (s.v[j + 1] - s.v[j]))
        .collect();
    let (bu, bv) = (cubic_basis(&u, &w), cubic_basis(&v, &w));
    let patch = &s.patches[i * (s.v.len() - 1) + j];
    let pwmax = patch.w.iter().copied().fold(0., f64::max);
    let mut d = vec![0.; 19];
    let mut q = [vec![0.; 19], vec![0.; 19], vec![0.; 19]];
    let mut magnitude = 0.;
    for x in 0..4 {
        for y in 0..4 {
            let coeff = product(&bu[x], &bv[y]);
            let r = patch.p[x * 4 + y].sub(p[0]);
            let r = [r.x, r.y, r.z];
            for k in 0..19 {
                let v = coeff[k] * patch.w[x * 4 + y] / pwmax;
                d[k] += v;
                for n in 0..3 {
                    q[n][k] += v * r[n];
                    magnitude += v.abs() * r[n].abs();
                }
            }
        }
    }
    let least =
        d.iter().copied().fold(f64::INFINITY, f64::min) * w.iter().copied().fold(1., f64::min);
    if least <= 0. {
        return Ok(f64::INFINITY);
    }
    let world: Vec<_> = p.iter().map(|q| q.sub(p[0])).collect();
    let mut r = vec![[0.; 3]; 22];
    for n in 0..3 {
        let a = product(&q[n], &w);
        let b = product(
            &d,
            &(0..4)
                .map(|k| w[k] * [world[k].x, world[k].y, world[k].z][n])
                .collect::<Vec<_>>(),
        );
        for k in 0..22 {
            r[k][n] = a[k] - b[k];
        }
    }
    let bound = r
        .iter()
        .map(|r| r[0].hypot(r[1]).hypot(r[2]))
        .fold(0., f64::max);
    Ok((bound + 8192. * f64::EPSILON * magnitude.max(1.)) / least)
}

fn polygon_domain(
    s: &SplineSurface,
    c: &Chart,
    body: &Body,
    face: &Face,
    tol: f64,
) -> Result<Vec<Vec<[f64; 2]>>> {
    if face.loops.is_empty() || face.loops.iter().map(Vec::len).sum::<usize>() > 2048 {
        return reject("spline trim loop count is empty or exceeds 2048 edges");
    }
    let mut loops = vec![];
    for lp in &face.loops {
        if lp.len() < 3 {
            return reject("polygonal spline trim needs at least three edges");
        }
        let mut points = vec![];
        for &(e, fwd) in lp {
            let curve = if fwd {
                body.edges[e].curve.clone()
            } else {
                body.edges[e].curve.reversed()
            };
            let (p, w) = crate::spline_curve::controls(&curve)?;
            crate::rational::validate(p, w)?;
            let uv = p.map(|p| c.point(p));
            let (a, b) = (uv[0], uv[3]);
            let d = [b[0] - a[0], b[1] - a[1]];
            let len = d[0].hypot(d[1]);
            if len < 1e-12
                || uv
                    .iter()
                    .flatten()
                    .any(|x| !x.is_finite() || *x < -1e-12 || *x > 1. + 1e-12)
            {
                return reject("spline trim is degenerate or outside the support domain");
            }
            let mut last = -1e-12;
            for q in uv {
                let t = ((q[0] - a[0]) * d[0] + (q[1] - a[1]) * d[1]) / (len * len);
                if cross(a, b, q).abs() > 1e-12 * len || t < last - 1e-12 || t > 1. + 1e-12 {
                    return reject("spline trim requires a monotone straight UV edge; curved UV trimming is not yet supported");
                }
                last = t;
            }
            let mut breaks = vec![0., 1.];
            for (axis, knots) in [&s.u, &s.v].iter().enumerate() {
                if d[axis].abs() < 1e-14 {
                    continue;
                }
                for &k in knots.iter().skip(1).take(knots.len() - 2) {
                    if k <= a[axis].min(b[axis]) || k >= a[axis].max(b[axis]) {
                        continue;
                    }
                    let (mut l, mut h) = (0., 1.);
                    for _ in 0..55 {
                        let m = (l + h) / 2.;
                        if (c.point(curve.eval(m).0)[axis] < k) == (d[axis] > 0.) {
                            l = m;
                        } else {
                            h = m;
                        }
                    }
                    breaks.push((l + h) / 2.);
                }
            }
            breaks.sort_by(f64::total_cmp);
            breaks.dedup();
            for pair in breaks.windows(2) {
                let m = c.point(curve.eval((pair[0] + pair[1]) / 2.).0);
                let i =
                    s.u.partition_point(|u| *u <= m[0])
                        .saturating_sub(1)
                        .min(s.u.len() - 2);
                let j =
                    s.v.partition_point(|v| *v <= m[1])
                        .saturating_sub(1)
                        .min(s.v.len() - 2);
                let part = crate::spline_curve::piece(&curve, pair[0], pair[1])?;
                let error = residual(s, c, &part, i, j)?;
                if !error.is_finite() || error > tol {
                    return Err(KernelError::Unsupported(format!("spline trim edge {e} support deviation bound {error:.9} inches exceeds limit {tol:.9}")));
                }
            }
            points.push(a);
        }
        loops.push(points);
    }
    polygon_checks(&loops, if face.sense { 1. } else { -1. })?;
    Ok(loops)
}

pub(crate) fn polygon_checks(loops: &[Vec<[f64; 2]>], sign: f64) -> Result<()> {
    if loops.is_empty()
        || loops.iter().any(|p| p.len() < 3)
        || loops.iter().map(Vec::len).sum::<usize>() > 8192
    {
        return reject("spline trim polygon is empty, degenerate, or exceeds 8192 points");
    }
    // Conservative geometric predicates reject touching or nearly unresolved
    // nonadjacent segments rather than assigning a guessed trim topology.
    let segments: Vec<_> = loops
        .iter()
        .enumerate()
        .flat_map(|(l, p)| (0..p.len()).map(move |i| (l, i, p[i], p[(i + 1) % p.len()])))
        .collect();
    for x in 0..segments.len() {
        for y in x + 1..segments.len() {
            let (l, i, a, b) = segments[x];
            let (m, j, c, d) = segments[y];
            if l == m && ((i + 1) % loops[l].len() == j || (j + 1) % loops[l].len() == i) {
                continue;
            }
            let separated = (0..2).any(|k| {
                a[k].max(b[k]) + 1e-12 < c[k].min(d[k]) || c[k].max(d[k]) + 1e-12 < a[k].min(b[k])
            });
            if separated {
                continue;
            }
            let (ab1, ab2, cd1, cd2) = (
                cross(a, b, c),
                cross(a, b, d),
                cross(c, d, a),
                cross(c, d, b),
            );
            if !(ab1 > 1e-12 && ab2 > 1e-12
                || ab1 < -1e-12 && ab2 < -1e-12
                || cd1 > 1e-12 && cd2 > 1e-12
                || cd1 < -1e-12 && cd2 < -1e-12)
            {
                return reject("spline trim loops intersect, touch, or cannot be separated");
            }
        }
    }
    for (i, p) in loops.iter().enumerate() {
        let area: f64 = (0..p.len())
            .map(|k| cross([0., 0.], p[k], p[(k + 1) % p.len()]))
            .sum();
        if area * sign * (if i == 0 { 1. } else { -1. }) <= 1e-12 {
            return reject("spline trim orientation is inconsistent with the face");
        }
        if i > 0
            && (!inside(p[0], &loops[0])
                || (1..loops.len()).any(|j| i != j && inside(p[0], &loops[j])))
        {
            return reject("spline trim hole lies outside its outer loop or is nested");
        }
    }
    Ok(())
}

fn polygon_measure(
    s: &SplineSurface,
    loops: &[Vec<[f64; 2]>],
    shift: Vec3,
    sign: f64,
    size: f64,
) -> Result<(f64, f64)> {
    let mut result = [0.; 2];
    for which in 0..2 {
        let failure = std::cell::RefCell::new(None);
        let count = std::cell::Cell::new(0usize);
        let budget = 1e-11 * size.powi(if which == 0 { 3 } else { 2 })
            / loops.iter().map(Vec::len).sum::<usize>() as f64;
        for lp in loops {
            for k in 0..lp.len() {
                let (a, b) = (lp[k], lp[(k + 1) % lp.len()]);
                let integrand = |t: f64| {
                    if failure.borrow().is_some() {
                        return 0.;
                    }
                    let (u, v) = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]));
                    let f = |x| {
                        count.set(count.get() + 1);
                        if count.get() > 2_000_000 {
                            *failure.borrow_mut() = Some(
                                "trimmed spline mass exceeds its evaluation budget".to_owned(),
                            );
                            return 0.;
                        }
                        let (q, du, dv) = s.eval(x, v);
                        let n = du.cross(dv);
                        if which == 0 {
                            q.add(shift).dot(n) / 3.
                        } else {
                            n.length() * sign
                        }
                    };
                    let mut total = 0.;
                    for span in s.u.windows(2) {
                        let hi = u.min(span[1]);
                        if hi <= span[0] {
                            continue;
                        }
                        match crate::mass::quadrature(&f, span[0], hi, budget) {
                            Ok(v) => total += v,
                            Err(e) => {
                                *failure.borrow_mut() = Some(e);
                                return 0.;
                            }
                        }
                    }
                    total * (b[1] - a[1])
                };
                let mut breaks = vec![0., 1.];
                if b[1] != a[1] {
                    for &v in &s.v {
                        let t = (v - a[1]) / (b[1] - a[1]);
                        if t > 0. && t < 1. {
                            breaks.push(t);
                        }
                    }
                }
                breaks.sort_by(f64::total_cmp);
                for p in breaks.windows(2) {
                    result[which] += crate::mass::quadrature(&integrand, p[0], p[1], budget)
                        .map_err(KernelError::Unsupported)?;
                }
            }
        }
        if let Some(e) = failure.into_inner() {
            return Err(KernelError::Unsupported(e));
        }
    }
    Ok((result[0], result[1]))
}

pub(crate) enum Domain {
    Polygon(Vec<Vec<[f64; 2]>>),
    Curved(Vec<Vec<Curve>>),
}
fn project(c: &Chart, curve: &Curve) -> Result<Curve> {
    let (p, w) = crate::spline_curve::controls(curve)?;
    Ok(Curve::RationalBezier {
        p: p.map(|p| {
            let uv = c.point(p);
            Vec3::new(uv[0], uv[1], 0.)
        }),
        w,
    })
}
fn polygonal(c: &Chart, curve: &Curve) -> bool {
    let Ok((p, _)) = crate::spline_curve::controls(curve) else {
        return false;
    };
    let uv = p.map(|p| c.point(p));
    let (a, b) = (uv[0], uv[3]);
    let d = [b[0] - a[0], b[1] - a[1]];
    let len = d[0].hypot(d[1]);
    if len < 1e-12 {
        return false;
    }
    uv.iter().all(|&p| cross(a, b, p).abs() <= 1e-12 * len)
}
/// Exact rational quadratic conic arcs, elevated in homogeneous coordinates.
/// Keeping each angular interval <= pi/2 gives strictly positive weights.
fn spans(curve: &Curve) -> Result<Vec<Curve>> {
    let (c, u, v, t0, t1) = match *curve {
        Curve::Arc { c, n, x, r, t0, t1 } => (c, x.scale(r), n.cross(x).scale(r), t0, t1),
        Curve::Conic { c, u, v, t0, t1 } => (c, u, v, t0, t1),
        Curve::Line { .. } | Curve::Bezier { .. } | Curve::RationalBezier { .. } => {
            return Ok(vec![curve.clone()])
        }
        _ => return reject("curved spline trims do not support cylinder-intersection curves"),
    };
    let n = ((t1 - t0) / std::f64::consts::FRAC_PI_2).ceil() as usize;
    if n == 0 || n > 4 {
        return reject("unusable conic angular domain for spline trimming");
    }
    let mut result = vec![];
    for i in 0..n {
        let a = t0 + (t1 - t0) * i as f64 / n as f64;
        let b = t0 + (t1 - t0) * (i + 1) as f64 / n as f64;
        let m = (a + b) / 2.;
        let weight = ((b - a) / 2.).cos();
        let q0 = u.scale(a.cos()).add(v.scale(a.sin()));
        let q1 = u.scale(m.cos()).add(v.scale(m.sin()));
        let q2 = u.scale(b.cos()).add(v.scale(b.sin()));
        let w = (1. + 2. * weight) / 3.;
        result.push(Curve::RationalBezier {
            p: [
                c.add(q0),
                c.add(q0.add(q1.scale(2.)).scale(1. / (3. * w))),
                c.add(q1.scale(2.).add(q2).scale(1. / (3. * w))),
                c.add(q2),
            ],
            w: [1., w, w, 1.],
        });
    }
    Ok(result)
}
fn split(c: &Curve) -> Result<[Curve; 2]> {
    Ok([
        crate::spline_curve::piece(c, 0., 0.5)?,
        crate::spline_curve::piece(c, 0.5, 1.)?,
    ])
}
/// Coordinate-plane roots split exact curves at every tensor knot. Derivative
/// critical points isolate all cubic roots, including tangencies and returns.
fn knot_breaks(s: &SplineSurface, uv: &Curve) -> Result<Vec<f64>> {
    let (p, w) = crate::spline_curve::controls(uv)?;
    let mut cuts = vec![0., 1.];
    for (axis, knots) in [&s.u, &s.v].iter().enumerate() {
        let q = p.map(|p| Vec3::new([p.x, p.y][axis], 0., 0.));
        let lo = q.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let hi = q.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        if hi - lo <= 1e-14 {
            continue;
        }
        let scalar = Curve::RationalBezier { p: q, w };
        for &k in &knots[1..knots.len() - 1] {
            if k < lo || k > hi {
                continue;
            }
            cuts.extend(crate::spline_curve::parameters(
                &scalar,
                Vec3::new(k, 0., 0.),
                1e-13,
            )?);
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    Ok(cuts)
}
fn support(
    s: &SplineSurface,
    c: &Chart,
    curve: &Curve,
    tol: f64,
    depth: usize,
    visits: &mut usize,
) -> Result<()> {
    *visits += 1;
    if *visits > 32768 {
        return reject("curved spline support certificate exceeds its subdivision budget");
    }
    let (p, w) = crate::spline_curve::controls(curve)?;
    crate::rational::validate(p, w)?;
    let uv = p.map(|p| c.point(p));
    let inside = uv.iter().flatten().all(|&x| x >= -1e-12 && x <= 1. + 1e-12);
    let lo: [f64; 2] =
        std::array::from_fn(|k| uv.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min));
    let hi: [f64; 2] =
        std::array::from_fn(|k| uv.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max));
    // Bound every patch touched by the positive-weight UV hull. Even numerical
    // knot-root error cannot hide a neighboring patch's geometry.
    let mut bound: f64 = 0.;
    let mut touched = 0;
    for i in 0..s.u.len() - 1 {
        if s.u[i] > hi[0] + 1e-13 || s.u[i + 1] < lo[0] - 1e-13 {
            continue;
        }
        for j in 0..s.v.len() - 1 {
            if s.v[j] > hi[1] + 1e-13 || s.v[j + 1] < lo[1] - 1e-13 {
                continue;
            }
            touched += 1;
            if touched > 4 {
                bound = f64::INFINITY;
                break;
            }
            bound = bound.max(residual(s, c, curve, i, j)?);
        }
        if touched > 4 {
            break;
        }
    }
    if touched == 0 {
        bound = f64::INFINITY;
    }
    if inside && bound <= tol {
        return Ok(());
    }
    if depth >= 24 {
        return Err(KernelError::Unsupported(format!("curved spline support remains uncertified at depth 24: bound {bound:.9} inches, limit {tol:.9}, UV hull inside domain: {inside}")));
    }
    for part in split(curve)? {
        support(s, c, &part, tol, depth + 1, visits)?;
    }
    Ok(())
}
#[derive(Clone)]
struct Leaf {
    curve: Curve,
    p: [Vec3; 4],
    error: f64,
}
fn point_segment(p: Vec3, a: Vec3, b: Vec3) -> f64 {
    let d = b.sub(a);
    let len = d.dot(d);
    if len == 0. {
        p.distance(a)
    } else {
        p.distance(a.add(d.scale((p.sub(a).dot(d) / len).clamp(0., 1.))))
    }
}
fn flatten(c: &Curve, out: &mut Vec<Leaf>, depth: usize) -> Result<()> {
    if depth > 24 || out.len() >= 4096 {
        return reject("curved spline trim subdivision exceeds its topology budget");
    }
    let (p, _) = crate::spline_curve::controls(c)?;
    let d = p[3].sub(p[0]);
    let len = d.dot(d);
    let mut previous = -1e-13;
    let mut monotone = len > 1e-24;
    for q in p {
        let t = q.sub(p[0]).dot(d) / len;
        if !t.is_finite() || t < previous - 1e-13 {
            monotone = false;
        }
        previous = t;
    }
    let error = p
        .iter()
        .map(|&q| point_segment(q, p[0], p[3]))
        .fold(0., f64::max)
        + 128. * f64::EPSILON;
    if monotone && error <= 1e-5 {
        out.push(Leaf {
            curve: c.clone(),
            p,
            error,
        });
        return Ok(());
    }
    for part in split(c)? {
        flatten(&part, out, depth + 1)?;
    }
    Ok(())
}
fn hull_separated(a: &Leaf, b: &Leaf, joint: Option<Vec3>) -> bool {
    let mut axes = vec![Vec3::new(1., 0., 0.), Vec3::new(0., 1., 0.)];
    for p in [a.p, b.p] {
        for i in 0..4 {
            for j in i + 1..4 {
                let d = p[j].sub(p[i]);
                axes.push(Vec3::new(-d.y, d.x, 0.));
            }
        }
    }
    let round = 512. * f64::EPSILON;
    for n in axes {
        let len = n.length();
        if len == 0. {
            continue;
        }
        let n = n.scale(1. / len);
        let origin = joint.unwrap_or(a.p[0]);
        let pa = a.p.map(|p| p.sub(origin).dot(n));
        let pb = b.p.map(|p| p.sub(origin).dot(n));
        for (x, y) in [(pa, pb), (pb, pa)] {
            let hi = x.into_iter().fold(f64::NEG_INFINITY, f64::max);
            let lo = y.into_iter().fold(f64::INFINITY, f64::min);
            if hi + round < lo {
                return true;
            }
            // Every rational Bernstein basis is positive in the open interval.
            // A strict control on either side keeps interiors off the joint line.
            if joint.is_some()
                && hi <= round
                && lo >= -round
                && (x.iter().any(|v| *v < -round) || y.iter().any(|v| *v > round))
            {
                return true;
            }
        }
    }
    false
}
fn separate(
    a: &Leaf,
    b: &Leaf,
    joint: Option<Vec3>,
    depth: usize,
    visits: &mut usize,
) -> Result<()> {
    *visits += 1;
    if *visits > 100000 || depth > 24 {
        return reject(
            "curved spline trim loops intersect, touch, or exceed the separation budget",
        );
    }
    if hull_separated(a, b, joint) {
        return Ok(());
    }
    let extent = |p: [Vec3; 4]| p.iter().map(|q| q.distance(p[0])).fold(0., f64::max);
    let (a, b) = if extent(a.p) >= extent(b.p) {
        (a, b)
    } else {
        (b, a)
    };
    for curve in split(&a.curve)? {
        let (p, _) = crate::spline_curve::controls(&curve)?;
        let retained = joint.filter(|j| p[0].distance(*j) <= 1e-13 || p[3].distance(*j) <= 1e-13);
        separate(
            &Leaf {
                curve,
                p,
                error: 0.,
            },
            b,
            retained,
            depth + 1,
            visits,
        )?;
    }
    Ok(())
}
fn contains_curve_loop(p: Vec3, leaves: &[Leaf]) -> Result<bool> {
    for leaf in leaves {
        if point_segment(p, leaf.p[0], leaf.p[3]) <= leaf.error + 1e-12 {
            return reject("curved trim containment is too close to certify");
        }
    }
    let polygon = leaves
        .iter()
        .map(|l| [l.p[0].x, l.p[0].y])
        .collect::<Vec<_>>();
    Ok(inside([p.x, p.y], &polygon))
}
pub(crate) fn domain(
    s: &SplineSurface,
    c: &Chart,
    body: &Body,
    face: &Face,
    tol: f64,
) -> Result<Domain> {
    if face
        .loops
        .iter()
        .flatten()
        .all(|&(e, _)| polygonal(c, &body.edges[e].curve))
    {
        return polygon_domain(s, c, body, face, tol).map(Domain::Polygon);
    }
    if face.loops.is_empty() || face.loops.iter().map(Vec::len).sum::<usize>() > 2048 {
        return reject("invalid curved spline trim loop count");
    }
    let mut curves = vec![];
    let mut flat: Vec<Vec<Leaf>> = vec![];
    let mut support_visits = 0;
    for lp in &face.loops {
        let mut loop_curves = vec![];
        let mut leaves = vec![];
        for &(e, fwd) in lp {
            let edge = if fwd {
                body.edges[e].curve.clone()
            } else {
                body.edges[e].curve.reversed()
            };
            for span in spans(&edge)? {
                let uv = project(c, &span)?;
                let cuts = knot_breaks(s, &uv)?;
                for interval in cuts.windows(2) {
                    let span = crate::spline_curve::piece(&span, interval[0], interval[1])?;
                    support(s, c, &span, tol, 0, &mut support_visits).map_err(|error| {
                        KernelError::Unsupported(format!("spline trim edge {e}: {error}"))
                    })?;
                    let uv = project(c, &span)?;
                    flatten(&uv, &mut leaves, 0)?;
                    loop_curves.push(uv);
                }
            }
        }
        if leaves.len() < 3 {
            return reject("curved trim loop is degenerate");
        }
        for i in 0..leaves.len() {
            if leaves[i].p[3].distance(leaves[(i + 1) % leaves.len()].p[0]) > 1e-12 {
                return reject("curved trim loop does not close in UV coordinates");
            }
        }
        curves.push(loop_curves);
        flat.push(leaves);
    }
    let all: Vec<_> = flat
        .iter()
        .enumerate()
        .flat_map(|(l, leaves)| leaves.iter().enumerate().map(move |(i, leaf)| (l, i, leaf)))
        .collect();
    if all.len() > 4096 {
        return reject("curved trim face exceeds 4096 topology spans");
    }
    let bounds: Vec<_> = all
        .iter()
        .map(|(_, _, l)| {
            let min = |k: usize| {
                l.p.iter()
                    .map(|p| [p.x, p.y][k])
                    .fold(f64::INFINITY, f64::min)
            };
            let max = |k: usize| {
                l.p.iter()
                    .map(|p| [p.x, p.y][k])
                    .fold(f64::NEG_INFINITY, f64::max)
            };
            ([min(0), min(1)], [max(0), max(1)])
        })
        .collect();
    let mut order: Vec<_> = (0..all.len()).collect();
    order.sort_by(|&i, &j| bounds[i].0[0].total_cmp(&bounds[j].0[0]));
    let mut visits = 0;
    let mut candidates = 0;
    for x in 0..order.len() {
        let a = order[x];
        for &b in &order[x + 1..] {
            if bounds[b].0[0] > bounds[a].1[0] + 1e-12 {
                break;
            }
            candidates += 1;
            if candidates > 2_000_000 {
                return reject("curved trim pair budget exceeded");
            }
            if bounds[a].1[1] + 1e-12 < bounds[b].0[1] || bounds[b].1[1] + 1e-12 < bounds[a].0[1] {
                continue;
            }
            let (l, i, first) = all[a];
            let (m, j, second) = all[b];
            let joint = if l == m && (i + 1) % flat[l].len() == j {
                Some(first.p[3])
            } else if l == m && (j + 1) % flat[l].len() == i {
                Some(second.p[3])
            } else {
                None
            };
            separate(first, second, joint, 0, &mut visits).map_err(|e| {
                KernelError::Unsupported(format!(
                    "trim loop {l} span {i} against loop {m} span {j}: {e}"
                ))
            })?;
        }
    }
    let sign = if face.sense { 1. } else { -1. };
    for (i, lp) in curves.iter().enumerate() {
        let area = crate::extrude::signed_area(lp, Vec3::new(0., 0., 1.));
        if !area.is_finite() || area * sign * (if i == 0 { 1. } else { -1. }) <= 1e-12 {
            return reject("curved trim winding disagrees with its face");
        }
        if i > 0 {
            let q = flat[i][0].p[0];
            if !contains_curve_loop(q, &flat[0])? {
                return reject("curved trim hole lies outside its outer loop");
            }
            for j in 1..flat.len() {
                if i != j && contains_curve_loop(q, &flat[j])? {
                    return reject("nested curved trim holes are not supported");
                }
            }
        }
    }
    Ok(Domain::Curved(curves))
}
pub(crate) fn measure(
    s: &SplineSurface,
    domain: &Domain,
    shift: Vec3,
    sign: f64,
    size: f64,
) -> Result<(f64, f64)> {
    let Domain::Curved(loops) = domain else {
        let Domain::Polygon(p) = domain else {
            unreachable!()
        };
        return polygon_measure(s, p, shift, sign, size);
    };
    let mut result = [0.; 2];
    for which in 0..2 {
        let failure = std::cell::RefCell::new(None);
        let count = std::cell::Cell::new(0usize);
        let budget = 1e-11 * size.powi(if which == 0 { 3 } else { 2 })
            / loops.iter().map(Vec::len).sum::<usize>() as f64;
        for curve in loops.iter().flatten() {
            let integrand = |t: f64| {
                if failure.borrow().is_some() {
                    return 0.;
                }
                let (uv, derivative) = curve.eval(t);
                let inner = |u| {
                    count.set(count.get() + 1);
                    if count.get() > 2_000_000 {
                        *failure.borrow_mut() =
                            Some("curved spline mass exceeds its evaluation budget".to_owned());
                        return 0.;
                    }
                    let (p, a, b) = s.eval(u, uv.y);
                    let n = a.cross(b);
                    if which == 0 {
                        p.add(shift).dot(n) / 3.
                    } else {
                        n.length() * sign
                    }
                };
                let mut total = 0.;
                for span in s.u.windows(2) {
                    let hi = uv.x.min(span[1]);
                    if hi <= span[0] {
                        continue;
                    }
                    match crate::mass::quadrature(&inner, span[0], hi, budget) {
                        Ok(v) => total += v,
                        Err(e) => {
                            *failure.borrow_mut() = Some(e);
                            return 0.;
                        }
                    }
                }
                total * derivative.y
            };
            result[which] += crate::mass::quadrature(&integrand, 0., 1., budget)
                .map_err(KernelError::Unsupported)?;
        }
        if let Some(e) = failure.into_inner() {
            return Err(KernelError::Unsupported(e));
        }
    }
    Ok((result[0], result[1]))
}

#[cfg(test)]
mod curved_tests {
    use super::*;
    #[test]
    fn conic_spans_retain_the_independent_ellipse_equation() {
        let curve = Curve::Conic {
            c: Vec3::new(3., -2., 1.),
            u: Vec3::new(2., 0., 0.),
            v: Vec3::new(0.5, 1., 0.),
            t0: -0.7,
            t1: 4.8,
        };
        let pieces = spans(&curve).unwrap();
        assert_eq!(pieces.len(), 4);
        assert!(pieces[0].start().distance(curve.start()) < 1e-14);
        assert!(pieces.last().unwrap().end().distance(curve.end()) < 1e-14);
        for p in pieces {
            for i in 0..=20 {
                let q = p.eval(i as f64 / 20.).0;
                let v = q.y + 2.;
                let u = (q.x - 3. - 0.5 * v) / 2.;
                assert!((u * u + v * v - 1.).abs() < 1e-13);
                assert!((q.z - 1.).abs() < 1e-14);
            }
        }
    }
    #[test]
    fn perpendicular_joints_are_separate_but_shared_overlaps_are_not() {
        let line = |a, b| {
            let curve = Curve::Line { a, b };
            let (p, _) = crate::spline_curve::controls(&curve).unwrap();
            Leaf {
                curve,
                p,
                error: 0.,
            }
        };
        let j = Vec3::new(0.5, 0.5, 0.);
        let a = line(Vec3::new(0.9, 0.5, 0.), j);
        let b = line(j, Vec3::new(0.5, 0.1, 0.));
        assert!(hull_separated(&a, &b, Some(j)));
        assert!(!hull_separated(&a, &b, None));
        let overlap = line(j, Vec3::new(0.8, 0.5, 0.));
        assert!(!hull_separated(&a, &overlap, Some(j)));
    }
}

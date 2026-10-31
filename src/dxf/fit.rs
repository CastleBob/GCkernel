//! Ellipses and splines as tangent-continuous lines and arcs within a tolerance (biarcs: two arcs meeting tangent,
//! matching the curve's end points and end tangents, equal tangent lengths; Bolton 1975, Meek and Walton 1992), plus
//! the B-spline evaluation and fit-point interpolation they need (Piegl and Tiller, The NURBS Book, A2.2/A2.3, 9.2.1).

use std::f64::consts::{PI, TAU};

use super::{Geometry, P2};
use crate::KernelError;

pub fn dot(a: P2, b: P2) -> f64 { a[0] * b[0] + a[1] * b[1] }
pub fn cross(a: P2, b: P2) -> f64 { a[0] * b[1] - a[1] * b[0] }
pub fn len(a: P2) -> f64 { dot(a, a).sqrt() }
pub fn dist(a: P2, b: P2) -> f64 { len(sub(b, a)) }
fn sub(a: P2, b: P2) -> P2 { [a[0] - b[0], a[1] - b[1]] }
fn add(a: P2, b: P2) -> P2 { [a[0] + b[0], a[1] + b[1]] }
fn mul(a: P2, s: f64) -> P2 { [a[0] * s, a[1] * s] }
fn unit(a: P2) -> Option<P2> {
    let l = len(a);
    (l > 0.0 && l.is_finite()).then(|| mul(a, 1.0 / l))
}

/// The knot span holding u for `count` control points of degree `p`.
fn find_span(count: usize, p: usize, u: f64, knots: &[f64]) -> usize {
    let n = count - 1;
    if u >= knots[n + 1] {
        return n;
    }
    if u <= knots[p] {
        return p;
    }
    let (mut low, mut high) = (p, n + 1);
    let mut mid = (low + high) / 2;
    while u < knots[mid] || u >= knots[mid + 1] {
        if u < knots[mid] { high = mid; } else { low = mid; }
        mid = (low + high) / 2;
    }
    mid
}

/// The p + 1 basis functions that are non-zero in `span` at u, and their first derivatives.
fn basis(span: usize, u: f64, p: usize, knots: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let values: Vec<f64> = (0..=p).map(|r| ndu[r][p]).collect();
    let mut ders = vec![0.0; p + 1];
    if p > 0 {
        for (r, d) in ders.iter_mut().enumerate() {
            let mut v = 0.0;
            if r >= 1 { v += ndu[r - 1][p - 1] / ndu[p][r - 1]; }
            if r < p { v -= ndu[r][p - 1] / ndu[p][r]; }
            *d = v * p as f64;
        }
    }
    (values, ders)
}

/// Point and first derivative of a (rational) B-spline at u, using the polynomial piece of `span`.
pub fn spline_at(degree: usize, knots: &[f64], points: &[P2], weights: &[f64], span: usize, u: f64) -> (P2, P2) {
    let (n, dn) = basis(span, u, degree, knots);
    let (mut a, mut da, mut w, mut dw) = ([0.0; 2], [0.0; 2], 0.0, 0.0);
    for j in 0..=degree {
        let i = span - degree + j;
        let wi = if weights.is_empty() { 1.0 } else { weights[i] };
        a = add(a, mul(points[i], n[j] * wi));
        da = add(da, mul(points[i], dn[j] * wi));
        w += n[j] * wi;
        dw += dn[j] * wi;
    }
    let c = mul(a, 1.0 / w);
    (c, mul(sub(da, mul(c, dw)), 1.0 / w))
}

/// The spline of degree <= `degree` through `through` (chord-length parameters, averaged knots), leaving along and
/// arriving along the `ends` tangents when given; None when the points repeat or are too many to solve.
pub fn interpolate(through: &[P2], degree: usize, ends: Option<(P2, P2)>) -> Option<Geometry> {
    let mut q: Vec<P2> = Vec::with_capacity(through.len());
    for p in through {
        if q.last().is_none_or(|l| dist(*l, *p) > 0.0) {
            q.push(*p);
        }
    }
    if q.len() < 2 || q.len() > 2000 {
        return None;
    }
    let n = q.len() - 1;
    let p = degree.clamp(1, n);
    let total: f64 = q.windows(2).map(|w| dist(w[0], w[1])).sum();
    let mut u = vec![0.0; n + 1];
    for k in 1..=n {
        u[k] = u[k - 1] + dist(q[k - 1], q[k]) / total;
    }
    u[n] = 1.0;
    let ends = ends.filter(|_| p >= 2).and_then(|(a, b)| Some((unit(a)?, unit(b)?)));
    // With end tangents two more control points carry them (The NURBS Book 9.2.2).
    let size = if ends.is_some() { n + 3 } else { n + 1 };
    let mut knots = vec![0.0; size + p + 1];
    for k in knots.iter_mut().skip(size) {
        *k = 1.0;
    }
    if ends.is_some() {
        for j in 1..=(n + 2 - p) {
            knots[j + p] = u[j - 1..j - 1 + p].iter().sum::<f64>() / p as f64;
        }
    } else {
        for j in 1..=(n - p) {
            knots[j + p] = u[j..j + p].iter().sum::<f64>() / p as f64;
        }
    }
    let mut a: Vec<Vec<f64>> = Vec::with_capacity(size);
    let mut b: Vec<P2> = Vec::with_capacity(size);
    let row_at = |t: f64, target: P2, a: &mut Vec<Vec<f64>>, b: &mut Vec<P2>| {
        let span = find_span(size, p, t, &knots);
        let (values, _) = basis(span, t, p, &knots);
        let mut row = vec![0.0; size];
        for j in 0..=p {
            row[span - p + j] = values[j];
        }
        a.push(row);
        b.push(target);
    };
    for k in 0..=n {
        row_at(u[k], q[k], &mut a, &mut b);
    }
    if let Some((t0, t1)) = ends {
        // In order beside the rows they neighbour, so the system stays banded.
        let mut start = vec![0.0; size];
        start[0] = -1.0;
        start[1] = 1.0;
        a.insert(1, start);
        b.insert(1, mul(t0, total * knots[p + 1] / p as f64));
        let mut end = vec![0.0; size];
        end[size - 2] = -1.0;
        end[size - 1] = 1.0;
        a.insert(size - 1, end);
        b.insert(size - 1, mul(t1, total * (1.0 - knots[size - 1]) / p as f64));
    }
    // Banded (half-width p + 1), solved with partial pivoting inside the band.
    let band = 2 * p + 3;
    for col in 0..size {
        let reach = (col + band).min(size);
        let pivot = (col..reach).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() < 1e-14 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..reach {
            let f = a[row][col] / a[col][col];
            if f != 0.0 {
                for k in col..(col + 2 * band).min(size) { a[row][k] -= f * a[col][k]; }
                b[row] = sub(b[row], mul(b[col], f));
            }
        }
    }
    let mut points = vec![[0.0; 2]; size];
    for row in (0..size).rev() {
        let mut s = b[row];
        for k in row + 1..(row + 2 * band).min(size) { s = sub(s, mul(points[k], a[row][k])); }
        points[row] = mul(s, 1.0 / a[row][row]);
    }
    Some(Geometry::Spline { degree: p, knots, points, weights: Vec::new() })
}

/// Lines, arcs, circles and polylines unchanged; ellipses and splines as tangent lines and arcs that stay within
/// `tol` of the curve.
pub fn lines_and_arcs(g: &Geometry, tol: f64) -> Result<Vec<Geometry>, KernelError> {
    if !(tol > 0.0 && tol.is_finite()) {
        return Err(KernelError::InvalidInput("the curve fitting tolerance must be positive".into()));
    }
    let mut out = Vec::new();
    let mut budget = MAX_SPANS;
    match g {
        Geometry::Line { .. } | Geometry::Arc { .. } | Geometry::Circle { .. } | Geometry::Polyline { .. } => {
            out.push(g.clone())
        }
        Geometry::Ellipse { c, u, v, t0, t1 } => {
            let size = len(*u).max(len(*v));
            let circular = (len(*u) - len(*v)).abs() <= 1e-9 * size && dot(*u, *v).abs() <= 1e-9 * size * size;
            if circular && t1 - t0 >= TAU - 1e-12 {
                out.push(Geometry::Circle { c: *c, r: len(*u) });
            } else if circular {
                let p0 = add(mul(*u, t0.cos()), mul(*v, t0.sin()));
                out.push(Geometry::Arc { c: *c, r: len(*u), start: p0[1].atan2(p0[0]),
                    sweep: (t1 - t0) * cross(*u, *v).signum() });
            } else {
                let f = |t: f64| (add(*c, add(mul(*u, t.cos()), mul(*v, t.sin()))), add(mul(*u, -t.sin()), mul(*v, t.cos())));
                let pieces = ((t1 - t0) / (PI / 2.0)).ceil().max(1.0) as usize;
                for k in 0..pieces {
                    let a = t0 + (t1 - t0) * k as f64 / pieces as f64;
                    let b = t0 + (t1 - t0) * (k + 1) as f64 / pieces as f64;
                    fit_span(&f, a, b, tol, &mut out, 0, &mut budget)?;
                }
            }
        }
        Geometry::Spline { degree, knots, points, weights } => {
            if *degree == 1 {
                for w in points.windows(2) {
                    if dist(w[0], w[1]) > 0.0 {
                        out.push(Geometry::Line { a: w[0], b: w[1] });
                    }
                }
            } else {
                for span in *degree..points.len() {
                    let (a, b) = (knots[span], knots[span + 1]);
                    if b > a {
                        let f = |t: f64| spline_at(*degree, knots, points, weights, span, t);
                        fit_span(&f, a, b, tol, &mut out, 0, &mut budget)?;
                    }
                }
            }
        }
    }
    Ok(out)
}

const SAMPLES: usize = 32;
const MAX_SPLITS: usize = 48;
/// Spans tried per curve before giving up (a damaged curve must not take forever).
const MAX_SPANS: usize = 50_000;

fn fit_span(f: &dyn Fn(f64) -> (P2, P2), a: f64, b: f64, tol: f64, out: &mut Vec<Geometry>, depth: usize,
    budget: &mut usize) -> Result<(), KernelError>
{
    if *budget == 0 {
        return Err(KernelError::Unsupported("a curve needs too many lines and arcs to fit".into()));
    }
    *budget -= 1;
    let samples: Vec<P2> = (0..=SAMPLES).map(|i| f(a + (b - a) * i as f64 / SAMPLES as f64).0).collect();
    let (p0, p1) = (samples[0], samples[SAMPLES]);
    let chord = dist(p0, p1);
    // Only a truly straight piece becomes a line; anything that bends gets a biarc, so every joint stays tangent.
    if chord > 0.0 && samples.iter().all(|s| to_segment(*s, p0, p1) <= 1e-9 * chord.max(tol)) {
        out.push(Geometry::Line { a: p0, b: p1 });
        return Ok(());
    }
    if chord == 0.0 && samples.iter().all(|s| dist(*s, p0) <= tol) {
        return Ok(()); // a piece no longer than the tolerance that starts and ends at one point
    }
    let t0 = tangent(f, a, b, true);
    let t1 = tangent(f, a, b, false);
    if let (Some(t0), Some(t1)) = (t0, t1) {
        if chord > 0.0 {
            if let Some(pair) = biarc(p0, t0, p1, t1) {
                // Checked at sample points, so aim inside the tolerance: the true deviation peaks between them.
                let aim = 0.8 * tol;
                let near = samples.iter().all(|s| pair.iter().map(|g| to_curve(*s, g)).fold(f64::MAX, f64::min) <= aim);
                let back = pair.iter().all(|g| along(g).iter().all(|q| to_polyline(*q, &samples) <= aim));
                if near && back {
                    out.extend(pair);
                    return Ok(());
                }
            }
        }
    }
    if depth >= MAX_SPLITS {
        return Err(KernelError::Unsupported(format!(
            "a curve could not be fitted with lines and arcs within {tol} (it turns too sharply)")));
    }
    let m = 0.5 * (a + b);
    fit_span(f, a, m, tol, out, depth + 1, budget)?;
    fit_span(f, m, b, tol, out, depth + 1, budget)
}

/// Unit tangent at the start (or end) of [a, b], from a nearby chord where the derivative vanishes (a cusp).
fn tangent(f: &dyn Fn(f64) -> (P2, P2), a: f64, b: f64, start: bool) -> Option<P2> {
    let (p, d) = f(if start { a } else { b });
    let h = (b - a) * 1e-4;
    let q = f(if start { a + h } else { b - h }).0;
    let chord = if start { sub(q, p) } else { sub(p, q) };
    match unit(d) {
        Some(t) if len(d) * h > 1e-9 * len(chord) => Some(t),
        _ => unit(chord),
    }
}

fn biarc(p0: P2, t0: P2, p1: P2, t1: P2) -> Option<[Geometry; 2]> {
    let v = sub(p1, p0);
    let vv = dot(v, v);
    let t = add(t0, t1);
    let denom = 2.0 * (1.0 - dot(t0, t1));
    let d = if denom.abs() < 1e-12 {
        let vt = dot(v, t1);
        if vt <= 0.0 { return None; }
        vv / (4.0 * vt)
    } else {
        let vt = dot(v, t);
        (-vt + (vt * vt + denom * vv).sqrt()) / denom
    };
    if !(d > 0.0 && d.is_finite()) {
        return None;
    }
    let q0 = add(p0, mul(t0, d));
    let q1 = sub(p1, mul(t1, d));
    let j = mul(add(q0, q1), 0.5);
    let tj = unit(sub(q1, q0)).unwrap_or(t0);
    Some([arc_from(p0, t0, j)?, arc_from(j, tj, p1)?])
}

/// The arc (or line) leaving p along unit tangent t and ending at q.
fn arc_from(p: P2, t: P2, q: P2) -> Option<Geometry> {
    let w = sub(q, p);
    let ww = dot(w, w);
    if ww == 0.0 {
        return None;
    }
    let turn = cross(t, w);
    if turn.abs() <= 1e-12 * ww {
        return (dot(t, w) > 0.0).then_some(Geometry::Line { a: p, b: q });
    }
    let signed = ww / (2.0 * turn);
    let c = add(p, mul([-t[1], t[0]], signed));
    let (a, b) = (sub(p, c), sub(q, c));
    let mut sweep = cross(a, b).atan2(dot(a, b));
    if turn > 0.0 && sweep < 0.0 { sweep += TAU; }
    if turn < 0.0 && sweep > 0.0 { sweep -= TAU; }
    Some(Geometry::Arc { c, r: signed.abs(), start: a[1].atan2(a[0]), sweep })
}

fn to_segment(p: P2, a: P2, b: P2) -> f64 {
    let ab = sub(b, a);
    let l = dot(ab, ab);
    if l == 0.0 {
        return dist(p, a);
    }
    let s = (dot(sub(p, a), ab) / l).clamp(0.0, 1.0);
    dist(p, add(a, mul(ab, s)))
}

fn to_polyline(p: P2, points: &[P2]) -> f64 {
    points.windows(2).map(|w| to_segment(p, w[0], w[1])).fold(f64::MAX, f64::min)
}

/// Distance from p to a line or arc piece.
fn to_curve(p: P2, g: &Geometry) -> f64 {
    match g {
        Geometry::Line { a, b } => to_segment(p, *a, *b),
        Geometry::Arc { c, r, start, sweep } => {
            let d = sub(p, *c);
            let mut phi = (d[1].atan2(d[0]) - start) * sweep.signum();
            phi = phi.rem_euclid(TAU);
            if phi <= sweep.abs() {
                (len(d) - r).abs()
            } else {
                let e = add(*c, mul([(start + sweep).cos(), (start + sweep).sin()], *r));
                let s = add(*c, mul([start.cos(), start.sin()], *r));
                dist(p, s).min(dist(p, e))
            }
        }
        _ => f64::MAX,
    }
}

/// A few points along a line or arc piece.
fn along(g: &Geometry) -> Vec<P2> {
    (1..8).map(|i| {
        let s = i as f64 / 8.0;
        match g {
            Geometry::Line { a, b } => add(*a, mul(sub(*b, *a), s)),
            Geometry::Arc { c, r, start, sweep } => {
                let t = start + sweep * s;
                add(*c, mul([t.cos(), t.sin()], *r))
            }
            _ => [0.0; 2],
        }
    }).collect()
}

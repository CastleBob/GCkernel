//! A painted subpath's pieces as drawing geometry. PDF has no circles: writers draw them as cubic Beziers with the
//! standard construction (a cubic from angle a to b on a circle of radius r has its handles along the tangents,
//! (4/3) tan((b - a)/4) r long; 0.5523 r for a quarter). Runs of such cubics on one circle become exact `Arc`s and
//! `Circle`s, four quarter cubics on conjugate semi-diameters an exact `Ellipse`, and cubics whose handles lie on
//! the chord become `Line`s. Everything else stays a cubic B-spline.
//!
//! Tolerances: 1e-6 of the size, plus the rounding the file's numbers carry (half the last written digit, mapped to
//! inches). Handle directions must match within that; handle lengths within a further 0.5 %, so the handle
//! constants writers use for a quarter (0.55228475, 0.5523, 0.551915, 0.55) all count. A 0.5 % longer or shorter
//! handle moves the cubic at most about 0.15 % of the radius off the true circle, which is what is output. A run of
//! arcs turning a full circle that misses its start by at most 0.1 % of the radius (single-precision drift in the
//! writer) is a circle.

use std::f64::consts::{PI, TAU};

use super::super::{fit, Geometry, P2};

#[derive(Clone, Copy, Debug)]
pub enum Seg {
    /// Ends and the rounding of the numbers that gave them (inches).
    Line(P2, P2, f64),
    Cubic([P2; 4], f64),
}

#[derive(Clone, Copy, Debug)]
struct ArcPiece {
    c: P2,
    r: f64,
    sweep: f64,
    a: P2,
    b: P2,
    q: f64,
}

enum Piece {
    Line(P2, P2),
    Arc(ArcPiece),
    Cubic([P2; 4]),
}

fn add(a: P2, b: P2) -> P2 { [a[0] + b[0], a[1] + b[1]] }
fn sub(a: P2, b: P2) -> P2 { [a[0] - b[0], a[1] - b[1]] }
fn mul(a: P2, s: f64) -> P2 { [a[0] * s, a[1] * s] }
fn rotate(a: P2, t: f64) -> P2 { [a[0] * t.cos() - a[1] * t.sin(), a[0] * t.sin() + a[1] * t.cos()] }

fn tiny(a: P2, b: P2) -> bool {
    fit::dist(a, b) <= 1e-12 * (1.0 + fit::len(a).max(fit::len(b)))
}

/// Handle length for a circular arc of sweep `theta` and radius r.
fn handle(theta: f64, r: f64) -> f64 {
    4.0 / 3.0 * (theta.abs() / 4.0).tan() * r
}

fn straight(p: &[P2; 4], q: f64) -> bool {
    let chord = sub(p[3], p[0]);
    let l2 = fit::dot(chord, chord);
    if l2 <= 0.0 {
        return false;
    }
    let l = l2.sqrt();
    [p[1], p[2]].iter().all(|&c| {
        let d = sub(c, p[0]);
        let t = fit::dot(d, chord) / l2;
        fit::cross(chord, d).abs() / l <= 1e-9 * l + q && (-1e-9..=1.0 + 1e-9).contains(&t)
    })
}

/// Whether a control point is the ideal handle `from + ideal` with its direction within `base` and its length within
/// `base` plus 0.5 %.
fn handle_ok(actual: P2, from: P2, ideal: P2, base: f64) -> bool {
    let h = fit::len(ideal);
    if h <= 0.0 {
        return false;
    }
    let e = sub(actual, add(from, ideal));
    let along = fit::dot(e, ideal) / h;
    let across = fit::cross(ideal, e) / h;
    along.abs() <= base + 5e-3 * h && across.abs() <= base
}

fn arc_of(p: &[P2; 4], q: f64) -> Option<ArcPiece> {
    let chord = sub(p[3], p[0]);
    let l = fit::len(chord);
    let d0 = sub(p[1], p[0]);
    let d3 = sub(p[3], p[2]);
    if l <= 0.0 || fit::len(d0) <= 0.0 || fit::len(d3) <= 0.0 {
        return None;
    }
    let a0 = fit::cross(chord, d0).atan2(fit::dot(chord, d0));
    let a3 = fit::cross(chord, d3).atan2(fit::dot(chord, d3));
    let half = (a3 - a0) / 2.0;
    if half.abs() < 1e-7 || half.abs() > PI / 2.0 + 1e-3 {
        return None;
    }
    let r = l / (2.0 * half.abs().sin());
    let unit = mul(chord, 1.0 / l);
    let n = [-unit[1], unit[0]];
    let c = add(mul(add(p[0], p[3]), 0.5), mul(n, r * half.cos() * half.signum()));
    let h = handle(2.0 * half, r);
    let p1 = add(p[0], mul(rotate(unit, -half), h));
    let p2 = sub(p[3], mul(rotate(unit, half), h));
    let base = 1e-6 * r + 4.0 * q;
    (handle_ok(p[1], p[0], sub(p1, p[0]), base) && handle_ok(p[2], p[3], sub(p2, p[3]), base))
        .then_some(ArcPiece { c, r, sweep: 2.0 * half, a: p[0], b: p[3], q })
}

/// The least-squares circle through points (exact for points on a circle), or None when they are collinear.
fn circle_fit(pts: &[P2]) -> Option<P2> {
    let n = pts.len() as f64;
    let m = pts.iter().fold([0.0, 0.0], |s, p| add(s, *p));
    let m = mul(m, 1.0 / n);
    let (mut suu, mut svv, mut suv, mut suuu, mut svvv, mut suvv, mut svuu) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for p in pts {
        let (u, v) = (p[0] - m[0], p[1] - m[1]);
        suu += u * u;
        svv += v * v;
        suv += u * v;
        suuu += u * u * u;
        svvv += v * v * v;
        suvv += u * v * v;
        svuu += v * u * u;
    }
    let det = suu * svv - suv * suv;
    if det.abs() <= 1e-12 * (suu * svv).max(1e-300) {
        return None;
    }
    let (b1, b2) = ((suuu + suvv) / 2.0, (svvv + svuu) / 2.0);
    let uc = (b1 * svv - b2 * suv) / det;
    let vc = (suu * b2 - suv * b1) / det;
    Some([m[0] + uc, m[1] + vc])
}

fn same_circle(x: &ArcPiece, y: &ArcPiece) -> bool {
    let tol = 1e-6 * x.r.max(y.r) + 16.0 * x.q.max(y.q);
    x.sweep.signum() == y.sweep.signum() && fit::dist(x.c, y.c) <= tol && (x.r - y.r).abs() <= tol
}

/// Four closed quarter cubics c + w_i -> c + w_i+1 (w = u, v, -u, -v) with handles k w_i+1 and k w_i.
fn ellipse(cubics: &[[P2; 4]], q: f64) -> Option<Geometry> {
    let e: Vec<P2> = cubics.iter().map(|p| p[0]).collect();
    let c = mul(e.iter().fold([0.0, 0.0], |s, p| add(s, *p)), 0.25);
    let (u, v) = (sub(e[0], c), sub(e[1], c));
    let size = fit::len(u).max(fit::len(v));
    if size <= 0.0 || fit::cross(u, v).abs() <= 1e-9 * size * size {
        return None;
    }
    let k = handle(PI / 2.0, 1.0);
    let tol = 1e-6 * size + 4.0 * q;
    if fit::dist(e[2], sub(c, u)) > tol || fit::dist(e[3], sub(c, v)) > tol || fit::dist(cubics[3][3], e[0]) > tol.max(1e-3 * size) {
        return None;
    }
    let w = [u, v, mul(u, -1.0), mul(v, -1.0)];
    for (i, p) in cubics.iter().enumerate() {
        let (wi, wj) = (w[i], w[(i + 1) % 4]);
        let (a, b) = (add(c, wi), add(c, wj));
        let end_tol = if i == 3 { tol.max(1e-3 * size) } else { tol };
        if !handle_ok(p[1], a, mul(wj, k), tol) || !handle_ok(p[2], b, mul(wi, k), tol) || fit::dist(p[3], b) > end_tol {
            return None;
        }
    }
    if (fit::len(u) - fit::len(v)).abs() <= tol && fit::dot(u, v).abs() <= tol * size {
        let r = e.iter().map(|p| fit::dist(*p, c)).sum::<f64>() / 4.0;
        return Some(Geometry::Circle { c, r });
    }
    Some(Geometry::Ellipse { c, u, v, t0: 0.0, t1: TAU })
}

fn run_geometry(run: &[ArcPiece]) -> Geometry {
    let total: f64 = run.iter().map(|a| a.sweep).sum();
    let first = run[0].a;
    let last = run[run.len() - 1].b;
    let q = run.iter().map(|a| a.q).fold(0.0, f64::max);
    let mut pts: Vec<P2> = std::iter::once(first).chain(run.iter().map(|a| a.b)).collect();
    // A full turn that misses its start by at most 0.1 % of the radius (single-precision drift) is a circle.
    let closed = fit::dist(first, last) <= (1e-6 * run[0].r + 4.0 * q).max(1e-3 * run[0].r) && total.abs() >= TAU - 2e-3;
    if closed {
        pts.pop();
    }
    let c = if pts.len() >= 3 { circle_fit(&pts) } else { None }.unwrap_or_else(|| {
        mul(run.iter().fold([0.0, 0.0], |s, a| add(s, a.c)), 1.0 / run.len() as f64)
    });
    let r = pts.iter().map(|p| fit::dist(*p, c)).sum::<f64>() / pts.len() as f64;
    if closed {
        return Geometry::Circle { c, r };
    }
    let start = (first[1] - c[1]).atan2(first[0] - c[0]);
    let end = (last[1] - c[1]).atan2(last[0] - c[0]);
    let mut sweep = end - start;
    if total > 0.0 {
        while sweep <= 0.0 { sweep += TAU; }
        while sweep >= TAU { sweep -= TAU; }
    } else {
        while sweep >= 0.0 { sweep -= TAU; }
        while sweep <= -TAU { sweep += TAU; }
    }
    // A sweep near a whole turn that does not close: keep the pieces' own total.
    if (sweep - total).abs() > 0.5 {
        sweep = total.clamp(-TAU + 1e-9, TAU - 1e-9);
    }
    Geometry::Arc { c, r, start, sweep }
}

fn spline(p: &[P2; 4]) -> Geometry {
    Geometry::Spline { degree: 3, knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], points: p.to_vec(), weights: Vec::new() }
}

/// The geometry of one subpath.
pub fn convert(segs: &[Seg]) -> Vec<Geometry> {
    let mut pieces: Vec<Piece> = Vec::with_capacity(segs.len());
    let mut cubics: Vec<[P2; 4]> = Vec::new();
    let mut q_max: f64 = 0.0;
    let mut lines = 0;
    for s in segs {
        match *s {
            Seg::Line(a, b, q) => {
                if !tiny(a, b) {
                    pieces.push(Piece::Line(a, b));
                    // A closing line no longer than the numbers' rounding does not stop a circle or ellipse.
                    if fit::dist(a, b) > 4.0 * q {
                        lines += 1;
                    }
                }
                q_max = q_max.max(q);
            }
            Seg::Cubic(p, q) => {
                q_max = q_max.max(q);
                if p.iter().all(|x| tiny(*x, p[0])) {
                    continue;
                }
                cubics.push(p);
                if straight(&p, q) {
                    pieces.push(Piece::Line(p[0], p[3]));
                } else if let Some(a) = arc_of(&p, q) {
                    pieces.push(Piece::Arc(a));
                } else {
                    pieces.push(Piece::Cubic(p));
                }
            }
        }
    }
    if lines == 0 && cubics.len() == 4 && pieces.len() - cubics.len() <= 1 {
        if let Some(g) = ellipse(&cubics, q_max) {
            let one_circle = pieces.iter().all(|p| matches!(p, Piece::Arc(_)));
            if !one_circle || matches!(g, Geometry::Circle { .. }) {
                return vec![g];
            }
        }
    }
    // A closed path whose first and last pieces are arcs of one circle: start it at a break so the arc stays whole.
    let is_closed = match (segs.first(), segs.last()) {
        (Some(f), Some(l)) => {
            let start = match f { Seg::Line(a, _, _) => *a, Seg::Cubic(p, _) => p[0] };
            let end = match l { Seg::Line(_, b, _) => *b, Seg::Cubic(p, _) => p[3] };
            fit::dist(start, end) <= 1e-9 * (1.0 + fit::len(start)) + 4.0 * q_max
        }
        _ => false,
    };
    if is_closed && pieces.len() > 1 {
        let joined = |x: &Piece, y: &Piece| matches!((x, y), (Piece::Arc(a), Piece::Arc(b)) if same_circle(a, b));
        if pieces.windows(2).any(|w| !joined(&w[0], &w[1])) {
            while joined(&pieces[pieces.len() - 1], &pieces[0]) {
                let last = pieces.pop().unwrap_or(Piece::Line([0.0; 2], [0.0; 2]));
                pieces.insert(0, last);
            }
        }
    }
    let mut out = Vec::with_capacity(pieces.len());
    let mut run: Vec<ArcPiece> = Vec::new();
    for p in pieces {
        match p {
            Piece::Arc(a) => {
                if run.last().is_some_and(|l| !same_circle(l, &a)) {
                    out.push(run_geometry(&run));
                    run.clear();
                }
                run.push(a);
            }
            other => {
                if !run.is_empty() {
                    out.push(run_geometry(&run));
                    run.clear();
                }
                match other {
                    Piece::Line(a, b) => out.push(Geometry::Line { a, b }),
                    Piece::Cubic(c) => out.push(spline(&c)),
                    Piece::Arc(_) => {}
                }
            }
        }
    }
    if !run.is_empty() {
        out.push(run_geometry(&run));
    }
    out
}

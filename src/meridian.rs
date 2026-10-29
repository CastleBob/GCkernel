//! Meridian curves of general surfaces of revolution (Claude 2026-10-29, Bob: the whole revolve subsystem): an ellipse
//! arc or a cubic Bezier curve in (radius, height). Exact evaluation and closed-form integrals (trigonometric
//! polynomials for the ellipse, polynomials for the Bezier curve) for volume and the divergence-theorem terms; the area
//! density r |M'| has no closed form and is integrated by adaptive Gauss-Legendre quadrature.

use crate::mass::quadrature;

/// A meridian in (radius, height).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Meridian {
    /// c + u cos t + v sin t (t any real; a piece uses an interval).
    Ellipse { c: (f64, f64), u: (f64, f64), v: (f64, f64) },
    /// Cubic Bezier with control points b[0..4], t in [0, 1].
    Bezier([(f64, f64); 4]),
}

/// A polynomial in t (power basis, lowest first).
#[derive(Clone, Debug)]
struct Poly(Vec<f64>);

impl Poly {
    fn mul(&self, o: &Poly) -> Poly {
        let mut out = vec![0.0; self.0.len() + o.0.len() - 1];
        for (i, a) in self.0.iter().enumerate() { for (j, b) in o.0.iter().enumerate() { out[i + j] += a * b; } }
        Poly(out)
    }
    fn add(&self, o: &Poly) -> Poly {
        let n = self.0.len().max(o.0.len());
        Poly((0..n).map(|i| self.0.get(i).unwrap_or(&0.0) + o.0.get(i).unwrap_or(&0.0)).collect())
    }
    fn scale(&self, k: f64) -> Poly { Poly(self.0.iter().map(|a| a * k).collect()) }
    fn deriv(&self) -> Poly { Poly(if self.0.len() <= 1 { vec![0.0] } else { (1..self.0.len()).map(|i| self.0[i] * i as f64).collect() }) }
    /// The antiderivative vanishing at 0, evaluated at t.
    fn integral_at(&self, t: f64) -> f64 { self.0.iter().enumerate().rev().fold(0.0, |acc, (i, a)| acc * t + a / (i + 1) as f64) * t }
}

/// A trigonometric polynomial: a0 + sum (a_k cos kt + b_k sin kt), stored as [a0, a1, b1, a2, b2, ...].
#[derive(Clone, Debug)]
struct Trig(Vec<f64>);

impl Trig {
    fn degree(&self) -> usize { (self.0.len() - 1) / 2 }
    fn coef(&self, k: usize) -> (f64, f64) {
        if k == 0 { (self.0[0], 0.0) } else if k <= self.degree() { (self.0[2 * k - 1], self.0[2 * k]) } else { (0.0, 0.0) }
    }
    fn from(deg: usize, f: impl Fn(usize) -> (f64, f64)) -> Trig {
        let mut v = vec![f(0).0];
        for k in 1..=deg { let (a, b) = f(k); v.push(a); v.push(b); }
        Trig(v)
    }
    fn mul(&self, o: &Trig) -> Trig {
        // cos a cos b = (cos(a-b) + cos(a+b)) / 2, sin a sin b = (cos(a-b) - cos(a+b)) / 2, sin a cos b = (sin(a+b) + sin(a-b)) / 2.
        let deg = self.degree() + o.degree();
        let mut acc = vec![(0.0, 0.0); deg + 1];
        for i in 0..=self.degree() {
            let (ai, bi) = self.coef(i);
            for j in 0..=o.degree() {
                let (aj, bj) = o.coef(j);
                let (s, d) = (i + j, i.abs_diff(j));
                let sign = if i >= j { 1.0 } else { -1.0 };
                let (ci, si) = if i == 0 { (ai, 0.0) } else { (ai, bi) };
                let (cj, sj) = if j == 0 { (aj, 0.0) } else { (aj, bj) };
                // cos*cos
                acc[s].0 += ci * cj / 2.0; acc[d].0 += ci * cj / 2.0;
                // sin*sin
                acc[d].0 += si * sj / 2.0; acc[s].0 -= si * sj / 2.0;
                // sin_i cos_j: (sin(i+j) + sin(i-j)) / 2
                acc[s].1 += si * cj / 2.0; acc[d].1 += sign * si * cj / 2.0;
                // cos_i sin_j: (sin(i+j) - sin(i-j)) / 2
                acc[s].1 += ci * sj / 2.0; acc[d].1 -= sign * ci * sj / 2.0;
            }
        }
        acc[0].1 = 0.0;
        Trig::from(deg, |k| acc[k])
    }
    fn add(&self, o: &Trig) -> Trig {
        let deg = self.degree().max(o.degree());
        Trig::from(deg, |k| { let (a, b) = self.coef(k); let (c, d) = o.coef(k); (a + c, b + d) })
    }
    fn scale(&self, s: f64) -> Trig { Trig(self.0.iter().map(|a| a * s).collect()) }
    fn deriv(&self) -> Trig { Trig::from(self.degree(), |k| { let (a, b) = self.coef(k); if k == 0 { (0.0, 0.0) } else { (b * k as f64, -a * k as f64) } }) }
    /// The antiderivative a0 t + sum (a_k sin kt - b_k cos kt) / k, evaluated at t.
    fn integral_at(&self, t: f64) -> f64 {
        let mut out = self.0[0] * t;
        for k in 1..=self.degree() {
            let (a, b) = self.coef(k);
            let kf = k as f64;
            out += (a * (kf * t).sin() - b * (kf * t).cos()) / kf;
        }
        out
    }
}

enum Coords { P(Poly, Poly), T(Trig, Trig) }

impl Meridian {
    fn coords(&self) -> Coords {
        match *self {
            Meridian::Ellipse { c, u, v } => Coords::T(Trig(vec![c.0, u.0, v.0]), Trig(vec![c.1, u.1, v.1])),
            Meridian::Bezier(b) => {
                // B(t) = b0 + 3 (b1 - b0) t + 3 (b0 - 2 b1 + b2) t^2 + (b3 - 3 b2 + 3 b1 - b0) t^3.
                let pol = |k: fn(&(f64, f64)) -> f64| Poly(vec![k(&b[0]), 3.0 * (k(&b[1]) - k(&b[0])),
                    3.0 * (k(&b[0]) - 2.0 * k(&b[1]) + k(&b[2])), k(&b[3]) - 3.0 * k(&b[2]) + 3.0 * k(&b[1]) - k(&b[0])]);
                Coords::P(pol(|p| p.0), pol(|p| p.1))
            }
        }
    }
    pub fn eval(&self, t: f64) -> (f64, f64) {
        match *self {
            Meridian::Ellipse { c, u, v } => (c.0 + u.0 * t.cos() + v.0 * t.sin(), c.1 + u.1 * t.cos() + v.1 * t.sin()),
            Meridian::Bezier(b) => {
                let s = 1.0 - t;
                let w = [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t];
                (0..4).fold((0.0, 0.0), |acc, k| (acc.0 + w[k] * b[k].0, acc.1 + w[k] * b[k].1))
            }
        }
    }
    pub fn deriv(&self, t: f64) -> (f64, f64) {
        match *self {
            Meridian::Ellipse { u, v, .. } => (-u.0 * t.sin() + v.0 * t.cos(), -u.1 * t.sin() + v.1 * t.cos()),
            Meridian::Bezier(b) => {
                let s = 1.0 - t;
                let w = [-3.0 * s * s, 3.0 * s * s - 6.0 * s * t, 6.0 * s * t - 3.0 * t * t, 3.0 * t * t];
                (0..4).fold((0.0, 0.0), |acc, k| (acc.0 + w[k] * b[k].0, acc.1 + w[k] * b[k].1))
            }
        }
    }
    /// The same curve run backwards over the negated interval (ellipse) or reversed control points (Bezier, same t range).
    pub fn reversed(&self) -> Meridian {
        match *self {
            Meridian::Ellipse { c, u, v } => Meridian::Ellipse { c, u, v: (-v.0, -v.1) },
            Meridian::Bezier(b) => Meridian::Bezier([b[3], b[2], b[1], b[0]]),
        }
    }
    /// Heights negated (a negative turn's mirror).
    pub fn mirrored(&self) -> Meridian {
        match *self {
            Meridian::Ellipse { c, u, v } => Meridian::Ellipse { c: (c.0, -c.1), u: (u.0, -u.1), v: (v.0, -v.1) },
            Meridian::Bezier(b) => Meridian::Bezier(b.map(|p| (p.0, -p.1))),
        }
    }
    /// Integral over [t0, t1] of (r z' - z r') (twice the signed area swept from the origin).
    pub fn area2(&self, t0: f64, t1: f64) -> f64 {
        match self.coords() {
            Coords::P(r, z) => { let f = r.mul(&z.deriv()).add(&z.mul(&r.deriv()).scale(-1.0)); f.integral_at(t1) - f.integral_at(t0) }
            Coords::T(r, z) => { let f = r.mul(&z.deriv()).add(&z.mul(&r.deriv()).scale(-1.0)); f.integral_at(t1) - f.integral_at(t0) }
        }
    }
    /// Integral over [t0, t1] of r^2 z' (Pappus).
    pub fn r2dz(&self, t0: f64, t1: f64) -> f64 {
        match self.coords() {
            Coords::P(r, z) => { let f = r.mul(&r).mul(&z.deriv()); f.integral_at(t1) - f.integral_at(t0) }
            Coords::T(r, z) => { let f = r.mul(&r).mul(&z.deriv()); f.integral_at(t1) - f.integral_at(t0) }
        }
    }
    /// The divergence-theorem antiderivatives at t for X = P + z(t) a + r(t) rho(theta), with P.a = pa:
    /// A^(t) = integral of r z', B^(t) = integral of (r^2 z' - r r' (z + pa)), both from t = 0 (mass.rs).
    pub fn flux_hats(&self, t: f64, pa: f64) -> (f64, f64) {
        match self.coords() {
            Coords::P(r, z) => {
                let a = r.mul(&z.deriv());
                let b = r.mul(&r).mul(&z.deriv()).add(&r.mul(&r.deriv()).mul(&z.add(&Poly(vec![pa]))).scale(-1.0));
                (a.integral_at(t), b.integral_at(t))
            }
            Coords::T(r, z) => {
                let a = r.mul(&z.deriv());
                let b = r.mul(&r).mul(&z.deriv()).add(&r.mul(&r.deriv()).mul(&z.add(&Trig(vec![pa]))).scale(-1.0));
                (a.integral_at(t), b.integral_at(t))
            }
        }
    }
    /// Integral from 0 to t of r |M'| (the area density), by quadrature.
    pub fn area_hat(&self, t: f64, size: f64) -> Result<f64, String> {
        quadrature(&|s| { let (r, _) = self.eval(s); let (dr, dz) = self.deriv(s); r * dr.hypot(dz) }, 0.0, t, 1e-14 * size * size)
    }
    /// The least radius strictly inside (t0, t1), if the curve turns there towards the axis (exact stationary points).
    pub fn inner_min_radius(&self, t0: f64, t1: f64) -> Option<f64> {
        let (lo, hi) = (t0.min(t1), t0.max(t1));
        let mut best: Option<f64> = None;
        let mut try_t = |t: f64| if t > lo && t < hi { let r = self.eval(t).0; best = Some(best.map_or(r, |b: f64| b.min(r))); };
        match self.coords() {
            Coords::T(r, _) => {
                // r' = -ur sin t + vr cos t = 0 at t = atan2(vr, ur) + k pi.
                let base = r.0[2].atan2(r.0[1]);
                for k in -4..=4 { try_t(base + k as f64 * std::f64::consts::PI); }
            }
            Coords::P(r, _) => {
                let d = r.deriv();
                let (c, b, a) = (d.0[0], *d.0.get(1).unwrap_or(&0.0), *d.0.get(2).unwrap_or(&0.0));
                if a.abs() > 1e-300 {
                    let disc = b * b - 4.0 * a * c;
                    if disc >= 0.0 { for sg in [-1.0, 1.0] { try_t((-b + sg * disc.sqrt()) / (2.0 * a)); } }
                } else if b.abs() > 1e-300 { try_t(-c / b); }
            }
        }
        best
    }
    /// The parameter of a point on (or near) the curve: exact for an ellipse (its unit-circle coordinates, on the branch
    /// nearest `near`), by Newton's method from the nearest of a few samples in [0, 1] for a Bezier curve.
    pub fn param_of(&self, q: (f64, f64), near: f64) -> f64 {
        match *self {
            Meridian::Ellipse { c, u, v } => {
                let det = u.0 * v.1 - u.1 * v.0;
                let (dx, dy) = (q.0 - c.0, q.1 - c.1);
                let (cs, sn) = ((dx * v.1 - dy * v.0) / det, (u.0 * dy - u.1 * dx) / det);
                let raw = sn.atan2(cs);
                near + (raw - near + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI
            }
            Meridian::Bezier(_) => {
                let d2 = |t: f64| { let p = self.eval(t); (p.0 - q.0).powi(2) + (p.1 - q.1).powi(2) };
                let mut t = (0..=32).map(|k| k as f64 / 32.0).min_by(|a, b| d2(*a).total_cmp(&d2(*b))).unwrap();
                for _ in 0..50 {
                    let (p, d) = (self.eval(t), self.deriv(t));
                    let (ex, ey) = (p.0 - q.0, p.1 - q.1);
                    // Gauss-Newton on |B(t) - q|^2 (second derivative omitted: the point lies on or very near the curve).
                    let g = ex * d.0 + ey * d.1;
                    let hss = d.0 * d.0 + d.1 * d.1;
                    if hss <= 0.0 { break; }
                    let step = g / hss;
                    t = (t - step).clamp(0.0, 1.0);
                    if step.abs() < 1e-15 { break; }
                }
                t
            }
        }
    }
    /// Distance from (rho, h) to the curve over [t0, t1] (display refinement; Newton from samples, ends included).
    pub fn distance(&self, q: (f64, f64), t0: f64, t1: f64) -> f64 {
        let d = |t: f64| { let p = self.eval(t); (p.0 - q.0).hypot(p.1 - q.1) };
        let (lo, hi) = (t0.min(t1), t0.max(t1));
        let mut t = (0..=64).map(|k| lo + (hi - lo) * k as f64 / 64.0).min_by(|a, b| d(*a).total_cmp(&d(*b))).unwrap();
        for _ in 0..40 {
            let (p, dv) = (self.eval(t), self.deriv(t));
            let g = (p.0 - q.0) * dv.0 + (p.1 - q.1) * dv.1;
            let hss = dv.0 * dv.0 + dv.1 * dv.1;
            if hss <= 0.0 { break; }
            let next = (t - g / hss).clamp(lo, hi);
            if (next - t).abs() < 1e-15 { t = next; break; }
            t = next;
        }
        d(t).min(d(lo)).min(d(hi))
    }
    /// A conservative box over [t0, t1]: the control polygon's for a Bezier piece (convex hull), the exact extremes for an
    /// ellipse arc.
    pub fn bounds(&self, t0: f64, t1: f64) -> (f64, f64, f64, f64) {
        let (lo, hi) = (t0.min(t1), t0.max(t1));
        match *self {
            Meridian::Bezier(_) => {
                let sub = self.sub_bezier(lo, hi);
                sub.iter().fold((f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
                    |b, p| (b.0.min(p.0), b.1.min(p.1), b.2.max(p.0), b.3.max(p.1)))
            }
            Meridian::Ellipse { u, v, .. } => {
                let mut ts = vec![lo, hi];
                for base in [v.0.atan2(u.0), v.1.atan2(u.1)] {
                    for k in -6..=6 { let t = base + k as f64 * std::f64::consts::PI; if t > lo && t < hi { ts.push(t); } }
                }
                ts.iter().map(|&t| self.eval(t)).fold((f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
                    |b, p| (b.0.min(p.0), b.1.min(p.1), b.2.max(p.0), b.3.max(p.1)))
            }
        }
    }
    /// The Bezier control points of the piece over [a, b] (de Casteljau twice).
    pub fn sub_bezier(&self, a: f64, b: f64) -> [(f64, f64); 4] {
        let Meridian::Bezier(p) = *self else { panic!("not a Bezier meridian") };
        let split = |p: [(f64, f64); 4], t: f64| -> ([(f64, f64); 4], [(f64, f64); 4]) {
            let l = |x: (f64, f64), y: (f64, f64)| (x.0 + (y.0 - x.0) * t, x.1 + (y.1 - x.1) * t);
            let (p01, p12, p23) = (l(p[0], p[1]), l(p[1], p[2]), l(p[2], p[3]));
            let (p012, p123) = (l(p01, p12), l(p12, p23));
            let m = l(p012, p123);
            ([p[0], p01, p012, m], [m, p123, p23, p[3]])
        };
        let right = split(p, a).1;
        if a >= 1.0 { return [p[3]; 4]; }
        split(right, ((b - a) / (1.0 - a)).clamp(0.0, 1.0)).0
    }
    /// The largest speed |M'| over [t0, t1] (bounded from the derivative's own form).
    pub fn max_speed(&self, t0: f64, t1: f64) -> f64 {
        (0..=32).map(|k| { let t = t0 + (t1 - t0) * k as f64 / 32.0; let d = self.deriv(t); d.0.hypot(d.1) }).fold(0.0, f64::max)
    }
}

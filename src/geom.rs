//! Points, vectors and the curve and surface kinds the GC Kernel handles (no freeform, Bob 2026-10-27).

use crate::KernelError;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self { Self { x, y, z } }
    pub fn add(self, o: Self) -> Self { Self::new(self.x + o.x, self.y + o.y, self.z + o.z) }
    pub fn sub(self, o: Self) -> Self { Self::new(self.x - o.x, self.y - o.y, self.z - o.z) }
    pub fn scale(self, k: f64) -> Self { Self::new(self.x * k, self.y * k, self.z * k) }
    pub fn dot(self, o: Self) -> f64 { self.x * o.x + self.y * o.y + self.z * o.z }
    pub fn cross(self, o: Self) -> Self {
        Self::new(self.y * o.z - self.z * o.y, self.z * o.x - self.x * o.z, self.x * o.y - self.y * o.x)
    }
    pub fn length(self) -> f64 { self.dot(self).sqrt() }
    pub fn is_finite(self) -> bool { self.x.is_finite() && self.y.is_finite() && self.z.is_finite() }
    /// The unit vector along this one; non-finite or zero-length directions are refused.
    pub fn unit(self) -> Result<Self, KernelError> {
        let l = self.length();
        if !self.is_finite() || !l.is_finite() || l < 1e-14 {
            Err(KernelError::InvalidInput(format!("not a usable direction: ({}, {}, {})", self.x, self.y, self.z)))
        } else {
            Ok(self.scale(1.0 / l))
        }
    }
    pub fn distance(self, o: Self) -> f64 { self.sub(o).length() }
    /// Any unit vector at right angles to this one (self must be a unit vector).
    pub fn perpendicular(self) -> Self {
        let helper = if self.x.abs() < 0.9 { Vec3::new(1.0, 0.0, 0.0) } else { Vec3::new(0.0, 1.0, 0.0) };
        let p = self.cross(helper);
        p.scale(1.0 / p.length())
    }
}

/// A curve as the neutral format carries it.
#[derive(Clone, Debug)]
pub enum Curve {
    /// Runs a -> b.
    Line { a: Vec3, b: Vec3 },
    /// Counter-clockwise about unit `n` from angle `t0` to `t1` (t1 > t0), angles from unit `x`
    /// (perpendicular to n) toward n x x.
    Arc { c: Vec3, n: Vec3, x: Vec3, r: f64, t0: f64, t1: f64 },
    /// An exact ellipse arc (a plane cutting a cylinder obliquely): `c + u cos t + v sin t` for t from `t0` to `t1`
    /// (t0 < t1 <= t0 + 2 pi). `u` and `v` are conjugate semi-diameters: finite and independent, not necessarily
    /// unit or at right angles (Astra 2026-10-27). Wire type 3.
    Conic { c: Vec3, u: Vec3, v: Vec3, t0: f64, t1: f64 },
    /// Where a cutter cylinder crosses a host cylinder (agreed with Astra 2026-10-27), wire type 4. t runs round the
    /// cutter (axis point `c`, unit axis `a`, unit `x` perpendicular to it, y = a x x, radius `rho`):
    /// Q(t) = c + rho (cos t x + sin t y); the point is Q(t) + h(t) a on the host (axis point `p`, unit axis `b`,
    /// radius `r`), h the root of |a'|^2 h^2 + 2 (D'.a') h + |D'|^2 - r^2 = 0 (D = Q - p; primes: the parts
    /// perpendicular to b) chosen by `branch` (+1 or -1). The discriminant is positive over [t0, t1].
    CylinderPair { c: Vec3, a: Vec3, x: Vec3, rho: f64, p: Vec3, b: Vec3, r: f64, branch: f64, t0: f64, t1: f64 },
}

/// The pieces of a CylinderPair at t: point, tangent, and the discriminant S(t).
#[allow(clippy::too_many_arguments)]
pub fn pair_eval(c: Vec3, a: Vec3, x: Vec3, rho: f64, p: Vec3, b: Vec3, r: f64, branch: f64, t: f64) -> (Vec3, Vec3, f64) {
    let y = a.cross(x);
    let perp = |w: Vec3| w.sub(b.scale(w.dot(b)));
    let q = c.add(x.scale(rho * t.cos())).add(y.scale(rho * t.sin()));
    let dq = x.scale(-rho * t.sin()).add(y.scale(rho * t.cos()));
    let (dp, ddp, ap) = (perp(q.sub(p)), perp(dq), perp(a));
    let alpha = ap.dot(ap);
    let (beta, dbeta) = (dp.dot(ap), ddp.dot(ap));
    let (gamma, dgamma) = (dp.dot(dp) - r * r, 2.0 * dp.dot(ddp));
    let s = beta * beta - alpha * gamma;
    let ds = 2.0 * beta * dbeta - alpha * dgamma;
    let root = s.max(0.0).sqrt();
    let h = (-beta + branch * root) / alpha;
    let dh = (-dbeta + branch * ds / (2.0 * root)) / alpha;
    (q.add(a.scale(h)), dq.add(a.scale(dh)), s)
}

/// The discriminant S(t) of a CylinderPair as a trigonometric polynomial s0 + s1 cos t + s2 sin t + s3 cos 2t +
/// s4 sin 2t (exact from eight samples: S has degree two).
#[allow(clippy::too_many_arguments)]
pub fn pair_disc_coefficients(c: Vec3, a: Vec3, x: Vec3, rho: f64, p: Vec3, b: Vec3, r: f64) -> [f64; 5] {
    let mut k = [0.0; 5];
    for i in 0..8 {
        let t = std::f64::consts::PI * i as f64 / 4.0;
        let s = pair_eval(c, a, x, rho, p, b, r, 1.0, t).2;
        k[0] += s / 8.0;
        k[1] += s * t.cos() / 4.0;
        k[2] += s * t.sin() / 4.0;
        k[3] += s * (2.0 * t).cos() / 4.0;
        k[4] += s * (2.0 * t).sin() / 4.0;
    }
    k
}

/// A proved lower bound of S over [t0, t1]: on each step S stays above the chord between its samples less
/// max |S''| step^2 / 8, so the least of n + 1 samples less that is a bound.
pub fn pair_disc_lower_bound(k: &[f64; 5], t0: f64, t1: f64, n: usize) -> f64 {
    let s = |t: f64| k[0] + k[1] * t.cos() + k[2] * t.sin() + k[3] * (2.0 * t).cos() + k[4] * (2.0 * t).sin();
    let bend = k[1].abs() + k[2].abs() + 4.0 * (k[3].abs() + k[4].abs());
    let step = (t1 - t0) / n as f64;
    let least = (0..=n).map(|i| s(t0 + step * i as f64)).fold(f64::INFINITY, f64::min);
    least - bend * step * step / 8.0
}

impl Curve {
    pub fn point_at_angle(c: Vec3, n: Vec3, x: Vec3, r: f64, t: f64) -> Vec3 {
        let y = n.cross(x);
        c.add(x.scale(r * t.cos())).add(y.scale(r * t.sin()))
    }
    pub fn conic_point(c: Vec3, u: Vec3, v: Vec3, t: f64) -> Vec3 { c.add(u.scale(t.cos())).add(v.scale(t.sin())) }
    /// The parameter range: 0..1 for a line, the angle range otherwise.
    pub fn domain(&self) -> (f64, f64) {
        match self {
            Curve::Line { .. } => (0.0, 1.0),
            Curve::Arc { t0, t1, .. } | Curve::Conic { t0, t1, .. } | Curve::CylinderPair { t0, t1, .. } => (*t0, *t1),
        }
    }
    /// Point and derivative at parameter t (see `domain`).
    pub fn eval(&self, t: f64) -> (Vec3, Vec3) {
        match *self {
            Curve::Line { a, b } => (a.add(b.sub(a).scale(t)), b.sub(a)),
            Curve::Arc { c, n, x, r, .. } => {
                let y = n.cross(x);
                (Self::point_at_angle(c, n, x, r, t), x.scale(-r * t.sin()).add(y.scale(r * t.cos())))
            }
            Curve::Conic { c, u, v, .. } => (Self::conic_point(c, u, v, t), u.scale(-t.sin()).add(v.scale(t.cos()))),
            Curve::CylinderPair { c, a, x, rho, p, b, r, branch, .. } => {
                let (q, d, _) = pair_eval(c, a, x, rho, p, b, r, branch, t);
                (q, d)
            }
        }
    }
    pub fn start(&self) -> Vec3 {
        match self {
            Curve::Line { a, .. } => *a,
            Curve::Arc { c, n, x, r, t0, .. } => Self::point_at_angle(*c, *n, *x, *r, *t0),
            Curve::Conic { c, u, v, t0, .. } => Self::conic_point(*c, *u, *v, *t0),
            Curve::CylinderPair { t0, .. } => self.eval(*t0).0,
        }
    }
    pub fn end(&self) -> Vec3 {
        match self {
            Curve::Line { b, .. } => *b,
            Curve::Arc { c, n, x, r, t1, .. } => Self::point_at_angle(*c, *n, *x, *r, *t1),
            Curve::Conic { c, u, v, t1, .. } => Self::conic_point(*c, *u, *v, *t1),
            Curve::CylinderPair { t1, .. } => self.eval(*t1).0,
        }
    }
    pub fn translated(&self, d: Vec3) -> Curve {
        match self {
            Curve::Line { a, b } => Curve::Line { a: a.add(d), b: b.add(d) },
            Curve::Arc { c, n, x, r, t0, t1 } => Curve::Arc { c: c.add(d), n: *n, x: *x, r: *r, t0: *t0, t1: *t1 },
            Curve::Conic { c, u, v, t0, t1 } => Curve::Conic { c: c.add(d), u: *u, v: *v, t0: *t0, t1: *t1 },
            Curve::CylinderPair { c, a, x, rho, p, b, r, branch, t0, t1 } => Curve::CylinderPair {
                c: c.add(d), a: *a, x: *x, rho: *rho, p: p.add(d), b: *b, r: *r, branch: *branch, t0: *t0, t1: *t1,
            },
        }
    }
    /// The same curve run the other way: an arc becomes counter-clockwise about -n from its old end.
    pub fn reversed(&self) -> Curve {
        match self {
            Curve::Line { a, b } => Curve::Line { a: *b, b: *a },
            Curve::Arc { c, n, r, t0, t1, .. } => {
                let x = self.end().sub(*c).scale(1.0 / r);
                Curve::Arc { c: *c, n: n.scale(-1.0), x, r: *r, t0: 0.0, t1: t1 - t0 }
            }
            // c + u cos(-s) + v sin(-s) = c + u cos s - v sin s, for s from -t1 to -t0.
            Curve::Conic { c, u, v, t0, t1 } => Curve::Conic { c: *c, u: *u, v: v.scale(-1.0), t0: -t1, t1: -t0 },
            // t -> -t flips y = a x x, so a flips; h a must stay, so h and with it the branch flip (Astra).
            Curve::CylinderPair { c, a, x, rho, p, b, r, branch, t0, t1 } => Curve::CylinderPair {
                c: *c, a: a.scale(-1.0), x: *x, rho: *rho, p: *p, b: *b, r: *r, branch: -branch, t0: -t1, t1: -t0,
            },
        }
    }
}

/// A surface as the neutral format carries it. Each has its own normal; a face's sense says whether
/// the face's outward normal agrees with it.
#[derive(Clone, Debug)]
pub enum Surface {
    /// Unit normal `n`, unit in-plane reference direction `x`.
    Plane { p: Vec3, n: Vec3, x: Vec3 },
    /// Right-handed about unit `a` (normal points away from the axis); `x` unit, perpendicular to `a`.
    Cylinder { p: Vec3, a: Vec3, x: Vec3, r: f64 },
}

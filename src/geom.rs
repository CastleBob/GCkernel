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
    /// Where two cylinders cross (agreed with Astra 2026-10-27), wire type 4. t runs round the **parameter** cylinder
    /// (axis point `c`, unit axis `a`, unit `x` perpendicular to it, y = a x x, radius `rho`):
    /// Q(t) = c + rho (cos t x + sin t y); the point is Q(t) + h(t) a on the **root** cylinder (axis point `p`, unit
    /// axis `b`, radius `r`), h the root of |a'|^2 h^2 + 2 (D'.a') h + |D'|^2 - r^2 = 0 (D = Q - p; primes: the parts
    /// perpendicular to b) chosen by `branch` (+1 or -1). The discriminant is positive over [t0, t1]. Which cylinder
    /// is the parameter one is the writer's choice of chart (the one smooth over the trim), not which is the tool.
    CylinderPair { c: Vec3, a: Vec3, x: Vec3, rho: f64, p: Vec3, b: Vec3, r: f64, branch: f64, t0: f64, t1: f64 },
}

/// The pieces of a CylinderPair at t: point less c (kept small, so callers can add c - o once for any origin o),
/// tangent, and the discriminant S(t).
#[allow(clippy::too_many_arguments)]
pub fn pair_eval(c: Vec3, a: Vec3, x: Vec3, rho: f64, p: Vec3, b: Vec3, r: f64, branch: f64, t: f64) -> (Vec3, Vec3, f64) {
    pair_eval_sided(c, a, x, rho, p, b, r, branch, t, 0.0)
}

/// `pair_eval` near a trim end where the two cylinders touch (S has a double root there, `pair_trim_touches`): with
/// S tiny, S' / (2 sqrt S) is 0/0 in floating point, so it takes its one-sided limit sign(t - t*) sqrt(S'' / 2),
/// `side` being that sign (+1 near t0, -1 near t1; 0 = never).
#[allow(clippy::too_many_arguments)]
fn pair_eval_sided(c: Vec3, a: Vec3, x: Vec3, rho: f64, p: Vec3, b: Vec3, r: f64, branch: f64, t: f64, side: f64)
    -> (Vec3, Vec3, f64) {
    let y = a.cross(x);
    let perp = |w: Vec3| w.sub(b.scale(w.dot(b)));
    let round = x.scale(rho * t.cos()).add(y.scale(rho * t.sin()));
    let dq = x.scale(-rho * t.sin()).add(y.scale(rho * t.cos()));
    // perp(Q - p) as perp(c - p) + perp(Q - c): the host's axis point may sit far along its axis, and the large
    // difference is then taken once, the same at every t, instead of afresh (and noisily) at each.
    let (dp, ddp, ap) = (perp(c.sub(p)).add(perp(round)), perp(dq), perp(a));
    let alpha = ap.dot(ap);
    let (beta, dbeta) = (dp.dot(ap), ddp.dot(ap));
    // S = beta^2 - alpha (|D'|^2 - r^2) = alpha r^2 - |D' x a'|^2 (Lagrange's identity): the cross product drops the
    // part of D' along a', so a parameter-cylinder anchor far along its axis does not cancel (Astra's probe:
    // c.z = 1e8 gave S = 0 and a point half an inch off the old way).
    let cross = dp.cross(ap);
    let s = alpha * r * r - cross.dot(cross);
    let ds = -2.0 * cross.dot(ddp.cross(ap));
    let root = s.max(0.0).sqrt();
    let h = (-beta + branch * root) / alpha;
    let slope = if side != 0.0 && s <= 1e-12 * alpha * r * r {
        // S'' = -2 (|C'|^2 + C . C''), C' = perp(Q') x a', C'' = perp(Q'') x a' = -perp(Q - c) x a'.
        let (dc, ddc) = (ddp.cross(ap), perp(round).scale(-1.0).cross(ap));
        side * (-(dc.dot(dc) + cross.dot(ddc))).max(0.0).sqrt()
    } else {
        ds / (2.0 * root)
    };
    let dh = (-dbeta + branch * slope) / alpha;
    (round.add(a.scale(h)), dq.add(a.scale(dh)), s)
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

/// Bob 2026-10-28: what agrees within the tolerance is the same. So a CylinderPair trim may END where its two
/// cylinders touch (S has a double root there: a hole's wall meeting a fillet at the fillet's tangent line), when
/// - S stays above -2 alpha r tol over the whole trim (alpha = |a'|^2): with S clamped at 0 no point is more than
///   tol off the root cylinder (|D' + h a'|^2 = r^2 - S / alpha); and
/// - S is proved above tol^2 everywhere except in end neighbourhoods, each at most a quarter of the trim, at an
///   end where S is within that bound of zero AND flat: S'^2 <= 4 floor max|S''|, a double root (a touch). A simple
///   root there is the parameter cylinder's ruling grazing the root cylinder while the surfaces cross: the curve's
///   speed is unbounded in this chart, so it is refused (the writer uses the other chart).
/// Nearing tangency inside the trim (a corner within the edge) is refused: the edge must be split there.
/// Returns whether the trim touches at t0 and at t1 (both false: the ordinary positive trim).
#[allow(clippy::too_many_arguments)]
pub fn pair_trim_touches(c: Vec3, a: Vec3, x: Vec3, rho: f64, p: Vec3, b: Vec3, r: f64, t0: f64, t1: f64, tol: f64)
    -> Option<(bool, bool)> {
    let k = pair_disc_coefficients(c, a, x, rho, p, b, r);
    if pair_disc_lower_bound(&k, t0, t1, 4096) > tol * tol { return Some((false, false)); }
    if !(tol > 0.0) { return None; }
    let ap = a.sub(b.scale(a.dot(b)));
    let floor = 2.0 * ap.dot(ap) * r * tol;
    if pair_disc_lower_bound(&k, t0, t1, 4096) < -floor { return None; }
    let s = |t: f64| k[0] + k[1] * t.cos() + k[2] * t.sin() + k[3] * (2.0 * t).cos() + k[4] * (2.0 * t).sin();
    let ds = |t: f64| -k[1] * t.sin() + k[2] * t.cos() - 2.0 * k[3] * (2.0 * t).sin() + 2.0 * k[4] * (2.0 * t).cos();
    let bend = k[1].abs() + k[2].abs() + 4.0 * (k[3].abs() + k[4].abs());
    let touches = |t: f64| s(t) <= floor && ds(t).powi(2) <= 4.0 * floor * bend;
    let (touch0, touch1) = (touches(t0), touches(t1));
    if !touch0 && !touch1 { return None; }
    // The narrowest end neighbourhoods outside which S is proved above tol^2.
    (2..=40).rev().map(|n| (t1 - t0) / 2f64.powi(n)).find_map(|d| {
        let (lo, hi) = (if touch0 { t0 + d } else { t0 }, if touch1 { t1 - d } else { t1 });
        (pair_disc_lower_bound(&k, lo, hi, 4096) > tol * tol).then_some((touch0, touch1))
    })
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
            Curve::CylinderPair { c, a, x, rho, p, b, r, branch, t0, t1 } => {
                let (offset, d, _) = pair_eval_sided(c, a, x, rho, p, b, r, branch, t, Self::pair_side(t, t0, t1));
                (c.add(offset), d)
            }
        }
    }
    /// `eval` with the point given from origin `o`. For a CylinderPair the large difference c - o is taken once,
    /// so a curve far from the world origin evaluates near `o` without noise (mass quadrature).
    pub fn eval_from(&self, o: Vec3, t: f64) -> (Vec3, Vec3) {
        match *self {
            Curve::CylinderPair { c, a, x, rho, p, b, r, branch, t0, t1 } => {
                let (offset, d, _) = pair_eval_sided(c, a, x, rho, p, b, r, branch, t, Self::pair_side(t, t0, t1));
                (c.sub(o).add(offset), d)
            }
            _ => { let (q, d) = self.eval(t); (q.sub(o), d) }
        }
    }
    /// Which way a touching end lies from t: +1 nearer t0 (the trim runs on past it), -1 nearer t1.
    fn pair_side(t: f64, t0: f64, t1: f64) -> f64 { if t - t0 <= t1 - t { 1.0 } else { -1.0 } }
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
    /// An elliptic cylinder (Claude 2026-10-29, Bob: exact ellipses): p + u cos t + v sin t + s a, unit axis `a`, `u` and
    /// `v` independent semi-diameters perpendicular to `a` (conjugate, as a Conic carries them) with (u x v).a > 0; wire
    /// type 3. Its normal (v cos t - u sin t) x a points away from the axis.
    EllipticCylinder { p: Vec3, a: Vec3, u: Vec3, v: Vec3 },
    /// A circular cone (Claude 2026-10-29, Bob: tapered Pole stock): p + s a + (r + k s)(x cos t + y sin t), y = a x x,
    /// unit axis `a`, unit `x` perpendicular to it, radius `r` at p and slope `k` (radius change per unit along a).
    /// Wire type 4. Its normal (x cos t + y sin t) - k a points away from the axis (as the cylinder's, which is k = 0).
    /// Faces stay where the radius is positive (no apex in a face).
    Cone { p: Vec3, a: Vec3, x: Vec3, r: f64, k: f64 },
}

impl Surface {
    /// The elliptic cylinder's parameters (t, s) of a point on (or near) it: s along the axis, t from solving
    /// q' = u cos t + v sin t in the plane across the axis (exact for points on the surface).
    pub fn elliptic_params(p: Vec3, a: Vec3, u: Vec3, v: Vec3, q: Vec3) -> (f64, f64) {
        let d = q.sub(p);
        let s = d.dot(a);
        let (c, n) = Self::elliptic_coordinates(a, u, v, d);
        (n.atan2(c), s)
    }

    /// Coordinates (c, n) of a vector d across the axis: d = c u + n v + (d.a) a. By Cramer's rule in the plane across
    /// a, with the triple products (w x v).a, (u x w).a over (u x v).a: no Gram-determinant cancellation when u and v
    /// are nearly parallel (Codex's acceptance case).
    pub fn elliptic_coordinates(a: Vec3, u: Vec3, v: Vec3, d: Vec3) -> (f64, f64) {
        let w = d.sub(a.scale(d.dot(a)));
        let det = u.cross(v).dot(a);
        (w.cross(v).dot(a) / det, u.cross(w).dot(a) / det)
    }
}

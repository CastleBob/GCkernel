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
}

impl Curve {
    pub fn point_at_angle(c: Vec3, n: Vec3, x: Vec3, r: f64, t: f64) -> Vec3 {
        let y = n.cross(x);
        c.add(x.scale(r * t.cos())).add(y.scale(r * t.sin()))
    }
    pub fn start(&self) -> Vec3 {
        match self {
            Curve::Line { a, .. } => *a,
            Curve::Arc { c, n, x, r, t0, .. } => Self::point_at_angle(*c, *n, *x, *r, *t0),
        }
    }
    pub fn end(&self) -> Vec3 {
        match self {
            Curve::Line { b, .. } => *b,
            Curve::Arc { c, n, x, r, t1, .. } => Self::point_at_angle(*c, *n, *x, *r, *t1),
        }
    }
    pub fn translated(&self, d: Vec3) -> Curve {
        match self {
            Curve::Line { a, b } => Curve::Line { a: a.add(d), b: b.add(d) },
            Curve::Arc { c, n, x, r, t0, t1 } => Curve::Arc { c: c.add(d), n: *n, x: *x, r: *r, t0: *t0, t1: *t1 },
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

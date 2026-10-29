//! Geometry primitives for GCmaker's plain-data connection (Bob 2026-10-28: toward Rust-owned geometry; contract with
//! Astra in GCm GCmaker/Handoff/Claude-GcKernel-FilletCuts-Request.md, op 7). Points and vectors in 2D and 3D, plane
//! frames and rigid motions, and the exact rigid motion of the kernel's curves, surfaces and bodies.
//!
//! Units are inches (protocol v1), angles radians, frames right-handed (z = x cross y). `ADMIT` is an admission check
//! on supplied unit / orthonormal data only: data outside it is refused, never normalized or repaired. It is not a
//! geometric tolerance; nothing is snapped. Non-finite values are refused everywhere.

use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;

/// How far supplied unit lengths, right angles and rotation matrices may be from exact before they are refused.
pub const ADMIT: f64 = 1e-12;

fn invalid<T>(message: impl Into<String>) -> Result<T, KernelError> { Err(KernelError::InvalidInput(message.into())) }

fn finite(v: Vec3, what: &str) -> Result<Vec3, KernelError> {
    if v.is_finite() { Ok(v) } else { invalid(format!("{what} is not finite")) }
}

/// A 2D displacement (sketch-local).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec2 { pub x: f64, pub y: f64 }

impl Vec2 {
    pub const fn new(x: f64, y: f64) -> Self { Self { x, y } }
    pub fn add(self, o: Self) -> Self { Self::new(self.x + o.x, self.y + o.y) }
    pub fn sub(self, o: Self) -> Self { Self::new(self.x - o.x, self.y - o.y) }
    pub fn scale(self, k: f64) -> Self { Self::new(self.x * k, self.y * k) }
    pub fn dot(self, o: Self) -> f64 { self.x * o.x + self.y * o.y }
    /// The z component of the 3D cross product: positive when `o` is counter-clockwise of `self`.
    pub fn cross(self, o: Self) -> f64 { self.x * o.y - self.y * o.x }
    pub fn length(self) -> f64 { self.dot(self).sqrt() }
    pub fn is_finite(self) -> bool { self.x.is_finite() && self.y.is_finite() }
    /// The unit vector along this one; non-finite or zero-length directions are refused.
    pub fn unit(self) -> Result<Self, KernelError> {
        let l = self.length();
        if !self.is_finite() || !l.is_finite() || l < 1e-14 {
            invalid(format!("not a usable direction: ({}, {})", self.x, self.y))
        } else {
            Ok(self.scale(1.0 / l))
        }
    }
}

/// A 2D position (sketch-local). Positions differ by displacements; a frame maps it with its origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point2 { pub x: f64, pub y: f64 }

impl Point2 {
    pub const fn new(x: f64, y: f64) -> Self { Self { x, y } }
    pub fn add(self, v: Vec2) -> Self { Self::new(self.x + v.x, self.y + v.y) }
    pub fn sub(self, o: Self) -> Vec2 { Vec2::new(self.x - o.x, self.y - o.y) }
    pub fn is_finite(self) -> bool { self.x.is_finite() && self.y.is_finite() }
}

/// A 3D position. The kernel's bodies store positions as `Vec3`; the conversions are explicit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point3 { pub x: f64, pub y: f64, pub z: f64 }

impl Point3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self { Self { x, y, z } }
    pub const fn from_vec3(v: Vec3) -> Self { Self::new(v.x, v.y, v.z) }
    pub const fn to_vec3(self) -> Vec3 { Vec3::new(self.x, self.y, self.z) }
    pub fn add(self, v: Vec3) -> Self { Self::from_vec3(self.to_vec3().add(v)) }
    pub fn sub(self, o: Self) -> Vec3 { self.to_vec3().sub(o.to_vec3()) }
    pub fn is_finite(self) -> bool { self.to_vec3().is_finite() }
}

/// A right-handed plane frame: origin `o`, unit axes `x` and `y` at right angles, `z = x cross y`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame { o: Vec3, x: Vec3, y: Vec3 }

impl Frame {
    /// Admits only a finite origin and a finite orthonormal pair (within `ADMIT`); the axes are kept as given.
    pub fn new(o: Point3, x: Vec3, y: Vec3) -> Result<Self, KernelError> {
        let o = finite(o.to_vec3(), "a frame origin")?;
        let x = finite(x, "a frame x axis")?;
        let y = finite(y, "a frame y axis")?;
        if (x.length() - 1.0).abs() > ADMIT { return invalid(format!("a frame x axis is not unit (length {})", x.length())); }
        if (y.length() - 1.0).abs() > ADMIT { return invalid(format!("a frame y axis is not unit (length {})", y.length())); }
        if x.dot(y).abs() > ADMIT { return invalid(format!("a frame's axes are not at right angles (x.y = {})", x.dot(y))); }
        Ok(Self { o, x, y })
    }
    pub fn origin(&self) -> Point3 { Point3::from_vec3(self.o) }
    pub fn x(&self) -> Vec3 { self.x }
    pub fn y(&self) -> Vec3 { self.y }
    pub fn z(&self) -> Vec3 { self.x.cross(self.y) }
    /// A frame-local position (x, y in the plane) as a world position.
    pub fn point2(&self, p: Point2) -> Point3 { Point3::from_vec3(self.o.add(self.x.scale(p.x)).add(self.y.scale(p.y))) }
    pub fn vector2(&self, v: Vec2) -> Vec3 { self.x.scale(v.x).add(self.y.scale(v.y)) }
    /// A frame-local position (x, y, z along the frame's axes) as a world position.
    pub fn point(&self, p: Point3) -> Point3 { Point3::from_vec3(self.o.add(self.vector(p.to_vec3()))) }
    pub fn vector(&self, v: Vec3) -> Vec3 { self.x.scale(v.x).add(self.y.scale(v.y)).add(self.z().scale(v.z)) }
    pub fn local_point(&self, p: Point3) -> Point3 { Point3::from_vec3(self.local_vector(p.to_vec3().sub(self.o))) }
    pub fn local_vector(&self, v: Vec3) -> Vec3 { Vec3::new(v.dot(self.x), v.dot(self.y), v.dot(self.z())) }
}

/// A proper rigid motion: world = m p + t for positions, m v for vectors (vectors never translate). `m` is a rotation
/// (orthonormal, determinant +1): no mirror, scale or shear is representable.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rigid { m: [[f64; 3]; 3], t: Vec3 }

impl Rigid {
    pub const IDENTITY: Rigid = Rigid { m: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], t: Vec3::new(0.0, 0.0, 0.0) };

    /// Admits a supplied matrix (row-major) and translation only when the matrix is a rotation within `ADMIT`.
    pub fn from_rows(m: [[f64; 3]; 3], t: Vec3) -> Result<Self, KernelError> {
        if m.iter().flatten().any(|v| !v.is_finite()) || !t.is_finite() { return invalid("a rigid motion is not finite"); }
        for i in 0..3 {
            for j in 0..3 {
                let d: f64 = (0..3).map(|k| m[k][i] * m[k][j]).sum::<f64>() - if i == j { 1.0 } else { 0.0 };
                if d.abs() > ADMIT { return invalid(format!("a rigid motion's matrix is not orthonormal (m^T m - I = {d} at {i}, {j})")); }
            }
        }
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        if det <= 0.0 { return invalid("a rigid motion's matrix is a reflection (determinant -1)"); }
        Ok(Self { m, t })
    }
    pub fn translation(d: Vec3) -> Result<Self, KernelError> {
        Ok(Self { t: finite(d, "a translation")?, ..Self::IDENTITY })
    }
    /// Rotation by `angle` radians counter-clockwise about the unit `axis` through `about` (Rodrigues).
    pub fn rotation(axis: Vec3, angle: f64, about: Point3) -> Result<Self, KernelError> {
        let k = finite(axis, "a rotation axis")?;
        if (k.length() - 1.0).abs() > ADMIT { return invalid(format!("a rotation axis is not unit (length {})", k.length())); }
        if !angle.is_finite() { return invalid("a rotation angle is not finite"); }
        let about = finite(about.to_vec3(), "a rotation centre")?;
        let (s, c) = angle.sin_cos();
        let v = 1.0 - c;
        let m = [
            [c + k.x * k.x * v, k.x * k.y * v - k.z * s, k.x * k.z * v + k.y * s],
            [k.y * k.x * v + k.z * s, c + k.y * k.y * v, k.y * k.z * v - k.x * s],
            [k.z * k.x * v - k.y * s, k.z * k.y * v + k.x * s, c + k.z * k.z * v],
        ];
        let r = Self { m, t: Vec3::new(0.0, 0.0, 0.0) };
        Ok(Self { m, t: about.sub(r.vector(about)) })
    }
    /// The motion taking frame `a` onto frame `b` (origin to origin, axis to axis).
    pub fn from_frames(a: &Frame, b: &Frame) -> Self {
        let (ax, ay, az) = (a.x(), a.y(), a.z());
        let (bx, by, bz) = (b.x(), b.y(), b.z());
        let col = |v: Vec3, i: usize| [v.x, v.y, v.z][i];
        // m = B A^T: row i, column j = sum over axes of B_i * A_j.
        let mut m = [[0.0; 3]; 3];
        for (i, row) in m.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = col(bx, i) * col(ax, j) + col(by, i) * col(ay, j) + col(bz, i) * col(az, j);
            }
        }
        let r = Self { m, t: Vec3::new(0.0, 0.0, 0.0) };
        Self { m, t: b.o.sub(r.vector(a.o)) }
    }
    pub fn rows(&self) -> ([[f64; 3]; 3], Vec3) { (self.m, self.t) }
    pub fn vector(&self, v: Vec3) -> Vec3 {
        let m = &self.m;
        Vec3::new(m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z, m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z)
    }
    pub fn point(&self, p: Point3) -> Point3 { Point3::from_vec3(self.vector(p.to_vec3()).add(self.t)) }
    fn position(&self, p: Vec3) -> Vec3 { self.vector(p).add(self.t) }
    /// `self` first, then `next`.
    pub fn then(&self, next: &Rigid) -> Rigid {
        let mut m = [[0.0; 3]; 3];
        for (i, row) in m.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = (0..3).map(|k| next.m[i][k] * self.m[k][j]).sum();
            }
        }
        Rigid { m, t: next.position(self.t) }
    }
    pub fn inverse(&self) -> Rigid {
        let m = &self.m;
        let mt = [[m[0][0], m[1][0], m[2][0]], [m[0][1], m[1][1], m[2][1]], [m[0][2], m[1][2], m[2][2]]];
        let r = Rigid { m: mt, t: Vec3::new(0.0, 0.0, 0.0) };
        Rigid { m: mt, t: r.vector(self.t).scale(-1.0) }
    }
    pub fn frame(&self, f: &Frame) -> Frame { Frame { o: self.position(f.o), x: self.vector(f.x), y: self.vector(f.y) } }
    /// The curve moved exactly: positions move, directions turn; radii, angles, parameters and branches are unchanged.
    pub fn curve(&self, c: &Curve) -> Curve {
        match *c {
            Curve::Line { a, b } => Curve::Line { a: self.position(a), b: self.position(b) },
            Curve::Arc { c, n, x, r, t0, t1 } => Curve::Arc { c: self.position(c), n: self.vector(n), x: self.vector(x), r, t0, t1 },
            Curve::Conic { c, u, v, t0, t1 } => Curve::Conic { c: self.position(c), u: self.vector(u), v: self.vector(v), t0, t1 },
            Curve::CylinderPair { c, a, x, rho, p, b, r, branch, t0, t1 } => Curve::CylinderPair {
                c: self.position(c), a: self.vector(a), x: self.vector(x), rho, p: self.position(p), b: self.vector(b), r, branch, t0, t1,
            },
        }
    }
    pub fn surface(&self, s: &Surface) -> Surface {
        match *s {
            Surface::Plane { p, n, x } => Surface::Plane { p: self.position(p), n: self.vector(n), x: self.vector(x) },
            Surface::Cylinder { p, a, x, r } => Surface::Cylinder { p: self.position(p), a: self.vector(a), x: self.vector(x), r },
            Surface::EllipticCylinder { p, a, u, v } => Surface::EllipticCylinder {
                p: self.position(p), a: self.vector(a), u: self.vector(u), v: self.vector(v),
            },
            Surface::Cone { p, a, x, r, k } => Surface::Cone { p: self.position(p), a: self.vector(a), x: self.vector(x), r, k },
        }
    }
    /// The body moved exactly; topology (indices, senses, loops) is untouched.
    pub fn body(&self, body: &Body) -> Body {
        Body {
            vertices: body.vertices.iter().map(|v| self.position(*v)).collect(),
            edges: body.edges.iter().map(|e| Edge { curve: self.curve(&e.curve), start: e.start, end: e.end }).collect(),
            faces: body.faces.iter().map(|f| Face { surface: self.surface(&f.surface), sense: f.sense, loops: f.loops.clone() }).collect(),
        }
    }
}

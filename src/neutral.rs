//! The neutral B-rep shared with OcctBridge.dll and GCmaker (format: GCm OcctBridge/src/occt_bridge.cpp).
//! Every value is a double; counts and indices are whole numbers. Reading is strict: non-finite values,
//! counts that are not whole or exceed the data present, bad references, degenerate frames, truncated and
//! trailing data are all errors. Nothing is repaired.

use crate::geom::{Curve, Surface, Vec3};
use crate::KernelError;

#[derive(Clone, Debug)]
pub struct Edge {
    pub curve: Curve,
    pub start: usize,
    pub end: usize,
}

/// A face's loops run as the face sits in the solid (outer counter-clockwise about the face's outward
/// normal, holes clockwise); each segment is (edge index, follows the edge's own direction).
#[derive(Clone, Debug)]
pub struct Face {
    pub surface: Surface,
    pub sense: bool,
    pub loops: Vec<Vec<(usize, bool)>>,
}

#[derive(Clone, Debug, Default)]
pub struct Body {
    pub vertices: Vec<Vec3>,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
}

/// A flat profile and the slab it is swept through: closed loops, the sweep direction `n`, and the
/// sweep from `low` to `high` along it.
#[derive(Clone, Debug)]
pub struct Profile {
    pub loops: Vec<Vec<Curve>>,
    pub n: Vec3,
    pub low: f64,
    pub high: f64,
}

const PERPENDICULAR: f64 = 1e-9;

fn invalid<T>(message: impl Into<String>) -> Result<T, KernelError> { Err(KernelError::InvalidInput(message.into())) }

pub struct Reader<'a> {
    data: &'a [f64],
    at: usize,
}

impl<'a> Reader<'a> {
    pub fn new(data: &'a [f64]) -> Self { Self { data, at: 0 } }

    pub fn next(&mut self) -> Result<f64, KernelError> {
        let Some(&v) = self.data.get(self.at) else { return invalid("neutral data ended early") };
        if !v.is_finite() { return invalid(format!("non-finite value at position {}", self.at)); }
        self.at += 1;
        Ok(v)
    }

    /// A whole-number count no larger than the values still unread (every counted item takes at least one).
    pub fn count(&mut self) -> Result<usize, KernelError> {
        let at = self.at;
        let v = self.next()?;
        let remaining = (self.data.len() - self.at) as f64;
        if v < 0.0 || v.fract() != 0.0 || v > remaining {
            return invalid(format!("value {v} at position {at} is not a count within the data"));
        }
        Ok(v as usize)
    }

    /// A whole-number index below `limit`.
    pub fn index(&mut self, limit: usize, what: &str) -> Result<usize, KernelError> {
        let at = self.at;
        let v = self.next()?;
        if v < 0.0 || v.fract() != 0.0 || v >= limit as f64 {
            return invalid(format!("{what} {v} at position {at} is not one of the {limit} available"));
        }
        Ok(v as usize)
    }

    pub fn flag(&mut self) -> Result<bool, KernelError> {
        match self.next()? {
            v if v == 0.0 => Ok(false),
            v if v == 1.0 => Ok(true),
            v => invalid(format!("expected 0 or 1, found {v}")),
        }
    }

    pub fn point(&mut self) -> Result<Vec3, KernelError> { Ok(Vec3::new(self.next()?, self.next()?, self.next()?)) }
    pub fn direction(&mut self) -> Result<Vec3, KernelError> { self.point()?.unit() }

    /// Unit `n` and unit `x` at right angles: an arc or cylinder frame.
    fn frame(&mut self) -> Result<(Vec3, Vec3), KernelError> {
        let n = self.direction()?;
        let x = self.direction()?;
        if n.dot(x).abs() > PERPENDICULAR { return invalid("a frame's reference direction is not perpendicular to its axis"); }
        Ok((n, x))
    }

    pub fn curve(&mut self) -> Result<Curve, KernelError> {
        match self.next()? {
            t if t == 0.0 => {
                let (a, b) = (self.point()?, self.point()?);
                if a.distance(b) < 1e-14 { return invalid("zero-length line"); }
                Ok(Curve::Line { a, b })
            }
            t if t == 1.0 => {
                let c = self.point()?;
                let (n, x) = self.frame()?;
                let (r, t0, t1) = (self.next()?, self.next()?, self.next()?);
                if r <= 0.0 { return invalid("an arc needs a positive radius"); }
                if t1 <= t0 || t1 - t0 > 2.0 * std::f64::consts::PI + 1e-12 {
                    return invalid("an arc needs t0 < t1 <= t0 + 2 pi");
                }
                Ok(Curve::Arc { c, n, x, r, t0, t1 })
            }
            t if t == 2.0 => Err(KernelError::Unsupported("B-spline curves are not supported by the GC Kernel yet".into())),
            t => invalid(format!("unknown curve type {t}")),
        }
    }

    pub fn surface(&mut self) -> Result<Surface, KernelError> {
        match self.next()? {
            t if t == 0.0 => {
                let p = self.point()?;
                let (n, x) = self.frame()?;
                Ok(Surface::Plane { p, n, x })
            }
            t if t == 1.0 => {
                let p = self.point()?;
                let (a, x) = self.frame()?;
                let r = self.next()?;
                if r <= 0.0 { return invalid("a cylinder needs a positive radius"); }
                Ok(Surface::Cylinder { p, a, x, r })
            }
            t if t == 2.0 => Err(KernelError::Unsupported("extruded (tabulated) surfaces are not supported by the GC Kernel yet".into())),
            t => invalid(format!("unknown surface type {t}")),
        }
    }

    pub fn profile(&mut self) -> Result<Profile, KernelError> {
        if self.next()? != 1.0 { return invalid("unsupported profile version"); }
        let loop_count = self.count()?;
        if loop_count == 0 { return invalid("a profile needs at least one loop"); }
        let mut loops = Vec::with_capacity(loop_count);
        for _ in 0..loop_count {
            let curves = self.count()?;
            if curves == 0 { return invalid("an empty profile loop"); }
            let mut lp = Vec::with_capacity(curves);
            for _ in 0..curves { lp.push(self.curve()?); }
            loops.push(lp);
        }
        let n = self.direction()?;
        let (low, high) = (self.next()?, self.next()?);
        if high <= low { return invalid("a profile sweep needs low < high"); }
        Ok(Profile { loops, n, low, high })
    }

    pub fn body(&mut self) -> Result<Body, KernelError> {
        if self.next()? != 1.0 { return invalid("unsupported body version"); }
        let vertex_count = self.count()?;
        let edge_count = self.count()?;
        let face_count = self.count()?;
        let mut body = Body::default();
        for _ in 0..vertex_count { body.vertices.push(self.point()?); }
        for _ in 0..edge_count {
            let curve = self.curve()?;
            let start = self.index(vertex_count, "start vertex")?;
            let end = self.index(vertex_count, "end vertex")?;
            body.edges.push(Edge { curve, start, end });
        }
        for _ in 0..face_count {
            let surface = self.surface()?;
            let sense = self.flag()?;
            let loop_count = self.count()?;
            let mut loops = Vec::with_capacity(loop_count);
            for _ in 0..loop_count {
                let segments = self.count()?;
                let mut lp = Vec::with_capacity(segments);
                for _ in 0..segments { lp.push((self.index(edge_count, "edge")?, self.flag()?)); }
                loops.push(lp);
            }
            body.faces.push(Face { surface, sense, loops });
        }
        Ok(body)
    }

    /// Everything must have been read.
    pub fn finish(self) -> Result<(), KernelError> {
        if self.at != self.data.len() { return invalid(format!("{} unread values after the data", self.data.len() - self.at)); }
        Ok(())
    }
}

fn put(out: &mut Vec<f64>, v: Vec3) { out.extend_from_slice(&[v.x, v.y, v.z]); }

pub fn put_curve(out: &mut Vec<f64>, c: &Curve) {
    match c {
        Curve::Line { a, b } => { out.push(0.0); put(out, *a); put(out, *b); }
        Curve::Arc { c, n, x, r, t0, t1 } => {
            out.push(1.0); put(out, *c); put(out, *n); put(out, *x); out.extend_from_slice(&[*r, *t0, *t1]);
        }
    }
}

pub fn write_body(body: &Body) -> Vec<f64> {
    let mut out = vec![1.0, body.vertices.len() as f64, body.edges.len() as f64, body.faces.len() as f64];
    for v in &body.vertices { put(&mut out, *v); }
    for e in &body.edges {
        put_curve(&mut out, &e.curve);
        out.push(e.start as f64);
        out.push(e.end as f64);
    }
    for f in &body.faces {
        match &f.surface {
            Surface::Plane { p, n, x } => { out.push(0.0); put(&mut out, *p); put(&mut out, *n); put(&mut out, *x); }
            Surface::Cylinder { p, a, x, r } => { out.push(1.0); put(&mut out, *p); put(&mut out, *a); put(&mut out, *x); out.push(*r); }
        }
        out.push(if f.sense { 1.0 } else { 0.0 });
        out.push(f.loops.len() as f64);
        for lp in &f.loops {
            out.push(lp.len() as f64);
            for (edge, forward) in lp {
                out.push(*edge as f64);
                out.push(if *forward { 1.0 } else { 0.0 });
            }
        }
    }
    out
}

pub fn write_profile(p: &Profile) -> Vec<f64> {
    let mut out = vec![1.0, p.loops.len() as f64];
    for lp in &p.loops {
        out.push(lp.len() as f64);
        for c in lp { put_curve(&mut out, c); }
    }
    put(&mut out, p.n);
    out.push(p.low);
    out.push(p.high);
    out
}

//! The neutral B-rep shared with OcctBridge.dll and GCmaker (format: GCm OcctBridge/src/occt_bridge.cpp).
//! Every value is a double; counts and indices are whole numbers. Reading is strict: non-finite values,
//! counts that are not whole or exceed the data present, bad references, degenerate frames, truncated and
//! trailing data are all errors. Nothing is repaired.

use crate::geom::{pair_trim_touches, Curve, Surface, Vec3};
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
    tolerance: f64,
}

impl<'a> Reader<'a> {
    /// Strict: every CylinderPair trim must have a proved positive discriminant.
    pub fn new(data: &'a [f64]) -> Self { Self { data, at: 0, tolerance: 0.0,
        } }
    /// The request's tolerance also admits a CylinderPair trim that ends where its cylinders touch
    /// (`geom::pair_trim_touches`, Bob 2026-10-28).
    pub fn with_tolerance(data: &'a [f64], tolerance: f64) -> Self { Self { data, at: 0, tolerance,
        } }

    pub fn next(&mut self) -> Result<f64, KernelError> {
        let Some(&v) = self.data.get(self.at) else { return invalid("neutral data ended early");
        };
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
            t if t == 2.0 => Err(KernelError::Unsupported("B-spline curves are not supported by the GC Kernel yet".into(),
            )),
            t if t == 3.0 => {
                let (c, u, v) = (self.point()?, self.point()?, self.point()?);
                let (t0, t1) = (self.next()?, self.next()?);
                let (lu, lv) = (u.length(), v.length());
                if !(lu > 1e-14 && lv > 1e-14) || u.cross(v).length() <= 1e-9 * lu * lv {
                    return invalid("a conic needs two independent semi-diameters");
                }
                if t1 <= t0 || t1 - t0 > 2.0 * std::f64::consts::PI + 1e-12 {
                    return invalid("a conic needs t0 < t1 <= t0 + 2 pi");
                }
                Ok(Curve::Conic { c, u, v, t0, t1 })
            }
            t if t == 4.0 => self.pair(),
            t if t == 7.0 => {
                let Curve::CylinderPair { c, a, x, rho, p, b, r, branch, t0, t1 } = self.pair()? else { unreachable!() };
                let s = self.next()?;
                if !(s > 0.0) { return invalid("a radial cylinder pair needs a positive radius"); }
                Ok(Curve::RadialPair { c, a, x, rho, p, b, r, branch, t0, t1, s })
            }
            t if t == 8.0 => {
                let c = self.point()?;
                let (a, x) = self.frame()?;
                let rho = self.next()?;
                let p = self.point()?;
                let b = self.direction()?;
                let (r, k, branch, t0, t1) = (self.next()?, self.next()?, self.next()?, self.next()?, self.next()?);
                if rho <= 0.0 || r <= 0.0 { return invalid("a cone pair needs a positive cylinder radius and cone radius at its point"); }
                if branch != 1.0 && branch != -1.0 { return invalid("a cone pair's branch must be 1 or -1"); }
                let ap = a.sub(b.scale(a.dot(b)));
                if ap.dot(ap) - k * k * a.dot(b).powi(2) < 1e-9 { return invalid("a cone pair's cylinder runs along a generator of its cone"); }
                if t1 <= t0 || t1 - t0 > 2.0 * std::f64::consts::PI + 1e-12 {
                    return invalid("a cone pair needs t0 < t1 <= t0 + 2 pi");
                }
                let q = crate::geom::cone_pair_disc_coefficients(c, a, x, rho, p, b, r, k);
                if !(crate::geom::pair_disc_lower_bound(&q, t0, t1, 4096) > 0.0) {
                    return invalid("a cone pair's discriminant is not proved positive over its interval");
                }
                Ok(Curve::ConePair { c, a, x, rho, p, b, r, k, branch, t0, t1 })
            }
            t if t == 5.0 => {
                let p = [self.point()?, self.point()?, self.point()?, self.point()?];
                if p.iter().any(|q| !q.is_finite()) { return invalid("a Bezier curve needs finite control points"); }
                if p[0].distance(p[3]) <= self.tolerance && p.iter().all(|q| q.distance(p[0]) <= self.tolerance) {
                    return invalid("a Bezier curve has no length");
                }
                Ok(Curve::Bezier { p })
            }
            t if t == 6.0 => {
                let p = [self.point()?, self.point()?, self.point()?, self.point()?];
                let w = [self.next()?, self.next()?, self.next()?, self.next()?];
                crate::rational::validate(p, w)?;
                if p.iter().all(|q| q.distance(p[0]) <= self.tolerance) {
                    return invalid("a rational Bezier curve has no length");
                }
                Ok(Curve::RationalBezier { p, w })
            }
            t => invalid(format!("unknown curve type {t}")),
        }
    }

    /// A CylinderPair's twenty values (wire type 4, and the first twenty of type 7).
    fn pair(&mut self) -> Result<Curve, KernelError> {
        let c = self.point()?;
        let (a, x) = self.frame()?;
        let rho = self.next()?;
        let p = self.point()?;
        let b = self.direction()?;
        let (r, branch, t0, t1) = (self.next()?, self.next()?, self.next()?, self.next()?);
        if rho <= 0.0 || r <= 0.0 { return invalid("a cylinder pair needs positive radii"); }
        if branch != 1.0 && branch != -1.0 { return invalid("a cylinder pair's branch must be 1 or -1"); }
        if a.sub(b.scale(a.dot(b))).length() < 1e-9 { return invalid("a cylinder pair needs axes that are not parallel"); }
        if t1 <= t0 || t1 - t0 > 2.0 * std::f64::consts::PI + 1e-12 {
            return invalid("a cylinder pair needs t0 < t1 <= t0 + 2 pi");
        }
        // The curve exists, and its two branches stay apart, only where the discriminant is positive (or, within
        // the tolerance, where the trim ends at a touch).
        if pair_trim_touches(c, a, x, rho, p, b, r, t0, t1, self.tolerance).is_none() {
            return invalid("a cylinder pair's discriminant is not proved positive over its interval");
        }
        Ok(Curve::CylinderPair { c, a, x, rho, p, b, r, branch, t0, t1 })
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
            t if t == 2.0 => Err(KernelError::Unsupported("extruded (tabulated) surfaces are not supported by the GC Kernel yet".into(),
            )),
            t if t == 3.0 => {
                let p = self.point()?;
                let (a, u, v) = (self.point()?, self.point()?, self.point()?);
                if (a.length() - 1.0).abs() > 1e-9 { return invalid("an elliptic cylinder needs a unit axis"); }
                let (lu, lv) = (u.length(), v.length());
                // Finite throughout: lengths and the triple product too (NaN or overflow would pass a comparison).
                let turn = u.cross(v).dot(a);
                if !(lu.is_finite() && lv.is_finite() && turn.is_finite() && (lu * lv).is_finite()) {
                    return invalid("an elliptic cylinder's semi-diameters must be finite");
                }
                if !(lu > 1e-14 && lv > 1e-14) || u.dot(a).abs() > 1e-9 * lu || v.dot(a).abs() > 1e-9 * lv {
                    return invalid("an elliptic cylinder needs two semi-diameters across its axis",
                    );
                }
                if u.cross(v).dot(a) <= 1e-9 * lu * lv {
                    return invalid("an elliptic cylinder's semi-diameters must be independent and turn about its axis");
                }
                Ok(Surface::EllipticCylinder { p, a, u, v })
            }
            t if t == 4.0 => {
                let p = self.point()?;
                let (a, x) = self.frame()?;
                let (r, k) = (self.next()?, self.next()?);
                if !(r.is_finite() && k.is_finite()) || r <= 0.0 { return invalid("a cone needs a positive radius at its point and a finite slope",
                    ); }
                Ok(Surface::Cone { p, a, x, r, k })
            }
            t if t == 5.0 => {
                let c = self.point()?;
                let (a, x) = self.frame()?;
                let r = self.next()?;
                if r <= 0.0 { return invalid("a sphere needs a positive radius"); }
                Ok(Surface::Sphere { c, a, x, r })
            }
            t if t == 6.0 => {
                let c = self.point()?;
                let (a, x) = self.frame()?;
                let (big, r) = (self.next()?, self.next()?);
                // Ring, horn or spindle: a face on it is admitted only where its radius from the axis stays positive.
                if !(r > 0.0 && big > 0.0) { return invalid("a torus needs a positive centre-circle radius and tube radius",
                    ); }
                Ok(Surface::Torus { c, a, x, big, r })
            }
            t if t == 7.0 => {
                let p = self.point()?;
                let (a, x) = self.frame()?;
                let kind = self.next()?;
                let m = if kind == 2.0 {
                    let v: Vec<f64> = (0..6).map(|_| self.next()).collect::<Result<_, _>>()?;
                    let det = v[2] * v[5] - v[3] * v[4];
                    if !det.is_finite() || det.abs() <= 1e-12 * (v[2].hypot(v[3]) * v[4].hypot(v[5])).max(1e-300) {
                        return invalid("a revolved ellipse needs independent semi-diameters");
                    }
                    crate::meridian::Meridian::Ellipse { c: (v[0], v[1]), u: (v[2], v[3]), v: (v[4], v[5]),
                    }
                } else if kind == 3.0 {
                    let v: Vec<f64> = (0..8).map(|_| self.next()).collect::<Result<_, _>>()?;
                    crate::meridian::Meridian::Bezier([(v[0], v[1]), (v[2], v[3]), (v[4], v[5]), (v[6], v[7]),
                    ])
                } else {
                    return invalid(format!("unknown meridian kind {kind}"));
                };
                Ok(Surface::Revolution { p, a, x, m })
            }
            t if t == 8.0 || t == 9.0 => {
                let (nu, nv) = (self.count()?, self.count()?);
                if nu == 0 || nv == 0 || nu > 4096 || nv > 4096 || nu * nv > 4096 {
                    return invalid("invalid spline patch counts");
                }
                let u = (0..=nu).map(|_| self.next()).collect::<Result<Vec<_>, _>>()?;
                let v = (0..=nv).map(|_| self.next()).collect::<Result<Vec<_>, _>>()?;
                let mut patches = Vec::with_capacity(nu * nv);
                for _ in 0..nu * nv {
                    let mut p = [Vec3::new(0., 0., 0.); 16];
                    let mut w = [0.; 16];
                    for q in &mut p { *q = self.point()?; }
                    for x in &mut w { *x = self.next()?; }
                    patches.push(crate::spline_surface::Patch { p, w });
                }
                let pcurves = if t == 9. {
                    let n = self.count()?;
                    if n == 0 || n > 2048 {
                        return invalid("invalid pcurve loop count");
                    }
                    let mut loops = Vec::new();
                    let mut total = 0;
                    for _ in 0..n {
                        let count = self.count()?;
                        total += count;
                        if count == 0 || total > 2048 {
                            return invalid("invalid pcurve edge count");
                        }
                        loops.push(
                            (0..count)
                                .map(|_| self.curve())
                                .collect::<Result<Vec<_>, _>>()?,
                        );
                    }
                    Some(loops)
                } else {
                    None
                };
                let data = crate::spline_surface::SplineSurface { u, v, patches,
                    pcurves,
                };
                data.validate()?;
                Ok(Surface::Spline { data: Box::new(data),
                })
            }
            t if t == 10.0 => {
                let c = self.point()?;
                let (a, x) = self.frame()?;
                let rho = self.next()?;
                let p = self.point()?;
                let b = self.direction()?;
                let (r, branch) = (self.next()?, self.next()?);
                if rho <= 0.0 || r <= 0.0 { return invalid("a conoid's cylinder pair needs positive radii"); }
                if branch != 1.0 && branch != -1.0 { return invalid("a conoid's branch must be 1 or -1"); }
                if a.sub(b.scale(a.dot(b))).length() < 1e-9 { return invalid("a conoid's cylinder pair needs axes that are not parallel"); }
                // Its height is defined all the way round: the pair's discriminant proved positive over the whole turn.
                let k = crate::geom::pair_disc_coefficients(c, a, x, rho, p, b, r);
                if !(crate::geom::pair_disc_lower_bound(&k, 0.0, 2.0 * std::f64::consts::PI, 4096) > 0.0) {
                    return invalid("a conoid's height is not proved defined all the way round");
                }
                Ok(Surface::Conoid { c, a, x, rho, p, b, r, branch })
            }
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
        Ok(Profile { loops, n, low, high,
        })
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
            body.faces.push(Face { surface, sense, loops,
            });
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
        Curve::Conic { c, u, v, t0, t1 } => {
            out.push(3.0); put(out, *c); put(out, *u); put(out, *v); out.extend_from_slice(&[*t0, *t1]);
        }
        Curve::CylinderPair { c, a, x, rho, p, b, r, branch, t0, t1,
        } => {
            out.push(4.0); put(out, *c); put(out, *a); put(out, *x); out.push(*rho);
            put(out, *p); put(out, *b); out.extend_from_slice(&[*r, *branch, *t0, *t1]);
        }
        Curve::RadialPair { c, a, x, rho, p, b, r, branch, t0, t1, s } => {
            out.push(7.0); put(out, *c); put(out, *a); put(out, *x); out.push(*rho);
            put(out, *p); put(out, *b); out.extend_from_slice(&[*r, *branch, *t0, *t1, *s]);
        }
        Curve::ConePair { c, a, x, rho, p, b, r, k, branch, t0, t1 } => {
            out.push(8.0); put(out, *c); put(out, *a); put(out, *x); out.push(*rho);
            put(out, *p); put(out, *b); out.extend_from_slice(&[*r, *k, *branch, *t0, *t1]);
        }
        Curve::Bezier { p } => { out.push(5.0); for q in p { put(out, *q); } }
        Curve::RationalBezier { p, w } => { out.push(6.0); for q in p { put(out, *q); } out.extend(w); }
    }
}

pub fn put_surface(out: &mut Vec<f64>, s: &Surface) {
    match s {
        Surface::Spline { data } => {
            out.extend([
                if data.pcurves.is_some() { 9. } else { 8. }, (data.u.len() - 1) as f64, (data.v.len() - 1) as f64,
            ]);
            out.extend(&data.u);
            out.extend(&data.v);
            for p in &data.patches {
                for q in p.p { out.extend([q.x, q.y, q.z]); }
                out.extend(p.w);
            }
            if let Some(loops) = &data.pcurves {
                out.push(loops.len() as f64);
                for lp in loops {
                    out.push(lp.len() as f64);
                    for c in lp {
                        put_curve(out, c);
                    }
                }
            }
        }
        Surface::Plane { p, n, x } => { out.push(0.0); put(out, *p); put(out, *n); put(out, *x); }
        Surface::Cylinder { p, a, x, r } => { out.push(1.0); put(out, *p); put(out, *a); put(out, *x); out.push(*r); }
        Surface::EllipticCylinder { p, a, u, v } => { out.push(3.0); put(out, *p); put(out, *a); put(out, *u); put(out, *v); }
        Surface::Cone { p, a, x, r, k } => { out.push(4.0); put(out, *p); put(out, *a); put(out, *x); out.extend_from_slice(&[*r, *k]); }
        Surface::Sphere { c, a, x, r } => { out.push(5.0); put(out, *c); put(out, *a); put(out, *x); out.push(*r); }
        Surface::Torus { c, a, x, big, r } => { out.push(6.0); put(out, *c); put(out, *a); put(out, *x); out.extend_from_slice(&[*big, *r]); }
        Surface::Conoid { c, a, x, rho, p, b, r, branch } => {
            out.push(10.0); put(out, *c); put(out, *a); put(out, *x); out.push(*rho); put(out, *p); put(out, *b);
            out.extend_from_slice(&[*r, *branch]);
        }
        Surface::Revolution { p, a, x, m } => {
            out.push(7.0); put(out, *p); put(out, *a); put(out, *x);
            match m {
                crate::meridian::Meridian::Ellipse { c, u, v } => {
                    out.extend_from_slice(&[2.0, c.0, c.1, u.0, u.1, v.0, v.1])
                }
                crate::meridian::Meridian::Bezier(b) => { out.push(3.0); for q in b { out.extend_from_slice(&[q.0, q.1]); } }
            }
        }
    }
}

pub fn write_body(body: &Body) -> Vec<f64> {
    let mut out = vec![1.0, body.vertices.len() as f64, body.edges.len() as f64, body.faces.len() as f64,
    ];
    for v in &body.vertices { put(&mut out, *v); }
    for e in &body.edges {
        put_curve(&mut out, &e.curve);
        out.push(e.start as f64);
        out.push(e.end as f64);
    }
    for f in &body.faces {
        put_surface(&mut out, &f.surface);
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

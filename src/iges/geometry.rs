use super::{bad, unsupported, Document, Reader, Result};
use crate::{
    frame::Rigid,
    geom::{Curve, Surface, Vec3},
    spline_surface::{Patch, SplineSurface},
};

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rational_tensor_elevation_preserves_independent_original_basis() {
        let mut values = vec![2., 2., 2., 2., 0., 0., 0., 0., 0.];
        values.extend([0., 0., 0., 1., 1., 1., 0., 0., 0., 1., 1., 1.]);
        let weights = (0..9).map(|i| 1. + 0.13 * i as f64).collect::<Vec<_>>();
        values.extend(&weights);
        let points = (0..9)
            .map(|k| {
                Vec3::new(
                    (k % 3) as f64,
                    (k / 3) as f64,
                    0.07 * ((k % 3) * (k / 3)) as f64,
                )
            })
            .collect::<Vec<_>>();
        for p in &points {
            values.extend([p.x, p.y, p.z]);
        }
        values.extend([0., 1., 0., 1.]);
        let file = super::super::syntax::File {
            scale: 1.,
            entities: std::collections::BTreeMap::from([(
                1,
                super::super::syntax::Entity {
                    kind: 128,
                    form: 0,
                    transform: 0,
                    label: String::new(),
                    values: values
                        .into_iter()
                        .map(super::super::syntax::Value::Number)
                        .collect(),
                },
            )]),
        };
        let doc = Document {
            file,
            tolerance: 1e-6,
            used: Default::default(),
        };
        let (Surface::Spline { data }, _) = doc.surface(1).unwrap() else {
            panic!()
        };
        for u in [0., 0.17, 0.61, 1.] {
            for v in [0., 0.27, 0.89, 1.] {
                let bu = [(1. - u) * (1. - u), 2. * u * (1. - u), u * u];
                let bv = [(1. - v) * (1. - v), 2. * v * (1. - v), v * v];
                let mut p = Vec3::new(0., 0., 0.);
                let mut w = 0.;
                for j in 0..3 {
                    for i in 0..3 {
                        let h = bu[i] * bv[j] * weights[j * 3 + i];
                        p = p.add(points[j * 3 + i].scale(h));
                        w += h;
                    }
                }
                assert!(data.eval(u, v).0.distance(p.scale(1. / w)) < 1e-12);
            }
        }
    }
}
fn elevate(p: &[Vec3], w: &[f64]) -> Result<([Vec3; 4], [f64; 4])> {
    if p.len() < 2 || p.len() > 4 || p.len() != w.len() || w.iter().any(|w| *w <= 0.) {
        return bad("invalid Bezier controls/weights");
    }
    let max = w.iter().copied().fold(0., f64::max);
    let mut weights = w.to_vec();
    if weights.iter().any(|w| *w / max < 1e-12) {
        return unsupported("weight ratio exceeds 1e12");
    }
    let origin = p[0];
    let mut q = p
        .iter()
        .zip(&weights)
        .map(|(p, w)| p.sub(origin).scale(*w))
        .collect::<Vec<_>>();
    while q.len() < 4 {
        let degree = q.len();
        let mut next = vec![q[0]];
        let mut nw = vec![weights[0]];
        for i in 1..degree {
            let a = i as f64 / degree as f64;
            next.push(q[i - 1].scale(a).add(q[i].scale(1. - a)));
            nw.push(weights[i - 1] * a + weights[i] * (1. - a));
        }
        next.push(*q.last().unwrap());
        nw.push(*weights.last().unwrap());
        q = next;
        weights = nw;
    }
    let p = std::array::from_fn(|i| origin.add(q[i].scale(1. / weights[i])));
    let w = std::array::from_fn(|i| weights[i]);
    crate::rational::validate(p, w)?;
    Ok((p, w))
}
fn single_knots(r: &mut Reader<'_>, n: usize, degree: usize) -> Result<(f64, f64)> {
    if degree < 1 || degree > 3 || n != degree + 1 {
        return unsupported("IGES spline intake currently requires one clamped degree-1..3 span");
    }
    let k = (0..n + degree + 1)
        .map(|_| r.next())
        .collect::<Result<Vec<_>>>()?;
    if k[..n].iter().any(|x| *x != k[0]) || k[n..].iter().any(|x| *x != k[n]) || k[0] >= k[n] {
        return unsupported("unclamped or multispan IGES spline");
    }
    Ok((k[0], k[n]))
}
impl Document {
    pub(super) fn motion(&self, id: usize) -> Result<Rigid> {
        if id == 0 {
            return Ok(Rigid::IDENTITY);
        }
        let e = self.entity(id)?;
        if e.form != 0 || e.transform != 0 {
            return unsupported("nested or non-Cartesian transformation");
        }
        let mut r = self.reader(id, 124)?;
        let mut m = [[0.; 3]; 3];
        let mut t = [0.; 3];
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] = r.next()?;
            }
            t[i] = r.next()? * self.file.scale;
        }
        r.finish()?;
        Rigid::from_rows(m, Vec3::new(t[0], t[1], t[2]))
    }
    fn position(&self, id: usize) -> Result<Vec3> {
        let e = self.entity(id)?;
        if e.form != 0 {
            return unsupported("point form must be zero");
        }
        let mut r = self.reader(id, 116)?;
        let p = r.point()?.scale(self.file.scale);
        if r.count()? != 0 {
            return unsupported("point display symbol reference");
        }
        r.finish()?;
        Ok(self
            .motion(e.transform)?
            .point(crate::frame::Point3::from_vec3(p))
            .to_vec3())
    }
    fn direction(&self, id: usize) -> Result<Vec3> {
        let e = self.entity(id)?;
        if e.form != 0 {
            return unsupported("direction form must be zero");
        }
        let mut r = self.reader(id, 123)?;
        let v = r.point()?;
        r.finish()?;
        self.motion(e.transform)?.vector(v).unit()
    }
    pub(super) fn curve(&self, id: usize, uv: bool) -> Result<Curve> {
        let e = self.entity(id)?;
        let mut r = self.reader(id, e.kind)?;
        let scale = if uv { 1. } else { self.file.scale };
        if uv && e.transform != 0 {
            return unsupported("transformed IGES parameter curve");
        }
        let curve = match e.kind {
            110 => {
                if e.form != 0 {
                    return unsupported("unbounded IGES line form");
                }
                Curve::Line {
                    a: r.point()?.scale(scale),
                    b: r.point()?.scale(scale),
                }
            }
            100 => {
                if e.form != 0 {
                    return unsupported("circular arc form must be zero");
                }
                if uv {
                    return unsupported("circular IGES parameter curves");
                }
                let z = r.next()? * scale;
                let c = Vec3::new(r.next()? * scale, r.next()? * scale, z);
                let a = Vec3::new(r.next()? * scale, r.next()? * scale, z).sub(c);
                let b = Vec3::new(r.next()? * scale, r.next()? * scale, z).sub(c);
                let radius = a.length();
                if radius <= 0. || (radius - b.length()).abs() > self.tolerance {
                    return bad("inconsistent circular arc endpoints");
                }
                let t0 = a.y.atan2(a.x);
                let mut sweep = (b.y.atan2(b.x) - t0).rem_euclid(std::f64::consts::TAU);
                if sweep == 0. {
                    sweep = std::f64::consts::TAU;
                }
                Curve::Arc {
                    c,
                    n: Vec3::new(0., 0., 1.),
                    x: Vec3::new(1., 0., 0.),
                    r: radius,
                    t0,
                    t1: t0 + sweep,
                }
            }
            126 => {
                if e.form != 0 {
                    return unsupported("specialized spline curve forms");
                }
                let n = r.count()? + 1;
                let degree = r.count()?;
                r.flag()?;
                let closed = r.flag()?;
                let polynomial = r.flag()?;
                if r.flag()? {
                    return unsupported("periodic IGES curves");
                }
                let (a, b) = single_knots(&mut r, n, degree)?;
                let mut w = (0..n).map(|_| r.next()).collect::<Result<Vec<_>>>()?;
                if polynomial && w.iter().any(|weight| *weight != w[0]) {
                    return bad("polynomial curve has unequal weights");
                }
                let max = w.iter().copied().fold(0., f64::max);
                for w in &mut w {
                    *w /= max;
                }
                let p = (0..n)
                    .map(|_| r.point().map(|p| p.scale(scale)))
                    .collect::<Result<Vec<_>>>()?;
                let lo = r.next()?;
                let hi = r.next()?;
                r.point()?;
                if lo < a || hi > b || lo >= hi {
                    return bad("invalid spline curve parameter bounds");
                }
                let (p, w) = elevate(&p, &w)?;
                let (p, w) = crate::rational::piece(p, w, (lo - a) / (b - a), (hi - a) / (b - a));
                if closed && p[0].distance(p[3]) > self.tolerance {
                    return bad("declared closed IGES curve does not close");
                }
                Curve::RationalBezier { p, w }
            }
            _ => return unsupported(format!("curve type {}", e.kind)),
        };
        r.finish()?;
        Ok(if uv {
            curve
        } else {
            self.motion(e.transform)?.curve(&curve)
        })
    }
    pub(super) fn surface(&self, id: usize) -> Result<(Surface, Option<[f64; 4]>)> {
        let e = self.entity(id)?;
        let mut r = self.reader(id, e.kind)?;
        let (surface, domain) = match e.kind {
            190 | 192 => {
                if e.form != 1 {
                    return unsupported("unparameterized analytic surface");
                }
                let p = self.position(r.count()?)?;
                let n = self.direction(r.count()?)?;
                let radius = if e.kind == 192 {
                    Some(r.next()? * self.file.scale)
                } else {
                    None
                };
                let x = self.direction(r.count()?)?;
                if n.dot(x).abs() > 1e-12 {
                    return bad("surface axes are not perpendicular");
                }
                (
                    if let Some(radius) = radius {
                        Surface::Cylinder {
                            p,
                            a: n,
                            x,
                            r: radius,
                        }
                    } else {
                        Surface::Plane { p, n, x }
                    },
                    None,
                )
            }
            128 => {
                if e.form != 0 {
                    return unsupported("specialized tensor surface forms");
                }
                let nu = r.count()? + 1;
                let nv = r.count()? + 1;
                let du = r.count()?;
                let dv = r.count()?;
                let closedu = r.flag()?;
                let closedv = r.flag()?;
                let polynomial = r.flag()?;
                let periodicu = r.flag()?;
                let periodicv = r.flag()?;
                if closedu || closedv || periodicu || periodicv {
                    return unsupported("closed/periodic IGES tensor surfaces");
                }
                let (ua, ub) = single_knots(&mut r, nu, du)?;
                let (va, vb) = single_knots(&mut r, nv, dv)?;
                let mut weights = (0..nu * nv).map(|_| r.next()).collect::<Result<Vec<_>>>()?;
                if polynomial && weights.iter().any(|w| *w != weights[0]) {
                    return bad("polynomial surface has unequal weights");
                }
                let max = weights.iter().copied().fold(0., f64::max);
                for w in &mut weights {
                    *w /= max;
                }
                let points = (0..nu * nv)
                    .map(|_| r.point().map(|p| p.scale(self.file.scale)))
                    .collect::<Result<Vec<_>>>()?;
                let (u0, u1, v0, v1) = (r.next()?, r.next()?, r.next()?, r.next()?);
                if u0 < ua || u1 > ub || u0 >= u1 || v0 < va || v1 > vb || v0 >= v1 {
                    return bad("invalid tensor parameter bounds");
                }
                let mut strips = Vec::new();
                for j in 0..nv {
                    strips.push(elevate(
                        &points[j * nu..(j + 1) * nu],
                        &weights[j * nu..(j + 1) * nu],
                    )?);
                }
                let mut p = [Vec3::new(0., 0., 0.); 16];
                let mut w = [0.; 16];
                for i in 0..4 {
                    let (q, weight) = elevate(
                        &strips.iter().map(|(p, _)| p[i]).collect::<Vec<_>>(),
                        &strips.iter().map(|(_, w)| w[i]).collect::<Vec<_>>(),
                    )?;
                    for j in 0..4 {
                        p[4 * i + j] = q[j];
                        w[4 * i + j] = weight[j];
                    }
                }
                let patch = Patch { p, w }.piece(
                    (u0 - ua) / (ub - ua),
                    (u1 - ua) / (ub - ua),
                    (v0 - va) / (vb - va),
                    (v1 - va) / (vb - va),
                );
                let data = SplineSurface {
                    u: vec![0., 1.],
                    v: vec![0., 1.],
                    patches: vec![patch],
                    pcurves: None,
                };
                data.validate()?;
                (
                    Surface::Spline {
                        data: Box::new(data),
                    },
                    Some([u0, u1, v0, v1]),
                )
            }
            _ => return unsupported(format!("surface type {}", e.kind)),
        };
        r.finish()?;
        Ok((self.motion(e.transform)?.surface(&surface), domain))
    }
}

//! Tensor-product knot insertion in homogeneous coordinates. Independently
//! written from STEP EXPRESS attributes and the tensor-product basis identity.
use super::{invalid, unsupported, Document, Result, Value};
use crate::{
    geom::{Surface, Vec3},
    spline_surface::{Patch, SplineSurface},
};
#[derive(Clone, Copy)]
struct H {
    q: Vec3,
    w: f64,
}
impl H {
    fn mix(self, b: Self, t: f64) -> Self {
        Self {
            q: self.q.scale(1. - t).add(b.q.scale(t)),
            w: self.w * (1. - t) + b.w * t,
        }
    }
}
fn integer(v: &Value) -> Result<usize> {
    let n = v.number()?;
    if n < 1. || n > 8192. || n.fract() != 0. {
        return invalid("invalid surface degree/multiplicity");
    }
    Ok(n as usize)
}
fn knots(m: &Value, k: &Value, d: usize, n: usize) -> Result<(Vec<usize>, Vec<f64>, Vec<f64>)> {
    let m = m.list()?.iter().map(integer).collect::<Result<Vec<_>>>()?;
    let k = k
        .list()?
        .iter()
        .map(Value::number)
        .collect::<Result<Vec<_>>>()?;
    if m.len() != k.len()
        || k.len() < 2
        || k.windows(2).any(|x| x[0] >= x[1])
        || m.iter().sum::<usize>() != n + d + 1
    {
        return invalid("inconsistent surface knot vector");
    }
    if m.iter().any(|&m| m > d + 1) {
        return invalid("surface knot multiplicity exceeds degree+1");
    }
    if m[1..m.len() - 1].iter().any(|&m| m > d) {
        return invalid("discontinuous spline surface");
    }
    let expanded = k
        .iter()
        .zip(&m)
        .flat_map(|(&k, &m)| std::iter::repeat_n(k, m))
        .collect::<Vec<_>>();
    let origin = expanded[d];
    let span = expanded[n] - origin;
    if !span.is_finite() || span <= 0. {
        return invalid("unusable surface parameter domain");
    }
    let k: Vec<_> = k.iter().map(|v| (v - origin) / span).collect();
    if k.windows(2).any(|x| x[0] >= x[1]) {
        return invalid("surface knots lose resolution");
    }
    let active = k.iter().copied().filter(|x| *x >= 0. && *x <= 1.).collect();
    Ok((m, k, active))
}
fn decompose(mut p: Vec<H>, d: usize, m: &[usize], k: &[f64]) -> Vec<[H; 4]> {
    let mut u: Vec<_> = k
        .iter()
        .zip(m)
        .flat_map(|(&k, &m)| std::iter::repeat_n(k, m))
        .collect();
    // Clamp both ends of the ACTIVE domain by exact Boehm insertion, then
    // discard controls outside that domain. This also handles periodic nets.
    for x in [0., 1.] {
        while u.iter().filter(|&&v| v == x).count() < d + 1 {
            let end = u.partition_point(|&v| v <= x) - 1;
            let count = u.iter().filter(|&&v| v == x).count();
            p = (0..=p.len())
                .map(|i| {
                    if i <= end - d {
                        p[i]
                    } else if i > end - count {
                        p[i - 1]
                    } else {
                        p[i - 1].mix(p[i], (x - u[i]) / (u[i + d] - u[i]))
                    }
                })
                .collect();
            u.insert(end + 1, x);
        }
    }
    let first = u.iter().position(|&x| x == 0.).unwrap();
    let last = u.iter().rposition(|&x| x == 1.).unwrap();
    p = p[first..first + last - first - d].to_vec();
    u = u[first..=last].to_vec();
    let k = k
        .iter()
        .copied()
        .filter(|x| *x >= 0. && *x <= 1.)
        .collect::<Vec<_>>();
    for &x in &k[1..k.len() - 1] {
        while u.iter().filter(|&&v| v == x).count() < d {
            let end = u.partition_point(|&v| v <= x) - 1;
            let count = u.iter().filter(|&&v| v == x).count();
            let mut next = Vec::with_capacity(p.len() + 1);
            for i in 0..=p.len() {
                next.push(if i <= end - d {
                    p[i]
                } else if i > end - count {
                    p[i - 1]
                } else {
                    p[i - 1].mix(p[i], (x - u[i]) / (u[i + d] - u[i]))
                });
            }
            p = next;
            u.insert(end + 1, x);
        }
    }
    (0..k.len() - 1)
        .map(|i| {
            let p = &p[i * d..=i * d + d];
            match d {
                1 => [p[0], p[0].mix(p[1], 1. / 3.), p[0].mix(p[1], 2. / 3.), p[1]],
                2 => [p[0], p[0].mix(p[1], 2. / 3.), p[1].mix(p[2], 1. / 3.), p[2]],
                _ => [p[0], p[1], p[2], p[3]],
            }
        })
        .collect()
}
impl Document {
    pub(super) fn spline_surface(&self, id: u64, scale: f64) -> Result<Surface> {
        let e = self.entity(id)?;
        let (base, k) = if e.has("B_SPLINE_SURFACE") {
            e.args("REPRESENTATION_ITEM", 1)?[0].text()?;
            (
                e.args("B_SPLINE_SURFACE", 7)?,
                e.args("B_SPLINE_SURFACE_WITH_KNOTS", 5)?,
            )
        } else {
            let a = e.args("B_SPLINE_SURFACE_WITH_KNOTS", 13)?;
            a[0].text()?;
            (&a[1..8], &a[8..13])
        };
        let (du, dv) = (integer(&base[0])?, integer(&base[1])?);
        if du > 3 || dv > 3 {
            return unsupported("spline surfaces support degree 1..3 in each direction");
        }
        if !matches!(
            base[3].enumeration()?,
            "PLANE_SURF"
                | "CYLINDRICAL_SURF"
                | "CONICAL_SURF"
                | "SPHERICAL_SURF"
                | "TOROIDAL_SURF"
                | "SURF_OF_REVOLUTION"
                | "RULED_SURF"
                | "GENERALISED_CONE"
                | "QUADRIC_SURF"
                | "SURF_OF_LINEAR_EXTRUSION"
                | "UNSPECIFIED"
        ) {
            return invalid("invalid spline surface_form");
        }
        for v in &base[4..7] {
            if !matches!(v.enumeration()?, "T" | "F" | "U") {
                return invalid("invalid surface logical attribute");
            }
        }
        if base[6].enumeration()? == "T" {
            return unsupported("self-intersecting spline surface");
        }
        if !matches!(
            k[4].enumeration()?,
            "UNIFORM_KNOTS" | "QUASI_UNIFORM_KNOTS" | "PIECEWISE_BEZIER_KNOTS" | "UNSPECIFIED"
        ) {
            return invalid("invalid surface knot_spec");
        }
        let rows = base[2].list()?;
        let nu = rows.len();
        if nu <= du || nu > 8192 {
            return invalid("invalid spline surface row count");
        }
        let nv = rows[0].list()?.len();
        if nv <= dv || nu * nv > 16384 {
            return invalid("invalid spline surface grid size");
        }
        let weights = if e.has("RATIONAL_B_SPLINE_SURFACE") {
            Some(e.args("RATIONAL_B_SPLINE_SURFACE", 1)?[0].list()?)
        } else {
            None
        };
        if weights.is_some_and(|w| w.len() != nu) {
            return invalid("surface weight row count mismatch");
        }
        let mut grid = Vec::new();
        for (i, row) in rows.iter().enumerate() {
            let row = row.list()?;
            if row.len() != nv {
                return invalid("ragged spline control grid");
            }
            let wr = weights.map(|w| w[i].list()).transpose()?;
            if wr.is_some_and(|w| w.len() != nv) {
                return invalid("ragged spline weight grid");
            }
            for (j, p) in row.iter().enumerate() {
                let p = self.point(p.reference()?)?.scale(scale);
                let w = wr.map(|w| w[j].number()).transpose()?.unwrap_or(1.);
                if !p.is_finite() || !w.is_finite() || w <= 0. {
                    return invalid("surface needs finite controls and positive weights");
                }
                grid.push(H { q: p, w });
            }
        }
        let max = grid.iter().map(|h| h.w).fold(0., f64::max);
        for h in &mut grid {
            h.w /= max;
            if h.w < 1e-12 {
                return unsupported("surface weight ratio exceeds 1e12");
            }
            h.q = h.q.scale(h.w);
        }
        let (um, uk, u) = knots(&k[0], &k[2], du, nu)?;
        let (vm, vk, v) = knots(&k[1], &k[3], dv, nv)?;
        if (u.len() - 1) * (v.len() - 1) > 4096 {
            return unsupported("surface exceeds 4096 patches");
        }
        let columns: Vec<_> = (0..nv)
            .map(|j| decompose((0..nu).map(|i| grid[i * nv + j]).collect(), du, &um, &uk))
            .collect();
        let mut patches = Vec::new();
        for a in 0..u.len() - 1 {
            let strips: Vec<_> = (0..4)
                .map(|i| decompose((0..nv).map(|j| columns[j][a][i]).collect(), dv, &vm, &vk))
                .collect();
            for b in 0..v.len() - 1 {
                let h: [H; 16] = std::array::from_fn(|n| strips[n / 4][b][n % 4]);
                patches.push(Patch {
                    p: h.map(|h| h.q.scale(1. / h.w)),
                    w: h.map(|h| h.w),
                });
            }
        }
        let data = SplineSurface {
            u,
            v,
            patches,
            pcurves: None,
        };
        data.validate()?;
        for axis in 0..2 {
            if base[4 + axis].enumeration()? != "T" {
                continue;
            }
            let nu = data.u.len() - 1;
            let nv = data.v.len() - 1;
            for i in 0..if axis == 0 { nv } else { nu } {
                let a = &data.patches[if axis == 0 { i } else { i * nv }];
                let b = &data.patches[if axis == 0 {
                    (nu - 1) * nv + i
                } else {
                    i * nv + nv - 1
                }];
                let edge = |p: &Patch, end: bool| crate::geom::Curve::RationalBezier {
                    p: std::array::from_fn(|k| {
                        p.p[if axis == 0 {
                            if end {
                                12 + k
                            } else {
                                k
                            }
                        } else {
                            4 * k + if end { 3 } else { 0 }
                        }]
                    }),
                    w: std::array::from_fn(|k| {
                        p.w[if axis == 0 {
                            if end {
                                12 + k
                            } else {
                                k
                            }
                        } else {
                            4 * k + if end { 3 } else { 0 }
                        }]
                    }),
                };
                if crate::spline_curve::difference(&edge(a, false), &edge(b, true))?
                    > self.tolerance
                {
                    return invalid(
                        "declared closed spline surface does not close within tolerance",
                    );
                }
            }
        }
        Ok(Surface::Spline {
            data: Box::new(data),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn unclamped_active_domain_matches_original_basis_including_endpoints() {
        for degree in 1..=3 {
            let n = degree + 5;
            let full = (0..n + degree + 1)
                .map(|i| (i as f64 - degree as f64) / 5.)
                .collect::<Vec<_>>();
            let original = (0..n)
                .map(|i| {
                    let w = 1. + 0.17 * (i % 3) as f64;
                    H {
                        q: Vec3::new(i as f64, 0.2 * (i * i) as f64, (i % 2) as f64).scale(w),
                        w,
                    }
                })
                .collect::<Vec<_>>();
            let patches = decompose(original.clone(), degree, &vec![1; full.len()], &full);
            assert_eq!(patches.len(), 5);
            for i in 0..=100 {
                let t = i as f64 / 100.;
                let at = t.min(1. - 1e-14);
                let mut q = Vec3::new(0., 0., 0.);
                let mut weight = 0.;
                for (j, p) in original.iter().enumerate() {
                    let value = basis(j, degree, at, &full);
                    q = q.add(p.q.scale(value));
                    weight += p.w * value;
                }
                let span = ((t * 5.).floor() as usize).min(4);
                let h = patches[span];
                let c = crate::geom::Curve::RationalBezier {
                    p: h.map(|h| h.q.scale(1. / h.w)),
                    w: h.map(|h| h.w),
                };
                assert!(
                    c.eval(t * 5. - span as f64)
                        .0
                        .distance(q.scale(1. / weight))
                        < 1e-10,
                    "degree {degree}, t {t}"
                );
            }
        }
    }
    fn basis(i: usize, d: usize, t: f64, k: &[f64]) -> f64 {
        if d == 0 {
            return if k[i] <= t && t < k[i + 1] { 1. } else { 0. };
        }
        let a = k[i + d] - k[i];
        let b = k[i + d + 1] - k[i + 1];
        (if a > 0. {
            (t - k[i]) / a * basis(i, d - 1, t, k)
        } else {
            0.
        }) + (if b > 0. {
            (k[i + d + 1] - t) / b * basis(i + 1, d - 1, t, k)
        } else {
            0.
        })
    }
    #[test]
    fn tensor_conversion_matches_independent_cox_de_boor_and_partial_derivatives() {
        for du in 1..=3 {
            for dv in 1..=3 {
                for rational in [false, true] {
                    let (nu, nv) = (du + 3, dv + 3);
                    let mut source = String::new();
                    let mut rows = vec![];
                    let mut weights = vec![];
                    let mut original = vec![];
                    for i in 0..nu {
                        let mut row = vec![];
                        let mut wr = vec![];
                        for j in 0..nv {
                            let id = 1 + i * nv + j;
                            let p = Vec3::new(i as f64, j as f64, 0.1 * ((i * j) % 3) as f64);
                            let w = if rational {
                                (1. + 0.1 * i as f64) * (1. + 0.15 * j as f64)
                            } else {
                                1.
                            };
                            original.push((p, w));
                            source +=
                                &format!("#{id}=CARTESIAN_POINT('',({},{},{}));", p.x, p.y, p.z);
                            row.push(format!("#{id}"));
                            wr.push(w.to_string());
                        }
                        rows.push(format!("({})", row.join(",")));
                        weights.push(format!("({})", wr.join(",")));
                    }
                    let base = format!("{du},{dv},({}),.UNSPECIFIED.,.F.,.F.,.F.", rows.join(","));
                    let knots = format!(
                        "({},1,1,{}),({},1,1,{}),(-7.,-1.6,7.4,11.),(2.,3.2,5.2,6.),.UNSPECIFIED.",
                        du + 1,
                        du + 1,
                        dv + 1,
                        dv + 1
                    );
                    source+=&format!("#1000=(BOUNDED_SURFACE()B_SPLINE_SURFACE({base})B_SPLINE_SURFACE_WITH_KNOTS({knots})GEOMETRIC_REPRESENTATION_ITEM(){}REPRESENTATION_ITEM('')SURFACE());",if rational {format!("RATIONAL_B_SPLINE_SURFACE(({}))",weights.join(","))}else{String::new()});
                    let file=format!("ISO-10303-21;HEADER;FILE_DESCRIPTION((''),'2;1');FILE_NAME('','',(''),(''),'','','');FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));ENDSEC;DATA;{source}ENDSEC;END-ISO-10303-21;");
                    let doc = Document {
                        exchange: super::super::syntax::parse(file.as_bytes()).unwrap(),
                        tolerance: 1e-6,
                    };
                    let Surface::Spline { data } = doc.spline_surface(1000, 1.).unwrap() else {
                        panic!()
                    };
                    let k = |d| {
                        std::iter::repeat_n(0., d + 1)
                            .chain([0.3, 0.8])
                            .chain(std::iter::repeat_n(1., d + 1))
                            .collect::<Vec<_>>()
                    };
                    let (uk, vk) = (k(du), k(dv));
                    let oracle = |u: f64, v: f64| {
                        let mut q = Vec3::new(0., 0., 0.);
                        let mut w = 0.;
                        for i in 0..nu {
                            for j in 0..nv {
                                let (p, weight) = original[i * nv + j];
                                let b = basis(i, du, u, &uk) * basis(j, dv, v, &vk) * weight;
                                q = q.add(p.scale(b));
                                w += b;
                            }
                        }
                        q.scale(1. / w)
                    };
                    for u in [0.04, 0.23, 0.45, 0.76, 0.91] {
                        for v in [0.07, 0.27, 0.51, 0.74, 0.97] {
                            let (p, a, b) = data.eval(u, v);
                            assert!(p.distance(oracle(u, v)) < 1e-10, "{du},{dv},{rational}");
                            let h = 1e-6;
                            assert!(
                                a.distance(oracle(u + h, v).sub(oracle(u - h, v)).scale(0.5 / h))
                                    < 1e-7
                            );
                            assert!(
                                b.distance(oracle(u, v + h).sub(oracle(u, v - h)).scale(0.5 / h))
                                    < 1e-7
                            );
                        }
                    }
                }
            }
        }
    }
}

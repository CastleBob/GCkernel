//! Positive-weight tensor cubic surfaces. Exact control nets remain authoritative;
//! interval and Bernstein bounds govern admission and display subdivision.
use crate::{
    geom::{Curve, Vec3},
    neutral::{Body, Face},
    KernelError,
};
type Result<T> = std::result::Result<T, KernelError>;
fn refuse<T>(s: &str) -> Result<T> {
    Err(KernelError::Unsupported(s.into()))
}
#[derive(Clone, Debug)]
pub struct Patch {
    pub p: [Vec3; 16],
    pub w: [f64; 16],
}
#[derive(Clone, Debug)]
pub struct SplineSurface {
    pub u: Vec<f64>,
    pub v: Vec<f64>,
    pub patches: Vec<Patch>,
}
fn basis(t: f64) -> [f64; 4] {
    let s = 1. - t;
    [s * s * s, 3. * s * s * t, 3. * s * t * t, t * t * t]
}
fn derivative(t: f64) -> [f64; 4] {
    let s = 1. - t;
    [
        -3. * s * s,
        3. * s * s - 6. * s * t,
        6. * s * t - 3. * t * t,
        3. * t * t,
    ]
}
impl Patch {
    pub fn eval(&self, u: f64, v: f64) -> (Vec3, Vec3, Vec3) {
        let (bu, bv, du, dv) = (basis(u), basis(v), derivative(u), derivative(v));
        let max = self.w.iter().copied().fold(0., f64::max);
        let o = self.p[0];
        let mut q = [Vec3::new(0., 0., 0.); 3];
        let mut w = [0.; 3];
        for i in 0..4 {
            for j in 0..4 {
                let k = i * 4 + j;
                let b = [bu[i] * bv[j], du[i] * bv[j], bu[i] * dv[j]];
                for n in 0..3 {
                    let a = b[n] * (self.w[k] / max);
                    q[n] = q[n].add(self.p[k].sub(o).scale(a));
                    w[n] += a;
                }
            }
        }
        let p = q[0].scale(1. / w[0]);
        (
            o.add(p),
            q[1].sub(p.scale(w[1])).scale(1. / w[0]),
            q[2].sub(p.scale(w[2])).scale(1. / w[0]),
        )
    }
    pub fn piece(&self, a: f64, b: f64, c: f64, d: f64) -> Self {
        let mut out = self.clone();
        for j in 0..4 {
            let p = std::array::from_fn(|i| self.p[i * 4 + j]);
            let w = std::array::from_fn(|i| self.w[i * 4 + j]);
            let (p, w) = crate::rational::piece(p, w, a, b);
            for i in 0..4 {
                out.p[i * 4 + j] = p[i];
                out.w[i * 4 + j] = w[i];
            }
        }
        for i in 0..4 {
            let p = std::array::from_fn(|j| out.p[i * 4 + j]);
            let w = std::array::from_fn(|j| out.w[i * 4 + j]);
            let (p, w) = crate::rational::piece(p, w, c, d);
            for j in 0..4 {
                out.p[i * 4 + j] = p[j];
                out.w[i * 4 + j] = w[j];
            }
        }
        out
    }
    fn regular(&self) -> bool {
        // Rational derivative numerators Qu W - Q Wu, Qv W - Q Wv.
        // Outward-rounded interval arithmetic encloses the entire control hull.
        let max = self.w.iter().copied().fold(0., f64::max);
        let w = self.w.map(|w| w / max);
        let o = self.p[0];
        let q = self.p.map(|p| p.sub(o));
        let wr = I::hull(&w);
        let mut cross_u = [I(0., 0.); 3];
        let mut cross_v = cross_u;
        let wu: Vec<_> = (0..3)
            .flat_map(|i| (0..4).map(move |j| 3. * (w[(i + 1) * 4 + j] - w[i * 4 + j])))
            .collect();
        let wv: Vec<_> = (0..4)
            .flat_map(|i| (0..3).map(move |j| 3. * (w[i * 4 + j + 1] - w[i * 4 + j])))
            .collect();
        for k in 0..3 {
            let x: [f64; 16] = std::array::from_fn(|i| [q[i].x, q[i].y, q[i].z][k] * w[i]);
            let xu: Vec<_> = (0..3)
                .flat_map(|i| (0..4).map(move |j| 3. * (x[(i + 1) * 4 + j] - x[i * 4 + j])))
                .collect();
            let xv: Vec<_> = (0..4)
                .flat_map(|i| (0..3).map(move |j| 3. * (x[i * 4 + j + 1] - x[i * 4 + j])))
                .collect();
            cross_u[k] = I::hull_error(
                &xu,
                128. * f64::EPSILON * x.iter().map(|v| v.abs()).fold(0., f64::max),
            )
            .mul(wr)
            .sub(I::hull(&x).mul(I::hull_error(&wu, 128. * f64::EPSILON)));
            cross_v[k] = I::hull_error(
                &xv,
                128. * f64::EPSILON * x.iter().map(|v| v.abs()).fold(0., f64::max),
            )
            .mul(wr)
            .sub(I::hull(&x).mul(I::hull_error(&wv, 128. * f64::EPSILON)));
        }
        (0..3).any(|i| {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            let x = cross_u[j].mul(cross_v[k]).sub(cross_u[k].mul(cross_v[j]));
            x.0.is_finite() && x.1.is_finite() && (x.0 > 0. || x.1 < 0.)
        })
    }
    pub fn certify_regular(&self) -> Result<()> {
        let mut stack = vec![(self.clone(), 0)];
        let mut visits = 0;
        while let Some((p, d)) = stack.pop() {
            visits += 1;
            if visits > 16384 {
                return refuse("spline surface regularity exceeds subdivision budget");
            }
            if p.regular() {
                continue;
            }
            if d >= 12 {
                return refuse("spline surface Jacobian is not certified nonzero");
            }
            for (a, b) in [(0., 0.5), (0.5, 1.)] {
                for (c, e) in [(0., 0.5), (0.5, 1.)] {
                    stack.push((p.piece(a, b, c, e), d + 1));
                }
            }
        }
        Ok(())
    }
}
#[derive(Clone, Copy)]
struct I(f64, f64);
impl I {
    fn hull(v: &[f64]) -> Self {
        if v.iter().any(|x| !x.is_finite()) {
            return Self(f64::NEG_INFINITY, f64::INFINITY);
        }
        let lo = v.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = v.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let e = 64. * f64::EPSILON * lo.abs().max(hi.abs()).max(1e-300);
        Self((lo - e).next_down(), (hi + e).next_up())
    }
    fn hull_error(v: &[f64], error: f64) -> Self {
        let r = Self::hull(v);
        Self((r.0 - error).next_down(), (r.1 + error).next_up())
    }
    fn sub(self, b: Self) -> Self {
        Self((self.0 - b.1).next_down(), (self.1 - b.0).next_up())
    }
    fn mul(self, b: Self) -> Self {
        Self::hull(&[self.0 * b.0, self.0 * b.1, self.1 * b.0, self.1 * b.1])
    }
}
impl SplineSurface {
    pub fn validate(&self) -> Result<()> {
        for k in [&self.u, &self.v] {
            if k.len() < 2
                || k.len() > 4097
                || k[0] != 0.
                || *k.last().unwrap() != 1.
                || k.iter().any(|x| !x.is_finite())
                || k.windows(2).any(|x| x[0] >= x[1])
            {
                return refuse("invalid normalized spline surface knots");
            }
        }
        if self.patches.len() != (self.u.len() - 1) * (self.v.len() - 1)
            || self.patches.len() > 4096
        {
            return refuse("invalid spline patch count");
        }
        for p in &self.patches {
            if p.p.iter().any(|q| !q.is_finite()) || p.w.iter().any(|w| !w.is_finite() || *w <= 0.)
            {
                return refuse("invalid spline surface controls/weights");
            }
            let max = p.w.iter().copied().fold(0., f64::max);
            if p.w.iter().any(|w| w / max < 1e-12) {
                return refuse("surface weight ratio exceeds 1e12");
            }
            if p.p.iter().any(|q| !q.sub(p.p[0]).length().is_finite()) {
                return refuse("spline control extent overflows");
            }
            p.certify_regular()?;
        }
        // Adjacent patches must agree over the entire common isoparametric edge.
        for i in 0..self.u.len() - 1 {
            for j in 0..self.v.len() - 1 {
                let p = &self.patches[i * (self.v.len() - 1) + j];
                for (next, side, other) in [
                    (
                        if i + 1 < self.u.len() - 1 {
                            Some((i + 1) * (self.v.len() - 1) + j)
                        } else {
                            None
                        },
                        1,
                        3,
                    ),
                    (
                        if j + 1 < self.v.len() - 1 {
                            Some(i * (self.v.len() - 1) + j + 1)
                        } else {
                            None
                        },
                        2,
                        0,
                    ),
                ] {
                    if let Some(n) = next {
                        if crate::spline_curve::difference(
                            &edge(p, side),
                            &edge(&self.patches[n], other).reversed(),
                        )? > 1e-9
                        {
                            return refuse("discontinuous spline patch join");
                        }
                    }
                }
            }
        }
        Ok(())
    }
    pub fn eval(&self, u: f64, v: f64) -> (Vec3, Vec3, Vec3) {
        let i = self
            .u
            .partition_point(|&k| k <= u)
            .saturating_sub(1)
            .min(self.u.len() - 2);
        let j = self
            .v
            .partition_point(|&k| k <= v)
            .saturating_sub(1)
            .min(self.v.len() - 2);
        let (du, dv) = (self.u[i + 1] - self.u[i], self.v[j + 1] - self.v[j]);
        let (p, a, b) = self.patches[i * (self.v.len() - 1) + j]
            .eval((u - self.u[i]) / du, (v - self.v[j]) / dv);
        (p, a.scale(1. / du), b.scale(1. / dv))
    }
    pub fn mapped(&self, f: impl Fn(Vec3) -> Vec3) -> Self {
        Self {
            u: self.u.clone(),
            v: self.v.clone(),
            patches: self
                .patches
                .iter()
                .map(|p| Patch {
                    p: p.p.map(&f),
                    w: p.w,
                })
                .collect(),
        }
    }
    fn boundary(&self) -> Vec<(Curve, f64, f64)> {
        let (nu, nv) = (self.u.len() - 1, self.v.len() - 1);
        let mut out = vec![];
        for i in 0..nu {
            out.push((edge(&self.patches[i * nv], 0), self.u[i], self.u[i + 1]));
        }
        for j in 0..nv {
            out.push((
                edge(&self.patches[(nu - 1) * nv + j], 1),
                1. + self.v[j],
                1. + self.v[j + 1],
            ));
        }
        for i in (0..nu).rev() {
            out.push((
                edge(&self.patches[i * nv + nv - 1], 2),
                3. - self.u[i + 1],
                3. - self.u[i],
            ));
        }
        for j in (0..nv).rev() {
            out.push((
                edge(&self.patches[j], 3),
                4. - self.v[j + 1],
                4. - self.v[j],
            ));
        }
        out
    }
    pub fn boundary_uv(&self, q: Vec3, tol: f64) -> Result<[f64; 2]> {
        let mut values = vec![];
        for (c, a, b) in self.boundary() {
            for t in crate::spline_curve::parameters(&c, q, tol)? {
                let uv = uv(a + (b - a) * t);
                if !values
                    .iter()
                    .any(|p: &[f64; 2]| (p[0] - uv[0]).abs() + (p[1] - uv[1]).abs() < 1e-10)
                {
                    values.push(uv);
                }
            }
        }
        if values.len() != 1 {
            return refuse("spline boundary point has missing or ambiguous UV coordinates");
        }
        Ok(values[0])
    }
    pub fn check_face(&self, body: &Body, face: &Face, tol: f64) -> Result<()> {
        if let Some(chart) = crate::spline_trim::chart(self) {
            crate::spline_trim::domain(self, &chart, body, face, tol)?;
            return Ok(());
        }
        if face.loops.len() != 1 {
            return refuse("spline surface requires one natural rectangular boundary; UV holes are unsupported");
        }
        let boundary = self.boundary();
        let mut runs = vec![];
        for &(e, fwd) in &face.loops[0] {
            let c = if fwd {
                body.edges[e].curve.clone()
            } else {
                body.edges[e].curve.reversed()
            };
            let mut found = vec![];
            for (iso, a, b) in &boundary {
                // Intersections of the two declared parameter intervals. Endpoints
                // of either can split the other; this also handles a straight edge
                // spanning several natural knot spans without fitting a new curve.
                let mut pairs = vec![];
                for t in crate::spline_curve::parameters(iso, c.start(), tol)? {
                    pairs.push((0., t));
                }
                for t in crate::spline_curve::parameters(iso, c.end(), tol)? {
                    pairs.push((1., t));
                }
                for t in crate::spline_curve::parameters(&c, iso.start(), tol)? {
                    pairs.push((t, 0.));
                }
                for t in crate::spline_curve::parameters(&c, iso.end(), tol)? {
                    pairs.push((t, 1.));
                }
                pairs.sort_by(|x, y| x.0.total_cmp(&y.0));
                pairs.dedup_by(|x, y| (x.0 - y.0).abs() < 1e-11 && (x.1 - y.1).abs() < 1e-11);
                if pairs.len() < 2 {
                    continue;
                }
                if pairs.len() != 2 {
                    return refuse("ambiguous spline boundary interval correspondence");
                }
                let ((s0, t0), (s1, t1)) = (pairs[0], pairs[1]);
                if s0 == s1 || t0 == t1 || (t1 > t0) != face.sense {
                    continue;
                }
                let part = crate::spline_curve::piece(iso, t0, t1)?;
                let actual = crate::spline_curve::piece(&c, s0, s1)?;
                let error = if matches!(c, Curve::Line { .. }) {
                    let (p, _) = crate::spline_curve::controls(&part)?;
                    let (start, end) = (actual.start(), actual.end());
                    let d = end.sub(start);
                    let len = d.dot(d);
                    p.iter()
                        .map(|q| {
                            q.distance(
                                start.add(d.scale((q.sub(start).dot(d) / len).clamp(0., 1.))),
                            )
                        })
                        .fold(0., f64::max)
                } else {
                    crate::spline_curve::difference(&part, &actual)?
                };
                if error.is_finite() && error <= tol {
                    found.push((s0, s1, a + (b - a) * t0, a + (b - a) * t1));
                }
            }
            found.sort_by(|x, y| x.0.total_cmp(&y.0));
            if found.is_empty() || found[0].0 != 0. || found.last().unwrap().1 != 1. {
                return refuse("spline edge cannot be certified on the natural surface boundary");
            }
            for pair in found.windows(2) {
                if (pair[0].1 - pair[1].0).abs() > 1e-11
                    || c.eval(pair[0].1).0.distance(c.eval(pair[1].0).0) > tol
                {
                    return refuse("spline boundary correspondence has a gap or overlap");
                }
            }
            runs.extend(found.into_iter().map(|(_, _, a, b)| (a, b)));
        }
        let length: f64 = runs.iter().map(|(a, b)| (b - a).abs()).sum();
        if (length - 4.).abs() > 1e-9 {
            return refuse("spline loop does not cover the complete natural boundary");
        }
        for i in 0..runs.len() {
            let gap = (runs[i].1 - runs[(i + 1) % runs.len()].0).abs();
            if gap.min((4. - gap).abs()) > 1e-9 {
                return refuse("spline loop has a gap or overlap in parameter space");
            }
        }
        Ok(())
    }
    pub fn measure(&self, shift: Vec3, sign: f64, size: f64) -> Result<(f64, f64)> {
        let mut volume = 0.;
        let mut area = 0.;
        for p in &self.patches {
            let p = Patch {
                p: p.p.map(|p| p.add(shift)),
                w: p.w,
            };
            for which in 0..2 {
                let abs =
                    1e-10 * size.powi(if which == 0 { 3 } else { 2 }) / self.patches.len() as f64;
                let failure = std::cell::RefCell::new(None);
                let samples = std::cell::Cell::new(0usize);
                let integral = |u| {
                    if failure.borrow().is_some() {
                        return 0.;
                    }
                    let integrand = |v| {
                        if failure.borrow().is_some() {
                            return 0.;
                        }
                        samples.set(samples.get() + 1);
                        let (q, a, b) = p.eval(u, v);
                        let normal = a.cross(b);
                        let value = if which == 0 {
                            q.dot(normal) * sign / 3.
                        } else {
                            normal.length()
                        };
                        if !value.is_finite() || samples.get() > 1_000_000 {
                            *failure.borrow_mut() = Some("spline mass integration is nonfinite or exceeds its evaluation budget".to_owned());
                            0.
                        } else {
                            value
                        }
                    };
                    match crate::mass::quadrature(&integrand, 0., 1., abs) {
                        Ok(x) => x,
                        Err(e) => {
                            *failure.borrow_mut() = Some(e);
                            0.
                        }
                    }
                };
                let value = crate::mass::quadrature(&integral, 0., 1., abs)
                    .map_err(KernelError::Unsupported)?;
                if let Some(error) = failure.into_inner() {
                    return Err(KernelError::Unsupported(error));
                }
                if which == 0 {
                    volume += value;
                } else {
                    area += value;
                }
            }
        }
        Ok((volume, area))
    }
    pub fn measure_face(
        &self,
        body: &Body,
        face: &Face,
        tol: f64,
        shift: Vec3,
        sign: f64,
        size: f64,
    ) -> Result<(f64, f64)> {
        if let Some(chart) = crate::spline_trim::chart(self) {
            let loops = crate::spline_trim::domain(self, &chart, body, face, tol)?;
            crate::spline_trim::measure(self, &loops, shift, sign, size)
        } else {
            self.measure(shift, sign, size)
        }
    }
    /// Bound S - L, where L is the triangle's affine interpolant, over its
    /// entire UV bounding rectangle. Q - W L has tensor Bernstein degree (4,4).
    /// Its positive denominator is at least the least restricted weight.
    fn triangle_bound(&self, t: [[f64; 2]; 3], p: [Vec3; 3]) -> f64 {
        let [a, b, c] = t;
        let det = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if det == 0. {
            return f64::INFINITY;
        }
        let linear = |u: f64, v: f64| {
            let s = ((u - a[0]) * (c[1] - a[1]) - (v - a[1]) * (c[0] - a[0])) / det;
            let t = ((b[0] - a[0]) * (v - a[1]) - (b[1] - a[1]) * (u - a[0])) / det;
            p[1].sub(p[0]).scale(s).add(p[2].sub(p[0]).scale(t))
        };
        let lo = [
            t.iter().map(|p| p[0]).fold(1., f64::min),
            t.iter().map(|p| p[1]).fold(1., f64::min),
        ];
        let hi = [
            t.iter().map(|p| p[0]).fold(0., f64::max),
            t.iter().map(|p| p[1]).fold(0., f64::max),
        ];
        let mut bound: f64 = 0.;
        let bin3 = [1., 3., 3., 1.];
        let bin4 = [1., 4., 6., 4., 1.];
        for (i, u) in self.u.windows(2).enumerate() {
            for (j, v) in self.v.windows(2).enumerate() {
                let (a, b, c, d) = (
                    lo[0].max(u[0]),
                    hi[0].min(u[1]),
                    lo[1].max(v[0]),
                    hi[1].min(v[1]),
                );
                if a >= b || c >= d {
                    continue;
                }
                let patch = self.patches[i * (self.v.len() - 1) + j].piece(
                    (a - u[0]) / (u[1] - u[0]),
                    (b - u[0]) / (u[1] - u[0]),
                    (c - v[0]) / (v[1] - v[0]),
                    (d - v[0]) / (v[1] - v[0]),
                );
                let l = [linear(a, c), linear(a, d), linear(b, c), linear(b, d)];
                let max = patch.w.iter().copied().fold(0., f64::max);
                let w = patch.w.map(|w| w / max);
                let mut local: f64 = 0.;
                for k in 0usize..5 {
                    for n in 0usize..5 {
                        let mut q = Vec3::new(0., 0., 0.);
                        for x in k.saturating_sub(1)..=k.min(3) {
                            for y in n.saturating_sub(1)..=n.min(3) {
                                let factor = bin3[x] * bin3[y] / (bin4[k] * bin4[n]);
                                q = q.add(
                                    patch.p[x * 4 + y]
                                        .sub(p[0])
                                        .sub(l[(k - x) * 2 + n - y])
                                        .scale(w[x * 4 + y] * factor),
                                );
                            }
                        }
                        local = local.max(q.length());
                    }
                }
                bound = bound.max(local / w.iter().copied().fold(1., f64::min));
            }
        }
        let round = 512. * f64::EPSILON * p.iter().map(|q| q.length()).fold(1., f64::max);
        bound + round
    }
    /// Refine only internal edges, splitting both incident triangles together.
    /// Canonical body-edge points are retained verbatim, so adjacent faces share
    /// the same display boundary. Unproved error bounds fail at a fixed budget.
    pub fn mesh(
        &self,
        mut points: Vec<Vec3>,
        loops: &[Vec<usize>],
        sign: f64,
        chord: f64,
    ) -> Result<(Vec<Vec3>, Vec<[usize; 3]>, Vec<Vec3>)> {
        let n = points.len();
        let (mut uv, mut tris) = if let Some(chart) = crate::spline_trim::chart(self) {
            let mut uv: Vec<_> = points.iter().map(|&p| chart.point(p)).collect();
            if loops.len() == 1
                && uv
                    .iter()
                    .all(|p| p.iter().any(|x| x.abs() < 1e-12 || (x - 1.).abs() < 1e-12))
            {
                uv.push([0.5, 0.5]);
                points.push(self.eval(0.5, 0.5).0);
                (uv, (0..n).map(|i| [i, (i + 1) % n, n]).collect())
            } else {
                let chart_points: Vec<_> = uv.iter().map(|p| [p[0], p[1] * sign]).collect();
                let tris = crate::tessellate::triangulate(&chart_points, loops, false)?;
                (uv, tris)
            }
        } else {
            if n < 4 {
                return refuse("spline face needs a rectangular display boundary");
            }
            let mut uv = points
                .iter()
                .map(|&p| self.boundary_uv(p, 1e-6))
                .collect::<Result<Vec<_>>>()?;
            uv.push([0.5, 0.5]);
            points.push(self.eval(0.5, 0.5).0);
            (uv, (0..n).map(|i| [i, (i + 1) % n, n]).collect())
        };
        let boundary_edges: std::collections::HashSet<_> = loops
            .iter()
            .flat_map(|lp| {
                (0..lp.len()).map(move |i| {
                    let (a, b) = (lp[i], lp[(i + 1) % lp.len()]);
                    (a.min(b), a.max(b))
                })
            })
            .collect();
        let boundary = |a: usize, b: usize| boundary_edges.contains(&(a.min(b), a.max(b)));
        let mut at = 0;
        while at < tris.len() {
            let tri = tris[at];
            let param = tri.map(|i| uv[i]);
            let world = tri.map(|i| points[i]);
            if self.triangle_bound(param, world) <= chord {
                at += 1;
                continue;
            }
            if tris.len() > 50000 {
                return refuse("spline display bound exceeds 50000 triangles");
            }
            let edge = (0..3)
                .map(|i| (tri[i], tri[(i + 1) % 3]))
                .filter(|&(a, b)| !boundary(a, b))
                .max_by(|&(a, b), &(c, d)| {
                    let len = |a: usize, b: usize| (uv[a][0] - uv[b][0]).hypot(uv[a][1] - uv[b][1]);
                    len(a, b).total_cmp(&len(c, d))
                });
            let Some((a, b)) = edge else {
                return refuse("spline boundary chord does not meet the display bound");
            };
            let t = [(uv[a][0] + uv[b][0]) / 2., (uv[a][1] + uv[b][1]) / 2.];
            if t == uv[a] || t == uv[b] {
                return refuse("spline display refinement loses parameter resolution");
            }
            let m = points.len();
            uv.push(t);
            points.push(self.eval(t[0], t[1]).0);
            let count = tris.len();
            let mut earliest = at;
            for i in 0..count {
                if let Some(k) = (0..3).find(|&k| {
                    let (x, y) = (tris[i][k], tris[i][(k + 1) % 3]);
                    (x == a && y == b) || (x == b && y == a)
                }) {
                    let (x, y, z) = (tris[i][k], tris[i][(k + 1) % 3], tris[i][(k + 2) % 3]);
                    tris[i] = [x, m, z];
                    tris.push([m, y, z]);
                    earliest = earliest.min(i);
                }
            }
            at = earliest;
        }
        let normals = uv
            .iter()
            .map(|t| {
                let (_, a, b) = self.eval(t[0], t[1]);
                a.cross(b).unit().map(|n| n.scale(sign))
            })
            .collect::<Result<Vec<_>>>()?;
        for t in &tris {
            let [a, b, c] = t.map(|i| points[i]);
            if b.sub(a).cross(c.sub(a)).dot(normals[t[0]]) <= 0. {
                return refuse("spline display triangle reverses the surface normal");
            }
        }
        Ok((points, tris, normals))
    }
}
fn edge(p: &Patch, side: usize) -> Curve {
    let ids: [usize; 4] = match side {
        0 => [0, 4, 8, 12],
        1 => [12, 13, 14, 15],
        2 => [15, 11, 7, 3],
        _ => [3, 2, 1, 0],
    };
    Curve::RationalBezier {
        p: ids.map(|i| p.p[i]),
        w: ids.map(|i| p.w[i]),
    }
}
fn uv(t: f64) -> [f64; 2] {
    if t <= 1. {
        [t, 0.]
    } else if t <= 2. {
        [1., t - 1.]
    } else if t <= 3. {
        [3. - t, 1.]
    } else {
        [0., 4. - t]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn singular_or_folded_parameterizations_are_not_admitted() {
        for x in [[0., 0., 0., 0.], [0., 2., -1., 1.]] {
            let p = Patch {
                p: std::array::from_fn(|i| Vec3::new(x[i / 4], (i % 4) as f64, 0.)),
                w: [1.; 16],
            };
            assert!(p.certify_regular().is_err());
        }
    }
    #[test]
    fn rational_restriction_and_whole_triangle_bound_cover_independent_samples() {
        let p = Patch {
            p: std::array::from_fn(|i| {
                Vec3::new(
                    (i / 4) as f64,
                    (i % 4) as f64,
                    if i == 5 || i == 10 { 1. } else { 0. },
                )
            }),
            w: std::array::from_fn(|i| 1. + 0.03 * i as f64),
        };
        let part = p.piece(0.2, 0.9, 0.1, 0.7);
        for i in 0..21 {
            for j in 0..21 {
                let (u, v) = (i as f64 / 20., j as f64 / 20.);
                assert!(
                    part.eval(u, v)
                        .0
                        .distance(p.eval(0.2 + 0.7 * u, 0.1 + 0.6 * v).0)
                        < 1e-12
                );
            }
        }
        let data = SplineSurface {
            u: vec![0., 1.],
            v: vec![0., 1.],
            patches: vec![p],
        };
        let uv = [[0.13, 0.22], [0.76, 0.28], [0.31, 0.84]];
        let q = uv.map(|p| data.eval(p[0], p[1]).0);
        let bound = data.triangle_bound(uv, q);
        for i in 0..21 {
            for j in 0..=20 - i {
                let (b, c) = (i as f64 / 20., j as f64 / 20.);
                let a = 1. - b - c;
                let (u, v) = (
                    a * uv[0][0] + b * uv[1][0] + c * uv[2][0],
                    a * uv[0][1] + b * uv[1][1] + c * uv[2][1],
                );
                let linear = q[0].scale(a).add(q[1].scale(b)).add(q[2].scale(c));
                assert!(data.eval(u, v).0.distance(linear) <= bound);
            }
        }
    }
}

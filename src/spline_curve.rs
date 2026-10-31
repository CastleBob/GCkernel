//! Bounded point inversion of positive-weight cubic curves, from Bernstein
//! polynomial identities. No iterative closest-point guess or geometry repair.
use crate::{
    geom::{Curve, Vec3},
    KernelError,
};

pub fn controls(c: &Curve) -> Result<([Vec3; 4], [f64; 4]), KernelError> {
    match *c {
        Curve::Bezier { p } => Ok((p, [1.; 4])),
        Curve::RationalBezier { p, w } => Ok((p, w)),
        Curve::Line { a, b } => Ok((
            [
                a,
                a.add(b.sub(a).scale(1. / 3.)),
                a.add(b.sub(a).scale(2. / 3.)),
                b,
            ],
            [1.; 4],
        )),
        _ => Err(KernelError::Unsupported(
            "spline boundary requires a line or cubic spline".into(),
        )),
    }
}

/// Find every intersection with the coordinate plane through q along the
/// largest control-net extent, then require actual 3D agreement. This is a
/// conservative point inversion, not a nearest-point recovery. Cubic derivative
/// roots partition the coordinate polynomial into monotone intervals.
pub fn parameters(c: &Curve, q: Vec3, tol: f64) -> Result<Vec<f64>, KernelError> {
    let (p, w) = controls(c)?;
    crate::rational::validate(p, w)?;
    let coord = |p: Vec3, k: usize| [p.x, p.y, p.z][k];
    let extent = |k| {
        p.iter()
            .map(|&v| coord(v, k))
            .fold(f64::NEG_INFINITY, f64::max)
            - p.iter().map(|&v| coord(v, k)).fold(f64::INFINITY, f64::min)
    };
    let axis = (0..3)
        .max_by(|&a, &b| extent(a).total_cmp(&extent(b)))
        .unwrap();
    if extent(axis) <= 1e-14 {
        return Err(KernelError::Unsupported(
            "point inversion has a collapsed control net".into(),
        ));
    }
    let max = w.iter().copied().fold(0., f64::max);
    let h: [f64; 4] = std::array::from_fn(|i| (coord(p[i], axis) - coord(q, axis)) * w[i] / max);
    let hs = h.iter().map(|x| x.abs()).fold(0., f64::max);
    if !hs.is_finite() || hs == 0. {
        return Err(KernelError::Unsupported(
            "point inversion coordinate is indeterminate".into(),
        ));
    }
    let h = h.map(|x| x / hs);
    let [d, linear, b, a] = [
        h[0],
        3. * (h[1] - h[0]),
        3. * (h[2] - 2. * h[1] + h[0]),
        h[3] - 3. * h[2] + 3. * h[1] - h[0],
    ];
    let at = |t: f64| {
        if t == 0. {
            h[0]
        } else if t == 1. {
            h[3]
        } else {
            ((a * t + b) * t + linear) * t + d
        }
    };
    let mut cuts = vec![0., 1.];
    let (aa, bb, cc) = (3. * a, 2. * b, linear);
    if aa == 0. {
        if bb != 0. {
            cuts.push(-cc / bb);
        }
    } else {
        let disc = bb * bb - 4. * aa * cc;
        if disc >= 0. {
            let z = -0.5 * (bb + disc.sqrt().copysign(bb));
            if z == 0. {
                cuts.push(-bb / (2. * aa));
            } else {
                cuts.push(z / aa);
                cuts.push(cc / z);
            }
        }
    }
    cuts.retain(|t| t.is_finite() && *t >= 0. && *t <= 1.);
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let eps = 128. * f64::EPSILON * h.iter().map(|x| x.abs()).fold(1e-300, f64::max);
    let mut roots = vec![];
    for pair in cuts.windows(2) {
        let (mut lo, mut hi) = (pair[0], pair[1]);
        if at(lo) * at(hi) >= 0. {
            continue;
        }
        let sign = at(lo).is_sign_positive();
        for _ in 0..55 {
            let mid = (lo + hi) / 2.;
            if at(mid).is_sign_positive() == sign {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        roots.push((lo + hi) / 2.);
    }
    roots.retain(|&t| c_eval(c, t).distance(q) <= tol);
    // Distinct isolated roots are NEVER merged by parameter proximity. An
    // endpoint merely within tolerance is not an extra root when its nearby
    // coordinate-plane intersection has already been isolated. Exact endpoint
    // roots and stationary candidates remain distinct (or cause a refusal).
    for &t in &cuts {
        let endpoint = t == 0. || t == 1.;
        if (at(t).abs() <= eps || endpoint) && c_eval(c, t).distance(q) <= tol {
            let numerical_endpoint_alias =
                endpoint && at(t) != 0. && roots.iter().any(|r| (*r - t).abs() < 1e-11);
            if !numerical_endpoint_alias {
                roots.push(t);
            }
        }
    }
    roots.sort_by(f64::total_cmp);
    roots.dedup();
    Ok(roots)
}
fn c_eval(c: &Curve, t: f64) -> Vec3 {
    c.eval(t).0
}

pub fn piece(c: &Curve, a: f64, b: f64) -> Result<Curve, KernelError> {
    let (p, w) = controls(c)?;
    let (p, w) = crate::rational::piece(p, w, a.min(b), a.max(b));
    let c = Curve::RationalBezier { p, w };
    Ok(if a > b { c.reversed() } else { c })
}

/// Entire-span rational difference bound, using degree-six Bernstein products.
pub fn difference(a: &Curve, b: &Curve) -> Result<f64, KernelError> {
    let (p, w) = controls(a)?;
    let (q, v) = controls(b)?;
    let bin3 = [1., 3., 3., 1.];
    let bin6 = [1., 6., 15., 20., 15., 6., 1.];
    let wm = w.iter().copied().fold(0., f64::max);
    let vm = v.iter().copied().fold(0., f64::max);
    let w = w.map(|x| x / wm);
    let v = v.map(|x| x / vm);
    let mut bound: f64 = 0.;
    for k in 0..7 {
        let mut r = Vec3::new(0., 0., 0.);
        for i in 0..4 {
            if k >= i && k - i < 4 {
                let j = k - i;
                r = r.add(
                    p[i].sub(q[j])
                        .scale(w[i] * v[j] * bin3[i] * bin3[j] / bin6[k]),
                );
            }
        }
        bound = bound.max(r.length());
    }
    let extent = p
        .iter()
        .chain(&q)
        .map(|x| x.distance(p[0]))
        .fold(1., f64::max);
    let round = 256. * f64::EPSILON * extent;
    Ok((bound + round)
        / (w.iter().copied().fold(1., f64::min) * v.iter().copied().fold(1., f64::min)))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_distinct_coordinate_roots_are_checked_in_three_dimensions() {
        // x=(t-.2)(t-.5)(t-.8): three solutions at the same point.
        let c = Curve::Bezier {
            p: [-0.08, 0.14, -0.14, 0.08].map(|x| Vec3::new(x, 0., 0.)),
        };
        let roots = parameters(&c, Vec3::new(0., 0., 0.), 1e-10).unwrap();
        assert_eq!(roots.len(), 3);
        for (a, b) in roots.iter().zip([0.2, 0.5, 0.8]) {
            assert!((a - b).abs() < 1e-12);
        }
        assert!(parameters(&c, Vec3::new(0., 0.01, 0.), 1e-6)
            .unwrap()
            .is_empty());
    }
    #[test]
    fn rational_point_inversion_and_subdivision_preserve_parameters() {
        let c = Curve::RationalBezier {
            p: [
                Vec3::new(-2., 0., 1.),
                Vec3::new(0., 1., 0.),
                Vec3::new(1., 2., 2.),
                Vec3::new(4., 3., 1.),
            ],
            w: [1., 0.3, 2., 1.],
        };
        for t in [0., 0.01, 0.2, 0.5, 0.91, 1.] {
            let roots = parameters(&c, c.eval(t).0, 1e-10).unwrap();
            assert_eq!(roots.len(), 1);
            assert!((roots[0] - t).abs() < 1e-11);
        }
        let p = piece(&c, 0.2, 0.8).unwrap();
        for i in 0..101 {
            let t = i as f64 / 100.;
            assert!(p.eval(t).0.distance(c.eval(0.2 + 0.6 * t).0) < 1e-12);
        }
    }
}

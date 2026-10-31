//! Positive-weight rational Bezier mathematics. Independent implementation of
//! homogeneous Bernstein evaluation and de Casteljau subdivision (MTU notes).
use crate::{geom::Vec3, KernelError};

pub fn validate(p: [Vec3; 4], w: [f64; 4]) -> Result<(), KernelError> {
    let max = w.iter().copied().fold(0.0, f64::max);
    if p.iter().any(|q| !q.is_finite()) || w.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(KernelError::InvalidInput(
            "rational Bezier needs finite points and positive finite weights".into(),
        ));
    }
    if w.iter().any(|v| v / max < 1e-12) {
        return Err(KernelError::Unsupported(
            "rational weight ratio exceeds 1e12".into(),
        ));
    }
    Ok(())
}

pub fn eval(p: [Vec3; 4], w: [f64; 4], t: f64) -> (Vec3, Vec3) {
    let max = w.iter().copied().fold(0.0, f64::max);
    let w = w.map(|v| v / max);
    let s = 1.0 - t;
    let b = [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t];
    let db = [
        -3.0 * s * s,
        3.0 * s * s - 6.0 * s * t,
        6.0 * s * t - 3.0 * t * t,
        3.0 * t * t,
    ];
    let den: f64 = (0..4).map(|i| b[i] * w[i]).sum();
    let dd: f64 = (0..4).map(|i| db[i] * w[i]).sum();
    let zero = Vec3::new(0.0, 0.0, 0.0);
    // Work relative to one control point to avoid cancellation after translation.
    let q = (0..4).fold(zero, |v, i| v.add(p[i].sub(p[0]).scale(b[i] * w[i] / den)));
    let d = (0..4).fold(zero, |v, i| v.add(p[i].sub(p[0]).scale(db[i] * w[i] / den)));
    (p[0].add(q), d.sub(q.scale(dd / den)))
}

/// Split in homogeneous coordinates, retaining positive weights and exact shape.
pub fn split(p: [Vec3; 4], w: [f64; 4], t: f64) -> [([Vec3; 4], [f64; 4]); 2] {
    let mut p = p;
    let mut w = w;
    let (mut lp, mut rp) = ([p[0]; 4], [p[3]; 4]);
    let (mut lw, mut rw) = ([w[0]; 4], [w[3]; 4]);
    for level in 1..4 {
        for i in 0..4 - level {
            let nw = (1.0 - t) * w[i] + t * w[i + 1];
            p[i] = p[i].add(p[i + 1].sub(p[i]).scale(t * w[i + 1] / nw));
            w[i] = nw;
        }
        lp[level] = p[0];
        lw[level] = w[0];
        rp[3 - level] = p[3 - level];
        rw[3 - level] = w[3 - level];
    }
    [(lp, lw), (rp, rw)]
}

pub fn piece(p: [Vec3; 4], w: [f64; 4], a: f64, b: f64) -> ([Vec3; 4], [f64; 4]) {
    if a >= 1.0 {
        return ([p[3]; 4], [w[3]; 4]);
    }
    let (p, w) = split(p, w, a)[1];
    split(p, w, (b - a) / (1.0 - a))[0]
}

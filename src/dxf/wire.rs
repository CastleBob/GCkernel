//! Ops 14 READ DXF (also SVG, told apart by content) and 15 WRITE DXF on the numeric wire. Text (a path, a layer name, a message) travels as its UTF-8
//! byte count followed by one value per byte.
//!
//! Op 14 request: [1, mode (0 exact, 1 lines and arcs, 2 header and layers only), fit tolerance (drawing units), path].
//! Reply: [1, $INSUNITS or -1, $MEASUREMENT or -1, $ACADVER text, layer count, per layer (colour, flags 1 frozen
//!        2 off 4 locked, name), curve count, per curve (kind, layer, colour, group, geometry), skipped count, per
//!        entry (type name, count), warning count, per warning (text)].
//! Op 15 request: [1, $INSUNITS or -1, fit tolerance, path, layer count, per layer (colour, flags, name), curve count,
//!        per curve (kind, layer, colour, geometry)]. Reply: [bytes written]. The file is replaced only once written.
//! Geometry by kind: 0 line (ax ay bx by), 1 arc (cx cy r start sweep), 2 circle (cx cy r), 3 ellipse (cx cy ux uy
//! vx vy t0 t1), 4 spline (degree, knot count, knots, point count, x y per point, weight count, weights),
//! 5 polyline (closed, count, x y bulge per vertex; written only).

use super::{fit, read, read_header, write_r12, Curve, Drawing, Geometry, Layer};
use crate::KernelError;

struct Cursor<'a> {
    data: &'a [f64],
    at: usize,
}

impl Cursor<'_> {
    fn next(&mut self) -> Result<f64, KernelError> {
        let v = *self.data.get(self.at).ok_or_else(|| KernelError::InvalidInput("the DXF request ended early".into()))?;
        self.at += 1;
        Ok(v)
    }

    fn count(&mut self, what: &str) -> Result<usize, KernelError> {
        let v = self.next()?;
        if v < 0.0 || v.fract() != 0.0 || v > (self.data.len() - self.at.min(self.data.len())) as f64 * 8.0 + 8.0 {
            return Err(KernelError::InvalidInput(format!("bad {what} count {v}")));
        }
        Ok(v as usize)
    }

    fn text(&mut self) -> Result<String, KernelError> {
        let n = self.count("text")?;
        let mut bytes = Vec::with_capacity(n);
        for _ in 0..n {
            let b = self.next()?;
            if !(0.0..=255.0).contains(&b) || b.fract() != 0.0 {
                return Err(KernelError::InvalidInput("text byte out of range".into()));
            }
            bytes.push(b as u8);
        }
        String::from_utf8(bytes).map_err(|_| KernelError::InvalidInput("text is not UTF-8".into()))
    }

    fn finish(&self) -> Result<(), KernelError> {
        if self.at == self.data.len() { Ok(()) } else {
            Err(KernelError::InvalidInput(format!("the DXF request has {} values left over", self.data.len() - self.at)))
        }
    }
}

fn put_text(out: &mut Vec<f64>, s: &str) {
    out.push(s.len() as f64);
    out.extend(s.bytes().map(|b| b as f64));
}

fn flags(l: &Layer) -> f64 {
    (u8::from(l.frozen) | (u8::from(l.off) << 1) | (u8::from(l.locked) << 2)) as f64
}

pub fn read_op(data: &[f64]) -> Result<Vec<f64>, KernelError> {
    let mut c = Cursor { data, at: 0 };
    if c.next()? != 1.0 {
        return Err(KernelError::Unsupported("DXF read request version".into()));
    }
    let mode = c.next()?;
    let tol = c.next()?;
    let path = c.text()?;
    c.finish()?;
    let bytes = std::fs::read(&path).map_err(|e| KernelError::InvalidInput(format!("could not read '{path}': {e}")))?;
    let mut drawing = if super::svg::looks_like(&bytes) {
        let mut d = super::svg::read(&bytes)?;
        if mode == 2.0 { d.curves.clear(); }
        d
    } else if mode == 2.0 {
        read_header(&bytes)?
    } else {
        read(&bytes)?
    };
    if mode == 1.0 {
        let mut fitted = Vec::with_capacity(drawing.curves.len());
        for curve in &drawing.curves {
            for g in fit::lines_and_arcs(&curve.geometry, tol)? {
                fitted.push(Curve { geometry: g, ..curve.clone() });
            }
        }
        drawing.curves = fitted;
    } else if mode != 0.0 && mode != 2.0 {
        return Err(KernelError::InvalidInput(format!("DXF read mode {mode}")));
    }

    let mut out = vec![1.0, drawing.insunits.map_or(-1.0, f64::from), drawing.measurement.map_or(-1.0, f64::from)];
    put_text(&mut out, &drawing.version);
    out.push(drawing.layers.len() as f64);
    for l in &drawing.layers {
        out.push(l.color as f64);
        out.push(flags(l));
        put_text(&mut out, &l.name);
    }
    out.push(drawing.curves.len() as f64);
    for curve in &drawing.curves {
        let (kind, values): (f64, Vec<f64>) = match &curve.geometry {
            Geometry::Line { a, b } => (0.0, vec![a[0], a[1], b[0], b[1]]),
            Geometry::Arc { c, r, start, sweep } => (1.0, vec![c[0], c[1], *r, *start, *sweep]),
            Geometry::Circle { c, r } => (2.0, vec![c[0], c[1], *r]),
            Geometry::Ellipse { c, u, v, t0, t1 } => (3.0, vec![c[0], c[1], u[0], u[1], v[0], v[1], *t0, *t1]),
            Geometry::Spline { degree, knots, points, weights } => {
                let mut v = vec![*degree as f64, knots.len() as f64];
                v.extend(knots);
                v.push(points.len() as f64);
                points.iter().for_each(|p| v.extend(p));
                v.push(weights.len() as f64);
                v.extend(weights);
                (4.0, v)
            }
            Geometry::Polyline { .. } => unreachable!("the reader gives polyline pieces"),
        };
        out.extend([kind, curve.layer as f64, curve.color as f64, curve.group as f64]);
        out.extend(values);
    }
    out.push(drawing.skipped.len() as f64);
    for (kind, n) in &drawing.skipped {
        put_text(&mut out, kind);
        out.push(*n as f64);
    }
    out.push(drawing.warnings.len() as f64);
    for w in &drawing.warnings {
        put_text(&mut out, w);
    }
    Ok(out)
}

pub fn write_op(data: &[f64]) -> Result<Vec<f64>, KernelError> {
    let mut c = Cursor { data, at: 0 };
    if c.next()? != 1.0 {
        return Err(KernelError::Unsupported("DXF write request version".into()));
    }
    let units = c.next()?;
    let tol = c.next()?;
    let path = c.text()?;
    let mut drawing = Drawing { insunits: (units >= 0.0).then_some(units as i32), ..Drawing::default() };
    for _ in 0..c.count("layer")? {
        let color = c.next()? as i32;
        let f = c.next()? as u32;
        let name = c.text()?;
        drawing.layers.push(Layer { name, color, frozen: f & 1 != 0, off: f & 2 != 0, locked: f & 4 != 0 });
    }
    for _ in 0..c.count("curve")? {
        let kind = c.next()?;
        let layer = c.count("layer index")?;
        let color = c.next()? as i32;
        let mut p = || -> Result<[f64; 2], KernelError> { Ok([c.next()?, c.next()?]) };
        let geometry = match kind as i32 {
            0 => Geometry::Line { a: p()?, b: p()? },
            1 => { let cc = p()?; Geometry::Arc { c: cc, r: c.next()?, start: c.next()?, sweep: c.next()? } }
            2 => { let cc = p()?; Geometry::Circle { c: cc, r: c.next()? } }
            3 => { let (cc, u, v) = (p()?, p()?, p()?); Geometry::Ellipse { c: cc, u, v, t0: c.next()?, t1: c.next()? } }
            4 => {
                let degree = c.count("degree")?;
                let knots = (0..c.count("knot")?).map(|_| c.next()).collect::<Result<Vec<_>, _>>()?;
                let n = c.count("point")?;
                let mut points = Vec::with_capacity(n);
                for _ in 0..n { points.push([c.next()?, c.next()?]); }
                let weights = (0..c.count("weight")?).map(|_| c.next()).collect::<Result<Vec<_>, _>>()?;
                Geometry::Spline { degree, knots, points, weights }
            }
            5 => {
                let closed = c.next()? != 0.0;
                let n = c.count("vertex")?;
                let (mut points, mut bulges) = (Vec::with_capacity(n), Vec::with_capacity(n));
                for _ in 0..n {
                    points.push([c.next()?, c.next()?]);
                    bulges.push(c.next()?);
                }
                Geometry::Polyline { points, bulges, closed }
            }
            k => return Err(KernelError::InvalidInput(format!("unknown DXF curve kind {k}"))),
        };
        drawing.curves.push(Curve { geometry, layer, color, group: 0 });
    }
    c.finish()?;
    if drawing.layers.is_empty() {
        drawing.layer_index("0");
    }
    let text = write_r12(&drawing, tol)?;
    let target = std::path::Path::new(&path);
    let temporary = target.with_extension("dxf.writing");
    std::fs::write(&temporary, text.as_bytes())
        .map_err(|e| KernelError::InvalidInput(format!("could not write '{}': {e}", temporary.display())))?;
    std::fs::rename(&temporary, target).map_err(|e| {
        let _ = std::fs::remove_file(&temporary);
        KernelError::InvalidInput(format!("could not replace '{path}': {e}"))
    })?;
    Ok(vec![text.len() as f64])
}

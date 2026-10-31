//! `Drawing` to an R12 (AC1009) ASCII DXF: HEADER (version, units, extents), TABLES (CONTINUOUS line type, the
//! layers), an empty BLOCKS section and the ENTITIES. Lines, arcs, circles and polylines are written exactly;
//! ellipses and splines, which R12 does not have, as tangent lines and arcs within the tolerance.

use std::f64::consts::TAU;
use std::fmt::Write as _;

use super::{fit, Drawing, Geometry, P2};
use crate::KernelError;

pub fn write_r12(d: &Drawing, tol: f64) -> Result<String, KernelError> {
    let mut pieces: Vec<(Geometry, usize, i32)> = Vec::new();
    for curve in &d.curves {
        if curve.layer >= d.layers.len() {
            return Err(KernelError::InvalidInput(format!("a curve names layer {} of {}", curve.layer, d.layers.len())));
        }
        check(&curve.geometry)?;
        for g in fit::lines_and_arcs(&curve.geometry, tol)? {
            pieces.push((g, curve.layer, curve.color));
        }
    }

    let mut out = String::new();
    let mut put = |code: i32, value: &str| {
        let _ = write!(out, "{code:>3}\r\n{value}\r\n");
    };

    put(0, "SECTION");
    put(2, "HEADER");
    put(9, "$ACADVER");
    put(1, "AC1009");
    if let Some(units) = d.insunits {
        put(9, "$INSUNITS");
        put(70, &units.to_string());
    }
    if let Some(m) = d.measurement {
        put(9, "$MEASUREMENT");
        put(70, &m.to_string());
    }
    let (lo, hi) = extents(pieces.iter().map(|p| &p.0));
    for (name, p) in [("$EXTMIN", lo), ("$EXTMAX", hi)] {
        put(9, name);
        put(10, &num(p[0]));
        put(20, &num(p[1]));
        put(30, "0.0");
    }
    put(0, "ENDSEC");

    put(0, "SECTION");
    put(2, "TABLES");
    put(0, "TABLE");
    put(2, "LTYPE");
    put(70, "1");
    put(0, "LTYPE");
    put(2, "CONTINUOUS");
    put(70, "0");
    put(3, "Solid line");
    put(72, "65");
    put(73, "0");
    put(40, "0.0");
    put(0, "ENDTAB");
    put(0, "TABLE");
    put(2, "LAYER");
    let has_zero = d.layers.iter().any(|l| l.name == "0");
    put(70, &(d.layers.len() + usize::from(!has_zero)).to_string());
    if !has_zero {
        layer(&mut put, "0", 7, 0);
    }
    for l in &d.layers {
        let color = l.color.clamp(1, 255);
        let flags = i32::from(l.frozen) | if l.locked { 4 } else { 0 };
        layer(&mut put, &name(&l.name), if l.off { -color } else { color }, flags);
    }
    put(0, "ENDTAB");
    put(0, "ENDSEC");

    put(0, "SECTION");
    put(2, "BLOCKS");
    put(0, "ENDSEC");

    put(0, "SECTION");
    put(2, "ENTITIES");
    for (g, layer, color) in &pieces {
        let layer = name(&d.layers[*layer].name);
        let common = |put: &mut dyn FnMut(i32, &str), kind: &str| {
            put(0, kind);
            put(8, &layer);
            if (1..=255).contains(color) {
                put(62, &color.to_string());
            }
        };
        match g {
            Geometry::Line { a, b } => {
                common(&mut put, "LINE");
                put(10, &num(a[0])); put(20, &num(a[1])); put(30, "0.0");
                put(11, &num(b[0])); put(21, &num(b[1])); put(31, "0.0");
            }
            Geometry::Circle { c, r } => {
                common(&mut put, "CIRCLE");
                put(10, &num(c[0])); put(20, &num(c[1])); put(30, "0.0");
                put(40, &num(*r));
            }
            Geometry::Arc { c, r, start, sweep } => {
                // DXF arcs run counter-clockwise from 50 to 51.
                let (from, to) = if *sweep > 0.0 { (*start, start + sweep) } else { (start + sweep, *start) };
                common(&mut put, "ARC");
                put(10, &num(c[0])); put(20, &num(c[1])); put(30, "0.0");
                put(40, &num(*r));
                put(50, &num(degrees(from)));
                put(51, &num(degrees(to)));
            }
            Geometry::Polyline { points, bulges, closed } => {
                common(&mut put, "POLYLINE");
                put(66, "1");
                put(10, "0.0"); put(20, "0.0"); put(30, "0.0");
                put(70, if *closed { "1" } else { "0" });
                for (i, p) in points.iter().enumerate() {
                    common(&mut put, "VERTEX");
                    put(10, &num(p[0])); put(20, &num(p[1])); put(30, "0.0");
                    let bulge = bulges.get(i).copied().unwrap_or(0.0);
                    if bulge != 0.0 {
                        put(42, &num(bulge));
                    }
                }
                common(&mut put, "SEQEND");
            }
            Geometry::Ellipse { .. } | Geometry::Spline { .. } => unreachable!("fitted above"),
        }
    }
    put(0, "ENDSEC");
    put(0, "EOF");
    Ok(out)
}

fn layer(put: &mut dyn FnMut(i32, &str), name: &str, color: i32, flags: i32) {
    put(0, "LAYER");
    put(2, name);
    put(70, &flags.to_string());
    put(62, &color.to_string());
    put(6, "CONTINUOUS");
}

/// R12 layer names: letters, digits, $ - _ only, at most 31 characters.
fn name(s: &str) -> String {
    let cleaned: String = s.chars()
        .map(|c| if c.is_ascii_alphanumeric() || "$-_".contains(c) { c.to_ascii_uppercase() } else { '_' })
        .take(31)
        .collect();
    if cleaned.is_empty() { "0".into() } else { cleaned }
}

/// Shortest text that reads back as the same number; never exponent notation, never "-0".
fn num(v: f64) -> String {
    let v = if v == 0.0 { 0.0 } else { v };
    let s = format!("{v}");
    if s.contains('.') { s } else { format!("{s}.0") }
}

fn degrees(radians: f64) -> f64 {
    let d = radians.rem_euclid(TAU).to_degrees();
    if d >= 360.0 { d - 360.0 } else { d }
}

fn check(g: &Geometry) -> Result<(), KernelError> {
    let finite = |ps: &[f64]| ps.iter().all(|v| v.is_finite());
    let ok = match g {
        Geometry::Line { a, b } => finite(&[a[0], a[1], b[0], b[1]]),
        Geometry::Arc { c, r, start, sweep } =>
            finite(&[c[0], c[1], *r, *start, *sweep]) && *r > 0.0 && *sweep != 0.0 && sweep.abs() < TAU,
        Geometry::Circle { c, r } => finite(&[c[0], c[1], *r]) && *r > 0.0,
        Geometry::Ellipse { c, u, v, t0, t1 } => finite(&[c[0], c[1], u[0], u[1], v[0], v[1], *t0, *t1]) && t1 > t0,
        Geometry::Spline { degree, knots, points, weights } =>
            *degree >= 1 && points.len() > *degree && knots.len() == points.len() + degree + 1
                && (weights.is_empty() || weights.len() == points.len())
                && finite(knots) && points.iter().all(|p| finite(p)) && weights.iter().all(|w| *w > 0.0),
        Geometry::Polyline { points, bulges, .. } =>
            points.len() >= 2 && (bulges.is_empty() || bulges.len() == points.len())
                && points.iter().all(|p| finite(p)) && finite(bulges),
    };
    if ok { Ok(()) } else { Err(KernelError::InvalidInput(format!("a curve to write is not valid: {g:?}"))) }
}

fn extents<'a>(pieces: impl Iterator<Item = &'a Geometry>) -> (P2, P2) {
    let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
    let mut add = |p: P2| {
        lo = [lo[0].min(p[0]), lo[1].min(p[1])];
        hi = [hi[0].max(p[0]), hi[1].max(p[1])];
    };
    for g in pieces {
        match g {
            Geometry::Line { a, b } => { add(*a); add(*b); }
            Geometry::Circle { c, r } => { add([c[0] - r, c[1] - r]); add([c[0] + r, c[1] + r]); }
            Geometry::Arc { c, r, start, sweep } => {
                let at = |t: f64| [c[0] + r * t.cos(), c[1] + r * t.sin()];
                add(at(*start));
                add(at(start + sweep));
                let (from, to) = if *sweep > 0.0 { (*start, start + sweep) } else { (start + sweep, *start) };
                for k in 0..8 {
                    let q = k as f64 * TAU / 4.0 - TAU;
                    if q > from && q < to { add(at(q)); }
                    if q + TAU > from && q + TAU < to { add(at(q + TAU)); }
                }
            }
            Geometry::Polyline { points, .. } => points.iter().for_each(|p| add(*p)),
            _ => {}
        }
    }
    if lo[0] > hi[0] { ([0.0; 2], [0.0; 2]) } else { (lo, hi) }
}

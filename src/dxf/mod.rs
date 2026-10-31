//! DXF drawings in and out (Bob 2026-10-31: "a robust DXF importer and possibly even an exporter in Rust", for the
//! flat plate nester). Written from Autodesk's public DXF Reference (group codes, section/table/entity layout, the
//! arbitrary axis algorithm, binary DXF); no third-party DXF source.
//!
//! Reading gives the drawing as seen from above (world XY): every model-space curve, blocks exploded, object
//! coordinate systems resolved, in drawing units. `fit` turns ellipses and splines into tangent lines and arcs within
//! a tolerance for consumers that take lines, arcs and circles only. Writing produces an R12 (AC1009) DXF, the
//! version every CAD and CAM program reads.

mod codes;
pub mod fit;
mod read;
pub mod svg;
pub mod wire;
mod write;

pub use read::{read, read_header};
pub use write::write_r12;

/// A point or vector in the drawing's XY plane.
pub type P2 = [f64; 2];

#[derive(Clone, Debug, PartialEq)]
pub enum Geometry {
    Line { a: P2, b: P2 },
    /// Circular arc: angles in radians, `sweep` signed (positive counter-clockwise), 0 < |sweep| < 2 pi.
    Arc { c: P2, r: f64, start: f64, sweep: f64 },
    /// Full circle, counter-clockwise from angle 0.
    Circle { c: P2, r: f64 },
    /// Elliptic arc c + u cos t + v sin t for t from t0 to t1 (t1 > t0; the whole ellipse when t1 - t0 is 2 pi).
    /// u and v are conjugate semi-diameters (perpendicular for an ELLIPSE entity, not after a skewing insert).
    Ellipse { c: P2, u: P2, v: P2, t0: f64, t1: f64 },
    /// B-spline (rational when `weights` is not empty) over [knots[degree], knots[len - degree - 1]].
    Spline { degree: usize, knots: Vec<f64>, points: Vec<P2>, weights: Vec<f64> },
    /// Vertices with the bulge (tan(sweep / 4), positive counter-clockwise) of the segment leaving each one. Written
    /// as one POLYLINE; the reader gives a polyline's pieces instead (sharing a `group`).
    Polyline { points: Vec<P2>, bulges: Vec<f64>, closed: bool },
}

#[derive(Clone, Debug, PartialEq)]
pub struct Curve {
    pub geometry: Geometry,
    /// Index into `Drawing::layers`.
    pub layer: usize,
    /// AutoCAD colour index as resolved for this curve: 1-255, or 256 BYLAYER.
    pub color: i32,
    /// Which top-level model-space entity the curve came from (a polyline's pieces or an insert's contents share it).
    pub group: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Layer {
    pub name: String,
    pub color: i32,
    pub frozen: bool,
    pub off: bool,
    pub locked: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Drawing {
    /// $ACADVER, e.g. "AC1009" or "AC1032"; empty when the file has no header.
    pub version: String,
    /// $INSUNITS (1 inch, 2 foot, 4 mm, 5 cm, 6 m, ...), when present.
    pub insunits: Option<i32>,
    /// $MEASUREMENT (0 imperial, 1 metric), when present.
    pub measurement: Option<i32>,
    pub layers: Vec<Layer>,
    pub curves: Vec<Curve>,
    /// Entities that are not curves (TEXT, DIMENSION, HATCH, ...) by DXF type name, with their counts.
    pub skipped: Vec<(String, usize)>,
    /// Things the reader noticed and handled, in words for the user (a missing block, curves off the XY plane, ...).
    pub warnings: Vec<String>,
}

impl Drawing {
    /// The index of the layer named `name` (DXF layer names ignore case), adding it when the table lacks it.
    pub fn layer_index(&mut self, name: &str) -> usize {
        let name = if name.trim().is_empty() { "0" } else { name.trim() };
        if let Some(i) = self.layers.iter().position(|l| l.name.eq_ignore_ascii_case(name)) {
            return i;
        }
        self.layers.push(Layer { name: name.to_string(), color: 7, frozen: false, off: false, locked: false });
        self.layers.len() - 1
    }

    pub(crate) fn skip(&mut self, kind: &str) {
        match self.skipped.iter_mut().find(|(k, _)| k == kind) {
            Some((_, n)) => *n += 1,
            None => self.skipped.push((kind.to_string(), 1)),
        }
    }

    pub(crate) fn warn(&mut self, message: String) {
        if !self.warnings.contains(&message) {
            self.warnings.push(message);
        }
    }
}

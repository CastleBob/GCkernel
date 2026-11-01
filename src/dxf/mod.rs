//! DXF drawings in and out (Bob 2026-10-31: "a robust DXF importer and possibly even an exporter in Rust", for the
//! flat plate nester). Written from Autodesk's public DXF Reference (group codes, section/table/entity layout, the
//! arbitrary axis algorithm, binary DXF); no third-party DXF source.
//!
//! Reading gives the drawing as seen from above (world XY): every model-space curve, blocks exploded, object
//! coordinate systems resolved, in drawing units. `fit` turns ellipses and splines into tangent lines and arcs within
//! a tolerance for consumers that take lines, arcs and circles only. Writing produces an R12 (AC1009) DXF, the
//! version every CAD and CAM program reads. `dwg` reads DWG files (R13 to 2018) into the same records as DXF; `pdf` reads PDF and PDF-compatible
//! Illustrator drawings.

mod codes;
pub mod dwg;
pub mod fit;
pub mod pdf;
mod read;
pub mod svg;
mod text;
pub mod wire;
mod write;

pub use read::{read, read_header};
pub use text::{mtext_plain, text_plain};
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

/// Where a text came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TextKind {
    /// TEXT: one line.
    Text,
    /// MTEXT: paragraphs, formatting codes stripped.
    MText,
    /// A visible ATTRIB of a block insert (its value).
    Attrib,
}

/// A text as seen from above (world XY), blocks and object coordinates resolved, in drawing units.
#[derive(Clone, Debug, PartialEq)]
pub struct Text {
    pub kind: TextKind,
    /// Plain text; lines (MTEXT paragraphs and line breaks) separated by a line feed. %%d, %%p and %%c are °, ± and Ø.
    pub text: String,
    /// The point the alignment places: the left end of the baseline for left/baseline text, the alignment point
    /// otherwise (TEXT 11 when 72/73 are not both 0), the attachment point of an MTEXT; the first point of aligned
    /// and fit text.
    pub at: P2,
    /// The second point of aligned (3) and fit (5) text; `at` otherwise.
    pub end: P2,
    /// Cap height measured square to the baseline.
    pub height: f64,
    /// Baseline direction, radians counter-clockwise from +X.
    pub rotation: f64,
    /// Relative character width (1 normal).
    pub width: f64,
    /// 0 left, 1 centre, 2 right, 3 aligned, 5 fit (TEXT 72; "middle" (4) is given as centre with `valign` 2).
    pub halign: i32,
    /// 0 baseline, 1 bottom, 2 middle, 3 top (TEXT 73, ATTRIB 74; MTEXT attachment rows are top, middle, bottom).
    pub valign: i32,
    /// MTEXT reference rectangle width (where lines wrap) in drawing units; 0 for no wrapping.
    pub box_width: f64,
    /// Text style name (code 7), "Standard" when the entity names none.
    pub style: String,
    /// The style's font: the TrueType family when the style gives one, else its font file ("txt", "romans.shx",
    /// "arial.ttf"); empty when the style is not in the file.
    pub font: String,
    /// The text was mirrored (a mirroring insert, a backward or upside-down TEXT, an object coordinate system seen
    /// from below): it is given readable instead, over the same place (baseline turned, left and right swapped).
    pub mirrored: bool,
    pub layer: usize,
    pub color: i32,
    pub group: usize,
}

/// A text style table entry.
#[derive(Clone, Debug, PartialEq)]
pub struct TextStyle {
    pub name: String,
    /// Primary font file (code 3).
    pub file: String,
    /// TrueType family name from the style's extended data (ACAD 1000), when present.
    pub family: String,
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
    /// Model-space texts: TEXT, MTEXT and visible insert ATTRIBs (DXF and DWG; not SVG or PDF).
    pub texts: Vec<Text>,
    /// The text style table.
    pub styles: Vec<TextStyle>,
    /// Entities that are not curves (TEXT, DIMENSION, HATCH, ...) by DXF type name, with their counts. TEXT and
    /// MTEXT are counted here as well as given in `texts`.
    pub skipped: Vec<(String, usize)>,
    /// Things the reader noticed and handled, in words for the user (a missing block, curves off the XY plane, ...).
    pub warnings: Vec<String>,
    /// The same about texts (mirrored, on edge, not decodable, ...), kept apart so that what a curves-only reader
    /// is told does not change.
    pub text_warnings: Vec<String>,
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

    /// The font a text of style `name` uses (family before file), empty when the style is not in the table.
    pub(crate) fn style_font(&self, name: &str) -> String {
        self.styles.iter().find(|s| s.name.eq_ignore_ascii_case(name.trim()))
            .map(|s| if s.family.trim().is_empty() { s.file.trim().to_string() } else { s.family.trim().to_string() })
            .unwrap_or_default()
    }

    pub(crate) fn warn_text(&mut self, message: String) {
        if !self.text_warnings.contains(&message) {
            self.text_warnings.push(message);
        }
    }

    pub(crate) fn warn(&mut self, message: String) {
        if !self.warnings.contains(&message) {
            self.warnings.push(message);
        }
    }
}

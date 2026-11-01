//! STEP import owns syntax-to-body conversion, never repair, modeling or document
//! publication. All solids must succeed before any result is returned. See README.
mod spline;
mod surface;
mod syntax;
pub mod diagnostic;
pub mod wire;
mod write;
pub use write::write;
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::{mass, neutral, topology, KernelError};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::{PI, TAU};
use std::io::Read;
use std::path::Path;
use syntax::{Entity, Exchange, Value};
type Result<T> = std::result::Result<T, KernelError>;
fn invalid<T>(s: impl Into<String>) -> Result<T> {
    Err(KernelError::InvalidInput(s.into()))
}
fn unsupported<T>(s: impl Into<String>) -> Result<T> {
    Err(KernelError::Unsupported(s.into()))
}
fn context(e: KernelError, id: u64) -> KernelError {
    match e {
        KernelError::InvalidInput(s) => KernelError::InvalidInput(format!("STEP #{id}: {s}")),
        KernelError::Unsupported(s) => KernelError::Unsupported(format!("STEP #{id}: {s}")),
        KernelError::Internal(s) => KernelError::Internal(format!("STEP #{id}: {s}")),
    }
}

#[derive(Clone, Debug)]
pub struct ImportedBody {
    /// Source identity only, not a persistent application topology name.
    pub entity_id: u64,
    pub name: String,
    pub source_length_to_inches: f64,
    pub body: Body,
    pub volume: f64,
}
#[derive(Clone, Debug)]
pub struct StepImport {
    pub bodies: Vec<ImportedBody>,
}

pub fn read_file(path: &Path, tolerance: f64) -> Result<StepImport> {
    let file = std::fs::File::open(path)
        .map_err(|e| KernelError::InvalidInput(format!("cannot read STEP: {e}")))?;
    let mut bytes = Vec::new();
    file.take(syntax::MAX_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| KernelError::InvalidInput(format!("cannot read STEP: {e}")))?;
    read(&bytes, tolerance)
}

/// First admission domain: explicit manifold solids in standalone shape
/// representations, analytic supports and exact line/circle/ellipse boundaries.
/// No assembly flattening, spline fitting, topology sewing or unit guessing.
pub fn read(bytes: &[u8], tolerance: f64) -> Result<StepImport> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return invalid("STEP tolerance must be finite and positive");
    }
    let exchange = syntax::parse(bytes)?;
    for schema in &exchange.schemas {
        let name = schema
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_ascii_uppercase();
        if !matches!(
            name.as_str(),
            "CONFIG_CONTROL_DESIGN"
                | "AUTOMOTIVE_DESIGN"
                | "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF"
        ) {
            return unsupported(format!("STEP schema {schema}"));
        }
    }
    let doc = Document {
        exchange,
        tolerance,
    };
    // These records can change placement, multiplicity or the material domain.
    // Refuse before selecting any body: no silently omitted assembly member.
    for (id, e) in &doc.exchange.entities {
        for kind in e.components.keys() {
            if kind.ends_with("SHAPE_REPRESENTATION")
                && !matches!(
                    kind.as_str(),
                    "SHAPE_REPRESENTATION" | "ADVANCED_BREP_SHAPE_REPRESENTATION"
                )
            {
                return unsupported(format!(
                    "STEP #{id} {kind} is outside the initial standalone-solid import domain"
                ));
            }
        }
        for forbidden in [
            "MAPPED_ITEM",
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE",
            "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION",
            "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION",
            "BREP_WITH_VOIDS",
            "FACETED_BREP",
            "SHELL_BASED_SURFACE_MODEL",
            "MANIFOLD_SURFACE_SHAPE_REPRESENTATION",
            "TESSELLATED_SHAPE_REPRESENTATION",
            "EXTERNALLY_DEFINED_REPRESENTATION",
            "CSG_SOLID",
            "SOLID_REPLICA",
        ] {
            if e.has(forbidden) {
                return unsupported(format!(
                    "STEP #{id} {forbidden} is outside the initial standalone-solid import domain"
                ));
            }
        }
    }
    let mut represented = BTreeMap::new();
    for (id, e) in &doc.exchange.entities {
        for kind in [
            "SHAPE_REPRESENTATION",
            "ADVANCED_BREP_SHAPE_REPRESENTATION",
            "MANIFOLD_SURFACE_SHAPE_REPRESENTATION",
        ] {
            if !e.has(kind) {
                continue;
            }
            let a = e.args(kind, 3).map_err(|err| context(err, *id))?;
            let (scale, angle) = doc.units(a[2].reference()?)?;
            for item in a[1].list()? {
                let item = item.reference()?;
                let entity = doc.entity(item)?;
                if entity.has("AXIS2_PLACEMENT_3D") {
                    doc.frame(item)?;
                    continue;
                }
                if !entity.has("MANIFOLD_SOLID_BREP") {
                    return unsupported(format!(
                        "STEP #{item} representation item {}",
                        entity.kind()
                    ));
                }
                if represented.insert(item, (scale, angle)).is_some() {
                    return unsupported(format!(
                        "STEP solid #{item} occurs in more than one representation"
                    ));
                }
            }
        }
    }
    let solids: Vec<u64> = doc
        .exchange
        .entities
        .iter()
        .filter(|(_, e)| e.has("MANIFOLD_SOLID_BREP"))
        .map(|(id, _)| *id)
        .collect();
    if solids.is_empty() {
        return invalid("STEP contains no supported closed solid");
    }
    let mut bodies = Vec::new();
    for id in solids {
        let Some(&(scale, angle)) = represented.get(&id) else {
            return invalid(format!(
                "STEP solid #{id} has no unambiguous representation/units"
            ));
        };
        // Borrow the same graph; unit context is per representation, not global.
        bodies.push(doc.body(id, scale, angle).map_err(|e| context(e, id))?);
    }
    Ok(StepImport { bodies })
}

struct Document {
    exchange: Exchange,
    tolerance: f64,
}
impl Document {
    fn entity(&self, id: u64) -> Result<&Entity> {
        self.exchange
            .entities
            .get(&id)
            .ok_or_else(|| KernelError::InvalidInput(format!("missing STEP #{id}")))
    }
    fn args(&self, id: u64, name: &str, n: usize) -> Result<&[Value]> {
        self.entity(id)?.args(name, n).map_err(|e| context(e, id))
    }
    fn triple(value: &Value) -> Result<Vec3> {
        let a = value.list()?;
        if a.len() != 3 {
            return invalid("expected three coordinates");
        }
        Ok(Vec3::new(a[0].number()?, a[1].number()?, a[2].number()?))
    }
    fn point(&self, id: u64) -> Result<Vec3> {
        Self::triple(&self.args(id, "CARTESIAN_POINT", 2)?[1])
    }
    fn direction(&self, id: u64) -> Result<Vec3> {
        Self::triple(&self.args(id, "DIRECTION", 2)?[1])?.unit()
    }
    fn vector(&self, id: u64) -> Result<Vec3> {
        let a = self.args(id, "VECTOR", 3)?;
        let n = a[2].number()?;
        if n <= 0.0 {
            return invalid(format!("STEP VECTOR #{id} needs positive magnitude"));
        }
        let v = self.direction(a[1].reference()?)?.scale(n);
        if !v.is_finite() {
            return invalid("vector overflow");
        }
        Ok(v)
    }
    fn frame(&self, id: u64) -> Result<(Vec3, Vec3, Vec3)> {
        let a = self.args(id, "AXIS2_PLACEMENT_3D", 4)?;
        let p = self.point(a[1].reference()?)?;
        let n = if matches!(a[2], Value::Omitted) {
            Vec3::new(0.0, 0.0, 1.0)
        } else {
            self.direction(a[2].reference()?)?
        };
        // ISO build_axes defaults and projection define the frame; not error recovery.
        let r = if matches!(a[3], Value::Omitted) {
            if n.cross(Vec3::new(1.0, 0.0, 0.0)).length() == 0.0 {
                Vec3::new(0.0, 1.0, 0.0)
            } else {
                Vec3::new(1.0, 0.0, 0.0)
            }
        } else {
            self.direction(a[3].reference()?)?
        };
        let x = r
            .sub(n.scale(r.dot(n)))
            .unit()
            .map_err(|e| context(e, id))?;
        Ok((p, n, x))
    }
    fn units(&self, context_id: u64) -> Result<(f64, f64)> {
        let e = self.entity(context_id)?;
        if e.args("GEOMETRIC_REPRESENTATION_CONTEXT", 1)?[0].number()? != 3.0 {
            return unsupported("STEP representation must be 3D");
        }
        let units = e.args("GLOBAL_UNIT_ASSIGNED_CONTEXT", 1)?[0].list()?;
        let (mut length, mut angle) = (None, None);
        for v in units {
            let id = v.reference()?;
            let unit = self.entity(id)?;
            if unit.has("LENGTH_UNIT") {
                if length.is_some() {
                    return invalid("multiple STEP length units in one context");
                }
                length = Some(self.unit_factor(id, "LENGTH_UNIT", &mut BTreeSet::new())? / 0.0254);
            } else if unit.has("PLANE_ANGLE_UNIT") {
                if angle.is_some() {
                    return invalid("multiple STEP angle units in one context");
                }
                angle = Some(self.unit_factor(id, "PLANE_ANGLE_UNIT", &mut BTreeSet::new())?);
            }
        }
        let (Some(l), Some(a)) = (length, angle) else {
            return invalid(format!(
                "STEP context #{context_id} lacks declared length/angle units"
            ));
        };
        if !l.is_finite() || l <= 0.0 || !a.is_finite() || a <= 0.0 {
            return invalid("invalid STEP unit factors");
        }
        Ok((l, a))
    }
    fn unit_factor(&self, id: u64, dimension: &str, visiting: &mut BTreeSet<u64>) -> Result<f64> {
        if visiting.len() >= 32 || !visiting.insert(id) {
            return invalid("cyclic or excessive STEP unit conversion chain");
        }
        let e = self.entity(id)?;
        e.args(dimension, 0)?;
        let result = if e.has("SI_UNIT") {
            let a = e.args("SI_UNIT", 2)?;
            let expected = if dimension == "LENGTH_UNIT" {
                "METRE"
            } else {
                "RADIAN"
            };
            if a[1].enumeration()? != expected {
                return invalid(format!("STEP {dimension} is not {expected}"));
            }
            let power = if matches!(a[0], Value::Omitted) {
                0
            } else {
                match a[0].enumeration()? {
                    "EXA" => 18,
                    "PETA" => 15,
                    "TERA" => 12,
                    "GIGA" => 9,
                    "MEGA" => 6,
                    "KILO" => 3,
                    "HECTO" => 2,
                    "DECA" => 1,
                    "DECI" => -1,
                    "CENTI" => -2,
                    "MILLI" => -3,
                    "MICRO" => -6,
                    "NANO" => -9,
                    "PICO" => -12,
                    "FEMTO" => -15,
                    "ATTO" => -18,
                    s => return unsupported(format!("STEP SI prefix {s}")),
                }
            };
            10_f64.powi(power)
        } else if e.has("CONVERSION_BASED_UNIT") {
            let a = e.args("CONVERSION_BASED_UNIT", 2)?;
            let measure_id = a[1].reference()?;
            let measure = self.entity(measure_id)?;
            let kind = if measure.has("LENGTH_MEASURE_WITH_UNIT") {
                "LENGTH_MEASURE_WITH_UNIT"
            } else if measure.has("PLANE_ANGLE_MEASURE_WITH_UNIT") {
                "PLANE_ANGLE_MEASURE_WITH_UNIT"
            } else {
                "MEASURE_WITH_UNIT"
            };
            let m = measure.args(kind, 2)?;
            let Value::Typed(tag, values) = &m[0] else {
                return invalid("unit conversion lacks typed measure");
            };
            let expected = if dimension == "LENGTH_UNIT" {
                "LENGTH_MEASURE"
            } else {
                "PLANE_ANGLE_MEASURE"
            };
            if tag != expected || values.len() != 1 {
                return invalid("unit conversion has wrong measure dimension");
            }
            values[0].number()? * self.unit_factor(m[1].reference()?, dimension, visiting)?
        } else {
            return unsupported(format!("STEP unit #{id}: {}", e.kind()));
        };
        visiting.remove(&id);
        if !result.is_finite() || result <= 0.0 {
            return invalid("nonpositive or overflowing unit conversion");
        }
        Ok(result)
    }
    fn conic(&self, id: u64, scale: f64) -> Result<(Vec3, Vec3, Vec3, f64, f64)> {
        let e = self.entity(id)?;
        let circle = e.has("CIRCLE");
        let a = e.args(
            if circle { "CIRCLE" } else { "ELLIPSE" },
            if circle { 3 } else { 4 },
        )?;
        let (c, n, x) = self.frame(a[1].reference()?)?;
        let r1 = a[2].number()? * scale;
        let r2 = if circle { r1 } else { a[3].number()? * scale };
        if !r1.is_finite() || !r2.is_finite() || r1 <= 0.0 || r2 <= 0.0 {
            return invalid("invalid STEP conic radii");
        }
        Ok((c.scale(scale), n, x, r1, r2))
    }
    fn edge_curve(
        &self,
        id: u64,
        from: Vec3,
        to: Vec3,
        closed: bool,
        sense: bool,
        scale: f64,
    ) -> Result<Curve> {
        let e = self.entity(id)?;
        if e.has("LINE") {
            let a = e.args("LINE", 3)?;
            let p = self.point(a[1].reference()?)?.scale(scale);
            let d = self.vector(a[2].reference()?)?.unit()?;
            if closed
                || from.distance(to) < 1e-14
                || from.sub(p).cross(d).length() > self.tolerance
                || to.sub(p).cross(d).length() > self.tolerance
            {
                return invalid(format!("STEP LINE #{id} disagrees with its vertices"));
            }
            if (to.sub(from).dot(d) > 0.0) != sense {
                return invalid(format!("STEP LINE #{id} disagrees with edge same_sense"));
            }
            return Ok(Curve::Line { a: from, b: to });
        }
        if e.has("CIRCLE") || e.has("ELLIPSE") {
            let (c, n, x, r1, r2) = self.conic(id, scale)?;
            let n = if sense { n } else { n.scale(-1.0) };
            let y = n.cross(x);
            let angle = |p: Vec3| {
                let w = p.sub(c);
                (w.dot(y) / r2).atan2(w.dot(x) / r1)
            };
            let t0 = angle(from);
            let sweep = if closed {
                TAU
            } else {
                (angle(to) - t0).rem_euclid(TAU)
            };
            if sweep == 0.0 {
                return invalid(format!("STEP conic #{id} has coincident distinct vertices"));
            }
            let curve = if e.has("CIRCLE") {
                Curve::Arc {
                    c,
                    n,
                    x,
                    r: r1,
                    t0,
                    t1: t0 + sweep,
                }
            } else {
                Curve::Conic {
                    c,
                    u: x.scale(r1),
                    v: y.scale(r2),
                    t0,
                    t1: t0 + sweep,
                }
            };
            if curve.start().distance(from) > self.tolerance
                || curve.end().distance(to) > self.tolerance
            {
                return invalid(format!("STEP conic #{id} disagrees with its vertices"));
            }
            return Ok(curve);
        }
        // Pcurves are not silently discarded: their reconciliation is a later domain.
        unsupported(format!("STEP curve #{id} {}", e.kind()))
    }
    fn surface(&self, id: u64, scale: f64, angle: f64) -> Result<(Surface, bool)> {
        let e = self.entity(id)?;
        if e.has("B_SPLINE_SURFACE_WITH_KNOTS") {
            return Ok((self.spline_surface(id, scale)?, true));
        }
        for kind in [
            "PLANE",
            "CYLINDRICAL_SURFACE",
            "CONICAL_SURFACE",
            "SPHERICAL_SURFACE",
            "TOROIDAL_SURFACE",
        ] {
            if !e.has(kind) {
                continue;
            }
            let count = match kind {
                "PLANE" => 2,
                "CONICAL_SURFACE" | "TOROIDAL_SURFACE" => 4,
                _ => 3,
            };
            let a = e.args(kind, count)?;
            let (p, n, x) = self.frame(a[1].reference()?)?;
            let p = p.scale(scale);
            let radius = if count > 2 {
                a[2].number()? * scale
            } else {
                0.0
            };
            let surface = match kind {
                "PLANE" => Surface::Plane { p, n, x },
                "CYLINDRICAL_SURFACE" => Surface::Cylinder {
                    p,
                    a: n,
                    x,
                    r: radius,
                },
                "SPHERICAL_SURFACE" => Surface::Sphere {
                    c: p,
                    a: n,
                    x,
                    r: radius,
                },
                "TOROIDAL_SURFACE" => Surface::Torus {
                    c: p,
                    a: n,
                    x,
                    big: radius,
                    r: a[3].number()? * scale,
                },
                _ => {
                    let semi_angle = a[3].number()? * angle;
                    if !semi_angle.is_finite() || semi_angle.abs() >= PI / 2.0 || semi_angle == 0.0
                    {
                        return invalid("invalid conical semi-angle");
                    }
                    Surface::Cone {
                        p,
                        a: n,
                        x,
                        r: radius,
                        k: semi_angle.tan(),
                    }
                }
            };
            return Ok((surface, true));
        }
        if e.has("SURFACE_OF_LINEAR_EXTRUSION") {
            let a = e.args("SURFACE_OF_LINEAR_EXTRUSION", 3)?;
            let curve_id = a[1].reference()?;
            let base = self.entity(curve_id)?;
            let axis = self.vector(a[2].reference()?)?.unit()?;
            if base.has("LINE") {
                let a = base.args("LINE", 3)?;
                let p = self.point(a[1].reference()?)?.scale(scale);
                let x = self.vector(a[2].reference()?)?.unit()?;
                let n = x.cross(axis).unit()?;
                return Ok((Surface::Plane { p, n, x }, true));
            }
            if base.has("CIRCLE") || base.has("ELLIPSE") {
                let (p, n, x, r1, r2) = self.conic(curve_id, scale)?;
                if n.cross(axis).length() > 1e-12 {
                    return unsupported(format!("STEP #{id} oblique conic extrusion"));
                }
                let same = n.dot(axis) > 0.0;
                return Ok((
                    if base.has("CIRCLE") {
                        Surface::Cylinder {
                            p,
                            a: axis,
                            x,
                            r: r1,
                        }
                    } else {
                        Surface::EllipticCylinder {
                            p,
                            a: axis,
                            u: x.scale(r1),
                            v: axis.cross(x).scale(r2),
                        }
                    },
                    same,
                ));
            }
        }
        if e.has("SURFACE_OF_REVOLUTION") {
            let a = e.args("SURFACE_OF_REVOLUTION", 3)?;
            let line_id = a[1].reference()?;
            let line = self.args(line_id, "LINE", 3)?;
            let q = self.point(line[1].reference()?)?.scale(scale);
            let d = self.vector(line[2].reference()?)?.unit()?;
            let axis = self.args(a[2].reference()?, "AXIS1_PLACEMENT", 3)?;
            let origin = self.point(axis[1].reference()?)?.scale(scale);
            let axis = if matches!(axis[2], Value::Omitted) {
                Vec3::new(0.0, 0.0, 1.0)
            } else {
                self.direction(axis[2].reference()?)?
            };
            let axial = d.dot(axis);
            let p = origin.add(axis.scale(q.sub(origin).dot(axis)));
            let radial = q.sub(p);
            let r = radial.length();
            if !r.is_finite() || r <= self.tolerance {
                return unsupported(format!("STEP #{id} line-revolution anchor on its axis"));
            }
            let x = radial.scale(1.0 / r);
            if d.dot(axis.cross(x)).abs() > 1e-12 {
                return unsupported(format!("STEP #{id} skew-line revolution"));
            }
            let dr = d.dot(x);
            if axial.abs() < 1e-14 {
                return Ok((
                    Surface::Plane {
                        p,
                        n: axis.scale(-dr.signum()),
                        x,
                    },
                    true,
                ));
            }
            // A straight meridian's radius is affine in height. This is an exact
            // change of analytic representation, not curve fitting or a retry.
            let support = if dr.abs() < 1e-14 {
                Surface::Cylinder { p, a: axis, x, r }
            } else {
                Surface::Cone {
                    p,
                    a: axis,
                    x,
                    r,
                    k: dr / axial,
                }
            };
            return Ok((support, axial > 0.0));
        }
        unsupported(format!("STEP surface #{id} {}", e.kind()))
    }
    fn body(&self, id: u64, scale: f64, angle: f64) -> Result<ImportedBody> {
        self.body_diagnostic(id,scale,angle,None)
    }
    fn body_diagnostic(&self, id:u64, scale:f64, angle:f64, report:Option<&mut Vec<diagnostic::Mismatch>>) -> Result<ImportedBody> {
        let a = self.args(id, "MANIFOLD_SOLID_BREP", 2)?;
        let name = a[0].text()?.to_owned();
        let shell_id = a[1].reference()?;
        let face_ids = self.args(shell_id, "CLOSED_SHELL", 2)?[1].list()?;
        if face_ids.is_empty() {
            return invalid("empty closed shell");
        }
        let mut body = Body::default();
        let mut vertices = BTreeMap::new();
        let mut edges: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        let mut seen_faces = BTreeSet::new();
        for v in face_ids {
            let face_id = v.reference()?;
            if !seen_faces.insert(face_id) {
                return invalid("duplicate face in closed shell");
            }
            let a = self.args(face_id, "ADVANCED_FACE", 4)?;
            let (surface, support_sense) = self.surface(a[2].reference()?, scale, angle)?;
            let mut face = Face {
                surface,
                sense: a[3].flag()? == support_sense,
                loops: Vec::new(),
            };
            let mut seen_bounds = BTreeSet::new();
            let mut outer = 0;
            for b in a[1].list()? {
                let bound_id = b.reference()?;
                if !seen_bounds.insert(bound_id) {
                    return invalid("duplicate bound in STEP face");
                }
                let bound = self.entity(bound_id)?;
                let is_outer = bound.has("FACE_OUTER_BOUND");
                if is_outer {
                    outer += 1;
                }
                let a = bound.args(
                    if is_outer {
                        "FACE_OUTER_BOUND"
                    } else {
                        "FACE_BOUND"
                    },
                    3,
                )?;
                let forward = a[2].flag()?;
                let loop_id = a[1].reference()?;
                let mut lp = Vec::new();
                for v in self.args(loop_id, "EDGE_LOOP", 2)?[1].list()? {
                    let a = self.args(v.reference()?, "ORIENTED_EDGE", 5)?;
                    if !matches!(a[1], Value::Derived) || !matches!(a[2], Value::Derived) {
                        return invalid("ORIENTED_EDGE endpoints must be derived");
                    }
                    let edge_id = a[3].reference()?;
                    let orientation = a[4].flag()?;
                    if !edges.contains_key(&edge_id) {
                        let a = self.args(edge_id, "EDGE_CURVE", 5)?;
                        let mut end = [0; 2];
                        for i in 0..2 {
                            let vertex_id = a[i + 1].reference()?;
                            end[i] = if let Some(&index) = vertices.get(&vertex_id) {
                                index
                            } else {
                                let v = self.args(vertex_id, "VERTEX_POINT", 2)?;
                                let p = self.point(v[1].reference()?)?.scale(scale);
                                let index = body.vertices.len();
                                body.vertices.push(p);
                                vertices.insert(vertex_id, index);
                                index
                            };
                        }
                        let support_id = a[3].reference()?;
                        let support = self.entity(support_id)?;
                        let curves = if support.has("B_SPLINE_CURVE_WITH_KNOTS")
                            || support.has("TRIMMED_CURVE")
                        {
                            self.spline_edge(
                                support_id,
                                body.vertices[end[0]],
                                body.vertices[end[1]],
                                a[4].flag()?,
                                scale,
                            )
                            .map_err(|e| context(e, support_id))?
                        } else {
                            vec![self.edge_curve(
                                a[3].reference()?,
                                body.vertices[end[0]],
                                body.vertices[end[1]],
                                end[0] == end[1],
                                a[4].flag()?,
                                scale,
                            )?]
                        };
                        if body.edges.len() + curves.len() > 100_000 {
                            return invalid("STEP body exceeds 100000 decomposed edges");
                        }
                        let mut indices = Vec::with_capacity(curves.len());
                        let count = curves.len();
                        let mut start = end[0];
                        for (i, curve) in curves.into_iter().enumerate() {
                            let finish = if i + 1 == count {
                                end[1]
                            } else {
                                let index = body.vertices.len();
                                body.vertices.push(curve.end());
                                index
                            };
                            indices.push(body.edges.len());
                            body.edges.push(Edge {
                                curve,
                                start,
                                end: finish,
                            });
                            start = finish;
                        }
                        edges.insert(edge_id, indices);
                    }
                    let indices = &edges[&edge_id];
                    if orientation {
                        lp.extend(indices.iter().map(|&i| (i, true)));
                    } else {
                        lp.extend(indices.iter().rev().map(|&i| (i, false)));
                    }
                }
                if !forward {
                    lp.reverse();
                    for (_, f) in &mut lp {
                        *f = !*f;
                    }
                }
                if is_outer {
                    face.loops.insert(0, lp);
                } else {
                    face.loops.push(lp);
                }
            }
            if outer != 1 {
                return unsupported(format!(
                    "STEP face #{face_id} requires exactly one explicit FACE_OUTER_BOUND"
                ));
            }
            body.faces.push(face);
        }
        // Existing strict neutral admission checks finite dimensions, radii, frames
        // and references before the topology validator accesses any indices.
        let data = neutral::write_body(&body);
        let mut reader = neutral::Reader::with_tolerance(&data, self.tolerance);
        let body = reader.body()?;
        reader.finish()?;
        if let Some(report)=report {
            report.extend(self.geometry_report(id,scale,&body,&edges,&face_ids.iter().map(Value::reference).collect::<Result<Vec<_>>>()?)?);
        }
        topology::check(&body, self.tolerance)?;
        // A CLOSED_SHELL is connected, even if disconnected shells share one file.
        let mut reached = BTreeSet::from([0]);
        let mut todo = vec![0];
        let mut owners: Vec<Vec<usize>> = vec![Vec::new(); body.edges.len()];
        for (f, face) in body.faces.iter().enumerate() {
            for &(e, _) in face.loops.iter().flatten() {
                owners[e].push(f);
            }
        }
        while let Some(f) = todo.pop() {
            for &(e, _) in body.faces[f].loops.iter().flatten() {
                for &other in &owners[e] {
                    if reached.insert(other) {
                        todo.push(other);
                    }
                }
            }
        }
        if reached.len() != body.faces.len() {
            return invalid("disconnected faces in one CLOSED_SHELL");
        }
        for (i, face) in body.faces.iter().enumerate() {
            if let Surface::Plane { n, .. } = face.surface {
                let outward = n.scale(if face.sense { 1.0 } else { -1.0 });
                for (j, lp) in face.loops.iter().enumerate() {
                    let curves: Vec<_> = lp
                        .iter()
                        .map(|&(e, f)| {
                            if f {
                                body.edges[e].curve.clone()
                            } else {
                                body.edges[e].curve.reversed()
                            }
                        })
                        .collect();
                    let area = crate::extrude::signed_area(&curves, outward);
                    if !area.is_finite() || (j == 0 && area <= 0.0) || (j != 0 && area >= 0.0) {
                        return invalid(format!("STEP face {i} bound {j} orientation disagrees with outward surface normal"));
                    }
                }
            } else {
                // Check the area sign on each curved face, not only the total
                // shell area: a reversed small wall can hide in a positive sum.
                // Compact this face to keep work proportional to its boundary.
                let mut local = Body::default();
                let mut vm = BTreeMap::new();
                let mut em = BTreeMap::new();
                let mut copy = face.clone();
                for lp in &mut copy.loops {
                    for (edge_id, _) in lp {
                        let source_id = *edge_id;
                        *edge_id = if let Some(&index) = em.get(&source_id) {
                            index
                        } else {
                            let mut edge = body.edges[source_id].clone();
                            for end in [&mut edge.start, &mut edge.end] {
                                let source = *end;
                                *end = *vm.entry(source).or_insert_with(|| {
                                    let index = local.vertices.len();
                                    local.vertices.push(body.vertices[source]);
                                    index
                                });
                            }
                            let index = local.edges.len();
                            local.edges.push(edge);
                            em.insert(source_id, index);
                            index
                        };
                    }
                }
                local.faces.push(copy);
                if mass::measure(&local, self.tolerance)?.area <= 0.0 {
                    return invalid(format!(
                        "STEP curved face {i} orientation disagrees with its boundary"
                    ));
                }
            }
        }
        let properties = mass::measure(&body, self.tolerance)?;
        if !properties.volume.is_finite()
            || properties.volume <= 0.0
            || !properties.area.is_finite()
            || properties.area <= 0.0
        {
            return invalid("STEP shell has nonpositive or non-finite volume/area");
        }
        Ok(ImportedBody {
            entity_id: id,
            name,
            source_length_to_inches: scale,
            body,
            volume: properties.volume,
        })
    }
}

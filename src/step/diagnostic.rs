//! Read-only diagnostics. An unadmitted body is never returned. The independent
//! source evaluator uses the original Cox-de Boor basis, not extracted Beziers.
use super::{invalid, syntax, Document, Result, Value};
use crate::{
    geom::{Surface, Vec3},
    neutral::Body,
    topology, KernelError,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Read,
    path::Path,
};
#[derive(Clone, Debug)]
pub struct Mismatch {
    pub solid: u64,
    pub face: u64,
    pub edge: u64,
    pub curve: u64,
    pub surface: u64,
    pub tolerance: f64,
    pub observed_departure: f64,
    pub source_observed_departure: Option<f64>,
    pub conversion_difference: Option<f64>,
    pub source_length_to_inches: f64,
}
pub struct Report {
    pub admitted: bool,
    pub admission_error: Option<String>,
    pub mismatches: Vec<Mismatch>,
    /// Geometry construction can fail before observations are available.
    pub body_errors: Vec<(u64, String)>,
    /// These reports currently inspect analytic supports only.
    pub analytic_supports_only: bool,
}
pub fn inspect_file(path: &Path, tolerance: f64) -> Result<Report> {
    let file = std::fs::File::open(path).map_err(|e| KernelError::InvalidInput(e.to_string()))?;
    let mut bytes = vec![];
    file.take(syntax::MAX_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| KernelError::InvalidInput(e.to_string()))?;
    inspect(&bytes, tolerance)
}
pub fn inspect(bytes: &[u8], tolerance: f64) -> Result<Report> {
    if !tolerance.is_finite() || tolerance <= 0. {
        return invalid("diagnostic tolerance must be positive");
    }
    let admission = super::read(bytes, tolerance);
    let doc = Document {
        exchange: syntax::parse(bytes)?,
        tolerance,
    };
    let mut represented = BTreeMap::new();
    for entity in doc.exchange.entities.values() {
        for kind in ["SHAPE_REPRESENTATION", "ADVANCED_BREP_SHAPE_REPRESENTATION"] {
            if !entity.has(kind) {
                continue;
            }
            let a = entity.args(kind, 3)?;
            let units = doc.units(a[2].reference()?)?;
            for item in a[1].list()? {
                let id = item.reference()?;
                if doc.entity(id)?.has("MANIFOLD_SOLID_BREP") {
                    if represented.insert(id, units).is_some() {
                        return invalid("diagnostic source solid has repeated representations");
                    }
                }
            }
        }
    }
    let mut mismatches = vec![];
    let mut body_errors = vec![];
    for (id, (scale, angle)) in represented {
        // Strict validation still runs and still fails. Only observations made
        // before that failure are retained, never a geometry result.
        if let Err(error) = doc.body_diagnostic(id, scale, angle, Some(&mut mismatches)) {
            body_errors.push((id, error.to_string()));
        }
    }
    Ok(Report {
        admitted: admission.is_ok(),
        admission_error: admission.err().map(|e| e.to_string()),
        mismatches,
        body_errors,
        analytic_supports_only: true,
    })
}
impl Document {
    pub(super) fn geometry_report(
        &self,
        solid: u64,
        scale: f64,
        body: &Body,
        edges: &BTreeMap<u64, Vec<usize>>,
        faces: &[u64],
    ) -> Result<Vec<Mismatch>> {
        let mut reverse = BTreeMap::new();
        for (&source, chain) in edges {
            for &e in chain {
                reverse.insert(e, source);
            }
        }
        let mut rows = vec![];
        let mut seen = BTreeSet::new();
        for (f, face) in body.faces.iter().enumerate() {
            if matches!(face.surface, Surface::Spline { .. }) {
                continue;
            }
            for &(e, _) in face.loops.iter().flatten() {
                let source = reverse[&e];
                if !seen.insert((f, source)) {
                    continue;
                }
                let chain = &edges[&source];
                if !chain.iter().any(|&e| {
                    topology::off_surface(&body.edges[e].curve, &face.surface, self.tolerance)
                        .is_some()
                }) {
                    continue;
                }
                let a = self.args(source, "EDGE_CURVE", 5)?;
                let curve = a[3].reference()?;
                let mut observed: f64 = 0.;
                for &e in chain {
                    let curve = &body.edges[e].curve;
                    let (a, b) = curve.domain();
                    for i in 0..=128 {
                        observed = observed.max(topology::distance(
                            curve.eval(a + (b - a) * i as f64 / 128.).0,
                            &face.surface,
                        ));
                    }
                }
                let oracle = self.source_spline_observations(
                    curve,
                    scale,
                    body,
                    chain,
                    a[4].flag()?,
                    &face.surface,
                )?;
                rows.push(Mismatch {
                    solid,
                    face: faces[f],
                    edge: source,
                    curve,
                    surface: self.args(faces[f], "ADVANCED_FACE", 4)?[2].reference()?,
                    tolerance: self.tolerance,
                    observed_departure: observed,
                    source_observed_departure: oracle.map(|p| p.0),
                    conversion_difference: oracle.map(|p| p.1),
                    source_length_to_inches: scale,
                });
            }
        }
        Ok(rows)
    }
    fn source_spline_observations(
        &self,
        id: u64,
        scale: f64,
        body: &Body,
        chain: &[usize],
        forward: bool,
        surface: &Surface,
    ) -> Result<Option<(f64, f64)>> {
        let e = self.entity(id)?;
        if !e.has("B_SPLINE_CURVE_WITH_KNOTS") {
            return Ok(None);
        }
        let (base, knots) = if e.has("B_SPLINE_CURVE") {
            (
                e.args("B_SPLINE_CURVE", 5)?,
                e.args("B_SPLINE_CURVE_WITH_KNOTS", 3)?,
            )
        } else {
            let a = e.args("B_SPLINE_CURVE_WITH_KNOTS", 9)?;
            (&a[1..6], &a[6..9])
        };
        let degree = base[0].number()? as usize;
        let p = base[1]
            .list()?
            .iter()
            .map(|v| self.point(v.reference()?).map(|p| p.scale(scale)))
            .collect::<Result<Vec<_>>>()?;
        let w = if e.has("RATIONAL_B_SPLINE_CURVE") {
            e.args("RATIONAL_B_SPLINE_CURVE", 1)?[0]
                .list()?
                .iter()
                .map(Value::number)
                .collect::<Result<Vec<_>>>()?
        } else {
            vec![1.; p.len()]
        };
        let distinct = knots[1]
            .list()?
            .iter()
            .map(Value::number)
            .collect::<Result<Vec<_>>>()?;
        let mult = knots[0]
            .list()?
            .iter()
            .map(Value::number)
            .collect::<Result<Vec<_>>>()?;
        if chain.len() + 1 != distinct.len() {
            return Ok(None);
        }
        // A bare spline can be trimmed by EDGE_CURVE vertices. Only compare
        // against the full original basis when its full endpoints are retained.
        let first = &body.edges[chain[0]].curve;
        let last = &body.edges[*chain.last().unwrap()].curve;
        let (a, b) = if forward {
            (first.start(), last.end())
        } else {
            (last.end(), first.start())
        };
        let round = 1024. * f64::EPSILON * p.iter().map(|q| q.length()).fold(1., f64::max);
        if a.distance(p[0]) > round || b.distance(*p.last().unwrap()) > round {
            return Ok(None);
        }
        let u: Vec<_> = distinct
            .iter()
            .zip(mult)
            .flat_map(|(&x, n)| std::iter::repeat_n(x, n as usize))
            .collect();
        let mut departure: f64 = 0.;
        let mut difference: f64 = 0.;
        for (span, k) in distinct.windows(2).enumerate() {
            for sample in 0..=128 {
                let t = sample as f64 / 128.;
                let at = k[0] + (k[1] - k[0]) * t;
                let point = if at == *u.last().unwrap() {
                    *p.last().unwrap()
                } else {
                    let mut basis: Vec<f64> = u
                        .windows(2)
                        .map(|k| if at >= k[0] && at < k[1] { 1. } else { 0. })
                        .collect();
                    for d in 1..=degree {
                        for i in 0..basis.len() - d {
                            let left = if u[i + d] > u[i] {
                                (at - u[i]) / (u[i + d] - u[i]) * basis[i]
                            } else {
                                0.
                            };
                            let right = if u[i + d + 1] > u[i + 1] {
                                (u[i + d + 1] - at) / (u[i + d + 1] - u[i + 1]) * basis[i + 1]
                            } else {
                                0.
                            };
                            basis[i] = left + right;
                        }
                    }
                    let mut q = Vec3::new(0., 0., 0.);
                    let mut weight = 0.;
                    for i in 0..p.len() {
                        let h = basis[i] * w[i];
                        q = q.add(p[i].sub(p[0]).scale(h));
                        weight += h;
                    }
                    p[0].add(q.scale(1. / weight))
                };
                let edge = if forward {
                    chain[span]
                } else {
                    chain[chain.len() - 1 - span]
                };
                let actual = body.edges[edge]
                    .curve
                    .eval(if forward { t } else { 1. - t })
                    .0;
                difference = difference.max(point.distance(actual));
                departure = departure.max(topology::distance(point, surface));
            }
        }
        Ok(Some((departure, difference)))
    }
}

//! Exact polynomial and positive-weight rational B-spline conversion. Homogeneous
//! knot insertion/degree elevation preserves shape, yielding cubic Bezier spans.
//! References: STEP EXPRESS b_spline_curve_with_knots; C.-K. Shene's public MTU
//! notes on B-spline knot insertion/subdivision. No third-party source code.
use super::{invalid, unsupported, Document, Result, Value};
use crate::geom::{Curve, Vec3};

impl Document {
    pub(super) fn spline_spans(&self, id: u64, scale: f64) -> Result<Vec<Curve>> {
        if self.entity(id)?.has("TRIMMED_CURVE") {
            return self.trimmed_spline(id, scale);
        }
        Ok(self.spline_data(id, scale)?.curves)
    }

    fn spline_data(&self, id: u64, scale: f64) -> Result<SplineData> {
        let e = self.entity(id)?;
        if !e.has("B_SPLINE_CURVE_WITH_KNOTS") {
            return unsupported(format!(
                "STEP trimmed basis #{id} is not a supported spline"
            ));
        }
        let rational = e.has("RATIONAL_B_SPLINE_CURVE");
        let (base, knots) = if e.has("B_SPLINE_CURVE") {
            e.args("REPRESENTATION_ITEM", 1)?[0].text()?;
            (
                e.args("B_SPLINE_CURVE", 5)?,
                e.args("B_SPLINE_CURVE_WITH_KNOTS", 3)?,
            )
        } else {
            let a = e.args("B_SPLINE_CURVE_WITH_KNOTS", 9)?;
            a[0].text()?;
            (&a[1..6], &a[6..9])
        };
        let degree = integer(&base[0], 1, 25)?;
        if degree > 3 {
            return unsupported(format!(
                "STEP spline #{id} degree {degree}; exact Bezier intake supports degrees 1..3"
            ));
        }
        match base[2].enumeration()? {
            "POLYLINE_FORM" | "CIRCULAR_ARC" | "ELLIPTIC_ARC" | "PARABOLIC_ARC"
            | "HYPERBOLIC_ARC" | "UNSPECIFIED" => (),
            _ => return invalid("invalid spline curve_form"),
        }
        for v in [&base[3], &base[4]] {
            if !matches!(v.enumeration()?, "T" | "F" | "U") {
                return invalid("invalid spline logical attribute");
            }
        }
        if base[4].enumeration()? == "T" {
            return unsupported("self-intersecting STEP spline boundary");
        }
        if !matches!(
            knots[2].enumeration()?,
            "UNIFORM_KNOTS" | "QUASI_UNIFORM_KNOTS" | "PIECEWISE_BEZIER_KNOTS" | "UNSPECIFIED"
        ) {
            return invalid("invalid spline knot_spec");
        }
        let refs = base[1].list()?;
        if refs.len() <= degree || refs.len() > 8192 {
            return invalid("spline needs degree+1..8192 control points");
        }
        let mut points: Vec<Vec3> = refs
            .iter()
            .map(|v| self.point(v.reference()?).map(|p| p.scale(scale)))
            .collect::<Result<_>>()?;
        if points.iter().any(|p| !p.is_finite()) {
            return invalid("nonfinite scaled spline control point");
        }
        let weights = if rational {
            e.args("RATIONAL_B_SPLINE_CURVE", 1)?[0]
                .list()?
                .iter()
                .map(Value::number)
                .collect::<Result<Vec<_>>>()?
        } else {
            vec![1.0; points.len()]
        };
        if weights.len() != points.len() || weights.iter().any(|w| !w.is_finite() || *w <= 0.0) {
            return invalid("rational spline needs one positive finite weight per control point");
        }
        let max = weights.iter().copied().fold(0.0, f64::max);
        if weights.iter().any(|w| w / max < 1e-12) {
            return unsupported("rational spline weight ratio exceeds 1e12");
        }
        let mut points: Vec<Homogeneous> = points
            .drain(..)
            .zip(weights)
            .map(|(p, w)| Homogeneous {
                q: p.scale(w / max),
                w: w / max,
            })
            .collect();
        let multiplicities = knots[0].list()?;
        let values = knots[1].list()?;
        if values.len() < 2
            || values.len() != multiplicities.len()
            || values.len() > refs.len() + degree + 1
        {
            return invalid("invalid spline knot/multiplicity counts");
        }
        let distinct: Vec<f64> = values.iter().map(Value::number).collect::<Result<_>>()?;
        if distinct.windows(2).any(|v| v[0] >= v[1]) {
            return invalid("spline knots must strictly increase");
        }
        let mult: Vec<usize> = multiplicities
            .iter()
            .map(|v| integer(v, 1, degree + 1))
            .collect::<Result<_>>()?;
        if mult.iter().sum::<usize>() != points.len() + degree + 1 {
            return invalid("spline knot sum disagrees with degree/control count");
        }
        if mult[0] != degree + 1 || mult[mult.len() - 1] != degree + 1 {
            return unsupported("unclamped/periodic STEP spline knot vector");
        }
        if mult[1..mult.len() - 1].iter().any(|&n| n > degree) {
            return invalid("discontinuous STEP spline boundary");
        }
        // Normalize the parameter domain before insertion. This changes only the
        // parameter, never control-point positions or shape.
        let parameters = distinct.clone();
        let lo = distinct[0];
        let span = distinct[distinct.len() - 1] - lo;
        if !span.is_finite() || span <= 0.0 {
            return invalid("unusable spline knot range");
        }
        let distinct: Vec<_> = distinct.iter().map(|u| (u - lo) / span).collect();
        if distinct.windows(2).any(|v| v[0] >= v[1]) {
            return invalid("spline knot range loses numerical resolution");
        }
        let mut expanded: Vec<f64> = distinct
            .iter()
            .zip(&mult)
            .flat_map(|(&u, &n)| std::iter::repeat_n(u, n))
            .collect();
        for (&u, &s) in distinct[1..distinct.len() - 1]
            .iter()
            .zip(&mult[1..mult.len() - 1])
        {
            for _ in s..degree {
                let k = expanded.partition_point(|&v| v <= u) - 1;
                let multiplicity = expanded.iter().filter(|&&v| v == u).count();
                let mut next = Vec::with_capacity(points.len() + 1);
                for i in 0..=points.len() {
                    next.push(if i <= k - degree {
                        points[i]
                    } else if i > k - multiplicity {
                        points[i - 1]
                    } else {
                        let alpha = (u - expanded[i]) / (expanded[i + degree] - expanded[i]);
                        points[i - 1].scale(1.0 - alpha).add(points[i].scale(alpha))
                    });
                }
                points = next;
                expanded.insert(k + 1, u);
            }
        }
        let mut out = Vec::with_capacity(distinct.len() - 1);
        for i in 0..distinct.len() - 1 {
            let p = &points[i * degree..=i * degree + degree];
            let p = match degree {
                1 => [
                    p[0],
                    p[0].scale(2.0 / 3.0).add(p[1].scale(1.0 / 3.0)),
                    p[0].scale(1.0 / 3.0).add(p[1].scale(2.0 / 3.0)),
                    p[1],
                ],
                2 => [
                    p[0],
                    p[0].scale(1.0 / 3.0).add(p[1].scale(2.0 / 3.0)),
                    p[1].scale(2.0 / 3.0).add(p[2].scale(1.0 / 3.0)),
                    p[2],
                ],
                _ => [p[0], p[1], p[2], p[3]],
            };
            if p.iter()
                .any(|p| !p.q.is_finite() || !p.w.is_finite() || p.w <= 0.0)
            {
                return invalid("spline decomposition overflow");
            }
            let weights = p.map(|h| h.w);
            let p = p.map(|h| h.q.scale(1.0 / h.w));
            if rational {
                crate::rational::validate(p, weights)?;
                out.push(Curve::RationalBezier { p, w: weights });
            } else {
                out.push(Curve::Bezier { p });
            }
        }
        Ok(SplineData {
            curves: out,
            parameters,
            closed: base[3].enumeration()? == "T",
        })
    }

    fn trimmed_spline(&self, id: u64, scale: f64) -> Result<Vec<Curve>> {
        let a = self.args(id, "TRIMMED_CURVE", 6)?;
        a[0].text()?;
        let basis = self.spline_data(a[1].reference()?, scale)?;
        let preference = a[5].enumeration()?;
        match preference {
            "PARAMETER" | "UNSPECIFIED" => (),
            "CARTESIAN" => {
                return unsupported(
                    "Cartesian-preferred STEP spline trimming requires point inversion",
                )
            }
            _ => return invalid("invalid STEP trimming preference"),
        }
        let trim = |value: &Value| -> Result<(f64, Option<Vec3>)> {
            let values = value.list()?;
            if values.is_empty() || values.len() > 2 {
                return invalid("trim set needs one or two selectors");
            }
            let (mut parameter, mut point) = (None, None);
            for value in values {
                match value {
                    Value::Typed(kind, args) if kind == "PARAMETER_VALUE" && args.len() == 1 => {
                        let u = args[0].number()?;
                        if !u.is_finite() || parameter.replace(u).is_some() {
                            return invalid("duplicate or nonfinite trim parameter");
                        }
                    }
                    Value::Ref(id) => {
                        let q = self.point(*id)?.scale(scale);
                        if !q.is_finite() || point.replace(q).is_some() {
                            return invalid("duplicate or nonfinite trim point");
                        }
                    }
                    _ => return invalid("invalid STEP trim selector"),
                }
            }
            let Some(u) = parameter else {
                return unsupported("Cartesian-only STEP spline trimming requires point inversion");
            };
            if point.is_some() && preference == "UNSPECIFIED" {
                return unsupported(
                    "dual STEP trim selectors require explicit parameter preference",
                );
            }
            Ok((u, point))
        };
        let (first, p1) = trim(&a[2])?;
        let (last, p2) = trim(&a[3])?;
        let forward = a[4].flag()?;
        let (lo, hi) = (basis.parameters[0], *basis.parameters.last().unwrap());
        if first < lo || first > hi || last < lo || last > hi {
            return invalid("STEP trim parameter is outside the source knot domain");
        }
        if first == last {
            return unsupported("equal STEP trim parameters do not specify an unambiguous span");
        }
        // Work in the basis's increasing direction, then reverse the COMPLETE chain.
        let (start, end) = if forward {
            (first, last)
        } else {
            (last, first)
        };
        let mut curves = if start < end {
            basis.interval(start, end)?
        } else {
            if !basis.closed {
                return invalid("open STEP spline trim order disagrees with sense_agreement");
            }
            let gap = basis.curves[0]
                .start()
                .distance(basis.curves.last().unwrap().end());
            if !gap.is_finite() || gap > self.tolerance {
                return invalid("closed STEP spline trim has a nonclosing seam");
            }
            let mut curves = basis.interval(start, hi)?;
            curves.extend(basis.interval(lo, end)?);
            if curves.is_empty() {
                return unsupported("STEP seam endpoints specify a zero-length trim");
            }
            curves
        };
        if !forward {
            curves.reverse();
            curves = curves.into_iter().map(|c| c.reversed()).collect();
        }
        for (point, actual) in [(p1, curves[0].start()), (p2, curves.last().unwrap().end())] {
            if let Some(point) = point {
                let gap = actual.distance(point);
                if !gap.is_finite() || gap > self.tolerance {
                    return invalid(format!(
                        "STEP parameter/Cartesian trim mismatch {gap} inches exceeds {}",
                        self.tolerance
                    ));
                }
            }
        }
        Ok(curves)
    }
}

struct SplineData {
    curves: Vec<Curve>,
    parameters: Vec<f64>,
    closed: bool,
}
impl SplineData {
    fn interval(&self, start: f64, end: f64) -> Result<Vec<Curve>> {
        let mut out = Vec::new();
        for (curve, knots) in self.curves.iter().zip(self.parameters.windows(2)) {
            let (lo, hi) = (start.max(knots[0]), end.min(knots[1]));
            if lo >= hi {
                continue;
            }
            if lo == knots[0] && hi == knots[1] {
                out.push(curve.clone());
                continue;
            }
            let a = (lo - knots[0]) / (knots[1] - knots[0]);
            let b = (hi - knots[0]) / (knots[1] - knots[0]);
            if !a.is_finite() || !b.is_finite() || a < 0.0 || b > 1.0 || a >= b {
                return invalid("STEP trim loses parameter resolution");
            }
            out.push(match *curve {
                Curve::Bezier { p } => Curve::Bezier {
                    p: crate::geom::bezier_piece(p, a, b),
                },
                Curve::RationalBezier { p, w } => {
                    let (p, w) = crate::rational::piece(p, w, a, b);
                    crate::rational::validate(p, w)?;
                    Curve::RationalBezier { p, w }
                }
                _ => return invalid("STEP spline decomposition produced a non-spline curve"),
            });
        }
        Ok(out)
    }
}

#[derive(Clone, Copy)]
struct Homogeneous {
    q: Vec3,
    w: f64,
}
impl Homogeneous {
    fn scale(self, a: f64) -> Self {
        Self {
            q: self.q.scale(a),
            w: self.w * a,
        }
    }
    fn add(self, b: Self) -> Self {
        Self {
            q: self.q.add(b.q),
            w: self.w + b.w,
        }
    }
}

fn integer(v: &Value, min: usize, max: usize) -> Result<usize> {
    let n = v.number()?;
    if !n.is_finite() || n.fract() != 0.0 || n < min as f64 || n > max as f64 {
        return invalid("invalid spline degree/multiplicity");
    }
    Ok(n as usize)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::syntax;
    fn document(
        degree: usize,
        multiplicity: usize,
        complex: bool,
    ) -> (Document, Vec<Vec3>, Vec<f64>) {
        let count = degree + 2 + multiplicity;
        let points: Vec<_> = (0..count)
            .map(|i| Vec3::new(i as f64, (i * i % 7) as f64, (i % 3) as f64))
            .collect();
        let refs = (0..count)
            .map(|i| format!("#{}", i + 1))
            .collect::<Vec<_>>()
            .join(",");
        let mut data = String::new();
        for (i, p) in points.iter().enumerate() {
            data += &format!("#{}=CARTESIAN_POINT('',({},{},{}));", i + 1, p.x, p.y, p.z);
        }
        let base = format!("{degree},({refs}),.UNSPECIFIED.,.F.,.F.");
        let knots = format!(
            "({},{multiplicity},1,{}),(0.,0.3,0.8,1.),.UNSPECIFIED.",
            degree + 1,
            degree + 1
        );
        data += &if complex {
            format!("#100=(BOUNDED_CURVE()B_SPLINE_CURVE({base})B_SPLINE_CURVE_WITH_KNOTS({knots})CURVE()GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM(''));" )
        } else {
            format!("#100=B_SPLINE_CURVE_WITH_KNOTS('',{base},{knots});")
        };
        let file=format!("ISO-10303-21;HEADER;FILE_DESCRIPTION((''),'2;1');FILE_NAME('','',(''),(''),'','','');FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));ENDSEC;DATA;{data}ENDSEC;END-ISO-10303-21;");
        let expanded = std::iter::repeat_n(0.0, degree + 1)
            .chain(std::iter::repeat_n(0.3, multiplicity))
            .chain([0.8])
            .chain(std::iter::repeat_n(1.0, degree + 1))
            .collect();
        (
            Document {
                exchange: syntax::parse(file.as_bytes()).unwrap(),
                tolerance: 1e-6,
            },
            points,
            expanded,
        )
    }
    // Independent Cox-de Boor basis oracle, not knot insertion/decomposition.
    fn basis(i: usize, p: usize, u: f64, k: &[f64]) -> f64 {
        if p == 0 {
            return if k[i] <= u && u < k[i + 1] { 1.0 } else { 0.0 };
        }
        let left = k[i + p] - k[i];
        let right = k[i + p + 1] - k[i + 1];
        (if left > 0.0 {
            (u - k[i]) / left * basis(i, p - 1, u, k)
        } else {
            0.0
        }) + (if right > 0.0 {
            (k[i + p + 1] - u) / right * basis(i + 1, p - 1, u, k)
        } else {
            0.0
        })
    }
    #[test]
    fn exact_spans_match_independent_basis_for_degrees_and_knot_multiplicities() {
        for degree in 1..=3 {
            for multiplicity in 1..=degree {
                for complex in [false, true] {
                    let (doc, points, knots) = document(degree, multiplicity, complex);
                    let curves = doc.spline_spans(100, 2.5).unwrap();
                    assert_eq!(curves.len(), 3);
                    for (span, range) in [0.0, 0.3, 0.8, 1.0].windows(2).enumerate() {
                        for j in 0..=100 {
                            let t = j as f64 / 100.0;
                            let u = range[0] + (range[1] - range[0]) * t;
                            let expected = if u == 1.0 {
                                *points.last().unwrap()
                            } else {
                                points
                                    .iter()
                                    .enumerate()
                                    .fold(Vec3::new(0., 0., 0.), |s, (i, p)| {
                                        s.add(p.scale(basis(i, degree, u, &knots)))
                                    })
                            }
                            .scale(2.5);
                            assert!(
                                curves[span].eval(t).0.distance(expected) < 1e-11,
                                "degree={degree} mult={multiplicity} span={span} t={t}"
                            );
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn rational_spans_match_independent_weighted_basis() {
        for degree in 1..=3 {
            for multiplicity in 1..=degree {
                let (mut doc, points, knots) = document(degree, multiplicity, true);
                let weights: Vec<f64> = (0..points.len())
                    .map(|i| 0.25 + (i % 4) as f64 * 0.7)
                    .collect();
                doc.exchange
                    .entities
                    .get_mut(&100)
                    .unwrap()
                    .components
                    .insert(
                        "RATIONAL_B_SPLINE_CURVE".into(),
                        vec![Value::List(
                            weights.iter().map(|&w| Value::Number(w)).collect(),
                        )],
                    );
                let curves = doc.spline_spans(100, 25.4).unwrap();
                for (span, range) in [0., 0.3, 0.8, 1.].windows(2).enumerate() {
                    for j in 0..100 {
                        let t = j as f64 / 100.;
                        let u = range[0] + (range[1] - range[0]) * t;
                        let bw: Vec<_> = (0..points.len())
                            .map(|i| weights[i] * basis(i, degree, u, &knots))
                            .collect();
                        let den: f64 = bw.iter().sum();
                        let expected = points
                            .iter()
                            .zip(&bw)
                            .fold(Vec3::new(0., 0., 0.), |q, (p, w)| {
                                q.add(p.scale(w / den * 25.4))
                            });
                        assert!(curves[span].eval(t).0.distance(expected) < 1e-10);
                    }
                }
            }
        }
    }
    fn add_trim(doc: &mut Document, first: f64, last: f64, forward: bool) {
        let parameter = |u| {
            Value::List(vec![Value::Typed(
                "PARAMETER_VALUE".into(),
                vec![Value::Number(u)],
            )])
        };
        doc.exchange.entities.insert(
            101,
            syntax::Entity {
                components: std::collections::BTreeMap::from([(
                    "TRIMMED_CURVE".into(),
                    vec![
                        Value::Text("".into()),
                        Value::Ref(100),
                        parameter(first),
                        parameter(last),
                        Value::Enum(if forward { "T" } else { "F" }.into()),
                        Value::Enum("PARAMETER".into()),
                    ],
                )]),
            },
        );
    }
    #[test]
    fn interior_trims_match_weighted_basis_for_all_supported_degrees_and_multiplicities() {
        for degree in 1..=3 {
            for mult in 1..=degree {
                for rational in [false, true] {
                    for reverse in [false, true] {
                        let (mut doc, points, knots) = document(degree, mult, true);
                        let weights: Vec<_> = (0..points.len())
                            .map(|i| if rational { 0.5 + (i % 3) as f64 } else { 1. })
                            .collect();
                        if rational {
                            doc.exchange
                                .entities
                                .get_mut(&100)
                                .unwrap()
                                .components
                                .insert(
                                    "RATIONAL_B_SPLINE_CURVE".into(),
                                    vec![Value::List(
                                        weights.iter().map(|&w| Value::Number(w)).collect(),
                                    )],
                                );
                        }
                        // Non-normalized source parameters; geometric units must not scale these.
                        let source = &mut doc
                            .exchange
                            .entities
                            .get_mut(&100)
                            .unwrap()
                            .components
                            .get_mut("B_SPLINE_CURVE_WITH_KNOTS")
                            .unwrap()[1];
                        *source = Value::List(
                            [0., 0.3, 0.8, 1.]
                                .iter()
                                .map(|&u| Value::Number(-7. + 20. * u))
                                .collect(),
                        );
                        for (lo, hi) in [(0.1, 0.9), (0.3, 0.8), (0.31, 0.7), (0., 1.)] {
                            let (first, last) = if reverse { (hi, lo) } else { (lo, hi) };
                            add_trim(&mut doc, -7. + 20. * first, -7. + 20. * last, !reverse);
                            let curves = doc.spline_spans(101, 25.4).unwrap();
                            let mut ranges: Vec<_> = [0., 0.3, 0.8, 1.]
                                .windows(2)
                                .filter_map(|r| {
                                    let (a, b) = (lo.max(r[0]), hi.min(r[1]));
                                    (a < b).then_some((a, b))
                                })
                                .collect();
                            if reverse {
                                ranges.reverse();
                            }
                            assert_eq!(curves.len(), ranges.len());
                            for (c, (a, b)) in curves.iter().zip(ranges) {
                                for i in 1..100 {
                                    let t = i as f64 / 100.;
                                    let u = if reverse {
                                        b - (b - a) * t
                                    } else {
                                        a + (b - a) * t
                                    };
                                    let bw: Vec<_> = (0..points.len())
                                        .map(|i| weights[i] * basis(i, degree, u, &knots))
                                        .collect();
                                    let den: f64 = bw.iter().sum();
                                    let expected = points
                                        .iter()
                                        .zip(&bw)
                                        .fold(Vec3::new(0., 0., 0.), |q, (p, w)| {
                                            q.add(p.scale(w / den * 25.4))
                                        });
                                    assert!(
                                        c.eval(t).0.distance(expected) < 1e-10,
                                        "degree={degree} reverse={reverse} range={lo}..{hi}"
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn malformed_and_unimplemented_splines_are_refused() {
        for (index, value) in [
            (1, Value::Number(3.5)),
            (1, Value::Number(4.0)),
            (4, Value::Enum("BAD".into())),
            (5, Value::Enum("T".into())),
            (6, Value::List(vec![Value::Number(4.0); 4])),
            (7, Value::List(vec![Value::Number(0.0); 4])),
            (8, Value::Enum("BAD".into())),
        ] {
            let (mut changed, _, _) = document(3, 2, false);
            changed
                .exchange
                .entities
                .get_mut(&100)
                .unwrap()
                .components
                .get_mut("B_SPLINE_CURVE_WITH_KNOTS")
                .unwrap()[index] = value;
            assert!(changed.spline_spans(100, 1.0).is_err());
        }
        let (mut rational, _, _) = document(3, 2, false);
        rational
            .exchange
            .entities
            .get_mut(&100)
            .unwrap()
            .components
            .insert("RATIONAL_B_SPLINE_CURVE".into(), vec![]);
        assert!(rational
            .spline_spans(100, 1.0)
            .unwrap_err()
            .to_string()
            .contains("RATIONAL_B_SPLINE_CURVE"));
    }
}

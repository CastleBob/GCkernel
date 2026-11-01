//! Explicit STEP pcurves. Contradictory associated geometry is never ignored.
use super::{invalid, unsupported, Document, Result, Value};
use crate::{
    geom::{Curve, Surface, Vec3},
    neutral::Face,
};
use std::collections::BTreeMap;
impl Document {
    pub(super) fn check_associated_supports(
        &self,
        edges: &BTreeMap<u64, std::collections::BTreeSet<u64>>,
    ) -> Result<()> {
        for (&edge, supports) in edges {
            let curve = self.args(edge, "EDGE_CURVE", 5)?[3].reference()?;
            if let Some(a) = self.surface_curve_args(curve)? {
                for reference in a[2].list()? {
                    let id = reference.reference()?;
                    let e = self.entity(id)?;
                    let support = if e.has("PCURVE") {
                        e.args("PCURVE", 3)?[1].reference()?
                    } else {
                        id
                    };
                    if !supports.contains(&support) {
                        return invalid(format!("surface curve #{curve} associates a support not incident to edge #{edge}"));
                    }
                }
            }
        }
        Ok(())
    }
    pub(super) fn point2(&self, id: u64) -> Result<Vec3> {
        let a = self.args(id, "CARTESIAN_POINT", 2)?[1].list()?;
        if a.len() != 2 {
            return invalid("pcurve control point must have dimension two");
        }
        Ok(Vec3::new(a[0].number()?, a[1].number()?, 0.))
    }
    fn surface_curve_args(&self, id: u64) -> Result<Option<&[Value]>> {
        let e = self.entity(id)?;
        for kind in ["SURFACE_CURVE", "SEAM_CURVE", "INTERSECTION_CURVE"] {
            if e.has(kind) {
                return Ok(Some(e.args(kind, 4)?));
            }
        }
        Ok(None)
    }
    pub(super) fn curve_3d(&self, id: u64) -> Result<u64> {
        let Some(a) = self.surface_curve_args(id)? else {
            return Ok(id);
        };
        let refs = a[2].list()?;
        if refs.is_empty() || refs.len() > 2 {
            return invalid("surface curve needs one or two associated geometries");
        }
        if !matches!(a[3].enumeration()?, "CURVE_3D" | "PCURVE_S1" | "PCURVE_S2") {
            return invalid("invalid surface curve master representation");
        }
        let curve = a[1].reference()?;
        if self.surface_curve_args(curve)?.is_some() || self.entity(curve)?.has("PCURVE") {
            return invalid("surface curve requires an independent 3D curve");
        }
        for (i, v) in refs.iter().enumerate() {
            let e = self.entity(v.reference()?)?;
            if (a[3].enumeration()? == "PCURVE_S1" && i == 0
                || a[3].enumeration()? == "PCURVE_S2" && i == 1)
                && !e.has("PCURVE")
            {
                return invalid("preferred surface curve representation is not a pcurve");
            }
        }
        if a[3].enumeration()? == "PCURVE_S2" && refs.len() != 2 {
            return invalid("missing preferred second pcurve");
        }
        Ok(curve)
    }
    pub(super) fn attach_pcurves(
        &self,
        face: &mut Face,
        surface_id: u64,
        edges: &BTreeMap<u64, Vec<usize>>,
        _scale: f64,
    ) -> Result<()> {
        let reverse: BTreeMap<usize, u64> = edges
            .iter()
            .flat_map(|(&id, chain)| chain.iter().map(move |&e| (e, id)))
            .collect();
        let mut candidates: BTreeMap<usize, Vec<Curve>> = BTreeMap::new();
        let mut explicit = false;
        for &(index, _) in face.loops.iter().flatten() {
            if candidates.contains_key(&index) {
                continue;
            }
            let id = reverse[&index];
            let edge = self.args(id, "EDGE_CURVE", 5)?;
            let Some(a) = self.surface_curve_args(edge[3].reference()?)? else {
                continue;
            };
            let mut associated = false;
            for value in a[2].list()? {
                let reference = value.reference()?;
                let e = self.entity(reference)?;
                if !e.has("PCURVE") {
                    if reference == surface_id {
                        associated = true;
                    } else {
                        // It must be an actual supported surface, even if it is
                        // associated with the other face of this shared edge.
                        self.surface(reference, _scale, 1.)?;
                    }
                    continue;
                }
                let pc = e.args("PCURVE", 3)?;
                if pc[1].reference()? != surface_id {
                    continue;
                }
                associated = true;
                explicit = true;
                if !matches!(face.surface, Surface::Spline { .. }) {
                    return unsupported(
                        "explicit pcurves on analytic supports are not yet reconciled",
                    );
                }
                let representation =
                    self.args(pc[2].reference()?, "DEFINITIONAL_REPRESENTATION", 3)?;
                let items = representation[1].list()?;
                if items.len() != 1 {
                    return invalid("pcurve definitional representation must contain one curve");
                }
                let context = self.entity(representation[2].reference()?)?;
                if context.args("GEOMETRIC_REPRESENTATION_CONTEXT", 1)?[0].number()? != 2. {
                    return invalid("pcurve context must have dimension two");
                }
                let mut curves = self.pcurve_spans(items[0].reference()?)?;
                let e = self.entity(surface_id)?;
                let (base, k) = if e.has("B_SPLINE_SURFACE") {
                    (
                        e.args("B_SPLINE_SURFACE", 7)?,
                        e.args("B_SPLINE_SURFACE_WITH_KNOTS", 5)?,
                    )
                } else {
                    let a = e.args("B_SPLINE_SURFACE_WITH_KNOTS", 13)?;
                    (&a[1..8], &a[8..13])
                };
                let rows = base[2].list()?;
                let counts = [rows.len(), rows[0].list()?.len()];
                let mut origin = [0.; 2];
                let mut extent = [0.; 2];
                for axis in 0..2 {
                    let degrees = base[axis].number()? as usize;
                    let expanded = k[axis + 2]
                        .list()?
                        .iter()
                        .zip(k[axis].list()?)
                        .map(|(v, m)| Ok((v.number()?, m.number()? as usize)))
                        .collect::<Result<Vec<_>>>()?
                        .into_iter()
                        .flat_map(|(v, m)| std::iter::repeat_n(v, m))
                        .collect::<Vec<_>>();
                    origin[axis] = expanded[degrees];
                    extent[axis] = expanded[counts[axis]] - origin[axis];
                }
                for curve in &mut curves {
                    let (p, w) = crate::spline_curve::controls(curve)?;
                    *curve = Curve::RationalBezier {
                        p: p.map(|p| {
                            Vec3::new(
                                (p.x - origin[0]) / extent[0],
                                (p.y - origin[1]) / extent[1],
                                0.,
                            )
                        }),
                        w,
                    };
                }
                if !edge[4].flag()? {
                    curves.reverse();
                    curves = curves.into_iter().map(|c| c.reversed()).collect();
                }
                let chain = &edges[&id];
                if curves.len() != chain.len() {
                    return unsupported("pcurve and physical edge have different span partitions");
                }
                for (&e, c) in chain.iter().zip(curves) {
                    candidates.entry(e).or_default().push(c);
                }
            }
            if !associated {
                return invalid("surface curve omits its incident face support");
            }
        }
        if !explicit {
            return Ok(());
        }
        let mut loops = Vec::new();
        for lp in &face.loops {
            let choices = lp
                .iter()
                .map(|&(e, fwd)| {
                    let Some(curves) = candidates.get(&e) else {
                        return unsupported("explicit spline trim requires pcurves for every edge");
                    };
                    Ok(curves
                        .iter()
                        .map(|c| if fwd { c.clone() } else { c.reversed() })
                        .collect::<Vec<_>>())
                })
                .collect::<Result<Vec<_>>>()?;
            // Seam pcurves are distinct occurrences. Select a unique closed UV
            // chain, with a fixed combinatorial budget; never guess a seam side.
            let mut paths = vec![Vec::<Curve>::new()];
            for options in choices {
                let mut next = Vec::new();
                for path in paths {
                    for curve in &options {
                        if path
                            .last()
                            .is_some_and(|p| p.end().distance(curve.start()) > 1e-12)
                        {
                            continue;
                        }
                        let mut p = path.clone();
                        p.push(curve.clone());
                        next.push(p);
                        if next.len() > 1024 {
                            return unsupported(
                                "pcurve seam correspondence exceeds 1024 candidates",
                            );
                        }
                    }
                }
                paths = next;
            }
            paths.retain(|p| p.last().unwrap().end().distance(p[0].start()) <= 1e-12);
            if paths.len() != 1 {
                return unsupported("pcurve loop has missing or ambiguous seam correspondence");
            }
            loops.push(paths.remove(0));
        }
        if let Surface::Spline { data } = &mut face.surface {
            data.pcurves = Some(loops);
        }
        Ok(())
    }
}

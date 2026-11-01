//! Independently authored AP203 B-rep serialization. Original curves, shared
//! topology and spline control nets are written directly; no meshing or fitting.
use super::{invalid, unsupported, Result};
use crate::{
    geom::{Curve, Surface, Vec3},
    mass,
    neutral::{self, Body},
    spline_surface::SplineSurface,
    topology,
};

fn real(x: f64) -> String {
    format!("{x:.17E}")
}
fn refs(ids: &[usize]) -> String {
    format!(
        "({})",
        ids.iter()
            .map(|i| format!("#{i}"))
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn reals(values: &[f64]) -> String {
    format!(
        "({})",
        values
            .iter()
            .map(|&v| real(v))
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn flag(v: bool) -> &'static str {
    if v {
        ".T."
    } else {
        ".F."
    }
}
fn string(value: &str) -> String {
    let mut out = String::from("'");
    for c in value.chars() {
        match c {
            '\'' => out.push_str("''"),
            '\\' => out.push_str("\\\\"),
            ' '..='~' => out.push(c),
            _ => out.push_str(&format!("\\X4\\{:08X}\\X0\\", c as u32)),
        }
    }
    out.push('\'');
    out
}
struct Writer {
    lines: Vec<String>,
    bytes: usize,
}
impl Writer {
    fn add(&mut self, value: String) -> Result<usize> {
        let id = self.lines.len() + 1;
        let line = format!("#{id}={value};\n");
        self.bytes += line.len();
        if self.bytes > super::syntax::MAX_BYTES - 4096 || id > 1_000_000 {
            return unsupported("STEP export exceeds its byte/entity budget");
        }
        self.lines.push(line);
        Ok(id)
    }
    fn point(&mut self, p: Vec3) -> Result<usize> {
        self.add(format!("CARTESIAN_POINT('',{})", reals(&[p.x, p.y, p.z])))
    }
    fn direction(&mut self, p: Vec3) -> Result<usize> {
        self.add(format!("DIRECTION('',{})", reals(&[p.x, p.y, p.z])))
    }
    fn axis(&mut self, p: Vec3, n: Vec3, x: Vec3) -> Result<usize> {
        let (p, n, x) = (self.point(p)?, self.direction(n)?, self.direction(x)?);
        self.add(format!("AXIS2_PLACEMENT_3D('',#{p},#{n},#{x})"))
    }
    fn vector(&mut self, p: Vec3) -> Result<usize> {
        let n = self.direction(p.unit()?)?;
        self.add(format!("VECTOR('',#{n},{})", real(p.length())))
    }
    fn ellipse(&mut self, c: Vec3, u: Vec3, v: Vec3) -> Result<usize> {
        // A rotation in parameter space orthogonalizes conjugate diameters;
        // determinant +1 preserves traversal. This is an exact conic identity.
        let angle = 0.5 * (2. * u.dot(v)).atan2(u.dot(u) - v.dot(v));
        let a = u.scale(angle.cos()).add(v.scale(angle.sin()));
        let b = v.scale(angle.cos()).sub(u.scale(angle.sin()));
        let axis = self.axis(c, a.cross(b).unit()?, a.unit()?)?;
        self.add(format!(
            "ELLIPSE('',#{axis},{},{})",
            real(a.length()),
            real(b.length())
        ))
    }
    fn curve(&mut self, c: &Curve) -> Result<usize> {
        match *c {
            Curve::Line { a, b } => {
                let p = self.point(a)?;
                let v = self.vector(b.sub(a))?;
                self.add(format!("LINE('',#{p},#{v})"))
            }
            Curve::Arc { c, n, x, r, .. } => {
                let axis = self.axis(c, n, x)?;
                self.add(format!("CIRCLE('',#{axis},{})", real(r)))
            }
            Curve::Conic { c, u, v, .. } => self.ellipse(c, u, v),
            Curve::Bezier { p } => {
                let ids = p
                    .iter()
                    .map(|&p| self.point(p))
                    .collect::<Result<Vec<_>>>()?;
                self.add(format!("B_SPLINE_CURVE_WITH_KNOTS('',3,{},.UNSPECIFIED.,{},.U.,(4,4),(0.,1.),.PIECEWISE_BEZIER_KNOTS.)",refs(&ids),flag(p[0]==p[3])))
            }
            Curve::RationalBezier { p, w } => {
                let ids = p
                    .iter()
                    .map(|&p| self.point(p))
                    .collect::<Result<Vec<_>>>()?;
                self.add(format!("(BOUNDED_CURVE() B_SPLINE_CURVE(3,{},.UNSPECIFIED.,{},.U.) B_SPLINE_CURVE_WITH_KNOTS((4,4),(0.,1.),.PIECEWISE_BEZIER_KNOTS.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE({}) REPRESENTATION_ITEM(''))",refs(&ids),flag(p[0]==p[3]),reals(&w)))
            }
            Curve::CylinderPair { .. } => {
                unsupported("STEP export of cylinder-intersection curves is not implemented")
            }
        }
    }
    fn spline(&mut self, s: &SplineSurface) -> Result<usize> {
        s.validate()?;
        let (nu, nv) = (s.u.len() - 1, s.v.len() - 1);
        let (rows, cols) = (nu * 3 + 1, nv * 3 + 1);
        if rows * cols > 16384 {
            return unsupported("STEP export spline grid exceeds 16384 controls");
        }
        let mut grid = vec![None::<(Vec3, f64)>; rows * cols];
        for i in 0..nu {
            for j in 0..nv {
                let p = &s.patches[i * nv + j];
                let scale = grid[i * 3 * cols + j * 3]
                    .map(|(_, w)| w / p.w[0])
                    .unwrap_or(1. / p.w[0]);
                for x in 0..4 {
                    for y in 0..4 {
                        let k = (i * 3 + x) * cols + j * 3 + y;
                        let value = (p.p[x * 4 + y], p.w[x * 4 + y] * scale);
                        if let Some((q, w)) = grid[k] {
                            let round =
                                256. * f64::EPSILON * q.length().max(value.0.length()).max(1.);
                            if q.distance(value.0) > round
                                || (w - value.1).abs()
                                    > 256. * f64::EPSILON * w.abs().max(value.1.abs())
                            {
                                return unsupported(
                                    "STEP export requires matching homogeneous spline patch seams",
                                );
                            }
                        } else {
                            grid[k] = Some(value);
                        }
                    }
                }
            }
        }
        let grid = grid.into_iter().map(Option::unwrap).collect::<Vec<_>>();
        let mut points = vec![];
        let mut weights = vec![];
        for row in grid.chunks(cols) {
            points.push(refs(
                &row.iter()
                    .map(|(p, _)| self.point(*p))
                    .collect::<Result<Vec<_>>>()?,
            ));
            weights.push(reals(&row.iter().map(|(_, w)| *w).collect::<Vec<_>>()));
        }
        let mult = |n: usize| {
            format!(
                "({})",
                (0..=n)
                    .map(|i| if i == 0 || i == n { "4" } else { "3" })
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        self.add(format!("(BOUNDED_SURFACE() B_SPLINE_SURFACE(3,3,({}),.UNSPECIFIED.,.U.,.U.,.U.) B_SPLINE_SURFACE_WITH_KNOTS({},{},{},{},.PIECEWISE_BEZIER_KNOTS.) GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(({})) REPRESENTATION_ITEM('') SURFACE())",points.join(","),mult(nu),mult(nv),reals(&s.u),reals(&s.v),weights.join(",")))
    }
    fn surface(&mut self, s: &Surface) -> Result<usize> {
        let (kind, p, a, x, values) = match *s {
            Surface::Plane { p, n, x } => ("PLANE", p, n, x, vec![]),
            Surface::Cylinder { p, a, x, r } => ("CYLINDRICAL_SURFACE", p, a, x, vec![r]),
            Surface::Cone { p, a, x, r, k } => {
                if k == 0. {
                    ("CYLINDRICAL_SURFACE", p, a, x, vec![r])
                } else {
                    ("CONICAL_SURFACE", p, a, x, vec![r, k.atan()])
                }
            }
            Surface::Sphere { c, a, x, r } => ("SPHERICAL_SURFACE", c, a, x, vec![r]),
            Surface::Torus { c, a, x, big, r } => ("TOROIDAL_SURFACE", c, a, x, vec![big, r]),
            Surface::EllipticCylinder { p, a, u, v } => {
                let c = self.ellipse(p, u, v)?;
                let vector = self.vector(a)?;
                return self.add(format!("SURFACE_OF_LINEAR_EXTRUSION('',#{c},#{vector})"));
            }
            Surface::Spline { ref data } => return self.spline(data),
            Surface::Revolution { .. } => {
                return unsupported(
                    "STEP export of general revolved spline surfaces is not implemented",
                )
            }
        };
        let axis = self.axis(p, a, x)?;
        let tail = values
            .iter()
            .map(|&v| format!(",{}", real(v)))
            .collect::<String>();
        self.add(format!("{kind}('',#{axis}{tail})"))
    }
    fn body(&mut self, name: &str, body: &Body) -> Result<usize> {
        let surfaces = body
            .faces
            .iter()
            .map(|f| self.surface(&f.surface))
            .collect::<Result<Vec<_>>>()?;
        let mut associations: Vec<Vec<usize>> = vec![Vec::new(); body.edges.len()];
        let mut owners: Vec<Vec<usize>> = vec![Vec::new(); body.edges.len()];
        let mut explicit = vec![false; body.edges.len()];
        for (fi, face) in body.faces.iter().enumerate() {
            for (li, lp) in face.loops.iter().enumerate() {
                for (ei, &(edge, forward)) in lp.iter().enumerate() {
                    let pcurves = match &face.surface {
                        Surface::Spline { data } => data.pcurves.as_ref(),
                        _ => None,
                    };
                    let association = if let Some(uv) = pcurves {
                        explicit[edge] = true;
                        let c = if forward {
                            uv[li][ei].clone()
                        } else {
                            uv[li][ei].reversed()
                        };
                        let (p, w) = crate::spline_curve::controls(&c)?;
                        let points = p
                            .iter()
                            .map(|p| {
                                self.add(format!("CARTESIAN_POINT('',{})", reals(&[p.x, p.y])))
                            })
                            .collect::<Result<Vec<_>>>()?;
                        let curve=self.add(format!("(BOUNDED_CURVE() B_SPLINE_CURVE(3,{},.UNSPECIFIED.,.F.,.U.) B_SPLINE_CURVE_WITH_KNOTS((4,4),(0.,1.),.PIECEWISE_BEZIER_KNOTS.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE({}) REPRESENTATION_ITEM(''))",refs(&points),reals(&w)))?;
                        let context = self.add(
                            "(GEOMETRIC_REPRESENTATION_CONTEXT(2) REPRESENTATION_CONTEXT('',''))"
                                .into(),
                        )?;
                        let rep = self.add(format!(
                            "DEFINITIONAL_REPRESENTATION('',(#{curve}),#{context})"
                        ))?;
                        self.add(format!("PCURVE('',#{},#{rep})", surfaces[fi]))?
                    } else {
                        surfaces[fi]
                    };
                    associations[edge].push(association);
                    owners[edge].push(fi);
                }
            }
        }
        let mut vertices = vec![];
        for &p in &body.vertices {
            let p = self.point(p)?;
            vertices.push(self.add(format!("VERTEX_POINT('',#{p})"))?);
        }
        let mut edges = vec![];
        for (i, e) in body.edges.iter().enumerate() {
            let mut c = self.curve(&e.curve)?;
            if explicit[i] {
                if associations[i].len() != 2 {
                    return invalid("pcurve edge must have two incident occurrences");
                }
                let kind = if owners[i][0] == owners[i][1] {
                    "SEAM_CURVE"
                } else {
                    "SURFACE_CURVE"
                };
                c = self.add(format!(
                    "{kind}('',#{c},{},.CURVE_3D.)",
                    refs(&associations[i])
                ))?;
            }
            edges.push(self.add(format!(
                "EDGE_CURVE('',#{},#{},#{c},.T.)",
                vertices[e.start], vertices[e.end]
            ))?);
        }
        let mut faces = vec![];
        for (fi, f) in body.faces.iter().enumerate() {
            let mut bounds = vec![];
            for (i, lp) in f.loops.iter().enumerate() {
                let mut oriented = vec![];
                for &(e, fwd) in lp {
                    oriented.push(self.add(format!(
                        "ORIENTED_EDGE('',*,*,#{},{})",
                        edges[e],
                        flag(fwd)
                    ))?);
                }
                let loop_id = self.add(format!("EDGE_LOOP('',{})", refs(&oriented)))?;
                bounds.push(self.add(format!(
                    "{}('',#{loop_id},.T.)",
                    if i == 0 {
                        "FACE_OUTER_BOUND"
                    } else {
                        "FACE_BOUND"
                    }
                ))?);
            }
            let surface = surfaces[fi];
            faces.push(self.add(format!(
                "ADVANCED_FACE('',{},#{surface},{})",
                refs(&bounds),
                flag(f.sense)
            ))?);
        }
        let shell = self.add(format!("CLOSED_SHELL('',{})", refs(&faces)))?;
        self.add(format!("MANIFOLD_SOLID_BREP({},#{shell})", string(name)))
    }
}

/// Return a complete deterministic AP203 file in inches. All bodies are checked
/// before serialization; errors return no partial file and perform no I/O.
pub fn write(bodies: &[(&str, &Body)], tolerance: f64) -> Result<String> {
    if !tolerance.is_finite() || tolerance <= 0. || bodies.is_empty() || bodies.len() > 1024 {
        return invalid("invalid STEP export tolerance/body count");
    }
    let mut checked = vec![];
    for &(name, b) in bodies {
        if name.len() > 32768
            || b.vertices.len() > 100000
            || b.edges.len() > 100000
            || b.faces.len() > 100000
        {
            return invalid("STEP export input exceeds its object budget");
        }
        let data = neutral::write_body(b);
        let mut reader = neutral::Reader::with_tolerance(&data, tolerance);
        let b = reader.body()?;
        reader.finish()?;
        topology::check(&b, tolerance)?;
        mass::measure(&b, tolerance)?;
        checked.push((name, b));
    }
    let mut w = Writer {
        lines: vec![],
        bytes: 0,
    };
    let mut solids = vec![];
    for (name, b) in &checked {
        solids.push(w.body(name, b)?);
    }
    let metre = w.add("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.))".into())?;
    let dims = w.add("DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.)".into())?;
    let factor = w.add(format!(
        "MEASURE_WITH_UNIT(LENGTH_MEASURE(0.0254),#{metre})"
    ))?;
    let inch = w.add(format!(
        "(CONVERSION_BASED_UNIT('inch',#{factor}) LENGTH_UNIT() NAMED_UNIT(#{dims}))"
    ))?;
    let radian = w.add("(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.))".into())?;
    let steradian = w.add("(NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT())".into())?;
    let uncertainty = w.add(format!(
        "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({}),#{inch},'distance_accuracy_value','')",
        real(tolerance)
    ))?;
    let context=w.add(format!("(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{uncertainty})) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{inch},#{radian},#{steradian})) REPRESENTATION_CONTEXT('',''))"))?;
    let rep = w.add(format!(
        "ADVANCED_BREP_SHAPE_REPRESENTATION('',{},#{context})",
        refs(&solids)
    ))?;
    let app=w.add("APPLICATION_CONTEXT('configuration controlled 3d designs of mechanical parts and assemblies')".into())?;
    w.add(format!("APPLICATION_PROTOCOL_DEFINITION('international standard','config_control_design',1994,#{app})"))?;
    let pc = w.add(format!("PRODUCT_CONTEXT('',#{app},'mechanical')"))?;
    let product = w.add(format!(
        "PRODUCT('GC Kernel export','GC Kernel export','',(#{pc}))"
    ))?;
    let formation = w.add(format!(
        "PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE('','',#{product},.NOT_KNOWN.)"
    ))?;
    let dc = w.add(format!(
        "PRODUCT_DEFINITION_CONTEXT('part definition',#{app},'design')"
    ))?;
    let definition = w.add(format!(
        "PRODUCT_DEFINITION('design','',#{formation},#{dc})"
    ))?;
    let shape = w.add(format!("PRODUCT_DEFINITION_SHAPE('','',#{definition})"))?;
    w.add(format!("SHAPE_DEFINITION_REPRESENTATION(#{shape},#{rep})"))?;
    let text = format!("ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('GC Kernel exact B-rep export'),'2;1');\nFILE_NAME('','',(''),(''),'GC Kernel','GC Kernel','');\nFILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));\nENDSEC;\nDATA;\n{}ENDSEC;\nEND-ISO-10303-21;\n",w.lines.concat());
    // Validate shell connectedness, individual bound orientation and the emitted
    // representation with the strict intake owner before returning any bytes.
    // This is an output check, not an independent correctness oracle.
    super::read(text.as_bytes(), tolerance)?;
    Ok(text)
}

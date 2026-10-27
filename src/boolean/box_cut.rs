//! Planar cutters on curved stock. Proven boxes with one globally through axis
//! use passage/intersection admission; other planar prisms and end slabs use a
//! compatible transverse cell chart selected before construction. No retries.
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Profile};
use crate::KernelError;
type Result<T> = std::result::Result<T, KernelError>;
fn unsupported(message: &str) -> KernelError {
    KernelError::Unsupported(message.into())
}
fn xyz(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    if !tol.is_finite() || tol <= 0. {
        return Err(KernelError::InvalidInput(
            "rectangular passage tolerance must be positive".into(),
        ));
    }
    // Oblique planar prisms enter the analytic section domain directly. No box
    // admission is attempted and failed operations never select another chart.
    if tool.vertices.len()!=8 || tool.edges.len()!=12 || tool.faces.len()!=6 || tool.faces.iter().any(|f|matches!(f.surface,Surface::Plane {n,..} if xyz(n).iter().filter(|v|**v!=0.).count()!=1)) {
        return super::endcut::planar_cross_cut(body,tool,tol);
    }
    let (lo, hi) = super::rectilinear::plain_box_bounds(tool, tol)?;
    if body.faces.len() > 5000 || body.edges.len() > 20000 || body.vertices.len() > 20000 {
        return Err(unsupported(
            "rectangular passage exceeds its input work limit",
        ));
    }
    let mut base_lo = [f64::INFINITY; 3];
    let mut base_hi = [f64::NEG_INFINITY; 3];
    for edge in &body.edges {
        let (l, h) = if let Curve::CylinderPair {
            c, a, rho, p, b, r, ..
        } = edge.curve
        {
            let a = xyz(a);
            let b = xyz(b);
            let c = xyz(c);
            let p = xyz(p);
            let nonzero: Vec<_> = (0..3).filter(|k| a[*k] != 0.).collect();
            if nonzero.len() != 1 || a[nonzero[0]].abs() != 1. || b[nonzero[0]] != 0. {
                return Err(unsupported(
                    "rectangular passage needs cardinal cylinder-pair axes",
                ));
            }
            let axis = nonzero[0];
            (
                std::array::from_fn(|k| if k == axis { p[k] - r } else { c[k] - rho }),
                std::array::from_fn(|k| if k == axis { p[k] + r } else { c[k] + rho }),
            )
        } else {
            super::passage::curve_bounds(&edge.curve)
        };
        for k in 0..3 {
            if !l[k].is_finite() || !h[k].is_finite() {
                return Err(KernelError::InvalidInput(
                    "nonfinite rectangular passage bounds".into(),
                ));
            }
            base_lo[k] = base_lo[k].min(l[k]);
            base_hi[k] = base_hi[k].max(h[k]);
        }
    }
    let axes: Vec<_> = (0..3)
        .filter(|k| lo[*k] <= base_lo[*k] && hi[*k] >= base_hi[*k])
        .collect();
    if axes.len() >= 2 {
        return super::endcut::planar_cross_cut(body, tool, tol);
    }
    if axes.is_empty() {
        return Err(unsupported(
            "rectangular tool must have one unambiguous axis spanning the complete stock bounds",
        ));
    }
    let axis = axes[0];
    let u = (axis + 1) % 3;
    let v = (axis + 2) % 3;
    let center = vector(std::array::from_fn(|k| (lo[k] + hi[k]) * 0.5));
    let radius = (hi[u] - lo[u]).hypot(hi[v] - lo[v]) * 0.5;
    if super::endcut::crosses_cylinder(body, axis, center, radius, lo[axis], hi[axis], tol) {
        return super::endcut::box_cross_cut(body, tool, axis, tol);
    }
    let passages =
        super::passage::round_passages(body, axis, center, radius, lo[axis], hi[axis], tol)?;
    let mut result = body.clone();
    for (bottom, top, bottom_face, top_face) in passages {
        let point = |s: f64, t: f64| {
            let mut q = [0.; 3];
            q[u] = s;
            q[v] = t;
            q[axis] = bottom;
            vector(q)
        };
        let corners = [
            point(lo[u], lo[v]),
            point(hi[u], lo[v]),
            point(hi[u], hi[v]),
            point(lo[u], hi[v]),
        ];
        let boundary = (0..4)
            .map(|i| Curve::Line {
                a: corners[i],
                b: corners[(i + 1) % 4],
            })
            .collect();
        let mut n = [0.; 3];
        n[axis] = 1.;
        let cavity = crate::extrude::extrude(
            &Profile {
                loops: vec![boundary],
                n: vector(n),
                low: 0.,
                high: top - bottom,
            },
            tol,
        )?;
        result = super::passage::splice_cavity(&result, cavity, axis, bottom_face, top_face, tol)?;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Surface;
    fn captured(name: &str) -> Body {
        let data = std::fs::read(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join(format!("tests/fixtures/rat_hole/{name}.f64")),
        )
        .unwrap();
        let values: Vec<_> = data
            .chunks_exact(8)
            .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
            .collect();
        crate::neutral::Reader::new(&values).body().unwrap()
    }
    fn tool(x0: f64, x1: f64, y0: f64, y1: f64, z0: f64, z1: f64) -> Body {
        let q = [
            Vec3::new(x0, y0, 0.),
            Vec3::new(x1, y0, 0.),
            Vec3::new(x1, y1, 0.),
            Vec3::new(x0, y1, 0.),
        ];
        crate::extrude::extrude(
            &Profile {
                loops: vec![(0..4)
                    .map(|i| Curve::Line {
                        a: q[i],
                        b: q[(i + 1) % 4],
                    })
                    .collect()],
                n: Vec3::new(0., 0., 1.),
                low: z0,
                high: z1,
            },
            1e-6,
        )
        .unwrap()
    }
    #[test]
    fn rectangular_web_slots_repeat_and_preserve_the_fillets() {
        let base = captured("base");
        let mut body = base.clone();
        for y in [5., 15., 25.] {
            body = crate::boolean::cut(&body, &tool(3., 7., y, y + 2., -6., 1.), 1e-6).unwrap();
            body = crate::neutral::Reader::new(&crate::neutral::write_body(&body))
                .body()
                .unwrap();
        }
        assert_eq!(crate::topology::check(&body, 1e-6).unwrap().genus(), 3.);
        assert!(
            (crate::mass::measure(&body, 1e-6).unwrap().volume
                - (155.58846611769187 - 3. * 8. * 0.23))
                .abs()
                < 1e-8
        );
        for (i, f) in base.faces.iter().enumerate() {
            if matches!(f.surface, Surface::Cylinder { .. }) {
                assert_eq!(format!("{f:?}"), format!("{:?}", body.faces[i]));
            }
        }
    }
    #[test]
    fn rectangular_punch_crosses_the_root_with_exact_conic_trims() {
        let body = crate::boolean::cut(&captured("base"), &tool(0.34, 1.14, 5., 7., -6., 1.), 1e-6)
            .unwrap();
        let integral = 0.5 * (0.4 * 0.3 + 0.25 * 0.8_f64.asin());
        let removed = 2. * (0.4 * 1.23 - 2. * integral + 0.4 * 0.23);
        assert_eq!(crate::topology::check(&body, 1e-6).unwrap().genus(), 1.);
        assert!(
            (crate::mass::measure(&body, 1e-6).unwrap().volume - (155.58846611769187 - removed))
                .abs()
                < 1e-8
        );
    }
    #[test]
    fn blind_box_cut_refuses_and_complete_removal_returns_empty() {
        assert!(matches!(
            crate::boolean::cut(&captured("base"), &tool(3., 7., 5., 7., -2., 1.), 1e-6),
            Err(KernelError::Unsupported(_))
        ));
        let empty =
            crate::boolean::cut(&captured("base"), &tool(-1., 11., -1., 37., -6., 1.), 1e-6)
                .unwrap();
        assert!(empty.vertices.is_empty() && empty.edges.is_empty() && empty.faces.is_empty());
    }
    #[test]
    fn rectangular_web_slot_follows_an_analytic_fillet_circle_through_protocol() {
        let n = Vec3::new(0., 0., 1.);
        let circle = crate::extrude::extrude(
            &Profile {
                loops: vec![vec![Curve::Arc {
                    c: Vec3::new(0.74, 5., 0.),
                    n,
                    x: Vec3::new(1., 0., 0.),
                    r: 0.4,
                    t0: 0.,
                    t1: std::f64::consts::TAU,
                }]],
                n,
                low: -6.,
                high: 1.,
            },
            1e-6,
        )
        .unwrap();
        let body = crate::boolean::cut(&captured("base"), &circle, 1e-6).unwrap();
        let mut data = crate::neutral::write_body(&body);
        data.extend(crate::neutral::write_body(&tool(3., 7., 14., 16., -6., 1.)));
        let request = crate::protocol::parse_request(&crate::protocol::encode_request(
            crate::protocol::OP_CUT,
            1e-6,
            &data,
        ))
        .unwrap();
        let out = crate::protocol::handle(&request).unwrap();
        let result = crate::neutral::Reader::new(&out).body().unwrap();
        assert_eq!(crate::topology::check(&result, 1e-6).unwrap().genus(), 2.);
        assert!(
            (crate::mass::measure(&result, 1e-6).unwrap().volume
                - (155.45067891355913 - 8. * 0.23))
                .abs()
                < 1e-8
        );
    }
    #[test]
    fn rectangular_root_punch_selects_each_coordinate_sweep_axis() {
        let rotate = |body: &Body| {
            let map = |v: Vec3| Vec3::new(v.y, v.z, v.x);
            let mut body = body.clone();
            for v in &mut body.vertices {
                *v = map(*v);
            }
            for e in &mut body.edges {
                match &mut e.curve {
                    Curve::Line { a, b } => {
                        *a = map(*a);
                        *b = map(*b);
                    }
                    Curve::Arc { c, n, x, .. } => {
                        *c = map(*c);
                        *n = map(*n);
                        *x = map(*x);
                    }
                    Curve::Conic { c, u, v, .. } => {
                        *c = map(*c);
                        *u = map(*u);
                        *v = map(*v);
                    }
                    Curve::CylinderPair { c, a, x, p, b, .. } => {
                        *c = map(*c);
                        *a = map(*a);
                        *x = map(*x);
                        *p = map(*p);
                        *b = map(*b);
                    }
                }
            }
            for f in &mut body.faces {
                match &mut f.surface {
                    Surface::Plane { p, n, x } => {
                        *p = map(*p);
                        *n = map(*n);
                        *x = map(*x);
                    }
                    Surface::Cylinder { p, a, x, .. } => {
                        *p = map(*p);
                        *a = map(*a);
                        *x = map(*x);
                    }
                }
            }
            body
        };
        let mut base = captured("base");
        let mut cutter = tool(0.34, 1.14, 5., 7., -6., 1.);
        let mut reference: Option<f64> = None;
        for _ in 0..3 {
            let result = crate::boolean::cut(&base, &cutter, 1e-6).unwrap();
            let volume = crate::mass::measure(&result, 1e-6).unwrap().volume;
            if let Some(v) = reference {
                assert!((volume - v).abs() < 1e-8);
            } else {
                reference = Some(volume);
            }
            base = rotate(&base);
            cutter = rotate(&cutter);
        }
    }
}

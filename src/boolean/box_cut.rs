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

// A polygonal prism can also punch a clear planar wall after earlier cuts have
// introduced several cylinder axes. Select by the finite cutter corridor before
// either local passage construction or transverse admission; never retry errors.
fn planar_profile_cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    if tool.faces.len() > 34 || tool.edges.len() > 96 || tool.vertices.len() > 64 {
        return Err(unsupported("planar profile cutter exceeds its work limit"));
    }
    crate::topology::check(tool, tol)?;
    if tool.faces.iter().any(|f| f.loops.len() != 1) {
        return super::endcut::planar_cross_cut(body, tool, tol);
    }
    let axes: Vec<_> = (0..3)
        .filter(|axis| {
            let mut caps = 0;
            tool.faces.iter().all(|f| match f.surface {
                Surface::Plane { n, .. } => {
                    let n = xyz(n);
                    if n[*axis].abs() == 1. && (0..3).all(|k| k == *axis || n[k] == 0.) {
                        caps += 1;
                        true
                    } else {
                        n[*axis] == 0.
                    }
                }
                _ => false,
            }) && caps == 2
        })
        .collect();
    if axes.len() != 1
        || tool
            .edges
            .iter()
            .any(|e| !matches!(e.curve, Curve::Line { .. }))
    {
        return super::endcut::planar_cross_cut(body, tool, tol);
    }
    let axis = axes[0];
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in &tool.vertices {
        for k in 0..3 {
            lo[k] = lo[k].min(xyz(*v)[k]);
            hi[k] = hi[k].max(xyz(*v)[k]);
        }
    }
    let center = vector(std::array::from_fn(|k| (lo[k] + hi[k]) * 0.5));
    let u = (axis + 1) % 3;
    let v = (axis + 2) % 3;
    let radius = (hi[u] - lo[u]).hypot(hi[v] - lo[v]) * 0.5;
    if super::endcut::crosses_cylinder(body, axis, center, radius, lo[axis], hi[axis], tol) {
        return super::endcut::planar_cross_cut(body, tool, tol);
    }
    let lower = tool
        .faces
        .iter()
        .find(|f| {
            matches!(f.surface, Surface::Plane { p, n, .. }
        if xyz(n)[axis].abs() == 1. && (xyz(p)[axis] - lo[axis]).abs() <= tol)
        })
        .ok_or_else(|| unsupported("planar profile has no lower cap"))?;
    if lower.loops.len() != 1 || lower.loops[0].len() < 3 || lower.loops[0].len() > 32 {
        return Err(unsupported(
            "planar profile requires one bounded convex polygon",
        ));
    }
    let outline: Vec<_> = lower.loops[0]
        .iter()
        .map(|(i, forward)| {
            if *forward {
                tool.edges[*i].curve.clone()
            } else {
                tool.edges[*i].curve.reversed()
            }
        })
        .collect();
    let mut normal = [0.; 3];
    normal[axis] = 1.;
    let n = vector(normal);
    let orientation = crate::extrude::signed_area(&outline, n).signum();
    for edge in &outline {
        let direction = edge.end().sub(edge.start());
        for point in outline.iter().map(Curve::start) {
            if direction.cross(point.sub(edge.start())).dot(n) * orientation < 0. {
                return Err(unsupported("planar profile must be convex"));
            }
        }
    }
    let expected = crate::extrude::extrude(
        &Profile {
            loops: vec![outline.clone()],
            n,
            low: 0.,
            high: hi[axis] - lo[axis],
        },
        tol,
    )?;
    super::capsule::match_prism(tool, &expected, tol)?;
    let body = super::passage::coalesce_caps(body, axis, tol)?;
    let passages =
        super::passage::round_passages(&body, axis, center, radius, lo[axis], hi[axis], tol)?;
    let mut result = body;
    for (bottom, top, bottom_face, top_face) in passages {
        let boundary = outline
            .iter()
            .map(|c| c.translated(n.scale(bottom - lo[axis])))
            .collect();
        let cavity = crate::extrude::extrude(
            &Profile {
                loops: vec![boundary],
                n,
                low: 0.,
                high: top - bottom,
            },
            tol,
        )?;
        result = super::passage::splice_cavity(&result, cavity, axis, bottom_face, top_face, tol)?;
    }
    Ok(result)
}

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    if !tol.is_finite() || tol <= 0. {
        return Err(KernelError::InvalidInput(
            "rectangular passage tolerance must be positive".into(),
        ));
    }
    // Select local polygon passage or transverse admission before construction.
    if tool.vertices.len() != 8
        || tool.edges.len() != 12
        || tool.faces.len() != 6
        || tool.faces.iter().any(|f| {
            matches!(f.surface,Surface::Plane {n,..}
            if xyz(n).iter().filter(|v| **v != 0.).count() != 1)
        })
    {
        return planar_profile_cut(body, tool, tol);
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
    #[test]
    fn convex_profile_passages_repeat_on_all_axes_and_refuse_incomplete_depth() {
        let n = Vec3::new(0., 0., 1.);
        let points = [
            Vec3::new(3., 5., 0.),
            Vec3::new(5., 5., 0.),
            Vec3::new(3.5, 7., 0.),
        ];
        let outline: Vec<_> = (0..3)
            .map(|i| Curve::Line {
                a: points[i],
                b: points[(i + 1) % 3],
            })
            .collect();
        let mut tool = crate::extrude::extrude(
            &Profile {
                loops: vec![outline.clone()],
                n,
                low: -6.,
                high: 1.,
            },
            1e-6,
        )
        .unwrap();
        let mut next = crate::extrude::extrude(
            &Profile {
                loops: vec![outline
                    .iter()
                    .map(|c| c.translated(Vec3::new(0., 10., 0.)))
                    .collect()],
                n,
                low: -6.,
                high: 1.,
            },
            1e-6,
        )
        .unwrap();
        let mut blind = crate::extrude::extrude(
            &Profile {
                loops: vec![outline],
                n,
                low: -2.,
                high: 1.,
            },
            1e-6,
        )
        .unwrap();
        let mut base = captured("base");
        let plate = self::tool(0., 10., 0., 20., -4., 0.);
        let punched = crate::boolean::cut(&plate, &tool, 1e-6).unwrap();
        assert!((crate::mass::measure(&punched, 1e-6).unwrap().volume - 792.).abs() < 1e-10);
        for _ in 0..3 {
            let original = crate::neutral::write_body(&base);
            let first = crate::boolean::cut(&base, &tool, 1e-6).unwrap();
            let second = crate::boolean::cut(&first, &next, 1e-6).unwrap();
            let expected = 155.58846611769187 - 2. * 2. * 0.23;
            assert!((crate::mass::measure(&second, 1e-6).unwrap().volume - expected).abs() < 1e-9);
            assert_eq!(crate::topology::check(&second, 1e-6).unwrap().genus(), 2.);
            assert!(matches!(
                crate::boolean::cut(&first, &tool, 1e-6),
                Err(KernelError::Unsupported(_))
            ));
            assert!(matches!(
                crate::boolean::cut(&base, &blind, 1e-6),
                Err(KernelError::Unsupported(_))
            ));
            assert_eq!(original, crate::neutral::write_body(&base));
            base = rotate(&base);
            tool = rotate(&tool);
            next = rotate(&next);
            blind = rotate(&blind);
        }
    }
    #[test]
    fn actual_rect_tube_fifth_profile_clears_partition_seams_and_retains_prior_cuts() {
        let load = |name: &str| {
            let bytes = std::fs::read(
                std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join(format!("tests/fixtures/parts/rect_tube/{name}.f64")),
            )
            .unwrap();
            let data: Vec<_> = bytes
                .chunks_exact(8)
                .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
                .collect();
            crate::neutral::Reader::new(&data).body().unwrap()
        };
        let mut before = load("base");
        for step in 1..=4 {
            let mut data = crate::neutral::write_body(&before);
            data.extend(crate::neutral::write_body(&load(&format!("step{step:02}"))));
            let request = crate::protocol::parse_request(&crate::protocol::encode_request(
                crate::protocol::OP_CUT,
                1e-6,
                &data,
            ))
            .unwrap();
            let out = crate::protocol::handle(&request).unwrap();
            before = crate::neutral::Reader::new(&out).body().unwrap();
        }
        let original = crate::neutral::write_body(&before);
        let merged = super::super::passage::coalesce_caps(&before, 0, 1e-6).unwrap();
        assert!(merged.faces.len() < before.faces.len());
        let old_mass = crate::mass::measure(&before, 1e-6).unwrap();
        let merged_mass = crate::mass::measure(&merged, 1e-6).unwrap();
        assert!((old_mass.volume - merged_mass.volume).abs() < 1e-10);
        assert!((old_mass.area - merged_mass.area).abs() < 1e-10);
        assert_eq!(
            crate::neutral::write_body(&merged),
            crate::neutral::write_body(
                &super::super::passage::coalesce_caps(&merged, 0, 1e-6).unwrap()
            )
        );
        let tool = load("step05");
        let result = crate::boolean::cut(&before, &tool, 1e-6).unwrap();
        assert!(matches!(
            crate::boolean::cut(&result, &tool, 1e-6),
            Err(KernelError::Unsupported(_))
        ));
        assert_eq!(original, crate::neutral::write_body(&before));
        assert_eq!(
            crate::topology::check(&result, 1e-6).unwrap().genus(),
            crate::topology::check(&before, 1e-6).unwrap().genus() + 1.
        );
        // Independent triangle area times the quarter-inch flat wall depth.
        let cap = tool
            .faces
            .iter()
            .find(|f| matches!(f.surface,Surface::Plane { n, .. } if n.x.abs() == 1.))
            .unwrap();
        let vertices: Vec<_> = cap.loops[0]
            .iter()
            .map(|(e, f)| {
                if *f {
                    tool.edges[*e].curve.start()
                } else {
                    tool.edges[*e].curve.end()
                }
            })
            .collect();
        assert_eq!(vertices.len(), 3);
        let area = vertices[1]
            .sub(vertices[0])
            .cross(vertices[2].sub(vertices[0]))
            .length()
            * 0.5;
        let volume = crate::mass::measure(&result, 1e-6).unwrap().volume;
        assert!((old_mass.volume - volume - area * 0.25).abs() < 1e-10);
        assert!((volume - 54.342254292).abs() < 0.001); // Claude's native mesh, 9ebdf7d.
        let mut turned_body = before.clone();
        let mut turned_tool = tool.clone();
        for _ in 0..2 {
            turned_body = rotate(&turned_body);
            turned_tool = rotate(&turned_tool);
            let turned = crate::boolean::cut(&turned_body, &turned_tool, 1e-6).unwrap();
            assert!((crate::mass::measure(&turned, 1e-6).unwrap().volume - volume).abs() < 1e-9);
        }
        let translated = |body: &Body| {
            let offset = Vec3::new(1e5, -1e5, 1e5);
            let mut body = body.clone();
            for p in &mut body.vertices {
                *p = p.add(offset);
            }
            for e in &mut body.edges {
                e.curve = e.curve.translated(offset);
            }
            for f in &mut body.faces {
                match &mut f.surface {
                    Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } => unreachable!("no elliptic-cylinder or revolved faces in these test bodies"),
                    Surface::Plane { p, .. } | Surface::Cylinder { p, .. } => *p = p.add(offset),
                }
            }
            body
        };
        let far = crate::boolean::cut(&translated(&before), &translated(&tool), 1e-6).unwrap();
        assert!((crate::mass::measure(&far, 1e-6).unwrap().volume - volume).abs() < 1e-7);
        for (i, face) in merged.faces.iter().enumerate() {
            if !matches!(face.surface, Surface::Plane { n, .. } if n.x.abs() == 1.) {
                assert_eq!(format!("{face:?}"), format!("{:?}", result.faces[i]));
            }
        }
        // Independent local ray test of the returned planar wall boundaries,
        // compared with the original flat wall minus the captured triangle.
        for ix in 0..17 {
            for iy in 0..31 {
                for iz in 0..29 {
                    let q = Vec3::new(
                        -0.093 + 0.031 * ix as f64,
                        2.613 + 0.083 * iy as f64,
                        -3.913 + 0.081 * iz as f64,
                    );
                    let mut sides = Vec::new();
                    for i in 0..3 {
                        sides.push(
                            vertices[(i + 1) % 3]
                                .sub(vertices[i])
                                .cross(q.sub(vertices[i]))
                                .x,
                        );
                    }
                    let in_triangle =
                        sides.iter().all(|s| *s > 0.) || sides.iter().all(|s| *s < 0.);
                    let expected = q.x > 0. && q.x < 0.25 && !in_triangle;
                    let mut winding = 0;
                    for face in &result.faces {
                        let Surface::Plane { p, n, .. } = face.surface else {
                            continue;
                        };
                        if n.x.abs() != 1. || p.x >= q.x {
                            continue;
                        }
                        let mut inside = false;
                        for &(edge, _) in face.loops.iter().flatten() {
                            let Curve::Line { a, b } = result.edges[edge].curve else {
                                panic!("flat side cap must have linear trims")
                            };
                            if (a.z > q.z) != (b.z > q.z)
                                && q.y < a.y + (q.z - a.z) * (b.y - a.y) / (b.z - a.z)
                            {
                                inside = !inside;
                            }
                        }
                        if inside {
                            winding -= (n.x * if face.sense { 1. } else { -1. }) as i32;
                        }
                    }
                    assert_eq!(winding == 1, expected, "side triangle probe {q:?}");
                }
            }
        }
    }
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
    fn rotate(body: &Body) -> Body {
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
                Curve::Bezier { .. } => unreachable!("no Bezier edges in these test bodies"),
            }
        }
        for f in &mut body.faces {
            match &mut f.surface {
                Surface::EllipticCylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } | Surface::Revolution { .. } => unreachable!("no elliptic-cylinder or revolved faces in these test bodies"),
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
    }
    #[test]
    fn rectangular_root_punch_selects_each_coordinate_sweep_axis() {
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

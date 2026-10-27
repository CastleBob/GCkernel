//! Bounded through-bore construction on a proven plain box or an admitted local
//! planar passage in a curved solid. A full cylindrical cutter must span every
//! intersected material interval and clear unrelated surfaces. The output retains existing faces, adds circular inner cap boundaries,
//! and introduces one inward-facing analytic cylinder. No mesh or offsets.
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;
use std::f64::consts::TAU;

type Result<T> = std::result::Result<T, KernelError>;
fn invalid(m: &str) -> KernelError {
    KernelError::InvalidInput(m.into())
}
fn unsupported(m: &str) -> KernelError {
    KernelError::Unsupported(m.into())
}
fn xyz(p: Vec3) -> [f64; 3] {
    [p.x, p.y, p.z]
}
fn vec(p: [f64; 3]) -> Vec3 {
    Vec3::new(p[0], p[1], p[2])
}

fn cardinal(v: Vec3) -> Result<(usize, f64)> {
    if !v.is_finite() {
        return Err(invalid("nonfinite bore direction"));
    }
    let p = xyz(v);
    let axes: Vec<usize> = (0..3).filter(|i| p[*i] != 0.).collect();
    if axes.len() != 1 {
        return Err(unsupported(
            "contained bore requires a coordinate-axis cylinder",
        ));
    }
    if p[axes[0]].abs() != 1. {
        return Err(invalid("bore direction is not unit length"));
    }
    Ok((axes[0], p[axes[0]]))
}

struct Cylinder {
    axis: usize,
    center: Vec3,
    radius: f64,
    low: f64,
    high: f64,
}

fn cylinder(tool: &Body, tolerance: f64) -> Result<Cylinder> {
    if tool.vertices.len() != 2 || tool.edges.len() != 3 || tool.faces.len() != 3 {
        return Err(unsupported(
            "bore cutter requires two full-circle cap edges and one seam",
        ));
    }
    if tool.vertices.iter().any(|v| !v.is_finite())
        || tool.edges.iter().any(|e| e.start >= 2 || e.end >= 2)
        || tool
            .faces
            .iter()
            .flat_map(|f| &f.loops)
            .flatten()
            .any(|(i, _)| *i >= 3)
    {
        return Err(invalid("invalid cylinder cutter topology"));
    }
    let cylinders: Vec<usize> = tool
        .faces
        .iter()
        .enumerate()
        .filter_map(|(i, f)| matches!(f.surface, Surface::Cylinder { .. }).then_some(i))
        .collect();
    if cylinders.len() != 1 {
        return Err(unsupported("bore cutter requires one cylindrical wall"));
    }
    let side = &tool.faces[cylinders[0]];
    let Surface::Cylinder { p, a, x, r } = side.surface else {
        unreachable!()
    };
    let (axis, _) = cardinal(a)?;
    if !p.is_finite()
        || !x.is_finite()
        || !r.is_finite()
        || r <= tolerance
        || xyz(x)[axis] != 0.
        || (x.length() - 1.).abs() > 1e-12
    {
        return Err(invalid("invalid cylinder frame or radius"));
    }
    if !side.sense || side.loops.len() != 1 || side.loops[0].len() != 4 {
        return Err(invalid("cylinder wall is not an outward closed boundary"));
    }
    let mut cap_data = Vec::new();
    let mut arc_indices = Vec::new();
    for (i, face) in tool.faces.iter().enumerate() {
        if i == cylinders[0] {
            continue;
        }
        let Surface::Plane { p: fp, n, x: fx } = face.surface else {
            return Err(unsupported("bore cutter requires planar caps"));
        };
        let (cap_axis, cap_direction) = cardinal(n)?;
        if cap_axis != axis
            || !fp.is_finite()
            || !fx.is_finite()
            || xyz(fx)[axis] != 0.
            || (fx.length() - 1.).abs() > 1e-12
            || face.loops.len() != 1
            || face.loops[0].len() != 1
        {
            return Err(invalid("invalid cylindrical cutter cap"));
        }
        let (index, forward) = face.loops[0][0];
        let edge = &tool.edges[index];
        let Curve::Arc {
            c,
            n: an,
            x: ax,
            r: ar,
            t0,
            t1,
        } = edge.curve
        else {
            return Err(unsupported("bore cap must be a full circle"));
        };
        let (arc_axis, arc_direction) = cardinal(an)?;
        if edge.start != edge.end
            || arc_axis != axis
            || !c.is_finite()
            || !ax.is_finite()
            || !t0.is_finite()
            || !t1.is_finite()
            || t1 - t0 != TAU
            || ar != r
            || xyz(ax)[axis] != 0.
            || (ax.length() - 1.).abs() > 1e-12
            || (0..3).any(|k| k != axis && xyz(c)[k] != xyz(p)[k])
            || xyz(fp)[axis] != xyz(c)[axis]
        {
            return Err(invalid(
                "cutter cap is not a circle on its cylindrical wall",
            ));
        }
        let outward = cap_direction * if face.sense { 1. } else { -1. };
        let travel = arc_direction * if forward { 1. } else { -1. };
        if outward != travel {
            return Err(invalid("cylinder cap loop has reversed orientation"));
        }
        cap_data.push((xyz(c)[axis], outward));
        arc_indices.push(index);
    }
    cap_data.sort_by(|a, b| a.0.total_cmp(&b.0));
    if cap_data[0].0 >= cap_data[1].0
        || cap_data[0].1 != -1.
        || cap_data[1].1 != 1.
        || arc_indices[0] == arc_indices[1]
    {
        return Err(invalid("cylinder cap ordering or orientation is invalid"));
    }
    let seam_index = (0..3).find(|i| !arc_indices.contains(i)).unwrap();
    let seam = &tool.edges[seam_index];
    let Curve::Line { a: sa, b: sb } = seam.curve else {
        return Err(invalid("cylinder seam must be straight"));
    };
    if !sa.is_finite()
        || !sb.is_finite()
        || seam.start == seam.end
        || (0..3).any(|i| i != axis && xyz(sa)[i] != xyz(sb)[i])
        || xyz(sa)[axis].min(xyz(sb)[axis]) != cap_data[0].0
        || xyz(sa)[axis].max(xyz(sb)[axis]) != cap_data[1].0
    {
        return Err(invalid("cylinder seam does not connect the cap circles"));
    }
    crate::topology::check(tool, tolerance)?;
    Ok(Cylinder {
        axis,
        center: p,
        radius: r,
        low: cap_data[0].0,
        high: cap_data[1].0,
    })
}

fn box_passage(body: &Body, cutter: &Cylinder, tolerance: f64) -> Result<(f64, f64, usize, usize)> {
    let (low, high) = super::rectilinear::plain_box_bounds(body, tolerance)?;
    let axis = cutter.axis;
    let center = xyz(cutter.center);
    let radius = cutter.radius;
    if cutter.low > low[axis] || cutter.high < high[axis] {
        return Err(unsupported(
            "only full through-bores are supported; cutter must span the box",
        ));
    }
    if (0..3).filter(|i| *i != axis).any(|i| {
        center[i] - radius - low[i] <= tolerance || high[i] - center[i] - radius <= tolerance
    }) {
        return Err(unsupported(
            "bore must remain strictly inside the box side walls",
        ));
    }
    let mut caps = [None, None];
    for (i, face) in body.faces.iter().enumerate() {
        let Surface::Plane { p, n, .. } = face.surface else {
            return Err(unsupported("bore base must be a plain box"));
        };
        let (face_axis, direction) = cardinal(n)?;
        if face_axis != axis {
            continue;
        }
        let sign = direction * if face.sense { 1. } else { -1. };
        let which = if sign < 0. { 0 } else { 1 };
        if xyz(p)[axis] != if which == 0 { low[axis] } else { high[axis] } || caps[which].is_some()
        {
            return Err(invalid("box cap orientation is inconsistent"));
        }
        caps[which] = Some(i);
    }
    let [Some(bottom_face), Some(top_face)] = caps else {
        return Err(invalid("box is missing a bore cap"));
    };
    Ok((low[axis], high[axis], bottom_face, top_face))
}

pub(super) fn cut(body: &Body, tool: &Body, tolerance: f64) -> Result<Body> {
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(invalid("bore tolerance must be finite and positive"));
    }
    let cutter = cylinder(tool, tolerance)?;
    let axis = cutter.axis;
    let center = xyz(cutter.center);
    let radius = cutter.radius;
    if super::endcut::crosses_cylinder(
        body,
        axis,
        cutter.center,
        radius,
        cutter.low,
        cutter.high,
        tolerance,
    ) {
        return super::endcut::cross_cut(body, tool, tolerance);
    }
    // Admission is selected by the operand representation before execution.
    // A failed box admission never retries against the general passage checker.
    let passages = if body.vertices.len() == 8 && body.edges.len() == 12 && body.faces.len() == 6 {
        vec![box_passage(body, &cutter, tolerance)?]
    } else {
        super::passage::round_passages(
            body,
            axis,
            cutter.center,
            radius,
            cutter.low,
            cutter.high,
            tolerance,
        )?
    };
    let mut result = body.clone();
    let mut normal = [0.; 3];
    normal[axis] = 1.;
    let normal = vec(normal);
    let mut reference = [0.; 3];
    reference[(axis + 1) % 3] = 1.;
    let reference = vec(reference);
    for (low, high, bottom_face, top_face) in passages {
        let mut bottom = center;
        bottom[axis] = low;
        let bottom = vec(bottom);
        let mut top = center;
        top[axis] = high;
        let top = vec(top);
        let vb = result.vertices.len();
        let vt = vb + 1;
        let pb = bottom.add(reference.scale(radius));
        let pt = top.add(reference.scale(radius));
        result.vertices.extend([pb, pt]);
        let eb = result.edges.len();
        let et = eb + 1;
        let seam = eb + 2;
        result.edges.push(Edge {
            curve: Curve::Arc {
                c: bottom,
                n: normal,
                x: reference,
                r: radius,
                t0: 0.,
                t1: TAU,
            },
            start: vb,
            end: vb,
        });
        result.edges.push(Edge {
            curve: Curve::Arc {
                c: top,
                n: normal,
                x: reference,
                r: radius,
                t0: 0.,
                t1: TAU,
            },
            start: vt,
            end: vt,
        });
        result.edges.push(Edge {
            curve: Curve::Line { a: pb, b: pt },
            start: vb,
            end: vt,
        });
        result.faces[bottom_face].loops.push(vec![(eb, true)]);
        result.faces[top_face].loops.push(vec![(et, false)]);
        result.faces.push(Face {
            surface: Surface::Cylinder {
                p: bottom,
                a: normal,
                x: reference,
                r: radius,
            },
            sense: false,
            loops: vec![vec![(eb, false), (seam, true), (et, true), (seam, false)]],
        });
    }
    crate::topology::check(&result, tolerance)
        .map_err(|e| KernelError::Internal(format!("bore result failed topology checks: {e}")))?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::neutral::Profile;

    // Independently integrate the RETURNED boundary via the divergence theorem.
    // Planar face contributions use their oriented loop areas; a full cylinder
    // contributes its actual radius and edge-derived axial span. No stock/cutter
    // recipe or boolean construction data is used.
    fn boundary_volume(body: &Body) -> f64 {
        let mut volume = 0.;
        for face in &body.faces {
            let sign = if face.sense { 1. } else { -1. };
            match face.surface {
                Surface::Plane { p, n, .. } => {
                    let outward = n.scale(sign);
                    let mut area = 0.;
                    for lp in &face.loops {
                        for (index, forward) in lp {
                            let travel = if *forward { 1. } else { -1. };
                            match body.edges[*index].curve {
                                Curve::Conic { .. } | Curve::CylinderPair { .. } => {
                                    panic!("round-bore fixture unexpectedly contains a conic")
                                }
                                Curve::Line { a, b } => {
                                    area += 0.5 * a.cross(b).dot(outward) * travel
                                }
                                Curve::Arc {
                                    n: an, r, t0, t1, ..
                                } => {
                                    assert_eq!(t1 - t0, TAU);
                                    area += std::f64::consts::PI * r * r * an.dot(outward) * travel;
                                }
                            }
                        }
                    }
                    assert!(area > 0.);
                    volume += p.dot(outward) * area / 3.;
                }
                Surface::Cylinder { a, r, .. } => {
                    let heights: Vec<f64> = face
                        .loops
                        .iter()
                        .flatten()
                        .flat_map(|(i, _)| {
                            let edge = &body.edges[*i];
                            [
                                body.vertices[edge.start].dot(a),
                                body.vertices[edge.end].dot(a),
                            ]
                        })
                        .collect();
                    let low = heights.iter().copied().fold(f64::INFINITY, f64::min);
                    let high = heights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                    volume += sign * TAU * r * r * (high - low) / 3.;
                }
            }
        }
        volume
    }

    fn fixture(axis: usize, center: [f64; 3], radius: f64, depth: (f64, f64)) -> (Body, Body) {
        let mut n = [0.; 3];
        n[axis] = 1.;
        let n = vec(n);
        let u = (axis + 1) % 3;
        let v = (axis + 2) % 3;
        let points: Vec<Vec3> = [(0., 0.), (10., 0.), (10., 8.), (0., 8.)]
            .iter()
            .map(|(a, b)| {
                let mut p = [0.; 3];
                p[u] = *a;
                p[v] = *b;
                vec(p)
            })
            .collect();
        let lines = (0..4)
            .map(|i| Curve::Line {
                a: points[i],
                b: points[(i + 1) % 4],
            })
            .collect();
        let body = crate::extrude::extrude(
            &Profile {
                loops: vec![lines],
                n,
                low: 0.,
                high: 4.,
            },
            1e-7,
        )
        .unwrap();
        let mut x = [0.; 3];
        x[u] = 1.;
        let circle = Curve::Arc {
            c: vec(center),
            n,
            x: vec(x),
            r: radius,
            t0: 0.,
            t1: TAU,
        };
        let tool = crate::extrude::extrude(
            &Profile {
                loops: vec![vec![circle]],
                n,
                low: depth.0,
                high: depth.1,
            },
            1e-7,
        )
        .unwrap();
        (body, tool)
    }

    #[test]
    fn round_through_hole_retains_exact_cylinder_and_cap_holes_on_every_axis() {
        for axis in 0..3 {
            let mut center = [0.; 3];
            center[(axis + 1) % 3] = 5.;
            center[(axis + 2) % 3] = 4.;
            let (body, tool) = fixture(axis, center, 1., (-1., 5.));
            let original = crate::neutral::write_body(&body);
            let result = crate::boolean::cut(&body, &tool, 1e-7).unwrap();
            assert_eq!(crate::neutral::write_body(&body), original);
            let report = crate::topology::check(&result, 1e-7).unwrap();
            assert_eq!(report.genus(), 1.);
            assert_eq!(result.faces.len(), 7);
            assert_eq!(
                result.faces.iter().filter(|f| f.loops.len() == 2).count(),
                2
            );
            let wall = result.faces.last().unwrap();
            assert!(!wall.sense);
            let Surface::Cylinder { p, a, r, .. } = wall.surface else {
                panic!("expected cylinder")
            };
            assert_eq!(r, 1.);
            assert_eq!(xyz(p)[axis], 0.);
            assert_eq!(xyz(a)[axis], 1.);
            assert_eq!(xyz(p), center);
            let retained = boundary_volume(&result);
            assert!((retained - (320. - 4. * std::f64::consts::PI)).abs() < 1e-10);
            let payload = crate::neutral::write_body(&result);
            let mut reader = crate::neutral::Reader::new(&payload);
            let returned = reader.body().unwrap();
            reader.finish().unwrap();
            crate::topology::check(&returned, 1e-7).unwrap();
        }
    }

    #[test]
    fn boundary_crossing_tangent_and_blind_bores_are_explicitly_unsupported() {
        for (x, depth) in [(0., (-1., 5.)), (1., (-1., 5.)), (5., (2., 5.))] {
            let (body, tool) = fixture(2, [x, 4., 0.], 1., depth);
            assert!(matches!(
                crate::boolean::cut(&body, &tool, 1e-7),
                Err(KernelError::Unsupported(_))
            ));
        }
    }

    #[test]
    fn bad_cap_radius_and_direction_are_refused() {
        let (body, mut tool) = fixture(2, [5., 4., 0.], 1., (-1., 5.));
        if let Curve::Arc { r, .. } = &mut tool.edges[0].curve {
            *r = 2.;
        }
        assert!(matches!(
            crate::boolean::cut(&body, &tool, 1e-7),
            Err(KernelError::InvalidInput(_))
        ));
        let (body, mut tool) = fixture(2, [5., 4., 0.], 1., (-1., 5.));
        tool.faces[0].sense = !tool.faces[0].sense;
        assert!(matches!(
            crate::boolean::cut(&body, &tool, 1e-7),
            Err(KernelError::InvalidInput(_))
        ));
    }

    #[test]
    fn contained_round_cut_succeeds_through_wire_protocol() {
        let (body, tool) = fixture(2, [5., 4., 0.], 1., (-1., 5.));
        let mut data = crate::neutral::write_body(&body);
        data.extend(crate::neutral::write_body(&tool));
        let request = crate::protocol::parse_request(&crate::protocol::encode_request(
            crate::protocol::OP_CUT,
            1e-7,
            &data,
        ))
        .unwrap();
        let result = crate::protocol::handle(&request).unwrap();
        let back = crate::neutral::Reader::new(&result).body().unwrap();
        assert_eq!(crate::topology::check(&back, 1e-7).unwrap().genus(), 1.);
        assert!((boundary_volume(&back) - (320. - 4. * std::f64::consts::PI)).abs() < 1e-9);
    }

    #[test]
    fn flush_axial_caps_work_but_recutting_the_same_hole_is_not_yet_admitted() {
        let (body, tool) = fixture(2, [5., 4., 0.], 1., (0., 4.));
        let result = crate::boolean::cut(&body, &tool, 1e-7).unwrap();
        assert!((boundary_volume(&result) - (320. - 4. * std::f64::consts::PI)).abs() < 1e-9);
        assert!(matches!(
            crate::boolean::cut(&result, &tool, 1e-7),
            Err(KernelError::Unsupported(_))
        ));
    }
}

//! Independent op-9 acceptance. The request is the documented public wire format,
//! not the revolve implementation's private constructors or integration routines.
//! No shared worker is rebuilt; expected values come from elementary solids and
//! independent Simpson integration of explicit meridians where area is not closed form.
use std::collections::{HashMap, HashSet};
use std::f64::consts::{PI, TAU};
use gc_kernel::geom::Vec3;
use gc_kernel::neutral::{put_surface, write_body, Body, Reader};
use gc_kernel::{mass, protocol, tessellate, topology};

const TOL: f64 = 1e-7;
const CHORD: f64 = 0.04;
type Loop = Vec<Vec<f64>>;

fn line(r: f64, z: f64) -> Vec<f64> { vec![0., r, z] }
fn arc(r: f64, z: f64, cr: f64, cz: f64, dir: f64) -> Vec<f64> {
    vec![1., r, z, cr, cz, dir]
}
fn ellipse(r: f64, z: f64, cr: f64, cz: f64, a: f64, b: f64) -> Vec<f64> {
    vec![2., r, z, cr, cz, a, 0., 0., b, 1.]
}
fn rectangle(r0: f64, r1: f64, z0: f64, z1: f64) -> Loop {
    vec![line(r0, z0), line(r1, z0), line(r1, z1), line(r0, z1)]
}

#[derive(Clone, Copy)]
struct Frame { p: Vec3, a: Vec3, x: Vec3 }
fn frame() -> Frame {
    Frame { p: Vec3::new(0., 0., 0.), a: Vec3::new(0., 1., 0.), x: Vec3::new(1., 0., 0.) }
}
fn request(loops: &[Loop], angle: f64, f: Frame) -> protocol::Request {
    let mut data = vec![f.p.x, f.p.y, f.p.z, f.a.x, f.a.y, f.a.z, f.x.x, f.x.y, f.x.z,
        angle, loops.len() as f64];
    for lp in loops {
        data.push(lp.len() as f64);
        for segment in lp { data.extend(segment); }
    }
    protocol::Request { op: 9, units: 1, tolerance: TOL, data }
}
fn near(actual: f64, expected: f64) {
    assert!(actual.is_finite() && (actual - expected).abs() <= 2e-8 * expected.abs().max(1.),
        "{actual} != independent {expected}");
}
fn simpson(f: impl Fn(f64) -> f64, lo: f64, hi: f64) -> f64 {
    const N: usize = 8192;
    let step = (hi - lo) / N as f64;
    let mut sum = f(lo) + f(hi);
    for i in 1..N { sum += (if i % 2 == 0 { 2. } else { 4. }) * f(lo + i as f64 * step); }
    sum * step / 3.
}

/// The callback measures distance from a point in the meridian frame (rho,z) to
/// its known lateral surface. Planar caps are checked against their own support.
fn mesh_checks(body: &Body, f: Frame, angle: f64, volume: f64, area: f64,
    lateral_distance: Option<&dyn Fn(f64, f64) -> f64>) {
    let before = write_body(body);
    let mesh = tessellate::tessellate(body, CHORD).expect("valid revolved body must tessellate");
    assert_eq!(before, write_body(body));
    assert!(!mesh.triangles.is_empty());
    assert_eq!(mesh.faces.len(), mesh.triangles.len());
    assert_eq!(mesh.normals.len(), mesh.positions.len());
    assert!(mesh.positions.iter().chain(&mesh.normals).all(|p| p.is_finite()));
    let bits = |v: f64| if v == 0. { 0 } else { v.to_bits() };
    let mut ids = HashMap::new();
    let welded: Vec<_> = mesh.positions.iter().map(|p| {
        let next = ids.len();
        *ids.entry((bits(p.x), bits(p.y), bits(p.z))).or_insert(next)
    }).collect();
    let mut edges: HashMap<(usize, usize), (usize, i32)> = HashMap::new();
    let mut links: HashMap<usize, HashMap<usize, Vec<usize>>> = HashMap::new();
    let origin = mesh.positions[0];
    let mut mesh_volume = 0.;
    for (tri, &face) in mesh.triangles.iter().zip(&mesh.faces) {
        assert!(tri.iter().all(|&i| i < mesh.positions.len()));
        assert!(face < body.faces.len());
        let [a, b, c] = tri.map(|i| mesh.positions[i]);
        let [u, v, w] = tri.map(|i| welded[i]);
        assert!(u != v && v != w && w != u, "collapsed triangle");
        let cross = b.sub(a).cross(c.sub(a));
        assert!(cross.length() > 1e-13, "zero-area triangle");
        let normal = tri.iter().fold(Vec3::new(0., 0., 0.), |sum, &i| sum.add(mesh.normals[i]));
        assert!(cross.dot(normal) > 0., "triangle runs against its outward normals");
        mesh_volume += a.sub(origin).dot(b.sub(origin).cross(c.sub(origin))) / 6.;
        for (s, e) in [(u, v), (v, w), (w, u)] {
            let value = edges.entry((s.min(e), s.max(e))).or_default();
            value.0 += 1; value.1 += if s < e { 1 } else { -1 };
        }
        for (at, left, right) in [(u, v, w), (v, w, u), (w, u, v)] {
            let link = links.entry(at).or_default();
            link.entry(left).or_default().push(right);
            link.entry(right).or_default().push(left);
        }
        let mut support = Vec::new();
        put_surface(&mut support, &body.faces[face].surface);
        for q in [a.add(b).scale(0.5), b.add(c).scale(0.5), c.add(a).scale(0.5), a.add(b).add(c).scale(1. / 3.)] {
            if support[0] == 0. {
                let p = Vec3::new(support[1], support[2], support[3]);
                let n = Vec3::new(support[4], support[5], support[6]);
                assert!(q.sub(p).dot(n).abs() <= TOL);
            } else if let Some(distance) = lateral_distance {
                let d = q.sub(f.p); let z = d.dot(f.a);
                let rho = d.sub(f.a.scale(z)).length();
                assert!(distance(rho, z) <= CHORD + 1e-8, "lateral mesh chord bound exceeded");
            }
        }
    }
    assert!(edges.values().all(|&(uses, winding)| uses == 2 && winding == 0),
        "mesh boundary is not closed and consistently oriented");
    for (vertex, link) in &links {
        assert!(link.values().all(|adjacent| adjacent.len() == 2), "nonmanifold vertex {vertex}");
        let mut pending = vec![*link.keys().next().unwrap()];
        let mut seen = HashSet::new();
        while let Some(v) = pending.pop() {
            if seen.insert(v) { pending.extend(link[&v].iter().copied()); }
        }
        assert_eq!(seen.len(), link.len(), "pinched/disconnected vertex fan {vertex}");
    }
    assert!(mesh_volume > 0.);
    assert!((mesh_volume - volume).abs() <= CHORD * area + 1e-8);
    if angle.abs() < TAU - 1e-10 {
        let y = f.a.cross(f.x);
        for p in &mesh.positions {
            let d = p.sub(f.p);
            let (x, y) = (d.dot(f.x), d.dot(y));
            if x.hypot(y) <= TOL { continue; }
            let t = (angle.signum() * y.atan2(x)).rem_euclid(TAU);
            assert!(t <= angle.abs() + 1e-8 || TAU - t <= 1e-8,
                "point lies outside signed sweep {angle}: azimuth {t}");
        }
    }
}

fn verify(loops: &[Loop], angle: f64, f: Frame, expected_volume: f64,
    expected_area: f64, exact: bool, distance: Option<&dyn Fn(f64, f64) -> f64>) {
    let reply = protocol::handle(&request(loops, angle, f)).expect("valid op-9 section must build");
    near(*reply.last().unwrap(), expected_volume);
    let wire = &reply[..reply.len() - 1];
    let mut reader = Reader::with_tolerance(wire, TOL);
    let body = reader.body().unwrap(); reader.finish().unwrap();
    topology::check(&body, TOL).unwrap();
    let m = mass::measure(&body, TOL).unwrap();
    near(m.volume, expected_volume); near(m.area, expected_area);
    assert_eq!(m.exact, exact, "mass closed-form/quadrature flag must be honest");
    let restored = write_body(&body);
    assert_eq!(wire.len(), restored.len());
    for (&a, &b) in wire.iter().zip(&restored) {
        assert!((a - b).abs() <= 1e-12 * a.abs().max(1.));
    }
    mesh_checks(&body, f, angle, expected_volume, expected_area, distance);
}

fn angles() -> [f64; 5] { [PI / 2., -PI / 2., PI, 3. * PI / 2., TAU] }
fn area_for_angle(full_area: f64, section_area: f64, angle: f64) -> f64 {
    full_area * angle.abs() / TAU + if angle.abs() < TAU { 2. * section_area } else { 0. }
}

#[test]
fn partial_frusta_have_signed_sweeps_and_closed_planar_section_caps() {
    let loops = vec![vec![line(0., 0.), line(3., 0.), line(1., 4.), line(0., 4.)]];
    for angle in angles() {
        let f = if angle == PI { Frame { p: Vec3::new(17., -23., 91.),
            a: Vec3::new(1., 2., 2.).unit().unwrap(), x: Vec3::new(2., -1., 0.).unit().unwrap() } } else { frame() };
        verify(&loops, angle, f, 52. * PI / 3. * angle.abs() / TAU,
            area_for_angle(PI * (10. + 4. * 20_f64.sqrt()), 8., angle), true,
            Some(&|rho, z| (rho - (3. - z / 2.)).abs() / 1.25_f64.sqrt()));
    }
}

#[test]
fn cone_sectors_share_the_apex_with_both_section_caps() {
    let loops = vec![vec![line(0., 0.), line(2., 0.), line(0., 3.)]];
    for angle in angles() {
        verify(&loops, angle, frame(), 4. * PI * angle.abs() / TAU,
            area_for_angle(PI * (4. + 2. * 13_f64.sqrt()), 3., angle), true,
            Some(&|rho, z| (rho - (2. - 2. * z / 3.)).abs() / (13_f64 / 9.).sqrt()));
    }
}

#[test]
fn sphere_sectors_keep_both_poles_manifold_and_have_analytic_mass() {
    let r = 2.;
    let loops = vec![vec![arc(0., -r, 0., 0., 1.), line(0., r)]];
    for angle in angles() {
        verify(&loops, angle, frame(), 4. * PI * r * r * r / 3. * angle.abs() / TAU,
            area_for_angle(4. * PI * r * r, PI * r * r / 2., angle), true,
            Some(&|rho, z| (rho.hypot(z) - r).abs()));
    }
}

#[test]
fn circular_tori_and_partial_elbows_keep_periodic_seams_and_caps_closed() {
    let (big, small) = (4., 1.);
    let loops = vec![vec![arc(big + small, 0., big, 0., 1.), arc(big - small, 0., big, 0., 1.)]];
    for angle in angles() {
        verify(&loops, angle, frame(), 2. * PI * PI * big * small * small * angle.abs() / TAU,
            area_for_angle(4. * PI * PI * big * small, PI * small * small, angle), true,
            Some(&|rho, z| ((rho - big).hypot(z) - small).abs()));
    }
}

#[test]
fn sphere_and_torus_hide_periodic_divisions_but_keep_partial_cap_boundaries() {
    let sphere = vec![vec![arc(0., -2., 0., 0., 1.), line(0., 2.)]];
    let torus = vec![vec![arc(5., 0., 4., 0., 1.), arc(3., 0., 4., 0., 1.)]];
    for (name, loops) in [("sphere", sphere), ("torus", torus)] {
        // Quarter sectors have two genuinely distinct cap supports; the half-turn
        // sphere's coplanar cap division would itself be an invisible internal edge.
        for angle in [TAU, PI / 2., -PI / 2.] {
            let reply = protocol::handle(&request(&loops, angle, frame())).unwrap();
            let mut reader = Reader::with_tolerance(&reply[..reply.len() - 1], TOL);
            let body = reader.body().unwrap(); reader.finish().unwrap();
            let mesh = tessellate::tessellate(&body, CHORD).unwrap();
            let expected: HashSet<usize> = body.faces.iter().filter(|face| {
                let mut surface = Vec::new();
                put_surface(&mut surface, &face.surface);
                surface[0] == 0.
            }).flat_map(|face| face.loops.iter().flatten().map(|&(edge, _)| edge)).collect();
            if angle == TAU {
                assert!(expected.is_empty(), "full {name} has no planar caps");
            } else {
                assert!(!expected.is_empty(), "partial {name} must retain cap boundaries");
            }
            let displayed: HashSet<_> = mesh.display_edges.iter().copied().collect();
            assert_eq!(displayed.len(), mesh.display_edges.len(), "duplicate display edge");
            assert_eq!(displayed, expected,
                "{name} sweep {angle}: smooth periodic divisions must be hidden and real cap edges retained");
        }
    }
}

#[test]
fn generic_revolution_hides_periodic_seams_but_keeps_cap_boundaries() {
    let torus = vec![vec![ellipse(5.2, 0., 4., 0., 1.2, 0.7),
        ellipse(2.8, 0., 4., 0., 1.2, 0.7)]];
    // A single cubic wall, with positive radii at both ends. Its two radial
    // closing segments produce real planar end caps even for a full turn.
    // This makes no claim about merging independently authored Bezier pieces.
    let cubic = vec![vec![vec![3., 1., 0., 1., 1., 2., 2.],
        line(2., 3.), line(0., 3.), line(0., 0.)]];
    for (name, loops) in [("elliptical torus", torus), ("single cubic wall", cubic)] {
        for angle in [TAU, PI / 2., -PI / 2.] {
            let reply = protocol::handle(&request(&loops, angle, frame())).unwrap();
            let mut reader = Reader::with_tolerance(&reply[..reply.len() - 1], TOL);
            let body = reader.body().unwrap(); reader.finish().unwrap();
            let mesh = tessellate::tessellate(&body, CHORD).unwrap();
            let expected: HashSet<usize> = body.faces.iter().filter(|face| {
                let mut surface = Vec::new();
                put_surface(&mut surface, &face.surface);
                surface[0] == 0.
            }).flat_map(|face| face.loops.iter().flatten().map(|&(edge, _)| edge)).collect();
            if name == "elliptical torus" && angle == TAU {
                assert!(expected.is_empty(), "full elliptical torus has no planar caps");
            } else {
                assert!(!expected.is_empty(), "real planar cap boundaries must be retained");
            }
            let displayed: HashSet<_> = mesh.display_edges.iter().copied().collect();
            assert_eq!(displayed.len(), mesh.display_edges.len(), "duplicate display edge");
            assert_eq!(displayed, expected,
                "{name} sweep {angle}: generic periodic seams must be hidden and real cap edges retained");
        }
    }
}

#[test]
fn valid_trimmed_spindle_arc_is_not_rejected_as_a_complete_self_intersecting_torus() {
    let h = 3_f64.sqrt();
    let loops = vec![vec![arc(0., -h, 1., 0., 1.), line(0., h)]];
    // Circle r=1+2*cos(u), clipped to r>=0. The planar region's radial first
    // moment is 8*pi/3 + 3*sqrt(3), so Pappus gives twice pi times that value.
    verify(&loops, TAU, frame(), 16. * PI * PI / 3. + 6. * PI * h,
        16. * PI * PI / 3. + 8. * PI * h, true,
        Some(&|rho, z| ((rho - 1.).hypot(z) - 2.).abs()));
}

#[test]
fn spheroid_and_elliptical_torus_preserve_conic_meridians() {
    let (a, c) = (2_f64, 3_f64);
    let eccentricity = (1. - a * a / (c * c)).sqrt();
    let full_area = 2. * PI * a * a * (1. + c / (a * eccentricity) * eccentricity.asin());
    let spheroid = vec![vec![ellipse(0., -c, 0., 0., a, c), line(0., c)]];
    for angle in [PI / 2., -PI, TAU] {
        verify(&spheroid, angle, frame(), 4. * PI * a * a * c / 3. * angle.abs() / TAU,
            area_for_angle(full_area, PI * a * c / 2., angle), false, None);
    }
    let (big, a, b) = (4., 1.2, 0.7);
    let circumference = simpson(|t| (a * t.sin()).hypot(b * t.cos()), 0., TAU);
    let torus = vec![vec![ellipse(big + a, 0., big, 0., a, b), ellipse(big - a, 0., big, 0., a, b)]];
    for angle in [PI / 2., -PI / 2., TAU] {
        verify(&torus, angle, frame(), 2. * PI * PI * big * a * b * angle.abs() / TAU,
            area_for_angle(TAU * big * circumference, PI * a * b, angle), false, None);
    }
}

#[test]
fn cubic_bezier_meridian_has_exact_polynomial_volume_and_independent_area() {
    // z(t)=3t, r(t)=1+3t^2-2t^3. Integral r^2 dz = 249/35;
    // section area = integral r dz = 9/2. Circular end caps add 5*pi.
    let loops = vec![vec![vec![3., 1., 0., 1., 1., 2., 2.], line(2., 3.), line(0., 3.), line(0., 0.)]];
    let lateral = TAU * simpson(|t| {
        let r = 1. + 3. * t * t - 2. * t * t * t;
        r * (6. * t - 6. * t * t).hypot(3.)
    }, 0., 1.);
    for angle in [PI / 2., -PI / 2., TAU] {
        verify(&loops, angle, frame(), 249. * PI / 35. * angle.abs() / TAU,
            area_for_angle(lateral + 5. * PI, 4.5, angle), false, None);
    }
}

#[test]
fn enclosed_profile_hole_subtracts_volume_and_adds_its_inner_boundary_area() {
    // Outer rectangle r=1..4,z=0..5; inner rectangular void r=2..3,z=1..4.
    // The supplied loops have the SAME orientation, which the contract permits.
    let loops = vec![rectangle(1., 4., 0., 5.), rectangle(2., 3., 1., 4.)];
    for angle in [PI / 2., -3. * PI / 2., TAU] {
        verify(&loops, angle, frame(), 60. * PI * angle.abs() / TAU,
            area_for_angle(120. * PI, 12., angle), true, None);
    }
}

#[test]
fn small_and_nearly_full_hollow_sectors_retain_distinct_planar_caps() {
    let loops = vec![rectangle(1., 3., 0., 4.)];
    for angle in [PI / 180., TAU - 1e-5] {
        let f = frame();
        verify(&loops, angle, f, 32. * PI * angle / TAU,
            area_for_angle(48. * PI, 8., angle), true, None);
        let reply = protocol::handle(&request(&loops, angle, f)).unwrap();
        let mut reader = Reader::with_tolerance(&reply[..reply.len() - 1], TOL);
        let body = reader.body().unwrap(); reader.finish().unwrap();
        // Axial end faces have normal parallel to a. The two additional section
        // caps have normals perpendicular to a; neither may disappear just because
        // their angular separation is much smaller than the display chord.
        let caps: Vec<_> = body.faces.iter().filter_map(|face| {
            let mut support = Vec::new();
            put_surface(&mut support, &face.surface);
            if support[0] != 0. { return None; }
            let normal = Vec3::new(support[4], support[5], support[6]);
            (normal.dot(f.a).abs() < 1e-10).then_some(normal)
        }).collect();
        assert_eq!(caps.len(), 2, "a non-full turn must retain both planar section caps");
        assert!(caps[0].cross(caps[1]).length() > angle.sin().abs() * 0.5,
            "distinct partial-turn cap planes were merged or made coincident");
    }
}

#[test]
fn negative_bezier_control_radius_does_not_reject_a_strictly_positive_meridian() {
    // z(t)=3t, r(t)=1-3.3t+9.6t^2-5.3t^3. Although control 1 has r=-0.1,
    // the curve minimum is r(11/53)=9568/14045 > 0; the other derivative root
    // is t=1. A control hull that crosses the axis is not proof the curve does.
    // Exact integrals: integral r^2 dz=18069/3500, integral r dz=147/40.
    let loops = vec![vec![vec![3., 1., 0., -0.1, 1., 2., 2.],
        line(2., 3.), line(0., 3.), line(0., 0.)]];
    let lateral = TAU * simpson(|t| {
        let r = 1. - 3.3 * t + 9.6 * t * t - 5.3 * t * t * t;
        r * (-3.3 + 19.2 * t - 15.9 * t * t).hypot(3.)
    }, 0., 1.);
    for angle in [PI / 2., TAU] {
        verify(&loops, angle, frame(), 18069. * PI / 3500. * angle / TAU,
            area_for_angle(lateral + 5. * PI, 147. / 40., angle), false, None);
    }
}

#[test]
fn sheared_conic_meridian_preserves_independent_conjugate_vectors() {
    // (r,z)=(6,0)+(1.2,0.4)*cos(t)+(0.3,0.9)*sin(t).
    // The conjugate vectors are independent but not perpendicular. The section
    // area is pi*det(u,v)=0.96*pi. Central symmetry places both its area centroid
    // and its perimeter centroid at radius 6, so Pappus gives both full measures.
    let (big, ur, uz, vr, vz) = (6_f64, 1.2_f64, 0.4_f64, 0.3_f64, 0.9_f64);
    let determinant = ur * vz - uz * vr;
    near(determinant, 0.96);
    let loops = vec![vec![
        vec![2., big + ur, uz, big, 0., ur, uz, vr, vz, 1.],
        vec![2., big - ur, -uz, big, 0., ur, uz, vr, vz, 1.],
    ]];
    let circumference = simpson(|t| {
        (-ur * t.sin() + vr * t.cos()).hypot(-uz * t.sin() + vz * t.cos())
    }, 0., TAU);
    for angle in [TAU, PI / 2., -PI / 2.] {
        verify(&loops, angle, frame(), 2. * PI * PI * big * determinant * angle.abs() / TAU,
            area_for_angle(TAU * big * circumference, PI * determinant, angle), false, None);
    }
}

#[test]
fn cubic_meridian_with_two_axis_endpoints_has_closed_manifold_poles() {
    // Cubic controls (0,0),(2,0),(2,3),(0,3): r=6t(1-t), z=9t^2-6t^3.
    // z'=18t(1-t); integral r^2 dz=648*Beta(4,4)=162/35 and
    // integral r dz=108*Beta(3,3)=18/5. Both poles close without circular caps.
    let loops = vec![vec![vec![3., 0., 0., 2., 0., 2., 3.], line(0., 3.)]];
    let lateral = TAU * simpson(|t| {
        let r = 6. * t * (1. - t);
        let dr = 6. - 12. * t;
        let dz = 18. * t * (1. - t);
        r * dr.hypot(dz)
    }, 0., 1.);
    for angle in [TAU, PI / 2., -PI / 2.] {
        verify(&loops, angle, frame(), 162. * PI / 35. * angle.abs() / TAU,
            area_for_angle(lateral, 18. / 5., angle), false, None);
    }
}

#[test]
fn malformed_and_nonmanifold_general_sections_fail_without_a_body() {
    let valid = vec![rectangle(1., 2., 0., 3.)];
    for angle in [0., f64::NAN, f64::INFINITY, -f64::INFINITY, TAU + 0.001] {
        assert!(protocol::handle(&request(&valid, angle, frame())).is_err());
    }
    let invalid = [
        vec![vec![line(1., 0.), line(3., 3.), line(3., 0.), line(1., 3.)]],
        vec![vec![line(0., 0.), line(2., -1.), line(2., 1.)]], // isolated axis pinch
        vec![vec![arc(2., 0., 1., 0., 1.), arc(0., 0., 1., 0., 1.)]], // full horn torus
        vec![vec![arc(2., 0., 0., 0., 1.), line(0., 3.)]], // mismatched arc endpoint
        vec![vec![vec![2., 2., 0., 0., 0., 2., 0., 1., 0., 1.], line(-2., 0.)]],
        vec![vec![vec![3., 1., 0., -5., 1., -5., 2.], line(1., 3.), line(0., 3.), line(0., 0.)]],
        vec![rectangle(1., 4., 0., 5.), rectangle(3., 5., 1., 4.)], // hole crosses outer
        vec![rectangle(1., 4., 0., 5.), rectangle(1., 2., 1., 4.)], // hole touches outer
    ];
    for loops in invalid {
        assert!(protocol::handle(&request(&loops, TAU, frame())).is_err(), "invalid profile admitted: {loops:?}");
    }
    let mut short = request(&valid, PI, frame()); short.data.pop();
    assert!(protocol::handle(&short).is_err());
    let mut trailing = request(&valid, PI, frame()); trailing.data.push(0.);
    assert!(protocol::handle(&trailing).is_err());
}

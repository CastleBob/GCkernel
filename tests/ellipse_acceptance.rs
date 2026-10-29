//! Independent acceptance for Bob's exact-ellipse milestone (Codex, 2026-10-29).
//! Public analytic mathematics only; see GCm/Handoff/GCKernel-ResearchWatch.md.
//! Prospective tests: the existing Conic profile must become an analytic solid,
//! never a polygonal substitute. No worker executable is rebuilt by these tests.

use std::collections::HashMap;
use std::f64::consts::{PI, TAU};
use gc_kernel::geom::{Curve, Surface, Vec3};
use gc_kernel::neutral::{Body, Profile, Reader, write_body};
use gc_kernel::{boolean, extrude, mass, tessellate, topology};

const TOL: f64 = 1e-7;
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn ellipse(c: Vec3, x: Vec3, y: Vec3, a: f64, b: f64, low: f64, high: f64) -> Profile {
    Profile { n: x.cross(y), low, high, loops: vec![vec![Curve::Conic {
        c, u: x.scale(a), v: y.scale(b), t0: 0.0, t1: TAU,
    }]] }
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 1e-9 * expected.abs().max(1.0),
        "actual {actual}, analytic {expected}");
}

fn checked_mesh(body: &Body, chord: f64) -> tessellate::Mesh {
    topology::check(body, TOL).unwrap();
    let original = write_body(body);
    let mut reader = Reader::with_tolerance(&original, TOL);
    let restored = reader.body().unwrap();
    reader.finish().unwrap();
    assert_eq!(write_body(&restored), original);
    let mesh = tessellate::tessellate(&restored, chord).unwrap();
    assert_eq!(write_body(body), original, "display changed exact geometry");
    assert!(mesh.positions.iter().chain(&mesh.normals).all(|p| p.is_finite()));
    assert!(!mesh.triangles.is_empty());
    // An edge used twice by one face is its periodic seam, not a selectable part edge.
    for (edge, _) in body.edges.iter().enumerate() {
        let owners: Vec<_> = body.faces.iter().enumerate().flat_map(|(face,f)|
            f.loops.iter().flatten().filter(move |&&(e,_)| e == edge).map(move |_| face)).collect();
        if owners.len() == 2 && owners[0] == owners[1] {
            assert!(!mesh.display_edges.contains(&edge), "periodic seam leaked into display edges");
        }
    }
    let bits = |n: f64| if n == 0.0 { 0 } else { n.to_bits() };
    let key = |p: Vec3| (bits(p.x), bits(p.y), bits(p.z));
    let mut edges = HashMap::<_, i32>::new();
    let origin = mesh.positions[0];
    let mut volume = 0.0;
    for tri in &mesh.triangles {
        let [a, b, c] = tri.map(|i| mesh.positions[i]);
        for (p, q) in [(a,b), (b,c), (c,a)] {
            if key(p) == key(q) { continue; }
            *edges.entry((key(p), key(q))).or_default() += 1;
            *edges.entry((key(q), key(p))).or_default() -= 1;
        }
        let cross = b.sub(a).cross(c.sub(a));
        if cross.length() > 1e-14 {
            let normal = tri.iter().fold(Vec3::new(0.0,0.0,0.0), |s, &i| s.add(mesh.normals[i]));
            assert!(cross.dot(normal) > 0.0, "triangle winds against its outward normals");
        }
        volume += a.sub(origin).dot(b.sub(origin).cross(c.sub(origin))) / 6.0;
    }
    assert!(edges.values().all(|&n| n == 0), "display mesh has an open or mismatched boundary");
    let measured = mass::measure(body, TOL).unwrap();
    assert!((volume - measured.volume).abs() <= chord * measured.area + 1e-9,
        "mesh volume {volume}, analytic volume {}", measured.volume);
    mesh
}

#[test]
fn exact_eccentric_prisms_roundtrip_and_display_in_rotated_planes() {
    for (a,b) in [(1.5,0.5), (0.5,1.5), (0.75,0.75), (2.0,0.1)] {
        for degrees in [0.0_f64, 37.0, 90.0] {
            for tilted in [false, true] {
                let angle = degrees.to_radians();
                let (base_x, base_y) = if tilted {
                    (Vec3::new(0.0,0.0,1.0), Vec3::new(1.0,0.0,0.0))
                } else { (Vec3::new(1.0,0.0,0.0), Vec3::new(0.0,1.0,0.0)) };
                let x = base_x.scale(angle.cos()).add(base_y.scale(angle.sin()));
                let y = base_y.scale(angle.cos()).sub(base_x.scale(angle.sin()));
                let c = Vec3::new(17.0,-23.0,5.0);
                let p = ellipse(c,x,y,a,b,-0.7,0.3);
                let body = extrude::extrude(&p,TOL).expect("exact ellipse extrusion is required");
                close(mass::measure(&body,TOL).unwrap().volume, PI*a*b);
                if a != b {
                    assert!(body.edges.iter().any(|e| matches!(e.curve, Curve::Conic { .. })),
                        "ellipse was replaced by non-elliptical edges");
                }
                let chord = 0.001;
                let mesh = checked_mesh(&body,chord);
                for (tri, &face) in mesh.triangles.iter().zip(&mesh.faces) {
                    if matches!(body.faces[face].surface, Surface::Plane { .. }) { continue; }
                    for &i in tri {
                        let d = mesh.positions[i].sub(c);
                        let (u,v,z) = (d.dot(x)/a,d.dot(y)/b,d.dot(p.n));
                        close(u*u+v*v,1.0);
                        assert!(z >= -0.7-TOL && z <= 0.3+TOL);
                        let outward = x.scale(u/a).add(y.scale(v/b)).unit().unwrap();
                        assert!(mesh.normals[i].dot(outward) > 1.0-1e-9);
                    }
                    for (i,j) in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])] {
                        let d = mesh.positions[i].sub(c);
                        let e = mesh.positions[j].sub(c);
                        let start = (d.dot(y)/b).atan2(d.dot(x)/a);
                        let end = (e.dot(y)/b).atan2(e.dot(x)/a);
                        let sweep = (end-start+PI).rem_euclid(TAU)-PI;
                        let t = start+sweep/2.0;
                        let z = (d.dot(p.n)+e.dot(p.n))/2.0;
                        let exact = c.add(x.scale(a*t.cos())).add(y.scale(b*t.sin())).add(p.n.scale(z));
                        let mid = mesh.positions[i].add(mesh.positions[j]).scale(0.5);
                        assert!(exact.distance(mid) <= chord+1e-9,
                            "ellipse side edge exceeds display chord: {}", exact.distance(mid));
                    }
                }
            }
        }
    }
}

#[test]
fn repeated_rotated_elliptical_through_cuts_have_analytic_removed_volume() {
    let axes = [Vec3::new(1.0,0.0,0.0),Vec3::new(0.0,1.0,0.0),Z];
    for axis in 0..3 {
        let (base_x,base_y) = (axes[axis],axes[(axis+1)%3]);
        let n = base_x.cross(base_y);
        let origin = Vec3::new(17.0,-23.0,5.0);
        let at = |x:f64,y:f64| origin.add(base_x.scale(x)).add(base_y.scale(y));
        let pts = [at(0.0,0.0),at(10.0,0.0),at(10.0,8.0),at(0.0,8.0)];
        let lines = (0..4).map(|i| Curve::Line { a:pts[i], b:pts[(i+1)%4] }).collect();
        let mut body = extrude::extrude(&Profile { loops:vec![lines],n,low:0.0,high:0.5 },TOL).unwrap();
        let mut expected = 40.0;
        for (cx,cy,a,b,degrees) in [(2.0,2.0,1.2,0.4,37.0_f64),(7.0,5.0,0.8,0.25,-23.0)] {
            let t = degrees.to_radians();
            let x = base_x.scale(t.cos()).add(base_y.scale(t.sin()));
            let y = n.cross(x);
            let tool = extrude::extrude(&ellipse(at(cx,cy),x,y,a,b,-0.1,0.6),TOL).unwrap();
            body = boolean::cut(&body,&tool,TOL).expect("an interior elliptical through-cut must complete");
            expected -= PI*a*b*0.5;
            close(mass::measure(&body,TOL).unwrap().volume,expected);
            checked_mesh(&body,0.001);
        }
    }
}

#[test]
fn elliptic_through_cut_at_plate_edge_keeps_a_separate_round_hole() {
    let pts = [Vec3::new(0.0,0.0,0.0),Vec3::new(10.0,0.0,0.0),
        Vec3::new(10.0,8.0,0.0),Vec3::new(0.0,8.0,0.0)];
    let lines = (0..4).map(|i| Curve::Line { a:pts[i], b:pts[(i+1)%4] }).collect();
    let plate = extrude::extrude(&Profile { loops:vec![lines],n:Z,low:0.0,high:0.5 },TOL).unwrap();
    let hole = extrude::extrude(&Profile { n:Z,low:-0.1,high:0.6,loops:vec![vec![Curve::Arc {
        c:Vec3::new(7.0,5.0,0.0),n:Z,x:Vec3::new(1.0,0.0,0.0),r:0.5,t0:0.0,t1:TAU,
    }]] },TOL).unwrap();
    let body = boolean::cut(&plate,&hole,TOL).unwrap();
    // The plate's x=0 boundary bisects the ellipse, so precisely half its area is removed.
    let cutter = extrude::extrude(&ellipse(Vec3::new(0.0,4.0,0.0),Vec3::new(1.0,0.0,0.0),
        Vec3::new(0.0,1.0,0.0),1.2,0.4,-0.1,0.6),TOL).unwrap();
    let cut = boolean::cut(&body,&cutter,TOL).expect("line-conic edge cut with disjoint round hole must complete");
    close(mass::measure(&cut,TOL).unwrap().volume,40.0-PI*0.25*0.5-0.5*PI*1.2*0.4*0.5);
    checked_mesh(&cut,0.001);
}

#[test]
fn singular_and_out_of_plane_ellipse_profiles_are_rejected() {
    let x = Vec3::new(1.0,0.0,0.0);
    for v in [x,Vec3::new(0.0,0.0,1.0)] {
        let bad = Profile { n:Z,low:0.0,high:1.0,loops:vec![vec![Curve::Conic {
            c:Vec3::new(0.0,0.0,0.0),u:x,v,t0:0.0,t1:TAU,
        }]] };
        assert!(extrude::extrude(&bad,TOL).is_err(), "invalid ellipse profile was accepted");
    }
}

#[test]
fn admitted_elliptic_wire_has_finite_stable_parameters() {
    // Near-parallel conjugates expose cancellation in a Gram-determinant inverse.
    // Explicit rejection is allowed; admitting a surface with NaN parameters is not.
    for height in [3.0, 1e-8] {
        let wire = [3.0, 0.0,0.0,0.0, 0.0,0.0,1.0, 1.0,0.0,0.0, 1.0,height,0.0];
        let mut reader = Reader::with_tolerance(&wire,TOL);
        let parsed = reader.surface();
        if height == 3.0 { assert!(parsed.is_ok()); }
        if let Ok(Surface::EllipticCylinder { p,a,u,v }) = parsed {
            reader.finish().unwrap();
            for t in [0.0_f64,0.7,2.3] {
                let q = p.add(u.scale(t.cos())).add(v.scale(t.sin())).add(a.scale(2.0));
                let (angle, axial) = Surface::elliptic_params(p,a,u,v,q);
                assert!(angle.is_finite() && axial.is_finite());
                close(angle,t);
                close(axial,2.0);
            }
        }
    }
    for (u,v) in [(Vec3::new(1.0,0.0,0.0),Vec3::new(2.0,0.0,0.0)),
        (Vec3::new(1e308,0.0,0.0),Vec3::new(0.0,1e308,0.0)),
        (Vec3::new(1e308,1e308,0.0),Vec3::new(1e308,1e308,0.0))] {
        let wire = [3.0,0.0,0.0,0.0,0.0,0.0,1.0,u.x,u.y,u.z,v.x,v.y,v.z];
        assert!(Reader::with_tolerance(&wire,TOL).surface().is_err(),
            "singular or overflowing elliptic wire was admitted");
    }
}

#[test]
fn crossing_ellipses_sharing_four_quarter_points_are_not_treated_as_coincident() {
    let x = Vec3::new(1.0,0.0,0.0);
    let y = Vec3::new(0.0,1.0,0.0);
    let c = Vec3::new(5.0,4.0,0.0);
    let stock = extrude::extrude(&ellipse(c,x,y,1.5,0.5,0.0,0.5),TOL).unwrap();
    // Both cutter conjugates land on the stock ellipse, but they are not conjugate
    // in its unit-circle frame. Four quarter-point samples cannot prove coincidence.
    let cutter = extrude::extrude(&Profile { n:Z,low:-0.1,high:0.6,loops:vec![vec![Curve::Conic {
        c,u:Vec3::new(1.5,0.0,0.0),v:Vec3::new(0.75,3.0_f64.sqrt()/4.0,0.0),t0:0.0,t1:TAU,
    }]] },TOL).unwrap();
    let error = boolean::cut(&stock,&cutter,TOL).err()
        .expect("unsupported crossing ellipses must be refused, not treated as one support");
    assert!(format!("{error:?}").contains("ellipse meeting"),
        "crossing supports must stop at admission, not later invalid geometry: {error:?}");
}

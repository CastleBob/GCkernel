//! Provenance tests (Claude 2026-10-30): every record is proved against the geometry it names, every request segment is
//! named exactly once per role, and the named ops' replies are the plain ops' replies followed by a well-formed envelope.

use std::f64::consts::PI;

use gc_kernel::extrude::{extrude, extrude_named};
use gc_kernel::geom::{Curve, Surface, Vec3};
use gc_kernel::neutral::{write_profile, Body, Profile};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_EXTRUDE, OP_EXTRUDE_NAMED};
use gc_kernel::provenance::*;

const TOL: f64 = 1e-7;
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn p(x: f64, y: f64) -> Vec3 { Vec3::new(x, y, 0.0) }
fn line(a: (f64, f64), b: (f64, f64)) -> Curve { Curve::Line { a: p(a.0, a.1), b: p(b.0, b.1) } }

/// Outer listed second and clockwise (with a clockwise arc), hole listed first and counter-clockwise: the extruder swaps
/// and reverses both, so request indices differ from built ones everywhere.
fn reordered() -> Profile {
    let hole = vec![line((1.0, 0.5), (3.0, 0.5)), line((3.0, 0.5), (3.0, 1.5)), line((3.0, 1.5), (1.0, 1.5)), line((1.0, 1.5), (1.0, 0.5))];
    let arc = Curve::Arc { c: p(2.0, 2.0), n: Z.scale(-1.0), x: Vec3::new(-1.0, 0.0, 0.0), r: 2.0, t0: 0.0, t1: PI };
    let outer = vec![line((0.0, 0.0), (0.0, 2.0)), arc, line((4.0, 2.0), (4.0, 0.0)), line((4.0, 0.0), (0.0, 0.0))];
    Profile { loops: vec![hole, outer], n: Z, low: -1.0, high: 3.0 }
}

fn plain() -> Profile {
    Profile {
        loops: vec![vec![line((0.0, 0.0), (5.0, 0.0)), line((5.0, 0.0), (5.0, 2.0)), line((5.0, 2.0), (0.0, 2.0)), line((0.0, 2.0), (0.0, 0.0))]],
        n: Z,
        low: 0.0,
        high: 1.5,
    }
}

fn near(a: Vec3, b: Vec3) -> bool { a.distance(b) <= TOL }

fn mid(c: &Curve) -> Vec3 {
    let (t0, t1) = c.domain();
    c.eval((t0 + t1) / 2.0).0
}

/// The edge is the request segment moved by d, in either direction.
fn same_curve(edge: &Curve, seg: &Curve, d: Vec3) -> bool {
    let (s, e) = (seg.start().add(d), seg.end().add(d));
    ((near(edge.start(), s) && near(edge.end(), e)) || (near(edge.start(), e) && near(edge.end(), s))) && near(mid(edge), mid(seg).add(d))
}

fn prove(profile: &Profile, body: &Body, prov: &Provenance) {
    assert_eq!((prov.faces.len(), prov.edges.len(), prov.vertices.len()), (body.faces.len(), body.edges.len(), body.vertices.len()));
    let (low, high) = (profile.n.scale(profile.low), profile.n.scale(profile.high));
    let seg = |r: &Record| &profile.loops[r.refs[0]][r.refs[1]];
    let edge_of = |role: u32, refs: &[usize]| prov.edges.iter().position(|r| r.role == role && r.refs == refs).expect("edge named");
    for (v, r) in prov.vertices.iter().enumerate() {
        assert_eq!(r.refs.len(), 2);
        let d = match r.role { LOW_POINT => low, HIGH_POINT => high, other => panic!("vertex role {other}") };
        assert!(near(body.vertices[v], seg(r).start().add(d)), "vertex {v} {r:?}");
    }
    for (e, r) in prov.edges.iter().enumerate() {
        assert_eq!(r.refs.len(), 2);
        let curve = &body.edges[e].curve;
        match r.role {
            LOW_SEGMENT => assert!(same_curve(curve, seg(r), low), "edge {e} {r:?}"),
            HIGH_SEGMENT => assert!(same_curve(curve, seg(r), high), "edge {e} {r:?}"),
            RULING => assert!(matches!(curve, Curve::Line { .. }) && near(curve.start(), seg(r).start().add(low)) && near(curve.end(), seg(r).start().add(high)), "edge {e} {r:?}"),
            other => panic!("edge role {other}"),
        }
    }
    for (f, r) in prov.faces.iter().enumerate() {
        let face = &body.faces[f];
        let edges: Vec<usize> = face.loops.iter().flatten().map(|&(e, _)| e).collect();
        match r.role {
            LOW_CAP | HIGH_CAP => {
                assert!(r.refs.is_empty());
                let Surface::Plane { n, .. } = face.surface else { panic!("cap {f} is not a plane") };
                let outward = if face.sense { n } else { n.scale(-1.0) };
                let want = if r.role == HIGH_CAP { profile.n } else { profile.n.scale(-1.0) };
                assert!(near(outward, want), "cap {f} {r:?}");
                let role = if r.role == HIGH_CAP { HIGH_SEGMENT } else { LOW_SEGMENT };
                assert!(edges.iter().all(|&e| prov.edges[e].role == role), "cap {f} bounded by its own segment copies");
            }
            SIDE => {
                assert_eq!(r.refs.len(), 2);
                assert!(edges.contains(&edge_of(LOW_SEGMENT, &r.refs)) && edges.contains(&edge_of(HIGH_SEGMENT, &r.refs)), "side {f} {r:?}");
            }
            other => panic!("face role {other}"),
        }
    }
    // Each request segment exactly once per role.
    for (l, lp) in profile.loops.iter().enumerate() {
        for s in 0..lp.len() {
            for (role, list) in [(SIDE, &prov.faces), (LOW_SEGMENT, &prov.edges), (HIGH_SEGMENT, &prov.edges), (RULING, &prov.edges),
                                 (LOW_POINT, &prov.vertices), (HIGH_POINT, &prov.vertices)] {
                assert_eq!(list.iter().filter(|r| r.role == role && r.refs == [l, s]).count(), 1, "role {role} loop {l} seg {s}");
            }
        }
    }
    for role in [LOW_CAP, HIGH_CAP] { assert_eq!(prov.faces.iter().filter(|r| r.role == role).count(), 1); }
}

#[test]
fn extrude_provenance_names_request_segments_through_reorder_and_reversal() {
    for profile in [reordered(), plain()] {
        let (body, prov) = extrude_named(&profile, TOL).unwrap();
        prove(&profile, &body, &prov);
    }
}

#[test]
fn named_extrude_body_is_the_plain_extrude_body() {
    for profile in [reordered(), plain()] {
        let (named, _) = extrude_named(&profile, TOL).unwrap();
        let plain = extrude(&profile, TOL).unwrap();
        assert_eq!(gc_kernel::neutral::write_body(&named), gc_kernel::neutral::write_body(&plain));
    }
}

/// The op 10 reply = the op 2 reply, then exactly one envelope that decodes to the in-process provenance.
#[test]
fn op10_reply_is_op2_reply_plus_envelope() {
    for profile in [reordered(), plain()] {
        let data = write_profile(&profile);
        let run = |op| handle(&parse_request(&encode_request(op, TOL, &data)).unwrap()).unwrap();
        let (two, ten) = (run(OP_EXTRUDE), run(OP_EXTRUDE_NAMED));
        assert_eq!(ten[..two.len()], two[..]);
        let env = &ten[two.len()..];
        assert_eq!((env[0], env[1]), (START, VERSION));
        let len = env[2] as usize;
        assert_eq!(env.len(), 3 + len + 1, "nothing after the end marker");
        assert_eq!(env[3 + len], END);
        let body = &env[3..3 + len];
        let (nf, ne, nv) = (body[0] as usize, body[1] as usize, body[2] as usize);
        let mut at = 3;
        let mut read = |n: usize| -> Vec<Record> {
            (0..n).map(|_| {
                let (role, k) = (body[at] as u32, body[at + 1] as usize);
                let r = Record { role, refs: body[at + 2..at + 2 + k].iter().map(|&v| v as usize).collect() };
                at += 2 + k;
                r
            }).collect()
        };
        let decoded = Provenance { faces: read(nf), edges: read(ne), vertices: read(nv) };
        assert_eq!(at, len);
        assert_eq!(decoded, extrude_named(&profile, TOL).unwrap().1);
    }
}

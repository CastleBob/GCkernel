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

// ---- op 11 named revolve ----

use gc_kernel::protocol::{OP_REVOLVE_GENERAL, OP_REVOLVE_NAMED};
use gc_kernel::revolve::{revolve_general, revolve_general_named, Kind, Revolution, Segment};

const TAU_: f64 = 2.0 * PI;

fn sg(r: f64, z: f64) -> Segment { Segment { r, z, kind: Kind::Line } }
fn sa(r: f64, z: f64, cr: f64, cz: f64, ccw: bool) -> Segment { Segment { r, z, kind: Kind::Arc { cr, cz, ccw } } }

fn frame(angle: f64, loops: Vec<Vec<Segment>>) -> Revolution {
    Revolution { p: Vec3::new(1.0, -2.0, 0.5), a: Vec3::new(0.0, 0.6, 0.8), x: Vec3::new(1.0, 0.0, 0.0), angle, loops }
}

/// Hollow ring: outer rectangle listed second and clockwise, a circular hole listed first counter-clockwise (both
/// swapped and reversed); a ball section given clockwise (reversed) on the axis, so faces split at the middle; a cup
/// with a flat disc bottom; full turns, small and large part turns, negative turns (mirrored).
fn revolutions() -> Vec<Revolution> {
    let ring = || vec![
        vec![sa(3.0, 1.5, 2.5, 1.5, true), sa(2.0, 1.5, 2.5, 1.5, true)],
        vec![sg(1.0, 0.0), sg(1.0, 3.0), sg(4.0, 3.0), sg(4.0, 0.0)],
    ];
    let ball = || vec![vec![sa(0.0, 1.0, 0.0, 0.0, false), sg(0.0, -1.0)]];
    let cup = || vec![vec![sg(0.0, 0.0), sg(2.0, 0.0), sg(2.0, 3.0), sg(1.5, 3.0), sg(1.5, 0.5), sg(0.0, 0.5)]];
    vec![frame(TAU_, ring()), frame(1.0, ring()), frame(-4.7, ring()), frame(TAU_, ball()), frame(-4.7, ball()), frame(2.0, ball()),
         frame(TAU_, cup()), frame(4.0, cup()), frame(-TAU_, cup())]
}

fn prove_revolve(rev: &Revolution, body: &Body, prov: &Provenance) {
    assert_eq!((prov.faces.len(), prov.edges.len(), prov.vertices.len()), (body.faces.len(), body.edges.len(), body.vertices.len()));
    let full = (rev.angle.abs() - TAU_).abs() <= 1e-12;
    let y = rev.a.cross(rev.x);
    let split_any = full || rev.angle.abs() > PI;
    let middle = if full { PI * rev.angle.signum() } else { rev.angle / 2.0 };
    // Request (radius, height) at signed angle t.
    let rho = |t: f64| rev.x.scale(t.cos()).add(y.scale(t.sin()));
    let at = |r: f64, z: f64, t: f64| rev.p.add(rev.a.scale(z)).add(rho(t).scale(r));
    let start = |l: usize, s: usize| { let q = rev.loops[l][s]; (q.r, q.z) };
    let next = |l: usize, s: usize| start(l, (s + 1) % rev.loops[l].len());
    let in_meridian = |q: Vec3, t: f64| { let w = q.sub(rev.p); w.dot(rev.a.cross(rho(t))).abs() <= TOL && w.dot(rho(t)) >= -TOL };
    let part_angles = |part: usize| match part {
        0 => if full { (0.0, 0.0) } else { (0.0, rev.angle) },
        1 => (0.0, middle),
        2 => (middle, if full { 0.0 } else { rev.angle }),
        _ => panic!("part {part}"),
    };
    for (v, r) in prov.vertices.iter().enumerate() {
        assert_eq!(r.refs.len(), 2);
        let (pr, pz) = start(r.refs[0], r.refs[1]);
        let want = match r.role {
            AXIS_POINT => { assert_eq!(pr, 0.0); at(0.0, pz, 0.0) }
            START_POINT => at(pr, pz, 0.0),
            END_POINT => { assert!(!full); at(pr, pz, rev.angle) }
            MID_POINT => { assert!(split_any); at(pr, pz, middle) }
            other => panic!("vertex role {other}"),
        };
        assert!(near(body.vertices[v], want), "vertex {v} {r:?}");
    }
    for (e, r) in prov.edges.iter().enumerate() {
        let edge = &body.edges[e];
        let (vs, ve) = (body.vertices[edge.start], body.vertices[edge.end]);
        let ends = |a: Vec3, b: Vec3| (near(vs, a) && near(ve, b)) || (near(vs, b) && near(ve, a));
        match r.role {
            PARALLEL => {
                assert_eq!(r.refs.len(), 3);
                let (pr, pz) = start(r.refs[0], r.refs[1]);
                let (t0, t1) = part_angles(r.refs[2]);
                assert!(matches!(edge.curve, Curve::Arc { r: ar, .. } if (ar - pr).abs() <= TOL), "edge {e} {r:?}");
                assert!(ends(at(pr, pz, t0), at(pr, pz, t1)), "edge {e} {r:?}");
                let tm = match (r.refs[2], full) { (0, true) => PI, (2, true) => 1.5 * middle, _ => (t0 + t1) / 2.0 };
                assert!(near(mid(&edge.curve), at(pr, pz, tm)), "edge {e} {r:?} midpoint");
            }
            START_MERIDIAN | END_MERIDIAN | MID_MERIDIAN => {
                assert_eq!(r.refs.len(), 2);
                let t = match r.role { START_MERIDIAN => 0.0, END_MERIDIAN => { assert!(!full); rev.angle } _ => middle };
                let ((r0, z0), (r1, z1)) = (start(r.refs[0], r.refs[1]), next(r.refs[0], r.refs[1]));
                assert!(ends(at(r0, z0, t), at(r1, z1, t)), "edge {e} {r:?}");
                assert!(in_meridian(mid(&edge.curve), t), "edge {e} {r:?} lies in its meridian half-plane");
            }
            AXIS_SEGMENT => {
                assert!(!full && r.refs.len() == 2);
                let ((r0, z0), (r1, z1)) = (start(r.refs[0], r.refs[1]), next(r.refs[0], r.refs[1]));
                assert!(r0 == 0.0 && r1 == 0.0 && ends(at(0.0, z0, 0.0), at(0.0, z1, 0.0)), "edge {e} {r:?}");
            }
            other => panic!("edge role {other}"),
        }
    }
    for (f, r) in prov.faces.iter().enumerate() {
        let face = &body.faces[f];
        let edges: Vec<&Record> = face.loops.iter().flatten().map(|&(e, _)| &prov.edges[e]).collect();
        match r.role {
            SWEPT => {
                assert_eq!(r.refs.len(), 3);
                let (l, s, part) = (r.refs[0], r.refs[1], r.refs[2]);
                let n = rev.loops[l].len();
                // Bounded only by this segment's meridians and the parallels of its two end points, of this part.
                for er in &edges {
                    let ok = match er.role {
                        START_MERIDIAN | END_MERIDIAN | MID_MERIDIAN => er.refs == [l, s],
                        PARALLEL => er.refs[0] == l && (er.refs[1] == s || er.refs[1] == (s + 1) % n) && (er.refs[2] == part || full && part == 0),
                        _ => false,
                    };
                    assert!(ok, "face {f} {r:?} has edge {er:?}");
                }
                let (t0, t1) = part_angles(part);
                let has = |role: u32| edges.iter().any(|er| er.role == role && er.refs == [l, s]);
                let role_at = |t: f64| if t == 0.0 { START_MERIDIAN } else if t == middle { MID_MERIDIAN } else { END_MERIDIAN };
                if edges.iter().any(|er| er.role != PARALLEL) {
                    assert!(has(role_at(t0)) && has(role_at(t1)), "face {f} {r:?} bounded by its meridians at {t0}, {t1}");
                } else {
                    assert!(full && part == 0, "a face of parallels only is a full-turn annulus or disc");
                }
            }
            START_CAP | END_CAP => {
                assert!(!full && r.refs.is_empty());
                let Surface::Plane { n, .. } = face.surface else { panic!("cap {f} is not a plane") };
                let outward = if face.sense { n } else { n.scale(-1.0) };
                let t = if r.role == START_CAP { 0.0 } else { rev.angle };
                // The solid lies on the inner side of the turn: outward = against the turn at the start, with it at the end.
                let turn = rev.a.cross(rho(t)).scale(rev.angle.signum());
                let want = if r.role == START_CAP { turn.scale(-1.0) } else { turn };
                assert!(near(outward, want), "cap {f} {r:?} outward {outward:?} want {want:?}");
                assert!(edges.iter().all(|er| matches!((r.role, er.role), (START_CAP, START_MERIDIAN) | (END_CAP, END_MERIDIAN) | (_, AXIS_SEGMENT))),
                        "cap {f} bounded by its meridians and axis segments");
            }
            other => panic!("face role {other}"),
        }
    }
    // No record twice; every request segment not on the axis sweeps a face; every request point is named.
    for list in [&prov.faces, &prov.edges, &prov.vertices] {
        for (i, a) in list.iter().enumerate() { assert!(!list[i + 1..].contains(a), "duplicate record {a:?}"); }
    }
    for (l, lp) in rev.loops.iter().enumerate() {
        for s in 0..lp.len() {
            let on_axis = start(l, s).0 == 0.0 && next(l, s).0 == 0.0 && matches!(lp[s].kind, Kind::Line);
            assert_eq!(prov.faces.iter().any(|r| r.role == SWEPT && r.refs[..2] == [l, s]), !on_axis, "segment {l} {s}");
            assert!(prov.vertices.iter().any(|r| r.refs == [l, s]), "point {l} {s}");
        }
    }
}

#[test]
fn revolve_provenance_names_request_segments_through_reorder_reversal_mirroring_and_splits() {
    for rev in revolutions() {
        let (body, prov) = revolve_general_named(&rev, TOL).unwrap();
        prove_revolve(&rev, &body, &prov);
        assert_eq!(gc_kernel::neutral::write_body(&body), gc_kernel::neutral::write_body(&revolve_general(&rev, TOL).unwrap()));
    }
}

fn write_revolution(rev: &Revolution) -> Vec<f64> {
    let mut d = vec![rev.p.x, rev.p.y, rev.p.z, rev.a.x, rev.a.y, rev.a.z, rev.x.x, rev.x.y, rev.x.z, rev.angle, rev.loops.len() as f64];
    for lp in &rev.loops {
        d.push(lp.len() as f64);
        for s in lp {
            match s.kind {
                Kind::Line => d.extend([0.0, s.r, s.z]),
                Kind::Arc { cr, cz, ccw } => d.extend([1.0, s.r, s.z, cr, cz, if ccw { 1.0 } else { -1.0 }]),
                _ => unreachable!(),
            }
        }
    }
    d
}

#[test]
fn op11_reply_is_op9_reply_plus_envelope() {
    for rev in revolutions() {
        let data = write_revolution(&rev);
        let run = |op| handle(&parse_request(&encode_request(op, TOL, &data)).unwrap()).unwrap();
        let (nine, eleven) = (run(OP_REVOLVE_GENERAL), run(OP_REVOLVE_NAMED));
        assert_eq!(eleven[..nine.len()], nine[..]);
        let mut want = Vec::new();
        write(&mut want, &revolve_general_named(&rev, TOL).unwrap().1);
        assert_eq!(eleven[nine.len()..], want[..]);
    }
}

// ---- op 12 named cut ----

use gc_kernel::protocol::{OP_CUT, OP_CUT_NAMED};

fn plate() -> Body {
    let r = |x0: f64, y0: f64, x1: f64, y1: f64| vec![line((x0, y0), (x1, y0)), line((x1, y0), (x1, y1)), line((x1, y1), (x0, y1)), line((x0, y1), (x0, y0))];
    extrude(&Profile { loops: vec![r(0.0, 0.0, 6.0, 4.0)], n: Z, low: 0.0, high: 0.5 }, TOL).unwrap()
}

/// A drill extruded along z through the plate (spanning it), centre (x, y).
fn drill(x: f64, y: f64, r: f64) -> Body {
    let circle = Curve::Arc { c: p(x, y), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r, t0: 0.0, t1: 2.0 * PI };
    extrude(&Profile { loops: vec![vec![circle]], n: Z, low: -1.0, high: 1.5 }, TOL).unwrap()
}

fn cut_reply(op: u32, body: &Body, tool: &Body) -> Result<Vec<f64>, gc_kernel::KernelError> {
    let data: Vec<f64> = gc_kernel::neutral::write_body(body).into_iter().chain(gc_kernel::neutral::write_body(tool)).collect();
    handle(&parse_request(&encode_request(op, TOL, &data)).unwrap())
}

fn on_plane(s: &Surface, q: Vec3) -> bool { let Surface::Plane { p, n, .. } = *s else { panic!("not a plane") }; q.sub(p).dot(n).abs() <= TOL }
fn on_cylinder(s: &Surface, q: Vec3) -> bool {
    let Surface::Cylinder { p, a, r, .. } = *s else { panic!("not a cylinder") };
    let w = q.sub(p);
    (w.sub(a.scale(w.dot(a))).length() - r).abs() <= TOL
}

fn prove_cut(body: &Body, tool: &Body, result: &Body, prov: &Provenance, outcomes: &[Outcome]) {
    assert_eq!((prov.faces.len(), prov.edges.len(), prov.vertices.len()), (result.faces.len(), result.edges.len(), result.vertices.len()));
    let arity = |role: u32| match role { FROM_BODY_VERTEX | FROM_TOOL_VERTEX => 1, _ => 2 };
    for r in prov.faces.iter().chain(&prov.edges).chain(&prov.vertices) { assert_eq!(r.refs.len(), arity(r.role), "{r:?}"); }
    for (v, r) in prov.vertices.iter().enumerate() {
        let q = result.vertices[v];
        match r.role {
            FROM_BODY_VERTEX => assert_eq!(q, body.vertices[r.refs[0]]),
            SEAM_POINT => assert!(on_plane(&body.faces[r.refs[0]].surface, q) && on_cylinder(&tool.faces[r.refs[1]].surface, q), "vertex {v} {r:?}"),
            other => panic!("vertex role {other}"),
        }
    }
    for (e, r) in prov.edges.iter().enumerate() {
        let edge = &result.edges[e];
        let (s, t) = (&prov.vertices[edge.start], &prov.vertices[edge.end]);
        match r.role {
            FROM_BODY_EDGE => {
                let source = &body.edges[r.refs[0]];
                assert_eq!((r.refs[1], format!("{:?}", edge.curve)), (0, format!("{:?}", source.curve)), "edge {e} {r:?}");
                assert!(s.role == FROM_BODY_VERTEX && s.refs == [source.start] && t.role == FROM_BODY_VERTEX && t.refs == [source.end], "edge {e} ends");
            }
            MEETING => {
                let (plane, wall) = (&body.faces[r.refs[0]].surface, &tool.faces[r.refs[1]].surface);
                assert!(matches!(edge.curve, Curve::Arc { .. }), "edge {e} {r:?}");
                let (t0, t1) = edge.curve.domain();
                for k in 0..8 {
                    let q = edge.curve.eval(t0 + (t1 - t0) * k as f64 / 8.0).0;
                    assert!(on_plane(plane, q) && on_cylinder(wall, q), "edge {e} {r:?} sample {k}");
                }
                assert!(s.role == SEAM_POINT && s.refs == r.refs && edge.start == edge.end, "edge {e} is a circle from its seam point");
            }
            TOOL_SEAM => {
                let wall = &tool.faces[r.refs[1]].surface;
                assert!(matches!(edge.curve, Curve::Line { .. }) && on_cylinder(wall, mid(&edge.curve)), "edge {e} {r:?}");
                assert!(s.role == SEAM_POINT && s.refs == r.refs && t.role == SEAM_POINT && t.refs[1] == r.refs[1], "edge {e} ends");
            }
            other => panic!("edge role {other}"),
        }
    }
    for (f, r) in prov.faces.iter().enumerate() {
        let face = &result.faces[f];
        let edges: Vec<&Record> = face.loops.iter().flatten().map(|&(e, _)| &prov.edges[e]).collect();
        match r.role {
            FROM_BODY_FACE => {
                let source = &body.faces[r.refs[0]];
                assert_eq!((format!("{:?}", face.surface), face.sense), (format!("{:?}", source.surface), source.sense), "face {f} support");
                assert_eq!(face.loops[..source.loops.len()], source.loops[..], "face {f} keeps its source's loops");
                let added: Vec<&Record> = face.loops[source.loops.len()..].iter().flatten().map(|&(e, _)| &prov.edges[e]).collect();
                assert_eq!(r.refs[1] == 1, !added.is_empty(), "face {f} trimmed flag");
                assert!(added.iter().all(|er| er.role == MEETING && er.refs[0] == r.refs[0]), "face {f} new loops are its meeting edges");
            }
            FROM_TOOL_FACE => {
                assert_eq!(r.refs[1], 1);
                let wall = &tool.faces[r.refs[0]].surface;
                let (Surface::Cylinder { p, a, r: ra, .. }, Surface::Cylinder { p: q, a: b, r: rb, .. }) = (face.surface.clone(), wall.clone()) else { panic!("face {f}") };
                assert!(ra == rb && a.cross(b).length() <= TOL && on_cylinder(wall, p.add(face_x(&face.surface).scale(ra))) && q.sub(p).cross(b).length() <= TOL, "face {f} lies on tool face {r:?}");
                assert!(edges.iter().all(|er| matches!(er.role, MEETING | TOOL_SEAM) && er.refs[1] == r.refs[0]), "face {f} bounded by its meetings and seam");
            }
            other => panic!("face role {other}"),
        }
    }
    // Outcomes: each input entity's descendants, read independently from the records.
    let inputs = [body.faces.len(), body.edges.len(), body.vertices.len(), tool.faces.len(), tool.edges.len(), tool.vertices.len()];
    assert_eq!(outcomes.len(), inputs.iter().sum::<usize>());
    let kinds = [(&prov.faces, FROM_BODY_FACE), (&prov.edges, FROM_BODY_EDGE), (&prov.vertices, FROM_BODY_VERTEX),
                 (&prov.faces, FROM_TOOL_FACE), (&prov.edges, FROM_TOOL_EDGE), (&prov.vertices, FROM_TOOL_VERTEX)];
    let mut o = outcomes.iter();
    for ((records, role), &n) in kinds.into_iter().zip(&inputs) {
        for i in 0..n {
            let out = o.next().unwrap();
            for &k in &out.outputs { assert!(records[k].role == role && records[k].refs[0] == i); }
            assert_eq!(out.outputs.len(), records.iter().filter(|r| r.role == role && r.refs[0] == i).count());
        }
    }
}

fn face_x(s: &Surface) -> Vec3 { match *s { Surface::Cylinder { x, .. } => x, _ => panic!() } }

#[test]
fn bore_history_names_the_plate_and_the_drill() {
    let (body, tool) = (plate(), drill(2.0, 1.5, 0.5));
    let (result, prov, outcomes) = gc_kernel::boolean::cut_named(&body, &tool, TOL).unwrap();
    prove_cut(&body, &tool, &result, &prov, &outcomes);
    assert_eq!(gc_kernel::neutral::write_body(&result), gc_kernel::neutral::write_body(&gc_kernel::boolean::cut(&body, &tool, TOL).unwrap()));
    // The two plate faces the hole passes through are modified, the rest unchanged; the drill's wall becomes the hole's
    // wall, everything else of the drill is consumed.
    let wall = tool.faces.iter().position(|f| matches!(f.surface, Surface::Cylinder { .. })).unwrap();
    let hole = prov.faces.iter().position(|r| r.role == FROM_TOOL_FACE).unwrap();
    let (nf, ne, nv) = (body.faces.len(), body.edges.len(), body.vertices.len());
    for (i, out) in outcomes[..nf].iter().enumerate() {
        let through = matches!(body.faces[i].surface, Surface::Plane { n, .. } if n.cross(Z).length() <= TOL);
        assert_eq!((out.code, &out.outputs[..]), (if through { MODIFIED } else { UNCHANGED }, &[i][..]), "plate face {i}");
    }
    assert!(outcomes[nf..nf + ne + nv].iter().enumerate().all(|(i, o)| o.code == UNCHANGED && o.outputs == [if i < ne { i } else { i - ne }]));
    for (j, out) in outcomes[nf + ne + nv..nf + ne + nv + tool.faces.len()].iter().enumerate() {
        assert_eq!((out.code, out.outputs.clone()), if j == wall { (MODIFIED, vec![hole]) } else { (CONSUMED, vec![]) }, "drill face {j}");
    }
    assert!(outcomes[nf + ne + nv + tool.faces.len()..].iter().all(|o| o.code == CONSUMED && o.outputs.is_empty()));
}

#[test]
fn op12_reply_is_op3_reply_plus_envelope() {
    let (body, tool) = (plate(), drill(4.0, 2.5, 0.75));
    let (three, twelve) = (cut_reply(OP_CUT, &body, &tool).unwrap(), cut_reply(OP_CUT_NAMED, &body, &tool).unwrap());
    assert_eq!(twelve[..three.len()], three[..]);
    let (_, prov, outcomes) = gc_kernel::boolean::cut_named(&body, &tool, TOL).unwrap();
    let mut want = Vec::new();
    let inputs = [body.faces.len(), body.edges.len(), body.vertices.len(), tool.faces.len(), tool.edges.len(), tool.vertices.len()];
    write_cut(&mut want, &prov, inputs, &outcomes);
    assert_eq!(twelve[three.len()..], want[..]);
    let envelope = &twelve[three.len()..];
    assert_eq!((envelope[0], envelope[1], *envelope.last().unwrap()), (START, VERSION, END));
    assert_eq!(envelope[2] as usize, envelope.len() - 4);
}

/// A second hole runs in the section domain, which records no history yet: op 12 refuses by name, op 3 still cuts.
#[test]
fn unrecorded_domain_refuses_the_named_cut() {
    let once = gc_kernel::boolean::cut(&plate(), &drill(2.0, 1.5, 0.5), TOL).unwrap();
    let second = drill(4.5, 2.5, 0.5);
    match gc_kernel::boolean::cut_named(&once, &second, TOL) {
        Err(gc_kernel::KernelError::Unsupported(m)) => assert!(m.contains("Section"), "{m}"),
        other => panic!("expected a refusal, got {:?}", other.map(|_| ())),
    }
    assert!(cut_reply(OP_CUT_NAMED, &once, &second).is_err());
    assert!(cut_reply(OP_CUT, &once, &second).is_ok());
}

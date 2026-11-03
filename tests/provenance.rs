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

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Curve> {
    vec![line((x0, y0), (x1, y0)), line((x1, y0), (x1, y1)), line((x1, y1), (x0, y1)), line((x0, y1), (x0, y0))]
}

fn plate() -> Body { extrude(&Profile { loops: vec![rect(0.0, 0.0, 6.0, 4.0)], n: Z, low: 0.0, high: 0.5 }, TOL).unwrap() }

/// A drill extruded along z through the plate (spanning it), centre (x, y).
fn drill(x: f64, y: f64, r: f64) -> Body {
    let circle = Curve::Arc { c: p(x, y), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r, t0: 0.0, t1: 2.0 * PI };
    extrude(&Profile { loops: vec![vec![circle]], n: Z, low: -1.0, high: 1.5 }, TOL).unwrap()
}

fn block(x0: f64, y0: f64, x1: f64, y1: f64) -> Body { extrude(&Profile { loops: vec![rect(x0, y0, x1, y1)], n: Z, low: -1.0, high: 1.5 }, TOL).unwrap() }

fn cut_reply(op: u32, body: &Body, tool: &Body) -> Result<Vec<f64>, gc_kernel::KernelError> {
    let data: Vec<f64> = gc_kernel::neutral::write_body(body).into_iter().chain(gc_kernel::neutral::write_body(tool)).collect();
    handle(&parse_request(&encode_request(op, TOL, &data)).unwrap())
}

fn on_surface(s: &Surface, q: Vec3) -> bool {
    match *s {
        Surface::Plane { p, n, .. } => q.sub(p).dot(n).abs() <= TOL,
        Surface::Cylinder { p, a, r, .. } => { let w = q.sub(p); (w.sub(a.scale(w.dot(a))).length() - r).abs() <= TOL }
        Surface::Cone { p, a, r, k, .. } => { let w = q.sub(p); let s = w.dot(a); (w.sub(a.scale(s)).length() - (r + k * s)).abs() <= TOL }
        _ => panic!("surface {s:?}"),
    }
}

fn on_curve(c: &Curve, q: Vec3) -> bool {
    match *c {
        Curve::Line { a, b } => { let d = b.sub(a); q.sub(a).cross(d).length() <= TOL * d.length() }
        Curve::Arc { c, n, r, .. } => q.sub(c).dot(n).abs() <= TOL && (q.sub(c).length() - r).abs() <= TOL,
        _ => panic!("curve {c:?}"),
    }
}

fn samples(c: &Curve) -> Vec<Vec3> { let (t0, t1) = c.domain(); (0..=8).map(|k| c.eval(t0 + (t1 - t0) * k as f64 / 8.0).0).collect() }

fn outward(s: &Surface, sense: bool) -> Option<Vec3> { match *s { Surface::Plane { n, .. } => Some(if sense { n } else { n.scale(-1.0) }), _ => None } }

fn prove_cut(body: &Body, tool: &Body, result: &Body, prov: &Provenance, outcomes: &[Outcome]) {
    assert_eq!((prov.faces.len(), prov.edges.len(), prov.vertices.len()), (result.faces.len(), result.edges.len(), result.vertices.len()));
    let arity = |role: u32| match role { FROM_BODY_VERTEX | FROM_TOOL_VERTEX => 1, MEETING | MEETING_POINT => 3, MEETING_END => 4, _ => 2 };
    for r in prov.faces.iter().chain(&prov.edges).chain(&prov.vertices) { assert_eq!(r.refs.len(), arity(r.role), "{r:?}"); }
    for list in [&prov.faces, &prov.edges, &prov.vertices] {
        for (i, a) in list.iter().enumerate() {
            if !matches!(a.role, FROM_BODY_FACE | FROM_TOOL_FACE | FROM_BODY_EDGE | FROM_TOOL_EDGE) { assert!(!list[i + 1..].contains(a), "generated name twice {a:?}"); }
        }
    }
    let bsurf = |i: usize| &body.faces[i].surface;
    let tsurf = |j: usize| &tool.faces[j].surface;
    for (v, r) in prov.vertices.iter().enumerate() {
        let q = result.vertices[v];
        let ok = match r.role {
            FROM_BODY_VERTEX => near(q, body.vertices[r.refs[0]]),
            SEAM_POINT | MEETING_POINT => on_surface(bsurf(r.refs[0]), q) && on_surface(tsurf(r.refs[1]), q),
            MEETING_END => on_surface(bsurf(r.refs[0]), q) && on_surface(tsurf(r.refs[1]), q) && on_surface(bsurf(r.refs[3]), q),
            TOOL_EDGE_END => on_curve(&tool.edges[r.refs[0]].curve, q) && on_surface(bsurf(r.refs[1]), q),
            other => panic!("vertex role {other}"),
        };
        assert!(ok, "vertex {v} {r:?}");
    }
    for (e, r) in prov.edges.iter().enumerate() {
        let edge = &result.edges[e];
        let (s, t) = (&prov.vertices[edge.start], &prov.vertices[edge.end]);
        match r.role {
            FROM_BODY_EDGE => {
                let source = &body.edges[r.refs[0]];
                assert!(samples(&edge.curve).iter().all(|&q| on_curve(&source.curve, q)), "edge {e} {r:?} lies on its source");
                if r.refs[1] == 0 {
                    let ends = [s.clone(), t.clone()];
                    let want = [Record::new(FROM_BODY_VERTEX, &[source.start]), Record::new(FROM_BODY_VERTEX, &[source.end])];
                    assert!(ends == want || ends == [want[1].clone(), want[0].clone()], "edge {e} {r:?} keeps its source's ends");
                    assert!(near(mid(&edge.curve), mid(&source.curve)), "edge {e} {r:?} is its whole source");
                }
            }
            FROM_TOOL_EDGE => {
                assert_eq!(r.refs[1], 1);
                assert!(samples(&edge.curve).iter().all(|&q| on_curve(&tool.edges[r.refs[0]].curve, q)), "edge {e} {r:?}");
                assert!([s, t].iter().all(|x| x.role == TOOL_EDGE_END && x.refs[0] == r.refs[0]), "edge {e} {r:?} ends");
            }
            MEETING => {
                assert!(samples(&edge.curve).iter().all(|&q| on_surface(bsurf(r.refs[0]), q) && on_surface(tsurf(r.refs[1]), q)), "edge {e} {r:?}");
                if edge.start == edge.end && s.role == SEAM_POINT { assert_eq!(s.refs, r.refs[..2], "edge {e} circle from its seam point"); }
                if edge.start == edge.end && s.role == MEETING_POINT { assert_eq!(s.refs, r.refs, "edge {e} loop from its own point"); }
                if s.role == MEETING_END && t.role == MEETING_END && s.refs[..3] == r.refs[..] { assert_eq!(t.refs[..3], r.refs[..], "edge {e} ruling ends"); }
            }
            TOOL_SEAM => {
                assert!(matches!(edge.curve, Curve::Line { .. }) && on_surface(tsurf(r.refs[1]), mid(&edge.curve)), "edge {e} {r:?}");
                let seam_end = |x: &Record| matches!(x.role, SEAM_POINT | MEETING_POINT) && x.refs[1] == r.refs[1];
                assert!(seam_end(s) && s.refs[0] == r.refs[0] && seam_end(t), "edge {e} ends");
            }
            other => panic!("edge role {other}"),
        }
    }
    for (f, r) in prov.faces.iter().enumerate() {
        let face = &result.faces[f];
        let bounding: Vec<usize> = face.loops.iter().flatten().map(|&(e, _)| e).collect();
        let source = match r.role { FROM_BODY_FACE => &body.faces[r.refs[0]], FROM_TOOL_FACE => &tool.faces[r.refs[0]], other => panic!("face role {other}") };
        for &e in &bounding { assert!(samples(&result.edges[e].curve).iter().all(|&q| on_surface(&source.surface, q)), "face {f} {r:?} edge {e} on its source surface"); }
        match r.role {
            FROM_BODY_FACE => {
                if let Some(n) = outward(&source.surface, source.sense) { assert!(near(outward(&face.surface, face.sense).unwrap(), n), "face {f} keeps its outward side"); }
                let mapped: Vec<&Record> = bounding.iter().map(|&e| &prov.edges[e]).collect();
                let whole = mapped.iter().all(|er| er.role == FROM_BODY_EDGE && er.refs[1] == 0)
                    && source.loops.iter().flatten().count() == bounding.len()
                    && source.loops.iter().flatten().all(|(se, _)| mapped.iter().any(|er| er.refs[0] == *se));
                assert_eq!(r.refs[1] == 0, whole, "face {f} {r:?} trimmed flag");
            }
            _ => {
                assert_eq!(r.refs[1], 1);
                for &e in &bounding {
                    let er = &prov.edges[e];
                    let ok = match er.role {
                        MEETING | TOOL_SEAM => er.refs[1] == r.refs[0],
                        FROM_TOOL_EDGE => tool.faces[r.refs[0]].loops.iter().flatten().any(|(te, _)| *te == er.refs[0]),
                        _ => false,
                    };
                    assert!(ok, "face {f} {r:?} bounded by {er:?}");
                }
            }
        }
    }
    // Outcomes: each input entity's descendants, read independently from the records.
    let inputs = [body.faces.len(), body.edges.len(), body.vertices.len(), tool.faces.len(), tool.edges.len(), tool.vertices.len()];
    assert_eq!(outcomes.len(), inputs.iter().sum::<usize>());
    let kinds = [(&prov.faces, FROM_BODY_FACE, DELETED), (&prov.edges, FROM_BODY_EDGE, DELETED), (&prov.vertices, FROM_BODY_VERTEX, DELETED),
                 (&prov.faces, FROM_TOOL_FACE, CONSUMED), (&prov.edges, FROM_TOOL_EDGE, CONSUMED), (&prov.vertices, FROM_TOOL_VERTEX, CONSUMED)];
    let mut o = outcomes.iter();
    for ((records, role, none), &n) in kinds.into_iter().zip(&inputs) {
        for i in 0..n {
            let out = o.next().unwrap();
            let want: Vec<usize> = (0..records.len()).filter(|&k| records[k].role == role && records[k].refs[0] == i).collect();
            assert_eq!(out.outputs, want);
            let code = match want.len() { 0 => none, 1 if records[want[0]].refs.get(1).map_or(true, |&t| t == 0) => UNCHANGED, 1 => MODIFIED, _ => SPLIT };
            assert_eq!(out.code, code);
        }
    }
}

fn named(body: &Body, tool: &Body) -> (Body, Provenance, Vec<Outcome>) {
    let (result, prov, outcomes) = gc_kernel::boolean::cut_named(body, tool, TOL).unwrap();
    prove_cut(body, tool, &result, &prov, &outcomes);
    assert_eq!(gc_kernel::neutral::write_body(&result), gc_kernel::neutral::write_body(&gc_kernel::boolean::cut(body, tool, TOL).unwrap()));
    (result, prov, outcomes)
}

fn count(prov: &Provenance, role: u32) -> usize { prov.faces.iter().chain(&prov.edges).chain(&prov.vertices).filter(|r| r.role == role).count() }

fn is_z_face(body: &Body, i: usize) -> bool { matches!(body.faces[i].surface, Surface::Plane { n, .. } if n.cross(Z).length() <= TOL) }

#[test]
fn bore_history_names_the_plate_and_the_drill() {
    let (body, tool) = (plate(), drill(2.0, 1.5, 0.5));
    let (_, prov, outcomes) = named(&body, &tool);
    // The two plate faces the hole passes through are modified, the rest unchanged; the drill's wall becomes the hole's
    // wall, everything else of the drill is consumed.
    let wall = tool.faces.iter().position(|f| matches!(f.surface, Surface::Cylinder { .. })).unwrap();
    let hole = prov.faces.iter().position(|r| r.role == FROM_TOOL_FACE).unwrap();
    let (nf, ne, nv) = (body.faces.len(), body.edges.len(), body.vertices.len());
    for (i, out) in outcomes[..nf].iter().enumerate() {
        assert_eq!((out.code, &out.outputs[..]), (if is_z_face(&body, i) { MODIFIED } else { UNCHANGED }, &[i][..]), "plate face {i}");
    }
    assert!(outcomes[nf..nf + ne + nv].iter().enumerate().all(|(i, o)| o.code == UNCHANGED && o.outputs == [if i < ne { i } else { i - ne }]));
    for (j, out) in outcomes[nf + ne + nv..nf + ne + nv + tool.faces.len()].iter().enumerate() {
        assert_eq!((out.code, out.outputs.clone()), if j == wall { (MODIFIED, vec![hole]) } else { (CONSUMED, vec![]) }, "drill face {j}");
    }
    assert!(outcomes[nf + ne + nv + tool.faces.len()..].iter().all(|o| o.code == CONSUMED && o.outputs.is_empty()));
}

/// A second hole runs in the section domain: the plate's sides and the first hole stay unchanged, its two faces are
/// modified, and the second drill's wall and seam become the new hole's.
#[test]
fn second_hole_history_keeps_the_first() {
    let once = gc_kernel::boolean::cut(&plate(), &drill(2.0, 1.5, 0.5), TOL).unwrap();
    let second = drill(4.5, 2.5, 0.75);
    let (_, prov, outcomes) = named(&once, &second);
    let nf = once.faces.len();
    for (i, out) in outcomes[..nf].iter().enumerate() {
        assert_eq!(out.code, if is_z_face(&once, i) { MODIFIED } else { UNCHANGED }, "face {i}");
    }
    let n = nf + once.edges.len() + once.vertices.len();
    assert!(outcomes[nf..n].iter().all(|o| o.code == UNCHANGED), "every edge and vertex of the holed plate is kept");
    let wall = second.faces.iter().position(|f| matches!(f.surface, Surface::Cylinder { .. })).unwrap();
    assert_eq!(outcomes[n + wall].code, MODIFIED);
    assert_eq!((count(&prov, FROM_TOOL_FACE), count(&prov, FROM_TOOL_EDGE), count(&prov, MEETING), count(&prov, TOOL_EDGE_END)), (1, 1, 2, 2));
}

/// A corner notch: two plate sides trimmed, the block's inner corner kept as a tool ruling, two generated meetings.
#[test]
fn corner_notch_history() {
    let (body, tool) = (plate(), block(5.0, 3.0, 7.0, 5.0));
    let (_, prov, outcomes) = named(&body, &tool);
    let nf = body.faces.len();
    let modified = outcomes[..nf].iter().filter(|o| o.code == MODIFIED).count();
    assert_eq!(modified, 4, "two caps and the two sides the notch crosses");
    assert_eq!(outcomes[..nf].iter().filter(|o| o.code == UNCHANGED).count(), 2);
    assert_eq!((count(&prov, FROM_TOOL_FACE), count(&prov, FROM_TOOL_EDGE), count(&prov, MEETING_END), count(&prov, TOOL_EDGE_END)), (2, 1, 4, 2));
    // The plate corner the notch removes is deleted, with its ruling.
    assert_eq!(outcomes[nf + body.edges.len()..nf + body.edges.len() + body.vertices.len()].iter().filter(|o| o.code == DELETED).count(), 2);
}

/// A slot across the plate's bottom side splits that side face in two; a drill across the side meets it twice,
/// numbered 0 and 1 along the side.
#[test]
fn split_side_and_numbered_meetings() {
    let (body, slot) = (plate(), block(2.0, -1.0, 3.0, 1.0));
    let (_, prov, outcomes) = named(&body, &slot);
    assert_eq!(outcomes[..body.faces.len()].iter().filter(|o| o.code == SPLIT).count(), 1, "the crossed side splits");
    assert_eq!(count(&prov, FROM_TOOL_FACE), 3);
    let (_, prov, outcomes) = named(&body, &drill(3.0, 0.2, 0.5));
    let side = (0..body.faces.len()).find(|&i| outcomes[i].code == SPLIT).expect("the drilled side splits");
    let mut n: Vec<usize> = prov.edges.iter().filter(|r| r.role == MEETING && r.refs[0] == side).map(|r| r.refs[2]).collect();
    n.sort();
    assert_eq!(n, [0, 1]);
}

#[test]
fn op12_reply_is_op3_reply_plus_envelope() {
    for (body, tool) in [(plate(), drill(4.0, 2.5, 0.75)), (plate(), block(5.0, 3.0, 7.0, 5.0))] {
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
}

/// A domain that records no history refuses op 12 by name; op 3 still cuts.
#[test]
fn unrecorded_domain_refuses_the_named_cut() {
    let (body, tool) = (plate(), extrude(&Profile { loops: vec![rect(1.0, 1.0, 2.0, 2.0)], n: Z, low: 0.2, high: 1.5 }, TOL).unwrap());
    match gc_kernel::boolean::cut_named(&body, &tool, TOL) {
        Err(gc_kernel::KernelError::Unsupported(m)) => assert!(m.contains("cut history is not recorded"), "{m}"),
        other => panic!("expected a refusal, got {:?}", other.map(|_| ())),
    }
    assert!(cut_reply(OP_CUT_NAMED, &body, &tool).is_err());
    assert!(cut_reply(OP_CUT, &body, &tool).is_ok());
}
/// Bob's NC1 plate p9999 corner-waste triangle (legs one ulp inside the plate faces; section-domain fixture).
#[test]
fn bobs_p9999_corner_cut_history() {
    let read = |name: &str| {
        let bytes = std::fs::read(format!("{}/tests/fixtures/section/p9999/{name}.f64", env!("CARGO_MANIFEST_DIR"))).unwrap();
        let values: Vec<f64> = bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
        gc_kernel::neutral::Reader::with_tolerance(&values, 1e-6).body().unwrap()
    };
    let (body, tool) = (read("body"), read("tool"));
    let (result, prov, outcomes) = gc_kernel::boolean::cut_named(&body, &tool, 1e-6).unwrap();
    assert_eq!(gc_kernel::neutral::write_body(&result), gc_kernel::neutral::write_body(&gc_kernel::boolean::cut(&body, &tool, 1e-6).unwrap()));
    prove_cut(&body, &tool, &result, &prov, &outcomes);
    // Two plate sides kept whole, the two the legs cross trimmed, the hypotenuse from the triangle.
    assert_eq!((count(&prov, FROM_TOOL_FACE), prov.faces.iter().filter(|r| r.role == FROM_BODY_FACE && r.refs[1] == 0).count()), (1, 2));
}

// ---- E3c: first interior sketch cuts whose arcs no passage domain takes run in the section domain ----

fn arc(cx: f64, cy: f64, r: f64, from: f64, sweep: f64) -> Curve {
    Curve::Arc { c: p(cx, cy), n: Z, x: Vec3::new(from.cos(), from.sin(), 0.0), r, t0: 0.0, t1: sweep }
}

/// Rounded rectangle (x0..x1, y0..y1), corner radius r, counter-clockwise: four arcs.
fn rounded(x0: f64, y0: f64, x1: f64, y1: f64, r: f64) -> Vec<Curve> {
    vec![
        line((x0 + r, y0), (x1 - r, y0)), arc(x1 - r, y0 + r, r, -PI / 2.0, PI / 2.0),
        line((x1, y0 + r), (x1, y1 - r)), arc(x1 - r, y1 - r, r, 0.0, PI / 2.0),
        line((x1 - r, y1), (x0 + r, y1)), arc(x0 + r, y1 - r, r, PI / 2.0, PI / 2.0),
        line((x0, y1 - r), (x0, y0 + r)), arc(x0 + r, y0 + r, r, PI, PI / 2.0),
    ]
}

/// A D: a flat bottom and a half circle on top, counter-clockwise: one arc.
fn d_shape(cx: f64, cy: f64, r: f64) -> Vec<Curve> { vec![line((cx - r, cy), (cx + r, cy)), arc(cx, cy, r, 0.0, PI)] }

fn through(loop_: Vec<Curve>) -> Body { extrude(&Profile { loops: vec![loop_], n: Z, low: -1.0, high: 1.5 }, TOL).unwrap() }

fn volume(body: &Body) -> f64 { gc_kernel::mass::measure(body, TOL).unwrap().volume }

/// Both refused outright before E3c; now cut exactly, recorded, and read back by op 13 as one hole of the tool's curves.
#[test]
fn first_interior_sketch_cuts_run_in_the_section_domain() {
    let r = 0.5;
    for (tool, area, curves) in [
        (through(rounded(1.5, 1.0, 4.5, 3.0, r)), 6.0 - (4.0 - PI) * r * r, 8),
        (through(d_shape(3.0, 2.0, 0.75)), PI * 0.75 * 0.75 / 2.0, 2),
    ] {
        let body = plate();
        let (result, prov, outcomes) = named(&body, &tool);
        let want = 0.5 * (24.0 - area);
        assert!((volume(&result) - want).abs() <= 1e-12 * want, "{} vs {want}", volume(&result));
        // The caps are modified, the four sides kept; every tool side becomes the hole's wall.
        let nf = body.faces.len();
        for (i, out) in outcomes[..nf].iter().enumerate() {
            assert_eq!(out.code, if is_z_face(&body, i) { MODIFIED } else { UNCHANGED }, "plate face {i}");
        }
        assert_eq!(count(&prov, FROM_TOOL_FACE), curves);
        let profile = gc_kernel::rebuild::prism_profile(&result, Z, TOL).unwrap();
        assert_eq!(profile.loops.iter().map(|l| l.len()).collect::<Vec<_>>(), [4, curves]);
    }
}

/// Either order of a drill and a rounded cut-out gives the same exact plate, every step recorded.
#[test]
fn rounded_cutout_and_drill_in_either_order() {
    let (rounded_tool, hole) = (through(rounded(1.5, 1.0, 4.0, 3.0, 0.5)), drill(5.0, 2.0, 0.4));
    let a = named(&named(&plate(), &rounded_tool).0, &hole).0;
    let b = named(&named(&plate(), &hole).0, &rounded_tool).0;
    let want = 0.5 * (24.0 - (5.0 - (4.0 - PI) * 0.25) - PI * 0.16);
    for v in [volume(&a), volume(&b)] { assert!((v - want).abs() <= 1e-12 * want, "{v} vs {want}"); }
}

fn slot_tool() -> Body { through(vec![line((2.0, 1.5), (4.0, 1.5)), arc(4.0, 2.0, 0.5, -PI / 2.0, PI), line((4.0, 2.5), (2.0, 2.5)), arc(2.0, 2.0, 0.5, PI / 2.0, PI)]) }
fn turned_tool() -> Body { through(vec![line((3.0, 1.0), (4.0, 2.0)), line((4.0, 2.0), (3.0, 3.0)), line((3.0, 3.0), (2.0, 2.0)), line((2.0, 2.0), (3.0, 1.0))]) }

/// The first interior rectangle, slot and turned square in a plain plate (2026-11-02): the box passage splices its
/// cavity (the rectangle no longer rebuilds cells), the slot rebuilds the box with the slot as a hole, the square
/// splices: each recorded, the caps modified, the four sides unchanged, exact volumes; op 12 body = op 3 body.
#[test]
fn first_interior_cuts_in_a_plate_record_history() {
    for (tool, area, walls) in [
        (block(1.0, 1.0, 2.0, 2.0), 1.0, 4),
        (slot_tool(), 2.0 + PI * 0.25, 4),
        (turned_tool(), 2.0, 4),
    ] {
        let body = plate();
        let (result, prov, outcomes) = named(&body, &tool);
        let want = 0.5 * (24.0 - area);
        assert!((volume(&result) - want).abs() <= 1e-12 * want, "{} vs {want}", volume(&result));
        let nf = body.faces.len();
        for (i, out) in outcomes[..nf].iter().enumerate() {
            assert_eq!(out.code, if is_z_face(&body, i) { MODIFIED } else { UNCHANGED }, "plate face {i}");
        }
        assert_eq!((count(&prov, FROM_TOOL_FACE), count(&prov, FROM_TOOL_EDGE), count(&prov, MEETING), count(&prov, TOOL_EDGE_END)),
            (walls, walls, 2 * walls, 2 * walls));
    }
}

// ---- Persistent names on every stock (Claude 2026-11-02): cuts through angle legs, tube walls, pipe and cone walls ----

fn p3(x: f64, y: f64, z: f64) -> Vec3 { Vec3::new(x, y, z) }
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);

/// A closed polygon in the XZ plane (x, z), swept along +Y from 0 to `length`: angle, tube and channel stock as the
/// shape factories lay them out.
fn stock(points: &[(f64, f64)], holes: &[&[(f64, f64)]], length: f64) -> Body {
    let polygon = |pts: &[(f64, f64)]| (0..pts.len()).map(|i| {
        let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
        Curve::Line { a: p3(a.0, 0.0, a.1), b: p3(b.0, 0.0, b.1) }
    }).collect::<Vec<_>>();
    let mut loops = vec![polygon(points)];
    loops.extend(holes.iter().map(|h| polygon(h)));
    extrude(&Profile { loops, n: Y, low: 0.0, high: length }, TOL).unwrap()
}

fn angle() -> Body { stock(&[(0.0, 0.0), (3.0, 0.0), (3.0, -0.25), (0.25, -0.25), (0.25, -3.0), (0.0, -3.0)], &[], 12.0) }

/// A cutter along z from `low` to `high` through (x, y) of the given section.
fn along_z(curves: Vec<Curve>, low: f64, high: f64) -> Body { extrude(&Profile { loops: vec![curves], n: Z, low, high }, TOL).unwrap() }

/// Hole, slot, rectangle and turned square through an angle's horizontal leg, one after another (each cutter spanning the
/// whole angle, as the app's Through All does): every cut recorded,
/// the angle's own faces that the cuts pass through modified, every other angle face unchanged at each step.
#[test]
fn cuts_through_an_angle_leg_keep_the_angle_names() {
    let circle = |x: f64, y: f64, r: f64| vec![Curve::Arc { c: p(x, y), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r, t0: 0.0, t1: 2.0 * PI }];
    let shift = |curves: Vec<Curve>, dy: f64| curves.into_iter().map(|c| c.translated(Vec3::new(0.0, dy, 0.0))).collect::<Vec<_>>();
    let slot = vec![line((1.0, 1.5), (2.0, 1.5)), arc(2.0, 2.0, 0.5, -PI / 2.0, PI), line((2.0, 2.5), (1.0, 2.5)), arc(1.0, 2.0, 0.5, PI / 2.0, PI)];
    let square = vec![line((1.5, 0.0), (2.0, 0.5)), line((2.0, 0.5), (1.5, 1.0)), line((1.5, 1.0), (1.0, 0.5)), line((1.0, 0.5), (1.5, 0.0))];
    let tools = [
        along_z(circle(1.5, 1.0, 0.25), -4.0, 1.0),
        along_z(shift(slot, 2.0), -4.0, 1.0),
        along_z(rect(1.0, 6.0, 2.0, 7.0), -4.0, 1.0),
        along_z(shift(square, 8.0), -4.0, 1.0),
    ];
    let mut body = angle();
    let start = volume(&body);
    let mut removed = 0.0;
    for (k, (tool, area)) in tools.iter().zip([PI * 0.0625, 1.0 + PI * 0.25, 1.0, 0.5]).enumerate() {
        let (result, prov, outcomes) = named(&body, tool);
        removed += 0.25 * area;
        assert!((start - volume(&result) - removed).abs() <= 1e-9, "cut {k}");
        let modified: Vec<usize> = (0..body.faces.len()).filter(|&i| outcomes[i].code == MODIFIED).collect();
        assert_eq!(modified.len(), 2, "cut {k}: the leg's top and bottom faces");
        assert!(modified.iter().all(|&i| is_z_face(&body, i)), "cut {k}");
        assert!(outcomes[..body.faces.len()].iter().all(|o| matches!(o.code, MODIFIED | UNCHANGED)), "cut {k}");
        assert!(count(&prov, FROM_TOOL_FACE) >= 1);
        body = result;
    }
}

/// A rectangular tube drilled through both walls, then cut through both walls by a rectangle: two passages each,
/// the four wall faces modified, every other face unchanged; generated names never repeat.
#[test]
fn cuts_through_both_walls_of_a_tube() {
    let tube = stock(&[(0.0, 0.0), (0.0, -2.0), (4.0, -2.0), (4.0, 0.0)], &[&[(0.125, -0.125), (3.875, -0.125), (3.875, -1.875), (0.125, -1.875)]], 12.0);
    let hole = along_z(vec![Curve::Arc { c: p(2.0, 3.0), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r: 0.25, t0: 0.0, t1: 2.0 * PI }], -3.0, 1.0);
    let (once, prov, outcomes) = named(&tube, &hole);
    assert_eq!(outcomes[..tube.faces.len()].iter().filter(|o| o.code == MODIFIED).count(), 4);
    assert_eq!((count(&prov, FROM_TOOL_FACE), count(&prov, MEETING), count(&prov, TOOL_SEAM)), (2, 4, 2));
    let (_, prov, outcomes) = named(&once, &along_z(rect(1.5, 6.0, 2.5, 7.5), -3.0, 1.0));
    assert_eq!(outcomes[..once.faces.len()].iter().filter(|o| o.code == MODIFIED).count(), 4);
    assert_eq!((count(&prov, FROM_TOOL_FACE), count(&prov, MEETING), count(&prov, TOOL_EDGE_END)), (8, 16, 16));
}

/// A pipe along +Y (outer ro, inner ri), as the shape factory places it.
fn pipe(ro: f64, ri: f64) -> Body {
    let ring = |r: f64| vec![Curve::Arc { c: p3(0.0, 0.0, 0.0), n: Vec3::new(0.0, -1.0, 0.0), x: Vec3::new(1.0, 0.0, 0.0), r, t0: 0.0, t1: 2.0 * PI }];
    extrude(&Profile { loops: vec![ring(ro), ring(ri)], n: Y, low: 0.0, high: 10.0 }, TOL).unwrap()
}

/// Radial bores through one wall and through both walls of a pipe, then of a tapered tube (cone_bore): the walls the
/// hole crosses are modified, everything else of the stock unchanged; each crossing a numbered meeting loop with its
/// own point (MEETING_POINT, role 42), each wall passage a seam and a wall from the drill.
#[test]
fn radial_and_conical_bores_record_history() {
    let drill_z = |y: f64, r: f64, low: f64, high: f64| along_z(vec![Curve::Arc { c: p(0.0, y), n: Z, x: Vec3::new(1.0, 0.0, 0.0), r, t0: 0.0, t1: 2.0 * PI }], low, high);
    let body = pipe(1.75, 1.5);
    for (tool, passages) in [(drill_z(4.0, 0.25, 1.0, 1.751), 1), (drill_z(6.0, 0.25, -2.5, 2.5), 2)] {
        let (_, prov, outcomes) = named(&body, &tool);
        assert_eq!(outcomes[..body.faces.len()].iter().filter(|o| o.code == MODIFIED).count(), 2, "outer and inner walls");
        assert_eq!((count(&prov, MEETING), count(&prov, MEETING_POINT), count(&prov, TOOL_SEAM), count(&prov, FROM_TOOL_FACE)),
            (2 * passages, 2 * passages, passages, passages));
        if passages == 2 {
            let mut n: Vec<usize> = prov.edges.iter().filter(|r| r.role == MEETING).map(|r| r.refs[2]).collect();
            n.sort();
            assert_eq!(n, [0, 0, 1, 1], "each wall met twice, numbered along the drill");
        }
    }
    use gc_kernel::revolve::{revolve_general, Kind, Revolution, Segment};
    let section = [(2.0, 0.0), (1.5, 10.0), (1.25, 10.0), (1.75, 0.0)];
    let cone = revolve_general(&Revolution { p: p3(0.0, 0.0, 0.0), a: Y, x: Vec3::new(1.0, 0.0, 0.0), angle: 2.0 * PI,
        loops: vec![section.iter().map(|&(r, z)| Segment { r, z, kind: Kind::Line }).collect()] }, TOL).unwrap();
    let (_, prov, outcomes) = named(&cone, &drill_z(5.0, 0.25, -3.0, 3.0));
    assert_eq!(outcomes[..cone.faces.len()].iter().filter(|o| o.code == MODIFIED).count(), 2, "both conical walls");
    assert_eq!((count(&prov, MEETING_POINT), count(&prov, TOOL_SEAM)), (4, 2));
}

/// Op 26: op 3's reply, then 1 and op 12's envelope for a recorded cut, or 0, n and the reason's n UTF-8 bytes.
#[test]
fn op26_reply_is_op3_reply_and_history_or_reason() {
    use gc_kernel::protocol::OP_CUT_RECORDED;
    let (body, tool) = (plate(), block(1.0, 1.0, 2.0, 2.0));
    let (three, recorded) = (cut_reply(OP_CUT, &body, &tool).unwrap(), cut_reply(OP_CUT_RECORDED, &body, &tool).unwrap());
    let twelve = cut_reply(OP_CUT_NAMED, &body, &tool).unwrap();
    assert_eq!(recorded[..three.len()], three[..]);
    assert_eq!(recorded[three.len()], 1.0);
    assert_eq!(recorded[three.len() + 1..], twelve[three.len()..]);
    // A blind pocket: same body as op 3, no history, and why.
    let pocket = extrude(&Profile { loops: vec![rect(1.0, 1.0, 2.0, 2.0)], n: Z, low: 0.2, high: 1.5 }, TOL).unwrap();
    let (three, recorded) = (cut_reply(OP_CUT, &body, &pocket).unwrap(), cut_reply(OP_CUT_RECORDED, &body, &pocket).unwrap());
    assert_eq!(recorded[..three.len()], three[..]);
    assert_eq!(recorded[three.len()], 0.0);
    let n = recorded[three.len() + 1] as usize;
    assert_eq!(recorded.len(), three.len() + 2 + n);
    let reason = String::from_utf8(recorded[three.len() + 2..].iter().map(|&b| b as u8).collect()).unwrap();
    assert!(reason.contains("cut history is not recorded"), "{reason}");
}

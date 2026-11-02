//! Union (ops 24/25, Claude 2026-11-02): SolidWorks-style bosses on plates, webs and flanges. Every result passes the
//! topology check and tessellates; every volume is the closed-form sum of the operands less their closed-form overlap.

use gc_kernel::boolean::{cut, union, union_named};
use gc_kernel::geom::{Curve, Surface, Vec3};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{write_body, Body, Profile, Reader};
use gc_kernel::provenance::{self, Outcome};
use gc_kernel::revolve::{Kind, Revolution, Segment};
use gc_kernel::topology::check;
use gc_kernel::{protocol, KernelError};
use std::f64::consts::{PI, TAU};

const TOL: f64 = 1e-6;

fn v(x: f64, y: f64, z: f64) -> Vec3 { Vec3::new(x, y, z) }

/// A polygon of lines in the plane z (counter-clockwise when the points are).
fn polygon(points: &[(f64, f64)], z: f64) -> Vec<Curve> {
    (0..points.len()).map(|i| {
        let (a, b) = (points[i], points[(i + 1) % points.len()]);
        Curve::Line { a: v(a.0, a.1, z), b: v(b.0, b.1, z) }
    }).collect()
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64, z: f64) -> Vec<Curve> { polygon(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)], z) }

fn circle(x: f64, y: f64, r: f64, z: f64) -> Vec<Curve> {
    vec![Curve::Arc { c: v(x, y, z), n: v(0., 0., 1.), x: v(1., 0., 0.), r, t0: 0., t1: TAU }]
}

fn prism(loops: Vec<Vec<Curve>>, n: Vec3, low: f64, high: f64) -> Body {
    gc_kernel::extrude::extrude(&Profile { loops, n, low, high }, TOL).unwrap()
}

fn up(loops: Vec<Vec<Curve>>, low: f64, high: f64) -> Body { prism(loops, v(0., 0., 1.), low, high) }

/// Bob's plate: 10 x 8, 0.5 thick, top at z 0.
fn plate() -> Body { up(vec![rect(0., 0., 10., 8., 0.)], -0.5, 0.) }

fn volume(b: &Body) -> f64 { measure(b, TOL).unwrap().volume }

/// The union, checked: valid, tessellates, its volume the expected one (to 1e-9 relative).
#[track_caller]
fn united(body: &Body, tool: &Body, expected: f64) -> Body {
    let result = match union(body, tool, TOL) { Ok(r) => r, Err(e) => panic!("{e}") };
    check(&result, TOL).unwrap();
    gc_kernel::tessellate::tessellate(&result, 0.01).unwrap();
    let got = volume(&result);
    assert!((got - expected).abs() <= 1e-9 * expected.max(1.), "union volume {got}, expected {expected}");
    result
}

#[track_caller]
fn refused(body: &Body, tool: &Body) -> String {
    match union(body, tool, TOL) {
        Err(KernelError::Unsupported(m)) => m,
        other => panic!("expected a named refusal, got {:?}", other.map(|b| b.faces.len())),
    }
}

fn planes(b: &Body) -> usize { b.faces.iter().filter(|f| matches!(f.surface, Surface::Plane { .. })).count() }

#[test]
fn rectangular_and_round_bosses_on_a_plate_keep_the_plate_entity_for_entity() {
    let p = plate();
    // 2 x 3 rectangle, 1 high: the plate's 6 faces kept, the top gains one inner loop, the boss adds 5 faces.
    let r = united(&p, &up(vec![rect(4., 2., 6., 5., 0.)], 0., 1.), 40. + 6.);
    assert_eq!(r.faces.len(), 6 + 5);
    for (i, f) in p.faces.iter().enumerate() {
        let extra = if f.loops.len() != r.faces[i].loops.len() { 1 } else { 0 };
        assert_eq!(r.faces[i].loops.len(), f.loops.len() + extra);
        assert_eq!(r.faces[i].loops[..f.loops.len()], f.loops[..]);
    }
    assert_eq!(r.faces.iter().filter(|f| f.loops.len() == 2).count(), 1, "only the top face gains a loop");
    assert_eq!(&r.vertices[..p.vertices.len()], &p.vertices[..]);
    // A round boss, r 1, 0.75 high.
    united(&p, &up(vec![circle(5., 4., 1., 0.)], 0., 0.75), 40. + PI * 0.75);
    // A ring boss (a tube standing on the plate): the plate inside the ring stays a face of its own.
    let ring = up(vec![circle(5., 4., 1.5, 0.), circle(5., 4., 1., 0.)], 0., 2.);
    let r = united(&p, &ring, 40. + PI * (2.25 - 1.) * 2.);
    assert_eq!(r.faces.len(), 6 + 3 + 1);
    // On the bottom face, downward; on the right side face, along +x.
    united(&p, &prism(vec![rect(1., 1., 3., 2., -0.5)], v(0., 0., -1.), 0., 0.25), 40. + 2. * 0.25);
    let side = prism(vec![polygon(&[(2., -0.4), (4., -0.4), (4., -0.1), (2., -0.1)], 0.)].into_iter()
        .map(|lp| lp.into_iter().map(|c| match c { Curve::Line { a, b } => Curve::Line { a: v(10., a.x, a.y), b: v(10., b.x, b.y) }, c => c }).collect())
        .collect(), v(1., 0., 0.), 0., 1.5);
    united(&p, &side, 40. + 2. * 0.3 * 1.5);
    // A second boss on the first boss's top, and one more on the plate beside it.
    let once = united(&p, &up(vec![rect(4., 2., 6., 5., 0.)], 0., 1.), 46.);
    let twice = united(&once, &up(vec![circle(5., 3.5, 0.5, 0.)], 1., 1.5), 46. + PI * 0.25 * 0.5);
    united(&twice, &up(vec![rect(8., 6., 9., 7., 0.)], 0., 0.2), 46. + PI * 0.125 + 0.2);
}

#[test]
fn revolved_bosses_on_a_plate() {
    let p = plate();
    // A dome: a quarter circle and the radius below it, about a vertical axis through (5, 4).
    let dome = Revolution {
        p: v(5., 4., 0.), a: v(0., 0., 1.), x: v(1., 0., 0.), angle: TAU,
        loops: vec![vec![
            Segment { r: 0., z: 0., kind: Kind::Line },
            Segment { r: 1.5, z: 0., kind: Kind::Arc { cr: 0., cz: 0., ccw: true } },
            Segment { r: 0., z: 1.5, kind: Kind::Line },
        ]],
    };
    let tool = gc_kernel::revolve::revolve_general(&dome, TOL).unwrap();
    united(&p, &tool, 40. + 2. / 3. * PI * 1.5f64.powi(3));
    // A post with a cone on top (a cylinder r 1 to 1, then up to the apex at 2).
    let post = Revolution {
        p: v(3., 3., 0.), a: v(0., 0., 1.), x: v(1., 0., 0.), angle: TAU,
        loops: vec![vec![
            Segment { r: 0., z: 0., kind: Kind::Line },
            Segment { r: 1., z: 0., kind: Kind::Line },
            Segment { r: 1., z: 1., kind: Kind::Line },
            Segment { r: 0., z: 2., kind: Kind::Line },
        ]],
    };
    let tool = gc_kernel::revolve::revolve_general(&post, TOL).unwrap();
    united(&p, &tool, 40. + PI + PI / 3.);
}

#[test]
fn bosses_over_and_flush_with_the_edge_and_bars_through_the_plate() {
    let p = plate();
    // Overhanging the right edge: 2 x 3 over x 9..11, 1 high: the overhang's underside is a face of its own.
    let r = united(&p, &up(vec![rect(9., 2., 11., 5., 0.)], 0., 1.), 46.);
    assert_eq!(r.faces.len(), 12);
    // Flush with the right edge: the plate's right side and the boss's are one face (an L), 10 faces.
    let r = united(&p, &up(vec![rect(8., 2., 10., 5., 0.)], 0., 1.), 46.);
    assert_eq!(r.faces.len(), 10);
    assert_eq!(planes(&r), 10);
    // Flush with the corner: two sides joined.
    let r = united(&p, &up(vec![rect(8., 0., 10., 3., 0.)], 0., 1.), 46.);
    assert_eq!(r.faces.len(), 9);
    // A bar through the plate, 1 x 1, from z -2 to 2: it stands out above and below.
    let r = united(&p, &up(vec![rect(5., 4., 6., 5., 0.)], -2., 2.), 40. + 1. * 3.5);
    assert_eq!(r.faces.len(), 6 + 5 + 5);
    // A round bar through a plate edge and past it.
    united(&p, &up(vec![circle(10., 4., 1., 0.)], -1., 1.), 40. + PI * 2. - PI / 2. * 0.5);
    // A boss into the plate, wholly inside it: the union is the plate.
    let r = united(&p, &up(vec![rect(1., 1., 2., 2., 0.)], -0.4, 0.), 40.);
    assert_eq!(r.faces.len(), 6);
    // The plate inside a bigger block: the block.
    united(&p, &up(vec![rect(-1., -1., 11., 9., 0.)], -1., 1.), 12. * 10. * 2.);
    // Two overlapping plates side by side (an L): 10 x 8 and 4 x 4 over x 8..12, y 6..10, same thickness.
    let r = united(&p, &up(vec![rect(8., 6., 12., 10., 0.)], -0.5, 0.), 40. + (16. - 4.) * 0.5);
    assert_eq!(r.faces.len(), 2 + 8);
}

#[test]
fn bosses_with_holes_and_round_corners_through_the_prism_domain() {
    // A plate with a hole; a boss overhanging the hole's edge (half the hole covered from above).
    let holed = up(vec![rect(0., 0., 10., 8., 0.), circle(5., 4., 1., 0.)], -0.5, 0.);
    let base = 40. - PI * 0.5;
    united(&holed, &up(vec![rect(5., 2., 8., 6., 0.)], 0., 1.), base + 12.);
    // A round boss overhanging the plate's corner: a quarter of it rests on the plate.
    united(&plate(), &up(vec![circle(10., 8., 1., 0.)], 0., 0.5), 40. + PI * 0.5);
    // A bar through the hole's wall.
    united(&holed, &up(vec![rect(5.5, 3.5, 7., 4.5, 0.)], -1., 1.), base + 1.5 * 2. - (1.5 - quarter_lens(1., 0.5)) * 0.5);
}

/// The part of the disc x^2 + y^2 < r^2 inside the strip x > a, |y| < a (a < r / sqrt 2), in closed form: full height
/// 2a out to x1 = sqrt(r^2 - a^2), then the circle's own height: 2a (x1 - a) + r^2 (pi / 2 - asin(x1 / r)) - x1 a.
fn quarter_lens(r: f64, a: f64) -> f64 {
    let x1 = (r * r - a * a).sqrt();
    2. * a * (x1 - a) + r * r * (PI / 2. - (x1 / r).asin()) - x1 * a
}

#[test]
fn refusals_are_named() {
    let p = plate();
    let far = up(vec![rect(20., 20., 21., 21., 0.)], 0., 1.);
    assert!(refused(&p, &far).contains("separate solids"));
    // Floating above the plate.
    assert!(refused(&p, &up(vec![rect(4., 2., 6., 5., 0.)], 1., 2.)).contains("separate solids"));
    // Touching only along an edge: a block on top, beyond the plate's right edge.
    assert!(refused(&p, &up(vec![rect(10., 2., 12., 5., 0.)], 0., 1.)).contains("non-manifold"));
    // Beside the plate without touching, on a body that is no prism along the boss (a filleted I-beam).
    let beam = ibeam();
    let boss = up(vec![circle(0., 10., 0.5, 0.)], 6., 7.);
    let r = union(&beam, &boss, TOL).unwrap();
    check(&r, TOL).unwrap();
    // A round boss over the flange's edge: the cylinder is no prism along the beam, so only the face domain could take
    // it, and its base crosses the face's edge.
    let over = up(vec![circle(3., 10., 0.5, 0.)], 6., 7.);
    let m = refused(&beam, &over);
    assert!(m.contains("touches or crosses the edge"), "{m}");
}

/// A filleted I-beam section in the xz plane (flanges 6 wide, 6 deep, 0.5 thick, web 0.3, root fillets 0.3), along y
/// 0..20.
fn ibeam() -> Body {
    let (w, d, tf, tw, r) = (3., 6., 0.5, 0.15, 0.3);
    let p = |x: f64, z: f64| v(x, 0., z);
    let line = |a: Vec3, b: Vec3| Curve::Line { a, b };
    // Root fillet about (cx, cz) from angle t0 by -pi/2 (clockwise in x-z, seen from -y) as an arc about -y... keep the
    // loop counter-clockwise about +y: x to z is clockwise about +y, so the loop is listed with that in mind.
    let arc = |cx: f64, cz: f64, t0: f64| Curve::Arc { c: p(cx, cz), n: v(0., -1., 0.), x: v(1., 0., 0.), r, t0, t1: t0 + PI / 2. };
    // In (x, z) with normal -y the angle runs from +x toward -y x x = ... use explicit points.
    let _ = arc;
    let mut lp = vec![];
    let pts = [(-w, 0.), (w, 0.), (w, tf), (tw + r, tf)];
    for i in 0..3 {
        lp.push(line(p(pts[i].0, pts[i].1), p(pts[i + 1].0, pts[i + 1].1)));
    }
    // Fillet from (tw + r, tf) to (tw, tf + r), centre (tw + r, tf + r).
    lp.push(fillet(p(tw + r, tf + r), p(tw + r, tf), p(tw, tf + r), r));
    lp.push(line(p(tw, tf + r), p(tw, d - tf - r)));
    lp.push(fillet(p(tw + r, d - tf - r), p(tw, d - tf - r), p(tw + r, d - tf), r));
    let pts = [(tw + r, d - tf), (w, d - tf), (w, d), (-w, d), (-w, d - tf), (-tw - r, d - tf)];
    for i in 0..5 {
        lp.push(line(p(pts[i].0, pts[i].1), p(pts[i + 1].0, pts[i + 1].1)));
    }
    lp.push(fillet(p(-tw - r, d - tf - r), p(-tw - r, d - tf), p(-tw, d - tf - r), r));
    lp.push(line(p(-tw, d - tf - r), p(-tw, tf + r)));
    lp.push(fillet(p(-tw - r, tf + r), p(-tw, tf + r), p(-tw - r, tf), r));
    let pts = [(-tw - r, tf), (-w, tf), (-w, 0.)];
    for i in 0..2 {
        lp.push(line(p(pts[i].0, pts[i].1), p(pts[i + 1].0, pts[i + 1].1)));
    }
    let body = prism(vec![lp], v(0., 1., 0.), 0., 20.);
    check(&body, TOL).unwrap();
    body
}

/// A quarter arc about c from a to b, whichever way round is the short one.
fn fillet(c: Vec3, a: Vec3, b: Vec3, r: f64) -> Curve {
    let (da, db) = (a.sub(c), b.sub(c));
    // Turning da to db: about n = da x db / |..|.
    let n = da.cross(db).scale(1. / (r * r));
    Curve::Arc { c, n, x: da.scale(1. / r), r, t0: 0., t1: PI / 2. }
}

#[test]
fn bosses_on_an_ibeam_flange_and_web() {
    let beam = ibeam();
    let vb = volume(&beam);
    // A round boss on the top flange (z 6), clear of the flange's edges.
    united(&beam, &up(vec![circle(1., 10., 0.5, 0.)], 6., 7.), vb + PI * 0.25);
    // A pad on the web's side (x 0.15), along +x, clear of the fillets: 1 x 2 in (y, z), 0.5 out... it would run
    // into the flanges if taller than the web room: here z 2..4, out 0.5 (flange edges at x 3).
    let pad = prism(vec![vec![
        Curve::Line { a: v(0.15, 9., 2.), b: v(0.15, 11., 2.) },
        Curve::Line { a: v(0.15, 11., 2.), b: v(0.15, 11., 4.) },
        Curve::Line { a: v(0.15, 11., 4.), b: v(0.15, 9., 4.) },
        Curve::Line { a: v(0.15, 9., 4.), b: v(0.15, 9., 2.) },
    ]], v(1., 0., 0.), 0., 0.5);
    united(&beam, &pad, vb + 2. * 2. * 0.5);
    // A block over the flange's edge is a prism along the beam like the beam itself: the prism domain unites it.
    united(&beam, &up(vec![rect(2.5, 9., 3.5, 11., 0.)], 6., 7.), vb + 2.);
    // So is a pad on the web reaching past the flanges' edges (it touches only the web).
    let tall = prism(vec![vec![
        Curve::Line { a: v(0.15, 9., 1.), b: v(0.15, 11., 1.) },
        Curve::Line { a: v(0.15, 11., 1.), b: v(0.15, 11., 5.) },
        Curve::Line { a: v(0.15, 11., 5.), b: v(0.15, 9., 5.) },
        Curve::Line { a: v(0.15, 9., 5.), b: v(0.15, 9., 1.) },
    ]], v(1., 0., 0.), 0., 3.5);
    united(&beam, &tall, vb + 2. * 4. * 3.5);
    // A round boss on the web, too big for the room between the flanges: refused by name, the flange is in the way.
    let disc = vec![Curve::Arc { c: v(0.15, 10., 3.), n: v(1., 0., 0.), x: v(0., 1., 0.), r: 2.1, t0: 0., t1: TAU }];
    let m = refused(&beam, &prism(vec![disc], v(1., 0., 0.), 0., 1.));
    assert!(m.contains("room"), "{m}");
}

fn decode(reply: &[f64]) -> (Body, Vec<f64>) {
    let mut r = Reader::new(reply);
    let body = r.body().unwrap();
    // The reader stops after the body; the rest is the envelope.
    let used = write_body(&body).len();
    (body, reply[used..].to_vec())
}

#[test]
fn ops_24_and_25_through_the_protocol() {
    let (p, boss) = (plate(), up(vec![rect(4., 2., 6., 5., 0.)], 0., 1.));
    let mut data = write_body(&p);
    data.extend(write_body(&boss));
    let run = |op: u32| protocol::handle(&protocol::parse_request(&protocol::encode_request(op, TOL, &data)).unwrap());
    let plain = run(24).unwrap();
    let named = run(25).unwrap();
    assert_eq!(&named[..plain.len()], &plain[..], "op 25's body is op 24's");
    let (body, env) = decode(&named);
    assert_eq!(env[0], provenance::START);
    assert_eq!(env[1], provenance::VERSION);
    assert_eq!(env[2] as usize, env.len() - 4);
    assert_eq!(*env.last().unwrap(), provenance::END);
    // The in-process history agrees, and the outcomes say: the plate's top modified, every other plate entity
    // unchanged; the boss's base consumed, its other faces unchanged.
    let (result, history, outcomes) = union_named(&p, &boss, TOL).unwrap();
    assert_eq!(write_body(&result), write_body(&body));
    let top = (0..p.faces.len()).find(|&i| matches!(p.faces[i].surface, Surface::Plane { n, .. } if n.z * if p.faces[i].sense { 1. } else { -1. } > 0.5)).unwrap();
    for (i, o) in outcomes[..p.faces.len()].iter().enumerate() {
        assert_eq!(o.code, if i == top { provenance::MODIFIED } else { provenance::UNCHANGED }, "plate face {i}");
    }
    let tool_faces = &outcomes[p.faces.len() + p.edges.len() + p.vertices.len()..][..boss.faces.len()];
    assert_eq!(tool_faces.iter().filter(|o| o.code == provenance::CONSUMED).count(), 1);
    assert_eq!(tool_faces.iter().filter(|o| o.code == provenance::UNCHANGED).count(), boss.faces.len() - 1);
    assert_eq!(history.faces.len(), result.faces.len());
    let _: &[Outcome] = &outcomes;
    // The stack domain refuses op 25 by name while op 24 builds.
    let over = up(vec![rect(9., 2., 11., 5., 0.)], 0., 1.);
    let mut data2 = write_body(&p);
    data2.extend(write_body(&over));
    let reply = protocol::handle(&protocol::parse_request(&protocol::encode_request(25, TOL, &data2)).unwrap());
    assert!(matches!(reply, Err(KernelError::Unsupported(ref m)) if m.contains("not recorded")), "{reply:?}");
    assert!(protocol::handle(&protocol::parse_request(&protocol::encode_request(24, TOL, &data2)).unwrap()).is_ok());
}

#[test]
fn a_hole_after_a_boss() {
    // The plate with a boss, then a hole through the plate beside it (op 3 on the union's result).
    let once = union(&plate(), &up(vec![rect(4., 2., 6., 5., 0.)], 0., 1.), TOL).unwrap();
    let drill = up(vec![circle(8., 4., 0.5, 0.)], -1., 0.5);
    match cut(&once, &drill, TOL) {
        Ok(r) => {
            check(&r, TOL).unwrap();
            assert!((volume(&r) - (46. - PI * 0.25 * 0.5)).abs() < 1e-9);
        }
        Err(e) => panic!("a hole beside a boss: {e}"),
    }
}

#[test]
fn later_bosses_on_a_stacked_body_overlap_overhang_and_rise_above_earlier_ones() {
    let p = plate();
    let once = united(&p, &up(vec![rect(4., 2., 6., 5., 0.)], 0., 1.), 46.);
    // A second boss across the first, as high: one plus-shaped top, its outline 8 walls.
    let plus = united(&once, &up(vec![rect(5., 3., 8., 4., 0.)], 0., 1.), 46. + 3. - 1.);
    assert_eq!(plus.faces.len(), 2 + 4 + 8 + 1);
    // Taller than the first: it rises out of it.
    united(&once, &up(vec![rect(5., 3., 8., 4., 0.)], 0., 2.), 46. + 6. - 1.);
    // Over the plate's edge after the first boss.
    let edge = united(&once, &up(vec![rect(9., 2., 11., 5., 0.)], 0., 1.), 52.);
    // And a third, joining the two bosses (a bridge between them, 1 high, standing on the plate).
    united(&edge, &up(vec![rect(5., 3., 10., 4., 0.)], 0., 1.), 52. + 5. - 1. - 1.);
    // Two overlapping round bosses r 1, centres 1 apart: the union of the discs (2 pi less the lens) times 0.5.
    let lens = 2. * (0.5f64).acos() - 0.5 * 3f64.sqrt();
    let round = united(&p, &up(vec![circle(4., 4., 1., 0.)], 0., 0.5), 40. + PI * 0.5);
    united(&round, &up(vec![circle(5., 4., 1., 0.)], 0., 0.5), 40. + (2. * PI - lens) * 0.5);
    // A pocketed plate takes a boss over the pocket's edge (the body is a stack of two prisms).
    let pocketed = cut(&p, &up(vec![rect(2., 2., 4., 4., 0.)], -0.2, 0.01), TOL).unwrap();
    united(&pocketed, &up(vec![rect(3., 3., 5., 5., 0.)], 0., 1.), 40. - 4. * 0.2 + 4.);
}

#[test]
fn a_boss_on_a_sloping_face_and_a_boss_over_a_hole() {
    // A bar along y whose section has a 45-degree face from (0, 0) to (2, 2) (x, z): 6 x 10 = 60.
    let section: Vec<Curve> = [(0., 0.), (4., 0.), (4., 2.), (2., 2.)]
        .iter()
        .enumerate()
        .map(|(i, _)| {
            let pts = [(0., 0.), (4., 0.), (4., 2.), (2., 2.)];
            let (a, b) = (pts[i], pts[(i + 1) % 4]);
            Curve::Line { a: v(a.0, 0., a.1), b: v(b.0, 0., b.1) }
        })
        .collect();
    let bar = prism(vec![section], v(0., 1., 0.), 0., 10.);
    assert!((volume(&bar) - 60.).abs() < 1e-9);
    // A 1 x 1 boss on the slope, centred at (1, 5, 1), 0.5 out along its outward normal.
    let (s, n) = (v(1., 0., 1.).scale(0.5f64.sqrt()), v(-1., 0., 1.).scale(0.5f64.sqrt()));
    let y = v(0., 1., 0.);
    let c = v(1., 5., 1.);
    let corner = |a: f64, b: f64| c.add(s.scale(a)).add(y.scale(b));
    let pts = [corner(-0.5, -0.5), corner(0.5, -0.5), corner(0.5, 0.5), corner(-0.5, 0.5)];
    let square: Vec<Curve> = (0..4).map(|i| Curve::Line { a: pts[i], b: pts[(i + 1) % 4] }).collect();
    united(&bar, &prism(vec![square], n, 0., 0.5), 60.5);
    // A boss over a plate's through hole: the hole becomes blind from below (the stack domain).
    let holed = up(vec![rect(0., 0., 10., 8., 0.), circle(5., 4., 0.5, 0.)], -0.5, 0.);
    let r = united(&holed, &up(vec![rect(4., 3., 6., 5., 0.)], 0., 1.), 40. - PI * 0.125 + 4.);
    assert_eq!(r.faces.len(), 6 + 1 + 1 + 5);
}

//! Normalized pipe ends (Claude 2026-10-31, Bob's "short pipe test": a 7 OD x 0.5 wall pipe, 2.5 long, whose right end
//! has Normalize End on over a 13 diameter Preserve Length cope at Cope Angle 90, Rotation 90). GCmaker's rule
//! (PipeCopeNormalizer, EndConditionSideCompositeBuilder): at each angle the deeper of the two walls' crossings with the
//! mating cylinder, carried radially across the wall. Every removed volume is checked against an independent value from
//! the dimensions (a quadrature of the section written here, not the kernel's B-rep or mass code), every result against
//! the topology check, and every display mesh against the exact volume.

use gc_kernel::boolean::cut;
use gc_kernel::geom::{Curve, Surface, Vec3};
use gc_kernel::mass::measure;
use gc_kernel::neutral::{Body, Profile};
use gc_kernel::normalize::{waste_tool, Envelope};
use gc_kernel::topology::check;
use gc_kernel::KernelError;
use std::f64::consts::{PI, TAU};

const TOL: f64 = 1e-6;

/// Adaptive Simpson, to `eps` absolute.
fn simpson(f: &dyn Fn(f64) -> f64, a: f64, b: f64, eps: f64) -> f64 {
    #[allow(clippy::too_many_arguments)]
    fn step(f: &dyn Fn(f64) -> f64, a: f64, b: f64, fa: f64, fm: f64, fb: f64, whole: f64, eps: f64, depth: u32) -> f64 {
        let m = (a + b) / 2.;
        let (lm, rm) = ((a + m) / 2., (m + b) / 2.);
        let (flm, frm) = (f(lm), f(rm));
        let left = (m - a) / 6. * (fa + 4. * flm + fm);
        let right = (b - m) / 6. * (fm + 4. * frm + fb);
        if depth == 0 || (left + right - whole).abs() <= 15. * eps {
            return left + right + (left + right - whole) / 15.;
        }
        step(f, a, m, fa, flm, fm, left, eps / 2., depth - 1) + step(f, m, b, fm, frm, fb, right, eps / 2., depth - 1)
    }
    let (fa, fm, fb) = (f(a), f((a + b) / 2.), f(b));
    step(f, a, b, fa, fm, fb, (b - a) / 6. * (fa + 4. * fm + fb), eps, 48)
}

/// A plain pipe along +Y from y 0 to `length`, as the shape factory places it.
fn pipe(ro: f64, ri: f64, length: f64) -> Body {
    let ring = |r: f64| vec![Curve::Arc { c: Vec3::new(0., 0., 0.), n: Vec3::new(0., -1., 0.), x: Vec3::new(1., 0., 0.), r, t0: 0., t1: TAU }];
    gc_kernel::extrude::extrude(&Profile { loops: vec![ring(ro), ring(ri)], n: Vec3::new(0., 1., 0.), low: 0., high: length }, TOL).unwrap()
}

const RO: f64 = 3.5;
const RI: f64 = 3.0;
const LENGTH: f64 = 2.5;
const RC: f64 = 6.5;

/// Bob's right-end cope: axis along X at y0 = -sqrt(Rc^2 - Ro^2) + 0.0005 (Preserve Length, its inboard nudge).
fn y0() -> f64 { -(RC * RC - RO * RO).sqrt() + 0.0005 }

fn envelope(o: Vec3, d: Vec3, axis_point: Vec3, axis: Vec3, rc: f64) -> Envelope {
    Envelope { o, d, x: Vec3::new(1., 0., 0.), ri: RI, ro: RO, p: axis_point, b: axis, rc }
}

fn right_end() -> Envelope { envelope(Vec3::new(0., 0., 0.), Vec3::new(0., 1., 0.), Vec3::new(0., y0(), 0.), Vec3::new(1., 0., 0.), RC) }

/// The removed volume for a square cope across X: at each angle the inner wall's crossing (its |z| is the smaller, so
/// it lies deeper), h(t) = y0 + sqrt(Rc^2 - Ri^2 sin^2 t) from the end; the waste over the annulus is (Ro^2 - Ri^2) / 2
/// times the integral of h round the turn.
fn square_cope_removed(y0: f64) -> f64 {
    let turn = 4. * simpson(&|t: f64| (RC * RC - RI * RI * t.sin().powi(2)).sqrt(), 0., PI / 2., 1e-15);
    (RO * RO - RI * RI) / 2. * (TAU * y0 + turn)
}

#[test]
fn bobs_short_pipe_right_end_normalized_over_a_square_cope() {
    let body = pipe(RO, RI, LENGTH);
    let (tool, tool_volume) = waste_tool(&right_end(), 1.5, 4.0, -4.5, TOL).unwrap();
    check(&tool, TOL).unwrap();
    // The tool's own volume: (sb^2 - sa^2) / 2 times the integral of (h - far) round the turn, against its B-rep.
    let measured = measure(&tool, TOL).unwrap();
    assert!(!measured.exact);
    assert!((measured.volume - tool_volume).abs() < 1e-9, "tool {} vs {}", measured.volume, tool_volume);
    let turn = 4. * simpson(&|t: f64| (RC * RC - RI * RI * t.sin().powi(2)).sqrt(), 0., PI / 2., 1e-15);
    let band = (16. - 2.25) / 2. * (TAU * (y0() + 4.5) + turn);
    assert!((tool_volume - band).abs() < 1e-9, "tool {tool_volume} vs independent {band}");
    let result = cut(&body, &tool, TOL).unwrap();
    let report = check(&result, TOL).unwrap();
    assert_eq!((report.faces, report.edges, report.vertices), (4, 7, 4));
    assert!(result.faces.iter().any(|f| matches!(f.surface, Surface::Conoid { .. })));
    let removed = measure(&body, TOL).unwrap().volume - measure(&result, TOL).unwrap().volume;
    let exact = square_cope_removed(y0());
    assert!((removed - exact).abs() < 1e-9, "removed {removed}, independent {exact}");
    // The end face's area: the integral over the annulus of sqrt(s^2 + h'^2) / s ds dt (h' from the closed form).
    let dh = |t: f64| -RI * RI * t.sin() * t.cos() / (RC * RC - RI * RI * t.sin().powi(2)).sqrt();
    let face = 4. * simpson(&|t: f64| {
        let g = dh(t);
        let a = |s: f64| s / 2. * (s * s + g * g).sqrt() + if g != 0. { g * g / 2. * (s / g.abs()).asinh() } else { 0. };
        a(RO) - a(RI)
    }, 0., PI / 2., 1e-15);
    let walls = |r: f64| r * (TAU * LENGTH - (TAU * y0() + 4. * simpson(&|t: f64| (RC * RC - RI * RI * t.sin().powi(2)).sqrt(), 0., PI / 2., 1e-15)));
    let far_cap = PI * (RO * RO - RI * RI);
    let area = measure(&result, TOL).unwrap().area;
    assert!((area - (face + walls(RO) + walls(RI) + far_cap)).abs() < 1e-9, "area {area}");
}

#[test]
fn normalized_end_display_mesh_encloses_the_exact_volume() {
    let body = pipe(RO, RI, LENGTH);
    let (tool, _) = waste_tool(&right_end(), 1.5, 4.0, -4.5, TOL).unwrap();
    let result = cut(&body, &tool, TOL).unwrap();
    let exact = measure(&result, TOL).unwrap().volume;
    for chord in [0.01, 0.001] {
        let mesh = gc_kernel::tessellate::tessellate(&result, chord).unwrap();
        let v: f64 = mesh.triangles.iter().map(|t| {
            let (a, b, c) = (mesh.positions[t[0]], mesh.positions[t[1]], mesh.positions[t[2]]);
            a.dot(b.cross(c)) / 6.
        }).sum();
        // Chords sag inward by at most `chord` over the surface: the mesh volume is within area x chord.
        let area = measure(&result, TOL).unwrap().area;
        assert!((v - exact).abs() < area * chord, "chord {chord}: mesh {v}, exact {exact}");
    }
    let mesh = gc_kernel::tessellate::tessellate(&tool, 0.01).unwrap();
    assert!(!mesh.triangles.is_empty());
}

/// One wall of a pipe (ro, ri) pierced square by a round hole of radius a: the integral of
/// 2 sqrt(a^2 - x^2) (sqrt(ro^2 - x^2) - sqrt(ri^2 - x^2)) over (-a, a), with x = a sin t.
fn one_wall(a: f64, ro: f64, ri: f64) -> f64 {
    simpson(&|t: f64| {
        let s = a * t.sin();
        2. * a * a * t.cos() * t.cos() * ((ro * ro - s * s).sqrt() - (ri * ri - s * s).sqrt())
    }, -PI / 2., PI / 2., 1e-15)
}

/// A round cutter of radius rc along unit d, its axis through q, from s0 to s1 along d.
fn cutter(q: Vec3, d: Vec3, rc: f64, s0: f64, s1: f64) -> Body {
    let x = d.perpendicular();
    let start = q.add(d.scale(s0));
    gc_kernel::extrude::extrude(
        &Profile { loops: vec![vec![Curve::Arc { c: start, n: d, x, r: rc, t0: 0., t1: TAU }]], n: d, low: 0., high: s1 - s0 },
        TOL,
    )
    .unwrap()
}

#[test]
fn bobs_top_hole_after_the_normalized_end_cuts_one_wall() {
    // Hole_1: 1 diameter, 1 deep from the top plane (z 3.5) at y 1.25, into the bore; after the end.
    let body = pipe(RO, RI, LENGTH);
    let (tool, _) = waste_tool(&right_end(), 1.5, 4.0, -4.5, TOL).unwrap();
    let ended = cut(&body, &tool, TOL).unwrap();
    let hole = cutter(Vec3::new(0., 1.25, 3.5), Vec3::new(0., 0., -1.), 0.5, -0.001, 1.0);
    let result = cut(&ended, &hole, TOL).unwrap();
    check(&result, TOL).unwrap();
    let removed = measure(&ended, TOL).unwrap().volume - measure(&result, TOL).unwrap().volume;
    let exact = one_wall(0.5, RO, RI);
    assert!((removed - exact).abs() < 1e-10, "removed {removed}, one wall {exact}");
    // GCmaker's final solid for Bob's part (sampled envelope, tessellated): 18.372154 in^3.
    let total = PI * (RO * RO - RI * RI) * LENGTH - square_cope_removed(y0()) - exact;
    assert!((measure(&result, TOL).unwrap().volume - total).abs() < 1e-9);
    gc_kernel::tessellate::tessellate(&result, 0.01).unwrap();
    // A hole reaching the end is refused by name, not cut.
    let low = cutter(Vec3::new(0., 0.5, 3.5), Vec3::new(0., 0., -1.), 0.5, -0.001, 1.0);
    match cut(&ended, &low, TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("normalized"), "{m}"),
        other => panic!("expected a named refusal, got {other:?}"),
    }
}

#[test]
fn left_end_normalized_mirrors_the_right() {
    // The same cope at the far end (y = 2.5, material toward -Y).
    let body = pipe(RO, RI, LENGTH);
    let e = envelope(Vec3::new(0., LENGTH, 0.), Vec3::new(0., -1., 0.), Vec3::new(0., LENGTH - y0(), 0.), Vec3::new(1., 0., 0.), RC);
    let (tool, _) = waste_tool(&e, 1.5, 4.0, -4.5, TOL).unwrap();
    let result = cut(&body, &tool, TOL).unwrap();
    check(&result, TOL).unwrap();
    let removed = measure(&body, TOL).unwrap().volume - measure(&result, TOL).unwrap().volume;
    assert!((removed - square_cope_removed(y0())).abs() < 1e-9, "removed {removed}");
}

/// The larger root y of the line (x, ., z) on the mating cylinder (axis point q, unit axis d, radius rc).
fn y_on(x: f64, z: f64, q: Vec3, d: Vec3, rc: f64) -> f64 {
    let w0 = Vec3::new(x - q.x, -q.y, z - q.z);
    let (a, b, c) = (1. - d.y * d.y, 2. * (w0.y - w0.dot(d) * d.y), w0.dot(w0) - w0.dot(d).powi(2) - rc * rc);
    (-b + (b * b - 4. * a * c).sqrt()) / (2. * a)
}

#[test]
fn turned_copes_follow_the_deeper_wall_and_tilted_ones_that_change_wall_are_refused() {
    // Cope Angle 90 turned 35 about the pipe, Preserve Length on (axis at y0): the removed volume is
    // (Ro^2 - Ri^2) / 2 times the integral round the turn of the deeper crossing max(h_in, h_out) (wall points at
    // (r cos t, -r sin t), as the waste band's frame has them), here always the inner one.
    let turn = 35_f64.to_radians();
    let axis = Vec3::new(turn.sin(), 0., turn.cos());
    let q = Vec3::new(0., y0(), 0.);
    let deeper = |q: Vec3, d: Vec3| move |t: f64| y_on(RO * t.cos(), -RO * t.sin(), q, d, RC).max(y_on(RI * t.cos(), -RI * t.sin(), q, d, RC));
    let e = envelope(Vec3::new(0., 0., 0.), Vec3::new(0., 1., 0.), q, axis, RC);
    let (tool, _) = waste_tool(&e, 1.5, 4.0, -8.0, TOL).unwrap();
    let body = pipe(RO, RI, 12.);
    let result = cut(&body, &tool, TOL).unwrap();
    check(&result, TOL).unwrap();
    let removed = measure(&body, TOL).unwrap().volume - measure(&result, TOL).unwrap().volume;
    let exact = (RO * RO - RI * RI) / 2. * simpson(&deeper(q, axis), 0., TAU, 1e-14);
    assert!((removed - exact).abs() < 1e-9, "removed {removed}, independent {exact}");
    // A turn changes nothing but the angles: the same volume as Bob's square cope.
    assert!((removed - square_cope_removed(y0())).abs() < 1e-9);
    // Off the pipe's centre line (0.4) the outer wall lies deeper where the cope's axis passes between the walls'
    // crossings: refused by name too.
    let off = Vec3::new(0.4 * turn.cos(), y0(), -0.4 * turn.sin());
    match waste_tool(&envelope(Vec3::new(0., 0., 0.), Vec3::new(0., 1., 0.), off, axis, RC), 1.5, 4.0, -8.0, TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("deeper wall changes"), "{m}"),
        other => panic!("expected a named refusal, got {:?}", other.map(|_| ())),
    }
    // Cope Angle 70 (tilt 20 about X), turned 35: the deeper wall changes round the pipe (both signs of h_in - h_out
    // occur, well beyond the tolerance), so the envelope has two pieces: refused by name.
    let tilt = 20_f64.to_radians();
    let axis = Vec3::new(tilt.cos() * turn.sin(), -tilt.sin(), tilt.cos() * turn.cos());
    let o = Vec3::new(0., 0., 0.);
    let gaps: Vec<f64> = (0..4096).map(|i| {
        let t = TAU * i as f64 / 4096.;
        y_on(RI * t.cos(), -RI * t.sin(), o, axis, RC) - y_on(RO * t.cos(), -RO * t.sin(), o, axis, RC)
    }).collect();
    assert!(gaps.iter().any(|g| *g > 1e-3) && gaps.iter().any(|g| *g < -1e-3));
    match waste_tool(&envelope(o, Vec3::new(0., 1., 0.), o, axis, RC), 1.5, 4.0, -8.0, TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("deeper wall changes"), "{m}"),
        other => panic!("expected a named refusal, got {:?}", other.map(|_| ())),
    }
}

#[test]
fn normalized_ends_outside_the_domain_are_refused_by_name() {
    let refused = |r: Result<(Body, f64), KernelError>| match r {
        Err(KernelError::Unsupported(m)) => m,
        other => panic!("expected a named refusal, got {:?}", other.map(|_| ())),
    };
    // A mating cylinder smaller than the pipe misses the walls at its sides.
    let small = envelope(Vec3::new(0., 0., 0.), Vec3::new(0., 1., 0.), Vec3::new(0., 0., 0.), Vec3::new(1., 0., 0.), 3.2);
    assert!(refused(waste_tool(&small, 1.5, 4.0, -4.5, TOL)).contains("all the way round"));
    // A mating cylinder along the pipe.
    let along = envelope(Vec3::new(0., 0., 0.), Vec3::new(0., 1., 0.), Vec3::new(0., 0., 0.), Vec3::new(0., 1., 0.), 8.);
    assert!(refused(waste_tool(&along, 1.5, 4.0, -4.5, TOL)).contains("along the pipe"));
    // An envelope reaching through the pipe: the cut domain refuses it (the end would not be one face).
    let body = pipe(RO, RI, 0.8);
    let (tool, _) = waste_tool(&right_end(), 1.5, 4.0, -4.5, TOL).unwrap();
    match cut(&body, &tool, TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("reaches a pipe end"), "{m}"),
        other => panic!("expected a named refusal, got {other:?}"),
    }
    // Any other cut on a body with a normalized end (here a box) is refused by name.
    let ended = cut(&pipe(RO, RI, LENGTH), &tool, TOL).unwrap();
    let slab = gc_kernel::extrude::extrude(&Profile {
        loops: vec![vec![
            Curve::Line { a: Vec3::new(-5., 2., -5.), b: Vec3::new(5., 2., -5.) },
            Curve::Line { a: Vec3::new(5., 2., -5.), b: Vec3::new(5., 2., 5.) },
            Curve::Line { a: Vec3::new(5., 2., 5.), b: Vec3::new(-5., 2., 5.) },
            Curve::Line { a: Vec3::new(-5., 2., 5.), b: Vec3::new(-5., 2., -5.) },
        ]],
        n: Vec3::new(0., 1., 0.), low: 0., high: 1.,
    }, TOL).unwrap();
    match cut(&ended, &slab, TOL) {
        Err(KernelError::Unsupported(m)) => assert!(m.contains("normalized pipe end"), "{m}"),
        other => panic!("expected a named refusal, got {other:?}"),
    }
}

#[test]
fn op_19_replies_the_waste_band_and_its_volume() {
    use gc_kernel::protocol::{encode_request, handle, parse_request, OP_NORMALIZED_END};
    let e = right_end();
    let data = [e.o.x, e.o.y, e.o.z, e.d.x, e.d.y, e.d.z, e.x.x, e.x.y, e.x.z, e.ri, e.ro, e.p.x, e.p.y, e.p.z, e.b.x, e.b.y, e.b.z, e.rc, 1.5, 4.0, -4.5];
    let reply = handle(&parse_request(&encode_request(OP_NORMALIZED_END, TOL, &data)).unwrap()).unwrap();
    let (body, volume) = (&reply[..reply.len() - 1], reply[reply.len() - 1]);
    let mut reader = gc_kernel::neutral::Reader::with_tolerance(body, TOL);
    let tool = reader.body().unwrap();
    reader.finish().unwrap();
    assert!((measure(&tool, TOL).unwrap().volume - volume).abs() < 1e-9);
    // The neutral round trip keeps the new curve and surface types exactly.
    let again = gc_kernel::neutral::write_body(&tool);
    assert_eq!(again, body.to_vec());
}


/// Preserve every exact curve/surface of the already-normalized cap, even as
/// topology indexes change. This checks geometry retention, not persistent names.
fn retained_cap(before: &Body, after: &Body) {
    let cap=before.faces.iter().find(|f|matches!(f.surface,Surface::Conoid{..})).unwrap();
    let signature=format!("{:?}",cap.surface);
    let kept=after.faces.iter().find(|f|format!("{:?}",f.surface)==signature).unwrap();
    assert_eq!(cap.sense,kept.sense);
    let curves=|body: &Body, face: &gc_kernel::neutral::Face| {
        face.loops.iter().map(|l|l.iter().map(|(e,f)|format!("{f}:{:?}",body.edges[*e].curve)).collect::<Vec<_>>()).collect::<Vec<_>>()
    };
    assert_eq!(curves(before,cap),curves(after,kept));
}
fn opposite_tools(length: f64, offset: f64, turn: f64) -> [Body;2] {
    let theta=turn.to_radians();
    let right=envelope(Vec3::new(0.,0.,0.),Vec3::new(0.,1.,0.),Vec3::new(0.,offset,0.),Vec3::new(1.,0.,0.),RC);
    let left=envelope(Vec3::new(0.,length,0.),Vec3::new(0.,-1.,0.),Vec3::new(0.,length-offset,0.),Vec3::new(theta.sin(),0.,theta.cos()),RC);
    [right,left].map(|e|waste_tool(&e,1.5,4.,-8.,TOL).unwrap().0)
}
fn closed_mesh(body: &Body, expected: f64) {
    let chord=0.002;
    let mesh=gc_kernel::tessellate::tessellate(body,chord).unwrap();
    let bits=|x:f64|if x==0. {0}else{x.to_bits()};
    let key=|p:Vec3|(bits(p.x),bits(p.y),bits(p.z));
    let mut boundary=std::collections::HashMap::<_,i32>::new();
    let mut volume=0.;
    for tri in &mesh.triangles {
        let [a,b,c]=tri.map(|i|mesh.positions[i]);
        volume+=a.dot(b.cross(c))/6.;
        for (p,q) in [(a,b),(b,c),(c,a)] {
            let (p,q)=(key(p),key(q));if p==q {continue;}
            let (edge,sign)=if p<q {((p,q),1)}else{((q,p),-1)};
            *boundary.entry(edge).or_default()+=sign;
        }
    }
    assert!(boundary.values().all(|v|*v==0),"open or inconsistently oriented display mesh");
    assert!(volume>0. && (volume-expected).abs()<measure(body,TOL).unwrap().area*chord);
}
#[test]
fn opposite_normalized_ends_keep_the_first_cap_in_either_order() {
    let length=20.;
    let body=pipe(RO,RI,length);let original=gc_kernel::neutral::write_body(&body);
    for offset in [y0(),0.] {for turn in [0.,35.,90.] {for order in [[0,1],[1,0]] {
        let tools=opposite_tools(length,offset,turn);
        let first=cut(&body,&tools[order[0]],TOL).unwrap();
        let snapshot=gc_kernel::neutral::write_body(&first);
        let result=cut(&first,&tools[order[1]],TOL).unwrap();
        check(&result,TOL).unwrap();
        assert_eq!(result.faces.iter().filter(|f|matches!(f.surface,Surface::Conoid{..})).count(),2);
        retained_cap(&first,&result);
        let expected=PI*(RO*RO-RI*RI)*length-2.*square_cope_removed(offset);
        assert!((measure(&result,TOL).unwrap().volume-expected).abs()<1e-9);
        let wire=gc_kernel::neutral::write_body(&result);
        let mut reader=gc_kernel::neutral::Reader::with_tolerance(&wire,TOL);
        let again=reader.body().unwrap();reader.finish().unwrap();
        // Reader::direction intentionally normalizes directions; a turned
        // mating axis can change by an ulp, so test geometry and structure.
        check(&again,TOL).unwrap();
        let back=gc_kernel::neutral::write_body(&again);
        assert_eq!(wire.len(),back.len());
        assert!(wire.iter().zip(&back).all(|(a,b)|(a-b).abs()<=16.*f64::EPSILON*a.abs().max(1.)));
        assert!((measure(&again,TOL).unwrap().volume-expected).abs()<1e-9);
        closed_mesh(&again,expected);
        assert_eq!(snapshot,gc_kernel::neutral::write_body(&first));
    }}}
    assert_eq!(original,gc_kernel::neutral::write_body(&body));
}
#[test]
fn opposite_normalized_end_rejects_touching_overlap_and_recut() {
    for length in [1.5,2.*(y0()+RC)] {
        let tools=opposite_tools(length,y0(),90.);
        let first=cut(&pipe(RO,RI,length),&tools[0],TOL).unwrap();
        let snapshot=gc_kernel::neutral::write_body(&first);
        let error=cut(&first,&tools[1],TOL).unwrap_err().to_string();
        assert!(error.contains("envelope reaches a pipe end"),"{error}");
        assert_eq!(snapshot,gc_kernel::neutral::write_body(&first));
    }
    let tools=opposite_tools(20.,y0(),35.);
    let first=cut(&pipe(RO,RI,20.),&tools[0],TOL).unwrap();
    assert!(cut(&first,&tools[0],TOL).is_err(),"same-end recut is outside this domain");
}
#[test]
fn radial_hole_after_two_normalized_ends_preserves_both_caps() {
    let tools=opposite_tools(20.,y0(),35.);
    let first=cut(&pipe(RO,RI,20.),&tools[0],TOL).unwrap();
    let ended=cut(&first,&tools[1],TOL).unwrap();
    let hole=cutter(Vec3::new(0.,10.,RO),Vec3::new(0.,0.,-1.),0.5,-0.001,1.);
    let result=cut(&ended,&hole,TOL).unwrap();check(&result,TOL).unwrap();
    let expected=PI*(RO*RO-RI*RI)*20.-2.*square_cope_removed(y0())-one_wall(0.5,RO,RI);
    assert!((measure(&result,TOL).unwrap().volume-expected).abs()<1e-9);
    closed_mesh(&result,expected);
}

#[test]
fn opposite_normalized_ends_support_translated_tilted_pipe_frames() {
    let o=Vec3::new(7.,11.,13.);let d=Vec3::new(0.6,0.8,0.);let x=Vec3::new(0.,0.,1.);
    let length=20.;
    let ring=|r|vec![Curve::Arc {c:o,n:d,x,r,t0:0.,t1:TAU}];
    let body=gc_kernel::extrude::extrude(&Profile {loops:vec![ring(RO),ring(RI)],n:d,low:0.,high:length},TOL).unwrap();
    let left_o=o.add(d.scale(length));
    let axes=[x,d.cross(x)];
    let tools=[Envelope{o,d,x,ri:RI,ro:RO,p:o.add(d.scale(y0())),b:axes[0],rc:RC},
       Envelope{o:left_o,d:d.scale(-1.),x,ri:RI,ro:RO,p:left_o.sub(d.scale(y0())),b:axes[1],rc:RC}]
       .map(|e|waste_tool(&e,1.5,4.,-8.,TOL).unwrap().0);
    for order in [[0,1],[1,0]] {
        let first=cut(&body,&tools[order[0]],TOL).unwrap();
        let result=cut(&first,&tools[order[1]],TOL).unwrap();
        retained_cap(&first,&result);check(&result,TOL).unwrap();
        let expected=PI*(RO*RO-RI*RI)*length-2.*square_cope_removed(y0());
        assert!((measure(&result,TOL).unwrap().volume-expected).abs()<1e-9);
        closed_mesh(&result,expected);
    }
}

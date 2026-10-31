//! DXF reader, fitter and writer (src/dxf). Hand-written files cover every entity the reader takes, block inserts
//! (nested, rotated, mirrored, arrayed), object coordinate systems, binary DXF and text encodings; the writer is
//! checked by reading its output back. GCKERNEL_DXF_CORPUS=<folder> also reads every .dxf under that folder.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use gc_kernel::dxf::{fit, read, write_r12, Curve, Drawing, Geometry, Layer, P2};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_DXF_READ, OP_DXF_WRITE};

const EPS: f64 = 1e-9;

/// An ASCII DXF from (code, value) pairs.
fn dxf(pairs: &[(i32, &str)]) -> Vec<u8> {
    let mut s = String::new();
    for (c, v) in pairs {
        s.push_str(&format!("{c}\n{v}\n"));
    }
    s.into_bytes()
}

fn entities(body: &[(i32, &str)]) -> Vec<u8> {
    let mut all = vec![(0, "SECTION"), (2, "ENTITIES")];
    all.extend_from_slice(body);
    all.extend_from_slice(&[(0, "ENDSEC"), (0, "EOF")]);
    dxf(&all)
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() <= EPS * (1.0 + a.abs().max(b.abs()))
}

fn near_p(a: P2, b: P2) -> bool {
    near(a[0], b[0]) && near(a[1], b[1])
}

fn arc_points(g: &Geometry) -> (P2, P2) {
    match g {
        Geometry::Arc { c, r, start, sweep } => (
            [c[0] + r * start.cos(), c[1] + r * start.sin()],
            [c[0] + r * (start + sweep).cos(), c[1] + r * (start + sweep).sin()],
        ),
        Geometry::Line { a, b } => (*a, *b),
        g => panic!("not a line or arc: {g:?}"),
    }
}

fn geometry(d: &Drawing) -> Vec<&Geometry> {
    d.curves.iter().map(|c| &c.geometry).collect()
}

#[test]
fn lines_arcs_and_circles_read_exactly() {
    let d = read(&entities(&[
        (0, "LINE"), (8, "CUT"), (10, "1"), (20, "2"), (30, "0"), (11, "4"), (21, "6"), (31, "0"),
        (0, "CIRCLE"), (8, "0"), (10, "10"), (20, "0"), (30, "0"), (40, "2.5"),
        (0, "ARC"), (10, "0"), (20, "0"), (40, "3"), (50, "350"), (51, "10"),
    ])).unwrap();
    let g = geometry(&d);
    assert_eq!(g.len(), 3);
    assert_eq!(*g[0], Geometry::Line { a: [1.0, 2.0], b: [4.0, 6.0] });
    assert_eq!(*g[1], Geometry::Circle { c: [10.0, 0.0], r: 2.5 });
    match g[2] {
        Geometry::Arc { c, r, start, sweep } => {
            assert!(near_p(*c, [0.0, 0.0]) && near(*r, 3.0));
            assert!(near(*sweep, 20f64.to_radians()), "arc across 0 degrees sweeps 20, got {}", sweep.to_degrees());
            assert!(near(start.rem_euclid(TAU), 350f64.to_radians()));
        }
        g => panic!("{g:?}"),
    }
    assert_eq!(d.layers[d.curves[0].layer].name, "CUT");
}

#[test]
fn polylines_give_their_lines_and_bulge_arcs() {
    // A 4 x 2 slot: two straight sides and two half-circle ends (bulge 1), closed.
    let d = read(&entities(&[
        (0, "LWPOLYLINE"), (90, "4"), (70, "1"),
        (10, "0"), (20, "0"), (10, "4"), (20, "0"), (42, "1"), (10, "4"), (20, "2"), (10, "0"), (20, "2"), (42, "1"),
    ])).unwrap();
    let g = geometry(&d);
    assert_eq!(g.len(), 4);
    assert_eq!(*g[0], Geometry::Line { a: [0.0, 0.0], b: [4.0, 0.0] });
    let (s, e) = arc_points(g[1]);
    assert!(near_p(s, [4.0, 0.0]) && near_p(e, [4.0, 2.0]));
    match g[1] {
        Geometry::Arc { c, r, sweep, .. } => {
            assert!(near_p(*c, [4.0, 1.0]) && near(*r, 1.0) && near(*sweep, PI), "{:?}", g[1]);
        }
        g => panic!("{g:?}"),
    }
    let (s, e) = arc_points(g[3]);
    assert!(near_p(s, [0.0, 2.0]) && near_p(e, [0.0, 0.0]));
    assert!(d.curves.iter().all(|c| c.group == 0));

    // Old-style POLYLINE with a clockwise (negative bulge) quarter arc, then a 3D polyline seen from above.
    let d = read(&entities(&[
        (0, "POLYLINE"), (66, "1"), (70, "0"), (10, "0"), (20, "0"), (30, "0"),
        (0, "VERTEX"), (10, "0"), (20, "0"), (42, "-0.41421356237309503"),
        (0, "VERTEX"), (10, "1"), (20, "1"),
        (0, "SEQEND"),
        (0, "POLYLINE"), (66, "1"), (70, "8"),
        (0, "VERTEX"), (10, "0"), (20, "0"), (30, "0"), (70, "32"),
        (0, "VERTEX"), (10, "3"), (20, "4"), (30, "5"), (70, "32"),
        (0, "SEQEND"),
    ])).unwrap();
    let g = geometry(&d);
    assert_eq!(g.len(), 2);
    match g[0] {
        Geometry::Arc { c, r, sweep, .. } => assert!(near_p(*c, [1.0, 0.0]) && near(*r, 1.0) && near(*sweep, -FRAC_PI_2)),
        g => panic!("{g:?}"),
    }
    assert_eq!(*g[1], Geometry::Line { a: [0.0, 0.0], b: [3.0, 4.0] });
    assert_eq!(d.curves[1].group, 1);
    assert!(d.warnings.iter().any(|w| w.contains("not in the XY plane")), "{:?}", d.warnings);
}

#[test]
fn object_coordinates_follow_the_arbitrary_axis_algorithm() {
    // Extrusion (0, 0, -1): object X is world -X, object Y is world +Y, so a counter-clockwise object arc turns
    // clockwise seen from above.
    let d = read(&entities(&[
        (0, "ARC"), (10, "5"), (20, "0"), (30, "0"), (40, "1"), (50, "0"), (51, "90"), (210, "0"), (220, "0"), (230, "-1"),
    ])).unwrap();
    match geometry(&d)[0] {
        Geometry::Arc { c, r, start, sweep } => {
            assert!(near_p(*c, [-5.0, 0.0]) && near(*r, 1.0));
            assert!(near(*sweep, -FRAC_PI_2), "sweep {}", sweep.to_degrees());
            assert!(near(start.rem_euclid(TAU), PI), "start {}", start.to_degrees());
        }
        g => panic!("{g:?}"),
    }
    assert!(d.warnings.is_empty(), "{:?}", d.warnings);
}

#[test]
fn ellipses_and_splines_read_exactly_and_fit_within_tolerance() {
    let d = read(&entities(&[
        (0, "ELLIPSE"), (10, "0"), (20, "0"), (30, "0"), (11, "4"), (21, "0"), (31, "0"), (40, "0.5"), (41, "0"), (42, "6.283185307179586"),
        (0, "SPLINE"), (70, "8"), (71, "3"), (72, "8"), (73, "4"),
        (40, "0"), (40, "0"), (40, "0"), (40, "0"), (40, "1"), (40, "1"), (40, "1"), (40, "1"),
        (10, "0"), (20, "0"), (30, "0"), (10, "1"), (20, "2"), (30, "0"), (10, "3"), (20, "-2"), (30, "0"), (10, "4"), (20, "0"), (30, "0"),
    ])).unwrap();
    let g = geometry(&d);
    assert_eq!(*g[0], Geometry::Ellipse { c: [0.0, 0.0], u: [4.0, 0.0], v: [0.0, 2.0], t0: 0.0, t1: TAU });
    assert!(matches!(g[1], Geometry::Spline { degree: 3, .. }));

    for tol in [1e-2, 1e-3, 1e-5] {
        for curve in &g {
            let pieces = fit::lines_and_arcs(curve, tol).unwrap();
            check_fit(curve, &pieces, tol);
        }
    }
}

/// The pieces chain end to end, turn tangent at every joint, and stay within tol of the curve both ways.
fn check_fit(curve: &Geometry, pieces: &[Geometry], tol: f64) {
    let exact = sample(curve, 20000);
    let mut along = Vec::new();
    for w in pieces.windows(2) {
        let (_, e) = arc_points(&w[0]);
        let (s, _) = arc_points(&w[1]);
        assert!(fit::dist(e, s) <= 1e-9, "pieces do not join: {e:?} {s:?}");
        let (t0, t1) = (end_tangent(&w[0], false), end_tangent(&w[1], true));
        assert!(fit::cross(t0, t1).abs() < 1e-6 && fit::dot(t0, t1) > 0.0, "corner between pieces {w:?}");
    }
    for p in pieces {
        along.extend(sample(p, 50));
    }
    let worst = |from: &[P2], to: &[P2]| from.iter()
        .map(|p| to.windows(2).map(|w| seg(*p, w[0], w[1])).fold(f64::MAX, f64::min))
        .fold(0.0, f64::max);
    let there = worst(&along, &exact);
    let back = worst(&exact, &along);
    assert!(there <= tol && back <= tol, "tol {tol}: off by {there} / {back} with {} pieces", pieces.len());
}

fn seg(p: P2, a: P2, b: P2) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let l = fit::dot(ab, ab);
    let s = if l == 0.0 { 0.0 } else { (fit::dot([p[0] - a[0], p[1] - a[1]], ab) / l).clamp(0.0, 1.0) };
    fit::dist(p, [a[0] + ab[0] * s, a[1] + ab[1] * s])
}

fn end_tangent(g: &Geometry, start: bool) -> P2 {
    match g {
        Geometry::Line { a, b } => { let l = fit::dist(*a, *b); [(b[0] - a[0]) / l, (b[1] - a[1]) / l] }
        Geometry::Arc { start: s, sweep, .. } => {
            let t = if start { *s } else { s + sweep };
            let k = sweep.signum();
            [-t.sin() * k, t.cos() * k]
        }
        g => panic!("{g:?}"),
    }
}

fn sample(g: &Geometry, n: usize) -> Vec<P2> {
    (0..=n).map(|i| {
        let s = i as f64 / n as f64;
        match g {
            Geometry::Line { a, b } => [a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s],
            Geometry::Arc { c, r, start, sweep } => { let t = start + sweep * s; [c[0] + r * t.cos(), c[1] + r * t.sin()] }
            Geometry::Circle { c, r } => { let t = TAU * s; [c[0] + r * t.cos(), c[1] + r * t.sin()] }
            Geometry::Ellipse { c, u, v, t0, t1 } => {
                let t = t0 + (t1 - t0) * s;
                [c[0] + u[0] * t.cos() + v[0] * t.sin(), c[1] + u[1] * t.cos() + v[1] * t.sin()]
            }
            Geometry::Spline { degree, knots, points, weights } => {
                let (a, b) = (knots[*degree], knots[knots.len() - degree - 1]);
                let u = a + (b - a) * s;
                let mut span = *degree;
                while span + 1 < points.len() && knots[span + 1] <= u { span += 1; }
                fit::spline_at(*degree, knots, points, weights, span, u).0
            }
            Geometry::Polyline { .. } => unreachable!(),
        }
    }).collect()
}

#[test]
fn a_rational_quarter_circle_spline_fits_one_arc() {
    let w = std::f64::consts::FRAC_1_SQRT_2.to_string();
    let d = read(&entities(&[
        (0, "SPLINE"), (70, "12"), (71, "2"), (72, "6"), (73, "3"),
        (40, "0"), (40, "0"), (40, "0"), (40, "1"), (40, "1"), (40, "1"),
        (41, "1"), (41, &w), (41, "1"),
        (10, "1"), (20, "0"), (30, "0"), (10, "1"), (20, "1"), (30, "0"), (10, "0"), (20, "1"), (30, "0"),
    ])).unwrap();
    let g = geometry(&d);
    let pieces = fit::lines_and_arcs(g[0], 1e-6).unwrap();
    check_fit(g[0], &pieces, 1e-6);
    for p in &pieces {
        if let Geometry::Arc { c, r, .. } = p {
            assert!(fit::dist(*c, [0.0, 0.0]) < 1e-5 && (r - 1.0).abs() < 1e-5, "{p:?}");
        }
    }
}

#[test]
fn fit_point_splines_pass_through_their_points() {
    let d = read(&entities(&[
        (0, "SPLINE"), (70, "8"), (71, "3"), (74, "5"),
        (11, "0"), (21, "0"), (31, "0"), (11, "1"), (21, "1"), (31, "0"), (11, "2"), (21, "0"), (31, "0"),
        (11, "3"), (21, "-1"), (31, "0"), (11, "4"), (21, "0"), (31, "0"),
    ])).unwrap();
    let g = geometry(&d);
    let along = sample(g[0], 4000);
    for q in [[0.0, 0.0], [1.0, 1.0], [2.0, 0.0], [3.0, -1.0], [4.0, 0.0]] {
        let best = along.iter().map(|p| fit::dist(*p, q)).fold(f64::MAX, f64::min);
        assert!(best < 2e-3, "fit point {q:?} missed by {best}");
    }
}

#[test]
fn blocks_nest_rotate_mirror_array_and_inherit_layer_and_colour() {
    let d = read(&dxf(&[
        (0, "SECTION"), (2, "BLOCKS"),
        (0, "BLOCK"), (2, "HOLE"), (70, "0"), (10, "1"), (20, "0"), (30, "0"),
        (0, "CIRCLE"), (8, "0"), (62, "0"), (10, "1"), (20, "0"), (30, "0"), (40, "0.25"),
        (0, "ARC"), (8, "KEEP"), (10, "1"), (20, "0"), (30, "0"), (40, "0.5"), (50, "0"), (51, "90"),
        (0, "ENDBLK"),
        (0, "BLOCK"), (2, "PAIR"), (70, "0"), (10, "0"), (20, "0"), (30, "0"),
        (0, "INSERT"), (8, "0"), (62, "0"), (2, "HOLE"), (10, "0"), (20, "0"), (30, "0"),
        (0, "INSERT"), (8, "0"), (62, "0"), (2, "HOLE"), (10, "10"), (20, "0"), (30, "0"), (41, "-1"),
        (0, "ENDBLK"),
        (0, "BLOCK"), (2, "LOOP"), (70, "0"), (10, "0"), (20, "0"), (30, "0"),
        (0, "INSERT"), (2, "LOOP"), (10, "0"), (20, "0"), (30, "0"),
        (0, "ENDBLK"),
        (0, "ENDSEC"),
        (0, "SECTION"), (2, "ENTITIES"),
        (0, "INSERT"), (8, "HOLES"), (62, "3"), (2, "PAIR"), (10, "100"), (20, "50"), (30, "0"), (50, "90"), (41, "2"), (42, "2"),
        (0, "INSERT"), (2, "HOLE"), (10, "0"), (20, "0"), (30, "0"), (70, "3"), (71, "2"), (44, "5"), (45, "7"),
        (0, "INSERT"), (2, "LOOP"), (10, "0"), (20, "0"), (30, "0"),
        (0, "INSERT"), (2, "MISSING"), (10, "0"), (20, "0"), (30, "0"),
        (0, "ENDSEC"), (0, "EOF"),
    ])).unwrap();
    // PAIR at (100, 50), rotated 90, scaled 2: HOLE's base (1, 0) sits on each inner insert point.
    let circles: Vec<_> = d.curves.iter().filter_map(|c| match c.geometry {
        Geometry::Circle { c: centre, r } => Some((centre, r, c.layer, c.color, c.group)),
        _ => None,
    }).collect();
    assert_eq!(circles.len(), 2 + 6, "{:?}", d.curves);
    assert!(near_p(circles[0].0, [100.0, 50.0]) && near(circles[0].1, 0.5));
    assert!(near_p(circles[1].0, [100.0, 70.0]) && near(circles[1].1, 0.5), "{:?}", circles[1]);
    let holes = d.layers.iter().position(|l| l.name == "HOLES").unwrap();
    assert_eq!((circles[0].2, circles[0].3), (holes, 3), "layer 0 and BYBLOCK take the insert's layer and colour");
    // The array: 3 columns 5 apart, 2 rows 7 apart, first hole at the origin.
    let mut grid: Vec<P2> = circles[2..].iter().map(|c| c.0).collect();
    grid.sort_by(|a, b| a[1].total_cmp(&b[1]).then(a[0].total_cmp(&b[0])));
    let expected = [[0.0, 0.0], [5.0, 0.0], [10.0, 0.0], [0.0, 7.0], [5.0, 7.0], [10.0, 7.0]];
    assert!(grid.iter().zip(expected).all(|(a, b)| near_p(*a, b)), "{grid:?}");
    assert!(circles[2..].iter().all(|c| c.4 == 1));

    // The mirrored HOLE's arc turns the other way; its own layer is kept.
    let arcs: Vec<_> = d.curves.iter().filter(|c| matches!(c.geometry, Geometry::Arc { .. })).collect();
    let keep = d.layers.iter().position(|l| l.name == "KEEP").unwrap();
    assert!(arcs.iter().all(|a| a.layer == keep));
    let sweeps: Vec<f64> = arcs[..2].iter().map(|a| match a.geometry { Geometry::Arc { sweep, .. } => sweep, _ => 0.0 }).collect();
    assert!(near(sweeps[0], FRAC_PI_2) && near(sweeps[1], -FRAC_PI_2), "{sweeps:?}");

    assert!(d.warnings.iter().any(|w| w.contains("'LOOP' contains itself")), "{:?}", d.warnings);
    assert!(d.warnings.iter().any(|w| w.contains("'MISSING' is inserted but not defined")), "{:?}", d.warnings);
}

#[test]
fn a_skewing_insert_turns_a_circle_into_an_ellipse() {
    let d = read(&dxf(&[
        (0, "SECTION"), (2, "BLOCKS"),
        (0, "BLOCK"), (2, "C"), (10, "0"), (20, "0"), (30, "0"),
        (0, "CIRCLE"), (10, "0"), (20, "0"), (30, "0"), (40, "1"),
        (0, "ENDBLK"), (0, "ENDSEC"),
        (0, "SECTION"), (2, "ENTITIES"),
        (0, "INSERT"), (2, "C"), (10, "0"), (20, "0"), (30, "0"), (41, "3"), (42, "1"),
        (0, "ENDSEC"), (0, "EOF"),
    ])).unwrap();
    assert_eq!(*geometry(&d)[0], Geometry::Ellipse { c: [0.0, 0.0], u: [3.0, 0.0], v: [0.0, 1.0], t0: 0.0, t1: TAU });
}

#[test]
fn header_layers_and_what_is_not_a_curve() {
    let d = read(&dxf(&[
        (999, "written by hand"),
        (0, "SECTION"), (2, "HEADER"), (9, "$ACADVER"), (1, "AC1015"), (9, "$INSUNITS"), (70, "4"), (9, "$MEASUREMENT"), (70, "1"),
        (0, "ENDSEC"),
        (0, "SECTION"), (2, "TABLES"), (0, "TABLE"), (2, "LAYER"), (70, "2"),
        (0, "LAYER"), (2, "Frozen"), (70, "1"), (62, "5"), (6, "CONTINUOUS"),
        (0, "LAYER"), (2, "Off"), (70, "0"), (62, "-2"), (6, "CONTINUOUS"),
        (0, "ENDTAB"), (0, "ENDSEC"),
        (0, "SECTION"), (2, "ENTITIES"),
        (0, "TEXT"), (8, "0"), (10, "0"), (20, "0"), (40, "1"), (1, "hello"),
        (0, "TEXT"), (8, "0"), (10, "0"), (20, "0"), (40, "1"), (1, "again"),
        (0, "HATCH"), (8, "0"),
        (0, "LINE"), (8, "frozen"), (10, "0"), (20, "0"), (11, "1"), (21, "0"),
        (0, "LINE"), (67, "1"), (10, "0"), (20, "0"), (11, "1"), (21, "0"),
        (0, "LINE"), (60, "1"), (10, "0"), (20, "0"), (11, "1"), (21, "0"),
        (0, "LINE"), (10, "2"), (20, "2"), (11, "2"), (21, "2"),
        (0, "ENDSEC"), (0, "EOF"),
    ])).unwrap();
    assert_eq!((d.version.as_str(), d.insunits, d.measurement), ("AC1015", Some(4), Some(1)));
    let frozen = &d.layers[d.curves[0].layer];
    assert!(frozen.name == "Frozen" && frozen.frozen && !frozen.off && frozen.color == 5, "layer names ignore case");
    assert!(d.layers.iter().any(|l| l.name == "Off" && l.off && l.color == 2));
    assert_eq!(d.curves.len(), 1, "paper space, invisible and zero-length lines are left out");
    assert_eq!(d.skipped, vec![("TEXT".to_string(), 2), ("HATCH".to_string(), 1)]);
    assert!(d.warnings.iter().any(|w| w.contains("paper-space")) && d.warnings.iter().any(|w| w.contains("zero-length")));
}

#[test]
fn line_endings_byte_order_mark_and_code_pages() {
    let lf = entities(&[(0, "LINE"), (8, "Schnitt\u{e4}"), (10, " 1.5 "), (20, "2"), (11, "3"), (21, "4E+00")]);
    let text = String::from_utf8(lf.clone()).unwrap();
    let crlf = text.replace('\n', "\r\n").into_bytes();
    let cr = text.replace('\n', "\r").into_bytes();
    let mut bom = b"\xEF\xBB\xBF".to_vec();
    bom.extend_from_slice(&crlf);
    // Windows-1252: a-umlaut is one byte (E4) and the euro sign is 80.
    let ansi = String::from_utf8(crlf.clone()).unwrap().replace("Schnitt\u{e4}", "SchnittX");
    let mut ansi = ansi.into_bytes();
    let at = ansi.windows(8).position(|w| w == b"SchnittX").unwrap() + 7;
    ansi[at] = 0xE4;
    let escaped = text.replace("Schnitt\u{e4}", "Schnitt\\U+00E4").into_bytes();
    for file in [lf, crlf, cr, bom, ansi, escaped] {
        let d = read(&file).unwrap();
        assert_eq!(*geometry(&d)[0], Geometry::Line { a: [1.5, 2.0], b: [3.0, 4.0] });
        assert_eq!(d.layers[d.curves[0].layer].name, "Schnitt\u{e4}");
    }
    assert!(read(b"0\nSECTION\n2\nENTITIES\nX\nLINE\n").unwrap_err().to_string().contains("line 5"));
    assert!(read(b"hello world").unwrap_err().to_string().contains("not a DXF"));
}

/// A binary DXF of the given pairs (R13+ two-byte group codes, or R12 one-byte codes).
fn binary(pairs: &[(i32, &str)], wide: bool) -> Vec<u8> {
    let mut out = b"AutoCAD Binary DXF\r\n\x1a\0".to_vec();
    for (code, value) in pairs {
        if wide { out.extend((*code as i16).to_le_bytes()); }
        else if *code < 255 { out.push(*code as u8); }
        else { out.push(255); out.extend((*code as i16).to_le_bytes()); }
        match code {
            10..=59 => out.extend(value.parse::<f64>().unwrap().to_le_bytes()),
            60..=79 | 170..=179 | 270..=289 => out.extend(value.parse::<i16>().unwrap().to_le_bytes()),
            90..=99 => out.extend(value.parse::<i32>().unwrap().to_le_bytes()),
            290..=299 => out.push(value.parse::<u8>().unwrap()),
            _ => { out.extend(value.as_bytes()); out.push(0); }
        }
    }
    out
}

#[test]
fn binary_dxf_reads_like_ascii() {
    let pairs = [
        (0, "SECTION"), (2, "HEADER"), (9, "$INSUNITS"), (70, "1"), (0, "ENDSEC"),
        (0, "SECTION"), (2, "ENTITIES"),
        (0, "LWPOLYLINE"), (8, "A"), (90, "3"), (70, "1"), (10, "0"), (20, "0"), (10, "2"), (20, "0"), (42, "0.5"), (10, "2"), (20, "2"),
        (0, "CIRCLE"), (8, "B"), (10, "5"), (20, "5"), (30, "0"), (40, "1"), (290, "1"),
        (0, "ENDSEC"), (0, "EOF"),
    ];
    let ascii = read(&dxf(&pairs)).unwrap();
    for wide in [true, false] {
        let b = read(&binary(&pairs, wide)).unwrap();
        assert_eq!(b.curves, ascii.curves);
        assert_eq!(b.insunits, Some(1));
    }
}

#[test]
fn what_is_written_reads_back() {
    let mut d = Drawing { insunits: Some(1), ..Drawing::default() };
    d.layers.push(Layer { name: "Cut Outside".into(), color: 1, frozen: false, off: false, locked: false });
    d.layers.push(Layer { name: "0".into(), color: 7, frozen: false, off: true, locked: true });
    let curves = [
        Geometry::Line { a: [0.0, 0.0], b: [10.125, 0.0] },
        Geometry::Arc { c: [5.0, 5.0], r: 2.0, start: 0.3, sweep: 1.2 },
        Geometry::Arc { c: [5.0, 5.0], r: 3.0, start: 2.0, sweep: -2.5 },
        Geometry::Circle { c: [1.0 / 3.0, -2.0], r: 0.0625 },
        Geometry::Polyline { points: vec![[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [0.0, 2.0]], bulges: vec![0.0, 1.0, 0.0, -0.5], closed: true },
        Geometry::Ellipse { c: [20.0, 0.0], u: [3.0, 1.0], v: [-0.5, 1.5], t0: 0.5, t1: 4.0 },
    ];
    for (i, g) in curves.iter().enumerate() {
        d.curves.push(Curve { geometry: g.clone(), layer: i % 2, color: if i == 0 { 5 } else { 256 }, group: 0 });
    }
    let tol = 1e-4;
    let text = write_r12(&d, tol).unwrap();
    assert!(text.starts_with("  0\r\nSECTION\r\n  2\r\nHEADER\r\n  9\r\n$ACADVER\r\n  1\r\nAC1009\r\n"));
    let back = read(text.as_bytes()).unwrap();
    assert_eq!((back.version.as_str(), back.insunits), ("AC1009", Some(1)));
    assert_eq!(back.layers.iter().map(|l| l.name.as_str()).collect::<Vec<_>>(), ["CUT_OUTSIDE", "0"]);
    assert!(back.layers[1].off && back.layers[1].locked && !back.layers[0].off);
    let g = geometry(&back);
    assert_eq!(*g[0], curves[0]);
    assert_eq!(back.curves[0].color, 5);
    assert!(same_arc(g[1], &curves[1]) && same_arc(g[2], &curves[2]), "{:?} {:?}", g[1], g[2]);
    assert_eq!(*g[3], curves[3], "numbers are written so they read back to the same value");
    // The polyline: line, half-circle bulge, line, clockwise bulge closing it.
    assert_eq!(*g[4], Geometry::Line { a: [0.0, 0.0], b: [4.0, 0.0] });
    let (s, e) = arc_points(g[5]);
    assert!(near_p(s, [4.0, 0.0]) && near_p(e, [4.0, 2.0]));
    let (s, e) = arc_points(g[7]);
    assert!(near_p(s, [0.0, 2.0]) && near_p(e, [0.0, 0.0]));
    assert!(matches!(g[7], Geometry::Arc { sweep, .. } if *sweep < 0.0));
    // The ellipse came back as its fitted lines and arcs.
    let fitted: Vec<Geometry> = g[8..].iter().map(|g| (*g).clone()).collect();
    check_fit(&curves[5], &fitted, tol);
}

fn same_arc(a: &Geometry, b: &Geometry) -> bool {
    let (sa, ea) = arc_points(a);
    let (sb, eb) = arc_points(b);
    match (a, b) {
        (Geometry::Arc { c: ca, r: ra, sweep: wa, .. }, Geometry::Arc { c: cb, r: rb, sweep: wb, .. }) => {
            let forward = near_p(sa, sb) && near_p(ea, eb) && near(*wa, *wb);
            let reversed = near_p(sa, eb) && near_p(ea, sb) && near(*wa, -*wb);
            near_p(*ca, *cb) && near(*ra, *rb) && (forward || reversed)
        }
        _ => false,
    }
}

fn text_values(s: &str) -> Vec<f64> {
    let mut v = vec![s.len() as f64];
    v.extend(s.bytes().map(f64::from));
    v
}

#[test]
fn worker_ops_read_and_write_files() {
    let dir = std::env::temp_dir().join(format!("gckernel-dxf-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("plate.dxf");
    let path_text = path.to_str().unwrap();

    // Write: one layer, a 10 x 5 closed polyline and a circle, inches.
    let mut req = vec![1.0, 1.0, 1e-4];
    req.extend(text_values(path_text));
    req.extend([1.0, 7.0, 0.0]);
    req.extend(text_values("PLATE"));
    req.extend([2.0, 5.0, 0.0, 256.0, 1.0, 4.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 10.0, 5.0, 0.0, 0.0, 5.0, 0.0]);
    req.extend([2.0, 0.0, 256.0, 5.0, 2.5, 1.0]);
    let reply = handle(&parse_request(&encode_request(OP_DXF_WRITE, 1e-6, &req)).unwrap()).unwrap();
    assert_eq!(reply[0] as usize, std::fs::metadata(&path).unwrap().len() as usize);
    assert!(!dir.join("plate.dxf.writing").exists());

    // Read it back as lines and arcs.
    let mut req = vec![1.0, 1.0, 1e-4];
    req.extend(text_values(path_text));
    let r = handle(&parse_request(&encode_request(OP_DXF_READ, 1e-6, &req)).unwrap()).unwrap();
    assert_eq!(&r[..3], &[1.0, 1.0, -1.0]);
    let version_len = r[3] as usize;
    let version: String = r[4..4 + version_len].iter().map(|b| *b as u8 as char).collect();
    assert_eq!(version, "AC1009");
    let mut at = 4 + version_len;
    let layers = r[at] as usize;
    at += 1;
    for _ in 0..layers {
        at += 2;
        at += 1 + r[at] as usize;
    }
    assert_eq!(r[at], 5.0, "four polyline sides and the circle");
    at += 1;
    let kinds: Vec<f64> = (0..5).map(|_| {
        let k = r[at];
        at += 4 + match k as i32 { 0 => 4, 2 => 3, _ => panic!("kind {k}") };
        k
    }).collect();
    assert_eq!(kinds, [0.0, 0.0, 0.0, 0.0, 2.0]);
    assert_eq!(&r[at..], &[0.0, 0.0], "nothing skipped, no warnings");

    let mut req = vec![1.0, 0.0, 1e-4];
    req.extend(text_values(dir.join("missing.dxf").to_str().unwrap()));
    let err = handle(&parse_request(&encode_request(OP_DXF_READ, 1e-6, &req)).unwrap()).unwrap_err();
    assert!(err.to_string().contains("could not read"), "{err}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn corpus_reads_and_fits() {
    let Ok(root) = std::env::var("GCKERNEL_DXF_CORPUS") else { return };
    let mut files = Vec::new();
    let mut stack = vec![std::path::PathBuf::from(root)];
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir).unwrap().flatten() {
            let p = entry.path();
            if p.is_dir() { stack.push(p); }
            else if p.extension().is_some_and(|e| e.eq_ignore_ascii_case("dxf")) { files.push(p); }
        }
    }
    files.sort();
    for f in &files {
        let d = read(&std::fs::read(f).unwrap()).unwrap_or_else(|e| panic!("{}: {e}", f.display()));
        let mut fitted = 0;
        for c in &d.curves {
            fitted += fit::lines_and_arcs(&c.geometry, 1e-3).unwrap_or_else(|e| panic!("{}: {e}", f.display())).len();
        }
        let mut kinds = std::collections::BTreeMap::new();
        for c in &d.curves {
            let k = match c.geometry { Geometry::Line { .. } => "line", Geometry::Arc { .. } => "arc",
                Geometry::Circle { .. } => "circle", Geometry::Ellipse { .. } => "ellipse",
                Geometry::Spline { .. } => "spline", Geometry::Polyline { .. } => "polyline" };
            *kinds.entry(k).or_insert(0) += 1;
        }
        println!("{} | {} | units {:?} | {kinds:?} -> {fitted} fitted | skipped {:?} | {:?}",
            f.display(), d.version, d.insunits, d.skipped, d.warnings);
    }
    assert!(!files.is_empty());
}

#[test]
fn fit_point_splines_leave_and_arrive_along_their_tangents() {
    let d = read(&entities(&[
        (0, "SPLINE"), (70, "8"), (71, "3"), (74, "3"),
        (12, "0"), (22, "1"), (32, "0"), (13, "1"), (23, "0"), (33, "0"),
        (11, "0"), (21, "0"), (31, "0"), (11, "2"), (21, "2"), (31, "0"), (11, "4"), (21, "1"), (31, "0"),
    ])).unwrap();
    let g = geometry(&d);
    let Geometry::Spline { degree, knots, points, weights } = g[0] else { panic!("{:?}", g[0]) };
    let (p0, d0) = fit::spline_at(*degree, knots, points, weights, *degree, knots[*degree]);
    let last = points.len() - 1;
    let (p1, d1) = fit::spline_at(*degree, knots, points, weights, last, knots[last + 1]);
    assert!(near_p(p0, [0.0, 0.0]) && near_p(p1, [4.0, 1.0]));
    assert!(d0[0].abs() < 1e-9 && d0[1] > 0.0, "leaves straight up: {d0:?}");
    assert!(d1[1].abs() < 1e-9 && d1[0] > 0.0, "arrives heading +X: {d1:?}");
    let along = sample(g[0], 4000);
    assert!(along.iter().map(|p| fit::dist(*p, [2.0, 2.0])).fold(f64::MAX, f64::min) < 2e-3);
}

#[test]
fn header_only_read_skips_the_drawing() {
    let file = dxf(&[
        (0, "SECTION"), (2, "HEADER"), (9, "$INSUNITS"), (70, "4"), (0, "ENDSEC"),
        (0, "SECTION"), (2, "ENTITIES"), (0, "LINE"), (10, "0"), (20, "0"), (11, "1"), (21, "0"), (0, "ENDSEC"), (0, "EOF"),
    ]);
    let d = gc_kernel::dxf::read_header(&file).unwrap();
    assert_eq!((d.insunits, d.curves.len()), (Some(4), 0));
}

/// Damaged files give an error or a drawing, never a panic: every byte of a real-looking file flipped, cut short or
/// repeated, ASCII and binary.
#[test]
fn damaged_files_never_panic() {
    let pairs = [
        (0, "SECTION"), (2, "BLOCKS"), (0, "BLOCK"), (2, "B"), (10, "0"), (20, "0"), (30, "0"),
        (0, "ARC"), (10, "1"), (20, "0"), (30, "0"), (40, "1"), (50, "0"), (51, "90"), (0, "ENDBLK"), (0, "ENDSEC"),
        (0, "SECTION"), (2, "ENTITIES"),
        (0, "LWPOLYLINE"), (90, "3"), (70, "1"), (10, "0"), (20, "0"), (42, "0.5"), (10, "2"), (20, "0"), (10, "2"), (20, "2"),
        (0, "SPLINE"), (71, "3"), (72, "8"), (73, "4"), (40, "0"), (40, "0"), (40, "0"), (40, "0"), (40, "1"), (40, "1"), (40, "1"), (40, "1"),
        (10, "0"), (20, "0"), (30, "0"), (10, "1"), (20, "2"), (30, "0"), (10, "3"), (20, "-2"), (30, "0"), (10, "4"), (20, "0"), (30, "0"),
        (0, "ELLIPSE"), (10, "0"), (20, "0"), (30, "0"), (11, "4"), (21, "0"), (31, "0"), (40, "0.5"),
        (0, "INSERT"), (2, "B"), (10, "0"), (20, "0"), (30, "0"), (41, "2"), (70, "2"), (71, "2"), (44, "3"), (45, "3"),
        (0, "ENDSEC"), (0, "EOF"),
    ];
    let mut seed = 0x2545F4914F6CDD1Du64;
    let mut next = || { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; seed };
    for base in [dxf(&pairs), binary(&pairs, true), binary(&pairs, false)] {
        for cut in (0..base.len()).step_by(3) {
            try_read(&base[..cut]);
        }
        for _ in 0..3000 {
            let mut f = base.clone();
            for _ in 0..1 + next() % 4 {
                let at = (next() as usize) % f.len();
                match next() % 4 {
                    0 => f[at] = next() as u8,
                    1 => f[at] = b"0123456789-.eE \n"[(next() % 16) as usize],
                    2 => { f.remove(at); }
                    _ => { let b = f[at]; f.insert(at, b); }
                }
            }
            try_read(&f);
        }
    }
}

fn try_read(bytes: &[u8]) {
    let owned = bytes.to_vec();
    let outcome = std::panic::catch_unwind(move || {
        if let Ok(d) = read(&owned) {
            for c in &d.curves {
                let _ = fit::lines_and_arcs(&c.geometry, 1e-3);
            }
            let _ = write_r12(&d, 1e-3);
        }
    });
    assert!(outcome.is_ok(), "panicked on {:?}", String::from_utf8_lossy(bytes));
}

#[test]
fn a_hundred_thousand_entities_read_in_well_under_a_second() {
    let mut s = String::from("0\nSECTION\n2\nENTITIES\n");
    for i in 0..100_000 {
        let x = i as f64 * 0.01;
        s.push_str(&format!("0\nLINE\n8\nCUT\n10\n{x}\n20\n0\n30\n0\n11\n{x}\n21\n1\n31\n0\n"));
    }
    s.push_str("0\nENDSEC\n0\nEOF\n");
    let start = std::time::Instant::now();
    let d = read(s.as_bytes()).unwrap();
    let took = start.elapsed();
    assert_eq!(d.curves.len(), 100_000);
    assert!(took.as_secs_f64() < 1.0, "took {took:?}");
}

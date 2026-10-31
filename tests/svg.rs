//! SVG reader (src/dxf/svg.rs): units, the y flip, every path command, shapes, transforms, use/symbol, hidden content.

use std::f64::consts::PI;

use gc_kernel::dxf::{fit, svg, Drawing, Geometry, P2};

fn read(s: &str) -> Drawing {
    svg::read(s.as_bytes()).unwrap_or_else(|e| panic!("{e}"))
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

fn near_p(a: P2, b: P2) -> bool {
    near(a[0], b[0]) && near(a[1], b[1])
}

fn ends(g: &Geometry) -> (P2, P2) {
    match g {
        Geometry::Line { a, b } => (*a, *b),
        Geometry::Arc { c, r, start, sweep } => (
            [c[0] + r * start.cos(), c[1] + r * start.sin()],
            [c[0] + r * (start + sweep).cos(), c[1] + r * (start + sweep).sin()],
        ),
        Geometry::Spline { points, .. } => (points[0], *points.last().unwrap()),
        Geometry::Ellipse { c, u, v, t0, t1 } => {
            let at = |t: f64| [c[0] + u[0] * t.cos() + v[0] * t.sin(), c[1] + u[1] * t.cos() + v[1] * t.sin()];
            (at(*t0), at(*t1))
        }
        g => panic!("{g:?}"),
    }
}

/// Consecutive pieces of one path meet end to start.
fn chained(gs: &[&Geometry]) -> bool {
    gs.windows(2).all(|w| fit::dist(ends(w[0]).1, ends(w[1]).0) < 1e-9)
}

#[test]
fn a_millimetre_plate_with_a_rounded_rectangle_and_holes() {
    let d = read(r#"<svg xmlns="http://www.w3.org/2000/svg" width="100mm" height="60mm" viewBox="0 0 100 60">
  <rect x="2" y="2" width="96" height="56" rx="8" fill="none" stroke="black"/>
  <circle cx="25" cy="30" r="5"/>
  <circle cx="75" cy="20" r="5"/>
</svg>"#);
    assert_eq!(d.insunits, Some(4));
    let g: Vec<&Geometry> = d.curves.iter().map(|c| &c.geometry).collect();
    assert_eq!(g.len(), 8 + 2);
    assert!(chained(&g[..8]) && fit::dist(ends(g[7]).1, ends(g[0]).0) < 1e-9, "the outline closes");
    let arcs: Vec<f64> = g[..8].iter().filter_map(|g| match g { Geometry::Arc { r, .. } => Some(*r), _ => None }).collect();
    assert_eq!(arcs.len(), 4);
    assert!(arcs.iter().all(|r| near(*r, 8.0)));
    // y turned up: SVG y 20 on a 60 tall page is 40.
    assert_eq!(*g[9], Geometry::Circle { c: [75.0, 40.0], r: 5.0 });
    let lowest = g[..8].iter().flat_map(|g| [ends(g).0[1], ends(g).1[1]]).fold(f64::MAX, f64::min);
    assert!(near(lowest, 2.0));
}

#[test]
fn every_path_command_in_compact_syntax() {
    // Absolute and relative, implicit repeats, "-" and "." running numbers together, arc flags without spaces.
    let d = read(r#"<svg width="10in" height="10in" viewBox="0 0 100 100"><path d="
        M10-0l10.5.5H40V10h-5v5
        C40 20 50 20 50 15S60 10 60 15
        Q65 20 70 15T80 15
        a5 5 0 00 10 0
        A5 5 0 1 1 90 30z
        m5 5 1 1"/></svg>"#);
    let g: Vec<&Geometry> = d.curves.iter().map(|c| &c.geometry).collect();
    let kinds: Vec<&str> = g.iter().map(|g| match g {
        Geometry::Line { .. } => "L", Geometry::Arc { .. } => "A", Geometry::Spline { degree: 3, .. } => "C",
        Geometry::Spline { degree: 2, .. } => "Q", _ => "?",
    }).collect();
    assert_eq!(kinds.join(""), "LLLLLCCQQAALL");
    assert!(chained(&g[..12]), "one continuous path until the second subpath");
    // Units: 100 user units across 10 in = 254 mm.
    let mm = 254.0 / 100.0;
    assert!(near_p(ends(g[0]).0, [10.0 * mm, 100.0 * mm]));
    assert!(near_p(ends(g[0]).1, [20.5 * mm, (100.0 - 0.5) * mm]), "l10.5.5 is 10.5 then 0.5: {:?}", ends(g[0]));
    // The smooth cubic reflects the previous control point.
    if let Geometry::Spline { points, .. } = g[6] {
        assert!(near_p(points[1], [50.0 * mm, (100.0 - 10.0) * mm]), "{points:?}");
    }
    // "a5 5 0 00 10 0": half circle from (80, 15) to (90, 15), centre (85, 15).
    match g[9] {
        Geometry::Arc { c, r, sweep, .. } => {
            assert!(near_p(*c, [85.0 * mm, 85.0 * mm]) && near(*r, 5.0 * mm), "{c:?} {r}");
            assert!(near(sweep.abs(), PI));
        }
        g => panic!("{g:?}"),
    }
    // The z closes the path back to (10, 0).
    assert!(near_p(ends(g[11]).1, [10.0 * mm, 100.0 * mm]));
}

#[test]
fn transforms_groups_use_symbol_and_hidden_content() {
    let d = read(r##"<svg xmlns:xlink="http://www.w3.org/1999/xlink" width="200mm" height="100mm" viewBox="0 0 200 100">
  <defs><circle id="hole" r="2"/></defs>
  <symbol id="pair"><use href="#hole" x="0"/><use href="#hole" x="10"/></symbol>
  <g transform="translate(50 50) rotate(90)"><line x1="0" y1="0" x2="10" y2="0"/></g>
  <g transform="matrix(2 0 0 2 0 0)"><use xlink:href="#hole" x="10" y="10"/></g>
  <use href="#pair" transform="translate(100,20) scale(2,1)"/>
  <g style="display:none"><circle cx="1" cy="1" r="1"/></g>
  <circle cx="1" cy="1" r="1" visibility="hidden"/>
  <use href="#missing"/>
  <text x="0" y="0">label</text>
</svg>"##);
    let g: Vec<&Geometry> = d.curves.iter().map(|c| &c.geometry).collect();
    assert_eq!(g.len(), 4, "{g:?}");
    // rotate(90) in y-down turns +x into +y (down the page), so the line runs from page (50, 50) to (50, 60).
    assert_eq!(*g[0], Geometry::Line { a: [50.0, 50.0], b: [50.0, 40.0] });
    assert_eq!(*g[1], Geometry::Circle { c: [20.0, 80.0], r: 4.0 });
    // scale(2, 1) stretches the pair's circles into ellipses 4 wide and 2 tall.
    for (geometry, x) in [(g[2], 100.0), (g[3], 120.0)] {
        match geometry {
            Geometry::Ellipse { c, u, v, .. } => {
                assert!(near_p(*c, [x, 80.0]), "{c:?}");
                assert!(near(fit::len(*u).max(fit::len(*v)), 4.0) && near(fit::len(*u).min(fit::len(*v)), 2.0));
            }
            g => panic!("{g:?}"),
        }
    }
    assert!(d.warnings.iter().any(|w| w.contains("#missing")));
    assert_eq!(d.skipped, vec![("TEXT".to_string(), 1)]);
}

#[test]
fn pixel_files_read_at_96_per_inch() {
    let d = read(r#"<?xml version="1.0"?><!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "x.dtd" [<!ENTITY a "b">]>
<!-- drawn in pixels --><svg width="288" height="192"><polygon points="0,0 96,0 96 48"/></svg>"#);
    assert_eq!(d.insunits, Some(1));
    assert!(d.warnings.iter().any(|w| w.contains("96 pixels per inch")));
    let g: Vec<&Geometry> = d.curves.iter().map(|c| &c.geometry).collect();
    assert_eq!(*g[0], Geometry::Line { a: [0.0, 2.0], b: [1.0, 2.0] });
    assert_eq!(g.len(), 3);
}

#[test]
fn view_box_aspect_is_kept_and_centred() {
    // A 100 x 100 box into a 200 x 100 mm page: scale 1, centred with 50 mm each side.
    let d = read(r#"<svg width="200mm" height="100mm" viewBox="0 0 100 100"><circle cx="50" cy="50" r="10"/></svg>"#);
    assert_eq!(d.curves[0].geometry, Geometry::Circle { c: [100.0, 50.0], r: 10.0 });
    let d = read(r#"<svg width="200mm" height="100mm" viewBox="0 0 100 100" preserveAspectRatio="none"><circle cx="50" cy="50" r="10"/></svg>"#);
    assert!(matches!(d.curves[0].geometry, Geometry::Ellipse { .. }));
}

#[test]
fn arcs_with_radii_too_small_grow_to_reach() {
    let d = read(r#"<svg width="100mm" height="100mm" viewBox="0 0 100 100"><path d="M0 50 A1 1 0 0 1 20 50"/></svg>"#);
    match d.curves[0].geometry {
        Geometry::Arc { c, r, sweep, .. } => assert!(near_p(c, [10.0, 50.0]) && near(r, 10.0) && near(sweep.abs(), PI)),
        ref g => panic!("{g:?}"),
    }
}

#[test]
fn damaged_svgs_never_panic() {
    let base = br##"<svg width="100mm" height="60mm" viewBox="0 0 100 60"><defs><circle id="h" r="2"/></defs>
<g transform="rotate(30 50 30) scale(1.5)"><rect x="2" y="2" width="96" height="56" rx="8"/>
<path d="M10 10 C20 20 30 0 40 10 S60 20 70 10 Q75 5 80 10 T90 10 A5 3 30 1 0 95 20 Z"/>
<use href="#h" x="5" y="5"/><polyline points="1 2 3 4 5 6"/><ellipse cx="5" cy="5" rx="3" ry="1"/></g></svg>"##;
    let mut seed = 0x9E3779B97F4A7C15u64;
    let mut next = || { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; seed };
    for cut in 0..base.len() {
        let _ = svg::read(&base[..cut]);
    }
    for _ in 0..5000 {
        let mut f = base.to_vec();
        for _ in 0..1 + next() % 4 {
            let at = (next() as usize) % f.len();
            match next() % 4 {
                0 => f[at] = next() as u8,
                1 => f[at] = b"0123456789-.eE ,<>/\"#"[(next() % 21) as usize],
                2 => { f.remove(at); }
                _ => { let b = f[at]; f.insert(at, b); }
            }
        }
        let outcome = std::panic::catch_unwind(|| {
            if let Ok(d) = svg::read(&f) {
                for c in &d.curves { let _ = fit::lines_and_arcs(&c.geometry, 1e-3); }
            }
        });
        assert!(outcome.is_ok(), "panicked on {}", String::from_utf8_lossy(&f));
    }
    assert!(svg::read(b"<html></html>").unwrap_err().to_string().contains("not an SVG"));
}

#[test]
fn corpus_svgs_read() {
    let Ok(root) = std::env::var("GCKERNEL_SVG_CORPUS") else { return };
    for entry in std::fs::read_dir(root).unwrap().flatten() {
        let p = entry.path();
        if !p.extension().is_some_and(|e| e.eq_ignore_ascii_case("svg")) { continue; }
        let d = svg::read(&std::fs::read(&p).unwrap()).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
        let mut fitted = 0;
        for c in &d.curves {
            fitted += fit::lines_and_arcs(&c.geometry, 1e-3).unwrap().len();
        }
        println!("{} | units {:?} | {} curves -> {fitted} fitted | layers {:?} | skipped {:?} | {:?}", p.display(),
            d.insunits, d.curves.len(), d.layers.iter().map(|l| l.name.as_str()).collect::<Vec<_>>(), d.skipped, d.warnings);
    }
}

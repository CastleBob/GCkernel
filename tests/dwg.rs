//! DWG reader (src/dxf/dwg). GCKERNEL_DWG_FIXTURES=<folder> compares every `*.dwg` there with its `.dxf` twin (the
//! same drawing saved by another program as DXF): the same curves, matched as a multiset by kind and geometry within
//! 1e-9 relative, on the same layers. The other tests read the sample drawings in tests/fixtures/dwg.

use std::f64::consts::TAU;
use std::path::{Path, PathBuf};

use gc_kernel::dxf::{dwg, fit, read, write_r12, Curve, Drawing, Geometry};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_DXF_READ};

fn near(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= 1e-9 * (1.0 + scale.max(a.abs()).max(b.abs()))
}

fn near_all(a: &[f64], b: &[f64], scale: f64) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| near(*x, *y, scale))
}

fn angle_near(a: f64, b: f64) -> bool {
    let d = (a - b).rem_euclid(TAU);
    d.min(TAU - d) <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

/// The same curve: equal kind and parameters within 1e-9 relative to the drawing's size.
fn same(a: &Geometry, b: &Geometry, scale: f64) -> bool {
    match (a, b) {
        (Geometry::Line { a: a0, b: a1 }, Geometry::Line { a: b0, b: b1 }) => near_all(&[a0[0], a0[1], a1[0], a1[1]], &[b0[0], b0[1], b1[0], b1[1]], scale),
        (Geometry::Circle { c: c0, r: r0 }, Geometry::Circle { c: c1, r: r1 }) => near_all(&[c0[0], c0[1], *r0], &[c1[0], c1[1], *r1], scale),
        (Geometry::Arc { c: c0, r: r0, start: s0, sweep: w0 }, Geometry::Arc { c: c1, r: r1, start: s1, sweep: w1 }) => {
            near_all(&[c0[0], c0[1], *r0], &[c1[0], c1[1], *r1], scale) && angle_near(*s0, *s1) && near(*w0, *w1, 1.0)
        }
        (Geometry::Ellipse { c: c0, u: u0, v: v0, t0: a0, t1: a1 }, Geometry::Ellipse { c: c1, u: u1, v: v1, t0: b0, t1: b1 }) => {
            near_all(&[c0[0], c0[1], u0[0], u0[1], v0[0], v0[1]], &[c1[0], c1[1], u1[0], u1[1], v1[0], v1[1]], scale)
                && angle_near(*a0, *b0) && near(a1 - a0, b1 - b0, 1.0)
        }
        (Geometry::Spline { degree: d0, knots: k0, points: p0, weights: w0 }, Geometry::Spline { degree: d1, knots: k1, points: p1, weights: w1 }) => {
            let flat = |p: &[[f64; 2]]| p.iter().flat_map(|q| [q[0], q[1]]).collect::<Vec<_>>();
            d0 == d1 && near_all(k0, k1, 1.0) && near_all(&flat(p0), &flat(p1), scale) && near_all(w0, w1, 1.0)
        }
        _ => false,
    }
}

fn kind(g: &Geometry) -> &'static str {
    match g {
        Geometry::Line { .. } => "line",
        Geometry::Arc { .. } => "arc",
        Geometry::Circle { .. } => "circle",
        Geometry::Ellipse { .. } => "ellipse",
        Geometry::Spline { .. } => "spline",
        Geometry::Polyline { .. } => "polyline",
    }
}

fn size(d: &Drawing) -> f64 {
    let mut s: f64 = 0.0;
    for c in &d.curves {
        let v: Vec<f64> = match &c.geometry {
            Geometry::Line { a, b } => vec![a[0], a[1], b[0], b[1]],
            Geometry::Arc { c, r, .. } | Geometry::Circle { c, r } => vec![c[0], c[1], *r],
            Geometry::Ellipse { c, u, v, .. } => vec![c[0], c[1], u[0], u[1], v[0], v[1]],
            Geometry::Spline { points, .. } => points.iter().flat_map(|p| [p[0], p[1]]).collect(),
            Geometry::Polyline { points, .. } => points.iter().flat_map(|p| [p[0], p[1]]).collect(),
        };
        s = v.iter().fold(s, |m, x| m.max(x.abs()));
    }
    s
}

fn layer(d: &Drawing, c: &Curve) -> String {
    d.layers.get(c.layer).map(|l| l.name.to_ascii_uppercase()).unwrap_or_default()
}

/// Differences between the DWG and DXF readings, in words; empty when the curves match.
fn compare(from_dwg: &Drawing, from_dxf: &Drawing) -> Vec<String> {
    let mut problems = Vec::new();
    let scale = size(from_dxf);
    let mut used = vec![false; from_dwg.curves.len()];
    let mut unmatched = Vec::new();
    for c in &from_dxf.curves {
        let found = from_dwg.curves.iter().enumerate()
            .position(|(i, d)| !used[i] && same(&d.geometry, &c.geometry, scale) && layer(from_dwg, d) == layer(from_dxf, c));
        match found {
            Some(i) => used[i] = true,
            None => unmatched.push(c),
        }
    }
    if from_dwg.curves.len() != from_dxf.curves.len() {
        problems.push(format!("{} curves from the DWG, {} from the DXF", from_dwg.curves.len(), from_dxf.curves.len()));
    }
    for c in unmatched.iter().take(5) {
        problems.push(format!("DXF curve not in the DWG: {} on '{}': {:?}", kind(&c.geometry), layer(from_dxf, c), c.geometry));
    }
    for (i, c) in from_dwg.curves.iter().enumerate().filter(|(i, _)| !used[*i]).take(5) {
        problems.push(format!("DWG curve {i} not in the DXF: {} on '{}': {:?}", kind(&c.geometry), layer(from_dwg, c), c.geometry));
    }
    // The layer table: layers in both files have the same visibility flags (a converter may add layers of
    // its own to the DXF, such as "Default", and drop names DXF does not allow).
    for m in &from_dwg.layers {
        match from_dxf.layers.iter().find(|l| l.name.eq_ignore_ascii_case(&m.name)) {
            None => {} // converters rename or drop some (e.g. *IDW_BlockReferenceLayer)
            Some(l) if (m.frozen, m.off, m.locked) != (l.frozen, l.off, l.locked) => problems.push(format!(
                "layer '{}': DWG frozen/off/locked {:?}, DXF {:?}", l.name, (m.frozen, m.off, m.locked), (l.frozen, l.off, l.locked))),
            _ => {}
        }
    }
    problems
}

/// Layer colours that differ (reported, not failed: a true-colour layer has no exact colour index).
fn colour_notes(from_dwg: &Drawing, from_dxf: &Drawing) -> Vec<String> {
    from_dxf.layers.iter().filter_map(|l| {
        let m = from_dwg.layers.iter().find(|m| m.name.eq_ignore_ascii_case(&l.name))?;
        (m.color != l.color).then(|| format!("{} {}/{}", l.name, m.color, l.color))
    }).collect()
}

fn fixtures() -> Option<Vec<(PathBuf, PathBuf)>> {
    let dir = std::env::var("GCKERNEL_DWG_FIXTURES").ok()?;
    let mut pairs: Vec<(PathBuf, PathBuf)> = std::fs::read_dir(&dir).ok()?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("dwg")))
        .filter_map(|p| { let twin = p.with_extension("dxf"); twin.exists().then_some((p, twin)) })
        .collect();
    pairs.sort();
    Some(pairs)
}

fn name(p: &Path) -> String {
    p.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default()
}

#[test]
fn every_dwg_reads_like_its_dxf_twin() {
    let Some(pairs) = fixtures() else { return };
    assert!(!pairs.is_empty(), "GCKERNEL_DWG_FIXTURES has no .dwg files with .dxf twins");
    let mut failed = Vec::new();
    for (dwg_path, dxf_path) in &pairs {
        let from_dxf = read(&std::fs::read(dxf_path).unwrap()).unwrap_or_else(|e| panic!("{}: {e}", name(dxf_path)));
        let from_dwg = match dwg::read(&std::fs::read(dwg_path).unwrap()) {
            Ok(d) => d,
            Err(e) => {
                println!("{}: FAIL {e}", name(dwg_path));
                failed.push(name(dwg_path));
                continue;
            }
        };
        let problems = compare(&from_dwg, &from_dxf);
        let skipped = |d: &Drawing| { let mut s = d.skipped.clone(); s.sort(); s };
        let mut note = if skipped(&from_dwg) != skipped(&from_dxf) {
            format!(" (skipped differ: DWG {:?}, DXF {:?})", skipped(&from_dwg), skipped(&from_dxf))
        } else {
            String::new()
        };
        let colours = colour_notes(&from_dwg, &from_dxf);
        if !colours.is_empty() {
            note += &format!(" (layer colours DWG/DXF differ: {})", colours.join(", "));
        }
        if problems.is_empty() {
            println!("{}: pass, {} {} curves{note}", name(dwg_path), from_dwg.version, from_dwg.curves.len());
        } else {
            println!("{}: FAIL {}{note}", name(dwg_path), from_dwg.version);
            for p in &problems {
                println!("    {p}");
            }
            failed.push(name(dwg_path));
        }
    }
    assert!(failed.is_empty(), "these DWGs do not match their DXF twins: {failed:?}");
}

#[test]
fn versions_are_reported() {
    let Some(pairs) = fixtures() else { return };
    for (dwg_path, _) in &pairs {
        let bytes = std::fs::read(dwg_path).unwrap();
        let d = dwg::read_header(&bytes).unwrap_or_else(|e| panic!("{}: {e}", name(dwg_path)));
        assert_eq!(d.version.as_bytes(), &bytes[0..6], "{}", name(dwg_path));
        let expected = [("R14", "AC1014"), ("2000", "AC1015"), ("2004", "AC1018"), ("2007", "AC1021"), ("2010", "AC1024"),
            ("2013", "AC1027"), ("2018", "AC1032")];
        for (tag, version) in expected {
            if name(dwg_path) == format!("fixture_{tag}.dwg") {
                assert_eq!(d.version, version);
            }
        }
        println!("{}: {} INSUNITS {:?} MEASUREMENT {:?} layers {}", name(dwg_path), d.version, d.insunits, d.measurement, d.layers.len());
    }
}

#[test]
fn dwg_headers_are_told_from_dxf_and_svg() {
    assert!(dwg::looks_like(b"AC1032\0\0\0\0\0\0"));
    assert!(dwg::looks_like(b"AC1015rest"));
    assert!(dwg::looks_like(b"AC1009"));
    assert!(!dwg::looks_like(b"  0\nSECTION\n  2\nHEADER\n"));
    assert!(!dwg::looks_like(b"0\r\nSECTION\r\n"));
    assert!(!dwg::looks_like(b"<?xml version=\"1.0\"?><svg/>"));
    assert!(!dwg::looks_like(b"AutoCAD Binary DXF\r\n\x1a\0"));
    assert!(!dwg::looks_like(b"AC10"));
    assert!(matches!(dwg::read(b"AC1009\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"), Err(gc_kernel::KernelError::Unsupported(_))));
    assert!(dwg::read(b"AC1032").is_err());
    assert!(dwg::read(b"").is_err());
}

/// The repository's sample DWGs (tests/fixtures/dwg): one drawing saved as R14, 2004, 2007 and 2018, so every file
/// organisation (R13-R15, R2004, R2007, R2010+) is covered without the fixtures folder.
fn samples() -> Vec<(String, Vec<u8>)> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/dwg");
    ["fixture_R14.dwg", "fixture_2004.dwg", "fixture_2007.dwg", "fixture_2018.dwg"].iter()
        .map(|n| (n.to_string(), std::fs::read(dir.join(n)).unwrap_or_else(|e| panic!("{n}: {e}"))))
        .collect()
}

/// Each sample: two layers' worth of lines, arcs, a circle, ellipses, a spline, a closed polyline and nested rotated
/// block inserts (18 curves), a TEXT left out, millimetres.
#[test]
fn the_sample_drawing_reads_in_every_version() {
    for (name, bytes) in samples() {
        let d = dwg::read(&bytes).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(d.curves.len(), 18, "{name}");
        assert_eq!(d.skipped, vec![("TEXT".to_string(), 1)], "{name}");
        let count = |k: &str| d.curves.iter().filter(|c| kind(&c.geometry) == k).count();
        assert_eq!((count("line"), count("arc"), count("circle"), count("ellipse"), count("spline")), (9, 2, 4, 2, 1), "{name}");
        for layer in ["CUT", "HOLES"] {
            assert!(d.layers.iter().any(|l| l.name.eq_ignore_ascii_case(layer)), "{name}: no layer {layer}");
        }
        assert!(d.warnings.is_empty(), "{name}: {:?}", d.warnings);
        assert_eq!(d.measurement, Some(1), "{name}");
        let expected_units = if name.contains("R14") { None } else { Some(6) };
        assert_eq!(d.insunits, expected_units, "{name}");
        let header = dwg::read_header(&bytes).unwrap();
        assert!(header.curves.is_empty() && header.layers.len() == d.layers.len(), "{name}");
    }
}

/// Damaged files give an error or a drawing, never a panic: truncations and random byte changes of each sample.
#[test]
fn damaged_files_never_panic() {
    let start = std::time::Instant::now();
    let mut seed = 0x9E3779B97F4A7C15u64;
    let mut next = || { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; seed };
    for (k, (name, base)) in samples().into_iter().enumerate() {
        let step = (base.len() / 300).max(1);
        for cut in (0..base.len()).step_by(step) {
            try_read(&base[..cut]);
        }
        // About six seconds per sample for random changes.
        let budget = 6.0 * (k + 1) as f64;
        let mut rounds = 0;
        while rounds < 3000 && start.elapsed().as_secs_f64() < budget {
            let mut f = base.clone();
            for _ in 0..1 + next() % 6 {
                let at = (next() as usize) % f.len();
                match next() % 3 {
                    0 => f[at] = next() as u8,
                    1 => f[at] ^= 1 << (next() % 8),
                    _ => f[at] = [0x00, 0xff, 0x7f, 0x80][(next() % 4) as usize],
                }
            }
            try_read(&f);
            rounds += 1;
        }
        println!("{name}: {rounds} damaged copies read without a panic");
    }
}

fn try_read(bytes: &[u8]) {
    let owned = bytes.to_vec();
    let outcome = std::panic::catch_unwind(move || {
        if let Ok(d) = dwg::read(&owned) {
            for c in &d.curves {
                let _ = fit::lines_and_arcs(&c.geometry, 1e-3);
            }
            let _ = write_r12(&d, 1e-3);
        }
        let _ = dwg::read_header(&owned);
    });
    assert!(outcome.is_ok(), "panicked on a damaged DWG of {} bytes", bytes.len());
}

/// Op 14 tells a DWG by its content: mode 1 gives lines and arcs, mode 2 the header and layers only.
#[test]
fn the_read_op_takes_dwg_files() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/dwg/fixture_2018.dwg");
    let path = path.to_str().unwrap();
    let request = |mode: f64| {
        let mut req = vec![1.0, mode, 1e-3, path.len() as f64];
        req.extend(path.bytes().map(f64::from));
        handle(&parse_request(&encode_request(OP_DXF_READ, 1e-6, &req)).unwrap()).unwrap()
    };
    for mode in [1.0, 2.0] {
        let r = request(mode);
        assert_eq!(&r[..3], &[1.0, 6.0, 1.0], "version 1, millimetres, metric");
        let version: String = r[4..4 + r[3] as usize].iter().map(|b| *b as u8 as char).collect();
        assert_eq!(version, "AC1032");
        let mut at = 4 + r[3] as usize;
        let layers = r[at] as usize;
        at += 1;
        for _ in 0..layers {
            at += 2;
            at += 1 + r[at] as usize;
        }
        if mode == 2.0 {
            assert_eq!(r[at], 0.0, "no curves in a header read");
        } else {
            assert!(r[at] >= 18.0, "ellipses and the spline become several lines and arcs");
        }
    }
}

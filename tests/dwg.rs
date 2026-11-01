//! DWG reader (src/dxf/dwg). GCKERNEL_DWG_FIXTURES=<folder> compares every `*.dwg` there with its `.dxf` twin (the
//! same drawing saved by another program as DXF): the same curves, matched as a multiset by kind and geometry within
//! 1e-9 relative, on the same layers with the same colours, and the same layer colours and flags. The other tests
//! read the sample drawings in tests/fixtures/dwg.

use std::f64::consts::TAU;
use std::path::{Path, PathBuf};

use gc_kernel::dxf::{dwg, fit, read, write_r12, Curve, Drawing, Geometry, TextKind};
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
            .position(|(i, d)| !used[i] && same(&d.geometry, &c.geometry, scale) && layer(from_dwg, d) == layer(from_dxf, c)
                && d.color == c.color);
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
    // The layer table: layers in both files have the same visibility flags and colour index (a true colour
    // becomes the nearest index in both; a converter may add layers of its own to the DXF, such as "Default", and
    // drop names DXF does not allow).
    for m in &from_dwg.layers {
        match from_dxf.layers.iter().find(|l| l.name.eq_ignore_ascii_case(&m.name)) {
            None => {} // converters rename or drop some (e.g. *IDW_BlockReferenceLayer)
            Some(l) if (m.frozen, m.off, m.locked) != (l.frozen, l.off, l.locked) => problems.push(format!(
                "layer '{}': DWG frozen/off/locked {:?}, DXF {:?}", l.name, (m.frozen, m.off, m.locked), (l.frozen, l.off, l.locked))),
            Some(l) if m.color != l.color => problems.push(format!(
                "layer '{}': DWG colour {}, DXF colour {}", l.name, m.color, l.color)),
            _ => {}
        }
    }
    problems
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
        let note = if skipped(&from_dwg) != skipped(&from_dxf) {
            format!(" (skipped differ: DWG {:?}, DXF {:?})", skipped(&from_dwg), skipped(&from_dxf))
        } else {
            String::new()
        };
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
    // R12 and older (AC1009 and before, down to AutoCAD 1.x and 2.x) are refused by name, whole files or not.
    for old in ["AC1009", "AC1006", "AC1004", "AC1003", "AC1002", "AC1001", "AC2.10", "AC1.50"] {
        let mut bytes = old.as_bytes().to_vec();
        assert!(dwg::looks_like(&bytes), "{old}");
        bytes.extend_from_slice(&[0u8; 0x400]);
        for result in [dwg::read(&bytes), dwg::read_header(&bytes), dwg::read(old.as_bytes())] {
            match result {
                Err(gc_kernel::KernelError::Unsupported(m)) => assert!(m.contains(old) && m.contains("R12 or earlier")
                    && m.contains("save the drawing as DWG R13 or later, or as DXF"), "{m}"),
                other => panic!("{old}: {other:?}"),
            }
        }
    }
    assert!(matches!(dwg::read(b"AC1040\0\0\0\0"), Err(gc_kernel::KernelError::Unsupported(m)) if m.contains("newer")));
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
/// GCKERNEL_DWG_FUZZ=<n> runs n times as long (and from another seed) for a deeper search.
#[test]
fn damaged_files_never_panic() {
    let start = std::time::Instant::now();
    let scale = std::env::var("GCKERNEL_DWG_FUZZ").ok().and_then(|v| v.parse::<u64>().ok()).unwrap_or(1).max(1);
    let mut seed = 0x9E3779B97F4A7C15u64 ^ scale.wrapping_mul(0x2545F4914F6CDD1D);
    let mut next = || { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; seed };
    for (k, (name, base)) in samples().into_iter().enumerate() {
        let step = (base.len() / 300).max(1);
        for cut in (0..base.len()).step_by(step) {
            try_read(&base[..cut]);
        }
        // About six seconds per sample for random changes.
        let budget = 6.0 * scale as f64 * (k + 1) as f64;
        let mut rounds = 0;
        while rounds < 3000 * scale && start.elapsed().as_secs_f64() < budget {
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

/// Where the R2007 sample (tests/fixtures/dwg/fixture_2007.dwg) keeps its Reed-Solomon code words: (file offset,
/// interleaved code words, parity bytes per word, what). Read off the file's header, page map and section map.
const R2007_HEADER: (usize, usize) = (0x80, 3);
const R2007_PAGE_MAP: (usize, usize) = (0x480, 4); // 9 copies of its 104 bytes
const R2007_SECTION_MAP: (usize, usize) = (0x4A60, 17); // 5 copies
const R2007_OBJECTS_PAGE: (usize, usize) = (0x1160, 28); // AcDb:AcDbObjects, (255, 251), encoding 4
const R2007_HEADER_PAGE: (usize, usize) = (0x4460, 6); // AcDb:Header, (255, 251), encoding 4

fn sample_2007() -> Vec<u8> {
    std::fs::read(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/dwg/fixture_2007.dwg")).unwrap()
}

/// Changes `errors` bytes (spread over the 255, parity included) of interleaved code word `word` of the `blocks`
/// words at `base`.
fn damage_word(f: &mut [u8], (base, blocks): (usize, usize), word: usize, errors: usize) {
    for m in 0..errors {
        let i = (word * 37 + m * 31 + 5) % 255;
        f[base + i * blocks + word] ^= 0x5A ^ m as u8;
    }
}

fn damage_all(f: &mut [u8], at: (usize, usize), errors: usize) {
    for word in 0..at.1 {
        damage_word(f, at, word, errors);
    }
}

/// The end-of-file copy of the file header, made unusable, so the header at 0x80 must be repaired.
fn spoil_header_copy(f: &mut [u8]) {
    let copy = (f.len() - 0x400, 3);
    damage_all(f, copy, 9);
}

/// Up to 8 byte errors per code word in the file header and system pages, 2 in data pages: the read equals the
/// undamaged read.
#[test]
fn r2007_byte_errors_are_corrected() {
    let base = sample_2007();
    let good = dwg::read(&base).unwrap();
    for errors in 1..=8 {
        let mut f = base.clone();
        spoil_header_copy(&mut f);
        damage_all(&mut f, R2007_HEADER, errors);
        assert_eq!(dwg::read(&f).unwrap(), good, "file header, {errors} errors per word");
        let mut f = base.clone();
        damage_all(&mut f, R2007_PAGE_MAP, errors);
        assert_eq!(dwg::read(&f).unwrap(), good, "page map, {errors} errors per word");
        let mut f = base.clone();
        damage_all(&mut f, R2007_SECTION_MAP, errors);
        assert_eq!(dwg::read(&f).unwrap(), good, "section map, {errors} errors per word");
    }
    for errors in 1..=2 {
        let mut f = base.clone();
        damage_all(&mut f, R2007_OBJECTS_PAGE, errors);
        damage_all(&mut f, R2007_HEADER_PAGE, errors);
        assert_eq!(dwg::read(&f).unwrap(), good, "data pages, {errors} errors per word");
    }
    // All at once: 8 per system word, 2 per data word.
    let mut f = base.clone();
    spoil_header_copy(&mut f);
    for at in [R2007_HEADER, R2007_PAGE_MAP, R2007_SECTION_MAP] {
        damage_all(&mut f, at, 8);
    }
    for at in [R2007_OBJECTS_PAGE, R2007_HEADER_PAGE] {
        damage_all(&mut f, at, 2);
    }
    assert_eq!(dwg::read(&f).unwrap(), good, "everywhere at once");
    // A copy beyond repair is passed over: the header's end-of-file copy, a later copy of the page map.
    let mut f = base.clone();
    damage_all(&mut f, R2007_HEADER, 9);
    damage_word(&mut f, R2007_PAGE_MAP, 0, 12);
    assert_eq!(dwg::read(&f).unwrap(), good, "second copies");
    // Random places and values: up to 8 errors in each system code word, up to 2 in each data code word.
    let mut seed = 0x6A09E667F3BCC909u64;
    let mut next = move || { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; seed };
    for round in 0..200 {
        let mut f = base.clone();
        spoil_header_copy(&mut f);
        for (at, most) in [(R2007_HEADER, 8), (R2007_PAGE_MAP, 8), (R2007_SECTION_MAP, 8), (R2007_OBJECTS_PAGE, 2),
            (R2007_HEADER_PAGE, 2)] {
            for word in 0..at.1 {
                let mut places: Vec<usize> = Vec::new();
                let errors = (next() % (most + 1)) as usize;
                while places.len() < errors {
                    let i = (next() % 255) as usize;
                    if !places.contains(&i) {
                        places.push(i);
                    }
                }
                for i in places {
                    f[at.0 + i * at.1 + word] ^= 1 + (next() % 255) as u8;
                }
            }
        }
        assert_eq!(dwg::read(&f).unwrap(), good, "random round {round}");
    }
}

/// More errors in a code word than its parity can correct (and no good copy): an error naming the damage.
#[test]
fn r2007_damage_beyond_repair_fails_by_name() {
    let base = sample_2007();
    let fails = |f: &[u8], what: &str| match dwg::read(f) {
        Err(e) => {
            let message = e.to_string();
            assert!(message.contains("the DWG 2007 file is damaged beyond repair at"), "{what}: {message}");
            assert!(message.contains(what), "{what}: {message}");
            println!("{message}");
        }
        Ok(_) => panic!("{what}: read a file damaged beyond repair"),
    };
    for errors in [9, 12, 40] {
        let mut f = base.clone();
        spoil_header_copy(&mut f);
        damage_word(&mut f, R2007_HEADER, 0, errors);
        fails(&f, "file header");
        let mut f = base.clone();
        damage_all(&mut f, R2007_PAGE_MAP, errors);
        fails(&f, "page map");
        let mut f = base.clone();
        damage_all(&mut f, R2007_SECTION_MAP, errors);
        fails(&f, "section map");
    }
    let mut f = base.clone();
    damage_all(&mut f, R2007_OBJECTS_PAGE, 3);
    fails(&f, "'AcDb:AcDbObjects' data page");
}

/// The error message for a read that must fail as damaged beyond repair (never a drawing, never another error).
fn beyond_repair(f: &[u8], what: &str) -> String {
    match dwg::read(f) {
        Err(e) => {
            let m = e.to_string();
            assert!(m.contains("the DWG 2007 file is damaged beyond repair at") && m.contains(what), "{what}: {m}");
            m
        }
        Ok(_) => panic!("{what}: a damaged file read as a drawing"),
    }
}

/// 3 or 4 byte errors in one data page code word: more than RS(255, 251) corrects. About half the time the decoder
/// takes such a word to another code word; the page CRC and checksum catch that, so every case fails by name.
#[test]
fn r2007_data_word_miscorrection_is_caught() {
    let base = sample_2007();
    let (mut refused, mut caught) = (0, 0);
    for (at, what) in [(R2007_OBJECTS_PAGE, "'AcDb:AcDbObjects' data page 8"), (R2007_HEADER_PAGE, "'AcDb:Header' data page 16")] {
        for word in 0..at.1 {
            for errors in [3, 4] {
                for pattern in 0..6usize {
                    let mut f = base.clone();
                    for m in 0..errors {
                        let i = (word * 37 + pattern * 53 + m * 71 + 11) % 255;
                        f[at.0 + i * at.1 + word] ^= (1 + (pattern * 29 + m * 83 + word) % 255) as u8;
                    }
                    let message = beyond_repair(&f, what);
                    if message.contains("error correction can repair") { refused += 1 } else { caught += 1 }
                }
            }
        }
    }
    println!("{refused} refused by the decoder, {caught} miscorrections caught by the page CRC");
    assert!(caught > 0 && refused > 0, "{refused} {caught}");
}

/// GF(2^8) product modulo `poly` (its x^8 term left out): 0x1D for the data page field x^8 + x^4 + x^3 + x^2 + 1,
/// 0x69 for the system page field x^8 + x^6 + x^5 + x^3 + 1.
fn gf_mul(mut a: u8, mut b: u8, poly: u8) -> u8 {
    let mut p = 0u8;
    while b != 0 {
        if b & 1 != 0 {
            p ^= a;
        }
        let high = a & 0x80 != 0;
        a <<= 1;
        if high {
            a ^= poly;
        }
        b >>= 1;
    }
    p
}

/// The code's generator polynomial (x - a)(x - a^2)..(x - a^parity), a = 2, lowest coefficient first: itself a
/// code word.
fn generator(parity: usize, poly: u8) -> Vec<u8> {
    let mut g = vec![1u8];
    let mut root = 1u8;
    for _ in 0..parity {
        root = gf_mul(root, 2, poly);
        let mut next = vec![0u8; g.len() + 1];
        for (i, &c) in g.iter().enumerate() {
            next[i + 1] ^= c;
            next[i] ^= gf_mul(c, root, poly);
        }
        g = next;
    }
    g
}

/// The same for system pages and the file header: a code word added to every word (every copy) passes Reed-Solomon
/// and fails their CRCs by name; where a copy is left whole, that copy is read.
#[test]
fn r2007_undetectable_system_changes_fail_their_crcs() {
    let g = generator(16, 0x69);
    assert_eq!(g.len(), 17);
    let base = sample_2007();
    let good = dwg::read(&base).unwrap();
    // The generator shifted to 0, 20, .. 220 and summed (still a code word): every copy in the word is touched.
    let add = |f: &mut [u8], (at, blocks): (usize, usize), word: usize| {
        for shift in (0..=220).step_by(20) {
            for (k, &c) in g.iter().enumerate() {
                f[at + (shift + k) * blocks + word] ^= c;
            }
        }
    };
    for (at, what) in [(R2007_PAGE_MAP, "page map"), (R2007_SECTION_MAP, "section map")] {
        let mut f = base.clone();
        for word in 0..at.1 {
            add(&mut f, at, word);
        }
        let message = beyond_repair(&f, what);
        assert!(message.contains("CRC does not match"), "{message}");
        // Only the first word: the copies held in later words are read.
        let mut f = base.clone();
        add(&mut f, at, 0);
        assert_eq!(dwg::read(&f).unwrap(), good, "{what}: the first word changed");
    }
    // The file header: both copies changed fails; the copy at 0x80 changed alone reads from the end of the file.
    let mut f = base.clone();
    let end = (f.len() - 0x400, 3);
    add(&mut f, R2007_HEADER, 0);
    add(&mut f, end, 0);
    // (Here the length field changes too, so the refusal names an impossible length rather than the CRC.)
    beyond_repair(&f, "file header");
    // The header's 0xC8-byte block is repeated 3 times in its 3 words (data bytes 0, 200, 400): a generator at data
    // bytes 40, 240 and 440 (word 0 position 40, word 1 positions 1 and 201) of both header copies changes every
    // copy's stored data and no length field.
    let mut f = base.clone();
    for (at, blocks) in [R2007_HEADER, end] {
        for (word, position) in [(0, 40), (1, 1), (1, 201)] {
            for (k, &c) in g.iter().enumerate() {
                f[at + (position + k) * blocks + word] ^= c;
            }
        }
    }
    let message = beyond_repair(&f, "file header");
    assert!(message.contains("CRC does not match"), "{message}");
    let mut f = base.clone();
    add(&mut f, R2007_HEADER, 0);
    assert_eq!(dwg::read(&f).unwrap(), good, "file header copy at 0x80 changed");
}

/// Adding a code word of the data page code (the generator (x - a)(x - a^2)(x - a^3)(x - a^4), a = 2, shifted)
/// leaves every syndrome zero, so Reed-Solomon sees nothing; only the page CRC can tell the data is wrong. A change
/// that only touches parity or the filler after the page's stored bytes changes nothing read, and reads as before.
#[test]
fn r2007_undetectable_word_changes_fail_the_page_crc() {
    let g = generator(4, 0x1D);
    assert_eq!(g.len(), 5);
    let base = sample_2007();
    let good = dwg::read(&base).unwrap();
    let (mut caught, mut unused) = (0, 0);
    // (page, its stored (compressed) byte count, name)
    for (at, stored, what) in [(R2007_OBJECTS_PAGE, 6942, "'AcDb:AcDbObjects' data page 8"),
        (R2007_HEADER_PAGE, 1404, "'AcDb:Header' data page 16")] {
        for word in 0..at.1 {
            for shift in [0usize, 100, 246, 250] {
                let mut f = base.clone();
                for (k, &c) in g.iter().enumerate() {
                    f[at.0 + ((shift + k) % 255) * at.1 + word] ^= c;
                }
                let used = (0..5).map(|k| (shift + k) % 255).any(|p| p < 251 && word * 251 + p < stored);
                if used {
                    let message = beyond_repair(&f, what);
                    assert!(message.contains("CRC does not match"), "{message}");
                    caught += 1;
                } else {
                    assert_eq!(dwg::read(&f).unwrap(), good, "{what} word {word} shift {shift}");
                    unused += 1;
                }
            }
        }
    }
    println!("{caught} wrong code words caught by the page CRC, {unused} in unused bytes read as before");
    assert!(caught > 100 && unused > 0);
}

/// GCKERNEL_DWG_SCAN=<folder>[;<folder>...] reads every `*.dwg` under those folders (read only) and reports each
/// one's version, curve count, texts (the first few) and warnings, or the reason it was refused; it fails only on a panic.
#[test]
fn scan_dwg_folders() {
    let Ok(dirs) = std::env::var("GCKERNEL_DWG_SCAN") else { return };
    let mut files = Vec::new();
    let mut pending: Vec<PathBuf> = dirs.split(';').filter(|d| !d.trim().is_empty()).map(PathBuf::from).collect();
    while let Some(dir) = pending.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else { continue };
        for path in entries.filter_map(|e| e.ok().map(|e| e.path())) {
            if path.is_dir() {
                pending.push(path);
            } else if path.extension().is_some_and(|e| e.eq_ignore_ascii_case("dwg")) {
                files.push(path);
            }
        }
    }
    files.sort();
    let (mut read_ok, mut refused) = (0, 0);
    for path in &files {
        let bytes = std::fs::read(path).unwrap();
        let owned = bytes.clone();
        let outcome = std::panic::catch_unwind(move || dwg::read(&owned));
        let version = String::from_utf8_lossy(&bytes[..bytes.len().min(6)]).to_string();
        match outcome {
            Ok(Ok(d)) => {
                read_ok += 1;
                let mut skipped = d.skipped.clone();
                skipped.sort();
                println!("{}: {version} {} curves, {} layers, {} texts, skipped {skipped:?}, warnings {:?}", path.display(),
                    d.curves.len(), d.layers.len(), d.texts.len(), d.warnings);
                for t in d.texts.iter().take(8) {
                    println!("    {:?} {:?} at ({:.4}, {:.4}) h {:.4} rot {:.2} deg align {}/{} style {} font {}", t.kind, t.text,
                        t.at[0], t.at[1], t.height, t.rotation.to_degrees(), t.halign, t.valign, t.style, t.font);
                }
            }
            Ok(Err(e)) => {
                refused += 1;
                println!("{}: {version} REFUSED {e}", path.display());
            }
            Err(_) => panic!("{}: panicked", path.display()),
        }
    }
    println!("{} DWGs: {read_ok} read, {refused} refused", files.len());
}

/// Each sample's one TEXT ("not a curve", at (0, -20), 5 high, style "Default" in Tahoma) is read as text in every
/// file organisation.
#[test]
fn the_sample_text_reads_in_every_version() {
    for (name, bytes) in samples() {
        let d = dwg::read(&bytes).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(d.texts.len(), 1, "{name}: {:?}", d.texts);
        let t = &d.texts[0];
        assert_eq!((t.kind, t.text.as_str()), (TextKind::Text, "not a curve"), "{name}");
        assert!(t.style.eq_ignore_ascii_case("Default"), "{name}: style '{}' (names ignore case; R14 stores DEFAULT)", t.style);
        assert!(near_all(&[t.at[0], t.at[1], t.height, t.rotation, t.width], &[0.0, -20.0, 5.0, 0.0, 1.0], 20.0),
            "{name}: {t:?}");
        assert_eq!((t.halign, t.valign, t.mirrored), (0, 0, false), "{name}");
        assert!(t.font.eq_ignore_ascii_case("Tahoma"), "{name}: font '{}', styles {:?}", t.font, d.styles);
    }
}

/// GCKERNEL_DWG_FIXTURES: every DWG's texts equal its DXF twin's (strings, places, heights, turns, alignment,
/// style and font), in order.
#[test]
fn every_dwg_text_equals_its_dxf_twin() {
    let Some(pairs) = fixtures() else { return };
    let mut failed = Vec::new();
    for (dwg_path, dxf_path) in &pairs {
        let from_dxf = read(&std::fs::read(dxf_path).unwrap()).unwrap();
        let Ok(from_dwg) = dwg::read(&std::fs::read(dwg_path).unwrap()) else { continue };
        let scale = size(&from_dxf).max(1.0);
        let same = from_dwg.texts.len() == from_dxf.texts.len() && from_dwg.texts.iter().zip(&from_dxf.texts).all(|(a, b)| {
            a.kind == b.kind && a.text == b.text && a.style.eq_ignore_ascii_case(&b.style) && a.font.eq_ignore_ascii_case(&b.font)
                && (a.halign, a.valign, a.mirrored) == (b.halign, b.valign, b.mirrored)
                && near_all(&[a.at[0], a.at[1], a.end[0], a.end[1], a.height, a.box_width],
                    &[b.at[0], b.at[1], b.end[0], b.end[1], b.height, b.box_width], scale)
                && angle_near(a.rotation, b.rotation) && near(a.width, b.width, 1.0)
        });
        println!("{}: {} texts from the DWG, {} from the DXF: {}", name(dwg_path), from_dwg.texts.len(), from_dxf.texts.len(),
            if same { "same" } else { "DIFFERENT" });
        if !same {
            println!("    DWG {:?}\n    DXF {:?}", from_dwg.texts, from_dxf.texts);
            failed.push(name(dwg_path));
        }
    }
    assert!(failed.is_empty(), "texts differ: {failed:?}");
}

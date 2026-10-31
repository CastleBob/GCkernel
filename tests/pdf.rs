//! PDF / Illustrator reader (src/dxf/pdf): geometry, transforms, pages, layers, filters, cross-reference forms,
//! damaged files, encryption, and the user's real drawings (GCKERNEL_PDF_CORPUS). The fixtures come from
//! tests/fixtures/pdf/make_fixtures.py (python -I, standard library only).

use std::f64::consts::{FRAC_PI_2, TAU};

use gc_kernel::dxf::{fit, pdf, Drawing, Geometry, P2};
use gc_kernel::protocol::{encode_request, handle, parse_request, OP_DXF_READ};

fn fixture(name: &str) -> Vec<u8> {
    std::fs::read(format!("{}/tests/fixtures/pdf/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap()
}

fn read(name: &str) -> Drawing {
    pdf::read(&fixture(name)).unwrap_or_else(|e| panic!("{name}: {e}"))
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

fn near_p(a: P2, b: P2) -> bool {
    near(a[0], b[0]) && near(a[1], b[1])
}

fn line(g: &Geometry) -> (P2, P2) {
    match g {
        Geometry::Line { a, b } => (*a, *b),
        g => panic!("not a line: {g:?}"),
    }
}

fn is_line(g: &Geometry, a: P2, b: P2) -> bool {
    matches!(g, Geometry::Line { a: x, b: y } if near_p(*x, a) && near_p(*y, b))
}

fn layer_of<'a>(d: &'a Drawing, i: usize) -> &'a str {
    &d.layers[d.curves[i].layer].name
}

fn kinds(d: &Drawing) -> String {
    let (mut l, mut a, mut c, mut e, mut s) = (0, 0, 0, 0, 0);
    for curve in &d.curves {
        match curve.geometry {
            Geometry::Line { .. } => l += 1,
            Geometry::Arc { .. } => a += 1,
            Geometry::Circle { .. } => c += 1,
            Geometry::Ellipse { .. } => e += 1,
            _ => s += 1,
        }
    }
    format!("{l} lines, {a} arcs, {c} circles, {e} ellipses, {s} splines")
}

#[test]
fn lines_curves_rectangles_circles_arcs_ellipses_forms_and_skipped_content() {
    let d = read("shapes.pdf");
    assert_eq!(d.insunits, Some(1));
    assert_eq!(d.version, "PDF 1.4");
    let g: Vec<&Geometry> = d.curves.iter().map(|c| &c.geometry).collect();
    assert_eq!(g.len(), 11, "{g:#?}");
    // PDF units are 1/72 inch.
    assert!(is_line(g[0], [1.0, 1.0], [2.0, 1.0]), "{:?}", g[0]);
    match g[1] {
        Geometry::Spline { degree: 3, knots, points, weights } => {
            assert_eq!(knots, &[0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]);
            assert!(weights.is_empty());
            let want = [[1.0, 2.0], [1.5, 2.5], [2.0, 1.5], [2.5, 2.0]];
            assert!(points.iter().zip(want).all(|(p, w)| near_p(*p, w)), "{points:?}");
        }
        g => panic!("{g:?}"),
    }
    // re: four lines around (3, 3) .. (4, 3.5).
    assert!(is_line(g[2], [3.0, 3.0], [4.0, 3.0]) && is_line(g[3], [4.0, 3.0], [4.0, 3.5]));
    assert!(is_line(g[4], [4.0, 3.5], [3.0, 3.5]) && is_line(g[5], [3.0, 3.5], [3.0, 3.0]));
    // Four 0.5523-handle quarters: an exact circle.
    match g[6] {
        Geometry::Circle { c, r } => assert!(near_p(*c, [5.0, 5.0]) && near(*r, 0.5), "{c:?} {r}"),
        g => panic!("{g:?}"),
    }
    // One quarter with the exact handle: an exact arc.
    match g[7] {
        Geometry::Arc { c, r, start, sweep } => {
            assert!(near_p(*c, [1.0, 8.0]) && near(*r, 1.0) && near(*start, 0.0) && near(*sweep, FRAC_PI_2), "{:?}", g[7]);
        }
        g => panic!("{g:?}"),
    }
    match g[8] {
        Geometry::Ellipse { c, u, v, t0, t1 } => {
            assert!(near_p(*c, [7.5, 2.5]) && near_p(*u, [1.0, 0.0]) && near_p(*v, [0.0, 0.5]), "{:?}", g[8]);
            assert!(near(*t0, 0.0) && near(*t1, TAU));
        }
        g => panic!("{g:?}"),
    }
    // Nested q/cm: translate 36 then scale 2.
    assert!(is_line(g[9], [1.0, 9.0], [2.0, 9.0]), "{:?}", g[9]);
    // The form's /Matrix lifts its line by 72.
    assert!(is_line(g[10], [4.0, 2.0], [5.0, 2.0]), "{:?}", g[10]);
    // "W n" (a clip with no painting) draws nothing; inline images (one holding " EI " in its data), the image
    // XObject, the text block and the shading are counted.
    let mut skipped = d.skipped.clone();
    skipped.sort();
    assert_eq!(skipped, vec![("IMAGE".to_string(), 3), ("SHADING".to_string(), 1), ("TEXT".to_string(), 1)]);
    assert!(d.warnings.is_empty(), "{:?}", d.warnings);
}

#[test]
fn crop_box_rotation_and_user_unit() {
    // Crop box [36 36 576 756] on a 612 x 792 page, /Rotate 90, /UserUnit 2: (72, 72) -> (36, 36) in the box,
    // turned a quarter clockwise in a 540 wide box -> (36, 504), at 2/72 inch per unit -> (1, 14).
    let d = read("rotated.pdf");
    assert_eq!(d.curves.len(), 1);
    assert!(is_line(&d.curves[0].geometry, [1.0, 14.0], [1.0, 12.0]), "{:?}", d.curves[0].geometry);
    assert_eq!(d.version, "PDF 1.6");
}

#[test]
fn first_page_in_tree_order_with_inherited_resources_and_split_contents() {
    // Page 1 sits in a nested Pages node that supplies its resources; the media box comes from the root; the
    // operand and operator of "/F0 Do" are in two content streams.
    let d = read("pages.pdf");
    assert_eq!(d.curves.len(), 1);
    assert!(is_line(&d.curves[0].geometry, [0.0, 0.0], [1.0, 1.0]));
    assert_eq!(d.warnings, vec!["page 1 of 3 was read".to_string()]);
}

#[test]
fn optional_content_groups_are_layers_and_off_groups_are_left_out() {
    let d = read("layers.pdf");
    let names: Vec<&str> = d.layers.iter().map(|l| l.name.as_str()).collect();
    assert_eq!(names, ["0", "Cut", "\u{c4}tch", "Hidden"]);
    assert!(d.layers[3].off && !d.layers[1].off);
    // Plain, Cut, Etch and the non-OC marked content; the OFF group, the all-on membership dictionary that includes
    // it, and the form whose /OC is the OFF group are left out.
    assert_eq!(d.curves.len(), 4, "{:#?}", d.curves);
    assert_eq!((layer_of(&d, 0), layer_of(&d, 1), layer_of(&d, 2), layer_of(&d, 3)), ("0", "Cut", "\u{c4}tch", "0"));
    assert!(is_line(&d.curves[3].geometry, [0.0, 4.0], [1.0, 4.0]));
    assert!(d.warnings.iter().any(|w| w == "layer 'Hidden' is off in the PDF; its content was left out"), "{:?}", d.warnings);
    // Header-only reading gives the layers without curves.
    let h = pdf::read_header(&fixture("layers.pdf")).unwrap();
    assert!(h.curves.is_empty() && h.layers.len() == 4);
}

#[test]
fn cross_reference_stream_object_stream_and_png_predictors() {
    // PDF 1.5: objects packed in an object stream, a cross-reference stream with the PNG Up predictor, and a
    // content stream whose rows use every PNG filter type (Predictor 15).
    let d = read("xrefstream.pdf");
    assert_eq!(d.curves.len(), 40);
    for (i, c) in d.curves.iter().enumerate() {
        let x = 10.0 * i as f64 / 72.0;
        assert!(is_line(&c.geometry, [x, 0.0], [x, 1.0]), "{i}: {:?}", c.geometry);
    }
    assert!(d.warnings.is_empty(), "{:?}", d.warnings);
}

#[test]
fn hybrid_reference_file_and_incremental_update() {
    let d = read("hybrid.pdf");
    assert_eq!(d.curves.len(), 1);
    assert!(is_line(&d.curves[0].geometry, [0.0, 0.0], [1.0, 0.0]));
    assert!(d.warnings.is_empty(), "{:?}", d.warnings);
    // The update's newer contents win over the original's.
    let d = read("updated.pdf");
    assert_eq!(d.curves.len(), 1);
    assert!(is_line(&d.curves[0].geometry, [0.0, 0.0], [2.0, 0.0]), "{:?}", d.curves[0].geometry);
}

#[test]
fn broken_startxref_and_wrong_stream_lengths_are_recovered() {
    let d = read("damaged.pdf");
    assert_eq!(d.curves.len(), 3, "{:?}", d.curves);
    for (i, c) in d.curves.iter().enumerate() {
        let y = i as f64;
        assert!(is_line(&c.geometry, [0.0, y], [1.0, y]), "{i}: {:?}", c.geometry);
    }
    assert!(d.warnings.iter().any(|w| w.contains("cross-reference table was missing or damaged")), "{:?}", d.warnings);
    assert!(d.warnings.iter().any(|w| w.contains("/Length") && w.contains("endstream")), "{:?}", d.warnings);
}

#[test]
fn lzw_ascii85_ascii_hex_and_run_length_filters() {
    let d = read("filters.pdf");
    // 300 lines through LZW (EarlyChange 1, with a clear code part way), 40 through LZW with EarlyChange 0, one
    // through [ASCII85 Flate] and one through [ASCIIHex RunLength]. Zero-length lines are dropped.
    let text = std::fs::read_to_string(format!("{}/tests/fixtures/pdf/filters_lines.txt", env!("CARGO_MANIFEST_DIR"))).unwrap();
    let mut expected = 0;
    for (i, l) in text.lines().enumerate() {
        let v: Vec<f64> = l.split_whitespace().filter_map(|t| t.parse().ok()).collect();
        let (a, b) = ([v[0] / 72.0, v[1] / 72.0], [v[2] / 72.0, v[3] / 72.0]);
        if fit::dist(a, b) > 0.0 {
            expected += 1;
            assert!(is_line(&d.curves[expected - 1].geometry, a, b), "line {i}: {l} vs {:?}", d.curves[expected - 1].geometry);
        }
    }
    assert_eq!(d.curves.len(), expected + 40 + 2, "{}", kinds(&d));
    let last = |n: usize| line(&d.curves[d.curves.len() - n].geometry);
    assert!(near(last(1).0[1], 10.0) && near(last(2).0[1], 710.0 / 72.0) && near(last(3).0[1], 700.0 / 72.0));
    assert!(d.warnings.is_empty(), "{:?}", d.warnings);
}

fn check_encrypted(name: &str, identity_stream: bool, method: &str) {
    let d = read(name);
    let names: Vec<&str> = d.layers.iter().map(|l| l.name.as_str()).collect();
    // The layer names are encrypted strings.
    assert_eq!(names, ["0", "Cut", "\u{c4}tch", "Hidden"], "{name}");
    assert_eq!(d.curves.len(), 4 + usize::from(identity_stream), "{name}: {:?}", d.curves);
    assert_eq!(layer_of(&d, 1), "Cut");
    if identity_stream {
        assert!(is_line(&d.curves[4].geometry, [0.0, 500.0 / 72.0], [1.0, 500.0 / 72.0]));
    }
    assert!(d.warnings.iter().any(|w| w.contains(method) && w.contains("opens without a password")), "{name}: {:?}", d.warnings);
}

#[test]
fn rc4_40_revision_2_with_no_user_password() {
    check_encrypted("rc4_40_r2.pdf", false, "RC4 40-bit");
}

#[test]
fn rc4_128_revision_3_with_no_user_password() {
    check_encrypted("rc4_r3.pdf", false, "RC4 128-bit");
}

#[test]
fn aes_128_revision_4_with_identity_filter_and_clear_metadata() {
    check_encrypted("aes128_r4.pdf", true, "AES 128-bit");
}

#[test]
fn aes_256_revision_6_with_object_streams() {
    check_encrypted("aes256_r6.pdf", true, "AES 256-bit");
}

#[test]
fn a_user_password_is_refused_by_name() {
    let e = pdf::read(&fixture("password.pdf")).unwrap_err();
    assert_eq!(e.to_string(), "unsupported: this PDF needs a password to open");
}

#[test]
fn illustrator_files_with_and_without_a_pdf_part() {
    let ai = fixture("illustrator.ai");
    assert!(pdf::looks_like(&ai) && pdf::looks_like(&fixture("shapes.pdf")));
    let d = pdf::read(&ai).unwrap();
    // The private Illustrator stream (holding "endstream" in its data) is never drawn.
    assert_eq!(d.curves.len(), 3);
    assert!(is_line(&d.curves[2].geometry, [1.0, 1.0], [0.0, 0.0]));
    let ps = b"%!PS-Adobe-3.0 EPSF-3.0\n%%Creator: Adobe Illustrator(R) 8.0\n%%BoundingBox: 0 0 100 100\n";
    assert!(!pdf::looks_like(ps) && pdf::looks_like_postscript(ps));
    let e = pdf::read(ps).unwrap_err().to_string();
    assert!(e.contains("this Illustrator file has no PDF part") && e.contains("re-save it with PDF compatibility"), "{e}");
    assert!(pdf::read(b"hello").unwrap_err().to_string().contains("not a PDF"));
}

/// Op 14 tells a PDF (or Illustrator file) by its content and refuses PostScript-only Illustrator files by name.
#[test]
fn the_read_op_takes_pdf_files() {
    let request = |path: &str, mode: f64| {
        let mut req = vec![1.0, mode, 1e-3, path.len() as f64];
        req.extend(path.bytes().map(f64::from));
        handle(&parse_request(&encode_request(OP_DXF_READ, 1e-6, &req)).unwrap())
    };
    let path = format!("{}/tests/fixtures/pdf/layers.pdf", env!("CARGO_MANIFEST_DIR"));
    for mode in [0.0, 1.0, 2.0] {
        let r = request(&path, mode).unwrap();
        assert_eq!(&r[..3], &[1.0, 1.0, -1.0], "version 1, inches, no measurement");
        let version: String = r[4..4 + r[3] as usize].iter().map(|b| *b as u8 as char).collect();
        assert_eq!(version, "PDF 1.5");
        let mut at = 4 + r[3] as usize;
        assert_eq!(r[at], 4.0, "four layers");
        at += 1;
        for _ in 0..4 {
            at += 2;
            at += 1 + r[at] as usize;
        }
        assert_eq!(r[at], if mode == 2.0 { 0.0 } else { 4.0 });
    }
    let dir = std::env::temp_dir().join(format!("gck_pdf_ps_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let ps = dir.join("old.ai");
    std::fs::write(&ps, b"%!PS-Adobe-2.0
%%Creator: Adobe Illustrator(TM) 5.0
").unwrap();
    let e = request(ps.to_str().unwrap(), 0.0).unwrap_err().to_string();
    let _ = std::fs::remove_dir_all(&dir);
    assert!(e.contains("this Illustrator file has no PDF part"), "{e}");
}

/// A page whose forms each draw the next twenty times, twelve deep: 20^11 uses.
fn hostile_forms() -> Vec<u8> {
    let mut body = String::from("%PDF-1.4\n");
    let mut offsets = Vec::new();
    let mut add = |body: &mut String, text: String| {
        offsets.push(body.len());
        body.push_str(&text);
    };
    add(&mut body, "1 0 obj <</Type/Catalog/Pages 2 0 R>> endobj\n".into());
    add(&mut body, "2 0 obj <</Type/Pages/Kids[3 0 R]/Count 1>> endobj\n".into());
    add(&mut body, "3 0 obj <</Type/Page/Parent 2 0 R/MediaBox[0 0 10 10]/Contents 4 0 R/Resources<</XObject<</X 5 0 R>>>>>> endobj\n".into());
    add(&mut body, "4 0 obj <</Length 5>>stream\n/X Do\nendstream endobj\n".into());
    let forms = 12;
    for i in 0..forms {
        let (draw, res) = if i + 1 < forms {
            ("/X Do ".repeat(20), format!("/Resources<</XObject<</X {} 0 R>>>>", 6 + i))
        } else {
            ("0 0 m 1 1 l S ".to_string(), String::new())
        };
        add(&mut body, format!("{} 0 obj <</Type/XObject/Subtype/Form/BBox[0 0 1 1]{res}/Length {}>>stream\n{draw}\nendstream endobj\n", 5 + i, draw.len()));
    }
    let xref = body.len();
    body.push_str(&format!("xref\n0 {}\n0000000000 65535 f \n", offsets.len() + 1));
    for o in &offsets {
        body.push_str(&format!("{o:010} 00000 n \n"));
    }
    body.push_str(&format!("trailer <</Size {}/Root 1 0 R>>\nstartxref\n{xref}\n%%EOF\n", offsets.len() + 1));
    body.into_bytes()
}

#[test]
fn damaged_files_never_panic() {
    let started = std::time::Instant::now();
    let mut seed = 0x9E3779B97F4A7C15u64;
    let mut next = || { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; seed };
    for name in ["shapes.pdf", "xrefstream.pdf", "layers.pdf", "filters.pdf", "rc4_r3.pdf", "aes128_r4.pdf", "damaged.pdf", "pages.pdf"] {
        let base = fixture(name);
        for cut in (0..base.len()).step_by(7) {
            try_read(&base[..cut]);
        }
        for _ in 0..1500 {
            let mut f = base.clone();
            for _ in 0..1 + next() % 4 {
                let at = (next() as usize) % f.len();
                match next() % 5 {
                    0 => f[at] = next() as u8,
                    1 => f[at] = b"0123456789-. /<>[]()R\nobjstreamxref"[(next() % 34) as usize],
                    2 => { f.remove(at); }
                    3 => { let b = f[at]; f.insert(at, b); }
                    _ => {
                        // A digit becomes a large number (lengths, counts, offsets, sizes).
                        if f[at].is_ascii_digit() {
                            f.splice(at..at + 1, b"99999999".iter().copied());
                        }
                    }
                }
            }
            try_read(&f);
        }
    }
    let fuzz = started.elapsed();
    let t = std::time::Instant::now();
    let e = pdf::read(&hostile_forms()).unwrap_err().to_string();
    assert!(e.contains("uses forms more than"), "{e}");
    println!("fuzz took {fuzz:?}; hostile forms refused in {:?}", t.elapsed());
}

fn try_read(bytes: &[u8]) {
    let owned = bytes.to_vec();
    let outcome = std::panic::catch_unwind(move || {
        if let Ok(d) = pdf::read(&owned) {
            for c in &d.curves {
                let _ = fit::lines_and_arcs(&c.geometry, 1e-3);
            }
        }
    });
    assert!(outcome.is_ok(), "panicked on {:?}", String::from_utf8_lossy(bytes));
}

#[test]
fn corpus_pdfs_read() {
    let Ok(list) = std::env::var("GCKERNEL_PDF_CORPUS") else { return };
    let text = std::fs::read_to_string(&list).unwrap();
    let mut failures = Vec::new();
    for line in text.lines().map(str::trim).filter(|l| !l.is_empty()) {
        // Paths may be written the Git Bash way (/c/...).
        let path = match line.strip_prefix("/c/") { Some(rest) => format!("C:/{rest}"), None => line.to_string() };
        let bytes = std::fs::read(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let start = std::time::Instant::now();
        match pdf::read(&bytes) {
            Ok(d) => {
                let mut fitted = 0;
                for c in &d.curves {
                    fitted += fit::lines_and_arcs(&c.geometry, 1e-3).unwrap().len();
                }
                println!("{path} | {} | {} | fitted {fitted} | layers {:?} | skipped {:?} | notes {:?} | {:?}", d.version, kinds(&d),
                    d.layers.iter().map(|l| l.name.as_str()).collect::<Vec<_>>(), d.skipped, d.warnings, start.elapsed());
            }
            Err(e) => {
                println!("{path} | FAILED {e}");
                failures.push(format!("{path}: {e}"));
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

use gc_kernel::{mass, neutral, protocol, step, tessellate, topology, KernelError};
use std::path::Path;
const TOL: f64 = 1e-6;

// Independently authored 2 x 3 x 4 box; source graph is intentionally unrelated
// to the importer implementation. Dimensions/volume are an analytic oracle.
fn box_data(mm: bool) -> String {
    let mut s = String::new();
    let scale = if mm { 25.4 } else { 1.0 };
    let vertices = [
        (0., 0., 0.),
        (2., 0., 0.),
        (2., 3., 0.),
        (0., 3., 0.),
        (0., 0., 4.),
        (2., 0., 4.),
        (2., 3., 4.),
        (0., 3., 4.),
    ];
    for (i, (x, y, z)) in vertices.into_iter().enumerate() {
        s += &format!(
            "#{}=CARTESIAN_POINT('',({},{},{}));#{}=VERTEX_POINT('',#{});\n",
            i + 1,
            x * scale,
            y * scale,
            z * scale,
            i + 11,
            i + 1
        );
    }
    let pairs = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];
    for (i, (a, b)) in pairs.into_iter().enumerate() {
        let n = 100 + i * 4;
        let d = (
            vertices[b].0 - vertices[a].0,
            vertices[b].1 - vertices[a].1,
            vertices[b].2 - vertices[a].2,
        );
        s += &format!("#{n}=DIRECTION('',({},{},{}));#{}=VECTOR('',#{n},1.);#{}=LINE('',#{},#{});#{}=EDGE_CURVE('',#{},#{},#{},.T.);\n",d.0,d.1,d.2,n+1,n+2,a+1,n+1,n+3,a+11,b+11,n+2);
    }
    let faces = [
        (
            [(3, false), (2, false), (1, false), (0, false)],
            (0, 0, -1),
            (1, 0, 0),
            1,
        ),
        (
            [(4, true), (5, true), (6, true), (7, true)],
            (0, 0, 1),
            (1, 0, 0),
            5,
        ),
        (
            [(0, true), (9, true), (4, false), (8, false)],
            (0, -1, 0),
            (1, 0, 0),
            1,
        ),
        (
            [(1, true), (10, true), (5, false), (9, false)],
            (1, 0, 0),
            (0, 1, 0),
            2,
        ),
        (
            [(2, true), (11, true), (6, false), (10, false)],
            (0, 1, 0),
            (1, 0, 0),
            3,
        ),
        (
            [(3, true), (8, true), (7, false), (11, false)],
            (-1, 0, 0),
            (0, 1, 0),
            1,
        ),
    ];
    for (i, (edges, n, x, p)) in faces.into_iter().enumerate() {
        let f = 200 + i * 20;
        for (j, (e, d)) in edges.into_iter().enumerate() {
            s += &format!(
                "#{}=ORIENTED_EDGE('',*,*,#{},.{}.);",
                f + j,
                103 + 4 * e,
                if d { "T" } else { "F" }
            );
        }
        s += &format!(
            "#{}=EDGE_LOOP('',(#{},#{},#{},#{}));#{}=FACE_OUTER_BOUND('',#{},.T.);",
            f + 4,
            f,
            f + 1,
            f + 2,
            f + 3,
            f + 5,
            f + 4
        );
        s += &format!(
            "#{}=DIRECTION('',({},{},{}));#{}=DIRECTION('',({},{},{}));",
            f + 6,
            n.0,
            n.1,
            n.2,
            f + 7,
            x.0,
            x.1,
            x.2
        );
        s += &format!("#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});#{}=PLANE('',#{});#{}=ADVANCED_FACE('',(#{}),#{},.T.);\n",f+8,p,f+6,f+7,f+9,f+8,f+10,f+5,f+9);
    }
    s += "#400=CLOSED_SHELL('',(#210,#230,#250,#270,#290,#310));#401=MANIFOLD_SOLID_BREP('Box',#400);";
    s + &units(mm, "#401")
}
fn units(mm: bool, items: &str) -> String {
    let length = if mm {
        "#500=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));".to_owned()
    } else {
        "#500=(CONVERSION_BASED_UNIT('inch',#503)LENGTH_UNIT()NAMED_UNIT(#504));#503=MEASURE_WITH_UNIT(LENGTH_MEASURE(0.0254),#505);#504=DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.);#505=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.));".into()
    };
    length + "#501=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));#502=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#500,#501))REPRESENTATION_CONTEXT('',''));"
        + &format!("#506=ADVANCED_BREP_SHAPE_REPRESENTATION('',({items}),#502);")
}
fn file(data: &str) -> String {
    format!("ISO-10303-21;HEADER;FILE_DESCRIPTION(('Test'),'2;1');FILE_NAME('test','',(''),(''),'GC test','','');FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));ENDSEC;DATA;{data}ENDSEC;END-ISO-10303-21;")
}
fn read(s: &str) -> step::StepImport {
    step::read(s.as_bytes(), TOL).unwrap()
}
fn refused(s: &str, fragment: &str) {
    let error = step::read(s.as_bytes(), TOL).unwrap_err().to_string();
    assert!(
        error.contains(fragment),
        "expected {fragment:?}, got {error}"
    );
}

#[test]
fn inches_and_mm_have_identical_dimensions_mass_topology_and_mesh() {
    let a = read(&file(&box_data(false)));
    let b = read(&file(&box_data(true)));
    for import in [&a, &b] {
        let b = &import.bodies[0];
        assert!((b.volume - 24.0).abs() < 1e-10);
        let m = mass::measure(&b.body, TOL).unwrap();
        assert!((m.area - 52.0).abs() < 1e-10);
        let report = topology::check(&b.body, TOL).unwrap();
        assert_eq!((report.vertices, report.edges, report.faces), (8, 12, 6));
        let mesh = tessellate::tessellate(&b.body, 0.01).unwrap();
        assert_eq!(mesh.triangles.len(), 12);
        let packet = neutral::write_body(&b.body);
        let mut r = neutral::Reader::new(&packet);
        let round = r.body().unwrap();
        r.finish().unwrap();
        assert_eq!(neutral::write_body(&round), packet);
    }
    for (a, b) in a.bodies[0]
        .body
        .vertices
        .iter()
        .zip(&b.bodies[0].body.vertices)
    {
        assert!(a.distance(*b) < 1e-12);
    }
}
#[test]
fn syntax_preserves_escaped_unicode_names_and_ignores_physical_line_breaks() {
    let data = box_data(false).replace(
        "'Box'",
        "'Bob''s \\X2\\03A9\\X0\\ \\X4\\0001F680\\X0\\ \\\\ model'",
    );
    let data = data.replace("MANIFOLD_SOLID_BREP", "MANIFOLD_SOL\r\nID_BREP");
    let b = read(&file(&data));
    assert_eq!(b.bodies[0].name, "Bob's Ω 🚀 \\ model");
    let b = read(&file(&box_data(false).replace("'Box'", "'Ω'")));
    assert_eq!(b.bodies[0].name, "Ω");
}
#[test]
fn ordering_and_comments_do_not_change_bodies() {
    let data = box_data(false);
    let original = read(&file(&data));
    let mut records: Vec<_> = data.split(';').filter(|s| !s.trim().is_empty()).collect();
    records.reverse();
    let shuffled = read(&file(&(records.join(";/* #999= ignored; */") + ";")));
    assert_eq!(
        neutral::write_body(&original.bodies[0].body),
        neutral::write_body(&shuffled.bodies[0].body)
    );
}
#[test]
fn malformed_and_ambiguous_input_is_refused() {
    let good = file(&box_data(false));
    let cases = [
        (
            good.replace("#1=CARTESIAN_POINT", "#2=CARTESIAN_POINT"),
            "duplicate STEP entity",
        ),
        (
            good.replace("#11=VERTEX_POINT('',#1)", "#11=VERTEX_POINT('',#99999)"),
            "missing #99999",
        ),
        (
            good.replace("LENGTH_MEASURE(0.0254)", "LENGTH_MEASURE(-1.)"),
            "nonpositive",
        ),
        (
            good.replace(
                "#503=MEASURE_WITH_UNIT(LENGTH_MEASURE(0.0254),#505)",
                "#503=MEASURE_WITH_UNIT(LENGTH_MEASURE(0.0254),#500)",
            ),
            "cyclic",
        ),
        (
            good.replace("SI_UNIT($,.RADIAN.)", "SI_UNIT($,.METRE.)"),
            "not RADIAN",
        ),
        (
            good.replace(
                "#103=EDGE_CURVE('',#11,#12,#102,.T.)",
                "#103=EDGE_CURVE('',#11,#12,#102,.F.)",
            ),
            "same_sense",
        ),
        (
            good.replace(
                "#205=FACE_OUTER_BOUND('',#204,.T.)",
                "#205=FACE_OUTER_BOUND('',#204,.F.)",
            ),
            "closed solid",
        ),
        (
            good.replace(
                "#210,#230,#250,#270,#290,#310",
                "#210,#230,#250,#270,#290,#310,#210",
            ),
            "duplicate face",
        ),
        (
            good.replace("LENGTH_MEASURE(0.0254)", "LENGTH_MEASURE(1.E999)"),
            "non-finite",
        ),
        (good.replace("END-ISO-10303-21;", ""), "expected a keyword"),
        (good.clone() + "garbage", "trailing"),
        (good.replace("'Box'", "'\\X2\\D800\\X0\\'"), "surrogate"),
        (
            good.replace("'Box'", "'unterminated"),
            "expected comma or closing parenthesis",
        ),
    ];
    for (text, reason) in cases {
        refused(&text, reason);
    }
    refused(&(good + "/*"), "unterminated comment");
}
#[test]
fn truncation_at_every_byte_returns_errors_without_panicking() {
    let good = file(&box_data(false));
    for i in 0..good.len() {
        assert!(
            step::read(&good.as_bytes()[..i], TOL).is_err(),
            "prefix {i}"
        );
    }
}
#[test]
fn all_bodies_are_kept_and_one_bad_body_refuses_the_whole_file() {
    let data = box_data(false).replace("#401),#502", "#401,#402),#502")
        + "#402=MANIFOLD_SOLID_BREP('Second',#400);";
    let imported = read(&file(&data));
    assert_eq!(imported.bodies.len(), 2);
    assert_eq!(imported.bodies[1].name, "Second");
    assert_eq!(imported.bodies[0].volume, imported.bodies[1].volume);
    refused(
        &file(&(data.replace("'Second',#400", "'Second',#403") + "#403=OPEN_SHELL('',(#210));")),
        "CLOSED_SHELL",
    );
    refused(
        &file(&(box_data(false) + "#600=MAPPED_ITEM('',#401,#401);")),
        "MAPPED_ITEM",
    );
    refused(
        &file(
            &(box_data(false)
                + "#600=GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION('',(),#502);"),
        ),
        "WIREFRAME",
    );
}
#[test]
fn missing_units_unrepresented_body_and_unknown_surface_never_get_guessed() {
    let data = box_data(false);
    refused(
        &file(&data.replace(
            "GLOBAL_UNIT_ASSIGNED_CONTEXT((#500,#501))",
            "GLOBAL_UNIT_ASSIGNED_CONTEXT((#501))",
        )),
        "lacks declared",
    );
    refused(
        &file(&data.replace("(#401),#502", "(),#502")),
        "no unambiguous",
    );
    refused(
        &file(&data.replace(
            "#209=PLANE('',#208)",
            "#209=B_SPLINE_SURFACE_WITH_KNOTS('',#208)",
        )),
        "B_SPLINE_SURFACE",
    );
    for t in [0., -1., f64::INFINITY, f64::NAN] {
        assert!(step::read(file(&data).as_bytes(), t).is_err());
    }
}
#[test]
fn bounded_parser_refuses_deep_nesting_and_unit_chain() {
    let nested = "(".repeat(70) + "1" + &")".repeat(70);
    refused(
        &file(&(box_data(false) + &format!("#600=OTHER({nested});"))),
        "nesting",
    );
    let huge_name = "x".repeat(32770);
    refused(
        &file(&box_data(false).replace("'Box'", &format!("'{huge_name}'"))),
        "string length",
    );
}

fn cylinder_data() -> String {
    let data="
#1=CARTESIAN_POINT('',(2.,0.,0.));#2=CARTESIAN_POINT('',(2.,0.,4.));
#3=VERTEX_POINT('',#1);#4=VERTEX_POINT('',#2);
#5=CARTESIAN_POINT('',(0.,0.,0.));#6=CARTESIAN_POINT('',(0.,0.,4.));
#7=DIRECTION('',(0.,0.,1.));#8=DIRECTION('',(1.,0.,0.));
#9=AXIS2_PLACEMENT_3D('',#5,#7,#8);#10=AXIS2_PLACEMENT_3D('',#6,#7,#8);
#11=CIRCLE('',#9,2.);#12=CIRCLE('',#10,2.);#13=VECTOR('',#7,1.);#14=LINE('',#1,#13);
#15=EDGE_CURVE('',#3,#3,#11,.T.);#16=EDGE_CURVE('',#4,#4,#12,.T.);#17=EDGE_CURVE('',#3,#4,#14,.T.);
#18=ORIENTED_EDGE('',*,*,#15,.F.);#19=EDGE_LOOP('',(#18));#20=FACE_OUTER_BOUND('',#19,.T.);
#21=PLANE('',#9);#22=ADVANCED_FACE('',(#20),#21,.F.);
#23=ORIENTED_EDGE('',*,*,#16,.T.);#24=EDGE_LOOP('',(#23));#25=FACE_OUTER_BOUND('',#24,.T.);
#26=PLANE('',#10);#27=ADVANCED_FACE('',(#25),#26,.T.);
#28=ORIENTED_EDGE('',*,*,#15,.T.);#29=ORIENTED_EDGE('',*,*,#17,.T.);#30=ORIENTED_EDGE('',*,*,#16,.F.);#31=ORIENTED_EDGE('',*,*,#17,.F.);
#32=EDGE_LOOP('',(#28,#29,#30,#31));#33=FACE_OUTER_BOUND('',#32,.T.);#34=CYLINDRICAL_SURFACE('',#9,2.);#35=ADVANCED_FACE('',(#33),#34,.T.);
#400=CLOSED_SHELL('',(#22,#27,#35));#401=MANIFOLD_SOLID_BREP('Cylinder',#400);
";
    data.to_owned() + &units(false, "#401")
}
#[test]
fn full_circle_seam_edges_remain_exact_and_allow_kernel_operations() {
    let imported = read(&file(&cylinder_data()));
    let body = &imported.bodies[0].body;
    assert!((imported.bodies[0].volume - 16. * std::f64::consts::PI).abs() < 1e-10);
    assert!(matches!(
        body.edges[0].curve,
        gc_kernel::geom::Curve::Arc { .. }
    ));
    assert!(tessellate::tessellate(body, 0.02).unwrap().triangles.len() > 20);
    let prism = gc_kernel::rebuild::rebuild_prism(body, TOL).unwrap();
    assert!((mass::measure(&prism, TOL).unwrap().volume - imported.bodies[0].volume).abs() < 1e-10);
    let reversed = cylinder_data()
        .replace(
            "#15=EDGE_CURVE('',#3,#3,#11,.T.)",
            "#15=EDGE_CURVE('',#3,#3,#11,.F.)",
        )
        .replace("#15,.F.", "#15,.SWAP.")
        .replace("#15,.T.", "#15,.F.")
        .replace("#15,.SWAP.", "#15,.T.");
    assert!((read(&file(&reversed)).bodies[0].volume - imported.bodies[0].volume).abs() < 1e-10);
}
#[test]
fn ellipse_extrusion_is_exact_not_a_polygon() {
    let d = cylinder_data()
        .replace("CIRCLE('',#9,2.)", "ELLIPSE('',#9,2.,1.)")
        .replace("CIRCLE('',#10,2.)", "ELLIPSE('',#10,2.,1.)")
        .replace(
            "CYLINDRICAL_SURFACE('',#9,2.)",
            "SURFACE_OF_LINEAR_EXTRUSION('',#11,#13)",
        );
    let imported = read(&file(&d));
    assert!((imported.bodies[0].volume - 8. * std::f64::consts::PI).abs() < 1e-9);
    assert!(matches!(
        imported.bodies[0].body.edges[0].curve,
        gc_kernel::geom::Curve::Conic { .. }
    ));
}

#[test]
fn a_small_reversed_curved_wall_cannot_hide_in_positive_total_area() {
    let bad = cylinder_data()
        .replace("(2.,0.,4.)", "(2.,0.,0.1)")
        .replace("(0.,0.,4.)", "(0.,0.,0.1)")
        .replace(
            "#35=ADVANCED_FACE('',(#33),#34,.T.)",
            "#35=ADVANCED_FACE('',(#33),#34,.F.)",
        );
    refused(&file(&bad), "curved face");
}

#[test]
fn cone_and_line_revolution_import_the_same_analytic_frustum() {
    let base = cylinder_data()
        .replace("(2.,0.,4.)", "(3.,0.,4.)")
        .replace("CIRCLE('',#10,2.)", "CIRCLE('',#10,3.)")
        .replace("#13=VECTOR('',#7,1.)", "#13=VECTOR('',#600,1.)")
        + "#600=DIRECTION('',(1.,0.,4.));#601=AXIS1_PLACEMENT('',#5,#7);";
    let cone = base.replace(
        "CYLINDRICAL_SURFACE('',#9,2.)",
        &format!("CONICAL_SURFACE('',#9,2.,{})", 0.25_f64.atan()),
    );
    let revolution = base.replace(
        "CYLINDRICAL_SURFACE('',#9,2.)",
        "SURFACE_OF_REVOLUTION('',#14,#601)",
    );
    for data in [cone, revolution] {
        let import = read(&file(&data));
        assert!((import.bodies[0].volume - 76.0 * std::f64::consts::PI / 3.0).abs() < 1e-9);
        assert!(!tessellate::tessellate(&import.bodies[0].body, 0.01)
            .unwrap()
            .triangles
            .is_empty());
    }
}

#[test]
fn deterministic_byte_mutations_never_panic_or_return_unchecked_topology() {
    let original = file(&box_data(false)).into_bytes();
    let mut seed = 47_u64;
    for _ in 0..1000 {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let index = (seed >> 16) as usize % original.len();
        let mut bytes = original.clone();
        bytes[index] = (seed >> 48) as u8;
        if let Ok(result) = step::read(&bytes, TOL) {
            for b in result.bodies {
                topology::check(&b.body, TOL).unwrap();
                assert!(b.volume > 0.0);
            }
        }
    }
}

#[test]
fn imported_box_accepts_an_existing_kernel_through_hole_without_conversion() {
    use gc_kernel::geom::{Curve, Vec3};
    let imported = read(&file(&box_data(false)));
    let tool = gc_kernel::extrude::extrude(
        &neutral::Profile {
            loops: vec![vec![Curve::Arc {
                c: Vec3::new(1.0, 1.5, 0.0),
                n: Vec3::new(0.0, 0.0, 1.0),
                x: Vec3::new(1.0, 0.0, 0.0),
                r: 0.25,
                t0: 0.0,
                t1: std::f64::consts::TAU,
            }]],
            n: Vec3::new(0.0, 0.0, 1.0),
            low: -1.0,
            high: 5.0,
        },
        TOL,
    )
    .unwrap();
    let cut = gc_kernel::boolean::cut(&imported.bodies[0].body, &tool, TOL).unwrap();
    topology::check(&cut, TOL).unwrap();
    assert!(
        (mass::measure(&cut, TOL).unwrap().volume - (24.0 - std::f64::consts::PI / 4.0)).abs()
            < 1e-9
    );
    assert!(!tessellate::tessellate(&cut, 0.01)
        .unwrap()
        .triangles
        .is_empty());
    assert_eq!(imported.bodies[0].volume, 24.0);
}
#[test]
fn worker_operation_reads_a_real_file_without_changing_it() {
    let path = std::env::temp_dir().join(format!("gck-step-{}-Ω.step", std::process::id()));
    let source = file(&box_data(false));
    std::fs::write(&path, &source).unwrap();
    let bytes = path.to_str().unwrap().as_bytes();
    let mut data = vec![1., bytes.len() as f64];
    data.extend(bytes.iter().map(|v| *v as f64));
    let request = protocol::encode_request(protocol::OP_STEP_READ, TOL, &data);
    let output = protocol::handle(&protocol::parse_request(&request).unwrap()).unwrap();
    assert_eq!(&output[..3], &[1., 1., 401.]);
    assert_eq!(std::fs::read_to_string(&path).unwrap(), source);
    std::fs::remove_file(path).unwrap();
    for bad in [
        vec![],
        vec![2., 0.],
        vec![1., 1., 0.],
        vec![1., 1., f64::NAN],
        vec![1., 2., 65.],
        vec![1., 1., 256.],
    ] {
        assert!(step::wire::read_op(&bad, TOL).is_err());
    }
}

#[test]
fn checked_in_real_step_corpus_has_explicit_results() {
    let Ok(root) = std::env::var("GC_STEP_CORPUS") else {
        return;
    };
    let admitted = [
        "233_a-s3_waveguide_bridge_hanger_brkt.STEP",
        "angle.step",
        "CChannel.step",
        "flag pole.step",
        "flat.step",
        "Hour Glass Flat.stp",
        "i-beam.step",
        "Nebraska.stp",
        "Pole.step",
        "Short Angle.step",
        "Simple Holes Pipe.step",
    ];
    let mut total = 0;
    for entry in std::fs::read_dir(&root).unwrap() {
        let path = entry.unwrap().path();
        if !matches!(
            path.extension()
                .and_then(|s| s.to_str())
                .map(str::to_ascii_lowercase)
                .as_deref(),
            Some("step" | "stp")
        ) {
            continue;
        }
        total += 1;
        let before = std::fs::read(&path).unwrap();
        let result = step::read_file(Path::new(&path), TOL);
        if admitted.contains(&path.file_name().unwrap().to_str().unwrap()) {
            let result = result.unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            assert_eq!(result.bodies.len(), 1);
            for item in result.bodies {
                topology::check(&item.body, TOL).unwrap();
                assert!(item.volume > 0.);
                // Respect source units even where historical fixtures declare metres.
                let mesh = tessellate::tessellate(&item.body, 0.02 * item.source_length_to_inches)
                    .unwrap();
                assert!(!mesh.triangles.is_empty());
            }
        } else {
            assert!(
                matches!(&result, Err(KernelError::Unsupported(_)))
                    || matches!(&result, Err(KernelError::InvalidInput(s)) if s.contains("Bezier span cannot be certified")),
                "{}: {result:?}",
                path.display()
            );
        }
        assert_eq!(std::fs::read(&path).unwrap(), before);
    }
    assert_eq!(total, 19);
}

fn spline_box(mm: bool, reverse: bool) -> String {
    let mut data = box_data(mm);
    let scale = if mm { 25.4 } else { 1.0 };
    let mut coords = [0.0, 0.2, 0.6, 1.0, 1.5, 2.0];
    if reverse {
        coords.reverse();
    }
    for (i, x) in coords.iter().enumerate() {
        data += &format!("#{}=CARTESIAN_POINT('',({},0.,0.));", 600 + i, x * scale);
    }
    data=data.replace("#102=LINE('',#1,#101);","#102=B_SPLINE_CURVE_WITH_KNOTS('',3,(#600,#601,#602,#603,#604,#605),.UNSPECIFIED.,.F.,.F.,(4,2,4),(0.,0.4,1.),.UNSPECIFIED.);");
    if reverse {
        data = data.replace(
            "#103=EDGE_CURVE('',#11,#12,#102,.T.);",
            "#103=EDGE_CURVE('',#11,#12,#102,.F.);",
        );
    }
    file(&data)
}

#[test]
fn spline_spans_share_topology_and_preserve_mass_units_and_orientation() {
    for mm in [false, true] {
        for reverse in [false, true] {
            let b = read(&spline_box(mm, reverse)).bodies.remove(0);
            assert_eq!(
                (
                    b.body.vertices.len(),
                    b.body.edges.len(),
                    b.body.faces.len()
                ),
                (9, 13, 6)
            );
            let report = topology::check(&b.body, TOL).unwrap();
            assert_eq!(report.genus(), 0.0);
            assert!((b.volume - 24.0).abs() < 1e-10);
            assert!((mass::measure(&b.body, TOL).unwrap().area - 52.0).abs() < 1e-10);
            assert_eq!(
                b.body
                    .edges
                    .iter()
                    .filter(|e| matches!(e.curve, gc_kernel::geom::Curve::Bezier { .. }))
                    .count(),
                2
            );
            let wire = neutral::write_body(&b.body);
            let mut decoder = neutral::Reader::with_tolerance(&wire, TOL);
            topology::check(&decoder.body().unwrap(), TOL).unwrap();
            decoder.finish().unwrap();
            assert!(!tessellate::tessellate(&b.body, 0.01)
                .unwrap()
                .triangles
                .is_empty());
        }
    }
}

#[test]
fn off_curve_or_malformed_spline_edges_fail_without_partial_body() {
    refused(
        &spline_box(false, false).replace("(0,0.,0.)", "(0.5,0.,0.)"),
        "point inversion found 0 distinct parameters",
    );
    refused(
        &spline_box(false, false).replace("(4,2,4)", "(4,1,4)"),
        "knot sum",
    );
}

// A source-defined polynomial ring within the EXISTING tolerance of its analytic
// cylinder. No conversion to circles: the source control points must survive.
fn spline_cylinder_data(spans: usize) -> String {
    let mut data = cylinder_data();
    for (curve, height, frame) in [(11, 0.0, 9), (12, 4.0, 10)] {
        let mut refs = Vec::new();
        let mut points = Vec::new();
        let step = std::f64::consts::TAU / spans as f64;
        let handle = 4.0 / 3.0 * (step / 4.0).tan();
        for i in 0..spans {
            let a = i as f64 * step;
            let b = (i + 1) as f64 * step;
            let start = (2.0 * a.cos(), 2.0 * a.sin());
            let end = (2.0 * b.cos(), 2.0 * b.sin());
            if i == 0 {
                points.push(start);
            }
            points.push((
                start.0 - 2.0 * handle * a.sin(),
                start.1 + 2.0 * handle * a.cos(),
            ));
            points.push((
                end.0 + 2.0 * handle * b.sin(),
                end.1 - 2.0 * handle * b.cos(),
            ));
            points.push(if i + 1 == spans { (2.0, 0.0) } else { end });
        }
        for (i, (x, y)) in points.iter().enumerate() {
            let id = 1000 + curve * 100 + i;
            refs.push(format!("#{id}"));
            data += &format!("#{id}=CARTESIAN_POINT('',({x},{y},{height}));");
        }
        let mult = (0..=spans)
            .map(|i| if i == 0 || i == spans { "4" } else { "3" })
            .collect::<Vec<_>>()
            .join(",");
        let knots = (0..=spans)
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join(",");
        data=data.replace(&format!("#{curve}=CIRCLE('',#{frame},2.);"),&format!("#{curve}=B_SPLINE_CURVE_WITH_KNOTS('',3,({}),.UNSPECIFIED.,.T.,.F.,({mult}),({knots}),.PIECEWISE_BEZIER_KNOTS.);",refs.join(",")));
    }
    data
}

#[test]
fn polynomial_spline_cylinder_keeps_source_curves_and_proves_whole_span_agreement() {
    let source = file(&spline_cylinder_data(16));
    let imported = read(&source);
    let body = &imported.bodies[0].body;
    assert_eq!(body.edges.len(), 33); // two rings of 16 spans and the shared seam
    assert_eq!(
        body.edges
            .iter()
            .filter(|e| matches!(e.curve, gc_kernel::geom::Curve::Bezier { .. }))
            .count(),
        32
    );
    topology::check(body, TOL).unwrap();
    let measured = mass::measure(body, TOL).unwrap();
    assert!(!measured.exact); // curved face integral uses quadrature, never claim closed form
    assert!(
        (measured.volume - 16.0 * std::f64::consts::PI).abs() < 2e-6,
        "{measured:?}"
    );
    assert!(
        (measured.area - 24.0 * std::f64::consts::PI).abs() < 2e-6,
        "{measured:?}"
    );
    assert!(!tessellate::tessellate(body, 0.01)
        .unwrap()
        .triangles
        .is_empty());
    // Exact cubic quadrant has a larger radial bulge, even though endpoints lie
    // on the cylinder. Do not accept it through endpoint-only or loose checks.
    refused(&file(&spline_cylinder_data(4)), "cannot be certified");
    // Only the tolerance explicitly requested by the caller may change admission.
    assert!(step::read(source.as_bytes(), 1e-10).is_err());
}

#[test]
fn nonuniform_parameter_trimmed_polynomial_box_preserves_units_and_shared_edge_chain() {
    for mm in [false, true] {
        for reversed in [false, true] {
            let scale = if mm { 25.4 } else { 1.0 };
            let mut data = box_data(mm);
            for (i, x) in [-1., 1., 3.].iter().enumerate() {
                data += &format!("#{}=CARTESIAN_POINT('',({},0.,0.));", 601 + i, x * scale);
            }
            data += "#600=B_SPLINE_CURVE_WITH_KNOTS('',1,(#601,#602,#603),.POLYLINE_FORM.,.F.,.F.,(2,1,2),(-10.,2.,30.),.UNSPECIFIED.);";
            let trim = if reversed {
                "#102=TRIMMED_CURVE('',#600,(PARAMETER_VALUE(16.),#2),(PARAMETER_VALUE(-4.),#1),.F.,.PARAMETER.);"
            } else {
                "#102=TRIMMED_CURVE('',#600,(PARAMETER_VALUE(-4.),#1),(PARAMETER_VALUE(16.),#2),.T.,.PARAMETER.);"
            };
            data = data.replace("#102=LINE('',#1,#101);", trim);
            if reversed {
                data = data.replace(
                    "#103=EDGE_CURVE('',#11,#12,#102,.T.);",
                    "#103=EDGE_CURVE('',#11,#12,#102,.F.);",
                );
            }
            let imported = read(&file(&data));
            let b = &imported.bodies[0].body;
            assert_eq!((b.vertices.len(), b.edges.len(), b.faces.len()), (9, 13, 6));
            assert_eq!(topology::check(b, TOL).unwrap().genus(), 0.);
            let m = mass::measure(b, TOL).unwrap();
            assert!((m.volume - 24.).abs() < 1e-10);
            assert!((m.area - 52.).abs() < 1e-10);
            assert!(!tessellate::tessellate(b, 0.01)
                .unwrap()
                .triangles
                .is_empty());
        }
    }
}

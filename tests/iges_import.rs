use gc_kernel::{
    extrude,
    geom::{Curve, Surface, Vec3},
    iges, mass,
    neutral::Profile,
    tessellate, topology,
};
const TOL: f64 = 1e-6;
struct Fixture {
    entries: Vec<(usize, usize, usize, String, String)>,
}
impl Fixture {
    fn add(&mut self, kind: usize, form: usize, transform: usize, data: String) -> usize {
        let id = self.entries.len() * 2 + 1;
        self.entries.push((
            kind,
            form,
            transform,
            format!("e{id}"),
            format!("{kind},{data};"),
        ));
        id
    }
    fn file(&self, mm: bool, model: f64) -> String {
        let record = |s: &str, section: char, n: usize| format!("{s:<72}{section}{n:7}\n");
        let mut out = record("Independent analytic IGES fixture", 'S', 1);
        let global=format!("1H,,1H;,8HGC, test,8Htest.igs,2HGC,2HGC,32,38,6,308,15,2HGC,{model},{},{},1,0.01,15H20261031.120000,0.000001,1000.,3HBob,2HGC,11,0,15H20261031.120000,;",if mm {2}else{1},if mm {"2HMM"}else{"2HIN"});
        let chunks = global.as_bytes().chunks(72).collect::<Vec<_>>();
        let ng = chunks.len();
        for (i, s) in chunks.iter().enumerate() {
            out += &record(std::str::from_utf8(s).unwrap(), 'G', i + 1);
        }
        let mut parameters = String::new();
        let mut np = 0;
        for (i, (kind, form, transform, label, data)) in self.entries.iter().enumerate() {
            let id = 2 * i + 1;
            let chunks = data.as_bytes().chunks(64).collect::<Vec<_>>();
            let first = format!(
                "{kind:8}{:8}{:8}{:8}{:8}{:8}{transform:8}{:8}{:08}",
                np + 1,
                0,
                0,
                0,
                0,
                0,
                if *kind == 186 {
                    0
                } else if *kind == 504 {
                    10001
                } else {
                    10000
                }
            );
            let second = format!(
                "{kind:8}{:8}{:8}{:8}{form:8}{:8}{:8}{label:<8}{:8}",
                0,
                0,
                chunks.len(),
                0,
                0,
                0
            );
            out += &record(&first, 'D', id);
            out += &record(&second, 'D', id + 1);
            for s in chunks {
                np += 1;
                parameters += &format!("{:<64}{id:8}P{np:7}\n", std::str::from_utf8(s).unwrap());
            }
        }
        out += &parameters;
        out += &record(
            &format!("S{:7}G{ng:7}D{:7}P{np:7}", 1, 2 * self.entries.len()),
            'T',
            1,
        );
        out
    }
}
fn box_fixture(mm: bool, model: f64, splines: bool, transform: bool) -> String {
    let corners = [
        Vec3::new(0., 0., 0.),
        Vec3::new(2., 0., 0.),
        Vec3::new(2., 3., 0.),
        Vec3::new(0., 3., 0.),
    ];
    let body = extrude::extrude(
        &Profile {
            loops: vec![(0..4)
                .map(|i| Curve::Line {
                    a: corners[i],
                    b: corners[(i + 1) % 4],
                })
                .collect()],
            n: Vec3::new(0., 0., 1.),
            low: 0.,
            high: 4.,
        },
        TOL,
    )
    .unwrap();
    let scale = if mm { 25.4 * model } else { model };
    let xyz = |p: Vec3| format!("{},{},{}", p.x * scale, p.y * scale, p.z * scale);
    let mut file = Fixture {
        entries: Vec::new(),
    };
    let vertices = file.add(
        502,
        1,
        0,
        format!(
            "{},{}",
            body.vertices.len(),
            body.vertices
                .iter()
                .map(|p| xyz(*p))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    let mut edges = Vec::new();
    for (i, e) in body.edges.iter().enumerate() {
        let c = if splines && i == 0 {
            let a = e.curve.start();
            let b = e.curve.end();
            let m = a.add(b).scale(0.5);
            file.add(
                126,
                0,
                0,
                format!(
                    "2,2,1,0,0,0,0.,0.,0.,1.,1.,1.,1.,0.7,1.,{},{},{},0.,1.,0.,0.,1.",
                    xyz(a),
                    xyz(m),
                    xyz(b)
                ),
            )
        } else {
            file.add(
                110,
                0,
                0,
                format!("{},{}", xyz(e.curve.start()), xyz(e.curve.end())),
            )
        };
        edges.push(format!(
            "{c},{vertices},{},{vertices},{}",
            e.start + 1,
            e.end + 1
        ));
    }
    let edge_list = file.add(504, 1, 0, format!("{},{}", edges.len(), edges.join(",")));
    let mut faces = Vec::new();
    for face in &body.faces {
        let Surface::Plane { p, n, x } = face.surface else {
            panic!()
        };
        let surface = if splines && p.z == 4. && n.z > 0.9 {
            let weights = (0..3)
                .flat_map(|j| {
                    (0..3)
                        .map(move |i| format!("{}", (1. + 0.2 * i as f64) * (1. + 0.1 * j as f64)))
                })
                .collect::<Vec<_>>()
                .join(",");
            let controls = (0..3)
                .flat_map(|j| (0..3).map(move |i| Vec3::new(i as f64, 1.5 * j as f64, 4.)))
                .map(xyz)
                .collect::<Vec<_>>()
                .join(",");
            file.add(128,0,0,format!("2,2,2,2,0,0,0,0,0,0.,0.,0.,1.,1.,1.,0.,0.,0.,1.,1.,1.,{weights},{controls},0.,1.,0.,1."))
        } else {
            let p = file.add(116, 0, 0, format!("{},0", xyz(p)));
            let n = file.add(123, 0, 0, format!("{},{},{}", n.x, n.y, n.z));
            let x = file.add(123, 0, 0, format!("{},{},{}", x.x, x.y, x.z));
            file.add(190, 1, 0, format!("{p},{n},{x}"))
        };
        let mut loops = Vec::new();
        for lp in &face.loops {
            let mut lp = lp.clone();
            if !face.sense {
                lp.reverse();
                for (_, forward) in &mut lp {
                    *forward = !*forward;
                }
            }
            loops.push(file.add(
                508,
                1,
                0,
                format!(
                        "{},{}",
                        lp.len(),
                        lp.iter()
                            .map(|(e, f)| format!("0,{edge_list},{},{},0", e + 1, usize::from(*f)))
                            .collect::<Vec<_>>()
                            .join(",")
                    ),
            ));
        }
        let f = file.add(
            510,
            1,
            0,
            format!(
                "{surface},{},1,{}",
                loops.len(),
                loops
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ),
        );
        faces.push(format!("{f},{}", usize::from(face.sense)));
    }
    let shell = file.add(514, 1, 0, format!("{},{}", faces.len(), faces.join(",")));
    let matrix = if transform {
        file.add(
            124,
            0,
            0,
            format!(
                "0.,-1.,0.,{},1.,0.,0.,{},0.,0.,1.,{}",
                10. * scale,
                20. * scale,
                30. * scale
            ),
        )
    } else {
        0
    };
    file.add(186, 0, matrix, format!("{shell},1,0"));
    file.file(mm, model)
}
#[test]
fn explicit_iges_solids_preserve_units_scale_splines_and_rigid_placement() {
    for (name, mm, scale, splines, transform) in [
        ("AnalyticBox.igs", false, 1., false, false),
        ("SplineBoxMm.igs", true, 0.5, true, false),
        ("SplineBoxMmUnitScale.igs", true, 1., true, false),
        ("PlacedBox.igs", false, 1., false, true),
    ] {
        let text = box_fixture(mm, scale, splines, transform);
        let items = iges::read(text.as_bytes(), TOL).unwrap();
        assert_eq!(items.len(), 1);
        let b = &items[0];
        assert!((b.volume - 24.).abs() < 1e-9);
        assert!((mass::measure(&b.body, TOL).unwrap().area - 52.).abs() < 1e-8);
        topology::check(&b.body, TOL).unwrap();
        assert!(!tessellate::tessellate(&b.body, 0.01)
            .unwrap()
            .triangles
            .is_empty());
        if transform {
            assert!(b
                .body
                .vertices
                .iter()
                .any(|p| p.distance(Vec3::new(10., 20., 30.)) < 1e-12));
        }
        if let Ok(dir) = std::env::var("GC_IGES_EXPORT_FIXTURES") {
            std::fs::write(std::path::Path::new(&dir).join(name), text).unwrap();
        }
    }
}
#[test]
fn iges_malformed_records_topology_and_unsupported_geometry_are_atomic_errors() {
    let text = box_fixture(false, 1., false, false);
    for count in [0, 1, 79, 80, text.len() - 1 - 81] {
        assert!(iges::read(&text.as_bytes()[..count], TOL).is_err());
    }
    let mut wrong = text.clone().into_bytes();
    wrong[79] = b'2';
    assert!(iges::read(&wrong, TOL).is_err());
    assert!(iges::read(text.replace("1H,,1H;", "1H,,1H,").as_bytes(), TOL).is_err());
    assert!(iges::read(text.replace("186,", "144,").as_bytes(), TOL).is_err());
    let mismatch = text.replace("110,0,0,0,2,0,0;", "110,1,0,0,2,0,0;");
    assert_ne!(mismatch, text);
    assert!(iges::read(mismatch.as_bytes(), TOL)
        .unwrap_err()
        .to_string()
        .contains("bounded curve endpoints"));
    assert!(iges::read(text.as_bytes(), f64::NAN).is_err());
    assert!(iges::wire::read_op(&[2., 1., 65.], TOL).is_err());
}

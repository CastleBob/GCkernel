use gc_kernel::{
    geom::Vec3,
    mass, neutral, step,
    stl::{self, Encoding, Facet, Mesh},
    tessellate, topology,
};
const TOL: f64 = 1e-6;
fn tetra() -> Mesh {
    let p = [
        Vec3::new(0., 0., 0.),
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    ];
    Mesh {
        encoding: Encoding::Ascii,
        facets: [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
            .map(|t| Facet {
                normal: Vec3::new(0., 0., 0.),
                vertices: t.map(|i| p[i]),
            })
            .to_vec(),
    }
}
fn ascii(m: &Mesh) -> Vec<u8> {
    let mut s = "solid independent tetrahedron\n".to_owned();
    for f in &m.facets {
        s += &format!(
            "facet normal {} {} {}\nouter loop\n",
            f.normal.x, f.normal.y, f.normal.z
        );
        for p in f.vertices {
            s += &format!("vertex {} {} {}\n", p.x, p.y, p.z);
        }
        s += "endloop\nendfacet\n";
    }
    s += "endsolid independent tetrahedron\n";
    s.into_bytes()
}
fn binary(m: &Mesh) -> Vec<u8> {
    let mut b = vec![0; 80];
    b[..5].copy_from_slice(b"solid");
    b.extend((m.facets.len() as u32).to_le_bytes());
    for f in &m.facets {
        for p in [f.normal, f.vertices[0], f.vertices[1], f.vertices[2]] {
            for v in [p.x, p.y, p.z] {
                b.extend((v as f32).to_le_bytes());
            }
        }
        b.extend([0, 0]);
    }
    b
}
#[test]
fn ascii_binary_explicit_units_and_exact_faceted_solid() {
    for (encoding, bytes) in [
        (Encoding::Ascii, ascii(&tetra())),
        (Encoding::Binary, binary(&tetra())),
    ] {
        for scale in [1., 1. / 25.4] {
            let m = stl::read(&bytes, scale).unwrap();
            assert_eq!(m.encoding, encoding);
            let b = stl::to_solid(&m, TOL).unwrap();
            assert_eq!(topology::check(&b, TOL).unwrap().genus(), 0.);
            assert!((mass::measure(&b, TOL).unwrap().volume - scale.powi(3) / 6.).abs() < 1e-12);
            tessellate::tessellate(&b, 0.001).unwrap();
            let exported = step::write(&[("STL tetrahedron", &b)], TOL).unwrap();
            assert_eq!(
                step::read(exported.as_bytes(), TOL).unwrap().bodies.len(),
                1
            );
            assert_eq!(b.vertices.len(), 4);
            assert_eq!(b.edges.len(), 6);
            assert_eq!(b.faces.len(), 4);
            if scale == 1. {
                if let Ok(dir) = std::env::var("GC_STL_EXPORT_FIXTURES") {
                    std::fs::write(
                        std::path::Path::new(&dir).join(if encoding == Encoding::Ascii {
                            "TetraAscii.stl"
                        } else {
                            "TetraBinary.stl"
                        }),
                        &bytes,
                    )
                    .unwrap();
                }
            }
        }
    }
}
#[test]
fn mesh_intake_never_silently_repairs_solid_failures() {
    let original = tetra();
    for kind in 0..5 {
        let mut m = original.clone();
        match kind {
            0 => {
                m.facets.pop();
            }
            1 => {
                m.facets[0].vertices.swap(1, 2);
            }
            2 => {
                for f in &mut m.facets {
                    f.vertices.swap(1, 2);
                }
            }
            3 => {
                m.facets[0].vertices[0].x += 1e-8;
            }
            _ => {
                let mut other = original.clone();
                for f in &mut other.facets {
                    for p in &mut f.vertices {
                        *p = p.add(Vec3::new(2., 0., 0.));
                    }
                }
                m.facets.extend(other.facets);
            }
        }
        let mesh = stl::read(&ascii(&m), 1.).unwrap();
        assert!(stl::to_solid(&mesh, TOL).is_err(), "kind {kind}");
    }
    assert!(stl::to_solid(&original, TOL).is_ok());
}
#[test]
fn malformed_numbers_records_attributes_and_truncation_are_errors() {
    let b = binary(&tetra());
    for n in 0..b.len() {
        assert!(stl::read(&b[..n], 1.).is_err());
    }
    for scale in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(stl::read(&b, scale).is_err());
    }
    let mut color = b.clone();
    color[132] = 1;
    assert!(stl::read(&color, 1.).is_err());
    let mut nan = b.clone();
    nan[96..100].copy_from_slice(&f32::NAN.to_le_bytes());
    assert!(stl::read(&nan, 1.).is_err());
    let text = String::from_utf8(ascii(&tetra())).unwrap();
    for bad in [
        text.replace("outer loop", "inner loop"),
        text.replace("endfacet", ""),
        text.replace("vertex 0 0 0", "vertex NaN 0 0"),
        text.replace("facet normal 0 0 0", "facet normal 0 0 2"),
        format!("{text}garbage"),
    ] {
        assert!(stl::read(bad.as_bytes(), 1.).is_err());
    }
}
#[test]
fn coplanar_shared_edges_and_large_translation_are_admitted() {
    let b = gc_kernel::step::read(
        include_bytes!("fixtures/vertex_trimmed_spline_box.step"),
        TOL,
    )
    .unwrap()
    .bodies
    .remove(0)
    .body;
    let mesh = tessellate::tessellate(&b, 0.01).unwrap();
    let facets = mesh
        .triangles
        .iter()
        .map(|t| Facet {
            normal: Vec3::new(0., 0., 0.),
            vertices: t.map(|i| mesh.positions[i].add(Vec3::new(1e6, -2e6, 3e6))),
        })
        .collect();
    let m = Mesh {
        encoding: Encoding::Ascii,
        facets,
    };
    let b = stl::to_solid(&m, TOL).unwrap();
    assert!((mass::measure(&b, TOL).unwrap().volume - 24.).abs() < 1e-8);
    assert!(!neutral::write_body(&b).is_empty());
}

#[test]
fn closed_but_crossing_and_pinched_shells_are_refused() {
    let points = [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(-1., 0., 0.),
        Vec3::new(0., -1., 0.),
        Vec3::new(0.8, 0.8, -0.5),
        Vec3::new(0., 0., -1.),
    ];
    let mut facets = vec![];
    for i in 0..4 {
        for ids in [[i, (i + 1) % 4, 4], [(i + 1) % 4, i, 5]] {
            facets.push(Facet {
                normal: Vec3::new(0., 0., 0.),
                vertices: ids.map(|k| points[k]),
            });
        }
    }
    let error = stl::to_solid(
        &Mesh {
            encoding: Encoding::Ascii,
            facets,
        },
        TOL,
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("facets"), "{error}");
    let mut pinch = tetra();
    let mut other = tetra();
    for f in &mut other.facets {
        for p in &mut f.vertices {
            *p = p.scale(-1.);
        }
        f.vertices.swap(1, 2);
    }
    pinch.facets.extend(other.facets);
    assert!(stl::to_solid(&pinch, TOL)
        .unwrap_err()
        .to_string()
        .contains("pinched"));
}


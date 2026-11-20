//! SU2 reader (D187): the fixtures of `tools/mesh/su2_test_meshes.py` load
//! through the format layer, cells are closed, volumes and boundary areas
//! are exact, markers name the boundaries, 2D files are extruded one layer.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::case::{MeshSection, build_mesh};
use foehncraft_mesh::format::{self, ReadOptions, su2};
use foehncraft_mesh::quality::QualityReport;

fn path(name: &str) -> String {
    format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn check(f: &format::MeshFile, volume: f64) -> QualityReport {
    f.mesh.validate().unwrap_or_else(|e| panic!("{e}"));
    let geo = Geometry::compute(&f.mesh);
    let (q, _) = QualityReport::compute(&f.mesh, &geo);
    assert!(q.ok(), "{}", q.to_markdown("su2"));
    assert!(
        ((q.total_volume - volume) / volume).abs() < 1e-12,
        "{}",
        q.total_volume
    );
    assert!(((q.boundary_volume - volume) / volume).abs() < 1e-12);
    q
}

fn area(q: &QualityReport, name: &str) -> f64 {
    q.boundaries
        .iter()
        .find(|b| b.0 == name)
        .unwrap_or_else(|| panic!("no boundary {name}"))
        .2
}

#[test]
fn mixed_3d_cells_and_markers() {
    let f = format::read(path("cubes_mixed_3d.su2"), None).unwrap();
    assert_eq!(f.reader, "su2");
    assert_eq!(f.format, "SU2 3D ASCII");
    let q = check(&f, 4.0);
    assert_eq!(f.mesh.n_cells, 20);
    for s in ["hexahedra", "prisms", "pyramids", "tetrahedra"] {
        assert!(q.shapes.iter().any(|x| x.0 == s), "{s}: {:?}", q.shapes);
    }
    let names: Vec<&str> = f.mesh.boundaries.iter().map(|b| b.name.as_str()).collect();
    assert_eq!(names, ["xmin", "xmax", "walls"]);
    assert!((area(&q, "xmin") - 1.0).abs() < 1e-12);
    assert!((area(&q, "xmax") - 1.0).abs() < 1e-12);
    assert!((area(&q, "walls") - 16.0).abs() < 1e-12);
    assert_eq!(f.mesh.zones, ["fluid"]);
    // the FFD section after the mesh is skipped with a note
    assert_eq!(f.warnings.len(), 1, "{:?}", f.warnings);
    assert!(f.warnings[0].contains("FFD_NBOX"), "{:?}", f.warnings);
}

#[test]
fn channel_2d_is_extruded_one_layer() {
    let f = format::read(path("channel_2d.su2"), None).unwrap();
    assert_eq!(f.format, "SU2 2D ASCII");
    // default thickness: the median edge length (0.5 here)
    let q = check(&f, 2.0 * 0.5);
    assert_eq!(f.mesh.n_cells, 4 + 8);
    assert!(q.shapes.iter().any(|x| x.0 == "hexahedra"));
    assert!(q.shapes.iter().any(|x| x.0 == "prisms"));
    let names: Vec<&str> = f.mesh.boundaries.iter().map(|b| b.name.as_str()).collect();
    assert_eq!(names, ["inlet", "outlet", "lower", "upper", "zmin", "zmax"]);
    assert!((area(&q, "zmin") - 2.0).abs() < 1e-12);
    assert!((area(&q, "inlet") - 0.5).abs() < 1e-12);
    assert!((area(&q, "lower") - 1.0).abs() < 1e-12);
    assert!(f.warnings.is_empty(), "{:?}", f.warnings);

    let opts = ReadOptions {
        thickness: Some(0.01),
    };
    let g = format::read_with(path("channel_2d.su2"), None, &opts).unwrap();
    check(&g, 2.0 * 0.01);
}

#[test]
fn case_reads_su2_with_renamed_markers_and_scale() {
    // the case path: format by extension, markers renamed to the BC names
    let sec = MeshSection {
        file: Some(path("channel_2d.su2").into()),
        thickness: Some(0.1),
        scale: 0.001,
        boundaries: [("lower", "wall"), ("upper", "top")]
            .into_iter()
            .map(|(a, b)| (a.to_string(), b.to_string()))
            .collect(),
        ..Default::default()
    };
    let (mesh, source, warnings) = build_mesh(&sec).unwrap();
    assert!(source.ends_with("(SU2 2D ASCII)"), "{source}");
    assert!(warnings.is_empty());
    let names: Vec<&str> = mesh.boundaries.iter().map(|b| b.name.as_str()).collect();
    assert_eq!(names, ["inlet", "outlet", "wall", "top", "zmin", "zmax"]);
    let geo = Geometry::compute(&mesh);
    let v: f64 = geo.cell_volume.iter().sum();
    assert!((v - 2.0 * 0.1 * 1e-9).abs() < 1e-21, "{v}");

    // the format forced by name on a file without the extension
    let dir = std::env::temp_dir().join(format!("su2-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let copy = dir.join("channel.mesh");
    std::fs::copy(path("channel_2d.su2"), &copy).unwrap();
    // by content
    assert_eq!(format::read(&copy, None).unwrap().reader, "su2");
    // by name
    assert_eq!(format::read(&copy, Some("su2")).unwrap().reader, "su2");
    let e = format::read(&copy, Some("xyz")).unwrap_err();
    assert!(
        e.contains("unknown mesh format 'xyz'") && e.contains("gmsh, su2, cgns, openfoam"),
        "{e}"
    );
    let e = format::read(&copy, Some("gmsh")).unwrap_err();
    assert!(e.contains("expected a section"), "{e}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn bad_files_give_clear_errors() {
    let o = ReadOptions::default();
    let e = |s: &str| su2::parse(s.as_bytes(), &o).unwrap_err();
    assert!(e("NDIME= 3\nNZONE= 2\n").contains("multi-zone"));
    assert!(e("NDIME= 4\n").contains("2 or 3"));
    assert!(e("NDIME= 3\nNELEM= 1\n24 0 1 2 3 4 5 6 7 8 9\n").contains("higher-order"));
    assert!(e("NDIME= 3\nNELEM= 1\n10 0 1 2\n").contains("needs 4 nodes"));
    let s = "NDIME= 3\nNELEM= 1\n10 0 1 2 7\nNPOIN= 4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\nNMARK= 0\n";
    assert!(e(s).contains("node 7"), "{}", e(s));
    let s = "NDIME= 2\nNELEM= 1\n10 0 1 2 3\nNPOIN= 4\n0 0\n1 0\n0 1\n1 1\nNMARK= 0\n";
    assert!(e(s).contains("not a 2D cell"), "{}", e(s));
    assert!(e("NDIME= 3\nhello\n").contains("line 2"));

    // a single tetrahedron without markers: its faces go to `default`
    let s = "NDIME= 3\nNELEM= 1\n10 0 1 2 3\nNPOIN= 4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n";
    let f = su2::parse(s.as_bytes(), &o).unwrap();
    assert_eq!(f.mesh.boundaries.len(), 1);
    assert_eq!(f.mesh.boundaries[0].name, "default");
    assert!(
        f.warnings.iter().any(|w| w.contains("no NMARK")),
        "{:?}",
        f.warnings
    );
    assert!(
        f.warnings.iter().any(|w| w.contains("belong to no marker")),
        "{:?}",
        f.warnings
    );
}

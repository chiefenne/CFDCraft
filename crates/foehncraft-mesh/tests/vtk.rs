//! VTK reader and writer (D192): legacy and XML files written by an
//! independent writer, every encoding, regions as boundaries and zones,
//! polyhedra and the round trip through [`write_vtu`].

use foehncraft_mesh::Geometry;
use foehncraft_mesh::format::vtk::write_vtu;
use foehncraft_mesh::format::{self, MeshFile};
use foehncraft_mesh::quality::QualityReport;

fn data(name: &str) -> String {
    format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn quality(f: &MeshFile) -> QualityReport {
    f.mesh.validate().unwrap_or_else(|e| panic!("{e}"));
    let geo = Geometry::compute(&f.mesh);
    let (q, _) = QualityReport::compute(&f.mesh, &geo);
    assert!(q.ok(), "{}", q.to_markdown("vtk"));
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
fn every_encoding_reads_the_same_mesh() {
    // the Gmsh file the VTK files were made from is the reference
    let g = format::read(data("channel_mixed_ascii.msh"), None).unwrap();
    let qg = quality(&g);
    for (file, how) in [
        ("channel_mixed_ascii.vtk", "VTK legacy ASCII"),
        ("channel_mixed_binary.vtk", "VTK legacy binary"),
        ("channel_mixed_ascii.vtu", "VTK XML (ascii)"),
        ("channel_mixed_binary.vtu", "VTK XML (binary)"),
        ("channel_mixed_zlib.vtu", "VTK XML (binary, zlib)"),
    ] {
        let f = format::read(data(&format!("vtk/{file}")), None).unwrap();
        assert_eq!(f.reader, "vtk", "{file}");
        assert_eq!(f.format, how, "{file}");
        let q = quality(&f);
        assert_eq!(f.mesh.n_cells, g.mesh.n_cells, "{file}");
        assert!((q.total_volume - qg.total_volume).abs() < 1e-12, "{file}");
        // the region numbers are the Gmsh physical tags
        let mut names: Vec<&str> = f.mesh.boundaries.iter().map(|b| b.name.as_str()).collect();
        names.sort();
        assert_eq!(
            names,
            [
                "boundary3",
                "boundary4",
                "boundary5",
                "boundary6",
                "boundary7"
            ],
            "{file}"
        );
        for (tag, name) in [
            (3, "back"),
            (4, "front"),
            (5, "inlet"),
            (6, "outlet"),
            (7, "walls"),
        ] {
            let a = area(&qg, name);
            assert!(
                (area(&q, &format!("boundary{tag}")) - a).abs() < 1e-12,
                "{file} {name}"
            );
        }
        assert_eq!(f.mesh.zones, ["zone1", "zone2"], "{file}");
        assert!(f.warnings.is_empty(), "{file}: {:?}", f.warnings);
    }
}

#[test]
fn round_trip_through_write_vtu() {
    let dir = std::env::temp_dir().join(format!("vtu-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    // hexahedra and prisms, and polyhedra, with their names
    for src in [
        data("channel_mixed_ascii.msh"),
        data("openfoam/dual_polyhedra"),
        data("openfoam/channel_ascii"),
    ] {
        let f = format::read(&src, None).unwrap();
        let qf = quality(&f);
        for binary in [false, true] {
            let p = dir.join(format!("m{binary}.vtu"));
            write_vtu(&f.mesh, &p, binary).unwrap();
            let r = format::read(&p, None).unwrap();
            let q = quality(&r);
            assert_eq!(r.mesh.n_cells, f.mesh.n_cells, "{src}");
            assert!((q.total_volume - qf.total_volume).abs() < 1e-12 * qf.total_volume);
            let names = |m: &MeshFile| {
                m.mesh
                    .boundaries
                    .iter()
                    .map(|b| (b.name.clone(), b.len))
                    .collect::<Vec<_>>()
            };
            assert_eq!(names(&r), names(&f), "{src}");
            assert_eq!(r.mesh.zones, f.mesh.zones, "{src}");
            for b in &f.mesh.boundaries {
                assert!((area(&q, &b.name) - area(&qf, &b.name)).abs() < 1e-12);
            }
            assert_eq!(r.mesh.n_faces(), f.mesh.n_faces(), "{src}");
        }
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn bad_files_give_clear_errors() {
    use format::MeshFormat;
    let e = |s: &str| {
        format::vtk::Vtk
            .parse(s.as_bytes(), &Default::default())
            .unwrap_err()
    };
    assert!(e("# vtk DataFile Version 3.0\nt\nASCII\nDATASET POLYDATA\n").contains("UNSTRUCTURED"));
    assert!(
        e("# vtk DataFile Version 3.0\nt\nASCII\nDATASET UNSTRUCTURED_GRID\nWHAT 1\n")
            .contains("WHAT")
    );
    assert!(e("<VTKFile type=\"ImageData\"></VTKFile>").contains("UnstructuredGrid"));
    assert!(
        e("<VTKFile type=\"UnstructuredGrid\" compressor=\"vtkLZ4DataCompressor\"></VTKFile>")
            .contains("LZ4")
    );
    // a grid of triangles only
    let s = "# vtk DataFile Version 3.0\nt\nASCII\nDATASET UNSTRUCTURED_GRID\n\
             POINTS 3 float\n0 0 0 1 0 0 0 1 0\nCELLS 1 4\n3 0 1 2\nCELL_TYPES 1\n5\n";
    assert!(e(s).contains("no 3D cells"), "{}", e(s));
}

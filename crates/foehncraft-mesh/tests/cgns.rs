#![cfg(feature = "hdf5")]
//! CGNS reader (D188): HDF5 files written by Gmsh (hexahedra, tetrahedra,
//! prisms) and a polyhedral NGON_n/NFACE_n file; Gmsh physical names become the
//! boundary names.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::format::{self, MeshFile};
use foehncraft_mesh::quality::QualityReport;

fn path(name: &str) -> String {
    format!("{}/tests/data/cgns/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn check(f: &MeshFile) -> QualityReport {
    f.mesh.validate().unwrap_or_else(|e| panic!("{e}"));
    let geo = Geometry::compute(&f.mesh);
    let (q, _) = QualityReport::compute(&f.mesh, &geo);
    assert!(q.ok(), "{}", q.to_markdown("cgns"));
    assert!((q.total_volume - 1.0).abs() < 1e-12, "{}", q.total_volume);
    assert!((q.boundary_volume - 1.0).abs() < 1e-12);
    for b in &q.boundaries {
        assert!((b.2 - 1.0).abs() < 1e-12, "{}: {}", b.0, b.2);
    }
    q
}

#[test]
fn gmsh_boxes() {
    for (file, cells, shape) in [
        ("box_hex.cgns", 8, "hexahedra"),
        ("box_tet.cgns", 184, "tetrahedra"),
        ("box_prism.cgns", 78, "prisms"),
    ] {
        let f = format::read(path(file), None).unwrap();
        assert_eq!(f.reader, "cgns");
        assert_eq!(f.format, "CGNS HDF5");
        let q = check(&f);
        assert_eq!(f.mesh.n_cells, cells, "{file}");
        assert!(
            q.shapes.iter().any(|x| x.0 == shape),
            "{file}: {:?}",
            q.shapes
        );
        let mut names: Vec<&str> = f.mesh.boundaries.iter().map(|b| b.name.as_str()).collect();
        names.sort();
        assert_eq!(
            names,
            ["xmax", "xmin", "ymax", "ymin", "zmax", "zmin"],
            "{file}"
        );
    }
}

#[test]
fn detected_by_content_and_name() {
    let p = path("box_hex.cgns");
    assert!(format::is_mesh_file(std::path::Path::new(&p)));
    assert_eq!(format::read(&p, Some("cgns")).unwrap().reader, "cgns");
}

#[test]
fn not_hdf5_is_refused() {
    let dir = std::env::temp_dir().join(format!("cgns-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("old.cgns");
    std::fs::write(&p, b"not hdf5 at all, an ADF file perhaps").unwrap();
    let e = format::read(&p, None).unwrap_err();
    assert!(e.contains("HDF5"), "{e}");
    std::fs::remove_dir_all(&dir).ok();
}

//! Multi-zone SU2 and CGNS files (D195): a rotor block and a stator block read
//! into one mesh, each a cell zone, the interface two boundaries (as the
//! sliding-mesh imprinting wants it), optionally merged into one connected mesh.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::format::{self, MeshFile, ReadOptions};
use foehncraft_mesh::quality::QualityReport;

fn path(name: &str) -> String {
    format!("{}/tests/data/multizone/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn area(f: &MeshFile, name: &str) -> f64 {
    let geo = Geometry::compute(&f.mesh);
    let (q, _) = QualityReport::compute(&f.mesh, &geo);
    q.boundaries
        .iter()
        .find(|b| b.0 == name)
        .unwrap_or_else(|| panic!("no boundary {name}"))
        .2
}

fn volume(f: &MeshFile) -> f64 {
    f.mesh.validate().unwrap_or_else(|e| panic!("{e}"));
    let geo = Geometry::compute(&f.mesh);
    let (q, _) = QualityReport::compute(&f.mesh, &geo);
    assert!(q.ok(), "{}", q.to_markdown("multizone"));
    q.total_volume
}

fn names(f: &MeshFile) -> Vec<String> {
    let mut v: Vec<String> = f.mesh.boundaries.iter().map(|b| b.name.clone()).collect();
    v.sort();
    v
}

fn kept_apart(f: &MeshFile, zones: [&str; 2], cells: [usize; 2]) {
    // two cell zones, no interior face between them: the interface is two boundaries
    assert!((volume(f) - 2.0).abs() < 1e-12);
    assert_eq!(
        f.mesh.zones.iter().map(|z| z.as_str()).collect::<Vec<_>>(),
        zones
    );
    assert_eq!(f.mesh.n_cells, cells[0] + cells[1]);
    assert_eq!(
        names(f),
        [
            "inlet",
            "outlet",
            "rotor_interface",
            "stator_interface",
            "wall"
        ]
    );
    assert!((area(f, "inlet") - 1.0).abs() < 1e-12);
    assert!((area(f, "outlet") - 1.0).abs() < 1e-12);
    assert!((area(f, "rotor_interface") - 1.0).abs() < 1e-12);
    assert!((area(f, "stator_interface") - 1.0).abs() < 1e-12);
    // the wall markers of both zones are one boundary: 4 sides of length 2
    assert!((area(f, "wall") - 8.0).abs() < 1e-12);
}

#[test]
fn su2_zones_kept_apart_non_matching() {
    let f = format::read(path("rotor_stator.su2"), None).unwrap();
    assert_eq!(f.reader, "su2");
    assert!(f.format.contains("2 zones"), "{}", f.format);
    kept_apart(&f, ["zone1", "zone2"], [8, 27]);
}

#[test]
fn su2_zones_merged_conforming() {
    let opts = ReadOptions {
        merge_nodes: true,
        ..Default::default()
    };
    let f = format::read_with(path("rotor_stator_conf.su2"), None, &opts).unwrap();
    // the interface is interior now; its two markers lie inside and are ignored
    assert!((volume(&f) - 2.0).abs() < 1e-12);
    assert_eq!(f.mesh.n_cells, 16);
    assert_eq!(names(&f), ["inlet", "outlet", "wall"]);
    assert!((area(&f, "wall") - 8.0).abs() < 1e-12);
    // 3 x 3 x 5 nodes of the joined block against 2 x 27 for two blocks
    assert_eq!(f.mesh.points.len(), 45);
    assert!(f.warnings.iter().any(|w| w.contains("merged")));

    // the same file without merging: 54 nodes, interface boundaries kept
    let g = format::read(path("rotor_stator_conf.su2"), None).unwrap();
    assert_eq!(g.mesh.points.len(), 54);
    assert!((volume(&g) - 2.0).abs() < 1e-12);
}

#[test]
fn su2_zone_count_and_errors() {
    let text = std::fs::read_to_string(path("rotor_stator_conf.su2")).unwrap();
    // a wrong NZONE is a warning, all zones are read
    let bad = text.replace("NZONE= 2", "NZONE= 3");
    let f = format::su2::parse(bad.as_bytes(), &ReadOptions::default()).unwrap();
    assert!(f.warnings.iter().any(|w| w.contains("NZONE= 3")));
    // a 2D multi-zone file is refused
    let two_d = text.replace("NDIME= 3", "NDIME= 2");
    assert!(format::su2::parse(two_d.as_bytes(), &ReadOptions::default()).is_err());
}

#[cfg(feature = "hdf5")]
#[test]
fn cgns_zones() {
    let f = format::read(path("rotor_stator.cgns"), None).unwrap();
    assert_eq!(f.reader, "cgns");
    kept_apart(&f, ["rotor", "stator"], [8, 8]);

    let opts = ReadOptions {
        merge_nodes: true,
        ..Default::default()
    };
    let g = format::read_with(path("rotor_stator.cgns"), None, &opts).unwrap();
    assert!((volume(&g) - 2.0).abs() < 1e-12);
    assert_eq!(g.mesh.n_cells, 16);
    assert_eq!(names(&g), ["inlet", "outlet", "wall"]);
    assert_eq!(g.mesh.points.len(), 45);
}

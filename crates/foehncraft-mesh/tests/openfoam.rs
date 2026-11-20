//! OpenFOAM reader (D188): polyMesh directories in ASCII and binary, polyhedral
//! cells, patches, cell zones, and clear errors.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::case::{MeshSection, build_mesh};
use foehncraft_mesh::format::{self, MeshFile, ReadOptions, openfoam};
use foehncraft_mesh::quality::QualityReport;

fn path(name: &str) -> String {
    format!("{}/tests/data/openfoam/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn check(f: &MeshFile, volume: f64) -> QualityReport {
    f.mesh.validate().unwrap_or_else(|e| panic!("{e}"));
    let geo = Geometry::compute(&f.mesh);
    let (q, _) = QualityReport::compute(&f.mesh, &geo);
    assert!(q.ok(), "{}", q.to_markdown("openfoam"));
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

fn names(f: &MeshFile) -> Vec<&str> {
    f.mesh.boundaries.iter().map(|b| b.name.as_str()).collect()
}

#[test]
fn ascii_channel() {
    let f = format::read(path("channel_ascii"), None).unwrap();
    assert_eq!(f.reader, "openfoam");
    assert!(
        f.format.starts_with("OpenFOAM polyMesh ASCII"),
        "{}",
        f.format
    );
    let q = check(&f, 2.0);
    assert_eq!(f.mesh.n_cells, 16);
    assert_eq!(names(&f), ["inlet", "outlet", "walls", "sides"]);
    assert!((area(&q, "inlet") - 1.0).abs() < 1e-12);
    assert!((area(&q, "walls") - 4.0).abs() < 1e-12);
    assert!((area(&q, "sides") - 4.0).abs() < 1e-12);
    assert_eq!(f.mesh.zones, ["right", "left"]);
    // empty patches are kept as boundaries, with a hint
    assert!(
        f.warnings.iter().any(|w| w.contains("empty")),
        "{:?}",
        f.warnings
    );
}

#[test]
fn binary_equals_ascii() {
    let a = format::read(path("channel_ascii"), None).unwrap();
    let b = format::read(path("channel_binary"), None).unwrap();
    assert!(b.format.contains("binary"), "{}", b.format);
    assert_eq!(a.mesh.n_cells, b.mesh.n_cells);
    assert_eq!(names(&a), names(&b));
    assert_eq!(a.mesh.zones, b.mesh.zones);
    let (ga, gb) = (Geometry::compute(&a.mesh), Geometry::compute(&b.mesh));
    for (x, y) in ga.cell_volume.iter().zip(&gb.cell_volume) {
        assert!((x - y).abs() < 1e-14);
    }
    check(&b, 2.0);
}

#[test]
fn polyhedral_dual_mesh() {
    let f = format::read(path("dual_polyhedra"), None).unwrap();
    let q = check(&f, 1.0);
    assert_eq!(f.mesh.n_cells, 81);
    assert_eq!(names(&f), ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]);
    for n in names(&f) {
        assert!((area(&q, n) - 1.0).abs() < 1e-12, "{n}");
    }
}

#[test]
fn directory_forms_and_case_path() {
    // the case directory, constant/ and polyMesh itself all work
    let base = path("channel_ascii");
    for p in [
        base.clone(),
        format!("{base}/constant"),
        format!("{base}/constant/polyMesh"),
    ] {
        assert!(format::is_mesh_file(std::path::Path::new(&p)), "{p}");
        assert_eq!(format::read(&p, None).unwrap().mesh.n_cells, 16, "{p}");
    }
    let sec = MeshSection {
        file: Some(base.into()),
        scale: 0.1,
        boundaries: [("walls", "wall")]
            .into_iter()
            .map(|(a, b)| (a.to_string(), b.to_string()))
            .collect(),
        ..Default::default()
    };
    let (mesh, source, _) = build_mesh(&sec).unwrap();
    assert!(source.contains("OpenFOAM"), "{source}");
    let n: Vec<&str> = mesh.boundaries.iter().map(|b| b.name.as_str()).collect();
    assert_eq!(n, ["inlet", "outlet", "wall", "sides"]);
    let v: f64 = Geometry::compute(&mesh).cell_volume.iter().sum();
    assert!((v - 2.0e-3).abs() < 1e-15, "{v}");
}

/// Copies a mesh into a fresh directory, changing one file.
fn variant(tag: &str, change: impl Fn(&str, String) -> String) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("foam-{tag}-{}", std::process::id()));
    let pm = dir.join("constant/polyMesh");
    std::fs::create_dir_all(&pm).unwrap();
    let src = format!("{}/constant/polyMesh", path("channel_ascii"));
    for n in [
        "points",
        "faces",
        "owner",
        "neighbour",
        "boundary",
        "cellZones",
    ] {
        let t = std::fs::read_to_string(format!("{src}/{n}")).unwrap();
        std::fs::write(pm.join(n), change(n, t)).unwrap();
    }
    dir
}

#[test]
fn bad_meshes_give_clear_errors() {
    let o = ReadOptions::default();
    let _ = &o;
    let e = |d: &std::path::Path| format::read(d, None).unwrap_err();

    let d = std::env::temp_dir().join(format!("foam-none-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    assert!(
        e(&d).contains("polyMesh") || e(&d).contains("mesh"),
        "{}",
        e(&d)
    );

    // a decomposed case
    let d = variant("proc", |n, t| {
        if n == "boundary" {
            t.replace("type            wall;", "type            processor;")
        } else {
            t
        }
    });
    assert!(e(&d).contains("processor"), "{}", e(&d));

    // patches that do not tile the faces
    let d = variant("start", |n, t| {
        if n == "boundary" {
            t.replacen("startFace", "startFace 3; //", 1)
        } else {
            t
        }
    });
    assert!(
        e(&d).contains("startFace") || e(&d).contains("patch"),
        "{}",
        e(&d)
    );

    // a compressed file
    let d = variant("gz", |_, t| t);
    let pm = d.join("constant/polyMesh");
    std::fs::rename(pm.join("points"), pm.join("points.gz")).unwrap();
    assert!(
        e(&d).contains("gz") || e(&d).contains("compressed"),
        "{}",
        e(&d)
    );

    let _ = openfoam::OpenFoam;
}

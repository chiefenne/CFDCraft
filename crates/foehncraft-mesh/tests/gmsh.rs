//! Stage-1 acceptance: meshes written by Gmsh load, ASCII and binary give
//! the same mesh, cells are closed and the volume is the domain volume.
//! The files are made by `tools/mesh/gmsh_test_meshes.py`.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::gmsh;
use foehncraft_mesh::quality::QualityReport;

fn load(name: &str) -> gmsh::GmshMesh {
    let p = format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"));
    gmsh::read(&p).unwrap_or_else(|e| panic!("{e}"))
}

fn check(name: &str, volume: f64) -> (gmsh::GmshMesh, QualityReport) {
    let g = load(name);
    g.mesh.validate().unwrap_or_else(|e| panic!("{name}: {e}"));
    let geo = Geometry::compute(&g.mesh);
    let (q, _) = QualityReport::compute(&g.mesh, &geo);
    assert!(q.ok(), "{name}\n{}", q.to_markdown(name));
    assert!(
        ((q.total_volume - volume) / volume).abs() < 1e-12,
        "{name}: {}",
        q.total_volume
    );
    assert!(
        ((q.boundary_volume - volume) / volume).abs() < 1e-12,
        "{name}"
    );
    eprintln!(
        "{name}: {} {:?} patches {:?} zones {:?} warnings {:?}",
        g.format,
        q.shapes,
        q.patches.iter().map(|p| (&p.0, p.1)).collect::<Vec<_>>(),
        g.mesh.zones,
        g.warnings
    );
    (g, q)
}

#[test]
fn box_tet_ascii_and_binary_are_equal() {
    let (a, qa) = check("box_tet_ascii.msh", 1.0);
    let (b, _) = check("box_tet_binary.msh", 1.0);
    assert_eq!(a.format, "MSH 4.1 ASCII");
    assert_eq!(b.format, "MSH 4.1 binary");
    // Gmsh writes 16 significant digits in ASCII, not always the exact double
    assert_eq!(a.mesh.points.len(), b.mesh.points.len());
    for (p, q) in a.mesh.points.iter().zip(&b.mesh.points) {
        assert!((0..3).all(|d| (p[d] - q[d]).abs() < 1e-14), "{p:?} {q:?}");
    }
    assert_eq!(a.mesh.face_nodes, b.mesh.face_nodes);
    assert_eq!(a.mesh.owner, b.mesh.owner);
    let names: Vec<&str> = a.mesh.patches.iter().map(|p| p.name.as_str()).collect();
    for n in ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"] {
        assert!(names.contains(&n), "{names:?}");
        let area = qa.patches.iter().find(|p| p.0 == n).unwrap().2;
        assert!((area - 1.0).abs() < 1e-12);
    }
    assert_eq!(a.mesh.zones, vec!["fluid".to_string()]);
    assert!(a.warnings.is_empty(), "{:?}", a.warnings);
}

#[test]
fn channel_with_zones_prisms_and_periodic_link() {
    let (a, q) = check("channel_mixed_ascii.msh", 1.0);
    let (b, _) = check("channel_mixed_binary.msh", 1.0);
    assert_eq!(a.mesh.face_nodes, b.mesh.face_nodes);
    assert_eq!(a.mesh.zones, vec!["left".to_string(), "right".to_string()]);
    assert!(q.shapes.iter().any(|s| s.0 == "hexahedra"));
    assert!(q.shapes.iter().any(|s| s.0 == "prisms"));
    assert_eq!(a.periodic.len(), b.periodic.len());
    assert!(a.periodic.iter().any(|p| p.dim == 2 && !p.nodes.is_empty()));
    let names: Vec<&str> = a.mesh.patches.iter().map(|p| p.name.as_str()).collect();
    assert_eq!(names, ["back", "front", "inlet", "outlet", "walls"]);
}

#[test]
fn pyramids_and_unnamed_surfaces() {
    let (g, q) = check("pyramids_ascii.msh", 2.0);
    assert!(q.shapes.iter().any(|s| s.0 == "pyramids"), "{:?}", q.shapes);
    assert!(
        g.mesh
            .patches
            .iter()
            .all(|p| p.name.starts_with("surface_"))
    );
}

#[test]
fn second_order_is_refused() {
    let txt = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 1 1 1\n3 1 0 1\n1\n0 0 0\n$EndNodes\n$Elements\n1 1 1 1\n3 1 11 1\n1 1 1 1 1 1 1 1 1 1 1\n$EndElements\n";
    let e = gmsh::parse(txt.as_bytes()).unwrap_err();
    assert!(e.contains("higher-order"), "{e}");
    let e = gmsh::parse(b"$MeshFormat\n2.2 0 8\n$EndMeshFormat\n").unwrap_err();
    assert!(e.contains("4.1"), "{e}");
}

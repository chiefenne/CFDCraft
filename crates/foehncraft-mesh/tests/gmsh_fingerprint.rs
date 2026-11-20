//! Regression guard of the mesh-format layer: the Gmsh fixtures must give
//! exactly the mesh they gave before the readers moved behind
//! `foehncraft_mesh::format` (node order, faces, owners, boundaries, zones).
//! The fingerprints were taken with the reader before the change.

use foehncraft_mesh::{Mesh, format, gmsh};

/// FNV-1a over the topology and the coordinate bits of a mesh.
fn fingerprint(m: &Mesh) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut eat = |x: u64| {
        for b in x.to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x100_0000_01b3);
        }
    };
    for p in &m.points {
        for x in p {
            eat(x.to_bits());
        }
    }
    for v in [
        &m.face_ptr,
        &m.face_nodes,
        &m.owner,
        &m.neighbour,
        &m.cell_zone,
    ] {
        eat(v.len() as u64);
        for &x in v.iter() {
            eat(x as u64);
        }
    }
    for b in &m.boundaries {
        for c in b.name.bytes() {
            eat(c as u64);
        }
        eat(b.start as u64);
        eat(b.len as u64);
    }
    for z in &m.zones {
        for c in z.bytes() {
            eat(c as u64);
        }
    }
    h
}

fn path(name: &str) -> String {
    format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))
}

const EXPECTED: [(&str, u64); 5] = [
    ("box_tet_ascii.msh", 0x93e0e8d1153954bd),
    ("box_tet_binary.msh", 0x128be9c6b4dcc0f0),
    ("channel_mixed_ascii.msh", 0xbd0c76de9d74014f),
    ("channel_mixed_binary.msh", 0x23578e665f19692a),
    ("pyramids_ascii.msh", 0x1aeec8ac6bd183e9),
];

#[test]
fn gmsh_meshes_are_unchanged() {
    let mut bad = Vec::new();
    for (name, want) in EXPECTED {
        let g = gmsh::read(path(name)).unwrap();
        let got = fingerprint(&g.mesh);
        eprintln!("(\"{name}\", 0x{got:016x}), warnings {:?}", g.warnings);
        if got != want {
            bad.push(name);
        }
        // the same mesh through the format registry (detected by extension)
        let f = format::read(path(name), None).unwrap();
        assert_eq!(f.reader, "gmsh");
        assert_eq!(f.format, g.format);
        assert_eq!(fingerprint(&f.mesh), got, "{name} via the registry");
        assert_eq!(f.warnings, g.warnings);
        let want_warnings: &[&str] = if name == "pyramids_ascii.msh" {
            &["9 surface elements lie inside the domain (internal surfaces) and are ignored"]
        } else {
            &[]
        };
        assert_eq!(g.warnings, want_warnings, "{name}");
    }
    assert!(bad.is_empty(), "changed: {bad:?}");
}

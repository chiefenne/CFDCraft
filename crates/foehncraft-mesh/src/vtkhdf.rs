//! VTKHDF output (D41).
//!
//! One HDF5 file in the VTKHDF format of VTK, read directly by ParaView
//! (open the `.vtkhdf` file; reader "VTKHDF"). The file is a
//! `PartitionedDataSetCollection`:
//!
//! - block `internal`: the cells as an `UnstructuredGrid`; tetrahedra,
//!   hexahedra, prisms and pyramids as VTK's linear cells, everything else
//!   as `VTK_POLYHEDRON` with its faces; cell fields under `CellData`;
//! - one `PolyData` block per boundary patch, named after the patch, with
//!   the patch index as a cell field;
//! - an assembly `internal`, `boundary/<patch>` for ParaView's block
//!   selection.
//!
//! Cells, nodes and boundary faces are written in the order of their global
//! identities (input order, R2), so the file does not depend on a
//! renumbering of the cells.
//!
//! The format is described in the VTK documentation, "VTK File Formats",
//! section "VTKHDF File Format". This writer is serial; the parallel writer
//! (all ranks into one file through MPI-IO) follows with the solver stages.

use crate::mesh::{CellShape, Mesh};
use anyhow::{Context, Result};
use h5::types::FixedAscii;
use hdf5_metno as h5;
use std::path::Path;

/// VTKHDF format version written.
pub const VERSION: [i64; 2] = [2, 4];

const VTK_TETRA: u8 = 10;
const VTK_HEXAHEDRON: u8 = 12;
const VTK_WEDGE: u8 = 13;
const VTK_PYRAMID: u8 = 14;
const VTK_POLYHEDRON: u8 = 42;

/// A cell field to write.
pub struct CellField<'a> {
    /// Name shown in ParaView.
    pub name: &'a str,
    /// One value per local cell.
    pub values: &'a [f64],
}

/// Writes the mesh and cell fields to `path`.
pub fn write(path: impl AsRef<Path>, mesh: &Mesh, fields: &[CellField]) -> Result<()> {
    let path = path.as_ref();
    for f in fields {
        anyhow::ensure!(
            f.values.len() == mesh.n_cells,
            "field {} has {} values for {} cells",
            f.name,
            f.values.len(),
            mesh.n_cells
        );
    }
    let file = h5::File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let root = group(&file, "VTKHDF")?;
    type_attr(&root, "PartitionedDataSetCollection")?;
    version_attr(&root)?;

    let internal = group(&root, "internal")?;
    write_internal(&internal, mesh, fields)?;
    index_attr(&internal, 0)?;

    let assembly = group(&root, "Assembly")?;
    assembly.link_soft("/VTKHDF/internal", "internal")?;
    let boundary = group(&assembly, "boundary")?;
    for (i, p) in mesh.patches.iter().enumerate() {
        let name = block_name(&p.name);
        let g = group(&root, &name)?;
        write_patch(&g, mesh, p.range(), i as f64)?;
        index_attr(&g, i as i64 + 1)?;
        boundary.link_soft(&format!("/VTKHDF/{name}"), &name)?;
    }
    Ok(())
}

/// Block names must not clash with the reserved group names.
fn block_name(patch: &str) -> String {
    match patch {
        "internal" | "Assembly" => format!("patch_{patch}"),
        _ => patch.replace('/', "_"),
    }
}

/// A group that tracks and indexes the creation order of its links: VTK
/// iterates over blocks and assembly nodes in creation order.
fn group(parent: &h5::Group, name: &str) -> Result<h5::Group> {
    use h5::plist::group_create::LinkCreationOrder;
    Ok(parent
        .create_group_builder()
        .with_gcpl(|p| p.link_creation_order(LinkCreationOrder::Indexed))
        .create(name)?)
}

fn type_attr(g: &h5::Group, t: &str) -> Result<()> {
    macro_rules! fixed {
        ($n:literal) => {{
            let v = FixedAscii::<$n>::from_ascii(t.as_bytes())?;
            g.new_attr::<FixedAscii<$n>>()
                .create("Type")?
                .write_scalar(&v)?;
        }};
    }
    match t.len() {
        8 => fixed!(8),
        16 => fixed!(16),
        28 => fixed!(28),
        _ => anyhow::bail!("unexpected type name {t}"),
    }
    Ok(())
}

fn version_attr(g: &h5::Group) -> Result<()> {
    g.new_attr::<i64>()
        .shape([2])
        .create("Version")?
        .write_raw(&VERSION)?;
    Ok(())
}

fn index_attr(g: &h5::Group, i: i64) -> Result<()> {
    g.new_attr::<i64>().create("Index")?.write_scalar(&i)?;
    Ok(())
}

fn ds<T: h5::H5Type>(g: &h5::Group, name: &str, data: &[T]) -> Result<()> {
    g.new_dataset::<T>()
        .shape([data.len()])
        .create(name)?
        .write_raw(data)?;
    Ok(())
}

fn points_ds(g: &h5::Group, pts: &[[f64; 3]]) -> Result<()> {
    let flat: Vec<f64> = pts.iter().flatten().copied().collect();
    g.new_dataset::<f64>()
        .shape([pts.len(), 3])
        .create("Points")?
        .write_raw(&flat)?;
    Ok(())
}

/// The nodes in the order of their identities, and the output position of
/// each local node.
fn points_in_id_order(mesh: &Mesh) -> (Vec<[f64; 3]>, Vec<u32>) {
    let mut order: Vec<usize> = (0..mesh.n_points()).collect();
    order.sort_by_key(|&v| mesh.point_id[v]);
    let mut pos = vec![0u32; mesh.n_points()];
    for (k, &v) in order.iter().enumerate() {
        pos[v] = k as u32;
    }
    (order.iter().map(|&v| mesh.points[v]).collect(), pos)
}

/// Nodes of face `f` ordered with the normal out of cell `c`.
fn outward(mesh: &Mesh, f: usize, is_nb: bool) -> Vec<i64> {
    let n = mesh.face(f);
    if is_nb {
        let mut v: Vec<i64> = Vec::with_capacity(n.len());
        v.push(n[0] as i64);
        v.extend(n[1..].iter().rev().map(|&x| x as i64));
        v
    } else {
        n.iter().map(|&x| x as i64).collect()
    }
}

/// VTK node order of a standard cell, or `None` for a polyhedron.
fn standard_cell(mesh: &Mesh, c: usize) -> Option<(u8, Vec<i64>)> {
    let shape = mesh.cell_shape(c);
    if shape == CellShape::Polyhedron {
        return None;
    }
    let faces: Vec<Vec<i64>> = mesh
        .cell_face_iter(c)
        .map(|(f, nb)| outward(mesh, f, nb))
        .collect();
    // node adjacent to `n` along an edge, not in `base`
    let up = |n: i64, base: &[i64]| -> i64 {
        for f in &faces {
            let k = f.len();
            for i in 0..k {
                if f[i] == n {
                    for m in [f[(i + 1) % k], f[(i + k - 1) % k]] {
                        if !base.contains(&m) {
                            return m;
                        }
                    }
                }
            }
        }
        unreachable!("standard cell without a side edge")
    };
    let first_with = |len: usize| faces.iter().find(|f| f.len() == len).unwrap().clone();
    Some(match shape {
        CellShape::Tetrahedron => {
            let t = &faces[0];
            let apex = up(t[0], t);
            (VTK_TETRA, vec![t[0], t[2], t[1], apex])
        }
        CellShape::Hexahedron => {
            let q = &faces[0];
            let base = vec![q[0], q[3], q[2], q[1]];
            let mut v = base.clone();
            v.extend(base.iter().map(|&n| up(n, &base)));
            (VTK_HEXAHEDRON, v)
        }
        CellShape::Prism => {
            // base normal towards the opposite triangle, as for the hexahedron
            let t = first_with(3);
            let base = vec![t[0], t[2], t[1]];
            let mut v = base.clone();
            v.extend(base.iter().map(|&n| up(n, &base)));
            (VTK_WEDGE, v)
        }
        CellShape::Pyramid => {
            let q = first_with(4);
            let base = vec![q[0], q[3], q[2], q[1]];
            let apex = up(base[0], &base);
            let mut v = base;
            v.push(apex);
            (VTK_PYRAMID, v)
        }
        CellShape::Polyhedron => unreachable!(),
    })
}

fn write_internal(g: &h5::Group, mesh: &Mesh, fields: &[CellField]) -> Result<()> {
    type_attr(g, "UnstructuredGrid")?;
    version_attr(g)?;
    let nc = mesh.n_cells;
    let (pts, pos) = points_in_id_order(mesh);
    let mut cells: Vec<usize> = (0..nc).collect();
    cells.sort_by_key(|&c| mesh.cell_id[c]);
    let mut types = Vec::with_capacity(nc);
    let mut offsets: Vec<i64> = Vec::with_capacity(nc + 1);
    let mut conn: Vec<i64> = Vec::with_capacity(8 * nc);
    let mut poly_offsets: Vec<i64> = Vec::with_capacity(nc + 1);
    let mut poly_faces: Vec<i64> = Vec::new();
    let mut face_offsets: Vec<i64> = vec![0];
    let mut face_conn: Vec<i64> = Vec::new();
    offsets.push(0);
    poly_offsets.push(0);
    let p = |n: i64| pos[n as usize] as i64;
    for &c in &cells {
        match standard_cell(mesh, c) {
            Some((t, v)) => {
                types.push(t);
                conn.extend(v.into_iter().map(p));
            }
            None => {
                types.push(VTK_POLYHEDRON);
                conn.extend(
                    mesh.cell_nodes(c)
                        .into_iter()
                        .map(|n| pos[n as usize] as i64),
                );
                for (f, nb) in mesh.cell_face_iter(c) {
                    poly_faces.push(face_offsets.len() as i64 - 1);
                    face_conn.extend(outward(mesh, f, nb).into_iter().map(p));
                    face_offsets.push(face_conn.len() as i64);
                }
            }
        }
        offsets.push(conn.len() as i64);
        poly_offsets.push(poly_faces.len() as i64);
    }
    ds(g, "NumberOfPoints", &[mesh.n_points() as i64])?;
    ds(g, "NumberOfCells", &[nc as i64])?;
    ds(g, "NumberOfConnectivityIds", &[conn.len() as i64])?;
    points_ds(g, &pts)?;
    ds(g, "Types", &types)?;
    ds(g, "Offsets", &offsets)?;
    ds(g, "Connectivity", &conn)?;
    ds(g, "NumberOfFaces", &[face_offsets.len() as i64 - 1])?;
    ds(g, "NumberOfFaceConnectivityIds", &[face_conn.len() as i64])?;
    ds(g, "NumberOfPolyhedronToFaceIds", &[poly_faces.len() as i64])?;
    ds(g, "FaceOffsets", &face_offsets)?;
    ds(g, "FaceConnectivity", &face_conn)?;
    ds(g, "PolyhedronOffsets", &poly_offsets)?;
    ds(g, "PolyhedronToFaces", &poly_faces)?;
    let cd = group(g, "CellData")?;
    for f in fields {
        let v: Vec<f64> = cells.iter().map(|&c| f.values[c]).collect();
        ds(&cd, f.name, &v)?;
    }
    let zone: Vec<f64> = cells.iter().map(|&c| mesh.cell_zone[c] as f64).collect();
    if !fields.iter().any(|f| f.name == "zone") {
        ds(&cd, "zone", &zone)?;
    }
    group(g, "PointData")?;
    group(g, "FieldData")?;
    Ok(())
}

fn write_patch(
    g: &h5::Group,
    mesh: &Mesh,
    faces: std::ops::Range<usize>,
    patch: f64,
) -> Result<()> {
    type_attr(g, "PolyData")?;
    version_attr(g)?;
    let mut map = std::collections::HashMap::new();
    let mut pts: Vec<[f64; 3]> = Vec::new();
    let mut conn: Vec<i64> = Vec::new();
    let mut offsets: Vec<i64> = vec![0];
    let mut order: Vec<usize> = faces.clone().collect();
    order.sort_by_key(|&f| mesh.face_id[f]);
    for f in order {
        for &n in mesh.face(f) {
            let id = *map.entry(n).or_insert_with(|| {
                pts.push(mesh.points[n as usize]);
                pts.len() as i64 - 1
            });
            conn.push(id);
        }
        offsets.push(conn.len() as i64);
    }
    ds(g, "NumberOfPoints", &[pts.len() as i64])?;
    points_ds(g, &pts)?;
    for (name, n, o, c) in [
        ("Vertices", 0, vec![0i64], vec![]),
        ("Lines", 0, vec![0], vec![]),
        ("Strips", 0, vec![0], vec![]),
        ("Polygons", faces.len(), offsets, conn),
    ] {
        let sg = group(g, name)?;
        ds(&sg, "NumberOfCells", &[n as i64])?;
        ds(&sg, "NumberOfConnectivityIds", &[c.len() as i64])?;
        ds(&sg, "Offsets", &o)?;
        ds(&sg, "Connectivity", &c)?;
    }
    let cd = group(g, "CellData")?;
    ds(&cd, "patch", &vec![patch; faces.len()])?;
    group(g, "PointData")?;
    group(g, "FieldData")?;
    Ok(())
}

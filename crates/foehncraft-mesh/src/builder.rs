//! Assembling a [`Mesh`] from faces or from cells.
//!
//! - [`MeshBuilder`]: faces are given with their cells (interior) or with
//!   their boundary (boundary); the builder orders them into the layout of
//!   [`Mesh`].
//! - [`CellSoup`]: cells are given as closed sets of polygons; faces shared
//!   by two cells are found by their node sets. A quadrilateral that meets two
//!   triangles of the cell on the other side (a hexahedron next to tetrahedra
//!   or prisms) is split into those two triangles, so the cell with the
//!   quadrilateral becomes a polyhedron with one more face. Faces left
//!   unmatched are boundary faces, named by a classification function.

use crate::geometry::{add, polygon_area_centroid, sub};
use crate::mesh::{Boundary, Mesh};
use std::collections::HashMap;

const NONE: u32 = u32::MAX;

/// Builds a mesh face by face.
#[derive(Debug, Default)]
pub struct MeshBuilder {
    points: Vec<[f64; 3]>,
    n_cells: usize,
    ptr: Vec<u32>,
    nodes: Vec<u32>,
    a: Vec<u32>,
    b: Vec<u32>,
    boundary: Vec<u32>,
    boundary_names: Vec<String>,
    cell_zone: Vec<u32>,
    zones: Vec<String>,
}

impl MeshBuilder {
    /// A builder for `n_cells` cells on the given nodes.
    pub fn new(points: Vec<[f64; 3]>, n_cells: usize) -> Self {
        Self {
            points,
            n_cells,
            ptr: vec![0],
            cell_zone: vec![0; n_cells],
            zones: vec!["fluid".to_string()],
            ..Default::default()
        }
    }

    /// Declares a boundary, so that boundaries appear in the order of declaration
    /// (also when empty). Returns its index.
    pub fn declare_boundary(&mut self, name: &str) -> u32 {
        if let Some(i) = self.boundary_names.iter().position(|p| p == name) {
            return i as u32;
        }
        self.boundary_names.push(name.to_string());
        (self.boundary_names.len() - 1) as u32
    }

    /// Sets the cell zones (names and one zone index per cell).
    pub fn set_zones(&mut self, zones: Vec<String>, cell_zone: Vec<u32>) {
        assert_eq!(cell_zone.len(), self.n_cells);
        self.zones = zones;
        self.cell_zone = cell_zone;
    }

    /// Interior face between cells `from` and `to`; the normal of `nodes`
    /// (right-hand rule) points from `from` to `to`.
    pub fn add_interior(&mut self, nodes: &[u32], from: u32, to: u32) {
        debug_assert!(from != to);
        self.nodes.extend_from_slice(nodes);
        self.ptr.push(self.nodes.len() as u32);
        self.a.push(from);
        self.b.push(to);
        self.boundary.push(NONE);
    }

    /// Boundary face of `cell` in boundary `name`; the normal points out of the
    /// domain.
    pub fn add_boundary(&mut self, nodes: &[u32], cell: u32, name: &str) {
        let p = self.declare_boundary(name);
        self.nodes.extend_from_slice(nodes);
        self.ptr.push(self.nodes.len() as u32);
        self.a.push(cell);
        self.b.push(NONE);
        self.boundary.push(p);
    }

    /// Orders the faces and builds the mesh. Identities are the positions
    /// in the new mesh.
    pub fn finish(self) -> Result<Mesh, String> {
        self.finish_with_order().map(|x| x.0)
    }

    /// As [`Self::finish`], also returning the insertion index of each face
    /// of the mesh.
    pub fn finish_with_order(self) -> Result<(Mesh, Vec<u32>), String> {
        let nf = self.a.len();
        // sort key: interior faces by (owner, neighbour), then boundary
        // faces by (boundary, owner); ties keep insertion order
        let mut order: Vec<u32> = (0..nf as u32).collect();
        let key = |f: usize| -> (u32, u32, u32) {
            if self.boundary[f] == NONE {
                let (p, n) = (self.a[f].min(self.b[f]), self.a[f].max(self.b[f]));
                (0, p, n)
            } else {
                (1 + self.boundary[f], self.a[f], 0)
            }
        };
        order.sort_by_key(|&f| (key(f as usize), f));
        let mut mesh = Mesh {
            points: self.points,
            n_cells: self.n_cells,
            cell_zone: self.cell_zone,
            zones: self.zones,
            ..Default::default()
        };
        mesh.face_id = (0..nf as u64).collect();
        mesh.n_inner_cells = self.n_cells;
        mesh.face_ptr.reserve(nf + 1);
        mesh.face_ptr.push(0);
        mesh.face_nodes.reserve(self.nodes.len());
        mesh.owner.reserve(nf);
        let mut boundary_count = vec![0usize; self.boundary_names.len()];
        for &f in &order {
            let f = f as usize;
            let nodes = &self.nodes[self.ptr[f] as usize..self.ptr[f + 1] as usize];
            let (a, b) = (self.a[f], self.b[f]);
            if a as usize >= self.n_cells || (b != NONE && b as usize >= self.n_cells) {
                return Err(format!("face with cell index out of range ({a}, {b})"));
            }
            if b == NONE {
                boundary_count[self.boundary[f] as usize] += 1;
                mesh.face_nodes.extend_from_slice(nodes);
                mesh.owner.push(a);
            } else if a < b {
                mesh.face_nodes.extend_from_slice(nodes);
                mesh.owner.push(a);
                mesh.neighbour.push(b);
            } else {
                // normal must point from the owner (smaller index): reverse
                mesh.face_nodes.push(nodes[0]);
                mesh.face_nodes.extend(nodes[1..].iter().rev());
                mesh.owner.push(b);
                mesh.neighbour.push(a);
            }
            mesh.face_ptr.push(mesh.face_nodes.len() as u32);
        }
        mesh.n_inner_faces = mesh.neighbour.len();
        let mut start = mesh.neighbour.len();
        for (name, len) in self.boundary_names.into_iter().zip(boundary_count) {
            mesh.boundaries.push(Boundary { name, start, len });
            start += len;
        }
        for w in mesh.owner[..mesh.n_inner_faces]
            .iter()
            .zip(&mesh.neighbour)
            .collect::<Vec<_>>()
            .windows(2)
        {
            if w[0] == w[1] {
                return Err(format!(
                    "cells {} and {} share more than one face",
                    w[0].0, w[0].1
                ));
            }
        }
        mesh.build_cell_faces();
        mesh.cell_id = (0..self.n_cells as u64).collect();
        mesh.point_id = (0..mesh.points.len() as u64).collect();
        Ok((mesh, order))
    }
}

/// Cells given as closed sets of polygons.
#[derive(Debug, Default)]
pub struct CellSoup {
    points: Vec<[f64; 3]>,
    cell_ptr: Vec<u32>,
    face_ptr: Vec<u32>,
    nodes: Vec<u32>,
    zone: Vec<u32>,
    zones: Vec<String>,
}

/// What [`CellSoup::finish`] needs to know about an unmatched face to name
/// its boundary.
#[derive(Debug, Clone, Copy)]
pub struct BoundaryFace<'a> {
    /// Nodes, ordered with the normal pointing out of the domain.
    pub nodes: &'a [u32],
    /// Area-weighted centroid.
    pub centroid: [f64; 3],
    /// Unit normal, out of the domain.
    pub normal: [f64; 3],
}

impl CellSoup {
    /// An empty soup on the given nodes.
    pub fn new(points: Vec<[f64; 3]>) -> Self {
        Self {
            points,
            cell_ptr: vec![0],
            face_ptr: vec![0],
            zones: vec!["fluid".to_string()],
            ..Default::default()
        }
    }

    /// Number of cells added.
    pub fn n_cells(&self) -> usize {
        self.cell_ptr.len() - 1
    }

    /// Nodes.
    pub fn points(&self) -> &[[f64; 3]] {
        &self.points
    }

    /// Sets the zone names; cells refer to them by index in
    /// [`Self::add_cell_in_zone`].
    pub fn set_zone_names(&mut self, zones: Vec<String>) {
        self.zones = zones;
    }

    /// Adds a cell from its faces. The orientation of the faces is made
    /// consistent (outward) from the sign of the cell volume, so either
    /// orientation may be given, as long as all faces of the cell agree.
    pub fn add_cell<'a>(&mut self, faces: impl IntoIterator<Item = &'a [u32]>) -> u32 {
        self.add_cell_in_zone(faces, 0)
    }

    /// As [`Self::add_cell`], in zone `zone`.
    pub fn add_cell_in_zone<'a>(
        &mut self,
        faces: impl IntoIterator<Item = &'a [u32]>,
        zone: u32,
    ) -> u32 {
        let first_face = self.face_ptr.len() - 1;
        for f in faces {
            self.nodes.extend_from_slice(f);
            self.face_ptr.push(self.nodes.len() as u32);
        }
        let last_face = self.face_ptr.len() - 1;
        // orientation from the signed volume about the node average
        let mut centre = [0.0; 3];
        let mut cnt = 0.0;
        for f in first_face..last_face {
            for &n in self.face_nodes(f) {
                centre = add(centre, self.points[n as usize]);
                cnt += 1.0;
            }
        }
        let centre = centre.map(|x| x / cnt);
        let mut vol = 0.0;
        for f in first_face..last_face {
            let (s, c) = polygon_area_centroid(&self.points, self.face_nodes(f));
            vol += crate::geometry::dot(s, sub(c, centre));
        }
        if vol < 0.0 {
            for f in first_face..last_face {
                let (a, b) = (self.face_ptr[f] as usize, self.face_ptr[f + 1] as usize);
                self.nodes[a + 1..b].reverse();
            }
        }
        self.cell_ptr.push(last_face as u32);
        self.zone.push(zone);
        (self.cell_ptr.len() - 2) as u32
    }

    fn face_nodes(&self, f: usize) -> &[u32] {
        &self.nodes[self.face_ptr[f] as usize..self.face_ptr[f + 1] as usize]
    }

    /// Matches the faces and builds the mesh. `boundaries` lists the boundary
    /// names in the order they should appear; `classify` names the boundary of
    /// each unmatched face.
    pub fn finish(
        self,
        boundaries: &[&str],
        classify: impl Fn(BoundaryFace) -> String,
    ) -> Result<Mesh, String> {
        let nf = self.face_ptr.len() - 1;
        let nc = self.n_cells();
        let mut face_cell = vec![0u32; nf];
        for c in 0..nc {
            for f in self.cell_ptr[c]..self.cell_ptr[c + 1] {
                face_cell[f as usize] = c as u32;
            }
        }
        // sorted node lists as matching keys
        let mut keys: Vec<u32> = Vec::with_capacity(self.nodes.len());
        let mut kptr = vec![0u32];
        for f in 0..nf {
            let start = keys.len();
            keys.extend_from_slice(self.face_nodes(f));
            keys[start..].sort_unstable();
            kptr.push(keys.len() as u32);
        }
        let key = |f: usize| &keys[kptr[f] as usize..kptr[f + 1] as usize];
        let mut order: Vec<u32> = (0..nf as u32).collect();
        order.sort_unstable_by(|&x, &y| {
            let (kx, ky) = (key(x as usize), key(y as usize));
            kx.len()
                .cmp(&ky.len())
                .then_with(|| kx.cmp(ky))
                .then(x.cmp(&y))
        });
        let mut partner = vec![NONE; nf];
        let mut i = 0;
        while i < nf {
            let mut j = i + 1;
            while j < nf && key(order[j] as usize) == key(order[i] as usize) {
                j += 1;
            }
            match j - i {
                1 => {}
                2 => {
                    let (x, y) = (order[i] as usize, order[i + 1] as usize);
                    if face_cell[x] == face_cell[y] {
                        return Err(format!("cell {} has a repeated face", face_cell[x]));
                    }
                    partner[x] = y as u32;
                    partner[y] = x as u32;
                }
                k => {
                    return Err(format!(
                        "{k} cells share the face with nodes {:?}",
                        key(order[i] as usize)
                    ));
                }
            }
            i = j;
        }

        // quadrilaterals against two triangles
        let mut tris: HashMap<[u32; 3], u32> = HashMap::new();
        for f in 0..nf {
            if partner[f] == NONE && kptr[f + 1] - kptr[f] == 3 {
                let k = key(f);
                tris.insert([k[0], k[1], k[2]], f as u32);
            }
        }
        let tri_key = |a: u32, b: u32, c: u32| {
            let mut k = [a, b, c];
            k.sort_unstable();
            k
        };
        // replacement faces: (nodes oriented out of `cell`, cell, partner face)
        let mut split: Vec<(Vec<u32>, u32, u32)> = Vec::new();
        let mut removed = vec![false; nf];
        if !tris.is_empty() {
            for f in 0..nf {
                if partner[f] != NONE || self.face_nodes(f).len() != 4 {
                    continue;
                }
                let q = self.face_nodes(f);
                let (a, b, c, d) = (q[0], q[1], q[2], q[3]);
                for (t1, t2) in [([a, b, c], [a, c, d]), ([a, b, d], [b, c, d])] {
                    let k1 = tri_key(t1[0], t1[1], t1[2]);
                    let k2 = tri_key(t2[0], t2[1], t2[2]);
                    if let (Some(&g1), Some(&g2)) = (tris.get(&k1), tris.get(&k2)) {
                        if partner[g1 as usize] != NONE || partner[g2 as usize] != NONE {
                            continue;
                        }
                        let cell = face_cell[f];
                        split.push((t1.to_vec(), cell, g1));
                        split.push((t2.to_vec(), cell, g2));
                        partner[g1 as usize] = f as u32;
                        partner[g2 as usize] = f as u32;
                        removed[f] = true;
                        break;
                    }
                }
            }
        }

        let mut b = MeshBuilder::new(self.points.clone(), nc);
        b.set_zones(self.zones.clone(), self.zone.clone());
        for p in boundaries {
            b.declare_boundary(p);
        }
        for f in 0..nf {
            if removed[f] {
                continue;
            }
            let p = partner[f];
            if p == NONE {
                let nodes = self.face_nodes(f);
                let (s, centroid) = polygon_area_centroid(&self.points, nodes);
                let n = crate::geometry::norm(s);
                let normal = s.map(|x| x / n);
                let name = classify(BoundaryFace {
                    nodes,
                    centroid,
                    normal,
                });
                b.add_boundary(nodes, face_cell[f], &name);
            } else if removed[p as usize] {
                // triangle matched to a split quadrilateral: added below
            } else if (f as u32) < p {
                b.add_interior(self.face_nodes(f), face_cell[f], face_cell[p as usize]);
            }
        }
        for (nodes, cell, g) in split {
            b.add_interior(&nodes, cell, face_cell[g as usize]);
        }
        b.finish()
    }
}

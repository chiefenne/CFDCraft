//! Synthetic face-based mesh for the kernel benchmarks.
//!
//! A periodic box of nx × ny × nz hexahedra, stored in the face-based form
//! of an unstructured solver (D23): per face the owner and neighbour cell,
//! per cell the list of its faces. Every face is an interior face, so every
//! cell has six faces. The polyhedral meshes of stage 1 have more faces per
//! cell; G1 is rerun on them (docs/benchmarks/g1.md).
//!
//! Cell numbering:
//! - `natural`: lexicographic, the locality of a well-ordered mesh
//!   (for example after reverse Cuthill–McKee renumbering);
//! - `shuffled`: a random permutation, the worst case for caches and for
//!   coalesced GPU loads. Real meshes lie between the two.

use foehncraft_core::index::MAX_LOCAL_INDEX;
use foehncraft_core::rng;
use foehncraft_kernels::{pack_face, unpack_face};
use serde::{Deserialize, Serialize};

/// Cell numbering of the synthetic mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "kebab-case")]
pub enum Ordering {
    /// Lexicographic numbering.
    #[default]
    Natural,
    /// Random permutation of the cells.
    Shuffled,
}

impl Ordering {
    /// Name used in case files and reports.
    pub fn name(self) -> &'static str {
        match self {
            Self::Natural => "natural",
            Self::Shuffled => "shuffled",
        }
    }
}

/// Face-based mesh connectivity.
#[derive(Debug, Clone)]
pub struct FaceMesh {
    /// Number of cells.
    pub n_cells: usize,
    /// Owner cell of each face (the smaller cell index).
    pub owner: Vec<u32>,
    /// Neighbour cell of each face.
    pub neighbour: Vec<u32>,
    /// Start of the face list of each cell in `cell_faces` (length
    /// `n_cells + 1`).
    pub cell_ptr: Vec<u32>,
    /// Faces of each cell, packed with their orientation
    /// ([`foehncraft_kernels::pack_face`]), in ascending face order.
    pub cell_faces: Vec<u32>,
}

impl FaceMesh {
    /// Number of faces.
    pub fn n_faces(&self) -> usize {
        self.owner.len()
    }

    /// Periodic box of `nx × ny × nz` cells.
    pub fn periodic_box(nx: usize, ny: usize, nz: usize, ordering: Ordering, seed: u64) -> Self {
        assert!(
            nx >= 3 && ny >= 3 && nz >= 3,
            "at least 3 cells per direction"
        );
        let n_cells = nx * ny * nz;
        assert!(
            3 * n_cells <= MAX_LOCAL_INDEX as usize,
            "mesh too large for 32-bit indices"
        );

        let label = match ordering {
            Ordering::Natural => (0..n_cells as u32).collect::<Vec<_>>(),
            Ordering::Shuffled => shuffled_labels(n_cells, seed),
        };
        let id = |i: usize, j: usize, k: usize| label[i + nx * (j + ny * k)];

        let mut faces: Vec<(u32, u32)> = Vec::with_capacity(3 * n_cells);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let c = id(i, j, k);
                    for d in [
                        id((i + 1) % nx, j, k),
                        id(i, (j + 1) % ny, k),
                        id(i, j, (k + 1) % nz),
                    ] {
                        faces.push((c.min(d), c.max(d)));
                    }
                }
            }
        }
        // Faces ordered by owner, then neighbour, as an unstructured mesh
        // reader would produce them.
        faces.sort_unstable();
        let (owner, neighbour): (Vec<u32>, Vec<u32>) = faces.into_iter().unzip();
        let (cell_ptr, cell_faces) = cell_face_lists(n_cells, &owner, &neighbour);
        let mesh = Self {
            n_cells,
            owner,
            neighbour,
            cell_ptr,
            cell_faces,
        };
        mesh.validate().expect("generated mesh is valid");
        mesh
    }

    /// Checks every index used by the kernels (the precondition of the
    /// unchecked loads in `foehncraft-kernels`).
    pub fn validate(&self) -> Result<(), String> {
        let nf = self.n_faces();
        if self.neighbour.len() != nf {
            return Err("owner and neighbour lengths differ".into());
        }
        for f in 0..nf {
            let (p, n) = (self.owner[f] as usize, self.neighbour[f] as usize);
            if p >= self.n_cells || n >= self.n_cells || p == n {
                return Err(format!("face {f}: invalid cells {p}, {n}"));
            }
        }
        if self.cell_ptr.len() != self.n_cells + 1 || self.cell_ptr[0] != 0 {
            return Err("cell_ptr has the wrong length or start".into());
        }
        if self.cell_ptr.windows(2).any(|w| w[0] > w[1])
            || *self.cell_ptr.last().unwrap() as usize != self.cell_faces.len()
        {
            return Err("cell_ptr is not a valid offset array".into());
        }
        for c in 0..self.n_cells {
            for k in self.cell_ptr[c] as usize..self.cell_ptr[c + 1] as usize {
                let (f, is_nb) = unpack_face(self.cell_faces[k]);
                if f >= nf {
                    return Err(format!("cell {c}: face {f} out of range"));
                }
                let expected = if is_nb {
                    self.neighbour[f]
                } else {
                    self.owner[f]
                };
                if expected as usize != c {
                    return Err(format!("cell {c}: face {f} has the wrong orientation"));
                }
            }
        }
        Ok(())
    }
}

impl FaceMesh {
    /// The interior faces of a stage-1 mesh (`foehncraft_mesh::Mesh`): the
    /// kernels of G1 run on the same arrays the solver will use. Boundary
    /// faces are left out (they add no gather work).
    pub fn from_mesh(m: &foehncraft_mesh::Mesh) -> Self {
        let owner = m.owner[..m.n_inner_faces].to_vec();
        let neighbour = m.neighbour[..m.n_inner_faces].to_vec();
        let (cell_ptr, cell_faces) = cell_face_lists(m.n_cells, &owner, &neighbour);
        let mesh = Self {
            n_cells: m.n_cells,
            owner,
            neighbour,
            cell_ptr,
            cell_faces,
        };
        mesh.validate().expect("stage-1 mesh is valid");
        mesh
    }
}

/// Cell-to-face lists from the face-to-cell arrays, faces in ascending order.
fn cell_face_lists(n_cells: usize, owner: &[u32], neighbour: &[u32]) -> (Vec<u32>, Vec<u32>) {
    let mut count = vec![0u32; n_cells + 1];
    for (&p, &n) in owner.iter().zip(neighbour) {
        count[p as usize + 1] += 1;
        count[n as usize + 1] += 1;
    }
    for c in 0..n_cells {
        count[c + 1] += count[c];
    }
    let cell_ptr = count.clone();
    let mut next = count;
    let mut cell_faces = vec![0u32; 2 * owner.len()];
    for (f, (&p, &n)) in owner.iter().zip(neighbour).enumerate() {
        cell_faces[next[p as usize] as usize] = pack_face(f as u32, false);
        next[p as usize] += 1;
        cell_faces[next[n as usize] as usize] = pack_face(f as u32, true);
        next[n as usize] += 1;
    }
    (cell_ptr, cell_faces)
}

/// Random permutation of `0..n` (Fisher–Yates with counter-based numbers).
fn shuffled_labels(n: usize, seed: u64) -> Vec<u32> {
    let mut p: Vec<u32> = (0..n as u32).collect();
    for i in (1..n).rev() {
        let j = (rng::keyed([foehncraft_core::rng::stream::MESH, seed, i as u64, 0x5eed])
            % (i as u64 + 1)) as usize;
        p.swap(i, j);
    }
    p
}

/// Box dimensions with about `cells` cells, as close to a cube as possible.
pub fn box_dims(cells: usize) -> (usize, usize, usize) {
    let n = (cells as f64).cbrt().round().max(3.0) as usize;
    let nz = (cells as f64 / (n * n) as f64).round().max(3.0) as usize;
    (n, n, nz)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn natural_and_shuffled_are_valid() {
        for ord in [Ordering::Natural, Ordering::Shuffled] {
            let m = FaceMesh::periodic_box(4, 5, 3, ord, 1);
            assert_eq!(m.n_cells, 60);
            assert_eq!(m.n_faces(), 180);
            assert!(m.cell_ptr.windows(2).all(|w| w[1] - w[0] == 6));
            m.validate().unwrap();
        }
    }

    #[test]
    fn box_dims_close_to_request() {
        let (a, b, c) = box_dims(1_000_000);
        assert_eq!((a, b, c), (100, 100, 100));
        let (a, b, c) = box_dims(10_000_000);
        let n = a * b * c;
        assert!((n as f64 / 1e7 - 1.0).abs() < 0.02);
    }
}

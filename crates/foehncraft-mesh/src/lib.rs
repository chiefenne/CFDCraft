//! The mesh of Foehncraft (stage 1).
//!
//! - [`mesh`]: the polyhedral face-based mesh (D23) in the layout the
//!   kernels of `foehncraft-kernels` use;
//! - [`builder`]: assembling a mesh from faces or from cells;
//! - [`geometry`]: face area vectors and centroids, cell volumes and
//!   centroids;
//! - [`quality`]: the mesh quality report;
//! - [`generate`]: the built-in generator for canonical test meshes (D34);
//! - [`cylinder`]: the built-in channel with a cylinder (Schäfer & Turek);
//! - [`wedge`]: axisymmetric wedge meshes revolved from a 2D mesh;
//! - [`gmsh`]: the Gmsh `.msh` 4.1 reader (D33);
//! - [`renumber`]: reverse Cuthill–McKee renumbering of the cells;
//! - [`locate`]: point location and cell-to-cell particle tracking (D47);
//! - [`partition`]: domain decomposition (D32);
//! - [`decompose`]: local meshes with halo cells and the halo exchange plan;
//! - [`vtkhdf`]: VTKHDF output (D41), with the feature `hdf5`.

pub mod builder;
pub mod case;
pub mod cylinder;
pub mod decompose;
pub mod generate;
pub mod geometry;
pub mod gmsh;
pub mod locate;
pub mod mesh;
pub mod partition;
pub mod quality;
pub mod renumber;
pub mod wedge;

pub use geometry::Geometry;
pub use mesh::{CellShape, Mesh, Patch};
#[cfg(feature = "hdf5")]
pub mod vtkhdf;

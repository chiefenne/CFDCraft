//! The geometry-query layer (D108): searches on surfaces and cells written
//! once, deterministic (fixed traversal rules, ties broken by global id,
//! D93).
//!
//! Part 1 (stage 6): a bounding volume hierarchy over the wall faces and the
//! closest-point query, for the exact wall distance of D48
//! ([`wall_distance`]), and the same distance for a mesh with a rigidly
//! turning part at any angle without a rebuild ([`TurningWallDistance`],
//! the sliding mesh of D182). Ray casting, the cell-box hierarchy and winding
//! numbers follow with their users (D108 parts 2 and 3).

pub mod bvh;
pub mod turning;

pub use bvh::{Bvh, Closest, SurfaceFaces, wall_distance};
pub use turning::TurningWallDistance;

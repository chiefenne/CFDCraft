//! The geometry-query layer (D108): searches on surfaces and cells written
//! once, deterministic (fixed traversal rules, ties broken by global id,
//! D93).
//!
//! Part 1 (stage 6): a bounding volume hierarchy over the wall faces and the
//! closest-point query, for the exact wall distance of D48
//! ([`wall_distance`]). Ray casting, the cell-box hierarchy and winding
//! numbers follow with their users (D108 parts 2 and 3).

pub mod bvh;

pub use bvh::{Bvh, Closest, SurfaceFaces, wall_distance};

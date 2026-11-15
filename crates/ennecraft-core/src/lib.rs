//! Core types shared by every Ennecraft crate.
//!
//! The crate is `no_std`, so that the same definitions can be compiled for
//! the CPU and for GPU targets (see `gpu/ptx-kernels`).
//!
//! - [`Real`]: the floating-point type is a generic parameter of the code
//!   (requirement D7). The reference path is `f64`.
//! - [`index`]: local (per rank) and global index types.
//! - [`repro`]: reductions with a fixed summation order (requirement D31).
//! - [`rng`]: counter-based random numbers keyed by identity, step and event,
//!   independent of thread or rank order (requirement D47).
#![no_std]

pub mod index;
pub mod real;
pub mod repro;
pub mod rng;

pub use index::{GlobalIndex, LocalIndex};
pub use real::Real;

/// Version of the Ennecraft crates, from the workspace `Cargo.toml`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

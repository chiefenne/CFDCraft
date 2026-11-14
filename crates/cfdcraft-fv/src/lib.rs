//! Finite-volume operators on a scalar transport equation (stage 2).
//!
//! - [`problem`]: the steady convection–diffusion problem, boundary
//!   conditions and manufactured solutions;
//! - [`disc`]: the host-side discretisation (geometry factors,
//!   least-squares weights, sliced-ELL structure);
//! - [`backend`]: where the kernel bodies of `cfdcraft_kernels` run: CPU
//!   threads or the NVIDIA GPU (feature `cuda`);
//! - [`solver`]: outer deferred-correction iterations and BiCGStab;
//! - [`verify`]: gradient exactness, limiter bounds and order-of-accuracy
//!   sweeps (the acceptance checks of stage 2);
//! - [`case`]: the `cfdcraft scalar` case file and run.

pub mod backend;
pub mod case;
pub mod disc;
pub mod problem;
pub mod solver;
pub mod verify;

pub use disc::Discretisation;
pub use problem::Problem;

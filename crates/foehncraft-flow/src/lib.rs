//! The incompressible flow solver of Foehncraft: steady laminar flow with the
//! SIMPLEC algorithm (stage 4) and transient flow with a PISO-type
//! algorithm (stage 5), with the Rhie–Chow interpolation on polyhedral
//! meshes, on the CPU or an NVIDIA GPU.
//!
//! - [`bc`]: boundary conditions (D29) and the check of their combination
//!   with guidance (D30);
//! - [`simplec`]: the outer iterations (D19, D27);
//! - [`analytic`]: exact solutions and reference data (Couette, Poiseuille,
//!   Ghia et al.);
//! - [`case`]: the `foehn flow` case file with validation (D39), the run
//!   with the log and CSV output (D43) and the verification checks;
//! - [`device`]: the devices of the flow loop (CPU, CUDA);
//! - [`transient`]: the time step (D28): BDF1, BDF2, Crank–Nicolson,
//!   PISO correctors, outer iterations;
//! - [`unsteady`]: the transient case driver with the adaptive time step;
//! - [`restart`]: restart files (D42);
//! - [`control`]: control files and re-reading the case file (D40);
//! - [`monitor`]: probes, force coefficients, histories, benchmarks (D43);
//! - [`plot`]: small SVG charts of the reports.

pub mod analytic;
pub mod bc;
pub mod case;
pub mod control;
pub mod device;
pub mod monitor;
pub mod plot;
pub mod restart;
pub mod simplec;
pub mod transient;
pub mod unsteady;

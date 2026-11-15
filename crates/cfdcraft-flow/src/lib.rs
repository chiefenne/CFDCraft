//! The incompressible flow solver of CFDCraft (stage 4): steady laminar
//! flow with the SIMPLEC algorithm and the Rhie–Chow interpolation on
//! polyhedral meshes.
//!
//! - [`bc`]: boundary conditions (D29) and the check of their combination
//!   with guidance (D30);
//! - [`simplec`]: the outer iterations (D19, D27);
//! - [`analytic`]: exact solutions and reference data (Couette, Poiseuille,
//!   Ghia et al.);
//! - [`case`]: the `cfdcraft flow` case file with validation (D39), the run
//!   with the log and CSV output (D43) and the verification checks;
//! - [`plot`]: small SVG charts of the reports.

pub mod analytic;
pub mod bc;
pub mod case;
pub mod plot;
pub mod simplec;

//! Gate G1 of the stage plan (docs/requirements.md, section 7): the kernel
//! benchmark that decides the GPU backend policy (D8) and confirms or
//! reopens the language decision (D11).
//!
//! Four kernels that dominate a cell-centred finite-volume solver — the two
//! passes of a face-based flux loop and the sparse matrix-vector product in
//! CSR and SELL-C-σ form — run on a synthetic face-based mesh of 10–50
//! million cells. The kernel bodies come from `inncraft-kernels`, one source
//! for the CPU and the GPU. They are timed against hand-written references:
//! CUDA C on the GPU (`kernels/g1_reference.cu`) and idiomatic Rust loops on
//! the CPU. The verdict: the single-source route passes if it is within 10 %
//! of the reference for every kernel.
//!
//! The procedure and the results are in `docs/benchmarks/g1.md`.

pub mod backend;
pub mod case;
pub mod matrix;
pub mod mesh;
pub mod report;
pub mod run;
pub mod sysinfo;
pub mod workload;

pub use case::G1Case;
pub use run::run_g1;

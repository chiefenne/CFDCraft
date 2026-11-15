//! Backends: who runs the kernel bodies, and where.

pub mod cpu;
#[cfg(feature = "cuda")]
pub mod cuda;

use crate::workload::{Kernel, Workload};

/// A way to run the G1 kernels.
pub trait Backend {
    /// Name used in case files and reports.
    fn name(&self) -> &'static str;
    /// Device description (CPU model and threads, GPU name).
    fn device(&self) -> String;
    /// Copies the workload to the device.
    fn setup(&mut self, w: &Workload) -> anyhow::Result<()>;
    /// Runs one kernel once, to completion, and returns its time in seconds.
    /// `w` is the workload given to [`Backend::setup`].
    fn run(&mut self, w: &Workload, k: Kernel) -> anyhow::Result<f64>;
    /// Output of the last run of `k`, on the host, in the original row
    /// order.
    fn output(&mut self, w: &Workload, k: Kernel) -> anyhow::Result<Vec<f64>>;
}

/// Names of the backends compiled into this build.
pub fn available() -> Vec<&'static str> {
    let mut v = vec!["cpu-serial", "cpu-threads", "cpu-native"];
    if cfg!(feature = "cuda") {
        v.extend(["cuda-reference", "cuda-rust"]);
    }
    v
}

/// Creates a backend by name.
pub fn create(name: &str, opts: &BackendOptions) -> anyhow::Result<Box<dyn Backend>> {
    match name {
        "cpu-serial" => Ok(Box::new(cpu::CpuBackend::new(
            cpu::Mode::Serial,
            opts.threads,
        )?)),
        "cpu-threads" => Ok(Box::new(cpu::CpuBackend::new(
            cpu::Mode::Threads,
            opts.threads,
        )?)),
        "cpu-native" => Ok(Box::new(cpu::CpuBackend::new(
            cpu::Mode::Native,
            opts.threads,
        )?)),
        #[cfg(feature = "cuda")]
        "cuda-reference" => Ok(Box::new(cuda::CudaBackend::reference(opts.device)?)),
        #[cfg(feature = "cuda")]
        "cuda-rust" => Ok(Box::new(cuda::CudaBackend::rust(
            opts.device,
            &opts.rust_ptx,
        )?)),
        #[cfg(not(feature = "cuda"))]
        "cuda-reference" | "cuda-rust" => anyhow::bail!(
            "backend `{name}` needs a build with the `cuda` feature: \
             cargo build --release -p foehncraft-cli --features cuda"
        ),
        _ => anyhow::bail!(
            "unknown backend `{name}`; available: {}",
            available().join(", ")
        ),
    }
}

/// Options that some backends need.
#[derive(Debug, Clone, Default)]
pub struct BackendOptions {
    /// Number of CPU threads (`None`: all cores).
    pub threads: Option<usize>,
    /// CUDA device ordinal.
    pub device: usize,
    /// Path of the PTX module built from `gpu/ptx-kernels`.
    pub rust_ptx: std::path::PathBuf,
}

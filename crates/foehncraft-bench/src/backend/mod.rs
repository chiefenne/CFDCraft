//! Backends: who runs the kernel bodies, and where. The GPU backends live
//! in `foehncraft-cuda` (D60) and register with [`register_gpu`].

pub mod cpu;

use crate::workload::{Kernel, Workload};
use std::sync::OnceLock;

/// Builds a GPU backend by name with the options.
pub type GpuFactory = fn(&str, &BackendOptions) -> anyhow::Result<Box<dyn Backend>>;

static GPU: OnceLock<(&'static [&'static str], GpuFactory)> = OnceLock::new();

/// Registers the GPU backends of this build (names and factory); called once
/// by the front end (`foehncraft_cuda::register`).
pub fn register_gpu(names: &'static [&'static str], f: GpuFactory) {
    let _ = GPU.set((names, f));
}

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
    if let Some((names, _)) = GPU.get() {
        v.extend(names.iter().copied());
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
        other => match GPU.get() {
            Some((names, f)) if names.contains(&other) => f(other, opts),
            _ if matches!(other, "cuda-reference" | "cuda-rust") => anyhow::bail!(
                "backend `{other}` needs a build with the `cuda` feature: \
                 cargo build --release -p foehncraft-cli --features cuda"
            ),
            _ => anyhow::bail!(
                "unknown backend `{other}`; available: {}",
                available().join(", ")
            ),
        },
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

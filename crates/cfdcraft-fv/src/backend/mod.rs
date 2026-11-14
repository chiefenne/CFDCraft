//! Backends: where the kernel bodies run.
//!
//! A backend owns the device-resident copy of the [`Discretisation`] and a
//! fixed set of named vectors ([`V`]); the solver drives it through the
//! operations of [`Backend`] and never touches the data directly. The CPU
//! backend runs the bodies of `cfdcraft_kernels` in (multi-threaded) loops,
//! the CUDA backend (feature `cuda`) launches the same bodies compiled to
//! PTX (`gpu/ptx-kernels`). Every operation writes each output from one
//! index only, and dot products use fixed blocks, so the results do not
//! depend on the number of threads, and CPU and GPU agree (D31).

pub mod cpu;
#[cfg(feature = "cuda")]
pub mod cuda;

use crate::disc::Discretisation;

/// The named vectors of a backend.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(missing_docs)]
pub enum V {
    /// Cell values φ (also the solution of the linear solve).
    Phi,
    /// Gradient components.
    Gx,
    Gy,
    Gz,
    /// Limiter value per cell.
    Psi,
    /// Face coefficients c_P, c_N, b_f.
    Cp,
    Cn,
    Bf,
    /// Right-hand side and inverse diagonal.
    Rhs,
    Dinv,
    /// BiCGStab work vectors.
    R,
    Rh,
    P,
    Ph,
    Vv,
    S,
    Sh,
    T,
}

impl V {
    /// Vectors with one value per cell.
    pub const CELL: [V; 15] = [
        V::Phi,
        V::Gx,
        V::Gy,
        V::Gz,
        V::Psi,
        V::Rhs,
        V::Dinv,
        V::R,
        V::Rh,
        V::P,
        V::Ph,
        V::Vv,
        V::S,
        V::Sh,
        V::T,
    ];
    /// Vectors with one value per face.
    pub const FACE: [V; 3] = [V::Cp, V::Cn, V::Bf];

    /// Index into a backend's vector table.
    pub fn index(self) -> usize {
        self as usize
    }

    /// Number of vectors.
    pub const COUNT: usize = 18;
}

/// Gradient limiter (D25).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Limiter {
    /// No limiter (ψ = 1).
    None,
    /// Venkatakrishnan (1995), ε² = (K h)³.
    Venkatakrishnan {
        /// K.
        k: f64,
    },
    /// Venkatakrishnan–Wang (Wang 2000), ε² = (K (φ_max − φ_min))².
    VenkatakrishnanWang {
        /// K.
        k: f64,
    },
}

/// Settings of the face coefficients.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaceSettings {
    /// 0 first-order upwind, 1 linear upwind.
    pub blend: f64,
    /// Apply the non-orthogonal correction.
    pub nonorth: bool,
    /// Limit of the non-orthogonal correction (negative: none).
    pub nonorth_limit: f64,
}

/// Operations of a backend.
pub trait Backend {
    /// Name (`cpu`, `cuda`).
    fn name(&self) -> String;
    /// Device description.
    fn device(&self) -> String;
    /// Number of cells.
    fn n_cells(&self) -> usize;
    /// Copies host values into a vector.
    fn upload(&mut self, v: V, data: &[f64]) -> anyhow::Result<()>;
    /// Copies a vector to the host.
    fn download(&mut self, v: V) -> anyhow::Result<Vec<f64>>;
    /// Waits for all queued work (for phase timers).
    fn sync(&mut self) -> anyhow::Result<()>;

    /// Least-squares gradient of `Phi` into `Gx`, `Gy`, `Gz`.
    fn gradient(&mut self) -> anyhow::Result<()>;
    /// Limiter of `Phi` and its gradient into `Psi`.
    fn limiter(&mut self, l: Limiter) -> anyhow::Result<()>;
    /// Face coefficients into `Cp`, `Cn`, `Bf`.
    fn faces(&mut self, s: FaceSettings) -> anyhow::Result<()>;
    /// Matrix rows, `Rhs` and `Dinv`, with under-relaxation `alpha`.
    fn assemble(&mut self, alpha: f64) -> anyhow::Result<()>;

    /// out = A x.
    fn spmv(&mut self, x: V, out: V) -> anyhow::Result<()>;
    /// out = b − A x.
    fn residual(&mut self, x: V, b: V, out: V) -> anyhow::Result<()>;
    /// (a, b) with the fixed block order.
    fn dot(&mut self, a: V, b: V) -> anyhow::Result<f64>;
    /// Minimum and maximum of a cell vector.
    fn min_max(&mut self, a: V) -> anyhow::Result<(f64, f64)>;
    /// P = R + β (P − ω Vv).
    fn update_p(&mut self, beta: f64, omega: f64) -> anyhow::Result<()>;
    /// out = x − α z.
    fn sub_scaled(&mut self, x: V, z: V, alpha: f64, out: V) -> anyhow::Result<()>;
    /// Phi += α Ph + ω Sh.
    fn update_x(&mut self, alpha: f64, omega: f64) -> anyhow::Result<()>;
    /// out = Dinv x.
    fn precondition(&mut self, x: V, out: V) -> anyhow::Result<()>;
    /// dst = src.
    fn copy(&mut self, src: V, dst: V) -> anyhow::Result<()>;
}

/// Builds a backend by name: `cpu` (all threads), `cpu-serial`, or `cuda`
/// (feature `cuda`; `ptx` is the module built by `gpu/ptx-kernels/build.sh`).
pub fn create(
    name: &str,
    disc: &Discretisation,
    threads: Option<usize>,
    ptx: &std::path::Path,
) -> anyhow::Result<Box<dyn Backend>> {
    let _ = ptx;
    match name {
        "cpu" => Ok(Box::new(cpu::CpuBackend::new(disc.clone(), threads)?)),
        "cpu-serial" => Ok(Box::new(cpu::CpuBackend::new(disc.clone(), Some(1))?)),
        #[cfg(feature = "cuda")]
        "cuda" => Ok(Box::new(cuda::CudaBackend::new(disc, 0, ptx)?)),
        #[cfg(not(feature = "cuda"))]
        "cuda" => anyhow::bail!("this build has no GPU backend; build with `--features cuda`"),
        other => anyhow::bail!("unknown backend `{other}` (cpu, cpu-serial, cuda)"),
    }
}

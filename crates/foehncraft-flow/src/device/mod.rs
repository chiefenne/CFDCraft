//! Where the loops of the flow solver run.
//!
//! The algorithms ([`crate::simplec`], [`crate::transient`]) are written once,
//! as a sequence of operations on device buffers ([`Buf`]); a
//! [`FlowDevice`] carries them out. The CPU device runs the single-source
//! bodies of `foehncraft_kernels::flow` in multi-threaded loops; a GPU device
//! (next step of stage 5) launches the same bodies compiled to PTX
//! (`gpu/ptx-kernels`) on device-resident data behind the same trait. Every
//! operation writes each output from one index only, and the reductions use
//! the fixed blocks of `foehncraft_core::repro`, so the result does not depend on
//! the number of threads, and CPU and GPU agree bitwise (D31).
//!
//! A device holds the static data of the run (mesh, geometry, least-squares
//! weights, boundary-condition arrays, matrix structure) and the linear
//! solvers (a backend of `foehncraft-fv`); the field buffers are allocated by
//! the algorithm and addressed by handle.

pub mod cpu;
pub mod kernel;

use crate::bc::BoundaryArrays;
use anyhow::Result;
use foehncraft_fv::amg::AmgSettings;
use foehncraft_fv::disc::{Discretisation, Soa3};
use foehncraft_fv::solver::Preconditioner;

/// Handle of a device buffer of f64 values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Buf(pub(crate) u32);

impl Buf {
    /// The handle of the `i`-th buffer of an executor (executors number
    /// their buffers from 0 in the order of allocation).
    pub fn new(i: usize) -> Self {
        Self(i as u32)
    }

    /// Index of the buffer in its executor.
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Three buffers, one per component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Buf3 {
    /// x components.
    pub x: Buf,
    /// y components.
    pub y: Buf,
    /// z components.
    pub z: Buf,
}

impl Buf3 {
    /// The three handles in order.
    pub fn all(self) -> [Buf; 3] {
        [self.x, self.y, self.z]
    }
}

/// Which least-squares weights a gradient uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Lsq {
    /// Weights of the velocity (outlet and symmetry stencil points
    /// projected).
    Velocity,
    /// Weights of the pressure (wall and symmetry points projected).
    Pressure,
}

/// ε² of the limiter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Eps2 {
    /// Per cell, Venkatakrishnan (stored on the device at creation).
    PerCell,
    /// One value, Venkatakrishnan–Wang.
    Single(f64),
}

/// Inputs and outputs of [`FlowDevice::boundary_values`].
#[derive(Debug, Clone, Copy)]
pub struct BoundaryArgs {
    /// Cell velocity.
    pub u: Buf3,
    /// Velocity gradients.
    pub grad_u: [Buf3; 3],
    /// Cell pressure.
    pub p: Buf,
    /// Pressure gradient.
    pub grad_p: Buf3,
    /// Face mass flux.
    pub mass_flux: Buf,
    /// Prescribed boundary velocity (faces with the fixed rule).
    pub ub_set: Buf3,
    /// New boundary velocity (output).
    pub ub: Buf3,
    /// New boundary pressure (output).
    pub pb: Buf,
}

/// Inputs of [`FlowDevice::momentum_faces`].
#[derive(Debug, Clone, Copy)]
pub struct MomentumFaceArgs {
    /// Convecting face mass flux.
    pub mass_flux: Buf,
    /// Cell velocity.
    pub u: Buf3,
    /// Velocity of the central correction (`u`, or extrapolated in time).
    pub uc: Buf3,
    /// Velocity gradients.
    pub grad_u: [Buf3; 3],
    /// Limiter per component.
    pub psi: Buf3,
    /// Boundary velocity.
    pub ub: Buf3,
    /// 1 linear upwind, 0 upwind.
    pub blend: f64,
    /// Fraction of the bounded central scheme.
    pub central: f64,
    /// β of its bounding.
    pub beta: f64,
    /// 1 with the non-orthogonal correction.
    pub nonorth: f64,
    /// 1 with convection.
    pub convection: f64,
    /// Output c_P.
    pub cp: Buf,
    /// Output c_N.
    pub cn: Buf,
    /// Output b_f per component.
    pub bf: Buf3,
}

/// Time terms of a momentum row.
#[derive(Debug, Clone, Copy)]
pub struct TimeTerms {
    /// ρ a₀ / Δt.
    pub tau: f64,
    /// ρ (−a₁) / Δt and ρ (−a₂) / Δt.
    pub c_old: [f64; 2],
    /// uⁿ.
    pub u_old: Buf3,
    /// uⁿ⁻¹ (used when `c_old[1]` is not zero).
    pub u_old2: Buf3,
    /// 1 − θ and the old spatial operator (Crank–Nicolson); weight 0: none.
    pub explicit: (f64, Buf3),
}

/// Inputs and outputs of [`FlowDevice::momentum_rows`].
#[derive(Debug, Clone, Copy)]
pub struct MomentumRowArgs {
    /// Face coefficients.
    pub cp: Buf,
    /// See `cp`.
    pub cn: Buf,
    /// Deferred parts.
    pub bf: Buf3,
    /// Pressure gradient.
    pub grad_p: Buf3,
    /// Body force per unit volume.
    pub body_force: [f64; 3],
    /// Current velocity (relaxation).
    pub u: Buf3,
    /// α_u.
    pub alpha: f64,
    /// θ.
    pub theta: f64,
    /// Time terms (None: steady).
    pub time: Option<TimeTerms>,
    /// Output: matrix values.
    pub val: Buf,
    /// Output: right-hand sides.
    pub rhs: Buf3,
    /// Output: D.
    pub d: Buf,
    /// Output: D^s.
    pub ds: Buf,
    /// Output: D^C.
    pub dc: Buf,
    /// Output: relaxed diagonal per component.
    pub diag: Buf3,
}

/// Inputs of [`FlowDevice::mass_flux`].
#[derive(Debug, Clone, Copy)]
pub struct FluxArgs {
    /// Cell velocity.
    pub u: Buf3,
    /// Velocity gradients.
    pub grad_u: [Buf3; 3],
    /// Cell pressure.
    pub p: Buf,
    /// Pressure gradient.
    pub grad_p: Buf3,
    /// Boundary pressure.
    pub pb: Buf,
    /// Boundary velocity.
    pub ub: Buf3,
    /// D per cell; None: the plain interpolated flux (no pressure terms).
    pub d: Option<Buf>,
    /// Relaxation term per face.
    pub relax: Buf,
    /// Time-step term: (D^s, Y, α_u).
    pub time: Option<(Buf, Buf, f64)>,
}

/// Static data of a device.
pub struct Setup<'a> {
    /// Discretisation (geometry, velocity weights, matrix structure).
    pub disc: &'a Discretisation,
    /// Least-squares weights of the pressure.
    pub lsq_p: &'a Soa3,
    /// Boundary arrays.
    pub bc: &'a BoundaryArrays,
    /// Density.
    pub rho: f64,
    /// Dynamic viscosity (constant in laminar flow).
    pub mu: f64,
    /// ε² of the Venkatakrishnan limiter per cell (empty: unused).
    pub eps2: &'a [f64],
    /// CPU threads.
    pub threads: Option<usize>,
}

/// The operations of the flow algorithms. Buffers are created by
/// [`alloc`](FlowDevice::alloc); an output buffer is never one of the
/// inputs of the same operation unless the operation says so.
pub trait FlowDevice {
    /// Backend name (`cpu`, `cuda`).
    fn name(&self) -> String;
    /// Device description.
    fn device(&self) -> String;

    /// A new buffer of `n` zeros.
    fn alloc(&mut self, n: usize) -> Result<Buf>;
    /// Copies host values into a buffer.
    fn upload(&mut self, b: Buf, v: &[f64]) -> Result<()>;
    /// Copies a buffer to the host.
    fn download(&mut self, b: Buf) -> Result<Vec<f64>>;
    /// dst = src.
    fn copy(&mut self, src: Buf, dst: Buf) -> Result<()>;
    /// b = v everywhere.
    fn fill(&mut self, b: Buf, v: f64) -> Result<()>;
    /// out = a (x − y); `out` may be `x` or `y`.
    fn scaled_difference(&mut self, out: Buf, a: f64, x: Buf, y: Buf) -> Result<()>;
    /// out = a x + b y; `out` may be `x` or `y`.
    fn lincomb(&mut self, out: Buf, a: f64, x: Buf, b: f64, y: Buf) -> Result<()>;
    /// Σ x² in the fixed block order.
    fn sum_sq(&mut self, x: Buf) -> Result<f64>;
    /// Σ x y in the fixed block order (x and y of equal length).
    fn dot(&mut self, x: Buf, y: Buf) -> Result<f64>;
    /// Minimum and maximum.
    fn min_max(&mut self, x: Buf) -> Result<(f64, f64)>;
    /// One value.
    fn value(&mut self, b: Buf, i: usize) -> Result<f64>;
    /// b −= s.
    fn subtract(&mut self, b: Buf, s: f64) -> Result<()>;
    /// `b[i] *= s`.
    fn scale_entry(&mut self, b: Buf, i: usize, s: f64) -> Result<()>;

    /// New boundary velocity and pressure from the current fields.
    fn boundary_values(&mut self, a: &BoundaryArgs) -> Result<()>;
    /// The per-boundary-face buffers of the convective outlets that the
    /// boundary values read ([`foehncraft_kernels::flow::UB_CONVECTIVE`]):
    /// the old-level part c (three components) and β, zero at creation.
    fn convective_buffers(&self) -> (Buf3, Buf);
    /// Backflow suppression after the predicted fluxes `m`: closes the
    /// pressure-outlet faces whose predicted flux enters the domain (their
    /// flux in `m` becomes 0, the closed state is kept on the device for
    /// the pressure correction of this iteration).
    fn close_backflow(&mut self, m: Buf) -> Result<()>;
    /// Least-squares gradient of `phi` with boundary values `value`, where
    /// `zero_gradient` (one of [`Lsq`]'s masks, by `mask`) is 1 the owner
    /// value is used instead.
    fn gradient(&mut self, w: Lsq, phi: Buf, mask: Mask, value: Buf, out: Buf3) -> Result<()>;
    /// Venkatakrishnan limiter of the velocity component `phi` with its
    /// gradient and boundary values.
    fn limiter(&mut self, phi: Buf, value: Buf, grad: Buf3, eps2: Eps2, out: Buf) -> Result<()>;
    /// Face mass fluxes (Rhie–Chow).
    fn mass_flux(&mut self, a: &FluxArgs, out: Buf) -> Result<()>;
    /// Momentum face coefficients.
    fn momentum_faces(&mut self, a: &MomentumFaceArgs) -> Result<()>;
    /// Momentum rows.
    fn momentum_rows(&mut self, a: &MomentumRowArgs) -> Result<()>;
    /// Spatial operator of the momentum equations with the coefficients
    /// `cp`, `cn`, `bf` and the velocity `u` (Crank–Nicolson).
    fn momentum_operator(&mut self, cp: Buf, cn: Buf, bf: Buf3, u: Buf3, out: Buf3) -> Result<()>;
    /// Momentum residual scale ‖a u‖ + ‖b − (a/α − a) u‖ over the three
    /// components (unrelaxed a and b).
    fn residual_scale(&mut self, u: Buf3, rhs: Buf3, diag: Buf3, alpha: f64) -> Result<f64>;
    /// Slot 0 of every row of `val` = `diag`.
    fn set_diagonal(&mut self, val: Buf, diag: Buf) -> Result<()>;
    /// Net outflow Σ ±ṁ per cell.
    fn net_outflow(&mut self, m: Buf, out: Buf) -> Result<()>;
    /// Throughflow Σ |ṁ| per cell.
    fn throughflow(&mut self, m: Buf, out: Buf) -> Result<()>;
    /// Courant number per cell, Δt Σ|ṁ| / (2 ρ V).
    fn courant(&mut self, m: Buf, dt: f64, out: Buf) -> Result<()>;
    /// Pressure-correction face coefficients from D^C.
    fn pressure_faces(&mut self, dc: Buf, out: Buf) -> Result<()>;
    /// Pressure-correction rows.
    fn pressure_rows(&mut self, cf: Buf, m: Buf, val: Buf, rhs: Buf) -> Result<()>;
    /// Corrected fluxes ṁ − c_f (p'_N − p'_P).
    fn correct_fluxes(&mut self, cf: Buf, m: Buf, pc: Buf, out: Buf) -> Result<()>;
    /// u_i −= D^C ∇p'_i.
    fn correct_velocity(&mut self, dc: Buf, grad: Buf3, u: Buf3) -> Result<()>;
    /// out = −Σ_N val_PN x_N / diag (the neighbour correction H(x)/A).
    fn neighbour_correction(&mut self, val: Buf, x: Buf, diag: Buf, out: Buf) -> Result<()>;
    /// Mass flux of a cell velocity increment.
    fn delta_flux(&mut self, delta: Buf3, out: Buf) -> Result<()>;
    /// out = Aᵀ in the storage of A = `val` (D59: the same solvers and the
    /// multigrid then work on Aᵀ).
    fn transpose(&mut self, val: Buf, out: Buf) -> Result<()>;
    /// out = A x for the matrix `val`.
    fn matrix_product(&mut self, val: Buf, x: Buf, out: Buf) -> Result<()>;
    /// out = Aᵀ x for the matrix `val`, without storing Aᵀ (D59).
    fn transposed_product(&mut self, val: Buf, x: Buf, out: Buf) -> Result<()>;

    /// Solves `val` x = `rhs` (x holds the start value): CG with `pc`
    /// when `cg`, else BiCGStab. Returns the iterations and ‖b − A x₀‖.
    #[allow(clippy::too_many_arguments)]
    fn solve(
        &mut self,
        val: Buf,
        rhs: Buf,
        x: Buf,
        pc: Preconditioner,
        cg: bool,
        tol: f64,
        max_iter: usize,
    ) -> Result<(usize, f64)>;
    /// Multigrid setup on the matrix `val`.
    fn amg_setup(&mut self, val: Buf, s: &AmgSettings) -> Result<()>;
    /// Waits for queued work (timers).
    fn sync(&mut self) -> Result<()>;
}

/// Boundary masks of the gradients: where the owner value replaces the face
/// value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mask {
    /// Zero normal gradient of the velocity (outlets).
    Velocity,
    /// Zero normal gradient of the pressure (walls, symmetry planes).
    Pressure,
    /// The pressure correction: zero on pressure boundaries, zero normal
    /// gradient elsewhere.
    PressureCorrection,
}

/// Builds a GPU flow device from the static data and the kernel module
/// (for CUDA the PTX file of `gpu/ptx-kernels/build.sh`).
pub type GpuFactory = fn(&Setup, &std::path::Path) -> Result<Box<dyn FlowDevice>>;

static GPU: std::sync::OnceLock<(&'static str, GpuFactory)> = std::sync::OnceLock::new();

/// Registers the GPU flow device of this build under `name` (`cuda`);
/// called once by the front end (`foehncraft_cuda::register`, D60: no
/// vendor code in this crate).
pub fn register_gpu(name: &'static str, f: GpuFactory) {
    let _ = GPU.set((name, f));
}

/// The flow backends of this build: `cpu`, and the GPU backend when one is
/// registered.
pub fn available() -> Vec<&'static str> {
    let mut v = vec!["cpu"];
    if let Some((n, _)) = GPU.get() {
        v.push(n);
    }
    v
}

/// Creates a device: `cpu`, or `cuda` when the build registered it (`ptx`
/// names its kernel module, from `gpu/ptx-kernels`).
pub fn create(name: &str, s: &Setup, ptx: &std::path::Path) -> Result<Box<dyn FlowDevice>> {
    match name {
        "cpu" => Ok(Box::new(cpu::CpuFlow::cpu(s)?)),
        "cuda" => match GPU.get() {
            Some((n, f)) if *n == name => f(s, ptx),
            _ => anyhow::bail!(
                "this build has no GPU backend; build with `--features cuda` \
                 (cargo build --release -p foehncraft-cli --features cuda)"
            ),
        },
        other => anyhow::bail!("unknown backend `{other}` (cpu, cuda)"),
    }
}

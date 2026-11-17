//! Launch entries of the flow solver (stage 5): one argument block and one
//! per-index body for every operation of the flow loop, shared by the CPU
//! loops (`foehncraft-flow`, `device::kernel`) and the PTX entry points
//! (`gpu/ptx-kernels/src/flow.rs`, entry names `flow_*`).
//!
//! An argument block is a `#[repr(C)]` struct of device addresses with
//! their lengths ([`Raw`]) and scalars. The host fills it with the
//! addresses of its buffers (host memory on the CPU, device memory on the
//! GPU) and passes it by value: as a kernel parameter on the GPU (the
//! NVPTX back end passes aggregates like CUDA C, `.param .align 8 .b8`), as
//! a reference into the loop body on the CPU. Both sides compile this one
//! definition, so the argument lists cannot drift apart, and both run the
//! same body for an index, so CPU and GPU agree bitwise (D31).
//!
//! The bodies call the kernel bodies of [`crate::flow`], [`crate::fv`] and
//! [`crate::blas`] and write their outputs; each output is written by one
//! index only. Inputs and outputs of one launch are distinct buffers, except
//! where an entry says that an output may be an input (element-wise
//! updates, which read index i before they write it, through raw pointers).

use crate::adjoint;
use crate::flow::{self as k, FaceGeometry, VelocityGradient};
use crate::fv::{self, Boundary, SellLayout, Topology, Vec3};
use crate::{blas, ld};
use core::marker::PhantomData;

/// A device array: address and length. The address is that of host memory
/// for the CPU loops and of GPU memory for the PTX entries.
#[repr(C)]
pub struct Raw<T> {
    /// Address of the first entry.
    pub ptr: u64,
    /// Number of entries.
    pub len: u64,
    _t: PhantomData<fn() -> T>,
}

impl<T> Clone for Raw<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Raw<T> {}

impl<T> Default for Raw<T> {
    fn default() -> Self {
        Self::EMPTY
    }
}

impl<T> core::fmt::Debug for Raw<T> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Raw({:#x}, {})", self.ptr, self.len)
    }
}

impl<T> Raw<T> {
    /// No array.
    pub const EMPTY: Self = Self {
        ptr: 0,
        len: 0,
        _t: PhantomData,
    };

    /// An array at `ptr` with `len` entries.
    pub const fn new(ptr: u64, len: usize) -> Self {
        Self {
            ptr,
            len: len as u64,
            _t: PhantomData,
        }
    }

    /// The entries as a slice (empty for [`Raw::EMPTY`]).
    ///
    /// # Safety
    /// `ptr` addresses `len` initialised entries that no thread writes while
    /// the slice is used.
    #[inline(always)]
    pub unsafe fn get<'a>(self) -> &'a [T] {
        if self.len == 0 {
            &[]
        } else {
            // SAFETY: guaranteed by the caller.
            unsafe { core::slice::from_raw_parts(self.ptr as usize as *const T, self.len as usize) }
        }
    }

    /// The address for writing.
    #[inline(always)]
    pub fn out(self) -> *mut T {
        self.ptr as usize as *mut T
    }
}

impl<T: Copy> Raw<T> {
    /// Entry `i` read through the raw address (for arrays that the same
    /// launch also writes, at the same index).
    ///
    /// # Safety
    /// `i < len`.
    #[inline(always)]
    pub unsafe fn read(self, i: usize) -> T {
        // SAFETY: guaranteed by the caller.
        unsafe { *(self.ptr as usize as *const T).add(i) }
    }
}

/// Three arrays, one per component.
#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
pub struct Raw3 {
    /// x components.
    pub x: Raw<f64>,
    /// y components.
    pub y: Raw<f64>,
    /// z components.
    pub z: Raw<f64>,
}

impl Raw3 {
    /// As a [`Vec3`].
    ///
    /// # Safety
    /// As [`Raw::get`] for all three.
    #[inline(always)]
    pub unsafe fn v<'a>(self) -> Vec3<'a, f64> {
        // SAFETY: guaranteed by the caller.
        unsafe {
            Vec3 {
                x: self.x.get(),
                y: self.y.get(),
                z: self.z.get(),
            }
        }
    }
}

/// The nine arrays of a velocity gradient.
#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
pub struct RawGrad {
    /// ∇u_x, ∇u_y, ∇u_z.
    pub g: [Raw3; 3],
}

impl RawGrad {
    /// As a [`VelocityGradient`].
    ///
    /// # Safety
    /// As [`Raw::get`] for all nine.
    #[inline(always)]
    pub unsafe fn v<'a>(self) -> VelocityGradient<'a, f64> {
        // SAFETY: guaranteed by the caller.
        unsafe {
            VelocityGradient {
                x: self.g[0].v(),
                y: self.g[1].v(),
                z: self.g[2].v(),
            }
        }
    }
}

/// Static data of a mesh: connectivity, geometry, matrix layout and the
/// boundary-condition arrays (uploaded once per run).
#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
pub struct Mesh {
    /// Number of local cells.
    pub n_cells: u64,
    /// Number of faces.
    pub n_faces: u64,
    /// Number of interior faces (boundary faces follow).
    pub n_interior: u64,
    /// Owner cell of each face.
    pub owner: Raw<u32>,
    /// Neighbour cell of each interior face.
    pub neighbour: Raw<u32>,
    /// Start of the face list of each cell.
    pub cell_ptr: Raw<u32>,
    /// Packed faces of each cell.
    pub cell_faces: Raw<u32>,
    /// Face area vectors.
    pub area: Raw3,
    /// Face centroids.
    pub face_centroid: Raw3,
    /// Cell centroids.
    pub cell_centroid: Raw3,
    /// |S|²/(d·S).
    pub dgeo: Raw<f64>,
    /// Non-orthogonal part k of the face area vector.
    pub k: Raw3,
    /// Owner weight of the face interpolation.
    pub weight: Raw<f64>,
    /// Cell volumes.
    pub volume: Raw<f64>,
    /// Least-squares weights of the velocity (per cell-face entry).
    pub lsq_u: Raw3,
    /// Least-squares weights of the pressure.
    pub lsq_p: Raw3,
    /// SELL slice height.
    pub chunk: u64,
    /// SELL slice starts.
    pub slice_ptr: Raw<u32>,
    /// Transpose map of the SELL storage (D59; see [`crate::adjoint`]).
    pub transpose: Raw<u32>,
    /// Momentum kind per boundary face.
    pub mom_kind: Raw<u32>,
    /// Flux kind per boundary face.
    pub flux_kind: Raw<u32>,
    /// Boundary velocity rule per boundary face.
    pub ub_rule: Raw<u32>,
    /// Boundary pressure rule per boundary face.
    pub pb_rule: Raw<u32>,
    /// Prescribed (static or total) pressure per boundary face.
    pub pset: Raw<f64>,
    /// Viscosity per face.
    pub mu_face: Raw<f64>,
    /// Closed state per boundary face of the backflow suppression (1
    /// closed, 0 open), written by [`CloseBackflow`] in every outer
    /// iteration.
    pub closed: Raw<f64>,
    /// Old-level part c of the convective outlets per boundary face, set
    /// per time step ([`k::UB_CONVECTIVE`]).
    pub conv: Raw3,
    /// β of the convective outlets per boundary face, set per time step.
    pub beta: Raw<f64>,
    /// First periodic face (`n_interior` without periodic faces; D29).
    pub periodic: u64,
    /// Shift of the neighbour centroid per periodic face.
    pub shift: Raw3,
}

impl Mesh {
    /// Connectivity.
    ///
    /// # Safety
    /// As [`Raw::get`].
    #[inline(always)]
    pub unsafe fn topo<'a>(&self) -> Topology<'a> {
        // SAFETY: guaranteed by the caller.
        unsafe {
            Topology {
                owner: self.owner.get(),
                neighbour: self.neighbour.get(),
                cell_ptr: self.cell_ptr.get(),
                cell_faces: self.cell_faces.get(),
                n_interior: self.n_interior as usize,
                periodic: self.periodic as usize,
            }
        }
    }

    /// Geometry.
    ///
    /// # Safety
    /// As [`Raw::get`].
    #[inline(always)]
    pub unsafe fn geo<'a>(&self) -> FaceGeometry<'a, f64> {
        // SAFETY: guaranteed by the caller.
        unsafe {
            FaceGeometry {
                area: self.area.v(),
                face_centroid: self.face_centroid.v(),
                cell_centroid: self.cell_centroid.v(),
                dgeo: self.dgeo.get(),
                k: self.k.v(),
                weight: self.weight.get(),
                shift: self.shift.v(),
            }
        }
    }

    /// Matrix layout.
    ///
    /// # Safety
    /// As [`Raw::get`].
    #[inline(always)]
    pub unsafe fn layout<'a>(&self) -> SellLayout<'a> {
        SellLayout {
            chunk: self.chunk as usize,
            // SAFETY: guaranteed by the caller.
            slice_ptr: unsafe { self.slice_ptr.get() },
        }
    }

    /// Number of boundary faces.
    #[inline(always)]
    pub fn n_boundary(&self) -> usize {
        (self.n_faces - self.n_interior) as usize
    }
}

/// An operation of the flow loop: its argument block and its body for one
/// index. The PTX entry `flow_<NAME>` runs [`Entry::run`] for every index
/// below the count given at launch.
pub trait Entry {
    /// Argument block.
    type Args: Copy + Send + Sync;
    /// Name of the operation (PTX entry `flow_<NAME>`).
    const NAME: &'static str;
    /// The body for index `i`.
    ///
    /// # Safety
    /// `i` is below the count of the launch, every address in `a` is valid
    /// for its length and the outputs are written by no other launch at the
    /// same time.
    unsafe fn run(i: usize, a: &Self::Args);
}

// SAFETY (of the impls below): the argument blocks hold plain addresses and
// numbers; sharing them between threads is what a launch does.
macro_rules! send_sync {
    ($($t:ty),*) => {$(
        // SAFETY: see above.
        unsafe impl Send for $t {}
        // SAFETY: see above.
        unsafe impl Sync for $t {}
    )*};
}

/// Arguments of the element-wise entries.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct LinArgs {
    /// a.
    pub a: f64,
    /// b.
    pub b: f64,
    /// x.
    pub x: Raw<f64>,
    /// y.
    pub y: Raw<f64>,
    /// Output (may be `x` or `y`).
    pub out: Raw<f64>,
}

/// out = a x + b y.
pub struct Lincomb;
impl Entry for Lincomb {
    type Args = LinArgs;
    const NAME: &'static str = "lincomb";
    #[inline(always)]
    unsafe fn run(i: usize, a: &LinArgs) {
        // SAFETY: guaranteed by the caller; index i only.
        unsafe { *a.out.out().add(i) = a.a * a.x.read(i) + a.b * a.y.read(i) }
    }
}

/// out = a (x − y).
pub struct ScaledDifference;
impl Entry for ScaledDifference {
    type Args = LinArgs;
    const NAME: &'static str = "scaled_difference";
    #[inline(always)]
    unsafe fn run(i: usize, a: &LinArgs) {
        // SAFETY: guaranteed by the caller; index i only.
        unsafe { *a.out.out().add(i) = a.a * (a.x.read(i) - a.y.read(i)) }
    }
}

/// out = a (every entry).
pub struct Fill;
impl Entry for Fill {
    type Args = LinArgs;
    const NAME: &'static str = "fill";
    #[inline(always)]
    unsafe fn run(i: usize, a: &LinArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(i) = a.a }
    }
}

/// out −= a.
pub struct Subtract;
impl Entry for Subtract {
    type Args = LinArgs;
    const NAME: &'static str = "subtract";
    #[inline(always)]
    unsafe fn run(i: usize, a: &LinArgs) {
        // SAFETY: guaranteed by the caller; index i only.
        unsafe { *a.out.out().add(i) = a.out.read(i) - a.a }
    }
}

/// out\[b as index\] *= a (launched with count 1).
pub struct ScaleEntry;
impl Entry for ScaleEntry {
    type Args = LinArgs;
    const NAME: &'static str = "scale_entry";
    #[inline(always)]
    unsafe fn run(_: usize, a: &LinArgs) {
        let j = a.b as usize;
        // SAFETY: guaranteed by the caller (j inside out).
        unsafe { *a.out.out().add(j) = a.out.read(j) * a.a }
    }
}

/// Arguments of the block reductions.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct BlockArgs {
    /// The array.
    pub x: Raw<f64>,
    /// Block sums (or minima), one per block of `foehncraft_core::repro::BLOCK`.
    pub out: Raw<f64>,
    /// Block maxima (min_max only).
    pub out2: Raw<f64>,
}

/// Σ x² per block (the leaves of the fixed tree of D31).
pub struct BlockSumSq;
impl Entry for BlockSumSq {
    type Args = BlockArgs;
    const NAME: &'static str = "block_sum_sq";
    #[inline(always)]
    unsafe fn run(b: usize, a: &BlockArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let x = a.x.get();
            *a.out.out().add(b) = blas::block_dot(b, x, x);
        }
    }
}

/// Arguments of [`BlockDot`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct DotArgs {
    /// First array.
    pub x: Raw<f64>,
    /// Second array (as long as `x`).
    pub y: Raw<f64>,
    /// Block sums, one per block of `foehncraft_core::repro::BLOCK`.
    pub out: Raw<f64>,
}

/// Σ x y per block (the leaves of the fixed tree of D31).
pub struct BlockDot;
impl Entry for BlockDot {
    type Args = DotArgs;
    const NAME: &'static str = "block_dot";
    #[inline(always)]
    unsafe fn run(b: usize, a: &DotArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(b) = blas::block_dot(b, a.x.get(), a.y.get());
        }
    }
}

/// Minimum and maximum per block.
pub struct BlockMinMax;
impl Entry for BlockMinMax {
    type Args = BlockArgs;
    const NAME: &'static str = "block_min_max";
    #[inline(always)]
    unsafe fn run(b: usize, a: &BlockArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let (mn, mx) = blas::block_min_max(b, a.x.get());
            *a.out.out().add(b) = mn;
            *a.out2.out().add(b) = mx;
        }
    }
}

/// Arguments of [`BoundaryConstants`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ConstArgs {
    /// Zero-gradient mask per boundary face.
    pub mask: Raw<f64>,
    /// Boundary values.
    pub value: Raw<f64>,
    /// Output per boundary face.
    pub out: Raw<f64>,
}

/// Boundary constants of a gradient ([`k::boundary_constant`]).
pub struct BoundaryConstants;
impl Entry for BoundaryConstants {
    type Args = ConstArgs;
    const NAME: &'static str = "boundary_constants";
    #[inline(always)]
    unsafe fn run(fb: usize, a: &ConstArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(fb) = k::boundary_constant(fb, a.mask.get(), a.value.get()) }
    }
}

/// Arguments of [`BoundaryValues`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct BoundaryValuesArgs {
    /// Mesh.
    pub m: Mesh,
    /// Density.
    pub rho: f64,
    /// Face mass flux.
    pub mass_flux: Raw<f64>,
    /// Cell velocity.
    pub u: Raw3,
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// Cell pressure.
    pub p: Raw<f64>,
    /// Pressure gradient.
    pub grad_p: Raw3,
    /// Prescribed boundary velocity.
    pub ub_set: Raw3,
    /// Output: boundary velocity.
    pub ub: Raw3,
    /// Output: boundary pressure.
    pub pb: Raw<f64>,
}

/// New boundary velocity and pressure of boundary face `fb`.
pub struct BoundaryValues;
impl Entry for BoundaryValues {
    type Args = BoundaryValuesArgs;
    const NAME: &'static str = "boundary_values";
    #[inline(always)]
    unsafe fn run(fb: usize, a: &BoundaryValuesArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = k::BoundaryInputs {
                topo: a.m.topo(),
                geo: a.m.geo(),
                rho: a.rho,
                mass_flux: a.mass_flux.get(),
            };
            let (u, g, ub0) = (a.u.v(), a.grad_u.v(), a.ub_set.v());
            let rule = ld(a.m.ub_rule.get(), fb);
            let conv = if rule == k::UB_CONVECTIVE {
                (a.m.conv.v().at(fb), ld(a.m.beta.get(), fb))
            } else {
                ([0.0; 3], 0.0)
            };
            let v = k::boundary_velocity(fb, &inp, rule, ub0.at(fb), conv, &u, &g);
            *a.ub.x.out().add(fb) = v[0];
            *a.ub.y.out().add(fb) = v[1];
            *a.ub.z.out().add(fb) = v[2];
            let (p, gp) = (a.p.get(), a.grad_p.v());
            *a.pb.out().add(fb) = k::boundary_pressure(
                fb,
                &inp,
                ld(a.m.pb_rule.get(), fb),
                ld(a.m.pset.get(), fb),
                p,
                &gp,
            );
        }
    }
}

/// Arguments of [`Gradient`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct GradientArgs {
    /// Mesh.
    pub m: Mesh,
    /// Least-squares weights.
    pub w: Raw3,
    /// Cell values.
    pub phi: Raw<f64>,
    /// Boundary coefficient a (the zero-gradient mask).
    pub bc_a: Raw<f64>,
    /// Boundary constant b.
    pub bc_b: Raw<f64>,
    /// Jump per periodic face (empty: plain periodic).
    pub jump: Raw<f64>,
    /// Output gradient.
    pub out: Raw3,
}

/// Least-squares gradient of cell `c` ([`fv::lsq_gradient`]).
pub struct Gradient;
impl Entry for Gradient {
    type Args = GradientArgs;
    const NAME: &'static str = "gradient";
    #[inline(always)]
    unsafe fn run(c: usize, a: &GradientArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let bc = Boundary {
                a: a.bc_a.get(),
                b: a.bc_b.get(),
                jump: a.jump.get(),
            };
            let g = fv::lsq_gradient(c, &a.m.topo(), &a.w.v(), a.phi.get(), &bc);
            *a.out.x.out().add(c) = g[0];
            *a.out.y.out().add(c) = g[1];
            *a.out.z.out().add(c) = g[2];
        }
    }
}

/// Arguments of [`Limiter`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct LimiterArgs {
    /// Mesh.
    pub m: Mesh,
    /// Cell values.
    pub phi: Raw<f64>,
    /// Boundary coefficient a.
    pub bc_a: Raw<f64>,
    /// Boundary constant b.
    pub bc_b: Raw<f64>,
    /// Jump per periodic face (empty: plain periodic).
    pub jump: Raw<f64>,
    /// Gradient.
    pub grad: Raw3,
    /// ε² per cell or one value.
    pub eps2: Raw<f64>,
    /// Output ψ.
    pub out: Raw<f64>,
}

/// Venkatakrishnan limiter of cell `c` ([`fv::venkat_limiter`]).
pub struct Limiter;
impl Entry for Limiter {
    type Args = LimiterArgs;
    const NAME: &'static str = "limiter";
    #[inline(always)]
    unsafe fn run(c: usize, a: &LimiterArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let bc = Boundary {
                a: a.bc_a.get(),
                b: a.bc_b.get(),
                jump: a.jump.get(),
            };
            *a.out.out().add(c) = fv::venkat_limiter(
                c,
                &a.m.topo(),
                a.phi.get(),
                &bc,
                &a.grad.v(),
                &a.m.cell_centroid.v(),
                &a.m.face_centroid.v(),
                a.eps2.get(),
            );
        }
    }
}

/// Arguments of [`MassFlux`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct MassFluxArgs {
    /// Mesh.
    pub m: Mesh,
    /// Density.
    pub rho: f64,
    /// Cell velocity.
    pub u: Raw3,
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// Cell pressure.
    pub p: Raw<f64>,
    /// Pressure gradient.
    pub grad_p: Raw3,
    /// Boundary pressure.
    pub pb: Raw<f64>,
    /// Boundary velocity.
    pub ub: Raw3,
    /// D per cell (zeros: the plain interpolated flux).
    pub d: Raw<f64>,
    /// Relaxation term per face.
    pub relax: Raw<f64>,
    /// 1 with the time-step term.
    pub transient: u32,
    /// D^s per cell.
    pub ds: Raw<f64>,
    /// Old-level part Y per face.
    pub time: Raw<f64>,
    /// α_u.
    pub alpha: f64,
    /// Output per face.
    pub out: Raw<f64>,
}

/// Rhie–Chow face mass flux ([`k::face_mass_flux`]).
pub struct MassFlux;
impl Entry for MassFlux {
    type Args = MassFluxArgs;
    const NAME: &'static str = "mass_flux";
    #[inline(always)]
    unsafe fn run(f: usize, a: &MassFluxArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = k::FluxInputs {
                topo: a.m.topo(),
                geo: a.m.geo(),
                rho: a.rho,
                u: a.u.v(),
                grad: a.grad_u.v(),
                p: a.p.get(),
                grad_p: a.grad_p.v(),
                pb: a.pb.get(),
                pset: a.m.pset.get(),
                ub: a.ub.v(),
                flux_kind: a.m.flux_kind.get(),
                mom_kind: a.m.mom_kind.get(),
                d: a.d.get(),
                relax: a.relax.get(),
                transient: a.transient != 0,
                ds: a.ds.get(),
                time: a.time.get(),
                alpha: a.alpha,
            };
            *a.out.out().add(f) = k::face_mass_flux(f, &inp);
        }
    }
}

/// Arguments of [`MomentumFaces`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct MomentumFaceArgs {
    /// Mesh.
    pub m: Mesh,
    /// Convecting mass flux.
    pub mass_flux: Raw<f64>,
    /// Cell velocity.
    pub u: Raw3,
    /// Velocity of the central correction.
    pub uc: Raw3,
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// Limiter per component.
    pub psi: Raw3,
    /// Boundary velocity.
    pub ub: Raw3,
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
    /// 1 with the transposed viscous stress (turbulent flow).
    pub transposed: f64,
    /// Output c_P.
    pub cp: Raw<f64>,
    /// Output c_N.
    pub cn: Raw<f64>,
    /// Output b_f.
    pub bf: Raw3,
}

/// Momentum face coefficients ([`k::momentum_face`]).
pub struct MomentumFaces;
impl Entry for MomentumFaces {
    type Args = MomentumFaceArgs;
    const NAME: &'static str = "momentum_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &MomentumFaceArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = k::MomentumInputs {
                topo: a.m.topo(),
                geo: a.m.geo(),
                mass_flux: a.mass_flux.get(),
                mu: a.m.mu_face.get(),
                u: a.u.v(),
                uc: a.uc.v(),
                grad: a.grad_u.v(),
                psi: a.psi.v(),
                ub: a.ub.v(),
                mom_kind: a.m.mom_kind.get(),
                blend: a.blend,
                central: a.central,
                beta: a.beta,
                nonorth: a.nonorth,
                convection: a.convection,
                transposed: a.transposed,
            };
            let (cp, cn, b) = k::momentum_face(f, &inp);
            *a.cp.out().add(f) = cp;
            *a.cn.out().add(f) = cn;
            *a.bf.x.out().add(f) = b[0];
            *a.bf.y.out().add(f) = b[1];
            *a.bf.z.out().add(f) = b[2];
        }
    }
}

/// Arguments of [`MomentumRows`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct MomentumRowArgs {
    /// Mesh.
    pub m: Mesh,
    /// Face coefficients.
    pub cp: Raw<f64>,
    /// See `cp`.
    pub cn: Raw<f64>,
    /// Deferred parts.
    pub bf: Raw3,
    /// Pressure gradient.
    pub grad_p: Raw3,
    /// Body force per unit volume.
    pub body_force: [f64; 3],
    /// Current velocity.
    pub u: Raw3,
    /// α_u.
    pub alpha: f64,
    /// θ.
    pub theta: f64,
    /// 1 in a time step.
    pub transient: u32,
    /// ρ a₀ / Δt.
    pub tau: f64,
    /// ρ (−a₁) / Δt, ρ (−a₂) / Δt.
    pub c_old: [f64; 2],
    /// uⁿ.
    pub u_old: Raw3,
    /// uⁿ⁻¹.
    pub u_old2: Raw3,
    /// 1 − θ of the explicit operator.
    pub explicit_weight: f64,
    /// Explicit operator.
    pub explicit: Raw3,
    /// Output: matrix values.
    pub val: Raw<f64>,
    /// Output: right-hand sides.
    pub rhs: Raw3,
    /// Output: D.
    pub d: Raw<f64>,
    /// Output: D^s.
    pub ds: Raw<f64>,
    /// Output: D^C.
    pub dc: Raw<f64>,
    /// Output: relaxed diagonals.
    pub diag: Raw3,
}

/// Momentum rows ([`k::momentum_row`]).
pub struct MomentumRows;
impl Entry for MomentumRows {
    type Args = MomentumRowArgs;
    const NAME: &'static str = "momentum_rows";
    #[inline(always)]
    unsafe fn run(c: usize, a: &MomentumRowArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = k::MomentumRowInputs {
                topo: a.m.topo(),
                cp: a.cp.get(),
                cn: a.cn.get(),
                bf: a.bf.v(),
                grad_p: a.grad_p.v(),
                body_force: a.body_force,
                volume: a.m.volume.get(),
                u: a.u.v(),
                alpha: a.alpha,
                layout: a.m.layout(),
                mom_kind: a.m.mom_kind.get(),
                face_area: a.m.area.v(),
                theta: a.theta,
                transient: a.transient != 0,
                tau: a.tau,
                c_old: a.c_old,
                u_old: a.u_old.v(),
                u_old2: a.u_old2.v(),
                explicit_weight: a.explicit_weight,
                explicit: a.explicit.v(),
            };
            let out = k::MomentumRowOut {
                val: a.val.out(),
                rhs: [a.rhs.x.out(), a.rhs.y.out(), a.rhs.z.out()],
                d: a.d.out(),
                ds: a.ds.out(),
                dc: a.dc.out(),
                diag: [a.diag.x.out(), a.diag.y.out(), a.diag.z.out()],
            };
            k::momentum_row(c, &inp, &out);
        }
    }
}

/// Arguments of [`MomentumOperator`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct OperatorArgs {
    /// Mesh.
    pub m: Mesh,
    /// Face coefficients.
    pub cp: Raw<f64>,
    /// See `cp`.
    pub cn: Raw<f64>,
    /// Deferred parts.
    pub bf: Raw3,
    /// Velocity.
    pub u: Raw3,
    /// Output.
    pub out: Raw3,
}

/// Spatial momentum operator ([`k::momentum_operator`]).
pub struct MomentumOperator;
impl Entry for MomentumOperator {
    type Args = OperatorArgs;
    const NAME: &'static str = "momentum_operator";
    #[inline(always)]
    unsafe fn run(c: usize, a: &OperatorArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = k::OperatorInputs {
                topo: a.m.topo(),
                cp: a.cp.get(),
                cn: a.cn.get(),
                bf: a.bf.v(),
                mom_kind: a.m.mom_kind.get(),
                face_area: a.m.area.v(),
                u: a.u.v(),
            };
            let r = k::momentum_operator(c, &inp);
            *a.out.x.out().add(c) = r[0];
            *a.out.y.out().add(c) = r[1];
            *a.out.z.out().add(c) = r[2];
        }
    }
}

/// Arguments of [`ResidualScale`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ResidualScaleArgs {
    /// Velocity.
    pub u: Raw3,
    /// Right-hand sides.
    pub rhs: Raw3,
    /// Relaxed diagonals.
    pub diag: Raw3,
    /// α_u.
    pub alpha: f64,
    /// Output a u (3 per cell, interleaved).
    pub v: Raw<f64>,
    /// Output b (3 per cell, interleaved).
    pub w: Raw<f64>,
}

/// Parts of the momentum residual scale ([`k::residual_scale`]).
pub struct ResidualScale;
impl Entry for ResidualScale {
    type Args = ResidualScaleArgs;
    const NAME: &'static str = "residual_scale";
    #[inline(always)]
    unsafe fn run(c: usize, a: &ResidualScaleArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let (v, w) = k::residual_scale(c, &a.u.v(), &a.rhs.v(), &a.diag.v(), a.alpha);
            let mut i = 0;
            while i < 3 {
                *a.v.out().add(3 * c + i) = v[i];
                *a.w.out().add(3 * c + i) = w[i];
                i += 1;
            }
        }
    }
}

/// Arguments of the cell and face entries with one input and one output.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct UnaryArgs {
    /// Mesh.
    pub m: Mesh,
    /// Scalar parameter.
    pub s: f64,
    /// Input.
    pub x: Raw<f64>,
    /// Output.
    pub out: Raw<f64>,
}

/// Slot 0 of row `c` of `out` (matrix values) = `x[c]`.
pub struct SetDiagonal;
impl Entry for SetDiagonal {
    type Args = UnaryArgs;
    const NAME: &'static str = "set_diagonal";
    #[inline(always)]
    unsafe fn run(c: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(a.m.layout().pos(c, 0)) = ld(a.x.get(), c) }
    }
}

/// `out[c]` = 1 / slot 0 of row `c` of `x` (the Jacobi diagonal of a
/// matrix assembled on the device).
pub struct DiagonalInverse;
impl Entry for DiagonalInverse {
    type Args = UnaryArgs;
    const NAME: &'static str = "diagonal_inverse";
    #[inline(always)]
    unsafe fn run(c: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(c) = 1.0 / ld(a.x.get(), a.m.layout().pos(c, 0)) }
    }
}

/// Net outflow Σ ±ṁ of cell `c`.
pub struct NetOutflow;
impl Entry for NetOutflow {
    type Args = UnaryArgs;
    const NAME: &'static str = "net_outflow";
    #[inline(always)]
    unsafe fn run(c: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(c) = k::net_outflow(c, &a.m.topo(), a.x.get()) }
    }
}

/// Throughflow Σ |ṁ| of cell `c`.
pub struct Throughflow;
impl Entry for Throughflow {
    type Args = UnaryArgs;
    const NAME: &'static str = "throughflow";
    #[inline(always)]
    unsafe fn run(c: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(c) = k::throughflow(c, &a.m.topo(), a.x.get()) }
    }
}

/// Courant number of cell `c` (`s` = Δt / (2ρ)).
pub struct Courant;
impl Entry for Courant {
    type Args = UnaryArgs;
    const NAME: &'static str = "courant";
    #[inline(always)]
    unsafe fn run(c: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(c) = k::courant(c, &a.m.topo(), a.x.get(), a.m.volume.get(), a.s)
        }
    }
}

/// Pressure-correction coefficient of face `f` from D^C (`s` = ρ).
pub struct PressureFaces;
impl Entry for PressureFaces {
    type Args = UnaryArgs;
    const NAME: &'static str = "pressure_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(f) = k::pressure_face(
                f,
                &a.m.topo(),
                &a.m.geo(),
                a.m.flux_kind.get(),
                a.s,
                a.x.get(),
                a.m.closed.get(),
            )
        }
    }
}

/// Arguments of [`CloseBackflow`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CloseArgs {
    /// Mesh (its `closed` array is written).
    pub m: Mesh,
    /// Predicted face mass flux, updated in place on the boundary faces.
    pub mass_flux: Raw<f64>,
}

/// Backflow suppression of boundary face `fb` ([`k::close_backflow`]):
/// writes the flux of the face (the output is the input, at the same
/// index) and its closed state.
pub struct CloseBackflow;
impl Entry for CloseBackflow {
    type Args = CloseArgs;
    const NAME: &'static str = "close_backflow";
    #[inline(always)]
    unsafe fn run(fb: usize, a: &CloseArgs) {
        // SAFETY: guaranteed by the caller; index fb only.
        unsafe {
            let f = fb + a.m.n_interior as usize;
            let kind = ld(a.m.flux_kind.get(), fb);
            if kind == k::FLUX_PRESSURE_NO_BACKFLOW {
                let (v, c) = k::close_backflow(kind, a.mass_flux.read(f));
                *a.mass_flux.out().add(f) = v;
                *a.m.closed.out().add(fb) = c;
            } else {
                *a.m.closed.out().add(fb) = 0.0;
            }
        }
    }
}

/// Arguments of [`PressureRows`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct PressureRowArgs {
    /// Mesh.
    pub m: Mesh,
    /// Face coefficients.
    pub cf: Raw<f64>,
    /// Mass flux.
    pub mass_flux: Raw<f64>,
    /// Output: matrix values.
    pub val: Raw<f64>,
    /// Output: right-hand side.
    pub rhs: Raw<f64>,
}

/// Pressure-correction rows ([`k::pressure_row`]).
pub struct PressureRows;
impl Entry for PressureRows {
    type Args = PressureRowArgs;
    const NAME: &'static str = "pressure_rows";
    #[inline(always)]
    unsafe fn run(c: usize, a: &PressureRowArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            k::pressure_row(
                c,
                &a.m.topo(),
                a.cf.get(),
                a.mass_flux.get(),
                &a.m.layout(),
                a.val.out(),
                a.rhs.out(),
            )
        }
    }
}

/// Arguments of [`CorrectFluxes`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CorrectFluxArgs {
    /// Mesh.
    pub m: Mesh,
    /// Face coefficients.
    pub cf: Raw<f64>,
    /// Predicted mass flux.
    pub mass_flux: Raw<f64>,
    /// Pressure correction.
    pub pc: Raw<f64>,
    /// Output.
    pub out: Raw<f64>,
}

/// Corrected face flux ([`k::correct_flux`]).
pub struct CorrectFluxes;
impl Entry for CorrectFluxes {
    type Args = CorrectFluxArgs;
    const NAME: &'static str = "correct_fluxes";
    #[inline(always)]
    unsafe fn run(f: usize, a: &CorrectFluxArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(f) =
                k::correct_flux(f, &a.m.topo(), a.cf.get(), a.mass_flux.get(), a.pc.get())
        }
    }
}

/// Arguments of [`CorrectVelocity`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CorrectVelocityArgs {
    /// D^C.
    pub dc: Raw<f64>,
    /// Gradient of the pressure correction.
    pub grad: Raw3,
    /// Velocity, updated in place.
    pub u: Raw3,
}

/// u_i −= D^C ∇p'_i of cell `c`.
pub struct CorrectVelocity;
impl Entry for CorrectVelocity {
    type Args = CorrectVelocityArgs;
    const NAME: &'static str = "correct_velocity";
    #[inline(always)]
    unsafe fn run(c: usize, a: &CorrectVelocityArgs) {
        // SAFETY: guaranteed by the caller; index c only.
        unsafe {
            let d = a.dc.read(c);
            *a.u.x.out().add(c) = a.u.x.read(c) - d * a.grad.x.read(c);
            *a.u.y.out().add(c) = a.u.y.read(c) - d * a.grad.y.read(c);
            *a.u.z.out().add(c) = a.u.z.read(c) - d * a.grad.z.read(c);
        }
    }
}

/// Arguments of [`NeighbourCorrection`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct NeighbourArgs {
    /// Mesh.
    pub m: Mesh,
    /// Matrix values.
    pub val: Raw<f64>,
    /// Cell values.
    pub x: Raw<f64>,
    /// Diagonal.
    pub diag: Raw<f64>,
    /// Output.
    pub out: Raw<f64>,
}

/// −Σ_N val_PN x_N / diag of cell `c`.
pub struct NeighbourCorrection;
impl Entry for NeighbourCorrection {
    type Args = NeighbourArgs;
    const NAME: &'static str = "neighbour_correction";
    #[inline(always)]
    unsafe fn run(c: usize, a: &NeighbourArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(c) =
                -k::neighbour_sum(c, &a.m.topo(), &a.m.layout(), a.val.get(), a.x.get())
                    / a.diag.read(c)
        }
    }
}

/// Slot `k` of Aᵀ in the storage of A: `out[k] = x[transpose[k]]`, for
/// every stored slot (D59).
pub struct Transpose;
impl Entry for Transpose {
    type Args = UnaryArgs;
    const NAME: &'static str = "transpose";
    #[inline(always)]
    unsafe fn run(k: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(k) = adjoint::transposed_value(k, a.m.transpose.get(), a.x.get())
        }
    }
}

/// Arguments of [`MatrixProduct`] and [`TransposedProduct`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ProductArgs {
    /// Mesh.
    pub m: Mesh,
    /// Matrix values.
    pub val: Raw<f64>,
    /// Cell values.
    pub x: Raw<f64>,
    /// Output.
    pub out: Raw<f64>,
}

/// Row `c` of A x.
pub struct MatrixProduct;
impl Entry for MatrixProduct {
    type Args = ProductArgs;
    const NAME: &'static str = "matrix_product";
    #[inline(always)]
    unsafe fn run(c: usize, a: &ProductArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(c) =
                adjoint::row_product(c, &a.m.topo(), &a.m.layout(), a.val.get(), a.x.get())
        }
    }
}

/// Row `c` of Aᵀ x, without storing Aᵀ (D59).
pub struct TransposedProduct;
impl Entry for TransposedProduct {
    type Args = ProductArgs;
    const NAME: &'static str = "transposed_product";
    #[inline(always)]
    unsafe fn run(c: usize, a: &ProductArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(c) = adjoint::transposed_row_product(
                c,
                &a.m.topo(),
                &a.m.layout(),
                a.m.transpose.get(),
                a.val.get(),
                a.x.get(),
            )
        }
    }
}

/// Arguments of [`DeltaFlux`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct DeltaFluxArgs {
    /// Mesh.
    pub m: Mesh,
    /// Density.
    pub rho: f64,
    /// Velocity increment.
    pub delta: Raw3,
    /// Output per face.
    pub out: Raw<f64>,
}

/// Mass flux of a velocity increment ([`k::delta_flux`]).
pub struct DeltaFlux;
impl Entry for DeltaFlux {
    type Args = DeltaFluxArgs;
    const NAME: &'static str = "delta_flux";
    #[inline(always)]
    unsafe fn run(f: usize, a: &DeltaFluxArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(f) = k::delta_flux(
                f,
                &a.m.topo(),
                &a.m.geo(),
                a.m.flux_kind.get(),
                a.rho,
                &a.delta.v(),
                a.m.closed.get(),
            )
        }
    }
}

/// Arguments of [`ScalarFaces`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ScalarFaceArgs {
    /// Mesh.
    pub m: Mesh,
    /// Boundary coefficient a of the scalar per boundary face.
    pub bc_a: Raw<f64>,
    /// Boundary constant b.
    pub bc_b: Raw<f64>,
    /// Jump per periodic face (empty: plain periodic).
    pub jump: Raw<f64>,
    /// Convecting face mass flux.
    pub mass_flux: Raw<f64>,
    /// ρΓ per face.
    pub gamma: Raw<f64>,
    /// Cell values.
    pub phi: Raw<f64>,
    /// Cell gradient.
    pub grad: Raw3,
    /// Limiter per cell.
    pub psi: Raw<f64>,
    /// 1 linear upwind, 0 upwind.
    pub blend: f64,
    /// Fraction of the bounded central scheme.
    pub central: f64,
    /// β of its bounding.
    pub beta: f64,
    /// 1 with the non-orthogonal correction.
    pub nonorth: f64,
    /// Limit of the non-orthogonal correction (< 0: none).
    pub nonorth_limit: f64,
    /// Output c_P.
    pub cp: Raw<f64>,
    /// Output c_N.
    pub cn: Raw<f64>,
    /// Output b_f.
    pub bf: Raw<f64>,
}

/// Face coefficients of a scalar ([`fv::face_coefficients`]).
pub struct ScalarFaces;
impl Entry for ScalarFaces {
    type Args = ScalarFaceArgs;
    const NAME: &'static str = "scalar_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &ScalarFaceArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = fv::FaceInputs {
                topo: a.m.topo(),
                bc: Boundary {
                    a: a.bc_a.get(),
                    b: a.bc_b.get(),
                    jump: a.jump.get(),
                },
                mass_flux: a.mass_flux.get(),
                gamma: a.gamma.get(),
                dgeo: a.m.dgeo.get(),
                k: a.m.k.v(),
                weight: a.m.weight.get(),
                face_centroid: a.m.face_centroid.v(),
                cell_centroid: a.m.cell_centroid.v(),
                shift: a.m.shift.v(),
                phi: a.phi.get(),
                grad: a.grad.v(),
                psi: a.psi.get(),
                blend: a.blend,
                central: a.central,
                beta: a.beta,
                nonorth: a.nonorth,
                nonorth_limit: a.nonorth_limit,
            };
            let (cp, cn, b) = fv::face_coefficients(f, &inp);
            *a.cp.out().add(f) = cp;
            *a.cn.out().add(f) = cn;
            *a.bf.out().add(f) = b;
        }
    }
}

/// Arguments of [`ScalarRows`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ScalarRowArgs {
    /// Mesh.
    pub m: Mesh,
    /// Face coefficients.
    pub cp: Raw<f64>,
    /// See `cp`.
    pub cn: Raw<f64>,
    /// Deferred parts.
    pub bf: Raw<f64>,
    /// Convecting face mass flux (bounded form).
    pub mass_flux: Raw<f64>,
    /// Current cell values.
    pub phi: Raw<f64>,
    /// Density.
    pub rho: f64,
    /// S_u per cell (empty: none).
    pub su: Raw<f64>,
    /// S_p per cell (empty: none).
    pub sp: Raw<f64>,
    /// α.
    pub alpha: f64,
    /// 1 in a time step.
    pub transient: u64,
    /// 1/Δt.
    pub inv_dt: f64,
    /// a₀, a₁, a₂.
    pub a: [f64; 3],
    /// φⁿ.
    pub old: Raw<f64>,
    /// φⁿ⁻¹.
    pub old2: Raw<f64>,
    /// Lower bound of the bounded BDF2.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
    /// 1 in the cells with a fixed value (empty: none).
    pub fix_mask: Raw<f64>,
    /// The fixed values.
    pub fix_val: Raw<f64>,
    /// 1: a positive variable (the deferred corrections split by sign).
    pub positive: u32,
    /// Output: matrix values.
    pub val: Raw<f64>,
    /// Output: right-hand side.
    pub rhs: Raw<f64>,
    /// Output: relaxed diagonal.
    pub diag: Raw<f64>,
    /// Output: 1 where the bounded BDF2 fell back to first order.
    pub order: Raw<f64>,
}

/// Row of a scalar ([`crate::scalar::scalar_row`]).
pub struct ScalarRows;
impl Entry for ScalarRows {
    type Args = ScalarRowArgs;
    const NAME: &'static str = "scalar_rows";
    #[inline(always)]
    unsafe fn run(c: usize, a: &ScalarRowArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = crate::scalar::ScalarRowInputs {
                topo: a.m.topo(),
                cp: a.cp.get(),
                cn: a.cn.get(),
                bf: a.bf.get(),
                mass_flux: a.mass_flux.get(),
                volume: a.m.volume.get(),
                phi: a.phi.get(),
                rho: a.rho,
                su: a.su.get(),
                sp: a.sp.get(),
                alpha: a.alpha,
                layout: a.m.layout(),
                transient: a.transient != 0,
                inv_dt: a.inv_dt,
                a: a.a,
                old: a.old.get(),
                old2: a.old2.get(),
                lo: a.lo,
                hi: a.hi,
                fix_mask: a.fix_mask.get(),
                fix_val: a.fix_val.get(),
                positive: a.positive != 0,
            };
            let r = crate::scalar::scalar_row(c, &inp, a.val.out());
            *a.rhs.out().add(c) = r.rhs;
            *a.diag.out().add(c) = r.diag;
            *a.order.out().add(c) = r.first_order;
        }
    }
}

/// Arguments of [`ScalarScale`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ScalarScaleArgs {
    /// Cell values.
    pub phi: Raw<f64>,
    /// Right-hand side.
    pub rhs: Raw<f64>,
    /// Relaxed diagonal.
    pub diag: Raw<f64>,
    /// α.
    pub alpha: f64,
    /// Output a φ.
    pub v: Raw<f64>,
    /// Output b − (a/α − a) φ.
    pub w: Raw<f64>,
}

/// Parts of the residual scale of a scalar
/// ([`crate::scalar::residual_parts`]).
pub struct ScalarScale;
impl Entry for ScalarScale {
    type Args = ScalarScaleArgs;
    const NAME: &'static str = "scalar_scale";
    #[inline(always)]
    unsafe fn run(c: usize, a: &ScalarScaleArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let (v, w) = crate::scalar::residual_parts(
                ld(a.phi.get(), c),
                ld(a.rhs.get(), c),
                ld(a.diag.get(), c),
                a.alpha,
            );
            *a.v.out().add(c) = v;
            *a.w.out().add(c) = w;
        }
    }
}

/// Arguments of [`SstCells`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct SstCellArgs {
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// k.
    pub k: Raw<f64>,
    /// ω.
    pub w: Raw<f64>,
    /// ∇k.
    pub grad_k: Raw3,
    /// ∇ω.
    pub grad_w: Raw3,
    /// Wall distance.
    pub d: Raw<f64>,
    /// ν.
    pub nu: f64,
    /// Constants.
    pub c: crate::turb::sst::SstConstants,
    /// Options.
    pub o: crate::turb::sst::SstOptions,
    /// Additional explicit sources of k and ω per cell (MMS, D104);
    /// empty: none.
    pub add_k: Raw<f64>,
    /// See `add_k`.
    pub add_w: Raw<f64>,
    /// Weight m of the wall production per cell (the automatic wall
    /// treatment, 0 off the near-wall cells); empty: none.
    pub wall_mask: Raw<f64>,
    /// The weighted wall production per cell, added to (1 − m) P̃_k.
    pub wall_pk: Raw<f64>,
    /// Output ν_t.
    pub nu_t: Raw<f64>,
    /// Output F₁.
    pub f1: Raw<f64>,
    /// Output P̃_k.
    pub pk: Raw<f64>,
    /// Output S_u of k (empty: ν_t only).
    pub su_k: Raw<f64>,
    /// Output S_p of k.
    pub sp_k: Raw<f64>,
    /// Output S_u of ω.
    pub su_w: Raw<f64>,
    /// Output S_p of ω.
    pub sp_w: Raw<f64>,
    /// Output σ_k ν_t.
    pub dk: Raw<f64>,
    /// Output σ_ω ν_t.
    pub dw: Raw<f64>,
    /// Output 1 where the cap of ν_t/ν acted.
    pub capped: Raw<f64>,
}

/// The k-ω SST model per cell ([`crate::turb::sst::sst_cell`]).
pub struct SstCells;
impl Entry for SstCells {
    type Args = SstCellArgs;
    const NAME: &'static str = "sst_cells";
    #[inline(always)]
    unsafe fn run(c: usize, a: &SstCellArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let g = a.grad_u.v();
            let gu = [g.x.at(c), g.y.at(c), g.z.at(c)];
            let inp = crate::turb::sst::SstCellInputs {
                k: ld(a.k.get(), c),
                w: ld(a.w.get(), c),
                grad_k: a.grad_k.v().at(c),
                grad_w: a.grad_w.v().at(c),
                grad_u: gu,
                d: ld(a.d.get(), c),
                nu: a.nu,
            };
            let r = crate::turb::sst::sst_cell(&inp, &a.c, &a.o);
            *a.nu_t.out().add(c) = r.nu_t;
            *a.f1.out().add(c) = r.f1;
            *a.pk.out().add(c) = r.pk;
            *a.dk.out().add(c) = r.dk;
            *a.dw.out().add(c) = r.dw;
            *a.capped.out().add(c) = r.capped;
            if a.su_k.len > 0 {
                let (ak, aw) = if a.add_k.len > 0 {
                    (ld(a.add_k.get(), c), ld(a.add_w.get(), c))
                } else {
                    (0.0, 0.0)
                };
                // an added source splits like the others: its positive
                // part explicit, its negative part implicit
                let (k, w) = (inp.k, inp.w);
                // near-wall cells of the automatic treatment: the
                // production (1 − m) P̃_k + m u_τ³/(κ y_P) (the weighted
                // wall part comes from [`WallCells`])
                let su_k = if a.wall_mask.len > 0 {
                    let m = ld(a.wall_mask.get(), c);
                    r.su_k - m * r.pk + ld(a.wall_pk.get(), c)
                } else {
                    r.su_k
                };
                *a.su_k.out().add(c) = su_k + ak.max(0.0);
                *a.sp_k.out().add(c) = r.sp_k + ak.min(0.0) / k.max(1e-300);
                *a.su_w.out().add(c) = r.su_w + aw.max(0.0);
                *a.sp_w.out().add(c) = r.sp_w + aw.min(0.0) / w.max(1e-300);
            }
        }
    }
}

/// Arguments of [`KepsCells`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct KepsCellArgs {
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// k.
    pub k: Raw<f64>,
    /// ε.
    pub e: Raw<f64>,
    /// ν.
    pub nu: f64,
    /// Constants.
    pub c: crate::turb::keps::KepsConstants,
    /// Options.
    pub o: crate::turb::keps::KepsOptions,
    /// Additional explicit sources of k and ε per cell (MMS, D104);
    /// empty: none.
    pub add_k: Raw<f64>,
    /// See `add_k`.
    pub add_e: Raw<f64>,
    /// Weight m_p of the wall production per cell (0 off the near-wall
    /// cells); empty: none.
    pub wall_mask: Raw<f64>,
    /// The weighted wall production per cell, added to (1 − m_p) P_k.
    pub wall_pk: Raw<f64>,
    /// Output ν_t.
    pub nu_t: Raw<f64>,
    /// Output C_μ.
    pub c_mu: Raw<f64>,
    /// Output P_k.
    pub pk: Raw<f64>,
    /// Output S_u of k (empty: ν_t only).
    pub su_k: Raw<f64>,
    /// Output S_p of k.
    pub sp_k: Raw<f64>,
    /// Output S_u of ε.
    pub su_e: Raw<f64>,
    /// Output S_p of ε.
    pub sp_e: Raw<f64>,
    /// Output ν_t/σ_k.
    pub dk: Raw<f64>,
    /// Output ν_t/σ_ε.
    pub de: Raw<f64>,
    /// Output 1 where the cap of ν_t/ν acted.
    pub capped: Raw<f64>,
}

/// The standard or realizable k-ε model per cell
/// ([`crate::turb::keps::keps_cell`]).
pub struct KepsCells;
impl Entry for KepsCells {
    type Args = KepsCellArgs;
    const NAME: &'static str = "keps_cells";
    #[inline(always)]
    unsafe fn run(c: usize, a: &KepsCellArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let g = a.grad_u.v();
            let inp = crate::turb::keps::KepsCellInputs {
                k: ld(a.k.get(), c),
                e: ld(a.e.get(), c),
                grad_u: [g.x.at(c), g.y.at(c), g.z.at(c)],
                nu: a.nu,
            };
            let r = crate::turb::keps::keps_cell(&inp, &a.c, &a.o);
            *a.nu_t.out().add(c) = r.nu_t;
            *a.c_mu.out().add(c) = r.c_mu;
            *a.pk.out().add(c) = r.pk;
            *a.dk.out().add(c) = r.dk;
            *a.de.out().add(c) = r.de;
            *a.capped.out().add(c) = r.capped;
            if a.su_k.len > 0 {
                let (ak, ae) = if a.add_k.len > 0 {
                    (ld(a.add_k.get(), c), ld(a.add_e.get(), c))
                } else {
                    (0.0, 0.0)
                };
                let (k, e) = (inp.k, inp.e);
                // near-wall cells: (1 − m_p) P_k + the wall production
                let su_k = if a.wall_mask.len > 0 {
                    let m = ld(a.wall_mask.get(), c);
                    r.su_k - m * r.pk + ld(a.wall_pk.get(), c)
                } else {
                    r.su_k
                };
                *a.su_k.out().add(c) = su_k + ak.max(0.0);
                *a.sp_k.out().add(c) = r.sp_k + ak.min(0.0) / k.max(1e-300);
                *a.su_e.out().add(c) = r.su_e + ae.max(0.0);
                *a.sp_e.out().add(c) = r.sp_e + ae.min(0.0) / e.max(1e-300);
            }
        }
    }
}

/// Arguments of [`TurbFaces`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct TurbFaceArgs {
    /// Mesh.
    pub m: Mesh,
    /// ν_t per cell.
    pub nu_t: Raw<f64>,
    /// σ_k ν_t per cell.
    pub dk: Raw<f64>,
    /// σ_ω ν_t per cell.
    pub dw: Raw<f64>,
    /// Laminar viscosity per face.
    pub mu: Raw<f64>,
    /// 1 on the wall faces (per boundary face): no turbulent diffusion
    /// through a smooth wall.
    pub wall: Raw<f64>,
    /// Density.
    pub rho: f64,
    /// Output μ + ρ ν_t per face.
    pub mu_eff: Raw<f64>,
    /// Output ρ(ν + σ_k ν_t) per face.
    pub gamma_k: Raw<f64>,
    /// Output ρ(ν + σ_ω ν_t) per face.
    pub gamma_w: Raw<f64>,
}

/// Face viscosities of the turbulent flow: linear interpolation of the
/// cell values on interior and periodic faces, the owner value on boundary
/// faces, none on walls.
pub struct TurbFaces;
impl Entry for TurbFaces {
    type Args = TurbFaceArgs;
    const NAME: &'static str = "turb_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &TurbFaceArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let ni = a.m.n_interior as usize;
            let p = ld(a.m.owner.get(), f) as usize;
            let mu = ld(a.mu.get(), f);
            let wall = f >= ni && ld(a.wall.get(), f - ni) > 0.0;
            let (nt, gk, gw) = if f < ni {
                let n = ld(a.m.neighbour.get(), f) as usize;
                let w = ld(a.m.weight.get(), f);
                let wn = 1.0 - w;
                let at = |b: Raw<f64>| w * ld(b.get(), p) + wn * ld(b.get(), n);
                (at(a.nu_t), at(a.dk), at(a.dw))
            } else {
                let s = 1.0 - ld(a.wall.get(), f - ni);
                (
                    s * ld(a.nu_t.get(), p),
                    s * ld(a.dk.get(), p),
                    s * ld(a.dw.get(), p),
                )
            };
            // wall faces: μ_w of the wall treatment ([`WallFaces`])
            if !wall {
                *a.mu_eff.out().add(f) = mu + a.rho * nt;
            }
            *a.gamma_k.out().add(f) = mu + a.rho * gk;
            *a.gamma_w.out().add(f) = mu + a.rho * gw;
        }
    }
}

/// Arguments of [`WallFaces`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct WallFaceArgs {
    /// Mesh.
    pub m: Mesh,
    /// The wall faces (global face index).
    pub faces: Raw<u32>,
    /// Unit normal of each wall face (outward).
    pub normal: Raw3,
    /// y_P of each wall face.
    pub y: Raw<f64>,
    /// k_s of each wall face.
    pub ks: Raw<f64>,
    /// Cell velocity.
    pub u: Raw3,
    /// Boundary face velocity (per boundary face).
    pub ub: Raw3,
    /// ν.
    pub nu: f64,
    /// Density.
    pub rho: f64,
    /// Treatment ([`crate::turb::wall::LOW_RE`], `AUTOMATIC`).
    pub treatment: u32,
    /// Constants.
    pub c: crate::turb::sst::SstConstants,
    /// Constants of the k-ε family (treatments `SCALABLE`, `STANDARD`,
    /// `COMPOUND`).
    pub ke: crate::turb::keps::KepsConstants,
    /// k per cell (the k-ε treatments).
    pub k: Raw<f64>,
    /// C_s of Cebeci & Bradshaw per wall face (the k-ε treatments).
    pub cs: Raw<f64>,
    /// u_τ of the last iteration per wall face (rough walls); may be the
    /// output `u_tau` (read before written, same index).
    pub u_tau_old: Raw<f64>,
    /// Output u_τ per wall face.
    pub u_tau: Raw<f64>,
    /// Output y⁺ per wall face.
    pub y_plus: Raw<f64>,
    /// Output ω per wall face.
    pub omega: Raw<f64>,
    /// Output: the viscous weight of the blending per wall face.
    pub blending: Raw<f64>,
    /// Output: the weighted wall production per wall face (SST automatic:
    /// (1 − b) u_τ³/(κ y_P)).
    pub pk: Raw<f64>,
    /// Output: the weight m_p of the wall production per wall face (the
    /// cell keeps (1 − m_p) of its own production).
    pub prod_w: Raw<f64>,
    /// Output (per face): μ_w on the wall faces.
    pub mu_eff: Raw<f64>,
    /// Output (per boundary face): the wall value of k (bc_b of k).
    pub k_b: Raw<f64>,
    /// Output (per boundary face): the wall value of ω (bc_b of ω).
    pub w_b: Raw<f64>,
}

/// The wall treatment per wall face ([`crate::turb::wall::wall_face`]).
pub struct WallFaces;
impl Entry for WallFaces {
    type Args = WallFaceArgs;
    const NAME: &'static str = "wall_faces";
    #[inline(always)]
    unsafe fn run(i: usize, a: &WallFaceArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let f = ld(a.faces.get(), i) as usize;
            let fb = f - a.m.n_interior as usize;
            let p = ld(a.m.owner.get(), f) as usize;
            let n = a.normal.v().at(i);
            let up = a.u.v().at(p);
            let uw = a.ub.v().at(fb);
            let du = [up[0] - uw[0], up[1] - uw[1], up[2] - uw[2]];
            let dn = du[0] * n[0] + du[1] * n[1] + du[2] * n[2];
            let t = [du[0] - dn * n[0], du[1] - dn * n[1], du[2] - dn * n[2]];
            let ut = foehncraft_core::math::sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
            if a.treatment >= crate::turb::wall::SCALABLE {
                // the k-ε family: ε fixed in the cell (weight 1), k and ε
                // with zero gradient at the wall
                let inp = crate::turb::wall::KepsWallInputs {
                    ut,
                    y: ld(a.y.get(), i),
                    nu: a.nu,
                    k: ld(a.k.get(), p),
                    ks: ld(a.ks.get(), i),
                    cs: ld(a.cs.get(), i),
                };
                let r = crate::turb::wall::keps_wall_face(&inp, &a.ke, a.treatment);
                *a.u_tau.out().add(i) = r.u_tau;
                *a.y_plus.out().add(i) = r.y_plus;
                *a.omega.out().add(i) = r.epsilon;
                *a.blending.out().add(i) = 0.0;
                *a.pk.out().add(i) = r.pk_wall;
                *a.prod_w.out().add(i) = r.prod_weight;
                *a.mu_eff.out().add(f) = a.rho * r.nu_wall;
                *a.k_b.out().add(fb) = 0.0;
                *a.w_b.out().add(fb) = 0.0;
                return;
            }
            let inp = crate::turb::wall::WallInputs {
                ut,
                y: ld(a.y.get(), i),
                nu: a.nu,
                ks: ld(a.ks.get(), i),
                u_tau_old: a.u_tau_old.read(i),
            };
            let r = crate::turb::wall::wall_face(&inp, &a.c, a.treatment);
            let m = 1.0 - r.blending;
            *a.u_tau.out().add(i) = r.u_tau;
            *a.y_plus.out().add(i) = r.y_plus;
            *a.omega.out().add(i) = r.omega;
            *a.blending.out().add(i) = r.blending;
            *a.pk.out().add(i) = m * r.pk_wall;
            *a.prod_w.out().add(i) = m;
            *a.mu_eff.out().add(f) = a.rho * r.nu_wall;
            *a.k_b.out().add(fb) = r.k_wall;
            *a.w_b.out().add(fb) = r.omega_wall;
        }
    }
}

/// Arguments of [`WallCells`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct WallCellArgs {
    /// The near-wall cells.
    pub cells: Raw<u32>,
    /// Start of the wall-face list of each near-wall cell (length + 1).
    pub ptr: Raw<u32>,
    /// Wall-face slots (indices into the wall-face arrays).
    pub slots: Raw<u32>,
    /// Area weight of each slot (the weights of a cell sum to 1).
    pub weight: Raw<f64>,
    /// ω per wall face.
    pub omega: Raw<f64>,
    /// Viscous weight of the blending per wall face (1: low-Re).
    pub blending: Raw<f64>,
    /// Weighted wall production per wall face.
    pub pk: Raw<f64>,
    /// Weight of the wall production per wall face.
    pub prod_w: Raw<f64>,
    /// Output: the area mean of `pk` per cell.
    pub wall_pk: Raw<f64>,
    /// Output: the area mean of `prod_w` per cell.
    pub prod_mask: Raw<f64>,
    /// Output: the fixed value per cell.
    pub fix_val: Raw<f64>,
    /// Output: the weight of the fixed value per cell, the area mean of
    /// 1 − blending over the smooth wall faces.
    pub fix_mask: Raw<f64>,
    /// The fixed variable, or [`Raw::EMPTY`]: if given, its value in the
    /// near-wall cells moves to the fixed one by the weight, the start
    /// value of the solve (the k-ε family: a start far from the fixed ε
    /// left the solve of the first iteration negative there).
    pub phi: Raw<f64>,
}

/// The area-weighted mean of the wall-face values of each near-wall cell,
/// in the stored order of its wall faces (stage 6 specification 4.1).
pub struct WallCells;
impl Entry for WallCells {
    type Args = WallCellArgs;
    const NAME: &'static str = "wall_cells";
    #[inline(always)]
    unsafe fn run(i: usize, a: &WallCellArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let c = ld(a.cells.get(), i) as usize;
            let lo = ld(a.ptr.get(), i) as usize;
            let hi = ld(a.ptr.get(), i + 1) as usize;
            let (mut v, mut m, mut p, mut q) = (0.0, 0.0, 0.0, 0.0);
            let mut j = lo;
            while j < hi {
                let s = ld(a.slots.get(), j) as usize;
                let w = ld(a.weight.get(), s);
                v += w * ld(a.omega.get(), s);
                m += w * (1.0 - ld(a.blending.get(), s));
                p += w * ld(a.pk.get(), s);
                q += w * ld(a.prod_w.get(), s);
                j += 1;
            }
            *a.fix_val.out().add(c) = v;
            *a.fix_mask.out().add(c) = m;
            if a.phi.len > 0 {
                let old = a.phi.read(c);
                *a.phi.out().add(c) = old + m * (v - old);
            }
            *a.wall_pk.out().add(c) = p;
            *a.prod_mask.out().add(c) = q;
        }
    }
}

/// Arguments of [`Floor`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct FloorArgs {
    /// Values (input and output).
    pub x: Raw<f64>,
    /// Lower bound.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
    /// Output: 1 where a bound acted, else 0.
    pub hit: Raw<f64>,
}

/// x ← min(max(x, lo), hi), and where a bound acted (the positivity floors
/// of the turbulence variables, stage 6 specification 6.2; counted in the
/// health section).
pub struct Floor;
impl Entry for Floor {
    type Args = FloorArgs;
    const NAME: &'static str = "floor";
    #[inline(always)]
    unsafe fn run(i: usize, a: &FloorArgs) {
        // SAFETY: guaranteed by the caller; x is read before it is written.
        unsafe {
            let v = a.x.read(i);
            let r = v.max(a.lo).min(a.hi);
            *a.hit.out().add(i) = if r != v { 1.0 } else { 0.0 };
            *a.x.out().add(i) = r;
        }
    }
}

send_sync!(
    LinArgs,
    BlockArgs,
    DotArgs,
    CloseArgs,
    ConstArgs,
    BoundaryValuesArgs,
    GradientArgs,
    LimiterArgs,
    MassFluxArgs,
    MomentumFaceArgs,
    MomentumRowArgs,
    OperatorArgs,
    ResidualScaleArgs,
    UnaryArgs,
    PressureRowArgs,
    CorrectFluxArgs,
    CorrectVelocityArgs,
    NeighbourArgs,
    DeltaFluxArgs,
    ProductArgs,
    ScalarFaceArgs,
    ScalarRowArgs,
    ScalarScaleArgs,
    SstCellArgs,
    KepsCellArgs,
    TurbFaceArgs,
    WallFaceArgs,
    WallCellArgs,
    FloorArgs
);

/// Names of all entries (the PTX module has `flow_<name>` for each).
pub const NAMES: [&str; 41] = [
    Lincomb::NAME,
    ScaledDifference::NAME,
    Fill::NAME,
    Subtract::NAME,
    ScaleEntry::NAME,
    BlockSumSq::NAME,
    BlockDot::NAME,
    BlockMinMax::NAME,
    BoundaryConstants::NAME,
    BoundaryValues::NAME,
    Gradient::NAME,
    Limiter::NAME,
    MassFlux::NAME,
    MomentumFaces::NAME,
    MomentumRows::NAME,
    MomentumOperator::NAME,
    ResidualScale::NAME,
    SetDiagonal::NAME,
    DiagonalInverse::NAME,
    NetOutflow::NAME,
    Throughflow::NAME,
    Courant::NAME,
    PressureFaces::NAME,
    CloseBackflow::NAME,
    PressureRows::NAME,
    CorrectFluxes::NAME,
    CorrectVelocity::NAME,
    NeighbourCorrection::NAME,
    DeltaFlux::NAME,
    Transpose::NAME,
    MatrixProduct::NAME,
    TransposedProduct::NAME,
    ScalarFaces::NAME,
    ScalarRows::NAME,
    ScalarScale::NAME,
    SstCells::NAME,
    KepsCells::NAME,
    TurbFaces::NAME,
    WallFaces::NAME,
    WallCells::NAME,
    Floor::NAME,
];

#[cfg(test)]
mod tests {
    use super::*;
    use core::mem::{align_of, size_of};

    /// CUDA limits the parameters of a kernel to 4096 bytes (Volta and
    /// newer accept more only with CUDA 12.1 and later); the count `n` takes
    /// another 8.
    #[test]
    fn argument_blocks_fit_the_kernel_parameter_limit() {
        let sizes = [
            size_of::<LinArgs>(),
            size_of::<BlockArgs>(),
            size_of::<DotArgs>(),
            size_of::<CloseArgs>(),
            size_of::<ConstArgs>(),
            size_of::<BoundaryValuesArgs>(),
            size_of::<GradientArgs>(),
            size_of::<LimiterArgs>(),
            size_of::<MassFluxArgs>(),
            size_of::<MomentumFaceArgs>(),
            size_of::<MomentumRowArgs>(),
            size_of::<OperatorArgs>(),
            size_of::<ResidualScaleArgs>(),
            size_of::<UnaryArgs>(),
            size_of::<PressureRowArgs>(),
            size_of::<CorrectFluxArgs>(),
            size_of::<CorrectVelocityArgs>(),
            size_of::<NeighbourArgs>(),
            size_of::<DeltaFluxArgs>(),
            size_of::<ProductArgs>(),
            size_of::<ScalarFaceArgs>(),
            size_of::<ScalarRowArgs>(),
            size_of::<ScalarScaleArgs>(),
            size_of::<SstCellArgs>(),
            size_of::<KepsCellArgs>(),
            size_of::<TurbFaceArgs>(),
            size_of::<WallFaceArgs>(),
            size_of::<WallCellArgs>(),
            size_of::<FloorArgs>(),
        ];
        for s in sizes {
            assert!(s + 8 <= 4096, "argument block of {s} bytes");
        }
    }

    /// The layout the PTX entries see: 64-bit addresses and lengths, 8-byte
    /// alignment, no padding inside an array descriptor.
    #[test]
    fn array_descriptors_are_two_words() {
        assert_eq!(size_of::<Raw<f64>>(), 16);
        assert_eq!(size_of::<Raw<u32>>(), 16);
        assert_eq!(align_of::<Raw<f64>>(), 8);
        assert_eq!(size_of::<Raw3>(), 48);
        assert_eq!(
            size_of::<Mesh>(),
            3 * 8
                + 4 * 16
                + 3 * 48
                + 16
                + 48
                + 16
                + 16
                + 2 * 48
                + 8
                + 2 * 16
                + 4 * 16
                + 16
                + 16
                + 16
                + 48
                + 16
                + 8
                + 48
        );
    }

    /// Every entry body of the element-wise operations computes what the
    /// documentation says, on host memory.
    #[test]
    fn element_wise_entries() {
        let x = [1.0, 2.0, 3.0];
        let mut y = [10.0, 20.0, 30.0];
        let r = |v: &[f64]| Raw::<f64>::new(v.as_ptr() as u64, v.len());
        let a = LinArgs {
            a: 2.0,
            b: -1.0,
            x: r(&x),
            y: r(&y),
            out: Raw::new(y.as_mut_ptr() as u64, 3),
        };
        for i in 0..3 {
            // SAFETY: index i of arrays of length 3; out is y (allowed).
            unsafe { Lincomb::run(i, &a) };
        }
        assert_eq!(y, [-8.0, -16.0, -24.0]);
        assert_eq!(NAMES.len(), 41);
    }
}

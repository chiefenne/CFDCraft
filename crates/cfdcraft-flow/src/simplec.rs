//! The steady incompressible flow solver: the SIMPLEC algorithm (Van
//! Doormaal & Raithby 1984), or SIMPLE (Patankar & Spalding 1972), on the
//! collocated, cell-centred grid with the Rhie–Chow face mass flux (D19,
//! D27).
//!
//! One outer iteration:
//!
//! 1. boundary values of velocity and pressure (extrapolated, symmetry,
//!    from the face flux or the total pressure) from the current fields;
//! 2. least-squares gradients of u, v, w and p (D24), limiters (D25);
//! 3. the relaxation term of the Rhie–Chow flux from the current fields,
//!    (1 − α_u)(ṁ_f − ρ ū_f·S);
//! 4. momentum: face coefficients, one matrix for the three components
//!    (relaxed with α_u), BiCGStab with Jacobi (D36) per component;
//! 5. predicted mass flux ṁ*_f by Rhie–Chow (relaxation independent: the
//!    converged result does not depend on α_u, Majumdar 1988);
//! 6. pressure correction: Σ c_f (p'_P − p'_N) = −Σ ±ṁ*_f with
//!    c_f = ρ D̄^C_f |S|²/(d·S), D^C = V / (a_P/α_u − Σ|a_N|) (SIMPLEC) or
//!    V / (a_P/α_u) (SIMPLE); CG with the smoothed-aggregation multigrid of
//!    stage 3 (D35), the hierarchy rebuilt every `amg_rebuild` iterations;
//! 7. corrections: ṁ_f −= c_f (p'_N − p'_P) (so the face fluxes satisfy
//!    continuity to the tolerance of the pressure solve), u −= D^C ∇p',
//!    p += α_p p'.
//!
//! The loops run on the host (CPU threads) and call the single-source bodies
//! of `cfdcraft_kernels::flow`; the linear solves run on a backend of
//! `cfdcraft-fv`. Each output is written by one index, so the result does
//! not depend on the number of threads (D31).
//!
//! Residuals (dimensionless, per outer iteration, before the solves):
//! momentum ‖b − A u_i‖₂ / ‖a_P u‖₂ per component with ‖a_P u‖ over all
//! three components, continuity ‖Σ_f ±ṁ*_f‖₂ / ‖Σ_f |ṁ_f|‖₂ (net outflow
//! of the predicted fluxes over the cell throughflow).

use crate::bc::BoundaryArrays;
use anyhow::{Result, ensure};
use cfdcraft_fv::amg::AmgSettings;
use cfdcraft_fv::backend::{Backend, V};
use cfdcraft_fv::disc::{Discretisation, Soa3};
use cfdcraft_fv::problem::{Boundary, BoundaryKind, Problem};
use cfdcraft_fv::solver::{self as linear, LimiterKind, Preconditioner, Scheme};
use cfdcraft_kernels::flow::{self as k, FaceGeometry, VelocityGradient};
use cfdcraft_kernels::fv::{self, SellLayout};
use cfdcraft_mesh::{Geometry, Mesh};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// Indices per task of the host loops.
const GRAIN: usize = 1 << 12;

/// Meshes with fewer cells run on one thread.
pub const SERIAL_BELOW: usize = 20_000;

/// Pressure–velocity algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "kebab-case")]
pub enum Algorithm {
    /// SIMPLEC (Van Doormaal & Raithby 1984).
    #[default]
    Simplec,
    /// SIMPLE (Patankar & Spalding 1972).
    Simple,
}

/// `[solver]` settings of a flow case.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Settings {
    /// Pressure–velocity algorithm.
    pub algorithm: Algorithm,
    /// Convection scheme of the momentum equations.
    pub scheme: Scheme,
    /// Limiter of the velocity gradients.
    pub limiter: LimiterKind,
    /// Constant K of the limiter.
    pub limiter_k: f64,
    /// Non-orthogonal correction of the viscous term (D26).
    pub nonorth_correction: bool,
    /// Under-relaxation of the momentum equations α_u.
    pub relaxation_velocity: f64,
    /// Relaxation of the pressure correction α_p (default 1 for SIMPLEC,
    /// 0.3 for SIMPLE).
    pub relaxation_pressure: Option<f64>,
    /// Largest number of outer iterations.
    pub max_iterations: usize,
    /// Convergence: every residual below this value.
    pub tolerance: f64,
    /// Relative tolerance of each momentum solve.
    pub momentum_tolerance: f64,
    /// Largest number of BiCGStab iterations per momentum solve.
    pub momentum_max_iter: usize,
    /// Relative tolerance of each pressure-correction solve.
    pub pressure_tolerance: f64,
    /// Largest number of CG iterations per pressure solve.
    pub pressure_max_iter: usize,
    /// Preconditioner of the pressure solve.
    pub pressure_preconditioner: Preconditioner,
    /// Multigrid settings.
    pub amg: AmgSettings,
    /// Rebuild the multigrid hierarchy every N outer iterations.
    pub amg_rebuild: usize,
    /// Log a line every N outer iterations.
    pub log_every: usize,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            algorithm: Algorithm::Simplec,
            scheme: Scheme::LinearUpwind,
            limiter: LimiterKind::None,
            limiter_k: 5.0,
            nonorth_correction: true,
            relaxation_velocity: 0.7,
            relaxation_pressure: None,
            max_iterations: 2000,
            tolerance: 1e-8,
            momentum_tolerance: 1e-2,
            momentum_max_iter: 50,
            pressure_tolerance: 1e-2,
            pressure_max_iter: 100,
            pressure_preconditioner: Preconditioner::Amg,
            amg: AmgSettings::default(),
            amg_rebuild: 10,
            log_every: 10,
        }
    }
}

impl Settings {
    /// α_p.
    pub fn alpha_p(&self) -> f64 {
        self.relaxation_pressure.unwrap_or(match self.algorithm {
            Algorithm::Simplec => 1.0,
            Algorithm::Simple => 0.3,
        })
    }

    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.relaxation_velocity > 0.0 && self.relaxation_velocity <= 1.0,
            "[solver] relaxation_velocity must lie in (0, 1]"
        );
        let ap = self.alpha_p();
        ensure!(
            ap > 0.0 && ap <= 1.0,
            "[solver] relaxation_pressure must lie in (0, 1]"
        );
        ensure!(
            !(self.algorithm == Algorithm::Simplec && self.relaxation_velocity >= 1.0),
            "[solver] SIMPLEC needs relaxation_velocity < 1 for a steady case \
             (D^C = V/(a_P/α − Σ|a_N|) is singular at α = 1); try 0.7"
        );
        ensure!(self.tolerance > 0.0, "[solver] tolerance must be > 0");
        ensure!(self.limiter_k >= 0.0, "[solver] limiter_k must be >= 0");
        for (name, t) in [
            ("momentum_tolerance", self.momentum_tolerance),
            ("pressure_tolerance", self.pressure_tolerance),
        ] {
            ensure!(t > 0.0 && t < 1.0, "[solver] {name} must lie in (0, 1)");
        }
        ensure!(self.amg_rebuild >= 1, "[solver] amg_rebuild must be >= 1");
        ensure!(self.log_every >= 1, "[solver] log_every must be >= 1");
        self.amg.validate()
    }
}

/// Physical inputs.
#[derive(Debug, Clone, Copy)]
pub struct Physics {
    /// Density ρ (kg/m³).
    pub rho: f64,
    /// Dynamic viscosity μ (Pa s).
    pub mu: f64,
    /// Body force per unit volume (N/m³).
    pub body_force: [f64; 3],
    /// Convection on (false: Stokes flow).
    pub convection: bool,
}

/// Residuals of one outer iteration.
#[derive(Debug, Clone, Copy, Default, Serialize)]
pub struct Residuals {
    /// Outer iteration (from 1).
    pub iteration: usize,
    /// Momentum residuals of u, v, w.
    pub momentum: [f64; 3],
    /// Continuity residual of the predicted fluxes.
    pub continuity: f64,
    /// BiCGStab iterations of the three momentum solves.
    pub momentum_iterations: [usize; 3],
    /// CG iterations of the pressure solve.
    pub pressure_iterations: usize,
}

impl Residuals {
    /// The largest residual.
    pub fn max(&self) -> f64 {
        self.momentum
            .iter()
            .copied()
            .fold(self.continuity, f64::max)
    }
}

/// Wall-clock seconds per phase (R6).
#[derive(Debug, Clone, Default, Serialize)]
pub struct Timings {
    /// Boundary values, gradients, limiters.
    pub gradients: f64,
    /// Momentum faces and rows.
    pub momentum_assembly: f64,
    /// Momentum solves.
    pub momentum_solve: f64,
    /// Rhie–Chow fluxes, pressure rows, corrections.
    pub pressure_assembly: f64,
    /// Pressure solves.
    pub pressure_solve: f64,
    /// Multigrid setups.
    pub amg_setup: f64,
}

impl Timings {
    /// Total seconds.
    pub fn total(&self) -> f64 {
        self.gradients
            + self.momentum_assembly
            + self.momentum_solve
            + self.pressure_assembly
            + self.pressure_solve
            + self.amg_setup
    }
}

#[derive(Clone, Copy)]
struct Ptr(*mut f64);
// SAFETY: every loop below writes disjoint entries, one index each.
unsafe impl Send for Ptr {}
// SAFETY: as above.
unsafe impl Sync for Ptr {}
impl Ptr {
    fn of(v: &mut [f64]) -> Self {
        Self(v.as_mut_ptr())
    }
    // methods (not field access) so that closures capture the whole `Ptr`
    fn get(&self) -> *mut f64 {
        self.0
    }
    unsafe fn at(&self, i: usize) -> *mut f64 {
        // SAFETY: guaranteed by the caller (i inside the slice).
        unsafe { self.0.add(i) }
    }
}

/// Runs `f(i)` for `i < n` on the pool, in chunks.
fn par_for(pool: &rayon::ThreadPool, n: usize, f: impl Fn(usize) + Sync) {
    pool.install(|| {
        (0..n.div_ceil(GRAIN)).into_par_iter().for_each(|ch| {
            for i in ch * GRAIN..((ch + 1) * GRAIN).min(n) {
                f(i);
            }
        })
    });
}

/// `out[i] = f(i)`.
fn par_vec(pool: &rayon::ThreadPool, out: &mut [f64], f: impl Fn(usize) -> f64 + Sync) {
    pool.install(|| {
        out.par_chunks_mut(GRAIN).enumerate().for_each(|(c, ch)| {
            for (j, o) in ch.iter_mut().enumerate() {
                *o = f(c * GRAIN + j);
            }
        })
    });
}

/// Reproducible sum of squares (fixed blocks, fixed combination).
fn sum_sq(v: &[f64]) -> f64 {
    let blocks: Vec<f64> = v
        .chunks(cfdcraft_core::repro::BLOCK)
        .map(|b| b.iter().map(|x| x * x).sum())
        .collect();
    cfdcraft_core::repro::combine_block_sums(&blocks)
}

/// The flow solver state on one partition.
pub struct Flow {
    /// Discretisation (geometry, least-squares weights, matrix structure).
    pub disc: Discretisation,
    /// Boundary arrays.
    pub bc: BoundaryArrays,
    /// Physics.
    pub phys: Physics,
    /// Settings.
    pub settings: Settings,
    /// Cell velocity.
    pub u: Soa3,
    /// Cell pressure.
    pub p: Vec<f64>,
    /// Face mass flux (kg/s) out of the owner.
    pub mass_flux: Vec<f64>,
    /// Boundary face velocity.
    pub ub: Soa3,
    /// Boundary face pressure.
    pub pb: Vec<f64>,
    /// Velocity gradients.
    pub grad_u: [Soa3; 3],
    /// Pressure gradient.
    pub grad_p: Soa3,
    /// Limiter of each velocity component.
    pub psi: Soa3,
    /// D_P = V/a_P of the last momentum assembly.
    pub d: Vec<f64>,
    /// D^C_P of the last momentum assembly.
    pub dc: Vec<f64>,
    /// Reference cell of the pressure (no pressure boundary).
    pub reference_cell: Option<usize>,
    /// Phase timings.
    pub timings: Timings,
    /// Multigrid setups done.
    pub amg_setups: usize,
    lsq_p: Soa3,
    mu_face: Vec<f64>,
    eps2: Vec<f64>,
    pool: rayon::ThreadPool,
    be: Box<dyn Backend>,
    iteration: usize,
}

fn disc_with(
    mesh: &Mesh,
    geo: &Geometry,
    points: &std::collections::BTreeMap<String, BoundaryKind>,
) -> Result<Discretisation> {
    let mut prob = Problem::default();
    for p in &mesh.patches {
        let kind = points
            .get(&p.name)
            .copied()
            .unwrap_or(BoundaryKind::FixedValue);
        let b = match kind {
            BoundaryKind::FixedValue => Boundary::value(0.0),
            _ => Boundary::zero_gradient(),
        };
        prob.boundary.insert(p.name.clone(), b);
    }
    Discretisation::new(mesh, geo, &prob)
}

/// The discretisations of the flow solver: the geometry, the matrix
/// structure and the least-squares weights of the velocity (boundary
/// stencil points at the face centroids, projected onto the face normal at
/// outlets and symmetry planes), and the least-squares weights of the
/// pressure (projected at walls and symmetry planes, where the pressure has
/// a zero normal gradient).
pub fn discretisation(
    mesh: &Mesh,
    geo: &Geometry,
    bc: &BoundaryArrays,
) -> Result<(Discretisation, Soa3)> {
    let d = disc_with(mesh, geo, &bc.u_points)?;
    let lsq_p = disc_with(mesh, geo, &bc.p_points)?.lsq;
    Ok((d, lsq_p))
}

impl Flow {
    /// A solver with the fields at `u0`, `p0` and the fluxes of `u0`;
    /// `reference_cell` fixes the pressure level when no boundary does.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        (disc, lsq_p): (Discretisation, Soa3),
        bc: BoundaryArrays,
        phys: Physics,
        settings: Settings,
        u0: [f64; 3],
        p0: f64,
        reference_cell: Option<usize>,
        threads: Option<usize>,
    ) -> Result<Self> {
        settings.validate()?;
        let nc = disc.n_cells;
        let nf = disc.n_faces;
        let nb = nf - disc.n_interior;
        // small meshes run on one thread: the latency of waking the pool
        // would dominate (the result is the same for any number of threads)
        let threads = if nc < SERIAL_BELOW { Some(1) } else { threads };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.unwrap_or(0))
            .build()?;
        let be = cfdcraft_fv::backend::create("cpu", &disc, threads, std::path::Path::new(""))?;
        let mut u = Soa3::zeros(nc);
        u.x.fill(u0[0]);
        u.y.fill(u0[1]);
        u.z.fill(u0[2]);
        let ub = Soa3 {
            x: bc.ub[0].clone(),
            y: bc.ub[1].clone(),
            z: bc.ub[2].clone(),
        };
        let pb: Vec<f64> = (0..nb)
            .map(|fb| {
                if bc.pb_rule[fb] == k::PB_EXTRAPOLATE {
                    p0
                } else {
                    bc.pset[fb]
                }
            })
            .collect();
        let eps2 = match settings.limiter {
            LimiterKind::Venkatakrishnan => disc.venkat_eps2(settings.limiter_k),
            _ => Vec::new(),
        };
        let reference_cell = if bc.has_pressure() {
            None
        } else {
            Some(reference_cell.unwrap_or(0))
        };
        let mut s = Self {
            mu_face: vec![phys.mu; nf],
            lsq_p,
            eps2,
            u,
            p: vec![p0; nc],
            mass_flux: vec![0.0; nf],
            ub,
            pb,
            grad_u: [Soa3::zeros(nc), Soa3::zeros(nc), Soa3::zeros(nc)],
            grad_p: Soa3::zeros(nc),
            psi: Soa3 {
                x: vec![1.0; nc],
                y: vec![1.0; nc],
                z: vec![1.0; nc],
            },
            d: vec![0.0; nc],
            dc: vec![0.0; nc],
            reference_cell,
            timings: Timings::default(),
            amg_setups: 0,
            pool,
            be,
            iteration: 0,
            disc,
            bc,
            phys,
            settings,
        };
        s.update_boundaries();
        s.gradients();
        // start fluxes: the interpolated velocity, without pressure terms
        let relax = vec![0.0; nf];
        let mut m = vec![0.0; nf];
        s.fluxes(&mut m, &relax, false);
        s.mass_flux = m;
        Ok(s)
    }

    /// Outer iterations done.
    pub fn iteration(&self) -> usize {
        self.iteration
    }

    fn geometry(&self) -> FaceGeometry<'_, f64> {
        let d = &self.disc;
        FaceGeometry {
            area: d.face_area.view(),
            face_centroid: d.face_centroid.view(),
            cell_centroid: d.cell_centroid.view(),
            dgeo: &d.dgeo,
            k: d.k.view(),
            weight: &d.weight,
        }
    }

    fn grad_view(&self) -> VelocityGradient<'_, f64> {
        VelocityGradient {
            x: self.grad_u[0].view(),
            y: self.grad_u[1].view(),
            z: self.grad_u[2].view(),
        }
    }

    /// Step 1: boundary values from the current fields.
    pub fn update_boundaries(&mut self) {
        let nb = self.disc.n_faces - self.disc.n_interior;
        let mut ux = vec![0.0; nb];
        let mut uy = vec![0.0; nb];
        let mut uz = vec![0.0; nb];
        let mut pb = vec![0.0; nb];
        {
            let inp = k::BoundaryInputs {
                topo: self.disc.topology(),
                geo: self.geometry(),
                rho: self.phys.rho,
                mass_flux: &self.mass_flux,
            };
            let (u, g, ub0) = (self.u.view(), self.grad_view(), self.ub.view());
            let (px, py, pz) = (Ptr::of(&mut ux), Ptr::of(&mut uy), Ptr::of(&mut uz));
            let bc = &self.bc;
            par_for(&self.pool, nb, |fb| {
                // SAFETY: fb < nb; the arrays were built for this mesh.
                unsafe {
                    let v = k::boundary_velocity(fb, &inp, bc.ub_rule[fb], ub0.at(fb), &u, &g);
                    *px.at(fb) = v[0];
                    *py.at(fb) = v[1];
                    *pz.at(fb) = v[2];
                }
            });
            let (p, gp) = (&self.p, self.grad_p.view());
            par_vec(&self.pool, &mut pb, |fb| {
                // SAFETY: as above.
                unsafe { k::boundary_pressure(fb, &inp, bc.pb_rule[fb], bc.pset[fb], p, &gp) }
            });
        }
        self.ub = Soa3 {
            x: ux,
            y: uy,
            z: uz,
        };
        self.pb = pb;
    }

    /// Least-squares gradient of `phi` with the weights `w`; boundary face
    /// values a φ_P + b with `a` = `zero_gradient` (1: zero normal gradient
    /// at a projected stencil point) and b = `value` where a = 0.
    fn gradient_of(
        &self,
        w: &Soa3,
        phi: &[f64],
        zero_gradient: &[f64],
        value: &[f64],
        out: &mut Soa3,
    ) {
        let topo = self.disc.topology();
        let w = w.view();
        let b: Vec<f64> = zero_gradient
            .iter()
            .zip(value)
            .map(|(&a, &v)| if a == 1.0 { 0.0 } else { v })
            .collect();
        let bc = fv::Boundary {
            a: zero_gradient,
            b: &b,
        };
        let (px, py, pz) = (
            Ptr::of(&mut out.x),
            Ptr::of(&mut out.y),
            Ptr::of(&mut out.z),
        );
        par_for(&self.pool, self.disc.n_cells, |c| {
            // SAFETY: c is a local cell; the weights were built for this mesh.
            unsafe {
                let g = fv::lsq_gradient(c, &topo, &w, phi, &bc);
                *px.at(c) = g[0];
                *py.at(c) = g[1];
                *pz.at(c) = g[2];
            }
        });
    }

    /// Step 2: gradients and limiters.
    pub fn gradients(&mut self) {
        let mut g = std::mem::take(&mut self.grad_u);
        let (lu, au) = (&self.disc.lsq, &self.bc.u_zero_gradient);
        self.gradient_of(lu, &self.u.x, au, &self.ub.x, &mut g[0]);
        self.gradient_of(lu, &self.u.y, au, &self.ub.y, &mut g[1]);
        self.gradient_of(lu, &self.u.z, au, &self.ub.z, &mut g[2]);
        self.grad_u = g;
        let mut gp = std::mem::take(&mut self.grad_p);
        self.gradient_of(
            &self.lsq_p,
            &self.p,
            &self.bc.p_zero_gradient,
            &self.pb,
            &mut gp,
        );
        self.grad_p = gp;
        if self.settings.scheme == Scheme::Upwind || self.settings.limiter == LimiterKind::None {
            return;
        }
        let mut psi = std::mem::take(&mut self.psi);
        for (i, out) in [&mut psi.x, &mut psi.y, &mut psi.z].into_iter().enumerate() {
            let (phi, b) = match i {
                0 => (&self.u.x, &self.ub.x),
                1 => (&self.u.y, &self.ub.y),
                _ => (&self.u.z, &self.ub.z),
            };
            let eps2 = match self.settings.limiter {
                LimiterKind::VenkatakrishnanWang => {
                    let mn = phi.iter().copied().fold(f64::INFINITY, f64::min);
                    let mx = phi.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                    let e = self.settings.limiter_k * (mx - mn);
                    vec![e * e]
                }
                _ => self.eps2.clone(),
            };
            let topo = self.disc.topology();
            let bm: Vec<f64> = self
                .bc
                .u_zero_gradient
                .iter()
                .zip(b)
                .map(|(&a, &v)| if a == 1.0 { 0.0 } else { v })
                .collect();
            let bc = fv::Boundary {
                a: &self.bc.u_zero_gradient,
                b: &bm,
            };
            let grad = self.grad_u[i].view();
            let (cc, fc) = (
                self.disc.cell_centroid.view(),
                self.disc.face_centroid.view(),
            );
            par_vec(&self.pool, out, |c| {
                // SAFETY: as in `gradient_of`.
                unsafe { fv::venkat_limiter(c, &topo, phi, &bc, &grad, &cc, &fc, &eps2) }
            });
        }
        self.psi = psi;
    }

    /// Face mass fluxes into `out`; with `pressure` the Rhie–Chow terms with
    /// the current D, otherwise the plain interpolated velocity flux.
    fn fluxes(&self, out: &mut [f64], relax: &[f64], pressure: bool) {
        let zero = vec![0.0; if pressure { 0 } else { self.disc.n_cells }];
        let inp = k::FluxInputs {
            topo: self.disc.topology(),
            geo: self.geometry(),
            rho: self.phys.rho,
            u: self.u.view(),
            grad: self.grad_view(),
            p: &self.p,
            grad_p: self.grad_p.view(),
            pb: &self.pb,
            ub: self.ub.view(),
            flux_kind: &self.bc.flux_kind,
            mom_kind: &self.bc.mom_kind,
            d: if pressure { &self.d } else { &zero },
            relax,
        };
        par_vec(&self.pool, out, |f| {
            // SAFETY: f is a face; all arrays were built for this mesh.
            unsafe { k::face_mass_flux(f, &inp) }
        });
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_linear(
        &mut self,
        val: &[f64],
        rhs: &[f64],
        x: &mut [f64],
        pc: Preconditioner,
        cg: bool,
        tol: f64,
        max_iter: usize,
    ) -> Result<(usize, f64)> {
        let be = self.be.as_mut();
        be.set_matrix(val)?;
        be.upload(V::Rhs, rhs)?;
        be.upload(V::Phi, x)?;
        be.residual(V::Phi, V::Rhs, V::R)?;
        let r0 = be.dot(V::R, V::R)?.sqrt();
        let it = if r0 == 0.0 {
            0
        } else if cg {
            linear::cg(be, pc, max_iter, tol, r0)?.0
        } else {
            linear::bicgstab(be, pc, max_iter, tol, r0)?
        };
        x.copy_from_slice(&be.download(V::Phi)?);
        Ok((it, r0))
    }

    /// One outer iteration; returns the residuals before the solves.
    pub fn iterate(&mut self) -> Result<Residuals> {
        self.iteration += 1;
        let mut res = Residuals {
            iteration: self.iteration,
            ..Default::default()
        };
        let nc = self.disc.n_cells;
        let nf = self.disc.n_faces;
        let s = self.settings.clone();
        let alpha = s.relaxation_velocity;
        let t = Instant::now();
        self.update_boundaries();
        self.gradients();
        self.timings.gradients += t.elapsed().as_secs_f64();

        // step 3: relaxation term of the Rhie–Chow flux
        let t = Instant::now();
        let mut relax = vec![0.0; nf];
        {
            let zero = vec![0.0; nf];
            self.fluxes(&mut relax, &zero, false);
            let m = &self.mass_flux;
            let r = Ptr::of(&mut relax);
            par_for(&self.pool, nf, |f| {
                // SAFETY: f < nf.
                unsafe { *r.at(f) = (1.0 - alpha) * (m[f] - *r.at(f)) }
            });
        }

        // step 4: momentum
        let mut cp = vec![0.0; nf];
        let mut cn = vec![0.0; nf];
        let mut bf = Soa3::zeros(nf);
        {
            let inp = k::MomentumInputs {
                topo: self.disc.topology(),
                geo: self.geometry(),
                mass_flux: &self.mass_flux,
                mu: &self.mu_face,
                u: self.u.view(),
                grad: self.grad_view(),
                psi: self.psi.view(),
                ub: self.ub.view(),
                mom_kind: &self.bc.mom_kind,
                blend: if s.scheme == Scheme::Upwind { 0.0 } else { 1.0 },
                nonorth: if s.nonorth_correction { 1.0 } else { 0.0 },
                convection: if self.phys.convection { 1.0 } else { 0.0 },
            };
            let (pcp, pcn) = (Ptr::of(&mut cp), Ptr::of(&mut cn));
            let (bx, by, bz) = (Ptr::of(&mut bf.x), Ptr::of(&mut bf.y), Ptr::of(&mut bf.z));
            par_for(&self.pool, nf, |f| {
                // SAFETY: f is a face of the mesh.
                unsafe {
                    let (a, b, c) = k::momentum_face(f, &inp);
                    *pcp.at(f) = a;
                    *pcn.at(f) = b;
                    *bx.at(f) = c[0];
                    *by.at(f) = c[1];
                    *bz.at(f) = c[2];
                }
            });
        }
        let stored = self.disc.sell.stored();
        let mut val = vec![0.0; stored];
        let mut rhs = Soa3::zeros(nc);
        let mut diag = Soa3::zeros(nc);
        {
            let layout = SellLayout {
                chunk: self.disc.sell.chunk,
                slice_ptr: &self.disc.sell.slice_ptr,
            };
            let inp = k::MomentumRowInputs {
                topo: self.disc.topology(),
                cp: &cp,
                cn: &cn,
                bf: bf.view(),
                grad_p: self.grad_p.view(),
                body_force: self.phys.body_force,
                volume: &self.disc.cell_volume,
                u: self.u.view(),
                alpha,
                layout,
                mom_kind: &self.bc.mom_kind,
                face_area: self.disc.face_area.view(),
            };
            let out = (
                Ptr::of(&mut val),
                Ptr::of(&mut rhs.x),
                Ptr::of(&mut rhs.y),
                Ptr::of(&mut rhs.z),
                Ptr::of(&mut self.d),
                Ptr::of(&mut self.dc),
                Ptr::of(&mut diag.x),
                Ptr::of(&mut diag.y),
                Ptr::of(&mut diag.z),
            );
            par_for(&self.pool, nc, |c| {
                let o = k::MomentumRowOut {
                    val: out.0.get(),
                    rhs: [out.1.get(), out.2.get(), out.3.get()],
                    d: out.4.get(),
                    dc: out.5.get(),
                    diag: [out.6.get(), out.7.get(), out.8.get()],
                };
                // SAFETY: c is a local cell; row c is written only here.
                unsafe { k::momentum_row(c, &inp, &o) }
            });
        }
        if s.algorithm == Algorithm::Simple {
            self.dc.copy_from_slice(&self.d);
        }
        // residual scale ‖a_P u‖ + ‖b‖ over the three components (unrelaxed
        // a_P and b), so that it is finite also from a fluid at rest
        let scale = {
            let mut v = vec![0.0; 3 * nc];
            let mut w = vec![0.0; 3 * nc];
            for c in 0..nc {
                let u = self.u.get(c);
                let b = [rhs.x[c], rhs.y[c], rhs.z[c]];
                let dg = [diag.x[c], diag.y[c], diag.z[c]];
                for i in 0..3 {
                    let dr = dg[i];
                    let a = dr * alpha;
                    v[3 * c + i] = a * u[i];
                    w[3 * c + i] = b[i] - (dr - a) * u[i];
                }
            }
            (sum_sq(&v).sqrt() + sum_sq(&w).sqrt()).max(f64::MIN_POSITIVE)
        };
        self.timings.momentum_assembly += t.elapsed().as_secs_f64();
        let t = Instant::now();
        let mut u = std::mem::take(&mut self.u);
        let sym = self.bc.mom_kind.contains(&k::MOM_SYMMETRY);
        for (i, ((x, b), dg)) in [(&mut u.x, &rhs.x), (&mut u.y, &rhs.y), (&mut u.z, &rhs.z)]
            .into_iter()
            .zip([&diag.x, &diag.y, &diag.z])
            .enumerate()
        {
            if sym {
                // the diagonal of this component (symmetry planes)
                let sell = &self.disc.sell;
                for c in 0..nc {
                    val[sell.pos(c, 0)] = dg[c];
                }
            }
            let (it, r0) = self.solve_linear(
                &val,
                b,
                x,
                Preconditioner::Jacobi,
                false,
                s.momentum_tolerance,
                s.momentum_max_iter,
            )?;
            res.momentum[i] = r0 / scale;
            res.momentum_iterations[i] = it;
        }
        self.u = u;
        self.timings.momentum_solve += t.elapsed().as_secs_f64();

        // step 5: predicted fluxes
        let t = Instant::now();
        let mut m = vec![0.0; nf];
        self.fluxes(&mut m, &relax, true);
        let topo = self.disc.topology();
        let mut defect = vec![0.0; nc];
        let mut through = vec![0.0; nc];
        // SAFETY: c is a local cell, m has one entry per face.
        par_vec(&self.pool, &mut defect, |c| unsafe {
            k::net_outflow(c, &topo, &m)
        });
        // SAFETY: as above.
        par_vec(&self.pool, &mut through, |c| unsafe {
            k::throughflow(c, &topo, &m)
        });
        res.continuity = sum_sq(&defect).sqrt() / sum_sq(&through).sqrt().max(f64::MIN_POSITIVE);

        // step 6: pressure correction
        let mut cf = vec![0.0; nf];
        {
            let geo = self.geometry();
            let (fk, dc, rho) = (&self.bc.flux_kind, &self.dc, self.phys.rho);
            // SAFETY: f is a face.
            par_vec(&self.pool, &mut cf, |f| unsafe {
                k::pressure_face(f, &topo, &geo, fk, rho, dc)
            });
        }
        let mut pval = vec![0.0; stored];
        let mut prhs = vec![0.0; nc];
        {
            let layout = SellLayout {
                chunk: self.disc.sell.chunk,
                slice_ptr: &self.disc.sell.slice_ptr,
            };
            let (pv, pr) = (Ptr::of(&mut pval), Ptr::of(&mut prhs));
            let (cf, m) = (&cf, &m);
            par_for(&self.pool, nc, |c| {
                // SAFETY: c is a local cell; row c is written only here.
                unsafe { k::pressure_row(c, &topo, cf, m, &layout, pv.get(), pr.get()) }
            });
        }
        if let Some(r) = self.reference_cell {
            let pos = self.disc.sell.pos(r, 0);
            pval[pos] *= 2.0;
        }
        self.timings.pressure_assembly += t.elapsed().as_secs_f64();
        let t = Instant::now();
        let use_amg = s.pressure_preconditioner == Preconditioner::Amg;
        if use_amg && (self.amg_setups == 0 || (self.iteration - 1).is_multiple_of(s.amg_rebuild)) {
            self.be.set_matrix(&pval)?;
            linear::setup_amg(self.be.as_mut(), &s.amg)?;
            self.amg_setups += 1;
            self.timings.amg_setup += t.elapsed().as_secs_f64();
        }
        let t = Instant::now();
        let mut pc = vec![0.0; nc];
        let (it, _) = self.solve_linear(
            &pval,
            &prhs,
            &mut pc,
            s.pressure_preconditioner,
            true,
            s.pressure_tolerance,
            s.pressure_max_iter,
        )?;
        res.pressure_iterations = it;
        self.timings.pressure_solve += t.elapsed().as_secs_f64();

        // step 7: corrections
        let t = Instant::now();
        let topo = self.disc.topology();
        if let Some(r) = self.reference_cell {
            let shift = pc[r];
            pc.iter_mut().for_each(|v| *v -= shift);
        }
        let mut mc = vec![0.0; nf];
        {
            let (cf, pcr, m) = (&cf, &pc, &m);
            // SAFETY: f is a face; cf, m have one entry per face, pc per cell.
            par_vec(&self.pool, &mut mc, |f| unsafe {
                k::correct_flux(f, &topo, cf, m, pcr)
            });
        }
        // gradient of p': zero on pressure boundaries, zero normal gradient
        // elsewhere
        let nb = nf - self.disc.n_interior;
        let pca: Vec<f64> = (0..nb)
            .map(|fb| {
                if self.bc.flux_kind[fb] == k::FLUX_PRESSURE {
                    0.0
                } else {
                    1.0
                }
            })
            .collect();
        let mut gpc = Soa3::zeros(nc);
        self.gradient_of(&self.lsq_p, &pc, &pca, &vec![0.0; nb], &mut gpc);
        let alpha_p = s.alpha_p();
        {
            let dc = &self.dc;
            let (ux, uy, uz) = (
                Ptr::of(&mut self.u.x),
                Ptr::of(&mut self.u.y),
                Ptr::of(&mut self.u.z),
            );
            let pp = Ptr::of(&mut self.p);
            let (g, pc) = (&gpc, &pc);
            par_for(&self.pool, nc, |c| {
                // SAFETY: c is a local cell, written once.
                unsafe {
                    *ux.at(c) -= dc[c] * g.x[c];
                    *uy.at(c) -= dc[c] * g.y[c];
                    *uz.at(c) -= dc[c] * g.z[c];
                    *pp.at(c) += alpha_p * pc[c];
                }
            });
        }
        self.timings.pressure_assembly += t.elapsed().as_secs_f64();
        anyhow::ensure!(
            res.max().is_finite(),
            "the solution diverged in outer iteration {} (residuals {:?}, continuity {:.3e}); \
             lower [solver] relaxation_velocity or start from upwind",
            self.iteration,
            res.momentum,
            res.continuity
        );
        self.mass_flux = mc;
        Ok(res)
    }
}

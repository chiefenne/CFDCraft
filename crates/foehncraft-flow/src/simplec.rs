//! The incompressible flow solver state and the steady algorithm: SIMPLEC
//! (Van Doormaal & Raithby 1984), or SIMPLE (Patankar & Spalding 1972), on
//! the collocated, cell-centred grid with the Rhie–Chow face mass flux (D19,
//! D27). The transient PISO algorithm on the same state is in
//! [`crate::transient`].
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
//! The algorithm is a sequence of operations of a [`FlowDevice`]: the CPU
//! device runs the single-source bodies of `foehncraft_kernels::flow` in
//! threaded loops, the CUDA device the same bodies on the GPU. Each output
//! is written by one index, so the result does not depend on the number of
//! threads (D31).
//!
//! Residuals (dimensionless, per outer iteration, before the solves):
//! momentum ‖b − A u_i‖₂ / ‖a_P u‖₂ per component with ‖a_P u‖ over all
//! three components, continuity ‖Σ_f ±ṁ*_f‖₂ / ‖Σ_f |ṁ_f|‖₂ (net outflow
//! of the predicted fluxes over the cell throughflow).

use crate::bc::BoundaryArrays;
use crate::device::{
    self, BoundaryArgs, Buf, Buf3, Eps2, FlowDevice, FluxArgs, Lsq, Mask, MomentumFaceArgs,
    MomentumRowArgs, Setup, TimeTerms,
};
use anyhow::Result;
use foehncraft_fv::amg::AmgSettings;
use foehncraft_fv::disc::{Discretisation, Soa3};
use foehncraft_fv::problem::{Boundary, BoundaryKind, Problem};
use foehncraft_fv::solver::{LimiterKind, Preconditioner, Scheme};
use foehncraft_kernels::flow as k;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::time::Instant;

pub use crate::device::cpu::SERIAL_BELOW;

/// Pressure–velocity algorithm.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Algorithm {
    /// SIMPLEC (Van Doormaal & Raithby 1984), steady.
    #[default]
    Simplec,
    /// SIMPLE (Patankar & Spalding 1972), steady.
    Simple,
    /// PISO (Issa 1986), transient: per time step a momentum predictor and
    /// `correctors` pressure corrections, optionally repeated in
    /// `outer_iterations` (a PIMPLE-type loop).
    Piso,
}

/// `[solver]` settings of a flow case.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Settings {
    /// Pressure–velocity algorithm (default SIMPLEC, PISO with `[time]`).
    pub algorithm: Option<Algorithm>,
    /// Convection scheme of the momentum equations.
    pub scheme: Scheme,
    /// `scheme = "central"`: fraction of the bounded central scheme on the
    /// interior faces, the rest linear upwind (0 to 1).
    pub central_blending: f64,
    /// `scheme = "central"`: β of the normalised-variable bounding; the
    /// scheme moves to upwind where the normalised upwind value leaves
    /// [β, 1) (0: plain central, unbounded).
    pub central_bounding: f64,
    /// Limiter of the velocity gradients.
    pub limiter: LimiterKind,
    /// Constant K of the limiter.
    pub limiter_k: f64,
    /// Non-orthogonal correction of the viscous term (D26).
    pub nonorth_correction: bool,
    /// Under-relaxation of the momentum equations α_u (default 0.7 steady,
    /// 1 transient).
    pub relaxation_velocity: Option<f64>,
    /// Relaxation of the pressure correction α_p (default 1 for SIMPLEC and
    /// PISO, 0.3 for SIMPLE).
    pub relaxation_pressure: Option<f64>,
    /// Largest number of outer iterations (steady).
    pub max_iterations: usize,
    /// Convergence: every residual below this value (steady).
    pub tolerance: f64,
    /// Relative tolerance of each momentum solve (default 1e-2 steady, 1e-8
    /// transient).
    pub momentum_tolerance: Option<f64>,
    /// Largest number of BiCGStab iterations per momentum solve.
    pub momentum_max_iter: usize,
    /// Relative tolerance of each pressure-correction solve (default 1e-2
    /// steady, 1e-8 transient).
    pub pressure_tolerance: Option<f64>,
    /// Largest number of CG iterations per pressure solve.
    pub pressure_max_iter: usize,
    /// Preconditioner of the pressure solve.
    pub pressure_preconditioner: Preconditioner,
    /// Multigrid settings.
    pub amg: AmgSettings,
    /// Rebuild the multigrid hierarchy every N outer iterations (steady) or
    /// N time steps (transient).
    pub amg_rebuild: usize,
    /// Log a line every N outer iterations or time steps.
    pub log_every: usize,
    /// PISO: pressure corrections per outer iteration.
    pub correctors: usize,
    /// PISO: outer iterations per time step (1: plain PISO).
    pub outer_iterations: usize,
    /// PISO: an outer iteration ends the time step when the momentum and
    /// continuity residuals of its start are below this value.
    pub outer_tolerance: f64,
    /// PISO: in the first outer iteration of a second-order step the
    /// explicit quantities are extrapolated from the two old levels to the
    /// new time: the convecting mass flux and the velocity gradients and
    /// limiters of the deferred corrections. Without it the lag of these
    /// terms makes plain PISO (one outer iteration) first order in time.
    pub extrapolation: bool,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            algorithm: None,
            scheme: Scheme::LinearUpwind,
            limiter: LimiterKind::None,
            limiter_k: 5.0,
            central_blending: 1.0,
            central_bounding: 0.1,
            nonorth_correction: true,
            relaxation_velocity: None,
            relaxation_pressure: None,
            max_iterations: 2000,
            tolerance: 1e-8,
            momentum_tolerance: None,
            momentum_max_iter: 50,
            pressure_tolerance: None,
            pressure_max_iter: 100,
            pressure_preconditioner: Preconditioner::Amg,
            amg: AmgSettings::default(),
            amg_rebuild: 10,
            log_every: 10,
            correctors: 2,
            outer_iterations: 1,
            outer_tolerance: 1e-6,
            extrapolation: true,
        }
    }
}

impl Settings {
    /// The algorithm in force: the one given, else PISO in a transient case
    /// and SIMPLEC in a steady one.
    pub fn algorithm(&self, transient: bool) -> Algorithm {
        self.algorithm.unwrap_or(if transient {
            Algorithm::Piso
        } else {
            Algorithm::Simplec
        })
    }

    /// Relative tolerance of the momentum solves. A time step needs
    /// accurate solves (their error is not iterated away), a steady outer
    /// iteration does not.
    pub fn momentum_tol(&self, transient: bool) -> f64 {
        self.momentum_tolerance
            .unwrap_or(if transient { 1e-8 } else { 1e-2 })
    }

    /// Relative tolerance of the pressure solves (see
    /// [`momentum_tol`](Self::momentum_tol)).
    pub fn pressure_tol(&self, transient: bool) -> f64 {
        self.pressure_tolerance
            .unwrap_or(if transient { 1e-8 } else { 1e-2 })
    }

    /// α_u.
    pub fn alpha_u(&self, transient: bool) -> f64 {
        self.relaxation_velocity
            .unwrap_or(if transient { 1.0 } else { 0.7 })
    }

    /// α_p.
    pub fn alpha_p(&self, transient: bool) -> f64 {
        self.relaxation_pressure
            .unwrap_or(match self.algorithm(transient) {
                Algorithm::Simplec | Algorithm::Piso => 1.0,
                Algorithm::Simple => 0.3,
            })
    }

    /// Checks the ranges; `transient` when the case has a `[time]` table.
    pub fn validate(&self, transient: bool) -> Result<()> {
        let alg = self.algorithm(transient);
        let au = self.alpha_u(transient);
        crate::need!(
            au > 0.0 && au <= 1.0,
            CASE_RANGE,
            "[solver] relaxation_velocity must lie in (0, 1]"
        );
        let ap = self.alpha_p(transient);
        crate::need!(
            ap > 0.0 && ap <= 1.0,
            CASE_RANGE,
            "[solver] relaxation_pressure must lie in (0, 1]"
        );
        if transient {
            crate::need!(
                alg == Algorithm::Piso,
                CASE_INCONSISTENT,
                "[solver] algorithm = \"{}\" is a steady algorithm; a case with a [time] table \
                 runs PISO (remove `algorithm` or set it to \"piso\")",
                match alg {
                    Algorithm::Simple => "simple",
                    _ => "simplec",
                }
            );
            crate::need!(
                self.correctors >= 1,
                CASE_RANGE,
                "[solver] correctors must be >= 1"
            );
            crate::need!(
                self.outer_iterations >= 1,
                CASE_RANGE,
                "[solver] outer_iterations must be >= 1"
            );
            crate::need!(
                self.outer_tolerance > 0.0,
                CASE_RANGE,
                "[solver] outer_tolerance must be > 0"
            );
            crate::need!(
                !(self.outer_iterations == 1 && au < 1.0),
                CASE_INCONSISTENT,
                "[solver] relaxation_velocity < 1 needs outer_iterations > 1 in a transient case \
                 (with one outer iteration the relaxation would lag the solution)"
            );
        } else {
            crate::need!(
                alg != Algorithm::Piso,
                CASE_INCONSISTENT,
                "[solver] algorithm = \"piso\" is transient and needs a [time] table"
            );
            crate::diag::check_fix(
                !(alg == Algorithm::Simplec && au >= 1.0),
                crate::diag::CASE_INCONSISTENT,
                &["solver.relaxation_velocity=0.7"],
                || {
                    String::from(
                        "[solver] SIMPLEC needs relaxation_velocity < 1 for a steady case \
                 (D^C = V/(a_P/α − Σ|a_N|) is singular at α = 1); try 0.7",
                    )
                },
            )?;
        }
        crate::need!(
            self.tolerance > 0.0,
            CASE_RANGE,
            "[solver] tolerance must be > 0"
        );
        crate::need!(
            self.limiter_k >= 0.0,
            CASE_RANGE,
            "[solver] limiter_k must be >= 0"
        );
        crate::need!(
            (0.0..=1.0).contains(&self.central_blending),
            CASE_RANGE,
            "[solver] central_blending must lie in [0, 1]"
        );
        crate::need!(
            (0.0..1.0).contains(&self.central_bounding),
            CASE_RANGE,
            "[solver] central_bounding must lie in [0, 1)"
        );
        for (name, t) in [
            ("momentum_tolerance", self.momentum_tol(transient)),
            ("pressure_tolerance", self.pressure_tol(transient)),
        ] {
            crate::need!(
                t > 0.0 && t < 1.0,
                CASE_RANGE,
                "[solver] {name} must lie in (0, 1)"
            );
        }
        crate::need!(
            self.amg_rebuild >= 1,
            CASE_RANGE,
            "[solver] amg_rebuild must be >= 1"
        );
        crate::need!(
            self.log_every >= 1,
            CASE_RANGE,
            "[solver] log_every must be >= 1"
        );
        self.amg.validate().map_err(|e| {
            let m = e.to_string();
            let mut d = crate::diag::Diagnostic::new(crate::diag::CASE_RANGE, m.clone());
            d.key = crate::diag::key_of(&m).map(|k| format!("solver.{k}"));
            crate::diag::err(d)
        })
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
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq)]
pub struct Residuals {
    /// Outer iteration (from 1).
    pub iteration: usize,
    /// Momentum residuals of u, v, w.
    pub momentum: [f64; 3],
    /// Continuity residual of the predicted fluxes.
    pub continuity: f64,
    /// BiCGStab iterations of the three momentum solves.
    pub momentum_iterations: [usize; 3],
    /// CG iterations of the pressure solve (all correctors).
    pub pressure_iterations: usize,
}

impl Residuals {
    /// Every residual is finite (a NaN does not hide behind
    /// [`max`](Self::max), which skips it).
    pub fn finite(&self) -> bool {
        self.continuity.is_finite() && self.momentum.iter().all(|r| r.is_finite())
    }

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

/// Host copy of the solution fields.
#[derive(Debug, Clone)]
pub struct State {
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
    /// Limiters of the velocity components.
    pub psi: Soa3,
}

/// A host copy of the fields with the discretisation and the physics, for
/// the output.
pub struct Solution<'a> {
    /// The fields.
    pub state: State,
    /// Discretisation.
    pub disc: &'a Discretisation,
    /// Physics.
    pub phys: Physics,
}

impl std::ops::Deref for Solution<'_> {
    type Target = State;
    fn deref(&self) -> &State {
        &self.state
    }
}

/// The field buffers that carry the solution from one iteration (or time
/// step) to the next.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Fields {
    pub u: Buf3,
    pub p: Buf,
    pub mass_flux: Buf,
    pub ub: Buf3,
    pub ub_set: Buf3,
    pub pb: Buf,
    pub grad_u: [Buf3; 3],
    pub grad_p: Buf3,
    pub psi: Buf3,
}

/// Work buffers of one iteration.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Work {
    pub zero_face: Buf,
    pub relax: Buf,
    pub cp: Buf,
    pub cn: Buf,
    pub bf: Buf3,
    pub val: Buf,
    pub rhs: Buf3,
    pub diag: Buf3,
    pub d: Buf,
    pub ds: Buf,
    pub dc: Buf,
    pub m: Buf,
    pub defect: Buf,
    pub through: Buf,
    pub cf: Buf,
    pub pval: Buf,
    pub prhs: Buf,
    pub pc: Buf,
    pub mc: Buf,
    pub gpc: Buf3,
    pub zero_boundary: Buf,
}

/// The flow solver state on one partition.
pub struct Flow {
    /// Host copy of the discretisation (geometry for the output).
    pub disc: Discretisation,
    /// Boundary arrays.
    pub bc: BoundaryArrays,
    /// Physics.
    pub phys: Physics,
    /// Settings.
    pub settings: Settings,
    /// Reference cell of the pressure (no pressure boundary).
    pub reference_cell: Option<usize>,
    /// Phase timings.
    pub timings: Timings,
    /// Multigrid setups done.
    pub amg_setups: usize,
    /// The pressure matrix of the last multigrid setup (for a restart that
    /// continues bit-identically).
    pub amg_matrix: Option<Vec<f64>>,
    pub(crate) dev: Box<dyn FlowDevice>,
    pub(crate) f: Fields,
    pub(crate) w: Work,
    pub(crate) transient: bool,
    pub(crate) iteration: usize,
    /// Device copies of the last valid u, p and ṁ and its iteration, for
    /// the divergence report (D63).
    valid: Option<(Buf3, Buf, Buf, usize)>,
    /// The residuals of the outer iteration or corrector that diverged.
    pub(crate) diverged: Option<Residuals>,
    /// The velocity of the central correction when it is not the current
    /// one (the extrapolated velocity of a time step, set by
    /// [`crate::transient`]).
    pub(crate) central_u: Option<Buf3>,
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

/// Where a flow runs.
#[derive(Debug, Clone)]
pub struct Target {
    /// Device name: `cpu` or `cuda`.
    pub device: String,
    /// CPU threads.
    pub threads: Option<usize>,
    /// PTX module of the CUDA device.
    pub ptx: std::path::PathBuf,
}

impl Default for Target {
    fn default() -> Self {
        Self {
            device: "cpu".into(),
            threads: None,
            ptx: std::path::PathBuf::new(),
        }
    }
}

impl Flow {
    /// A solver with the fields at `u0`, `p0` and the fluxes of `u0`;
    /// `reference_cell` fixes the pressure level when no boundary does.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        disc: (Discretisation, Soa3),
        bc: BoundaryArrays,
        phys: Physics,
        settings: Settings,
        u0: [f64; 3],
        p0: f64,
        reference_cell: Option<usize>,
        target: &Target,
    ) -> Result<Self> {
        Self::on(disc, bc, phys, settings, false, reference_cell, target).and_then(|mut s| {
            let nc = s.disc.n_cells;
            s.initialise(
                &Soa3 {
                    x: vec![u0[0]; nc],
                    y: vec![u0[1]; nc],
                    z: vec![u0[2]; nc],
                },
                &vec![p0; nc],
            )?;
            Ok(s)
        })
    }

    /// A solver on `target` with zero fields; set them with
    /// [`initialise`](Self::initialise) or [`load`](Self::load).
    pub fn on(
        (disc, lsq_p): (Discretisation, Soa3),
        bc: BoundaryArrays,
        phys: Physics,
        settings: Settings,
        transient: bool,
        reference_cell: Option<usize>,
        target: &Target,
    ) -> Result<Self> {
        settings.validate(transient)?;
        let nc = disc.n_cells;
        let nf = disc.n_faces;
        let nb = nf - disc.n_interior;
        let eps2 = match settings.limiter {
            LimiterKind::Venkatakrishnan => disc.venkat_eps2(settings.limiter_k),
            _ => Vec::new(),
        };
        let mut dev = device::create(
            &target.device,
            &Setup {
                disc: &disc,
                lsq_p: &lsq_p,
                bc: &bc,
                rho: phys.rho,
                mu: phys.mu,
                eps2: &eps2,
                threads: target.threads,
            },
            &target.ptx,
        )?;
        let stored = disc.sell.stored();
        let cell3 = |dev: &mut Box<dyn FlowDevice>| -> Result<Buf3> {
            Ok(Buf3 {
                x: dev.alloc(nc)?,
                y: dev.alloc(nc)?,
                z: dev.alloc(nc)?,
            })
        };
        let u = cell3(&mut dev)?;
        let grad_u = [cell3(&mut dev)?, cell3(&mut dev)?, cell3(&mut dev)?];
        let grad_p = cell3(&mut dev)?;
        let psi = cell3(&mut dev)?;
        for b in psi.all() {
            dev.fill(b, 1.0)?;
        }
        let ub_set = Buf3 {
            x: dev.alloc(nb)?,
            y: dev.alloc(nb)?,
            z: dev.alloc(nb)?,
        };
        dev.upload(ub_set.x, &bc.ub[0])?;
        dev.upload(ub_set.y, &bc.ub[1])?;
        dev.upload(ub_set.z, &bc.ub[2])?;
        let ub = Buf3 {
            x: dev.alloc(nb)?,
            y: dev.alloc(nb)?,
            z: dev.alloc(nb)?,
        };
        let f = Fields {
            u,
            p: dev.alloc(nc)?,
            mass_flux: dev.alloc(nf)?,
            ub,
            ub_set,
            pb: dev.alloc(nb)?,
            grad_u,
            grad_p,
            psi,
        };
        let w = Work {
            zero_face: dev.alloc(nf)?,
            relax: dev.alloc(nf)?,
            cp: dev.alloc(nf)?,
            cn: dev.alloc(nf)?,
            bf: Buf3 {
                x: dev.alloc(nf)?,
                y: dev.alloc(nf)?,
                z: dev.alloc(nf)?,
            },
            val: dev.alloc(stored)?,
            rhs: cell3(&mut dev)?,
            diag: cell3(&mut dev)?,
            d: dev.alloc(nc)?,
            ds: dev.alloc(nc)?,
            dc: dev.alloc(nc)?,
            m: dev.alloc(nf)?,
            defect: dev.alloc(nc)?,
            through: dev.alloc(nc)?,
            cf: dev.alloc(nf)?,
            pval: dev.alloc(stored)?,
            prhs: dev.alloc(nc)?,
            pc: dev.alloc(nc)?,
            mc: dev.alloc(nf)?,
            gpc: cell3(&mut dev)?,
            zero_boundary: dev.alloc(nb)?,
        };
        let reference_cell = if bc.has_pressure() {
            None
        } else {
            Some(reference_cell.unwrap_or(0))
        };
        Ok(Self {
            disc,
            bc,
            phys,
            settings,
            reference_cell,
            timings: Timings::default(),
            amg_setups: 0,
            amg_matrix: None,
            dev,
            f,
            w,
            transient,
            iteration: 0,
            valid: None,
            diverged: None,
            central_u: None,
        })
    }

    /// Keeps a device copy of the current u, p and ṁ as the last valid
    /// state (D63): three cell fields and one face field, copied on the
    /// device, against the dozens of passes of an iteration.
    pub fn keep_valid(&mut self) -> Result<()> {
        let (nc, nf) = (self.disc.n_cells, self.disc.n_faces);
        let v = match self.valid {
            Some(v) => v,
            None => {
                let u = Buf3 {
                    x: self.dev.alloc(nc)?,
                    y: self.dev.alloc(nc)?,
                    z: self.dev.alloc(nc)?,
                };
                (u, self.dev.alloc(nc)?, self.dev.alloc(nf)?, 0)
            }
        };
        for (a, b) in self.f.u.all().into_iter().zip(v.0.all()) {
            self.dev.copy(a, b)?;
        }
        self.dev.copy(self.f.p, v.1)?;
        self.dev.copy(self.f.mass_flux, v.2)?;
        self.valid = Some((v.0, v.1, v.2, self.iteration));
        Ok(())
    }

    /// Puts the last valid state of [`keep_valid`](Self::keep_valid) back
    /// with its boundary values and gradients; false when none was kept.
    pub fn rewind(&mut self) -> Result<bool> {
        let Some((u, p, m, it)) = self.valid else {
            return Ok(false);
        };
        for (a, b) in u.all().into_iter().zip(self.f.u.all()) {
            self.dev.copy(a, b)?;
        }
        self.dev.copy(p, self.f.p)?;
        self.dev.copy(m, self.f.mass_flux)?;
        self.iteration = it;
        self.update_boundaries()?;
        self.gradients()?;
        Ok(true)
    }

    /// Sets velocity and pressure and starts from the fluxes of the
    /// interpolated velocity (no pressure terms).
    pub fn initialise(&mut self, u: &Soa3, p: &[f64]) -> Result<()> {
        self.dev.upload(self.f.u.x, &u.x)?;
        self.dev.upload(self.f.u.y, &u.y)?;
        self.dev.upload(self.f.u.z, &u.z)?;
        self.dev.upload(self.f.p, p)?;
        self.update_boundaries()?;
        self.gradients()?;
        let a = self.flux_args(None, self.w.zero_face, None);
        self.dev.mass_flux(&a, self.f.mass_flux)
    }

    /// Potential-flow start (D52): the velocity of the potential φ that
    /// carries the prescribed boundary fluxes through the domain, ∇²φ = 0
    /// with ∂φ/∂n from the velocity inlets, ∂φ/∂n = 0 on walls and symmetry
    /// planes and φ = 0 on pressure boundaries. It is one pressure-correction
    /// solve from zero cell velocities with D = 1, so the AMG of the
    /// pressure solves it with the same boundary conditions; the face fluxes
    /// are continuous to the solver tolerance. The cell velocity comes from
    /// the face fluxes, u_P = Σ_f ±ṁ_f (x_f − x_P)/(ρ V_P), which is exact
    /// for a uniform flow and, unlike the gradient of φ, also holds next to
    /// the inlets (where the correction gradient assumes ∂φ/∂n = 0). The
    /// pressure stays at its start value, the multigrid of the run is built
    /// anew. Returns the CG iterations.
    pub fn potential_flow(&mut self, mesh: &Mesh, geo: &Geometry) -> Result<usize> {
        let (f, w) = (self.f, self.w);
        for b in f.u.all() {
            self.dev.fill(b, 0.0)?;
        }
        let kept = self.settings.clone();
        self.settings.pressure_tolerance = Some(1e-10);
        self.settings.pressure_max_iter = self.settings.pressure_max_iter.max(500);
        self.update_boundaries()?;
        let a = self.flux_args(None, w.zero_face, None);
        self.dev.mass_flux(&a, f.mass_flux)?;
        self.dev.fill(w.dc, 1.0)?;
        // the boundary fluxes, and per pass the non-orthogonal part of the
        // face flux −ρ (∇φ)_f·k_f of the potential of the pass before
        // (deferred, Jasak 1996), so that the potential is that of the full
        // Laplacian on non-orthogonal meshes
        let m0 = self.dev.download(f.mass_flux)?;
        let rho = self.phys.rho;
        let mut it = 0;
        let mut r = Ok(());
        for pass in 0..=POTENTIAL_NONORTH_PASSES {
            if pass > 0 {
                self.dev.gradient(
                    Lsq::Pressure,
                    w.pc,
                    Mask::PressureCorrection,
                    w.zero_boundary,
                    w.gpc,
                )?;
                let g = [
                    self.dev.download(w.gpc.x)?,
                    self.dev.download(w.gpc.y)?,
                    self.dev.download(w.gpc.z)?,
                ];
                let d = &self.disc;
                let mut ms = m0.clone();
                for (fi, m) in ms.iter_mut().enumerate().take(d.n_interior) {
                    let (p, n, wt) = (d.owner[fi] as usize, d.neighbour[fi] as usize, d.weight[fi]);
                    let k = [d.k.x[fi], d.k.y[fi], d.k.z[fi]];
                    let gk: f64 = (0..3)
                        .map(|i| (wt * g[i][p] + (1.0 - wt) * g[i][n]) * k[i])
                        .sum();
                    *m -= rho * gk;
                }
                self.dev.upload(f.mass_flux, &ms)?;
            }
            match self.pressure_correction(f.mass_flux, pass == 0) {
                Ok(n) => it += n,
                Err(e) => {
                    r = Err(e);
                    break;
                }
            }
        }
        self.settings = kept;
        r?;
        self.dev.correct_fluxes(w.cf, f.mass_flux, w.pc, w.mc)?;
        self.dev.copy(w.mc, f.mass_flux)?;
        let m = self.dev.download(f.mass_flux)?;
        let u = velocity_of_fluxes(mesh, geo, &m, self.phys.rho, self.disc.n_cells);
        self.dev.upload(f.u.x, &u.x)?;
        self.dev.upload(f.u.y, &u.y)?;
        self.dev.upload(f.u.z, &u.z)?;
        self.amg_setups = 0;
        self.amg_matrix = None;
        self.update_boundaries()?;
        self.gradients()?;
        self.dev.sync()?;
        Ok(it)
    }

    /// Outer iterations done.
    pub fn iteration(&self) -> usize {
        self.iteration
    }

    /// Device name and description.
    pub fn device(&self) -> String {
        self.dev.device()
    }

    /// Host copy of the fields.
    pub fn state(&mut self) -> Result<State> {
        let f = self.f;
        let mut v3 = |b: Buf3| -> Result<Soa3> {
            Ok(Soa3 {
                x: self.dev.download(b.x)?,
                y: self.dev.download(b.y)?,
                z: self.dev.download(b.z)?,
            })
        };
        Ok(State {
            u: v3(f.u)?,
            ub: v3(f.ub)?,
            grad_u: [v3(f.grad_u[0])?, v3(f.grad_u[1])?, v3(f.grad_u[2])?],
            grad_p: v3(f.grad_p)?,
            psi: v3(f.psi)?,
            p: self.dev.download(f.p)?,
            mass_flux: self.dev.download(f.mass_flux)?,
            pb: self.dev.download(f.pb)?,
        })
    }

    /// Replaces the fields by `s` (a restart).
    pub fn load(&mut self, s: &State) -> Result<()> {
        let f = self.f;
        let mut put = |b: Buf3, v: &Soa3| -> Result<()> {
            self.dev.upload(b.x, &v.x)?;
            self.dev.upload(b.y, &v.y)?;
            self.dev.upload(b.z, &v.z)
        };
        put(f.u, &s.u)?;
        put(f.ub, &s.ub)?;
        for i in 0..3 {
            put(f.grad_u[i], &s.grad_u[i])?;
        }
        put(f.grad_p, &s.grad_p)?;
        put(f.psi, &s.psi)?;
        self.dev.upload(f.p, &s.p)?;
        self.dev.upload(f.mass_flux, &s.mass_flux)?;
        self.dev.upload(f.pb, &s.pb)
    }

    /// New prescribed boundary velocity (time-dependent conditions); faces
    /// with another update rule ignore it.
    pub fn set_boundary_velocity(&mut self, ub: &[Vec<f64>; 3]) -> Result<()> {
        let s = self.f.ub_set;
        self.dev.upload(s.x, &ub[0])?;
        self.dev.upload(s.y, &ub[1])?;
        self.dev.upload(s.z, &ub[2])
    }

    pub(crate) fn flux_args(
        &self,
        d: Option<Buf>,
        relax: Buf,
        time: Option<(Buf, Buf, f64)>,
    ) -> FluxArgs {
        let f = &self.f;
        FluxArgs {
            u: f.u,
            grad_u: f.grad_u,
            p: f.p,
            grad_p: f.grad_p,
            pb: f.pb,
            ub: f.ub,
            d,
            relax,
            time,
        }
    }

    /// Step 1: boundary values from the current fields.
    pub fn update_boundaries(&mut self) -> Result<()> {
        let f = &self.f;
        self.dev.boundary_values(&BoundaryArgs {
            u: f.u,
            grad_u: f.grad_u,
            p: f.p,
            grad_p: f.grad_p,
            mass_flux: f.mass_flux,
            ub_set: f.ub_set,
            ub: f.ub,
            pb: f.pb,
        })
    }

    /// Step 2: gradients and limiters.
    pub fn gradients(&mut self) -> Result<()> {
        let f = self.f;
        for i in 0..3 {
            let (u, ub) = (f.u.all()[i], f.ub.all()[i]);
            self.dev
                .gradient(Lsq::Velocity, u, Mask::Velocity, ub, f.grad_u[i])?;
        }
        self.dev
            .gradient(Lsq::Pressure, f.p, Mask::Pressure, f.pb, f.grad_p)?;
        let s = &self.settings;
        if s.scheme == Scheme::Upwind || s.limiter == LimiterKind::None {
            return Ok(());
        }
        let (lim, kk) = (s.limiter, s.limiter_k);
        for i in 0..3 {
            let (u, ub, psi) = (f.u.all()[i], f.ub.all()[i], f.psi.all()[i]);
            let eps2 = match lim {
                LimiterKind::VenkatakrishnanWang => {
                    let (mn, mx) = self.dev.min_max(u)?;
                    let e = kk * (mx - mn);
                    Eps2::Single(e * e)
                }
                _ => Eps2::PerCell,
            };
            self.dev.limiter(u, ub, f.grad_u[i], eps2, psi)?;
        }
        Ok(())
    }

    /// Step 3: the relaxation term (1 − α)(ṁ − ρ ū·S) of the current fields
    /// into `w.relax` (and returns nothing else).
    pub(crate) fn relaxation_term(&mut self, alpha: f64) -> Result<()> {
        let a = self.flux_args(None, self.w.zero_face, None);
        self.dev.mass_flux(&a, self.w.relax)?;
        self.dev
            .scaled_difference(self.w.relax, 1.0 - alpha, self.f.mass_flux, self.w.relax)
    }

    /// Step 4 (assembly): momentum faces with the convecting flux `conv`
    /// and the rows; returns the residual scale.
    pub(crate) fn momentum_assembly(
        &mut self,
        conv: Buf,
        alpha: f64,
        theta: f64,
        time: Option<TimeTerms>,
    ) -> Result<f64> {
        self.momentum_assembly_faces(conv)?;
        let (f, w) = (self.f, self.w);
        self.dev.momentum_rows(&MomentumRowArgs {
            cp: w.cp,
            cn: w.cn,
            bf: w.bf,
            grad_p: f.grad_p,
            body_force: self.phys.body_force,
            u: f.u,
            alpha,
            theta,
            time,
            val: w.val,
            rhs: w.rhs,
            d: w.d,
            ds: w.ds,
            dc: w.dc,
            diag: w.diag,
        })?;
        if self.settings.algorithm(self.transient) != Algorithm::Simplec {
            self.dev.copy(w.d, w.dc)?;
        }
        // residual scale ‖a_P u‖ + ‖b‖ over the three components (unrelaxed
        // a_P and b), so that it is finite also from a fluid at rest
        self.dev.residual_scale(f.u, w.rhs, w.diag, alpha)
    }

    /// The momentum face coefficients of the current fields convected by
    /// `conv` into `w.cp`, `w.cn`, `w.bf`.
    pub(crate) fn momentum_assembly_faces(&mut self, conv: Buf) -> Result<()> {
        let s = &self.settings;
        let (f, w) = (self.f, self.w);
        let uc = self.central_u.unwrap_or(f.u);
        self.dev.momentum_faces(&MomentumFaceArgs {
            mass_flux: conv,
            u: f.u,
            uc,
            grad_u: f.grad_u,
            psi: f.psi,
            ub: f.ub,
            blend: if s.scheme == Scheme::Upwind { 0.0 } else { 1.0 },
            central: if s.scheme == Scheme::Central {
                s.central_blending
            } else {
                0.0
            },
            beta: s.central_bounding,
            nonorth: if s.nonorth_correction { 1.0 } else { 0.0 },
            convection: if self.phys.convection { 1.0 } else { 0.0 },
            cp: w.cp,
            cn: w.cn,
            bf: w.bf,
        })
    }

    /// Step 4 (solve): the three momentum components; residuals over
    /// `scale` and iterations into `res`.
    pub(crate) fn momentum_solve(&mut self, scale: f64, res: &mut Residuals) -> Result<()> {
        let (f, w) = (self.f, self.w);
        let sym = self.bc.mom_kind.contains(&k::MOM_SYMMETRY);
        let s = &self.settings;
        let (tol, max) = (s.momentum_tol(self.transient), s.momentum_max_iter);
        for i in 0..3 {
            if sym {
                // the diagonal of this component (symmetry planes)
                self.dev.set_diagonal(w.val, w.diag.all()[i])?;
            }
            let (it, r0) = self.dev.solve(
                w.val,
                w.rhs.all()[i],
                f.u.all()[i],
                Preconditioner::Jacobi,
                false,
                tol,
                max,
            )?;
            res.momentum[i] = r0 / scale;
            res.momentum_iterations[i] = it;
        }
        Ok(())
    }

    /// Continuity residual of the fluxes `m`.
    pub(crate) fn continuity(&mut self, m: Buf) -> Result<f64> {
        let w = self.w;
        self.dev.net_outflow(m, w.defect)?;
        self.dev.throughflow(m, w.through)?;
        let d = self.dev.sum_sq(w.defect)?.sqrt();
        let t = self.dev.sum_sq(w.through)?.sqrt();
        Ok(d / t.max(f64::MIN_POSITIVE))
    }

    /// Step 6: the pressure correction of the fluxes `m` into `w.pc`
    /// (shifted to zero in the reference cell); `rebuild` sets the
    /// multigrid up first. Returns the CG iterations.
    pub(crate) fn pressure_correction(&mut self, m: Buf, rebuild: bool) -> Result<usize> {
        let t = Instant::now();
        let w = self.w;
        self.dev.pressure_faces(w.dc, w.cf)?;
        self.dev.pressure_rows(w.cf, m, w.pval, w.prhs)?;
        if let Some(r) = self.reference_cell {
            let pos = self.disc.sell.pos(r, 0);
            self.dev.scale_entry(w.pval, pos, 2.0)?;
        }
        self.dev.sync()?;
        self.timings.pressure_assembly += t.elapsed().as_secs_f64();
        let t = Instant::now();
        let s = self.settings.clone();
        if s.pressure_preconditioner == Preconditioner::Amg && (self.amg_setups == 0 || rebuild) {
            self.dev.amg_setup(w.pval, &s.amg)?;
            self.amg_matrix = Some(self.dev.download(w.pval)?);
            self.amg_setups += 1;
            self.timings.amg_setup += t.elapsed().as_secs_f64();
        }
        let t = Instant::now();
        self.dev.fill(w.pc, 0.0)?;
        let (it, _) = self.dev.solve(
            w.pval,
            w.prhs,
            w.pc,
            s.pressure_preconditioner,
            true,
            s.pressure_tol(self.transient),
            s.pressure_max_iter,
        )?;
        if let Some(r) = self.reference_cell {
            let shift = self.dev.value(w.pc, r)?;
            self.dev.subtract(w.pc, shift)?;
        }
        self.dev.sync()?;
        self.timings.pressure_solve += t.elapsed().as_secs_f64();
        Ok(it)
    }

    /// Step 7: corrected fluxes `m` → `w.mc`, u −= D^C ∇p', p += α_p p'.
    pub(crate) fn corrections(&mut self, m: Buf, alpha_p: f64) -> Result<()> {
        let (f, w) = (self.f, self.w);
        self.dev.correct_fluxes(w.cf, m, w.pc, w.mc)?;
        // gradient of p': zero on pressure boundaries, zero normal gradient
        // elsewhere
        self.dev.gradient(
            Lsq::Pressure,
            w.pc,
            Mask::PressureCorrection,
            w.zero_boundary,
            w.gpc,
        )?;
        self.dev.correct_velocity(w.dc, w.gpc, f.u)?;
        self.dev.lincomb(f.p, 1.0, f.p, alpha_p, w.pc)
    }

    /// Rebuilds the AMG hierarchy of a restart from the stored pressure
    /// matrix, so that the run continues bit-identically.
    pub fn restore_amg(&mut self, m: &[f64], setups: usize) -> Result<()> {
        let w = self.w;
        self.dev.upload(w.pval, m)?;
        let amg = self.settings.amg.clone();
        self.dev.amg_setup(w.pval, &amg)?;
        self.amg_setups = setups;
        self.amg_matrix = Some(m.to_vec());
        Ok(())
    }

    /// Re-validated settings that may change during a run (D40).
    pub fn set_settings(&mut self, s: Settings) -> Result<()> {
        s.validate(self.transient)?;
        self.settings = s;
        Ok(())
    }

    /// One steady outer iteration; returns the residuals before the solves.
    pub fn iterate(&mut self) -> Result<Residuals> {
        self.iteration += 1;
        let mut res = Residuals {
            iteration: self.iteration,
            ..Default::default()
        };
        let s = self.settings.clone();
        let alpha = s.alpha_u(false);
        let t = Instant::now();
        self.update_boundaries()?;
        self.gradients()?;
        self.dev.sync()?;
        self.timings.gradients += t.elapsed().as_secs_f64();

        // step 3: relaxation term of the Rhie–Chow flux
        let t = Instant::now();
        self.relaxation_term(alpha)?;
        // step 4: momentum
        let scale = self.momentum_assembly(self.f.mass_flux, alpha, 1.0, None)?;
        self.dev.sync()?;
        self.timings.momentum_assembly += t.elapsed().as_secs_f64();
        let t = Instant::now();
        self.momentum_solve(scale, &mut res)?;
        self.dev.sync()?;
        self.timings.momentum_solve += t.elapsed().as_secs_f64();

        // step 5: predicted fluxes
        let t = Instant::now();
        let a = self.flux_args(Some(self.w.d), self.w.relax, None);
        self.dev.mass_flux(&a, self.w.m)?;
        res.continuity = self.continuity(self.w.m)?;
        self.dev.sync()?;
        self.timings.pressure_assembly += t.elapsed().as_secs_f64();

        // step 6: pressure correction
        let rebuild = (self.iteration - 1).is_multiple_of(s.amg_rebuild);
        res.pressure_iterations = self.pressure_correction(self.w.m, rebuild)?;

        // step 7: corrections
        let t = Instant::now();
        self.corrections(self.w.m, s.alpha_p(false))?;
        self.dev.sync()?;
        self.timings.pressure_assembly += t.elapsed().as_secs_f64();
        if !res.finite() {
            self.diverged = Some(res);
            return Err(crate::diag::err(crate::diag::Diagnostic::new(
                crate::diag::RUN_DIVERGED,
                format!(
                    "the solution diverged in outer iteration {} (residuals {:?}, continuity {:.3e})",
                    self.iteration, res.momentum, res.continuity
                ),
            )));
        }
        std::mem::swap(&mut self.f.mass_flux, &mut self.w.mc);
        Ok(res)
    }
}

/// Cell velocities from face mass fluxes, u_P = Σ_f ±ṁ_f (x_f − x_P)/(ρ V_P):
/// exact for a uniform velocity, since Σ_f (x_f − x_P) S_fᵀ = V_P I for a
/// closed cell (the start of D52).
fn velocity_of_fluxes(mesh: &Mesh, geo: &Geometry, m: &[f64], rho: f64, nc: usize) -> Soa3 {
    let ni = mesh.n_interior_faces();
    let mut u = Soa3::zeros(nc);
    for (fi, &mf) in m.iter().enumerate() {
        let xf = geo.face_centroid[fi];
        let sides = [(mesh.owner[fi] as usize, 1.0), (usize::MAX, 0.0)];
        let sides = if fi < ni {
            [sides[0], (mesh.neighbour[fi] as usize, -1.0)]
        } else {
            sides
        };
        for (c, s) in sides {
            if c < nc {
                let xp = geo.cell_centroid[c];
                u.x[c] += s * mf * (xf[0] - xp[0]);
                u.y[c] += s * mf * (xf[1] - xp[1]);
                u.z[c] += s * mf * (xf[2] - xp[2]);
            }
        }
    }
    for c in 0..nc {
        let r = 1.0 / (rho * geo.cell_volume[c]);
        u.x[c] *= r;
        u.y[c] *= r;
        u.z[c] *= r;
    }
    u
}

/// Passes of the non-orthogonal correction of the potential-flow start
/// (D52) after the first solve.
const POTENTIAL_NONORTH_PASSES: usize = 4;

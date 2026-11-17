//! RANS turbulence models on the flow solver (stage 6, D129–D139): the
//! k-ω SST model, 2003 version (D130, stage 6 PR 6b), with the production
//! options of D138, the low-Re and automatic wall treatments and rough
//! walls (D134), the inlet values (D137), the modified pressure
//! p̂ = p + ⅔ρk and the exact or Poisson wall distance (D48).
//!
//! Case file:
//!
//! ```toml
//! [turbulence]
//! model = "k-omega-sst"
//! wall_treatment = "automatic"        # or "low-re"
//! production_limiter = "menter"       # "kato-launder", "durbin", "none"
//! # ambient = { k = 1e-6, omega = 1.0 }  # sustaining terms
//! # wall_distance = "pde"             # Poisson distance instead of the exact one
//! # constants = { a1 = 0.31 }         # TURB007
//!
//! [initial]
//! turbulence_intensity = 0.01          # or k and omega (numbers or expressions)
//! viscosity_ratio = 1.0
//!
//! [boundary.inlet]
//! type = "velocity-inlet"
//! velocity = [10.0, 0.0, 0.0]
//! turbulence_intensity = 0.05          # with turbulence_length_scale or viscosity_ratio,
//! viscosity_ratio = 10.0               # or k and omega
//!
//! [solver.turbulence]
//! relaxation = 0.8
//! ```
//!
//! Coupling (stage 6 specification 6.4): after the pressure correction of
//! every outer iteration (and of every step of a transient run), with the
//! velocity gradient of the corrected velocity:
//!
//! 1. inlet values; gradients and limiters of k and ω;
//! 2. the SST sources per cell (F₁, F₂, CD_kω, P̃_k and the split of
//!    section 6.2), the wall treatment per wall face and the ω of the
//!    near-wall cells;
//! 3. the face diffusivities, then the k equation and the ω equation (the
//!    near-wall rows replaced by the fixed ω), the positivity floors;
//! 4. ν_t from the new fields, μ_eff = μ + ρν_t on the faces (the wall
//!    faces take μ_w of the wall treatment), the pressure boundaries shifted
//!    to p̂.
//!
//! The momentum equation then takes the face viscosity μ_eff and the
//! transposed part ∇·[μ_eff (∇u)ᵀ] of the stress as a deferred term; the
//! pressure of the solver is p̂, the output and the forces take the static
//! pressure p = p̂ − ⅔ρk ([`crate::simplec::Flow::state`]).

pub mod inlet;
pub mod wall;

use crate::bc::{BoundaryArrays, BoundaryCondition};
use crate::device::{
    Buf, Buf3, Eps2, FlowDevice, Lsq, Mask, ScalarFaceArgs, ScalarRowArgs, ScalarTime, SstCellOps,
    TurbFaceOps, WallCellOps, WallFaceOps,
};
use crate::diag::{self, Diagnostic};
use crate::expr::{self, Env, Scalar, Vars};
use anyhow::Result;
use foehncraft_fv::problem::BoundaryKind;
use foehncraft_fv::solver::{LimiterKind, Preconditioner, Scheme};
use foehncraft_kernels::turb::sst::{self as ksst, SstConstants, SstOptions};
use foehncraft_kernels::turb::wall as kwall;
use foehncraft_mesh::{Geometry, Mesh};
use inlet::InletSpec;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use wall::Walls;

/// `turbulence.model`.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum TurbulenceModel {
    /// No turbulence model.
    #[default]
    Laminar,
    /// k-ω SST, 2003 version (Menter, Kuntz & Langtry 2003).
    KOmegaSst,
}

/// `turbulence.wall_treatment` (D134).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum WallTreatment {
    /// The y⁺-insensitive blend of the viscous and the logarithmic layer
    /// (the default of SST).
    Automatic,
    /// The integrated form to the wall: k_w = 0, ω of the viscous sublayer
    /// in the first cell (first-cell y⁺ ≲ 1).
    LowRe,
}

/// `turbulence.production_limiter` (D138).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum ProductionLimiter {
    /// P̃_k = min(ν_t S², 10 β* k ω) (2003).
    Menter,
    /// Kato & Launder (1993): P_k = ν_t S Ω, then the 2003 limiter.
    KatoLaunder,
    /// Durbin (1996): ν_t ≤ C_lim k/(√6 |S|), then the 2003 limiter.
    Durbin,
    /// P_k = ν_t S², unlimited.
    None,
}

/// `turbulence.wall_distance` (D48).
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum WallDistanceMethod {
    /// The exact distance to the nearest wall face (bounding-volume
    /// hierarchy).
    #[default]
    Exact,
    /// From the Poisson equation ∇²φ = −1 (Tucker 2003).
    Pde,
}

/// `solver.turbulence.destruction`.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Destruction {
    /// −βω² ≈ βω*² − 2βω* ω.
    #[default]
    Newton,
    /// −βω² ≈ −βω* ω.
    Picard,
}

/// `turbulence.ambient`: the sustaining terms of Spalart & Rumsey (2007).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Ambient {
    /// Ambient k (m²/s²).
    pub k: f64,
    /// Ambient ω (1/s).
    pub omega: f64,
}

/// `[turbulence]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct TurbulenceSection {
    /// Turbulence model (D12, D130).
    pub model: TurbulenceModel,
    /// Momentum wall treatment (D134; SST default `automatic`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wall_treatment: Option<WallTreatment>,
    /// Production form of SST (default `menter`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub production_limiter: Option<ProductionLimiter>,
    /// C_lim of the Durbin bound.
    pub durbin_coefficient: f64,
    /// Sustaining terms (Spalart & Rumsey 2007): ambient k and ω.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ambient: Option<Ambient>,
    /// Cap of ν_t/ν (counted, TURB008).
    pub max_viscosity_ratio: f64,
    /// Smooth min and max in the primal solution (specification 2.6).
    pub smooth_limits: bool,
    /// Wall distance (D48).
    pub wall_distance: WallDistanceMethod,
    /// Constant overrides by their literature name (`beta_star`, `a1`,
    /// `kappa`, `gamma1`, `gamma2`, `beta1`, `beta2`, `sigma_k1`,
    /// `sigma_k2`, `sigma_omega1`, `sigma_omega2`); TURB007.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub constants: BTreeMap<String, f64>,
}

impl Default for TurbulenceSection {
    fn default() -> Self {
        Self {
            model: TurbulenceModel::Laminar,
            wall_treatment: None,
            production_limiter: None,
            durbin_coefficient: 1.0,
            ambient: None,
            max_viscosity_ratio: 1e5,
            smooth_limits: false,
            wall_distance: WallDistanceMethod::Exact,
            constants: BTreeMap::new(),
        }
    }
}

/// The names of the SST constants that `[turbulence.constants]` takes.
pub const SST_CONSTANTS: &[&str] = &[
    "beta_star",
    "a1",
    "kappa",
    "gamma1",
    "gamma2",
    "beta1",
    "beta2",
    "sigma_k1",
    "sigma_k2",
    "sigma_omega1",
    "sigma_omega2",
];

/// Relative width of the smooth switches with `smooth_limits = true`.
pub const SMOOTH_WIDTH: f64 = 1e-2;

impl TurbulenceSection {
    /// The default table (laminar flow; left out of the written case).
    pub fn is_default(&self) -> bool {
        *self == Self::default()
    }

    /// No turbulence model.
    pub fn is_laminar(&self) -> bool {
        self.model == TurbulenceModel::Laminar
    }

    /// The model constants with the overrides.
    pub fn sst_constants(&self) -> Result<SstConstants> {
        let mut c = SstConstants::default();
        for (name, &v) in &self.constants {
            let slot = match name.as_str() {
                "beta_star" => &mut c.beta_star,
                "a1" => &mut c.a1,
                "kappa" => &mut c.kappa,
                "gamma1" => &mut c.gamma1,
                "gamma2" => &mut c.gamma2,
                "beta1" => &mut c.beta1,
                "beta2" => &mut c.beta2,
                "sigma_k1" => &mut c.sigma_k1,
                "sigma_k2" => &mut c.sigma_k2,
                "sigma_omega1" => &mut c.sigma_w1,
                "sigma_omega2" => &mut c.sigma_w2,
                _ => {
                    let hint = crate::bc::did_you_mean(name, SST_CONSTANTS.iter().copied())
                        .map(|h| format!(" (did you mean `{h}`?)"))
                        .unwrap_or_default();
                    return Err(diag::err(
                        Diagnostic::new(
                            diag::CASE_UNKNOWN_KEY,
                            format!(
                                "[turbulence.constants] has no constant `{name}`{hint}; the SST \
                                 constants are {}",
                                SST_CONSTANTS.join(", ")
                            ),
                        )
                        .key(format!("turbulence.constants.{name}")),
                    ));
                }
            };
            crate::need!(
                v > 0.0 && v.is_finite(),
                CASE_RANGE,
                "[turbulence.constants] {name} must be > 0, is {v}"
            );
            *slot = v;
        }
        Ok(c)
    }

    /// The run-time options of the kernel.
    pub fn sst_options(&self, s: &TurbulenceSettings) -> SstOptions {
        let (k_amb, w_amb) = self.ambient.map_or((0.0, 0.0), |a| (a.k, a.omega));
        SstOptions {
            limiter: match self.production_limiter.unwrap_or(ProductionLimiter::Menter) {
                ProductionLimiter::Menter => ksst::LIMITER_MENTER,
                ProductionLimiter::KatoLaunder => ksst::LIMITER_KATO_LAUNDER,
                ProductionLimiter::Durbin => ksst::LIMITER_DURBIN,
                ProductionLimiter::None => ksst::LIMITER_NONE,
            },
            newton: u32::from(s.destruction == Destruction::Newton),
            c_lim: self.durbin_coefficient,
            smooth: if self.smooth_limits {
                SMOOTH_WIDTH
            } else {
                0.0
            },
            k_amb,
            w_amb,
            max_ratio: self.max_viscosity_ratio,
        }
    }

    /// The wall treatment (the model's default when not set).
    pub fn treatment(&self) -> WallTreatment {
        self.wall_treatment.unwrap_or(WallTreatment::Automatic)
    }

    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        if self.is_laminar() {
            return Ok(());
        }
        self.sst_constants()?;
        crate::need!(
            self.durbin_coefficient > 0.0 && self.durbin_coefficient.is_finite(),
            CASE_RANGE,
            "turbulence.durbin_coefficient must be > 0, is {}",
            self.durbin_coefficient
        );
        crate::need!(
            self.max_viscosity_ratio > 0.0,
            CASE_RANGE,
            "turbulence.max_viscosity_ratio must be > 0, is {}",
            self.max_viscosity_ratio
        );
        if let Some(a) = self.ambient {
            crate::need!(
                a.k > 0.0 && a.omega > 0.0 && a.k.is_finite() && a.omega.is_finite(),
                CASE_RANGE,
                "turbulence.ambient: k and omega must be > 0"
            );
        }
        Ok(())
    }
}

/// `[solver.turbulence]`: the transport equations of the turbulence model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct TurbulenceSettings {
    /// Convection scheme (`central` only in transient runs, TURB006).
    pub scheme: Scheme,
    /// Limiter of the gradients (`venkatakrishnan` and
    /// `venkatakrishnan-wang` are the same for k and ω: the threshold is
    /// scaled with the variable's range).
    pub limiter: LimiterKind,
    /// Relative threshold K of the limiter, ε = K (φ_max − φ_min) over the
    /// domain (Wang 2000). The default 0 (taken as 10⁻¹²) is the bounded
    /// Venkatakrishnan function: the dimensional threshold (K h)³ of the
    /// momentum equations is far above the differences of k (10⁻⁸ to
    /// 10⁻³ m²/s² on a flat plate) and would switch the limiter off.
    pub limiter_k: f64,
    /// Under-relaxation of k and ω (default 0.8 steady, 1 transient).
    pub relaxation: Option<f64>,
    /// Relative tolerance per solve (default 1e-2 steady, 1e-8 transient).
    pub tolerance: Option<f64>,
    /// Largest number of BiCGStab iterations per solve.
    pub max_iter: usize,
    /// Linearisation of the ω destruction.
    pub destruction: Destruction,
    /// Relaxation α_νt of ν_t.
    pub viscosity_relaxation: f64,
}

impl Default for TurbulenceSettings {
    fn default() -> Self {
        Self {
            scheme: Scheme::LinearUpwind,
            limiter: LimiterKind::Venkatakrishnan,
            limiter_k: 0.0,
            relaxation: None,
            tolerance: None,
            max_iter: 50,
            destruction: Destruction::Newton,
            viscosity_relaxation: 1.0,
        }
    }
}

impl TurbulenceSettings {
    /// α of k and ω.
    pub fn alpha(&self, transient: bool) -> f64 {
        self.relaxation.unwrap_or(if transient { 1.0 } else { 0.8 })
    }

    /// Relative tolerance of a solve.
    pub fn tol(&self, transient: bool) -> f64 {
        self.tolerance
            .unwrap_or(if transient { 1e-8 } else { 1e-2 })
    }

    /// Checks the ranges.
    pub fn validate(&self, transient: bool) -> Result<()> {
        let a = self.alpha(transient);
        crate::need!(
            a > 0.0 && a <= 1.0,
            CASE_RANGE,
            "[solver.turbulence] relaxation must be in (0, 1], is {a}"
        );
        crate::need!(
            self.viscosity_relaxation > 0.0 && self.viscosity_relaxation <= 1.0,
            CASE_RANGE,
            "[solver.turbulence] viscosity_relaxation must be in (0, 1], is {}",
            self.viscosity_relaxation
        );
        crate::need!(
            self.tol(transient) > 0.0,
            CASE_RANGE,
            "[solver.turbulence] tolerance must be > 0"
        );
        crate::need!(
            self.max_iter >= 1,
            CASE_RANGE,
            "[solver.turbulence] max_iter must be >= 1"
        );
        crate::need!(
            self.limiter_k >= 0.0 && self.limiter_k.is_finite(),
            CASE_RANGE,
            "[solver.turbulence] limiter_k must be >= 0"
        );
        crate::need!(
            transient || self.scheme != Scheme::Central,
            TURB_UNBOUNDED_SCHEME,
            "[solver.turbulence] scheme = \"central\" is not bounded for k and ω in a steady run; \
             use \"linear-upwind\""
        );
        Ok(())
    }
}

/// The start values of the turbulence (`[initial]`, specification 7.3).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InitialTurbulence {
    /// Intensity.
    pub turbulence_intensity: Option<f64>,
    /// ν_t/ν.
    pub viscosity_ratio: Option<f64>,
    /// k directly (number or expression of x, y, z).
    pub k: Option<Scalar>,
    /// ω directly.
    pub omega: Option<Scalar>,
}

/// Checks the turbulence keys of a case.
pub fn validate(case: &crate::case::FlowCase) -> Result<()> {
    let t = &case.turbulence;
    t.validate()?;
    if t.is_laminar() {
        let i = &case.initial;
        crate::need!(
            i.turbulence_intensity.is_none()
                && i.viscosity_ratio.is_none()
                && i.k.is_none()
                && i.omega.is_none(),
            CASE_INCONSISTENT,
            "[initial] turbulence values in a laminar case (turbulence.model = \"laminar\")"
        );
        return Ok(());
    }
    let transient = case.time.is_some();
    case.solver.turbulence.validate(transient)?;
    let i = &case.initial;
    crate::need!(
        i.k.is_some() == i.omega.is_some(),
        CASE_INCONSISTENT,
        "[initial] k and omega go together"
    );
    crate::need!(
        i.k.is_none() || (i.turbulence_intensity.is_none() && i.viscosity_ratio.is_none()),
        CASE_INCONSISTENT,
        "[initial]: k and omega, or turbulence_intensity and viscosity_ratio, not both"
    );
    let env = case.env()?;
    for (key, v) in [("initial.k", &i.k), ("initial.omega", &i.omega)] {
        if let Some(v) = v {
            expr::value(v, &env, Vars::SPACE, key)?;
        }
    }
    for (name, bc) in &case.boundary {
        for (error, code, msg) in inlet::check(bc) {
            if error {
                return Err(diag::err(
                    Diagnostic::new(code, format!("[boundary.{name}]: {msg}"))
                        .key(format!("boundary.{name}")),
                ));
            }
        }
        if inlet::is_inlet(bc) || matches!(bc, BoundaryCondition::Wall { .. }) {
            InletSpec::of(
                &bc.turbulence(),
                &env,
                transient,
                &format!("boundary.{name}"),
            )?;
        }
        if let BoundaryCondition::Wall { .. } = bc {
            let ks = bc.roughness();
            crate::need!(
                ks >= 0.0 && ks.is_finite(),
                CASE_RANGE,
                "[boundary.{name}] roughness must be >= 0, is {ks}"
            );
        }
    }
    Ok(())
}

/// The run-time findings of the setup (notes and warnings for the log).
pub fn setup_findings(case: &crate::case::FlowCase) -> Vec<crate::bc::Finding> {
    use crate::bc::{Finding, Level};
    let mut out = Vec::new();
    let t = &case.turbulence;
    if t.is_laminar() {
        return out;
    }
    let mut push = |level, code, message: String| {
        out.push(Finding {
            level,
            code,
            boundary: None,
            message,
        })
    };
    if !t.constants.is_empty() {
        push(
            Level::Warning,
            diag::TURB_CONSTANTS,
            format!(
                "SST constants changed in [turbulence.constants]: {}",
                t.constants
                    .iter()
                    .map(|(k, v)| format!("{k} = {v}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        );
    }
    if !case
        .boundary
        .values()
        .any(|b| matches!(b, BoundaryCondition::Wall { .. }))
    {
        push(
            Level::Note,
            diag::TURB_NO_WALL,
            "no wall boundary: the wall distance is infinite (F₁ = 0, F₂ = 0)".into(),
        );
    }
    if case
        .time
        .as_ref()
        .is_some_and(|t| t.scheme == crate::transient::TimeScheme::CrankNicolson)
    {
        push(
            Level::Note,
            diag::TURB_UNBOUNDED_SCHEME,
            "Crank–Nicolson is not bounded for k and ω: they take BDF2".into(),
        );
    }
    out
}

/// Results of the turbulence solves of one outer iteration.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct TurbResiduals {
    /// Scaled residuals of k and ω (before the solves).
    pub residual: Vec<f64>,
    /// BiCGStab iterations of k and ω.
    pub iterations: Vec<usize>,
    /// Cells where the positivity floor acted, k and ω.
    pub floors: Vec<usize>,
    /// Cells where the ν_t cap acted.
    pub capped: usize,
}

impl TurbResiduals {
    /// No turbulence model.
    pub fn is_empty(&self) -> bool {
        self.residual.is_empty()
    }
}

/// The names of the equations, in the order of the residuals.
pub const EQUATIONS: [&str; 2] = ["k", "omega"];

/// Host copy of the turbulence state (restart, output).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TurbState {
    /// k per cell.
    pub k: Vec<f64>,
    /// ω per cell.
    pub omega: Vec<f64>,
    /// ν_t per cell.
    pub nu_t: Vec<f64>,
    /// Transient: k and ω of the level before.
    pub k_old: Option<Vec<f64>>,
    /// See `k_old`.
    pub omega_old: Option<Vec<f64>>,
    /// The face viscosity μ_eff of the momentum equation.
    pub mu_eff: Vec<f64>,
    /// The boundary values of k and ω per boundary face (inlets and rough
    /// walls).
    pub k_b: Vec<f64>,
    /// See `k_b`.
    pub omega_b: Vec<f64>,
    /// u_τ per boundary face (0 off the walls).
    pub u_tau: Vec<f64>,
    /// First-cell y⁺ per boundary face (0 off the walls).
    pub y_plus: Vec<f64>,
    /// The modified pressure p̂ = p + ⅔ρk of the solver, per cell (the
    /// state's `p` is the static pressure).
    pub p_modified: Vec<f64>,
    /// The modified pressure per boundary face.
    pub pb_modified: Vec<f64>,
    /// The wall distance per cell.
    pub wall_distance: Vec<f64>,
}

/// Time terms of a step.
pub type StepTime = crate::scalar::StepTime;

/// One transported turbulence variable on the device.
#[derive(Debug, Clone, Copy)]
struct Var {
    phi: Buf,
    grad: Buf3,
    psi: Buf,
    bc_a: Buf,
    bc_b: Buf,
    /// Inlet values per boundary face (0 elsewhere).
    b_in: Buf,
    /// Wall values per boundary face (rough walls; 0 elsewhere).
    b_wall: Buf,
    lsq: Buf3,
    old: Option<(Buf, Buf, Buf)>,
    valid: Option<Buf>,
    floor: f64,
}

/// The wall treatment on the device.
#[derive(Debug, Clone, Copy)]
struct WallDev {
    normal: Buf3,
    y: Buf,
    ks: Buf,
    u_tau: Buf,
    y_plus: Buf,
    omega: Buf,
    blending: Buf,
    pk: Buf,
    weight: Buf,
}

/// The inlet faces.
#[derive(Debug, Clone)]
struct Inlets {
    specs: Vec<InletSpec>,
    /// Per inlet face: boundary face index, spec, centroid.
    faces: Vec<(u32, usize, [f64; 3])>,
    speed: bool,
    time: bool,
}

/// The k-ω SST model of a flow on its device.
pub struct Turbulence {
    /// The model settings.
    pub section: TurbulenceSection,
    /// The solver settings.
    pub settings: TurbulenceSettings,
    c: SstConstants,
    o: SstOptions,
    treatment: u32,
    nu: f64,
    rho: f64,
    nc: usize,
    ni: usize,
    nb: usize,
    k: Var,
    w: Var,
    nu_t: Buf,
    nu_t_prev: Option<Buf>,
    f1: Buf,
    pk: Buf,
    dk: Buf,
    dw: Buf,
    capped: Buf,
    src: [Buf; 4],
    d: Buf,
    gu: [Buf3; 3],
    mu_lam: Buf,
    gamma_k: Buf,
    gamma_w: Buf,
    mask_b: Buf,
    walls: Walls,
    wd: Option<WallDev>,
    fix_mask: Buf,
    fix_val: Buf,
    wall_pk: Buf,
    zero_face: Buf,
    cp: Buf,
    cn: Buf,
    bf: Buf,
    val: Buf,
    rhs: Buf,
    diag: Buf,
    order: Buf,
    work: Buf,
    inlets: Inlets,
    /// Host copies of a and of the inlet b of k (for p̂ at the pressure
    /// boundaries).
    ka_host: Vec<f64>,
    kin_host: Vec<f64>,
    /// The pressure faces (boundary face index) and their owners.
    pressure_faces: Vec<u32>,
    pressure_owners: Vec<u32>,
    added: Option<(Buf, Buf)>,
    t_set: f64,
    /// The wall distance per cell (host).
    pub wall_distance: Vec<f64>,
    /// Health (D69): per equation (k, ω) the iterations in which the floor
    /// acted on more than 0.1 % of the cells, and the most cells floored in
    /// one iteration; the most cells with a capped ν_t.
    pub health: TurbHealth,
}

/// A wall's y⁺ statistics: boundary, min, mean, max, area fraction above 2.
pub type YPlusRow = (String, f64, f64, f64, f64);

/// Floor and cap counts of a run (D69, TURB002).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct TurbHealth {
    /// Iterations with floors on more than 0.1 % of the cells, k and ω.
    pub floor_iterations: [usize; 2],
    /// Of those, the ones in the final tenth of the solves so far (the
    /// inexact linear solves of the first iterations may undershoot; a
    /// converging run does not need the floor at its end).
    pub floor_iterations_late: [usize; 2],
    /// Per solve, a bit per equation: the floor acted on more than 0.1 %.
    #[serde(skip)]
    flags: Vec<u8>,
    /// Most cells floored in one iteration, k and ω.
    pub floor_cells_max: [usize; 2],
    /// Most cells with ν_t capped by `max_viscosity_ratio` in one iteration.
    pub capped_max: usize,
}

/// A new buffer with the values `v`.
fn up(dev: &mut dyn FlowDevice, v: &[f64]) -> Result<Buf> {
    let b = dev.alloc(v.len())?;
    dev.upload(b, v)?;
    Ok(b)
}

/// The boundary kinds of k and ω per boundary (for the gradient weights)
/// and the coefficient a per boundary face.
fn kinds(
    mesh: &Mesh,
    bcs: &BTreeMap<String, BoundaryCondition>,
    treatment: WallTreatment,
    omega: bool,
) -> (BTreeMap<String, BoundaryKind>, Vec<f64>) {
    let ni = mesh.n_interior_faces();
    let mut a = vec![1.0; mesh.n_boundary_faces()];
    let mut points = BTreeMap::new();
    for p in &mesh.boundaries {
        let bc = bcs.get(&p.name);
        let fixed = match bc {
            Some(b) if inlet::is_inlet(b) => true,
            Some(b @ BoundaryCondition::Wall { .. }) => {
                omega || b.roughness() > 0.0 || treatment == WallTreatment::LowRe
            }
            _ => false,
        };
        let kind = if fixed {
            BoundaryKind::FixedValue
        } else {
            BoundaryKind::ZeroGradient
        };
        points.insert(p.name.clone(), kind);
        if fixed {
            for f in p.range() {
                a[f - ni] = 0.0;
            }
        }
    }
    (points, a)
}

impl Turbulence {
    /// The SST model of a case on the device `dev` of a flow.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        dev: &mut dyn FlowDevice,
        case: &crate::case::FlowCase,
        mesh: &Mesh,
        geo: &Geometry,
        disc: &foehncraft_fv::disc::Discretisation,
        bc_arrays: &BoundaryArrays,
        rho: f64,
        mu: f64,
        face_viscosity: Buf,
    ) -> Result<Self> {
        let section = case.turbulence.clone();
        let settings = case.solver.turbulence.clone();
        let transient = case.time.is_some();
        let env = case.env()?;
        let (nc, nf, ni) = (disc.n_cells, disc.n_faces, disc.n_interior);
        let nb = nf - ni;
        let stored = disc.sell.stored();
        let c = section.sst_constants()?;
        let o = section.sst_options(&settings);
        let treatment = section.treatment();
        let nu = mu / rho;
        // the inlets
        let mut specs = Vec::new();
        let mut faces = Vec::new();
        let mut spec_of = BTreeMap::new();
        for p in &mesh.boundaries {
            let Some(bc) = case.boundary.get(&p.name) else {
                continue;
            };
            if !inlet::is_inlet(bc) {
                continue;
            }
            let s = InletSpec::of(
                &bc.turbulence(),
                &env,
                transient,
                &format!("boundary.{}", p.name),
            )?;
            spec_of.insert(p.name.clone(), specs.len());
            specs.push(s);
            for f in p.range() {
                faces.push(((f - ni) as u32, specs.len() - 1, geo.face_centroid[f]));
            }
        }
        let inlets = Inlets {
            speed: specs.iter().any(InletSpec::speed_dependent),
            time: specs.iter().any(InletSpec::time_dependent),
            specs,
            faces,
        };
        // the walls
        let walls = Walls::new(mesh, geo, &case.boundary);
        let var = |dev: &mut dyn FlowDevice, omega: bool| -> Result<Var> {
            let (points, a) = kinds(mesh, &case.boundary, treatment, omega);
            let lsq = crate::simplec::disc_with(mesh, geo, &points)?.lsq;
            let lsq = Buf3 {
                x: up(dev, &lsq.x)?,
                y: up(dev, &lsq.y)?,
                z: up(dev, &lsq.z)?,
            };
            let bc_a = up(dev, &a)?;
            let psi = dev.alloc(nc)?;
            dev.fill(psi, 1.0)?;
            Ok(Var {
                phi: dev.alloc(nc)?,
                grad: Buf3 {
                    x: dev.alloc(nc)?,
                    y: dev.alloc(nc)?,
                    z: dev.alloc(nc)?,
                },
                psi,
                bc_a,
                bc_b: dev.alloc(nb)?,
                b_in: dev.alloc(nb)?,
                b_wall: dev.alloc(nb)?,
                lsq,
                old: if transient {
                    Some((dev.alloc(nc)?, dev.alloc(nc)?, dev.alloc(nc)?))
                } else {
                    None
                },
                valid: None,
                floor: 0.0,
            })
        };
        let k = var(dev, false)?;
        let w = var(dev, true)?;
        let ka_host = kinds(mesh, &case.boundary, treatment, false).1;
        let a3 = |dev: &mut dyn FlowDevice| -> Result<Buf3> {
            Ok(Buf3 {
                x: dev.alloc(nc)?,
                y: dev.alloc(nc)?,
                z: dev.alloc(nc)?,
            })
        };
        let gu = [a3(dev)?, a3(dev)?, a3(dev)?];
        let mu_lam = dev.alloc(nf)?;
        dev.copy(face_viscosity, mu_lam)?;
        let mask_b = {
            let b = dev.alloc(nb)?;
            dev.upload(b, &walls.mask_b)?;
            b
        };
        let wd = if walls.is_empty() {
            None
        } else {
            let nw = walls.len();
            let wd = WallDev {
                normal: Buf3 {
                    x: up(dev, &walls.normal[0])?,
                    y: up(dev, &walls.normal[1])?,
                    z: up(dev, &walls.normal[2])?,
                },
                y: up(dev, &walls.y)?,
                ks: up(dev, &walls.ks)?,
                weight: up(dev, &walls.weight)?,
                u_tau: dev.alloc(nw)?,
                y_plus: dev.alloc(nw)?,
                omega: dev.alloc(nw)?,
                blending: dev.alloc(nw)?,
                pk: dev.alloc(nw)?,
            };
            dev.wall_setup(&walls.stat)?;
            Some(wd)
        };
        let fix_mask = up(dev, &vec![0.0; nc])?;
        let pressure: Vec<u32> = bc_arrays
            .pb_rule
            .iter()
            .enumerate()
            .filter(|(_, r)| **r != foehncraft_kernels::flow::PB_EXTRAPOLATE)
            .map(|(fb, _)| fb as u32)
            .collect();
        let pressure_owners = pressure
            .iter()
            .map(|&fb| mesh.owner[ni + fb as usize])
            .collect();
        let mut t = Self {
            section,
            settings,
            c,
            o,
            treatment: match treatment {
                WallTreatment::LowRe => kwall::LOW_RE,
                WallTreatment::Automatic => kwall::AUTOMATIC,
            },
            nu,
            rho,
            nc,
            ni,
            nb,
            k,
            w,
            nu_t: dev.alloc(nc)?,
            nu_t_prev: None,
            f1: dev.alloc(nc)?,
            pk: dev.alloc(nc)?,
            dk: dev.alloc(nc)?,
            dw: dev.alloc(nc)?,
            capped: dev.alloc(nc)?,
            src: [
                dev.alloc(nc)?,
                dev.alloc(nc)?,
                dev.alloc(nc)?,
                dev.alloc(nc)?,
            ],
            d: dev.alloc(nc)?,
            gu,
            mu_lam,
            gamma_k: dev.alloc(nf)?,
            gamma_w: dev.alloc(nf)?,
            mask_b,
            walls,
            wd,
            fix_mask,
            fix_val: dev.alloc(nc)?,
            wall_pk: up(dev, &vec![0.0; nc])?,
            zero_face: dev.alloc(nf)?,
            cp: dev.alloc(nf)?,
            cn: dev.alloc(nf)?,
            bf: dev.alloc(nf)?,
            val: dev.alloc(stored)?,
            rhs: dev.alloc(nc)?,
            diag: dev.alloc(nc)?,
            order: dev.alloc(nc)?,
            work: dev.alloc(nc)?,
            inlets,
            ka_host,
            kin_host: vec![0.0; nb],
            pressure_faces: pressure,
            pressure_owners,
            added: None,
            t_set: f64::NAN,
            wall_distance: Vec::new(),
            health: TurbHealth::default(),
        };
        if t.settings.viscosity_relaxation < 1.0 {
            t.nu_t_prev = Some(dev.alloc(nc)?);
        }
        // the wall distance (D48)
        let d = if t.walls.is_empty() {
            vec![1e10; nc]
        } else {
            match t.section.wall_distance {
                WallDistanceMethod::Exact => crate::derived::wall_distance(case, mesh, geo)?,
                WallDistanceMethod::Pde => t.poisson_distance(dev, mesh, geo, disc, case)?,
            }
        };
        dev.upload(t.d, &d)?;
        t.wall_distance = d;
        // the start values (specification 7.3)
        let (k0, w0) = t.initial_values(dev, case, geo, bc_arrays, &env)?;
        dev.upload(t.k.phi, &k0)?;
        dev.upload(t.w.phi, &w0)?;
        let kref = k0
            .iter()
            .copied()
            .fold(0.0, f64::max)
            .max(t.inlet_k_max(bc_arrays));
        let kref = if kref > 0.0 { kref } else { 1e-6 };
        t.k.floor = 1e-10 * kref;
        t.w.floor = t.k.floor / (1e-2 * nu);
        Ok(t)
    }

    /// The largest inlet k (at the prescribed boundary velocity).
    fn inlet_k_max(&self, bc: &BoundaryArrays) -> f64 {
        self.inlets
            .faces
            .iter()
            .map(|&(fb, s, x)| {
                let fb = fb as usize;
                let u = [bc.ub[0][fb], bc.ub[1][fb], bc.ub[2][fb]];
                let sp = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
                self.inlets.specs[s]
                    .at(sp, x, 0.0, self.nu, self.c.beta_star)
                    .0
            })
            .fold(0.0, f64::max)
    }

    /// The start values: `[initial]`, else the values of the first inlet
    /// (by name) at its mean speed, else I = 0.01 and r = 1 with the
    /// initial speed (1 m/s when that is 0).
    fn initial_values(
        &self,
        _dev: &mut dyn FlowDevice,
        case: &crate::case::FlowCase,
        geo: &Geometry,
        bc: &BoundaryArrays,
        env: &Env,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        let nc = self.nc;
        let i = &case.initial;
        let u0 = i.velocity;
        let speed0 = (u0[0] * u0[0] + u0[1] * u0[1] + u0[2] * u0[2]).sqrt();
        let speed0 = if speed0 > 0.0 { speed0 } else { 1.0 };
        let cells = &geo.cell_centroid[..nc];
        if let (Some(k), Some(w)) = (&i.k, &i.omega) {
            let k = expr::value(k, env, Vars::SPACE, "initial.k")?;
            let w = expr::value(w, env, Vars::SPACE, "initial.omega")?;
            let kv: Vec<f64> = cells.iter().map(|&x| k.at(x, 0.0)).collect();
            let wv: Vec<f64> = cells.iter().map(|&x| w.at(x, 0.0)).collect();
            return Ok((kv, wv));
        }
        let from = |ii: f64, r: f64, sp: f64| {
            let k = 1.5 * (ii * sp) * (ii * sp);
            (k, k / (r * self.nu))
        };
        let (k, w) = if i.turbulence_intensity.is_some() || i.viscosity_ratio.is_some() {
            from(
                i.turbulence_intensity.unwrap_or(0.01),
                i.viscosity_ratio.unwrap_or(1.0),
                speed0,
            )
        } else if let Some(&(_, s, _)) = self.inlets.faces.first() {
            // the first inlet: area-free mean of its face values
            let (mut ks, mut ws, mut n) = (0.0, 0.0, 0.0);
            for &(fb, si, x) in &self.inlets.faces {
                if si != s {
                    continue;
                }
                let fb = fb as usize;
                let u = [bc.ub[0][fb], bc.ub[1][fb], bc.ub[2][fb]];
                let sp = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
                let sp = if sp > 0.0 { sp } else { speed0 };
                let (kk, ww) = self.inlets.specs[s].at(sp, x, 0.0, self.nu, self.c.beta_star);
                ks += kk;
                ws += ww;
                n += 1.0;
            }
            (ks / n, ws / n)
        } else {
            from(0.01, 1.0, speed0)
        };
        Ok((vec![k; nc], vec![w; nc]))
    }

    /// The Poisson wall distance (Tucker 2003): ∇²φ = −1, φ = 0 on the
    /// walls, zero normal gradient elsewhere, solved with the scalar
    /// machinery (non-orthogonal corrections in sweeps).
    fn poisson_distance(
        &mut self,
        dev: &mut dyn FlowDevice,
        mesh: &Mesh,
        geo: &Geometry,
        disc: &foehncraft_fv::disc::Discretisation,
        case: &crate::case::FlowCase,
    ) -> Result<Vec<f64>> {
        let nc = self.nc;
        let ni = disc.n_interior;
        let mut points = BTreeMap::new();
        let mut a = vec![1.0; self.nb];
        for p in &mesh.boundaries {
            let wall = matches!(
                case.boundary.get(&p.name),
                Some(BoundaryCondition::Wall { .. })
            );
            points.insert(
                p.name.clone(),
                if wall {
                    BoundaryKind::FixedValue
                } else {
                    BoundaryKind::ZeroGradient
                },
            );
            if wall {
                for f in p.range() {
                    a[f - ni] = 0.0;
                }
            }
        }
        let lsq = crate::simplec::disc_with(mesh, geo, &points)?.lsq;
        let lsq = Buf3 {
            x: up(dev, &lsq.x)?,
            y: up(dev, &lsq.y)?,
            z: up(dev, &lsq.z)?,
        };
        let bc_a = up(dev, &a)?;
        let bc_b = up(dev, &vec![0.0; self.nb])?;
        let gamma = up(dev, &vec![1.0; disc.n_faces])?;
        let su = up(dev, &vec![1.0 / self.rho; nc])?;
        let phi = dev.alloc(nc)?;
        let psi = dev.alloc(nc)?;
        dev.fill(psi, 1.0)?;
        let g = Buf3 {
            x: dev.alloc(nc)?,
            y: dev.alloc(nc)?,
            z: dev.alloc(nc)?,
        };
        for _ in 0..8 {
            dev.gradient_with(lsq, phi, bc_a, bc_b, None, g)?;
            dev.scalar_faces(&ScalarFaceArgs {
                bc_a,
                bc_b,
                jump: None,
                mass_flux: self.zero_face,
                gamma,
                phi,
                grad: g,
                psi,
                blend: 0.0,
                central: 0.0,
                beta: 0.1,
                nonorth: 1.0,
                cp: self.cp,
                cn: self.cn,
                bf: self.bf,
            })?;
            dev.scalar_rows(&ScalarRowArgs {
                cp: self.cp,
                cn: self.cn,
                bf: self.bf,
                mass_flux: self.zero_face,
                phi,
                su: Some(su),
                sp: None,
                alpha: 1.0,
                time: None,
                fix: None,
                positive: false,
                val: self.val,
                rhs: self.rhs,
                diag: self.diag,
                order: self.order,
            })?;
            dev.solve(
                self.val,
                self.rhs,
                phi,
                Preconditioner::Jacobi,
                false,
                1e-12,
                2000,
            )?;
        }
        dev.gradient_with(lsq, phi, bc_a, bc_b, None, g)?;
        let p = dev.download(phi)?;
        let (gx, gy, gz) = (dev.download(g.x)?, dev.download(g.y)?, dev.download(g.z)?);
        Ok(wall::poisson_distance(&p, [&gx, &gy, &gz]))
    }

    /// Sets the wall distance (verification: a manufactured solution with
    /// its own distance field).
    pub fn set_wall_distance(&mut self, dev: &mut dyn FlowDevice, d: &[f64]) -> Result<()> {
        anyhow::ensure!(d.len() == self.nc, "wall distance: one value per cell");
        dev.upload(self.d, d)?;
        self.wall_distance = d.to_vec();
        Ok(())
    }

    /// Adds sources of k and ω per unit mass (m²/s³ and 1/s²; a
    /// manufactured solution), split by sign into the rows.
    pub fn set_added_sources(
        &mut self,
        dev: &mut dyn FlowDevice,
        sk: &[f64],
        sw: &[f64],
    ) -> Result<()> {
        let (a, b) = match self.added {
            Some(x) => x,
            None => (dev.alloc(self.nc)?, dev.alloc(self.nc)?),
        };
        dev.upload(a, sk)?;
        dev.upload(b, sw)?;
        self.added = Some((a, b));
        Ok(())
    }

    /// Sets k and ω (verification and tests).
    pub fn set_fields(&mut self, dev: &mut dyn FlowDevice, k: &[f64], w: &[f64]) -> Result<()> {
        dev.upload(self.k.phi, k)?;
        dev.upload(self.w.phi, w)
    }

    /// The inlet values of k and ω for the boundary velocity `ub` at time
    /// `t` (only when they depend on the speed or on t, unless `all`).
    fn set_inlets(&mut self, dev: &mut dyn FlowDevice, ub: Buf3, t: f64, all: bool) -> Result<()> {
        if self.inlets.faces.is_empty() {
            return Ok(());
        }
        if !all && !self.inlets.speed && !(self.inlets.time && t != self.t_set) {
            return Ok(());
        }
        let idx: Vec<u32> = self.inlets.faces.iter().map(|f| f.0).collect();
        let (ux, uy, uz) = if self.inlets.speed {
            (
                dev.gather(ub.x, &idx)?,
                dev.gather(ub.y, &idx)?,
                dev.gather(ub.z, &idx)?,
            )
        } else {
            let z = vec![0.0; idx.len()];
            (z.clone(), z.clone(), z)
        };
        let mut kb = vec![0.0; self.nb];
        let mut wb = vec![0.0; self.nb];
        for (i, &(fb, s, x)) in self.inlets.faces.iter().enumerate() {
            let sp = (ux[i] * ux[i] + uy[i] * uy[i] + uz[i] * uz[i]).sqrt();
            let (k, w) = self.inlets.specs[s].at(sp, x, t, self.nu, self.c.beta_star);
            kb[fb as usize] = k.max(self.k.floor);
            wb[fb as usize] = w.max(self.w.floor);
        }
        dev.upload(self.k.b_in, &kb)?;
        dev.upload(self.w.b_in, &wb)?;
        self.kin_host = kb;
        self.t_set = t;
        Ok(())
    }

    /// bc_b = inlet values + wall values.
    fn boundary_values(&self, dev: &mut dyn FlowDevice) -> Result<()> {
        for v in [&self.k, &self.w] {
            dev.lincomb(v.bc_b, 1.0, v.b_in, 1.0, v.b_wall)?;
        }
        Ok(())
    }

    /// Gradients and limiters of k and ω.
    fn gradients(&self, dev: &mut dyn FlowDevice) -> Result<()> {
        let s = &self.settings;
        for v in [&self.k, &self.w] {
            dev.gradient_with(v.lsq, v.phi, v.bc_a, v.bc_b, None, v.grad)?;
            if s.scheme != Scheme::Upwind && s.limiter != LimiterKind::None {
                let (mn, mx) = dev.min_max(v.phi)?;
                let e = s.limiter_k.max(1e-12) * (mx - mn);
                let eps2 = Eps2::Single((e * e).max(1e-300));
                dev.limiter_with(v.phi, v.bc_a, v.bc_b, None, v.grad, eps2, v.psi)?;
            }
        }
        Ok(())
    }

    /// The model per cell: with `sources`, also the source split.
    fn cells(&self, dev: &mut dyn FlowDevice, sources: bool) -> Result<()> {
        dev.sst_cells(&SstCellOps {
            grad_u: self.gu,
            k: self.k.phi,
            w: self.w.phi,
            grad_k: self.k.grad,
            grad_w: self.w.grad,
            d: self.d,
            nu: self.nu,
            c: self.c,
            o: self.o,
            add: if sources { self.added } else { None },
            wall: (sources && self.treatment == kwall::AUTOMATIC && self.wd.is_some())
                .then_some((self.fix_mask, self.wall_pk)),
            nu_t: self.nu_t,
            f1: self.f1,
            pk: self.pk,
            sources: sources.then_some(self.src),
            dk: self.dk,
            dw: self.dw,
            capped: self.capped,
        })
    }

    /// The wall treatment: u_τ, y⁺, μ_w, the wall values of rough walls and
    /// the ω of the near-wall cells.
    fn walls(&self, dev: &mut dyn FlowDevice, u: Buf3, ub: Buf3) -> Result<()> {
        let Some(wd) = self.wd else { return Ok(()) };
        dev.wall_faces(&WallFaceOps {
            normal: wd.normal,
            y: wd.y,
            ks: wd.ks,
            u,
            ub,
            nu: self.nu,
            treatment: self.treatment,
            c: self.c,
            u_tau: wd.u_tau,
            y_plus: wd.y_plus,
            omega: wd.omega,
            blending: wd.blending,
            pk: wd.pk,
            mu_eff: dev.face_viscosity(),
            k_b: self.k.b_wall,
            w_b: self.w.b_wall,
        })?;
        dev.wall_cells(&WallCellOps {
            omega: wd.omega,
            weight: wd.weight,
            blending: wd.blending,
            pk: wd.pk,
            wall_pk: self.wall_pk,
            fix_val: self.fix_val,
            fix_mask: self.fix_mask,
        })
    }

    /// μ_eff on the faces off the walls and the diffusivities of k and ω.
    fn faces(&self, dev: &mut dyn FlowDevice) -> Result<()> {
        dev.turb_faces(&TurbFaceOps {
            nu_t: self.nu_t,
            dk: self.dk,
            dw: self.dw,
            mu: self.mu_lam,
            wall: self.mask_b,
            mu_eff: dev.face_viscosity(),
            gamma_k: self.gamma_k,
            gamma_w: self.gamma_w,
        })
    }

    /// The velocity gradient of `u` with the boundary values `ub`.
    fn velocity_gradient(&self, dev: &mut dyn FlowDevice, u: Buf3, ub: Buf3) -> Result<()> {
        for ((ui, bi), g) in u.all().into_iter().zip(ub.all()).zip(self.gu) {
            dev.gradient(Lsq::Velocity, ui, Mask::Velocity, bi, g)?;
        }
        Ok(())
    }

    /// The face viscosity of the start (before the first iteration and after
    /// a restart without turbulence state).
    pub(crate) fn start(
        &mut self,
        dev: &mut dyn FlowDevice,
        u: Buf3,
        ub: Buf3,
        bc: &BoundaryArrays,
    ) -> Result<()> {
        self.set_inlets(dev, ub, 0.0, true)?;
        self.velocity_gradient(dev, u, ub)?;
        self.walls(dev, u, ub)?;
        self.boundary_values(dev)?;
        self.gradients(dev)?;
        self.cells(dev, false)?;
        self.faces(dev)?;
        self.modified_pressure(dev, bc)
    }

    /// The pressure boundaries take p̂ = p + ⅔ρk with k of the face.
    fn modified_pressure(&self, dev: &mut dyn FlowDevice, bc: &BoundaryArrays) -> Result<()> {
        if self.pressure_faces.is_empty() {
            return Ok(());
        }
        let kp = dev.gather(self.k.phi, &self.pressure_owners)?;
        let mut v = bc.pset.clone();
        for (i, &fb) in self.pressure_faces.iter().enumerate() {
            let fb = fb as usize;
            let kf = self.ka_host[fb] * kp[i] + self.kin_host[fb];
            v[fb] += 2.0 / 3.0 * self.rho * kf;
        }
        dev.set_pressure_values(&v)
    }

    /// The turbulence solves of one outer iteration with the corrected
    /// velocity `u` and the fluxes `mass_flux`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve(
        &mut self,
        dev: &mut dyn FlowDevice,
        u: Buf3,
        ub: Buf3,
        mass_flux: Buf,
        bc: &BoundaryArrays,
        step: Option<StepTime>,
    ) -> Result<TurbResiduals> {
        let tr = step.is_some();
        let t = step.map_or(0.0, |s| s.time);
        self.set_inlets(dev, ub, t, false)?;
        // 1. velocity gradient of the corrected velocity, the wall
        // treatment, gradients of k and ω
        self.velocity_gradient(dev, u, ub)?;
        self.walls(dev, u, ub)?;
        self.boundary_values(dev)?;
        self.gradients(dev)?;
        // 2. the sources and the diffusivities from the fields at the start
        self.cells(dev, true)?;
        self.faces(dev)?;
        let s = self.settings.clone();
        let alpha = s.alpha(tr);
        let mut out = TurbResiduals::default();
        for (e, (v, gamma)) in [(self.k, self.gamma_k), (self.w, self.gamma_w)]
            .into_iter()
            .enumerate()
        {
            dev.scalar_faces(&ScalarFaceArgs {
                bc_a: v.bc_a,
                bc_b: v.bc_b,
                jump: None,
                mass_flux,
                gamma,
                phi: v.phi,
                grad: v.grad,
                psi: v.psi,
                blend: if s.scheme == Scheme::Upwind { 0.0 } else { 1.0 },
                central: if s.scheme == Scheme::Central {
                    1.0
                } else {
                    0.0
                },
                beta: 0.1,
                nonorth: 1.0,
                cp: self.cp,
                cn: self.cn,
                bf: self.bf,
            })?;
            let time = match (step, v.old) {
                (Some(st), Some((o1, o2, _))) => Some(ScalarTime {
                    inv_dt: 1.0 / st.dt,
                    a: st.a,
                    old: o1,
                    old2: o2,
                    bounds: (0.0, f64::INFINITY),
                }),
                _ => None,
            };
            let (su, sp) = (self.src[2 * e], self.src[2 * e + 1]);
            dev.scalar_rows(&ScalarRowArgs {
                cp: self.cp,
                cn: self.cn,
                bf: self.bf,
                mass_flux,
                phi: v.phi,
                su: Some(su),
                sp: Some(sp),
                alpha,
                time,
                fix: (e == 1
                    && self.treatment == kwall::AUTOMATIC
                    && self.wd.is_some()
                    && !self.walls.stat.cells.is_empty())
                .then_some((self.fix_mask, self.fix_val)),
                positive: true,
                val: self.val,
                rhs: self.rhs,
                diag: self.diag,
                order: self.order,
            })?;
            let scale = dev.scalar_scale(v.phi, self.rhs, self.diag, alpha)?;
            let (it, r0) = dev.solve(
                self.val,
                self.rhs,
                v.phi,
                Preconditioner::Jacobi,
                false,
                s.tol(tr),
                s.max_iter,
            )?;
            out.residual.push(if scale > 0.0 { r0 / scale } else { r0 });
            out.iterations.push(it);
            out.floors
                .push(dev.floor(v.phi, v.floor, f64::INFINITY, self.work)?);
        }
        // 3. ν_t from the new fields
        if let Some(prev) = self.nu_t_prev {
            dev.copy(self.nu_t, prev)?;
        }
        self.cells(dev, false)?;
        if let Some(prev) = self.nu_t_prev {
            let a = s.viscosity_relaxation;
            dev.lincomb(self.nu_t, a, self.nu_t, 1.0 - a, prev)?;
        }
        out.capped = dev.sum_sq(self.capped)?.round() as usize;
        let nc = self.n_cells();
        let mut flag = 0u8;
        for (e, &n) in out.floors.iter().enumerate() {
            if n * 1000 > nc {
                self.health.floor_iterations[e] += 1;
                flag |= 1 << e;
            }
            self.health.floor_cells_max[e] = self.health.floor_cells_max[e].max(n);
        }
        let h = &mut self.health;
        h.flags.push(flag);
        let tail = h.flags.len().div_ceil(10);
        for e in 0..2 {
            h.floor_iterations_late[e] = h.flags[h.flags.len() - tail..]
                .iter()
                .filter(|&&f| f & (1 << e) != 0)
                .count();
        }
        self.health.capped_max = self.health.capped_max.max(out.capped);
        self.faces(dev)?;
        self.modified_pressure(dev, bc)?;
        Ok(out)
    }

    /// Start of a time step (as the scalars: kⁿ⁻¹ ← kⁿ ← k, with the
    /// extrapolated start value of BDF2).
    pub fn step_start(
        &mut self,
        dev: &mut dyn FlowDevice,
        has_old: bool,
        extrapolate: Option<f64>,
    ) -> Result<()> {
        for v in [&self.k, &self.w] {
            let Some((o1, o2, keep)) = v.old else {
                continue;
            };
            dev.copy(o2, keep)?;
            if has_old {
                dev.copy(o1, o2)?;
            }
            dev.copy(v.phi, o1)?;
            if has_old && let Some(w) = extrapolate {
                dev.lincomb(v.phi, 1.0 + w, o1, -w, o2)?;
                dev.floor(v.phi, v.floor, f64::INFINITY, self.work)?;
            }
        }
        Ok(())
    }

    /// Takes the last step back exactly.
    pub fn step_back(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for v in [&self.k, &self.w] {
            let Some((o1, o2, keep)) = v.old else {
                continue;
            };
            dev.copy(o1, v.phi)?;
            dev.copy(o2, o1)?;
            dev.copy(keep, o2)?;
        }
        Ok(())
    }

    /// Keeps a device copy of k and ω as the last valid state (D63).
    pub fn keep_valid(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for v in [&mut self.k, &mut self.w] {
            let b = match v.valid {
                Some(b) => b,
                None => dev.alloc(self.nc)?,
            };
            dev.copy(v.phi, b)?;
            v.valid = Some(b);
        }
        Ok(())
    }

    /// Puts the last valid k and ω back.
    pub fn rewind(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for v in [&self.k, &self.w] {
            if let Some(b) = v.valid {
                dev.copy(b, v.phi)?;
            }
        }
        Ok(())
    }

    /// After a rewind of a time step: the old level of n−1 back into kⁿ.
    pub fn rewind_levels(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for v in [&self.k, &self.w] {
            if let Some((o1, o2, _)) = v.old {
                dev.copy(o2, o1)?;
            }
        }
        Ok(())
    }

    /// Host copy of the state; `p`, `pb` are the solver's p̂ (cells and
    /// boundary faces).
    pub(crate) fn state(
        &self,
        dev: &mut dyn FlowDevice,
        p: &[f64],
        pb: &[f64],
    ) -> Result<TurbState> {
        let ni = self.ni;
        let k = dev.download(self.k.phi)?;
        let omega = dev.download(self.w.phi)?;
        let k_b = dev.download(self.k.bc_b)?;
        let omega_b = dev.download(self.w.bc_b)?;
        let mut u_tau = vec![0.0; self.nb];
        let mut y_plus = vec![0.0; self.nb];
        if let Some(wd) = self.wd {
            let ut = dev.download(wd.u_tau)?;
            let yp = dev.download(wd.y_plus)?;
            for (i, &f) in self.walls.stat.faces.iter().enumerate() {
                u_tau[f as usize - ni] = ut[i];
                y_plus[f as usize - ni] = yp[i];
            }
        }
        Ok(TurbState {
            nu_t: dev.download(self.nu_t)?,
            k_old: match self.k.old {
                Some((o1, _, _)) => Some(dev.download(o1)?),
                None => None,
            },
            omega_old: match self.w.old {
                Some((o1, _, _)) => Some(dev.download(o1)?),
                None => None,
            },
            mu_eff: dev.download(dev.face_viscosity())?,
            k_b,
            omega_b,
            u_tau,
            y_plus,
            p_modified: p.to_vec(),
            pb_modified: pb.to_vec(),
            wall_distance: self.wall_distance.clone(),
            k,
            omega,
        })
    }

    /// k on the boundary faces, a k_P + b.
    pub(crate) fn k_face(&self, k: &[f64], k_b: &[f64], owner: &[u32]) -> Vec<f64> {
        (0..self.nb)
            .map(|fb| self.ka_host[fb] * k[owner[self.ni + fb] as usize] + k_b[fb])
            .collect()
    }

    /// Loads a host state (a restart).
    pub(crate) fn load(&mut self, dev: &mut dyn FlowDevice, s: &TurbState) -> Result<()> {
        dev.upload(self.k.phi, &s.k)?;
        dev.upload(self.w.phi, &s.omega)?;
        dev.upload(self.nu_t, &s.nu_t)?;
        if let (Some((o1, _, _)), Some(o)) = (self.k.old, &s.k_old) {
            dev.upload(o1, o)?;
        }
        if let (Some((o1, _, _)), Some(o)) = (self.w.old, &s.omega_old) {
            dev.upload(o1, o)?;
        }
        dev.upload(dev.face_viscosity(), &s.mu_eff)?;
        // inlet and wall parts of the boundary values
        let mut kin = vec![0.0; self.nb];
        let mut win = vec![0.0; self.nb];
        let mut kw = vec![0.0; self.nb];
        let mut ww = vec![0.0; self.nb];
        for &(fb, _, _) in &self.inlets.faces {
            let fb = fb as usize;
            kin[fb] = s.k_b[fb];
            win[fb] = s.omega_b[fb];
        }
        for fb in 0..self.nb {
            if self.walls.mask_b[fb] > 0.0 {
                kw[fb] = s.k_b[fb];
                ww[fb] = s.omega_b[fb];
            }
        }
        dev.upload(self.k.b_in, &kin)?;
        dev.upload(self.w.b_in, &win)?;
        dev.upload(self.k.b_wall, &kw)?;
        dev.upload(self.w.b_wall, &ww)?;
        self.kin_host = kin;
        self.boundary_values(dev)?;
        if let Some(wd) = self.wd {
            let ut: Vec<f64> = self
                .walls
                .stat
                .faces
                .iter()
                .map(|&f| s.u_tau[f as usize - self.ni])
                .collect();
            dev.upload(wd.u_tau, &ut)?;
        }
        Ok(())
    }

    /// The modified pressures to load into the solver.
    pub(crate) fn pressure_of(s: &TurbState) -> (&[f64], &[f64]) {
        (&s.p_modified, &s.pb_modified)
    }

    /// Per wall boundary: the area-weighted min, mean and max of y⁺ and
    /// the area fraction outside the claimed band (low-Re: above 2;
    /// automatic: the buffer layer 3 to 20), for the report and TURB001.
    pub fn y_plus_summary(&self, dev: &mut dyn FlowDevice, mesh: &Mesh) -> Result<Vec<YPlusRow>> {
        let Some(wd) = self.wd else {
            return Ok(Vec::new());
        };
        let yp = dev.download(wd.y_plus)?;
        let mut per: BTreeMap<usize, (f64, f64, f64, f64, f64)> = BTreeMap::new();
        for (i, &b) in self.walls.boundary.iter().enumerate() {
            let a = self.walls.area[i];
            let e = per
                .entry(b)
                .or_insert((f64::INFINITY, 0.0, f64::NEG_INFINITY, 0.0, 0.0));
            e.0 = e.0.min(yp[i]);
            e.1 += a * yp[i];
            e.2 = e.2.max(yp[i]);
            e.3 += a;
            // area outside the claimed band
            // (automatic: the buffer layer, where the blended wall
            // values overpredict C_f of the flat plate by up to 16 %)
            let out = match self.treatment {
                kwall::LOW_RE => yp[i] > 2.0,
                _ => yp[i] > 3.0 && yp[i] < 20.0,
            };
            if out {
                e.4 += a;
            }
        }
        Ok(per
            .into_iter()
            .map(|(b, e)| {
                (
                    mesh.boundaries[b].name.clone(),
                    e.0,
                    e.1 / e.3,
                    e.2,
                    e.4 / e.3,
                )
            })
            .collect())
    }

    /// The model's name for the log.
    pub fn label(&self) -> String {
        format!(
            "k-omega SST (2003), {} wall treatment, {} production",
            match self.treatment {
                kwall::LOW_RE => "low-Re",
                _ => "automatic",
            },
            match self.o.limiter {
                ksst::LIMITER_KATO_LAUNDER => "Kato-Launder",
                ksst::LIMITER_DURBIN => "Durbin-bounded",
                ksst::LIMITER_NONE => "unlimited",
                _ => "Menter-limited",
            }
        )
    }

    /// Number of cells.
    pub fn n_cells(&self) -> usize {
        self.nc
    }

    /// The ν_t cap that acted in the last iteration, ν_t/ν max and the
    /// cell (for the extremes block).
    pub fn ratio_max(&self, dev: &mut dyn FlowDevice) -> Result<(f64, usize)> {
        let v = dev.download(self.nu_t)?;
        let (mut m, mut at) = (0.0, 0);
        for (c, &x) in v.iter().enumerate() {
            if x > m {
                m = x;
                at = c;
            }
        }
        Ok((m / self.nu, at))
    }
}

/// The first-cell y⁺ of one wall boundary.
#[derive(Debug, Clone, Serialize)]
pub struct WallYPlus {
    /// Boundary.
    pub boundary: String,
    /// Smallest.
    pub min: f64,
    /// Area-weighted mean.
    pub mean: f64,
    /// Largest.
    pub max: f64,
}

/// The turbulence model in the report of a solve.
#[derive(Debug, Clone, Serialize)]
pub struct TurbulenceSummary {
    /// Model and options.
    pub model: String,
    /// First-cell y⁺ per wall.
    pub walls: Vec<WallYPlus>,
    /// Largest ν_t/ν and its cell.
    pub viscosity_ratio_max: f64,
    /// The cell of `viscosity_ratio_max`.
    pub viscosity_ratio_cell: usize,
    /// Floor and cap counts of the run.
    pub health: TurbHealth,
}

/// The summary of the turbulence model of `fl`, the findings it gives
/// (TURB001) and the log lines.
pub fn summary(
    fl: &mut crate::simplec::Flow,
    mesh: &Mesh,
) -> Result<Option<(TurbulenceSummary, Vec<Diagnostic>)>> {
    let Some(t) = fl.turbulence.as_ref() else {
        return Ok(None);
    };
    let dev = fl.dev.as_mut();
    let walls = t.y_plus_summary(dev, mesh)?;
    let (rmax, cell) = t.ratio_max(dev)?;
    let mut finds = Vec::new();
    for (e, name) in EQUATIONS.iter().enumerate() {
        let n = t.health.floor_iterations_late[e];
        if n > 0 {
            finds.push(Diagnostic::new(
                diag::TURB_FLOOR,
                format!(
                    "the positivity floor of {name} acted on more than 0.1 % of the cells in \
                     {n} of the final tenth of the solves, at most {} cells in one (a \
                     converging run does not need it at its end: check the mesh near the \
                     boundary-layer edge and the inlet values)",
                    t.health.floor_cells_max[e]
                ),
            ));
        }
    }
    for (b, _, _, _, out) in &walls {
        if *out > 0.1 {
            let text = if t.treatment == kwall::LOW_RE {
                format!(
                    "wall `{b}`: first-cell y⁺ above 2 on {:.0} % of the area with the low-Re \
                     treatment (it needs y⁺ ≲ 1; the automatic treatment takes any y⁺)",
                    100.0 * out
                )
            } else {
                format!(
                    "wall `{b}`: first-cell y⁺ between 3 and 20 on {:.0} % of the area; the \
                     automatic treatment overpredicts the wall shear stress there by up to 16 % \
                     (flat plate, gallery 68): put the first cell at y⁺ ≲ 1 or y⁺ ≳ 25",
                    100.0 * out
                )
            };
            finds.push(Diagnostic::new(diag::TURB_Y_PLUS, text));
        }
    }
    Ok(Some((
        TurbulenceSummary {
            model: t.label(),
            walls: walls
                .into_iter()
                .map(|(boundary, min, mean, max, _)| WallYPlus {
                    boundary,
                    min,
                    mean,
                    max,
                })
                .collect(),
            viscosity_ratio_max: rmax,
            viscosity_ratio_cell: cell,
            health: t.health.clone(),
        },
        finds,
    )))
}

/// Adds the turbulence model of `case` to the flow solver `fl` (no-op for a
/// laminar case); call after [`crate::scalar::attach`] and before loading a
/// restart.
pub fn attach(
    fl: &mut crate::simplec::Flow,
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
) -> Result<()> {
    if case.turbulence.is_laminar() {
        return Ok(());
    }
    let fv = fl.dev.face_viscosity();
    let t = Turbulence::new(
        fl.dev.as_mut(),
        case,
        mesh,
        geo,
        &fl.disc,
        &fl.bc,
        fl.phys.rho,
        fl.phys.mu,
        fv,
    )?;
    fl.add_turbulence(t)
}

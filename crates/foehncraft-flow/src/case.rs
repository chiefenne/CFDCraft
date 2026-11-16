//! The `foehn flow` case: steady incompressible laminar flow on one mesh,
//! or an order-of-accuracy sweep over refinement levels.
//!
//! Case file (TOML, D39). Every key is optional unless noted, SI units
//! (D44), unknown keys and values are errors with a "did you mean" hint:
//!
//! ```toml
//! [mesh]                         # as in `foehn mesh` (file or generate)
//! [mesh.generate]
//! type = "hex"
//! dimension = 2
//! cells = [64, 64, 1]
//! [fluid]
//! density = 1.0                  # ρ (kg/m³)
//! viscosity = 0.01               # μ (Pa s), or an expression of x, y, z (D53)
//! # material = "water"           # a material instead (D94-D102, crate::materials)
//! # temperature = 293.15         # K, where the material is evaluated
//! [flow]
//! convection = true              # false: Stokes flow
//! body_force = [0.0, 0.0, 0.0]   # N/m³, numbers or expressions of x, y, z, t
//! [constants]                    # named numbers of the expressions (D53)
//! u_max = 1.5
//! [tables.ramp]                  # a table, called as ramp(t) (crate::expr)
//! x = [0.0, 1.0]
//! y = [0.0, 1.0]
//! [initial]
//! velocity = [0.0, 0.0, 0.0]
//! pressure = 0.0
//! [pressure]                     # only used when no boundary fixes the pressure
//! reference_point = [0.5, 0.5, 0.0]
//! reference_value = 0.0
//! [boundary.ymax]                # one table per patch (required), crate::bc
//! type = "wall"
//! velocity = [1.0, 0.0, 0.0]
//! [solver]                       # crate::simplec::Settings
//! relaxation_velocity = 0.9
//! [analytic]                     # optional exact solution, crate::analytic
//! type = "poiseuille"
//! [verification]
//! levels = [0, 1, 2]             # refine levels of [mesh.generate]: a sweep
//! expected_order = [1.8, 2.2]
//! max_error = 1e-9               # largest relative L2 velocity error
//! reference = "ghia-re100"       # centre-line data of the cavity
//! reference_tolerance = 0.02
//! [output]
//! vtkhdf = true
//! history_every = 1              # residual history: every N iterations
//! wall_csv = true                # one CSV of wall quantities per wall
//! derived = ["vorticity", "q-criterion"]  # derived fields (D72), crate::derived
//! [output.reference]             # reference values of pressure-coefficient
//! velocity = 1.0
//! [[output.lines]]               # velocity and pressure along a line
//! name = "vertical"
//! from = [0.5, 0.0, 0.0]
//! to = [0.5, 1.0, 0.0]
//! ```
//!
//! Output (D43): `flow.log` (header with what is written and how large),
//! `residuals.csv`, `patches.csv` (mass flow, mean pressure and force per
//! patch), `wall_<patch>.csv`, `line_<name>.csv`, `fields.vtkhdf`,
//! `case.resolved.toml` (the case with every default filled in), `flow.json`
//! and `flow.md` (the report), SVG charts.

use crate::analytic::{Analytic, Benchmark, Reference};
use crate::balance::{self, BalanceSummary};
use crate::bc::{self, BoundaryCondition, Level};
use crate::control::{self, Control};
use crate::diag;
use crate::divergence;
use crate::engine::Steer;
use crate::expr::{Env, Scalar, Value, Vars};
use crate::plot::{self, Series};
use crate::restart::{self, Restart};
use crate::simplec::{self, Flow, Physics, Residuals, Settings, Solution, Timings};
use crate::transient::TimeScheme;
use anyhow::{Context, Result, bail};
use foehncraft_mesh::case::{MeshSection, build_mesh};
use foehncraft_mesh::geometry::{dot, norm, sub};
use foehncraft_mesh::locate::Locator;
use foehncraft_mesh::renumber::renumber;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// `[fluid]`: the fluid's properties, as numbers (the shorthand of one
/// constant material, D100) or from a material (D94–D102).
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct FluidSection {
    /// Density ρ (kg/m³); default 1, or the material's. One number: a
    /// density that varies comes with the energy equation (L1).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub density: Option<Scalar>,
    /// Dynamic viscosity μ (Pa s): a number or an expression of x, y, z
    /// (D53); default 0.001, or the material's.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub viscosity: Option<Scalar>,
    /// A material of `[materials]` or of the database (`foehn materials`);
    /// `density` and `viscosity` given here replace its values.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub material: Option<String>,
    /// The temperature (K) at which the material's properties are
    /// evaluated (isothermal flow); default 293.15 K.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperature: Option<f64>,
    /// The operating pressure (Pa) of the ideal-gas density; default
    /// 101 325 Pa.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub operating_pressure: Option<f64>,
}

impl FluidSection {
    /// The shorthand of one constant fluid.
    pub fn constant(density: f64, viscosity: f64) -> Self {
        Self {
            density: Some(Scalar::Number(density)),
            viscosity: Some(Scalar::Number(viscosity)),
            ..Self::default()
        }
    }

    /// The conditions at which the properties are evaluated.
    pub fn conditions(&self) -> crate::materials::Conditions {
        crate::materials::Conditions {
            temperature: self
                .temperature
                .unwrap_or(crate::materials::DEFAULT_TEMPERATURE),
            pressure: self
                .operating_pressure
                .unwrap_or(crate::materials::DEFAULT_PRESSURE),
        }
    }
}

/// `[flow]`
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct FlowSection {
    /// Convection on (false: Stokes flow).
    pub convection: bool,
    /// Body force per unit volume (N/m³): three numbers, or expressions of
    /// x, y, z and, in a transient case, t (D53).
    pub body_force: crate::expr::Vector,
}

impl Default for FlowSection {
    fn default() -> Self {
        Self {
            convection: true,
            body_force: crate::expr::vector([0.0; 3]),
        }
    }
}

/// `[initial]`
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct InitialSection {
    /// Start velocity (m/s).
    pub velocity: [f64; 3],
    /// Start pressure (Pa).
    pub pressure: f64,
    /// Start from potential flow (D52) instead of `velocity`: the
    /// irrotational velocity that carries the inflow through the domain to
    /// the pressure boundaries, from one Laplace solve with the pressure
    /// multigrid.
    pub potential_flow: bool,
}

/// `[pressure]`
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct PressureSection {
    /// Point whose cell holds the reference pressure (default: the cell
    /// nearest to the lower corner of the bounding box).
    pub reference_point: Option<[f64; 3]>,
    /// Pressure in the reference cell (Pa).
    pub reference_value: f64,
}

/// `[run]`
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct RunSection {
    /// CPU threads (default: all cores; meshes below 20 000 cells run on
    /// one).
    pub threads: Option<usize>,
    /// Where the flow loop runs: `cpu` or `cuda` (a build with the `cuda`
    /// feature; the results are bitwise those of `cpu`, D31).
    pub backend: String,
    /// PTX module of the `cuda` backend.
    pub ptx: PathBuf,
    /// Restart file to continue from (D42).
    pub restart: Option<PathBuf>,
    /// Write `restart.h5` every N iterations or time steps (0: never).
    pub restart_every: usize,
    /// Also write a numbered `restart_<N>.h5` at these iterations or steps.
    pub restart_at: Vec<usize>,
    /// Transient: write `restart.h5` every this much simulated time (s), at
    /// the multiples of it; the time step is shortened to land on them
    /// (D66).
    pub restart_interval: Option<f64>,
    /// Transient: also write a numbered `restart_<N>.h5` at these simulated
    /// times (s), landed on exactly (D66).
    pub restart_times: Vec<f64>,
    /// Write `restart.h5` at the end of the run.
    pub restart_final: bool,
    /// Re-read the case file every N iterations or time steps (0: only on a
    /// `.read_toml_N` control file), D40.
    pub reread_every: usize,
    /// Limits of the run (D83).
    pub limits: Limits,
}

/// `[run.limits]` (D83): when one is reached the run writes a restart file
/// and the final output and ends with exit code 3 (stopped, D64).
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Limits {
    /// Largest wall time of the run (s), setup included.
    pub max_wall_time: Option<f64>,
    /// Largest number of time steps of the run (transient; a continued run
    /// counts from its restart).
    pub max_steps: Option<usize>,
    /// Largest number of outer iterations of the run (steady).
    pub max_iterations: Option<usize>,
    /// Largest volume of the output directory (bytes).
    pub max_output: Option<f64>,
}

impl Limits {
    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        if let Some(t) = self.max_wall_time {
            crate::need!(
                t > 0.0,
                CASE_RANGE,
                "[run.limits] max_wall_time must be > 0"
            );
        }
        if let Some(n) = self.max_steps {
            crate::need!(n >= 1, CASE_RANGE, "[run.limits] max_steps must be >= 1");
        }
        if let Some(n) = self.max_iterations {
            crate::need!(
                n >= 1,
                CASE_RANGE,
                "[run.limits] max_iterations must be >= 1"
            );
        }
        if let Some(b) = self.max_output {
            crate::need!(b > 0.0, CASE_RANGE, "[run.limits] max_output must be > 0");
        }
        Ok(())
    }

    /// True when no limit is set.
    pub fn is_empty(&self) -> bool {
        *self == Self::default()
    }
}

impl Default for RunSection {
    fn default() -> Self {
        Self {
            threads: None,
            backend: "cpu".into(),
            ptx: PathBuf::from("target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx"),
            restart: None,
            restart_every: 0,
            restart_at: Vec::new(),
            restart_interval: None,
            restart_times: Vec::new(),
            restart_final: true,
            reread_every: 0,
            limits: Limits::default(),
        }
    }
}

/// `[time]`: a transient case (D28).
#[derive(Debug, Clone, Copy, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct TimeSection {
    /// Time integration.
    pub scheme: TimeScheme,
    /// Time step Δt (s); with `adaptive` the first one.
    pub time_step: f64,
    /// End time (s).
    pub end_time: Option<f64>,
    /// Number of time steps (instead of, or as well as, `end_time`).
    pub steps: Option<usize>,
    /// Time of the initial fields (s).
    pub start_time: f64,
    /// Adapt Δt to the largest Courant number `max_cfl`.
    pub adaptive: bool,
    /// Target of the largest Courant number with `adaptive`.
    pub max_cfl: f64,
    /// Largest ratio of two successive time steps with `adaptive`.
    pub max_growth: f64,
    /// Smallest time step with `adaptive` (s).
    pub min_time_step: f64,
    /// Largest time step with `adaptive` (s).
    pub max_time_step: Option<f64>,
    /// Relative tolerance of the local time-integration error with
    /// `adaptive` (D65): the step follows from the error estimate of BDF2
    /// or Crank–Nicolson, bounded by `max_cfl` and `max_growth`, and a
    /// step whose error exceeds the tolerance is repeated with a shorter
    /// one. Off when absent.
    pub error_tolerance: Option<f64>,
    /// Absolute tolerance of the local time-integration error (m/s), added
    /// to `error_tolerance` × the RMS velocity.
    pub error_absolute: f64,
}

impl Default for TimeSection {
    fn default() -> Self {
        Self {
            scheme: TimeScheme::Bdf2,
            time_step: 0.0,
            end_time: None,
            steps: None,
            start_time: 0.0,
            adaptive: false,
            max_cfl: 1.0,
            max_growth: 1.2,
            min_time_step: 0.0,
            max_time_step: None,
            error_tolerance: None,
            error_absolute: 1e-6,
        }
    }
}

impl TimeSection {
    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        crate::need!(
            self.time_step > 0.0,
            CASE_RANGE,
            "[time] time_step must be > 0 (the first step with `adaptive`)"
        );
        crate::need!(
            self.end_time.is_some() || self.steps.is_some(),
            CASE_MISSING_KEY,
            "[time] needs `end_time` or `steps`"
        );
        if let Some(e) = self.end_time {
            crate::need!(
                e > self.start_time,
                CASE_RANGE,
                "[time] end_time must be later than start_time"
            );
        }
        crate::need!(self.max_cfl > 0.0, CASE_RANGE, "[time] max_cfl must be > 0");
        crate::diag::check_fix(
            self.max_growth > 1.0 && self.max_growth < 1.0 + std::f64::consts::SQRT_2,
            crate::diag::CASE_RANGE,
            &["time.max_growth=1.2"],
            || {
                String::from(
                    "[time] max_growth must lie in (1, 2.414): BDF2 with variable steps is zero-stable \
             only for step ratios below 1 + √2",
                )
            },
        )?;
        crate::need!(
            self.min_time_step >= 0.0,
            CASE_RANGE,
            "[time] min_time_step must be >= 0"
        );
        if let Some(m) = self.max_time_step {
            crate::need!(
                m >= self.min_time_step && m > 0.0,
                CASE_RANGE,
                "[time] max_time_step must be > 0 and >= min_time_step"
            );
        }
        if let Some(e) = self.error_tolerance {
            crate::need!(
                e > 0.0 && e < 1.0,
                CASE_RANGE,
                "[time] error_tolerance must lie in (0, 1)"
            );
            crate::need!(
                self.adaptive,
                CASE_RANGE,
                "[time] error_tolerance needs adaptive = true"
            );
            crate::need!(
                self.scheme != TimeScheme::Bdf1,
                CASE_RANGE,
                "[time] error_tolerance needs scheme bdf2 or crank-nicolson (the estimate is that \
                 of a second-order method)"
            );
        }
        crate::need!(
            self.error_absolute >= 0.0,
            CASE_RANGE,
            "[time] error_absolute must be >= 0"
        );
        Ok(())
    }
}

/// A point probe, `[[output.probes]]`.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Probe {
    /// Name (columns `<name>_u` … in `probes.csv`).
    pub name: String,
    /// The point.
    pub point: [f64; 3],
}

/// Force coefficients, `[output.coefficients]`: C = 2 F / (ρ U² A) of the
/// sum of the forces on `patches`.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Coefficients {
    /// Patches whose forces are summed.
    pub patches: Vec<String>,
    /// Reference velocity U (m/s).
    pub velocity: f64,
    /// Reference area A (m²).
    pub area: f64,
    /// Drag direction (unit vector).
    #[serde(default = "x_axis")]
    pub drag_direction: [f64; 3],
    /// Lift direction (unit vector).
    #[serde(default = "y_axis")]
    pub lift_direction: [f64; 3],
}

fn x_axis() -> [f64; 3] {
    [1.0, 0.0, 0.0]
}

fn y_axis() -> [f64; 3] {
    [0.0, 1.0, 0.0]
}

/// `[verification]`
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct VerificationSection {
    /// Refinement levels of `[mesh.generate]`; two or more give an
    /// order-of-accuracy sweep (needs `[analytic]`).
    pub levels: Vec<u32>,
    /// Accepted range of the observed order on the finest pair.
    pub expected_order: Option<[f64; 2]>,
    /// Largest accepted relative L2 velocity error (needs `[analytic]`).
    pub max_error: Option<f64>,
    /// Reference data of the lid-driven cavity.
    pub reference: Option<Reference>,
    /// Largest accepted deviation from the reference, in units of the lid
    /// speed.
    pub reference_tolerance: Option<f64>,
    /// Lid speed and cavity size for the reference comparison.
    pub lid_velocity: Option<f64>,
    /// Fails the case unless it converges.
    pub require_convergence: Option<bool>,
    /// Time steps of a time-order sweep on one mesh (transient, at least
    /// three, largest first): the order follows from the differences of
    /// successive solutions.
    pub time_steps: Vec<f64>,
    /// Accepted range of the observed time order on the smallest steps.
    pub expected_time_order: Option<[f64; 2]>,
    /// Benchmark data (the cylinder of Schäfer & Turek).
    pub benchmark: Option<Benchmark>,
    /// Patch of the cylinder of the benchmark.
    pub benchmark_patch: Option<String>,
    /// Largest accepted distance from the reference ranges, relative to
    /// the midpoint of each range.
    pub benchmark_tolerance: Option<f64>,
    /// Start of the evaluation window of a periodic benchmark (s; default:
    /// the second half of the run).
    pub evaluation_start: Option<f64>,
}

/// A line along which velocity and pressure are sampled.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Line {
    /// Name (file `line_<name>.csv`).
    pub name: String,
    /// Start point.
    pub from: [f64; 3],
    /// End point.
    pub to: [f64; 3],
    /// Number of points.
    #[serde(default = "default_points")]
    pub points: usize,
}

fn default_points() -> usize {
    101
}

/// `[output]`
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct OutputSection {
    /// Output directory (default `output/<case name>`).
    pub directory: Option<PathBuf>,
    /// Write `fields.vtkhdf` (velocity, pressure, error).
    pub vtkhdf: bool,
    /// Write a row of `residuals.csv` every N outer iterations.
    pub history_every: usize,
    /// Write `wall_<patch>.csv` for every wall.
    pub wall_csv: bool,
    /// Sample lines.
    pub lines: Vec<Line>,
    /// Transient: write `fields_<step>.vtkhdf` every N steps (0: none; the
    /// final fields are `fields.vtkhdf`).
    pub fields_every: usize,
    /// Transient: also write the fields at these steps.
    pub fields_at: Vec<usize>,
    /// Transient: write `fields_<step>.vtkhdf` every this much simulated
    /// time (s), at the multiples of it; the time step is shortened to
    /// land on them (D66).
    pub fields_interval: Option<f64>,
    /// Transient: also write the fields at these simulated times (s),
    /// landed on exactly (D66).
    pub fields_times: Vec<f64>,
    /// A row of `forces.csv` and `probes.csv` every N steps (transient).
    pub monitor_every: usize,
    /// Transient: a row of `forces.csv` and `probes.csv` every this much
    /// simulated time (s) instead of every `monitor_every` steps, landed
    /// on exactly (D66): equally spaced samples for spectra.
    pub monitor_interval: Option<f64>,
    /// Point probes: velocity and pressure per step (`probes.csv`).
    pub probes: Vec<Probe>,
    /// Force coefficients of some patches (`forces.csv`, the report).
    pub coefficients: Option<Coefficients>,
    /// Derived fields written with the fields (D72), from the cell
    /// gradients: vorticity, q-criterion, lambda2, enstrophy, strain-rate,
    /// pressure-coefficient, total-pressure, wall-shear-stress, y-plus.
    pub derived: Vec<crate::derived::Derived>,
    /// Reference values of the pressure coefficient.
    pub reference: crate::derived::ReferenceValues,
}

impl Default for OutputSection {
    fn default() -> Self {
        Self {
            directory: None,
            vtkhdf: true,
            history_every: 1,
            wall_csv: true,
            lines: Vec::new(),
            fields_every: 0,
            fields_at: Vec::new(),
            fields_interval: None,
            fields_times: Vec::new(),
            monitor_every: 1,
            monitor_interval: None,
            probes: Vec::new(),
            coefficients: None,
            derived: Vec::new(),
            reference: Default::default(),
        }
    }
}

/// A flow case.
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct FlowCase {
    /// The mesh.
    pub mesh: MeshSection,
    /// Fluid properties.
    pub fluid: FluidSection,
    /// Materials of the case (D94): new ones, or changes of database
    /// materials with `base`.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub materials: BTreeMap<String, crate::materials::MaterialSection>,
    /// Regions of cell zones with a kind and a material (D112); zones in no
    /// region are fluid with the `[fluid]` material.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub regions: BTreeMap<String, crate::materials::RegionSection>,
    /// Named constants of the expressions (D53): numbers, or expressions
    /// of other constants.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub constants: BTreeMap<String, Scalar>,
    /// Tables of the expressions (D53), piecewise linear in one argument.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub tables: BTreeMap<String, crate::expr::TableSection>,
    /// Flow model.
    pub flow: FlowSection,
    /// Start fields.
    pub initial: InitialSection,
    /// Pressure level.
    pub pressure: PressureSection,
    /// Boundary conditions per patch.
    pub boundary: BTreeMap<String, BoundaryCondition>,
    /// Solver settings.
    pub solver: Settings,
    /// Exact solution.
    pub analytic: Option<Analytic>,
    /// Checks.
    pub verification: VerificationSection,
    /// Threads, device, restart and re-read schedule.
    pub run: RunSection,
    /// Output.
    pub output: OutputSection,
    /// Time stepping; present in a transient case.
    pub time: Option<TimeSection>,
    /// Conservation ledger: report schedule and tolerances (D88, D89).
    pub balance: crate::balance::BalanceSection,
    /// The case file (for re-reading it during a run).
    #[serde(skip)]
    pub source: Option<PathBuf>,
    /// The `--set key=value` overrides applied to the case file.
    #[serde(skip)]
    pub overrides: Vec<String>,
}

/// Adds "did you mean" to the unknown-field and unknown-variant messages of
/// the TOML parser.
pub fn explain_toml_error(msg: &str) -> String {
    let mut out = msg.to_string();
    for key in ["unknown field `", "unknown variant `"] {
        let Some(i) = msg.find(key) else { continue };
        let rest = &msg[i + key.len()..];
        let Some(j) = rest.find('`') else { continue };
        let word = &rest[..j];
        let after = &rest[j + 1..];
        let end = after.find('\n').unwrap_or(after.len());
        let cands: Vec<&str> = after[..end].split('`').skip(1).step_by(2).collect();
        if let Some(c) = bc::did_you_mean(word, cands.iter().copied()) {
            out.push_str(&format!("\n  did you mean `{c}`?"));
        }
    }
    out
}

/// `run.restart = "latest"` (`--restart latest`): the newest restart file
/// of the output directory, by the step or iteration it holds (D64).
pub const LATEST: &str = "latest";

/// The diagnostic of a TOML error: the code from its kind, the location
/// and the key path from the span in `text` (D78).
pub fn toml_diagnostic(text: &str, e: &toml::de::Error) -> diag::Diagnostic {
    let msg = e.message();
    let code = if msg.contains("unknown field") {
        diag::CASE_UNKNOWN_KEY
    } else if msg.contains("missing field") {
        diag::CASE_MISSING_KEY
    } else if msg.contains("unknown variant")
        || msg.contains("invalid type")
        || msg.contains("invalid value")
        || msg.contains("invalid length")
    {
        diag::CASE_BAD_VALUE
    } else {
        diag::CASE_SYNTAX
    };
    let mut d = diag::Diagnostic::new(code, explain_toml_error(&e.to_string()));
    if let Some(span) = e.span()
        && span.start <= text.len()
    {
        let before = &text[..span.start];
        let line = before.matches('\n').count() + 1;
        let col = before.rsplit('\n').next().map_or(0, |l| l.chars().count()) + 1;
        d.location = Some([line, col]);
        // the table header above the span and the key at it
        let table = before
            .lines()
            .rev()
            .map(str::trim)
            .find(|l| l.starts_with('[') && !l.starts_with("[["))
            .map(|l| l.trim_matches(|c| c == '[' || c == ']').trim().to_string());
        let at = text[span.clone()].trim();
        let key: String = at
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric() || *c == '_' || *c == '-' || *c == '.')
            .collect();
        let line_text = before.rsplit('\n').next().unwrap_or("").trim();
        let is_key = !key.is_empty() && (line_text.is_empty() || !line_text.contains('='));
        d.key = match (table, is_key) {
            (Some(t), true) => Some(format!("{t}.{key}")),
            (None, true) => Some(key),
            (Some(t), false) => {
                let k = line_text.split('=').next().unwrap_or("").trim();
                Some(if k.is_empty() { t } else { format!("{t}.{k}") })
            }
            (None, false) => {
                let k = line_text.split('=').next().unwrap_or("").trim();
                (!k.is_empty()).then(|| k.to_string())
            }
        };
    }
    d
}

/// Applies one `key=value` override (`--set`, D62) to a case table: the key
/// is a dotted path (`fluid.viscosity`, `boundary.inlet.velocity`), the
/// value a TOML value (`1e-3`, `[1.0, 0.0, 0.0]`, `true`, `"rcm"`); a value
/// that is not TOML is taken as a string (`run.backend=cuda`).
pub fn apply_set(v: &mut toml::Table, set: &str) -> Result<()> {
    let bad = |m: String| diag::err(diag::Diagnostic::new(diag::CASE_SET, m));
    let Some((key, value)) = set.split_once('=') else {
        return Err(bad(format!("--set {set}: expected key=value")));
    };
    let key = key.trim();
    let value = value.trim();
    if key.is_empty() || key.split('.').any(str::is_empty) {
        return Err(bad(format!(
            "--set {set}: the key `{key}` is not a dotted path"
        )));
    }
    let parsed = toml::from_str::<toml::Table>(&format!("v = {value}"))
        .ok()
        .and_then(|mut t| t.remove("v"))
        .unwrap_or_else(|| toml::Value::String(value.to_string()));
    let mut parts: Vec<&str> = key.split('.').collect();
    let last = parts.pop().unwrap_or_default();
    let mut cur = v;
    let mut path = String::new();
    for p in parts {
        if !path.is_empty() {
            path.push('.');
        }
        path.push_str(p);
        let next = cur
            .entry(p.to_string())
            .or_insert_with(|| toml::Value::Table(Default::default()));
        let toml::Value::Table(t) = next else {
            return Err(bad(format!(
                "--set {set}: `{path}` is a value, not a table"
            )));
        };
        cur = t;
    }
    cur.insert(last.to_string(), parsed);
    Ok(())
}

impl FlowCase {
    /// Reads and validates a case file; relative paths are taken from its
    /// directory.
    pub fn from_file(path: &Path) -> Result<Self> {
        Self::from_file_with(path, &[])
    }

    /// As [`from_file`](Self::from_file), with `key=value` overrides
    /// applied before the case is validated (`--set`, D62). A re-read of
    /// the case file during the run applies them again.
    pub fn from_file_with(path: &Path, sets: &[String]) -> Result<Self> {
        let mut c = Self::from_file_unloaded(path, sets)?;
        c.load_tables()?;
        Ok(c)
    }

    /// As [`from_file_with`](Self::from_file_with), without reading the
    /// files the case refers to (the profile tables): the caller checks
    /// their paths first (read roots, D82) and then calls
    /// [`load_tables`](Self::load_tables).
    pub fn from_file_unloaded(path: &Path, sets: &[String]) -> Result<Self> {
        let text = std::fs::read_to_string(path).map_err(|e| {
            diag::err(diag::Diagnostic::new(
                diag::CASE_FILE,
                format!("reading {}: {e}", path.display()),
            ))
        })?;
        let mut c = Self::from_toml_parts(&text, path.parent().unwrap_or(Path::new(".")), sets)
            .with_context(|| format!("in case file {}", path.display()))?;
        c.source = Some(path.to_path_buf());
        Ok(c)
    }

    /// Reads the profile tables of the boundary conditions (their paths are
    /// resolved against the case file already).
    pub fn load_tables(&mut self) -> Result<()> {
        for b in self.boundary.values_mut() {
            b.load_table(Path::new(""))?;
        }
        for (name, t) in &mut self.tables {
            t.load(name)?;
        }
        // the expressions can be checked once the tables are read
        self.validate_expressions()?;
        Ok(())
    }

    /// Parses and validates a case from TOML text.
    pub fn from_toml(text: &str, dir: &Path) -> Result<Self> {
        Self::from_toml_with(text, dir, &[])
    }

    /// Parses a case from TOML text, applies the `key=value` overrides and
    /// validates it.
    pub fn from_toml_with(text: &str, dir: &Path, sets: &[String]) -> Result<Self> {
        let mut c = Self::from_toml_parts(text, dir, sets)?;
        c.load_tables()?;
        Ok(c)
    }

    /// As [`from_toml_with`](Self::from_toml_with), without reading the
    /// profile tables.
    fn from_toml_parts(text: &str, dir: &Path, sets: &[String]) -> Result<Self> {
        let mut c: Self = if sets.is_empty() {
            toml::from_str(text).map_err(|e| diag::err(toml_diagnostic(text, &e)))?
        } else {
            let mut v: toml::Table =
                toml::from_str(text).map_err(|e| diag::err(toml_diagnostic(text, &e)))?;
            for s in sets {
                apply_set(&mut v, s)?;
            }
            Self::deserialize(toml::Value::Table(v)).map_err(|e| {
                let mut d = toml_diagnostic("", &e);
                d.message = format!("after --set: {}", d.message);
                diag::err(d)
            })?
        };
        c.overrides = sets.to_vec();
        if let Some(f) = &c.mesh.file {
            c.mesh.file = Some(dir.join(f));
        }
        if let Some(r) = &c.run.restart
            && r.as_os_str() != LATEST
        {
            c.run.restart = Some(dir.join(r));
        }
        for b in c.boundary.values_mut() {
            b.resolve_paths(dir);
        }
        crate::expr::resolve_paths(&mut c.tables, dir);
        c.validate()?;
        Ok(c)
    }

    /// The constants and tables of the expressions.
    pub fn env(&self) -> Result<Env> {
        Env::new(&self.constants, &self.tables)
    }

    /// The resolved property program (D94 (3)): the `[fluid]` material and
    /// the regions, evaluated at the case conditions.
    pub fn properties(&self) -> Result<crate::materials::PropertySet> {
        use crate::materials::{RegionKind, resolve_material};
        let env = self.env()?;
        let cond = self.fluid.conditions();
        let mut warnings = Vec::new();
        let fluid = resolve_material(
            self.fluid.material.as_deref(),
            &self.materials,
            self.fluid.density.as_ref(),
            self.fluid.viscosity.as_ref(),
            &env,
            cond,
            &mut warnings,
        )?;
        let mut regions = Vec::new();
        let mut seen: BTreeMap<&str, &str> = BTreeMap::new();
        for (name, r) in &self.regions {
            crate::need!(
                r.kind == RegionKind::Fluid,
                REG_NOT_FLUID,
                "[regions.{name}] kind = \"{}\": only fluid regions are solved until the energy \
                 equation (L1) brings solids and porous zones (D112)",
                match r.kind {
                    RegionKind::Solid => "solid",
                    RegionKind::Porous => "porous",
                    RegionKind::Fluid => "fluid",
                }
            );
            crate::need!(
                !r.zones.is_empty(),
                REG_ZONE,
                "[regions.{name}] zones must name at least one cell zone"
            );
            for z in &r.zones {
                if let Some(other) = seen.insert(z.as_str(), name.as_str()) {
                    crate::need!(
                        false,
                        REG_ZONE,
                        "[regions.{name}] zone `{z}` is also in region `{other}`"
                    );
                }
            }
            let m = match &r.material {
                Some(mat) => resolve_material(
                    Some(mat),
                    &self.materials,
                    None,
                    None,
                    &env,
                    cond,
                    &mut warnings,
                )?,
                None => fluid.clone(),
            };
            let (a, b) = (m.density.value.constant(), fluid.density.value.constant());
            crate::need!(
                a.map(f64::to_bits) == b.map(f64::to_bits),
                MAT_FORM,
                "[regions.{name}] material `{}` has the density {:?} kg/m³, the fluid {:?}: \
                 the density is one number for the whole fluid until the energy equation (L1)",
                r.material.as_deref().unwrap_or(""),
                a.unwrap_or(f64::NAN),
                b.unwrap_or(f64::NAN)
            );
            regions.push((name.clone(), r.zones.clone(), m));
        }
        Ok(crate::materials::PropertySet {
            conditions: cond,
            fluid,
            regions,
            warnings,
        })
    }

    /// The density ρ (kg/m³).
    pub fn rho(&self) -> f64 {
        self.properties().map_or(f64::NAN, |p| p.density())
    }

    /// The reference viscosity μ (Pa s): the `[fluid]` value (Reynolds
    /// numbers, exact solutions).
    pub fn mu(&self) -> f64 {
        self.properties()
            .map_or(f64::NAN, |p| p.reference_viscosity())
    }

    /// The body force as compiled values.
    pub fn body_force_values(&self) -> Result<[Value; 3]> {
        let allowed = if self.time.is_some() {
            Vars::ALL
        } else {
            Vars::SPACE
        };
        crate::expr::vector_value(
            &self.flow.body_force,
            &self.env()?,
            allowed,
            "flow.body_force",
        )
    }

    /// The compiled boundary values (D53).
    pub fn compiled_boundaries(&self) -> Result<bc::Compiled> {
        bc::compile_all(&self.boundary, &self.env()?, self.time.is_some())
    }

    /// The uniform body force, when all three components are numbers.
    pub fn uniform_body_force(&self) -> Option<[f64; 3]> {
        crate::expr::numbers(&self.flow.body_force)
    }

    /// Checks every expression of the case: syntax, names, the variables
    /// each key allows (D53), and the property program (D94).
    pub fn validate_expressions(&self) -> Result<()> {
        let env = self.env()?;
        self.body_force_values()?;
        let transient = self.time.is_some();
        for (name, b) in &self.boundary {
            b.compile(name, &env, transient)?;
        }
        let p = self.properties()?;
        let rho = p.density();
        crate::need!(
            rho > 0.0 && rho.is_finite(),
            CASE_RANGE,
            "[fluid] density must be > 0"
        );
        if let Some(mu) = p.fluid.viscosity.value.constant() {
            crate::need!(
                mu > 0.0 && mu.is_finite(),
                CASE_RANGE,
                "[fluid] viscosity must be > 0"
            );
        }
        Ok(())
    }

    /// The resolved case (`case.resolved.toml`): every default filled in,
    /// the fluid's density and viscosity as resolved, and the resolved
    /// property program as comments (D94).
    pub fn resolved_toml(&self) -> Result<String> {
        let mut c = self.clone();
        let p = self.properties()?;
        if c.fluid.density.is_none() {
            c.fluid.density = Some(Scalar::Number(p.density()));
        }
        if c.fluid.viscosity.is_none()
            && let Some(mu) = p.fluid.viscosity.value.constant()
        {
            c.fluid.viscosity = Some(Scalar::Number(mu));
        }
        let mut text = toml::to_string_pretty(&c)?;
        text.push('\n');
        text.push_str(&p.text());
        Ok(text)
    }

    /// Ranges and consistency of the settings (the expressions and the
    /// file tables are checked by [`load_tables`](Self::load_tables)).
    pub fn validate(&self) -> Result<()> {
        if let Some(d) = &self.fluid.density
            && let Some(v) = d.number()
        {
            crate::need!(v > 0.0, CASE_RANGE, "[fluid] density must be > 0");
        }
        if let Some(d) = &self.fluid.viscosity
            && let Some(v) = d.number()
        {
            crate::need!(v > 0.0, CASE_RANGE, "[fluid] viscosity must be > 0");
        }
        if let Some(t) = self.fluid.temperature {
            crate::need!(
                t > 0.0 && t.is_finite(),
                CASE_RANGE,
                "[fluid] temperature must be > 0 (K)"
            );
        }
        if let Some(p) = self.fluid.operating_pressure {
            crate::need!(
                p > 0.0 && p.is_finite(),
                CASE_RANGE,
                "[fluid] operating_pressure must be > 0 (Pa)"
            );
        }
        if self
            .tables
            .values()
            .all(|t| t.file.is_none() || t.data.is_some())
        {
            self.validate_expressions()?;
        }
        let transient = self.time.is_some();
        self.solver.validate(transient)?;
        self.run.limits.validate()?;
        self.balance.validate()?;
        if let Some(t) = &self.time {
            t.validate()?;
        }
        crate::need!(
            self.output.history_every >= 1,
            CASE_RANGE,
            "[output] history_every must be >= 1"
        );
        crate::need!(
            self.output.monitor_every >= 1,
            CASE_RANGE,
            "[output] monitor_every must be >= 1"
        );
        // D66: output in simulated time
        let (o, r) = (&self.output, &self.run);
        for (key, v) in [
            ("output.fields_interval", o.fields_interval),
            ("output.monitor_interval", o.monitor_interval),
            ("run.restart_interval", r.restart_interval),
        ] {
            if let Some(v) = v {
                crate::need!(
                    v > 0.0 && v.is_finite(),
                    CASE_RANGE,
                    "{key} must be > 0 (s)"
                );
                crate::need!(
                    transient,
                    CASE_INCONSISTENT,
                    "{key} schedules output in simulated time and needs a [time] table"
                );
            }
        }
        for (key, v) in [
            ("output.fields_times", &o.fields_times),
            ("run.restart_times", &r.restart_times),
        ] {
            if !v.is_empty() {
                crate::need!(
                    v.iter().all(|x| x.is_finite()),
                    CASE_RANGE,
                    "{key} must hold finite times (s)"
                );
                crate::need!(
                    transient,
                    CASE_INCONSISTENT,
                    "{key} schedules output in simulated time and needs a [time] table"
                );
            }
        }
        if let Some(a) = &self.analytic
            && a.is_transient()
            && !transient
        {
            anyhow::bail!(
                "[analytic] type = \"taylor-green\" is transient and needs a [time] table"
            );
        }
        if let Some(c) = &self.output.coefficients {
            crate::need!(
                c.velocity > 0.0 && c.area > 0.0,
                CASE_RANGE,
                "[output.coefficients] velocity and area must be > 0"
            );
            crate::need!(
                !c.patches.is_empty(),
                CASE_INCONSISTENT,
                "[output.coefficients] patches must name at least one patch"
            );
        }
        crate::derived::validate(self)?;
        let v = &self.verification;
        if !v.time_steps.is_empty() {
            crate::need!(
                transient,
                CASE_INCONSISTENT,
                "[verification] time_steps need a [time] table"
            );
            crate::need!(
                v.time_steps.len() >= 3,
                CASE_INCONSISTENT,
                "[verification] time_steps need at least three steps for an observed order"
            );
            crate::need!(
                v.levels.len() < 2,
                CASE_INCONSISTENT,
                "[verification] a sweep is either over levels or over time_steps"
            );
            crate::diag::check_fix(
                self.time.is_some_and(|t| !t.adaptive),
                crate::diag::CASE_INCONSISTENT,
                &["time.adaptive=false"],
                || {
                    "[verification] time_steps need a fixed time step ([time] adaptive = false)"
                        .into()
                },
            )?;
        }
        if let Some(b) = v.benchmark {
            let re = self.rho() * b.mean_velocity() * Benchmark::DIAMETER / self.mu();
            crate::need!(
                (re - b.reynolds()).abs() < 1e-6 * b.reynolds(),
                CASE_INCONSISTENT,
                "[verification] benchmark {} needs Re = ρ Ū D / μ = {} (Ū = {} m/s, D = 0.1 m), \
                 the case has {re}",
                b.label(),
                b.reynolds(),
                b.mean_velocity()
            );
            if b == Benchmark::SchaferTurek2d2 {
                crate::need!(
                    transient,
                    CASE_INCONSISTENT,
                    "[verification] benchmark schafer-turek-2d2 is periodic and needs a [time] table"
                );
            }
        }
        let v = &self.verification;
        if v.levels.len() >= 2 {
            crate::need!(
                self.mesh.generate.is_some() || self.mesh.cylinder.is_some(),
                CASE_INCONSISTENT,
                "[verification] levels need a [mesh.generate] or [mesh.cylinder] table"
            );
            crate::need!(
                self.analytic.is_some() || v.benchmark.is_some(),
                CASE_INCONSISTENT,
                "[verification] levels need an [analytic] solution or a benchmark"
            );
        }
        if v.max_error.is_some() {
            crate::need!(
                self.analytic.is_some(),
                CASE_INCONSISTENT,
                "[verification] max_error needs an [analytic] solution"
            );
        }
        if let Some(reference) = v.reference {
            let re = self.rho() * v.lid_velocity.unwrap_or(1.0) / self.mu();
            let want = reference.reynolds();
            crate::need!(
                (re - want).abs() < 1e-6 * want,
                CASE_INCONSISTENT,
                "[verification] reference {:?} is for Re = {want}, the case has Re = ρ U / μ = {re} \
                 (unit cavity, lid_velocity = {})",
                reference,
                v.lid_velocity.unwrap_or(1.0)
            );
        }
        for l in &self.output.lines {
            crate::need!(
                l.points >= 2,
                CASE_RANGE,
                "[[output.lines]] {}: points must be >= 2",
                l.name
            );
        }
        Ok(())
    }

    pub(crate) fn analytic(&self) -> Option<Analytic> {
        self.analytic.map(|a| a.with_fluid(self.rho(), self.mu()))
    }

    /// What the boundary-condition check needs to know about the case.
    pub fn check_context(&self) -> bc::CheckContext {
        bc::CheckContext {
            laminar: true,
            transient: self.time.is_some(),
            potential_flow: self.initial.potential_flow,
        }
    }

    /// The exact solution at time `t` (transient solutions decay in time).
    pub fn analytic_at(&self, t: f64) -> Option<Analytic> {
        self.analytic().map(|a| a.at_time(t))
    }

    /// Reynolds number ρ U L / μ: the lid speed of the cavity reference,
    /// the scales of the exact solution or the mean velocity and the
    /// cylinder diameter of the benchmark.
    pub fn reynolds(&self) -> Option<f64> {
        let rho = self.rho();
        let mu = self.mu();
        if let Some(b) = self.verification.benchmark {
            return Some(rho * b.mean_velocity() * crate::analytic::Benchmark::DIAMETER / mu);
        }
        self.verification
            .reference
            .map(|_| rho * self.verification.lid_velocity.unwrap_or(1.0) / mu)
            .or_else(|| {
                self.analytic()
                    .map(|a| rho * a.velocity_scale() * a.length_scale() / mu)
            })
    }
}

/// Error norms against the exact solution.
#[derive(Debug, Clone, Serialize)]
pub struct ErrorNorms {
    /// Cells.
    pub cells: usize,
    /// Mesh size h (cube root of the mean cell volume; square root of the
    /// mean cell area in 2D).
    pub h: f64,
    /// L2 velocity error over the velocity scale.
    pub l2: f64,
    /// Largest velocity error over the velocity scale.
    pub linf: f64,
    /// Largest pressure error over the pressure range of the exact solution
    /// (absolute when the exact pressure is constant).
    pub pressure_linf: f64,
}

/// Comparison with the reference data.
#[derive(Debug, Clone, Serialize)]
pub struct ReferenceComparison {
    /// The data.
    pub reference: String,
    /// Largest |u − u_ref| / U on the vertical centre line.
    pub u_max_dev: f64,
    /// Largest |v − v_ref| / U on the horizontal centre line.
    pub v_max_dev: f64,
    /// (y, u/U) computed at the reference points.
    pub u: Vec<(f64, f64)>,
    /// (x, v/U) computed at the reference points.
    pub v: Vec<(f64, f64)>,
    /// (y, u/U) on 201 points of the vertical centre line, for the chart.
    #[serde(skip)]
    pub u_line: Vec<(f64, f64)>,
    /// (x, v/U) on 201 points of the horizontal centre line.
    #[serde(skip)]
    pub v_line: Vec<(f64, f64)>,
}

/// Forces and fluxes of one patch.
#[derive(Debug, Clone, Serialize)]
pub struct PatchSummary {
    /// Name.
    pub name: String,
    /// Condition.
    pub kind: String,
    /// Area (m²).
    pub area: f64,
    /// Mass flow out of the domain (kg/s).
    pub mass_flow: f64,
    /// Area-weighted mean pressure (Pa).
    pub mean_pressure: f64,
    /// Force of the fluid on the patch, pressure and viscous (N).
    pub force: [f64; 3],
}

/// One solve (one level).
#[derive(Debug, Clone, Serialize)]
pub struct Solve {
    /// Mesh description.
    pub mesh: String,
    /// SHA-256 of the mesh (D85).
    pub mesh_hash: String,
    /// Where it ran.
    pub device: String,
    /// Cells.
    pub cells: usize,
    /// Outer iterations.
    pub iterations: usize,
    /// Converged to the tolerance.
    pub converged: bool,
    /// Last residuals.
    pub residuals: Residuals,
    /// Phase timings.
    pub timings: Timings,
    /// Multigrid setups.
    pub amg_setups: usize,
    /// Nanoseconds per cell and outer iteration (R6).
    pub ns_per_cell_iteration: f64,
    /// Error norms.
    pub errors: Option<ErrorNorms>,
    /// Reference comparison.
    pub reference: Option<ReferenceComparison>,
    /// Patches.
    pub patches: Vec<PatchSummary>,
    /// Net mass flow out of the domain over the largest patch mass flow (in
    /// a closed domain over ρ |u|max times the largest patch area).
    pub mass_imbalance: f64,
    /// Force coefficients (drag, lift) of `[output.coefficients]` at the
    /// end.
    pub coefficients: Option<[f64; 2]>,
    /// Transient run.
    pub transient: Option<crate::unsteady::TransientSummary>,
    /// Benchmark comparison.
    pub benchmark: Option<crate::monitor::BenchmarkResult>,
    /// Conservation ledger (D88, D89).
    pub balance: Option<BalanceSummary>,
}

/// Report of a flow case.
#[derive(Debug, Clone, Serialize)]
pub struct FlowReport {
    /// Case name.
    pub name: String,
    /// Reynolds number ρ U L / μ when the case defines one (lid speed or
    /// analytic velocity scale and height), else `None`.
    pub reynolds: Option<f64>,
    /// Findings of the boundary-condition check.
    pub findings: Vec<bc::Finding>,
    /// Solves per level.
    pub solves: Vec<Solve>,
    /// Observed orders of the sweep.
    pub orders: Vec<f64>,
    /// Observed time orders of a time-step sweep.
    pub time_orders: Vec<f64>,
    /// The run was stopped early: `.stop`, a signal or a limit.
    pub stopped: bool,
    /// Why it stopped early.
    pub stop: Option<crate::engine::Stop>,
    /// Warnings of the run with their codes (re-reads, restarts).
    pub warnings: Vec<diag::Diagnostic>,
    /// Problems found by the checks.
    pub failures: Vec<String>,
    /// Seconds.
    pub seconds: f64,
}

pub(crate) fn mesh_for(case: &FlowCase, level: Option<u32>) -> Result<(Mesh, Geometry, String)> {
    let mut sec = case.mesh.clone();
    if let (Some(l), Some(g)) = (level, sec.generate.as_mut()) {
        g.refine = l;
    }
    if let (Some(l), Some(g)) = (level, sec.cylinder.as_mut()) {
        g.refine = l;
    }
    let (mut mesh, source, _) = build_mesh(&sec)?;
    renumber(&mut mesh, sec.renumber, sec.seed).map_err(anyhow::Error::msg)?;
    let geo = Geometry::compute(&mesh);
    let desc = format!("{source}: {} cells", mesh.n_cells);
    Ok((mesh, geo, desc))
}

pub(crate) fn dimension(case: &FlowCase) -> u8 {
    if case.mesh.cylinder.is_some() {
        return 2;
    }
    case.mesh.generate.as_ref().map_or(3, |g| g.dimension)
}

/// Writes log lines to the terminal callback and to `flow.log`.
pub(crate) struct Log<'a, F: FnMut(&str)> {
    pub(crate) out: &'a mut F,
    pub(crate) text: String,
}

impl<F: FnMut(&str)> Log<'_, F> {
    pub(crate) fn line(&mut self, l: &str) {
        (self.out)(l);
        self.text.push_str(l);
        self.text.push('\n');
    }
}

/// Velocity and pressure at `x` in cell `c`, reconstructed with the cell
/// gradients (second order).
pub(crate) fn value_at(fl: &Solution, geo: &Geometry, c: usize, x: [f64; 3]) -> ([f64; 3], f64) {
    let r = sub(x, geo.cell_centroid[c]);
    let u = fl.u.get(c);
    let g = [
        fl.grad_u[0].get(c),
        fl.grad_u[1].get(c),
        fl.grad_u[2].get(c),
    ];
    (
        [
            u[0] + dot(g[0], r),
            u[1] + dot(g[1], r),
            u[2] + dot(g[2], r),
        ],
        fl.p[c] + dot(fl.grad_p.get(c), r),
    )
}

pub(crate) fn locate(loc: &Locator, geo: &Geometry, x: [f64; 3]) -> usize {
    loc.locate_search(x).unwrap_or_else(|| {
        (0..geo.cell_volume.len())
            .min_by(|&a, &b| {
                let da = norm(sub(geo.cell_centroid[a], x));
                let db = norm(sub(geo.cell_centroid[b], x));
                da.total_cmp(&db)
            })
            .unwrap_or(0)
    })
}

/// Runs a case; `log` receives progress lines; files go to `out`.
pub fn run(case: &FlowCase, name: &str, out: &Path, log: impl FnMut(&str)) -> Result<FlowReport> {
    let mut steer = Steer::quiet(out);
    run_steered(case, name, out, &mut steer, log)
}

/// The steady restart file of the current state.
fn steady_restart(fl: &mut Flow, mesh: &Mesh, name: &str) -> Result<Restart> {
    Ok(Restart {
        case: name.into(),
        cell_id: mesh.cell_id[..mesh.n_cells].to_vec(),
        face_id: mesh.face_id.clone(),
        state: fl.state()?,
        iteration: fl.iteration(),
        amg_setups: fl.amg_setups,
        amg_matrix: fl.amg_matrix.clone(),
        transient: None,
        geometry: None,
    })
}

/// The physics of a case without a mesh: the density, the reference
/// viscosity and the uniform body force (zero when it is given by
/// expressions, which become a momentum source).
pub fn physics(case: &FlowCase) -> Result<Physics> {
    let p = case.properties()?;
    Ok(Physics::constant(
        p.density(),
        p.reference_viscosity(),
        case.uniform_body_force().unwrap_or([0.0; 3]),
        case.flow.convection,
    ))
}

/// `phys` on a mesh: the viscosity per face of the property program
/// (D94), when it is not one number.
pub fn physics_on(case: &FlowCase, phys: &Physics, mesh: &Mesh, geo: &Geometry) -> Result<Physics> {
    for (name, r) in &case.regions {
        for z in &r.zones {
            if !mesh.zones.contains(z) {
                let hint = match bc::did_you_mean(z, mesh.zones.iter().map(String::as_str)) {
                    Some(c) => format!(" (did you mean `{c}`?)"),
                    None => String::new(),
                };
                return Err(diag::err(
                    diag::Diagnostic::new(
                        diag::REG_ZONE,
                        format!(
                            "[regions.{name}] zone `{z}` is not a cell zone of the mesh{hint}; \
                             the mesh has {}",
                            mesh.zones
                                .iter()
                                .map(|z| format!("`{z}`"))
                                .collect::<Vec<_>>()
                                .join(", ")
                        ),
                    )
                    .key(format!("regions.{name}.zones")),
                ));
            }
        }
    }
    let mut p = phys.clone();
    p.mu_face = case.properties()?.face_viscosity(mesh, geo);
    Ok(p)
}

/// The momentum source per unit volume (N/m³) of a body force given by
/// expressions (D53), at the cell centroids at time `t`; `None` for a
/// uniform body force.
pub fn momentum_source(
    case: &FlowCase,
    geo: &Geometry,
    n_cells: usize,
    t: f64,
) -> Result<Option<[Vec<f64>; 3]>> {
    if case.uniform_body_force().is_some() {
        return Ok(None);
    }
    let v = case.body_force_values()?;
    let out: [Vec<f64>; 3] = std::array::from_fn(|i| {
        (0..n_cells)
            .map(|c| v[i].at(geo.cell_centroid[c], t))
            .collect()
    });
    if let Some((i, c)) = out
        .iter()
        .enumerate()
        .find_map(|(i, s)| s.iter().position(|x| !x.is_finite()).map(|c| (i, c)))
    {
        return Err(diag::err(
            diag::Diagnostic::new(
                diag::EXPR_VALUE,
                format!(
                    "flow.body_force[{i}] is not finite at the centroid {:?} of cell {c} (t = {t} s)",
                    geo.cell_centroid[c]
                ),
            )
            .key(format!("flow.body_force[{i}]")),
        ));
    }
    Ok(Some(out))
}

/// True when the body force or a boundary value depends on t.
pub fn time_dependent(case: &FlowCase, values: &bc::Compiled) -> Result<(bool, bool)> {
    let source = case.body_force_values()?.iter().any(|v| v.uses().time());
    let boundary = values.values().any(|v| v.uses().time());
    Ok((source, boundary))
}

/// As [`run`], with the requests, limits and run record of `steer`
/// (crate::engine).
pub fn run_steered(
    case: &FlowCase,
    name: &str,
    out: &Path,
    steer: &mut Steer,
    log: impl FnMut(&str),
) -> Result<FlowReport> {
    if case.time.is_some() {
        return crate::unsteady::run(case, name, out, steer, log);
    }
    let t0 = Instant::now();
    std::fs::create_dir_all(out)?;
    let control = Control::new(out)?;
    let mut case = case.clone();
    let mut log = log;
    let mut lg = Log {
        out: &mut log,
        text: String::new(),
    };
    let levels: Vec<Option<u32>> = if case.verification.levels.is_empty() {
        vec![None]
    } else {
        case.verification.levels.iter().map(|&l| Some(l)).collect()
    };
    let analytic = case.analytic();
    let phys = physics(&case)?;
    let values = case.compiled_boundaries()?;
    let at = bc::At {
        analytic: analytic.as_ref(),
        values: &values,
        t: 0.0,
    };
    let reynolds = case.reynolds();
    let mut rep = FlowReport {
        name: name.into(),
        reynolds,
        findings: Vec::new(),
        solves: Vec::new(),
        orders: Vec::new(),
        time_orders: Vec::new(),
        stopped: false,
        stop: None,
        warnings: Vec::new(),
        failures: Vec::new(),
        seconds: 0.0,
    };
    std::fs::write(out.join("case.resolved.toml"), case.resolved_toml()?)?;
    let last = levels.len() - 1;
    let single = levels.len() == 1;
    if !single && case.run.restart.is_some() {
        bail!("[run] restart continues a single run, not a sweep over levels");
    }
    let mut histories: Vec<(String, Vec<(f64, f64)>)> = Vec::new();
    let mut stopped = false;
    for (li, level) in levels.iter().enumerate() {
        let (mesh, geo, desc) = mesh_for(&case, *level)?;
        let mesh_hash = crate::engine::mesh_hash(&mesh);
        // boundary-condition check (D30)
        let findings = bc::check(
            &mesh,
            &geo,
            &case.boundary,
            phys.rho,
            at,
            case.check_context(),
        );
        if li == 0 {
            for f in &findings {
                lg.line(&f.to_string());
            }
            rep.findings = findings.clone();
        }
        if findings.iter().any(|f| f.level == Level::Error) {
            let msg: Vec<String> = findings
                .iter()
                .filter(|f| f.level == Level::Error)
                .map(|f| f.to_string())
                .collect();
            std::fs::write(out.join("flow.log"), &lg.text)?;
            bail!(
                "the boundary conditions cannot be run:\n  {}",
                msg.join("\n  ")
            );
        }
        let nc = mesh.n_cells;
        if li == 0 || li == last {
            header(&mut lg, &case, &mesh, &desc, li == last);
        }
        lg.line(&format!("mesh  {desc}"));
        let arrays = bc::arrays(&mesh, &geo, &case.boundary, phys.rho, at)?;
        let disc = simplec::discretisation(&mesh, &geo, &arrays)?;
        let phys = physics_on(&case, &phys, &mesh, &geo)?;
        let loc = Locator::new(&mesh, &geo);
        let reference_cell = if arrays.has_pressure() {
            None
        } else {
            let x = case.pressure.reference_point.unwrap_or_else(|| {
                let (lo, hi) = mesh.bounding_box();
                std::array::from_fn(|i| lo[i] + 1e-6 * (hi[i] - lo[i]))
            });
            Some(locate(&loc, &geo, x))
        };
        let p0 = if reference_cell.is_some() {
            case.pressure.reference_value
        } else {
            case.initial.pressure
        };
        let mut fl = Flow::new(
            disc,
            arrays,
            phys,
            case.solver.clone(),
            case.initial.velocity,
            p0,
            reference_cell,
            &crate::simplec::Target {
                device: case.run.backend.clone(),
                threads: case.run.threads,
                ptx: case.run.ptx.clone(),
            },
        )?;
        if let Some(src) = momentum_source(&case, &geo, nc, 0.0)? {
            fl.set_momentum_source(src)?;
        }
        if let Some(path) = case.run.restart.clone() {
            let r = restart::read(&path)?;
            let (r, same) =
                r.to_local(&mesh.cell_id[..nc], &mesh.face_id, mesh.n_interior_faces())?;
            if r.transient.is_some() {
                bail!(
                    "{} is the restart file of a transient run; this case is steady (no [time] table)",
                    path.display()
                );
            }
            fl.load(&r.state)?;
            fl.iteration = r.iteration;
            if let Some(m) = &r.amg_matrix {
                fl.restore_amg(m, r.amg_setups)?;
            }
            if !same {
                let m = "note: the restart file has another cell or face order; the fields were \
                         mapped by identity and the multigrid rebuilt (the continuation is not \
                         bit-identical)";
                lg.line(m);
                rep.warnings
                    .push(diag::Diagnostic::new(diag::RUN_RESTART_MAPPED, m));
            }
            lg.line(&format!(
                "restart from {} at iteration {}",
                path.display(),
                r.iteration
            ));
        } else if case.initial.potential_flow {
            let it = fl.potential_flow(&mesh, &geo)?;
            lg.line(&format!("start  potential flow (D52): {it} CG iterations"));
        }
        let first = fl.iteration();
        let mut csv = csv_preamble(
            "residuals",
            name,
            None,
            "iteration -, momentum_x -, momentum_y -, momentum_z -, continuity - (scaled residuals), \
             bicgstab_x -, bicgstab_y -, bicgstab_z -, cg - (iterations)",
        );
        csv.push_str(
            "iteration,momentum_x,momentum_y,momentum_z,continuity,bicgstab_x,bicgstab_y,bicgstab_z,cg\n",
        );
        let (bunits, bhead) = balance::csv_head(false);
        let mut bcsv = csv_preamble("balance", name, None, &bunits);
        bcsv.push_str(&bhead);
        bcsv.push('\n');
        let nonorth = case.solver.nonorth_correction;
        let mut hist = Vec::new();
        let mut converged = false;
        let mut res = Residuals::default();
        let ts = Instant::now();
        loop {
            if fl.iteration() >= case.solver.max_iterations {
                break;
            }
            let previous = res;
            fl.keep_valid()?;
            res = match fl.iterate() {
                Ok(r) => r,
                Err(e) if divergence::is_divergence(&e, &mut fl) => {
                    // D63: the last valid state, the snapshot and the report
                    lg.line(&format!("      {e:#}"));
                    let at = fl.iteration();
                    let ctx = divergence::Context {
                        name,
                        case: &case,
                        at,
                        done: at.saturating_sub(first + 1),
                        time: None,
                        dt: None,
                        residuals: fl.diverged.unwrap_or_default(),
                        last_residuals: (at > first + 1).then_some(previous),
                        first_courant: None,
                    };
                    let prepared = (|| -> Result<(simplec::State, Vec<&'static str>)> {
                        let failed = fl.state()?;
                        anyhow::ensure!(fl.rewind()?, "no valid state was kept");
                        let mut files = Vec::new();
                        if cfg!(feature = "hdf5") && case.output.vtkhdf {
                            let sol = Solution {
                                state: fl.state()?,
                                disc: &fl.disc,
                                phys: fl.phys.clone(),
                            };
                            write_fields(
                                &case,
                                &out.join("fields_last_valid.vtkhdf"),
                                &mesh,
                                &geo,
                                &sol,
                                analytic.as_ref(),
                                Some((fl.iteration(), None)),
                            )?;
                            files.push("fields_last_valid.vtkhdf");
                        }
                        if single && restart::AVAILABLE {
                            let rs = steady_restart(&mut fl, &mesh, name)?;
                            restart::write(&out.join("restart.h5"), &rs)?;
                            files.push("restart.h5");
                        }
                        Ok((failed, files))
                    })();
                    let e =
                        divergence::finish(&mut fl, &ctx, &mesh, &geo, out, prepared, &mut |l| {
                            lg.line(l)
                        });
                    steer.record.flush();
                    std::fs::write(out.join("flow.log"), &lg.text)?;
                    return Err(e);
                }
                Err(e) => return Err(e),
            };
            let it = fl.iteration();
            let done = res.max() < case.solver.tolerance;
            if it % case.output.history_every == 0 || done {
                let _ = writeln!(
                    csv,
                    "{},{:e},{:e},{:e},{:e},{},{},{},{}",
                    it,
                    res.momentum[0],
                    res.momentum[1],
                    res.momentum[2],
                    res.continuity,
                    res.momentum_iterations[0],
                    res.momentum_iterations[1],
                    res.momentum_iterations[2],
                    res.pressure_iterations
                );
                hist.push((it as f64, res.max()));
                if steer.record.on() {
                    steer.record.row(&serde_json::json!({
                        "type": "iteration",
                        "iteration": it,
                        "momentum": res.momentum,
                        "continuity": res.continuity,
                        "bicgstab": res.momentum_iterations,
                        "cg": res.pressure_iterations,
                    }));
                }
            }
            if it == first + 1 || it % case.solver.log_every == 0 || done {
                steer.record.flush();
                lg.line(&format!(
                    "iter {it:6}  momentum {:.3e} {:.3e} {:.3e}  continuity {:.3e}  bicgstab {:?}  cg {}",
                    res.momentum[0],
                    res.momentum[1],
                    res.momentum[2],
                    res.continuity,
                    res.momentum_iterations,
                    res.pressure_iterations
                ));
            }
            if done {
                converged = true;
                break;
            }
            // the ledger of the current state (D88, D89); the next
            // iteration recomputes the boundary values and gradients
            if case.balance.every > 0 && it.is_multiple_of(case.balance.every) {
                fl.update_boundaries()?;
                fl.gradients()?;
                let sol = Solution {
                    state: fl.state()?,
                    disc: &fl.disc,
                    phys: fl.phys.clone(),
                };
                let r = balance::rates(
                    &sol,
                    &mesh,
                    &geo,
                    &case,
                    &fl.bc.mom_kind,
                    &geo.cell_volume,
                    nonorth,
                );
                let b = balance::steady(&r, it);
                bcsv.push_str(&balance::csv_row(&b));
                bcsv.push('\n');
                if steer.record.on() {
                    steer.record.row(&balance::record_row(&b));
                }
                lg.line(
                    balance::table(&format!("      balance of iteration {it}"), &b, false)
                        .trim_end(),
                );
            }
            // control files (D40), signals and limits (D64, D83)
            let req = control.check(it, case.run.reread_every);
            if req.result && li == last {
                fl.update_boundaries()?;
                fl.gradients()?;
                let sol = Solution {
                    state: fl.state()?,
                    disc: &fl.disc,
                    phys: fl.phys.clone(),
                };
                write_fields(
                    &case,
                    &out.join(format!("fields_{it:06}.vtkhdf")),
                    &mesh,
                    &geo,
                    &sol,
                    analytic.as_ref(),
                    Some((it, None)),
                )?;
                lg.line(&format!(
                    "      fields of iteration {it} written (fields_{it:06}.vtkhdf)"
                ));
            }
            let r = &case.run;
            let restart_now = single
                && ((r.restart_every > 0 && it.is_multiple_of(r.restart_every))
                    || r.restart_at.contains(&it)
                    || steer.restart_requested());
            if restart_now && restart::AVAILABLE {
                let rs = steady_restart(&mut fl, &mesh, name)?;
                restart::write(&out.join("restart.h5"), &rs)?;
                if case.run.restart_at.contains(&it) {
                    restart::write(&out.join(format!("restart_{it:06}.h5")), &rs)?;
                }
                lg.line(&format!("      iteration {it}: restart.h5 written"));
            }
            if req.reread {
                reread_case(
                    &mut case,
                    &mut fl,
                    &mut lg,
                    steer,
                    &mut rep.warnings,
                    it,
                    out,
                )?;
            }
            if let Some(s) = steer.poll(&case.run.limits, it - first, false, req.stop) {
                lg.line(&format!("      iteration {it}: {}", s.message));
                stopped = true;
                break;
            }
        }
        steer.record.flush();
        let secs = ts.elapsed().as_secs_f64();
        let iters = fl.iteration();
        lg.line(&format!(
            "      {} in {iters} outer iterations, {secs:.2} s, {:.0} ns/cell/iteration",
            if converged {
                "converged"
            } else if stopped {
                "stopped"
            } else {
                "NOT converged"
            },
            1e9 * fl.timings.total() / (nc as f64 * (iters - first).max(1) as f64)
        ));
        let suffix = level.map_or(String::new(), |l| format!("_level{l}"));
        std::fs::write(out.join(format!("residuals{suffix}.csv")), &csv)?;
        histories.push((
            level.map_or("residual".into(), |l| format!("level {l}")),
            hist,
        ));
        if !converged && !stopped && case.verification.require_convergence.unwrap_or(true) {
            rep.failures.push(format!(
                "{desc}: not converged in {iters} outer iterations (largest residual {:.2e})",
                res.max()
            ));
        }
        // final gradients for the output
        fl.update_boundaries()?;
        fl.gradients()?;
        let fl_host = Solution {
            state: fl.state()?,
            disc: &fl.disc,
            phys: fl.phys.clone(),
        };
        let sol = &fl_host;
        let errors = analytic.map(|a| error_norms(sol, &mesh, &geo, &a, dimension(&case)));
        if let Some(e) = &errors {
            lg.line(&format!(
                "      velocity error L2 {:.3e}, Linf {:.3e} (relative), pressure Linf {:.3e}",
                e.l2, e.linf, e.pressure_linf
            ));
        }
        let reference = case
            .verification
            .reference
            .map(|r| compare_reference(sol, &mesh, &geo, &loc, r, &case));
        if let Some(r) = &reference {
            lg.line(&format!(
                "      {}: largest deviation u {:.4}, v {:.4} (lid speeds)",
                r.reference, r.u_max_dev, r.v_max_dev
            ));
        }
        let patches = patch_summaries(sol, &mesh, &geo, &case);
        let total: f64 = patches.iter().map(|p| p.mass_flow).sum();
        // scale: the largest patch mass flow, or in a closed domain the
        // flow ρ |u|max A through the largest patch
        let umax = (0..nc).map(|c| norm(sol.u.get(c))).fold(0.0f64, f64::max);
        let big = patches
            .iter()
            .map(|p| p.mass_flow.abs().max(case.rho() * umax * p.area))
            .fold(0.0f64, f64::max);
        let mass_imbalance = if big > 0.0 { total / big } else { 0.0 };
        let coefficients = crate::monitor::coefficients(&case, &patches);
        if let Some([cd, cl]) = coefficients {
            lg.line(&format!(
                "      force coefficients C_D {cd:.5}, C_L {cl:.5}"
            ));
        }
        let benchmark =
            crate::unsteady::steady_benchmark_of(&case, sol, &mesh, &geo, &loc, &patches);
        if let Some(b) = &benchmark {
            crate::unsteady::log_benchmark(b, |l| lg.line(l));
        }
        // the ledger of the final state (D88, D89, D90 part 4)
        let bal = balance::steady(
            &balance::rates(
                sol,
                &mesh,
                &geo,
                &case,
                &fl.bc.mom_kind,
                &geo.cell_volume,
                nonorth,
            ),
            iters,
        );
        lg.line(balance::table("      balance of the final state", &bal, false).trim_end());
        bcsv.push_str(&balance::csv_row(&bal));
        bcsv.push('\n');
        if steer.record.on() {
            steer.record.row(&balance::record_row(&bal));
        }
        std::fs::write(out.join(format!("balance{suffix}.csv")), &bcsv)?;
        let summary = BalanceSummary {
            largest_mass_relative: bal.mass.relative.abs(),
            last: bal,
            ..Default::default()
        };
        for d in balance::check(&summary, &case.balance) {
            lg.line(&format!("      warning {}: {}", d.code, d.message));
            steer.event(d.code, iters, &d.message);
            rep.warnings.push(d);
        }
        if li == last || stopped {
            std::fs::write(
                out.join("balance.svg"),
                balance::chart(name, &summary.last, false),
            )?;
            write_outputs(
                &case,
                out,
                &mesh,
                &geo,
                sol,
                analytic.as_ref(),
                &patches,
                &loc,
                Some((iters, None)),
            )?;
            if let Some(r) = &reference {
                write_reference_svg(out, name, r)?;
            }
        }
        drop(fl_host);
        if single && (case.run.restart_final || stopped) && restart::AVAILABLE {
            let rs = steady_restart(&mut fl, &mesh, name)?;
            restart::write(&out.join("restart.h5"), &rs)?;
        }
        rep.solves.push(Solve {
            mesh: desc,
            mesh_hash,
            device: fl.device(),
            cells: nc,
            iterations: iters,
            converged,
            residuals: res,
            ns_per_cell_iteration: 1e9 * fl.timings.total()
                / (nc as f64 * (iters - first).max(1) as f64),
            timings: fl.timings.clone(),
            amg_setups: fl.amg_setups,
            errors,
            reference,
            patches,
            mass_imbalance,
            coefficients,
            transient: None,
            benchmark,
            balance: Some(summary),
        });
        if stopped {
            rep.stopped = true;
            rep.stop = steer.stop.clone();
            break;
        }
    }
    // checks
    let v = &case.verification;
    for s in &rep.solves {
        if let (Some(tol), Some(b)) = (v.benchmark_tolerance, &s.benchmark)
            && b.max_deviation() > tol
        {
            rep.failures.push(format!(
                "{}: largest distance from the {} ranges {:.4} above {tol}",
                s.mesh,
                b.benchmark,
                b.max_deviation()
            ));
        }
    }
    if levels.len() >= 2 && rep.solves.iter().all(|s| s.errors.is_some()) {
        let e: Vec<&ErrorNorms> = rep
            .solves
            .iter()
            .filter_map(|s| s.errors.as_ref())
            .collect();
        rep.orders = e
            .windows(2)
            .map(|w| (w[0].l2 / w[1].l2).ln() / (w[0].h / w[1].h).ln())
            .collect();
        lg.line(&format!("observed orders {:.3?}", rep.orders));
        if let (Some([lo, hi]), Some(&o)) = (v.expected_order, rep.orders.last())
            && !(lo..=hi).contains(&o)
        {
            rep.failures.push(format!(
                "observed order {o:.3} on the finest pair, outside [{lo}, {hi}]"
            ));
        }
        let pts: Vec<(f64, f64)> = e.iter().map(|e| (e.h, e.l2)).collect();
        let slope: Vec<(f64, f64)> = match (pts.first(), pts.last()) {
            (Some(&a), Some(&b)) => vec![(a.0, b.1 * (a.0 / b.0).powi(2)), b],
            _ => Vec::new(),
        };
        std::fs::write(
            out.join("convergence.svg"),
            plot::chart(
                &format!("{name}: velocity error over mesh size"),
                "h (m)",
                "L2 error / U",
                &[
                    Series {
                        label: "Foehncraft",
                        points: &pts,
                        markers: true,
                    },
                    Series {
                        label: "slope 2",
                        points: &slope,
                        markers: false,
                    },
                ],
                true,
                true,
            ),
        )?;
    }
    for s in &rep.solves {
        if let (Some(max), Some(e)) = (v.max_error, &s.errors)
            && e.l2 > max
        {
            rep.failures.push(format!(
                "{}: relative L2 velocity error {:.3e} above {max:e}",
                s.mesh, e.l2
            ));
        }
        if let (Some(tol), Some(r)) = (v.reference_tolerance, &s.reference)
            && r.u_max_dev.max(r.v_max_dev) > tol
        {
            rep.failures.push(format!(
                "{}: deviation from {} of {:.4} above {tol}",
                s.mesh,
                r.reference,
                r.u_max_dev.max(r.v_max_dev)
            ));
        }
        if s.mass_imbalance.abs() > 1e-6 {
            rep.failures.push(format!(
                "{}: net mass flow out of the domain {:.2e} of the largest boundary flow",
                s.mesh, s.mass_imbalance
            ));
        }
    }
    let series: Vec<Series> = histories
        .iter()
        .map(|(l, h)| Series {
            label: l,
            points: h,
            markers: false,
        })
        .collect();
    std::fs::write(
        out.join("residuals.svg"),
        plot::chart(
            &format!("{name}: largest residual"),
            "outer iteration",
            "residual",
            &series,
            false,
            true,
        ),
    )?;
    rep.seconds = t0.elapsed().as_secs_f64();
    for f in &rep.failures {
        lg.line(&format!("FAILED: {f}"));
    }
    std::fs::write(out.join("flow.log"), &lg.text)?;
    std::fs::write(out.join("flow.json"), serde_json::to_string_pretty(&rep)?)?;
    std::fs::write(out.join("flow.md"), rep.to_markdown())?;
    Ok(rep)
}

/// Re-reads the case file at iteration or step `n` (D40) and applies what
/// may change during a run; returns the case as it was when something was
/// applied.
pub(crate) fn reread_case<F: FnMut(&str)>(
    case: &mut FlowCase,
    fl: &mut Flow,
    lg: &mut Log<F>,
    steer: &mut Steer,
    warnings: &mut Vec<diag::Diagnostic>,
    n: usize,
    out: &Path,
) -> Result<Option<FlowCase>> {
    let what = if case.time.is_some() {
        "step"
    } else {
        "iteration"
    };
    let Some(src) = case.source.clone() else {
        lg.line(&format!("      {what} {n}: no case file to re-read"));
        return Ok(None);
    };
    match control::reread(&src, case) {
        Ok(rr) => {
            if rr.applied.is_empty() && rr.not_applied.is_empty() {
                lg.line(&format!(
                    "      {what} {n}: case file re-read, nothing changed"
                ));
            }
            for l in &rr.applied {
                lg.line(&format!("      {what} {n}: changed {l}"));
            }
            for l in &rr.not_applied {
                let m = format!("{what} {n}: not applied (cannot change during a run): {l}");
                lg.line(&format!("      {m}"));
                steer.event(diag::RUN_NOT_APPLIED, n, &m);
                warnings.push(diag::Diagnostic::new(diag::RUN_NOT_APPLIED, m));
            }
            let old = if rr.applied.is_empty() {
                None
            } else {
                fl.set_settings(rr.case.solver.clone())?;
                Some(std::mem::replace(case, rr.case))
            };
            std::fs::write(
                out.join(format!("case.read_toml_{n}.toml")),
                toml::to_string_pretty(&*case)?,
            )?;
            Ok(old)
        }
        Err(e) => {
            let m =
                format!("{what} {n}: the case file was not applied, it does not validate: {e:#}");
            lg.line(&format!("      {m}"));
            steer.event(diag::RUN_REREAD_INVALID, n, &m);
            warnings.push(diag::Diagnostic::new(diag::RUN_REREAD_INVALID, m));
            Ok(None)
        }
    }
}

pub(crate) fn header<F: FnMut(&str)>(
    lg: &mut Log<F>,
    case: &FlowCase,
    mesh: &Mesh,
    desc: &str,
    last: bool,
) {
    let nc = mesh.n_cells;
    let walls = case
        .boundary
        .iter()
        .filter(|(_, b)| matches!(b, BoundaryCondition::Wall { .. }))
        .count();
    let rows = case.solver.max_iterations / case.output.history_every;
    lg.line(&format!(
        "Foehncraft {} flow solver (SIMPLEC, steady, laminar)",
        env!("CARGO_PKG_VERSION")
    ));
    lg.line(&format!(
        "fluid  ρ = {} kg/m³, μ = {} Pa s; convection {}",
        case.rho(),
        case.mu(),
        if case.flow.convection {
            "on"
        } else {
            "off (Stokes)"
        }
    ));
    if last {
        lg.line(&format!("output for {desc}:"));
        lg.line(&format!(
            "  residuals.csv: up to {rows} rows (every {} iterations), about {} kB",
            case.output.history_every,
            rows * 90 / 1000 + 1
        ));
        if case.output.vtkhdf {
            lg.line(&format!(
                "  fields.vtkhdf: velocity, pressure{} on {nc} cells, about {} MB",
                if case.analytic.is_some() {
                    ", error"
                } else {
                    ""
                },
                (nc * (8 * 6) + mesh.face_nodes.len() * 8 + mesh.points.len() * 24) / 1_000_000 + 1
            ));
        }
        if case.output.wall_csv && walls > 0 {
            lg.line(&format!(
                "  wall_<patch>.csv for {walls} wall(s): pressure and wall shear stress per face"
            ));
        }
        for l in &case.output.lines {
            lg.line(&format!("  line_{}.csv: {} points", l.name, l.points));
        }
        lg.line("  balance.csv, balance.svg, patches.csv, flow.json, flow.md, residuals.svg, case.resolved.toml");
    }
}

pub(crate) fn error_norms(
    fl: &Solution,
    mesh: &Mesh,
    geo: &Geometry,
    a: &Analytic,
    dim: u8,
) -> ErrorNorms {
    let (mut e2, mut einf, mut vol, mut pinf) = (0.0, 0.0f64, 0.0, 0.0f64);
    let (mut pmin, mut pmax) = (f64::INFINITY, f64::NEG_INFINITY);
    // transient solutions: the pressure level is free (reference cell), so
    // the volume means are compared
    let (mut pm, mut pem) = (0.0, 0.0);
    if a.is_transient() {
        let mut v = 0.0;
        for c in 0..mesh.n_cells {
            let w = geo.cell_volume[c];
            pm += fl.p[c] * w;
            pem += a.pressure(geo.cell_centroid[c]) * w;
            v += w;
        }
        pm /= v;
        pem /= v;
    }
    for c in 0..mesh.n_cells {
        let x = geo.cell_centroid[c];
        let ue = a.velocity(x);
        let e = norm(sub(fl.u.get(c), ue));
        e2 += e * e * geo.cell_volume[c];
        einf = einf.max(e);
        vol += geo.cell_volume[c];
        let pe = a.pressure(x);
        pmin = pmin.min(pe);
        pmax = pmax.max(pe);
        pinf = pinf.max(((fl.p[c] - pm) - (pe - pem)).abs());
    }
    let n = mesh.n_cells as f64;
    let h = if dim == 2 {
        let (lo, hi) = mesh.bounding_box();
        (vol / (hi[2] - lo[2]) / n).sqrt()
    } else {
        (vol / n).cbrt()
    };
    let us = a.velocity_scale().max(f64::MIN_POSITIVE);
    let prange = pmax - pmin;
    ErrorNorms {
        cells: mesh.n_cells,
        h,
        l2: (e2 / vol).sqrt() / us,
        linf: einf / us,
        pressure_linf: if prange > 0.0 { pinf / prange } else { pinf },
    }
}

pub(crate) fn compare_reference(
    fl: &Solution,
    mesh: &Mesh,
    geo: &Geometry,
    loc: &Locator,
    r: Reference,
    case: &FlowCase,
) -> ReferenceComparison {
    let (ys, us, xs, vs) = r.centre_lines();
    let lid = case.verification.lid_velocity.unwrap_or(1.0);
    let (lo, hi) = mesh.bounding_box();
    let zc = 0.5 * (lo[2] + hi[2]);
    let span = [hi[0] - lo[0], hi[1] - lo[1]];
    let mut u = Vec::new();
    let mut v = Vec::new();
    let (mut du, mut dv) = (0.0f64, 0.0f64);
    for i in 0..ys.len() {
        // keep the probe points inside the domain
        let y = ys[i].clamp(1e-9, 1.0 - 1e-9);
        let x = [lo[0] + 0.5 * span[0], lo[1] + y * span[1], zc];
        let c = locate(loc, geo, x);
        let val = value_at(fl, geo, c, x).0[0] / lid;
        let val = if ys[i] == 1.0 {
            1.0
        } else if ys[i] == 0.0 {
            0.0
        } else {
            val
        };
        du = du.max((val - us[i]).abs());
        u.push((ys[i], val));
        let xx = xs[i].clamp(1e-9, 1.0 - 1e-9);
        let x = [lo[0] + xx * span[0], lo[1] + 0.5 * span[1], zc];
        let c = locate(loc, geo, x);
        let val = if xs[i] == 1.0 || xs[i] == 0.0 {
            0.0
        } else {
            value_at(fl, geo, c, x).0[1] / lid
        };
        dv = dv.max((val - vs[i]).abs());
        v.push((xs[i], val));
    }
    let n = 200;
    let (mut u_line, mut v_line) = (Vec::new(), Vec::new());
    for i in 0..=n {
        let t = i as f64 / n as f64;
        let tt = t.clamp(1e-9, 1.0 - 1e-9);
        let x = [lo[0] + 0.5 * span[0], lo[1] + tt * span[1], zc];
        let uu = if i == n {
            1.0
        } else if i == 0 {
            0.0
        } else {
            value_at(fl, geo, locate(loc, geo, x), x).0[0] / lid
        };
        u_line.push((t, uu));
        let x = [lo[0] + tt * span[0], lo[1] + 0.5 * span[1], zc];
        let vv = if i == 0 || i == n {
            0.0
        } else {
            value_at(fl, geo, locate(loc, geo, x), x).0[1] / lid
        };
        v_line.push((t, vv));
    }
    ReferenceComparison {
        reference: r.label().into(),
        u_max_dev: du,
        v_max_dev: dv,
        u,
        v,
        u_line,
        v_line,
    }
}

pub(crate) fn wall_face_values(
    fl: &Solution,
    geo: &Geometry,
    f: usize,
    ni: usize,
) -> (f64, [f64; 3], [f64; 3]) {
    // pressure at the face (zero normal gradient: tangential extrapolation),
    // viscous force of the fluid on the face, wall shear stress
    let fb = f - ni;
    let c = fl.disc.owner[f] as usize;
    let s = geo.face_area[f];
    let r = sub(geo.face_centroid[f], geo.cell_centroid[c]);
    let ss = dot(s, s);
    let rn = dot(r, s) / ss;
    let rt = [r[0] - rn * s[0], r[1] - rn * s[1], r[2] - rn * s[2]];
    let p = fl.p[c] + dot(fl.grad_p.get(c), rt);
    let ub = fl.ub.get(fb);
    let up = fl.u.get(c);
    let k = fl.disc.k.get(f);
    let dg = fl.disc.dgeo[f];
    let mu = fl.phys.mu_at(f);
    let mut fv = [0.0; 3];
    for i in 0..3 {
        fv[i] = -mu * (dg * (ub[i] - up[i]) + dot(k, fl.grad_u[i].get(c)));
    }
    let a = ss.sqrt();
    let n = [s[0] / a, s[1] / a, s[2] / a];
    // shear: tangential part of the viscous force per area, on the wall
    let fn_ = dot(fv, n);
    let tau = [
        (fv[0] - fn_ * n[0]) / a,
        (fv[1] - fn_ * n[1]) / a,
        (fv[2] - fn_ * n[2]) / a,
    ];
    (p, fv, tau)
}

pub(crate) fn patch_summaries(
    fl: &Solution,
    mesh: &Mesh,
    geo: &Geometry,
    case: &FlowCase,
) -> Vec<PatchSummary> {
    let ni = mesh.n_interior_faces();
    mesh.patches
        .iter()
        .map(|p| {
            let bc = &case.boundary[&p.name];
            let (mut area, mut q, mut pa) = (0.0, 0.0, 0.0);
            let mut force = [0.0; 3];
            for f in p.range() {
                let s = geo.face_area[f];
                let a = norm(s);
                area += a;
                q += fl.mass_flux[f];
                let (pf, fv) = if matches!(bc, BoundaryCondition::Wall { .. }) {
                    let (pf, fv, _) = wall_face_values(fl, geo, f, ni);
                    (pf, fv)
                } else {
                    (fl.pb[f - ni], [0.0; 3])
                };
                pa += pf * a;
                for i in 0..3 {
                    force[i] += pf * s[i] + fv[i];
                }
            }
            PatchSummary {
                name: p.name.clone(),
                kind: bc.type_name().into(),
                area,
                mass_flow: q,
                mean_pressure: if area > 0.0 { pa / area } else { 0.0 },
                force,
            }
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn write_outputs(
    case: &FlowCase,
    out: &Path,
    mesh: &Mesh,
    geo: &Geometry,
    fl: &Solution,
    analytic: Option<&Analytic>,
    patches: &[PatchSummary],
    loc: &Locator,
    at: Option<(usize, Option<f64>)>,
) -> Result<()> {
    let name = case_name(case);
    let mut csv = csv_preamble(
        "patches",
        &name,
        at,
        "patch -, type -, area m², mass_flow kg/s (out of the domain), mean_pressure Pa, \
         force_x force_y force_z N",
    );
    csv.push_str("patch,type,area,mass_flow,mean_pressure,force_x,force_y,force_z\n");
    for p in patches {
        let _ = writeln!(
            csv,
            "{},{},{:e},{:e},{:e},{:e},{:e},{:e}",
            p.name,
            p.kind,
            p.area,
            p.mass_flow,
            p.mean_pressure,
            p.force[0],
            p.force[1],
            p.force[2]
        );
    }
    std::fs::write(out.join("patches.csv"), csv)?;
    let ni = mesh.n_interior_faces();
    if case.output.wall_csv {
        for p in &mesh.patches {
            if !matches!(case.boundary[&p.name], BoundaryCondition::Wall { .. }) {
                continue;
            }
            let mut csv = csv_preamble(
                "wall",
                &name,
                at,
                "x y z m (face centroid), pressure Pa, tau_x tau_y tau_z Pa (wall shear stress)",
            );
            csv.push_str("x,y,z,pressure,tau_x,tau_y,tau_z\n");
            let mut rows: Vec<(u64, String)> = p
                .range()
                .map(|f| {
                    let x = geo.face_centroid[f];
                    let (pf, _, tau) = wall_face_values(fl, geo, f, ni);
                    (
                        mesh.face_id[f],
                        format!(
                            "{:e},{:e},{:e},{:e},{:e},{:e},{:e}",
                            x[0], x[1], x[2], pf, tau[0], tau[1], tau[2]
                        ),
                    )
                })
                .collect();
            // global face order: independent of the renumbering (R2)
            rows.sort_by_key(|r| r.0);
            for (_, r) in rows {
                csv.push_str(&r);
                csv.push('\n');
            }
            std::fs::write(out.join(format!("wall_{}.csv", p.name)), csv)?;
        }
    }
    for l in &case.output.lines {
        let mut csv = csv_preamble(
            "line",
            &name,
            at,
            "s m (distance from `from`), x y z m, u v w m/s, p Pa",
        );
        csv.push_str("s,x,y,z,u,v,w,p\n");
        let d = sub(l.to, l.from);
        let len = norm(d);
        for i in 0..l.points {
            let t = i as f64 / (l.points - 1) as f64;
            let x = std::array::from_fn(|k| l.from[k] + t * d[k]);
            let c = locate(loc, geo, x);
            let (u, p) = value_at(fl, geo, c, x);
            let _ = writeln!(
                csv,
                "{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e}",
                t * len,
                x[0],
                x[1],
                x[2],
                u[0],
                u[1],
                u[2],
                p
            );
        }
        std::fs::write(out.join(format!("line_{}.csv", l.name)), csv)?;
    }
    write_fields(
        case,
        &out.join("fields.vtkhdf"),
        mesh,
        geo,
        fl,
        analytic,
        at,
    )
}

/// The case name (the stem of the case file).
pub(crate) fn case_name(case: &FlowCase) -> String {
    case.source
        .as_ref()
        .and_then(|p| p.file_stem())
        .map_or("flow".into(), |s| s.to_string_lossy().into_owned())
}

/// The comment lines at the top of a CSV file (D75): format and version,
/// the Foehncraft version, the case, the iteration or step and time of
/// the state it holds, and the units of the columns (`-`: dimensionless).
pub(crate) fn csv_preamble(
    kind: &str,
    case: &str,
    at: Option<(usize, Option<f64>)>,
    units: &str,
) -> String {
    let when = match at {
        Some((n, Some(t))) => format!("; step {n}; time {t:e} s"),
        Some((n, None)) => format!("; iteration {n}"),
        None => String::new(),
    };
    format!(
        "# format foehncraft-{kind} 1; Foehncraft {}; case {case}{when}\n# units: {units}\n",
        foehncraft_core::VERSION
    )
}

/// Units of the fields of `fields*.vtkhdf` (attribute `units`, D75).
pub const FIELD_UNITS: &str = "coordinates m; velocity_x m/s; velocity_y m/s; velocity_z m/s; \
     speed m/s; pressure Pa; velocity_error m/s";

/// Writes the fields; `at` is the iteration or step and the time (D75).
#[cfg(feature = "hdf5")]
pub(crate) fn write_fields(
    case: &FlowCase,
    file: &Path,
    mesh: &Mesh,
    geo: &Geometry,
    fl: &Solution,
    analytic: Option<&Analytic>,
    at: Option<(usize, Option<f64>)>,
) -> Result<()> {
    use foehncraft_mesh::vtkhdf::{Attr, CellField, FaceField, write_all};
    if !case.output.vtkhdf {
        return Ok(());
    }
    let derived = crate::derived::compute(case, &case.output.derived, mesh, geo, fl)?;
    let speed: Vec<f64> = (0..mesh.n_cells).map(|c| norm(fl.u.get(c))).collect();
    let mut fields = vec![
        CellField {
            name: "velocity_x",
            values: &fl.u.x,
        },
        CellField {
            name: "velocity_y",
            values: &fl.u.y,
        },
        CellField {
            name: "velocity_z",
            values: &fl.u.z,
        },
        CellField {
            name: "speed",
            values: &speed,
        },
        CellField {
            name: "pressure",
            values: &fl.p,
        },
    ];
    let err: Vec<f64>;
    if let Some(a) = analytic {
        err = (0..mesh.n_cells)
            .map(|c| norm(sub(fl.u.get(c), a.velocity(geo.cell_centroid[c]))))
            .collect();
        fields.push(CellField {
            name: "velocity_error",
            values: &err,
        });
    }
    for (name, values) in &derived.cells {
        fields.push(CellField { name, values });
    }
    let faces: Vec<FaceField> = derived
        .faces
        .iter()
        .map(|(name, values)| FaceField { name, values })
        .collect();
    let name = case_name(case);
    let units = format!(
        "{FIELD_UNITS}{}",
        crate::derived::Fields::units(&case.output.derived)
    );
    let mut attrs = vec![
        ("foehncraft_format", Attr::Str("foehncraft-fields")),
        ("foehncraft_format_version", Attr::U64(1)),
        ("foehncraft_version", Attr::Str(foehncraft_core::VERSION)),
        ("case", Attr::Str(&name)),
        ("units", Attr::Str(&units)),
    ];
    match at {
        Some((n, Some(t))) => {
            attrs.push(("step", Attr::U64(n as u64)));
            attrs.push(("time", Attr::F64(t)));
        }
        Some((n, None)) => attrs.push(("iteration", Attr::U64(n as u64))),
        None => {}
    }
    write_all(file, mesh, &fields, &faces, &attrs)
}

#[cfg(not(feature = "hdf5"))]
pub(crate) fn write_fields(
    _: &FlowCase,
    _: &Path,
    _: &Mesh,
    _: &Geometry,
    _: &Solution,
    _: Option<&Analytic>,
    _: Option<(usize, Option<f64>)>,
) -> Result<()> {
    Ok(())
}

pub(crate) fn write_reference_svg(out: &Path, name: &str, r: &ReferenceComparison) -> Result<()> {
    let (ys, us, xs, vs) = match r.reference.contains("1000") {
        true => Reference::GhiaRe1000.centre_lines(),
        false => Reference::GhiaRe100.centre_lines(),
    };
    let ref_u: Vec<(f64, f64)> = ys.iter().zip(us.iter()).map(|(&y, &u)| (u, y)).collect();
    let ref_v: Vec<(f64, f64)> = xs.iter().zip(vs.iter()).map(|(&x, &v)| (x, v)).collect();
    let cu: Vec<(f64, f64)> = r.u_line.iter().map(|&(y, u)| (u, y)).collect();
    let cv = r.v_line.clone();
    std::fs::write(
        out.join("centreline_u.svg"),
        plot::chart(
            &format!("{name}: u on the vertical centre line"),
            "u / U",
            "y / L",
            &[
                Series {
                    label: "Foehncraft",
                    points: &cu,
                    markers: false,
                },
                Series {
                    label: "Ghia et al. 1982",
                    points: &ref_u,
                    markers: true,
                },
            ],
            false,
            false,
        ),
    )?;
    std::fs::write(
        out.join("centreline_v.svg"),
        plot::chart(
            &format!("{name}: v on the horizontal centre line"),
            "x / L",
            "v / U",
            &[
                Series {
                    label: "Foehncraft",
                    points: &cv,
                    markers: false,
                },
                Series {
                    label: "Ghia et al. 1982",
                    points: &ref_v,
                    markers: true,
                },
            ],
            false,
            false,
        ),
    )?;
    Ok(())
}

/// Output directory of a case: `[output] directory`, `--out`, or
/// `output/<case name>`.
/// Checks the case without solving it (`foehn flow --check`): builds the
/// mesh of the first level and prints the findings of the boundary-condition
/// check with their guidance; fails on an error.
pub fn check_case(case: &FlowCase, path: &Path, mut log: impl FnMut(&str)) -> Result<()> {
    let level = case.verification.levels.first().copied();
    let (mesh, geo, desc) = mesh_for(case, level)?;
    log(&format!("{}: {desc}", path.display()));
    // the property program on this mesh: regions name existing zones (D112)
    physics_on(case, &physics(case)?, &mesh, &geo)?;
    let findings = bc::check(
        &mesh,
        &geo,
        &case.boundary,
        case.rho(),
        bc::At {
            analytic: case.analytic().as_ref(),
            values: &case.compiled_boundaries()?,
            t: case.time.map_or(0.0, |t| t.start_time),
        },
        case.check_context(),
    );
    for f in &findings {
        log(&f.to_string());
    }
    let errors = findings.iter().filter(|f| f.level == Level::Error).count();
    if errors > 0 {
        bail!("{errors} error(s) in the boundary conditions");
    }
    log("the case can be run");
    Ok(())
}

/// The case name and output directory.
pub fn output_dir(case: &FlowCase, path: &Path, out: Option<PathBuf>) -> (String, PathBuf) {
    let name = path
        .file_stem()
        .map_or("flow".into(), |s| s.to_string_lossy().into_owned());
    let dir = out
        .or_else(|| case.output.directory.clone())
        .unwrap_or_else(|| PathBuf::from("output").join(&name));
    (name, dir)
}

/// Fails when the report lists failures.
pub fn check(rep: &FlowReport) -> Result<()> {
    if rep.failures.is_empty() {
        Ok(())
    } else {
        bail!(
            "{} check(s) failed:\n  {}",
            rep.failures.len(),
            rep.failures.join("\n  ")
        )
    }
}

impl FlowReport {
    /// Markdown report.
    pub fn to_markdown(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "# Flow case `{}`\n", self.name);
        if let Some(re) = self.reynolds {
            let _ = writeln!(s, "Reynolds number {re}\n");
        }
        if !self.findings.is_empty() {
            let _ = writeln!(s, "Boundary-condition check:\n");
            for f in &self.findings {
                let _ = writeln!(s, "- {f}");
            }
            s.push('\n');
        }
        let _ = writeln!(
            s,
            "| Mesh | Iterations | Converged | Largest residual | ns/cell/iteration | L2 error | Linf error | p error | Reference deviation |\n|---|---|---|---|---|---|---|---|---|"
        );
        for v in &self.solves {
            let _ = writeln!(
                s,
                "| {} | {} | {} | {:.2e} | {:.0} | {} | {} | {} | {} |",
                v.mesh,
                v.iterations,
                if v.converged { "yes" } else { "no" },
                v.residuals.max(),
                v.ns_per_cell_iteration,
                v.errors
                    .as_ref()
                    .map_or("–".into(), |e| format!("{:.3e}", e.l2)),
                v.errors
                    .as_ref()
                    .map_or("–".into(), |e| format!("{:.3e}", e.linf)),
                v.errors
                    .as_ref()
                    .map_or("–".into(), |e| format!("{:.3e}", e.pressure_linf)),
                v.reference.as_ref().map_or("–".into(), |r| format!(
                    "u {:.4}, v {:.4}",
                    r.u_max_dev, r.v_max_dev
                )),
            );
        }
        if !self.orders.is_empty() {
            let _ = writeln!(s, "\nObserved orders: {:.3?}", self.orders);
        }
        if !self.time_orders.is_empty() {
            let _ = writeln!(s, "\nObserved time orders: {:.3?}", self.time_orders);
        }
        if self.solves.iter().any(|v| v.transient.is_some()) {
            let _ = writeln!(
                s,
                "\n| Mesh | Scheme | Steps | End time (s) | Δt (s) | Largest Courant | Outer per step | CG per step | ns/cell/step |\n|---|---|---|---|---|---|---|---|---|"
            );
            for v in &self.solves {
                if let Some(t) = &v.transient {
                    let _ = writeln!(
                        s,
                        "| {} | {} | {} | {:.6e} | {:.3e}–{:.3e} | {:.3} | {:.2} | {:.1} | {:.0} |",
                        v.cells,
                        t.scheme,
                        t.steps,
                        t.time,
                        t.dt_min,
                        t.dt_max,
                        t.cfl_max,
                        t.outer_mean,
                        t.cg_mean,
                        t.ns_per_cell_step
                    );
                }
            }
        }
        if self
            .solves
            .iter()
            .any(|v| v.benchmark.is_some() || v.coefficients.is_some())
        {
            let _ = writeln!(s);
            for v in &self.solves {
                if let Some([cd, cl]) = v.coefficients {
                    let _ = writeln!(
                        s,
                        "{} cells: force coefficients C_D {cd:.5}, C_L {cl:.5}\n",
                        v.cells
                    );
                }
                if let Some(b) = &v.benchmark {
                    let _ = writeln!(
                        s,
                        "{} cells, {}{}:\n\n| Quantity | Computed | Reference | Distance |\n|---|---|---|---|",
                        v.cells,
                        b.benchmark,
                        if b.periods > 0 {
                            format!(" ({} periods)", b.periods)
                        } else {
                            String::new()
                        }
                    );
                    for x in &b.values {
                        let _ = writeln!(
                            s,
                            "| {} | {:.5} | {} – {} | {:.4} |",
                            x.range.name, x.value, x.range.lower, x.range.upper, x.deviation
                        );
                    }
                    s.push('\n');
                }
            }
        }
        for v in &self.solves {
            if let Some(b) = &v.balance {
                s.push_str(&balance::markdown(b, v.cells));
            }
        }
        if self.stopped {
            let _ = writeln!(s, "\nThe run was stopped by `.stop`.");
        }
        if let Some(v) = self.solves.last() {
            let _ = writeln!(
                s,
                "\n| Patch | Type | Mass flow (kg/s) | Mean pressure (Pa) | Force (N) |\n|---|---|---|---|---|"
            );
            for p in &v.patches {
                let _ = writeln!(
                    s,
                    "| {} | {} | {:.6e} | {:.6e} | [{:.4e}, {:.4e}, {:.4e}] |",
                    p.name,
                    p.kind,
                    p.mass_flow,
                    p.mean_pressure,
                    p.force[0],
                    p.force[1],
                    p.force[2]
                );
            }
            let t = &v.timings;
            let _ = writeln!(
                s,
                "\nSeconds per phase (R6): gradients {:.2}, momentum assembly {:.2}, momentum solve {:.2}, pressure assembly {:.2}, pressure solve {:.2}, multigrid setup {:.2} ({} setups)",
                t.gradients,
                t.momentum_assembly,
                t.momentum_solve,
                t.pressure_assembly,
                t.pressure_solve,
                t.amg_setup,
                v.amg_setups
            );
        }
        if self.failures.is_empty() {
            let _ = writeln!(s, "\nAll checks passed ({:.1} s).", self.seconds);
        } else {
            let _ = writeln!(s, "\nFailed checks:\n");
            for f in &self.failures {
                let _ = writeln!(s, "- {f}");
            }
        }
        s
    }
}

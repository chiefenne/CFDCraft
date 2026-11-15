//! The `ennecraft flow` case: steady incompressible laminar flow on one mesh,
//! or an order-of-accuracy sweep over refinement levels.
//!
//! Case file (TOML, D39). Every key is optional unless noted, SI units
//! (D44), unknown keys and values are errors with a "did you mean" hint:
//!
//! ```toml
//! [mesh]                         # as in `ennecraft mesh` (file or generate)
//! [mesh.generate]
//! type = "hex"
//! dimension = 2
//! cells = [64, 64, 1]
//! [fluid]
//! density = 1.0                  # ρ (kg/m³)
//! viscosity = 0.01               # μ (Pa s)
//! [flow]
//! convection = true              # false: Stokes flow
//! body_force = [0.0, 0.0, 0.0]   # N/m³
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
use crate::bc::{self, BoundaryCondition, Level};
use crate::plot::{self, Series};
use crate::simplec::{self, Flow, Physics, Residuals, Settings, Solution, Timings};
use crate::transient::TimeScheme;
use anyhow::{Context, Result, bail};
use ennecraft_mesh::case::{MeshSection, build_mesh};
use ennecraft_mesh::geometry::{dot, norm, sub};
use ennecraft_mesh::locate::Locator;
use ennecraft_mesh::renumber::renumber;
use ennecraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// `[fluid]`
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct FluidSection {
    /// Density ρ (kg/m³).
    pub density: f64,
    /// Dynamic viscosity μ (Pa s).
    pub viscosity: f64,
}

impl Default for FluidSection {
    fn default() -> Self {
        Self {
            density: 1.0,
            viscosity: 1e-3,
        }
    }
}

/// `[flow]`
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct FlowSection {
    /// Convection on (false: Stokes flow).
    pub convection: bool,
    /// Body force per unit volume (N/m³).
    pub body_force: [f64; 3],
}

impl Default for FlowSection {
    fn default() -> Self {
        Self {
            convection: true,
            body_force: [0.0; 3],
        }
    }
}

/// `[initial]`
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct InitialSection {
    /// Start velocity (m/s).
    pub velocity: [f64; 3],
    /// Start pressure (Pa).
    pub pressure: f64,
}

/// `[pressure]`
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PressureSection {
    /// Point whose cell holds the reference pressure (default: the cell
    /// nearest to the lower corner of the bounding box).
    pub reference_point: Option<[f64; 3]>,
    /// Pressure in the reference cell (Pa).
    pub reference_value: f64,
}

/// `[run]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RunSection {
    /// CPU threads (default: all cores; meshes below 20 000 cells run on
    /// one).
    pub threads: Option<usize>,
    /// Where the flow loop runs: `cpu` (the `cuda` backend of the flow loop
    /// follows with the GPU part of stage 5).
    pub backend: String,
    /// PTX module of the `cuda` backend.
    pub ptx: PathBuf,
    /// Restart file to continue from (D42).
    pub restart: Option<PathBuf>,
    /// Write `restart.h5` every N iterations or time steps (0: never).
    pub restart_every: usize,
    /// Also write a numbered `restart_<N>.h5` at these iterations or steps.
    pub restart_at: Vec<usize>,
    /// Write `restart.h5` at the end of the run.
    pub restart_final: bool,
    /// Re-read the case file every N iterations or time steps (0: only on a
    /// `.read_toml_N` control file), D40.
    pub reread_every: usize,
}

impl Default for RunSection {
    fn default() -> Self {
        Self {
            threads: None,
            backend: "cpu".into(),
            ptx: PathBuf::from("target/nvptx64-nvidia-cuda/release/ennecraft_ptx_kernels.ptx"),
            restart: None,
            restart_every: 0,
            restart_at: Vec::new(),
            restart_final: true,
            reread_every: 0,
        }
    }
}

/// `[time]`: a transient case (D28).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
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
        }
    }
}

impl TimeSection {
    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        anyhow::ensure!(
            self.time_step > 0.0,
            "[time] time_step must be > 0 (the first step with `adaptive`)"
        );
        anyhow::ensure!(
            self.end_time.is_some() || self.steps.is_some(),
            "[time] needs `end_time` or `steps`"
        );
        if let Some(e) = self.end_time {
            anyhow::ensure!(
                e > self.start_time,
                "[time] end_time must be later than start_time"
            );
        }
        anyhow::ensure!(self.max_cfl > 0.0, "[time] max_cfl must be > 0");
        anyhow::ensure!(
            self.max_growth > 1.0 && self.max_growth < 1.0 + std::f64::consts::SQRT_2,
            "[time] max_growth must lie in (1, 2.414): BDF2 with variable steps is zero-stable \
             only for step ratios below 1 + √2"
        );
        anyhow::ensure!(
            self.min_time_step >= 0.0,
            "[time] min_time_step must be >= 0"
        );
        if let Some(m) = self.max_time_step {
            anyhow::ensure!(
                m >= self.min_time_step && m > 0.0,
                "[time] max_time_step must be > 0 and >= min_time_step"
            );
        }
        Ok(())
    }
}

/// A point probe, `[[output.probes]]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Probe {
    /// Name (columns `<name>_u` … in `probes.csv`).
    pub name: String,
    /// The point.
    pub point: [f64; 3],
}

/// Force coefficients, `[output.coefficients]`: C = 2 F / (ρ U² A) of the
/// sum of the forces on `patches`.
#[derive(Debug, Clone, Serialize, Deserialize)]
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
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
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
#[derive(Debug, Clone, Serialize, Deserialize)]
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
#[derive(Debug, Clone, Serialize, Deserialize)]
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
    /// A row of `forces.csv` and `probes.csv` every N steps (transient).
    pub monitor_every: usize,
    /// Point probes: velocity and pressure per step (`probes.csv`).
    pub probes: Vec<Probe>,
    /// Force coefficients of some patches (`forces.csv`, the report).
    pub coefficients: Option<Coefficients>,
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
            monitor_every: 1,
            probes: Vec::new(),
            coefficients: None,
        }
    }
}

/// A flow case.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct FlowCase {
    /// The mesh.
    pub mesh: MeshSection,
    /// Fluid properties.
    pub fluid: FluidSection,
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
    /// The case file (for re-reading it during a run).
    #[serde(skip)]
    pub source: Option<PathBuf>,
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

impl FlowCase {
    /// Reads and validates a case file; relative paths are taken from its
    /// directory.
    pub fn from_file(path: &Path) -> Result<Self> {
        let text =
            std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        let mut c = Self::from_toml(&text, path.parent().unwrap_or(Path::new(".")))
            .with_context(|| format!("in case file {}", path.display()))?;
        c.source = Some(path.to_path_buf());
        Ok(c)
    }

    /// Parses and validates a case from TOML text.
    pub fn from_toml(text: &str, dir: &Path) -> Result<Self> {
        let mut c: Self = toml::from_str(text)
            .map_err(|e| anyhow::anyhow!(explain_toml_error(&e.to_string())))?;
        if let Some(f) = &c.mesh.file {
            c.mesh.file = Some(dir.join(f));
        }
        if let Some(r) = &c.run.restart {
            c.run.restart = Some(dir.join(r));
        }
        c.validate()?;
        Ok(c)
    }

    /// Ranges and consistency of the settings.
    pub fn validate(&self) -> Result<()> {
        anyhow::ensure!(self.fluid.density > 0.0, "[fluid] density must be > 0");
        anyhow::ensure!(self.fluid.viscosity > 0.0, "[fluid] viscosity must be > 0");
        let transient = self.time.is_some();
        self.solver.validate(transient)?;
        if let Some(t) = &self.time {
            t.validate()?;
        }
        anyhow::ensure!(
            self.output.history_every >= 1,
            "[output] history_every must be >= 1"
        );
        anyhow::ensure!(
            self.output.monitor_every >= 1,
            "[output] monitor_every must be >= 1"
        );
        if let Some(a) = &self.analytic
            && a.is_transient()
            && !transient
        {
            anyhow::bail!(
                "[analytic] type = \"taylor-green\" is transient and needs a [time] table"
            );
        }
        if let Some(c) = &self.output.coefficients {
            anyhow::ensure!(
                c.velocity > 0.0 && c.area > 0.0,
                "[output.coefficients] velocity and area must be > 0"
            );
            anyhow::ensure!(
                !c.patches.is_empty(),
                "[output.coefficients] patches must name at least one patch"
            );
        }
        let v = &self.verification;
        if !v.time_steps.is_empty() {
            anyhow::ensure!(transient, "[verification] time_steps need a [time] table");
            anyhow::ensure!(
                v.time_steps.len() >= 3,
                "[verification] time_steps need at least three steps for an observed order"
            );
            anyhow::ensure!(
                v.levels.len() < 2,
                "[verification] a sweep is either over levels or over time_steps"
            );
            anyhow::ensure!(
                self.time.is_some_and(|t| !t.adaptive),
                "[verification] time_steps need a fixed time step ([time] adaptive = false)"
            );
        }
        if let Some(b) = v.benchmark {
            let re =
                self.fluid.density * b.mean_velocity() * Benchmark::DIAMETER / self.fluid.viscosity;
            anyhow::ensure!(
                (re - b.reynolds()).abs() < 1e-6 * b.reynolds(),
                "[verification] benchmark {} needs Re = ρ Ū D / μ = {} (Ū = {} m/s, D = 0.1 m), \
                 the case has {re}",
                b.label(),
                b.reynolds(),
                b.mean_velocity()
            );
            if b == Benchmark::SchaferTurek2d2 {
                anyhow::ensure!(
                    transient,
                    "[verification] benchmark schafer-turek-2d2 is periodic and needs a [time] table"
                );
            }
        }
        let v = &self.verification;
        if v.levels.len() >= 2 {
            anyhow::ensure!(
                self.mesh.generate.is_some() || self.mesh.cylinder.is_some(),
                "[verification] levels need a [mesh.generate] or [mesh.cylinder] table"
            );
            anyhow::ensure!(
                self.analytic.is_some() || v.benchmark.is_some(),
                "[verification] levels need an [analytic] solution or a benchmark"
            );
        }
        if v.max_error.is_some() {
            anyhow::ensure!(
                self.analytic.is_some(),
                "[verification] max_error needs an [analytic] solution"
            );
        }
        if let Some(reference) = v.reference {
            let re = self.fluid.density * v.lid_velocity.unwrap_or(1.0) / self.fluid.viscosity;
            let want = reference.reynolds();
            anyhow::ensure!(
                (re - want).abs() < 1e-6 * want,
                "[verification] reference {:?} is for Re = {want}, the case has Re = ρ U / μ = {re} \
                 (unit cavity, lid_velocity = {})",
                reference,
                v.lid_velocity.unwrap_or(1.0)
            );
        }
        for l in &self.output.lines {
            anyhow::ensure!(
                l.points >= 2,
                "[[output.lines]] {}: points must be >= 2",
                l.name
            );
        }
        Ok(())
    }

    fn analytic(&self) -> Option<Analytic> {
        self.analytic
            .map(|a| a.with_fluid(self.fluid.density, self.fluid.viscosity))
    }

    /// The exact solution at time `t` (transient solutions decay in time).
    pub fn analytic_at(&self, t: f64) -> Option<Analytic> {
        self.analytic().map(|a| a.at_time(t))
    }

    /// Reynolds number ρ U L / μ: the lid speed of the cavity reference,
    /// the scales of the exact solution or the mean velocity and the
    /// cylinder diameter of the benchmark.
    pub fn reynolds(&self) -> Option<f64> {
        let rho = self.fluid.density;
        let mu = self.fluid.viscosity;
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
    /// The run was stopped by a `.stop` control file.
    pub stopped: bool,
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
    if case.time.is_some() {
        return crate::unsteady::run(case, name, out, log);
    }
    let t0 = Instant::now();
    std::fs::create_dir_all(out)?;
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
    let phys = Physics {
        rho: case.fluid.density,
        mu: case.fluid.viscosity,
        body_force: case.flow.body_force,
        convection: case.flow.convection,
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
        failures: Vec::new(),
        seconds: 0.0,
    };
    std::fs::write(
        out.join("case.resolved.toml"),
        toml::to_string_pretty(case)?,
    )?;
    let last = levels.len() - 1;
    let mut histories: Vec<(String, Vec<(f64, f64)>)> = Vec::new();
    for (li, level) in levels.iter().enumerate() {
        let (mesh, geo, desc) = mesh_for(case, *level)?;
        // boundary-condition check (D30)
        let findings = bc::check(
            &mesh,
            &geo,
            &case.boundary,
            phys.rho,
            analytic.as_ref(),
            true,
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
            header(&mut lg, case, &mesh, &desc, li == last);
        }
        lg.line(&format!("mesh  {desc}"));
        let arrays = bc::arrays(&mesh, &geo, &case.boundary, phys.rho, analytic.as_ref())?;
        let disc = simplec::discretisation(&mesh, &geo, &arrays)?;
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
            case.run.threads,
        )?;
        let mut csv = String::from(
            "iteration,momentum_x,momentum_y,momentum_z,continuity,bicgstab_x,bicgstab_y,bicgstab_z,cg\n",
        );
        let mut hist = Vec::new();
        let mut converged = false;
        let mut res = Residuals::default();
        let ts = Instant::now();
        for it in 1..=case.solver.max_iterations {
            res = fl.iterate()?;
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
            }
            if it == 1 || it % case.solver.log_every == 0 || done {
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
        }
        let secs = ts.elapsed().as_secs_f64();
        let iters = fl.iteration();
        lg.line(&format!(
            "      {} in {iters} outer iterations, {secs:.2} s, {:.0} ns/cell/iteration",
            if converged {
                "converged"
            } else {
                "NOT converged"
            },
            1e9 * fl.timings.total() / (nc as f64 * iters.max(1) as f64)
        ));
        let suffix = level.map_or(String::new(), |l| format!("_level{l}"));
        std::fs::write(out.join(format!("residuals{suffix}.csv")), &csv)?;
        histories.push((
            level.map_or("residual".into(), |l| format!("level {l}")),
            hist,
        ));
        if !converged && case.verification.require_convergence.unwrap_or(true) {
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
            phys: fl.phys,
        };
        let sol = &fl_host;
        let errors = analytic.map(|a| error_norms(sol, &mesh, &geo, &a, dimension(case)));
        if let Some(e) = &errors {
            lg.line(&format!(
                "      velocity error L2 {:.3e}, Linf {:.3e} (relative), pressure Linf {:.3e}",
                e.l2, e.linf, e.pressure_linf
            ));
        }
        let reference = case
            .verification
            .reference
            .map(|r| compare_reference(sol, &mesh, &geo, &loc, r, case));
        if let Some(r) = &reference {
            lg.line(&format!(
                "      {}: largest deviation u {:.4}, v {:.4} (lid speeds)",
                r.reference, r.u_max_dev, r.v_max_dev
            ));
        }
        let patches = patch_summaries(sol, &mesh, &geo, case);
        let total: f64 = patches.iter().map(|p| p.mass_flow).sum();
        // scale: the largest patch mass flow, or in a closed domain the
        // flow ρ |u|max A through the largest patch
        let umax = (0..nc).map(|c| norm(sol.u.get(c))).fold(0.0f64, f64::max);
        let big = patches
            .iter()
            .map(|p| p.mass_flow.abs().max(case.fluid.density * umax * p.area))
            .fold(0.0f64, f64::max);
        let mass_imbalance = if big > 0.0 { total / big } else { 0.0 };
        let coefficients = crate::monitor::coefficients(case, &patches);
        if let Some([cd, cl]) = coefficients {
            lg.line(&format!(
                "      force coefficients C_D {cd:.5}, C_L {cl:.5}"
            ));
        }
        let benchmark =
            crate::unsteady::steady_benchmark_of(case, sol, &mesh, &geo, &loc, &patches);
        if let Some(b) = &benchmark {
            crate::unsteady::log_benchmark(b, |l| lg.line(l));
        }
        if li == last {
            write_outputs(
                case,
                out,
                &mesh,
                &geo,
                sol,
                analytic.as_ref(),
                &patches,
                &loc,
            )?;
            if let Some(r) = &reference {
                write_reference_svg(out, name, r)?;
            }
        }
        rep.solves.push(Solve {
            mesh: desc,
            cells: nc,
            iterations: iters,
            converged,
            residuals: res,
            ns_per_cell_iteration: 1e9 * fl.timings.total() / (nc as f64 * iters.max(1) as f64),
            timings: fl.timings.clone(),
            amg_setups: fl.amg_setups,
            errors,
            reference,
            patches,
            mass_imbalance,
            coefficients,
            transient: None,
            benchmark,
        });
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
                        label: "Ennecraft",
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
        "Ennecraft {} flow solver (SIMPLEC, steady, laminar)",
        env!("CARGO_PKG_VERSION")
    ));
    lg.line(&format!(
        "fluid  ρ = {} kg/m³, μ = {} Pa s; convection {}",
        case.fluid.density,
        case.fluid.viscosity,
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
        lg.line("  patches.csv, flow.json, flow.md, residuals.svg, case.resolved.toml");
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
    let mu = fl.phys.mu;
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
) -> Result<()> {
    let mut csv = String::from("patch,type,area,mass_flow,mean_pressure,force_x,force_y,force_z\n");
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
            let mut csv = String::from("x,y,z,pressure,tau_x,tau_y,tau_z\n");
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
        let mut csv = String::from("s,x,y,z,u,v,w,p\n");
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
    write_fields(case, &out.join("fields.vtkhdf"), mesh, geo, fl, analytic)
}

#[cfg(feature = "hdf5")]
pub(crate) fn write_fields(
    case: &FlowCase,
    file: &Path,
    mesh: &Mesh,
    geo: &Geometry,
    fl: &Solution,
    analytic: Option<&Analytic>,
) -> Result<()> {
    use ennecraft_mesh::vtkhdf::{CellField, write};
    if !case.output.vtkhdf {
        return Ok(());
    }
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
    write(file, mesh, &fields)
}

#[cfg(not(feature = "hdf5"))]
pub(crate) fn write_fields(
    _: &FlowCase,
    _: &Path,
    _: &Mesh,
    _: &Geometry,
    _: &Solution,
    _: Option<&Analytic>,
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
                    label: "Ennecraft",
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
                    label: "Ennecraft",
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
/// Checks the case without solving it (`ennecraft flow --check`): builds the
/// mesh of the first level and prints the findings of the boundary-condition
/// check with their guidance; fails on an error.
pub fn check_case(case: &FlowCase, path: &Path, mut log: impl FnMut(&str)) -> Result<()> {
    let level = case.verification.levels.first().copied();
    let (mesh, geo, desc) = mesh_for(case, level)?;
    log(&format!("{}: {desc}", path.display()));
    let findings = bc::check(
        &mesh,
        &geo,
        &case.boundary,
        case.fluid.density,
        case.analytic().as_ref(),
        true,
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

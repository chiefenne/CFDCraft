//! The `cfdcraft poisson` case: a study of the pressure solver (stage 3)
//! on a Poisson problem −∇·(Γ∇φ) = s, the model of the pressure equation.
//!
//! On every refinement level the system is assembled once (orthogonal
//! diffusion, no convection; the matrix is symmetric positive definite) and
//! solved from φ = 0 to ‖b − Aφ‖ ≤ tol ‖b‖ by every solver of the list on
//! every backend:
//!
//! - `amg`: flexible CG preconditioned by the smoothed-aggregation V-cycle
//!   with the settings of `[linear.amg]` (default: mixed precision);
//! - `amg-double`: the same hierarchy in float64;
//! - `amg-plain`: plain (unsmoothed) aggregation;
//! - `jacobi`: CG with Jacobi preconditioning;
//! - `hypre`: hypre's PCG with BoomerAMG (feature `hypre`, D20), once per
//!   level on its own device.
//!
//! The report gives iterations, setup and solve times, ns per cell and CG
//! iteration (R6), the hierarchy, the true final residual and the largest
//! difference to the reference solution (hypre when it is in the list, else
//! the first solve of the level). Checks: iteration growth from the
//! coarsest to the finest level, agreement with the reference, bitwise
//! equality of each solver across backends.
//!
//! ```toml
//! [mesh.generate]
//! type = "poly"
//! cells = [4, 4, 4]
//! [problem]                     # velocity must be zero
//! diffusivity = 1.0
//! manufactured = {}
//! [problem.boundary.xmin]       # ... one table per patch
//! type = "fixed-value"
//! profile = "manufactured"
//! [linear]
//! tolerance = 1e-10
//! [linear.amg]                  # crate::amg::AmgSettings
//! smoother = "chebyshev"
//! [study]
//! levels = [1, 2, 3, 4]
//! solvers = ["amg", "amg-double", "jacobi", "hypre"]
//! backends = ["cpu", "cuda"]
//! max_iteration_ratio = 1.5     # amg solvers: finest / coarsest level
//! [hypre]                       # BoomerAMG settings
//! ```

use crate::amg::{Aggregation, AmgSettings, HierarchyInfo, Precision};
use crate::backend::{self, Backend, FaceSettings, V};
use crate::disc::Discretisation;
use crate::problem::Problem;
use crate::solver::{self, Preconditioner};
use anyhow::{Context, Result, bail, ensure};
use cfdcraft_mesh::case::{MeshSection, build_mesh};
use cfdcraft_mesh::renumber::renumber;
use cfdcraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// A Poisson study case.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PoissonCase {
    /// Mesh (file or generated).
    pub mesh: MeshSection,
    /// The problem (velocity zero).
    pub problem: Problem,
    /// Linear solver settings.
    pub linear: LinearSection,
    /// The study.
    pub study: StudySection,
    /// hypre's settings.
    pub hypre: HypreSettings,
    /// Output.
    pub output: PoissonOutput,
}

/// `[linear]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct LinearSection {
    /// Relative residual tolerance ‖b − Aφ‖ / ‖b‖.
    pub tolerance: f64,
    /// Largest number of CG iterations.
    pub max_iterations: usize,
    /// Multigrid settings.
    pub amg: AmgSettings,
}

impl Default for LinearSection {
    fn default() -> Self {
        Self {
            tolerance: 1e-10,
            max_iterations: 2000,
            amg: AmgSettings::default(),
        }
    }
}

/// A solver of the study.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum StudySolver {
    /// CG + AMG with `[linear.amg]`.
    Amg,
    /// CG + AMG in float64.
    AmgDouble,
    /// CG + plain aggregation.
    AmgPlain,
    /// CG + Jacobi.
    Jacobi,
    /// hypre PCG + BoomerAMG.
    Hypre,
}

impl StudySolver {
    /// Name in reports.
    pub fn name(self) -> &'static str {
        match self {
            StudySolver::Amg => "amg",
            StudySolver::AmgDouble => "amg-double",
            StudySolver::AmgPlain => "amg-plain",
            StudySolver::Jacobi => "jacobi",
            StudySolver::Hypre => "hypre",
        }
    }
}

impl std::str::FromStr for StudySolver {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self> {
        Ok(match s {
            "amg" => StudySolver::Amg,
            "amg-double" => StudySolver::AmgDouble,
            "amg-plain" => StudySolver::AmgPlain,
            "jacobi" => StudySolver::Jacobi,
            "hypre" => StudySolver::Hypre,
            _ => bail!("unknown solver `{s}` (amg, amg-double, amg-plain, jacobi, hypre)"),
        })
    }
}

/// `[study]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct StudySection {
    /// Refinement levels of `[mesh.generate]` (empty: the mesh as given).
    pub levels: Vec<u32>,
    /// Solvers.
    pub solvers: Vec<StudySolver>,
    /// Backends of the in-house solvers: `cpu`, `cpu-serial`, `cuda`.
    pub backends: Vec<String>,
    /// CPU threads (default: all cores).
    pub threads: Option<usize>,
    /// PTX module of the GPU backend.
    pub rust_ptx: PathBuf,
    /// Largest ratio of the iterations on the finest and the coarsest level
    /// for the solvers in `check_iterations`.
    pub max_iteration_ratio: Option<f64>,
    /// Solvers whose iteration growth is checked (default: `amg`,
    /// `amg-double`).
    pub check_iterations: Vec<StudySolver>,
    /// Largest max-norm difference to the reference solution, relative to
    /// the max norm of the reference.
    pub max_difference: f64,
    /// Repetitions of each solve (the fastest is reported).
    pub repeat: usize,
}

impl Default for StudySection {
    fn default() -> Self {
        Self {
            levels: Vec::new(),
            solvers: vec![StudySolver::Amg, StudySolver::Jacobi],
            backends: vec!["cpu".into()],
            threads: None,
            rust_ptx: "target/nvptx64-nvidia-cuda/release/cfdcraft_ptx_kernels.ptx".into(),
            max_iteration_ratio: None,
            check_iterations: vec![StudySolver::Amg, StudySolver::AmgDouble],
            max_difference: 1e-6,
            repeat: 1,
        }
    }
}

/// `[output]`
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PoissonOutput {
    /// Output directory (default `output/<case name>`).
    pub directory: Option<PathBuf>,
    /// Write the solution of the finest level (VTKHDF).
    pub vtkhdf: bool,
}

/// BoomerAMG settings of the case file (`[hypre]`). The defaults are a
/// common choice for the pressure equation that runs on the CPU and the GPU:
/// HMIS coarsening, extended+i interpolation with at most 4 entries per
/// row, one level of aggressive coarsening, symmetric hybrid Gauss–Seidel on
/// the CPU (l1-Jacobi, `relax_type = 18`, on the GPU).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct HypreSettings {
    /// Strength threshold (0.25 is hypre's default; 0.5 is often advised in
    /// 3D).
    pub strong_threshold: f64,
    /// Coarsening (8 PMIS, 10 HMIS).
    pub coarsen_type: i32,
    /// Interpolation (6 extended+i, 14 extended).
    pub interp_type: i32,
    /// Relaxation; `None`: 6 (symmetric hybrid Gauss–Seidel) on the CPU, 18
    /// (l1-Jacobi) on the GPU.
    pub relax_type: Option<i32>,
    /// Sweeps per level.
    pub num_sweeps: i32,
    /// Levels of aggressive coarsening.
    pub agg_num_levels: i32,
    /// Largest number of interpolation entries per row.
    pub pmax: i32,
    /// Largest coarsest level.
    pub max_coarse: i32,
    /// hypre's print level (0 silent).
    pub print_level: i32,
}

impl Default for HypreSettings {
    fn default() -> Self {
        Self {
            strong_threshold: 0.5,
            coarsen_type: 10,
            interp_type: 6,
            relax_type: None,
            num_sweeps: 1,
            agg_num_levels: 1,
            pmax: 4,
            max_coarse: 256,
            print_level: 0,
        }
    }
}

impl PoissonCase {
    /// Reads a case file; relative paths are taken from its directory.
    pub fn from_file(path: &Path) -> Result<Self> {
        let text =
            std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        let mut c: Self =
            toml::from_str(&text).with_context(|| format!("in case file {}", path.display()))?;
        let dir = path.parent().unwrap_or(Path::new("."));
        if let Some(f) = &c.mesh.file {
            c.mesh.file = Some(dir.join(f));
        }
        c.problem.link();
        c.validate()?;
        Ok(c)
    }

    /// Checks the settings.
    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.problem.velocity == [0.0; 3],
            "[problem] velocity must be zero in a Poisson study (a symmetric matrix for CG)"
        );
        ensure!(
            self.linear.tolerance > 0.0 && self.linear.tolerance < 1.0,
            "[linear] tolerance must lie in (0, 1)"
        );
        self.linear.amg.validate()?;
        ensure!(!self.study.solvers.is_empty(), "[study] solvers is empty");
        ensure!(!self.study.backends.is_empty(), "[study] backends is empty");
        ensure!(self.study.repeat >= 1, "[study] repeat must be >= 1");
        if self.study.solvers.contains(&StudySolver::Hypre) && !cfg!(feature = "hypre") {
            bail!(
                "[study] solvers lists `hypre`, but this build has no hypre; build with \
                 `--features cfdcraft-fv/hypre` and HYPRE_DIR set (recipe 17)"
            );
        }
        if self.study.levels.len() >= 2 {
            ensure!(
                self.mesh.generate.is_some(),
                "[study] levels need a [mesh.generate] table"
            );
        }
        Ok(())
    }

    fn amg_settings(&self, s: StudySolver) -> AmgSettings {
        let mut a = self.linear.amg.clone();
        match s {
            StudySolver::AmgDouble => a.precision = Precision::Double,
            StudySolver::AmgPlain => a.aggregation = Aggregation::Plain,
            _ => {}
        }
        a
    }
}

/// One solve of the study.
#[derive(Debug, Clone, Serialize)]
pub struct PoissonSolve {
    /// Refinement level.
    pub level: Option<u32>,
    /// Cells.
    pub cells: usize,
    /// Solver.
    pub solver: String,
    /// Backend (`hypre` for hypre).
    pub backend: String,
    /// Device.
    pub device: String,
    /// Iterations.
    pub iterations: usize,
    /// ‖b − Aφ‖ / ‖b‖ recomputed in float64 after the solve.
    pub true_residual: f64,
    /// Seconds of the setup (host for the in-house AMG).
    pub setup_seconds: f64,
    /// Seconds of the solve.
    pub solve_seconds: f64,
    /// Figure of merit (R6): ns per cell and iteration of the solve.
    pub ns_per_cell_iteration: f64,
    /// Hierarchy of the in-house AMG.
    pub hierarchy: Option<HierarchyInfo>,
    /// max |φ − φ_ref| / max |φ_ref|.
    pub diff_to_reference: Option<f64>,
    /// Bitwise equal to the same solver on the first backend.
    pub bitwise_to_first_backend: Option<bool>,
}

/// Report of a study.
#[derive(Debug, Clone, Default, Serialize)]
pub struct PoissonReport {
    /// Case name.
    pub name: String,
    /// Mesh per level.
    pub meshes: Vec<String>,
    /// Reference of the differences.
    pub reference: String,
    /// All solves.
    pub solves: Vec<PoissonSolve>,
    /// Iteration ratio finest / coarsest per checked solver.
    pub iteration_ratios: Vec<(String, f64)>,
    /// Failed checks.
    pub failures: Vec<String>,
    /// Seconds for everything.
    pub seconds: f64,
}

fn mesh_for(case: &PoissonCase, level: Option<u32>) -> Result<(Mesh, Geometry, String)> {
    let mut sec = case.mesh.clone();
    if let (Some(l), Some(g)) = (level, sec.generate.as_mut()) {
        g.refine = l;
    }
    let (mut mesh, source, _) = build_mesh(&sec)?;
    renumber(&mut mesh, sec.renumber, sec.seed).map_err(anyhow::Error::msg)?;
    let geo = Geometry::compute(&mesh);
    let desc = format!("{source}: {} cells", mesh.n_cells);
    Ok((mesh, geo, desc))
}

/// Assembles the Poisson system in a backend: orthogonal diffusion at
/// φ = 0, no relaxation.
pub fn assemble(be: &mut dyn Backend) -> Result<()> {
    let n = be.n_cells();
    be.upload(V::Phi, &vec![0.0; n])?;
    be.gradient()?;
    be.faces(FaceSettings {
        blend: 0.0,
        nonorth: false,
        nonorth_limit: -1.0,
    })?;
    be.assemble(1.0)?;
    be.sync()
}

/// Solves the assembled system from φ = 0 with an in-house solver; returns
/// (iterations, setup seconds, solve seconds, hierarchy).
pub fn solve_inhouse(
    be: &mut dyn Backend,
    pc: Preconditioner,
    amg: &AmgSettings,
    tol: f64,
    max_iter: usize,
) -> Result<(usize, f64, f64, Option<HierarchyInfo>)> {
    let n = be.n_cells();
    let (info, setup) = if pc == Preconditioner::Amg {
        let t = Instant::now();
        let i = solver::setup_amg(be, amg)?;
        (Some(i), t.elapsed().as_secs_f64())
    } else {
        (None, 0.0)
    };
    be.upload(V::Phi, &vec![0.0; n])?;
    be.copy(V::Rhs, V::R)?;
    let r0 = be.dot(V::R, V::R)?.sqrt();
    be.sync()?;
    let t = Instant::now();
    let (it, _) = solver::cg(be, pc, max_iter, tol, r0)?;
    be.sync()?;
    Ok((it, setup, t.elapsed().as_secs_f64(), info))
}

fn true_residual(be: &mut dyn Backend) -> Result<f64> {
    be.residual(V::Phi, V::Rhs, V::R)?;
    let r = be.dot(V::R, V::R)?.sqrt();
    let b = be.dot(V::Rhs, V::Rhs)?.sqrt();
    Ok(r / b.max(f64::MIN_POSITIVE))
}

fn max_rel_diff(a: &[f64], reference: &[f64]) -> f64 {
    let scale = reference
        .iter()
        .fold(0.0f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    a.iter()
        .zip(reference)
        .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()))
        / scale
}

/// Runs a study; `log` receives progress lines; files go to `out`.
pub fn run(
    case: &PoissonCase,
    name: &str,
    out: &Path,
    mut log: impl FnMut(&str),
) -> Result<PoissonReport> {
    let t0 = Instant::now();
    std::fs::create_dir_all(out)?;
    let levels: Vec<Option<u32>> = if case.study.levels.is_empty() {
        vec![None]
    } else {
        case.study.levels.iter().map(|&l| Some(l)).collect()
    };
    let has_hypre = case.study.solvers.contains(&StudySolver::Hypre);
    let mut rep = PoissonReport {
        name: name.into(),
        reference: if has_hypre {
            "hypre".into()
        } else {
            format!(
                "{} on {}",
                case.study.solvers[0].name(),
                case.study.backends[0]
            )
        },
        ..PoissonReport::default()
    };
    let tol = case.linear.tolerance;
    let maxit = case.linear.max_iterations;
    for (li, &level) in levels.iter().enumerate() {
        let (mesh, geo, desc) = mesh_for(case, level)?;
        log(&format!("level {}: {desc}", level.unwrap_or(0)));
        rep.meshes.push(desc);
        let disc = Discretisation::new(&mesh, &geo, &case.problem)?;
        let n = mesh.n_cells;
        // solutions of this level: (solver, backend) -> φ
        let mut sols: Vec<(StudySolver, String, Vec<f64>)> = Vec::new();
        let mut first_solves = Vec::new();
        let mut matrix = None;
        for be_name in &case.study.backends {
            let mut be = backend::create(be_name, &disc, case.study.threads, &case.study.rust_ptx)
                .with_context(|| format!("backend {be_name}"))?;
            assemble(be.as_mut())?;
            if matrix.is_none() && has_hypre {
                let val = be.matrix_values()?;
                matrix = Some((disc.matrix_csr(&val), be.download(V::Rhs)?));
            }
            let device = be.device();
            for &s in &case.study.solvers {
                if s == StudySolver::Hypre {
                    continue;
                }
                let pc = if s == StudySolver::Jacobi {
                    Preconditioner::Jacobi
                } else {
                    Preconditioner::Amg
                };
                let amg = case.amg_settings(s);
                let mut best: Option<(usize, f64, f64, Option<HierarchyInfo>)> = None;
                for _ in 0..case.study.repeat {
                    let r = solve_inhouse(be.as_mut(), pc, &amg, tol, maxit)?;
                    best = Some(match best {
                        Some(b) if b.1 + b.2 <= r.1 + r.2 => b,
                        _ => r,
                    });
                }
                let (it, setup, solve, info) = best.expect("repeat >= 1");
                let tr = true_residual(be.as_mut())?;
                let phi = be.download(V::Phi)?;
                let bitwise = sols
                    .iter()
                    .find(|(ss, _, _)| *ss == s)
                    .map(|(_, _, p)| p.iter().zip(&phi).all(|(a, b)| a.to_bits() == b.to_bits()));
                if bitwise == Some(false) {
                    rep.failures.push(format!(
                        "level {}: {} on {be_name} is not bitwise equal to {}",
                        level.unwrap_or(0),
                        s.name(),
                        case.study.backends[0]
                    ));
                }
                if it >= maxit {
                    rep.failures.push(format!(
                        "level {}: {} on {be_name} did not converge in {maxit} iterations",
                        level.unwrap_or(0),
                        s.name()
                    ));
                }
                log(&format!(
                    "  {:10} {:10} {:4} it, true residual {tr:.2e}, setup {setup:.3} s, solve {solve:.3} s, {:.1} ns/cell/it{}",
                    s.name(),
                    be_name,
                    it,
                    1e9 * solve / (n as f64 * it.max(1) as f64),
                    match bitwise {
                        Some(true) => ", bitwise equal to the first backend",
                        Some(false) => ", NOT bitwise equal to the first backend",
                        None => "",
                    }
                ));
                rep.solves.push(PoissonSolve {
                    level,
                    cells: n,
                    solver: s.name().into(),
                    backend: be_name.clone(),
                    device: device.clone(),
                    iterations: it,
                    true_residual: tr,
                    setup_seconds: setup,
                    solve_seconds: solve,
                    ns_per_cell_iteration: 1e9 * solve / (n as f64 * it.max(1) as f64),
                    hierarchy: info,
                    diff_to_reference: None,
                    bitwise_to_first_backend: bitwise,
                });
                first_solves.push(rep.solves.len() - 1);
                sols.push((s, be_name.clone(), phi));
            }
        }
        let mut reference = None;
        if let Some((a, b)) = matrix.take() {
            let (sv, x) = run_hypre(case, level, &a, &b, &mut log)?;
            if sv.iterations >= maxit {
                rep.failures.push(format!(
                    "level {}: hypre did not converge in {maxit} iterations",
                    level.unwrap_or(0)
                ));
            }
            rep.solves.push(sv);
            reference = Some(x);
        }
        let reference = reference.unwrap_or_else(|| sols[0].2.clone());
        for (k, (s, b, phi)) in first_solves.iter().zip(&sols) {
            let d = max_rel_diff(phi, &reference);
            rep.solves[*k].diff_to_reference = Some(d);
            if d > case.study.max_difference {
                rep.failures.push(format!(
                    "level {}: {} on {b} differs from the reference ({}) by {d:.2e}",
                    level.unwrap_or(0),
                    s.name(),
                    rep.reference
                ));
            }
        }
        if li + 1 == levels.len() && case.output.vtkhdf {
            write_solution(out, &mesh, &sols[0].2)?;
        }
    }
    // iteration growth
    if levels.len() >= 2 {
        for &s in &case.study.check_iterations {
            if !case.study.solvers.contains(&s) {
                continue;
            }
            let its: Vec<usize> = rep
                .solves
                .iter()
                .filter(|v| v.solver == s.name() && v.backend == case.study.backends[0])
                .map(|v| v.iterations)
                .collect();
            if let (Some(&a), Some(&b)) = (its.first(), its.last()) {
                let ratio = b as f64 / a.max(1) as f64;
                rep.iteration_ratios.push((s.name().into(), ratio));
                if let Some(mx) = case.study.max_iteration_ratio
                    && ratio > mx
                {
                    rep.failures.push(format!(
                        "{}: iterations grow from {a} to {b} (ratio {ratio:.2} > {mx})",
                        s.name()
                    ));
                }
            }
        }
        std::fs::write(
            out.join("iterations.svg"),
            iterations_svg(name, &rep, &case.study.backends[0]),
        )?;
    }
    rep.seconds = t0.elapsed().as_secs_f64();
    std::fs::write(
        out.join("poisson.json"),
        serde_json::to_string_pretty(&rep)?,
    )?;
    std::fs::write(out.join("poisson.md"), rep.to_markdown())?;
    Ok(rep)
}

/// Solves the system once with hypre.
#[cfg(feature = "hypre")]
fn run_hypre(
    case: &PoissonCase,
    level: Option<u32>,
    a: &crate::amg::Csr,
    b: &[f64],
    log: &mut impl FnMut(&str),
) -> Result<(PoissonSolve, Vec<f64>)> {
    let n = a.n_rows;
    let h = crate::hypre::solve(
        a,
        b,
        case.linear.tolerance,
        case.linear.max_iterations,
        &case.hypre,
    )?;
    let ax = a.spmv(&h.x);
    let rr: f64 = b.iter().zip(&ax).map(|(b, y)| (b - y) * (b - y)).sum();
    let bb: f64 = b.iter().map(|v| v * v).sum();
    let tr = (rr / bb.max(f64::MIN_POSITIVE)).sqrt();
    let ns = 1e9 * h.solve_seconds / (n as f64 * h.iterations.max(1) as f64);
    log(&format!(
        "  {:10} {:10} {:4} it, true residual {tr:.2e}, setup {:.3} s, solve {:.3} s, {ns:.1} ns/cell/it",
        "hypre", "hypre", h.iterations, h.setup_seconds, h.solve_seconds,
    ));
    Ok((
        PoissonSolve {
            level,
            cells: n,
            solver: "hypre".into(),
            backend: "hypre".into(),
            device: crate::hypre::device(),
            iterations: h.iterations,
            true_residual: tr,
            setup_seconds: h.setup_seconds,
            solve_seconds: h.solve_seconds,
            ns_per_cell_iteration: ns,
            hierarchy: None,
            diff_to_reference: Some(0.0),
            bitwise_to_first_backend: None,
        },
        h.x,
    ))
}

#[cfg(not(feature = "hypre"))]
fn run_hypre(
    _: &PoissonCase,
    _: Option<u32>,
    _: &crate::amg::Csr,
    _: &[f64],
    _: &mut impl FnMut(&str),
) -> Result<(PoissonSolve, Vec<f64>)> {
    bail!("this build has no hypre (feature `hypre`)")
}

#[cfg(feature = "hdf5")]
fn write_solution(out: &Path, mesh: &Mesh, phi: &[f64]) -> Result<()> {
    use cfdcraft_mesh::vtkhdf::{CellField, write};
    write(
        out.join("fields.vtkhdf"),
        mesh,
        &[CellField {
            name: "phi",
            values: phi,
        }],
    )
}

#[cfg(not(feature = "hdf5"))]
fn write_solution(_out: &Path, _mesh: &Mesh, _phi: &[f64]) -> Result<()> {
    bail!("this build has no HDF5 (VTKHDF output); set [output] vtkhdf = false")
}

impl PoissonReport {
    /// Markdown report.
    pub fn to_markdown(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "# Poisson study `{}`\n", self.name);
        for m in &self.meshes {
            let _ = writeln!(s, "- mesh: {m}");
        }
        let _ = writeln!(s, "- reference of the differences: {}", self.reference);
        for (n, r) in &self.iteration_ratios {
            let _ = writeln!(s, "- iterations finest / coarsest level, {n}: {r:.2}");
        }
        let _ = writeln!(s, "- total time: {:.2} s\n", self.seconds);
        let _ = writeln!(
            s,
            "| Cells | Solver | Backend | Device | Levels | Op. complexity | Iterations | True residual | Setup (s) | Solve (s) | ns/cell/it. | vs reference | Bitwise vs first backend |"
        );
        let _ = writeln!(
            s,
            "|---:|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"
        );
        for v in &self.solves {
            let (lv, oc) = v.hierarchy.as_ref().map_or(("–".into(), "–".into()), |h| {
                (
                    h.rows.len().to_string(),
                    format!("{:.2}", h.operator_complexity),
                )
            });
            let _ = writeln!(
                s,
                "| {} | {} | {} | {} | {lv} | {oc} | {} | {:.2e} | {:.3} | {:.3} | {:.1} | {} | {} |",
                v.cells,
                v.solver,
                v.backend,
                v.device,
                v.iterations,
                v.true_residual,
                v.setup_seconds,
                v.solve_seconds,
                v.ns_per_cell_iteration,
                v.diff_to_reference
                    .map_or("–".into(), |d| format!("{d:.1e}")),
                match v.bitwise_to_first_backend {
                    Some(true) => "yes",
                    Some(false) => "**no**",
                    None => "–",
                }
            );
        }
        if self.failures.is_empty() {
            let _ = writeln!(s, "\nAll checks passed.");
        } else {
            let _ = writeln!(s, "\n**Failed checks:**\n");
            for f in &self.failures {
                let _ = writeln!(s, "- {f}");
            }
        }
        s
    }
}

const SVG_HEAD: &str = r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 440" width="640" height="440" font-family="Helvetica, Arial, sans-serif" font-size="12"><rect width="640" height="440" fill="#ffffff"/>"##;
const COLOURS: [&str; 5] = ["#2a62c9", "#d1495b", "#2e933c", "#e09f3e", "#6c4ab6"];

/// Iterations over the number of cells (log-log), one curve per solver.
fn iterations_svg(name: &str, rep: &PoissonReport, backend: &str) -> String {
    let (x0, x1, y0, y1) = (70.0, 600.0, 380.0, 50.0);
    let mut series: Vec<(String, Vec<(f64, f64)>)> = Vec::new();
    for v in &rep.solves {
        if v.backend != backend && v.backend != "hypre" {
            continue;
        }
        let p = (
            (v.cells as f64).log10(),
            (v.iterations.max(1) as f64).log10(),
        );
        match series.iter_mut().find(|(n, _)| *n == v.solver) {
            Some((_, pts)) => pts.push(p),
            None => series.push((v.solver.clone(), vec![p])),
        }
    }
    let all = series.iter().flat_map(|(_, p)| p.iter());
    let (mut cmin, mut cmax, mut imax) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64);
    for &(c, i) in all {
        cmin = cmin.min(c);
        cmax = cmax.max(c);
        imax = imax.max(i);
    }
    let (cmin, cmax) = (cmin.floor(), cmax.ceil().max(cmin.floor() + 1.0));
    let (imin, imax) = (0.0, imax.ceil().max(1.0));
    let px = |v: f64| x0 + (v - cmin) / (cmax - cmin) * (x1 - x0);
    let py = |v: f64| y0 - (v - imin) / (imax - imin) * (y0 - y1);
    let mut s = String::from(SVG_HEAD);
    let _ = write!(
        s,
        r##"<text x="320" y="28" text-anchor="middle" font-size="15">{name}: CG iterations to the tolerance</text><rect x="{x0}" y="{y1}" width="{}" height="{}" fill="none" stroke="#444"/>"##,
        x1 - x0,
        y0 - y1
    );
    let mut d = cmin;
    while d <= cmax + 1e-9 {
        let _ = write!(
            s,
            r##"<text x="{:.1}" y="{}" text-anchor="middle">1e{}</text>"##,
            px(d),
            y0 + 18.0,
            d as i32
        );
        d += 1.0;
    }
    for e in 0..=(imax as i32) {
        for m in [1.0f64, 2.0, 5.0] {
            let v = (m * 10f64.powi(e)).log10();
            if v > imax + 1e-9 {
                continue;
            }
            let _ = write!(
                s,
                r##"<line x1="{x0}" x2="{x1}" y1="{y:.1}" y2="{y:.1}" stroke="#ddd"/><text x="{}" y="{:.1}" text-anchor="end">{}</text>"##,
                x0 - 6.0,
                py(v) + 4.0,
                m * 10f64.powi(e),
                y = py(v)
            );
        }
    }
    for (k, (n, pts)) in series.iter().enumerate() {
        let c = COLOURS[k % COLOURS.len()];
        let line: Vec<String> = pts
            .iter()
            .map(|(a, b)| format!("{:.1},{:.1}", px(*a), py(*b)))
            .collect();
        let _ = write!(
            s,
            r##"<polyline points="{}" fill="none" stroke="{c}" stroke-width="2"/>"##,
            line.join(" ")
        );
        for (a, b) in pts {
            let _ = write!(
                s,
                r##"<circle cx="{:.1}" cy="{:.1}" r="4" fill="{c}"/>"##,
                px(*a),
                py(*b)
            );
        }
        let _ = write!(
            s,
            r##"<rect x="{}" y="{}" width="14" height="4" fill="{c}"/><text x="{}" y="{}">{n}</text>"##,
            x0 + 12.0,
            y1 + 14.0 + 18.0 * k as f64,
            x0 + 32.0,
            y1 + 19.0 + 18.0 * k as f64
        );
    }
    let _ = write!(
        s,
        r##"<text x="320" y="420" text-anchor="middle">cells (log scale); in-house solvers on the {backend} backend</text></svg>"##
    );
    s
}

/// Case name and output directory of a case file.
pub fn output_dir(case: &PoissonCase, path: &Path, out: Option<PathBuf>) -> (String, PathBuf) {
    let name = path
        .file_stem()
        .map_or("poisson".into(), |s| s.to_string_lossy().into_owned());
    let dir = out
        .or_else(|| case.output.directory.clone())
        .unwrap_or_else(|| PathBuf::from("output").join(&name));
    (name, dir)
}

/// Fails with the list of failed checks.
pub fn check(rep: &PoissonReport) -> Result<()> {
    if rep.failures.is_empty() {
        Ok(())
    } else {
        bail!("checks failed:\n  {}", rep.failures.join("\n  "))
    }
}

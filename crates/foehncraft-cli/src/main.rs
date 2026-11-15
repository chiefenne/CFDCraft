//! `foehn`: the command line of Foehncraft (D4).
//!
//! Commands: `info` (build and platform), `mesh` (build or read a mesh,
//! check, partition and write it; stage 1), `scalar` (steady scalar
//! transport, verification sweeps, CPU and GPU; stage 2), `poisson` (the
//! pressure solver study: CG with algebraic multigrid against hypre; stage
//! 3), `flow` (steady incompressible laminar flow with SIMPLEC; stage 4)
//! and `bench g1` (the GPU kernel gate of stage 0).

use anyhow::Result;
use clap::{Parser, Subcommand};
use foehncraft_bench::G1Case;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "foehn",
    version,
    about = "Foehncraft – general purpose CFD solver"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Print the version, the enabled features and the platform.
    Info,
    /// Build or read a mesh, check it, partition it and write the quality
    /// report and a VTKHDF file for ParaView.
    Mesh(MeshArgs),
    /// Solve a steady scalar convection–diffusion case (stage 2 operators):
    /// one mesh, an order-of-accuracy sweep, or several solver variants, on
    /// the CPU and the GPU.
    Scalar(ScalarArgs),
    /// Study the pressure solver on a Poisson problem (stage 3): CG with the
    /// smoothed-aggregation multigrid, Jacobi and hypre, over refinement
    /// levels, on the CPU and the GPU.
    Poisson(PoissonArgs),
    /// Solve an incompressible laminar flow case: steady with SIMPLEC or
    /// SIMPLE (stage 4), transient with PISO and BDF2 or Crank–Nicolson
    /// when the case has a `[time]` table (stage 5); on one mesh, an
    /// order-of-accuracy sweep in space or a sweep over time steps.
    /// `--check` only checks the case file and the boundary conditions.
    /// Control files in the output directory steer a running case: `.stop`,
    /// `.result_N`, `.read_toml_N`.
    Flow(FlowArgs),
    /// Benchmarks.
    Bench {
        #[command(subcommand)]
        bench: Bench,
    },
}

#[derive(Subcommand)]
enum Bench {
    /// Gate G1: kernel benchmark, single-source route against hand-written
    /// references (docs/benchmarks/g1.md).
    G1(G1Args),
    /// List the backends compiled into this build.
    Backends,
}

#[derive(clap::Args)]
struct MeshArgs {
    /// Mesh case file (TOML) or a Gmsh `.msh` file.
    input: PathBuf,
    /// Number of parts (overrides `partition.parts`).
    #[arg(long)]
    parts: Option<usize>,
    /// Partitioning method: metis or rcb (overrides `partition.method`).
    #[arg(long)]
    method: Option<String>,
    /// Cell order: none, rcm or shuffle (overrides `mesh.renumber`).
    #[arg(long)]
    renumber: Option<String>,
    /// Random points of the point-location test (overrides
    /// `check.locate_points`).
    #[arg(long)]
    locate: Option<usize>,
    /// Do not write the VTKHDF file.
    #[arg(long)]
    no_vtkhdf: bool,
    /// Output directory (overrides `output.directory`).
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(clap::Args)]
struct ScalarArgs {
    /// Scalar case file (TOML).
    case: PathBuf,
    /// Comma-separated backends: cpu, cpu-serial, cuda (overrides
    /// `run.backends`).
    #[arg(long, value_delimiter = ',')]
    backends: Option<Vec<String>>,
    /// CPU threads (overrides `run.threads`).
    #[arg(long)]
    threads: Option<usize>,
    /// PTX module of the GPU backend (overrides `run.rust_ptx`).
    #[arg(long)]
    rust_ptx: Option<PathBuf>,
    /// Refinement levels of `[mesh.generate]`, comma-separated (overrides
    /// `verification.levels`).
    #[arg(long, value_delimiter = ',')]
    levels: Option<Vec<u32>>,
    /// Output directory (overrides `output.directory`).
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(clap::Args)]
struct PoissonArgs {
    /// Poisson case file (TOML).
    case: PathBuf,
    /// Comma-separated backends: cpu, cpu-serial, cuda (overrides
    /// `study.backends`).
    #[arg(long, value_delimiter = ',')]
    backends: Option<Vec<String>>,
    /// Comma-separated solvers: amg, amg-double, amg-plain, jacobi, hypre
    /// (overrides `study.solvers`).
    #[arg(long, value_delimiter = ',')]
    solvers: Option<Vec<String>>,
    /// CPU threads (overrides `study.threads`).
    #[arg(long)]
    threads: Option<usize>,
    /// PTX module of the GPU backend (overrides `study.rust_ptx`).
    #[arg(long)]
    rust_ptx: Option<PathBuf>,
    /// Refinement levels of `[mesh.generate]`, comma-separated (overrides
    /// `study.levels`).
    #[arg(long, value_delimiter = ',')]
    levels: Option<Vec<u32>>,
    /// Repetitions of each solve (overrides `study.repeat`).
    #[arg(long)]
    repeat: Option<usize>,
    /// Output directory (overrides `output.directory`).
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(clap::Args)]
struct FlowArgs {
    /// Flow case file (TOML).
    case: PathBuf,
    /// Check the case file and the boundary conditions, then stop.
    #[arg(long)]
    check: bool,
    /// CPU threads (overrides `run.threads`).
    #[arg(long)]
    threads: Option<usize>,
    /// Refinement levels of `[mesh.generate]`, comma-separated (overrides
    /// `verification.levels`).
    #[arg(long, value_delimiter = ',')]
    levels: Option<Vec<u32>>,
    /// Largest number of outer iterations (overrides
    /// `solver.max_iterations`).
    #[arg(long)]
    max_iterations: Option<usize>,
    /// Do not write the VTKHDF file.
    #[arg(long)]
    no_vtkhdf: bool,
    /// Continue from this restart file (overrides `run.restart`).
    #[arg(long)]
    restart: Option<PathBuf>,
    /// End time of a transient case in s (replaces `time.end_time` and
    /// `time.steps`).
    #[arg(long)]
    end_time: Option<f64>,
    /// Where the flow loop runs: cpu or cuda (overrides `run.backend`).
    #[arg(long)]
    backend: Option<String>,
    /// PTX module of the cuda backend (overrides `run.ptx`).
    #[arg(long)]
    ptx: Option<PathBuf>,
    /// Output directory (overrides `output.directory`).
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(clap::Args)]
struct G1Args {
    /// Case file (TOML); every key is optional.
    case: Option<PathBuf>,
    /// Approximate number of cells (overrides `mesh.cells`).
    #[arg(long)]
    cells: Option<usize>,
    /// Cell numbering: natural or shuffled (overrides `mesh.ordering`).
    #[arg(long)]
    ordering: Option<String>,
    /// Comma-separated backends (overrides `run.backends`).
    #[arg(long, value_delimiter = ',')]
    backends: Option<Vec<String>>,
    /// Timed repetitions (overrides `run.repeat`).
    #[arg(long)]
    repeat: Option<usize>,
    /// CPU threads (overrides `run.threads`).
    #[arg(long)]
    threads: Option<usize>,
    /// PTX module of the single-source GPU route (overrides `run.rust_ptx`).
    #[arg(long)]
    rust_ptx: Option<PathBuf>,
    /// Output directory (overrides `output.directory`).
    #[arg(long)]
    out: Option<PathBuf>,
}

fn main() -> Result<()> {
    // GPU backends (feature `cuda`; nothing without it)
    foehncraft_cuda::register();
    match Cli::parse().command {
        Command::Info => info(),
        Command::Mesh(a) => mesh(a),
        Command::Scalar(a) => scalar(a),
        Command::Poisson(a) => poisson(a),
        Command::Flow(a) => flow(a),
        Command::Bench {
            bench: Bench::Backends,
        } => {
            for b in foehncraft_bench::backend::available() {
                println!("{b}");
            }
            Ok(())
        }
        Command::Bench {
            bench: Bench::G1(a),
        } => g1(a),
    }
}

fn info() -> Result<()> {
    println!("Foehncraft {}", foehncraft_core::VERSION);
    println!("platform     {}", foehncraft_bench::sysinfo::platform());
    println!("cpu          {}", foehncraft_bench::sysinfo::cpu_model());
    println!("threads      {}", rayon::current_num_threads());
    println!("precision    f64 (reference path, D7)");
    {
        use foehncraft_parallel::Communicator as _;
        let c = foehncraft_parallel::SerialComm;
        println!(
            "ranks        {} ({} communicator; MPI runs start with mpirun)",
            c.size(),
            c.name()
        );
    }
    let features: Vec<&str> = [
        ("cuda", cfg!(feature = "cuda")),
        ("mpi", cfg!(feature = "mpi")),
        ("checked-indexing", cfg!(feature = "checked-indexing")),
        ("hdf5", cfg!(feature = "hdf5")),
        ("metis", cfg!(feature = "metis")),
        ("hypre", cfg!(feature = "hypre")),
    ]
    .into_iter()
    .filter_map(|(n, on)| on.then_some(n))
    .collect();
    println!(
        "features     {}",
        if features.is_empty() {
            "none".to_string()
        } else {
            features.join(", ")
        }
    );
    println!(
        "backends     {}",
        foehncraft_bench::backend::available().join(", ")
    );
    Ok(())
}

fn mesh(a: MeshArgs) -> Result<()> {
    use foehncraft_mesh::case::MeshCase;
    let is_msh = a.input.extension().is_some_and(|e| e == "msh");
    let mut case = if is_msh {
        MeshCase::for_file(&a.input)
    } else {
        MeshCase::from_file(&a.input)?
    };
    let name = a
        .input
        .file_stem()
        .map_or("mesh".into(), |s| s.to_string_lossy().into_owned());
    if let Some(p) = a.parts {
        case.partition.parts = p;
    }
    if let Some(m) = a.method {
        case.partition.method = match m.as_str() {
            "metis" => foehncraft_mesh::partition::Method::Metis,
            "rcb" => foehncraft_mesh::partition::Method::Rcb,
            _ => anyhow::bail!("--method: metis or rcb"),
        };
    }
    if let Some(r) = a.renumber {
        use foehncraft_mesh::renumber::Ordering;
        case.mesh.renumber = match r.as_str() {
            "none" => Ordering::None,
            "rcm" => Ordering::Rcm,
            "shuffle" => Ordering::Shuffle,
            _ => anyhow::bail!("--renumber: none, rcm or shuffle"),
        };
    }
    if let Some(n) = a.locate {
        case.check.locate_points = n;
    }
    if a.no_vtkhdf {
        case.output.vtkhdf = false;
    }
    if a.out.is_some() {
        case.output.directory = a.out;
    }
    let (report, files) = foehncraft_mesh::case::run(&case, &name, &mut |l| println!("{l}"))?;
    for f in files {
        println!("wrote {}", f.display());
    }
    anyhow::ensure!(report.quality.ok(), "the mesh is invalid (see the report)");
    Ok(())
}

fn scalar(a: ScalarArgs) -> Result<()> {
    use foehncraft_fv::case::{ScalarCase, check, output_dir, run};
    let mut case = ScalarCase::from_file(&a.case)?;
    if let Some(b) = a.backends {
        case.run.backends = b;
    }
    if a.threads.is_some() {
        case.run.threads = a.threads;
    }
    if let Some(p) = a.rust_ptx {
        case.run.rust_ptx = p;
    }
    if let Some(l) = a.levels {
        case.verification.levels = l;
    }
    let (name, dir) = output_dir(&case, &a.case, a.out);
    let rep = run(&case, &name, &dir, |l| println!("{l}"))?;
    if !rep.orders.is_empty() {
        println!(
            "observed order: {}",
            rep.orders
                .iter()
                .map(|o| format!("{o:.3}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    println!("wrote {}", dir.join("scalar.md").display());
    check(&rep)
}

fn flow(a: FlowArgs) -> Result<()> {
    use foehncraft_flow::case::{FlowCase, check, check_case, output_dir, run};
    let mut case = FlowCase::from_file(&a.case)?;
    if a.threads.is_some() {
        case.run.threads = a.threads;
    }
    if let Some(l) = a.levels {
        case.verification.levels = l;
    }
    if let Some(n) = a.max_iterations {
        case.solver.max_iterations = n;
    }
    if a.no_vtkhdf {
        case.output.vtkhdf = false;
    }
    if let Some(b) = a.backend {
        case.run.backend = b;
    }
    if let Some(p) = a.ptx {
        case.run.ptx = p;
    }
    if a.restart.is_some() {
        case.run.restart = a.restart;
    }
    if let Some(e) = a.end_time {
        let Some(t) = case.time.as_mut() else {
            anyhow::bail!("--end-time needs a transient case (a [time] table)");
        };
        t.end_time = Some(e);
        t.steps = None;
        t.validate()?;
    }
    if a.check {
        return check_case(&case, &a.case, |l| println!("{l}"));
    }
    let (name, dir) = output_dir(&case, &a.case, a.out);
    let rep = run(&case, &name, &dir, |l| println!("{l}"))?;
    println!("wrote {}", dir.join("flow.md").display());
    check(&rep)
}

fn poisson(a: PoissonArgs) -> Result<()> {
    use foehncraft_fv::poisson::{PoissonCase, StudySolver, check, output_dir, run};
    let mut case = PoissonCase::from_file(&a.case)?;
    if let Some(b) = a.backends {
        case.study.backends = b;
    }
    if let Some(s) = a.solvers {
        case.study.solvers = s
            .iter()
            .map(|n| n.parse::<StudySolver>())
            .collect::<Result<_>>()?;
    }
    if a.threads.is_some() {
        case.study.threads = a.threads;
    }
    if let Some(p) = a.rust_ptx {
        case.study.rust_ptx = p;
    }
    if let Some(l) = a.levels {
        case.study.levels = l;
    }
    if let Some(r) = a.repeat {
        case.study.repeat = r;
    }
    case.validate()?;
    let (name, dir) = output_dir(&case, &a.case, a.out);
    let rep = run(&case, &name, &dir, |l| println!("{l}"))?;
    for (n, r) in &rep.iteration_ratios {
        println!("iterations finest / coarsest level, {n}: {r:.2}");
    }
    println!("wrote {}", dir.join("poisson.md").display());
    check(&rep)
}

fn g1(a: G1Args) -> Result<()> {
    let mut case = match &a.case {
        Some(p) => G1Case::from_file(p)?,
        None => G1Case::default(),
    };
    if let Some(n) = a.cells {
        case.mesh.cells = n;
        case.mesh.nx = None;
    }
    if let Some(o) = a.ordering {
        case.mesh.ordering = match o.as_str() {
            "natural" => foehncraft_bench::mesh::Ordering::Natural,
            "shuffled" => foehncraft_bench::mesh::Ordering::Shuffled,
            _ => anyhow::bail!("--ordering: natural or shuffled"),
        };
    }
    if let Some(b) = a.backends {
        case.run.backends = b;
    }
    if let Some(r) = a.repeat {
        case.run.repeat = r;
    }
    if a.threads.is_some() {
        case.run.threads = a.threads;
    }
    if let Some(p) = a.rust_ptx {
        case.run.rust_ptx = p;
    }
    if let Some(o) = a.out {
        case.output.directory = o;
    }
    let report = foehncraft_bench::run_g1(&case, &mut |line| println!("{line}"))?;
    for v in &report.verdicts {
        println!(
            "G1 {} against {}: {}",
            v.candidate,
            v.reference,
            if v.pass { "PASS" } else { "FAIL" }
        );
    }
    let files = report.write(&case.output.directory, &case.output.name)?;
    for f in files {
        println!("wrote {}", f.display());
    }
    let all_correct = report
        .results
        .iter()
        .all(|b| b.kernels.iter().all(|k| k.within_tolerance));
    anyhow::ensure!(
        all_correct,
        "a backend did not reproduce the reference results"
    );
    Ok(())
}

//! `cfdcraft`: the command line of CFDCraft (D4).
//!
//! Commands: `info` (build and platform), `mesh` (build or read a mesh,
//! check, partition and write it; stage 1), `scalar` (steady scalar
//! transport, verification sweeps, CPU and GPU; stage 2) and `bench g1` (the
//! GPU kernel gate of stage 0). The flow solver commands follow with the
//! later stages.

use anyhow::Result;
use cfdcraft_bench::G1Case;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "cfdcraft",
    version,
    about = "CFDCraft – general purpose CFD solver"
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
    match Cli::parse().command {
        Command::Info => info(),
        Command::Mesh(a) => mesh(a),
        Command::Scalar(a) => scalar(a),
        Command::Bench {
            bench: Bench::Backends,
        } => {
            for b in cfdcraft_bench::backend::available() {
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
    println!("CFDCraft {}", cfdcraft_core::VERSION);
    println!("platform     {}", cfdcraft_bench::sysinfo::platform());
    println!("cpu          {}", cfdcraft_bench::sysinfo::cpu_model());
    println!("threads      {}", rayon::current_num_threads());
    println!("precision    f64 (reference path, D7)");
    {
        use cfdcraft_parallel::Communicator as _;
        let c = cfdcraft_parallel::SerialComm;
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
        cfdcraft_bench::backend::available().join(", ")
    );
    Ok(())
}

fn mesh(a: MeshArgs) -> Result<()> {
    use cfdcraft_mesh::case::MeshCase;
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
            "metis" => cfdcraft_mesh::partition::Method::Metis,
            "rcb" => cfdcraft_mesh::partition::Method::Rcb,
            _ => anyhow::bail!("--method: metis or rcb"),
        };
    }
    if let Some(r) = a.renumber {
        use cfdcraft_mesh::renumber::Ordering;
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
    let (report, files) = cfdcraft_mesh::case::run(&case, &name, &mut |l| println!("{l}"))?;
    for f in files {
        println!("wrote {}", f.display());
    }
    anyhow::ensure!(report.quality.ok(), "the mesh is invalid (see the report)");
    Ok(())
}

fn scalar(a: ScalarArgs) -> Result<()> {
    use cfdcraft_fv::case::{ScalarCase, check, output_dir, run};
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
            "natural" => cfdcraft_bench::mesh::Ordering::Natural,
            "shuffled" => cfdcraft_bench::mesh::Ordering::Shuffled,
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
    let report = cfdcraft_bench::run_g1(&case, &mut |line| println!("{line}"))?;
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

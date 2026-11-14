//! `cfdcraft`: the command line of CFDCraft (D4).
//!
//! Stage 0 has two commands: `info` (build and platform) and `bench g1`
//! (the GPU kernel gate). The solver commands follow with the later stages.

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

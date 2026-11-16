<p align="center">
  <img src="assets/logo/foehncraft.svg" alt="Foehncraft logo: a peregrine falcon" width="180">
</p>

<h1 align="center">Foehncraft</h1>

<p align="center">
  An open-source, industry-strength, general purpose CFD solver, written in Rust and built for modern GPUs.
</p>

---

> **Status: stage 4 done (laminar steady flow, release 0.4).** The
> requirements are in [docs/requirements.md](docs/requirements.md). The code
> so far is the workspace, continuous integration, the GPU kernel benchmark
> that decided the GPU backend ([G1](docs/benchmarks/g1.md)), the
> polyhedral mesh (Gmsh reader, test mesher, quality report, partitioning,
> VTKHDF output), the finite-volume operators on a scalar transport
> equation on the CPU and the GPU ([stage 2](docs/verification/stage2.md)),
> the pressure solver: aggregation AMG with CG in mixed precision, with
> hypre as the reference ([stage 3](docs/verification/stage3.md)), and the
> first flow solver: steady incompressible laminar flow with SIMPLEC, the
> Rhie–Chow interpolation and checked boundary conditions, verified on
> Couette, Poiseuille and the lid-driven cavity
> ([stage 4](docs/verification/stage4.md)). Stage 5 is in progress: the
> transient core is in (PISO with BDF2 and Crank–Nicolson, adaptive time
> step, restart, control files, verified on the Taylor–Green vortex and the
> cylinder benchmark of Schäfer & Turek, [stage 5](docs/verification/stage5.md)),
> the whole flow loop runs on NVIDIA GPUs with results bitwise equal to
> the CPU, the kernels are ready for a discrete adjoint (dual-number
> derivatives, transposed matrices, binomial checkpointing), every run
> books its mass and kinetic-energy balance, and the numerics gained a
> bounded central scheme, a time step chosen by the time-integration
> error, output at exact times and a potential-flow start.

[![CI](https://github.com/chiefenne/Foehncraft/actions/workflows/ci.yml/badge.svg)](https://github.com/chiefenne/Foehncraft/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/chiefenne/Foehncraft)](https://github.com/chiefenne/Foehncraft/releases)
[![Rust 1.89+](https://img.shields.io/badge/rust-1.89%2B-orange.svg)](Cargo.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Goal

Foehncraft aims to be a general purpose computational fluid dynamics solver of
industry strength: the class of tool that engineers rely on every day for
fluid flow, heat and mass transfer, species transport and reacting flows, on real
geometries and at real problem sizes. It is free and open source under the MIT
licence.

The guiding targets are:

- **General purpose.** Arbitrary polyhedral meshes in 3D, with a physics
  roadmap that reaches from incompressible flow to conjugate heat transfer,
  reacting flow, multiphase flow, moving meshes and scale-resolving
  turbulence.
- **Performance first.** Modern parallel GPUs are the primary platform, with
  a design target of more than 100 million cells on several GPUs and on
  clusters. The same solver runs on any CPU.
- **Accuracy that can be checked.** Double precision, numerics taken from the
  published literature, and every development stage verified against analytic
  solutions and published reference data.
- **Solver first.** A case file and a command line; meshes from established
  open tools, results to established open post-processors.

## First milestone

Incompressible flow on unstructured meshes, with two or three RANS turbulence
models, run from a text case file. Everything after that follows the roadmap
in the [requirements](docs/requirements.md).

## Build and run

Needs Rust 1.89 or newer, a C compiler and CMake (HDF5 and METIS are built
from source on the first build; recipe 1 has the details).

```
cargo test --workspace
cargo build --release -p foehncraft-cli
target/release/foehn info
target/release/foehn bench g1 examples/cases/bench/g1_ci_small.toml
target/release/foehn mesh examples/cases/mesh/poly_partition.toml
target/release/foehn scalar examples/cases/scalar/mms_poly_box.toml
target/release/foehn poisson examples/cases/poisson/poisson_poly_box.toml
target/release/foehn flow examples/cases/flow/cavity_re100.toml
target/release/foehn flow examples/cases/flow/taylor_green_bdf2.toml
```

## Documentation

| Kind | What it is for | Where |
|---|---|---|
| Tutorial | a first mesh, end to end, in about ten minutes | [docs/tutorial.md](docs/tutorial.md) |
| How-to | steps for one task: build, run, check, partition, rerun G1, order of accuracy, pressure solver and hypre, GPU runs, flow cases and boundary conditions, transient runs, restart and control files, scripts and agents, balances and divergence reports | [docs/recipes.md](docs/recipes.md) |
| Reference | every case file, its keys and defaults; the machine interfaces (JSON formats, exit codes, signals, diagnostic codes, balance and divergence files) | [docs/cases.md](docs/cases.md), [docs/keys.md](docs/keys.md), [docs/interfaces.md](docs/interfaces.md) |
| Examples | what each feature computes, with pictures and the commands | [docs/gallery.md](docs/gallery.md) |
| Explanation | requirements, decisions, benchmark and verification results and their reading | [docs/requirements.md](docs/requirements.md), [G1](docs/benchmarks/g1.md), [stage 2](docs/verification/stage2.md), [stage 3](docs/verification/stage3.md), [stage 4](docs/verification/stage4.md), [stage 5](docs/verification/stage5.md) |

The [changelog](CHANGELOG.md) lists every release.

| Crate | Content |
|---|---|
| `foehncraft-core` | generic floating-point type, index types, reproducible sums, counter-based random numbers (`no_std`) |
| `foehncraft-kernels` | compute kernels, one source for CPU and GPU (`no_std`): benchmark, scalar operators, multigrid, flow |
| `foehncraft-parallel` | communicator (serial, MPI), reproducible threaded reductions |
| `foehncraft-mesh` | polyhedral mesh, geometry, quality report, Gmsh reader, test mesher (also the cylinder benchmark mesh), partitioning and halo, VTKHDF output, point location |
| `foehncraft-fv` | finite-volume operators on a scalar (gradients, limiters, convection, diffusion, assembly), BiCGStab, and the pressure solver: CG with smoothed-aggregation multigrid in mixed precision, optional hypre reference; CPU backend (the GPU backend is in `foehncraft-cuda`) |
| `foehncraft-flow` | incompressible flow: SIMPLEC and SIMPLE (steady), PISO with BDF2 and Crank–Nicolson (transient), the Rhie–Chow interpolation, boundary conditions and their check, restart and control files, probes and force coefficients, the flow device written once over a CPU or GPU executor, the `foehn flow` case, exact solutions, reference and benchmark data |
| `foehncraft-bench` | gate G1: kernel benchmark, CPU backends |
| `foehncraft-cuda` | all NVIDIA GPU code (D60): the flow executor, the scalar and multigrid backend, the G1 backends with the CUDA C reference |
| `foehncraft-cli` | the `foehn` command |
| `gpu/ptx-kernels` | the kernels compiled to NVIDIA PTX (nightly toolchain) |

## How it is built

Foehncraft follows the paradigm of the
[ArtCraft Crafting Apps](https://getartcraft.com/apps/craft-anything): open
source software in Rust, one engine with thin front ends on top, and
development carried out together with a coding agent in a loop of exploring,
building, checking and refining. The falcon joins the animals of that family.

What the agent builds is decided by a person. The development is driven by the
expert knowledge of the maintainer, from 40 years in computational fluid
dynamics: the requirements, the choice of methods, the test cases and the
judgement of every result come from that experience. The agent writes and
checks code; it does not decide what a CFD solver has to be.

Three rules follow from this:

1. **Requirements before code.** Each decision is recorded with its reason in
   [docs/requirements.md](docs/requirements.md) before anything is implemented.
2. **Every stage ends with results that can be tested.** No stage is done
   until its verification cases pass.
3. **Clean-room numerics.** Methods come from the published literature, with
   references; no code is taken from other solvers.

## Why a falcon

The peregrine falcon is the fastest animal there is, and its shape is the work
of the flow around it.

## Why the name

**Foehn** (German *Föhn*) is pronounced [føːn], with a long "ö" as in French *feu*, followed by "n". English speakers say /fɜːn/ ("fern") or /feɪn/. Foehncraft is therefore [ˈføːnkraft].

The foehn is the warm, dry wind that falls down the lee side of a mountain
range. Innsbruck, where the maintainer grew up and where the name comes from,
is one of the best-known foehn places in the world. With a southerly flow, air crosses the main Alpine ridge
at the Brenner Pass, the lowest gap at about 1370 m, runs down the Wipp valley
past Mount Patscherkofel and breaks into the Inn valley at Innsbruck, often in
strong gusts. Sinking without cloud, the air warms by about 1 °C per 100 m of
descent; when the foehn reaches the valley floor at about 580 m, the
temperature in town can rise by 5 to 10 °C within an hour, and the cold air
pooled in the valley is swept away.

Gap flow, a downslope windstorm, stratified flow over complex terrain and a
moving front between warm and cold air: the foehn is computational fluid
dynamics in the open air, and it has been studied at the University of
Innsbruck since the early twentieth century. The command is `foehn`.

## Contributing

The project is at stage 4 and has one maintainer. Comments on the
requirements are welcome as GitHub issues. Contribution guidelines will
follow with the first solver code.

## Licence

Distributed under the MIT licence. See [LICENSE](LICENSE) for more information.

<!-- year:start -->2026<!-- year:end --> Andreas Ennemoser – andreas.ennemoser@aon.at

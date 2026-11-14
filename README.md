<p align="center">
  <img src="assets/logo/cfdcraft.svg" alt="CFDCraft logo: a peregrine falcon" width="180">
</p>

<h1 align="center">CFDCraft</h1>

<p align="center">
  An open-source, industry-strength, general purpose CFD solver, written in Rust and built for modern GPUs.
</p>

---

> **Status: stage 2 done (operators on a scalar equation, release 0.2).** The
> requirements are in [docs/requirements.md](docs/requirements.md). The code
> so far is the workspace, continuous integration, the GPU kernel benchmark
> that decided the GPU backend ([G1](docs/benchmarks/g1.md)), the
> polyhedral mesh (Gmsh reader, test mesher, quality report, partitioning,
> VTKHDF output) and the finite-volume operators on a scalar transport
> equation, on the CPU and the GPU: least-squares gradients, limiters,
> convection and diffusion, reproducible assembly and BiCGStab
> ([stage 2](docs/verification/stage2.md)). The flow solver starts with
> stage 3 (pressure solver).

[![CI](https://github.com/chiefenne/CFDCraft/actions/workflows/ci.yml/badge.svg)](https://github.com/chiefenne/CFDCraft/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/chiefenne/CFDCraft)](https://github.com/chiefenne/CFDCraft/releases)
[![Rust 1.89+](https://img.shields.io/badge/rust-1.89%2B-orange.svg)](Cargo.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Goal

CFDCraft aims to be a general purpose computational fluid dynamics solver of
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
cargo build --release -p cfdcraft-cli
target/release/cfdcraft info
target/release/cfdcraft bench g1 examples/cases/bench/g1_ci_small.toml
target/release/cfdcraft mesh examples/cases/mesh/poly_partition.toml
target/release/cfdcraft scalar examples/cases/scalar/mms_poly_box.toml
```

## Documentation

| Kind | What it is for | Where |
|---|---|---|
| Tutorial | a first mesh, end to end, in about ten minutes | [docs/tutorial.md](docs/tutorial.md) |
| How-to | steps for one task: build, run, check, partition, rerun G1, order of accuracy, GPU runs | [docs/recipes.md](docs/recipes.md) |
| Reference | every case file, its keys and defaults | [docs/cases.md](docs/cases.md) |
| Examples | what each feature computes, with pictures and the commands | [docs/gallery.md](docs/gallery.md) |
| Explanation | requirements, decisions, benchmark and verification results and their reading | [docs/requirements.md](docs/requirements.md), [G1](docs/benchmarks/g1.md), [stage 2](docs/verification/stage2.md) |

The [changelog](CHANGELOG.md) lists every release.

| Crate | Content |
|---|---|
| `cfdcraft-core` | generic floating-point type, index types, reproducible sums, counter-based random numbers (`no_std`) |
| `cfdcraft-kernels` | compute kernels, one source for CPU and GPU (`no_std`) |
| `cfdcraft-parallel` | communicator (serial, MPI), reproducible threaded reductions |
| `cfdcraft-mesh` | polyhedral mesh, geometry, quality report, Gmsh reader, test mesher, partitioning and halo, VTKHDF output, point location |
| `cfdcraft-bench` | gate G1: kernel benchmark, CPU and CUDA backends |
| `cfdcraft-cli` | the `cfdcraft` command |
| `gpu/ptx-kernels` | the kernels compiled to NVIDIA PTX (nightly toolchain) |

## How it is built

CFDCraft follows the paradigm of the
[ArtCraft Crafting Apps](https://getartcraft.com/apps/craft-anything): open
source software in Rust, one engine with thin front ends on top, and
development carried out together with a coding agent in a loop of exploring,
building, checking and refining. The falcon joins the animals of that family.

What the agent builds is decided by a person. The development is driven by the
expert knowledge of the maintainer, Andreas Ennemoser, from 40 years in
computational fluid dynamics: the requirements, the choice of methods, the
test cases and the judgement of every result come from that experience. The
agent writes and checks code; it does not decide what a CFD solver has to be.

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

## Contributing

The project is at stage 1 and has one maintainer. Comments on the
requirements are welcome as GitHub issues. Contribution guidelines will
follow with the first solver code.

## Licence

Distributed under the MIT licence. See [LICENSE](LICENSE) for more information.

<!-- year:start -->2026<!-- year:end --> Andreas Ennemoser – andreas.ennemoser@aon.at

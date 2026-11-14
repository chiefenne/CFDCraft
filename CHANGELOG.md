# Changelog

## Summary

| Version | Stage | Highlights |
|---|---|---|
| Unreleased | 0 Foundations and gate G1 | Cargo workspace, CI on Linux, macOS and Windows with MPI runs on 4 ranks, the G1 kernel benchmark (CPU measured; GPU run pending the maintainer's GPU) |

The tag `v0.0` is set when both verdicts of gate G1 are in: the CPU verdict
from the maintainer's machine and the GPU verdict from a free Google Colab
T4 (O13).

---

## Unreleased – stage 0: foundations and gate G1

Workspace (`Cargo.toml`, Rust 2024 edition, minimum Rust 1.89):

- `cfdcraft-core` (`no_std`): the generic floating-point type `Real`
  (f64 reference path, f32 possible, D7), 32-bit local and 64-bit global
  indices, pairwise summation in a fixed order (D31), counter-based random
  numbers keyed by identity, step and event (D47).
- `cfdcraft-kernels` (`no_std`): the kernel bodies, one source for CPU and
  GPU: face flux, cell gather, SpMV in CSR and SELL-C-σ. Gather-based, no
  atomics. Feature `checked-indexing` for bounds-checked loads.
- `cfdcraft-parallel`: `Communicator` trait with `SerialComm` and `MpiComm`
  (feature `mpi`, system MPI library, D9); reproducible all-rank sums;
  multi-threaded sums bitwise equal to the sequential sum for any number of
  threads.
- `cfdcraft-bench`: gate G1. Synthetic periodic face-based mesh (natural or
  shuffled numbering), CSR and SELL-C-σ matrices, backends `cpu-serial`,
  `cpu-threads`, `cpu-native`, and with feature `cuda` `cuda-reference`
  (hand-written CUDA C, compiled by NVRTC at run time) and `cuda-rust` (the
  kernel bodies compiled to PTX). TOML case files, JSON, Markdown and SVG
  reports, G1 verdict per candidate/reference pair.
- `cfdcraft-cli`: the `cfdcraft` command with `info`, `bench g1` and
  `bench backends`.
- `gpu/ptx-kernels`: the kernel bodies as PTX entry points, built with the
  nightly toolchain (`build.sh`).

Continuous integration (`.github/workflows/ci.yml`, D22): format, clippy and
docs; tests in debug and release on Linux, macOS and Windows; minimum Rust
version; MPI smoke test with 2 and 4 ranks on Linux and macOS (O12); PTX
build of the single-source kernels; G1 smoke run on every platform.

Cloud GPU runs: `tools/g1_cloud.sh` (installs Rust, builds PTX and CLI for
the GPU found, runs G1 twice) and `tools/g1_colab.ipynb` (the same on a free
Google Colab T4).

Documentation: [gallery](docs/gallery.md), [case catalogue](docs/cases.md),
[recipes](docs/recipes.md), [G1 procedure and results](docs/benchmarks/g1.md).

Findings:

- The rustc NVPTX target supports compute capability 7.0 (Volta) and newer
  only; G1 needs such a GPU for the single-source route (O13). The
  maintainer's GTX 960 (compute capability 5.2) runs the CPU verdict and the
  CUDA reference (`g1_gpu_reference_4m.toml`, recipe 2); the GPU verdict is
  taken on a free Colab T4 (recipe 3).
- On the CPU the single-source bodies match hand-written loops once the
  parallel loop uses large enough work packets (64 Ki outputs): with 4 Ki the
  face-flux kernel was about 30 % slower on 10 million cells.
- A branching max/min in the upwind flux cost bandwidth on random mass-flux
  signs; `Real::max_of`/`min_of` are now branch-free.

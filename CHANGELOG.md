# Changelog

## Summary

| Version | Stage | Highlights |
|---|---|---|
| 0.3 (2026-11-14) | 3 Pressure solver | Pressure linear solver: smoothed-aggregation algebraic multigrid as the preconditioner of CG (Vaněk, Mandel, Brezina 1996), symmetric V-cycle with Chebyshev or Jacobi smoothing, float32 hierarchy with float64 Krylov iteration (mixed precision, D37), one kernel source for the CPU and the NVIDIA GPU, optional hypre 3.2.1 reference (BoomerAMG-PCG, built from source). **10–13 iterations from 512 to 274 625 cells on all mesh types; same solution as hypre to 1.6e-9; on a Tesla T4 6.5 ns per cell and iteration against 10.5 for hypre at 912 673 cells, all GPU results bitwise equal to the CPU.** New command `cfdcraft poisson`, four Poisson case files, gallery 10–11, recipes 15–17 |
| 0.2 (2026-11-14) | 2 Operators on a scalar equation | Finite-volume operators on a scalar transport equation, one source for the CPU and the NVIDIA GPU: least-squares gradients (exact for linear fields, error 6e-15 on all mesh types), linear upwind with Venkatakrishnan and Venkatakrishnan–Wang limiters (no new extrema at K = 0), over-relaxed non-orthogonal diffusion, reproducible gather-based assembly into sliced ELL, BiCGStab. Observed order 1.87–2.08 on hex, tet, poly and mixed meshes in 2D and 3D. **GPU and CPU results bitwise equal on a Tesla T4** (912 673 polyhedra, 311 ns per cell and iteration). New command `cfdcraft scalar`, seven case files, gallery 8–9, recipes 12–14 |
| 0.1 (2026-11-14) | 1 Mesh | Polyhedral face-based mesh with exact geometry (every cell closed to 1e-15, volumes exact), Gmsh 4.1 reader, built-in test mesher (hex, tet, poly, mixed; 2D and 3D; graded, refined, perturbed), quality report, VTKHDF output for ParaView, METIS and RCB partitioning with a solver-ready halo layout and multi-field halo exchange (tested over MPI on 2–4 ranks), particle point location on all mesh types, golden-value test matrix, `cfdcraft mesh` command. G1 rerun on polyhedral and shuffled meshes: bitwise equal; quiet-machine verdict open (O15) |
| 0.0 (2026-11-14) | 0 Foundations and gate G1 | Cargo workspace, CI on Linux, macOS and Windows with MPI runs on 4 ranks, the G1 kernel benchmark. **G1 GPU verdict passed on a Tesla T4: single-source Rust kernels within 1.0 % of hand-written CUDA C, bitwise equal results**; **CPU verdict passed on an i9-9900K: single-source bodies within 1.8 % of hand-written loops, in two runs** |

Release 0.3 closes stage 3: all acceptance tests of the stage plan pass
(iteration count nearly independent of the mesh size, same converged result
as hypre, GPU time per cell compared with hypre and reported). The Cargo
version of the crates is now 0.3.0; the release is the Git tag `v0.3` (D46).

Release 0.2 closes stage 2: all acceptance tests of the stage plan pass
(exact gradients, no new extrema, observed order 1.8–2.2 on all mesh types,
bitwise reproducible results, GPU and CPU bitwise equal on a Tesla T4). The
Cargo version of the crates is now 0.2.0; the release is the Git tag `v0.2`
(D46).

Release 0.1 closes stage 1: all acceptance tests of the stage plan pass
(closed cells, exact volumes, Gmsh files from Gmsh load, VTKHDF files load
with VTK's reader, point location on all mesh types). The Cargo version of
the crates is now 0.1.0; the release is the Git tag `v0.1` (D46).

Release 0.0 closes stage 0: both verdicts of gate G1 are in (GPU: Tesla T4,
CPU: i9-9900K, 2026-11-14, both PASS; O13 closed). The Cargo version of the
crates stays 0.0.0; the release is the Git tag `v0.0` (D46).

---

## 0.3 (2026-11-14) – stage 3: pressure solver

New module `cfdcraft-fv::amg`, kernel bodies `cfdcraft-kernels::amg` (one
source for the CPU and the GPU, float32 and float64) and the command
`cfdcraft poisson`:

- **Smoothed-aggregation multigrid (D35).** Strength of connection,
  three-phase greedy aggregation in the cell order, filtered and
  Jacobi-smoothed prolongator, Galerkin coarse operators, explicit
  restriction (gather, no atomics), dense coarsest level; deterministic host
  setup in float64. Plain aggregation as an option.
- **Symmetric V-cycle** with a Chebyshev (default) or damped Jacobi smoother
  in D⁻¹A, λ_max from a fixed power iteration.
- **Mixed precision (D37).** The hierarchy (including a float32 copy of the
  fine matrix) in float32, CG in float64; `precision = "double"` for
  comparison. R10 allows f32 only in the multigrid kernels and rejects
  contracted, approximate and flush-to-zero f32 arithmetic.
- **Flexible CG** (Notay 2000). `[solver] linear_solver = "cg"` and
  `preconditioner = "amg"` in `cfdcraft scalar` (also AMG with BiCGStab);
  the report shows the linear iterations and the multigrid setup time.
- **GPU.** PTX entries `amg_*_f32`, `amg_*_f64`, `blas_axpy`, `blas_xpby`;
  the CUDA backend runs the V-cycle with all vectors resident.
- **hypre reference (D20).** Cargo feature `hypre`: a C shim over hypre's
  public interface (BoomerAMG with PCG), `tools/build_hypre.sh` builds hypre
  from source (CPU or CUDA); a CI job builds it and compares.
- **`cfdcraft poisson`.** Mesh sweeps with several solvers and backends;
  report of iterations, setup and solve times, ns per cell and iteration,
  hierarchy, differences to hypre, bitwise comparison across backends;
  checks of iteration growth and agreement.
- **Verification.** Stage-3 tests (iterations nearly independent of the
  mesh size on all mesh types, mixed precision equals float64,
  SPD V-cycle, bitwise reproducibility, agreement with hypre), five Poisson
  case files and one scalar case, gallery 10–11, recipes 15–17, the
  explanation page `docs/verification/stage3.md`, `tools/stage3_colab.ipynb`.

## 0.2 (2026-11-14) – stage 2: operators on a scalar equation

New crate `cfdcraft-fv` (features `hdf5`, default on, and `cuda`) and new
kernel bodies in `cfdcraft-kernels` (`fv`, `blas`), one source for the CPU
and the GPU:

- **Least-squares gradient (D24).** Weighted least squares over the face
  neighbours; weights per cell-face entry precomputed on the host;
  fixed-value, fixed-gradient and symmetry boundaries keep linear fields
  exact on every mesh type.
- **Convection (D25).** First-order upwind and linear upwind as a deferred
  correction; Venkatakrishnan and Venkatakrishnan–Wang limiters.
- **Diffusion (D26).** Over-relaxed non-orthogonal correction, optionally
  limited.
- **Assembly (D31, D37).** Face kernel for the flux coefficients, cell
  kernel that gathers them into a sliced-ELL row (C = 32, σ = 1), the
  right-hand side and the Jacobi diagonal; implicit under-relaxation. No
  atomics; results independent of the thread count.
- **BiCGStab (D36)** with Jacobi preconditioning; dot products in fixed
  blocks of 256, the same order on the CPU and the GPU.
- **Backends.** CPU (rayon) and NVIDIA GPU (`cuda`: the bodies as PTX entry
  points `fv_*`, `blas_*` in `gpu/ptx-kernels`, all fields resident on the
  GPU).
- **`cfdcraft scalar`.** Case file with `[problem]` (manufactured solutions,
  boundary profiles), `[solver]`, `[run]` (backends compared with each
  other), `[verification]` (refinement sweeps, expected order, bounds),
  `[[variants]]` and `[output]` (VTKHDF fields, sampling line). Reports in
  Markdown and JSON with time per phase and ns per cell and iteration (R6),
  convergence and profile plots as SVG.
- **Verification.** Tests for exact gradients, no new extrema, order of
  accuracy 1.8–2.2 on hex, tet, poly and mixed meshes in 2D and 3D
  (release builds, R7), bitwise reproducibility; seven scalar case files,
  gallery examples 8 and 9, recipes 12–14, the explanation page
  `docs/verification/stage2.md`.
- **CI.** Every scalar case except the GPU one runs as a test (R5);
  `tools/check_kernels.sh` refuses atomics, f32 and fused or approximate f64
  arithmetic in the kernels and the PTX (R10).
- **GPU check.** `tools/stage2_colab.ipynb` and `tools/stage2_cloud.sh`
  (CPU against two GPU runs on 0.9 million polyhedra). Result on a Tesla
  T4 (Colab): both GPU runs, the four oblique-step schemes and the sweep
  are bitwise equal to the CPU; 7.1 s against 158.8 s on the 2-thread
  Colab VM.

---

## 0.1 (2026-11-14) – stage 1: mesh

- Release workflow (`.github/workflows/release.yml`): creates the tag
  `v<version>` on `main` and a GitHub release with the changelog section of
  that version as notes; dry-run by default (recipe 6).

New crate `cfdcraft-mesh` (features `hdf5` and `metis`, both on by default;
`mpi` for the MPI halo test):

- **Mesh structure (D23).** Faces as node lists, `owner` and `neighbour`
  per face, the faces of each cell with their orientation; patches and cell
  zones. Order: inner faces, processor faces, boundary faces by patch;
  cells inner, owned cells touching the halo, halo cells by source rank,
  each group in reverse Cuthill-McKee order (R1). u64 global ids for every
  cell, face and node in input order (R2). `validate()` checks index ranges,
  ordering and topological closure (every edge used twice).
- **Geometry.** Face area vectors and centroids from a triangle fan about
  the node average, cell volumes and centroids from tetrahedra on the same
  fan: Σ S = 0 to round-off for every cell, volumes add up to the domain
  volume.
- **Quality report.** Closure, total against enclosed volume, patch areas,
  non-orthogonality, skewness, face warp, aspect ratio, volume ratio,
  inverted faces, non-positive volumes; Markdown and JSON.
- **Gmsh reader (D33).** MSH 4.1 ASCII and binary (either endianness),
  tetrahedra, hexahedra, prisms, pyramids; physical surfaces as patches,
  physical volumes as cell zones, periodic links; refuses higher-order
  elements and other versions with a clear message.
- **Test mesher (D34).** Box meshes of hexahedra, tetrahedra (Kuhn split),
  polyhedra (dual of the tetrahedra), mixed cells; 3D or 2D one cell thick
  (D3); exponential grading per direction, one- or two-sided; nested
  refinement levels; node perturbation by the counter-based generator,
  boundary nodes kept on their faces. Conversion of any tetrahedral mesh into
  polyhedra with the exact volume (also on curved walls).
- **Renumbering.** Reverse Cuthill-McKee and shuffle, with the ids kept.
- **Partitioning (D32).** METIS (vendored) and recursive coordinate
  bisection; imbalance, cut faces and neighbours in the report.
- **Decomposition and halo (D38, R3).** Local meshes with halo depth as a
  parameter; per-neighbour send lists and contiguous receive ranges;
  several fields per message; host-staged transport now, GPU-aware selected
  at run time later on the same plan. Halo geometry filled by exchange.
- **VTKHDF writer (D41).** VTKHDF 2.4 `PartitionedDataSetCollection`:
  block `internal` with cell fields and one block per patch; standard VTK
  cells and `VTK_POLYHEDRON`; global-id order. HDF5 built from source and
  linked statically (needs CMake).
- **Point location and tracking (D47).** Tetrahedral-decomposition walk
  from cell to cell through the faces, with a bin search where the walk
  stops at a non-convex boundary; reports where a track leaves the mesh or
  the rank.

Tests: closure, volume and patch areas of 16 generated variants; Gmsh ASCII
against binary, zones, prisms, pyramids, periodic links, refusals; VTKHDF of
all types; point location on all types; decomposition with METIS and RCB on
1, 2, 3 and 5 parts and 1 and 2 halo layers; **golden values (R4)** for 32
cases of mesh type × dimension × perturbation × refinement × 1, 2 and 4
ranks, checking that owned volumes sum to the total and that every halo
cell has exactly one owner; the MPI halo exchange on 2 and 4 ranks in CI.

CLI: `cfdcraft mesh CASE.toml | FILE.msh` builds or reads, renumbers,
checks, partitions and writes the report and `mesh.vtkhdf`. G1:
`mesh.type = "poly" | "tet" | "hex" | "mixed"` runs the benchmark on
stage-1 meshes.

Tools: `tools/check_vtkhdf.py` (checks files with VTK's reader),
`tools/render_mesh.py` (gallery pictures), `tools/gmsh_test_meshes.py`
(the Gmsh test files), `examples/meshes/tjunction.geo`.

Documentation: gallery examples 3–7, mesh cases and keys in the case
catalogue, recipes 7–11, the G1 stage-1 rerun, implementation notes in
D23, D32, D38, D41 and D47, new open item O15.

Findings:

- Polyhedral boundary faces merged across a curved wall lose volume
  (0.34 % on the T-junction); keeping one planar piece per triangle there
  makes the dual exactly as large as the tetrahedral mesh.
- VTK computes the volume of a polyhedron as that of the convex hull of its
  nodes, so the VTKHDF check computes it from VTK's faces instead.
- **G1 on polyhedra (shared cloud machine, 2 million cells):** bitwise
  equal; the single-source bodies are 25–42 % faster on cell gather, equal
  on SpMV-CSR, 4–30 % slower on face flux, as on a hexahedral box of the
  same size. Shuffled order costs a factor 3 to 4. Quiet-machine and GPU
  reruns are open (O15).

---

## 0.0 (2026-11-14) – stage 0: foundations and gate G1

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

CLI: `cfdcraft info` reports the communicator and the rank count.

Documentation: [gallery](docs/gallery.md), [case catalogue](docs/cases.md),
[recipes](docs/recipes.md), [G1 procedure and results](docs/benchmarks/g1.md).

Findings:

- **G1 GPU verdict, Tesla T4 (Google Colab), 10 million cells: PASS.** The
  single-source Rust kernels, compiled to PTX, run within 1.0 % of the
  hand-written CUDA C kernels on all four kernels (220–255 GB/s of 320 GB/s),
  in two runs, with results bitwise equal to the CPU. D8: one kernel source;
  D11 stands. The single-source GPU route needs nightly Rust (O14).
- **G1 CPU verdict, Intel i9-9900K (macOS), 10 million cells: PASS** in two
  runs: largest time ratio single-source against hand-written Rust loops
  1.018, results bitwise equal, 20–26 GB/s.
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

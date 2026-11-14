# CFDCraft – Requirements

Living document, updated after each requirement step. Form taken from the
[PyAeroFlow requirements](https://github.com/chiefenne/pyaeroflow/blob/main/docs/requirements.md).

- **D** – decided
- **L** – planned for later (decided, not in the first milestone)
- **P** – possible later (mentioned, not decided)
- **X** – not planned
- **O** – open requirement step, with the current proposal

Status: requirements engineering in progress. Steps 1 (scope) and 2 (language,
hardware) answered on 2026-11-14; step 3 open. No code yet.


## 1. Purpose and scope

- General purpose CFD solver for fluid dynamics, heat and mass transfer,
  species transport and combustion. Long-term reference class: the
  industry-standard general purpose CFD solvers.
- Open source (MIT), developed in the manner of the
  [ArtCraft Crafting Apps](https://getartcraft.com/apps/craft-anything): Rust,
  one engine with thin front ends on top (first a CLI), built step by step
  with a coding agent in an explore, build, check, refine loop.
- First milestone: incompressible flow on unstructured 3D meshes with two or
  three turbulence models, run from a case file by an expert user.
- Top goal: performance on modern parallel GPUs, at 100+ million cells. Runs
  on any CPU as well.
- Step-by-step development; every stage ends with results that can be tested
  against analytic solutions or published reference data.


## 2. Decisions

| ID  | Step | Topic | Decision | Notes |
|-----|------|-------|----------|-------|
| D1  | –  | Repository | github.com/chiefenne/CFDCraft, public from the start, MIT licence | |
| D2  | –  | Project form | Open source; one maintainer at the start, with unlimited time; the project must stay open to contributors | Contributor rules: later step |
| D3  | 1  | Dimension | 3D only, arbitrary polyhedral cells from day one. 2D cases run as a 3D mesh one cell layer thick, with symmetry conditions on the two side planes | No separate 2D core |
| D4  | 1  | Product boundary, milestone 1 | Solver only: case file plus CLI; meshes from Gmsh or other tools; results to ParaView | Desktop app: L7. Mesher: L8 |
| D5  | 1  | Reference user, milestone 1 | An expert who writes a text case file | Guidance for users coming from GUI codes is not a milestone-1 requirement |
| D6  | 1  | Problem size and parallelism | Design target 100+ million cells. GPU first; several GPUs and clusters with MPI are required, so domain decomposition is part of the data structures from the start | 100 million cells at float64 need roughly 100–200 GB for a segregated solver, more than one GPU holds: multi-GPU is not optional at the target size |
| D7  | 2  | Precision | float64 is mandatory for the reference solution path on GPU and CPU. The floating-point type is a generic parameter of the code, so a float32 or mixed-precision mode costs no second code base | Consequence for Apple GPUs: O6 |
| D8  | 2  | GPU backend policy | Top speed decides. One kernel source for all targets is preferred, but only if it loses no performance against a native implementation on the primary GPU platform; otherwise two backends | Measured in stage 0 (gate G1, section 7) |
| D9  | 2  | Foreign libraries | Linking C libraries is allowed where it costs no performance. MPI is linked as a system library | A call into a C library at solver granularity has no measurable overhead |
| D10 | 2  | Solver ecosystem | Sparse linear algebra, Krylov solvers and algebraic multigrid are written in this project, from the published literature, provided this proves feasible | Feasibility: the methods are well published (aggregation and classical AMG, GPU variants). This is the largest single work item of the project. Proposal for a safety net: O8 |
| D11 | 2  | Language | Rust (stable toolchain, cargo workspace), conditional on gate G1 | Alternatives considered: C++ with Kokkos, modern Fortran with OpenACC. See section 2a |


## 2a. Language assessment (step 2)

Facts that shaped D11:

- With float64 mandatory (D7), Apple GPUs cannot run the reference path in
  any language, so Metal no longer argues for or against Rust.
- On NVIDIA, a GPU kernel compiles to the same device code whatever language
  hosts it. Rust can launch CUDA kernels with no overhead, so Rust costs no
  GPU speed. The open point is tooling: whether kernels can be written once
  in Rust for CUDA and CPU (and other vendors) without a penalty. Gate G1
  measures that before any solver code is written.
- What C++ with Kokkos would give for free and Rust does not: a proven
  performance-portability layer and direct use of hypre, PETSc and Trilinos.
  D10 accepts building the solver part in-house.
- What Rust gives: compile-time memory and data-race safety for the large
  host side (mesh, partitioning, I/O, configuration, orchestration), one
  build tool, and a good fit for agent-driven development, where the compiler
  catches whole classes of errors.
- Risk accepted: fewer CFD contributors know Rust than C++ or Fortran.

If G1 fails (no route within 10 % of a hand-written CUDA reference), the
fallback is two backends (D8) with CUDA kernels in CUDA C launched from Rust;
only if that also proves unworkable is the language reopened.


## 3. Planned for later

In the order given in step 1.

| ID | Feature | Notes |
|----|---------|-------|
| L1 | Energy equation and conjugate heat transfer | Multi-region meshes (fluid and solid) must not be blocked by the mesh data structure |
| L2 | Species transport and combustion | Needs variable density (low-Mach formulation) and stiff chemistry integration |
| L3 | Multiphase: VOF and Euler–Euler | Both required |
| L4 | Moving meshes | |
| L5 | Scale-resolving turbulence: LES, DES | Needs low-dissipation convection schemes and second-order time integration |
| L6 | Clusters at scale | MPI is in the design from the start (D6); tuning and scaling tests on real clusters come after milestone 1 |
| L7 | Desktop app, in the style of the Crafting Apps | On the same engine as the CLI |
| L8 | Built-in mesher | Much later; solver first |
| L9 | Logo in the ArtCraft style (engraved animal) | After the requirements; animal open |


## 4. Possible later (not decided)

| ID | Feature | Notes |
|----|---------|-------|
| P1 | Fully compressible flow (transonic, supersonic) | Not named in step 1; to be confirmed (see O9) |
| P2 | float32 or mixed-precision GPU mode | Enabled by D7; would let Apple GPUs take part (O6) |
| P3 | Web build (WebAssembly), MCP server | Part of the ArtCraft form; no requirement yet |


## 5. Not planned

None yet.


## 6. Open requirement steps

| ID | Step | Topic | Current proposal |
|----|------|-------|------------------|
| O3 | 3 | Pressure–velocity coupling | Pressure-based, collocated finite volume with Rhie–Chow interpolation; segregated SIMPLE-family algorithm (SIMPLEC steady, PISO-type transient), assembled sparse matrices, AMG for the pressure equation. A pressure-based coupled solver later. Reason: it is the route of the industry-standard solvers, extends to energy, species, combustion and multiphase, and its cost concentrates in the linear solvers of D10. PyAeroFlow's coupled Newton–Krylov with direct LU does not reach 100 million cells |
| O4 | 3 | Turbulence models and wall treatment, milestone 1 | k-ω SST (2003), Spalart–Allmaras (SA-neg), realizable k-ε. Wall treatment valid for any y⁺ (wall-resolved and wall functions), because general purpose meshes cannot guarantee y⁺ ≈ 1 |
| O5 | 3 | Steady and transient | Both in milestone 1 |
| O6 | 2 | Apple platforms | float64 on the CPU of a Mac (multi-threaded, SIMD) is the supported path. Metal GPU only through P2, as a declared reduced-precision mode, after milestone 1 |
| O7 | 2 | GPU vendors | NVIDIA (CUDA) is the primary platform. AMD and Intel GPUs: required or not? Several large European clusters are AMD |
| O8 | 2 | Reference linear solver | hypre as an optional, linkable reference behind the same interface as the in-house solvers (D10), to verify them and to measure them against the state of the art |
| O9 | 1 | Compressible flow | Low-Mach variable density with L2; fully compressible as P1 |
| O10 | – | Clean-room rule | Numerics from the published literature only; no code taken from other solvers, so the project stays MIT (as in PyAeroFlow) |
| O11 | – | Continuous integration | CPU tests on GitHub for every pull request (unlike PyAeroFlow X5), since contributors are wanted; GPU tests run locally |

Steps still to come: discretization (gradients, convection and diffusion
schemes, limiters, non-orthogonal correction), time integration, linear
solvers in detail, mesh formats and partitioning, boundary conditions, case
file format, output and restart, verification and validation cases, stage
plan.


## 7. Stage plan

Not yet defined (last requirement step). One item is fixed already:

| Stage | Content | Tests, acceptance |
|-------|---------|-------------------|
| **0 Gate G1: kernel benchmark** | Two representative kernels (a face-based flux loop on an unstructured mesh, a sparse matrix-vector product), float64, at 10–50 million cells: hand-written CUDA reference against the candidate single-source Rust routes; the same kernels on CPU | A single-source route within 10 % of the CUDA reference confirms D11 with one kernel source; otherwise two backends (D8) |


## References

To be filled per step.

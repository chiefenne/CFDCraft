# CFDCraft – Requirements

Living document, updated after each requirement step. Form taken from the
[PyAeroFlow requirements](https://github.com/chiefenne/pyaeroflow/blob/main/docs/requirements.md).

- **D** – decided
- **L** – planned for later (decided, not in the first milestone)
- **P** – possible later (mentioned, not decided)
- **X** – not planned
- **O** – open requirement step, with the current proposal

Status: requirements engineering in progress. Steps 1 (scope) and 2 (language,
hardware) answered on 2026-11-14; step 3 mostly answered on 2026-11-14 (open:
O3 pressure–velocity coupling, O8, O10, O11). No code yet.

Terminology: the project says *reacting flows*, not "combustion", everywhere.


## 1. Purpose and scope

- General purpose CFD solver for fluid dynamics, heat and mass transfer,
  species transport and reacting flows. Long-term reference class: the
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
| D12 | 3  | Turbulence models, milestone 1 | k-ω SST (Menter, Kuntz & Langtry 2003); standard k-ε (Launder & Spalding 1974); the modern k-ε variants realizable (Shih et al. 1995) and RNG (Yakhot et al. 1992) | The k-ε variants are there to reproduce literature data. All models through one generic transport-equation interface; the three k-ε variants share one implementation and differ in constants and source terms. Spalart–Allmaras: P4 |
| D13 | 3  | Wall treatment, milestone 1 | Wall functions are required: standard (log-law) and scalable wall functions for the k-ε family; an automatic, y⁺-insensitive wall treatment for SST that blends the viscous sublayer and the log law | General purpose meshes cannot guarantee y⁺ ≈ 1. Low-Re, wall-resolved k-ε: P5 |
| D14 | 3  | Steady and transient | Both in milestone 1 | |
| D15 | 3  | GPU vendor | NVIDIA (CUDA) only | AMD and Intel GPUs: P6. The single-source question of D8 and gate G1 is therefore NVIDIA against CPU |
| D16 | 3  | Apple platforms | A Mac runs the float64 path on its CPU (multi-threaded, SIMD). Apple GPUs only through P2, as a declared reduced-precision mode | Follows from D7: Apple GPUs have no float64 hardware and the Metal shading language has no `double` |
| D17 | 3  | Compressibility | Low-Mach variable density first, together with reacting flows (L2); fully compressible high-Mach flow later (L10) | |
| D18 | –  | Documentation per feature | Every implemented feature ships with (1) a gallery example (`docs/gallery.md`: one picture per case with the commands that reproduce it), (2) one or more commented case files in `examples/cases/` listed in the case catalogue (`docs/cases.md`: what the case is, the mesh it needs, where it is discussed), (3) a recipe (`docs/recipes.md`: short command sequence for the task, naming the case-file keys it uses). A feature is not done without all three | Owner rule of 2026-11-14, same form as PyAeroFlow |


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


## 2b. Segregated or coupled (step 3, open)

Is the segregated approach of the industry standards legacy, or an
advantage? Both, depending on the problem:

- **Memory.** A coupled solver for u, v, w, p stores a 4×4 block per matrix
  entry, 16 times the entries of one scalar matrix, and its multigrid
  hierarchy on top. At 100 million polyhedral cells (about 15 entries per row)
  the coupled matrix alone is about 190 GB in float64; one scalar matrix is
  about 12 GB and is reused for every equation. On GPUs, where memory is the
  scarce resource, this is a real advantage of the segregated approach.
- **Transient flow.** With small time steps (PISO type, CFL around 1) a
  segregated step needs few outer corrections and is cheap. The coupled
  solver's main strength, fewer outer iterations, matters little there.
- **Further physics.** Energy, species, turbulence and volume fractions are
  added as further scalar equations. Segregated solvers extend naturally;
  coupled solvers usually couple only momentum and pressure and solve the rest
  segregated anyway.
- **Steady, stiff flows.** Here coupled wins: the SIMPLE family needs
  under-relaxation, and its convergence degrades on high-aspect-ratio
  boundary-layer meshes and with strong body forces. A coupled solver
  typically converges in several times fewer iterations and less wall time,
  at roughly 1.5–2 times the memory. This is why the industry standards offer
  a coupled pressure-based solver next to the segregated one, and recommend
  it for steady external aerodynamics.
- **GPU fit.** Both fit: block-sparse products have a higher arithmetic
  intensity than scalar ones, so the coupled solver uses a GPU better per
  byte, but needs a block algebraic multigrid, which is harder to build
  in-house (D10) than a scalar one.

So the segregated default is not only legacy: it is the better choice for
transient, multi-physics and memory-bound runs, and the coupled solver is the
better choice for steady, stiff ones. The proposal is to build both on one
discretization, segregated first (it is needed for every later physics and
needs only scalar AMG) and coupled second.


## 3. Planned for later

In the order given in step 1.

| ID | Feature | Notes |
|----|---------|-------|
| L1 | Energy equation and conjugate heat transfer | Multi-region meshes (fluid and solid) must not be blocked by the mesh data structure |
| L2 | Species transport and reacting flows | Low-Mach variable density (D17) and stiff chemistry integration |
| L3 | Multiphase: VOF and Euler–Euler | Both required |
| L4 | Moving meshes | |
| L5 | Scale-resolving turbulence: LES, DES | Needs low-dissipation convection schemes and second-order time integration |
| L6 | Clusters at scale | MPI is in the design from the start (D6); tuning and scaling tests on real clusters come after milestone 1 |
| L7 | Desktop app, in the style of the Crafting Apps | On the same engine as the CLI |
| L8 | Built-in mesher | Much later; solver first |
| L9 | Logo in the ArtCraft style | Peregrine falcon, first version in the repository (`assets/logo/cfdcraft.svg`); variants under review |
| L10 | Fully compressible high-Mach flow (transonic, supersonic) | After low-Mach variable density (D17) |


## 4. Possible later (not decided)

| ID | Feature | Notes |
|----|---------|-------|
| P2 | float32 or mixed-precision GPU mode | Enabled by D7; would let Apple GPUs take part (O6) |
| P3 | Web build (WebAssembly), MCP server | Part of the ArtCraft form; no requirement yet |
| P4 | Spalart–Allmaras | Not in milestone 1 (D12); cheap to add through the transport-equation interface |
| P5 | Low-Re, wall-resolved k-ε variants | D13 |
| P6 | AMD and Intel GPUs | D15 |


## 5. Not planned

None yet.


## 6. Open requirement steps

| ID | Step | Topic | Current proposal |
|----|------|-------|------------------|
| O3 | 3 | Pressure–velocity coupling | Pressure-based, collocated finite volume with Rhie–Chow interpolation. Segregated SIMPLE-family algorithm first (SIMPLEC steady, PISO-type transient), coupled pressure-based solver (momentum and pressure in one block system) as a planned second algorithm. Both work on the same discretization: the assembly produces per-face coefficients from which either the scalar matrices or the block matrix are built. Assessment in section 2b |
| O8 | 2 | Reference linear solver | hypre (open-source solver library of Lawrence Livermore National Laboratory, MIT or Apache-2.0 licence; best known for the BoomerAMG algebraic multigrid; MPI and CUDA) as an optional, linkable reference behind the same interface as the in-house solvers (D10), to verify them and to measure them against the state of the art |
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

- F. R. Menter, M. Kuntz, R. Langtry, *Ten years of industrial experience with
  the SST turbulence model*, Turbulence, Heat and Mass Transfer 4 (2003)
- B. E. Launder, D. B. Spalding, *The numerical computation of turbulent
  flows*, Comput. Methods Appl. Mech. Eng. 3, 269–289 (1974)
- T.-H. Shih, W. W. Liou, A. Shabbir, Z. Yang, J. Zhu, *A new k-ε eddy
  viscosity model for high Reynolds number turbulent flows*, Computers &
  Fluids 24, 227–238 (1995)
- V. Yakhot, S. A. Orszag, S. Thangam, T. B. Gatski, C. G. Speziale,
  *Development of turbulence models for shear flows by a double expansion
  technique*, Phys. Fluids A 4, 1510–1520 (1992)
- C. M. Rhie, W. L. Chow, *Numerical study of the turbulent flow past an
  airfoil with trailing edge separation*, AIAA J. 21, 1525–1532 (1983)
- R. I. Issa, *Solution of the implicitly discretised fluid flow equations
  by operator-splitting*, J. Comput. Phys. 62, 40–65 (1986)
- J. P. Van Doormaal, G. D. Raithby, *Enhancements of the SIMPLE method for
  predicting incompressible fluid flows*, Numer. Heat Transfer 7, 147–163
  (1984)
- NASA Turbulence Modeling Resource: https://tmbwg.github.io/turbmodels/
- hypre: https://github.com/hypre-space/hypre

# CFDCraft – Requirements

Living document, updated after each requirement step. Form taken from the
[PyAeroFlow requirements](https://github.com/chiefenne/pyaeroflow/blob/main/docs/requirements.md).

- **D** – decided
- **L** – planned for later (decided, not in the first milestone)
- **P** – possible later (mentioned, not decided)
- **X** – not planned
- **O** – open requirement step, with the current proposal

Status: requirements engineering in progress. Steps 1 (scope) and 2 (language,
hardware) answered on 2026-11-14; steps 3 to 6 (turbulence and platform, discretization, mesh input and solvers, case file and output) answered on 2026-11-14. No code yet.

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
| D19 | 3  | Pressure–velocity coupling | Pressure-based, collocated finite volume with Rhie–Chow interpolation, on one shared discretization (the assembly produces per-face coefficients from which scalar or block matrices are built). Segregated SIMPLE-family algorithm first (SIMPLEC steady, PISO-type transient); a coupled pressure-based solver (momentum and pressure in one block system) as the planned second algorithm (L11) | Segregated first: every later physics needs it, and it needs only scalar AMG. Assessment in section 2b |
| D20 | 3  | Reference linear solver | hypre (open-source solver library of Lawrence Livermore National Laboratory, MIT or Apache-2.0 licence; BoomerAMG, MPI, CUDA) as an optional, linkable reference behind the same interface as the in-house solvers (D10), used to verify them and to benchmark them against the state of the art. Never required to build or run CFDCraft | |
| D21 | 3  | Clean-room rule | Methods come from the published literature, with references. Code from other solvers is never copied. Open-source solvers may be read for understanding where their licence allows it | Keeps the project under the MIT licence (D1) |
| D22 | 3  | Continuous integration | GitHub Actions runs the CPU tests (build, unit tests, small verification cases) on every pull request and every push to `main`, on GitHub's free hosted runners. Free for public repositories. GPU tests are run locally by the maintainer once suitable hardware is available; they are marked separately and are not part of the hosted checks | The maintainer has no GPU hardware yet (2026-11-14), so GPU verification starts when it is available |
| D23 | 4  | Spatial discretization | Cell-centred finite volume on arbitrary polyhedral cells, face-based data structure (owner and neighbour cell per face) | Every later physics fits it (multiphase, moving meshes, conjugate heat transfer) |
| D24 | 4  | Gradient reconstruction | Weighted least squares on the face-neighbour stencil as the default; Green-Gauss as a cheaper option | Plain Green-Gauss is only zero-order accurate on skewed polyhedral meshes (Syrakos et al. 2017, as in PyAeroFlow D22) |
| D25 | 4  | Convection schemes | Second-order linear upwind with a limiter as the default; first-order upwind for start-up; a blended central scheme as an option (needed for LES, L5). Limiters in milestone 1: Venkatakrishnan and Venkatakrishnan–Wang | Further TVD limiters later through the same interface |
| D26 | 4  | Diffusion on non-orthogonal meshes | Over-relaxed non-orthogonal correction, applied as a deferred correction with a limiter on strongly skewed cells | |
| D27 | 4  | Rhie–Chow interpolation | Time-step-independent form: steady results do not depend on under-relaxation or time step | |
| D28 | 4  | Time integration | Second-order backward (BDF2) as the default, first order for the first step; Crank–Nicolson as an option; optional adaptive time step controlled by CFL | |
| D29 | 4  | Boundary conditions, milestone 1 | Velocity inlet, mass-flow inlet, static-pressure inlet, total-pressure inlet, static-pressure outlet, total-pressure outlet, no-slip wall (also moving), symmetry, periodic. Turbulence values at inlets as turbulence intensity and length scale, or viscosity ratio. Static and total pressure conditions can be combined, for example static pressure at an inlet with total pressure at an outlet | The compatible combinations are enforced by D30 |
| D30 | 4  | Boundary-condition check and guidance | Before the run, a check of the combination of boundary conditions on the whole domain: a pressure level must be fixed (a pressure boundary or a reference cell), the net mass balance must be possible, over-specified combinations are errors (for example velocity inlet together with pressure inlet on the same flow path), and unused or unsuitable turbulence values are reported. Errors and warnings name the offending boundary and say what to change. The documentation carries a guidance table of sensible combinations and the cases each one suits | Same idea as PyAeroFlow D39; a recipe (D18) shows typical combinations |
| D31 | 4  | Reproducible results | Results are reproducible from run to run on the same hardware and decomposition: gather-based assembly of face contributions with no floating-point atomics, and a fixed summation order in reductions | Mandatory; costs possibly a few percent of speed. Reproducibility across different GPU or rank counts is not promised, because the summation order changes with the decomposition |
| D32 | 4  | Mesh partitioning | METIS or Scotch (linked C library, allowed by D9) for domain decomposition | Decision between the two at implementation |
| D33 | 5  | Mesh input, milestone 1 | Gmsh `.msh` version 4 (ASCII and binary), with physical groups as the named boundary-face sets | Further formats (CGNS, other open formats) later through the same reader interface |
| D34 | 5  | Built-in mesh generator for canonical test cases | In code: block-structured hexahedral meshes with varying cell size and boundary-layer grading (2D as one cell layer thick with symmetry side planes, D3, and full 3D); tetrahedral meshes (hexahedra split); polyhedral meshes (dual of a tetrahedral mesh); mixed meshes combining these. Optional random node perturbation and refinement levels for order-of-accuracy studies | Used by the verification cases; same role as PyAeroFlow D42. Not the general mesher of L8 |
| D35 | 5  | Pressure linear solver | Aggregation-based algebraic multigrid as preconditioner of CG, V-cycle, Jacobi or Chebyshev smoothers (GPU friendly). ILU-type smoothers excluded on GPU | In-house (D10), verified against hypre (D20) |
| D36 | 5  | Momentum and scalar linear solvers | BiCGStab with Jacobi preconditioning, a few iterations to a loose relative tolerance per outer iteration | Standard practice for segregated solvers |
| D37 | 5  | Matrix storage and precision in the linear solvers | Matrix stored once in a GPU-friendly sparse format (ELL or sliced ELL; final choice from gate G1). The multigrid hierarchy is stored and applied in float32; the outer Krylov iteration, residuals and solution stay in float64 (D7) | Mixed precision inside the preconditioner halves memory and bandwidth of the hierarchy without changing the converged result, because the outer float64 iteration corrects it |
| D38 | 5  | Parallel communication | Halo exchange overlapped with interior computation; GPU-aware MPI. A multi-node run is a milestone-1 acceptance test | Needs access to a multi-node machine (O12) |
| D39 | 6  | Case file | TOML (standard-library style parser in Rust), literature names and established acronyms, every key with a sensible default optional, validation before the run (unknown keys with "did you mean", types, ranges, consistency), resolved configuration with defaults written to the output directory. Same principles as PyAeroFlow D26–D30 | |
| D40 | 6  | Re-reading the case file during a run | The case file can be re-read while a run is going on, on a schedule the user sets: every iteration, every time step, every N iterations or time steps, or never (default). Only settings that are safe to change at run time are applied (output, restart and monitoring settings, stop criteria, relaxation factors, time step); every change is logged, and a change to anything else is reported as not applied. In addition the control files of PyAeroFlow L7 are kept: `.stop` (write a restart and the results, end cleanly), `.result_N` (write results at iteration or step N), `.read_toml_N` (read the case file at iteration or step N) | PyAeroFlow L7 (5) and (6) |
| D41 | 6  | Field output | VTKHDF: one parallel HDF5 file per write (one file per time step or a transient series in one file), read directly by ParaView. Linked HDF5 library (D9) | Per-rank VTU files are not workable at 100+ million cells |
| D42 | 6  | Restart | Same HDF5 file layout, written in parallel by all ranks, readable with a different number of ranks or GPUs. A restart continues the run (counters, time, solver state), as PyAeroFlow D56 and D58 | |
| D43 | 6  | Monitoring and user-configurable output | A human-readable `.log` file, CSV histories (residuals, forces, integrals), one CSV of wall quantities per boundary, point probes. Which fields, monitors, probes and wall quantities are written, and how often (per iteration, per time step, every N, at listed steps), is chosen by the user in the case file; defaults are modest. The log header states what is written and how large the output will be | Same idea as PyAeroFlow D59 |
| D44 | 6  | Units | SI on input and output and dimensional inside the solver | PyAeroFlow's non-dimensional form (its D51) is not adopted: reacting flows have many reference scales |


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


## 2b. Segregated or coupled (step 3, decided in D19)

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
| L11 | Coupled pressure-based solver | Second pressure–velocity algorithm on the shared discretization (D19); needs a block AMG |
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
| O12 | 5 | Multi-node acceptance test | D38 makes a multi-node run a milestone-1 test. The maintainer has an old PC with an NVIDIA GPU (model to be determined) and does not want to rent hardware. Proposal: the milestone-1 tests are multi-rank MPI runs on one machine (CPU ranks in CI, D22, and ranks sharing the one GPU locally); a true multi-node run is accepted when a cluster allocation is available (university or research-computing access, a community contribution) and is then added as a test. Until then the multi-node path is designed for and tested with several ranks on one machine, not claimed as verified |

Steps still to come: verification and validation cases, stage plan.


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

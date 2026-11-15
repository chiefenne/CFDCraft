# CFDCraft – Requirements

Living document, updated after each requirement step. Form taken from the
[PyAeroFlow requirements](https://github.com/chiefenne/pyaeroflow/blob/main/docs/requirements.md).

- **D** – decided
- **L** – planned for later (decided, not in the first milestone)
- **P** – possible later (mentioned, not decided)
- **X** – not planned
- **O** – open requirement step, with the current proposal

Status: requirements engineering in progress. Steps 1 (scope) and 2 (language,
hardware) answered on 2026-11-14; steps 3 to 7 (turbulence and platform, discretization, mesh input and solvers, case file and output, test cases) answered on 2026-11-14; the stage plan was confirmed in order by the maintainer on 2026-11-14 (stage 7 changed to internal flows) and sprays (L12) were added. Stage 0 completed on 2026-11-14 (tag `v0.0`, gate G1 passed; see `docs/benchmarks/g1.md`). Stage 1 (mesh) completed on 2026-11-14 (release 0.1): all acceptance tests pass; the G1 reruns on polyhedral meshes on a quiet machine and on the GPU passed (O15 closed). Stage 2 (operators on a scalar equation) completed on 2026-11-14 (release 0.2): all CPU acceptance tests pass; the GPU run against the CPU passed on a Tesla T4 (Colab), bitwise equal (docs/verification/stage2.md). From the MFC assessment the maintainer adopted R5, R6, R7 and R10 for stage 2 (2026-11-14).

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
| D8  | 2  | GPU backend policy | Top speed decides. One kernel source for all targets is preferred, but only if it loses no performance against a native implementation on the primary GPU platform; otherwise two backends | Measured in stage 0 (gate G1, section 7). Result of 2026-11-14: on a Tesla T4 the single-source Rust kernels are within 1.0 % of the hand-written CUDA C kernels, on an i9-9900K within 1.8 % of hand-written Rust loops (limit 10 %), bitwise equal results, so **one kernel source** is used for the GPU and the CPU on this evidence. Reservation: memory-bound kernels on a well-ordered mesh; rerun on shuffled and polyhedral meshes (stage 1) and on the first compute-heavy kernels (docs/benchmarks/g1.md) |
| D9  | 2  | Foreign libraries | Linking C libraries is allowed where it costs no performance. MPI is linked as a system library | A call into a C library at solver granularity has no measurable overhead |
| D10 | 2  | Solver ecosystem | Sparse linear algebra, Krylov solvers and algebraic multigrid are written in this project, from the published literature, provided this proves feasible | Feasibility: the methods are well published (aggregation and classical AMG, GPU variants). This is the largest single work item of the project. Proposal for a safety net: O8 |
| D11 | 2  | Language | Rust (stable toolchain, cargo workspace), conditional on gate G1 | Alternatives considered: C++ with Kokkos, modern Fortran with OpenACC. See section 2a. **G1 passed on 2026-11-14 (GPU: T4, CPU: i9-9900K), D11 stands.** Deviation found: the single-source GPU route (`gpu/ptx-kernels`, rustc target `nvptx64-nvidia-cuda`) needs the nightly toolchain, so the "stable toolchain" holds for the host code only (O14) |
| D12 | 3  | Turbulence models, milestone 1 | k-ω SST (Menter, Kuntz & Langtry 2003); standard k-ε (Launder & Spalding 1974); realizable k-ε (Shih et al. 1995); k-ε-ζ-f (Hanjalić, Popovac & Hadžiabdić 2004), the elliptic-relaxation eddy-viscosity model, **in place of RNG k-ε** (decision of 2026-11-14, assessment `k-zeta-f-assessment.md`) | The k-ε variants are there to reproduce literature data; ζ-f is the near-wall model for separation and heat transfer (best of the linear eddy-viscosity models on the asymmetric diffuser and, with L1, on impingement and near-wall heat transfer). All models through one generic transport-equation interface, which also takes the steady elliptic-relaxation equation of ζ-f (Helmholtz operator, no time derivative, solved with the scalar path of stage 2, no AMG); the two k-ε variants share one implementation. The time and length scales are bounded as in the 2004 paper (realizability limits). Inlets carry an optional ζ (default 2/3). RNG k-ε is not planned (X1). Spalart–Allmaras: P4. PANS on ζ-f: P7 |
| D13 | 3  | Wall treatment, milestone 1 | Wall functions are required: standard (log-law) and scalable wall functions for the k-ε family; an automatic, y⁺-insensitive wall treatment for SST that blends the viscous sublayer and the log law; the compound wall treatment (Popovac & Hanjalić 2007; Popovac 2023) for ζ-f, blending the integrated and the wall-function limits, including the thermal field once L1 exists | General purpose meshes cannot guarantee y⁺ ≈ 1. The wall-treatment interface takes first-cell y⁺, wall distance of the first cell and a blending function per wall face and is not specific to SST, so that the compound treatment and the thermal wall function of L1 use the same path. Low-Re, wall-resolved k-ε: P5 |
| D14 | 3  | Steady and transient | Both in milestone 1 | |
| D15 | 3  | GPU vendor | NVIDIA (CUDA) only | AMD and Intel GPUs: P6. The single-source question of D8 and gate G1 is therefore NVIDIA against CPU |
| D16 | 3  | Apple platforms | A Mac runs the float64 path on its CPU (multi-threaded, SIMD). Apple GPUs only through P2, as a declared reduced-precision mode | Follows from D7: Apple GPUs have no float64 hardware and the Metal shading language has no `double` |
| D17 | 3  | Compressibility | Low-Mach variable density first, together with reacting flows (L2); fully compressible high-Mach flow later (L10) | |
| D18 | –  | Documentation per feature | Every implemented feature ships with (1) a gallery example (`docs/gallery.md`: one picture per case with the commands that reproduce it), (2) one or more commented case files in `examples/cases/` listed in the case catalogue (`docs/cases.md`: what the case is, the mesh it needs, where it is discussed), (3) a recipe (`docs/recipes.md`: short command sequence for the task, naming the case-file keys it uses). A feature is not done without all three | Owner rule of 2026-11-14, same form as PyAeroFlow |
| D19 | 3  | Pressure–velocity coupling | Pressure-based, collocated finite volume with Rhie–Chow interpolation, on one shared discretization (the assembly produces per-face coefficients from which scalar or block matrices are built). Segregated SIMPLE-family algorithm first (SIMPLEC steady, PISO-type transient); a coupled pressure-based solver (momentum and pressure in one block system) as the planned second algorithm (L11) | Segregated first: every later physics needs it, and it needs only scalar AMG. Assessment in section 2b |
| D20 | 3  | Reference linear solver | hypre (open-source solver library of Lawrence Livermore National Laboratory, MIT or Apache-2.0 licence; BoomerAMG, MPI, CUDA) as an optional, linkable reference behind the same interface as the in-house solvers (D10), used to verify them and to benchmark them against the state of the art. Never required to build or run CFDCraft | **Stage 3:** cargo feature `hypre`, a C shim over hypre's public interface (IJ, ParCSR PCG, BoomerAMG), hypre built from source by `tools/build_hypre.sh` (sequential, shared; CPU in CI, CUDA on Colab); nothing redistributed, no code copied (D21) |
| D21 | 3  | Clean-room rule | Methods come from the published literature, with references. Code from other solvers is never copied. Open-source solvers may be read for understanding where their licence allows it | Keeps the project under the MIT licence (D1) |
| D22 | 3  | Continuous integration | GitHub Actions runs the CPU tests (build, unit tests, small verification cases) on every pull request and every push to `main`, on GitHub's free hosted runners. Free for public repositories. GPU tests are run locally by the maintainer once suitable hardware is available; they are marked separately and are not part of the hosted checks | The maintainer has no GPU hardware yet (2026-11-14), so GPU verification starts when it is available |
| D23 | 4  | Spatial discretization | Cell-centred finite volume on arbitrary polyhedral cells, face-based data structure (owner and neighbour cell per face) | Every later physics fits it (multiphase, moving meshes, conjugate heat transfer). **Implementation (stage 1, `cfdcraft-mesh`):** faces stored as node lists (CSR), `owner` and `neighbour` per face (u32 local indices), the faces of each cell as a CSR list with the orientation in the top bit. Local order: cells inner, then owned cells touching the halo, then halo cells by source rank and layer, each group in reverse Cuthill-McKee order with the permutation kept; faces inner, then processor faces, then boundary faces by patch (R1). Every cell, face and node carries a u64 global id in input order, so output and restart are written in global order independent of the partition (R2). Geometry: face area vectors and centroids from a triangle fan about the node average, cell volumes and centroids from tetrahedra on the same fan, so Σ S = 0 holds to round-off for every cell and the volumes add up to the domain volume exactly |
| D24 | 4  | Gradient reconstruction | Weighted least squares on the face-neighbour stencil as the default; Green-Gauss as a cheaper option | Plain Green-Gauss is only zero-order accurate on skewed polyhedral meshes (Syrakos et al. 2017, as in PyAeroFlow D22) **Implementation (stage 2, `cfdcraft-kernels::fv`):** weights 1/|d|², the weight vectors per cell-face entry precomputed on the host, so the kernel is a gather without roots or inversions; fixed-gradient and symmetry faces enter at the normal projection of the face centroid, so linear fields stay exact on every mesh type (error ≤ 6e-15, docs/verification/stage2.md). Green-Gauss not yet implemented |
| D25 | 4  | Convection schemes | Second-order linear upwind with a limiter as the default; first-order upwind for start-up; a blended central scheme as an option (needed for LES, L5). Limiters in milestone 1: Venkatakrishnan and Venkatakrishnan–Wang | Further TVD limiters later through the same interface **Implementation (stage 2):** linear upwind as a deferred correction on implicit upwind; ε² = (K h)³ (Venkatakrishnan) or (K Δφ_domain)² (Venkatakrishnan–Wang), ψ clipped to 1; with K = 0 no new extrema (tested), K > 0 trades small overshoots for smoothness. Blended central scheme not yet implemented |
| D26 | 4  | Diffusion on non-orthogonal meshes | Over-relaxed non-orthogonal correction, applied as a deferred correction with a limiter on strongly skewed cells | **Implementation (stage 2):** over-relaxed decomposition S = Δ + k (Jasak 1996), k·∇φ_f deferred, optional limit as a multiple of the orthogonal part (`solver.nonorth_limit`) |
| D27 | 4  | Rhie–Chow interpolation | Time-step-independent form: steady results do not depend on under-relaxation or time step | |
| D28 | 4  | Time integration | Second-order backward (BDF2) as the default, first order for the first step; Crank–Nicolson as an option; optional adaptive time step controlled by CFL | |
| D29 | 4  | Boundary conditions, milestone 1 | Velocity inlet, mass-flow inlet, static-pressure inlet, total-pressure inlet, static-pressure outlet, total-pressure outlet, no-slip wall (also moving), symmetry, periodic. Turbulence values at inlets as turbulence intensity and length scale, or viscosity ratio. Static and total pressure conditions can be combined, for example static pressure at an inlet with total pressure at an outlet | The compatible combinations are enforced by D30 |
| D30 | 4  | Boundary-condition check and guidance | Before the run, a check of the combination of boundary conditions on the whole domain: a pressure level must be fixed (a pressure boundary or a reference cell), the net mass balance must be possible, over-specified combinations are errors (for example velocity inlet together with pressure inlet on the same flow path), and unused or unsuitable turbulence values are reported. Errors and warnings name the offending boundary and say what to change. The documentation carries a guidance table of sensible combinations and the cases each one suits | Same idea as PyAeroFlow D39; a recipe (D18) shows typical combinations |
| D31 | 4  | Reproducible results | Results are reproducible from run to run on the same hardware and decomposition: gather-based assembly of face contributions with no floating-point atomics, and a fixed summation order in reductions | Mandatory; costs possibly a few percent of speed. Reproducibility across different GPU or rank counts is not promised, because the summation order changes with the decomposition **Stage 2:** face kernel computes per-face coefficients, cell kernel gathers them in stored order into the matrix row; dot products in blocks of 256 combined in a fixed tree, the same blocks on CPU and GPU; tested bitwise for 1–4 threads and run to run |
| D32 | 4  | Mesh partitioning | METIS or Scotch (linked C library, allowed by D9) for domain decomposition | Decision between the two at implementation. **Decided in stage 1: METIS** (vendored, built from source, k-way above 8 parts, recursive bisection up to 8) as the default, plus a built-in recursive coordinate bisection (RCB) that needs no C library and gives the same parts on every machine; Scotch not used. Partition quality (imbalance, cut faces, neighbours) is in the mesh report |
| D33 | 5  | Mesh input, milestone 1 | Gmsh `.msh` version 4 (ASCII and binary), with physical groups as the named boundary-face sets | Further formats (CGNS, other open formats) later through the same reader interface |
| D34 | 5  | Built-in mesh generator for canonical test cases | In code: block-structured hexahedral meshes with varying cell size and boundary-layer grading (2D as one cell layer thick with symmetry side planes, D3, and full 3D); tetrahedral meshes (hexahedra split); polyhedral meshes (dual of a tetrahedral mesh); mixed meshes combining these. Optional random node perturbation and refinement levels for order-of-accuracy studies | Used by the verification cases; same role as PyAeroFlow D42. Not the general mesher of L8 |
| D35 | 5  | Pressure linear solver | Aggregation-based algebraic multigrid as preconditioner of CG, V-cycle, Jacobi or Chebyshev smoothers (GPU friendly). ILU-type smoothers excluded on GPU | In-house (D10), verified against hypre (D20) **Implementation (stage 3):** smoothed aggregation (Vaněk, Mandel, Brezina 1996) with a deterministic host setup, explicit restriction R = Pᵀ (gather, no atomics), dense coarsest level; symmetric V-cycle with Chebyshev (default, degree 2) or damped Jacobi smoothing; flexible CG (Notay 2000). Strength threshold 0.04 instead of Vaněk's 0.08, because the unstructured cells have 10–20 neighbours (docs/verification/stage3.md). Setup on the host in stage 3; GPU setup and reuse of the aggregates across outer iterations with stages 4 and 7 |
| D36 | 5  | Momentum and scalar linear solvers | BiCGStab with Jacobi preconditioning, a few iterations to a loose relative tolerance per outer iteration | Standard practice for segregated solvers **Implementation (stage 2):** right-preconditioned BiCGStab (van der Vorst 1992), default relative tolerance 1e-3 per outer iteration |
| D37 | 5  | Matrix storage and precision in the linear solvers | Matrix stored once in a GPU-friendly sparse format (ELL or sliced ELL; final choice from gate G1). The multigrid hierarchy is stored and applied in float32; the outer Krylov iteration, residuals and solution stay in float64 (D7) | Mixed precision inside the preconditioner halves memory and bandwidth of the hierarchy without changing the converged result, because the outer float64 iteration corrects it **Stage 2:** sliced ELL (SELL-C-σ, C = 32, σ = 1): rows stay in cell order, so the assembly writes them in place; slot 0 the diagonal, then the interior faces of the cell. CSR and SELL give bitwise the same products (G1). The float32 hierarchy follows with the AMG of stage 3 **Stage 3:** every matrix, diagonal and vector of the V-cycle in float32 (including a float32 copy of the fine matrix for the fine-level smoother); the residual is rounded on entry, the correction widened on exit. Same iterations and solution as the float64 hierarchy on all mesh types; R10 allows f32 only in the multigrid kernels |
| D38 | 5  | Parallel communication | Halo exchange overlapped with interior computation; GPU-aware MPI. A multi-node run is a milestone-1 acceptance test | Needs access to a multi-node machine (O12). **Implementation note (stage 1, R3):** each rank has a halo plan per neighbour: a send list of local cells and a contiguous receive range in the halo block (halo ordered by source rank, D23). The exchange packs all fields of one call into one contiguous buffer per neighbour and receives directly into the halo range, one message per neighbour and call. Halo depth is a parameter (default 1, two layers tested). The transport is chosen at run time: host-staged now, GPU-aware MPI from stage 2 on the same plan (no change to the layout). Overlap with interior computation follows from the cell order: inner cells need no halo data. Tested with 1 to 5 parts in process and over MPI with 2, 3 and 4 ranks (`mpi_halo` example, CI) |
| D39 | 6  | Case file | TOML (standard-library style parser in Rust), literature names and established acronyms, every key with a sensible default optional, validation before the run (unknown keys with "did you mean", types, ranges, consistency), resolved configuration with defaults written to the output directory. Same principles as PyAeroFlow D26–D30 | |
| D40 | 6  | Re-reading the case file during a run | The case file can be re-read while a run is going on, on a schedule the user sets: every iteration, every time step, every N iterations or time steps, or never (default). Only settings that are safe to change at run time are applied (output, restart and monitoring settings, stop criteria, relaxation factors, time step); every change is logged, and a change to anything else is reported as not applied. In addition the control files of PyAeroFlow L7 are kept: `.stop` (write a restart and the results, end cleanly), `.result_N` (write results at iteration or step N), `.read_toml_N` (read the case file at iteration or step N) | PyAeroFlow L7 (5) and (6) |
| D41 | 6  | Field output | VTKHDF: one parallel HDF5 file per write (one file per time step or a transient series in one file), read directly by ParaView. Linked HDF5 library (D9) | Per-rank VTU files are not workable at 100+ million cells. **Implementation (stage 1):** VTKHDF 2.4 `PartitionedDataSetCollection` with the block `internal` and one block per patch, an assembly by name, standard cells as VTK cells and all others as `VTK_POLYHEDRON`, cells and points in global-id order. HDF5 is built from source and linked statically (`hdf5-metno`, feature `static`; needs CMake), so no system HDF5 is needed; serial writer in stage 1, the parallel write (one file, each rank its hyperslab in global order) with the solver output of stage 6. Files are checked with VTK's own reader (`tools/check_vtkhdf.py`); ParaView 5.13 or newer reads them (to be confirmed on the maintainer's ParaView) |
| D42 | 6  | Restart | Same HDF5 file layout, written in parallel by all ranks, readable with a different number of ranks or GPUs. A restart continues the run (counters, time, solver state), as PyAeroFlow D56 and D58 | |
| D43 | 6  | Monitoring and user-configurable output | A human-readable `.log` file, CSV histories (residuals, forces, integrals), one CSV of wall quantities per boundary, point probes. Which fields, monitors, probes and wall quantities are written, and how often (per iteration, per time step, every N, at listed steps), is chosen by the user in the case file; defaults are modest. The log header states what is written and how large the output will be | Same idea as PyAeroFlow D59 |
| D44 | 6  | Units | SI on input and output and dimensional inside the solver | PyAeroFlow's non-dimensional form (its D51) is not adopted: reacting flows have many reference scales |
| D45 | 7  | Verification and validation cases | Canonical CFD cases. Verification: plane Couette and Poiseuille flow, 2D Taylor–Green vortex (time and space order), manufactured solutions on hexahedral, tetrahedral, polyhedral and mixed meshes (D34). Laminar benchmarks: lid-driven cavity 2D and 3D (Ghia et al. 1982), flow around a cylinder Re 20 and Re 100 (Schäfer & Turek 1996), backward-facing step (Armaly et al. 1983), Rayleigh–Bénard convection. Turbulent validation: zero-pressure-gradient flat plate (NASA TMR), channel flow against DNS (Lee & Moser 2015), backward-facing step (Driver & Seegmiller 1985), and general-purpose internal flows: fully developed turbulent pipe flow (DNS of El Khoury et al. 2013), asymmetric plane diffuser (Buice & Eaton 1997), U-bend (Monson et al. 1990), a 3D T-junction mixer on a polyhedral mesh (grid convergence and qualitative mixing). External aerodynamics (NACA 0012, NASA hump) is not part of milestone 1 | Each ships with gallery example, case file and recipe (D18). The maintainer's own reference data may be added later |
| D46 | 7  | Versions, tags and changelog | Versions 0.x per stage with a Git tag `v0.<stage>`, 1.0 at the end of milestone 1. The assistant does the version bump and the tag. `CHANGELOG.md` carries a high-level summary on top of the full, detailed changelog for every release | Owner rule of 2026-11-14 |
| D47 | 8  | Design constraints from the discrete-phase (spray) feature, milestone 1 | The face-based mesh (D23) provides what Lagrangian tracking needs: cell adjacency through faces for tracking a particle from cell to cell, a robust point-in-cell search on polyhedra (also after a mesh partition), and halo cells wide enough for particles crossing rank boundaries. Particle data are kept in structure-of-arrays form so that kernels run on GPU. Reproducibility (D31) extends to particles: random numbers come from a counter-based generator keyed by particle identity, step and event, never from a stream that depends on rank or thread order | Costs little now and avoids a redesign when L12 starts; verified in stage 1 by the point-location test on all mesh types. **Stage 1:** cell-to-cell tracking by a tetrahedral-decomposition walk through the faces (`locate.rs`, after Macpherson, Nordin & Weller 2009) with a bin search as fallback at non-convex boundaries; it reports the face where a track leaves the mesh or the rank. The point-location test passes on all eight generated mesh types and on the Gmsh T-junction (tet and poly), 0 points in a wrong cell; halo depth is a parameter (D38). The counter-based generator exists since stage 0 and drives the mesh perturbation |


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

Order of work after milestone 1: L1 (with the compound wall treatment of D13 for heat transfer), L2, L12 (sprays, after energy and species because evaporation needs them), L3, L4, L5; L6 to L11 as the need arises.

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
| L13 | Immersed boundary method | Bodies given as surfaces (for example STL) on a background mesh that is not fitted to them, as an alternative to body-fitted meshes: cut cells or ghost cells on the polyhedral face-based mesh (D23), no change to the solver data structures; moving and deforming bodies together with L4. Needs a wall treatment for the cut boundary that reuses the interface of D13 and boundary conditions on the immersed surface through the D29 set. Literature: Peskin 2002, Mittal & Iaccarino 2005 (survey), Cartesian cut-cell methods. Added by the maintainer 2026-11-14; scheme, accuracy targets and test cases are decided when the item starts |
| L12 | Discrete droplet sprays (Lagrangian) | Droplets are tracked as parcels, each standing for many identical droplets (tracking every droplet is impossible at 100+ million cells). Models in this order of priority: breakup (for example KH-RT), evaporation and heat transfer (with L1 and L2), wall films and wall interaction, collision; injection models, drag and other forces, turbulent dispersion and one-way and two-way coupling to the gas phase come with the first version. Must be parallel: parcels migrate between ranks, the load is balanced when parcels cluster, kernels run on GPU (D47). Complements the Euler–Euler and VOF multiphase of L3 |
| L10 | Fully compressible high-Mach flow (transonic, supersonic) | After low-Mach variable density (D17) |


## 4. Possible later (not decided)

| ID | Feature | Notes |
|----|---------|-------|
| P2 | float32 or mixed-precision GPU mode | Enabled by D7; would let Apple GPUs take part (O6) |
| P3 | Web build (WebAssembly), MCP server | Part of the ArtCraft form; no requirement yet |
| P4 | Spalart–Allmaras | Not in milestone 1 (D12); cheap to add through the transport-equation interface |
| P5 | Low-Re, wall-resolved k-ε variants | D13 |
| P6 | AMD and Intel GPUs | D15 |
| P7 | PANS on ζ-f (partially averaged Navier–Stokes) as a scale-resolving model, with the resolution parameter f_k computed dynamically | Kept for later by the maintainer (2026-11-14), second hybrid route after SST-based DES within L5. Needs ζ-f of D12 and low-dissipation convection (D25). Literature: Basara et al. 2011, Basara 2015, Basara, Pavlovic & Girimaji 2018 and 2019 |
| P8 | Curvature and rotation correction for SST and ζ-f | Motivated by the U-bend case of D45, which neither model captures without it |


## 5. Not planned

| ID | Feature | Notes |
|----|---------|-------|
| X1 | RNG k-ε (Yakhot et al. 1992) | Dropped by the maintainer on 2026-11-14: replaced by k-ε-ζ-f in D12; the standard and realizable k-ε remain for literature data |


## 6. Open requirement steps

| ID | Step | Topic | Current proposal |
|----|------|-------|------------------|
| O13 | 0 | GPU for gate G1 | **Closed (2026-11-14): GPU verdict on a Tesla T4 and CPU verdict on an i9-9900K iMac, both PASS (docs/benchmarks/g1.md).** The rustc NVPTX target supports NVIDIA compute capability 7.0 (Volta) and newer only. The maintainer's PC (i7-4790, 16 GB, GTX 960 with 2 GB, compute capability 5.2, Windows) is older, so the GPU verdict was taken on a free Google Colab T4 with `tools/g1_colab.ipynb` and passed (docs/benchmarks/g1.md). Optional, no longer needed for stage 0: the Windows run on the i7-4790 (`g1_cpu_10m.toml`) and the check of the CUDA reference on the GTX 960 (`g1_gpu_reference_4m.toml`). A GPU of compute capability 7.5 or newer of the maintainer's own is wanted for the local GPU tests of later stages (D22), not for stage 0 |
| O14 | 0 | Nightly toolchain for the GPU route | The single-source GPU kernels need nightly Rust (`abi_ptx`, `stdarch_nvptx`) and may break when nightly changes. Proposal: pin the nightly date in `gpu/ptx-kernels/rust-toolchain.toml` once stage 2 starts and update it deliberately; CI builds the PTX on every pull request, so a break is seen at once. Fallback if the route is dropped: the kernels as CUDA C behind the same entry points (D8's second option), which the G1 harness already carries as `cuda-reference` |
| O15 | 1 | G1 on polyhedral meshes, quiet machine and GPU | D8 asks for G1 on shuffled and polyhedral meshes. Cloud run (shared 4-core machine, 2 million cells, docs/benchmarks/g1.md): all results bitwise equal; single-source bodies 25–42 % faster on cell gather, equal on SpMV-CSR, 4–30 % slower on face flux, also on a hexahedral box of the same size, so a short-kernel effect of the shared machine rather than of polyhedra. Proposal: the maintainer runs `g1_cpu_poly_2m`, `g1_cpu_poly_2m_shuffled` and `g1_cpu_tet_2m` twice on the iMac (recipe 10) and `g1_gpu_poly_4m` on Colab; if face flux stays above the 10 % limit, tune the parallel loop of the CPU backend (packet size, loop form) in stage 2, the kernel bodies stay single-source **CPU part closed (2026-11-14):** the quiet i9-9900K (iMac) passes all six runs, largest time ratio 1.09, bitwise equal (docs/benchmarks/g1.md). **GPU part closed (2026-11-14):** Tesla T4, 4.07 M polyhedra, two runs PASS, largest time ratio 1.08, bitwise equal (docs/benchmarks/g1.md). **Closed.**|
| O12 | 5 | Multi-node acceptance test | D38 makes a multi-node run a milestone-1 test. The maintainer has an old PC with an NVIDIA GPU (model to be determined) and does not want to rent hardware. Proposal: the milestone-1 tests are multi-rank MPI runs on one machine (CPU ranks in CI, D22, and ranks sharing the one GPU locally); a true multi-node run is accepted when a cluster allocation is available (university or research-computing access, a community contribution) and is then added as a test. Until then the multi-node path is designed for and tested with several ranks on one machine, not claimed as verified |

Steps still to come: verification and validation cases, stage plan.


## 7. Stage plan

Milestone 1 is version 1.0: incompressible flow, steady and transient, on
unstructured 3D meshes (polyhedra), with k-ω SST and the k-ε variants, on
NVIDIA GPUs and CPUs, with MPI. Every stage ends with:

- passing tests (CPU tests on GitHub, D22; GPU tests run locally by the
  maintainer, marked `gpu`; long runs marked `slow`);
- a gallery example, case files and a recipe for each new feature (D18);
- updated documentation and `CHANGELOG.md` (D46);
- a Git tag `v0.<stage>`.

| Stage | Content | Tests, acceptance |
|-------|---------|-------------------|
| **0 Foundations and gate G1** | Cargo workspace, CI (D22), documentation skeleton (gallery, cases, recipes). Gate G1: kernel benchmark, a face-based flux loop and a sparse matrix-vector product at 10–50 million cells, hand-written CUDA reference against the single-source Rust route, on the maintainer's GPU and on CPU | G1: single-source route within 10 % of the CUDA reference, or two backends (D8). Language decision D11 confirmed or reopened |
| **1 Mesh** | Polyhedral face-based mesh structure (D23), geometry (face areas and normals, volumes, centroids, non-orthogonality and skewness), Gmsh reader (D33), built-in test mesher (D34), partitioning and halo layout (D32), VTKHDF writer (D41), mesh quality report | Every cell closed (Σ n·dS = 0); total volume equals domain volume; Gmsh files written by Gmsh load; VTKHDF files load in ParaView. **Done 2026-11-14 (release 0.1):** closure below 1e-15 and exact volumes on all generated types and the Gmsh T-junction; Gmsh 4.1 ASCII and binary files from Gmsh load; VTKHDF checked with VTK's reader; point location on all types; golden values for mesh type × perturbation × refinement × 1, 2, 4 ranks (R4); G1 rerun on polyhedral and shuffled meshes (O15) |
| **2 Operators on a scalar equation, GPU** | Least-squares gradients (D24), linear upwind with Venkatakrishnan and Venkatakrishnan–Wang limiters (D25), non-orthogonal diffusion correction (D26), sparse matrix format (D37), BiCGStab (D36), gather-based reproducible assembly (D31), CPU and GPU kernels | Gradients exact for linear fields on all mesh types; limiter creates no new extrema; convection–diffusion with an exact solution, observed order 1.8–2.2 on hexahedral, tetrahedral, polyhedral and mixed meshes; bitwise identical results from run to run; GPU and CPU agree to round-off tolerance  **Implemented 2026-11-14 (`cfdcraft-fv`, `cfdcraft scalar`):** gradient error ≤ 6e-15 on all types; limited face values inside the neighbour range and oblique step in [0, 1] with K = 0; observed order 1.87–2.08 on the finest pair of hex, tet, poly and mixed meshes in 2D and 3D (perturbed); bitwise equal for 1–4 threads; **GPU against CPU on a Tesla T4 (Colab, `tools/stage2_colab.ipynb`): bitwise equal on 912 673 polyhedra with the limiter, on the oblique step and on the sweep; 311 ns per cell and iteration, 21 times faster than the 2-thread Colab CPU.** Order sweeps are CI tests (R7), gallery cases run in CI (R5), ns per cell and iteration per phase in every report (R6), kernel lint in CI (R10) |
| **3 Pressure solver** | Aggregation AMG with CG (D35), mixed-precision hierarchy (D37), hypre reference (D20) | Poisson problem on all mesh types: iteration count nearly independent of mesh size; same converged result as hypre; GPU time per cell compared with hypre and reported  **Implemented 2026-11-14 (`cfdcraft-fv::amg`, `cfdcraft poisson`):** smoothed aggregation (θ = 0.04), Chebyshev V-cycle, flexible CG; 10–13 iterations to 1e-10 on hex, tet, poly and mixed meshes from 512 to 274 625 cells (growth 1.1–1.3; Jacobi-CG 4–7 times); same solution as hypre 3.2.1 BoomerAMG-PCG to 1.6e-9; float32 hierarchy: same iterations and solution as float64; bitwise equal on 1 and 4 threads (docs/verification/stage3.md). GPU on a Tesla T4 (Colab): bitwise equal to the CPU; 6.5 ns per cell and iteration for the in-house AMG against 10.5 for hypre at 912 673 cells, and 10 to 11 iterations from 729 to 274 625 cells (docs/verification/stage3.md) |
| **4 Laminar steady flow** | SIMPLEC, time-step-independent Rhie–Chow (D27, D19), boundary conditions with combination check and guidance (D29, D30), TOML case file and validation (D39), log and CSV output (D43) | Couette: round-off error. Poiseuille: observed order 1.8–2.2 on all mesh types. Lid-driven cavity 2D and 3D, Re 100 and 1000, against Ghia et al.  **Implemented 2026-11-15 (`cfdcraft-flow`, `cfdcraft flow`), CPU (GPU flow loop with stage 5, owner decision 2026-11-14):** Couette error 4e-11 to 2e-10 (at the iteration tolerance) on hex, tet, poly and mixed cells, perturbed; Poiseuille observed order 1.80–2.00 on the finest pair of all four types; cavity against Ghia et al. on 2D hexahedra 0.4 % / 1.0 % (Re 100) and 0.4 % / 1.4 % (Re 1000) of the lid speed, on 3D tetrahedra and polyhedra (slab with symmetry planes) 0.7–2.5 %; converged result independent of α_u to 3e-11 (D27); boundary-condition check with guidance and "did you mean" (D30, D39); bitwise equal on 1 and 4 threads. Periodic boundaries follow with stage 6 (face pairing); on perturbed prisms the Poiseuille error stalls on the finest level (open item, docs/verification/stage4.md) |
| **5 Transient flow, restart, control** | PISO-type algorithm, BDF2 and Crank–Nicolson (D28), adaptive time step, restart (D42), re-reading the case file and control files (D40), configurable output and probes (D43) | Taylor–Green: time order 2, space order 2. Cylinder Re 20 (C_D, C_L, Δp) and Re 100 (C_D,max, C_L,max, Strouhal number) against Schäfer & Turek. Restart continues a run bit-identically on the same decomposition **Part 1 implemented 2026-11-15 (`cfdcraft-flow::transient`, `unsteady`, `restart`, `control`, `monitor`; docs/verification/stage5.md):** PISO with two correctors and optional outer iterations, BDF2 with variable steps and Crank–Nicolson (first step BDF1), extrapolated convecting flux and deferred corrections; Taylor–Green time order 2.21 (BDF2) and 2.00 (Crank–Nicolson), space order 2.09; cylinder Re 20 inside all reference ranges on 98 304 cells; Re 100 on 6 144 cells C_D,max, St and Δp inside, C_L,max 1.9 % low; restart bit-identical (BDF2 adaptive, Crank–Nicolson with outer iterations). The GPU flow loop and the further stage-5 items follow in their own pull requests |
| **6 Turbulence and wall functions** | k-ω SST (2003), standard and realizable k-ε, k-ε-ζ-f, wall functions, the automatic wall treatment and the compound wall treatment (D12, D13); generic transport-equation interface including the elliptic-relaxation equation | Flat plate: C_f and velocity profile against the NASA TMR values; channel flow against DNS (ζ-f: near-wall profiles of k, ε and ζ with y⁺ ≈ 1, and unchanged wall shear for y⁺ from 1 to 100 with the compound treatment); backward-facing step: reattachment length against experiment |
| **7 Parallel performance, internal flows, milestone 1** | Multi-rank MPI with overlapped halo exchange (D38), GPU-aware MPI, strong and weak scaling measurements as far as the hardware allows, the general-purpose internal-flow cases of D45 (pipe, plane diffuser, U-bend, T-junction on a polyhedral mesh), user documentation, tag `v1.0` | Same result on 1, 2 and 4 ranks within solver tolerance; scaling report states what was measured and where, with the multi-node path labelled "designed for" until a cluster allocation exists (O12); pipe-flow wall shear and mean profile against DNS; diffuser separation and pressure recovery against Buice & Eaton; U-bend pressure and velocity profiles against Monson et al.; T-junction grid convergence on tetrahedral, polyhedral and mixed meshes; all cases of D45 in the gallery |

## References

- F. R. Menter, M. Kuntz, R. Langtry, *Ten years of industrial experience with
  the SST turbulence model*, Turbulence, Heat and Mass Transfer 4 (2003)
- B. E. Launder, D. B. Spalding, *The numerical computation of turbulent
  flows*, Comput. Methods Appl. Mech. Eng. 3, 269–289 (1974)
- T.-H. Shih, W. W. Liou, A. Shabbir, Z. Yang, J. Zhu, *A new k-ε eddy
  viscosity model for high Reynolds number turbulent flows*, Computers &
  Fluids 24, 227–238 (1995)
- K. Hanjalić, M. Popovac, M. Hadžiabdić, *A robust near-wall
  elliptic-relaxation eddy-viscosity turbulence model for CFD*, Int. J. Heat
  Fluid Flow 25, 1047–1051 (2004)
- M. Popovac, K. Hanjalić, *Compound wall treatment for RANS computation of
  complex turbulent flows and heat transfer*, FTaC (Springer journal) 78,
  177–202 (2007)
- M. Popovac, *Robust eddy viscosity turbulence modeling with elliptic
  relaxation and compound wall treatment*, Energies 16(9), 3685 (2023)
- K. Hanjalić, B. Launder, *Modelling Turbulence in Engineering and the
  Environment*, Cambridge University Press (2011)
- B. Basara, S. Krajnović, S. Girimaji, Z. Pavlovic, *Near-wall formulation
  of the partially averaged Navier–Stokes turbulence model*, AIAA J. 49(12),
  2627–2636 (2011)
- B. Basara, *PANS method as a computational framework from an industrial
  perspective*, in Progress in Hybrid RANS-LES Modelling, NNFM 130,
  3–17 (2015)
- R. Mittal, G. Iaccarino, *Immersed boundary methods*, Annu. Rev. Fluid
  Mech. 37, 239–261 (2005)
- C. M. Rhie, W. L. Chow, *Numerical study of the turbulent flow past an
  airfoil with trailing edge separation*, AIAA J. 21, 1525–1532 (1983)
- R. I. Issa, *Solution of the implicitly discretised fluid flow equations
  by operator-splitting*, J. Comput. Phys. 62, 40–65 (1986)
- J. P. Van Doormaal, G. D. Raithby, *Enhancements of the SIMPLE method for
  predicting incompressible fluid flows*, Numer. Heat Transfer 7, 147–163
  (1984)
- NASA Turbulence Modeling Resource: https://tmbwg.github.io/turbmodels/
- U. Ghia, K. N. Ghia, C. T. Shin, *High-Re solutions for incompressible flow
  using the Navier–Stokes equations and a multigrid method*, J. Comput. Phys.
  48, 387–411 (1982)
- M. Schäfer, S. Turek, *Benchmark computations of laminar flow around a
  cylinder*, Notes on Numerical Fluid Mechanics 52, 547–566 (1996)
- B. F. Armaly, F. Durst, J. C. F. Pereira, B. Schönung, *Experimental and
  theoretical investigation of backward-facing step flow*, J. Fluid Mech. 127,
  473–496 (1983)
- D. M. Driver, H. L. Seegmiller, *Features of a reattaching turbulent shear
  layer in divergent channel flow*, AIAA J. 23, 163–171 (1985)
- M. Lee, R. D. Moser, *Direct numerical simulation of turbulent channel flow
  up to Re_τ ≈ 5200*, J. Fluid Mech. 774, 395–415 (2015)
- G. K. El Khoury et al., *Direct numerical simulation of turbulent pipe flow
  at moderately high Reynolds numbers*, Flow Turbul. Combust. 91, 475–495
  (2013)
- C. U. Buice, J. K. Eaton, *Experimental investigation of flow through an
  asymmetric plane diffuser*, Stanford Report TSD-107 (1997)
- D. J. Monson, H. L. Seegmiller, P. K. McConnaughey, Y. S. Chen,
  *Comparison of experiment with calculations using curvature-corrected
  zero and two equation turbulence models for a two-dimensional U-duct*,
  AIAA 90-1484 (1990)
- hypre: https://github.com/hypre-space/hypre
- G. B. Macpherson, N. Nordin, H. G. Weller, *Particle tracking in
  unstructured, arbitrary polyhedral meshes for use in CFD and molecular
  dynamics*, Commun. Numer. Methods Eng. 25, 263–273 (2009)
- G. Karypis, V. Kumar, *A fast and high quality multilevel scheme for
  partitioning irregular graphs*, SIAM J. Sci. Comput. 20, 359–392 (1998)
- E. Cuthill, J. McKee, *Reducing the bandwidth of sparse symmetric
  matrices*, Proc. 24th ACM National Conference, 157–172 (1969)

# Changelog

## Unreleased

- **CI: shorter waits with parallel work.** Pull requests run on Linux
  only (about 5 minutes); macOS and Windows join when a pull request touches
  build files, the toolchain, the workflow, the command line, MPI, VTKHDF,
  restart, run control or the console. Pushes to main run only what the
  gate did not (macOS and Windows tests and quick cases, MPI on macOS, the
  long cases, debug tests), skip docs-only merges, and a newer merge
  cancels the older run. The full tier runs nightly, by hand and with the
  label `full-ci`. HDF5 file locking is off in CI (flaky on macOS).
- **Outer iterations in the log.** Every logged time step shows the
  residuals and the linear iterations of each outer iteration under its
  box. New `docs/sayings.md` (*Residuen-gaffen*), linked from the README.
- **Log look (D125).** The log opens with "Foehncraft" in block letters
  with depth, on one row. Every block is a box of asterisks as wide as its
  content and centred on the page, which has no fixed width any more
  (the banner is 141 columns). Summaries have blank rows, a heading and
  dashed rules; a logged time step starts with a box (step, time, dt,
  Courant number, outer iterations). Messages read
  `%FOEHN-I-CASE, text`: severity and the capitalised name of the source
  file that writes them. Only text changes; results, CSV and JSON are
  bitwise the same. Gallery 27 updated.
- **Cell-quality remediation (D50).** Steady runs on tetrahedra with skewed
  faces (Gmsh meshes without `Mesh.OptimizeNetgen = 1`) diverged with the
  least-squares pressure gradient. The solver now measures the face
  skewness at start-up; when faces exceed `solver.skewness_limit` (0.9),
  `solver.remediation = "auto"` (the default) switches to the Gauss
  pressure gradient (new key `solver.pressure_gradient`) and the
  Venkatakrishnan limiter, and reports it (`MESH001`, `MESH002`,
  `flow.log`, `flow.md`, `flow.json`). Meshes without such faces run
  bitwise as before (all catalogue cases). The raw Gmsh T-junction that
  diverged in iteration 87 converges; gallery 26, recipe 36, case
  `flow/pipe_gmsh_raw.toml` with its mesh in the repository.
- **Divergence report.** A pressure multigrid that breaks down in the
  middle of a run (the pressure-correction matrix no longer positive
  definite) is reported as a divergence (`RUN014`, `divergence.md`)
  instead of "the pressure level is not fixed"; `DIV003` names skewed
  faces anywhere in the mesh when the remediation is off.
- The mesh quality report counts the faces with a skewness above 0.9.
- `examples/meshes/tjunction_netgen.geo`: the T-junction with the Netgen
  optimisation; `tjunction.geo` loses a refinement field that Gmsh 4.15
  rejects (it never applied, the mesh is unchanged).
- `tools/plot/residuals_chart.py`: residual histories of several runs in
  one chart.

**Log output (D125–D127).** The log of `foehn flow` is rewritten so that a log passed on alone shows the run: a banner with version and build, framed blocks (input and output files, run control, mesh information with quality, fluid properties, solution settings as a matrix per equation, boundary conditions), tables with repeating header rows for iterations and time steps, summaries (convergence, verification, results per boundary, ledger, run time, warnings) and a closing frame with the exit code. Frames use asterisk rules; `[log] style = "box"` draws them with box characters. Console and `flow.log` are identical. Gallery 27, recipe 37.

**Configurable log table (D125, part 2).** New `[log] columns` (a list) chooses the columns of the iteration and time-step table: the residuals, the linear iterations, values at probes (`probe.<name>.<quantity>`) or in one cell (`cell.<index>.<quantity>`), the largest and smallest value of a quantity with the cell where it occurs (`max.velocity_magnitude`, `min.pressure`), the mass flow of a boundary and the time quantities; `[log] digits` sets the significant digits. Without `[log] columns` the table is the one of before, and no result changes with the columns (CPU and GPU), so `[log]` is not part of the case hash. `run.jsonl` and the machine formats go to 1.1 (additive): a `columns` row describes the table, a `table_row` row carries the values of every table line, and a `section` row carries every framed block of the log with its label, value and unit, so a monitor draws what the terminal shows. Gallery 28, recipe 38, case `flow/cavity_re100_columns.toml`.

**Case-file format 0.6 (D126).** New key `case_format` (absent = 0.5); old keys are mapped on reading and listed in the log; `foehn case upgrade [--write|--check]` rewrites files and keeps their comments.

**Renamed: "patch" is now "boundary" (D126).** Keys `[boundary.<boundary>]` stay as they were; the names that change are the file `patches.csv` → `boundaries.csv`, the field `patches` → `boundaries` in `result.json`, `flow.json` and the ledger, the quantity names `patches.*`, the VTKHDF array `patch` → `boundary`, and the words in the log and the reports. Scripts that read the old names must follow.

**Solver settings per equation (D127).** `[solver.momentum]` and `[solver.pressure]` replace the flat keys: `relaxation_velocity` → `momentum.relaxation`, `relaxation_pressure` → `pressure.relaxation`, `momentum_tolerance` → `momentum.tolerance`, `momentum_max_iter` → `momentum.max_iter`, `pressure_tolerance` → `pressure.tolerance`, `pressure_max_iter` → `pressure.max_iter`, `pressure_preconditioner` → `pressure.preconditioner`, `amg` → `pressure.amg`, and `scheme`, `limiter`, `limiter_k`, `central_blending`, `central_bounding` → `momentum.*`. All 73 example case files are upgraded. Old files still run (format 0.5).

**Settings that change during a run (D128).** New table `[[schedule]]`: entries with `at_iteration` (steady), `at_step` or `at_time` (transient) and `set = { key = value }` change settings from a prepared point of the run, so a whole run is set up before it starts (under-relaxation, tolerances, iteration limits, time step, output and log cadences). Schedulable are the keys a re-read may change (D40) without `run.limits.*` and `run.reread_every`; everything else (mesh, boundaries, schemes, fluid) is refused when the case is read, with the reason (`SCHED001`). Every change is a line in `flow.log`, a row in the "Run control" frame and an `event` row `RUN015` in `run.jsonl` (`key`, `old`, `new`); entries that have fired stay in force over a re-read, a restart applies the entries before its start. Case format stays 0.6 (a new optional table). Without `[[schedule]]` results, files and logs are bitwise unchanged. Gallery 29, recipe 39, case `flow/cavity_re100_schedule.toml`.

CPU results are bitwise unchanged against 0.5 (checked with `tools/check/compare_runs.py --across-versions` on the steady and the transient cavity; only the renamed names differ).

## Summary

| Version | Stage | Highlights |
|---|---|---|
| 0.5 (2026-11-16) | 5 Transient flow, run control, Python | Transient incompressible flow: PISO with variable-step BDF2 and Crank–Nicolson, outer correctors for large time steps, a time step from the integration error, output at exact times, restart files and control files; the whole flow loop on NVIDIA GPUs; adjoint-ready kernels (dual numbers, transposed matrices, checkpointing); the run contract for scripts and agents (`--json`, `--set`, `result.json`, run record, exit codes and signals, case schema, diagnostic codes, estimate, limits, provenance); the conservation ledger and the divergence report; the bounded central scheme and the boundary conditions of D29 (slip, wedge, tables, backflow, convective outlet); derived fields; expressions, the material database and regions; the Python package `foehncraft` with studies; the grid convergence index of every sweep and a generated verification manual; all scripts in Python. **Taylor–Green time order 2.0–2.4 and space order 2.1; Schäfer & Turek 2D-1 inside all reference ranges, C_D extrapolated to 5.5797 ± 0.0004 from three meshes; GPU flow loop bitwise equal to the CPU on a Tesla T4 (722 ns per cell and iteration on 970 225 polyhedra); a case built in Python is bitwise its TOML file.** Commands `foehn materials`, `foehn post derive`, `foehn schema`; gallery 14–25, recipes 21–35 |
| 0.4 (2026-11-15) | 4 Laminar steady flow | First flow solver: steady incompressible laminar flow with SIMPLEC (and SIMPLE) on polyhedral meshes, relaxation-independent Rhie–Chow interpolation (Majumdar 1988), CG + AMG for the pressure correction, boundary conditions of D29 (walls also moving, velocity and mass-flow inlets, static and total pressure inlets and outlets, symmetry) with a check of their combination and guidance (D30), TOML case file with "did you mean" (D39), log, CSV, wall and line output (D43). **Couette exact to the iteration tolerance on all mesh types; Poiseuille observed order 1.80–2.00; lid-driven cavity within 0.4–1.4 % of Ghia et al. at Re 100 and 1000 in 2D and within 2.5 % on 3D tetrahedra and polyhedra.** New crate `cfdcraft-flow`, command `cfdcraft flow`, thirteen case files, gallery 12–13, recipes 18–20. CPU only; the GPU flow loop follows with stage 5 |
| 0.3 (2026-11-14) | 3 Pressure solver | Pressure linear solver: smoothed-aggregation algebraic multigrid as the preconditioner of CG (Vaněk, Mandel, Brezina 1996), symmetric V-cycle with Chebyshev or Jacobi smoothing, float32 hierarchy with float64 Krylov iteration (mixed precision, D37), one kernel source for the CPU and the NVIDIA GPU, optional hypre 3.2.1 reference (BoomerAMG-PCG, built from source). **10–13 iterations from 512 to 274 625 cells on all mesh types; same solution as hypre to 1.6e-9; on a Tesla T4 6.5 ns per cell and iteration against 10.5 for hypre at 912 673 cells, all GPU results bitwise equal to the CPU.** New command `cfdcraft poisson`, four Poisson case files, gallery 10–11, recipes 15–17 |
| 0.2 (2026-11-14) | 2 Operators on a scalar equation | Finite-volume operators on a scalar transport equation, one source for the CPU and the NVIDIA GPU: least-squares gradients (exact for linear fields, error 6e-15 on all mesh types), linear upwind with Venkatakrishnan and Venkatakrishnan–Wang limiters (no new extrema at K = 0), over-relaxed non-orthogonal diffusion, reproducible gather-based assembly into sliced ELL, BiCGStab. Observed order 1.87–2.08 on hex, tet, poly and mixed meshes in 2D and 3D. **GPU and CPU results bitwise equal on a Tesla T4** (912 673 polyhedra, 311 ns per cell and iteration). New command `cfdcraft scalar`, seven case files, gallery 8–9, recipes 12–14 |
| 0.1 (2026-11-14) | 1 Mesh | Polyhedral face-based mesh with exact geometry (every cell closed to 1e-15, volumes exact), Gmsh 4.1 reader, built-in test mesher (hex, tet, poly, mixed; 2D and 3D; graded, refined, perturbed), quality report, VTKHDF output for ParaView, METIS and RCB partitioning with a solver-ready halo layout and multi-field halo exchange (tested over MPI on 2–4 ranks), particle point location on all mesh types, golden-value test matrix, `cfdcraft mesh` command. G1 rerun on polyhedral and shuffled meshes: bitwise equal; quiet-machine verdict open (O15) |
| 0.0 (2026-11-14) | 0 Foundations and gate G1 | Cargo workspace, CI on Linux, macOS and Windows with MPI runs on 4 ranks, the G1 kernel benchmark. **G1 GPU verdict passed on a Tesla T4: single-source Rust kernels within 1.0 % of hand-written CUDA C, bitwise equal results**; **CPU verdict passed on an i9-9900K: single-source bodies within 1.8 % of hand-written loops, in two runs** |

Release 0.5 closes stage 5: the acceptance checks of the stage plan pass
(Taylor–Green time and space order 2, also at Courant number 9 with outer
correctors; the cylinder benchmark of Schäfer & Turek, at Re 20 inside all
reference ranges, at Re 100 with C_L,max 0.5 % below its range on 24 576
cells; restart bit-identical; expressions from TOML and from Python
give the same run; the energy budget of the cavity closes; the GPU flow
loop gives bitwise the CPU results). The generated
[verification manual](docs/verification/manual.md) lists every case with
its pass bands. The Cargo version of the crates is now 0.5.0; the release
is the Git tag `v0.5` (D46).

Release 0.4 closes stage 4: all acceptance tests of the stage plan pass
(Couette at round-off, that is at the iteration tolerance; Poiseuille order
1.8–2.2 on all mesh types; lid-driven cavity in 2D and 3D at Re 100 and 1000
against Ghia et al.). The Cargo version of the crates is now 0.4.0; the
release is the Git tag `v0.4` (D46).

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

## 0.5 (2026-11-16) – stage 5: transient laminar flow, run control and Python

**Highlights.**

- **Transient flow.** PISO with variable-step BDF2 and Crank–Nicolson,
  outer correctors (Courant number 9 at time order 2), an adaptive time
  step from the integration error with exact step rejection, output at
  exact simulated times, potential-flow start, restart files that
  continue a run bit for bit, control files and re-reading the case file
  during a run, probes and force histories. Verified on the Taylor–Green
  vortex (time order 2.0–2.4, space order 2.1) and the cylinder benchmark
  of Schäfer & Turek (Re 20 and Re 100).
- **GPU.** The whole flow loop runs on NVIDIA GPUs from the same kernel
  source as the CPU, bitwise equal to the CPU on a Tesla T4.
- **Adjoint-ready kernels.** Dual-number derivatives of every kernel,
  transposed matrices as gathers, binomial checkpointing (D59).
- **Run contract and health.** `--json`, `--set`, `result.json`, the run
  record, exit codes and signals, the case schema, diagnostic codes, the
  estimate before a run, limits, provenance; the mass and kinetic-energy
  ledger and a divergence report with ranked causes and fixes.
- **Numerics and boundary conditions.** The bounded central scheme, slip
  walls, wedges, inlet profiles from tables, backflow treatment and
  convective outlets; derived fields (vorticity, Q, λ₂, wall shear stress,
  y⁺) with the output and from restart files.
- **Expressions and materials.** Boundary values, body forces and the
  viscosity as expressions of x, y, z and t; a material database with
  sources; regions.
- **Python.** The package `foehncraft`: cases, expression objects, runs in
  process with callbacks, results, materials and parametric studies; all
  scripts of the repository are Python now.
- **Verification.** The grid convergence index of every sweep of three or
  more meshes or time steps, and the generated verification manual.

The parts in detail, newest first:

**Scripts in Python (stage 5 part 6d, D61 (2)).** The shell scripts are
ported to Python and `tools/` is split into groups that run as modules
from the repository root: `tools/check/` (`ci_cases`, `compare_runs`,
`backend`, `c_free`, `kernels`, `vtkhdf`), `tools/bench/` (`build_hypre`
and the cloud GPU drivers with their shared steps in `cloud`),
`tools/mesh/`, `tools/plot/` (the charts), `tools/release/`
(`verification_manual`) and `tools/notebooks/` (the Colab notebooks).
`tools/ci_cases.sh SET` is now `python3 -m tools.check.ci_cases SET`,
`gpu/ptx-kernels/build.sh` is `python3 gpu/ptx-kernels/build.py`; CI,
docs, examples and error messages are updated, and CI lints the scripts
(ruff, syntax and names). See `tools/README.md`.

**Discretisation uncertainty and the verification manual (stage 5 part
6c, D67, D68).** A sweep of three or more meshes or time steps reports
for every monitored quantity (patch mass flows, mean pressures and forces,
force coefficients, benchmark values, viscous and numerical dissipation)
the observed order by the fixed-point iteration of Celik et al. (2008),
the extrapolated value, the band φ1 ± U, the grid convergence index and a
verdict (monotonic, oscillatory, divergent, unchanged): a table in
`flow.md` and `uncertainty` in `flow.json` and `result.json` (an
addition; `Result.uncertainty` in Python). Level sweeps no longer need an
exact solution or a benchmark; `expected_order` without `[analytic]` is
now refused (CASE005) instead of being ignored. New case
`cylinder_re20_gci.toml` (the benchmark 2D-1 on 6 144 to 98 304 cells,
`long` set): C_D extrapolates to 5.5797 ± 0.0004 inside the reference
range. New script `tools/release/verification_manual.py` writes
`docs/verification/manual.md`, `results.json` and a chart of the observed
orders from the case files and a catalogue run. Gallery 25, recipes 34
and 35 on new Verification pages; `tools/gci_chart.py`. Runs are bitwise
unchanged.

**The Python package (stage 5 part 6b, D61, D62).** `python/` holds the
package `foehncraft` with a native module (`crates/foehncraft-py`, PyO3,
built with maturin, abi3 for Python 3.10+): `Case` (read, change with
dotted keys, write; the TOML round trip holds for every example case;
relative paths kept when saved elsewhere), expression objects with Python
operators that render the strings of the case file (`foehncraft.expr`,
evaluated with the solver's own compiler), `run` in process or through
`foehn` with callbacks on every iteration or time step that can stop the
run, `Result` (summary, CSV tables, VTKHDF fields), the material database
and material builders, and `foehncraft.study` (grid, list, Latin
hypercube, Sobol; one case file per variant, `study.csv`, resume). The
engine gains `engine::Observer` (sees every row of the run record) and
`RunOptions::observer`; runs are bitwise unchanged. New CI job "Python
package" (ruff, pytest). New scripts `examples/python/poiseuille_inflow.py`
and `examples/studies/couette_viscosity_study.py`; gallery 24, recipe 33
on the new Python pages, case catalogue page "Python scripts".

**Expressions, materials and regions (stage 5 part 6a, D53, D102, D112,
D121, D124 W3/W4).** Boundary values (velocity, normal velocity, pressure,
total pressure, mass flow), the body force and the viscosity may be
strings of a small expression language in x, y, z and t, with
`[constants]` and piecewise linear `[tables]` (inline or from a file); the
strings are checked when the case is read (EXPR001–EXPR005, with the
place and a "did you mean") and compiled to a stack-machine program that
the solver evaluates with new in-house elementary functions
(`foehncraft_core::math`, D93 (6)), the same on every platform. Fluid
properties from a built-in material database (air, water, engine oil,
ethylene glycol, glycerol, each with its source) with `[fluid] material`,
`temperature` and `operating_pressure`; `[materials.<name>]` with `base`
and the forms constant, polynomial, Sutherland, power law, table (linear
or monotone cubic) and expression; the new command `foehn materials`
lists and evaluates them (also `--json`). `[regions.<name>]` assigns cell
zones to materials; non-fluid regions and unknown zones are refused
(REG001, REG002); different viscosities per region give a viscosity per
face. Walls take a `mass_flux` (blowing and suction), checked for closed
volumes (BC026), and expose per-face exchange coefficients
(`derived::wall_exchange`). The resolved property program is written into
`case.resolved.toml`; `properties_hash` and `materials_database` join the
provenance. The random-number key gains a stream word (`rng::keyed`;
meshes unchanged). `foehn flow --check` also checks the regions against
the mesh. New exact solutions `couette-viscosity` and
`couette-transpiration`; new cases `flow/couette_viscosity.toml`,
`couette_transpiration.toml`, `poiseuille_water.toml`,
`poiseuille_expression.toml`, `taylor_green_expression.toml` and
`channel_ramp.toml`; gallery 23, recipes 31 and 32 on the new
"Expressions and materials" pages; `tools/property_chart.py`. Cases
without the new keys run bitwise as before.

**Derived fields of stage 5 (part 5c, D72).** `[output] derived` adds
fields computed from the solution to every `fields*.vtkhdf`: `vorticity`
(components and magnitude), `q-criterion`, `lambda2`, `enstrophy`,
`strain-rate`, `pressure-coefficient` (reference state in
`[output.reference]`: velocity, pressure, density), `total-pressure`, and
on the patch blocks of the file `wall-shear-stress` and `y-plus`. The new
command `foehn post derive CASE` writes the same file afterwards from a
restart file (`--restart`, default `latest`; `--fields`, `--file`), equal
byte for byte to the file of the run. The VTKHDF writer gains boundary
face fields on the patch blocks; the `units` attribute lists the derived
datasets; the estimate counts them; `output.derived` and
`output.reference` may change during a run. `tools/render_mesh.py` gains
`--diverging`. New case `flow/couette_derived.toml` (exact on any mesh);
gallery 22, recipe 30 on the new Post-processing pages. Without `derived`
every output is as before.

**Boundary conditions of stage 5 (part 5b, D29).** New condition types
`slip-wall`, `wedge` and `convective-outlet` (∂u/∂t + U_c ∂u/∂n = 0 in
transient runs, the outlet flow balanced with the inflow, U_c by default
the mean outflow speed); inlet profiles from a table (`profile = "table"`,
`table`, `table_coordinate` x, y, z or the radius about an axis, the
table's sha256 in the case hash); backflow treatment at static and total
pressure outlets (`backflow = "total-pressure"`, the new default,
`"zero-gradient"`, the treatment before, or `"suppress"`, which closes the
faces where the flow would enter). Axisymmetric wedge meshes from any 2D
mesh of the box generator (`mesh.generate.wedge`, the axis without a
patch), and the exact solutions `uniform` and `pipe-poiseuille`. Seven new
diagnostic codes BC019–BC025. The restart file of a run with a
convective outlet gains the datasets `ub_old_x|y|z`. New cases
`flow/channel_slip.toml`, `poiseuille_table.toml` (with
`poiseuille_profile.csv`), `pipe_wedge.toml`, `cylinder_re100_short.toml`
and `cylinder_re100_convective.toml`; gallery 21, recipe 29. Runs whose
outlets see no backflow are bitwise as before; the Poiseuille cases on
perturbed prisms, polygons and mixed cells, where a few outlet faces see
inflow in the first iterations, change by 2e-6 relative
(`backflow = "zero-gradient"` gives the earlier results bitwise).

**Wall film requirements (D123, D124).** Requirements for the wall film of the Lagrangian sprays (L12): an Eulerian thin film on the wall faces with parcels for splash, rebound, separation and stripping; impingement regimes (Bai–Gosman, Kuhnke), capillarity and dry patches in the first version, CHT coupling of the heated film, ledger and reproducibility rules, staging F1 to F5 and validation cases; Milestone 1 keeps a surface mesh and surface halo of the wall faces, per-face wall boundary values including a wall mass flux, and gravity as a case vector. Documentation only.

**Numerics of stage 5 (part 5a).** A bounded central convection scheme
(`solver.scheme = "central"`, D25) with the bound β and a blend with linear
upwind; the time step from the error estimate of BDF2 and Crank–Nicolson
(`time.error_tolerance`, D65), with steps above the tolerance taken back
exactly and repeated shorter (`solves[].transient.rejected` in
`result.json`, an `error` value in the step rows of `run.jsonl`); output in
simulated time that the steps land on (`output.fields_interval`,
`fields_times`, `monitor_interval`, `run.restart_interval`,
`restart_times`, D66); and a potential-flow start
(`initial.potential_flow`, D52). The restart file gains the attribute
`dt_prev2` and, with the error estimate, the datasets `u_old2_x|y|z`;
older files still read. New cases `flow/taylor_green_error.toml`,
`flow/cylinder_potential.toml` and the variant `bounded-central` of
`scalar/oblique_step_2d.toml`; gallery 20, recipe 28,
`tools/time_step_chart.py`. Runs without the new keys are bitwise as
before.

**Material property requirements (D94–D103, L23–L25).** Requirements for a
material database and the handling of properties: database, case materials
and a resolved, hashed property program; property forms (constant,
polynomial, NASA-7 and NASA-9, Sutherland, tables with linear or monotone
cubic interpolation, polynomial c_p expressions, mixing rules, kinetic
theory, equations of state); checks; GPU evaluation; SI only; no property
changes at run time; a built-in starter set with IAPWS water and steam;
importers, table generators and a generic interface for external property
sources. Implementation starts in stage 5 PR 6 (μ and ρ) and continues with L1.

**Renamed from CFDCraft (and the short-lived Ennecraft) to Foehncraft**, after
the foehn wind of Innsbruck (see the README), so that the name cannot be
confused with other projects. Crates are now `foehncraft-*`, the command is `foehn`, the
environment variables are `FOEHNCRAFT_*`, and the repository is
`chiefenne/Foehncraft`. Entries below keep the names of their time.

**Faster CI.** Pull requests get a fast gate (about 5 minutes instead of
about 25): a `ci` build profile (release optimisation without LTO, results
bitwise equal), the case catalogue split over parallel jobs
(`tools/ci_cases.sh`), and a VTKHDF check that no longer computes unused
volumes of field files (11.5 minutes to seconds). Push to main, a nightly run
at 02:07 Vienna time and pull requests with the label `full-ci` run the full
tier: `release` profile, debug tests, every case including the cylinder at
Re 100, and the slow tests. One summary check, "CI OK". Main warms the `ci`-profile cache that pull
requests read (without it every PR compiled its dependencies from scratch).

**Run contract for scripts and agents (stage 5, part 4a).** `foehn flow`
can be driven by a program: `--set key=value` overrides (TOML values),
global `--json` (one JSON object on standard output, human text on
standard error) and `--quiet`, exit codes 0 finished, 1 error, 2 failed
checks, 3 stopped; SIGUSR1 writes a restart file, SIGTERM, SIGINT and
SIGUSR2 stop cleanly, `--restart latest` continues bitwise; budgets in
`[run.limits]` (wall time, steps or iterations, output volume). Every run
writes `result.json` (status, overrides, provenance with commit and case
hash, diagnostics, results, SHA-256 of the outputs; byte-identical between
identical runs) and the run record `run.jsonl`. 41 stable diagnostic
codes with key paths, locations and `--set` fixes; `foehn schema` prints
the JSON Schema of the case file, `docs/keys.md` is generated from it;
`--check` estimates memory, output and wall time; `--read-root` and
`--write-root` confine the files; CSV, VTKHDF and restart files describe
themselves; HDF5 files have no timestamps. The CLI is a thin layer on
`foehncraft_flow::engine`; `foehncraft-flow` without default features is
checked C-free in CI. Formats and codes: `docs/interfaces.md`; recipe 26,
gallery 18, case `cavity_re100_budget.toml`. No result changes: the
catalogue gives bitwise the fields of main.

**Health checks: conservation ledger and divergence report (stage 5,
part 4b).** Every flow run books mass and kinetic energy (D88): storage,
the terms of each patch (viscous work, convected energy, pressure work,
mass flow), the body force, the viscous dissipation of the discrete
operator, and the remainder, the numerical dissipation, which the second
law wants positive (D90 part 4), split into its pressure part (the
momentum interpolation) and the part of the convection and time
schemes. Steady runs give the
table every `[balance] every` iterations and at the end; transient runs a
compact line with every progress line and a time-integrated ledger of the
run (trapezoidal rule). `balance.csv`, `balance.svg`, `balance` rows in
`run.jsonl`, `solves[].balance` in `result.json`, tolerances with the
codes `BAL001`–`BAL003` (the new key `[balance] tolerance_energy` for
the negative numerical dissipation) and an optional stop. Measured: on the lid-driven
cavity the numerical dissipation is positive, about 1.5 % of the lid
power and falling with the mesh; on the Taylor–Green vortex it falls with
order 2.5, its pressure part with order 3. A run that diverges (D63) now writes its last valid state and
restart file, a VTKHDF snapshot of the failing cells and their
neighbours, and a report (`divergence.md`, `divergence.json`) with the
likely causes ranked by their evidence and a `--set` fix for each, and
ends with `RUN014` (exit code 1). A linear-solver breakdown counts as a
divergence, and a NaN in the continuity residual no longer slips past the
check. The restart file can hold point coordinates and cell volumes
(D91). Recipe 27, gallery 19, case `cavity_re100_energy.toml`,
`tools/balance_chart.py`, `tools/energy_budget_chart.py`. No result
changes: the catalogue gives bitwise the fields of part 4a.

### Adjoint-ready kernels (part 3)

Preparation for a discrete adjoint solver (D59; the solver itself is L17).
No result changes: all flow cases give bitwise the files of part 2.

- **Dual numbers:** `foehncraft_core::Dual<N>` implements `Real`, so every
  kernel body also runs on dual numbers and returns the exact derivatives
  of its results, from the same source. Switches (upwind direction,
  limiter branch, clipping) follow the values, i.e. they are frozen at the
  linearisation state; the limiter array ψ can be held constant.
  **Jacobians of 11 face, cell and boundary kernels agree with central
  differences in every input array (320 580 derivatives, largest relative
  difference 5e-7).**
- **Shape derivatives:** the geometry of `foehncraft-mesh` (face areas and
  centroids, cell volumes and centroids) is generic over the coordinate
  type, so derivatives with respect to the node coordinates come from the
  same code (165 024 derivatives against central differences, 2.8e-7).
- **Transposed matrices:** a transpose map of the matrix storage gives Aᵀ
  in the storage of A and the product Aᵀx without storing it, both gathers
  without atomics (`foehncraft_kernels::adjoint`; flow-device calls
  `transpose`, `matrix_product`, `transposed_product`; PTX entries
  `flow_transpose`, `flow_matrix_product`, `flow_transposed_product`).
  ⟨Av, w⟩ = ⟨v, Aᵀw⟩ to 1e-14 on the matrices of a SIMPLEC solve; CG with
  the multigrid and BiCGStab solve the transposed systems.
- **Binomial checkpointing** (`foehncraft_flow::checkpoint`, Griewank &
  Walther 2000): an optimal schedule (checked against a dynamic program),
  a driver for any primal run and the in-memory snapshot of the flow
  solver; the reverse sweep of a BDF2 run meets every level bitwise.
- Gallery 17, recipe 25, `tools/checkpoint_chart.py`, report
  [stage 5, part 3](docs/verification/stage5.md#adjoint-ready-kernels-part-3).

### GPU flow loop (part 2)

- **The flow loop on the GPU:** `run.backend = "cuda"` or
  `foehn flow --backend cuda` runs every flow case, steady and transient,
  on an NVIDIA GPU (compute capability 7.0 or newer), with all fields,
  matrices and old levels in GPU memory; results bitwise equal to the CPU
  (D31). **Tesla T4 (Colab): all 8 cases bitwise equal to the CPU; 722 ns
  per cell and SIMPLEC iteration on 970 225 polyhedra (7.1 times the two
  CPU threads of the VM), 1 049 ns per cell and PISO step on 98 304 cells
  (9.2 times).** Gallery 16.
- **One source down to the argument lists:** the 27 operations of the flow
  loop are entries of `foehncraft_kernels::entry` (a `#[repr(C)]` argument
  block and a body per index). The CPU loops and the PTX entries `flow_*`
  run the same body with the same block; the flow device is written once
  over an executor (`device::kernel`).
- **CUDA in one crate (D60):** new crate `foehncraft-cuda` holds all CUDA
  code (flow executor, the scalar and multigrid backend, the G1 backends
  and the CUDA C reference, moved from `foehncraft-fv` and
  `foehncraft-bench`); the solver crates find it through registered
  factories. CI lint `tools/check_backend.py`. The feature is now
  `foehncraft-cuda/cuda` (the CLI's `cuda` feature is unchanged).
- New tools `tools/compare_runs.py` (bitwise comparison of two runs),
  `tools/stage5_cloud.py` and `tools/stage5_colab.ipynb` (CPU against GPU
  on eight cases up to one million polyhedra); new cases
  `cavity_poly3d_gpu.toml`, `cylinder_re100_gpu.toml`; recipe 24.
- The CPU results of all flow cases are bitwise those of part 1. A missing
  NVIDIA driver now gives a message instead of a panic.

### Transient laminar flow (part 1)

Stage 5 comes in several pull requests and is released as 0.5 when it is
complete. The first part, the transient core and its benchmarks
([stage 5](docs/verification/stage5.md)):

- **PISO (Issa 1986)** with variable-step BDF2, Crank–Nicolson (first step
  BDF1, the start of Rannacher 1984) and BDF1; two correctors with the
  neighbour correction between them; optional outer iterations per step
  (PIMPLE-type) with a residual tolerance; the convecting flux and the
  deferred corrections extrapolated to the new time, which keeps plain
  PISO second order. A case is transient when it has a `[time]` table.
- **Adaptive time step** from the largest Courant number, with a growth
  limit inside the zero-stability range of BDF2.
- **Rhie–Chow time-step term (D27)** whose steady state is exactly the
  steady solution, for any Δt (cavity: equal to SIMPLEC to 2e-12).
- **Restart files (D42, HDF5):** a continued run is bit-identical to the
  run without the interruption. **Control files and re-reading (D40):**
  `.stop`, `.result_N`, `.read_toml_N`, `run.reread_every`. **Monitors
  (D43):** `history.csv`, `forces.csv`, `probes.csv`, field series,
  force coefficients.
- **Built-in cylinder mesh** `[mesh.cylinder]` (O-grid and wake block) and
  the parabolic inlet profile; the Schäfer & Turek benchmarks 2D-1 and
  2D-2 with their evaluation (C_D, C_L, St, Δp).
- **Results:** Taylor–Green time order 2.21 (BDF2) and 2.00
  (Crank–Nicolson) on the finest pair, space order 2.09, order 2 at
  Courant numbers up to 9 with outer iterations; cylinder Re 20 inside all
  reference ranges on 98 304 cells; Re 100 on 6 144 cells with C_D,max,
  St and Δp inside the ranges, on 24 576 cells C_L,max 0.5 % below its
  range and the others inside.
- `cfdcraft flow --restart FILE --end-time T`; seven new case files,
  gallery 14–15, recipes 21–23. The Release workflow now refreshes the
  cached README badges (recipe 6).

## 0.4 (2026-11-15) – stage 4: laminar steady flow

New crate `cfdcraft-flow`, kernel bodies `cfdcraft-kernels::flow` (one
source for the CPU and the GPU) and the command `cfdcraft flow`:

- **SIMPLEC and SIMPLE (D19).** Segregated outer iterations on the
  collocated, cell-centred polyhedral grid: momentum with implicit upwind
  and deferred linear-upwind convection, orthogonal diffusion implicit and
  the non-orthogonal and skewness corrections deferred, implicit
  under-relaxation; pressure correction with the SIMPLEC factor
  D^C = V/(a_P/α − Σ|a_N|), solved with CG and the stage-3 multigrid
  (hierarchy rebuilt every 10 outer iterations); reference cell when no
  boundary fixes the pressure.
- **Rhie–Chow interpolation (D27)** with the relaxation term of Majumdar
  (1988): the converged result does not depend on the under-relaxation
  (lid force equal to 3e-11 for α_u = 0.7 and 0.9).
- **Boundary conditions (D29):** wall (also moving), velocity inlet (vector,
  normal speed or the exact profile), mass-flow inlet, static- and
  total-pressure inlet and outlet (freely combined), symmetry (implicit per
  component). Periodic boundaries follow with stage 6.
- **Boundary-condition check (D30)** before every run and with
  `cfdcraft flow --check`: missing and unknown patches with "did you mean",
  net mass flow without a pressure boundary, over-specified flow paths,
  wall normal velocity, unused turbulence values, non-planar symmetry,
  reference-cell note; the guidance table of sensible combinations in the
  case catalogue.
- **Case file (D39):** `[fluid]`, `[flow]`, `[initial]`, `[pressure]`,
  `[boundary.<patch>]`, `[solver]`, `[analytic]`, `[verification]`, `[run]`,
  `[output]`; unknown keys and values with "did you mean"; the resolved case
  with all defaults in `case.resolved.toml`.
- **Output (D43):** `flow.log` with a header of what is written and how
  large, `residuals.csv`, `patches.csv` (mass flow, mean pressure, force),
  `wall_<patch>.csv` (pressure and wall shear stress per face, in global
  face order), `line_<name>.csv`, `fields.vtkhdf`, `flow.md`, `flow.json`
  and SVG charts (residuals, convergence, centre lines against the
  reference).
- **Verification:** exact solutions (Couette, Poiseuille) and the
  centre-line data of Ghia et al. (1982) for Re 100 and 1000; thirteen case
  files (CI runs all of them on Linux and the quick ones on macOS and
  Windows); tests for Couette, Poiseuille order, the cavity, relaxation
  independence, quick convergence with symmetry planes, the check, "did you
  mean" and bitwise equality on 1 and 4 threads.
- `Backend::set_matrix` (new matrix values on the same pattern) on the CPU
  and CUDA backends; `tools/render_mesh.py --no-edges`.
- Documentation: [stage 4 verification](docs/verification/stage4.md),
  gallery 12–13, recipes 18–20, flow cases, keys, boundary conditions and
  guidance in the case catalogue.

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

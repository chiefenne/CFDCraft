# Run control, health and adjoint readiness

Part of the [Example gallery](../gallery.md): adjoint-ready kernels, runs that a program drives, balances and the divergence report.

| # | Example | Stage | Feature |
|---|---|---|---|
| 17 | [Adjoint-ready kernels: derivatives and checkpointing](#adjoint-ready) | 5 | dual numbers through the kernel bodies, transposed matrices, binomial checkpointing |
| 18 | [A run that a program can drive](#run-contract) | 5 | `--json`, exit codes, `--set`, limits, signals, `--restart latest`, diagnostic codes, estimate before the run, byte-identical files |
| 19 | [Balances and the divergence report](#health) | 5 | mass and kinetic-energy ledger, numerical dissipation of the schemes (second law), divergence report with ranked causes |
| 26 | [A log that can be circulated](#log) | 5 | frames, blocks of settings, tables with headers, summaries, closing frame with exit code (D125) |

---

<a id="adjoint-ready"></a>
## 17. Adjoint-ready kernels: derivatives and checkpointing

![Primal steps per time step of a reverse sweep against the snapshots in memory](../figures/gallery/checkpoint_steps.svg)

Preparation for a discrete adjoint solver (D59, L17). The kernel bodies are
generic over the floating-point type, so they also run on dual numbers and
return the exact derivatives of their results; matrices can be applied
transposed as a gather; and a transient run can be reversed with few
states in memory.

The chart shows the cost of reversing a transient run with binomial
checkpointing (Griewank & Walther 2000): with 20 snapshots, the reverse
sweep of 10 000 time steps runs each step 3.8 times on average and never
more than 4 times; with 40 snapshots, 100 000 steps cost 3.9 primal steps
each (values above 12 drawn at 12). The schedule is optimal: its count of
recomputed steps equals that of a dynamic program over all positions of
the snapshots.

| Check | Result |
|---|---|
| Jacobians of 11 face, cell and boundary kernels from dual numbers, every input array, against central differences | 320 580 derivatives, largest relative difference 5.0e-7 |
| Shape derivatives of areas, centroids and volumes from the node coordinates | 165 024 derivatives, 2.8e-7; total volume unchanged by interior nodes |
| ⟨Av, w⟩ = ⟨v, Aᵀw⟩ on the momentum and pressure matrices of a SIMPLEC solve | to 10⁻¹⁴; Aᵀ without storing it equals the stored Aᵀ bit for bit, on 1 and 4 threads |
| CG + multigrid and BiCGStab on Aᵀ | converged to 10⁻¹⁰ |
| Reverse sweep of 13 BDF2 steps with 2 snapshots | every level bitwise equal to the forward run, 24 primal steps (the optimum) |

| | |
|---|---|
| Code | `foehncraft_core::Dual`, `foehncraft_kernels::adjoint`, `foehncraft_flow::checkpoint`, `foehncraft_mesh::geometry` (generic over `Coord`) |
| Tests | `crates/foehncraft-flow/tests/adjoint.rs`, `crates/foehncraft-flow/src/adjoint_tests.rs`, `checkpoint.rs`, `geometry.rs` |
| Report | [stage 5, part 3](../verification/stage5.md#adjoint-ready-kernels-part-3) |

```
cargo test --release -p foehncraft-flow --test adjoint -- --nocapture   # Jacobians, one line per kernel
cargo test --release -p foehncraft-flow --lib -- adjoint_tests checkpoint   # transposed matrices, checkpointing
python3 -m tools.plot.checkpoint_chart checkpoint_steps.svg
```

[↑ back to top](#run-control-health-and-adjoint-readiness)

---

<a id="run-contract"></a>
## 18. A run that a program can drive

![Memory and file sizes estimated by foehn flow --check against eleven runs](../figures/gallery/estimate.svg)

A parametric study, a batch job or an agent drives `foehn` without
reading its log: it sets keys with `--set`, reads one JSON object, acts on
the exit code and continues a stopped run with `--restart latest`. Before
it starts a run, `--check` tells it what the run will cost. The chart
shows that estimate against eleven runs of the catalogue, from 128 to
28 227 cells, hexahedra to polyhedra: peak memory within 6 %, field and
restart files within 9 %. The memory model is a base of 14 MB plus
1960 bytes per interior face and 710 per boundary face.

```
$ foehn --json flow cavity_re100_budget.toml --set run.limits.max_iterations=300; echo "exit $?"
{
  "format": "foehncraft.result",
  "format_version": "1.0",
  "status": "stopped",
  "exit_code": 3,
  "stop": { "code": "RUN004", "kind": "iterations", "message": "stopped at the limit of 300 outer iterations ([run.limits] max_iterations)" },
  ...
exit 3
$ foehn --json flow cavity_re100_budget.toml --restart latest; echo "exit $?"
  ...
  "status": "finished",
exit 0
```

| Check | Result |
|---|---|
| Stop at a limit or on SIGTERM, then `--restart latest` | fields and restart file bitwise equal to the run without the stop (steady and transient) |
| Two identical runs | `result.json`, `fields.vtkhdf`, `restart.h5` and every CSV byte for byte the same; only the log, the reports and `run.jsonl` (timings) differ |
| Exit codes | 0 finished, 1 error (as a JSON object with code and key), 2 failed check (`RUN006`, `RUN007`), 3 stopped (`RUN001`–`RUN005`) |
| A typo in the case file | `CASE002` with the key `mesh.generate.cels`, line 4, column 1 and "did you mean `cells`" |
| The 20 flow cases of the catalogue, this build against main | every field, CSV and JSON value bitwise equal; the resolved case gains the empty `[run.limits]` table |

| | |
|---|---|
| Code | `foehncraft_flow::engine` (run, check, estimate, limits, result, run record), `diag` (codes), `schema` (JSON Schema), `foehncraft-cli` |
| Tests | `crates/foehncraft-cli/tests/run_contract.rs`, unit tests in `engine.rs`, `diag.rs`, `schema.rs` |
| Formats | [interfaces.md](../interfaces.md), [keys.md](../keys.md) |
| Case | `examples/cases/flow/cavity_re100_budget.toml` |
| Recipe | [26](../recipes/runs.md#run-contract) |

```
cargo test --release -p foehncraft-cli --test run_contract
python3 -m tools.plot.estimate_chart estimate.svg
```

[↑ back to top](#run-control-health-and-adjoint-readiness)

---

<a id="health"></a>
## 19. Balances and the divergence report

![Numerical dissipation of kinetic energy over the mesh size, lid-driven cavity and Taylor–Green vortex](../figures/gallery/energy_budget.svg)

Every run books mass and kinetic energy: what flows and works through
each boundary, what the body force adds, what viscosity dissipates, and what
is left, the dissipation of the numerical schemes. The second law asks
that this rest is not negative. In the steady lid-driven cavity the lid
works on the fluid and everything ends up as dissipation; the schemes take
about 1.5 % of it, falling slowly with the mesh because the two corners
where the lid meets the walls are singular. The ledger splits the rest in
two: the pressure part, which the momentum interpolation of the face
fluxes dissipates, and the part of the convection and time schemes. In
the cavity the pressure part is nearly all of it and stays of the same
size at the singular corners, while the linear-upwind correction gives
back a little, less with every refinement. On the smooth Taylor–Green
vortex both parts are positive: the numerical dissipation falls with
order 2.5, its pressure part with order 3.

```
$ foehn flow cavity_re100_energy.toml
      balance of the final state
        kinetic energy
          ymax viscous work in                         +2.919438e-3 W
          viscous dissipation                          +2.877098e-3 W
          storage                                       +0.000000e0 W
          numerical dissipation (+1.45e-2)             +4.234029e-5 W
            of it pressure                             +5.115046e-5 W
            of it schemes                              -8.810171e-6 W
```

When a run diverges it keeps its last valid state and says why:

```
$ foehn flow cavity_re1000.toml --set 'solver.algorithm="simple"' --set solver.momentum.relaxation=0.99 --set solver.pressure.relaxation=1.0
      DIVERGED in outer iteration 16: 1 failing cells, velocity first; report in divergence.md
        DIV002 under-relaxation (0.80): relaxation factors velocity 0.99, pressure 1 (SIMPLE); try --set solver.momentum.relaxation=0.69 --set solver.pressure.relaxation=0.31
        DIV006 linear solve (0.60): momentum BiCGStab at its limit of 50 iterations; try --set solver.momentum.max_iter=100
Error RUN014: the solution diverged in outer iteration 16 ...
```

| Check | Result |
|---|---|
| Cavity Re 100, 16² to 128² | numerical dissipation 2.05 %, 1.59 %, 1.45 %, 1.33 % of the lid power, positive; pressure part 7.6e-5 to 5.0e-5 W, scheme part −3.6e-5 to −4.4e-6 W; mass ledger zero |
| Taylor–Green, 16² to 64², fixed Courant number | numerical dissipation 4.1e-4, 7.4e-5, 1.2e-5 J over 0.4 s (pressure part 1.3e-4, 1.6e-5, 1.4e-6 J); mass imbalance at round-off |
| Steady divergence (SIMPLE, no relaxation) | `RUN014`, the first cause is the relaxation; its fix runs |
| Transient divergence (time step 100 times too large) | the first cause is the time step, with a smaller one as the fix |
| The catalogue against part 4a | fields, CSV files and report values bitwise equal; one new warning, `BAL003` on the Courant-10 Taylor–Green sweep |

| | |
|---|---|
| Code | `foehncraft_flow::balance` (ledger, tables, CSV, chart), `divergence` (analysis, snapshot, report), `simplec::Flow::keep_valid`, `rewind` |
| Tests | `crates/foehncraft-flow/tests/health.rs`, unit tests in `balance.rs`, `divergence.rs` |
| Formats | [interfaces.md](../interfaces.md#conservation-ledger-d88-d89-d90-part-4) |
| Case | `examples/cases/flow/cavity_re100_energy.toml` |
| Recipe | [27](../recipes/runs.md#health) |

```
cargo test --release -p foehncraft-flow --test health
python3 -m tools.plot.energy_budget_chart energy_budget.svg
```

[↑ back to top](#run-control-health-and-adjoint-readiness)

---

<a id="log"></a>
## 26. A log that can be circulated

Part of what a run says about itself (D125): the log shows what the case file set, so a log mailed on its own is enough to see the run. Every block is framed by asterisk rules (80 columns, no colour); the console and `flow.log` are identical. Excerpt of the transient cavity (Re 100, 64 x 64, BDF2), blocks between the banner and the table left out:

```
********************************************************************************
  Foehncraft  -  general purpose CFD solver
  Version                       0.5.0  (build a1b2c3d)
  Date and time                 2026-11-16 12:00:00 UTC
  Host name                     vm
  Platform                      linux-x86_64, Intel(R) Xeon(R) Processor @
                                2.10GHz
  Features                      hdf5, metis
  Threads / backend             4 / cpu
  ...
********************************************************************************
  Time steps (PISO, BDF2)

         -- time ------------------- -- scaled residuals ----------------- -- solvers --
    Step    Time [s]   Dt [s] Cour.       ResU     ResV     ResW     ResM   Outer    CG 
       1  1.0000E-01 1.00E-01 4.645   1.00E+00 0.00E+00 0.00E+00 6.33E-02       1    22
      50  5.0000E+00 1.00E-01 6.089   1.03E-04 9.33E-05 0.00E+00 1.20E-05       1    22
  ...
********************************************************************************
  Run finished: every check passed
  Exit code                     0
  Elapsed time                  0:00:25
  Results in                    output/cavity_re100_transient
  Report                        flow.md, result.json
********************************************************************************

```

Before the table, the log lists the input and output files, the run control, the mesh information with its quality, the fluid properties, the solution settings as a matrix per equation and the boundary conditions with their values. After it come the summaries: convergence or time stepping, verification, results per boundary, the ledger, the run time and the warnings. Messages are tagged `[I]`, `[W]` or `[E]`; a failed run ends with `Run FAILED: <code>`, the reason and the exit code. `[log] style = "box"` draws the frames with box characters.

| | |
|---|---|
| Code | `crates/foehncraft-flow/src/console.rs`, `logblocks.rs`, `case.rs`, `unsteady.rs` |
| Tests | `crates/foehncraft-flow/src/console.rs` (frames, tables), `tests/health.rs`, `tests/stage5.rs` |
| Case files | every file of `examples/cases/flow/` |

```
target/release/foehn flow examples/cases/flow/cavity_re100_transient.toml
less output/cavity_re100_transient/flow.log
```

[↑ back to top](#run-control-health-and-adjoint-readiness)

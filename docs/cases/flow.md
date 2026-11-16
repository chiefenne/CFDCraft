# Flow cases (`examples/cases/flow/`)

Part of the [Case catalogue](../cases.md): incompressible laminar flow, steady and transient, grouped by physics; what a flow run writes.

On this page:

- [Steady laminar flow](#steady-laminar-flow)
- [Transient flow](#transient-flow)
- [Runs driven by programs, balances](#runs-driven-by-programs-balances)
- [The flow loop on a GPU](#the-flow-loop-on-a-gpu)

Stage 4: steady incompressible laminar flow, SIMPLEC with the Rhie–Chow
interpolation ([explanation](../verification/stage4.md)). Stage 5: transient
flow with PISO, BDF2 or Crank–Nicolson, when the case has a `[time]` table
([explanation](../verification/stage5.md)). Run with
`foehn flow FILE`; `foehn flow FILE --check` only checks the case file
and the boundary conditions (D30); command-line options override the file
(`foehn flow --help`). Each run writes `flow.log` (the header says what
is written and how large it is), `residuals.csv` (one row per outer
iteration; `residuals_level<N>.csv` per level of a sweep), `patches.csv`
(per patch: area, mass flow, mean pressure, force), `wall_<patch>.csv` (per
wall face: position, pressure, wall shear stress), `line_<name>.csv` for
each sampling line, `fields.vtkhdf` (velocity, speed, pressure and, with an
exact solution, the velocity error), `case.resolved.toml` (the case with
every default filled in), `flow.md` and `flow.json` (the report),
`result.json` (the result summary for programs, without timings) and
`run.jsonl` (the run record, one JSON object per row; both in
[interfaces.md](../interfaces.md)), and the
charts `residuals.svg`, `convergence.svg` (sweeps) and `centreline_u.svg`,
`centreline_v.svg` (cavity against the reference data). The run fails when
a check of the case fails (observed order, error, reference deviation,
convergence, mass balance). CI runs every case on Linux and the quick ones
on macOS and Windows (R5); the sets are in `tools/ci_cases.sh`, and the
long cylinder at Re 100 runs on main, nightly and with the label `full-ci`.

A transient case also writes `history.csv` (one row per time step: time,
Δt, Courant number, residuals, outer and CG iterations; the residuals of a
steady case stay in `residuals.csv`), `forces.csv` (force per wall and the
coefficients, every `output.monitor_every` steps or `monitor_interval`
seconds), `probes.csv`, `fields_<step>.vtkhdf` on the schedule of
`output.fields_every`, `fields_at`, `fields_interval` and `fields_times`
(the last two in simulated time, landed on exactly), `restart.h5` (D42), `coefficients.svg` (force coefficients;
for the periodic cylinder the evaluation window) and, for a time-order
sweep, `time_convergence.svg`. Control files in the output directory steer
a running case (D40): `.stop` ends it with its restart file and results,
`.result_N` writes the fields at step (or iteration) N, `.read_toml_N`
re-reads the case file at N and applies what may change during a run.

## Steady laminar flow

Couette and Poiseuille flow, the lid-driven cavity and the steady cylinder.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/couette_hex.toml` | generated hex 2D, 128 cells | plane Couette flow with convection: error at the solver tolerance (≤ 1e-9 required) | [recipe 18](../recipes/flow.md#flow-verification) |
| `flow/couette_tet.toml`, `couette_poly.toml`, `couette_mixed.toml` | generated 2D, perturbed by 15 %, 128 to 256 cells | the same Stokes flow on perturbed tetrahedra (prisms), polyhedra and mixed cells | [recipe 18](../recipes/flow.md#flow-verification) |
| `flow/poiseuille_hex.toml` | generated hex 2D, perturbed, 512 to 8 192 cells (3 levels) | plane Poiseuille flow, observed order of the velocity error 1.8–2.2 required | [gallery 13](../gallery/flow.md#poiseuille), [recipe 18](../recipes/flow.md#flow-verification) |
| `flow/poiseuille_tet.toml`, `poiseuille_poly.toml`, `poiseuille_mixed.toml` | generated 2D, 3 levels | the same on prisms, perturbed polyhedra and mixed cells | [gallery 13](../gallery/flow.md#poiseuille), [recipe 18](../recipes/flow.md#flow-verification) |
| `flow/cavity_re100.toml` | generated hex 2D, 64 × 64 | lid-driven cavity, Re = 100, against Ghia et al. (1982), deviation ≤ 2 % of the lid speed | [gallery 12](../gallery/flow.md#cavity), [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cavity_re1000.toml` | generated hex 2D, 128 × 128 | the same at Re = 1000 | [gallery 12](../gallery/flow.md#cavity), [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cavity_re100_tet3d.toml` | generated tet 3D, 24 576 tetrahedra in a slab with symmetry planes | the cavity on fully three-dimensional tetrahedra, deviation ≤ 3 % | [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cavity_re100_poly3d.toml` | generated poly 3D, 5 445 polyhedra | the same on polyhedra | [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cavity_re1000_poly3d.toml` | generated poly 3D, 28 227 polyhedra | Re = 1000 on polyhedra, deviation ≤ 3 % | [gallery 12](../gallery/flow.md#cavity), [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cylinder_re20.toml` | built-in cylinder mesh, refine 1: 24 576 hexahedra | Schäfer & Turek 2D-1, steady, Re = 20: C_D, C_L, Δp within 0.5 % of the reference ranges | [gallery 15](../gallery/flow.md#cylinder), [recipe 23](../recipes/flow.md#cylinder) |
| `flow/cylinder_potential.toml` | built-in cylinder mesh, refine 1: 24 576 hexahedra | the potential-flow start (D52) around the cylinder, written without outer iterations | [gallery 20](../gallery/flow.md#time-control), [recipe 28](../recipes/flow.md#time-control) |

## Transient flow

Time-stepped cavity, Taylor–Green vortex, periodic cylinder, error-based time step.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/cavity_re100_transient.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 time-stepped (PISO, BDF2, Courant number 6) to its steady state: the SIMPLEC solution, deviation from Ghia ≤ 2 % | [recipe 21](../recipes/flow.md#transient-flow) |
| `flow/taylor_green_bdf2.toml` | generated hex 2D, 64 × 64 | Taylor–Green vortex, plain PISO with BDF2 over five time steps: time order 1.8–2.3 required | [gallery 14](../gallery/flow.md#taylor-green), [recipe 21](../recipes/flow.md#transient-flow) |
| `flow/taylor_green_cn.toml` | generated hex 2D, 64 × 64 | the same with Crank–Nicolson | [gallery 14](../gallery/flow.md#taylor-green), [recipe 21](../recipes/flow.md#transient-flow) |
| `flow/taylor_green_space.toml` | generated hex 2D, 16² to 128² (4 levels) | Taylor–Green vortex, space order 1.8–2.2 required | [gallery 14](../gallery/flow.md#taylor-green), [recipe 21](../recipes/flow.md#transient-flow) |
| `flow/taylor_green_cfl10.toml` | generated hex 2D, 64 × 64 | Taylor–Green vortex at Courant numbers 1.3 to 9 with outer iterations: time order 2 | [gallery 14](../gallery/flow.md#taylor-green), [recipe 21](../recipes/flow.md#transient-flow) |
| `flow/taylor_green_error.toml` | generated hex 2D, 64 × 64 | the time step from the BDF2 error estimate with rejected steps (D65), probe rows every 0.05 s and fields every 0.1 s landed on exactly (D66) | [gallery 20](../gallery/flow.md#time-control), [recipe 28](../recipes/flow.md#time-control) |
| `flow/cylinder_re100.toml` | built-in cylinder mesh, refine 0: 6 144 hexahedra | Schäfer & Turek 2D-2, periodic, Re = 100: C_D,max, C_L,max, St, Δp within 2.5 % (about 8 minutes; Linux, on main, nightly and with the label `full-ci`) | [gallery 15](../gallery/flow.md#cylinder), [recipe 23](../recipes/flow.md#cylinder) |

## Runs driven by programs, balances

Limits, exit codes, restart and the energy ledger.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/cavity_re100_budget.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 with `[run.limits]` (wall time, iterations, output volume), run the way a script or an agent runs it; a tighter limit stops it with exit code 3 and `--restart latest` continues bitwise | [gallery 18](../gallery/runs.md#run-contract), [recipe 26](../recipes/runs.md#run-contract) |
| `flow/cavity_re100_energy.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 with the balance table every 100 iterations: the lid power goes into viscous and numerical dissipation, the numerical part positive (D90 part 4) | [gallery 19](../gallery/runs.md#health), [recipe 27](../recipes/runs.md#health) |

## The flow loop on a GPU

CPU = GPU bitwise, not run by CI.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/cavity_poly3d_gpu.toml` | generated poly 3D, 970 225 polyhedra | the steady flow loop on the GPU: 30 SIMPLEC iterations of the cavity at Re = 100 on the CPU and the GPU, bitwise equal, time per cell and iteration (GPU; not run by CI) | [gallery 16](../gallery/benchmarks.md#gpu-flow), [recipe 24](../recipes/gpu.md#flow-gpu) |
| `flow/cylinder_re100_gpu.toml` | built-in cylinder mesh, refine 2: 98 304 hexahedra | the transient flow loop on the GPU: 200 adaptive BDF2 steps of the Re = 100 cylinder on the CPU and the GPU, bitwise equal (GPU; not run by CI) | [gallery 16](../gallery/benchmarks.md#gpu-flow), [recipe 24](../recipes/gpu.md#flow-gpu) |

The adjoint-ready kernels of stage 5 (D59) have no case file: their
checks are tests on generated meshes (gallery 17, recipe 25).

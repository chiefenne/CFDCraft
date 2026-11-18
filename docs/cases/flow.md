# Flow cases (`examples/cases/flow/`)

Part of the [Case catalogue](../cases.md): incompressible laminar flow, steady and transient, grouped by physics; what a flow run writes.

On this page:

- [Steady laminar flow](#steady-laminar-flow)
- [Transient flow](#transient-flow)
- [Runs driven by programs, balances](#runs-driven-by-programs-balances)
- [The flow loop on a GPU](#the-flow-loop-on-a-gpu)
- [Boundary conditions](#boundary-conditions)
- [Derived fields](#derived-fields)
- [Expressions and materials](#expressions-and-materials)
- [Passive scalars](#passive-scalars)

Stage 4: steady incompressible laminar flow, SIMPLEC with the Rhie–Chow
interpolation ([explanation](../verification/stage4.md)). Stage 5: transient
flow with PISO, BDF2 or Crank–Nicolson, when the case has a `[time]` table
([explanation](../verification/stage5.md)). Run with
`foehn flow FILE`; `foehn flow FILE --check` only checks the case file
and the boundary conditions (D30); command-line options override the file
(`foehn flow --help`). Each run writes `flow.log` (the header says what
is written and how large it is), `residuals.csv` (one row per outer
iteration; `residuals_level<N>.csv` per level of a sweep), `boundaries.csv`
(per boundary: area, mass flow, mean pressure, force), `wall_<boundary>.csv` (per
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
on macOS and Windows (R5); the sets are in `tools/check/ci_cases.py`, and the
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
| `flow/cavity_re1000_coupled.toml` | generated hex 2D, 128 × 128 | the cavity at Re = 1000 with the coupled algorithm and the automatic pseudo-time step: about 450 iterations instead of 1 397, same result | [gallery 78](../gallery/flow.md#coupled), [gallery 80](../gallery/flow.md#momentum-multigrid), [recipe 78](../recipes/flow.md#coupled) |
| `flow/cavity_re100_tet3d.toml` | generated tet 3D, 24 576 tetrahedra in a slab with symmetry planes | the cavity on fully three-dimensional tetrahedra, deviation ≤ 3 % | [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cavity_re100_poly3d.toml` | generated poly 3D, 5 445 polyhedra | the same on polyhedra | [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cavity_re1000_poly3d.toml` | generated poly 3D, 28 227 polyhedra | Re = 1000 on polyhedra, deviation ≤ 3 % | [gallery 12](../gallery/flow.md#cavity), [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cylinder_re20.toml` | built-in cylinder mesh, refine 1: 24 576 hexahedra | Schäfer & Turek 2D-1, steady, Re = 20: C_D, C_L, Δp within 0.5 % of the reference ranges | [gallery 15](../gallery/flow.md#cylinder), [recipe 23](../recipes/flow.md#cylinder) |
| `flow/cylinder_potential.toml` | built-in cylinder mesh, refine 1: 24 576 hexahedra | the potential-flow start (D52) around the cylinder, written without outer iterations | [gallery 20](../gallery/flow.md#time-control), [recipe 28](../recipes/flow.md#time-control) |
| `flow/cylinder_re20_gci.toml` | built-in cylinder mesh, refine 0 to 2: 6 144 to 98 304 hexahedra (3 levels) | the benchmark 2D-1 on three meshes: observed order, extrapolated value and grid convergence index of C_D, C_L and Δp (D67); every mesh within 8 % of the reference ranges (about 16 minutes, in the `long` set) | [gallery 25](../gallery/verification.md#gci), [recipe 34](../recipes/verification.md#gci) |
| `flow/pipe_gmsh_raw.toml` | Gmsh tetrahedra without optimisation (`examples/meshes/pipe_gmsh_raw.msh`), 1 572 cells, 28 faces with a skewness above 0.9 | developing pipe flow at Re = 100 on a poor mesh: the cell-quality remediation (D50) switches to the Gauss pressure gradient and the limiter, reports it and converges; with `--set solver.remediation=off` the run diverges with a report that names the skewed faces | [gallery 26](../gallery/meshes.md#skewed-tets), [recipe 36](../recipes/meshes.md#skewed-tets) |

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
| `flow/cavity_re100_columns.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 with `[log] columns`: residuals, a probe in the lower-left vortex, one cell, the largest velocity and lowest pressure with their cells, the lid mass flow; bitwise the results of `cavity_re100.toml` | [gallery 28](../gallery/runs.md#log-columns), [recipe 38](../recipes/runs.md#log-columns) |
| `flow/cavity_re100_schedule.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 started with momentum under-relaxation 0.5; `[[schedule]]` raises it to 0.9 from iteration 100 and thins the log from iteration 200; every change is in the log and in `run.jsonl` | [gallery 29](../gallery/runs.md#schedule), [recipe 39](../recipes/runs.md#schedule) |
| `flow/cavity_re100_energy.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 with the balance table every 100 iterations: the lid power goes into viscous and numerical dissipation, the numerical part positive (D90 part 4) | [gallery 19](../gallery/runs.md#health), [recipe 27](../recipes/runs.md#health) |

## The flow loop on a GPU

CPU = GPU bitwise, not run by CI.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/cavity_poly3d_gpu.toml` | generated poly 3D, 970 225 polyhedra | the steady flow loop on the GPU: 30 SIMPLEC iterations of the cavity at Re = 100 on the CPU and the GPU, bitwise equal, time per cell and iteration (GPU; not run by CI) | [gallery 16](../gallery/benchmarks.md#gpu-flow), [recipe 24](../recipes/gpu.md#flow-gpu) |
| `flow/cylinder_re100_gpu.toml` | built-in cylinder mesh, refine 2: 98 304 hexahedra | the transient flow loop on the GPU: 200 adaptive BDF2 steps of the Re = 100 cylinder on the CPU and the GPU, bitwise equal (GPU; not run by CI) | [gallery 16](../gallery/benchmarks.md#gpu-flow), [recipe 24](../recipes/gpu.md#flow-gpu) |

The adjoint-ready kernels of stage 5 (D59) have no case file: their
checks are tests on generated meshes (gallery 17, recipe 25).

## Boundary conditions

Slip walls, axisymmetric wedges, inlet profiles from tables, backflow
treatment at pressure outlets and the convective outlet (D29, stage 5
part 5b), periodic pairs with the heat-flux jump of the scalars (D29,
D109 (1), stage 6 part 6p); the keys are on the [key page](flow-keys.md#boundary-conditions-d29).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/cylinder_re100.toml` | built-in cylinder mesh, refine 0: 6 144 hexahedra | Schäfer & Turek 2D-2, periodic, Re = 100: C_D,max, C_L,max, St, Δp within 2.5 % (about 8 minutes; Linux, on main, nightly and with the label `full-ci`) | gallery 15, recipe 23 |
| `flow/channel_slip.toml` | generated mixed 2D, perturbed by 20 %, 432 cells | uniform flow between slip walls, exact (error ≤ 1e-9 required) | [gallery 21](../gallery/boundary-conditions.md#boundary-conditions), [recipe 29](../recipes/boundary-conditions.md#boundary-conditions) |
| `flow/poiseuille_table.toml` | generated hex 2D, 512 cells | the inlet profile read from `poiseuille_profile.csv`: the run with the exact parabola | [gallery 21](../gallery/boundary-conditions.md#boundary-conditions), [recipe 29](../recipes/boundary-conditions.md#boundary-conditions) |
| `flow/pipe_wedge.toml` | generated hex 2D revolved into a 0.5° wedge, perturbed, 256 to 4 096 cells (3 levels) | Hagen–Poiseuille flow on the axisymmetric wedge, observed order 1.8–2.2 required | [gallery 21](../gallery/boundary-conditions.md#boundary-conditions), [recipe 29](../recipes/boundary-conditions.md#boundary-conditions) |
| `flow/cylinder_re100_short.toml` | built-in cylinder mesh with the outlet 5.5 D behind the cylinder, 4 928 hexahedra | the Re = 100 vortex street through a static pressure outlet with backflow treatment (Linux, on main, nightly and with the label `full-ci`) | [gallery 21](../gallery/boundary-conditions.md#boundary-conditions), [recipe 29](../recipes/boundary-conditions.md#boundary-conditions) |
| `flow/cylinder_re100_convective.toml` | as `cylinder_re100_short.toml` | the same with a convective outlet (Linux, on main, nightly and with the label `full-ci`) | [gallery 21](../gallery/boundary-conditions.md#boundary-conditions), [recipe 29](../recipes/boundary-conditions.md#boundary-conditions) |
| `flow/poiseuille_periodic.toml` | generated hex 2D, perturbed by 15 %, 512 to 8 192 cells (3 levels) | plane Poiseuille flow on one periodic pair driven by its mass flow 0.1 kg/s: observed order 1.8–2.2 required (1.98, 1.94), the driving force converges to 12 μU/h² = 1.2 Pa/m (+0.09 % on the finest mesh) | [gallery 67](../gallery/boundary-conditions.md#periodic), [recipe 67](../recipes/boundary-conditions.md#periodic) |
| `flow/channel_heated_periodic.toml` | as `poiseuille_periodic.toml` | the laminar heated channel on one periodic module with the heat-flux jump: observed order of the temperature error 1.8–2.2 required, Nu = 8.2460 against 140/17 = 8.2353 on the finest mesh | [gallery 67](../gallery/boundary-conditions.md#periodic), [recipe 67](../recipes/boundary-conditions.md#periodic) |
| `flow/pipe_heated_periodic.toml` | generated hex 2D revolved into a 0.5° wedge, perturbed, 64 to 1 024 cells (3 levels) | the short periodic heated pipe (one diameter, D111): Hagen–Poiseuille flow at a fixed mass flow and the heat-flux jump; orders 1.8–2.2 (velocity) and 1.8–2.6 (temperature) required, Nu = 4.3691 against 48/11 = 4.3636 on the finest mesh | [gallery 67](../gallery/boundary-conditions.md#periodic), [recipe 67](../recipes/boundary-conditions.md#periodic) |

## Turbulence models

RANS turbulence models with their wall treatments (stage 6, D129–D139);
the keys are on the [key page](flow-keys.md#turbulence-d129d139).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/tmr_plate_sst.toml` | generated hex 2D on the TMR grid family (8 + 24 cells along x, plate split from the symmetry plane at x = 0, wall grading 38 000), 3 072 cells (`verification.levels = [0, 1, 2, 3]`: 768 to 49 152) | the TMR zero-pressure-gradient flat plate at Re = 5·10⁶/m with k-ω SST and `low-re`: C_f(0.97) against the TMR value 2.6907e-3 (+0.37 % on the finest mesh), the wall shear along the plate, the profile at x = 0.97, the D166 features | [gallery 68](../gallery/turbulence.md#sst), [recipe 68](../recipes/turbulence.md#sst) |
| `flow/channel_sst.toml` | generated hex 2D, 2 × 96 cells, graded towards both walls (ratio 60), periodic in x | fully developed channel flow at Re_τ = 395 driven by the DNS bulk velocity: U_b⁺ = 17.43 against 17.55 (Moser, Kim & Mansour), −0.7 % | [gallery 68](../gallery/turbulence.md#sst), [recipe 68](../recipes/turbulence.md#sst) |
| `flow/tmr_plate_keps.toml` | generated hex 2D, 16 + 48 cells along x, 48 across with wall grading 361 (first-cell y⁺ 30 at x = 0.97), 3 072 cells | the TMR flat plate with the standard k-ε model, `scalable` wall functions and the sustaining terms: C_f(0.97) = 2.829e-3, −1.6 % from White's correlation | [gallery 71](../gallery/turbulence.md#k-epsilon), [recipe 71](../recipes/turbulence.md#k-epsilon) |
| `flow/channel_keps.toml` | generated hex 2D, 2 × 12 uniform cells (first-cell y⁺ 32), periodic in x | channel flow at Re_τ = 395 with the standard k-ε and `scalable`: U_b⁺ = 18.31, +4.3 % from the DNS | [gallery 71](../gallery/turbulence.md#k-epsilon), [recipe 71](../recipes/turbulence.md#k-epsilon) |
| `flow/pipe_rough_keps.toml` | generated hex 2D revolved into a 0.5° wedge, 2 × 40 uniform cells, periodic in x | fully developed rough-pipe flow at Re = 10⁶, k_s/D = 10⁻³, Cebeci–Bradshaw shift: f = 0.02081, +4.4 % from Colebrook; base case of the Moody study | [gallery 71](../gallery/turbulence.md#k-epsilon), [recipe 71](../recipes/turbulence.md#k-epsilon) |
| `flow/channel_zeta_f.toml` | generated hex 2D, 2 × 96 cells graded 60 towards both walls (first-cell y⁺ 0.27), periodic in x | channel flow at Re_τ = 395 with k-ε-ζ-f and the `compound` treatment: U_b⁺ = 17.89, +1.9 % from the DNS | [gallery 73](../gallery/turbulence.md#zeta-f), [recipe 73](../recipes/turbulence.md#zeta-f) |
| `flow/tmr_plate_zeta_f.toml` | generated hex 2D, 16 + 48 cells along x, 48 across with wall grading 17 555 (first-cell y⁺ 1 at x = 0.97), 3 072 cells | the TMR flat plate with k-ε-ζ-f, `compound`, a turbulent start and the sustaining terms: C_f(0.97) = 3.023e-3, +5.1 % from White's correlation (finding) | [gallery 73](../gallery/turbulence.md#zeta-f), [recipe 73](../recipes/turbulence.md#zeta-f) |
| `flow/channel_zeta_f0.toml` | generated hex 2D, 2 × 48 cells graded 27.3 towards both walls (first-cell y⁺ 1), periodic in x | channel flow at Re_τ = 395 with k-ε-ζ-f₀ (Popovac 2023) and `compound-2023`: U_b⁺ = 17.79, +1.4 % from the DNS | [gallery 74](../gallery/turbulence.md#zeta-f0), [recipe 74](../recipes/turbulence.md#zeta-f0) |
| `flow/tmr_plate_zeta_f0.toml` | the mesh of `tmr_plate_zeta_f.toml` (first-cell y⁺ 1 at x = 0.97), 3 072 cells | the TMR flat plate with k-ε-ζ-f₀ from the TMR free stream: turns turbulent on its own, every field converged, C_f(0.97) = 2.802e-3, −2.5 % from White's correlation | [gallery 74](../gallery/turbulence.md#zeta-f0), [recipe 74](../recipes/turbulence.md#zeta-f0) |
| `flow/backstep_zeta_f0.toml` | generated hex 2D, 125 × 80 cells with the step cut out (`remove`), 8 200 cells, first-cell y⁺ about 1.5 | the backward-facing step of Driver & Seegmiller (Re_H = 36 000) with k-ε-ζ-f₀: reattachment x_r/H = 6.38, +2 % from the experiment | [gallery 74](../gallery/turbulence.md#zeta-f0), [recipe 74](../recipes/turbulence.md#zeta-f0) |
| `flow/backstep_zeta_f0_coupled.toml` | the mesh of `backstep_zeta_f0.toml`, 8 200 cells | the same step with the coupled algorithm and pseudo-time steps, at most 150 iterations: x_r within 1 % after about 40 iterations instead of 550 to 600, the same x_r/H = 6.385 | [gallery 78](../gallery/flow.md#coupled), [gallery 80](../gallery/flow.md#momentum-multigrid), [gallery 81](../gallery/flow.md#turbulence-multigrid), [recipe 78](../recipes/flow.md#coupled), [recipe 80](../recipes/flow.md#momentum-multigrid), [recipe 81](../recipes/flow.md#turbulence-multigrid) |
| `studies/backstep/keps.toml`, `zf0.toml`, `sst.toml` (in `examples/studies/`) | the mesh of `backstep_zeta_f0.toml` with each interval halved, 250 × 160 cells with the step cut out, 32 800 cells | the same step with k-ε and scalable wall functions, k-ε-ζ-f₀ and k-ω SST, with field output: x_r/H = 5.36, 6.35 and 6.80 against the experiment's 6.26 after 3 000 iterations; the velocity and k figures of the gallery | [gallery 74](../gallery/turbulence.md#zeta-f0), [recipe 74](../recipes/turbulence.md#zeta-f0) |
| `flow/channel_eb.toml` | the mesh of `channel_zeta_f.toml`, 2 × 96 cells graded 60 towards both walls (first-cell y⁺ 0.27), periodic in x | channel flow at Re_τ = 395 with the elliptic-blending k-ε (Billard & Laurence 2012): C_f at 99.0 % of the DNS (the paper: 98 %), U_b⁺ = 17.64, +0.5 % | [gallery 88](../gallery/turbulence.md#elliptic-blending), [recipe 88](../recipes/turbulence.md#elliptic-blending) |
| `flow/tmr_plate_eb.toml` | the mesh of `tmr_plate_zeta_f.toml` (first-cell y⁺ 1 at x = 0.97), 3 072 cells | the TMR flat plate with the elliptic-blending k-ε from a turbulent start: C_f(0.97) = 2.889e-3, +0.5 % from White's correlation | [gallery 88](../gallery/turbulence.md#elliptic-blending), [recipe 88](../recipes/turbulence.md#elliptic-blending) |
| `flow/backstep_eb.toml` | the mesh of `backstep_zeta_f0.toml`, 8 200 cells | the backward-facing step with the elliptic-blending k-ε (SIMPLEC): x_r/H = 7.59, +21 % from the experiment (finding: the variable C*_ε2; 6.24 with it constant) | [gallery 88](../gallery/turbulence.md#elliptic-blending), [recipe 88](../recipes/turbulence.md#elliptic-blending) |

## Derived fields

Fields derived from the solution with the output or afterwards from a
restart file (D72, stage 5 part 5c); the keys are on the
[key page](flow-keys.md#derived-fields-d72).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/couette_derived.toml` | generated poly 2D, perturbed by 15 %, 153 cells | every derived field of plane Couette flow, exact (vorticity −1/s, Q = λ₂ = 0, wall shear stress 0.1 Pa); `foehn post derive` writes the same file from `restart.h5` | [gallery 22](../gallery/post-processing.md#derived-fields), [recipe 30](../recipes/post-processing.md#derived-fields) |
| `flow/cylinder_re100_short.toml` | built-in cylinder mesh, 4 928 hexahedra | vorticity and Q criterion of the vortex street at t = 8 s, from the restart file with `foehn post derive` | [gallery 22](../gallery/post-processing.md#derived-fields), [recipe 30](../recipes/post-processing.md#derived-fields) |

## Expressions and materials

Values as expressions of x, y, z and t, constants and tables (D53), the
material database and regions (D102, D112), walls with a mass flux
(D124 W3), stage 5 part 6a; the keys are on the
[key page](flow-keys.md#expressions-materials-and-regions-d53-d102-d112).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/couette_viscosity.toml` | generated hex 2D, perturbed by 10 %, 128 to 2 048 cells (3 levels) | Couette flow with the viscosity μ0 (1 + 4 y) from an expression: u = U ln(1 + a y)/ln(1 + a), observed order 1.8–2.2 required | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 31](../recipes/properties.md#expressions) |
| `flow/couette_transpiration.toml` | generated hex 2D, perturbed by 10 %, 128 to 2 048 cells (3 levels) | Couette flow with blowing and suction through the walls (`mass_flux`), Re = 4, observed order 1.8–2.2 required | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 31](../recipes/properties.md#expressions) |
| `flow/poiseuille_water.toml` | generated hex 2D, perturbed by 10 %, 2 048 cells | water at 313.15 K from the material database, in a region of the zone `fluid`; Poiseuille flow in a 1 mm channel, error ≤ 1e-3 required | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 32](../recipes/properties.md#materials) |
| `flow/poiseuille_expression.toml` | as `poiseuille_hex.toml`, 3 levels | the inlet parabola and half of the driving pressure gradient as expressions: the errors of `poiseuille_hex.toml` | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 31](../recipes/properties.md#expressions) |
| `flow/taylor_green_expression.toml` | generated hex 2D, 4 096 cells | the decaying vortex prescribed on the sides as expressions of x, y and t, error ≤ 2e-3 at t = 0.4 s required | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 31](../recipes/properties.md#expressions) |
| `flow/channel_ramp.toml` | generated hex 2D, 512 cells | start-up of a channel with the inflow ramped by a table in t, BDF2 with an adaptive step; steady Poiseuille flow at t = 4 s (error ≤ 2e-3) | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 31](../recipes/properties.md#expressions) |

## Passive scalars

Scalars carried by the flow (D51) with the generic scalar transport
equation of stage 6 (D116), the bounded BDF2 (D28 (1)) and the exact wall
distance (D48), stage 6 part 6a; the keys are on the
[key page](flow-keys.md#passive-scalars-d51-d116).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/channel_heated.toml` | generated hex 2D, perturbed by 15 %, 512 to 8 192 cells (3 levels) | laminar heated channel: a passive temperature with uniform wall flux in Poiseuille flow, Pe = 1000; observed order of the temperature error 1.8–2.2 required, Nu = 8.2398 against 140/17 = 8.2353 on the finest mesh | [gallery 63](../gallery/scalars.md#passive-scalars), [recipe 63](../recipes/scalars.md#passive-scalars) |
| `flow/scalar_front_bdf2.toml` | generated hex 2D, 40 cells | a front of a scalar bounded by [0, 1] in a uniform flow, upwind, BDF2 at Courant number 0.8: the bounded BDF2 keeps it inside the bounds (plain BDF2 overshoots to 1.00033) | [gallery 63](../gallery/scalars.md#passive-scalars), [recipe 63](../recipes/scalars.md#passive-scalars) |

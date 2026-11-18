# Incompressible flow

Part of the [Recipes](../recipes.md): Couette, Poiseuille, cavity and cylinder runs, boundary conditions, transient runs, time-step control.

| # | Recipe |
|---|---|
| 18 | [Verify the flow solver on Couette and Poiseuille flow](#flow-verification) |
| 19 | [Run the lid-driven cavity and compare with Ghia et al.](#flow-cavity) |
| 20 | [Set up the boundary conditions of a flow case](#flow-boundaries) |
| 21 | [Run a transient flow case and check its time order](#transient-flow) |
| 23 | [Run the cylinder benchmark of Schäfer & Turek](#cylinder) |
| 28 | [Let the error choose the time step, write output at exact times, start from potential flow](#time-control) |
| 78 | [Converge a steady case in fewer iterations with the coupled algorithm](#coupled) |
| 80 | [Choose and tune the momentum sweep of the coupled algorithm](#momentum-multigrid) |
| 81 | [Use the multigrid of the turbulence equations and its larger steps](#turbulence-multigrid) |

---

<a id="flow-verification"></a>
## 18. Verify the flow solver on Couette and Poiseuille flow

The exact solutions of stage 4 ([explanation](../verification/stage4.md)):
Couette flow must come out at the iteration error on every mesh type,
Poiseuille flow with an observed order of 1.8–2.2. Each case checks its own
result and fails otherwise.

```
cargo build --release -p foehncraft-cli
for m in hex tet poly mixed; do
  target/release/foehn flow examples/cases/flow/couette_${m}.toml
done
for m in hex tet poly mixed; do
  target/release/foehn flow examples/cases/flow/poiseuille_${m}.toml
done
cat output/poiseuille_poly/flow.md          # errors per level, observed order, boundaries
# plots: output/poiseuille_poly/convergence.svg, residuals.svg
```

Keys used: `[analytic]` (`type = "couette"` or `"poiseuille"`),
`profile = "analytic"` at the velocity inlet, `[verification]` `levels`,
`expected_order`, `max_error`. More levels or other settings from the
command line:

```
target/release/foehn flow examples/cases/flow/poiseuille_hex.toml --levels 0,1,2,3
target/release/foehn flow examples/cases/flow/poiseuille_hex.toml --threads 1 --max-iterations 8000
```

The checks as tests (the solving ones in release builds):

```
cargo test --release -p foehncraft-flow --test stage4
```

[↑ back to top](#incompressible-flow)

---

<a id="flow-cavity"></a>
## 19. Run the lid-driven cavity and compare with Ghia et al.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cavity_re100.toml       # about 5 s
target/release/foehn flow examples/cases/flow/cavity_re1000.toml      # about 1 min
target/release/foehn flow examples/cases/flow/cavity_re100_poly3d.toml
cat output/cavity_re1000/flow.md
# centre lines against Ghia et al.: output/cavity_re1000/centreline_u.svg, centreline_v.svg
# velocity and pressure along the lines: line_vertical.csv, line_horizontal.csv
# wall pressure and shear stress: wall_ymax.csv (lid) and the other walls
pip install vtk
python3 -m tools.mesh.render_mesh output/cavity_re1000/fields.vtkhdf cavity.png --field speed --top --no-edges
```

Keys used: `[verification]` `reference = "ghia-re100"` or `"ghia-re1000"`,
`reference_tolerance` (lid speeds), `lid_velocity`; `[[output.lines]]`
for the sampling lines. The Reynolds number is ρ U L / μ with the unit
cavity; the case is refused if it does not match the reference.

Convergence: the outer iterations grow with the mesh and the Reynolds
number. A larger `relaxation_velocity` (0.9 to 0.95) needs fewer; lower it
(0.7, or `scheme = "upwind"` for the first iterations) if a run diverges:

```toml
[solver]
relaxation_velocity = 0.9
tolerance = 1e-8
max_iterations = 10000
```

[↑ back to top](#incompressible-flow)

---

<a id="flow-boundaries"></a>
## 20. Set up the boundary conditions of a flow case

Every boundary of the mesh needs a `[boundary.<boundary>]` table. Check a case
before running it: the check (D30) names each boundary that is missing,
unknown (with "did you mean"), over-specified or unsuitable, and says what
to change; `--check` stops after it.

```
target/release/foehn mesh examples/cases/mesh/hex_2d.toml     # boundary names: output/hex_2d/mesh_report.md
target/release/foehn flow my_case.toml --check
```

A channel with a known flow rate:

```toml
[boundary.inlet]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]        # or normal_velocity = 1.0, or mass-flow-inlet
[boundary.outlet]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.walls]
type = "wall"
```

A pressure-driven flow (the flow rate is a result):

```toml
[boundary.inlet]
type = "total-pressure-inlet"
total_pressure = 10.0
[boundary.outlet]
type = "static-pressure-outlet"
pressure = 0.0
```

An enclosed flow (no pressure boundary: the reference cell holds the level):

```toml
[boundary.lid]
type = "wall"
velocity = [1.0, 0.0, 0.0]
[pressure]
reference_point = [0.5, 0.5, 0.0]
reference_value = 0.0
```

The table of sensible combinations and the cases each one suits is in the
[case catalogue](../cases/flow-keys.md#sensible-combinations-d30); every condition and
key is listed in [boundary conditions](../cases/flow-keys.md#boundary-conditions-d29).

[↑ back to top](#incompressible-flow)

---

<a id="transient-flow"></a>
## 21. Run a transient flow case and check its time order

A case is transient when it has a `[time]` table; it then runs PISO with
BDF2 (or Crank–Nicolson) instead of the steady SIMPLEC iterations.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/taylor_green_bdf2.toml    # about 40 s
target/release/foehn flow examples/cases/flow/taylor_green_cn.toml
target/release/foehn flow examples/cases/flow/taylor_green_space.toml   # space order, about 30 s
target/release/foehn flow examples/cases/flow/taylor_green_cfl10.toml   # Courant 1 to 9, outer iterations
cat output/taylor_green_bdf2/flow.md
# plot: output/taylor_green_bdf2/time_convergence.svg, per step: history.csv
```

The smallest transient case file adds this to a steady one:

```toml
[time]
scheme = "bdf2"          # or "crank-nicolson", "bdf1"
time_step = 0.01         # s
end_time = 2.0           # or steps = 200
adaptive = true          # optional: Δt from the largest Courant number
max_cfl = 0.8
```

Plain PISO (one outer iteration, two correctors) is second order in time
up to Courant numbers of about 1. For larger steps add outer iterations,
which end when the residuals of a step are below `outer_tolerance`:

```toml
[solver]
outer_iterations = 50
outer_tolerance = 1e-8
relaxation_velocity = 0.9   # below 1 only with outer_iterations > 1
```

A time-order sweep runs the case once per step of
`verification.time_steps` (largest first) on the same mesh and reports
the order from the differences of successive solutions; with
`expected_time_order` it fails outside that range:

```toml
[verification]
time_steps = [0.01, 0.005, 0.0025, 0.00125]
expected_time_order = [1.8, 2.3]
```

Keys: [transient cases](../cases/flow-keys.md#flow-case-file-keys); the method and the
results: [stage 5](../verification/stage5.md).

[↑ back to top](#incompressible-flow)

---

<a id="cylinder"></a>
## 23. Run the cylinder benchmark of Schäfer & Turek

```
target/release/foehn flow examples/cases/flow/cylinder_re20.toml     # steady, Re 20, about 1 min
target/release/foehn flow examples/cases/flow/cylinder_re100.toml    # periodic, Re 100, about 8 min
cat output/cylinder_re100/flow.md
# C_D, C_L per step: output/cylinder_re100/forces.csv, chart: coefficients.svg
# probes in front of and behind the cylinder: probes.csv
pip install vtk
python3 -m tools.mesh.render_mesh output/cylinder_re20/fields.vtkhdf cylinder.png --field speed --top --no-edges --size 1800x420 --zoom 2.9
```

The mesh is built in (`[mesh.cylinder]`); `refine = 1` and `2` give 24 576
and 98 304 cells. Re 100 with `refine = 1` takes about 40 minutes on two
cores and brings C_L,max to 0.5 % of its range; with a finer mesh lower the
first `time_step` (the adaptive step starts from it). On 98 304 cells the steady case needs
`relaxation_velocity = 0.5` and about 12 000 outer iterations:

```toml
[mesh.cylinder]
refine = 2
[solver]
relaxation_velocity = 0.5
max_iterations = 12000
```

Keys: `[verification]` `benchmark = "schafer-turek-2d1"` or
`"schafer-turek-2d2"`, `benchmark_tolerance` (relative to the midpoints of
the reference ranges), `evaluation_start` (2D-2). Force coefficients of
any case: `[output.coefficients]` with `boundaries`, `velocity`, `area`;
point histories: `[[output.probes]]` with `name` and `point`.

[↑ back to top](#incompressible-flow)

---

<a id="time-control"></a>
## 28. Let the error choose the time step, write output at exact times, start from potential flow

**The time step from the error** (D65). With `adaptive = true` and
`error_tolerance`, BDF2 and Crank–Nicolson estimate the local error of
each step from the difference between the step and the quadratic
extrapolation of the last three levels, and take the next step that meets
the tolerance; `max_cfl` and `max_growth` stay upper bounds. A step whose
error is more than the tolerance is taken back exactly and repeated
shorter; the log names each one, `history.csv` holds the steps taken, and
`result.json` counts the repeats under `solves[].transient.rejected`.

```
target/release/foehn flow examples/cases/flow/taylor_green_error.toml --out output/tg_error
grep -E "rejected|error" output/tg_error/flow.log | head
target/release/foehn flow examples/cases/flow/taylor_green_error.toml --set time.error_tolerance=1e-6 --out output/tg_error6
python3 -m tools.plot.time_step_chart output/time_step_error.svg
```

```toml
[time]
scheme = "bdf2"            # or "crank-nicolson"
time_step = 0.001          # the first step; the estimate starts with the fourth
adaptive = true
max_cfl = 5.0              # upper bound
max_growth = 1.2           # upper bound of the step ratio
error_tolerance = 1e-4     # relative to the RMS velocity, per step
error_absolute = 1e-6      # m/s
```

The tolerance is per step: the error at the end of a run is up to the
number of steps times larger. Keep it at least ten times above the
tolerances of the linear solvers (`solver.momentum.tolerance`,
`pressure_tolerance`), whose noise the estimate would otherwise chase. With
plain PISO at Courant numbers above about 2 the error grows faster than
the step, and the estimate cuts the step there (the steep drops in the
gallery chart); outer iterations (`solver.outer_iterations`) let it take
longer steps. A restart continues such a run bitwise.

**Output at exact times** (D66). In simulated time instead of steps:

```toml
[output]
fields_interval = 0.1      # fields_<step>.vtkhdf at t = 0.1, 0.2, ... s
fields_times = [0.25]      # and at these times
monitor_interval = 0.05    # forces.csv and probes.csv rows every 0.05 s
[run]
restart_interval = 1.0     # restart.h5 at t = 1, 2, ... s
restart_times = [2.5]      # restart_<step>.h5 at these times
```

The step before a scheduled time is shortened to land on it, and when the
time lies within two steps, the two steps share the distance, so no step
falls below half the planned one. The files keep the step number in their
names; the time is inside (the root attribute `time` of the VTKHDF file,
the `time` column of the CSV rows).

**Potential-flow start** (D52). `initial.potential_flow = true` replaces
the uniform start velocity by the irrotational flow that carries the
inflow to the pressure boundaries: one Laplace solve with the pressure
multigrid and four passes of the non-orthogonal correction; walls and
symmetry planes let it slip. To see the start field itself, run no outer
iteration:

```
target/release/foehn flow examples/cases/flow/cylinder_potential.toml --out output/potential
python3 -m tools.mesh.render_mesh output/potential/fields.vtkhdf potential.png --field speed --top --no-edges --size 1800x420 --zoom 2.9
target/release/foehn flow examples/cases/flow/cylinder_re20.toml --set initial.potential_flow=true --out output/cyl_pf
```

The last run converges in 1527 outer iterations instead of 1570 from rest.
A restart file, when given, wins over the start.

[↑ back to top](#incompressible-flow)

---

<a id="coupled"></a>
## 78. Converge a steady case in fewer iterations with the coupled algorithm

One key switches a steady case from SIMPLEC to the coupled
pressure–velocity algorithm with pseudo-time steps
([gallery 78](../gallery/flow.md#coupled)):

```toml
[solver]
algorithm = "coupled"
```

or on the command line, for any steady case file:

```
target/release/foehn flow examples/cases/flow/channel_zeta_f0.toml --set solver.algorithm=coupled
```

Remove a `[solver.momentum] relaxation` and a relaxation `[[schedule]]`
that the case set for SIMPLEC: with the coupled algorithm the relaxation
is only a cap on the local Courant number α/(1 − α), and its default is
0.98. The result is the one of SIMPLEC; only the number of iterations
changes.

The step case of gallery 78 and the cavity:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml
target/release/foehn flow examples/cases/flow/cavity_re1000_coupled.toml
```

**Reading the log.** The table has one more column of linear
iterations: `K`, the Krylov iterations of the coupled solve of that outer
iteration (`linear.coupled` in `[log] columns`; `krylov` in
`residuals.csv` and the run record). `U V W` are then the multigrid
cycles of all momentum sweeps of the preconditioner (the BiCGStab
iterations with `momentum_sweep = "bicgstab"`, recipe 80) and `P`
the CG iterations of its pressure solves and of the final pressure
correction. A `K` at its limit (`solver.coupled.max_iter`, 30) in every
iteration means that the step or the cap is too large for the
preconditioner: lower `solver.pseudo_time.factor` or the relaxation.

**The pseudo-time step.** Default: automatic, Δτ = L/U with L the volume
of the domain over the area of its boundaries (without symmetry planes)
and U the larger of the largest prescribed boundary speed (inlets, moving
walls) and the root mean square of the cell speeds. The summary of the
run shows the step of the last iteration (`Pseudo-time step`), the run
record has it per iteration (`pseudo_step`).

```toml
[solver.pseudo_time]
factor = 4.0        # a larger step: fewer iterations, dearer linear solves
# step = 5.0        # or the step itself, in s
# momentum = false  # relaxation only (the cap alone sets the step)
# turbulence = 0.5  # c_T of the turbulence step min(Δτ, c_T T); 0: off
```

A larger `factor` pays on diffusion-dominated laminar cases (cavity at
Re = 400: 292, 130 and 84 iterations with 1, 4 and 20). On the turbulent
step the count does not change between 5 and 29 H/U. These keys and those
of `[solver.coupled]` may change during a run (re-read, `[[schedule]]`).
Iteration counts repeat exactly for the same input, but they are good to
about 10 % only: a change in the last digit of the step can move them
(cavity at Re = 1000: 423 iterations with `step = 0.25`, 452 with the
automatic step, which is the same to 16 digits). The result is the same.

**The cap.** `solver.momentum.relaxation` (0.98) caps the local Courant
number of the momentum equations at 49, `solver.turbulence.relaxation`
(0.98 with a pseudo-time step) that of the turbulence equations. A higher
momentum cap (0.99, 0.999) saves outer iterations on fine meshes and
costs more in each; lower it (0.95) if a run diverges at the start.
Leave the turbulence value: 0.9 makes the turbulence lag the momentum
equations and costs several times the iterations, 1 can diverge from a
uniform start.

**With SIMPLEC.** The turbulence step alone also helps the segregated
algorithm; it is off by default there, so that existing results do not
change:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0.toml \
    --set solver.pseudo_time.turbulence=0.5
```

The coarse step then settles in about 230 iterations at α_u = 0.8
instead of 525.

**Has the result settled?** As with SIMPLEC the residual of the step
case stalls near 10⁻⁴, so judge by a monitored value: write fields at
chosen iterations (`touch output/<case>/.result_50 .result_100 …` before
the run, recipe 74) and compare, or watch a probe in the log
(`[log] columns`, recipe 38). The figure of gallery 78 was drawn from
such histories:

```
python3 -m tools.plot.settling_chart docs/figures/gallery/coupled_step.svg "x_r / H" \
    "k-ε-ζ-f₀:SIMPLEC=simplec_zf0.csv" "k-ε-ζ-f₀:coupled=coupled_zf0.csv" …
```

with one two-column file (iteration, value) per run.

[↑ back to top](#incompressible-flow)

---

<a id="momentum-multigrid"></a>
## 80. Choose and tune the momentum sweep of the coupled algorithm

The preconditioner of the coupled algorithm (recipe 78) applies one
V-cycle of a reduction multigrid to each velocity component
([gallery 80](../gallery/flow.md#momentum-multigrid)). That is the
default; a case needs no key for it:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml
target/release/foehn flow examples/cases/flow/cavity_re1000_coupled.toml
```

**Reading the log and the report.** The columns `U V W` count the
multigrid cycles of the momentum sweeps of an outer iteration: one per
Krylov iteration (`K`) and component with the default settings, none for
a component without a residual (the third one of a planar case). The
report of the run (`flow.md`) has the time of the cycles under
`momentum solve` and that of the setups under `multigrid setup`, with
the number of hierarchies built for the momentum sweeps:

```
multigrid setup 1.78 (45 setups, 28 of the momentum sweeps)
```

**Compare with the old sweep.** BiCGStab with a diagonal preconditioner,
to `momentum_tolerance` or `momentum_max_iter`:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml \
    --set solver.coupled.momentum_sweep=bicgstab
```

The converged result is the same; the path to it and the counts differ
a little (the fine step settles in 175 instead of 150 iterations).

**The keys.**

```toml
[solver.coupled]
momentum_sweep = "multigrid"   # or "bicgstab"
momentum_cycles = 1            # cycles per sweep
momentum_reuse = 0.02          # keep a hierarchy up to this change of the matrix

[solver.coupled.multigrid]
strength = 0.3                 # strong couplings: above this part of the largest of a row
restriction_strength = 0.05    # the same for the entries that enter the restriction
drop_tolerance = 1e-4          # coarse entries below this part of the diagonal are dropped
max_coarse = 64                # unknowns of the coarsest level (dense inverse)
max_levels = 30
precision = "mixed"            # hierarchy in float32; "double": float64
```

- `momentum_reuse`: a hierarchy is kept while the sum of the absolute
  changes of the matrix coefficients since its setup stays below this
  part of the sum of their absolute values. 0 builds one in every outer
  iteration, which costs more than it gains (a setup takes as long as
  twenty to thirty cycles). Above about 0.05 a kept hierarchy stops
  converging; the guard then rebuilds it at once (a kept hierarchy must
  halve the residual in the first sweep of an outer iteration), so a
  large value costs time but not the run.
- `momentum_cycles = 2`: for a run whose `K` stays high although the
  step and the cap are the defaults. On the cases of the gallery it does
  not pay.
- `[solver.coupled.multigrid]`: leave the defaults. A larger
  `drop_tolerance` (10⁻³) gives a smaller hierarchy and is enough with
  the cap 0.98; without the cap it loses the convergence on closed
  streamlines.

The three keys of the sweep may change during a run (re-read,
`[[schedule]]`); those of `[solver.coupled.multigrid]` are fixed for a
run.

**If a setup fails** (a zero on the diagonal, a coarsening that stalls),
the run goes on: BiCGStab does the momentum sweeps of that outer
iteration, and the log carries one warning,
`the multigrid setup of the momentum sweep failed (…); BiCGStab does the
sweeps of such iterations`. Send such a case in: it should not happen.

**Restart.** The restart file holds the matrices of the kept
hierarchies (`momentum_amg_matrix_0` to `_2`); the continued run
rebuilds them and stays bit-identical. A restart file from before this
change works too (the hierarchies are built anew).

**The figure of gallery 80.** The Krylov iterations per outer
iteration come from `residuals.csv`:

```
python3 -m tools.plot.krylov_chart docs/figures/gallery/momentum_multigrid.svg \
    "8 200 cells:BiCGStab sweep=output/step_bicgstab" "8 200 cells:multigrid sweep=output/step_multigrid" \
    "32 800 cells:BiCGStab sweep=output/fine_bicgstab" "32 800 cells:multigrid sweep=output/fine_multigrid"
```

[↑ back to top](#incompressible-flow)

---

<a id="turbulence-multigrid"></a>
## 81. Use the multigrid of the turbulence equations and its larger steps

With the coupled algorithm (recipe 78) in a steady run, the transport
equations of k, ε and ζ of the k-ε family (k-ε, realizable k-ε,
k-ε-ζ-f, k-ε-ζ-f₀) take one cycle of the reduction multigrid of
recipe 80 as the preconditioner of BiCGStab, and the relaxations of the
momentum and turbulence equations default to 0.995 instead of 0.98
([gallery 81](../gallery/flow.md#turbulence-multigrid)). A case needs no
key for it:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml
```

The step settles in about 40 outer iterations and converges to its
tolerance in about 110 (120 with the pressure tolerance 0.1 of the
coupled preconditioner before #82).

**Reading the report.** `flow.md` counts the hierarchies built for the
turbulence equations after those of the momentum sweeps:

```
multigrid setup 1.17 (4 setups, 76 of the momentum sweeps, 114 of the turbulence equations)
```

While the recirculation zone forms, most solves rebuild their
hierarchy; once the flow has settled they are kept. The BiCGStab counts
of k, ε and ζ are in `residuals.csv` (columns `bicgstab_k`, `bicgstab_epsilon`, `bicgstab_zeta`): 1 to 2
per solve with the multigrid, up to the limit 50 with Jacobi at 0.995.

**Back to the defaults of part 6f** (Jacobi and 0.98), for a comparison
or a case that does not take the larger steps:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml \
    --set solver.turbulence.preconditioner=jacobi \
    --set solver.turbulence.relaxation=0.98 --set solver.momentum.relaxation=0.98
```

Setting `preconditioner = "jacobi"` alone keeps the relaxations at 0.98
too; setting a relaxation keeps the others at their defaults. A
relaxation above 0.99 with Jacobi runs the solves into their limit and
the result jitters: do not combine them.

**The keys.**

```toml
[solver.turbulence]
preconditioner = "reduction"   # default with the coupled algorithm and the k-ε family;
                               # "jacobi" elsewhere; "amg" is the aggregation multigrid
                               # of recipe 73 (not with the coupled algorithm: its
                               # coarsest level is not positive definite there)
reuse = 0.02                   # keep a hierarchy up to this change of the matrix
relaxation = 0.995             # default with "reduction", else 0.98

[solver.momentum]
relaxation = 0.995             # the cap on the local step; default with "reduction", else 0.98
```

- `reuse`: as `momentum_reuse` of recipe 80, the sum of the absolute
  changes of the matrix since the setup over the sum of its absolute
  values; a hierarchy whose solve needs more than 8 BiCGStab iterations
  is rebuilt in the next solve anyway. 0 builds one in every solve.
- The multigrid settings are those of `[solver.coupled.multigrid]`.

**k-ω SST** keeps Jacobi and 0.98: with accurate solves of k and ω its
steady state on the step is not held (gallery 81). `--set
solver.turbulence.preconditioner=reduction` on an SST case is accepted
but not recommended.

**If a setup fails**, the solves of that equation take Jacobi until a
setup succeeds, and the log carries one warning, `the multigrid setup
of a turbulence equation failed (…); its solves take the Jacobi
preconditioner until a setup succeeds`.

**Restart.** The restart file holds the matrices of the kept
hierarchies (`turb_amg_matrix_<slot>`); the continued run stays
bit-identical. A restart file from before this change works too (the
hierarchies are built anew).

**The figure of gallery 81.** Write x_r every 5 iterations (recipe 78)
into a two-column file (iteration, value) per run, then:

```
python3 -m tools.plot.settling_chart docs/figures/gallery/turbulence_multigrid.svg "x_r / H" \
    "8 200 cells:part 6f (relaxations 0.98)=coarse_6f.csv" "8 200 cells:part 6g (0.995, multigrid)=coarse_6g.csv" \
    "32 800 cells:part 6f (relaxations 0.98)=fine_6f.csv" "32 800 cells:part 6g (0.995, multigrid)=fine_6g.csv"
```

[↑ back to top](#incompressible-flow)

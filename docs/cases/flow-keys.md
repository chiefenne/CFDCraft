# Flow case-file keys

Part of the [Case catalogue](../cases.md): every key of a flow case, the boundary conditions (D29) and their sensible combinations (D30).

`[mesh]`, `[mesh.generate]` as in the [mesh cases](meshes.md#mesh-case-file-keys).
Unknown keys and values are errors with a "did you mean" hint (D39). SI
units. The complete list, generated from the code with types, defaults,
units and whether a key may change during a run, is
[keys.md](../keys.md); `--set key=value` overrides any key from the command
line (recipe 26).

| Key | Default | Meaning |
|---|---|---|
| `fluid.density` | 1 | ρ (kg/m³), or from `fluid.material` |
| `fluid.viscosity` | 1e-3 | dynamic viscosity μ (Pa s), or from `fluid.material`; an expression of x, y, z gives a viscosity that varies in space |
| `fluid.material` | none | a material of `[materials]` or of the database (`foehn materials`), below |
| `fluid.temperature` | 293.15 | the temperature at which the material properties are evaluated (K; isothermal until the energy equation) |
| `fluid.operating_pressure` | 101325 | the pressure at which they are evaluated (Pa) |
| `flow.convection` | true | false: Stokes flow |
| `flow.body_force` | `[0, 0, 0]` | body force per unit volume (N/m³); each component a number or an expression of x, y, z, t |
| `initial.velocity` | `[0, 0, 0]` | start velocity (m/s) |
| `initial.pressure` | 0 | start pressure (Pa) |
| `initial.potential_flow` | false | start from potential flow instead of `initial.velocity` (D52): the irrotational velocity that carries the inflow to the pressure boundaries, one Laplace solve with the pressure multigrid and four non-orthogonal correction passes |
| `pressure.reference_point` | lower corner of the bounding box | the cell containing it holds the pressure level when no boundary fixes it |
| `pressure.reference_value` | 0 | pressure in that cell (Pa) |
| `boundary.<patch>` | required for every patch | condition, table below |
| `solver.algorithm` | `"simplec"` | `simplec` (Van Doormaal & Raithby 1984) or `simple` |
| `solver.scheme` | `"linear-upwind"` | convection of momentum: `upwind`, `linear-upwind` or `central` (bounded central, D25) |
| `solver.central_blending`, `solver.central_bounding` | 1, 0.1 | with `central`: share of the central correction (rest linear upwind) and the bound β (as for the scalar cases) |
| `solver.limiter` | `"none"` | gradient limiter: `none`, `venkatakrishnan`, `venkatakrishnan-wang` |
| `solver.limiter_k` | 5 | constant K of the limiter |
| `solver.nonorth_correction` | true | deferred non-orthogonal and skewness correction of the viscous term |
| `solver.relaxation_velocity` | 0.7 | under-relaxation α_u of the momentum equations (implicit, Patankar 1980) |
| `solver.relaxation_pressure` | 1 (SIMPLEC), 0.3 (SIMPLE) | relaxation α_p of the pressure correction |
| `solver.max_iterations` | 2000 | largest number of outer iterations |
| `solver.tolerance` | 1e-8 | converged when every scaled residual (momentum, continuity) is below it |
| `solver.momentum_tolerance` | 1e-2 steady, 1e-8 transient | relative tolerance of each momentum solve (BiCGStab, Jacobi) |
| `solver.momentum_max_iter` | 50 | largest number of BiCGStab iterations per solve |
| `solver.pressure_tolerance` | 1e-2 steady, 1e-8 transient | relative tolerance of each pressure-correction solve |
| `solver.pressure_max_iter` | 100 | largest number of CG iterations per solve |
| `solver.pressure_preconditioner` | `"amg"` | `amg` (stage 3) or `jacobi` |
| `solver.amg.*` | as `linear.amg.*` of the [Poisson cases](poisson.md#poisson-case-file-keys) | multigrid settings |
| `solver.amg_rebuild` | 10 | rebuild the hierarchy every N outer iterations |
| `solver.log_every` | 10 | a log line every N outer iterations |
| `analytic.type` | none | exact solution for the checks: `uniform`, `couette`, `poiseuille` or `pipe-poiseuille` |
| `analytic.lower`, `analytic.upper` | 0, 1 | channel walls (y) |
| `analytic.wall_velocity` | 1 | Couette: speed of the upper wall |
| `analytic.mean_velocity` | 1 | Poiseuille: mean velocity |
| `analytic.x_outlet`, `analytic.pressure` / `pressure_outlet` | 1, 0 | Poiseuille: outlet position and pressure; Couette: the constant pressure |
| `verification.levels` | none | refinement levels of `[mesh.generate]`: an order-of-accuracy sweep |
| `verification.expected_order` | none | accepted range of the observed order on the finest pair |
| `verification.max_error` | none | largest accepted relative L2 velocity error |
| `verification.reference` | none | `ghia-re100` or `ghia-re1000`: centre-line data of the unit cavity |
| `verification.reference_tolerance` | none | largest accepted deviation, in lid speeds |
| `verification.lid_velocity` | 1 | lid speed for the reference comparison and the Reynolds number |
| `verification.require_convergence` | none | fail unless the outer iterations converge |
| `run.threads` | all cores | CPU threads (meshes below 20 000 cells run on one; the results are bitwise the same) |
| `run.backend` | `cpu` | where the flow loop runs: `cpu`, or `cuda` in a build with the `cuda` feature (NVIDIA GPU of compute capability 7.0 or newer); the results are bitwise those of `cpu` (D31). `foehn flow --backend` overrides it |
| `run.ptx` | `target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx` | PTX module of the `cuda` backend (`gpu/ptx-kernels/build.sh`); `foehn flow --ptx` overrides it |
| `output.directory` | `output/<case name>` | output directory |
| `output.vtkhdf` | true | write `fields.vtkhdf` |
| `output.history_every` | 1 | a row of `residuals.csv` every N outer iterations |
| `output.wall_csv` | true | write `wall_<patch>.csv` for every wall |
| `output.lines` | none | sampling lines `[[output.lines]]` with `name`, `from`, `to`, `points` (101) |
| `output.fields_every` | 0 | transient: write `fields_<step>.vtkhdf` every N steps (0: only the final `fields.vtkhdf`) |
| `output.fields_at` | none | transient: also write the fields at these steps |
| `output.fields_interval`, `output.fields_times` | none | transient: also write the fields at the multiples of this simulated time (s) and at these times; the steps land on them exactly (D66) |
| `output.monitor_every` | 1 | a row of `forces.csv` and `probes.csv` every N steps |
| `output.monitor_interval` | none | transient: a row of `forces.csv` and `probes.csv` at the multiples of this simulated time (s) instead of every `monitor_every` steps, landed on exactly (equally spaced samples for spectra, D66) |
| `output.probes` | none | point probes `[[output.probes]]` with `name` and `point`: velocity and pressure in `probes.csv` |
| `output.coefficients` | none | force coefficients C = 2F/(ρU²A) of the sum over `patches`, with `velocity` U, `area` A, `drag_direction` (x), `lift_direction` (y) |
| `output.derived` | none | derived fields in every `fields*.vtkhdf` (D72): `vorticity`, `q-criterion`, `lambda2`, `enstrophy`, `strain-rate`, `pressure-coefficient`, `total-pressure`, `wall-shear-stress`, `y-plus`; see [Derived fields](#derived-fields-d72) |
| `output.reference.velocity`, `.pressure`, `.density` | `[output.coefficients] velocity`, 0, `[fluid] density` | reference state U_ref (m/s), p_ref (Pa), ρ_ref (kg/m³) of `pressure-coefficient` |

Transient cases (stage 5):

| Key | Default | Meaning |
|---|---|---|
| `time.scheme` | `"bdf2"` | `bdf1` (backward Euler), `bdf2` (variable steps) or `crank-nicolson`; the first step of BDF2 and Crank–Nicolson is a BDF1 step |
| `time.time_step` | required | Δt (s); with `adaptive` the first step |
| `time.end_time`, `time.steps` | one of them required | end time (s), number of steps |
| `time.start_time` | 0 | time of the initial fields (s) |
| `time.adaptive` | false | adapt Δt to the largest Courant number `max_cfl` |
| `time.max_cfl` | 1 | target Courant number with `adaptive` |
| `time.max_growth` | 1.2 | largest ratio of successive steps, below 1 + √2 (zero stability of BDF2) |
| `time.min_time_step`, `time.max_time_step` | 0, none | bounds of the adaptive Δt (s) |
| `time.error_tolerance` | none | with `adaptive` and `bdf2` or `crank-nicolson`: relative tolerance of the local time-integration error per step (D65); the step follows from the error estimate, inside `max_cfl` and `max_growth`, and a step above the tolerance is repeated shorter (at most 8 times, not below `min_time_step`). Keep it at least ten times above the linear-solver tolerances |
| `time.error_absolute` | 1e-6 | absolute tolerance of the error estimate (m/s), added to `error_tolerance` × the RMS velocity |
| `solver.algorithm` | `"piso"` with `[time]` | PISO (Issa 1986) |
| `solver.correctors` | 2 | pressure corrections per outer iteration (one diverges on meshes with high diffusion numbers, such as the cylinder) |
| `solver.outer_iterations` | 1 | outer iterations per time step (1: plain PISO; more: PIMPLE-type, for large steps) |
| `solver.outer_tolerance` | 1e-6 | ends the outer iterations of a step when the residuals at the start of one are below it |
| `solver.extrapolation` | true | extrapolate the convecting flux, velocity gradients and limiters to the new time in the first outer iteration (keeps plain PISO second order) |
| `solver.relaxation_velocity` | 1 transient | below 1 needs `outer_iterations` > 1 |
| `run.restart` | none | restart file to continue from (D42) |
| `run.restart_every` | 0 | write `restart.h5` every N steps or iterations |
| `run.restart_at` | none | also write `restart_<N>.h5` at these steps |
| `run.restart_interval`, `run.restart_times` | none | transient: write `restart.h5` at the multiples of this simulated time (s), and a numbered `restart_<N>.h5` at these times; landed on exactly (D66) |
| `run.restart_final` | true | write `restart.h5` at the end |
| `run.reread_every` | 0 | re-read the case file every N steps or iterations (0: only on `.read_toml_N`), D40 |
| `run.limits.max_wall_time` | none | stop after this wall time (s) with the restart file and results, exit code 3 (D83) |
| `run.limits.max_steps`, `max_iterations` | none | stop after this many time steps (transient) or outer iterations (steady) |
| `run.limits.max_output` | none | stop when the output directory holds this many bytes |
| `balance.every` | 0 | the full balance table (mass, kinetic energy per patch and term) in the log every N iterations or steps; 0: at the end only (D89) |
| `balance.tolerance_step` | 1e-8 | relative mass imbalance of a time step or of the final steady state before the warning `BAL001` |
| `balance.tolerance_run` | 1e-4 | relative mass imbalance over the run (`BAL002`) |
| `balance.tolerance_energy` | 1e-2 | negative numerical dissipation of kinetic energy, relative to the largest term (`BAL003`) |
| `balance.action` | `"warn"` | `"stop"`: a transient run stops at the step with exit code 3 |
| `verification.time_steps` | none | a time-order sweep on one mesh: three or more steps, largest first |
| `verification.expected_time_order` | none | accepted range of the observed time order on the smallest steps |
| `verification.benchmark` | none | `schafer-turek-2d1` (steady, Re 20) or `schafer-turek-2d2` (periodic, Re 100) |
| `verification.benchmark_patch` | `"cylinder"` | patch of the cylinder |
| `verification.benchmark_tolerance` | none | largest accepted distance from the reference ranges, relative to their midpoints |
| `verification.evaluation_start` | half the run | start of the evaluation window of 2D-2 (s) |
| `analytic.type = "taylor-green"` | | the decaying Taylor–Green vortex with `velocity` U (1), `length` L (1, the period) and `pressure` p₀ (0) |
| `analytic.type = "uniform"` | | uniform flow with `velocity` ([1, 0, 0]), `pressure` (0) and `length` (1, for the Reynolds number): exact between slip walls |
| `analytic.type = "pipe-poiseuille"` | | Hagen–Poiseuille flow along the x axis with `radius` R (1), `mean_velocity` (1), `x_outlet` (1) and `pressure_outlet` (0) |

Built-in cylinder mesh `[mesh.cylinder]` (D34; patches `inlet`, `outlet`,
`walls`, `cylinder`, `zmin`, `zmax`): an O-grid around the cylinder that
fills the square 0 ≤ x, y ≤ height and a graded wake block to the outlet.

| Key | Default | Meaning |
|---|---|---|
| `length`, `height`, `depth` | 2.2, 0.41, 0.01 | channel length, height and the layer thickness in z (m) |
| `center`, `diameter` | [0.2, 0.2], 0.1 | cylinder centre (x, y) and diameter (m) |
| `n` | 32 | intervals per side of the O-grid (4n around the cylinder) |
| `radial`, `radial_ratio` | 32, 6 | intervals from the cylinder to the square, ratio of the last to the first |
| `wake`, `wake_ratio` | 64, 8 | intervals in x of the wake block, ratio of the last to the first |
| `refine` | 0 | halve every interval this many times (refine 0: 6 144 cells, each level four times more) |

# Derived fields (D72)

`[output] derived` adds fields computed from the solution to every
`fields*.vtkhdf`; `foehn post derive CASE` writes the same file afterwards
from a restart file ([recipe 30](../recipes/post-processing.md#derived-fields)).
G is the cell velocity gradient (Gᵢⱼ = ∂uᵢ/∂xⱼ, the weighted
least-squares gradient of the solver), S = ½ (G + Gᵀ), Ω = ½ (G − Gᵀ),
‖·‖ the Frobenius norm.

| Name | Datasets | Definition | Unit |
|---|---|---|---|
| `vorticity` | `vorticity_x`, `_y`, `_z`, `vorticity_magnitude` | ω = ∇ × u | 1/s |
| `q-criterion` | `q_criterion` | ½ (‖Ω‖² − ‖S‖²), vortex cores where Q > 0 | 1/s² |
| `lambda2` | `lambda2` | middle eigenvalue of S² + Ω², vortex cores where λ₂ < 0 | 1/s² |
| `enstrophy` | `enstrophy` | ½ ω·ω | 1/s² |
| `strain-rate` | `strain_rate` | √(2 S:S) | 1/s |
| `pressure-coefficient` | `pressure_coefficient` | (p − p_ref)/(½ ρ_ref U_ref²); needs a reference velocity | - |
| `total-pressure` | `total_pressure` | p + ½ ρ u·u | Pa |
| `wall-shear-stress` | `wall_shear_stress_x`, `_y`, `_z`, `wall_shear_stress_magnitude`, on the patch blocks | tangential viscous stress of the fluid on the wall (as `wall_<patch>.csv`), zero off walls | Pa |
| `y-plus` | `y_plus`, on the patch blocks | d √(τ_w/ρ)/ν, d the wall distance of the cell centroid, zero off walls | - |

The datasets are written in the order of this table, whatever the order of
the list. `output.derived` and `output.reference` may change during a run.

# Boundary conditions (D29)

| `type` | Keys | Velocity | Pressure | Mass flux |
|---|---|---|---|---|
| `wall` | `velocity` (default 0; only the tangential part is kept), `mass_flux` (kg/(m² s) into the domain through the wall, D124 W3: blowing > 0, suction < 0; a number or an expression) | prescribed (with the transpiration velocity m''/ρ along −n) | extrapolated | 0, or −m''·\|S\| |
| `slip-wall` | none | no shear stress, zero normal component (as `symmetry`, on any surface) | extrapolated | 0 |
| `velocity-inlet` | one of `velocity`, `normal_velocity` (into the domain), `profile = "analytic"` (the exact solution, in a transient case at each time level), `profile = "parabolic"` with `max_velocity` and `profile_axes` (developed laminar inflow, for example `["y"]`), `profile = "table"` with `table` and `table_coordinate` (a measured or computed profile, below) | prescribed | extrapolated | ρ u_b·S |
| `mass-flow-inlet` | `mass_flow` (kg/s, spread by area) | normal, from the mass flow | extrapolated | prescribed |
| `static-pressure-inlet` | `pressure` | normal to the face, from the flux | prescribed | Rhie–Chow |
| `total-pressure-inlet` | `total_pressure` | normal to the face, from the flux | p₀ − ½ρu_n² | Rhie–Chow |
| `static-pressure-outlet` | `pressure`, `backflow` | extrapolated; on backflow by `backflow` | prescribed; on backflow by `backflow` | Rhie–Chow |
| `total-pressure-outlet` | `total_pressure`, `backflow` | extrapolated; on backflow by `backflow` | p₀ − ½ρu_n² | Rhie–Chow |
| `convective-outlet` | `convection_velocity` U_c (default: the mean outflow speed) | ∂u/∂t + U_c ∂u/∂n = 0 (transient cases) | extrapolated | ρ u_b·S, balanced with the inflow |
| `symmetry` | none | no shear stress, zero normal component | extrapolated | 0 |
| `wedge` | none | as `symmetry`: the two side planes of an axisymmetric wedge mesh ([`mesh.generate.wedge`](meshes.md)) | extrapolated | 0 |
| `periodic` | `partner` | comes with stage 6 (the face pairing); an error until then | | |

**Backflow at pressure outlets** (`backflow`). Where the flow enters
through an outlet face (a vortex crossing the outlet, a recirculation
reaching it):

| `backflow` | Velocity of the entering fluid | Pressure | Suits |
|---|---|---|---|
| `"total-pressure"` (default) | normal to the face | the prescribed value as the total pressure, p − ½ρu_n² | outlets that vortices or a recirculation cross; converges like an outflow |
| `"zero-gradient"` | extrapolated | the prescribed value as the static pressure | the treatment before stage 5; converges slowly where the flow enters |
| `"suppress"` | the face is closed for the iteration: a slip wall with the extrapolated pressure | | outlets where nothing may enter; a closed face opens when the prescribed pressure drives the flow out |

**Convective outlet.** The velocity leaves the domain by the convection
equation ∂u/∂t + U_c ∂u/∂n = 0, discretised with the backward difference
of the run between the face and its cell, so vortices pass the outlet with
little reflection and the outlet can sit close behind a body. It needs a
transient case. The outlet flow is balanced with the prescribed inflow in
every iteration; without a pressure boundary the reference cell fixes the
pressure level. U_c defaults to the mean outflow speed (the prescribed
inflow over the outlet area); give `convection_velocity` when nothing
prescribes an inflow.

**Inlet profile from a table** (`profile = "table"`). `table` is a file
(relative to the case file) of rows of a coordinate and either the normal
speed into the domain (two columns) or the three velocity components (four
columns); `#` starts a comment, columns are separated by commas,
semicolons or blanks, a first line that is not numeric is a header, and
the coordinate increases. Between the rows the values are linear, beyond
the ends constant (the check warns, BC020). `table_coordinate` is `"x"`,
`"y"`, `"z"` or `"radius"`, the distance from the axis through
`table_origin` (default the origin) along `table_direction` (needed). The
table's sha256 is part of the case hash (D85), and the table must lie
inside the read root of the run (D82).

**Wedge** (axisymmetric flow). `mesh.generate.wedge = θ` (degrees, up to
90) revolves a 2D mesh of the half plane y ≥ 0 about the x axis (y is the
radius) into a wedge one cell thick, symmetric about the x–y plane; nodes
on the axis merge, so the cells next to it become prisms or polyhedra with
an edge on the axis and there is no axis patch. `zmin` and `zmax` are the
wedge planes and take `type = "wedge"`. The faces are planar, so a wedge
of θ is a slice of a pipe whose section is a polygon of 360°/θ sides: the
result differs from the axisymmetric one by about θ²/8 (1e-3 at 5°); the
check warns above 10° (BC024).

Inlets also take `turbulence_intensity`, `turbulence_length_scale` and
`viscosity_ratio` for the turbulence models of stage 6, and so do pressure
outlets for the entering fluid; in laminar flow the check reports them as
unused.

# Sensible combinations (D30)

The check runs before every solve and with `--check`. Errors stop the run
and name the boundary and what to change (a missing or unknown patch with
"did you mean", a net mass flow with no pressure boundary to take it, a
flow path whose rate is prescribed twice, a wall with a normal velocity, a
negative mass flow); warnings and notes say what the run will do (no
pressure boundary: the reference cell holds the level; nothing drives a
flow; unused turbulence values; a symmetry patch that is not planar; a
table that is missing or malformed; a convective outlet in a steady case
or without a convection velocity; wedge planes that do not form one
wedge).

| Combination | Suits | Pressure level |
|---|---|---|
| velocity-inlet + static-pressure-outlet (+ walls, symmetry) | internal flows with a known flow rate: channels, pipes, ducts, the Poiseuille case | outlet pressure |
| mass-flow-inlet + static-pressure-outlet | as above when the mass flow is known rather than a velocity profile | outlet pressure |
| total-pressure-inlet + static-pressure-outlet | pressure-driven flows: the flow rate is a result (fans, valves, plenum to plenum) | both |
| static-pressure-inlet + static-pressure-outlet | pressure-driven flows with a known static pressure drop | both |
| static-pressure-inlet + total-pressure-outlet | flows that exhaust into a large volume at a known total pressure | both |
| walls only, one wall moving (+ symmetry) | enclosed flows: the lid-driven cavity | reference cell |
| velocity-inlet + velocity-inlet + static-pressure-outlet | mixing of two streams: the T-junction | outlet pressure |
| velocity-inlet + convective-outlet (transient) | wakes and vortex streets with the outlet close behind the body | reference cell |
| velocity-inlet + static-pressure-outlet + slip walls | external flow in a finite box: the slip walls stand for the far field | outlet pressure |
| velocity-inlet + static-pressure-outlet + wedges (+ walls) | axisymmetric flow: pipes, nozzles, jets | outlet pressure |

# Expressions, materials and regions (D53, D102, D112)

Stage 5 part 6a ([gallery 23](../gallery/properties.md#expressions-materials),
[recipe 31](../recipes/properties.md#expressions),
[recipe 32](../recipes/properties.md#materials)).

**Expressions.** Wherever the table below says so, a value may be a string
of the expression language instead of a number. The variables are `x`,
`y`, `z` (m, the face or cell centroid) and `t` (s, simulated time); the
constants `pi`, `e` and those of `[constants]`; the operators
`+ - * / % ^` (also `**`), `< <= > >= == !=`, `and or not` (true is 1);
the functions `sin cos tan asin acos atan sinh cosh tanh exp log log10
sqrt abs floor ceil round sign` (one argument; `log` is the natural
logarithm), `atan2 pow min max hypot step` (two; `step(a, b)` is 1 where
a ≥ b), `clamp(a, lo, hi)`, `if(condition, a, b)` and every
`[tables.<name>]` as a function of one argument. Each string is checked
and compiled when the case is read; the solver evaluates the compiled
program, on the CPU and the GPU alike, with the in-house elementary
functions of D93 (6), so a result does not depend on the platform's
mathematics library. A number stays a number: a case without expressions
runs bit for bit as before.

| Key | Takes | Meaning |
|---|---|---|
| `constants.<name>` | number or expression of other constants | a named constant; the order does not matter, cycles are errors |
| `tables.<name>.x`, `.y` | numbers | a piecewise linear table, clamped at the ends, called as `name(arg)` |
| `tables.<name>.file` | path | the same from a two-column file (CSV or blanks, `#` comments, header allowed) |
| `boundary.<patch>.velocity`, `normal_velocity`, `mass_flux` | x, y, z, t | inlet and wall values per face |
| `boundary.<patch>.pressure`, `total_pressure` | x, y, z, t | pressure boundaries per face |
| `boundary.<patch>.mass_flow` | t | one value for the patch |
| `flow.body_force` | x, y, z, t | per cell; an expression is added to the momentum source, a uniform number keeps the kernel path |
| `fluid.viscosity`, `materials.<name>.viscosity` | x, y, z | per face (properties do not change during a run, D100) |

`t` is allowed only in a transient case; time-dependent values are
evaluated at the new time level of every step.

**Materials.** `fluid.material` names a material of the case's
`[materials]` or of the built-in database; a `[materials.<name>]` table
with `base = "<database name>"` overrides single properties of a database
entry. Each property is a number, an expression or a table:

| Key | Meaning |
|---|---|
| `materials.<name>.base` | the material it starts from |
| `materials.<name>.kind` | `gas`, `liquid` or `solid` |
| `materials.<name>.temperature_range` | `[T_min, T_max]` (K): outside it the properties are clamped, with warning MAT003 |
| `materials.<name>.molar_mass` | kg/mol (the ideal-gas density) |
| `materials.<name>.density`, `.viscosity` | a number, an expression, or a table with `form` = `constant` (`value`), `polynomial` (`coefficients` a₀, a₁, … in T), `sutherland` (`mu_ref`, `t_ref`, `s`), `power-law` (`value_ref`, `t_ref`, `exponent`), `table` (`temperature` and `values`, or `file`; `interpolation` = `linear` or `monotone-cubic`), `ideal-gas` (density p M / (R T)), each with `source` and optionally `range`; the complete list is in [keys.md](../keys.md) |

The density must be constant in the domain (a density that varies comes
with the energy equation, MAT002). The run writes the resolved property
program into `case.resolved.toml`, and its hash into the provenance of
every output file (`properties_hash`, `materials_database`).

**Regions.** `[regions.<name>]` assigns cell zones to a material:

| Key | Meaning |
|---|---|
| `regions.<name>.kind` | `fluid`; `solid` and `porous` are refused until the energy equation (REG001) |
| `regions.<name>.zones` | cell zones of the mesh (each zone in one region; unknown zones are REG002) |
| `regions.<name>.material` | its material (default: that of `[fluid]`) |

Regions with different viscosities give a viscosity per face, the harmonic
mean across a zone boundary.


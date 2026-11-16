# Case catalogue

The case files in `examples/cases/` with what each one is, the mesh it needs
and where it is discussed. Pictures are in the [Example gallery](gallery.md),
command sequences in the [Recipes](recipes.md).

How to read the **Mesh** column:

* **synthetic** – the mesh is generated in memory from the `[mesh]` section
  of the case file; no mesh file is needed.
* **generated** – the built-in test mesher makes the mesh from the
  `[mesh.generate]` section (D34).
* **Gmsh** – the mesh is read from a Gmsh MSH 4.1 file (D33); the `.geo`
  source is in `examples/meshes/`, the `.msh` is made from it with Gmsh
  and not stored in the repository.

## Meshes (`examples/cases/mesh/`)

Stage 1: build or read a mesh, check it, renumber and partition it, write
the quality report (`mesh_report.md`, `mesh_report.json`) and `mesh.vtkhdf`
for ParaView. Run with `foehn mesh FILE`; a Gmsh file can be given
directly instead of a case file (`foehn mesh pipe.msh`), command-line
options override the file (`foehn mesh --help`).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `mesh/hex_box.toml` | generated, 576 hexahedra, graded, perturbed | hexahedral test mesh, quality report | gallery 3 |
| `mesh/tet_box.toml` | generated, 3 456 tetrahedra | tetrahedral test mesh | gallery 3 |
| `mesh/poly_box.toml` | generated, 819 polyhedra | polyhedral test mesh (dual of tetrahedra) | gallery 3 |
| `mesh/mixed_box.toml` | generated, 1 728 mixed cells | hexahedra, prisms, tetrahedra and transitional polyhedra in one mesh | gallery 3 |
| `mesh/hex_2d.toml`, `tet_2d.toml`, `poly_2d.toml`, `mixed_2d.toml` | generated, 2D (one cell thick) | the four types in 2D, patch renaming | gallery 5 |
| `mesh/poly_partition.toml` | generated, 2 197 polyhedra | METIS partition into 4 parts, halo layout, point location | gallery 4, recipe 9 |
| `mesh/tjunction_tet.toml` | Gmsh, `meshes/tjunction.geo`, 77 776 tetrahedra | Gmsh reader, physical groups as patches and zones, RCM renumbering, 4 parts | gallery 6, recipe 8 |
| `mesh/tjunction_poly.toml` | Gmsh, the same, converted into 15 927 polyhedra | conversion of tetrahedra into polyhedra with the exact volume | gallery 6, recipe 8 |

## Scalar cases (`examples/cases/scalar/`)

Stage 2: steady convection–diffusion of a scalar, ∇·(ρuφ) − ∇·(Γ∇φ) = s,
with the operators of stage 2 ([explanation](verification/stage2.md)). Run
with `foehn scalar FILE`; command-line options override the file
(`foehn scalar --help`). Each run writes `scalar.md` and `scalar.json`
(per solve: iterations, residual, error norms, φ range, time per phase, ns
per cell and iteration, comparison with the first backend), `fields.vtkhdf`
(φ, limiter, error), `convergence.svg` for a sweep and `profile.svg` for a
sampling line. The run fails when a check of the case fails (observed
order, bounds, convergence, backends differing by more than 1e-10). CI runs
every case except `*_gpu.toml` (R5).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `scalar/mms_poly_box.toml` | generated poly 3D, 729 to 35 937 cells (3 levels) | manufactured solution, order of accuracy (1.8–2.2 required) | gallery 8, recipe 12 |
| `scalar/mms_hex_box.toml` | generated hex 3D, 512 to 32 768 cells | the same on perturbed (warped) hexahedra | recipe 12 |
| `scalar/mms_tet_box.toml` | generated tet 3D, 3 072 to 196 608 cells | the same on tetrahedra | recipe 12 |
| `scalar/mms_mixed_box.toml` | generated mixed 3D, 1 344 to 95 232 cells | the same on hexahedra, prisms, tetrahedra and polyhedra | recipe 12 |
| `scalar/mms_poly_2d.toml` | generated poly 2D, 289 to 16 641 cells (4 levels) | planar solution, zero-gradient (symmetry) z patches | recipe 12 |
| `scalar/oblique_step_2d.toml` | generated poly 2D, 2 401 cells | upwind, linear upwind, both limiters and the bounded central scheme (stage 5, D25) as `[[variants]]`; bounds check | gallery 9, recipe 13 |
| `scalar/mms_poly_gpu.toml` | generated poly 3D, 912 673 cells | CPU against two GPU runs: GPU equals CPU and GPU run to run | recipe 14, `tools/stage2_colab.ipynb` |
| `scalar/poisson_amg_poly_box.toml` | generated poly 3D, 729 to 35 937 cells | Poisson manufactured solution with CG + AMG (stage 3) against BiCGStab + Jacobi as `[[variants]]`; order of accuracy | gallery 11, recipe 15 |

### Scalar case-file keys

`[mesh]` and `[mesh.generate]` as in the [mesh cases](#mesh-case-file-keys).
A boundary condition is required for every patch.

| Key | Default | Meaning |
|---|---|---|
| `problem.density` | 1.0 | ρ (kg/m³) |
| `problem.diffusivity` | 0.05 | Γ (kg/(m s)) |
| `problem.velocity` | `[1.0, 0.6, 0.3]` | u (m/s), constant |
| `problem.source` | 0 | constant volumetric source |
| `problem.manufactured` | none | `{}` for the default manufactured solution, or `{ offset, amplitude, wavenumber, phase }`; sets the source, the `manufactured` boundary values and the error norms |
| `problem.boundary.<patch>.type` | – (required) | `fixed-value`, `fixed-gradient` (outward normal gradient) or `zero-gradient` (also symmetry and outflow) |
| `problem.boundary.<patch>.profile` | `"constant"` | `constant` (`value`), `step` (`axis`, `at`, `below`, `above`), `linear` (`value`, `gradient`), `manufactured` |
| `solver.scheme` | `"linear-upwind"` | `upwind`, `linear-upwind` or `central` (bounded central, D25) |
| `solver.central_blending` | 1 | with `central`: share c of the bounded central correction; the rest, 1 − c, is linear upwind |
| `solver.central_bounding` | 0.1 | with `central`: β of the normalised-variable bound (central where φ̃_U ≥ β, blended towards upwind below, upwind outside (0, 1)); 0 for plain central |
| `solver.limiter` | `"none"` | `none`, `venkatakrishnan`, `venkatakrishnan-wang` |
| `solver.limiter_k` | 5.0 | K: ε² = (K h)³ (Venkatakrishnan) or ε = K (φ_max − φ_min) (Venkatakrishnan–Wang); 0 for strictly bounded |
| `solver.nonorth_correction` | true | non-orthogonal diffusion correction |
| `solver.nonorth_limit` | −1 | limit of that correction as a multiple of the orthogonal part (negative: none) |
| `solver.relaxation` | 1.0 | implicit under-relaxation factor |
| `solver.max_outer` | 200 | outer (deferred-correction) iterations |
| `solver.tolerance` | 1e-10 | outer tolerance on ‖b − Aφ‖/‖b‖ |
| `solver.linear_max_iter` | 200 | linear iterations per outer iteration |
| `solver.linear_tolerance` | 1e-3 | relative tolerance of each linear solve |
| `solver.linear_solver` | `"bicgstab"` | `bicgstab` (any matrix) or `cg` (flexible CG, symmetric positive definite matrices: no convection) |
| `solver.preconditioner` | `"jacobi"` | `jacobi` or `amg` (one V-cycle of the smoothed-aggregation multigrid, stage 3) |
| `solver.amg` | defaults | multigrid settings, the keys of [`linear.amg`](#poisson-case-file-keys) |
| `run.backends` | `["cpu"]` | `cpu`, `cpu-serial`, `cuda`; later backends are compared with the first |
| `run.threads` | all cores | CPU threads |
| `run.rust_ptx` | `target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx` | PTX module of the GPU backend |
| `verification.levels` | none | refinement levels of `[mesh.generate]`; two or more: order-of-accuracy sweep |
| `verification.expected_order` | none | `[lo, hi]`: fail if the order on the finest pair is outside |
| `verification.bounds` | none | `[lo, hi]`: fail if φ leaves this range (by more than 1e-6) |
| `[[variants]]` | none | `name`, `solver = { ... }`, optional `bounds`: the case once per variant |
| `output.directory` | `output/<case name>` | output directory |
| `output.vtkhdf` | true | write `fields.vtkhdf` |
| `output.profile` | none | `{ from, to, points }`: sample φ along a line (`profile.csv`, `profile.svg`) |

## Poisson studies (`examples/cases/poisson/`)

Stage 3: the pressure solver on the Poisson equation −∇·(Γ∇φ) = s
([explanation](verification/stage3.md)). Run with `foehn poisson FILE`;
command-line options override the file (`foehn poisson --help`). Each
level is assembled once and solved from φ = 0 by every solver of
`study.solvers` on every backend of `study.backends` (hypre once per
level, on its own device). Each run writes `poisson.md` and `poisson.json`
(per solve: iterations, true final residual, setup and solve time, ns per
cell and iteration, hierarchy sizes and complexities, difference to the
reference, bitwise comparison with the first backend) and
`iterations.svg` for a sweep. The run fails when the iterations of the
checked solvers grow by more than `study.max_iteration_ratio`, when a
solution differs from the reference by more than `study.max_difference`,
when a solver differs between backends, or when a solve does not converge.
CI runs every case except `*_gpu.toml`, and again with hypre (R5).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `poisson/poisson_poly_box.toml` | generated poly 3D, 729 to 274 625 cells (4 levels) | CG + AMG (mixed and float64) and Jacobi-CG, iteration growth ≤ 1.5, 4 threads against 1 | gallery 10, recipe 15 |
| `poisson/poisson_hex_box.toml` | generated hex 3D, 512 to 262 144 cells | the same on perturbed hexahedra | recipe 15 |
| `poisson/poisson_tet_box.toml` | generated tet 3D, 3 072 to 196 608 cells (3 levels) | the same on tetrahedra | recipe 15 |
| `poisson/poisson_mixed_box.toml` | generated mixed 3D, 1 344 to 95 232 cells (3 levels) | the same on hexahedra, prisms, tetrahedra and polyhedra | recipe 15 |
| `poisson/poisson_poly_gpu.toml` | generated poly 3D, 912 673 cells | CPU and GPU against hypre on the GPU: times per cell, GPU equals CPU | recipe 16, `tools/stage3_colab.ipynb` |

### Poisson case-file keys

`[mesh]`, `[mesh.generate]` as in the [mesh cases](#mesh-case-file-keys),
`[problem]` as in the [scalar cases](#scalar-case-file-keys) (the velocity
must be zero).

| Key | Default | Meaning |
|---|---|---|
| `linear.tolerance` | 1e-10 | stop when ‖b − Aφ‖ ≤ tolerance ‖b‖ |
| `linear.max_iterations` | 2000 | largest number of CG iterations |
| `linear.amg.aggregation` | `"smoothed"` | `smoothed` (Vaněk, Mandel, Brezina 1996) or `plain` (no prolongator smoothing) |
| `linear.amg.strength` | 0.04 | strength threshold θ on the finest level (halved per level) |
| `linear.amg.smoother` | `"chebyshev"` | `chebyshev` or `jacobi` (damped) |
| `linear.amg.degree` | 2 | Chebyshev degree or Jacobi sweeps, before and after the coarse correction |
| `linear.amg.chebyshev_interval` | `[0.1, 1.1]` | Chebyshev interval as fractions of the estimated λ_max(D⁻¹A) |
| `linear.amg.jacobi_weight` | 4/3 | Jacobi damping as a multiple of 1/λ_max |
| `linear.amg.prolongation_weight` | 4/3 | prolongator smoothing ω as a multiple of 1/λ_max |
| `linear.amg.max_coarse` | 256 | largest coarsest level (dense inverse) |
| `linear.amg.max_levels` | 20 | largest number of levels |
| `linear.amg.power_iterations` | 20 | power iterations for λ_max |
| `linear.amg.precision` | `"mixed"` | `mixed` (float32 hierarchy, float64 CG; D37) or `double` |
| `study.levels` | none | refinement levels of `[mesh.generate]` |
| `study.solvers` | `["amg", "jacobi"]` | `amg` (`linear.amg`), `amg-double`, `amg-plain`, `jacobi`, `hypre` (feature `hypre`) |
| `study.backends` | `["cpu"]` | `cpu`, `cpu-serial`, `cuda` for the in-house solvers |
| `study.threads` | all cores | CPU threads |
| `study.rust_ptx` | `target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx` | PTX module of the GPU backend |
| `study.max_iteration_ratio` | none | fail if iterations on the finest / coarsest level exceed it |
| `study.check_iterations` | `["amg", "amg-double"]` | solvers whose iteration growth is checked |
| `study.max_difference` | 1e-6 | fail if max \|φ − φ_ref\| / max \|φ_ref\| exceeds it (reference: hypre if listed, else the first solve) |
| `study.repeat` | 1 | repetitions of each solve; the fastest is reported |
| `hypre.strong_threshold` | 0.5 | BoomerAMG strength threshold |
| `hypre.coarsen_type` | 10 | 8 PMIS, 10 HMIS |
| `hypre.interp_type` | 6 | 6 extended+i, 14 extended |
| `hypre.relax_type` | 6 on the CPU, 18 on the GPU | 6 symmetric hybrid Gauss–Seidel, 18 l1-Jacobi |
| `hypre.num_sweeps` | 1 | sweeps per level |
| `hypre.agg_num_levels` | 1 | levels of aggressive coarsening |
| `hypre.pmax` | 4 | interpolation entries per row |
| `hypre.max_coarse` | 256 | largest coarsest level |
| `hypre.print_level` | 0 | hypre's own output |
| `output.directory` | `output/<case name>` | output directory |
| `output.vtkhdf` | false | write φ of the finest level (`fields.vtkhdf`) |

## Flow cases (`examples/cases/flow/`)

Stage 4: steady incompressible laminar flow, SIMPLEC with the Rhie–Chow
interpolation ([explanation](verification/stage4.md)). Stage 5: transient
flow with PISO, BDF2 or Crank–Nicolson, when the case has a `[time]` table
([explanation](verification/stage5.md)). Run with
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
[interfaces.md](interfaces.md)), and the
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

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/couette_hex.toml` | generated hex 2D, 128 cells | plane Couette flow with convection: error at the solver tolerance (≤ 1e-9 required) | recipe 18 |
| `flow/couette_tet.toml`, `couette_poly.toml`, `couette_mixed.toml` | generated 2D, perturbed by 15 %, 128 to 256 cells | the same Stokes flow on perturbed tetrahedra (prisms), polyhedra and mixed cells | recipe 18 |
| `flow/poiseuille_hex.toml` | generated hex 2D, perturbed, 512 to 8 192 cells (3 levels) | plane Poiseuille flow, observed order of the velocity error 1.8–2.2 required | gallery 13, recipe 18 |
| `flow/poiseuille_tet.toml`, `poiseuille_poly.toml`, `poiseuille_mixed.toml` | generated 2D, 3 levels | the same on prisms, perturbed polyhedra and mixed cells | gallery 13, recipe 18 |
| `flow/cavity_re100.toml` | generated hex 2D, 64 × 64 | lid-driven cavity, Re = 100, against Ghia et al. (1982), deviation ≤ 2 % of the lid speed | gallery 12, recipe 19 |
| `flow/cavity_re1000.toml` | generated hex 2D, 128 × 128 | the same at Re = 1000 | gallery 12, recipe 19 |
| `flow/cavity_re100_tet3d.toml` | generated tet 3D, 24 576 tetrahedra in a slab with symmetry planes | the cavity on fully three-dimensional tetrahedra, deviation ≤ 3 % | recipe 19 |
| `flow/cavity_re100_poly3d.toml` | generated poly 3D, 5 445 polyhedra | the same on polyhedra | recipe 19 |
| `flow/cavity_re1000_poly3d.toml` | generated poly 3D, 28 227 polyhedra | Re = 1000 on polyhedra, deviation ≤ 3 % | gallery 12, recipe 19 |
| `flow/cavity_re100_budget.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 with `[run.limits]` (wall time, iterations, output volume), run the way a script or an agent runs it; a tighter limit stops it with exit code 3 and `--restart latest` continues bitwise | gallery 18, recipe 26 |
| `flow/cavity_re100_energy.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 with the balance table every 100 iterations: the lid power goes into viscous and numerical dissipation, the numerical part positive (D90 part 4) | gallery 19, recipe 27 |
| `flow/cavity_re100_transient.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 time-stepped (PISO, BDF2, Courant number 6) to its steady state: the SIMPLEC solution, deviation from Ghia ≤ 2 % | recipe 21 |
| `flow/taylor_green_bdf2.toml` | generated hex 2D, 64 × 64 | Taylor–Green vortex, plain PISO with BDF2 over five time steps: time order 1.8–2.3 required | gallery 14, recipe 21 |
| `flow/taylor_green_cn.toml` | generated hex 2D, 64 × 64 | the same with Crank–Nicolson | gallery 14, recipe 21 |
| `flow/taylor_green_space.toml` | generated hex 2D, 16² to 128² (4 levels) | Taylor–Green vortex, space order 1.8–2.2 required | gallery 14, recipe 21 |
| `flow/taylor_green_cfl10.toml` | generated hex 2D, 64 × 64 | Taylor–Green vortex at Courant numbers 1.3 to 9 with outer iterations: time order 2 | gallery 14, recipe 21 |
| `flow/taylor_green_error.toml` | generated hex 2D, 64 × 64 | the time step from the BDF2 error estimate with rejected steps (D65), probe rows every 0.05 s and fields every 0.1 s landed on exactly (D66) | gallery 20, recipe 28 |
| `flow/cylinder_re20.toml` | built-in cylinder mesh, refine 1: 24 576 hexahedra | Schäfer & Turek 2D-1, steady, Re = 20: C_D, C_L, Δp within 0.5 % of the reference ranges | gallery 15, recipe 23 |
| `flow/cylinder_potential.toml` | built-in cylinder mesh, refine 1: 24 576 hexahedra | the potential-flow start (D52) around the cylinder, written without outer iterations | gallery 20, recipe 28 |
| `flow/cylinder_re100.toml` | built-in cylinder mesh, refine 0: 6 144 hexahedra | Schäfer & Turek 2D-2, periodic, Re = 100: C_D,max, C_L,max, St, Δp within 2.5 % (about 8 minutes; Linux, on main, nightly and with the label `full-ci`) | gallery 15, recipe 23 |
| `flow/channel_slip.toml` | generated mixed 2D, perturbed by 20 %, 432 cells | uniform flow between slip walls, exact (error ≤ 1e-9 required) | gallery 21, recipe 29 |
| `flow/poiseuille_table.toml` | generated hex 2D, 512 cells | the inlet profile read from `poiseuille_profile.csv`: the run with the exact parabola | gallery 21, recipe 29 |
| `flow/pipe_wedge.toml` | generated hex 2D revolved into a 0.5° wedge, perturbed, 256 to 4 096 cells (3 levels) | Hagen–Poiseuille flow on the axisymmetric wedge, observed order 1.8–2.2 required | gallery 21, recipe 29 |
| `flow/cylinder_re100_short.toml` | built-in cylinder mesh with the outlet 5.5 D behind the cylinder, 4 928 hexahedra | the Re = 100 vortex street through a static pressure outlet with backflow treatment (Linux, on main, nightly and with the label `full-ci`) | gallery 21, recipe 29 |
| `flow/cylinder_re100_convective.toml` | as `cylinder_re100_short.toml` | the same with a convective outlet (Linux, on main, nightly and with the label `full-ci`) | gallery 21, recipe 29 |
| `flow/cavity_poly3d_gpu.toml` | generated poly 3D, 970 225 polyhedra | the steady flow loop on the GPU: 30 SIMPLEC iterations of the cavity at Re = 100 on the CPU and the GPU, bitwise equal, time per cell and iteration (GPU; not run by CI) | gallery 16, recipe 24 |
| `flow/cylinder_re100_gpu.toml` | built-in cylinder mesh, refine 2: 98 304 hexahedra | the transient flow loop on the GPU: 200 adaptive BDF2 steps of the Re = 100 cylinder on the CPU and the GPU, bitwise equal (GPU; not run by CI) | gallery 16, recipe 24 |

The adjoint-ready kernels of stage 5 (D59) have no case file: their
checks are tests on generated meshes (gallery 17, recipe 25).

### Flow case-file keys

`[mesh]`, `[mesh.generate]` as in the [mesh cases](#mesh-case-file-keys).
Unknown keys and values are errors with a "did you mean" hint (D39). SI
units. The complete list, generated from the code with types, defaults,
units and whether a key may change during a run, is
[keys.md](keys.md); `--set key=value` overrides any key from the command
line (recipe 26).

| Key | Default | Meaning |
|---|---|---|
| `fluid.density` | 1 | ρ (kg/m³) |
| `fluid.viscosity` | 1e-3 | dynamic viscosity μ (Pa s) |
| `flow.convection` | true | false: Stokes flow |
| `flow.body_force` | `[0, 0, 0]` | body force per unit volume (N/m³) |
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
| `solver.amg.*` | as `linear.amg.*` of the [Poisson cases](#poisson-case-file-keys) | multigrid settings |
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

### Boundary conditions (D29)

| `type` | Keys | Velocity | Pressure | Mass flux |
|---|---|---|---|---|
| `wall` | `velocity` (default 0; only the tangential part is kept) | prescribed | extrapolated | 0 |
| `slip-wall` | none | no shear stress, zero normal component (as `symmetry`, on any surface) | extrapolated | 0 |
| `velocity-inlet` | one of `velocity`, `normal_velocity` (into the domain), `profile = "analytic"` (the exact solution, in a transient case at each time level), `profile = "parabolic"` with `max_velocity` and `profile_axes` (developed laminar inflow, for example `["y"]`), `profile = "table"` with `table` and `table_coordinate` (a measured or computed profile, below) | prescribed | extrapolated | ρ u_b·S |
| `mass-flow-inlet` | `mass_flow` (kg/s, spread by area) | normal, from the mass flow | extrapolated | prescribed |
| `static-pressure-inlet` | `pressure` | normal to the face, from the flux | prescribed | Rhie–Chow |
| `total-pressure-inlet` | `total_pressure` | normal to the face, from the flux | p₀ − ½ρu_n² | Rhie–Chow |
| `static-pressure-outlet` | `pressure`, `backflow` | extrapolated; on backflow by `backflow` | prescribed; on backflow by `backflow` | Rhie–Chow |
| `total-pressure-outlet` | `total_pressure`, `backflow` | extrapolated; on backflow by `backflow` | p₀ − ½ρu_n² | Rhie–Chow |
| `convective-outlet` | `convection_velocity` U_c (default: the mean outflow speed) | ∂u/∂t + U_c ∂u/∂n = 0 (transient cases) | extrapolated | ρ u_b·S, balanced with the inflow |
| `symmetry` | none | no shear stress, zero normal component | extrapolated | 0 |
| `wedge` | none | as `symmetry`: the two side planes of an axisymmetric wedge mesh (`mesh.generate.wedge`) | extrapolated | 0 |
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

### Sensible combinations (D30)

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

## Benchmarks (`examples/cases/bench/`)

Gate G1, the kernel benchmark of stage 0 ([G1](benchmarks/g1.md)). Run with
`foehn bench g1 FILE`; every key is optional, command-line options
override the file (`foehn bench g1 --help`).

| Case file | Backends | Mesh | Shows | Where |
|---|---|---|---|---|
| `bench/g1_ci_small.toml` | CPU | synthetic, 0.2 M cells, shuffled | smoke test of the harness: all backends bitwise equal; run by CI on every platform | CI, recipe 1 |
| `bench/g1_cpu_10m.toml` | CPU | synthetic, 10 M cells, natural | the CPU half of G1: single-source bodies against hand-written loops | gallery 1, G1 results |
| `bench/g1_cpu_10m_shuffled.toml` | CPU | synthetic, 10 M cells, shuffled | worst-case memory access, the cost of a badly numbered mesh | gallery 1 |
| `bench/g1_gpu_10m.toml` | GPU and CPU | synthetic, 10 M cells, natural | **the G1 verdict**: single-source Rust PTX against hand-written CUDA C | recipe 3, G1 procedure |
| `bench/g1_gpu_reference_4m.toml` | GPU (reference only) and CPU | synthetic, 4 M cells, natural | the CUDA path on a small or pre-Volta GPU: GPU results equal the CPU bitwise; no verdict | recipe 2, G1 maintainer's machine |
| `bench/g1_gpu_50m.toml` | GPU | synthetic, 50 M cells, natural | G1 at the upper end of the stage-plan range (17 GB GPU memory) | recipe 3 |
| `bench/g1_cpu_poly_2m.toml` | CPU | generated, 2 M polyhedra, natural | G1 rerun on polyhedra (stage 1): uneven rows, 14 faces per cell | gallery 7, recipe 10 |
| `bench/g1_cpu_poly_2m_shuffled.toml` | CPU | generated, 2 M polyhedra, shuffled | the same, worst-case ordering | gallery 7, recipe 10 |
| `bench/g1_cpu_tet_2m.toml` | CPU | generated, 2 M tetrahedra, natural | G1 rerun on tetrahedra | gallery 7, recipe 10 |
| `bench/g1_gpu_poly_4m.toml` | GPU and CPU | generated, 4 M polyhedra, natural | the GPU verdict rerun on polyhedra (fits a Colab T4) | recipe 10 |

### Case-file keys

| Key | Default | Meaning |
|---|---|---|
| `mesh.type` | `"periodic-box"` | `periodic-box` (stage 0), or a generated mesh of the unit cube: `hex`, `tet`, `poly`, `mixed` (inner faces only, no boundary) |
| `mesh.cells` | 10 000 000 | approximate cell count (of the near-cubic box, or of the generated mesh) |
| `mesh.nx`, `mesh.ny`, `mesh.nz` | – | exact box dimensions (all three, at least 3 each) |
| `mesh.ordering` | `"natural"` | `"natural"` or `"shuffled"` cell numbering |
| `mesh.seed` | 2026 | seed of the fields and of the shuffle |
| `run.backends` | the three CPU backends | see `foehn bench backends` |
| `run.kernels` | all four | `face-flux`, `cell-gather`, `spmv-csr`, `spmv-sell` |
| `run.repeat` | 10 | timed runs per kernel (median and minimum reported) |
| `run.warmup` | 2 | untimed runs first |
| `run.threads` | all cores | CPU threads |
| `run.device` | 0 | CUDA device ordinal |
| `run.rust_ptx` | `target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx` | PTX of the single-source GPU kernels |
| `run.tolerance` | 1e-12 | allowed relative difference to `cpu-serial` |
| `sell.chunk` | 32 | SELL slice height C |
| `sell.sigma` | 1 | SELL sorting window σ |
| `output.directory` | `"output/g1"` | report directory |
| `output.name` | `"g1"` | base name of `NAME.json`, `NAME.md`, `NAME.svg` |

### Mesh case-file keys

Every key is optional; a case needs either `mesh.file` or a
`[mesh.generate]` table. Paths are relative to the case file.

| Key | Default | Meaning |
|---|---|---|
| `mesh.file` | – | Gmsh MSH 4.1 file, ASCII or binary (D33) |
| `mesh.scale` | 1.0 | factor on the coordinates of the file (0.001: millimetres to metres) |
| `mesh.polyhedral` | false | convert a tetrahedral mesh into polyhedra (dual mesh, one cell per node) |
| `mesh.renumber` | `"none"` | cell order: `none` (as read), `rcm` (reverse Cuthill-McKee), `shuffle` |
| `mesh.seed` | 2026 | seed of `shuffle` |
| `[mesh.patches]` | – | patch renaming, `old = "new"` |
| `mesh.generate.type` | `"hex"` | `hex`, `tet`, `poly` or `mixed` |
| `mesh.generate.dimension` | 3 | 2: one cell layer thick in z (D3) |
| `mesh.generate.lower`, `upper` | `[0,0,0]`, `[1,1,1]` | corners of the box |
| `mesh.generate.cells` | `[8, 8, 8]` | intervals of the base hexahedral grid per direction (z ignored in 2D) |
| `mesh.generate.grading.x` (`.y`, `.z`) | `{ ratio = 1 }` | `ratio`: last to first interval; with `symmetric = true`: middle to end intervals |
| `mesh.generate.refine` | 0 | halve every interval this many times (nested meshes for grid-convergence tests) |
| `mesh.generate.perturb` | 0 | random node displacement as a fraction of the local interval, 0 to 0.3; boundary nodes stay on their faces |
| `mesh.generate.wedge` | none | revolve the 2D mesh (`dimension = 2`, `lower[1] ≥ 0`) into an axisymmetric wedge of this angle (degrees, up to 90) about the x axis; `zmin` and `zmax` become the wedge planes, `ymin` disappears on the axis |
| `mesh.generate.seed` | 2026 | seed of the perturbation (counter-based, the same on every machine) |
| `partition.parts` | 1 | number of parts (1: no decomposition) |
| `partition.method` | `"metis"` | `metis` or `rcb` (recursive coordinate bisection) |
| `partition.halo_layers` | 1 | halo depth |
| `check.locate_points` | 1000 | random points of the point-location test (0: none) |
| `check.seed` | 2026 | seed of the points |
| `output.directory` | `output/<case name>` | report directory |
| `output.vtkhdf` | true | write `mesh.vtkhdf` (fields: volume, non-orthogonality, skewness, aspect ratio, volume ratio, closure, rank) |

The generator names the six sides `xmin`, `xmax`, `ymin`, `ymax`, `zmin`,
`zmax`. A Gmsh file gives its surface physical groups as patches (or
`surface_<tag>` for surfaces without one) and its volume physical groups as
cell zones.

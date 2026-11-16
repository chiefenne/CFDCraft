# Poisson studies (`examples/cases/poisson/`)

Part of the [Case catalogue](../cases.md): the pressure solver on the Poisson equation, with hypre as the reference; Poisson keys.

On this page:

- [Poisson case-file keys](#poisson-case-file-keys)

Stage 3: the pressure solver on the Poisson equation −∇·(Γ∇φ) = s
([explanation](../verification/stage3.md)). Run with `foehn poisson FILE`;
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
| `poisson/poisson_poly_box.toml` | generated poly 3D, 729 to 274 625 cells (4 levels) | CG + AMG (mixed and float64) and Jacobi-CG, iteration growth ≤ 1.5, 4 threads against 1 | [gallery 10](../gallery/scalar-pressure.md#poisson-iterations), [recipe 15](../recipes/scalar-pressure.md#poisson-study) |
| `poisson/poisson_hex_box.toml` | generated hex 3D, 512 to 262 144 cells | the same on perturbed hexahedra | [recipe 15](../recipes/scalar-pressure.md#poisson-study) |
| `poisson/poisson_tet_box.toml` | generated tet 3D, 3 072 to 196 608 cells (3 levels) | the same on tetrahedra | [recipe 15](../recipes/scalar-pressure.md#poisson-study) |
| `poisson/poisson_mixed_box.toml` | generated mixed 3D, 1 344 to 95 232 cells (3 levels) | the same on hexahedra, prisms, tetrahedra and polyhedra | [recipe 15](../recipes/scalar-pressure.md#poisson-study) |
| `poisson/poisson_poly_gpu.toml` | generated poly 3D, 912 673 cells | CPU and GPU against hypre on the GPU: times per cell, GPU equals CPU | [recipe 16](../recipes/gpu.md#poisson-gpu), `tools/stage3_colab.ipynb` |

## Poisson case-file keys

`[mesh]`, `[mesh.generate]` as in the [mesh cases](meshes.md#mesh-case-file-keys),
`[problem]` as in the [scalar cases](scalar.md#scalar-case-file-keys) (the velocity
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

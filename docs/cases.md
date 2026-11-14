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
for ParaView. Run with `cfdcraft mesh FILE`; a Gmsh file can be given
directly instead of a case file (`cfdcraft mesh pipe.msh`), command-line
options override the file (`cfdcraft mesh --help`).

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
with `cfdcraft scalar FILE`; command-line options override the file
(`cfdcraft scalar --help`). Each run writes `scalar.md` and `scalar.json`
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
| `scalar/oblique_step_2d.toml` | generated poly 2D, 2 401 cells | upwind, linear upwind and both limiters as `[[variants]]`; bounds check | gallery 9, recipe 13 |
| `scalar/mms_poly_gpu.toml` | generated poly 3D, 912 673 cells | CPU against two GPU runs: GPU equals CPU and GPU run to run | recipe 14, `tools/stage2_colab.ipynb` |

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
| `solver.scheme` | `"linear-upwind"` | `upwind` or `linear-upwind` |
| `solver.limiter` | `"none"` | `none`, `venkatakrishnan`, `venkatakrishnan-wang` |
| `solver.limiter_k` | 5.0 | K: ε² = (K h)³ (Venkatakrishnan) or ε = K (φ_max − φ_min) (Venkatakrishnan–Wang); 0 for strictly bounded |
| `solver.nonorth_correction` | true | non-orthogonal diffusion correction |
| `solver.nonorth_limit` | −1 | limit of that correction as a multiple of the orthogonal part (negative: none) |
| `solver.relaxation` | 1.0 | implicit under-relaxation factor |
| `solver.max_outer` | 200 | outer (deferred-correction) iterations |
| `solver.tolerance` | 1e-10 | outer tolerance on ‖b − Aφ‖/‖b‖ |
| `solver.linear_max_iter` | 200 | BiCGStab iterations per outer iteration |
| `solver.linear_tolerance` | 1e-3 | relative tolerance of each linear solve |
| `run.backends` | `["cpu"]` | `cpu`, `cpu-serial`, `cuda`; later backends are compared with the first |
| `run.threads` | all cores | CPU threads |
| `run.rust_ptx` | `target/nvptx64-nvidia-cuda/release/cfdcraft_ptx_kernels.ptx` | PTX module of the GPU backend |
| `verification.levels` | none | refinement levels of `[mesh.generate]`; two or more: order-of-accuracy sweep |
| `verification.expected_order` | none | `[lo, hi]`: fail if the order on the finest pair is outside |
| `verification.bounds` | none | `[lo, hi]`: fail if φ leaves this range (by more than 1e-6) |
| `[[variants]]` | none | `name`, `solver = { ... }`, optional `bounds`: the case once per variant |
| `output.directory` | `output/<case name>` | output directory |
| `output.vtkhdf` | true | write `fields.vtkhdf` |
| `output.profile` | none | `{ from, to, points }`: sample φ along a line (`profile.csv`, `profile.svg`) |

## Benchmarks (`examples/cases/bench/`)

Gate G1, the kernel benchmark of stage 0 ([G1](benchmarks/g1.md)). Run with
`cfdcraft bench g1 FILE`; every key is optional, command-line options
override the file (`cfdcraft bench g1 --help`).

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
| `run.backends` | the three CPU backends | see `cfdcraft bench backends` |
| `run.kernels` | all four | `face-flux`, `cell-gather`, `spmv-csr`, `spmv-sell` |
| `run.repeat` | 10 | timed runs per kernel (median and minimum reported) |
| `run.warmup` | 2 | untimed runs first |
| `run.threads` | all cores | CPU threads |
| `run.device` | 0 | CUDA device ordinal |
| `run.rust_ptx` | `target/nvptx64-nvidia-cuda/release/cfdcraft_ptx_kernels.ptx` | PTX of the single-source GPU kernels |
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

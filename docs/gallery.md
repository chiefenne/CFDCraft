# Example gallery

What Foehncraft computes, one picture per case, with the commands that
reproduce it. Command sequences for everyday tasks are in
[Recipes](recipes.md); all case files are listed in the
[Case catalogue](cases.md).

All commands run from the repository root. Case files are in
`examples/cases/`.

| # | Example | Stage | Feature |
|---|---|---|---|
| 1 | [Kernel benchmark G1 on the CPU](#g1-cpu) | 0 | gate G1, reproducible kernels |
| 2 | [Kernel benchmark G1 on a Tesla T4](#g1-gpu) | 0 | gate G1, single-source Rust on the GPU |
| 3 | [The four test-mesh types in 3D](#mesh-types-3d) | 1 | built-in mesher, grading, perturbation, quality report |
| 4 | [A polyhedral mesh in four parts](#poly-partition) | 1 | partitioning, halo layout, point location |
| 5 | [The four test-mesh types in 2D](#mesh-types-2d) | 1 | one-cell-thick 2D meshes, patch renaming |
| 6 | [T-junction from Gmsh](#tjunction) | 1 | Gmsh reader, conversion to polyhedra, renumbering |
| 7 | [G1 on polyhedral and tetrahedral meshes](#g1-stage1) | 1 | G1 rerun on stage-1 meshes, natural and shuffled order |
| 8 | [Convection–diffusion, manufactured solution](#scalar-mms) | 2 | least-squares gradients, linear upwind, non-orthogonal diffusion, BiCGStab, order of accuracy |
| 9 | [Oblique step: five convection schemes](#oblique-step) | 2, 5 | Venkatakrishnan and Venkatakrishnan–Wang limiters, bounded central scheme, bounded solutions |
| 10 | [Pressure solver: iterations over the mesh size](#poisson-iterations) | 3 | CG with smoothed-aggregation AMG, mixed precision, hypre reference |
| 11 | [Poisson manufactured solution with the multigrid](#poisson-amg-mms) | 3 | CG + AMG inside the outer iterations, non-orthogonal correction |
| 12 | [Lid-driven cavity against Ghia et al.](#cavity) | 4 | SIMPLEC, Rhie–Chow, reference cell, centre-line comparison, 2D and 3D cells |
| 13 | [Poiseuille flow: order of accuracy of the flow solver](#poiseuille) | 4 | second order on all mesh types, velocity inlet with the exact profile, pressure outlet |
| 14 | [Taylor–Green vortex: second order in time and space](#taylor-green) | 5 | PISO, variable-step BDF2, Crank–Nicolson, extrapolation, outer iterations at large steps |
| 15 | [Cylinder in a channel: the Schäfer & Turek benchmark](#cylinder) | 5 | built-in O-grid, steady Re 20, vortex street at Re 100, adaptive time step, force coefficients |
| 16 | [The flow loop on a GPU](#gpu-flow) | 5 | single-source flow kernels, CPU = GPU bitwise, time per cell on a Tesla T4 |
| 17 | [Adjoint-ready kernels: derivatives and checkpointing](#adjoint-ready) | 5 | dual numbers through the kernel bodies, transposed matrices, binomial checkpointing |
| 18 | [A run that a program can drive](#run-contract) | 5 | `--json`, exit codes, `--set`, limits, signals, `--restart latest`, diagnostic codes, estimate before the run, byte-identical files |
| 19 | [Balances and the divergence report](#health) | 5 | mass and kinetic-energy ledger, numerical dissipation of the schemes (second law), divergence report with ranked causes |
| 20 | [Time step from the error, output at exact times, potential-flow start](#time-control) | 5 | BDF2 error estimate with rejected steps, output in simulated time, potential flow from the pressure multigrid |
| 21 | [Boundary conditions: slip walls, wedges, tables, backflow, convective outlets](#boundary-conditions) | 5 | slip wall, axisymmetric wedge mesh, inlet profile from a table, backflow treatment at pressure outlets, convective outlet |

---

<a id="g1-cpu"></a>
## 1. Kernel benchmark G1 on the CPU

![G1 effective bandwidth of the four kernels on 10 million cells](figures/gallery/g1_cpu_10m.svg)

The four kernels that dominate a finite-volume solver — face flux, cell
gather and the sparse matrix-vector product in CSR and SELL-C-σ form — on a
periodic box of 10 million hexahedra in face-based form. The single-source
kernel bodies (`cpu-threads`) run as fast as hand-written Rust loops
(`cpu-native`) and give bitwise the same results as the single-thread run
(`cpu-serial`), for any number of threads. The bars are effective bandwidth:
bytes that must move over the median time. Procedure, numbers and reading:
[G1](benchmarks/g1.md).

| | |
|---|---|
| Mesh | periodic box 215 × 215 × 216 = 9 984 600 cells, 29 953 800 faces, natural numbering |
| Machine | chart: Intel Xeon @ 2.10 GHz, 4 virtual cores (shared cloud machine); the verdict run on a quiet i9-9900K is in [G1](benchmarks/g1.md#cpu-verdict-imac-i9-9900k-2026-11-14) |
| Case file | `bench/g1_cpu_10m.toml` |

```
cargo build --release -p foehncraft-cli
target/release/foehn bench g1 examples/cases/bench/g1_cpu_10m.toml
# picture: output/g1_cpu_10m/g1.svg, numbers: g1.md and g1.json
```

The same run with a randomly numbered mesh (`bench/g1_cpu_10m_shuffled.toml`)
shows what a badly ordered mesh costs: three to five times less bandwidth.

![G1 with shuffled cell numbering](figures/gallery/g1_cpu_10m_shuffled.svg)

[↑ back to top](#example-gallery)

---

<a id="g1-gpu"></a>
## 2. Kernel benchmark G1 on a Tesla T4

![G1 effective bandwidth on a Tesla T4, 10 million cells](figures/gallery/g1_gpu_t4_10m.svg)

The same four kernels on a Tesla T4 (Google Colab, free tier). Blue: the
hand-written CUDA C kernels; orange: the single-source Rust kernels compiled
to PTX. Both reach 220–255 GB/s of the T4's 320 GB/s, the Rust kernels are
within 1.0 % of the CUDA C kernels (gate limit: 10 %), and all results are
bitwise equal to the CPU. The two short bars per group are the CPU backends
of the Colab machine (2 virtual cores). Reading and limits:
[G1](benchmarks/g1.md#gpu-verdict-tesla-t4-google-colab-2026-11-14).

| | |
|---|---|
| Mesh | periodic box 215 × 215 × 216 = 9 984 600 cells, natural numbering |
| Machine | Tesla T4 (sm_75, 16 GB), Intel Xeon @ 2.00 GHz, 2 threads |
| Case file | `bench/g1_gpu_10m.toml` |

```
bash tools/g1_cloud.sh      # on any Linux machine with an NVIDIA GPU of compute capability 7.0+
# or open tools/g1_colab.ipynb in Google Colab (recipe 3)
```

[↑ back to top](#example-gallery)

---

<a id="mesh-types-3d"></a>
## 3. The four test-mesh types in 3D

![Hexahedral, tetrahedral, polyhedral and mixed test meshes of a box, cut open](figures/gallery/mesh_types_3d.png)

The built-in test mesher (D34) makes four mesh types of the same box
2 × 1 × 1: hexahedra, tetrahedra (six per hexahedron), polyhedra (the dual
of the tetrahedral mesh, one cell per node) and a mixed mesh of hexahedra,
prisms and tetrahedra, with polyhedra where a hexahedron meets the
triangles of its neighbours. All four are clustered towards the walls
y = 0 and y = 1 (middle-to-wall ratio 4) and have perturbed nodes (15 % of
the local interval), so the geometry is never trivially orthogonal. Every
run writes the quality report (`mesh_report.md` and `.json`): closure of
every cell, total volume against the enclosed volume, non-orthogonality,
skewness, face warp, aspect and volume ratio, inverted faces and negative
volumes. All four meshes have the exact volume 2 and closure below 1e-15.

| | |
|---|---|
| Meshes | hex 576, tet 3 456, poly 819, mixed 1 728 cells |
| Case files | `mesh/hex_box.toml`, `mesh/tet_box.toml`, `mesh/poly_box.toml`, `mesh/mixed_box.toml` |

```
cargo build --release -p foehncraft-cli
for t in hex tet poly mixed; do
  target/release/foehn mesh examples/cases/mesh/${t}_box.toml
done
# output/<case>/mesh_report.md, mesh_report.json, mesh.vtkhdf (open in ParaView)
python3 tools/render_mesh.py output/poly_box/mesh.vtkhdf poly_box.png --cut
```

[↑ back to top](#example-gallery)

---

<a id="poly-partition"></a>
## 4. A polyhedral mesh in four parts

![Polyhedral unit cube in four METIS parts](figures/gallery/poly_partition.png)

A polyhedral mesh of the unit cube (2 197 cells, nodes perturbed by 20 %)
split into four parts with METIS (D32), one colour per rank. Each part is
laid out as solver ranks will see it: inner cells, then cells touching the
halo, then the halo cells by source rank, each group in reverse
Cuthill-McKee order; faces inner, processor, then boundary by patch. The
report lists the cells, halo cells, processor faces and neighbours of each
rank (imbalance 1.001, 1 230 cut faces). The point-location test of D47
throws 5 000 random points and finds each one by walking from a random cell
through the faces; all 5 000 are found in the right cell.

| | |
|---|---|
| Mesh | poly, 12 × 12 × 12 base grid, 2 197 cells |
| Case file | `mesh/poly_partition.toml` |

```
target/release/foehn mesh examples/cases/mesh/poly_partition.toml
target/release/foehn mesh examples/cases/mesh/poly_partition.toml --method rcb --parts 6
python3 tools/render_mesh.py output/poly_partition/mesh.vtkhdf part.png --cut --field rank --categorical
```

[↑ back to top](#example-gallery)

---

<a id="mesh-types-2d"></a>
## 5. The four test-mesh types in 2D

![Hexahedral, tetrahedral, polyhedral and mixed 2D meshes](figures/gallery/mesh_types_2d.png)

2D meshes are one cell thick in z (D3): quadrilaterals become hexahedra,
triangles become prisms, polygons become prismatic polyhedra. The mixed
mesh combines both. Clustering towards x = 0 and towards the walls, nodes
perturbed in x and y. The case files rename the generator's patches
(`xmin` to `inlet` and so on) with a `[mesh.patches]` table; `zmin` and
`zmax` take the symmetry condition in a 2D case.

| | |
|---|---|
| Meshes | hex 160, tet 320, poly 187, mixed 240 cells |
| Case files | `mesh/hex_2d.toml`, `mesh/tet_2d.toml`, `mesh/poly_2d.toml`, `mesh/mixed_2d.toml` |

```
for t in hex tet poly mixed; do
  target/release/foehn mesh examples/cases/mesh/${t}_2d.toml
done
python3 tools/render_mesh.py output/poly_2d/mesh.vtkhdf poly_2d.png --top
```

[↑ back to top](#example-gallery)

---

<a id="tjunction"></a>
## 6. T-junction from Gmsh

![T-junction, tetrahedral mesh from Gmsh, cut at the mid plane, four ranks](figures/gallery/tjunction_tet.png)

![The same T-junction converted into polyhedra](figures/gallery/tjunction_poly.png)

The T-junction mixer of stage 7 (D45) meshed with Gmsh
(`examples/meshes/tjunction.geo`, main pipe 100 mm diameter, branch 50 mm)
and read by the Gmsh v4 reader (D33). Physical groups become patches
(`inlet_main`, `inlet_branch`, `outlet`, `wall`) and the cell zone `fluid`.
The second case converts the 77 776 tetrahedra into 15 927 polyhedra, one
per node; on the curved wall the boundary faces follow the triangles of
the tetrahedral mesh, so both meshes have the same volume to the last
digit (8.310580753620e-3 m³). Both are renumbered by reverse
Cuthill-McKee and split into four parts. The point-location walk stops at
the curved, non-convex wall for some points (176 of 1 195); the bin search
finds those, none lands in a wrong cell.

| | |
|---|---|
| Meshes | tet 77 776 cells, poly 15 927 cells (7 to 23 faces, mean 15) |
| Case files | `mesh/tjunction_tet.toml`, `mesh/tjunction_poly.toml` |

```
gmsh -3 examples/meshes/tjunction.geo -format msh41 -o examples/meshes/tjunction.msh
target/release/foehn mesh examples/cases/mesh/tjunction_tet.toml
target/release/foehn mesh examples/cases/mesh/tjunction_poly.toml
python3 tools/render_mesh.py output/tjunction_poly/mesh.vtkhdf tj.png \
    --field rank --categorical --half --top --crop 0.2 0.8 -0.06 0.31
```

[↑ back to top](#example-gallery)

---

<a id="g1-stage1"></a>
## 7. G1 on polyhedral and tetrahedral meshes

![G1 on 2 million polyhedra, natural order](figures/gallery/g1_cpu_poly_2m.svg)

The G1 kernels rerun on the stage-1 meshes, as D8 asks: 2 million
polyhedra (14 faces per interior cell, uneven rows) in the generator's
order and shuffled, and 2 million tetrahedra. All results are bitwise
equal between the backends. On this shared 4-core cloud machine the
single-source bodies beat the hand-written loops on cell gather (by 25 to
40 %), tie on SpMV-CSR and lose 4 to 30 % on face flux, where one kernel
call takes only 2 to 10 ms; the same gap appears on a 2 million periodic
box, so it is a small-kernel effect, not a polyhedral one. The verdict
needs a quiet machine (recipe 10); the numbers are in
[G1](benchmarks/g1.md#stage-1-rerun-polyhedral-and-tetrahedral-meshes-2026-11-14).

![G1 on 2 million polyhedra, shuffled order](figures/gallery/g1_cpu_poly_2m_shuffled.svg)

A shuffled polyhedral mesh costs a factor 3 to 4 in bandwidth: the reason
the mesh is renumbered by reverse Cuthill-McKee on input (example 6).

![G1 on 2 million tetrahedra](figures/gallery/g1_cpu_tet_2m.svg)

| | |
|---|---|
| Meshes | poly 2 048 383 cells, 14 145 894 faces; tet 1 999 620 cells |
| Machine | Intel Xeon @ 2.30 GHz, 4 virtual cores (shared cloud machine) |
| Case files | `bench/g1_cpu_poly_2m.toml`, `bench/g1_cpu_poly_2m_shuffled.toml`, `bench/g1_cpu_tet_2m.toml`, `bench/g1_gpu_poly_4m.toml` |

```
target/release/foehn bench g1 examples/cases/bench/g1_cpu_poly_2m.toml
target/release/foehn bench g1 examples/cases/bench/g1_cpu_poly_2m_shuffled.toml
target/release/foehn bench g1 examples/cases/bench/g1_cpu_tet_2m.toml
```

[↑ back to top](#example-gallery)

---

<a id="scalar-mms"></a>
## 8. Convection–diffusion, manufactured solution

![Error of φ on 35 937 polyhedra, one octant cut away](figures/gallery/mms_poly_error.png)

Steady convection–diffusion of a scalar with a manufactured solution
(Roache 2002) on perturbed polyhedral meshes: φ = 1 + sin(πx + 0.5)
cos(πy + 0.3) cos(πz/2 + 0.2), u = (1, 0.6, 0.3) m/s, Γ = 0.2 kg/(m s),
fixed values of the exact solution on all six patches. The operators of
stage 2 at work: weighted least-squares gradients, linear upwind as a
deferred correction, diffusion with the over-relaxed non-orthogonal
correction, gather-based assembly into a sliced-ELL matrix, BiCGStab with
Jacobi. The picture is the error φ − φ_exact on the finest level; it is
largest in the interior, where the solution has its largest curvature.

![Convergence: L2 error over mesh size, observed order 1.83 and 1.98](figures/gallery/mms_poly_convergence.svg)

The same case on hexahedra, tetrahedra, mixed cells and 2D polygons gives
observed orders between 1.87 and 2.08 on the finest pair; table and
discussion in [Stage 2 verification](verification/stage2.md#order-of-accuracy).
On an NVIDIA GPU the operators give bitwise the CPU result: the 912 673-cell
case `mms_poly_gpu` on a Tesla T4 runs at 311 ns per cell and iteration
([GPU against CPU](verification/stage2.md#gpu-against-cpu-tesla-t4-2026-11-14),
recipe 14).

| | |
|---|---|
| Meshes | generated poly 3D, base 4³ intervals, refine 1–3: 729, 4 913, 35 937 cells, perturb 0.15 |
| Machine | Intel Xeon @ 2.10 GHz, 4 virtual cores; 2.8 s for the three levels |
| Case files | `scalar/mms_poly_box.toml`, `mms_hex_box.toml`, `mms_tet_box.toml`, `mms_mixed_box.toml`, `mms_poly_2d.toml` |

```
cargo build --release -p foehncraft-cli
target/release/foehn scalar examples/cases/scalar/mms_poly_box.toml
# report: output/mms_poly_box/scalar.md, plot: convergence.svg, fields: fields.vtkhdf
python3 tools/render_mesh.py output/mms_poly_box/fields.vtkhdf mms_poly_error.png --field error --cut
```

[↑ back to top](#example-gallery)

---

<a id="oblique-step"></a>
## 9. Oblique step: five convection schemes

| first-order upwind | linear upwind, no limiter | linear upwind, Venkatakrishnan (K = 0) | bounded central (β = 0.5) |
|---|---|---|---|
| ![upwind](figures/gallery/oblique_step_upwind.png) | ![linear upwind](figures/gallery/oblique_step_linear-upwind.png) | ![Venkatakrishnan](figures/gallery/oblique_step_venkatakrishnan.png) | ![bounded central](figures/gallery/oblique_step_bounded-central.png) |

![φ along x = 0.8 for the five schemes](figures/gallery/oblique_step_profile.svg)

Pure convection of a step at 30° to the x axis across 2 401 polygonal
prisms. Upwind is bounded but smears the step over many cells. Linear
upwind without a limiter is sharp but over- and undershoots by 7 % and
5 %. With the Venkatakrishnan limiter (and with Venkatakrishnan–Wang) at
K = 0 the step is as sharp as with linear upwind and φ stays in [0, 1].
The bounded central scheme of stage 5 (D25) takes the central face value
where the normalised variable of the upwind cell shows a smooth profile,
blends towards upwind near extrema and is upwind at them; on these
distorted polygons it needs the limiter of the upwind cell as well. It
stays in [0, 1] with a 10–90 % width of the step of 0.090 m (linear
upwind with the limiter 0.082 m, upwind 0.172 m). The run checks the
bounds of every limited variant and fails if one leaves [0, 1].

| | |
|---|---|
| Mesh | generated poly 2D, 48 × 48 intervals, perturb 0.1: 2 401 cells |
| Case file | `scalar/oblique_step_2d.toml` (five `[[variants]]`) |

```
target/release/foehn scalar examples/cases/scalar/oblique_step_2d.toml
# report: output/oblique_step_2d/scalar.md, profiles: profile.svg,
# fields per variant: output/oblique_step_2d/<variant>/fields.vtkhdf
python3 tools/render_mesh.py output/oblique_step_2d/venkatakrishnan/fields.vtkhdf step.png --field phi --top --range -0.08 1.08
```

[↑ back to top](#example-gallery)

---

<a id="poisson-iterations"></a>
## 10. Pressure solver: iterations over the mesh size

![CG iterations to 1e-10 over the number of cells on polyhedral meshes](figures/gallery/poisson_poly_iterations.svg)

The Poisson equation −∇·(Γ∇φ) = s, the model of the pressure equation, on
perturbed polyhedral meshes from 729 to 274 625 cells, each solved from
φ = 0 to ‖b − Aφ‖ ≤ 1e-10 ‖b‖. CG with the smoothed-aggregation multigrid
of stage 3 needs 10 to 11 iterations on every level, in mixed precision
(float32 hierarchy) as in float64; hypre's BoomerAMG-PCG needs 8 to 14;
plain aggregation grows from 10 to 34 and Jacobi-CG from 33 to 223, the
growth of a method without a coarse-grid correction. All solutions agree
with hypre's to 1.3e-9 (relative), and each in-house solver gives bitwise
the same solution on 4 threads and on 1. The same holds on hexahedra,
tetrahedra and mixed meshes (10 to 13 iterations; table in
[Stage 3 verification](verification/stage3.md#iterations-over-the-mesh-size)).

| | |
|---|---|
| Meshes | generated poly 3D, base 4³ intervals, refine 1–4: 729, 4 913, 35 937, 274 625 cells, perturb 0.15 |
| Machine | Intel Xeon @ 2.30 GHz, 4 virtual cores; about 25 s with hypre and both backends |
| Case files | `poisson/poisson_poly_box.toml`, `poisson_hex_box.toml`, `poisson_tet_box.toml`, `poisson_mixed_box.toml` |

```
cargo build --release -p foehncraft-cli
target/release/foehn poisson examples/cases/poisson/poisson_poly_box.toml
# report: output/poisson_poly_box/poisson.md, plot: iterations.svg
# with hypre as the reference (recipe 17):
tools/build_hypre.sh target/hypre
HYPRE_DIR=$PWD/target/hypre cargo build --release -p foehncraft-cli --features hypre
target/release/foehn poisson examples/cases/poisson/poisson_poly_box.toml \
    --solvers amg,amg-double,amg-plain,jacobi,hypre
```

[↑ back to top](#example-gallery)

---

<a id="poisson-amg-mms"></a>
## 11. Poisson manufactured solution with the multigrid

The pressure solver inside the outer iterations of `foehn scalar`: the
Poisson equation with the manufactured solution of gallery 8, Γ = 1, on
polyhedra with the non-orthogonal correction as a deferred correction,
solved with `linear_solver = "cg"` and `preconditioner = "amg"` and, as a
second variant, with BiCGStab and Jacobi. Both reach the same solution and
the same observed order (1.75, 1.90); on 35 937 polyhedra the multigrid
needs 110 linear iterations in 23 outer iterations against 944 for
BiCGStab with Jacobi, and the time per outer iteration drops from 2 357 to
933 ns per cell.

| Variant | Cells | Outer it. | Linear it. | L2 error | ns/cell/it. |
|---|---:|---:|---:|---:|---:|
| cg-amg | 729 | 17 | 64 | 3.772e-3 | 4 511 |
| bicgstab-jacobi | 729 | 17 | 149 | 3.772e-3 | 5 807 |
| cg-amg | 4 913 | 21 | 99 | 1.239e-3 | 2 238 |
| bicgstab-jacobi | 4 913 | 21 | 409 | 1.239e-3 | 3 250 |
| cg-amg | 35 937 | 23 | 110 | 3.517e-4 | 933 |
| bicgstab-jacobi | 35 937 | 23 | 944 | 3.517e-4 | 2 357 |

| | |
|---|---|
| Meshes | generated poly 3D, base 4³ intervals, refine 1–3, perturb 0.15 |
| Case file | `scalar/poisson_amg_poly_box.toml` (two `[[variants]]`) |

```
target/release/foehn scalar examples/cases/scalar/poisson_amg_poly_box.toml
# report: output/poisson_amg_poly_box/scalar.md (linear iterations, time per phase incl. multigrid setup)
```

[↑ back to top](#example-gallery)

---

<a id="cavity"></a>
## 12. Lid-driven cavity against Ghia et al.

![Speed in the lid-driven cavity at Re = 1000 on 128 x 128 hexahedra](figures/gallery/cavity_re1000_speed.png)

The first flow of Foehncraft: the unit cavity with the lid sliding at 1 m/s,
steady and laminar, solved with SIMPLEC and the relaxation-independent
Rhie–Chow interpolation of stage 4. No boundary fixes the pressure, so a
reference cell holds its level. At Re = 1000 the primary vortex moves
towards the centre and the two corner vortices at the bottom grow; the
velocity on the two centre lines matches the tabulated multigrid solution
of Ghia, Ghia & Shin (1982) to 0.4 % (u) and 1.4 % (v) of the lid speed.

![u on the vertical centre line, Re = 1000](figures/gallery/cavity_re1000_u.svg)
![v on the horizontal centre line, Re = 1000](figures/gallery/cavity_re1000_v.svg)

The same flow on fully three-dimensional cells, a slab with symmetry planes
front and back: 0.7 % / 1.3 % at Re = 100 on 5 445 polyhedra, 2.0 % / 2.4 %
at Re = 1000 on 28 227 polyhedra, 2.5 % / 2.2 % at Re = 100 on 24 576
tetrahedra. Re = 100 on 64 × 64 hexahedra: 0.44 % / 0.95 %
([u](figures/gallery/cavity_re100_u.svg), [v](figures/gallery/cavity_re100_v.svg)).

![Speed at Re = 100 on polyhedra](figures/gallery/cavity_re100_poly3d_speed.png)

| | |
|---|---|
| Mesh | generated hex 2D, 128 × 128 cells, one layer with symmetry planes |
| Physics | ρ = 1 kg/m³, μ = 0.001 Pa s, lid 1 m/s: Re = 1 000 |
| Solver | SIMPLEC, linear upwind, α_u = 0.9, CG + AMG for the pressure correction; 1 397 outer iterations to 1e-8, 54 s on a shared 4-core cloud machine |
| Case files | `flow/cavity_re1000.toml`, `cavity_re100.toml`, `cavity_re100_poly3d.toml`, `cavity_re100_tet3d.toml`, `cavity_re1000_poly3d.toml` |

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cavity_re1000.toml
# report: output/cavity_re1000/flow.md, charts: centreline_u.svg, centreline_v.svg, residuals.svg
python3 tools/render_mesh.py output/cavity_re1000/fields.vtkhdf cavity.png --field speed --top --no-edges
```

[↑ back to top](#example-gallery)

---

<a id="poiseuille"></a>
## 13. Poiseuille flow: order of accuracy of the flow solver

![L2 velocity error over the mesh size on perturbed polyhedra](figures/gallery/poiseuille_poly_convergence.svg)

Plane channel flow between walls at rest, the exact parabola at the inlet
(`profile = "analytic"`), a static-pressure outlet, Re = 10, on three nested
meshes per cell type. The velocity error falls with the square of the mesh
size on every type, the check of stage 4 (1.8–2.2 on the finest pair);
Couette flow (`couette_*.toml`) comes out at the iteration error, as a
linear profile must.

| Mesh | Cells (finest) | L2 error (finest) | Observed order |
|---|---:|---:|---:|
| hex, perturbed 15 % | 8 192 | 1.10e-4 | 1.98, 1.89 |
| tet (prisms) | 16 384 | 6.70e-5 | 1.93, 1.80 |
| poly, perturbed 15 % | 8 385 | 7.74e-5 | 2.01, 2.00 |
| mixed | 12 288 | 9.05e-5 | 1.94, 1.88 |

| | |
|---|---|
| Meshes | generated 2D, 32 × 16 to 128 × 64 intervals in a 2 × 1 channel |
| Case files | `flow/poiseuille_hex.toml`, `poiseuille_tet.toml`, `poiseuille_poly.toml`, `poiseuille_mixed.toml` |

```
target/release/foehn flow examples/cases/flow/poiseuille_poly.toml
# report: output/poiseuille_poly/flow.md, plot: convergence.svg
```

[↑ back to top](#example-gallery)

---

<a id="taylor-green"></a>
## 14. Taylor–Green vortex: second order in time and space

![Change of the solution over the time step, BDF2](figures/gallery/taylor_green_bdf2_time.svg)

The decaying vortex on the unit square, Re = 100, the exact velocity on
the four sides at each time level. Five runs on the same 64 × 64 mesh
with time steps from 0.01 to 0.000625 s: the change of the solution when
the step is halved falls with Δt², for BDF2 and for Crank–Nicolson, with
one outer iteration per step (plain PISO with two correctors). That needs
the convecting flux and the deferred corrections extrapolated to the new
time; without it the lag of these terms makes PISO first order.

| Scheme | Observed time order (pairs of steps) |
|---|---|
| BDF2 | 2.35, 2.42, 2.21 |
| Crank–Nicolson | 2.40, 2.41, 2.00 |
| BDF2, Courant 1.3 to 9, outer iterations to 1e-9 | 1.97, 2.00 |

![Velocity error over the mesh size](figures/gallery/taylor_green_space.svg)

In space, at a fixed small step, the error falls from 2.6e-2 on 16² to
1.3e-4 on 128² cells, order 2.09 on the finest pair.

| | |
|---|---|
| Mesh | generated hex 2D, 64 × 64 (time), 16² to 128² (space) |
| Physics | ρ = 1 kg/m³, μ = 0.01 Pa s, U = 1 m/s, L = 1 m: Re = 100, t = 0 to 0.4 s |
| Solver | PISO, two correctors, BDF2 or Crank–Nicolson, linear solves to 1e-8 |
| Case files | `flow/taylor_green_bdf2.toml`, `taylor_green_cn.toml`, `taylor_green_space.toml`, `taylor_green_cfl10.toml` |

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/taylor_green_bdf2.toml
# report: output/taylor_green_bdf2/flow.md, plot: time_convergence.svg
target/release/foehn flow examples/cases/flow/taylor_green_space.toml
# report: output/taylor_green_space/flow.md, plot: convergence.svg
```

[↑ back to top](#example-gallery)

---

<a id="cylinder"></a>
## 15. Cylinder in a channel: the Schäfer & Turek benchmark

![Speed around the cylinder at Re = 20, 98 304 cells](figures/gallery/cylinder_re20_speed.png)

The laminar cylinder benchmark of Schäfer & Turek (1996): a cylinder of
0.1 m in a 2.2 × 0.41 m channel, slightly off the centre line, with a
parabolic inflow. At Re = 20 the flow is steady with two attached
vortices; on 98 304 cells the drag and lift coefficients and the pressure
difference across the cylinder lie inside the reference ranges
(C_D 5.5800, C_L 0.01060, Δp 0.11751 Pa).

![Vortex street at Re = 100, 24 576 cells, t = 8 s](figures/gallery/cylinder_re100_speed.png)

At Re = 100 a vortex street forms. The time step follows a Courant number
of 0.8; from t = 5 s the run evaluates the largest drag and lift, the
Strouhal number of the lift and the pressure difference half a period
after the lift maximum.

![Drag and lift coefficients over time, Re = 100](figures/gallery/cylinder_re100_coefficients.svg)

| Case | Cells | C_D (max) | C_L (max) | St | Δp (Pa) |
|---|---:|---:|---:|---:|---:|
| Re 20, reference | | 5.57–5.59 | 0.0104–0.0110 | | 0.1172–0.1176 |
| Re 20 | 24 576 | 5.5808 | 0.01036 | | 0.11750 |
| Re 20 | 98 304 | 5.5800 | 0.01060 | | 0.11751 |
| Re 100, reference | | 3.22–3.24 | 0.99–1.01 | 0.295–0.305 | 2.46–2.50 |
| Re 100 | 6 144 | 3.2319 | 0.9712 | 0.2955 | 2.4819 |
| Re 100 | 24 576 | 3.2305 | 0.9847 | 0.3008 | 2.4854 |

![The O-grid around the cylinder, refine 0](figures/gallery/cylinder_mesh.png)

| | |
|---|---|
| Mesh | built in (`[mesh.cylinder]`): O-grid around the cylinder, graded wake block; refine 0, 1, 2: 6 144, 24 576, 98 304 hexahedra in one layer |
| Physics | ρ = 1 kg/m³, μ = 0.001 Pa s, Ū = 0.2 m/s (Re 20) or 1 m/s (Re 100), D = 0.1 m |
| Solver | Re 20 steady SIMPLEC; Re 100 PISO with BDF2, adaptive Δt (Courant 0.8), t = 0 to 8 s |
| Case files | `flow/cylinder_re20.toml`, `cylinder_re100.toml` |

```
target/release/foehn flow examples/cases/flow/cylinder_re20.toml
target/release/foehn flow examples/cases/flow/cylinder_re100.toml
# report: output/cylinder_re100/flow.md, chart: coefficients.svg, histories: forces.csv, probes.csv
python3 tools/render_mesh.py output/cylinder_re20/fields.vtkhdf cylinder.png --field speed --top --no-edges --size 1800x420 --zoom 2.9
```

[↑ back to top](#example-gallery)

---

<a id="gpu-flow"></a>
## 16. The flow loop on a GPU

![Time per cell and iteration (or step) of eight flow cases, CPU against a Tesla T4](figures/gallery/gpu_flow_times.svg)

The same flow cases run on the CPU and on an NVIDIA GPU
(`run.backend = "cuda"`). The operations of the flow loop are written once
in Rust and run as CPU loops and as PTX; fields, matrices and old time
levels stay in GPU memory, and the pressure and momentum solves run on the
GPU backend. On a free Colab Tesla T4 (sm_75) against the two CPU threads
of the VM, **all eight cases give bitwise the same fields, residual
histories, forces, reports and restart files** (D31: gather-based
assembly, no fused multiply-adds, no atomics, fixed reduction tree).

| Case | Cells | CPU ns/cell | GPU ns/cell | CPU s | GPU s |
|---|---:|---:|---:|---:|---:|
| couette_poly | 153 | 12 358 | 20 594 | 0.2 | 1.5 |
| poiseuille_mixed | 12 288 | 2 925 | 855 | 43.1 | 15.3 |
| cavity_re100 | 4 096 | 4 011 | 1 389 | 7.3 | 2.9 |
| cavity_re100_poly3d | 5 445 | 5 848 | 1 559 | 6.4 | 2.3 |
| cavity_re100_transient | 4 096 | 14 085 | 6 140 | 28.9 | 13.0 |
| taylor_green_cn | 4 096 | 8 678 | 2 621 | 46.4 | 15.1 |
| cylinder_re100_gpu | 98 304 | 9 692 | 1 049 | 191.8 | 22.2 |
| cavity_poly3d_gpu | 970 225 | 5 116 | 722 | 200.1 | 136.4 |

Time per cell and SIMPLEC iteration (steady cases) or PISO step
(transient cases); the wall times include mesh generation and output, which
is why the million-cell case gains less in wall time (1.5 times, 30
iterations) than per iteration. Small meshes do not fill the GPU: the
153-cell case is slower on the GPU. The Krylov solvers still read a dot
product back per iteration (D76 follows later).

| | |
|---|---|
| Mesh | generated cavity (hexahedra, polyhedra, 3D polyhedra up to 970 225 cells), built-in cylinder mesh, refine 2 |
| Solver | SIMPLEC (steady) and PISO with BDF2 or Crank–Nicolson (transient), CG + AMG and BiCGStab on the GPU |
| Case files | `flow/cavity_poly3d_gpu.toml`, `flow/cylinder_re100_gpu.toml`, and the catalogue cases of the table |
| Report | `docs/benchmarks/stage5_gpu_t4.md` |

```
cargo build --release --features cuda   # nightly Rust, PTX: gpu/ptx-kernels/build.sh
python3 tools/stage5_cloud.py           # eight cases on cpu and cuda, bitwise comparison, report.md
target/release/foehn flow examples/cases/flow/cylinder_re100_gpu.toml --backend cuda
python3 tools/gpu_chart.py output/stage5_gpu/report.md gpu_flow_times.svg
```

Without a GPU of your own: `tools/stage5_colab.ipynb` on a free Colab T4
(recipe 24).

[↑ back to top](#example-gallery)

---

<a id="adjoint-ready"></a>
## 17. Adjoint-ready kernels: derivatives and checkpointing

![Primal steps per time step of a reverse sweep against the snapshots in memory](figures/gallery/checkpoint_steps.svg)

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
| Report | [stage 5, part 3](verification/stage5.md#adjoint-ready-kernels-part-3) |

```
cargo test --release -p foehncraft-flow --test adjoint -- --nocapture   # Jacobians, one line per kernel
cargo test --release -p foehncraft-flow --lib -- adjoint_tests checkpoint   # transposed matrices, checkpointing
python3 tools/checkpoint_chart.py checkpoint_steps.svg
```

[↑ back to top](#example-gallery)

---

<a id="run-contract"></a>
## 18. A run that a program can drive

![Memory and file sizes estimated by foehn flow --check against eleven runs](figures/gallery/estimate.svg)

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
| Formats | [interfaces.md](interfaces.md), [keys.md](keys.md) |
| Case | `examples/cases/flow/cavity_re100_budget.toml` |
| Recipe | [26](recipes.md#run-contract) |

```
cargo test --release -p foehncraft-cli --test run_contract
python3 tools/estimate_chart.py estimate.svg
```

[↑ back to top](#example-gallery)

<a id="health"></a>
## 19. Balances and the divergence report

![Numerical dissipation of kinetic energy over the mesh size, lid-driven cavity and Taylor–Green vortex](figures/gallery/energy_budget.svg)

Every run books mass and kinetic energy: what flows and works through
each patch, what the body force adds, what viscosity dissipates, and what
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
$ foehn flow cavity_re1000.toml --set 'solver.algorithm="simple"' --set solver.relaxation_velocity=0.99 --set solver.relaxation_pressure=1.0
      DIVERGED in outer iteration 16: 1 failing cells, velocity first; report in divergence.md
        DIV002 under-relaxation (0.80): relaxation factors velocity 0.99, pressure 1 (SIMPLE); try --set solver.relaxation_velocity=0.69 --set solver.relaxation_pressure=0.31
        DIV006 linear solve (0.60): momentum BiCGStab at its limit of 50 iterations; try --set solver.momentum_max_iter=100
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
| Formats | [interfaces.md](interfaces.md#conservation-ledger-d88-d89-d90-part-4) |
| Case | `examples/cases/flow/cavity_re100_energy.toml` |
| Recipe | [27](recipes.md#health) |

```
cargo test --release -p foehncraft-flow --test health
python3 tools/energy_budget_chart.py energy_budget.svg
```

[↑ back to top](#example-gallery)

---

<a id="time-control"></a>
## 20. Time step from the error, output at exact times, potential-flow start

![Time step over time for three error tolerances and for the Courant number alone, Taylor–Green vortex](figures/gallery/time_step_error.svg)

The Taylor–Green vortex with the step chosen by the error of BDF2 (D65):
from the fourth step on, the local error is estimated from the difference
between the step and the quadratic extrapolation of the last three levels,
and the next step is the one that meets the tolerance, at most 1.2 times
the last and below Courant number 5. With 1e-4 the error never limits and
the run follows the growth limit like the Courant control alone; with
1e-5 and 1e-6 the estimate cuts the step where plain PISO loses accuracy,
at Courant numbers about 2, and repeats the steps that missed (the
crosses). The velocity error at the end falls from 1.6e-3 to 6.3e-4, near
the spatial error of the mesh (5.5e-4), for twice the steps. A rejected
step is taken back exactly: the repeated step is bitwise the step a run
with the shorter Δt takes, and a restart continues bitwise. At the end
every run lands on t = 0.4 s in two equal steps; the same holds for
field, monitor and restart output in simulated time (D66).

![Potential-flow start around the cylinder, 24 576 cells](figures/gallery/cylinder_potential_speed.png)

The start field of `initial.potential_flow` (D52) for the cylinder of
gallery 15: the irrotational flow that carries the parabolic inflow
around the cylinder to the outlet, with the stagnation points on the axis
and the fastest flow above and below the cylinder. It is one Laplace solve
with the pressure multigrid (the pressure-correction equation with unit
coefficients) and four passes of the non-orthogonal correction; the cell
velocity is rebuilt from the face fluxes, which are continuous to the
solver tolerance. The steady run at Re 20 from it converges in 1527
instead of 1570 outer iterations.

| Check | Result |
|---|---|
| Taylor–Green, 64², BDF2, tolerance 1e-4 / 1e-5 / 1e-6 | 26 / 37 / 50 steps, 0 / 6 / 9 rejected, velocity error L2 1.53e-3 / 8.36e-4 / 6.25e-4 (Courant control alone: 25 steps, 1.64e-3) |
| Rejected step, BDF2 and Crank–Nicolson | the repeated step and the next one bitwise equal to the run that took the short step at once |
| Restart of an error-controlled run with rejections | bit-identical continuation |
| Output at exact times | monitor rows at multiples of 0.05 s to 1e-12 s, restart at its time, no step below half the one before |
| Potential start in a channel with uniform inflow | uniform flow to 1e-11 on hexahedra, RMS deviation 1.2 % on polygons perturbed by 15 % (12 % without the non-orthogonal correction) |
| The catalogue against part 4b | 22 flow cases and the scalar cases bitwise equal without the new keys |

| | |
|---|---|
| Code | `foehncraft_flow::transient` (`ErrorControl`, `Transient::reject`), `unsteady` (rejection, `Schedule`), `simplec::Flow::potential_flow` |
| Tests | `adjoint_tests.rs` (`rejected_step_repeats_bitwise`, `potential_flow_start_in_a_channel`), `tests/stage5.rs` (`error_controlled_steps_restart_bit_identically`, `output_lands_on_scheduled_times`) |
| Case files | `flow/taylor_green_error.toml`, `flow/cylinder_potential.toml` |
| Recipe | [28](recipes.md#time-control) |

```
target/release/foehn flow examples/cases/flow/taylor_green_error.toml
python3 tools/time_step_chart.py time_step_error.svg
target/release/foehn flow examples/cases/flow/cylinder_potential.toml
python3 tools/render_mesh.py output/cylinder_potential/fields.vtkhdf potential.png --field speed --top --no-edges --size 1800x420 --zoom 2.9
```

[↑ back to top](#example-gallery)

---

<a id="boundary-conditions"></a>
## 21. Boundary conditions: slip walls, wedges, tables, backflow, convective outlets

![Hagen–Poiseuille flow on a wedge of 20°, coloured by speed](figures/gallery/pipe_wedge_speed.png)

**Axisymmetric flow on a wedge.** A 2D mesh of the half plane y ≥ 0 is
revolved about the x axis into a wedge one cell thick
(`mesh.generate.wedge`); the nodes on the axis merge, so the cells next to
it are prisms with an edge on the axis and the mesh has no axis patch. The
picture is the pipe at 20° to show the shape; the verification runs at
0.5°, where Hagen–Poiseuille flow converges at second order (2.07 and 2.04
on 256 to 4 096 cells). Every cell type of the box generator revolves,
perturbed or not. The planar faces make a wedge of θ a slice of a
polygonal pipe of 360°/θ sides; its result differs from the round pipe by
about θ²/8 (9e-4 at 5°), so applications take 2–5° and a convergence
study 0.5–1°.

![Lift coefficient over the last second at Re = 100: long channel, short channel with a pressure outlet, short channel with a convective outlet](figures/gallery/cylinder_outlets.svg)

**The outlet close behind a body.** The cylinder of gallery 15 at Re = 100
with the outlet 5.5 D behind the cylinder instead of 20 D: the vortices
cross the outlet, so fluid enters through it in parts of every period.
A static pressure outlet with the default backflow treatment (the
prescribed pressure taken as the total pressure of the entering fluid)
gives the lift of the long channel on the same near-wake cells; a
convective outlet (∂u/∂t + U_c ∂u/∂n = 0, U_c the mean outflow speed)
gives C_L,max 1.6 % lower. Drag, Strouhal number and Δp stay inside the
reference ranges with both.

| Case | Cells | C_D,max | C_L,max | St | Δp (Pa) |
|---|---:|---:|---:|---:|---:|
| reference | | 3.22–3.24 | 0.99–1.01 | 0.295–0.305 | 2.46–2.50 |
| long channel (gallery 15) | 6 144 | 3.2319 | 0.9712 | 0.2955 | 2.4819 |
| short, pressure outlet, `backflow = "total-pressure"` | 4 928 | 3.2336 | 0.9713 | 0.2951 | 2.4825 |
| short, convective outlet | 4 928 | 3.2261 | 0.9555 | 0.2952 | 2.4779 |

| Check | Result |
|---|---|
| Slip walls, uniform flow on hex, tet, poly and mixed cells perturbed by 20 % | velocity error 1.5e-11, pressure 1.0e-11 (exact to the solver tolerance) |
| Inlet profile from a table (rows at the face centres) | the run with the exact parabola, velocity error 1.592e-3 for both |
| Reversed flow through a pressure outlet, `backflow = "total-pressure"` | uniform flow and the static pressure p − ½ρU² exact to 3e-11 |
| A side outlet above the channel pressure, `backflow = "suppress"` | the faces close, act as a slip wall, no mass crosses (exactly 0), uniform flow exact to 1e-11 |
| Convective outlet, no pressure boundary | outlet flow = inflow to round-off in every step; the start-up disturbance leaves (error 4e-5 after two flow-through times) |
| Convective outlet, restart and rejected steps | bit-identical continuation, the old boundary levels in the restart file |

| | |
|---|---|
| Code | `foehncraft_mesh::wedge`, `foehncraft_flow::bc` (`Backflow`, `ProfileTable`, the new condition types and checks), `simplec::Convective`, kernels `flow::close_backflow`, `UB_OUTLET`, `UB_CONVECTIVE`, `UB_NO_BACKFLOW`, `PB_OUTLET`, `PB_NO_BACKFLOW` |
| Tests | `tests/bcs.rs`, `tests/stage5.rs` (`convective_outlet_restarts_bit_identically`), `wedge.rs` (`wedges_of_every_type`) |
| Case files | `flow/channel_slip.toml`, `poiseuille_table.toml`, `pipe_wedge.toml`, `cylinder_re100_short.toml`, `cylinder_re100_convective.toml` |
| Recipe | [29](recipes.md#boundary-conditions) |

```
target/release/foehn flow examples/cases/flow/pipe_wedge.toml
target/release/foehn flow examples/cases/flow/pipe_wedge.toml --set mesh.generate.wedge=20.0 --set 'verification.levels=[1]' --set mesh.generate.perturb=0.0 --set output.vtkhdf=true --out output/wedge20
python3 tools/render_mesh.py output/wedge20/fields.vtkhdf pipe_wedge.png --field speed --size 1400x600 --zoom 1.6
target/release/foehn flow examples/cases/flow/cylinder_re100.toml
target/release/foehn flow examples/cases/flow/cylinder_re100_short.toml
target/release/foehn flow examples/cases/flow/cylinder_re100_convective.toml
python3 tools/outlet_chart.py cylinder_outlets.svg
```

[↑ back to top](#example-gallery)

# Stage 2 – operators on a scalar equation

Stage 2 of the [stage plan](../requirements.md#7-stage-plan) builds the
finite-volume operators every later equation uses. They are tested on the
steady convection–diffusion equation of a scalar φ,

∇·(ρ u φ) − ∇·(Γ ∇φ) = s,

on all four mesh types of the test mesher (D34), on the CPU and on an
NVIDIA GPU. This page explains the methods, the checks and the results.
The commands are in [recipes 12–14](../recipes/scalar-pressure.md#scalar-mms), the pictures
in [gallery 8–9](../gallery/scalar-pressure.md#scalar-mms), the case files in the
[case catalogue](../cases/scalar.md#scalar-cases-examplescasesscalar).

## Acceptance checks and results

| Check (stage plan row 2) | Result | Where |
|---|---|---|
| Gradients exact for linear fields on all mesh types | largest relative error 6e-15 on hex, tet, poly and mixed meshes, 2D and 3D, perturbed and unperturbed, with fixed-value and fixed-gradient boundaries | test `gradient_of_a_linear_field_is_exact_on_every_mesh_type` |
| The limiter creates no new extrema | with K = 0 every limited face value lies inside the range of the cell and its neighbours (random field, all mesh types); the converged oblique step stays in [0, 1] to 1e-7 with both limiters | tests `limited_face_values_create_no_new_extrema`, `oblique_step_stays_bounded_with_the_limiters`; case `oblique_step_2d` |
| Convection–diffusion with an exact solution, observed order 1.8–2.2 on hexahedral, tetrahedral, polyhedral and mixed meshes | 1.87 to 2.08 on the finest pair of all eight meshes (table below) | test `order_of_accuracy_on_all_mesh_types` (release builds, CI); cases `mms_*` |
| Bitwise identical results from run to run | bitwise equal for 1, 2, 3 and 4 threads and from run to run, with the limiter on | test `results_are_bitwise_reproducible` |
| GPU and CPU agree to round-off | **bitwise equal** on a Tesla T4 (Colab, 2026-11-14): 912 673 polyhedra with the Venkatakrishnan limiter, both GPU runs equal the CPU run bit for bit; the four oblique-step schemes and the polyhedral sweep are bitwise equal too (see below) | case `mms_poly_gpu`, `tools/notebooks/stage2_colab.ipynb` |

### Order of accuracy

Manufactured solution φ = 1 + sin(πx + 0.5) cos(πy + 0.3) cos(πz/2 + 0.2)
on the unit box (2D: without the z factor), u = (1, 0.6, 0.3) m/s
(2D: (1, 0.6, 0)), ρ = 1 kg/m³, Γ = 0.2 kg/(m s), fixed values of the exact
solution on all boundaries (2D: zero gradient on the two z boundaries). Nodes
perturbed by 15 % of the local interval, a new random perturbation on every
level. Linear upwind without a limiter, non-orthogonal correction on. Error:
volume-weighted L2 norm of φ − φ_exact at the cell centroids; h = (V/N)^(1/3)
(2D: (A/N)^(1/2)).

| Mesh | Cells (levels) | L2 error (finest) | Observed order | Outer it. (finest) | ns/cell/it. (finest) |
|---|---|---:|---|---:|---:|
| hex 3D | 512, 4096, 32768 | 2.809e-4 | 1.70, 1.87 | 10 | 2172 |
| hex 2D | 256, 1024, 4096 | 1.506e-4 | 1.80, 1.94 | 10 | 5704 |
| tet 3D | 3072, 24576, 196608 | 2.228e-4 | 1.92, 1.97 | 23 | 2042 |
| tet 2D | 512, 2048, 8192 | 8.680e-5 | 1.86, 2.08 | 18 | 6295 |
| poly 3D | 729, 4913, 35937 | 5.115e-4 | 1.83, 1.98 | 23 | 2040 |
| poly 2D | 289, 1089, 4225 | 1.826e-4 | 1.79, 1.91 | 16 | 10752 |
| mixed 3D | 1344, 11776, 95232 | 2.640e-4 | 1.68, 1.94 | 21 | 1963 |
| mixed 2D | 384, 1536, 6144 | 1.508e-4 | 1.82, 1.97 | 17 | 7662 |

Machine: Intel Xeon @ 2.10 GHz, 4 virtual cores (shared cloud machine).
Reproduce with `cargo run --release -p foehncraft-fv --example sweep_table`;
one level finer with `cargo test --release -p foehncraft-fv -- --ignored`.

**Why Γ = 0.2.** Convection alone (Γ → 0) and diffusion alone (u = 0) each
converge with order 2.0 on the uniform hexahedral mesh. Combined at
Γ = 0.05 (Péclet number about 24) the observed order jumps between 1.4 and
2.3 up to 32 768 cells: the leading convective and diffusive truncation
errors have opposite signs and partly cancel on the coarse meshes, so the
error is not yet in its asymptotic range. At Γ = 0.2 (Péclet about 6) and
above the sweep is clean, at Γ = 0.01 it is above 2 (2.3 to 2.6). The
scheme is second order in both terms; the choice of Γ only moves the
asymptotic range into meshes small enough for CI.

**Why the coarse pairs are lower.** Each level gets its own random
perturbation, so the meshes are not smooth refinements of each other; the
observed order settles towards 2 from about 4 000 cells on.

### GPU against CPU (Tesla T4, 2026-11-14)

`tools/notebooks/stage2_colab.ipynb` on a free Colab T4 (driver 580.82, CUDA 13.0,
the CPU side is the Colab VM: Intel Xeon @ 2.0 GHz, 2 threads). Case
`mms_poly_gpu`: 912 673 polyhedra, linear upwind with the Venkatakrishnan
limiter (K = 5), 21 outer iterations, final residual 8.0e-11, L2 error
7.241e-5, solved on the CPU and twice on the GPU:

| Backend | Time | ns/cell/it. | vs CPU | gradient / limiter / faces / assembly / linear solver (s) |
|---|---:|---:|---|---|
| CPU, 2 threads | 158.8 s | 8 238 | – | 0.95 / 4.08 / 3.31 / 1.84 / 147.7 |
| Tesla T4, run 1 | 7.6 s | 311 | bitwise equal | 0.23 / 0.18 / 0.07 / 0.21 / 5.26 |
| Tesla T4, run 2 | 7.1 s | 317 | bitwise equal | 0.23 / 0.17 / 0.07 / 0.21 / 5.39 |

The GPU result is bit for bit the CPU result and the same from run to run
(D31). The GPU run is 21 to 22 times faster than the two-thread VM; that
figure only shows that the GPU path works and is not a performance claim,
because the CPU side is a small shared machine and the linear solver
(Jacobi BiCGStab, 97 iterations per outer iteration at the start) dominates
both: 93 % of the time on the CPU, 74 % on the GPU. Per outer iteration
the operator kernels (gradient, limiter, faces, assembly) take 0.03 s on
the GPU for 0.9 million cells. The oblique step (four schemes, 2 401 cells)
is bitwise equal on the GPU for all four variants, and the polyhedral sweep
on the GPU alone gives the observed orders 1.830 and 1.976, the same as the
CPU. Small meshes are launch-bound on the GPU (2 800 ns/cell/it. at 2 401
cells against 311 at 0.9 million).

### Limiters on the oblique step

Pure convection (Γ = 10⁻⁶) of a step at 30° to the x axis on 2 401
polygonal prisms, profiles along x = 0.8 in
[gallery 9](../gallery/scalar-pressure.md#oblique-step):

| Scheme | φ range | Outer iterations |
|---|---|---:|
| first-order upwind | [0, 1] | 4 |
| linear upwind, no limiter | [−0.047, 1.072] | 32 |
| linear upwind, Venkatakrishnan, K = 0 | [0, 1] | 62 |
| linear upwind, Venkatakrishnan–Wang, K = 0 | [0, 1] | 62 |

With K > 0 the Venkatakrishnan limiters allow small over- and undershoots
in exchange for smooth convergence and unlimited gradients in smooth
regions; that is the purpose of ε (Venkatakrishnan 1995). On the same
problem K = 0.3 gives 0.1 to 0.6 % overshoot on hex, tet and poly meshes,
K = 5 about 5 %. Use K = 0 where bounds matter (volume fractions, mass
fractions later), a moderate K where convergence matters.

## Methods

All operators are cell-centred finite volumes on the face-based mesh of
stage 1 (D23). Every kernel is a *body* that computes one output from its
index, compiled once for the CPU and once for the GPU (D8, gate G1).

**Least-squares gradient (D24).** Weighted least squares over the face
neighbours, weights 1/|d|² (Syrakos et al. 2017). The weight vectors
w_k = G⁻¹ ω_k d_k, G = Σ ω_k d_k d_kᵀ are computed once per cell-face entry
on the host, so the kernel is a plain gather ∇φ_c = Σ w_k (φ_k − φ_c), with
no square roots or matrix inversions on the GPU. Boundary faces enter as
stencil points: a fixed-value face at its centroid; a fixed-gradient,
zero-gradient or symmetry face at the projection of the face centroid onto
the face normal through the cell centre, with the value φ_c + g d_n. Both
keep the gradient of a linear field exact, also on the one-cell-thick 2D
meshes, where the two symmetry faces supply the z direction.

**Linear upwind with limiters (D25).** The face value of the upwind cell U
is φ_U + ψ_U ∇φ_U·(x_f − x_U), applied as a deferred correction on top of
implicit first-order upwind. ψ_U is Venkatakrishnan's limiter
(Venkatakrishnan 1995), the minimum over the faces of
(Δ₊² + ε² + 2Δ₋Δ₊)/(Δ₊² + 2Δ₋² + Δ₋Δ₊ + ε²) with Δ₋ the unlimited
increment and Δ₊ the distance to the neighbourhood maximum or minimum,
clipped to 1. ε² = (K h)³ with h³ the cell volume (Venkatakrishnan), or
ε² = (K (φ_max − φ_min))² over the domain (Venkatakrishnan–Wang; Wang 2000),
which makes the limiter independent of the mesh size.

**Diffusion with non-orthogonal correction (D26).** The face area vector is
split as S = Δ + k with Δ = d |S|²/(d·S), the over-relaxed decomposition
(Jasak 1996): the part along d (cell centres, or cell centre to face
centroid on a boundary) is implicit, Γ |S|²/(d·S) (φ_N − φ_P); the rest,
Γ k·∇φ_f with the gradient interpolated to the face, is a deferred
correction. The correction can be limited to a multiple of the orthogonal
part on strongly skewed faces (`nonorth_limit`; off by default).

**Assembly (D31, D37).** One face kernel computes the flux of each face as
F_f = c_P φ_P + c_N φ_N + b_f. One cell kernel gathers the faces of its
cell in stored order into the matrix row, the right-hand side and the
inverse diagonal; no face writes into a cell, so there are no atomics and
the result is independent of the thread count. The matrix is stored once
in sliced ELL (SELL-C-σ, Kreutzer et al. 2014, with C = 32 and σ = 1): rows
in cell order, slot 0 the diagonal, then the interior faces of the cell;
32 rows per slice padded to the longest, so the 32 threads of a GPU warp
load consecutive addresses. Implicit under-relaxation after Patankar (1980)
is part of the row kernel (`relaxation`, default 1).

**BiCGStab (D36).** Van der Vorst's BiCGStab (1992) with right Jacobi
preconditioning, a loose relative tolerance per outer iteration (default
10⁻³). Dot products are summed in blocks of 256 values in sequential order
and the block sums are combined in a fixed binary tree (D31), the same on
the CPU and on the GPU. Outer iterations repeat gradient, limiter, face
coefficients, assembly and linear solve until the normalised residual
‖b − Aφ‖/‖b‖ of the full system falls below `tolerance` (default 10⁻¹⁰).

**GPU.** The same bodies are compiled by rustc to PTX (`gpu/ptx-kernels`,
entry points `fv_*` and `blas_*`) and launched through the CUDA driver; all
fields stay on the GPU, only the block sums of the dot products (one value
per 256 cells) and the limiter range come back. The PTX has no fused
multiply-add on f64, IEEE-rounded division and no atomic instructions,
which CI checks on every pull request (`python3 -m tools.check.kernels`, R10), so
the GPU results are expected to equal the CPU results bitwise.

**Figure of merit (R6).** Every run reports wall-clock seconds per phase
(gradient, limiter, faces, assembly, linear solver), measured after device
synchronisation, and ns per cell and outer iteration. The linear solver
dominates (about 85 % on the CPU) because Jacobi-preconditioned BiCGStab
needs 20–90 iterations per outer iteration; the aggregation AMG of stage 3
is the remedy for the pressure equation, and a loose tolerance is enough for
transported scalars.

## Limits of stage 2

- One partition: the operators run on a mesh without halo cells; the
  multi-rank operators (gradients and limiters of halo cells, global dot
  products) come with stage 7. The kernels already take processor faces as
  interior faces with the neighbour in the halo block.
- Constant ρ, Γ and velocity: the mass flux of a constant velocity is
  conservative to round-off on every mesh. Variable fields come with the
  flow solver (stage 4).
- Steady problems only; time integration is stage 5.

## References

- J. H. Ferziger, M. Perić, R. L. Street, *Computational Methods for Fluid
  Dynamics*, 4th ed., Springer (2020)
- H. Jasak, *Error analysis and estimation for the finite volume method
  with applications to fluid flows*, PhD thesis, Imperial College London
  (1996)
- M. Kreutzer, G. Hager, G. Wellein, H. Fehske, A. R. Bishop, *A unified
  sparse matrix data format for efficient general sparse matrix-vector
  multiplication on modern processors with wide SIMD units*, SIAM J. Sci.
  Comput. 36, C401–C423 (2014)
- B. P. Leonard, *A stable and accurate convective modelling procedure
  based on quadratic upstream interpolation*, Comput. Methods Appl. Mech.
  Eng. 19, 59–98 (1979)
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere
  (1980)
- P. J. Roache, *Code verification by the method of manufactured
  solutions*, J. Fluids Eng. 124, 4–10 (2002)
- A. Syrakos, S. Varchanis, Y. Dimakopoulos, A. Goulas, J. Tsamopoulos,
  *A critical analysis of some popular methods for the discretisation of
  the gradient operator in finite volume methods*, Phys. Fluids 29, 127103
  (2017)
- H. A. van der Vorst, *Bi-CGSTAB: a fast and smoothly converging variant
  of Bi-CG for the solution of nonsymmetric linear systems*, SIAM J. Sci.
  Stat. Comput. 13, 631–644 (1992)
- V. Venkatakrishnan, *Convergence to steady state solutions of the Euler
  equations on unstructured grids with limiters*, J. Comput. Phys. 118,
  120–130 (1995)
- Z. J. Wang, *A fast nested multi-grid viscous flow solver for adaptive
  Cartesian/quad grids*, Int. J. Numer. Methods Fluids 33, 657–680 (2000)

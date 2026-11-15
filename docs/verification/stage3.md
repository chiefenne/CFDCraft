# Stage 3 – the pressure solver

Stage 3 of the [stage plan](../requirements.md#7-stage-plan) builds the
linear solver of the pressure equation, the largest single cost of a
pressure-based flow solver: CG preconditioned by an aggregation-based
algebraic multigrid (D35), with the multigrid hierarchy stored and applied
in float32 while CG stays in float64 (D37), checked against hypre's
BoomerAMG (D20). The model problem is the Poisson equation

−∇·(Γ ∇φ) = s,

discretised exactly as the diffusion part of stage 2 (orthogonal part
implicit, non-orthogonal correction deferred), so the matrix is symmetric
positive definite and has the structure of the pressure-correction
equation of SIMPLE. This page explains the method, the checks and the
results. The commands are in [recipes 15–17](../recipes.md#poisson-study),
the pictures in [gallery 10–11](../gallery.md#poisson-iterations), the case
files in the [case catalogue](../cases.md#poisson-studies-examplescasespoisson).

## Acceptance checks and results

| Check (stage plan row 3) | Result | Where |
|---|---|---|
| Poisson problem on all mesh types: iteration count nearly independent of mesh size | 10 to 13 CG iterations to ‖r‖ ≤ 1e-10 ‖b‖ on hex, tet, poly and mixed meshes from 512 to 274 625 cells; the count grows by a factor of 1.1 to 1.3 while the number of cells grows 64 to 512 times. Jacobi-CG grows 4 to 7 times over the same range (table below) | cases `poisson_*_box` (fail above 1.5), test `amg_iterations_nearly_independent_of_mesh_size` |
| Same converged result as hypre | largest difference to hypre's solution 1.6e-9 of max \|φ\| on all mesh types and levels, both solved to 1e-10 | test `same_converged_result_as_hypre` (feature `hypre`), CI job "Pressure solver against hypre" |
| GPU time per cell compared with hypre and reported | done on a Tesla T4 (Colab): 6.5 against 10.5 ns per cell and iteration at 912 673 cells, bitwise equal to the CPU (section below) | case `poisson_poly_gpu`, recipe 16 |
| Mixed precision does not change the converged result (D37) | float32 and float64 hierarchies need the same number of iterations on every mesh and level; solutions agree to 1e-8 of max \|φ\| (to 1e-11 tolerance); one V-cycle differs at the float32 rounding (about 3e-7) | test `mixed_precision_gives_the_float64_result` |
| Reproducible (D31) | setup and solve bitwise equal on 1 and 4 threads, all mesh types | test `amg_solve_is_bitwise_reproducible`, column "Bitwise vs first backend" of every study |

### Iterations over the mesh size

Problem: Γ = 1, φ = 1 + sin(πx + 0.5) cos(πy + 0.3) cos(πz/2 + 0.2) for the
source and the fixed boundary values on all six patches, nodes perturbed by
15 %. Each level is assembled once and solved from φ = 0 to
‖b − Aφ‖ ≤ 1e-10 ‖b‖. In-house AMG with the defaults (smoothed
aggregation, θ = 0.04, Chebyshev degree 2, coarsest level ≤ 256
unknowns). hypre 3.2.1: PCG with one BoomerAMG V-cycle (HMIS coarsening,
extended+i interpolation with at most 4 entries per row, one level of
aggressive coarsening, θ = 0.5, symmetric hybrid Gauss–Seidel), on one
CPU thread. Times are for one thread on both sides (backend `cpu-serial`);
setup / solve in seconds.

| Mesh | Cells | AMG mixed | AMG float64 | plain aggr. | Jacobi | hypre | max diff. to hypre | AMG setup / solve | hypre setup / solve | ns/cell/it. AMG / hypre | op. complexity |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|
| hex | 512 | 10 | 10 | 11 | 38 | 8 | 1.5e-10 | 0.002 / 0.005 | 0.000 / 0.000 | 1043 / 65 | 1.41 |
| hex | 4 096 | 11 | 11 | 17 | 76 | 11 | 2.6e-10 | 0.008 / 0.015 | 0.003 / 0.005 | 332 / 110 | 1.51 |
| hex | 32 768 | 11 | 11 | 26 | 147 | 13 | 3.7e-10 | 0.072 / 0.051 | 0.034 / 0.038 | 141 / 89 | 1.56 |
| hex | 262 144 | 12 | 12 | 39 | 272 | 14 | 1.1e-9 | 0.390 / 0.355 | 0.199 / 0.330 | 113 / 90 | 1.58 |
| tet | 3 072 | 11 | 11 | 18 | 97 | 13 | 2.7e-10 | 0.005 / 0.006 | 0.002 / 0.005 | 177 / 121 | 1.75 |
| tet | 24 576 | 12 | 12 | 27 | 191 | 16 | 2.8e-10 | 0.045 / 0.024 | 0.018 / 0.024 | 82 / 60 | 1.79 |
| tet | 196 608 | 13 | 13 | 41 | 357 | 18 | 1.6e-9 | 0.278 / 0.237 | 0.167 / 0.268 | 93 / 76 | 1.80 |
| poly | 729 | 10 | 10 | 10 | 33 | 8 | 1.1e-10 | 0.002 / 0.003 | 0.001 / 0.001 | 428 / 98 | 1.23 |
| poly | 4 913 | 10 | 10 | 15 | 63 | 10 | 5.6e-10 | 0.012 / 0.009 | 0.005 / 0.009 | 184 / 180 | 1.23 |
| poly | 35 937 | 11 | 11 | 22 | 114 | 13 | 1.0e-9 | 0.075 / 0.063 | 0.037 / 0.049 | 159 / 105 | 1.22 |
| poly | 274 625 | 11 | 11 | 34 | 223 | 14 | 1.3e-9 | 0.554 / 0.603 | 0.374 / 0.799 | 200 / 208 | 1.21 |
| mixed | 1 344 | 10 | 10 | 15 | 72 | 11 | 2.5e-10 | 0.011 / 0.005 | 0.001 / 0.007 | 359 / 450 | 1.58 |
| mixed | 11 776 | 12 | 12 | 24 | 150 | 14 | 4.5e-10 | 0.025 / 0.018 | 0.010 / 0.016 | 130 / 94 | 1.69 |
| mixed | 95 232 | 13 | 13 | 36 | 292 | 17 | 1.0e-9 | 0.179 / 0.135 | 0.087 / 0.144 | 109 / 89 | 1.73 |

Machine: Intel Xeon @ 2.30 GHz, 4 virtual cores (shared cloud machine;
single runs, so times scatter by ±20 %). Reproduce with recipe 17.

**Reading.** The smoothed-aggregation V-cycle keeps CG at 10 to 13
iterations on every mesh type; hypre's classical AMG needs 8 to 18 and
grows a little faster on these meshes. Plain aggregation (no prolongator
smoothing) is cheaper to set up and has a lower operator complexity, but
with a V-cycle its iterations grow like those of a one-level method (11
to 41), which is why smoothing is the default. Per iteration the in-house
cycle (two Chebyshev steps before and after, float32) costs about as much
as hypre's (one hybrid Gauss–Seidel sweep, float64); the solve times are
close. The host setup is 1.5 to 2 times slower than hypre's; it is not
tuned yet (see limits).

**Mixed against float64.** The iteration counts are identical and the
solutions agree to the tolerance, as D37 expects: the float32 rounding of
the hierarchy perturbs the preconditioner by about 3e-7, far less than
the factor by which one V-cycle reduces the error. With four threads the
float32 hierarchy solves 1.2 to 1.4 times faster than the float64 one
(half the bytes of every matrix and vector inside the cycle); on a GPU the
memory traffic dominates more and the gain is expected to be larger.

**The strength threshold.** Vaněk, Mandel and Brezina use θ = 0.08,
halved on each coarser level. The unstructured cells of the test meshes
have 10 to 20 neighbours, each coupling correspondingly weaker, so with
θ = 0.08 many couplings count as weak and the aggregates stay small
(coarsening ratio 6 on polyhedra). θ = 0.04 gives a coarsening ratio of
about 15, lowers the operator complexity on polyhedra from 1.64 to 1.21
and the iteration growth from 1.56 to 1.10; it is the default. Values
above about 0.07 on polyhedra make every coupling weak and the setup
stops with a message.

### Inside the outer iterations

`inncraft scalar` uses the same solver with `[solver] linear_solver = "cg"`
and `preconditioner = "amg"`. On the Poisson manufactured solution with the
non-orthogonal correction (case `poisson_amg_poly_box`) the outer
deferred-correction iterations, the solution and the observed order
(1.90) are the same as with BiCGStab and Jacobi; the linear iterations
drop from 944 to 110 in total on 35 937 polyhedra (to a relative tolerance
of 1e-4 per outer iteration), and the time per outer iteration from 2 357
to 933 ns per cell. The hierarchy is built once, because the matrix of a
linear problem does not change between outer iterations.

## Methods

**Smoothed aggregation** (Vaněk, Mandel, Brezina 1996), set up on the host
in float64, in a fixed order, so the hierarchy depends only on the matrix:

1. *Strength of connection.* j is strongly coupled to i when
   |a_ij| ≥ θ_l √(|a_ii a_jj|), θ_l = θ · 2⁻ˡ on level l.
2. *Aggregation*, three greedy phases over the rows in stored order (the
   RCM cell order of the mesh): (1) a node whose strong neighbours are all
   free forms an aggregate with them; (2) a free node joins the aggregate
   of its most strongly coupled neighbour aggregated in phase 1; (3) the
   remaining nodes form aggregates with their free strong neighbours. Nodes
   without any strong coupling stay out of every aggregate; the smoother
   alone treats them.
3. *Tentative prolongator* P̃: the constant vector (the near-null space of
   the pressure operator) on each aggregate, normalised.
4. *Prolongator smoothing* P = (I − ω D_F⁻¹ A_F) P̃ with the filtered matrix
   A_F (weak couplings lumped onto the diagonal, so that A_F keeps the row
   sums of A) and ω = (4/3) / λ_max(D_F⁻¹ A_F).
5. *Galerkin coarse operator* A_c = Pᵀ A P, computed row by row
   (Gustavson's algorithm, entries summed in a fixed order). The
   restriction R = Pᵀ is stored explicitly, so that restriction, like
   every other operation of the cycle, is a gather (D31: no atomics).
6. *Coarsest level*: at most `max_coarse` unknowns (256), solved by its
   dense inverse from a Cholesky factorisation. A singular coarsest level
   (pressure level not fixed) stops the setup with a message that points to
   D30.

λ_max(D⁻¹A) comes from a fixed number of power iterations (20) from a
deterministic start vector, with the Rayleigh quotient (x, Ax)/(x, Dx),
capped by the Gershgorin bound.

**V-cycle.** The same polynomial smoother before and after the coarse
correction: Chebyshev of degree 2 in D⁻¹A on the interval
[0.1 λ_max, 1.1 λ_max] (Adams, Brezina, Hu, Tuminaro 2003) by the
three-term recurrence (Saad 2003, algorithm 12.1), or damped Jacobi with
ω = (4/3)/λ_max. Each step is d ← c₁ d + c₂ D⁻¹(b − A x), x ← x + d with
host-computed coefficients (c₁, c₂); the first step from x = 0 needs no
product with A. A polynomial smoother is self-adjoint in the A inner
product, so the cycle is a symmetric positive definite preconditioner
(test `vcycle_is_symmetric_positive_definite`: (M⁻¹x, y) = (x, M⁻¹y) to
1e-12 in float64). Polynomial smoothers need only matrix-vector products,
which is why D35 chose them over Gauss–Seidel and ILU for the GPU.

**Mixed precision (D37).** In the default mode every matrix of the
hierarchy, including a float32 copy of the fine matrix used by the
fine-level smoother, the diagonals, the dense coarsest inverse and the
cycle vectors are float32. The residual enters the cycle rounded to
float32 and the correction leaves it widened to float64; CG, its residual,
its search directions and the solution are float64. The outer float64
iteration corrects the rounding of the preconditioner, so the converged
result is that of float64 (to the tolerance).

**Flexible CG** (Notay 2000): β = (z_{k+1}, r_{k+1} − r_k) / (z_k, r_k)
(the Polak–Ribière form). For an exactly symmetric preconditioner it is
the classical CG; it stays robust when the preconditioner is symmetric only
up to the float32 rounding. Four dot products per iteration, all in the
fixed block order of stage 2, so CPU and GPU give the same scalars.

**Kernels.** The cycle has eight kernel bodies in `inncraft-kernels::amg`
(smoother step, first smoother step, update, residual, product with a
rectangular matrix, prolongation and correction, dense row, precision
conversion), generic over the floating-point type and instantiated as PTX
entry points `amg_*_f32` and `amg_*_f64`. The R10 lint now allows f32 only
in the multigrid modules and rejects fused multiply-adds, approximate
division and flush-to-zero in f32 as in f64, so the CPU and GPU results
are expected to be bitwise equal in both precisions.

**hypre (D20).** Optional (cargo feature `hypre`): a C shim
(`crates/inncraft-fv/hypre/shim.c`) calls only hypre's documented user
interface (IJ matrix and vectors, ParCSR PCG, BoomerAMG) on the CSR matrix
that the in-house setup also starts from. hypre is licensed under MIT or
Apache-2.0, so linking it is compatible with Inncraft's MIT licence; no
hypre code is copied (D21), and a build without the feature does not need
it. `tools/build_hypre.sh` builds it from source (sequential, shared, CPU
or CUDA) for CI, Colab and local machines; no binary is redistributed.

## GPU against CPU and hypre (Tesla T4)

Run on a Google Colab Tesla T4 (sm_75, 15 GB, driver 580.82.07) with
`tools/stage3_colab.ipynb`. hypre 3.2.1 was built with CUDA for sm_75 and
without Umpire (the Colab image has no Umpire; hypre allocates with
`cudaMalloc`). The CPU column is 2 threads on the same Colab machine.

**912 673-cell polyhedral box, tolerance 1e-10:**

| Solver | Backend | Levels | Iterations | Setup (s) | Solve (s) | ns/cell/it. | Bitwise vs CPU |
|---|---|---:|---:|---:|---:|---:|---|
| amg | CPU, 2 threads | 4 | 12 | 4.408 | 2.151 | 196.4 | – |
| amg | CUDA, T4 | 4 | 12 | 4.415 | 0.072 | 6.5 | yes |
| amg-double | CUDA, T4 | 4 | 12 | 4.418 | 0.101 | 9.2 | yes |
| jacobi | CUDA, T4 | – | 329 | 0.000 | 0.603 | 2.0 | yes |
| hypre (BoomerAMG-PCG) | GPU | – | 23 | 1.430 | 0.221 | 10.5 | – |

The GPU solve is 30 times faster than the CPU solve. Per cell and iteration
the in-house AMG is about 1.6 times faster than hypre on the GPU. hypre needs
about twice the iterations and a third of the setup time, because its setup
runs on the GPU and ours still runs on the host (4.4 s, the same on CPU and
GPU).

**Polyhedral sweep on the T4 (`poisson_poly_box`, cells from 729 to 274 625):**

| Cells | amg it. | amg ns/cell/it. | amg-double ns/cell/it. | hypre it. | hypre ns/cell/it. |
|---:|---:|---:|---:|---:|---:|
| 729 | 10 | 290.1 | 316.7 | 17 | 1281.0 |
| 4 913 | 10 | 52.1 | 64.4 | 18 | 129.1 |
| 35 937 | 11 | 21.5 | 27.3 | 21 | 29.9 |
| 274 625 | 11 | 10.1 | 13.4 | 22 | 14.0 |

The in-house AMG keeps 10 to 11 iterations from 729 to 274 625 cells on
the GPU, as on the CPU. hypre grows from 17 to 22. The per-cell numbers
at 729 and 4 913 cells are latency bound (the kernels do not fill the GPU
at these sizes), so only the two largest sizes say something about
throughput. All GPU results are bitwise equal to the CPU results.

## Limits of stage 3

- The multigrid setup runs on the host (multi-threaded, float64) and is
  uploaded once; only the cycle runs on the GPU. A GPU setup, and the
  reuse of the aggregates when only the values of the matrix change (every
  outer iteration of SIMPLEC), follow with stage 4 and the performance
  work of stage 7.
- One partition only; the parallel aggregation across ranks (aggregates
  inside each rank, halo couplings in the coarse operators) follows with
  the multi-rank solver of stage 7.
- The coarsest level is a dense inverse (up to 4 000 unknowns). The pure
  Neumann pressure problem needs a reference value (D30) before the setup;
  a singular coarsest level is reported, not regularised.
- hypre is used sequentially (no MPI); its GPU build is checked on Colab
  only, and without Umpire, so its GPU setup time is not its best.

## References

- P. Vaněk, J. Mandel, M. Brezina, *Algebraic multigrid by smoothed
  aggregation for second and fourth order elliptic problems*, Computing 56,
  179–196 (1996)
- M. Adams, M. Brezina, J. Hu, R. Tuminaro, *Parallel multigrid smoothing:
  polynomial versus Gauss–Seidel*, J. Comput. Phys. 188, 593–610 (2003)
- Y. Notay, *Flexible conjugate gradients*, SIAM J. Sci. Comput. 22,
  1444–1460 (2000)
- Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM
  (2003)
- F. G. Gustavson, *Two fast algorithms for sparse matrices: multiplication
  and permuted transposition*, ACM Trans. Math. Softw. 4, 250–269 (1978)
- V. E. Henson, U. M. Yang, *BoomerAMG: a parallel algebraic multigrid
  solver and preconditioner*, Appl. Numer. Math. 41, 155–177 (2002)
- hypre: https://github.com/hypre-space/hypre

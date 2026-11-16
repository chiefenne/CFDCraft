# Stage 4 – laminar steady flow

Stage 4 of the [stage plan](../requirements.md#7-stage-plan) is the first
flow solver of Foehncraft: steady incompressible laminar flow with the
SIMPLEC algorithm (Van Doormaal & Raithby 1984) on the collocated,
cell-centred polyhedral grid, with the Rhie–Chow interpolation in its
relaxation-independent form (D19, D27), the boundary conditions of D29 with
the check of their combination (D30), the `foehn flow` case file with
"did you mean" validation (D39) and the log and CSV output (D43). The
outer iterations run on the CPU; every face and cell body is a single-source
kernel of `foehncraft-kernels` (`flow.rs`), so the GPU port of stage 5 reuses
them unchanged. This page explains the method, the checks and the results.
The commands are in [recipes 18–20](../recipes/flow.md#flow-verification), the
pictures in [gallery 12–13](../gallery/flow.md#cavity), the case files in the
[case catalogue](../cases/flow.md#flow-cases-examplescasesflow).

## Acceptance checks and results

| Check (stage plan row 4) | Result | Where |
|---|---|---|
| Couette: round-off error | relative L2 velocity error 4e-11 to 2e-10 on hexahedra (with convection) and on tetrahedra (prisms), polyhedra and mixed cells perturbed by 15 % (Stokes flow), with the outer iterations converged to 1e-11; the error follows the tolerance, not the mesh | cases `couette_*` (fail above 1e-9), test `couette_is_exact_on_all_mesh_types` |
| Poiseuille: observed order 1.8–2.2 on all mesh types | finest pair: hex (perturbed) 1.89, tet 1.80, poly (perturbed) 2.00, mixed 1.88 (table below) | cases `poiseuille_*` (fail outside 1.8–2.2), test `poiseuille_is_second_order_on_hexahedra` |
| Lid-driven cavity 2D, Re 100 and 1000, against Ghia et al. (1982) | largest deviation on the centre lines: Re 100 on 64 × 64 hexahedra 0.44 % (u) and 0.95 % (v) of the lid speed; Re 1000 on 128 × 128 0.40 % and 1.4 % | cases `cavity_re100`, `cavity_re1000` (fail above 2 %), test `cavity_re100_agrees_with_ghia` |
| Lid-driven cavity 3D, Re 100 and 1000 | the same flow on fully three-dimensional cells (a slab with symmetry planes front and back): Re 100 on 24 576 tetrahedra 2.5 % / 2.2 %, on 5 445 polyhedra 0.70 % / 1.3 %; Re 1000 on 28 227 polyhedra 2.0 % / 2.4 % | cases `cavity_re100_tet3d`, `cavity_re100_poly3d`, `cavity_re1000_poly3d` (fail above 3 %) |
| Rhie–Chow independent of the under-relaxation (D27) | lid force of the cavity with α_u = 0.7 and 0.9: equal to 3e-11 (relative), both converged to 1e-11 | test `converged_result_does_not_depend_on_relaxation` |
| Boundary-condition check and guidance (D30) | errors, warnings and notes with the boundary name and what to change; "did you mean" for boundary names, keys and values | `foehn flow --check`, tests `boundary_check_*`, `unknown_keys_and_values_get_a_hint`; guidance table in the [case catalogue](../cases/flow-keys.md#sensible-combinations-d30) |
| Reproducible (D31) | residuals and boundary forces bitwise equal on 1 and 4 threads (22 500 cells) | test `results_are_bitwise_equal_on_1_and_4_threads` |

### Couette and Poiseuille flow

Channel 0 < y < 1 of length 2, one cell layer in z with symmetry planes,
ρ = 1 kg/m³, μ = 0.1 Pa s (Re = 10), the exact profile at the inlet
(`profile = "analytic"`) and a static-pressure outlet. Couette: upper wall
moving at 1 m/s, linear profile, constant pressure; a linear field is
reproduced exactly by every part of the discretisation, so the error is
the iteration error. Poiseuille: both walls at rest, parabolic profile with
mean velocity 1 m/s, linear pressure drop 12 μ U / h². Errors relative to
the velocity scale, h the square root of the mean cell area.

| Mesh | Cells | Outer it. | L2 error | Order | p error (max, rel.) |
|---|---:|---:|---:|---:|---:|
| hex, perturbed 15 % | 512 | 170 | 1.60e-3 | | 1.6e-2 |
| | 2 048 | 169 | 4.07e-4 | 1.98 | 5.1e-3 |
| | 8 192 | 321 | 1.10e-4 | 1.89 | 2.2e-3 |
| tet (prisms) | 1 024 | 223 | 8.92e-4 | | 8.3e-3 |
| | 4 096 | 363 | 2.34e-4 | 1.93 | 3.5e-3 |
| | 16 384 | 1 191 | 6.70e-5 | 1.80 | 1.6e-3 |
| poly, perturbed 15 % | 561 | 175 | 1.16e-3 | | 7.6e-3 |
| | 2 145 | 177 | 3.03e-4 | 2.01 | 3.5e-3 |
| | 8 385 | 332 | 7.74e-5 | 2.00 | 1.9e-3 |
| mixed | 768 | 226 | 1.28e-3 | | 1.3e-2 |
| | 3 072 | 336 | 3.33e-4 | 1.94 | 5.8e-3 |
| | 12 288 | 1 077 | 9.05e-5 | 1.88 | 2.8e-3 |

Outer iterations to 1e-9 with α_u = 0.9. The prisms of the `tet` and
`mixed` meshes are not perturbed in these cases; on prisms perturbed by
15 % the order falls from 1.7 to 0.5 on the finest pair (see the
[limits](#limits-of-stage-4)).

### Lid-driven cavity

Unit square, lid (`ymax`) at 1 m/s, the other walls at rest, ρ = 1 kg/m³,
μ = 0.01 and 0.001 Pa s. No boundary fixes the pressure, so a reference
cell holds the level (the check says so in a note). Linear upwind for
momentum, no limiter, α_u = 0.9, converged to 1e-8. The deviation is the
largest difference to the tabulated centre-line velocities of Ghia et al.
(1982), in lid speeds, with the computed values interpolated linearly
(cell value plus gradient) at the reference points.

| Case | Cells | Outer it. | Time | max \|u − u_Ghia\| | max \|v − v_Ghia\| |
|---|---:|---:|---:|---:|---:|
| Re 100, hex 2D 64 × 64 | 4 096 | 439 | 5 s | 0.0044 | 0.0095 |
| Re 1000, hex 2D 128 × 128 | 16 384 | 1 397 | 54 s | 0.0040 | 0.0136 |
| Re 100, tet 3D 32 × 32 × 4 | 24 576 | 753 | 40 s | 0.0254 | 0.0224 |
| Re 100, poly 3D 32 × 32 × 4 | 5 445 | 192 | 4 s | 0.0070 | 0.0133 |
| Re 1000, poly 3D 96 × 96 × 2 | 28 227 | 1 240 | 129 s | 0.0198 | 0.0239 |

Times on one or four cores of a shared cloud machine (Intel Xeon, 4
virtual cores), 2 000 to 4 000 ns per cell and outer iteration.

## Methods

**Momentum** (per component, one matrix for the three components): implicit
first-order upwind convection and the orthogonal part of the diffusion,
with the linear-upwind correction ṁ ψ ∇u_U·(x_f − x_U) and the
non-orthogonal part −μ k·∇u_f deferred, as for the scalar of stage 2. The
deferred part also carries a skewness correction of the face gradient,
−μ D (∇u_N − ∇u_P)·(x_f − x_m), x_m the midpoint of the line P–N: the
difference u_N − u_P gives the exact derivative of a quadratic at x_m, not
at the face centroid. The pressure gradient enters as the source
−∇p_P V_P with the least-squares cell gradient. Under-relaxation is
implicit, a_P/α_u with (1 − α_u) a_P/α_u u_P on the right-hand side
(Patankar 1980). Each component is solved with BiCGStab and Jacobi to a
relative tolerance of 1e-2.

**Rhie–Chow face flux** (Rhie & Chow 1983) in the relaxation-independent
form of Majumdar (1988):

ṁ_f = ρ ū_f·S − ρ D̄_f (d·S)⁻¹|S|² [(p_N − p_P) − ∇p̄_f·(x_N − x_P)]
     + (1 − α_u)(ṁ_f^old − ρ ū_f^old·S),

with D_P = V_P/(a_P/α_u), ū_f the linear interpolation to the point of the
line P–N nearest to the face centroid plus the interpolated gradient
times the remaining offset, and ∇p̄_f the interpolated cell gradient. The
last term removes the dependence of the converged flux on α_u (D27); the
test above shows it.

**SIMPLEC** (Van Doormaal & Raithby 1984). Per outer iteration: (1)
boundary values; (2) gradients (and limiters if chosen); (3) the
relaxation term of the flux; (4) momentum assembly and solves; (5) the
predicted fluxes and the continuity defect; (6) the pressure correction
∇·(ρ D^C ∇p') = ∇·ṁ*, D^C_P = V_P/(a_P/α_u − Σ|a_N|), compact two-point
coefficients, solved with CG and the smoothed-aggregation multigrid of
stage 3 to 1e-2 (hierarchy rebuilt every 10 outer iterations); (7) the
corrections ṁ = ṁ* − c_f (p'_N − p'_P), u_P −= D^C_P ∇p'_P, p += α_p p'
(α_p = 1 for SIMPLEC). `algorithm = "simple"` uses D in place of D^C and
α_p = 0.3 by default. Without a pressure boundary the diagonal of the
reference cell is doubled; because the sum of the right-hand sides is zero
in a closed domain, the solution then has p' = 0 in that cell and
satisfies every continuity equation, and p' is shifted so that the
reference value is kept.

**Boundary conditions.** Each boundary face has a momentum kind (prescribed
velocity: walls and velocity inlets; outflow: zero-gradient outlets with
the inflowing part taking the boundary value; symmetry), a flux kind
(prescribed ρ u_b·S, or Rhie–Chow with the boundary pressure for pressure
boundaries) and update rules for the boundary velocity and pressure
(extrapolation along the face for outlets, normal velocity from the flux
for pressure inlets, total pressure p₀ − ½ ρ u_n²). On a symmetry plane the
diffusive flux is μ D (u_P·n̂) n̂: the part μ D n̂_i² u_P,i is implicit in
the diagonal of component i, the cross-component part is deferred. An
earlier version took the whole flux from the lagged face value
u_P − (u_P·n̂) n̂; on the thin slabs of the 2D cases that acted as a strong
extra under-relaxation of the tangential components and the cavity needed
about 1 700 outer iterations on any mesh, against 93 now for the same
Stokes case. The test `stokes_cavity_converges_quickly` keeps it fixed.

**Gradients.** Weighted least squares with per-variable boundary stencils:
the velocity uses the boundary values on walls and inlets and is
projected (zero normal gradient) on outlets and symmetry planes; the
pressure is projected on walls and symmetry planes and uses the boundary
value on pressure boundaries. Least squares is exact for linear fields;
on triangles and prisms next to a wall it is only first-order for a
quadratic profile, which is why the face gradient of the skewness
correction is interpolated rather than taken from one cell.

**Residuals.** Momentum: ‖b − A u‖ of each component before its solve,
over ‖a_P u‖ + ‖b‖ (finite from a fluid at rest). Continuity: the norm of
the net outflow of the predicted fluxes over the norm of the throughflow
of every cell. The case converges when all of them are below
`solver.tolerance`.

**Threads.** Below 20 000 cells the loops run on one thread (the thread
pool costs more than it saves); the results are bitwise the same either
way (gather assembly, sums in a fixed order).

## Outer iterations

SIMPLE-family algorithms need more outer iterations on finer meshes: in the
tables above the count grows about threefold per refinement on the
channel and with the Reynolds number in the cavity. It also falls with a
larger α_u (cavity Re 100 on 64²: 971 outer iterations with α_u = 0.8, 439
with 0.9, 218 with 0.95). This is the known behaviour of segregated
pressure–velocity coupling, not a solver defect; the coupled solver (L11)
is the planned remedy. A tighter momentum or pressure tolerance than 1e-2
does not change the count.

## Limits of stage 4

- CPU only. The kernels are single-source and checked by the R10 lint;
  the GPU flow loop and its Colab check follow in stage 5 together with
  PISO (owner decision, 2026-11-14).
- On prisms perturbed by 15 % (the `tet` and `mixed` 2D meshes with
  `perturb = 0.15`) the Poiseuille error stalls on the finest level
  (16 384 cells: L2 2.3e-4 instead of about 7e-5, pressure error 1 %),
  while perturbed hexahedra and polyhedra keep second order. Exact
  velocity gradients do not remove it, so it comes from the geometry
  (skewed triangles next to the walls and the outlet), not from the
  gradient. The acceptance cases use unperturbed prisms; this is an open
  item for the mesh-quality work of stage 7.
- The 3D cavity is the two-dimensional flow on three-dimensional cells
  (slab with symmetry planes), compared with Ghia et al. A cubic cavity
  with three-dimensional reference data belongs to the validation of
  stage 7 (D45).
- Periodic boundaries need the periodic face pairing and come with
  stage 6 (D29); the check reports them as an error until then.
- Turbulence values at inlets are accepted and reported as unused.
- One partition only; the multi-rank flow solver follows in stage 7.

## References

- C. M. Rhie, W. L. Chow, *Numerical study of the turbulent flow past an
  airfoil with trailing edge separation*, AIAA J. 21, 1525–1532 (1983)
- S. Majumdar, *Role of underrelaxation in momentum interpolation for
  calculation of flow with nonstaggered grids*, Numer. Heat Transfer 13,
  125–132 (1988)
- J. P. Van Doormaal, G. D. Raithby, *Enhancements of the SIMPLE method for
  predicting incompressible fluid flows*, Numer. Heat Transfer 7, 147–163
  (1984)
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere
  (1980)
- U. Ghia, K. N. Ghia, C. T. Shin, *High-Re solutions for incompressible
  flow using the Navier–Stokes equations and a multigrid method*,
  J. Comput. Phys. 48, 387–411 (1982)
- J. H. Ferziger, M. Perić, R. L. Street, *Computational Methods for Fluid
  Dynamics*, 4th ed., Springer (2020)
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in
  Computational Fluid Dynamics*, Springer (2016)

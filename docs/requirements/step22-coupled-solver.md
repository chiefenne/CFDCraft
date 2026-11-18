# Step 22 – Coupled pressure–velocity algorithm and pseudo-time steps: specification

Requirements step 22 (2026-11-18). The iteration budget of the
backward-facing step (gallery 74) showed that the segregated SIMPLEC needs
550 to 600 outer iterations on 8 200 cells and 1 500 to 1 750 on 32 800,
and that the count roughly triples when the spacing halves. The maintainer
chose the coupled pressure–velocity solver (L11, D19) as the next step,
with the targets "a case like this needs at most about 500 iterations" and
"we want to beat the industry standards: performance, speed, precision,
robustness", and with one constraint: the kinetic-energy conservation
needed later for LES (L39) must not be broken. This page is the
specification and the record of the measurements behind it. The
implementation is stage 6 part 6e
([verification](../verification/stage6.md#coupled-algorithm-and-pseudo-time-steps-part-6e)).

Clean-room rule (D21): methods from the cited papers only; other codes are
not named.

Contents

1. [Short answer](#1-short-answer)
2. [Why SIMPLEC needs so many iterations](#2-why-simplec-needs-so-many-iterations)
3. [The coupled algorithm (D176)](#3-the-coupled-algorithm-d176)
4. [Pseudo-time steps (D177, D178)](#4-pseudo-time-steps-d177-d178)
5. [Kinetic energy and LES](#5-kinetic-energy-and-les)
6. [Case keys](#6-case-keys)
7. [Measurements](#7-measurements)
8. [State of the art and what comes next](#8-state-of-the-art-and-what-comes-next)
9. [Items and their stages](#9-items-and-their-stages)
10. [References](#10-references)

## 1. Short answer

- **Two changes, not one.** The coupled solve alone halves the count
  (300 instead of 550 to 600 on the coarse step). The larger gain comes
  from replacing under-relaxation, which is a local time step that shrinks
  with the cell, by pseudo-time steps on physical scales: one step for all
  cells in the momentum equations and the time scale of the turbulence
  model in the turbulence equations. Together: 65 to 85 iterations on the
  coarse step for the three models instead of 525 to 675, 170 to 180 on
  the fine mesh instead of 1 750, 452 instead of 1 397 on the cavity at
  Re = 1000.
- **Same equations, same result.** The coupled algorithm solves the
  discrete equations of SIMPLEC in increment form. The converged result is
  the SIMPLEC one and depends neither on the step nor on the cap on the
  local step. The face fluxes satisfy continuity after every outer
  iteration.
- **No block matrix.** The coupled operator is applied matrix-free with
  the kernels that exist; the preconditioner is a SIMPLE-type sweep with
  the scalar multigrid of the pressure. The memory of the 4 × 4 block
  matrix and its hierarchy (section 2b of the requirements) is not needed.
- **What limits it today** is the momentum sweep of the preconditioner
  (BiCGStab with a diagonal preconditioner): it sets the cap on the local
  step, and the cap sets the count on fine meshes. It also makes a
  coupled iteration cost five to seven SIMPLEC iterations, so the run
  time of the step case gains a factor of about 1.5 where the iteration
  count gains 8 to 10 (section 7.5). A multigrid for
  convection-dominated matrices removes it (L44).

## 2. Why SIMPLEC needs so many iterations

Implicit under-relaxation with the factor α adds (1 − α)/α a_P to the
diagonal of a row. With a_P ≈ ρV/Δt_P of the cell (Δt_P the time a signal
needs to cross the cell by convection or diffusion) this is a pseudo-time
step

    Δτ_P = α/(1 − α) Δt_P,

a local Courant number of α/(1 − α): 2.3 at α = 0.7, 4 at 0.8, 9 at 0.9.
The step shrinks with the cell. The flow of the step case needs a fixed
physical time to establish its recirculation zone (some tens of H/U), so
the number of iterations grows as the mesh is refined, as measured. A
larger α helps only up to the stability limit of the segregated algorithm
(0.9 stalled the residual at 10⁻⁴ on the step).

The two-equation turbulence models add their own limit: the source terms
are stiff on the model time scale T (k/ε, 1/(β* ω)), and the lagged
coupling between the momentum equations and ν_t makes relaxation factors
above 0.9 oscillate.

## 3. The coupled algorithm (D176)

One outer iteration, from the fields u, p and the face fluxes ṁ:

1. Boundary values, gradients, limiters and the momentum rows A as in
   SIMPLEC, with the pseudo-time term of section 4.
2. The residuals of the current fields: r_u = b − A u − V ∇p per component
   and r_c = −Σ_f ±ṁ*_f, with ṁ* the Rhie–Chow fluxes of u and p (D27).
3. The increments δu, δp from one linear system,

       | A  G | |δu|   |r_u|        G δp = V ∇(δp),
       | D  L | |δp| = |r_c|        D δu = Σ_f ± ρ (δu)_f·S,

   with L δp = Σ_f c_f (δp_P − δp_N), c_f = ρ D̄_f |S|²/(d·S) the
   Rhie–Chow coefficient: the linearisation of the same flux formula. The
   interpolated pressure gradient of the Rhie–Chow term and the skewness
   corrections stay at the old fields (they are in r_c), so the system is
   an approximate Jacobian; the residuals are exact.
4. u += δu, p += δp, ṁ = ṁ* + ρ (δu)_f·S − c_f (δp_N − δp_P), then one
   pressure correction on the remaining continuity defect (the linear
   solve is loose, tolerance 0.1): the face fluxes satisfy continuity to
   the tolerance of `[solver.pressure]` after every outer iteration, as
   with SIMPLEC.
5. Turbulence and scalars as in SIMPLEC, with the fluxes of step 4.

**Same result as SIMPLEC.** The residuals are those of the segregated
algorithm, and the face flux depends neither on α (Majumdar 1988) nor on
the pseudo-time step (the time-step term of D27 with the current iterate
as the old level). Acceptance: the lid force of a cavity agrees with
SIMPLEC to 10⁻⁸ for the automatic step, a small step, a large step and
relaxation only (`tests/coupled.rs`).

**Linear solve.** Flexible GMRES (Saad 1993) with right preconditioning,
from a zero start, to a relative residual of 0.1. The operator is applied
with the operations of the device: the scalar momentum matrix per
component, the gradient with the mask of a pressure correction, the flux
of a velocity increment, the pressure matrix. The preconditioner is one
SIMPLE-type sweep on the increments (Patankar & Spalding 1972; as a
preconditioner Klaij & Vuik 2013):

    z_u = Ã⁻¹ r_u,   Ŝ z_p = r_p − D z_u,   z_u −= D̂ ∇z_p.

- Ã⁻¹: BiCGStab with a diagonal preconditioner per component, to 10⁻² or
  10 iterations.
- D̂ = min(D^C, V/((1 − α) a_P/α + ρV/(α Δτ))): the SIMPLEC diagonal of
  Van Doormaal & Raithby (1984), bounded by its value for fluxes that
  satisfy continuity (D^C is singular in cells where they do not yet).
- Ŝ: the pressure matrix with the face coefficient
  ρ (D̄ + D̂̄)_f |S|²/(d·S), the Rhie–Chow part L plus the compact form of
  −D D̂ G; CG with the multigrid of the pressure (D35), to 0.1.
- The continuity row is scaled with Σ a_P/(ρ Σ|S|) to the units of the
  momentum rows.

Findings behind these choices (measured on the matrices of the step case,
section 7.3): with an exact momentum inverse the SIMPLE-type Schur
complement needs 3 Krylov iterations for the first decade and the
pressure-convection–diffusion form (Kay, Loghin & Wathen 2002) more on
the stretched mesh; the momentum sweep must reach 10⁻² to 10⁻³, a looser
or truncated sweep makes the Krylov count and the outer iteration
erratic.

**Memory.** No block matrix and no block hierarchy: next to the buffers
of SIMPLEC, one more pressure matrix, 20 values per cell of work vectors
and 8 values per cell for every Krylov iteration (the basis and the
preconditioned basis, allocated as the solves need them; typically 5 to
10 iterations). At 100 million cells and 7 Krylov iterations these are
about 60 GB in float64, against 190 GB for the block matrix alone
(section 2b of the requirements) and its hierarchy on top. The Krylov
basis can later be held in single precision (D37), which halves its
part.

**Scope.** Steady runs. A transient coupled algorithm is not planned:
with time steps near CFL 1 the segregated step needs few corrections
(section 2b).

## 4. Pseudo-time steps (D177, D178)

### 4.1 Momentum: one step for all cells (D177)

The momentum rows get the term ρV/Δτ (u − uᵏ) with uᵏ the current
iterate: it vanishes at convergence. It enters through the time term of
the transient solver (D28), and the face flux through the
time-step-independent Rhie–Chow form of D27, so nothing new is
discretised.

The step is the same in every cell. Default (automatic):

    Δτ = factor × L / U,     L = ΣV / ΣA_b,

with A_b the areas of the boundary faces without symmetry planes
(periodic faces are interior faces) and U the larger of the largest
prescribed boundary speed (inlets, moving walls) and the root mean square
of the cell speeds, evaluated in every outer iteration. Outlets and
symmetry planes take their velocity from the cells and do not count: on a
mesh of one cell layer the symmetry planes carry the whole field, and its
largest velocity components shortened the step of the cavity by a fifth
(567 instead of 452 iterations). L is the half height of a plane channel, a quarter
of the diameter of a pipe, 7.3 H on the step case, and does not depend on
the thickness of a mesh of one cell layer, unlike the cube root of the
volume that is common practice in three dimensions. On the step the count
does not depend on the step between 5 and 29 H/U; the cost of the linear
solves grows with it (section 7.2), so the smaller length is the default
and `factor` scales it.

α_u stays as a cap on the local Courant number, α/(1 − α): default 0.98
(49). The cap is what the momentum sweep of the preconditioner can
afford: at 0.999 the count falls further (cavity: 166 instead of 423) but
BiCGStab with a diagonal preconditioner then costs more than it saves.

### 4.2 Turbulence: the time scale of the model (D178)

The transport equations of the turbulence model take the step

    Δτ_t = min(Δτ, c_T T)   per cell,

with T the time scale of the model: k/ε (k-ε family), 1/(β* ω) (SST) and
the bounded scale of k-ε-ζ-f; c_T = 0.5. The term ρV/Δτ_t (φ − φᵏ) is
added to the row before the relaxation; the residual scale does not
contain it. The relaxation of the transport equations is then 0.98 (a cap
as for the momentum equations).

Measured on the step: a global step for the turbulence equations is
limited to 0.5 H/U (1 and 2 oscillate); relaxation alone above 0.9
oscillates; with the time-scale rule c_T from 0.25 to 1 works (75, 55, 45
to 50 iterations at 0.25, 0.5, 1); no relaxation at all lets k-ε diverge
in the first iterations from a uniform start, and 0.9 makes the
turbulence lag the momentum equations (several times the iterations). The
rule follows the observation in the literature that the source terms and
the positivity of two-equation models limit the step (Mor-Yossef & Levy
2009).

With SIMPLEC the rule is off by default (results unchanged) and can be
switched on with `[solver.pseudo_time] turbulence = 0.5`: the coarse step
then needs about 230 iterations at α_u = 0.8 instead of 550 to 600.

### 4.3 No start-up phase

The pseudo-time steps act from the first iteration, from a uniform start.
The relaxation ramp of the segregated cases (0.5, then 0.8) is not
needed.

## 5. Kinetic energy and LES

The maintainer's constraint: the kinetic-energy conservation needed for
LES (L39) must hold with any new algorithm. What the coupled algorithm
changes and what it does not:

- **Not the discrete equations.** Convection scheme, gradients,
  Rhie–Chow flux and boundary conditions are those of SIMPLEC and PISO.
  What the scheme conserves it conserves with either algorithm; the
  `kep` scheme of L39 needs no change for it.
- **Not the continuity of the fluxes.** After every outer iteration the
  face fluxes satisfy continuity to the tolerance of the pressure solve
  (step 4 of section 3), the condition for a skew-symmetric convection
  term.
- **The pressure in increment form**, as in the PISO (D27): the
  kinetic-energy error of the Rhie–Chow term is that of the segregated
  solver.

Acceptance check with the ledger (D88): on the cavity at Re = 1000 the
numerical dissipation of kinetic energy is 0.637 % of the lid power with
both algorithms, and its pressure part (the Rhie–Chow term) agrees to
seven digits (section 7.4). A transient coupled algorithm, if one is ever
wanted for LES, has to pass the inviscid periodic test of L39 with the
same ledger.

## 6. Case keys

```toml
[solver]
algorithm = "coupled"        # steady; default "simplec"

[solver.momentum]
relaxation = 0.98            # default with "coupled": the cap α/(1 − α) = 49

[solver.coupled]             # the linear solve
tolerance = 0.1              # relative residual of the coupled system
max_iter = 30                # Krylov iterations per outer iteration
momentum_tolerance = 0.01    # momentum sweep of the preconditioner
momentum_max_iter = 10
pressure_tolerance = 0.1     # pressure solve of the preconditioner
pressure_max_iter = 20

[solver.pseudo_time]
# step = 5.0                 # Δτ in s; default: automatic
factor = 1.0                 # on the automatic step
momentum = true              # false: relaxation only
# turbulence = 0.5           # c_T; default 0.5 with "coupled", 0 (off) otherwise
```

The log shows the Krylov iterations of the coupled solve in the column
`K` (`linear.coupled`), next to the totals of the momentum sweeps and the
pressure solves; `residuals.csv` and the run record carry them as
`krylov`.

## 7. Measurements

Backward-facing step of Driver & Seegmiller (gallery 74), uniform start,
iterations until the reattachment length x_r stays within 1 % of its
final value. Coarse: 8 200 cells; fine: 32 800 cells.

### 7.1 Outer iterations

| Case | SIMPLEC (0.5, then 0.8) | Coupled, relaxation only | Coupled with pseudo-time steps |
|---|---|---|---|
| k-ε-ζ-f₀, coarse | 525 | 300 | 65 |
| k-ε, coarse | 600 | | 75 |
| k-ω SST, coarse | 675 | | 85 |
| k-ε-ζ-f₀, fine | 1 750 | 700 to 800 | 170 to 180 |
| Channel Re_τ = 395, k-ε-ζ-f₀, residual 10⁻⁶ | about 2 100 | | 323 |
| Cavity Re = 1000, 128², residual 10⁻⁸ | 1 397 (α = 0.9) | | 452 |

The final x_r is the one of SIMPLEC (6.385, 5.287, 6.81 H on the coarse
mesh), the bulk velocity of the channel too (U_b⁺ = 17.79). The fine mesh
needs more than twice the count of the coarse one because of the cap on
the local step (section 7.2). The counts are good to about 10 %
(section 7.6).

Small laminar cases are the exception: the cavity at Re = 400 on 48²
cells takes 255 iterations with SIMPLEC at α = 0.9 and 292 with the
coupled algorithm and the automatic step; `factor` 4 and 20 give 130 and
84 iterations at a higher cost per iteration (at Re = 1000 on 64²: 584
iterations in 8 s with SIMPLEC; 625 in 13 s, 276 in 24 s, 178 in 43 s
with `factor` 1, 4, 20). The automatic step is a careful one, and the
momentum sweep makes larger steps expensive (L44).

### 7.2 The step and the cap

k-ε-ζ-f₀, coarse mesh, cap 0.98, 150 iterations each:

| Δτ U/H | Iterations to 1 % | Krylov iterations per outer iteration |
|---|---|---|
| 2 | 75 | 4.6 |
| 4.7 | 65 to 70 | 7.2 |
| 7.3 (automatic) | 65 to 70 | 8.0 |
| 29 | 65 to 70 | 12.5 |
| none (relaxation 0.98 only) | 65 to 70 | 17.9 |

The cap at Δτ = 5 H/U: 0.95 needs 85 iterations (3.3 Krylov iterations
each), 0.98 65 to 70 (7), 0.999 60 to 65 (13 and more). On the fine mesh
0.95 needs about 200 and 0.98 145. On the cavity (Δτ = 0.25 s given): 0.95
1 012, 0.98 423, and 0.999 with Δτ = 1 s 166 iterations, at four times
the run time of 0.98; a tighter coupled solve (0.01) changes nothing
(167).

### 7.3 The momentum sweep

Momentum matrix of the step at Δτ = 5 H/U without a cap (diagonal
dominance 1.03), BiCGStab iterations to 10⁻² / 10⁻³:

| Preconditioner | Iterations |
|---|---|
| Diagonal | 113 to 137 |
| ILU(0) | 43 to 59 |
| Symmetric Gauss–Seidel | 21 to 24 |
| Aggregation multigrid of the symmetric part | 90 to 194 |
| Approximate ideal restriction multigrid (Manteuffel, Ruge & Southworth 2018) | 1 |

The aggregation multigrid of the pressure does not work for
convection-dominated matrices; the reduction-based multigrid does, at an
operator complexity of 4. This is L44.

### 7.4 Kinetic-energy ledger

Cavity at Re = 1000, 128² cells, converged to 10⁻⁸ (W per m of depth
times the cell thickness; lid power 2.7628 × 10⁻⁴):

| | SIMPLEC | Coupled |
|---|---|---|
| Viscous dissipation | 2.74524 × 10⁻⁴ | 2.74525 × 10⁻⁴ |
| Numerical dissipation | 1.7605 × 10⁻⁶ (0.637 %) | 1.7593 × 10⁻⁶ (0.637 %) |
| of it pressure (Rhie–Chow) | 3.996223 × 10⁻⁶ | 3.996224 × 10⁻⁶ |
| of it schemes | −2.2357 × 10⁻⁶ | −2.2369 × 10⁻⁶ |
| Deviation from Ghia et al. (u, v) | 0.0040, 0.0136 | 0.0040, 0.0136 |

### 7.5 Run time

k-ε-ζ-f₀ on the step, one thread, each run to the iteration from which
x_r stays within 1 %:

| Mesh | SIMPLEC (0.5, then 0.8) | Coupled | Time per outer iteration |
|---|---|---|---|
| Coarse, 8 200 cells | 525 iterations, 31 s | 65 iterations, 21 s | 0.06 s and 0.32 s |
| Fine, 32 800 cells | 1 750 iterations, 342 s | 175 iterations, 235 s | 0.20 s and 1.35 s |

A coupled iteration costs five to seven SIMPLEC iterations: it takes 10
(coarse) and 14 (fine) Krylov iterations, each with a momentum sweep that
runs into its limit of 10 iterations and with a pressure solve. So the
run time gains a factor of about 1.5 where the iteration count gains 8 to
10. The sweep is where the run time is to be won: it decides the number
of Krylov iterations and the cap on the local step (sections 7.2 and
7.3), and L44 replaces it.

### 7.6 Scatter of the counts

A run repeats bit by bit for the same input, on any number of threads
and across a restart. But the count reacts to the last digit of an
input. The cavity needs 423 iterations with Δτ = 0.25 s given and 452
with the automatic step, which is 0.25 s less one unit in the 17th
digit; the fine step needed 145 iterations with a slightly shorter step
(a velocity scale that counted the symmetry planes) and needs 170 to 180
at 7.3 H/U. The cause: the coupled solve stops at a tolerance of 0.1,
its momentum sweep is cut off after 10 iterations and the multigrid
cycle of its pressure solve runs in single precision, so two runs
separate in the fifth digit within one outer iteration and reach the
criterion on different paths. The converged result is the same
(section 7.4). Counts are therefore quoted to about 10 %; a smooth
dependence of the count on the step needs a converging sweep (L44) and
the step control of L45.

## 8. State of the art and what comes next

Published reference points for pressure-based coupled solvers:

- Fully coupled block systems with a block algebraic multigrid (Darwish,
  Sraj & Moukalled 2009; Mangani, Buchmayr & Darwish 2014): on a
  backward-facing step with 12 000, 48 000 and 195 000 cells the coupled
  solver was 7.6, 16.3 and 42 times faster than the segregated one, with a
  pseudo-time step that does not change with the mesh; the turbulence
  equations stay segregated.
- SIMPLE-type block preconditioners with flexible GMRES for
  Reynolds-averaged flows at high Reynolds numbers (Klaij & Vuik 2013):
  sub-solves to 10⁻², some 20 Krylov iterations.
- Common practice of the industry standards for the pseudo-time step: one
  automatic step from a length of the domain over a velocity scale, a
  factor near 1 (3 to 10 can speed a run up), the turbulence equations on
  the same step.

Where Foehncraft stands after this step: the iteration counts are in the
published range on the coarse mesh (7 to 9 times fewer iterations than
SIMPLEC), the count still grows with the mesh (170 to 180 on the fine
mesh against 65), and the gain in run time is far smaller than the gain in
iterations (a factor of about 1.5 on the step, section 7.5, where the
published coupled solvers with a block multigrid report 8 to 40) because
of the momentum sweeps. The turbulence step from the model time scale is our own
addition; the published solvers leave the turbulence equations on the
global step.

What closes the gap, in this order:

- **L44, multigrid for convection-dominated matrices** (approximate ideal
  restriction): for the momentum sweep and the turbulence equations. It
  removes the cap on the local step, so the count no longer grows with
  the mesh, and makes each Krylov iteration cheap. Clean-room from
  Manteuffel, Ruge & Southworth (2018); CPU and GPU kernels.
- **L45, automatic step control and a stop on monitored values**: the
  step grows as the residual falls (switched evolution relaxation, Mulder
  & van Leer 1985; Kelley & Keyes 1998) and falls back after a failed
  iteration; the run ends when monitored values (forces, x_r, a probe)
  have stopped changing, because the residual is no stop measure on such
  cases (stage 6 verification, part 6d2).
- **L46, Newton–Krylov**: the exact linearisation of the Rhie–Chow flux
  and of the convection term by matrix-free products, for quadratic
  convergence near the solution.

## 9. Items and their stages

Rows in [`requirements.md`](../requirements.md).

| ID | Item | Stage |
|---|---|---|
| D176 | Coupled pressure–velocity algorithm (steady), matrix-free, flexible GMRES with a SIMPLE-type preconditioner; `algorithm = "coupled"`, `[solver.coupled]` | stage 6 part 6e |
| D177 | Pseudo-time step of the momentum equations, automatic from ΣV/ΣA_b over a velocity scale; α_u as the cap on the local Courant number; `[solver.pseudo_time]` | stage 6 part 6e |
| D178 | Pseudo-time step of the turbulence equations from the model time scale, min(Δτ, c_T T) | stage 6 part 6e |
| L44 | Multigrid for convection-dominated matrices (approximate ideal restriction) for the momentum sweep and the turbulence equations | the next solver step |
| L45 | Automatic step control (switched evolution relaxation) and a stop on monitored values | after L44 |
| L46 | Newton–Krylov on the coupled system | after L45 |

## 10. References

- S. V. Patankar, D. B. Spalding, *A calculation procedure for heat, mass
  and momentum transfer in three-dimensional parabolic flows*, Int. J.
  Heat Mass Transfer 15, 1787–1806 (1972)
- J. P. Van Doormaal, G. D. Raithby, *Enhancements of the SIMPLE method
  for predicting incompressible fluid flows*, Numer. Heat Transfer 7,
  147–163 (1984)
- C. M. Rhie, W. L. Chow, *Numerical study of the turbulent flow past an
  airfoil with trailing edge separation*, AIAA J. 21, 1525–1532 (1983)
- S. Majumdar, *Role of underrelaxation in momentum interpolation for
  calculation of flow with nonstaggered grids*, Numer. Heat Transfer 13,
  125–132 (1988)
- Y. Saad, *A flexible inner-outer preconditioned GMRES algorithm*, SIAM
  J. Sci. Comput. 14, 461–469 (1993)
- D. Kay, D. Loghin, A. Wathen, *A preconditioner for the steady-state
  Navier–Stokes equations*, SIAM J. Sci. Comput. 24, 237–256 (2002)
- C. M. Klaij, C. Vuik, *SIMPLE-type preconditioners for cell-centered,
  colocated finite volume discretization of incompressible
  Reynolds-averaged Navier–Stokes equations*, Int. J. Numer. Meth. Fluids
  71, 830–849 (2013)
- M. Darwish, I. Sraj, F. Moukalled, *A coupled finite volume solver for
  the solution of incompressible flows on unstructured grids*, J. Comput.
  Phys. 228, 180–201 (2009)
- L. Mangani, M. Buchmayr, M. Darwish, a fully coupled pressure-based
  solver for steady incompressible turbulent flows with a block algebraic
  multigrid, Numer. Heat Transfer B 66, 1–20 (2014)
- T. A. Manteuffel, J. Ruge, B. S. Southworth, *Nonsymmetric algebraic
  multigrid based on local approximate ideal restriction (ℓAIR)*, SIAM J.
  Sci. Comput. 40, A4105–A4130 (2018)
- W. A. Mulder, B. van Leer, *Experiments with implicit upwind methods
  for the Euler equations*, J. Comput. Phys. 59, 232–246 (1985)
- C. T. Kelley, D. E. Keyes, *Convergence analysis of pseudo-transient
  continuation*, SIAM J. Numer. Anal. 35, 508–523 (1998)
- Y. Mor-Yossef, Y. Levy, *The unconditionally positive-convergent
  implicit time integration scheme for two-equation turbulence models:
  revisited*, Computers & Fluids 38, 1984–1994 (2009)
- D. M. Driver, H. L. Seegmiller, *Features of a reattaching turbulent
  shear layer in divergent channel flow*, AIAA J. 23, 163–171 (1985)

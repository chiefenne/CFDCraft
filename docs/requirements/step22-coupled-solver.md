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
([verification](../verification/stage6.md#coupled-algorithm-and-pseudo-time-steps-part-6e));
the multigrid of the momentum sweep (D179, sections 3.1 and 7.7) is
part 6f
([verification](../verification/stage6.md#multigrid-for-the-momentum-sweep-part-6f)).

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
  model in the turbulence equations. Together: 65 to 80 iterations on the
  coarse step for the three models instead of 525 to 675, 175 on
  the fine mesh instead of 1 600, 452 instead of 1 397 on the cavity at
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
- **The momentum sweep of the preconditioner** was what limited part 6e
  (BiCGStab with a diagonal preconditioner; the run time of the step
  case gained a factor of 1.4 to 1.5 where the iteration count gained 8
  to 9). Since part 6f it is one cycle of a reduction multigrid for
  convection-dominated matrices (D179, sections 3.1 and 7.7): half the
  Krylov iterations, 150 instead of 175 iterations on the fine step, and
  a run time 2.7 (coarse) and 3.3 (fine) times shorter than SIMPLEC.
- **What limits it today**: the turbulence model and the pressure
  solves take 60 % of a coupled run, the setups of the momentum
  hierarchies cost more than their cycles, and the cap on the local
  step has to stay because of the Schur complement approximation
  (section 8).

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

- Ã⁻¹: one cycle of a reduction multigrid per component (D179,
  section 3.1). Part 6e used BiCGStab with a diagonal preconditioner, to
  10⁻² or 10 iterations; it remains as `momentum_sweep = "bicgstab"`.
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
10 iterations), and the hierarchies of the momentum sweeps in float32
(section 3.1; one to three, each two to five times the entries of the
momentum matrix). At 100 million cells and 7 Krylov iterations these are
about 60 GB in float64, against 190 GB for the block matrix alone
(section 2b of the requirements) and its hierarchy on top. The Krylov
basis can later be held in single precision (D37), which halves its
part.

**Scope.** Steady runs. A transient coupled algorithm is not planned:
with time steps near CFL 1 the segregated step needs few corrections
(section 2b).

### 3.1 The momentum sweep: a reduction multigrid (D179)

With large pseudo-time steps the momentum rows are far from diagonally
dominant, and they are not symmetric: upwind convection makes them
nearly triangular along the flow. Point smoothers and the aggregation
multigrid of the pressure rely on a symmetric positive definite matrix
and fail on them (section 7.3). A reduction multigrid uses the
triangular structure instead (Ries, Trottenberg & Winter 1983;
Manteuffel, Ruge & Southworth 2018):

- The unknowns are split into coarse points C and fine points F such
  that the block A_FF is strongly diagonally dominant. With the exact
  Schur complement on the C points, one Jacobi sweep on the F points
  would solve the system.
- The Schur complement is approximated by the Galerkin operator
  A_c = R A P, with the restriction R = [Z I],
  Z = −A_CF (2 D_F⁻¹ − D_F⁻¹ A_FF D_F⁻¹): the ideal restriction
  −A_CF A_FF⁻¹, which makes the coarse residual equation independent of
  the error on the F points, with two terms of the Neumann series of
  A_FF⁻¹, on the C–F entries above 0.05 of the largest entry of their
  row. The interpolation P = [W; I] is the direct interpolation of
  classical multigrid (Stüben 2001); it keeps the row sums, so the same
  hierarchy serves where diffusion dominates.
- The same again on A_c, down to at most 64 unknowns, which a dense
  inverse solves.

**Setup** (host, float64, the same on any number of threads, D31), per
level: strong dependencies by the classical measure, j is strong for i
when |a_ij| ≥ θ max_{k≠i} |a_ik|, θ = 0.3; the first pass of the
Ruge–Stüben splitting (the undecided point on which most points depend
becomes C, the points that depend on it become F; ties go to the
smallest index); R, P and A_c, of which the entries below 10⁻⁴ of the
diagonal of their row are dropped (not lumped).

**Cycle.** The V-cycle and the level operations of the scalar multigrid
(D35, D37) are used unchanged, on the CPU and on the GPU, with the
hierarchy in float32: the hierarchy stores an inverse diagonal that is
zero on the C points, so the Jacobi step of the cycle relaxes the F
points only, once before and once after the coarse correction. One
cycle from a zero start is a momentum sweep
(`momentum_cycles`). Its cost is 5 to 8 products with the momentum
matrix.

**What the studies behind the design showed** (a prototype on matrices
written out of step and cavity runs, then the Rust code; section 7.7):

- The number of F sweeps per level must be even with the two-term
  series (an odd number leaves error modes with a factor near −1): one
  before and one after the coarse correction.
- Direct interpolation instead of the one-point interpolation of the
  paper: 0.32 instead of 0.63 per cycle where diffusion dominates.
- A drop tolerance of 10⁻³ is enough with the cap on the local step and
  not without it on closed streamlines (0.74 to 1.01 per cycle on the
  rotating-flow test matrices; 0.40 to 0.53 with 10⁻⁴).
- Relaxing the C points as well, lumping the dropped entries, other
  strength thresholds, a parallel independent-set splitting, a clean-up
  of weakly dominant F rows and restrictions from local dense solves
  gave nothing for their cost.

**Reuse.** A setup costs as much as twenty to thirty cycles, and the momentum
matrix changes little from one outer iteration to the next once the
flow has formed. A hierarchy is kept while
Σ|a_ij − a_ij⁰| ≤ 0.02 Σ|a_ij⁰| (`momentum_reuse`; a⁰ the matrix of its
setup): measured on the step case, a hierarchy from a matrix 1.3 %
away converges with 0.27 per cycle, 2.7 % with 0.34, 5.4 % with 0.6, and
from 6.4 % on it does not converge. As a guard, a kept hierarchy whose
first cycle of an outer iteration leaves more than half of the residual
is rebuilt at once. A restart file holds the matrices of the kept
hierarchies, so a continued run builds the same ones and stays
bit-identical.

**Details.** Each velocity component has its own hierarchy when the
matrices differ (the diagonals differ at symmetry planes); components
with the same matrix share one. A component whose part of a Krylov
vector is zero up to rounding (the third component of a planar case)
takes no sweep and needs no hierarchy. If a setup fails (a zero
diagonal, a coarsening that stalls), BiCGStab does the sweeps of that
outer iteration and the log carries one warning.

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
momentum_sweep = "multigrid" # or "bicgstab"
momentum_cycles = 1          # multigrid cycles per sweep
momentum_reuse = 0.02        # keep a hierarchy up to this matrix change; 0: never
momentum_tolerance = 0.01    # a BiCGStab sweep
momentum_max_iter = 10
pressure_tolerance = 0.1     # pressure solve of the preconditioner
pressure_max_iter = 20

[solver.coupled.multigrid]   # the multigrid of the momentum sweeps (D179)
strength = 0.3               # θ of the splitting and the interpolation
restriction_strength = 0.05  # θ of the entries that enter the restriction
drop_tolerance = 1e-4        # coarse entries below this part of the diagonal go
max_coarse = 64              # unknowns of the coarsest level (dense inverse)
max_levels = 30
precision = "mixed"          # float32 hierarchy; "double": float64

[solver.pseudo_time]
# step = 5.0                 # Δτ in s; default: automatic
factor = 1.0                 # on the automatic step
momentum = true              # false: relaxation only
# turbulence = 0.5           # c_T; default 0.5 with "coupled", 0 (off) otherwise
```

The log shows the Krylov iterations of the coupled solve in the column
`K` (`linear.coupled`), next to the totals of the momentum sweeps (the
multigrid cycles, or the BiCGStab iterations) and the pressure solves;
`residuals.csv` and the run record carry them as `krylov`. The report of
a run counts the setups of the momentum hierarchies.

## 7. Measurements

Backward-facing step of Driver & Seegmiller (gallery 74), uniform start,
iterations until the reattachment length x_r stays within 1 % of its
final value. Coarse: 8 200 cells; fine: 32 800 cells.

Sections 7.1 to 7.6 are the record of part 6e: they were measured with
the BiCGStab momentum sweep (`momentum_sweep = "bicgstab"` today).
Section 7.7 repeats the main cases with the multigrid sweep of D179, the
default since part 6f, and holds the run times of the current code.

### 7.1 Outer iterations

| Case | SIMPLEC (0.5, then 0.8) | Coupled, relaxation only | Coupled with pseudo-time steps |
|---|---|---|---|
| k-ε-ζ-f₀, coarse | 525 | 300 | 65 |
| k-ε, coarse | 600 | | 75 |
| k-ω SST, coarse | 675 | | 80 |
| k-ε-ζ-f₀, fine | 1 600 | 700 to 800 | 175 |
| Channel Re_τ = 395, k-ε-ζ-f₀, residual 10⁻⁶ | about 2 100 | | 323 |
| Cavity Re = 1000, 128², residual 10⁻⁸ | 1 397 (α = 0.9) | | 452 |

The final x_r is the one of SIMPLEC (6.385, 5.287, 6.81 H on the coarse
mesh, 6.352 H on the fine one), the bulk velocity of the channel too (U_b⁺ = 17.79). The fine mesh
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
| Fine, 32 800 cells | 1 600 iterations, 327 s | 175 iterations, 235 s | 0.20 s and 1.35 s |

A coupled iteration costs five to seven SIMPLEC iterations: it takes 10
(coarse) and 14 (fine) Krylov iterations, each with a momentum sweep that
runs into its limit of 10 iterations and with a pressure solve. So the
run time gains a factor of 1.4 to 1.5 where the iteration count gains 8 to
9. The sweep is where the run time is to be won: it decides the number
of Krylov iterations and the cap on the local step (sections 7.2 and
7.3), and L44 replaces it.

These are the times of part 6e. Section 7.7 has those of the current
code, which are shorter for both algorithms.

### 7.6 Scatter of the counts

A run repeats bit by bit for the same input, on any number of threads
and across a restart. But the count reacts to the last digit of an
input. The cavity needs 423 iterations with Δτ = 0.25 s given and 452
with the automatic step, which is 0.25 s less one unit in the 17th
digit; the fine step needed 145 iterations with a slightly shorter step
(a velocity scale that counted the symmetry planes) and needs 175
at 7.3 H/U. The cause: the coupled solve stops at a tolerance of 0.1,
its momentum sweep is cut off after 10 iterations and the multigrid
cycle of its pressure solve runs in single precision, so two runs
separate in the fifth digit within one outer iteration and reach the
criterion on different paths. The converged result is the same
(section 7.4). Counts are therefore quoted to about 10 %; a smooth
dependence of the count on the step needs a converging sweep (L44) and
the step control of L45.

### 7.7 The multigrid sweep (D179)

**The cycle.** Reduction of the residual per cycle (stationary
iteration, asymptotic) and operator complexity (the entries of all
levels over those of the matrix):

| Matrices | Factor per cycle | Operator complexity |
|---|---|---|
| Upwind convection–diffusion in a solid-body rotation, 64² to 256² cells, cell Péclet number 1 to 300 at the rim, diagonal raised by 2 % (the cap 0.98) | 0.07 to 0.22 | 2.4 to 3.5 |
| The same without the raised diagonal (closed streamlines, no cap), Péclet number 5 and 30 | 0.25 on 64², 0.40 and 0.53 on 128² | 3.2 to 3.5 |
| Diffusion only, 32² | 0.26 | 2.4 |
| Scalar transport on three-dimensional meshes of hexahedra, tetrahedra, polyhedra and a mix, cell Péclet number 0.5 and 20 | 0.04 to 0.28 | 2.2 to 2.7 on polyhedra, 3.5 to 5.4 on the others |
| Momentum matrices of the step case (8 200 and 32 800 cells) and of the cavity, from runs | 0.13 to 0.45 | 2.1 to 2.8 |

With the cap the factor does not grow with the mesh (0.15 on 128², 0.16
on 256² at a Péclet number of 30). Without it, on closed streamlines, it
does (0.25 and 0.40 on 64² and 128²): only the diffusion across the
streamlines damps an error along them. The operator complexity on
hexahedra and tetrahedra in three dimensions is the weak point: the
coarse operators of a reduction method get wide stencils there, and
dropping more entries costs convergence. (`amg::reduction` tests and
`tests/reduction.rs`; the last row from the prototype, with a drop
tolerance of 10⁻³.)

**The runs.** One thread; iterations until x_r stays within 1 % of its
final value (step) or to the tolerance (cavity 10⁻⁸); Krylov iterations
per outer iteration; run time to that iteration:

| Case | SIMPLEC | Coupled, BiCGStab sweep | Coupled, multigrid sweep |
|---|---|---|---|
| Step, k-ε-ζ-f₀, 8 200 cells | 525; –; 18.7 s | 65; 10.0; 8.3 s | 65; 5.3; 7.0 s |
| Step, k-ε-ζ-f₀, 32 800 cells | 1 600; –; 243 s | 175; 14.2; 131 s | 150; 6.0; 74 s |
| Cavity, Re = 1000, 128² | 1 397; –; 28.3 s | 452; 2.2; 18.8 s | 449; 2.1; 12.1 s |

With the multigrid sweep k-ε and k-ω SST settle on the coarse step in
75 and 80 iterations as before (2.2 and 5.1 Krylov iterations per outer
iteration over 160 iterations), the channel converges in 324 iterations
(323). The final values are those of section 7.1 (x_r = 6.385, 5.287,
6.81 H; 6.35 H on the fine mesh; the cavity deviates from Ghia et al.
by 0.0040 and 0.0136 as before).

Where the time goes, in seconds (the phases of the run report; "other"
is the turbulence model, the start and the output):

| Run | Momentum sweeps | Setups of their hierarchies | Pressure solves | Assembly and the rest of the coupled solve | Other |
|---|---|---|---|---|---|
| Coarse, 65 iterations, BiCGStab sweep | 1.8 | – | 2.6 | 1.1 | 2.8 |
| Coarse, 65 iterations, multigrid sweep | 0.3 | 1.4 (96) | 1.6 | 1.1 | 2.6 |
| Fine, 175 iterations, BiCGStab sweep | 33.9 | – | 39.5 | 14.3 | 42.9 |
| Fine, 150 iterations, multigrid sweep | 4.1 | 7.3 (140) | 18.5 | 12.8 | 31.4 |

"Other" split up in a repeat of the two multigrid runs (7.1 s and
70.5 s; the report now gives the time of the turbulence model):

| Run | Turbulence model | of it: linear solves | of it: assembly, gradients, sources, wall treatment | Start, output |
|---|---|---|---|---|
| Coarse, 65 iterations | 2.3 | 1.3 | 1.0 | 0.4 |
| Fine, 150 iterations | 24.5 | 14.8 | 9.8 | 4.8 |

The linear solves of the fine run: f with CG and the aggregation
multigrid 18 iterations per solve and 5.8 s; k, ε and ζ with
Jacobi-BiCGStab 33, 24 and 39 iterations per solve and 3.1 s, 2.2 s and
3.7 s (coarse: 14 for f, 25, 15 and 31).

Findings:

1. One cycle is as good a momentum inverse for the Krylov solve as the
   BiCGStab sweep was meant to be: the count per outer iteration halves
   (10 → 5 on the coarse step, 14 → 6 on the fine one), and with it the
   pressure solves. A sweep costs a third of a BiCGStab sweep that
   runs into its limit of 10 iterations.
2. The fine mesh settles in 150 iterations instead of 175 and in 0.57
   of the time; against SIMPLEC the run time gains a factor of 3.3
   on the fine mesh and 2.7 on the coarse one (2.3 on the cavity).
3. The setups cost more than the cycles they serve: a setup takes as
   long as twenty to thirty cycles, and while the flow forms the matrix
   changes by more than 2 % in most outer iterations (96 setups for two
   components in 65 iterations on the coarse step, 140 in 150 on the
   fine one; the cavity needs 28 in 449). On the coarse mesh setups
   and cycles together cost what the BiCGStab sweeps did, and the gain
   is the halved Krylov count. A cheaper setup (the host code builds
   its matrices row by row with many small allocations) and a refresh
   of the values on a kept splitting are the next steps here.
4. The cap on the local step stays. Without it (relaxation 0.999) the
   fine step settles in about 120 iterations, but the Krylov solve runs
   into its limit of 30 in every outer iteration. With exact sub-solves
   the SIMPLE-type preconditioner still needs 10 Krylov iterations on
   these matrices (4 to 6 with the cap), and 51 on the cavity with a
   step of 20 times the automatic one: the approximation of the Schur
   complement limits then, not the sweep alone (section 8).
5. The momentum sweep is no longer what limits a coupled run, and no
   single part is. On the fine mesh the turbulence model takes 35 % of
   the time (its four linear solves 21 %, the rest of it 14 %), the
   pressure solves 25 %, the assembly and the rest of the coupled solve
   16 %, the setups of the hierarchies 11 %.
6. Large steps on laminar cases are no longer expensive. Cavity at
   Re = 1000 on 64²: SIMPLEC at α = 0.9 584 iterations in 3.0 s; coupled
   with `factor` 1, 4, 20: 624, 275, 170 iterations in 3.7 s, 2.8 s,
   2.6 s (with the BiCGStab sweep the same three runs took 1, 1.9 and
   3.3 times the time of `factor` 1, section 7.1). On 128² the cap
   limits: 449, 358, 290 iterations in 12.1 s, 15.9 s, 13.5 s, so the
   automatic step stays the default.

**One thread without a pool.** These measurements showed a loss that
has nothing to do with the algorithm: on a pool of one thread (meshes
below 20 000 cells always run on one) every kernel call was still
handed to the worker thread and waited for, about 2 000 hand-offs per
outer iteration of the cavity. A pool of one thread and loops of
one task (the coarse levels of a multigrid cycle) now run on the
calling thread, with the same values. That alone took the coupled
coarse step from 21 s to 8.3 s and the SIMPLEC one from 31 s to 18.7 s;
the times of sections 7.1 and 7.5 are from before.

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

Where Foehncraft stands after part 6f: the iteration counts are in the
published range (8 to 11 times fewer iterations than SIMPLEC), the count
still grows with the mesh (150 on the fine mesh against 65), and the
gain in run time is a factor of 2.7 on the coarse step and 3.3 on the
fine one (section 7.7), where the published coupled solvers with a
block multigrid report 8 to 40 on meshes of 12 000 to 195 000 cells.
The gain grows with the mesh here too, but the comparison is not
closed. The turbulence step from the model time scale is our own
addition; the published solvers leave the turbulence equations on the
global step.

What closes the gap, in this order:

- **L44, the rest.** A multigrid with a coarse–fine splitting for
  the pressure (2 to 5 iterations per solve in the prototype where the
  aggregation multigrid needs 5 to 23; the pressure solves are a quarter
  of a coupled run); the multigrid of D179 as the preconditioner of the
  turbulence equations (their linear solves are a fifth of a coupled
  run, with 24 to 39 BiCGStab iterations per solve); a cheaper setup.
  By estimate each of them removes 5 to 15 % of the run time, none of
  them more.
- **A better Schur complement** for steps without the cap on the local
  step: the SIMPLE-type approximation needs 4 to 6 Krylov iterations
  with the cap 0.98 and 10 to 50 without it, with exact sub-solves. With
  it the count would no longer grow with the mesh. By estimate this is
  the largest lever: the fine step needs 150 outer iterations where the
  coarse one needs 65.
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
| D179 | Reduction multigrid (approximate ideal restriction) as the momentum sweep of the coupled preconditioner; `momentum_sweep`, `momentum_cycles`, `momentum_reuse`, `[solver.coupled.multigrid]` | stage 6 part 6f |
| L44 | The same multigrid for the turbulence equations, a cheaper setup, a splitting-based multigrid for the pressure, a Schur complement without the cap | the next solver steps |
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
- T. A. Manteuffel, S. Münzenmaier, J. Ruge, B. S. Southworth,
  *Nonsymmetric reduction-based algebraic multigrid*, SIAM J. Sci.
  Comput. 41, S242–S268 (2019)
- M. Ries, U. Trottenberg, G. Winter, *A note on MGR methods*, Linear
  Algebra Appl. 49, 1–26 (1983)
- J. W. Ruge, K. Stüben, *Algebraic multigrid*, in S. F. McCormick
  (ed.), Multigrid Methods, SIAM, 73–130 (1987)
- K. Stüben, *A review of algebraic multigrid*, J. Comput. Appl. Math.
  128, 281–309 (2001)
- W. A. Mulder, B. van Leer, *Experiments with implicit upwind methods
  for the Euler equations*, J. Comput. Phys. 59, 232–246 (1985)
- C. T. Kelley, D. E. Keyes, *Convergence analysis of pseudo-transient
  continuation*, SIAM J. Numer. Anal. 35, 508–523 (1998)
- Y. Mor-Yossef, Y. Levy, *The unconditionally positive-convergent
  implicit time integration scheme for two-equation turbulence models:
  revisited*, Computers & Fluids 38, 1984–1994 (2009)
- D. M. Driver, H. L. Seegmiller, *Features of a reattaching turbulent
  shear layer in divergent channel flow*, AIAA J. 23, 163–171 (1985)

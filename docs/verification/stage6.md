# Stage 6 – RANS turbulence and wall treatments

Stage 6 of the [stage plan](../requirements.md#7-stage-plan) brings the
RANS turbulence models and their wall treatments
([specification](../requirements/stage6-turbulence.md)). This page grows
with the pull requests of the stage. Part 6a gives every later model its
transport equation, its wall distance and its solver path
([below](#passive-scalars-and-the-wall-distance-part-6a)). The commands
are in [recipe 63](../recipes/scalars.md#passive-scalars), the pictures in
[gallery 63](../gallery/scalars.md#passive-scalars), the generated tables
of every case in the [verification manual](manual.md).

## Passive scalars and the wall distance (part 6a)

**The equation.** Every scalar of stage 6 and later (the passive scalars
of D51 now, k, ω, ε, ζ, f of the turbulence models, the temperature of
stage 7) is one instance of the generic transport equation of D116,

ρ ∂φ/∂t + ∇·(ṁ φ) − φ ∇·ṁ − ∇·(ρΓ ∇φ) = ρ (S_u + S_p φ),  S_p ≤ 0,

written once as single-source kernels (`foehncraft-kernels/src/scalar.rs`,
CPU loops and PTX entries `scalar_faces`, `scalar_rows`, `scalar_scale`,
gathers in the stored face order, so CPU and GPU give the same rows bit
for bit, D31). The face coefficients are those of the stage 2 convection
and diffusion operators: upwind with the linear-upwind or bounded central
deferred correction under the Venkatakrishnan limiter, over-relaxed
non-orthogonal diffusion with the least-squares gradient of the scalar.
Three parts are new:

- **Bounded form.** The term −φ_P Σ_f ṁ_f, the net mass outflow of the
  cell from the current fluxes, is implicit. The upwind operator then has
  zero row sums and stays an M-matrix while continuity is still
  converging, so a scalar between two bounds stays between them in the
  outer iterations of a steady run (Patankar 1980; the form of the
  industry standards).
- **Sources.** S_u explicit, S_p ≤ 0 implicit in the diagonal; a
  positive S_p is taken as zero, so the diagonal can only grow (the
  turbulence destruction terms use it in 6b).
- **Bounded BDF2** (D28 (1)). BDF2 is backward Euler from the old value
  φ* = −(a₁ φⁿ + a₂ φⁿ⁻¹)/a₀ with the step Δt/a₀. Where φ* leaves the
  `bounds` of the scalar, the cell takes the coefficients of BDF1 for this
  step: a select per cell, no branch beyond it, and the history counts the
  cells. Crank–Nicolson runs integrate the scalars with BDF2. At the start
  of a BDF2 step the scalar is extrapolated, φ ← φⁿ + ω (φⁿ − φⁿ⁻¹) with
  ω = Δt/Δt_prev, before its deferred correction is evaluated: with one
  outer iteration per step (PISO) the correction would otherwise lag a
  whole step and the scheme would be first order in time (observed 1.23
  before, 1.98 after).

**Boundary conditions.** The face value is φ_f = a φ_P + b per boundary
face: a fixed value (a = 0, b = φ_b) at inlets and walls, a fixed
diffusive flux (a = 1, b = q d_n/Γ, D116) at walls, slip walls, symmetry
planes and wedges, zero gradient at outlets. The tables
`[boundary.<name>.scalar.<scalar>]` are read before the flow conditions
are parsed, so the flow keys of every boundary kind stay checked as
before; `SCAL001` reports a table that does not fit its boundary.

**Exact wall distance** (D48, D108 part 1). A bounding volume hierarchy
over the triangles of the wall faces (the fans around the node average that the face geometry uses; binned
surface area heuristic with 16 bins, at most 4 triangles per leaf) answers
closest-point queries (Ericson 2005, section 5.1.5, for the point on a
triangle). The traversal visits the nearer child first and also every box
at the same distance as the best found, and a tie between faces goes to
the smallest global id, so the result does not depend on the threads, the
partition or the order of the faces (D93). The tests compare it with a
brute-force search on random surfaces and with min(y, 1 − y) in a channel
on hex, poly and tet cells (exact to 1e-12).

**Smooth limits** (spec 2.6). `foehncraft_core::math` has the hard and
smooth min, max and abs that the turbulence limiters of 6b to 6e use;
`Switches::HARD` gives the published models exactly, `Switches::SMOOTH`
the C¹ versions for the adjoint (D104).

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Laminar heated channel (`channel_heated.toml`), 512 / 2 048 / 8 192 cells | observed order of the temperature error | 1.8 – 2.2 | 2.01, 2.14 |
| | Nu on the finest mesh against 140/17 = 8.2353 | ± 0.5 % | 8.2398 (+0.055 %) |
| Manufactured solution 1 + cos(πx/2) sin(πy) in Poiseuille flow, Pe = 50 (`tests/scalars.rs`) | observed order on hex cells | 1.8 – 2.3 | 2.07, 2.15 |
| | on polyhedra | 1.6 – 2.6 | 2.50, 2.34 |
| | on triangular prisms | 1.6 – 2.6 | 2.25, 1.91 |
| Taylor–Green vortex carrying a scalar, Δt = 0.04 … 0.005 s | observed time order, BDF2 | 1.8 – 2.6 | 2.02, 1.98 |
| | observed time order, Crank–Nicolson | 1.8 – 2.6 | 2.02, 2.02 |
| Bounded front (`scalar_front_bdf2.toml`) | largest value with `bounds = [0, 1]` | ≤ 1 | 1 − 7e-11 (plain BDF2: 1.00033) |
| Restart of a transient run with a scalar after 8 of 16 steps | the scalar and its old level at the end | bit-identical | bit-identical |
| Wall distance in a channel, hex, poly, tet | against min(y, 1 − y) | 1e-12 | passed |

The heated channel is the open variant of the 9.2 case: the developed
profile at the inlet and a zero-gradient outlet, whose missing axial
diffusion changes the temperature by about Γ G / U = 2e-6 K in the last
cell. The periodic variant with the heat-flux jump of D109 follows with
the periodic boundaries (part 6p).

## Periodic boundaries (part 6p)

Work in progress (pull request of part 6p): translational periodic pairs
for the flow and the passive scalars, driven by a mass flow or a mean
pressure gradient, with the heat-flux jump of D109 (1) for the scalars.

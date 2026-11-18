# Stage 6 – RANS turbulence and wall treatments

Stage 6 of the [stage plan](../requirements.md#7-stage-plan) brings the
RANS turbulence models and their wall treatments
([specification](../requirements/stage6-turbulence.md)). This page grows
with the pull requests of the stage. Part 6a gives every later model its
transport equation, its wall distance and its solver path
([below](#passive-scalars-and-the-wall-distance-part-6a)), part 6p the
periodic boundaries that the channel cases need
([below](#periodic-boundaries-part-6p)), part 6b the k-ω SST model
([below](#k-ω-sst-part-6b)), part 6c the k-ε family with its wall
functions ([below](#k-ε-family-part-6c)), part 6d the k-ε-ζ-f model and
the multigrid of the turbulence equations
([below](#k-ε-ζ-f-and-the-scalar-multigrid-part-6d)), part 6d2 the
k-ε-ζ-f₀ model of Popovac (2023) with its wall treatment
([below](#k-ε-ζ-f₀-and-the-compound-treatment-of-2023-part-6d2)), part 6e the
coupled pressure–velocity algorithm and the pseudo-time steps that cut the
iterations of the steady turbulent cases
([below](#coupled-algorithm-and-pseudo-time-steps-part-6e)), part 6f the
multigrid of its momentum sweep
([below](#multigrid-for-the-momentum-sweep-part-6f)), part 6g the
multigrid of the turbulence equations and the larger steps with it
([below](#multigrid-for-the-turbulence-equations-part-6g)), part 6h the
pressure solve of the coupled preconditioner
([below](#the-pressure-solve-of-the-coupled-preconditioner-part-6h)). The commands are in
[recipe 63](../recipes/scalars.md#passive-scalars),
[recipe 67](../recipes/boundary-conditions.md#periodic),
[recipe 68](../recipes/turbulence.md#sst),
[recipe 71](../recipes/turbulence.md#k-epsilon),
[recipe 73](../recipes/turbulence.md#zeta-f),
[recipe 74](../recipes/turbulence.md#zeta-f0),
[recipe 78](../recipes/flow.md#coupled),
[recipe 80](../recipes/flow.md#momentum-multigrid),
[recipe 81](../recipes/flow.md#turbulence-multigrid) and
[recipe 83](../recipes/flow.md#pressure-hierarchy), the pictures in
[gallery 63](../gallery/scalars.md#passive-scalars),
[gallery 67](../gallery/boundary-conditions.md#periodic),
[gallery 68](../gallery/turbulence.md#sst),
[gallery 71](../gallery/turbulence.md#k-epsilon),
[gallery 73](../gallery/turbulence.md#zeta-f),
[gallery 74](../gallery/turbulence.md#zeta-f0),
[gallery 78](../gallery/flow.md#coupled),
[gallery 80](../gallery/flow.md#momentum-multigrid),
[gallery 81](../gallery/flow.md#turbulence-multigrid) and
[gallery 83](../gallery/flow.md#pressure-hierarchy), the generated
tables of every case in the [verification manual](manual.md).

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

**The pairs.** Two boundaries of `type = "periodic"` that name each other
form a translational pair (D29). `foehncraft_mesh::periodic::pair`
matches each face of one side with the face of the other whose centroid
lies at the translated position (a hash grid of the face centroids, within
1e-6 of the face size; the translation from the area-weighted centroids or
from `translation`) and checks that the two face areas are opposite. The
paired faces move to a block between the inner and the boundary faces,
side `a` and side `b` of each pair one after the other. Each side keeps
its own face with its own geometry, owned by its cell and listed only in
that cell: the neighbour of the face is the owner of its partner, seen at
its centroid plus the shift (−T on side `a`, +T on side `b`). The same
upper-triangular face and matrix layout then serves every kernel: the
face loops take the shifted centroid for the distance, the interpolation
weight and the least-squares stencil (`fv::neighbour_centroid`), the
transposed matrix pairs the two sides through the partner. A pair carries
the flux twice, once from each side, with the same values from the same
inputs, so mass is conserved across it to round-off. The pairing is the
last mesh step, after the renumbering; the mesh needs two cells across the
period and runs on one partition for now. Generated boxes give the nodes
of the two sides the same perturbation when the case pairs them.

**The driving.** The pressure is periodic, its level held in the reference
cell; a uniform body force β t̂ along the translation stands for the mean
pressure gradient. With `pressure_gradient = G` it is fixed, β = G. With
`mass_flow` it is adjusted after every pressure correction: a force
increment δβ moves the cell velocities by δβ D^C t̂ (the response of the
SIMPLEC rows) and the fluxes by the flux of that increment, so
δβ = (ṁ_target − ṁ)/K with K the flux of a unit increment through side
`b`; after the correction the mass flow through the pair is the target to
round-off. A uniform start velocity carries the target flow from the first
iteration. β is in the log, the report (`periodic`) and the restart file,
and it is part of the body force of the kinetic-energy balance.

**The scalars** (D109 (1)). A scalar crosses a pair unchanged
(`periodic = "periodic"`, the default) or with the jump of the heat-flux
mode, φ_N seen = φ_N ± Δ (minus on side `a`, plus on side `b`), in the
face coefficients (the jump goes into the source with the implicit
neighbour coefficient), the gradients and the limiters. Δ = Q/ṁ with Q
ρ times the integrals of the wall fluxes and the source and ṁ the mass
flow through the pair of the current fluxes, so the balance of the module
closes in every iteration. In a steady run the level is fixed after each
solve: the bulk value of the plane, the mass-flux mean of the face values
reconstructed with the gradients from both sides, is shifted to `bulk`. A
transient run keeps the level of its start. Fixed values on walls or
inlets are refused with the heat-flux mode (D109 (4)).

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Periodic Poiseuille flow (`poiseuille_periodic.toml`), mass flow 0.1 kg/s, 512 / 2 048 / 8 192 cells | observed order of the velocity error | 1.8 – 2.2 | 1.98, 1.94 (open channel 1.98, 1.89) |
| | L2 error against the open channel on the same meshes | factor 2 | 2.7e-3, 6.7e-4, 1.8e-4 (open 1.6e-3, 4.1e-4, 1.1e-4) |
| | driving force against 12 μU/h² = 1.2 Pa/m | | −0.62 %, −0.06 %, +0.09 % |
| | the same driven by 1.2 Pa/m on uniform meshes: mass flow against 0.1 kg/s | order 1.8 – 2.2 | +0.78 %, +0.20 %, +0.05 % |
| Periodic heated channel (`channel_heated_periodic.toml`), heat-flux jump | observed order of the temperature error | 1.8 – 2.2 | 2.40, 2.02 |
| | Nu on the finest mesh against 140/17 = 8.2353 | ± 0.5 % | 8.2460 (+0.13 %; open channel +0.055 %) |
| | jump Δ = Q/ṁ and the bulk value at the plane | 1e-12 | 4e-3 K, −17/70 K |
| Short periodic heated pipe (`pipe_heated_periodic.toml`), one diameter, wedge of 0.5°, 64 / 256 / 1 024 cells | observed order of the velocity error | 1.8 – 2.2 | 2.02, 1.98 |
| | observed order of the temperature error | 1.8 – 2.6 | 2.30, 2.28 |
| | Nu on the finest mesh against 48/11 = 4.3636 | ± 0.5 % | 4.3691 (+0.13 %) |
| | driving force against 8 μU/R² = 3.2 Pa/m | | −1.06 %, −0.28 %, −0.09 % |
| Pairing of perturbed hex, poly and tet boxes in x and z | every face matched, cell geometry | round-off | passed |
| Transient periodic run with the heat-flux mode, restart after 6 of 12 steps | fields, scalar and β at the end | bit-identical | bit-identical |
| Every case without periodic boundaries | results | bit-identical | unchanged |

The temperature of the pipe runs without a limiter: the Venkatakrishnan
limiter clips the smooth minimum on the axis in the prisms next to it, and
on the perturbed meshes the observed order falls to 1.65 and 1.42 (on a
fourth mesh); without it the order is 2.3 and 2.2, and on unperturbed
meshes 2.0 either way. The first pair of the channel sweep converges
faster than second order (2.40); the finest pair, which the check uses,
gives 2.02.

## k-ω SST (part 6b)

**The model** (D130). SST 2003 of Menter, Kuntz & Langtry: k and ω on the
generic scalar path of 6a with the F₁ blend of the constants, the cross
diffusion (1 − F₁) 2σ_ω2/ω ∇k·∇ω split by its sign (the positive part
explicit, the negative part in the diagonal), the eddy viscosity
ν_t = a₁k/max(a₁ω, S F₂), the production P̃ = min(ν_t S², 10 β* k ω)
(or Kato–Launder, Durbin, none), the sustaining terms of Spalart & Rumsey
and a cap of ν_t/ν. The destruction is linearised point-implicitly
(β*ω in the diagonal of k, 2βω with the Newton option for ω), so both
stay positive in a converging run; the floors (10⁻¹⁰ of the largest start or
inlet k, and the ω that gives ν_t/ν = 100 with it) are counted (TURB002 when they act in the final tenth of
the solves). Without walls the wall distance is infinite and F₁ = 0
(TURB005). The gradients of k and ω use the Venkatakrishnan limiter with
a threshold relative to the range of the variable, ε = K (φ_max − φ_min),
K = 0 by default: the dimensional threshold of the momentum equations is
far above the differences of k on a flat plate and would switch the
limiter off.

**The walls** (D134). `low-re`: k = 0 and ω_w = 60ν/(β₁ d₁²) as Dirichlet
values on the wall face, d₁ = 2 y_P the height of the wall cell, the value
of the TMR codes (the first-cell value 6ν/(β₁ y_P²) of the draft converged
at first order to +1.1 %; specification 4.3). `automatic`:
u_τ = (u_vis⁴ + u_log⁴)^¼ with the log law floored at y⁺ = 11.06, the
first-cell ω_P = (ω_vis² + ω_log²)^½ imposed with the weight m = 1 − b,
b = u_vis⁴/u_τ⁴, through a fractional fix of the row,
(1 − m)·row + m a_P (φ − ω_P), and the production of the first cell
blended likewise, (1 − m) P̃ + m u_τ³/(κ y_P). Rough walls take the wall
values of Aupoix (2015). The first-cell height y_P is floored at
10⁻³ V/A of the wall cell.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| TMR flat plate (`tmr_plate_sst.toml`), `low-re`, 768 / 3 072 / 12 288 / 49 152 cells, first-cell y⁺ 0.09 on the finest | C_f(0.970084) on the finest mesh against the TMR value 2.6907e-3 | ±1 % + U_num | 2.7007e-3 (+0.37 %) |
| | observed order and GCI_fine of C_f (three finest meshes) | p 1.5 – 2.5, GCI ≤ 1 % | p = 0.43, GCI 1.87 %, extrapolated 2.741e-3: **finding** |
| | C_D (one side) against the TMR value 2.8488e-3 | ±1 % + U_num | 2.8606e-3 (+0.41 %); p = 0.74, GCI 0.84 % |
| | C_f(Re_x) against White (2006), 10⁶ ≤ Re_x ≤ 10⁷ | ±5 % | −8.7 % to −5.3 % (the TMR value of SST itself is −6.4 % at x = 0.97): **finding**, the model's difference from the correlation |
| | u⁺(y⁺) at x = 0.97 against Spalding's law, 1 ≤ y⁺ ≤ 300 | ±3 % | 2.3 % largest |
| The plate with `automatic`, 3 072 cells, first-cell y⁺ 1.7 / 4.6 / 9.3 / 27 / 92 | C_f(0.97) against the TMR value | ±2 % | +2.3 %, +15.6 %, +14.8 %, −0.2 %, −0.7 %: **finding** in the buffer layer, TURB001 warns there |
| Channel at Re_τ = 395 (`channel_sst.toml`), periodic module, mass flow of the DNS bulk velocity, `low-re`, 96 cells across, first-cell y⁺ 0.28 | U_b⁺ against the DNS of Moser, Kim & Mansour (17.55) | ±3 % | 17.43 (−0.70 %; Dean's correlation 17.24: +1.1 %), Re_τ = 397.8 |
| | U_c⁺ against the DNS (about 20.1) | ±4 % | 19.64 (about −2.4 %) |
| | observed order and GCI_fine of U_b⁺ on the three meshes of 96, 192 and 384 cells across (U_b⁺ 17.427, 17.381, 17.346) | p 1.8 – 2.2, GCI ≤ 0.5 % | p = 0.40, extrapolated 17.237, GCI_fine 0.78 %: **finding** (the fourth mesh, 768 cells across, did not reach its tolerance in 30 000 iterations and is left out) |
| Homogeneous decay in a periodic box, k₀ = 1, ω₀ = 10, t = 0.5, 10 / 20 / 40 steps | BDF2 time order of k and ω against the exact solution (F₁ = 0) | 1.8 – 2.2 | k 1.92, 1.96 and ω 2.07, 2.03 (two pairs of meshes) |
| Steady plate, restart after half the iterations | fields at the end | bit-identical | bit-identical |
| Plate on 1 and 4 threads | fields after 8 iterations | bit-identical | bit-identical |
| Poisson wall distance on a plane wall | d against the exact distance | round-off | passed |

**Findings.** (0) U_b⁺ of the channel approaches 17.24 (the extrapolated value, Dean's correlation) at the low order 0.40 with the first-cell y⁺ falling from 0.28 to 0.07, like the plate (0.43): the wall value of ω, which grows as 1/d₁² with the mesh, converges slowly in both; a sweep of the wall value is an item for 6c. (1) The C_f of the plate converges monotonically but at a
low observed order on these meshes (0.43); the finest value is within
0.4 % of the TMR, which reports its values after five grids, and the
differences between the three finest meshes are below 0.7 %. A fifth mesh
(512 × 384) or the TMR's own family will show whether the order rises;
the band stays as it is. (2) The automatic treatment overpredicts C_f by
up to 16 % with the first cell in the buffer layer: the blended u_τ of the
original formulation and the log-layer production meet badly there. The
check TURB001 warns for y⁺ between 3 and 20; an improved blending is an
item for the wall-function work of 6c. (3) White's correlation is a fit
to experiments that lies 5 to 9 % above the SST result along the whole
plate, as does the TMR's SST result: the band of the specification for
this row was set for the model, and SST does not meet it.

**Not in this part.** The manufactured solution of SST (moved to 6g with
the other model MMS), the D166 features as named variables of the
expression language (with their first consumer, D104), turbulence values
of the backflow at outlets (the cell's values are used; BC028 says so),
and curvature correction and QCR (D170, part 6h).

## k-ε family (part 6c)

**The models** (D131). The standard k-ε of Launder & Spalding (C_μ 0.09,
C_ε1 1.44, C_ε2 1.92, σ_k 1.0, σ_ε 1.3) and the realizable k-ε of Shih
et al. (1995): variable C_μ = 1/(A₀ + A_s U* k/ε) with A₀ = 4.0, the
destruction C₂ ε²/(k + √(νε)), C₁ = max(0.43, η/(η + 5)), σ_ε 1.2. One
kernel (`turb::keps`) for both, generic over the number type, so the
dual-number Jacobian sees the closures; the arccos of A_s is a helper of
`foehncraft-core::math` with its derivative. The standard model takes the
Durbin bound by default (inactive in equilibrium, it keeps the impulsive
start of near-wall cells bounded); Kato–Launder, the clip
min(P_k, 10ε) and `none` are options. The ε destruction is linearised
with the Picard form by default: the Newton form removes at most half
of ε per iteration and let k collapse in the decaying free stream of the
plate. Inlets and initial values take ε directly or ω (converted with
ε = C_μ k ω); a restart file of the other family is converted the same
way.

**The walls** (D134). All treatments take u* = C_μ^¼ k_P^½ and
y* = u* y_P/ν, with no iteration. `scalable`: the log law at
ỹ* = max(y*, 11.06), ε_P = u*³/(κ ỹ), P_k,P = (τ_w/ρ) u*/(κ ỹ).
`standard`: the log law from y* = 11.06, the viscous sublayer below it.
`compound`: the blend of Popovac & Hanjalić with Γ = 0.01 y*⁴/(1 + 5y*)
for u⁺, ε_P and P_k,P. ε_P is fixed in the near-wall cell (weight 1) and
written into the field before each solve; k and ε have zero gradient at
the wall. Rough walls shift the log law by ΔB of Cebeci & Bradshaw with
C_s = 0.5.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| TMR flat plate (`tmr_plate_keps.toml`), standard k-ε, `scalable`, 3 072 cells, sustaining terms at the inlet values, first-cell y⁺ 1 / 2 / 5 / 10 / 30 / 60 / 100 at x = 0.97 | C_f(0.97) against White, 2.875e-3 | ±5 % | −2.6, −1.7, −1.0, −0.3, −1.6, −1.5, −0.8 % |
| | y⁺ insensitivity against the y⁺ ≈ 1 result | ±5 % (scalable from y⁺ 30) | +1.1 to +2.4 % (every y⁺ within 2.4 %) |
| The plate, standard k-ε, `compound` | C_f(0.97) for y⁺ 30 / 60 / 100 | ±5 % | −2.3, −3.6, −3.4 % |
| | C_f(0.97) for y⁺ 1.3 / 2.6 / 6 / 11 | not claimed | +58, +52, +40, +21 %: **finding**, the high-Re model has no damping in the buffer layer; TURB001 below y⁺ 30 |
| The plate, realizable k-ε, `scalable` and `compound`, y⁺ 30 / 60 / 100 | C_f(0.97) against White | ±5 % | scalable −6.8, −5.2, −3.3 %; compound −7.4, −7.1, −5.9 %: **finding** |
| Channel at Re_τ = 395 (`channel_keps.toml`), standard k-ε, `scalable`, 12 uniform cells, first-cell y⁺ 32 | U_b⁺ against the DNS (17.55) | ±5 % | 18.31 (+4.3 %) |
| The channel on 4 to 40 cells (y⁺ 93 to 9), `scalable` | U_b⁺ | ±5 % from y⁺ 30 | +6.0 % (y⁺ 93), +4.3 to +4.9 % (y⁺ 47 to 9) |
| The channel, `compound`, y⁺ 93 to 11; graded mesh at y⁺ 0.4 | U_b⁺ | ±3 % at y⁺ ≈ 1 | +5.9 % (93) to +1.6 % (24), −4.4 % (16), −12.9 % (11); −32 % at y⁺ 0.4 (τ_w 2.2 times too high): **finding**, resolved walls go to ζ-f and EB (owner decision 2026-11-17) |
| The channel, realizable, `scalable`, y⁺ 93 to 9 | U_b⁺ | ±5 % | +6.0 to +10.7 %: **finding** |
| Rough pipe (`pipe_rough_keps.toml`), Re = 10⁶, k_s/D = 10⁻³ | f against Colebrook (0.01994) | ±5 % | 0.02081 (+4.4 %) |
| Moody study (`examples/studies/moody_study.py`), Re 3·10⁴ to 10⁷, k_s/D 0 to 10⁻² | f against Colebrook | ±5 % | smooth −5.0 to −1.7 %; fully rough +8 to +11 % (C_s = 0.4: +1 to +3 %); transitional −17 to +8 %: **finding** |
| Homogeneous decay in a box, k₀ = ε₀ = 1, t = 0.5, 10 / 20 / 40 steps | BDF2 time order of k and ε against the exact solution | 1.8 – 2.2 | asserted in `tests/keps.rs` |
| Sustaining terms in a box at rest, both models | k and ε constant | round-off | asserted |
| Restart of an SST run with the k-ε model | ε = C_μ k ω on the way in, ε in the new file | – | asserted |
| Plate with `compound` on 1 and 4 threads | residuals, wall shear after 8 iterations | bit-identical | asserted |
| Kernels | realizable C_μ in simple shear, the realizability bound ν_t ≤ k/(3λ_max), the sustaining terms, the source signs, dual-number derivatives; the wall treatments in the log layer and the sublayer, ΔB continuous at k_s⁺ 2.25 and 90 | – | kernel tests of `turb::keps` and `turb::wall` |

**Findings.** (1) The compound treatment with the high-Re k-ε models
overpredicts the wall shear below y⁺ 30 (up to 2.2 times at y⁺ ≈ 1): the
models have no damping in the buffer layer. The owner chose to keep it
as a finding and leave resolved walls to ζ-f (6d) and EB (6e); TURB001
warns. (2) The realizable model gives 3 to 7 % less wall friction than
the correlations on the plate and 6 to 7 % more U_b⁺ in the channel, in
both cases in the same direction; the log-layer consistency of its
constants holds (κ² = σ_ε(C₂√C_μ − C₁) = 0.168), so the near-wall cell,
where its C_μ responds to the poorly resolved velocity gradient, is the
first suspect; open for 6g. (3) With the default C_s = 0.5 the fully
rough friction factor lies 8 to 11 % above Colebrook; C_s = 0.4 puts it
within 3 %. In the transitional range the Nikuradse-type ΔB cannot follow
Colebrook's monotonic curve (up to 17 % below). (4) With the TMR inlet
values and no sustaining terms the free stream of the k-ε family decays
(the realizable k to zero): the plate runs with the sustaining terms.

**Not in this part.** The improved blending of the SST `automatic`
treatment and a sweep of the SST wall value of ω, both left open by 6b,
move to 6g with the model comparison; the GCI of the k-ε plate and
channel also follows there.

## k-ε-ζ-f and the scalar multigrid (part 6d)

**The model** (D132). Hanjalić, Popovac & Hadžiabdić (2004):
ν_t = C_μ ζ k T, C_μ = 0.22, C_ε1 = 1.4(1 + 0.012/ζ), C_ε2 = 1.9,
C₁ = 1.4, C₂′ = 0.65, σ_k 1, σ_ε 1.3, σ_ζ 1.2, C_T 6, C_L 0.36, C_η 85,
the realizability coefficient 0.6 (constants checked against the page of
the turbulence modelling benchmark working group, TMBWG). The bounds of T
and L take the strain magnitude S = (2 S_ij S_ij)^½, as the authors'
implementation does (6d amendment of specification 2.1, which read
(S_ij S_ij)^½: that form gave U_b⁺ +2.9 % in the channel against +1.9 %).
C_ε1 P is evaluated as 1.4(P + 0.012 P/ζ) with P/ζ = C_μ k T S², finite
as ζ → 0. One kernel (`turb::zeta_f`) gives T, L², ν_t, the source split
of k, ε, ζ and the right-hand side of f, generic over the number type
(dual-number test). Solve order per outer iteration: k and ε; the cell
values again with the new k and ε, then f (CG); the cell values again,
then ζ with the new f, bounded to [10⁻¹⁰, 2]; ν_t at the end. A transport
solve may lower k, ε or ζ by at most a factor 10 (`ZETA_F_DROP`; 6d
amendment): an unconverged BiCGStab iterate of the loose steady tolerance
can go negative where the solution is small, and the absolute floor
alone put ε near zero, which blew up T, L and P/ε in the source of f
(the plate at y⁺ 0.5 diverged in iteration 36).

**The walls** (D134). `low-re`: k = ζ = 0 (Dirichlet), ε_w = 2νk₁/d₁²
and f_w = −2νζ₁/d₁² on the wall faces. `compound` (default): k and ζ by
zero gradient (dk/dy = dζ/dy = 0 at the wall), ε and P_k of the first
cell fixed by the blend of Popovac & Hanjalić (2007) with C_μ = 0.09 in
the log-law values, f_w = e^(−Γ)(−2νζ₁/d₁²). Roughness with `low-re` is
refused.

**The multigrid** (D139). The backends hold hierarchies in slots
(0 pressure, 1 f, 2 k, 3 ε or ω, 4 ζ). The symmetric positive definite f
equation runs CG with its own hierarchy by default; the transport
equations take one with `solver.turbulence.preconditioner = "amg"`,
built from the symmetric part ½(A + Aᵀ) and preconditioning BiCGStab.
All rebuild on the schedule of the pressure (`solver.amg_rebuild`), and
the matrices of their last setups go into the restart file
(`turb_amg_matrix_<slot>`), so a continued run builds the same
hierarchies. A rebuild in every solve (the first version) cost 21 s of
the 56 s of the channel: on a small mesh the coarsest level holds all
cells and its dense inverse dominates.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Channel at Re_τ = 395 (`channel_zeta_f.toml`), 96 cells graded 60, first-cell y⁺ 0.27, `compound` | U_b⁺ against the DNS (17.55) | ±3 % | 17.89 (+1.9 %), 3 700 iterations |
| The channel, `low-re` | U_b⁺ | ±3 % | 17.79 (+1.4 %) |
| TMR flat plate (`tmr_plate_zeta_f.toml`), `compound`, 64 × 48 cells, first-cell y⁺ 1, turbulent start | C_f(0.97) against White, 2.875e-3 | ±3 % | 3.023e-3 (+5.1 %): **finding** |
| The plate on 64 × 96 cells (y⁺ 1) | C_f(0.97) | ±3 % | +4.7 %: mesh-independent, **finding** |
| The plate with an inlet of Tu 0.8 % (k = 10⁻⁴, ν_t/ν = 10), 48 cells | C_f(0.97) | – | +5.0 % |
| The plate from the TMR free-stream values (ν_t/ν = 0.009) | transition | – | laminar over the whole plate, C_f = 3.0e-4: **finding** |
| The plate at y⁺ 0.5 | stability | – | diverged without the relative bound; converges with it |
| Run time of the channel (192 cells, 3 000 iterations) | multigrid rebuild per solve against the pressure schedule | – | 56 s against 26 s (SST: 16 s) |
| Homogeneous decay, k₀ = ε₀ = 1, ζ₀ = 0.3, t = 0.5, 10 / 20 / 40 steps | BDF2 order of k, ε and of the return of ζ to 2/3 | 1.8 – 2.2 | asserted in `tests/zeta_f.rs` |
| Sustaining terms in a box at rest | k, ε constant, ζ = 2/3, f = 0 | 10⁻¹⁰ | asserted |
| Transient restart with the multigrid of all four equations | fields at the end | bit-identical | asserted |
| The plate on 1 and 4 threads (20 480 cells, multigrid of all equations) | residuals, wall shear after 8 iterations | bit-identical | asserted |
| A second hierarchy slot | the V-cycle of the pressure slot | bit-identical | asserted (`tests/stage3.rs`) |
| Kernels | T, L and ν_t in simple shear, the decay and sustaining limits, the source signs, the low-Re wall values, dual-number derivatives | – | kernel tests of `turb::zeta_f` |

**Findings.** (1) On the TMR plate the model gives C_f 5 % above White's
correlation, independent of the mesh (48 or 96 cells across), outside
the ±3 % band; in the channel it lies within 2 %. (2) From the TMR free
stream (ν_t/ν = 0.009) the boundary layer stays laminar over the whole
plate: the length scale of f in the weak free stream,
C_L C_η (ν³/ε)^¼ ≈ 7 mm, exceeds the laminar layer, and f drives ζ, and
with it ν_t, to zero across it; the case starts from a turbulent field.
(3) In the laminar part near the leading edge the published sink −ε of
the k equation has no positive steady solution once T sits at its
Kolmogorov bound (ε then no longer falls with k): the floor holds k in
about 120 cells and its residual stays at 7·10⁻⁶ while every other field
converges to round-off (TURB002). The ζ-f₀ variant of Popovac (2023),
with f₀ = 0 at the wall, follows as the next part for robustness in y⁺.

## k-ε-ζ-f₀ and the compound treatment of 2023 (part 6d2)

**The model** (D132, specification 3.3.1). Popovac (2023): f₀ = f +
(N_ζ − 1)ζ/T in place of f, N_ζ = 2, so f₀,w = 0;
ν_t = (C_μ/ζ_eq) ζ k T with C_μ = 0.09, ζ_eq = 0.41;
C_ε1 = 1.4(1 + 0.012/√ζ); the other constants of 2004. The same kernel
as ζ-f (`turb::zeta_f`, `FORM_F0`). The author's code is not available
(maintainer, 2026-11-17), so the default was chosen on the validation:

| Variant | Channel U_b⁺, y⁺ 1 | Plate C_f, y⁺ 1 | Step x_r/H |
|---|---|---|---|
| As printed (sign of (N_ζ − 1) in eq. 7 corrected, k sink −ε) | −11 % | +21 % (TMR), floor holds k | 7.34 (Gmsh mesh) |
| C_ε1 with 1/ζ (2004) | +8.5 %, not converged | – | – |
| k sink k/T with the realizability bound | +30 % (relaminarizes) | – | – |
| k sink k/T_K, a = 0.6 | +1.3 % | −4.0 % (White) | 7.67 |
| **Default: k sink k/T_K, a = 0.6 √2** | **+1.4 %** | **−2.5 % (White)** | **6.38** |

**The wall treatment** (D134, `compound-2023`). u* = C_μ^¼ k_P^½,
y* = u* y_P/ν; u⁺ = [y*^(−n) + u⁺_log^(−n)]^(−1/n), n = 10/3,
u⁺_log = ln(8.9 max(y*, 1))/(κψ); ψ = 1 − C_U max(y_P − y_v/2, 0)/(U_t κ u*)
with y_v = 11ν/u* and C_U = t·[(U·∇)U + ∇p/ρ] of the wall cell (steady
part), bounded to [0.5, 2] (the paper gives no bounds); ν_w = ν y*/u⁺.
Γ = du⁺/dy* from the blend, ε_P = Γ 2νk_P/y_P² + (1 − Γ) u*³/(κ y_P)
and P_k,P = u*² (τ_w/ρν) Γ, fixed in the cell; k and ζ by zero flux.
It also serves ζ-f (option).

**The mesher.** `mesh.generate.remove` takes the cells whose centroid lies
in a box out of a hexahedral, tetrahedral or mixed mesh and names the
exposed faces; the step case uses it.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Channel at Re_τ = 395 (`channel_zeta_f0.toml`), 48 cells, first-cell y⁺ 1 | U_b⁺ against the DNS (17.55) | ±3 % | 17.79 (+1.4 %), converged |
| The channel at y⁺ 0.3 / 10 / 30 | U_b⁺ | – | +2.7 / −4.2 / −1.5 % (ζ-f, 2007 treatment: +1.9 / −7.1 / +4.8 %) |
| TMR plate (`tmr_plate_zeta_f0.toml`) from the TMR free stream, y⁺ 1 | C_f(0.97) against White, 2.875e-3; convergence | ±3 % | 2.802e-3 (−2.5 %), converged to 10⁻⁷, no floor |
| The plate at y⁺ 10 / 30 | C_f(0.97) against White | ±3 % | +0.9 / −2.0 %, converged (ζ-f: not converged) |
| Backward-facing step (`backstep_zeta_f0.toml`), 8 200 cells | x_r/H against 6.26 | ±10 % | 6.38 (+2 %) |
| The step on 32 800 cells | x_r/H | ±10 % | 6.35 (+1.5 %), within 0.5 % of the coarse mesh |
| The step with ζ-f (2004, `compound`) | x_r/H | ±10 % | 7.82 (+25 %; 6.80 with `a = 0.849`) |
| The step with SST (`automatic`) on 8 200 / 32 800 cells | x_r/H | ±10 % | 6.81 / 6.80 (+9 %); the 5.79 reported first came from 4 000 iterations at α_u = 0.5, before the bubble had settled |
| The step with k-ε (`scalable`), 32 800 cells | x_r/H | – | 5.36 (−14 %) |
| The plate on 1 and 4 threads (20 480 cells) | residuals, wall shear after 8 iterations | bit-identical | asserted |
| Kernels | f₀ form against the 2004 source, the limits of the blend (y⁺ 30 to 300), Γ as the derivative of the blend, ψ, dual numbers | – | kernel tests of `turb::zeta_f` and `turb::wall` |
| Mesher | step volume, exposed area, unused nodes, taken names | exact | `removed_box_builds_a_step` |

**Findings.** (1) The realizability bound of T with a = 0.6 and
S = (2 S_ij S_ij)^½ is active in the separated shear layer of the step
and halves ν_t there; the bubble then grows to 7.4 to 7.7 H with both
ζ-f models. (2) The published k sink −ε has no positive steady k in the
laminar part of a boundary layer (the 6d finding); k/T_K removes it.
(3) C_p, C_f and the profiles of the step against the experiment and a
GCI follow with the stage 6 validation. (4) Iterations of the step
(measured after 6d2, fields every 50 to 250 iterations): with the
momentum under-relaxation 0.5 throughout, x_r is within 1 % of its final
value after 950 (ζ-f₀), 1 250 (k-ε) and 1 500 (SST) iterations on 8 200
cells, and after about 2 900 (ζ-f₀) and more than 4 000 (SST) on 32 800
cells. With 0.5 for the first 100 iterations and 0.8 after them it takes
550 to 600 and about 1 750; 0.7 from the start works as well, 0.9 stalls
the residual at 10⁻⁴, and a turbulence relaxation of 0.95 oscillates. The
count grows about three times when the mesh spacing halves: the
relaxation sets a local pseudo time step of a few cell transit times, and
the bubble needs a fixed time to form. The scaled residual is 5·10⁻⁵ to
2·10⁻⁴ when x_r has settled; a residual tolerance does not show when the
step result is done.

## Coupled algorithm and pseudo-time steps (part 6e)

Specification: [requirements step 22](../requirements/step22-coupled-solver.md)
(D176 to D178), which also holds the measurements behind the choices.

**The algorithm** (D176, `foehncraft-flow::coupled`). One outer
iteration solves

    | A  G | |δu|   |r_u|
    | D  L | |δp| = |r_c|

for the increments, with the momentum rows A, the residuals and the
Rhie–Chow fluxes of SIMPLEC; then u += δu, p += δp, the fluxes follow
from the same increments, and one pressure correction removes the
continuity defect that the loose linear solve (0.1) leaves. The operator
is applied with the operations of the device (no block matrix); flexible
GMRES, right-preconditioned with a SIMPLE-type sweep: BiCGStab on the
momentum rows to 10⁻², the pressure matrix with ρ(D̄ + D̂̄) by CG and the
multigrid of the pressure, the velocity correction with
D̂ = min(D^C, V/((1 − α) a_P/α + ρV/(αΔτ))) (kernel `schur_diagonal`).

**The pseudo-time steps** (D177, D178). Momentum: ρV/Δτ (u − uᵏ) on the
rows through the time term of the transient solver and the
time-step-independent Rhie–Chow flux (D27) with the current iterate as
the old level; Δτ = `factor` × ΣV/ΣA_b / U, one value for all cells;
α_u = 0.98 as the cap on the local Courant number. Turbulence:
ρV max(1/Δτ, 1/(c_T T)) (φ − φᵏ) on the scalar rows (kernel
`pseudo_coefficient`), c_T = 0.5, relaxation 0.98; the residual scale
does not contain the term.

### Results

Backward-facing step of part 6d2, uniform start, iterations until x_r
stays within 1 % of its final value (a result every 5 iterations for the
coupled runs, every 25 for SIMPLEC, every 50 on the fine mesh; the fine
runs went on to 400 and 3 000 iterations for the final value):

| Test | Check | SIMPLEC (0.5, then 0.8 from 100) | Coupled |
|---|---|---|---|
| Step, k-ε-ζ-f₀, 8 200 cells (`backstep_zeta_f0_coupled.toml`) | iterations to 1 %; final x_r/H | 525; 6.384 | 65; 6.385 |
| Step, k-ε (`scalable`), 8 200 cells | the same | 600; 5.287 | 75; 5.287 |
| Step, k-ω SST, 8 200 cells | the same | 675; 6.82 | 80; 6.81 |
| Step, k-ε-ζ-f₀, 32 800 cells | the same | 1 600; 6.352 | 175; 6.352 |
| Channel at Re_τ = 395, k-ε-ζ-f₀ (`channel_zeta_f0.toml`) | iterations to 10⁻⁶; U_b⁺ | about 2 100; 17.79 | 323; 17.79 |
| Cavity at Re = 1000, 128² (`cavity_re1000_coupled.toml`) | iterations to 10⁻⁸; deviation from Ghia (u, v) | 1 397 (α = 0.9); 0.0040, 0.0136 | 452; 0.0040, 0.0136 |
| The same cavity, kinetic-energy ledger (D88) | numerical dissipation over lid power; its pressure part (W) | 0.637 %; 3.996223 × 10⁻⁶ | 0.637 %; 3.996224 × 10⁻⁶ |

| Test | Check | Band | Result |
|---|---|---|---|
| Cavity 16², Re = 100, converged to 10⁻¹¹ | lid force against SIMPLEC, for the automatic step, Δτ = 0.02 s, `factor = 20` with α = 0.9, and relaxation only | 10⁻⁸ | asserted (`coupled_result_is_the_simplec_result`) |
| Cavity 48², Re = 400 | iterations: SIMPLEC at α = 0.9; coupled with `factor` 1, 4, 20 | factor 4 at most 0.7 of SIMPLEC, 20 fewer than 4 | 255; 292, 130, 84 |
| Cavity 16², 10 iterations, lid speed 1 and 2 m/s | the automatic step (L = 0.25 m over the lid speed; the symmetry planes of the one cell layer must not enter); `factor`, `step`, `momentum = false` | 0.25 s and 0.125 s to 10⁻¹² | asserted (`automatic_step_is_the_length_over_the_prescribed_speed`) |
| Cavity 16², 24 iterations, restart at 10 | velocity, pressure, fluxes, multigrid matrix of the continued run | bit-identical | asserted |
| Cavity 150², 8 iterations, 1 and 4 threads | residuals, Krylov count, boundary forces | bit-identical | asserted |
| Input | `[time]` with `coupled`, relaxation 1, ranges of `[solver.coupled]` and `[solver.pseudo_time]`, the defaults | refused with CASE004 / CASE005 | asserted |
| Kernels | pseudo-time term of the scalar row and its residual scale, `pseudo_coefficient`, `schur_diagonal` | – | kernel tests |

Run time of the step with k-ε-ζ-f₀ until x_r has settled, one thread:
21 s against 31 s on the coarse mesh (65 against 525 iterations), about
235 s against 327 s on the fine one (175 against 1 600). A coupled
iteration costs five to seven SIMPLEC iterations (10 and 14 Krylov
iterations, each with a momentum sweep of 10 iterations and a pressure
solve), so the run time gains a factor of 1.4 to 1.5 where the count
gains 8 to 9.

**Findings.** (1) The coupled solve alone, with relaxation as before,
halves the count (300 on the coarse step, 700 to 800 on the fine one):
the larger part of the gain is the pseudo-time step. (2) On the step the
count does not depend on Δτ between 5 and 29 H/U, nor on having a global
step at all once the cap is 0.98; the global step makes the linear
solves cheaper (8 instead of 18 Krylov iterations per outer iteration).
(3) The cap on the local step sets the count on fine meshes and on the
cavity (Δτ given; 0.95: 1 012, 0.98: 423, 0.999: 166 iterations); a tighter
coupled solve does not (167). At 0.999 the momentum sweeps cost four
times the run time, so 0.98 is the default. (4) The momentum sweep must
reach 10⁻² to 10⁻³; truncated sweeps make the Krylov count and x_r
erratic. On the momentum matrix of the step at Δτ = 5 H/U, BiCGStab
needs 113 to 137 iterations with a diagonal preconditioner, 21 to 24
with symmetric Gauss–Seidel, 90 to 194 with the aggregation multigrid of
the symmetric part and 1 with a multigrid of approximate ideal
restriction (L44). (5) Turbulence: a global step is limited to 0.5 H/U;
with the time-scale rule c_T = 0.25, 0.5, 1 gives 75, 55, 45 to 50
iterations (cap 0.999); relaxation 1 lets k-ε diverge at the start, 0.9
makes the turbulence lag the momentum equations (x_r overshoots to 11 H
and has not settled after 150 iterations). (6) With SIMPLEC the
time-scale rule alone (`solver.pseudo_time.turbulence = 0.5`) gives 230
iterations at α_u = 0.8; it is off by default there. (7) On small
laminar cavities the coupled algorithm is not faster than SIMPLEC at
α = 0.9: the automatic step (0.25 L/U) is a careful one, and larger
steps cost more per iteration than they save (cavity at Re = 1000 on
64²: 584 iterations in 8 s with SIMPLEC; 625 in 13 s, 276 in 24 s and
178 in 43 s with `factor` 1, 4 and 20). The momentum sweep is the limit
in every one of these findings. (8) A run repeats bit by bit for the
same input, but its count reacts to the last digit of an input: the
cavity needs 423 iterations with Δτ = 0.25 s given and 452 with the
automatic step, which differs from 0.25 s in the 17th digit; the fine
step needed 145 iterations with a slightly shorter step (a velocity
scale that counted the symmetry planes) and needs 175 at
7.3 H/U. The
coupled solve stops at 0.1, its momentum sweep is cut off after 10
iterations and the multigrid cycle runs in single precision, so two runs
separate in the fifth digit within one outer iteration. The converged
result does not change; the counts of this section are good to about
10 %. The histories of the figure of gallery 78 come from intermediate
results every 5 iterations (every 25 with SIMPLEC). An intermediate
result leaves a run unchanged (stage 5, control files; test
`intermediate_results_and_balances_do_not_change_the_run`), so they are
the histories of the runs without them: the fine step gives x_r/H =
6.27, 6.34, 6.37 after 170, 180, 190 iterations either way. Before that
correction a result file shifted the path within this scatter (the SST
count read 85 instead of 80).

## Multigrid for the momentum sweep (part 6f)

Specification: [requirements step 22](../requirements/step22-coupled-solver.md)
sections 3.1 and 7.7 (D179).

**The method** (`foehncraft-fv::amg::reduction`). The momentum sweep of
the coupled preconditioner is one V-cycle of a reduction multigrid per
velocity component instead of BiCGStab. Setup on the host, per level:
strong dependencies (|a_ij| ≥ 0.3 max_{k≠i} |a_ik|); the first pass of
the Ruge–Stüben splitting into coarse points C and fine points F; the
restriction R = [Z I] with Z = −A_CF (2 D_F⁻¹ − D_F⁻¹ A_FF D_F⁻¹), an
approximate ideal restriction; direct interpolation; A_c = R A P with
the entries below 10⁻⁴ of the diagonal dropped; a dense inverse at 64
unknowns or fewer. The cycle is the one of the scalar multigrid on the
same level operations (CPU and GPU, float32 hierarchy), with an inverse
diagonal that is zero on the C points: Jacobi on the F points, once
before and once after the coarse correction.

**In the coupled algorithm** (`foehncraft-flow::coupled`). One hierarchy
per velocity component, shared by components with the same matrix and
built when a sweep first needs it; kept while the momentum matrix has
changed by at most 2 % since its setup (`momentum_reuse`) and while its
cycle still halves the residual of the first sweep of an outer
iteration; stored with its matrix in the restart file. A component
whose part of a Krylov vector is zero up to rounding takes no sweep. A
failed setup leaves the sweeps of that outer iteration to BiCGStab, with
a warning in the log.

**One thread without a pool** (`backend::cpu`, `device::cpu`). A pool of
one thread and loops of one task run on the calling thread instead of
being handed to a worker (same values); found with the measurements
below.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Splitting of a 5-point Laplacian and of a convection matrix | C and F points; every F point has a strong C neighbour; coarsening along the flow | – | asserted (`splitting_of_a_laplacian`, `splitting_follows_the_flow`) |
| Convection–diffusion in a solid-body rotation, 32² to 256², cell Péclet number 0 to 300 | factor per cycle with the diagonal raised by 2 % (the cap 0.98); without; operator complexity | 0.35; 0.5 to 0.7; 4 | 0.07 to 0.22 (0.15 and 0.16 on 128² and 256²); 0.25 to 0.53; 2.4 to 3.5 |
| Scalar transport on hexahedra, tetrahedra, polyhedra and a mixed mesh in 3D, cell Péclet number 0.5 and 20 (`tests/reduction.rs`, on the backend, float32 hierarchy) | residual after 8 cycles; operator complexity | 10⁻³; 5.5 | passes; 2.2 to 5.4 |
| The same hierarchy in two multigrid slots, float32 and float64, repeated | cycles | bit-identical repeat; the precisions agree to 10⁻⁵ | asserted (`precision_and_slots`) |
| Setup on 1 and 4 threads | hierarchy | identical | asserted (`setup_does_not_depend_on_threads`) |
| Matrix without couplings; dense inverse of a non-symmetric matrix | cycle solves; A A⁻¹ = I | 10⁻¹² | asserted |
| Cavity 24², converged to 10⁻¹¹ | lid force with the multigrid sweep, the BiCGStab sweep, two cycles, `momentum_reuse = 0`, a float64 hierarchy | 10⁻⁸ | asserted (`momentum_sweeps_give_the_same_result`) |
| Cavity 12 × 12 × 6 with a symmetry plane (the components have different matrices) | result against the BiCGStab sweep; a hierarchy per component | 10⁻⁷ | asserted (`multigrid_sweep_with_symmetry_planes_in_3d`) |
| Cavity 16² (restart at 10 of 24) and 24² (restart at 70 of 90, with kept hierarchies) | velocity, pressure, fluxes, kept hierarchies of the continued run | bit-identical | asserted (`coupled_restart_continues_bit_identically`) |
| Cavity 150², 8 iterations, 1 and 4 threads | residuals, Krylov count, boundary forces | bit-identical | asserted |
| Cavity 70² with `max_levels = 1` (the coarsest level would exceed the dense limit) | the run goes on with BiCGStab sweeps; one warning in `flow.log` | – | asserted (`failed_multigrid_setup_falls_back_to_bicgstab`) |
| Input | ranges of `momentum_cycles`, `momentum_reuse`, `[solver.coupled.multigrid]` | refused with CASE005 | asserted |

Runs on one thread, to the iteration from which x_r stays within 1 % of
its final value (step) or to 10⁻⁸ (cavity); Krylov iterations per outer
iteration in brackets:

| Case | SIMPLEC | Coupled, BiCGStab sweep | Coupled, multigrid sweep |
|---|---|---|---|
| Step, k-ε-ζ-f₀, 8 200 cells (`backstep_zeta_f0_coupled.toml`) | 525 iterations, 18.7 s | 65 (10.0), 8.3 s | 65 (5.3), 7.0 s |
| Step, k-ε-ζ-f₀, 32 800 cells | 1 600 iterations, 243 s | 175 (14.2), 131 s | 150 (6.0), 74 s |
| Cavity at Re = 1000, 128² (`cavity_re1000_coupled.toml`) | 1 397 iterations, 28.3 s | 452 (2.2), 18.8 s | 449 (2.1), 12.1 s |
| Step, k-ε (`scalable`) and k-ω SST, 8 200 cells | 600 and 675 iterations | 75 and 80 | 75 and 80 (2.2 and 5.1) |
| Channel at Re_τ = 395, k-ε-ζ-f₀ | about 2 100 iterations | 323 | 324 |

Final values as in part 6e: x_r/H = 6.385 (k-ε-ζ-f₀), 5.287 (k-ε), 6.81
(SST), 6.35 on the fine mesh; the cavity deviates from Ghia et al. by
0.0040 and 0.0136 with both sweeps.

**Findings.** (1) One cycle per component replaces a BiCGStab sweep of
10 iterations at a third of its cost, and the Krylov count per outer
iteration halves on the step (so do the pressure solves). (2) The fine
step settles in 150 iterations instead of 175 and in 74 s instead of
131 s; against SIMPLEC the run time gains a factor of 2.7 on the coarse
step, 3.3 on the fine one and 2.3 on the cavity. (3) A setup costs as
much as twenty to thirty cycles, and while the flow forms the matrix
changes by more than 2 % in most iterations: 96 setups in the 65
iterations of the coarse step (1.4 s of 7.0 s, against 0.3 s for all
cycles), 140 in 150 on the fine one, 28 in 449 on the cavity. A setup
in every outer iteration made the cavity slower than with BiCGStab,
which is why hierarchies are kept. Measured on step matrices, a
hierarchy from a matrix 1.3 %, 2.7 % and 5.4 % away converges with 0.27,
0.34 and 0.6 per cycle; from 6.4 % on it does not. (4) Without the cap
on the local step (relaxation 0.999) the fine step settles in about 120
iterations, but every outer iteration runs into the limit of 30 Krylov
iterations: the cap stays at 0.98 until the Schur complement
approximation is better. (5) The BiCGStab sweeps took a quarter of the
run time on the fine mesh; the cycles and their setups take 16 %, and
the turbulence model 35 % (its linear solves 21 %) and the pressure
solves 25 %. (6) Thread hand-off: with a pool of one
thread every kernel call went to the worker thread and back, about
2 000 times per outer iteration of the cavity. Running such loops on
the calling thread took the coarse step with SIMPLEC from 31 s to
18.7 s and with the coupled algorithm and the BiCGStab sweep from 21 s
to 8.3 s, without changing a value. The times of part 6e above are from
before that change. (7) Larger pseudo-time steps on small laminar
cases, which cost two to three times the default with the BiCGStab
sweep (finding 7 of part 6e), are now as cheap as it: cavity at
Re = 1000 on 64², SIMPLEC 584 iterations in 3.0 s, coupled with
`factor` 1, 4, 20: 624, 275, 170 iterations in 3.7 s, 2.8 s, 2.6 s.

## Multigrid for the turbulence equations (part 6g)

Specification: [requirements step 22](../requirements/step22-coupled-solver.md)
sections 4.2 and 7.8 (D180).

**The method** (`foehncraft-flow::turbulence`). With the coupled
algorithm and pseudo-time steps (steady runs), the transport equations
of the k-ε family (k, ε and ζ of k-ε, realizable k-ε, k-ε-ζ-f and
k-ε-ζ-f₀; multigrid slots 2 to 4) take one V-cycle of the reduction
multigrid of part 6f as the preconditioner of BiCGStab instead of the
diagonal (`[solver.turbulence] preconditioner = "reduction"`, the
default there). A hierarchy is kept while the matrix of its equation
has changed by at most 2 % since its setup (`reuse`, the measure of
`momentum_reuse`), and dropped when a solve with it needs more than 8
BiCGStab iterations. A failed setup leaves the solves of that equation
to Jacobi until a setup succeeds, with one warning in the log. The
restart file holds the matrices of the kept hierarchies
(`turb_amg_matrix_<slot>`). The elliptic f equation keeps CG with the
aggregation multigrid.

**The steps.** With the multigrid the relaxations default to 0.995
(`solver.turbulence.relaxation`, `solver.momentum.relaxation`) instead
of 0.98. A relaxation α acts as a local pseudo step α/(1−α) ρV/a_P,
about 25 h²/ν_eff in a flat cell at the wall at 0.98: it shrinks with
the cells, the rule min(Δτ, 0.5 T) of D178 does not. k-ω SST keeps
Jacobi and 0.98, laminar cases 0.98, transient runs and SIMPLEC are
unchanged; a key set in the case wins.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Channel at Re_τ = 395, k-ε-ζ-f₀, coupled, to 10⁻¹⁰ | driving force with the multigrid at 0.995 against Jacobi at 0.98; setups counted and fewer than three per iteration; fewer iterations; SST keeps the old defaults | 10⁻⁶ | asserted (`turbulence_multigrid_gives_the_same_result`) |
| The same channel, restart at 30 of 60 iterations (hdf5) | velocity, pressure, turbulence fields, kept hierarchies in the restart file | bit-identical | asserted (`turbulence_multigrid_restart_continues_bit_identically`) |
| Input | `reuse` outside [0, 0.2] | refused with CASE005 | asserted |

Runs on one thread, quiet machine, to the iteration from which x_r stays
within 1 % of its final value (step) or to the tolerance (channel);
part 6f in brackets:

| Case | Iterations | Run time |
|---|---|---|
| Step, k-ε-ζ-f₀, 8 200 cells (`backstep_zeta_f0_coupled.toml`) | 40 (65) | 7.5 s (7.0 s) |
| Step, k-ε-ζ-f₀, 32 800 cells | 65 (150) | 56 s (74 s) |
| Step, k-ε (`scalable`), 8 200 cells | 45 (75) | 5.7 s (4.5 s) |
| Channel at Re_τ = 395, k-ε-ζ-f₀ | 231 (323) | – |
| Step, k-ω SST; cavity at Re = 1000 | unchanged: 80; 449 | unchanged |

Final values: x_r/H = 6.385 and 6.353 (k-ε-ζ-f₀, coarse and fine), 5.287
(k-ε), as before. The step case file converges to its tolerance in 117
iterations with every check passed.

**Findings.** (1) The turbulence relaxation 0.98 made the outer count
grow with the mesh (65 to 150 iterations per four times the cells): on
the fine step, momentum caps of 0.99 and 0.995 gave 150 and 145
iterations, `pseudo_time.factor = 2` 140, `pseudo_time.turbulence = 1`
125, the turbulence relaxation 0.99 95. (2) With Jacobi a relaxation of
0.99 or more runs the BiCGStab solves of k, ε and ζ into their limit of
50 in every outer iteration and x_r jitters or runs away; with the
multigrid they need 1 to 2 iterations. (3) Both relaxations at 0.995
together give 65 iterations on the fine step, the count of the coarse
step before: the count grows by 1.6 per four times the cells instead of
2.3. (4) An iteration costs more: the momentum cap 0.995 doubles the
Krylov iterations (9 instead of 5 per outer iteration) and with them
the pressure solves, and while the flow forms the hierarchies of the
turbulence equations are rebuilt in most solves (114 setups in 120
solves on the coarse step, 172 in 195 on the fine one; a setup costs
about as much as twenty cycles). The coarse steps take about the same
time as before; the fine step is 1.3 times faster, and the gain grows
with the mesh. Against SIMPLEC the fine step is now 4.3 times faster
(243 s). (5) SST on the step: with accurate solves of k and ω (the
multigrid, or Jacobi to 10⁻⁴) the residuals grow again after about 160
iterations and x_r wanders from 6.81 to 6.85 H at any relaxation from
0.95 to 0.995; the loose Jacobi solves of today damp this. A relaxation
of ν_t of 0.2 cures it at the caps 0.98 (206 iterations to the
tolerance, x_r = 6.809 H) but not at the larger ones, so SST keeps the
old defaults. (6) The aggregation multigrid of the scalar equations
(`preconditioner = "amg"`, part 6d) fails with pseudo steps: its
coarsest level is not positive definite.

## The pressure solve of the coupled preconditioner (part 6h)

Specification: [requirements step 22](../requirements/step22-coupled-solver.md)
sections 7.9 and 7.10 (D181).

**The method** (`foehncraft-flow::coupled`). The pressure solve inside
the SIMPLE-type preconditioner of the coupled algorithm stops at the
relative tolerance 0.3 instead of 0.1 (`pressure_tolerance`): the
Krylov iteration corrects what it leaves. The aggregation multigrid of
its matrix Ŝ is built anew when Ŝ has changed by more than 5 % since
the setup (`pressure_reuse`, the measure of `momentum_reuse`), besides
every `amg_rebuild` iterations as before. The restart file and the
state kept across an intermediate result hold the matrix of the last
setup, as before.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Cavity 24², Re = 100, coupled, to 10⁻¹¹ | a new hierarchy whenever Ŝ changes (`pressure_reuse = 0`) against the periodic rebuild alone (`1`): more setups, the periodic count exactly every 10 iterations, the same lid force | 10⁻⁸ | asserted (`pressure_hierarchy_follows_its_matrix`) |
| Coupled tests of parts 6e to 6g, restart included (hdf5) | results with the new defaults | as before | asserted |
| Input | `pressure_reuse` outside [0, 1] | refused | asserted |

Runs on one thread, quiet machine, to the iteration from which x_r stays
within 1 % of its final value; part 6g with the tolerance 0.1 in
brackets:

| Case | Iterations | Run time | Pressure solves |
|---|---|---|---|
| Step, k-ε-ζ-f₀, 8 200 cells (`backstep_zeta_f0_coupled.toml`) | 40 (40) | 6.2 s (8.0 s) | 0.69 s (2.02 s) |
| Step, k-ε-ζ-f₀, 32 800 cells | 65 (65) | 39.8 s (56.9 s) | 5.81 s (20.0 s) |
| Step, k-ε (`scalable`), 8 200 cells | 45 (45) | 4.7 s (6.4 s) | 0.85 s (2.03 s) |
| Step, k-ω SST, 8 200 cells | 80 (80) | 5.4 s | 0.82 s |
| Cavity at Re = 1000 | 442 to the tolerance (449) | – | – |

The step case file converges to its tolerance in 115 iterations with
every check passed.

**Findings.** (1) Two published Schur complements were measured against
the SIMPLE-type one and do not pay on the step: SIMPLER (a pressure
prediction first) saves a few Krylov iterations for a second pressure
solve, the least-squares commutator runs into the Krylov limit near the
inlet, outlet and walls (section 7.9). (2) The pressure solves, not the
Schur complement, were the cost: a hierarchy up to nine outer
iterations old needed 9 to 20 CG iterations per solve where a fresh one
needs two. (3) The tolerance 0.3 alone gains 8 % on the coarse meshes
and nothing on the fine one (its solves stop at the limit of 20 either
way); with the hierarchy following its matrix the fine step needs 71 CG
iterations per outer iteration instead of 160 and runs 1.4 times
faster, 6.1 times faster than SIMPLEC (243 s). (4) Laminar cases do not
trigger the rule: the pseudo-time term holds D̂, and the cavity builds
the same 45 hierarchies either way.

## Thermal wall function and turbulent Prandtl number (part 6f)

Specification: [requirements stage 6](../requirements/stage6-turbulence.md)
section 5 (D135), `foehncraft-kernels::turb::thermal`.

**The method.** A passive scalar of a turbulent flow diffuses with
ρ(Γ + ν_t/Pr_t) on every face: Pr_t is constant (`turbulent_prandtl`,
0.85) or Kays & Crawford's correlation of ν_t/ν and Pr. On a wall face
the diffusivity is Γ_w = μ y⁺/T⁺, so the face conductance is
h_w = ρ u\*/T⁺ with T⁺ of Kader (1981; all y⁺, any Pr; the default) or
Jayatilleke (1969; log region, blended to Pr y⁺ by Kader's Γ). Fixed-flux
walls take b = q d_n/Γ_w, so the wall value is φ_P + q/h_w. y⁺ is that of
the momentum wall treatment (u\* of the k-ε family, the blended u_τ of
SST); with the scalable treatment T⁺ is taken at ỹ⁺ = max(y⁺, 11.06),
the virtual wall distance of its log law (part 6f, pull request 87:
before, the thermal law saw the true y⁺ and the channel Nu was +28 % on a
y⁺ ≈ 8 mesh). The default is `kader` where the specification wrote `none` with
`low-re`: the two agree for y⁺ below 1, and `kader` is also right at
larger y⁺, so no switch is needed. Both kernels are single-source
(GPU and CPU). Cases without a turbulence model do not run them.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| T⁺ of Kader for Pr 0.01, 0.71, 7, 100 at y⁺ 1, 10, 100; Jayatilleke; Kays–Crawford | formulas evaluated apart in Python | 1e-12 relative | asserted (`turb::thermal`) |
| Limits | T⁺ → Pr y⁺ at the wall, log slope 2.12, Pr_t → 2 Pr_t∞ at ν_t → 0 | 1e-3 | asserted |
| Heated pipe, Re = 10⁵, Pr = 0.71, k-ε scalable, 36 cells | Nu against Gnielinski (1976), 180.2 | ±10 % | +1.3 % (Kader), +0.6 % (Kays–Crawford), +4.7 % (Jayatilleke); `none` −39 % (`heated_pipe_meets_gnielinski`) |
| Isothermal flat plate, Re = 5·10⁶ per unit length, 0.25 ≤ x ≤ 1 | St against 0.0287 Re_x^−0.2 Pr^−0.4 | ±10 % | +1.7 % at 0.29, +7.5 % at 0.97 (`heated_plate_meets_the_stanton_correlation`) |
| Input | `turbulent_prandtl` ≤ 0, unknown wall function | refused | asserted |
| Heated channel, Re_τ = 395, Pr = 0.71, 12 cells, y⁺ 32 (`channel_heated_keps.toml`) | Nu on 2h against 36.39 of the DNS of Abe, Kawamura & Matsuo (JAXA DNS database, Re_b = 13 849) | ±7 % | −1.0 % (constant Pr_t), −2.0 % (Kays–Crawford) |
| Heated channel, Re_τ 180 … 1020, Pr 0.71 and 0.025 (table below) | Nu on 2h against the DNS | Pr 0.71: ±7 % both models; Pr 0.025: ±8 % Kays–Crawford | asserted for all eight points (`heated_channel_meets_the_dns_nusselt_number`) |
| Heated channel, Re_τ = 395, Pr = 0.71, 48 cells, y⁺ 7.8 | thermal law at the scalable floor | ±13 % | +10.5 % (+28 % with T⁺ at the true y⁺) (`the_thermal_law_follows_the_scalable_floor`) |

**The channel at all DNS points** (pull request 87). The case of
`channel_heated_keps.toml` with ν = 1/Re_τ, Γ = ν/Pr, the mass flow of
each DNS bulk Reynolds number (the Re_b of the file header of each
thermal case) and N uniform cells across 2h, chosen for a first-cell
y⁺ ≈ 30, the claimed band of the scalable wall functions. Reference: the
Nu of each file header, except Re_τ = 1020, Pr = 0.025 (header 278.71,
see below), where the Nu of the database's own profiles is used (9.81).
y⁺ and Re_τ are those of the run (u_τ from the driving force; standard
k-ε predicts a friction 3–7 % low at the fixed mass flow).

| Re_τ | Pr | Re_b | N | y⁺ | Nu DNS | Nu, Pr_t 0.85 | Nu, Kays–Crawford |
|---|---|---|---|---|---|---|---|
| 180 | 0.71 | 5 546 | 6 | 28.1 | 18.25 | 18.43 (+1.0 %) | 18.22 (−0.2 %) |
| 395 | 0.71 | 13 849 | 12 | 31.5 | 36.39 | 35.98 (−1.1 %) | 35.62 (−2.1 %) |
| 640 | 0.71 | 24 095 | 20 | 30.9 | 55.57 | 54.62 (−1.7 %) | 54.10 (−2.6 %) |
| 1020 | 0.71 | 41 441 | 32 | 31.3 | 83.00 | 82.87 (−0.2 %) | 82.15 (−1.0 %) |
| 180 | 0.025 | 5 731 | 6 | 28.9 | 5.48 | 6.16 (+12.4 %) | 5.80 (+5.9 %) |
| 395 | 0.025 | 14 147 | 12 | 32.1 | 6.60 | 7.88 (+19.4 %) | 7.07 (+7.2 %) |
| 640 | 0.025 | 24 428 | 20 | 31.3 | 7.89 | 9.57 (+21.3 %) | 8.42 (+6.7 %) |
| 1020 | 0.025 | 41 441 | 32 | 31.3 | 9.81 | 11.98 (+22.1 %) | 10.49 (+6.9 %) |

At Pr 0.71 both Pr_t models are within 2.6 % at every Re_τ. At Pr 0.025
(a liquid metal) conduction dominates: Pe_t = (ν_t/ν) Pr stays below
about 2, where the real Pr_t rises well above 0.85. The constant Pr_t
then over-predicts the turbulent heat flux and Nu by 12–22 %, growing
with Re_τ; Kays & Crawford, whose Pr_t tends to 2 Pr_t∞ at small Pe_t,
brings it to +6 … +7 %. The rest is not the mesh: at Re_τ 395, Pr 0.025,
Nu is 6.94, 7.07, 7.13 and 7.08 on 8, 12, 24 and 48 cells. At Pr 0.71
the k-ε model without damping sets the limit below y⁺ 30: Nu 35.52,
35.62, 38.21, 40.20 on 8, 12, 24, 48 cells (y⁺ 47, 32, 16, 8). Meshes
finer than y⁺ ≈ 30 near heated walls take ζ-f or EB with the compound
treatment, not the k-ε family.

The Re_τ = 180 channel on 6 cells is exactly one-dimensional; its pressure
residual (about 10⁻⁴⁰) lies below the float32 range of the multigrid
hierarchy, so the preconditioned residual vanished and CG reported a
breakdown. CG now stops there, as BiCGStab did already.

The DNS data belong to JAXA (all rights reserved; citation of Abe et al.
requested), so the repository holds only the published Nusselt number.
`examples/studies/channel_heat_dns.py` downloads the tables to a
git-ignored folder, recomputes U_b⁺, Θ_b⁺ and Nu for the eight cases
(Re_τ 180, 395, 640, 1020 at Pr 0.71 and 0.025; the integrated Nu
agrees with the file headers to 0.4 % except `ch1020pr0025.dat`, whose
header repeats the Re_b of the Pr = 0.71 case and gives Nu = 278.71
where the profiles give 9.81) and compares our run.

Not in this part: the heated pipe with the other models, and the plate for x below 0.2, where
the k-ε start of the turbulent layer, not the thermal law, sets St
(50 % low).

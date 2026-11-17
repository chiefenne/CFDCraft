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
([below](#k-ε-ζ-f-and-the-scalar-multigrid-part-6d)). The commands are in
[recipe 63](../recipes/scalars.md#passive-scalars),
[recipe 67](../recipes/boundary-conditions.md#periodic),
[recipe 68](../recipes/turbulence.md#sst),
[recipe 71](../recipes/turbulence.md#k-epsilon) and
[recipe 73](../recipes/turbulence.md#zeta-f), the pictures in
[gallery 63](../gallery/scalars.md#passive-scalars),
[gallery 67](../gallery/boundary-conditions.md#periodic),
[gallery 68](../gallery/turbulence.md#sst),
[gallery 71](../gallery/turbulence.md#k-epsilon) and
[gallery 73](../gallery/turbulence.md#zeta-f), the generated
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

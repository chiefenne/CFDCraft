# Stage 6 – RANS turbulence and wall treatments

Stage 6 of the [stage plan](../requirements.md#7-stage-plan) brings the
RANS turbulence models and their wall treatments
([specification](../requirements/stage6-turbulence.md)). This page grows
with the pull requests of the stage. Part 6a gives every later model its
transport equation, its wall distance and its solver path
([below](#passive-scalars-and-the-wall-distance-part-6a)), part 6p the
periodic boundaries that the channel cases need
([below](#periodic-boundaries-part-6p)), part 6b the k-ω SST model
([below](#k-ω-sst-part-6b)). The commands are in
[recipe 63](../recipes/scalars.md#passive-scalars),
[recipe 67](../recipes/boundary-conditions.md#periodic) and
[recipe 68](../recipes/turbulence.md#sst), the pictures in
[gallery 63](../gallery/scalars.md#passive-scalars),
[gallery 67](../gallery/boundary-conditions.md#periodic) and
[gallery 68](../gallery/turbulence.md#sst), the generated
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
| | GCI_fine of U_b⁺ (four meshes, 96 to 768 cells across) | ≤ 0.5 % | CHANNEL_GCI |
| Homogeneous decay in a periodic box, k₀ = 1, ω₀ = 10, t = 0.5, 10 / 20 / 40 steps | BDF2 time order of k and ω against the exact solution (F₁ = 0) | 1.8 – 2.2 | DECAY_ORDERS |
| Steady plate, restart after half the iterations | fields at the end | bit-identical | bit-identical |
| Plate on 1 and 4 threads | fields after 8 iterations | bit-identical | bit-identical |
| Poisson wall distance on a plane wall | d against the exact distance | round-off | passed |

**Findings.** (1) The C_f of the plate converges monotonically but at a
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

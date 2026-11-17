# Turbulence models

Part of the [Example gallery](../gallery.md): the RANS turbulence models of stage 6 with their wall treatments (D129–D139), checked against the NASA Turbulence Modeling Resource (TMR), DNS and correlations.

| # | Example | Stage | Feature |
|---|---|---|---|
| 68 | [k-ω SST on the TMR flat plate and in a channel](#sst) | 6 | `[turbulence] model = "k-omega-sst"`, `wall_treatment`, inlet turbulence values, `roughness`, derived field `turbulence-features` |
| 71 | [k-ε family with wall functions: flat plate, channel, Moody chart](#k-epsilon) | 6 | `model = "k-epsilon"`, `"realizable-k-epsilon"`, `wall_treatment = "scalable"`, `"standard"`, `"compound"`, `epsilon` at inlets, `roughness_constant` |

---

<a id="sst"></a>
## 68. k-ω SST on the TMR flat plate and in a channel

![Left: the skin friction along the TMR flat plate on 256 x 192 cells against White's correlation and the TMR value of SST at x = 0.97. Middle: the velocity at x = 0.97 and in the channel at Re_τ = 395 in wall units against Spalding's law. Right: C_f at x = 0.97 with the automatic wall treatment for first-cell y⁺ from 2 to 100 against the TMR value, with the ±2 % band](../figures/gallery/sst.svg)

**The model.** The k-ω SST model of Menter, Kuntz & Langtry (2003): the
blend F₁ between k-ω near walls and k-ε away from them, the eddy
viscosity ν_t = a₁k/max(a₁ω, S F₂) and the production limiter
P̃ = min(ν_t S², 10 β* k ω). k and ω are two transport equations of the
generic scalar path of stage 6a, solved after each pressure correction
with a point-implicit (Newton) linearisation of the destruction, so k and ω
stay positive without clipping in a converging run. Two wall treatments:

```toml
[turbulence]
model = "k-omega-sst"
wall_treatment = "low-re"   # or "automatic" (the default)
```

`low-re` integrates to the wall (y⁺ ≲ 1): k = 0 and ω = 60ν/(β₁ d₁²) on
the wall face, d₁ the height of the wall cell, as in the TMR's own codes.
`automatic` takes any first-cell y⁺: the wall shear stress and ω in the
first cell are blended between their viscous-sublayer and log-law values,
and so is the production of k there.

**The TMR flat plate** (`tmr_plate_sst.toml`). The 2D zero-pressure-gradient
plate at Re = 5·10⁶ per metre with the TMR's free-stream values
(k = 2.25·10⁻⁷ U², ν_t/ν = 0.009), on the built-in mesher's copy
of the TMR grid family (point counts and stretching, not the TMR files):
four meshes from 32 × 24 to 256 × 192 cells, first-cell y⁺ 0.09 on the
finest. The skin friction at x = 0.97 converges monotonically:

| Cells | C_f(0.97) | against TMR | C_D (one side) | against TMR |
|---|---|---|---|---|
| 768 | 2.6937e-3 | +0.11 % | – | – |
| 3 072 | 2.6682e-3 | −0.83 % | 2.8259e-3 | −0.80 % |
| 12 288 | 2.6868e-3 | −0.14 % | 2.8477e-3 | −0.04 % |
| 49 152 | 2.7007e-3 | +0.37 % | 2.8606e-3 | +0.41 % |
| TMR, SST, finest grid | 2.6907e-3 | | 2.8488e-3 | |

The validation band of the specification, |C_f − C_f,TMR| ≤ 1 % plus
the numerical uncertainty, holds on every mesh but the first. The
verification band does not: the observed order of the three finest meshes
is 0.43 for C_f (GCI 1.87 %, extrapolated 2.741e-3) and 0.74 for C_D
(GCI 0.84 %), where the specification asks for p ≥ 1.5 and GCI ≤ 1 %.
This is recorded as a finding: on the coarsest mesh the stretched grid
puts only a few cells into the laminar leading-edge region, and the
convergence is not yet in its asymptotic range; the TMR itself reports
its finest-grid values after five grids. White's correlation lies 5 to 9 %
above the computed C_f along the plate, the TMR value of SST 6.4 % at
x = 0.97, so the model, not the code, sets that difference. u⁺ at x = 0.97
follows Spalding's law (κ = 0.41, B = 5.0) within 2.3 % for
1 ≤ y⁺ ≤ 300.

**The automatic treatment** on the second mesh with the wall grading
changed so that the first cell sits at y⁺ ≈ 2 to 100:

| First-cell y⁺ at x = 0.97 | 1.7 | 4.6 | 9.3 | 27 | 92 |
|---|---|---|---|---|---|
| C_f against the TMR value | +2.3 % | +15.6 % | +14.8 % | −0.2 % | −0.7 % |

In the sublayer and in the log layer the result holds; in the buffer
layer (3 ≲ y⁺ ≲ 20) the blended wall values give a wall shear stress up
to 16 % high, outside the ±2 % band of the specification. This is the
known weakness of the blending of the original formulation; until it is
improved, TURB001 warns when more than a tenth of a wall lies in that
range, and meshes should put the first cell at y⁺ ≲ 1 or y⁺ ≳ 25.

**The channel** (`channel_sst.toml`). Fully developed flow between two
walls at Re_τ = 395 on one periodic module (stage 6p) driven by the bulk
velocity of the DNS of Moser, Kim & Mansour (1999), 96 cells across,
first-cell y⁺ 0.28. The driving force gives u_τ and with it
Re_τ = 397.8 and U_b⁺ = 17.43: −0.7 % from the DNS (17.55) and +1.1 %
from Dean's correlation C_f = 0.073 Re_b^−¼ (17.24), inside the ±3 %
band. On three meshes (96, 192 and 384 cells across) U_b⁺ is 17.427, 17.381 and 17.346: the observed order is 0.40, the extrapolated value 17.237 (Dean: 17.24) and GCI_fine 0.78 %, outside the verification band of the specification (finding, as for the plate: the convergence is slow). The fourth mesh did not reach its tolerance in 30 000 iterations and is left out. The centreline velocity U_c⁺ = 19.64 lies about 2.4 % below the DNS
(about 20.1), inside ±4 %.

| Check | Result |
|---|---|
| Homogeneous decay of k and ω (periodic box, no walls) against the exact ODE solution | BDF2 time order of k 1.92 and 1.96, of ω 2.07 and 2.03 (band 1.8–2.2) |
| Restart in the middle of a steady SST run | bit-identical to the run without restart |
| 1 against 4 threads on the flat plate | bit-identical |
| Wall distance by the Poisson PDE on a plane wall | exact to round-off |
| Bad input (unknown model, constants, central convection in a steady run, k without ω) | refused with CASE002/CASE004/CASE005/TURB006 |
| Inlet turbulence values out of range, defaulted, or given at an outlet | BC029, TURB009, BC028 |

| | |
|---|---|
| Code | kernels `turb::sst` (closure, blending, production), `turb::wall` (low-re, automatic and Aupoix wall values), `turb::invariants` (the D166 features), `scalar` (the fractional first-cell fix); `foehncraft_flow::turbulence` (`Turbulence`, settings, health counts, `summary`), `turbulence::wall` (thickness, Poisson distance), `turbulence::inlet` (intensity, length scale, viscosity ratio), `derived` (`turbulence-features`) |
| Tests | `tests/turbulence.rs`, `turbulence/inlet.rs`, `turbulence/wall.rs`, kernel tests of `turb::sst` and `turb::wall` |
| Case files | `flow/tmr_plate_sst.toml`, `flow/channel_sst.toml` |
| Recipe | [68](../recipes/turbulence.md#sst) |

```
target/release/foehn flow examples/cases/flow/tmr_plate_sst.toml --out output/plate \
    --set 'verification.levels=[0, 1, 2, 3]'
target/release/foehn flow examples/cases/flow/channel_sst.toml --out output/channel
for r in "2 7983.6" "5 2859.1" "10 1298.9" "30 361.1" "100 83.3"; do set -- $r
  target/release/foehn flow examples/cases/flow/tmr_plate_sst.toml --out output/yp$1 \
      --set 'turbulence.wall_treatment="automatic"' --set "mesh.generate.grading.y.ratio=$2"
done
python3 -m tools.plot.sst_chart docs/figures/gallery/sst.svg output/plate output/channel \
    output/yp2 output/yp5 output/yp10 output/yp30 output/yp100
```

[↑ back to top](#turbulence-models)

---

<a id="k-epsilon"></a>
## 71. k-ε family with wall functions: flat plate, channel, Moody chart

![Left: C_f at x = 0.97 on the TMR flat plate with the standard k-ε model against the first-cell y⁺ from 1 to 100, scalable and compound wall treatments, with White's correlation and its ±5 % band. Right: the Darcy friction factor of rough pipes from Re 3·10⁴ to 10⁷ and k_s/D from 0 to 10⁻² against Colebrook's equation](../figures/gallery/keps.svg)

**The models.** The standard k-ε of Launder & Spalding (1974) and the
realizable k-ε of Shih et al. (1995), with its variable C_μ and the
destruction C₂ε²/(k + √(νε)), share one kernel. The standard model is
bounded by Durbin's realizability limit by default, which is inactive in
an equilibrium boundary layer and keeps the impulsive start of the
near-wall cells in check; the ε destruction is linearised in the Picard
form, which lets ε follow k when the free stream decays. Inlets take k
with ε or with ω (ε = C_μ k ω), and a run can start from an SST restart
file.

```toml
[turbulence]
model = "k-epsilon"            # or "realizable-k-epsilon"
wall_treatment = "scalable"    # "standard", "compound"
```

**The wall treatments.** All take the velocity scale u* = C_μ^¼ k^½ of
the first cell, so they need no iteration. `scalable` evaluates the log
law at max(y*, 11.06) and fixes ε and the production of the first cell
from it; `standard` falls back to the viscous sublayer below 11.06;
`compound` blends sublayer and log law with the function of Popovac &
Hanjalić. Rough walls shift the log law by ΔB of Cebeci & Bradshaw.

**The flat plate** (`tmr_plate_keps.toml`). The TMR plate of gallery 68
on 64 × 48 cells, its free-stream values held by the sustaining terms
(without them the free stream of the k-ε family decays before the plate).
With the wall grading changed so that the first cell sits at y⁺ 1 to 100
at x = 0.97, `scalable` gives C_f within −2.6 to −0.3 % of White's
correlation at every y⁺ (left). `compound` agrees from y⁺ 30
(−2.3 to −3.6 %) but overpredicts below, by 21 % at y⁺ 11 and 58 % at
y⁺ 1: the high-Re models have no damping in the buffer layer. This is a
finding, and TURB001 warns below y⁺ 30; walls resolved to y⁺ ≈ 1 are for
the ζ-f and EB models. The realizable model gives 3 to 7 % less friction
(finding).

**The channel** (`channel_keps.toml`). Re_τ = 395 on 12 uniform cells
(first-cell y⁺ 32): U_b⁺ = 18.31, +4.3 % from the DNS, inside the ±5 %
band of wall functions; +4.3 to +4.9 % for every y⁺ from 9 to 47.

**Rough pipes** (`pipe_rough_keps.toml`, `examples/studies/moody_study.py`).
A wedge pipe of one periodic module; the friction factor follows from the
driving force (right). Smooth pipes lie 2 to 5 % below Colebrook's
equation. Fully rough pipes lie 8 to 11 % above with the default
C_s = 0.5 and within 3 % with C_s = 0.4 (finding). In the transitional
range the Cebeci–Bradshaw shift follows Nikuradse's sand-grain dip and
lies below Colebrook's monotonic curve for commercial pipes.

| Check | Result |
|---|---|
| Plate, `scalable`, y⁺ 1 to 100: C_f(0.97) against White | −2.6 % to −0.3 % (band ±5 %) |
| Plate, `compound`, y⁺ 30 / 60 / 100 | −2.3 / −3.6 / −3.4 %; below y⁺ 30 up to +58 % (finding, TURB001) |
| Channel Re_τ 395, `scalable`, y⁺ 32: U_b⁺ | 18.31, +4.3 % (band ±5 %) |
| Realizable k-ε, plate and channel | 3 to 7 % low friction (finding) |
| Rough pipe Re 10⁶, k_s/D 10⁻³ | f 0.02081, +4.4 % from Colebrook |
| Moody study, 16 runs | smooth −5.0 to −1.7 %, fully rough +8 to +11 % (C_s 0.4: +1 to +3 %) |
| Homogeneous decay of k and ε, BDF2 | second order (asserted 1.8 to 2.2) |
| SST restart continued with k-ε; 1 against 4 threads | converted; bit-identical |
| Bad input (SST treatment with k-ε, realizable with a limiter, ω ambient, ω and ε at one inlet) | CASE005, BC029 |

| | |
|---|---|
| Code | kernels `turb::keps` (both closures, sources, Durbin bound), `turb::wall` (`keps_wall_face`, `roughness_shift`, `compound_gamma`), `entry` (`KepsCells`, the k-ε branch of `WallFaces`, the field update of `WallCells`); `foehncraft_flow::turbulence` (model and treatment keys, ω–ε conversion of inlets, initial values and restart files) |
| Tests | `tests/keps.rs`, kernel tests of `turb::keps` and `turb::wall`, `turbulence/inlet.rs` |
| Case files | `flow/tmr_plate_keps.toml`, `flow/channel_keps.toml`, `flow/pipe_rough_keps.toml`; study `examples/studies/moody_study.py` |
| Recipe | [71](../recipes/turbulence.md#k-epsilon) |

```
target/release/foehn flow examples/cases/flow/tmr_plate_keps.toml
target/release/foehn flow examples/cases/flow/channel_keps.toml
target/release/foehn flow examples/cases/flow/pipe_rough_keps.toml
python3 examples/studies/moody_study.py output/moody --jobs 4 --chart output/moody.svg
```

[↑ back to top](#turbulence-models)

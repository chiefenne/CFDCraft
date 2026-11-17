# Stage 6 specification: RANS turbulence and wall treatments

Specification of stage 6 of the [stage plan](../requirements.md#7-stage-plan),
drafted on 2026-11-17 from the maintainer's stage-6 brief, before any code
(rule: requirements before code). It details the decisions **D129–D139** of
[`docs/requirements.md`](../requirements.md#2-decisions) and builds on the
earlier stage-6 decisions D12, D13, D28 (1), D29, D36, D48, D49, D51, D55,
D104, D108, D109 (1), D111, D116 and D124 (W4). Where this page and a D row
differ, the D row wins until this page is amended. *Amended 2026-11-17:*
the optional curvature and rotation correction and the quadratic
constitutive relation of SST (**D170**, section 3.6) and the explicit
algebraic stress model as a later item (**D171**, section 3.7).

Rules that hold for every item below:

- **Clean room** (D21): equations, constants, limiters and wall treatments
  come from the cited papers. No code of other solvers is read for
  implementation or copied. Every constant in the code carries the paper,
  and where possible the equation or table number, in a doc comment. Constants
  marked ◇ below are written from the literature as remembered at drafting
  time. They must be checked against the paper itself when the model is
  implemented (checklist item of the pull request), and the check is recorded
  in the doc comment.
- float64 on the reference path (D7). Gather-based assembly, no
  floating-point atomics, fixed summation order in global face order
  (D31, D93). The CPU and GPU run the same single-source kernel bodies
  (D8), bitwise equal.
- Adjoint-ready kernels (D59): pure functions of their inputs, generic over
  the `Real` type (dual numbers included), with smooth or frozen variants of
  every non-differentiable switch (section 2.6).
- The project says *reacting flows*, never the other word, and names no
  other CFD code.

Contents:

1. [Scope and reconciliation with earlier decisions](#1-scope-and-reconciliation-with-earlier-decisions)
2. [Common framework (D129)](#2-common-framework-d129)
3. [Turbulence models (D130–D133)](#3-turbulence-models-d130d133)
4. [Momentum wall treatment (D134)](#4-momentum-wall-treatment-d134)
5. [Thermal wall treatment and turbulent Prandtl number (D135)](#5-thermal-wall-treatment-and-turbulent-prandtl-number-d135)
6. [Numerics, GPU kernels and stability (D136)](#6-numerics-gpu-kernels-and-stability-d136)
7. [Boundary conditions, initial values and checks (D137)](#7-boundary-conditions-initial-values-and-checks-d137)
8. [Case-file keys (D138)](#8-case-file-keys-d138)
9. [Verification and validation suite, GCI pass bands (D139)](#9-verification-and-validation-suite-gci-pass-bands-d139)
10. [Implementation plan](#10-implementation-plan)
11. [Questions for the maintainer](#11-questions-for-the-maintainer)
12. [References](#12-references)


## 1. Scope and reconciliation with earlier decisions

### 1.1 What stage 6 delivers

| Item | Decision | Source of the requirement |
|---|---|---|
| Generic transport equation for scalars, passive scalars | D116, D51, D129 | earlier stage 6 scope |
| Exact wall distance from the wall-face hierarchy, PDE form as option | D48, D108 (1) | earlier stage 6 scope |
| k-ω SST (Menter, Kuntz & Langtry 2003) with the production limiter of the paper, Kato–Launder and Durbin options | D12, D130 | M1 decision 2026-11-14; brief 2026-11-17 adds the limiters |
| Standard k-ε (Launder & Spalding 1974) and realizable k-ε (Shih et al. 1995) | D12, D131 | M1 decision |
| k-ε-ζ-f (Hanjalić, Popovac & Hadžiabdić 2004) | D12, D132 | M1 decision (in place of RNG k-ε) |
| Elliptic-blending k-ε (Billard & Laurence 2012) | D133 (section 3.4) | brief 2026-11-17; maintainer decision 2026-11-17: joins ζ-f |
| Lag elliptic-blending k-ε (Lardeau & Billard 2016) | P11 (**proposed**) | brief, "optionally" |
| All-y⁺ momentum wall treatments: automatic (SST), compound (ζ-f, EB, optional for k-ε), scalable and standard wall functions, rough walls | D13, D134 | M1 decision; brief |
| Thermal law of the wall for any Pr (Kader 1981; Jayatilleke P-function), constant and Kays–Crawford turbulent Prandtl number | D135 | brief; needed by D111 |
| Turbulence numerics: bounded schemes, point-implicit sources, positivity, coupling order, bounded BDF2, AMG for scalars | D28 (1), D36, D136 | earlier scope; brief |
| Boundary conditions, inlet turbulence, checks and diagnostic codes | D29, D30, D137 | earlier scope; brief |
| Case-file keys | D138 | brief |
| Verification suite and GCI pass bands | D139 | brief; replaces the acceptance text of the stage-6 row where they differ |
| Moving reference frame, driven and rotational periodic boundaries with the affine map, study layer, statistics, derived fields, robustness controls | D55, D29, D109 (1), D62, D71, D72, D49 | earlier scope, **kept** in stage 6 unchanged (section 10) |

### 1.2 Reconciliation

| Topic | Earlier decision | Brief of 2026-11-17 | This specification |
|---|---|---|---|
| Near-wall model for separation and heat transfer | ζ-f (Hanjalić 2004) with the compound wall treatment, in M1 (2026-11-14); PANS on ζ-f later (P7); literature line Hanjalić, newer Basara | Elliptic-blending k-ε (Billard & Laurence 2012), optionally the lag variant; compare with ζ-f | ζ-f stays as decided. EB k-ε **joins** ζ-f, after it in stage 6 (maintainer decision 2026-11-17, D133; section 3.5 compares the two). Nothing of ζ-f is dropped. |
| SST variant | SST 2003 | SST 2003 with Kato–Launder production limiter or Durbin realizability constraint | Kept. The 2003 limiter is the default; Kato–Launder and Durbin are options (D130). The ω production is written with the limited P̃_k, as corrected for the 2003 paper (section 3.1). |
| Realizable k-ε | M1 | "fast-converging industrial baseline" | Kept, unchanged (D131). |
| Standard k-ε | M1, for literature data | not mentioned | Kept (D12); its known stagnation-point anomaly becomes a test of the limiters (section 9). |
| Wall treatment | D13: standard and scalable WF for k-ε, automatic for SST, compound for ζ-f, thermal field with L1, roughness | All-y⁺ momentum and thermal treatments, CWT for ε models | Kept and detailed (D134, D135). The compound treatment is also offered for the k-ε family as an option. |
| Thermal wall functions and turbulent heat flux | "once L1 exists" (D13); energy equation is L1 (after M1) | stage 6: thermal law of the wall and Pr_t; stage 7: energy equation, GGDH, buoyancy | The thermal wall function and Pr_t are specified **and verified in stage 6 on a temperature carried as a constant-property passive scalar** (forced convection with constant properties is one-way coupled, so the passive scalar is exact for it). The energy equation, GGDH, buoyancy production and CHT belong to the stage-7 specification (D150 on), which takes the interfaces of D129 and D135 (agreed with the stage-7 thread). |
| Linear solver of the elliptic equations | D12: f of ζ-f "with the scalar path of stage 2, no AMG" | "leverage the AMG + CG solver for α" | Amended: f (ζ-f) and α (EB) are symmetric positive definite Helmholtz problems; they are solved with CG and the aggregation AMG of stage 3 by default, Jacobi-CG as option (D136). This uses D36's AMG for scalars. |
| Stage numbering | Stage 7 = parallel performance, internal flows, v1.0 | Stage 7 = energy, buoyancy, CHT | Not decided here. The stage-7 specification thread handles it. Stage 6 is the same stage in both plans. |
| Keys file | `docs/keys.md` is generated from the code (`foehn schema --markdown`, D77) | "TOML schema additions (`docs/keys.md`)" | The new keys are specified in section 8. `docs/keys.md` is regenerated in each implementation pull request, as the generator requires. |
| Crate names | `foehncraft-kernels` (no_std), `foehncraft-fv`, `foehncraft-cuda`, `foehncraft-flow` | same | Names match; new modules listed in section 10. |
| Wall-film and periodic constraints | D124 W4 (per-face exchange coefficients), D109 (1) (affine map on periodic pairs) | – | The wall treatment fills the shear, heat and mass slots of `derived::wall_exchange` (D124 W4). All turbulence variables cross periodic pairs through the affine map with a = 1, b = 0 (D109 (1)). |


## 2. Common framework (D129)

### 2.1 Notation

Incompressible, constant density ρ in stage 6 (variable density with
stage 7). ν = μ/ρ, ν_t = μ_t/ρ. Velocity gradient ∇u with components
∂u_i/∂x_j, from the least-squares gradient of D24.

- Strain rate S_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i), rotation
  W_ij = ½(∂u_i/∂x_j − ∂u_j/∂x_i).
- Invariants: S = (2 S_ij S_ij)^½, Ω = (2 W_ij W_ij)^½, and
  |S| = (S_ij S_ij)^½ = S/√2 where a paper uses that form (ζ-f, EB,
  realizable k-ε). Each model section states which one it uses.
- Production P_k = ν_t S² per unit mass (incompressible, Boussinesq),
  with the model's limiter where it has one.
- y: wall distance from D48 (exact closest-point distance by default).
  y_P: normal distance of the centroid of a near-wall cell from its wall
  face (section 4.1), used by the wall treatment.

### 2.2 Generic transport equation

Every turbulence variable is solved with the generic equation of D116,

  ∂(φ)/∂t + ∇·(F φ) − φ ∇·F − ∇·(Γ_φ ∇φ) = S_u + S_p φ,  S_u ≥ 0, S_p ≤ 0,

with F the face volume flux (mass flux over ρ), Γ_φ = ν + ν_t/σ_φ, and the
term −φ∇·F making the convection bounded while continuity is not yet
converged (section 6.1). The elliptic equations of ζ-f and EB use the
same interface with no time derivative and no convection:

  −∇·(L² ∇ψ) + ψ = S_ψ   (ψ = f or α).

### 2.3 Fields

Registered in the field registry (D115), written to the VTKHDF output on
request and to the restart file always:

| Field | Unit | Models | Notes |
|---|---|---|---|
| `k` | m²/s² | all | turbulent kinetic energy |
| `epsilon` | m²/s³ | k-ε family, ζ-f, EB; derived ε = β* k ω for SST | |
| `omega` | 1/s | SST | specific dissipation rate |
| `nu_t` | m²/s | all | kinematic eddy viscosity |
| `zeta` | – | ζ-f | v²/k |
| `f` | 1/s | ζ-f | elliptic relaxation function |
| `phi` | – | EB | v²/k (same physical quantity as ζ) |
| `alpha` | – | EB | elliptic blending parameter, 0 at walls, 1 far away |
| `wall_distance` | m | all | D48 |
| `turbulent_time_scale` | s | all | T of the model (k/ε, or bounded T of ζ-f and EB; 1/(β*ω) for SST) |
| `v2` | m²/s² | ζ-f (ζk), EB (φk) | wall-normal variance, for GGDH in stage 7 |
| `wall_normal` | – | EB | n = ∇α/\|∇α\| (zero where \|∇α\| vanishes), for GGDH in stage 7 |
| `reynolds_stress` | m²/s² | all | estimate of u'_i u'_j, symmetric, 6 components: Boussinesq (2/3)k δ_ij − 2ν_t S_ij for every model (stage 7 may refine it for ζ-f and EB) |

Wall fields (per wall face, boundary blocks of the output, D124 W6):
`y_plus`, `y_star`, `u_tau`, `wall_shear_stress`, `blending` (the Γ of the
compound treatment, or the automatic weight), and per scalar
`heat_transfer_coefficient` (h_w) and `t_plus`.

### 2.4 Coupling to momentum

- Momentum diffusion uses μ_eff = μ + ρν_t in the full stress
  ∇·[μ_eff(∇u + (∇u)ᵀ)]; the transposed part is a deferred explicit term
  (it vanishes for constant μ_eff and incompressible flow, not for a
  varying ν_t).
- The isotropic part (2/3)ρk of the Reynolds stress is absorbed into a
  modified pressure p̂ = p + (2/3)ρk, which the pressure equation solves.
  Pressure boundary values given by the user are static pressures p and
  are converted with the boundary k. The output `pressure` is the static
  pressure p = p̂ − (2/3)ρk; `modified_pressure` is written on request.
  Forces and the pressure coefficient use p.
- The momentum wall shear at wall faces comes from the wall treatment
  (section 4) as a per-face effective viscosity, so the momentum kernel of
  stage 4 is unchanged in structure.

### 2.5 Hooks for later stages and correction fields

- **Model-correction fields (D104):** each turbulence equation takes the
  optional per-cell multiplier β on its production term and an additive
  source. Default off.
- **Buoyancy production (for stage 7):** the k, ε and ω equations each take
  a per-cell explicit part S_u ≥ 0 and implicit part S_p ≤ 0 named
  `buoyancy_production`, zero in stage 6. Stage 7 specifies G_b and its
  coefficients (D150 on).
- **Sustaining terms** (Spalart & Rumsey 2007) for the free-stream decay:
  optional ambient values (section 3.1 and 3.2).

### 2.6 Smooth and published switches

The models contain min, max, |·| and clipping. The primal solution uses the
**published hard form** by default, so that results can be compared with
the literature and the NASA Turbulence Modeling Resource (TMR) values. All
such operations go through one helper in `foehncraft-core::math`:

- `pmin(a, b)`, `pmax(a, b)`, `pabs(a)`: hard form (select instructions on
  the GPU, no branch);
- the same with a smooth variant, `smin(a, b; w) = ½(a + b − √((a − b)² + (w·m)²))`
  with m = ½(|a| + |b|) and relative width w (default 10⁻³), and
  `smax` likewise. It is C¹ and its error is bounded by w·m/2;
- the linearisation and the adjoint (D59 (3)) always use the smooth
  variant or the frozen switch; `turbulence.smooth_limits = true` makes
  the primal use it as well. A test checks that the flat-plate C_f with
  smooth limits differs from the hard form by less than 0.1 %.

### 2.7 Software structure

- The model is an `enum TurbulenceModel` chosen on the host. Kernels are
  monomorphised per model (const generic or one function per model); no
  trait objects and no runtime branch on the model inside a kernel.
- A host-side trait `TurbulenceModule` in `foehncraft-flow` gives each
  model's fields, equations and their solve order, the wall-treatment kind
  and the log columns; it plugs into the module and step schedule of D118.
- Kernels (no_std, single source): `foehncraft-kernels::turb`
  (`sst`, `keps`, `zeta_f`, `eb`, `wall`, `thermal`, `invariants`). The
  device trait `FlowDevice` of stage 4 gets the matching operations; the
  CUDA device runs the same bodies (D8).


## 3. Turbulence models (D130–D133)

Equations are written per unit mass for constant ρ, with D/Dt = ∂/∂t + u·∇.

### 3.1 k-ω SST, 2003 version (D130)

Menter, Kuntz & Langtry (2003); constants of Menter (1994) where the 2003
paper keeps them; the TMR's SST-2003 page records the differences and the
correction of the ω production.

  Dk/Dt = P̃_k − β* k ω + ∇·[(ν + σ_k ν_t)∇k]

  Dω/Dt = (γ/ν_t) P̃_k − β ω² + ∇·[(ν + σ_ω ν_t)∇ω] + 2(1 − F₁) σ_ω2 (1/ω) ∇k·∇ω

- ν_t = a₁ k / max(a₁ ω, S F₂), a₁ = 0.31, with the strain-rate invariant
  S (2003 change from the vorticity of 1994).
- P_k = ν_t S², P̃_k = min(P_k, 10 β* k ω) (2003: factor 10, applied in
  both equations because the ω production uses P̃_k).
- F₁ = tanh(arg₁⁴), arg₁ = min[max(√k/(β* ω y), 500ν/(y² ω)), 4 σ_ω2 k/(CD_kω y²)],
  CD_kω = max(2 σ_ω2 (1/ω) ∇k·∇ω, 10⁻¹⁰) (2003 floor; in units of the
  per-unit-mass form, documented as such in the code).
- F₂ = tanh(arg₂²), arg₂ = max(2√k/(β* ω y), 500ν/(y² ω)).
- Blending of each constant: c = F₁ c₁ + (1 − F₁) c₂.
- Constants: β* = 0.09, κ = 0.41; set 1: γ₁ = 5/9, β₁ = 0.075, σ_k1 = 0.85,
  σ_ω1 = 0.5; set 2: γ₂ = 0.44, β₂ = 0.0828, σ_k2 = 1.0, σ_ω2 = 0.856.

**Production options** (`turbulence.production_limiter`, D138):

- `"menter"` (default): P̃_k as above.
- `"kato-launder"` (Kato & Launder 1993): P_k = ν_t S Ω instead of ν_t S²,
  then the 2003 clip. In a pure shear flow S = Ω and nothing changes; at a
  stagnation point Ω ≈ 0 and the spurious production vanishes. Ω comes
  from the absolute-frame velocity gradient, which is the solved variable
  in a rotating zone (D55); a note (TURB004) says so when both are active.
- `"durbin"` (Durbin 1996): the eddy viscosity is bounded by the
  realizability limit ν_t ≤ C_lim k/(√6 |S|), applied as
  ν_t = min(ν_t,SST, C_lim k/(√6 |S|)), C_lim = 1 (Durbin: "≤ 1").
  Derivation recorded in the code: the normal stresses
  (2/3)k − 2ν_t λ_i of the Boussinesq form stay non-negative when
  ν_t ≤ k/(3 λ_max), and λ_max ≤ (2/3)^½ |S| for a trace-free S_ij.
  Followed by the 2003 clip.

**Sustaining terms** (Spalart & Rumsey 2007; TMR "SST-sust"), optional
(`turbulence.ambient`): + β* ω_amb k_amb in the k equation and
+ β ω_amb² in the ω equation, so that free-stream values do not decay
between inlet and body. Off by default. The flat plate of section 9 uses
the TMR settings with them off and records the decay.

### 3.2 Standard and realizable k-ε (D131)

One implementation with two closures (D12).

  Dk/Dt = P_k − ε + ∇·[(ν + ν_t/σ_k)∇k]

**Standard** (Launder & Spalding 1974): ν_t = C_μ k²/ε,

  Dε/Dt = (ε/k)(C_ε1 P_k − C_ε2 ε) + ∇·[(ν + ν_t/σ_ε)∇ε],

C_μ = 0.09, C_ε1 = 1.44, C_ε2 = 1.92, σ_k = 1.0, σ_ε = 1.3. The
Kato–Launder option of 3.1 applies (P_k = ν_t S Ω), which is the form
Kato & Launder proposed for this model.

**Realizable** (Shih, Liou, Shabbir, Yang & Zhu 1995):

  Dε/Dt = C₁ S ε − C₂ ε²/(k + √(νε)) + ∇·[(ν + ν_t/σ_ε)∇ε],

C₁ = max(0.43, η/(η + 5)), η = S k/ε, C₂ = 1.9, σ_k = 1.0, σ_ε = 1.2;
ν_t = C_μ k²/ε with the variable

  C_μ = 1/(A₀ + A_s U* k/ε), U* = (S_ij S_ij + W_ij W_ij)^½,
  A_s = √6 cos φ_s, φ_s = ⅓ arccos(√6 W̃), W̃ = S_ij S_jk S_ki / |S|³,

W̃ clipped to [−1/√6, 1/√6] before the arccos (smooth clip in the
linearisation). A₀ = 4.0 ◇ (a value 4.04 is also in use; the paper's value
is taken). In a rotating zone the rotation term uses the rotation relative to
the frame, as in the paper; stated in the doc comment.

**Sustaining terms** for both, the analogue of 3.1 (Spalart & Rumsey
2007), optional with ambient values k_amb and ε_amb: the k equation gains
+ ε_amb and the ε equation + C_ε2 ε_amb²/k_amb (realizable:
+ C₂ ε_amb²/(k_amb + √(ν ε_amb))), so that the ambient state is a steady
solution of the decay.

### 3.3 k-ε-ζ-f (D132)

Hanjalić, Popovac & Hadžiabdić (2004); constants of the paper ◇.

  ν_t = C_μ ζ k T

  Dk/Dt = P_k − ε + ∇·[(ν + ν_t/σ_k)∇k]

  Dε/Dt = (C_ε1 P_k − C_ε2 ε)/T + ∇·[(ν + ν_t/σ_ε)∇ε]

  Dζ/Dt = f − (ζ/k) P_k + ∇·[(ν + ν_t/σ_ζ)∇ζ]

  L²∇²f − f = (1/T)(c₁ + C₂′ P_k/ε)(ζ − 2/3)

- T = max[min(k/ε, a/(√6 C_μ |S| ζ)), C_T (ν/ε)^½],
  L = C_L max[min(k^{3/2}/ε, k^{1/2}/(√6 C_μ |S| ζ)), C_η ν^{3/4}/ε^{1/4}].
- C_μ = 0.22, C_ε1 = 1.4(1 + 0.012/ζ), C_ε2 = 1.9, c₁ = 0.4, C₂′ = 0.65,
  σ_k = 1, σ_ε = 1.3, σ_ζ = 1.2, a = 0.6, C_T = 6, C_L = 0.36, C_η = 85 ◇.
- Wall values (integration to the wall): k_w = 0, ζ_w = 0,
  ε_w = 2ν k_P/y_P², f_w = −2ν ζ_P/y_P²; with the compound treatment these
  become blended values (section 4.4).
- Bounds: 0 ≤ ζ ≤ 2 (v² ≤ 2k), enforced by the source split of section
  6.2, not by clipping; a clip remains as the last guard and is counted.
- Curvature correction for the U-bend (P8) and PANS on ζ-f (P7) stay later.
  The later literature of Basara and co-workers (P7) builds on this form.

### 3.4 Elliptic-blending k-ε (D133)

Billard & Laurence (2012), the "BL-v²/k" model. Variables k, ε, φ = v²/k and
the blending parameter α:

  ν_t = C_μ φ k T,  T = max(k/ε, C_T (ν/ε)^½),
  L = C_L max(k^{3/2}/ε, C_η ν^{3/4}/ε^{1/4})

  α − L²∇²α = 1,  α_w = 0

  Dk/Dt = P_k − ε + ∇·[(ν + ν_t/σ_k)∇k]

  Dε/Dt = (C′_ε1 P_k − C*_ε2 ε)/T + E + ∇·[(ν + ν_t/σ_ε)∇ε]

  Dφ/Dt = (1 − α³) f_w + α³ f_h − (φ/k) P_k + (2/k)(ν_t/σ_k) ∇φ·∇k + ∇·[(ν + ν_t/σ_φ)∇φ]

with the near-wall part f_w = −(ε/2)(φ/k) and the homogeneous part
f_h = −(1/T)(C₁ − 1 + C₂ P_k/ε)(φ − 2/3), a coefficient C′_ε1 that grows
towards the wall through (1 − α³), the variable C*_ε2 and the extra
near-wall term E of the paper. The structure above is ◇ as well. The exact forms of C′_ε1, C*_ε2 and E and
**all constants (C_μ, C_ε1, C_ε2, C_ε3, A₁, C₁, C₂, C_T, C_L, C_η, σ_k,
σ_ε, σ_φ) are taken from Billard & Laurence (2012) when the model is
implemented** and are not written here from memory. Acceptance requires
that the channel at Re_τ = 395 reproduces the paper's own channel result
(section 9) before any other case.

- Wall values: k_w = 0, φ_w = 0, α_w = 0, ε_w = 2ν k_P/y_P² (or the
  compound value of section 4.4).
- Bounds: 0 ≤ α ≤ 1 by the discrete maximum principle (section 3.5),
  0 ≤ φ ≤ 2.
- The wall-normal direction n = ∇α/|∇α| is exported (section 2.3); it is
  what GGDH and algebraic flux models of stage 7 use near walls.

**Lag EB k-ε (P11, proposed as possible later).** Lardeau & Billard (2016)
add a transport equation for the lag between stress and strain to the EB
model (a fifth equation, conference paper, AIAA 2016-1600). Proposal: keep
it as P11 and decide after the EB verification of stage 6, because the
gain is in non-equilibrium flows that are not in the stage-6 cases.

### 3.5 EB k-ε compared with ζ-f, and recommendation

| Aspect | ζ-f (Hanjalić et al. 2004) | EB k-ε (Billard & Laurence 2012) |
|---|---|---|
| Elliptic operator | −∇·(L²∇f) + f | −∇·(L²∇α) + α |
| Matrix | symmetric, positive definite; M-matrix on orthogonal meshes, non-orthogonal part deferred (D26) | same |
| Condition number | grows like 1 + (L/h)²_max; L spans C_η η (near walls) to k^{3/2}/ε (outer flow), so the coefficient jumps by orders of magnitude: AMG needed | same operator, same scaling, same need for AMG |
| Right-hand side | dimensional (1/s), depends on P_k/ε, T and ζ − 2/3; changes sign | the constant 1 |
| Wall boundary value | f_w = −2νζ_P/y_P²: depends on the first cell (ζ_P and y_P²), couples back to ζ, and changes from face to face on polyhedral and prismatic wall layers with uneven y_P | α_w = 0: independent of the first cell |
| Bounds of the solution | none for f; f enters the ζ source with either sign | 0 ≤ α ≤ 1 by the discrete maximum principle (with an M-matrix and deferred correction limited, section 6.2) |
| Lagging | f must follow ζ every outer iteration | α depends only on L (k, ε); it can be solved every N outer iterations (`solver.elliptic.every`) |
| Near-wall source in the v²/k equation | through f (elliptic, global) | (1 − α³) f_w + α³ f_h: local, algebraic |
| Literature line | Hanjalić; PANS-ζ-f of Basara and co-workers (P7) builds on it | Billard & Laurence; elliptic blending also in later Reynolds-stress and algebraic heat-flux models |

The two matrices are equally conditioned; the difference is in the wall
value and the right-hand side. On arbitrary polyhedral meshes with uneven
first-cell distances the α equation is the better behaved one: its
boundary value is exact and its solution is bounded, whereas the f wall
value scales with 1/y_P² and links neighbouring wall faces of different
height through ζ. The compound treatment of ζ-f (section 4.4) removes most
of that sensitivity by blending f_w, which is why the 2004 and 2007 papers
pair them.

**Decision (maintainer, 2026-11-17): EB k-ε joins ζ-f, delivered after it in stage 6**
(pull request 6e, section 10). It reuses the k and ε equations, the v²/k
transport equation and the elliptic solver of ζ-f, so the extra code is one
source kernel and its tests. ζ-f stays the decided near-wall model and the
base of P7. After the stage-6 validation (channel, backward-facing step,
impinging jet; section 9) the maintainer decides which of the two is the
recommended near-wall model for heat transfer in the user documentation.
Alternatives: EB replaces ζ-f (loses the base of P7 and the maintainer's
literature line), or EB follows later as a P item (stage 6 stays smaller).


### 3.6 Optional corrections of SST: curvature and rotation, QCR (D170)

Two options of the SST model of 3.1, both off by default, so the published
SST 2003 stays the reference path. Both work on the velocity gradient that
the solver already has, so they add no transport equation.

**Curvature and rotation correction** (`turbulence.curvature_correction`):
the correction of Spalart & Shur (1997) applied to SST by Smirnov & Menter
(2009). The production P̃_k in the k and the ω equation is multiplied by

  f_r1 = max[0, min(f_rotation, 1.25)] ◇,

  f_rotation = (1 + c_r1) · 2 r*/(1 + r*) · [1 − c_r3 atan(c_r2 r̃)] − c_r1,

  r* = S/Ω,  r̃ = 2 Ω_ik S_jk [DS_ij/Dt + (ε_imn S_jn + ε_jmn S_in) Ω_m^rot] / (Ω D³),

  D² = max(S², 0.09 ω²),  c_r1 = 1.0, c_r2 = 2.0, c_r3 = 1.0 ◇,

with S and Ω the strain-rate and vorticity invariants of 3.1, Ω_ik the
rotation tensor in the frame of the zone and Ω^rot the angular velocity of
a rotating zone (D55; zero elsewhere). In a pure shear flow r* = 1 and
r̃ = 0, so f_r1 = 1 and the model is SST exactly; this is the first check.
DS_ij/Dt is the material derivative of the strain-rate tensor: its
convective part u·∇S_ij from the least-squares gradient of the six
components (one extra gradient of a symmetric tensor per outer iteration),
its time part from the BDF levels in a transient run. The upper limit 1.25
and the clip at 0 are those of the paper; the smooth limits of 2.6 replace
them in `smooth_limits` runs. A note (TURB code) says when the correction
is combined with Kato–Launder production, since both act on stagnation
and swirl regions.

**Quadratic constitutive relation** (`turbulence.qcr`): the Boussinesq
stress τ_ij = 2 μ_t S*_ij − (2/3) ρ k δ_ij (S* trace-free) is extended by
the QCR of Spalart (2000),

  τ_ij,QCR = τ_ij − C_cr1 (O_ik τ_jk + O_jk τ_ik),
  O_ik = 2 W_ik / √(∂u_m/∂x_n ∂u_m/∂x_n),  C_cr1 = 0.3 ◇,

with W_ik the rotation tensor, and, for `"qcr2013"`, the isotropic term of
Mani, Babcock, Winkler & Spalart (2013),

  − C_cr2 μ_t √(2 S*_mn S*_mn) δ_ij,  C_cr2 = 2.5 ◇.

For a model that carries k, the TMR page on QCR records how the 2013 term
relates to the (2/3) ρ k δ_ij term; which form is used is decided at
implementation against the paper and the TMR page (◇). The non-linear part
enters the momentum equation as the divergence of an explicit stress, a
deferred correction in the stored face order (D26, D31); the implicit
diffusion stays μ + μ_t. Its effect is the anisotropy of the normal
stresses that drives the secondary flow in non-circular ducts and corners,
which a linear model cannot produce. The production P_k of the turbulence
equations keeps the linear form (as the papers do).

### 3.7 Explicit algebraic Reynolds-stress model, later item (D171)

The explicit algebraic Reynolds-stress model (EARSM) of Wallin & Johansson
(2000) on the k-ω platform recalibrated for it by Hellsten (2005) is
planned after stage 6 (not in v0.6): the Reynolds-stress anisotropy is an
explicit tensor polynomial in the normalised strain and rotation tensors,

  a_ij = β₁ S_ij + β₃ (Ω_ik Ω_kj − II_Ω δ_ij/3) + β₄ (S_ik Ω_kj − Ω_ik S_kj)
         + β₆ (…) + β₉ (…)  (the full basis of the paper; in two dimensions
         β₁, β₃, β₄ remain) ◇,

with the β coefficients from the self-consistent solution of the paper
(a cubic for the production-to-dissipation ratio, solved in closed form, so
the model stays explicit and has no iteration). The k-ω equations are those
of Hellsten (2005) with its constants and cross diffusion (not SST's) ◇.
It reuses the explicit-stress path of the QCR (3.6) and the invariant
features of D166, and gives the stress anisotropy that the internal flows
of the validation (U-bend, square-duct secondary flow, swirl) need beyond
QCR. Its pull request, constants check and verification (rotating channel,
square duct, U-bend) are planned with the stage after v0.6; until then
D171 fixes only the model and the interfaces it needs, which stage 6 keeps
open: an explicit stress term in the momentum assembly (used first by the
QCR) and the strain and rotation tensors of the invariants module.


## 4. Momentum wall treatment (D134)

### 4.1 Wall-face interface

For every wall face f with near-wall cell P (the owner cell of the boundary
face):

- y_P = (x_P − x_f)·n_f, the normal distance of the cell centroid from the
  face plane, with n_f the outward unit normal; a floor of 10⁻³ of the
  cell's characteristic length guards degenerate cells (counted in health).
- U_t = |u_P − u_w − ((u_P − u_w)·n_f) n_f|, the tangential velocity
  relative to the wall velocity (D29 moving walls), and its direction t_f.
- Inputs: k_P, ν (face value), roughness k_s (section 4.5), the model's
  constants, and for scalars Pr (or Sc) and Pr_t.
- Outputs per face: u_τ, y⁺ = u_τ y_P/ν, y* = C_μ^{1/4} k_P^{1/2} y_P/ν,
  the wall shear stress τ_w (vector along −t_f), the effective wall
  viscosity μ_w = τ_w y_P / U_t used by the momentum kernel, the blending
  weight, and the first-cell values the model needs (ω_P, ε_P, the
  production P_k,P, f_w or α-related values), plus the per-face exchange
  coefficients of D124 W4 (shear; heat and mass from section 5).

The interface is model-independent: one wall-face kernel per wall
treatment, selected on the host. The same face list serves the thermal
treatment (section 5) and later the wall film (D123, D124).

**Near-wall cells with several wall faces** (corner cells, prism layers in
concave corners, polyhedra cut by two walls): the cell values that the
treatment imposes (ω_P, ε_P, P_k,P) are the area-weighted means over the
cell's wall faces, gathered in global face order from a per-cell list of
wall faces (CSR, built once on the host, D93 (1)), so the result is
bitwise reproducible.

### 4.2 Laws of the wall

- Log law u⁺ = (1/κ) ln(E y⁺), κ = 0.41, E = 9.793 (B = 5.2 ◇).
- Spalding (1961), one formula for all y⁺:
  y⁺ = u⁺ + e^{−κB}[e^{κu⁺} − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6],
  solved for u⁺ by Newton iteration with a fixed number of steps (no
  data-dependent loop on the GPU: 8 steps from the log-law guess, then a
  residual check counted in health).
- Viscous sublayer u⁺ = y⁺.

### 4.3 Treatments per model

| Treatment | Models | Wall shear | First-cell turbulence values | Claimed y⁺ range |
|---|---|---|---|---|
| `low-re` (integration to the wall) | SST, ζ-f, EB | τ_w = μ U_t/y_P | SST: ω_P = 6ν/(β₁ y_P²) fixed in the first cell (Menter 1994; the limit y⁺ → 0 of the automatic treatment); ζ-f, EB: ε_w, f_w or α_w of section 3 | y⁺ ≲ 1 |
| `automatic` (Menter, Kuntz & Langtry 2003) | SST (default) | u_τ = (u_τ,vis⁴ + u_τ,log⁴)^{1/4}, u_τ,vis = (ν U_t/y_P)^½, u_τ,log = κ U_t/ln(E y⁺) (fixed number of Newton steps, as for Spalding's law) | ω_P = (ω_vis² + ω_log²)^½, ω_vis = 6ν/(β₁ y_P²), ω_log = u_τ/(√β* κ y_P) | any; validated 0.5 ≤ y⁺ ≤ 100 |
| `compound` (Popovac & Hanjalić 2007) | ζ-f, EB (default), k-ε family (option) | τ_w from the blended profile: u⁺ = e^{−Γ} y⁺ + e^{−1/Γ} u⁺_log, Γ = 0.01 (y⁺)⁴/(1 + 5 y⁺) ◇ | ε_P = e^{−Γ} 2νk_P/y_P² + e^{−1/Γ} C_μ^{3/4} k_P^{3/2}/(κ y_P); P_k,P blended likewise between the integrated and the log-law value τ_w u_τ/(κ y_P); f_w (ζ-f) blended between −2νζ_P/y_P² and 0 | any; validated 0.5 ≤ y⁺ ≤ 100 |
| `scalable` (Grotjans & Menter 1998) | k-ε family (default) | log law evaluated at ỹ* = max(y*, 11.06) | ε_P = C_μ^{3/4} k_P^{3/2}/(κ y_P) (with ỹ*), P_k,P = τ_w u_τ/(κ y_P) | y⁺ ≥ 1 without failure; accuracy claimed for y⁺ ≥ 30 |
| `standard` (Launder & Spalding 1974) | k-ε family | log law with u_τ = C_μ^{1/4} k_P^{1/2} | as scalable, without the floor | 30 ≤ y⁺ ≤ 300 |

The exponentials e^{−Γ} and e^{−1/Γ} are evaluated for every face (no
branch on the regime); Γ is floored at 10⁻³⁰⁰ before 1/Γ, and the in-house
exp (D93 (6)) returns exactly 0 below its underflow limit. For SST the
compound blending of Popovac & Hanjalić is not used: the automatic
treatment is the SST paper's own.

Popovac (2023) is cited for the later refinements of the compound treatment;
which of them enter is decided with the ζ-f pull request (6d) after reading
the paper (◇).

### 4.4 Imposing first-cell values

ω_P (SST) and ε_P (wall-function treatments) are fixed values in the
near-wall cell: the matrix row of that cell is replaced by a_P = 1, b = value,
off-diagonals 0, through a per-cell mask applied with select instructions
in the row-assembly kernel (no branch); the neighbour rows keep their
coupling. For P_k,P the cell's production is replaced by the area-weighted
wall-face value. The compound treatment of ζ-f also blends the f wall
value. These replacements are part of the assembly, so the dual-number
Jacobian (D59) sees them.

### 4.5 Rough walls

Equivalent sand-grain roughness k_s per wall boundary (`roughness`, m),
k_s⁺ = u_τ k_s/ν:

- **Log-law shift** (all wall-function and compound treatments): u⁺ is
  reduced by ΔB(k_s⁺) of Cebeci & Bradshaw (1977): ΔB = 0 for k_s⁺ < 2.25;
  transitional ΔB = (1/κ) ln[(k_s⁺ − 2.25)/87.75 + C_s k_s⁺]·sin[0.4258 (ln k_s⁺ − 0.811)]
  for 2.25 ≤ k_s⁺ < 90; fully rough ΔB = (1/κ) ln(1 + C_s k_s⁺) for
  k_s⁺ ≥ 90; C_s = 0.5 (`roughness_constant`) ◇. The pieces join
  continuously; the smooth helper of section 2.6 joins them C¹ for the
  linearisation.
- **SST** (`low-re` and `automatic`): rough-wall ω boundary value of
  Aupoix (2015), ω_w = (u_τ²/ν) ω_w⁺(k_s⁺) with the correlation of that
  paper ◇; with the automatic treatment ω_vis is replaced by this value.
- **ζ-f, EB**: the same ΔB in the compound profile.
- A warning (BC027) when k_s exceeds y_P on more than 10 % of a wall's
  area (the roughness elements then reach beyond the first cell centre).


## 5. Thermal wall treatment and turbulent Prandtl number (D135)

Stage 6 applies this to passive scalars (D51) and to a temperature carried
as a constant-property passive scalar (molecular Prandtl number
Pr = ν/a with a the thermal diffusivity, given as the scalar's
`diffusivity`). Stage 7 applies the same interface to the energy equation
(D92): the wall heat flux is then q_w = h_w (T_w − T_P) with the h_w below.

### 5.1 Turbulent diffusivity

Γ_φ = ν/Pr + ν_t/Pr_t (Sc and Sc_t for a mass-like scalar). Options per
scalar (`turbulent_prandtl_model`):

- `"constant"` (default): Pr_t = 0.85 (`turbulent_prandtl`).
- `"kays-crawford"` (Kays 1994, after Kays & Crawford):
  Pr_t = 1 / {1/(2 Pr_t∞) + C Pe_t/√Pr_t∞ − (C Pe_t)² [1 − exp(−1/(C Pe_t √Pr_t∞))]},
  Pe_t = (ν_t/ν) Pr, Pr_t∞ = 0.85, C = 0.3. It raises Pr_t near walls and
  at low Pr (liquid metals).

GGDH and other anisotropic heat-flux closures are stage 7 (D150 on); they
use the fields of D129.

### 5.2 Thermal law of the wall

Wall-face interface as in 4.1, with Pr and Pr_t added. Outputs: T⁺ (or φ⁺),
the face conductance and the heat-transfer coefficient
h_w = ρ c_p u*/T⁺ (for a passive scalar: the exchange coefficient
ρ u*/φ⁺), with u* = C_μ^{1/4} k_P^{1/2} for the k-ε family and ζ-f/EB,
and the blended u_τ of 4.3 for SST.

- **`"kader"`** (default, all y⁺, any Pr; Kader 1981, wall form with the
  outer-layer factor dropped, y/δ → 0):
  T⁺ = Pr y⁺ e^{−Γ_t} + [2.12 ln(1 + y⁺) + β(Pr)] e^{−1/Γ_t},
  β(Pr) = (3.85 Pr^{1/3} − 1.3)² + 2.12 ln Pr,
  Γ_t = 0.01 (Pr y⁺)⁴/(1 + 5 Pr³ y⁺).
  Kader's paper covers temperature and concentration, so the same form
  serves a mass-like scalar with Sc.
- **`"jayatilleke"`** (log region only): T⁺ = Pr_t [u⁺ + P(Pr/Pr_t)],
  P = 9.24 [(Pr/Pr_t)^{3/4} − 1][1 + 0.28 e^{−0.007 Pr/Pr_t}]
  (Jayatilleke 1969). Blended with the conduction limit T⁺ = Pr y⁺ by the
  Γ_t of the Kader form, so that it is also defined for small y⁺.
- **Integration to the wall** (`"none"`): with a `low-re` momentum
  treatment the scalar flux is the molecular conduction ν/Pr (φ_P − φ_w)/y_P.
- With rough walls the thermal roughness correction is not part of stage 6;
  a note says that the smooth-wall T⁺ is used.

Boundary kinds of a scalar at a wall: fixed value, fixed flux (the flux
enters the first cell directly; the wall value is reported from T⁺), and
zero flux (default). Robin conditions and conjugate interfaces are stage 7.


## 6. Numerics, GPU kernels and stability (D136)

### 6.1 Discretisation

- **Convection:** linear upwind as a deferred correction on implicit
  upwind (D25) with the Venkatakrishnan limiter on by default for the
  turbulence variables (K = 5); first-order upwind as an option and for
  the first N iterations by the ramp of D49 (2). The bounded central
  scheme is refused for turbulence variables in steady runs (TURB006) and
  allowed in transient runs with β ≥ 0.1. The deferred correction is
  limited so that it cannot make the right-hand side negative (6.2).
- **Bounded form:** the term −φ_P ∇·F (with ∇·F the net volume outflow
  of the cell from the current fluxes) is added implicitly, so that the
  convection operator is an M-matrix whatever the continuity residual.
- **Diffusion:** over-relaxed non-orthogonal correction, deferred (D26);
  face diffusivity Γ_f by linear interpolation of ν_t (harmonic mean as
  option for strongly varying ν_t). The non-orthogonal correction of the
  elliptic equations is limited like the convective deferred correction.
- **Gradients:** least squares (D24) of k and ω for the SST cross
  diffusion; the same gradients serve the limiter.
- **Time:** BDF2 with the bounded fall-back of D28 (1), first order in the
  first step; Crank–Nicolson is refused for turbulence variables (not
  bounded), with a note that BDF2 is used for them.

### 6.2 Point-implicit linearisation and positivity

Each source is split into S_u ≥ 0 and S_p ≤ 0 with values φ* of the
previous outer iteration (Patankar 1980), so that the matrix keeps a
positive diagonal and the solution stays positive for positive boundary
values (M-matrix property). Every split uses select operations
(`max(x, 0)`, `min(x, 0)`), evaluated in all cells: no data-dependent
branch, so no warp divergence on the GPU.

| Equation | Term | S_u | S_p |
|---|---|---|---|
| k (all) | production P̃_k | P̃_k | 0 |
| k (SST) | −β* k ω | 0 | −β* ω* |
| k (ε models) | −ε | 0 | −ε*/k* |
| ω | (γ/ν_t) P̃_k | γ P̃_k/ν_t* (ν_t* floored, 6.2 floors) | 0 |
| ω | −β ω² | β ω*² (Newton) | −2β ω* |
| ω | cross diffusion CD = 2(1 − F₁) σ_ω2 ∇k·∇ω/ω | max(CD, 0) | min(CD, 0)/ω* |
| ε (standard, ζ-f) | (C_ε1 P_k − C_ε2 ε)/T | C_ε1 P_k/T | −C_ε2/T |
| ε (realizable) | C₁ S ε − C₂ ε²/(k + √(νε)) | C₁ S ε* | −C₂ ε*/(k* + √(νε*)) |
| ζ | f − (ζ/k) P_k | max(f, 0) | min(f, 0)/ζ* − P_k/k* |
| φ (EB) | sources of 3.4 | positive parts | negative parts divided by φ* |
| f, α | elliptic | right-hand side (α: 1) | the unit coefficient (+1 on the diagonal) |
| all | deferred corrections D | max(D, 0) | min(D, 0)/φ* |
| k, ε, ω | buoyancy hook (stage 7) | its S_u | its S_p |

The divisions by φ* use φ* ≥ φ_floor. The Newton form of the ω
destruction makes the diagonal twice as strong as the Picard form and
converges faster; the Picard form (S_p = −β ω*) is an option.

**Floors** after each solve: k ≥ k_min, ε ≥ ε_min, ω ≥ ω_min with
k_min = 10⁻¹⁰ k_ref and ω_min, ε_min from k_min and the reference
viscosity ratio 10⁻², where k_ref is the largest inlet k (or the initial
k). With the splits above the floors should never act on a converging
run. Each activation is counted per equation in the health section
(D69) and gives TURB002 when more than 0.1 % of the cells needed it in an
iteration. ζ, φ are bounded to [0, 2] and α to [0, 1] the same way.

**Eddy-viscosity floor and bound:** ν_t ≥ 0 by construction; an optional
maximum viscosity ratio (`turbulence.max_viscosity_ratio`, default 10⁵)
caps ν_t/ν, counted in health.

### 6.3 Linear solvers

- k, ε, ω, ζ, φ: right-preconditioned BiCGStab with Jacobi (D36);
  AMG on the equation's own matrix as an option (`preconditioner = "amg"`)
  for stiff high-aspect-ratio layers (D36 *2026-11-15*).
- f, α: CG with the aggregation AMG of stage 3 (the matrices are
  symmetric positive definite, which the stage-3 AMG requires); Jacobi-CG
  as option. The hierarchy is rebuilt with `solver.amg_rebuild` like the
  pressure, and L² from the lagged k, ε.
- Tolerances per outer iteration: 10⁻² steady, 10⁻⁸ transient (as
  momentum), maximum iterations 50 (BiCGStab) and 100 (CG).

### 6.4 Coupling order

**SIMPLEC / SIMPLE outer iteration** (steps 1–7 of `simplec.rs` unchanged,
with μ_eff from the current ν_t and p̂ for the pressure):

8. velocity gradient of the corrected velocity; invariants S, Ω, |S|, W̃
   per cell (one invariants kernel);
9. wall treatment per wall face (section 4), then the per-cell gather of
   first-cell values;
10. the model's equations in this order, each assembled with the newest
    values of the ones before it:
    - SST: k, then ω (F₁, F₂ and CD_kω from the fields at the start of
      step 10);
    - standard and realizable k-ε: k, then ε;
    - ζ-f: k, ε, then f, then ζ;
    - EB: α (every `solver.elliptic.every` iterations, default 1), then
      k, ε, φ;
11. ν_t from the new fields, with optional relaxation α_νt (default 1);
    then the passive scalars (D51) with the new ν_t.

Under-relaxation defaults (steady): k, ε, ω, ζ, φ 0.8; f, α 1 (no time
derivative; they follow the lagged fields); transient 1.

**PISO** (plain): steps 8–11 once per time step after the pressure
correctors. **Iterative transient** (D28 (2), outer iterations per step):
once per outer iteration, so that a converged step is converged in all
variables. The residual of each turbulence equation (scaled as the
momentum residual) is part of the convergence test (`solver.tolerance` and
`outer_tolerance`).

**Wall distance** (D48) is computed once before the first iteration and
again only when the mesh moves (D56, L4).

### 6.5 GPU kernel design

- One **invariants** cell kernel (S, Ω, |S|, W̃ from ∇u), one **sources**
  cell kernel per model that writes S_u, S_p of all its equations and ν_t
  candidates, one **wall-face** kernel per treatment over the wall-face
  range, one **first-cell gather** kernel over the near-wall cells (CSR of
  wall faces per cell), and the generic face and row kernels of stage 2
  for convection and diffusion.
- No runtime branch on the model or regime inside a kernel: the model is a
  monomorphised instantiation; regimes are blended (compound, automatic)
  or selected with min, max and select.
- Transcendental functions (exp, ln, tanh, pow, cos, arccos, cube root as
  pow) come from the in-house routines of `foehncraft-core::math` on CPU
  and GPU alike (D93 (6)); no FMA contraction left to the compiler.
- All reductions (residual norms, y⁺ statistics, floor counts) are exact
  and order-independent (D31, D93 (4)); counts are integers.
- Memory per cell for SST: k, ω, ν_t, two old time levels of k and ω,
  S_u and S_p reused across equations, gradients of k and ω: about 20
  float64 values per cell on top of the flow solver; ζ-f and EB about 30.
  Measured and reported in the pull requests (D81 estimate updated).

### 6.6 Reproducibility, parity and adjoint

- Bitwise equal CPU and GPU results for each model on the flat plate and
  the channel (Colab T4 notebook, as for stage 5), and bitwise equal
  results on 1, 2 and 4 threads (CI).
- Dual-number Jacobian of every new kernel (sources, wall faces, Kader,
  Spalding Newton steps with fixed count) against central finite
  differences, as the stage-5 tests do (D59).
- Restart continues bitwise (D42), including the turbulence fields, the
  old time levels and the wall-treatment state.


## 7. Boundary conditions, initial values and checks (D137)

### 7.1 Inlets and backflow

Turbulence at velocity, mass-flow, static- and total-pressure inlets and for
the backflow of pressure outlets (keys already reserved in stage 5), one of:

- `turbulence_intensity` I and `turbulence_length_scale` l:
  k = (3/2)(I |u|)², ε = C_μ^{3/4} k^{3/2}/l, ω = k^{1/2}/(C_μ^{1/4} l);
- `turbulence_intensity` and `viscosity_ratio` r = ν_t/ν:
  ε = C_μ k²/(r ν), ω = k/(r ν);
- `k` with `epsilon` or `omega` directly (numbers or expressions, D53),
  as the TMR cases specify them;
- defaults when nothing is given: I = 0.05 and r = 10, with a note
  (TURB009) naming the boundary.

`zeta` (v²/k, ζ-f and EB, default 2/3) at inlets (D12). With a table
profile (D29) k and ε or ω can be columns of the table, so that a profile
from a periodic pipe or channel run (D111) feeds a jet or a duct.

### 7.2 Walls, symmetry, periodic, wedge

- Walls: section 4 (k_w = 0 in the integrated limit; zero normal gradient of
  k with wall functions), `roughness`, moving walls through U_t.
- Symmetry, slip wall, wedge: zero normal gradient of all turbulence
  scalars; f, α likewise (α is only 0 on no-slip walls).
- Periodic and rotational periodic: every turbulence variable through the
  paired-face map with a = 1, b = 0 (D109 (1)).
- Convective outlet: the turbulence scalars take the same convective
  condition as the velocity.

### 7.3 Initial values

`initial.turbulence_intensity` and `initial.viscosity_ratio`, or
`initial.k` with `initial.epsilon`/`initial.omega` (numbers or expressions);
by default the values of the first inlet in the order of the case file,
else I = 0.01 and r = 1 with the initial speed (or 1 m/s when that is 0).
ζ, φ start at 2/3, α at 1 with the elliptic equation solved before the
first iteration, f at 0.

### 7.4 Checks and diagnostic codes

The boundary-condition check (D30) and the validation of D39 gain these
codes (registry of D78, documented in `docs/interfaces.md` when
implemented):

| Code | Severity | Title |
|---|---|---|
| `BC027` | warning | A wall roughness k_s larger than the first-cell distance on more than 10 % of the wall's area. |
| `BC028` | warning | A turbulence value that the selected model does not use (`zeta` with SST or k-ε, `omega` with an ε model); BC009 stays for laminar cases. |
| `BC029` | error | Two incompatible turbulence specifications on one inlet (for example `turbulence_length_scale` and `viscosity_ratio`, or `k` with intensity). |
| `TURB001` | warning | First-cell y⁺ outside the claimed range of the wall treatment (section 4.3) on more than 10 % of a wall's area (reported per boundary with the area-weighted min, mean and max). |
| `TURB002` | warning | A positivity floor or bound acted on more than 0.1 % of the cells in an iteration (equation and cell count). |
| `TURB003` | warning | Inlet turbulence outside the usual range: I > 0.5, or a viscosity ratio below 10⁻² or above 10⁵. |
| `TURB004` | note | Kato–Launder production in a rotating zone: Ω is the absolute-frame vorticity. |
| `TURB005` | note | No wall boundary: the wall distance is infinite (F₁ = 0, F₂ = 0; the ζ-f and EB near-wall terms vanish). |
| `TURB006` | error | A convection scheme or time scheme that is not bounded for a turbulence variable (central in a steady run, Crank–Nicolson). |
| `TURB007` | warning | Model constants changed in `[turbulence.constants]`. |
| `TURB008` | warning | The ν_t cap `max_viscosity_ratio` acted (cell count). |
| `TURB009` | note | Inlet turbulence defaulted to I = 0.05, r = 10 (boundary named). |

The guidance table of D30 gets the recommended wall treatment per model
and y⁺ band and the inlet settings of the canonical cases.


## 8. Case-file keys (D138)

Additions to the case file (format version bumped with the first
implementation pull request, D126). `docs/keys.md` is generated from the
code; the table below is the target.

| Key | Type | Default | Unit | Run time | Description |
|---|---|---|---|---|---|
| `turbulence.model` | "laminar" \| "k-omega-sst" \| "k-epsilon" \| "realizable-k-epsilon" \| "k-epsilon-zeta-f" \| "elliptic-blending-k-epsilon" | `"laminar"` | | | Turbulence model (D12, D130–D133). |
| `turbulence.wall_treatment` | "automatic" \| "compound" \| "low-re" \| "scalable" \| "standard" | per model: SST `automatic`, ζ-f and EB `compound`, k-ε family `scalable` | | | Momentum wall treatment (D134). |
| `turbulence.production_limiter` | "menter" \| "kato-launder" \| "durbin" \| "none" | SST `menter`, k-ε `none` | | | Production form of SST and standard k-ε (3.1, 3.2). |
| `turbulence.durbin_coefficient` | number | `1.0` | | | C_lim of the Durbin bound. |
| `turbulence.ambient` | { k, omega \| epsilon } | none | SI | | Sustaining terms (Spalart & Rumsey 2007). |
| `turbulence.curvature_correction` | boolean | `false` | | | SST: the curvature and rotation correction of Smirnov & Menter (2009) (3.6, D170); constants `c_r1`, `c_r2`, `c_r3`. |
| `turbulence.qcr` | "none" \| "qcr2000" \| "qcr2013" | `"none"` | | | Quadratic constitutive relation (3.6, D170); constants `c_cr1`, `c_cr2`. |
| `turbulence.max_viscosity_ratio` | number | `1e5` | | yes | Cap of ν_t/ν. |
| `turbulence.smooth_limits` | boolean | `false` | | | Smooth min and max in the primal solution (2.6). |
| `turbulence.wall_distance` | "exact" \| "pde" | `"exact"` | | | D48. |
| `turbulence.constants.<name>` | number | the paper's | | | Constant override by its literature name (`beta_star`, `a1`, `c_mu`, …); TURB007. |
| `solver.turbulence.scheme` | "upwind" \| "linear-upwind" \| "central" | `"linear-upwind"` | | | Convection scheme of the turbulence equations. |
| `solver.turbulence.limiter` | "none" \| "venkatakrishnan" \| "venkatakrishnan-wang" | `"venkatakrishnan"` | | | Limiter. |
| `solver.turbulence.limiter_k` | number | `5.0` | | | Constant K of the limiter. |
| `solver.turbulence.relaxation` | number | 0.8 steady, 1 transient | | yes | Under-relaxation of k, ε, ω, ζ, φ. |
| `solver.turbulence.tolerance` | number | 1e-2 steady, 1e-8 transient | | yes | Relative tolerance per solve. |
| `solver.turbulence.max_iter` | integer | `50` | | yes | BiCGStab iterations per solve. |
| `solver.turbulence.preconditioner` | "jacobi" \| "amg" | `"jacobi"` | | | Preconditioner (D36). |
| `solver.turbulence.destruction` | "newton" \| "picard" | `"newton"` | | | Linearisation of −βω² and similar terms (6.2). |
| `solver.turbulence.viscosity_relaxation` | number | `1.0` | | yes | Relaxation α_νt of ν_t. |
| `solver.elliptic.tolerance` | number | 1e-2 steady, 1e-8 transient | | yes | f and α solves. |
| `solver.elliptic.max_iter` | integer | `100` | | yes | CG iterations. |
| `solver.elliptic.preconditioner` | "amg" \| "jacobi" | `"amg"` | | | Preconditioner of CG. |
| `solver.elliptic.every` | integer | `1` | | yes | Solve α (EB) every N outer iterations. |
| `boundary.<name> (type = "wall").roughness` | number | `0.0` | m | | Equivalent sand-grain roughness k_s. |
| `boundary.<name> (type = "wall").roughness_constant` | number | `0.5` | | | C_s of Cebeci & Bradshaw. |
| `boundary.<name> (inlets).k`, `.epsilon`, `.omega` | number \| string | none | SI | | Direct turbulence values (7.1). |
| `boundary.<name> (inlets).zeta` | number | `0.6667` | | | v²/k for ζ-f and EB. |
| `initial.turbulence_intensity`, `initial.viscosity_ratio`, `initial.k`, `initial.epsilon`, `initial.omega` | number \| string | from the first inlet | SI | | Start values (7.3). |
| `scalars.<name>.diffusivity` | number \| string | required | m²/s | | Molecular diffusivity (ν/Pr for a temperature). |
| `scalars.<name>.turbulent_prandtl` | number | `0.85` | | | Pr_t or Sc_t (D135). |
| `scalars.<name>.turbulent_prandtl_model` | "constant" \| "kays-crawford" | `"constant"` | | | D135. |
| `scalars.<name>.wall_function` | "kader" \| "jayatilleke" \| "none" | `"kader"` (`"none"` with `low-re`) | | | Thermal or mass law of the wall. |
| `scalars.<name>.initial` | number \| string | `0.0` | | | Start value. |
| `boundary.<name>.scalar.<scalar>.value`, `.flux` | number \| string | zero flux at walls, 0 at inlets | | | Per-physics sub-table of D116. |
| `output.fields` | adds `k`, `epsilon`, `omega`, `nu_t`, `zeta`, `f`, `phi`, `alpha`, `wall_distance`, `turbulent_time_scale`, `v2`, `wall_normal`, `reynolds_stress`, `modified_pressure`, the scalars | | | | Section 2.3. |
| `output.wall` | adds `y_plus`, `y_star`, `u_tau`, `blending`, `heat_transfer_coefficient`, `t_plus` | | | | Wall fields. |

Log (D125): residual columns `k`, `epsilon`, `omega`, `zeta`, `f`, `phi`,
`alpha` and one per scalar; the solution-settings matrix gets the
turbulence and elliptic rows (D127); the extremes block gets ν_t/ν max with
its cell and the wall y⁺ min, mean and max per wall boundary.


## 9. Verification and validation suite, GCI pass bands (D139)

### 9.1 How a verdict is formed

- **Discretisation uncertainty** from D67 (Celik et al. 2008) on the three
  finest of at least four systematically refined levels (refinement ratio
  r ≥ 1.3 per level, the same in each direction): observed order p,
  extrapolated value φ_ext and U_num = GCI_fine = 1.25 |e₂₁|/(r^p − 1)
  (relative). Oscillatory or divergent convergence is reported by D67
  and counts as a failed verification.
- **Verification verdict:** p inside the stated band and GCI_fine below
  the stated limit.
- **Validation verdict:** |φ_fine − φ_ref| ≤ B_model + U_num + U_ref, with
  B_model the model band in the tables below, U_ref the uncertainty of the
  reference (experimental uncertainty, or the spread of the DNS or TMR
  values), all relative. Model bands are set before the runs and are not
  changed to make a model pass; a model that misses its band is reported
  as a finding with the published result of the same model next to it.
- Every case runs on the hexahedral mesh of the built-in mesher and on
  its polyhedral variant (`mesh.generate.type = "poly"`) where the case
  allows; both must pass. Every case ships with gallery entry, case
  files and recipe (D18), and its numbers enter the verification manual
  (D68).

### 9.2 Code verification (exact solutions)

| Case | Quantity | Reference | Pass |
|---|---|---|---|
| Manufactured solutions for each model (construction after Eça, Hoekstra, Hay & Pelletier 2007): smooth k, ε or ω, ζ/φ, f/α and ν_t fields, all source terms of section 3 active, Dirichlet boundaries | L2 and L∞ errors of every variable | exact | p of the L2 error in [1.8, 2.2] on hexahedra, [1.6, 2.4] on polyhedra and tetrahedra; the same with the dual-number Jacobian test green |
| Homogeneous decay in a periodic box (no walls, uniform start) | k(t), ε(t) or ω(t) | exact ODE solutions: k-ε k = k₀(1 + (C_ε2 − 1)ε₀t/k₀)^{−1/(C_ε2 − 1)}; k-ω (SST with F₁ fixed at 1 or 0) ω = ω₀/(1 + βω₀t), k = k₀(1 + βω₀t)^{−β*/β} | time order of BDF2 in [1.9, 2.1]; error at the finest Δt below 10⁻⁶ relative; bounded BDF2 never leaves positivity |
| Sustaining terms | k, ω constant in time at the ambient values | exact | to round-off |
| Wall distance (D48) | d on planes, cylinders, spheres, and the PDE form | exact | exact (BVH) to round-off; PDE form within 1 % at the first cell layer |
| Wall functions as functions | u⁺(y⁺) of Spalding and the blended profile, T⁺ of Kader for Pr 0.01, 0.71, 7, 100, P-function | published formulas evaluated in high precision | relative error below 10⁻¹² (Newton steps converge) |
| Laminar heated channel, passive scalar: in 6a open (developed profile at the inlet, zero-gradient outlet), periodic with the heat-flux jump of D109 (1) once 6p is in | Nu | 140/17 = 8.235 | p in [1.8, 2.2], Nu within 0.5 % (D111) |

### 9.3 Zero-pressure-gradient flat plate

Set-up of the TMR flat plate verification case (2D, incompressible form of
its M = 0.2 case): Re = 5·10⁶ per unit length, plate from x = 0 to
x = 2 with a symmetry stretch ahead of it, the TMR's free-stream turbulence
values. The built-in mesher reproduces the TMR family of grids (the same
point counts and stretching, not the TMR files). Temperature as a passive
scalar with Pr = 0.71, isothermal plate.

| Quantity | Reference | Verification | Validation band B_model |
|---|---|---|---|
| C_f at x = 0.97 (SST, low-re treatment) | TMR finest-grid value of SST, stored with its retrieval date in the case's reference file | p in [1.5, 2.5], GCI_fine ≤ 1 % | ≤ 1 % (code-to-code: SST 2003 and the 1994 SST differ by S against Ω in ν_t, which is the same in a boundary layer to well below 1 %) |
| C_f(Re_x), 10⁶ ≤ Re_x ≤ 10⁷, every model with its default treatment | White (2006): C_f = 0.455/ln²(0.06 Re_x) | GCI_fine ≤ 2 % at x = 0.97 | ±5 % |
| u⁺(y⁺) at Re_θ ≈ 10⁴ (y⁺ ≈ 1 meshes) | Spalding's law for 1 ≤ y⁺ ≤ 30, log law for 30 ≤ y⁺ ≤ 300 | – | ±3 % (SST, ζ-f, EB) |
| St_x(Re_x), 5·10⁵ ≤ Re_x ≤ 5·10⁶ | St = 0.0287 Re_x^{−0.2} Pr^{−0.4} (Kays, Crawford & Weigand 2005) | GCI_fine ≤ 2 % | ±10 % |
| y⁺ insensitivity: first-cell y⁺ of 0.5, 1, 2, 5, 10, 30, 60, 100 at x = 0.97 | the model's own y⁺ ≈ 1 result (low-re) | – | C_f within ±2 % (SST automatic), ±3 % (ζ-f and EB compound), ±5 % (k-ε compound, and scalable for y⁺ ≥ 30); St within ±5 % (Kader) |
| Smooth limits against hard limits (2.6) | the hard-limit result | – | C_f within 0.1 % |

### 9.4 Fully developed channel flow

Periodic plane channel driven by a mean pressure gradient that fixes
Re_τ (D29 *2026-11-15*), one cell column (2D, symmetric or full height),
refined in y only (four levels, r = 2); for the low-re models the finest
y⁺ of the first cell is 0.25. Re_τ = 180, 395 and 590 against the DNS of
Moser, Kim & Mansour (1999). Passive temperature with uniform heat flux
(the jump of D109 (1)), Pr = 0.71, against the DNS of Kasagi, Tomita &
Kuroda (1992, Re_τ = 150) and Kawamura et al. (1998, Re_τ = 180; Re_τ = 395 from the later DNS of the same group, chosen at implementation).

| Quantity | Models | Verification | Validation band |
|---|---|---|---|
| Bulk velocity U_b⁺ (equivalently C_f) | all | GCI_fine ≤ 0.5 % | ±3 % (SST, ζ-f, EB, realizable and standard k-ε with compound at y⁺ ≈ 1); ±5 % (k-ε family with scalable wall functions at y⁺ ≈ 30, Re_τ 395 and 590 only) |
| Centreline velocity U_c⁺ | all | GCI_fine ≤ 0.5 % | ±4 % |
| u⁺(y⁺), 5 ≤ y⁺ ≤ Re_τ | SST, ζ-f, EB | – | maximum deviation ≤ 5 % |
| k⁺ peak value and position; v²⁺ = ζk⁺ or φk⁺ near-wall profile | ζ-f, EB | – | peak position within ±20 %; reported against the model papers' own channel curves (the paper result is the acceptance target for EB at Re_τ = 395, section 3.4) |
| Nu (heat-flux jump), constant Pr_t and Kays–Crawford | all, Re_τ 150/180/395 | GCI_fine ≤ 1 % | ±7 % |
| Wall-function insensitivity: first-cell y⁺ from 1 to 100 | compound and automatic | – | wall shear within ±3 % of the y⁺ ≈ 1 result (stage-6 row of the plan) |

### 9.5 Backward-facing step (Driver & Seegmiller 1985)

Expansion ratio 1.125, Re_H ≈ 36 000 on the step height H, inflow boundary
layer matched to the experiment's (δ ≈ 1.5 H at x/H = −4 ◇), as in the TMR's
2D backward-facing step case (its grid family reproduced by the built-in
mesher, three to four levels). Reattachment length x_R/H = 6.26 ± 0.10 from
the experiment.

| Quantity | Models | Verification | Validation band |
|---|---|---|---|
| x_R/H | SST, ζ-f, EB | GCI_fine ≤ 2 %, p in [1, 3] | ±10 % of 6.26 |
| x_R/H | standard and realizable k-ε | GCI_fine ≤ 2 % | ±5 % of a published result of the same model on a comparable grid (selected and recorded in the case before the run; the standard k-ε is known to reattach early) |
| C_p and C_f along the lower wall at the experiment's stations | all | – | reported against the experiment; RMS difference of C_p ≤ 0.02 for SST, ζ-f, EB |
| Velocity profiles at x/H = 1, 4, 6, 10 | all | – | reported (plots in the gallery) |

### 9.6 Axisymmetric impinging jet, H/D = 2

Fully developed pipe jet (profile from a periodic pipe run, D111, fed in as a
table inlet with k and ε or ω columns, section 7.1) on a flat plate at
H/D = 2, axisymmetric on a wedge mesh (D29 *Part 5b*). Heat transfer:
Re = 23 750 with uniform wall heat flux, Baughn & Shimizu (1989) and Baughn
et al. (1991); velocity: Re = 23 000, Cooper, Jackson, Launder & Liao (1993).
Temperature as a passive scalar with Pr = 0.71 and uniform flux. The
digitised reference data (Nu(r/D); mean velocity and turbulence profiles at
the published r/D stations) are stored in the case with their source,
figure or table, and retrieval note.

| Quantity | Models | Verification | Validation band |
|---|---|---|---|
| Stagnation Nusselt number Nu₀ | ζ-f, EB | GCI_fine ≤ 3 % | ±15 % |
| Nu₀ | SST with `kato-launder` or `durbin` | GCI_fine ≤ 3 % | ±20 % |
| Secondary maximum of Nu(r/D) | ζ-f, EB | – | present, located at r/D within ±0.5 of the experiment's |
| Nu(r/D), 0 ≤ r/D ≤ 4 | ζ-f, EB, SST with limiter | – | RMS difference ≤ 15 % |
| Stagnation-point anomaly | standard k-ε without and with Kato–Launder; SST `none` against `menter`, `kato-launder`, `durbin` | – | without limiter Nu₀ over-predicted by more than 50 % (the anomaly is reproduced); the limiter removes at least 70 % of the excess k at the stagnation point. This tests the limiters, not the models |
| Wall-jet velocity profiles at r/D = 0.5, 1, 2, 2.5, 3 | all | – | reported |

### 9.7 Cases kept from the earlier stage-6 scope

Unchanged from the stage plan: rough-pipe friction factor against
Colebrook (periodic pipe, fully rough and transitional, ±5 %), the study
over Re and k_s/D that reproduces the Moody chart (D62), laminar
Taylor–Couette flow in a rotating frame (D55) against the analytic
solution, the turbulent heated pipe as a passive scalar against Gnielinski
(1976) (±10 %, 10⁴ ≤ Re ≤ 10⁶, D111) with each model and both Pr_t
options, and the CPU–GPU and thread-count parity of section 6.6.

### 9.8 The SST options of 3.6 (D170)

Bands fixed before the runs, as everywhere in this section; the ◇ bands
are proposals to be confirmed with the papers.

| Case | Quantity | Reference | Band |
|---|---|---|---|
| Flat plate and channel at Re_τ = 395 with `curvature_correction` | C_f, u⁺ profile | SST without the option (pure shear: f_r1 = 1) | relative difference ≤ 1e-6 |
| Channel in spanwise system rotation, Re_τ = 194 (rotating zone, D55; periodic, 6p) | wall shear stress on the unstable and the stable side, mean-velocity profile | DNS of Kristoffersen & Andersson (1993) ◇ | ratio of the two wall shear stresses within ±15 % ◇ |
| Fully developed square duct, Re_b ≈ 3 500 (periodic, 6p) with `qcr = "qcr2013"` | largest secondary velocity, its position | DNS of Pinelli, Uhlmann, Sekimoto & Kawahara (2010) ◇ | within ±30 % and the corner vortex pair present; linear SST gives none (a check that the option acts) |
| Manufactured solution with the corrections switched on | observed order | — | 1.8 – 2.2 (the corrections must not lower the order) |


## 10. Implementation plan

Pull requests in order, each with tests, gallery entry, case files, recipe
(D18), changelog fragment and the verification-manual update; each starts
only on the maintainer's go. Branch names from the task (project rule).

| PR | Content | Main files | Acceptance |
|---|---|---|---|
| 6a | Generic scalar transport (D116), passive scalars and their boundary sub-tables (D51), bounded BDF2 (D28 (1)), smooth-limit helpers (2.6), exact wall distance with the BVH (D48, D108 (1)) and the derived field `wall-distance` (done, PR #63) | `foehncraft-kernels/src/scalar.rs`, `foehncraft-mesh/src/query/bvh.rs`, `foehncraft-flow/src/scalar.rs`, `foehncraft-core/src/math.rs` | 9.2 laminar heated channel (developed profile at the inlet, uniform wall flux, zero-gradient outlet; Nu within 0.5 %), wall distance, scalar MMS on hex, poly and tet cells, scalar time order |
| 6p | Periodic boundaries for the flow and the scalars (translational pairs; the affine temperature map of D109–D111 with the energy stage) | `foehncraft-mesh` periodic pairing, `foehncraft-flow/src/bc.rs` | periodic Poiseuille and heated channel equal to the open ones; needed by the 9.4 channel at Re_τ 180/395/590 |
| 6b | k-ω SST 2003 with the three production options, low-re and automatic treatments, rough walls (Aupoix), inlet turbulence, TURB codes, log columns, p̂; the optional PDE wall distance with the key `turbulence.wall_distance` (the exact distance of 6a is the default); the invariant features of D166 as named expression variables (same definitions as the D106 export) | `foehncraft-kernels/src/turb/{mod,invariants,sst,wall}.rs`, `foehncraft-flow/src/turbulence/{mod,wall,inlet}.rs` | 9.2 MMS and decay (SST), 9.3 flat plate, 9.4 channel (SST) |
| 6c | Standard and realizable k-ε, scalable and standard wall functions, compound treatment for the k-ε family, Cebeci–Bradshaw roughness | `turb/keps.rs` | 9.2, 9.3, 9.4 for the k-ε family, rough pipe and Moody study |
| 6d | k-ε-ζ-f with the compound treatment (Popovac & Hanjalić 2007; Popovac 2023 refinements as decided then), elliptic solver path (CG + AMG), the AMG preconditioner for the scalar equations (D36; 6a solves them with Jacobi-preconditioned BiCGStab, enough for the laminar cases) | `turb/zeta_f.rs`, `foehncraft-flow/src/turbulence/elliptic.rs` | 9.2, 9.3, 9.4 incl. near-wall profiles |
| 6e | EB k-ε (D133) | `turb/eb.rs` | 9.4 at Re_τ = 395 against the paper first, then 9.2–9.4 |
| 6f | Thermal wall treatment and Pr_t (D135) on passive scalars | `turb/thermal.rs` | 9.2 wall-function functions, St on the flat plate, channel Nu, heated pipe against Gnielinski |
| 6g | Validation cases: backward-facing step, impinging jet; MRF (D55) with Taylor–Couette; statistics (D71) and derived fields (D72) of the turbulence models; GPU parity notebook; GCI reports | cases, gallery, `tools/notebooks/stage6_colab.ipynb` | 9.5, 9.6, 9.7, 6.6 |
| 6h | SST options of D170: curvature and rotation correction (Smirnov & Menter 2009), QCR2000 and QCR2013 with the explicit-stress path in the momentum assembly that D171 reuses | `turb/sst.rs`, `turb/invariants.rs`, `foehncraft-flow/src/turbulence/stress.rs` | 9.8 |
| release | v0.6: version bump, changelog summary, README, verification manual, Release workflow (dry run first) | | all of section 9 green or findings recorded |

Order rationale: 6a gives every later model its transport, wall distance
and solver path; 6p brings the periodic boundaries that the channel cases
of 9.4 need, so it comes before the first model is verified there; SST first because it is the general-purpose default and
the flat plate and channel can be checked against the TMR at once; ζ-f
before EB because EB reuses its elliptic path; the thermal wall function
after the models because it needs u* from each of them; the SST options
(6h) after the validation cases, because their checks need the rotating
zone and the periodic duct. The EARSM of D171 is not in v0.6.


## 11. Questions for the maintainer

- **Q1 (closed 2026-11-17):** EB k-ε joins ζ-f in stage 6 (D133).
- **Q2:** Lag EB k-ε as P11 (possible later), decided after the EB
  verification (recommended), or in stage 6?
- **Q3:** The model bands of section 9 are a proposal; change any before
  the first run (they are fixed once a run has been made).
- **Q4 (open, not confirmed):** Add the k-ω model of Wilcox (2006) as a
  cheap variant on the SST code path (its stress limiter and cross
  diffusion; one more closure of the same kernels)? Not part of the plan
  until the maintainer confirms it.


## 12. References

- B. Aupoix, *Roughness corrections for the k–ω shear stress transport
  model: status and proposals*, J. Fluids Eng. 137, 021202 (2015)
- J. W. Baughn, S. Shimizu, *Heat transfer measurements from a surface
  with uniform heat flux and an impinging jet*, J. Heat Transfer 111,
  1096–1098 (1989)
- J. W. Baughn, A. E. Hechanova, X. Yan, *An experimental study of
  entrainment effects on the heat transfer from a flat surface to a heated
  circular impinging jet*, J. Heat Transfer 113, 1023–1025 (1991)
- F. Billard, D. Laurence, *A robust k-ε-v²/k elliptic blending turbulence
  model applied to near-wall, separated and buoyant flows*, Int. J. Heat
  Fluid Flow 33, 45–58 (2012)
- T. Cebeci, P. Bradshaw, *Momentum transfer in boundary layers*,
  Hemisphere (1977)
- I. B. Celik et al., *Procedure for estimation and reporting of
  uncertainty due to discretization in CFD applications*, J. Fluids Eng.
  130, 078001 (2008)
- D. Cooper, D. C. Jackson, B. E. Launder, G. X. Liao, *Impinging jet
  studies for turbulence model assessment — I. Flow-field experiments*,
  Int. J. Heat Mass Transfer 36, 2675–2684 (1993)
- D. M. Driver, H. L. Seegmiller, *Features of a reattaching turbulent
  shear layer in divergent channel flow*, AIAA J. 23, 163–171 (1985)
- P. A. Durbin, *On the k-ε stagnation point anomaly*, Int. J. Heat Fluid
  Flow 17, 89–90 (1996)
- L. Eça, M. Hoekstra, A. Hay, D. Pelletier, *On the construction of
  manufactured solutions for one and two-equation eddy-viscosity models*,
  Int. J. Numer. Methods Fluids 54, 119–154 (2007)
- V. Gnielinski, *New equations for heat and mass transfer in turbulent
  pipe and channel flow*, Int. Chem. Eng. 16, 359–368 (1976)
- G. Grotjans, F. R. Menter, *Wall functions for general application CFD
  codes*, ECCOMAS 98, 1112–1117 (1998)
- K. Hanjalić, M. Popovac, M. Hadžiabdić, *A robust near-wall
  elliptic-relaxation eddy-viscosity turbulence model for CFD*, Int. J. Heat
  Fluid Flow 25, 1047–1051 (2004)
- A. Hellsten, *New advanced k-ω turbulence model for high-lift
  aerodynamics*, AIAA J. 43, 1857–1869 (2005)
- C. L. V. Jayatilleke, *The influence of Prandtl number and surface
  roughness on the resistance of the laminar sub-layer to momentum and
  heat transfer*, Prog. Heat Mass Transfer 1, 193–329 (1969)
- B. A. Kader, *Temperature and concentration profiles in fully turbulent
  boundary layers*, Int. J. Heat Mass Transfer 24, 1541–1544 (1981)
- N. Kasagi, Y. Tomita, A. Kuroda, *Direct numerical simulation of passive
  scalar field in a turbulent channel flow*, J. Heat Transfer 114, 598–606
  (1992)
- M. Kato, B. E. Launder, *The modelling of turbulent flow around stationary
  and vibrating square cylinders*, 9th Symp. Turbulent Shear Flows, Kyoto,
  10-4 (1993)
- H. Kawamura, K. Ohsaka, H. Abe, K. Yamamoto, *DNS of turbulent heat
  transfer in channel flow with low to medium-high Prandtl number fluid*,
  Int. J. Heat Fluid Flow 19, 482–491 (1998)
- W. M. Kays, *Turbulent Prandtl number — where are we?*, J. Heat Transfer
  116, 284–295 (1994)
- W. M. Kays, M. E. Crawford, B. Weigand, *Convective heat and mass
  transfer*, 4th ed., McGraw-Hill (2005)
- R. Kristoffersen, H. I. Andersson, *Direct simulations of low-Reynolds-number
  turbulent flow in a rotating channel*, J. Fluid Mech. 256, 163–197 (1993)
- S. Lardeau, F. Billard, *Development of an elliptic-blending lag model
  for industrial applications*, AIAA 2016-1600 (2016)
- B. E. Launder, D. B. Spalding, *The numerical computation of turbulent
  flows*, Comput. Methods Appl. Mech. Eng. 3, 269–289 (1974)
- M. Mani, D. A. Babcock, C. M. Winkler, P. R. Spalart, *Predictions of a
  supersonic turbulent flow in a square duct*, AIAA 2013-0860 (2013) ◇
- F. R. Menter, *Two-equation eddy-viscosity turbulence models for
  engineering applications*, AIAA J. 32, 1598–1605 (1994)
- F. R. Menter, M. Kuntz, R. Langtry, *Ten years of industrial experience
  with the SST turbulence model*, Turbulence, Heat and Mass Transfer 4,
  625–632 (2003)
- R. D. Moser, J. Kim, N. N. Mansour, *Direct numerical simulation of
  turbulent channel flow up to Re_τ = 590*, Phys. Fluids 11, 943–945 (1999)
- NASA Langley Turbulence Modeling Resource (TMR): SST and SST-2003 pages,
  2D zero-pressure-gradient flat plate and 2D backward-facing step cases
  (values retrieved at implementation, retrieval date stored with the data)
- S. V. Patankar, *Numerical heat transfer and fluid flow*, Hemisphere (1980)
- A. Pinelli, M. Uhlmann, A. Sekimoto, G. Kawahara, *Reynolds number
  dependence of mean flow structure in square duct turbulence*, J. Fluid
  Mech. 644, 107–122 (2010)
- M. Popovac, K. Hanjalić, *Compound wall treatment for RANS computation of
  complex turbulent flows and heat transfer*, Flow Turbul. Combust. 78,
  177–202 (2007)
- M. Popovac, Energies 16, 3685 (2023) (refinements of the compound wall
  treatment; read before PR 6d)
- T.-H. Shih, W. W. Liou, A. Shabbir, Z. Yang, J. Zhu, *A new k-ε eddy
  viscosity model for high Reynolds number turbulent flows*, Comput. Fluids
  24, 227–238 (1995)
- P. E. Smirnov, F. R. Menter, *Sensitization of the SST turbulence model
  to rotation and curvature by applying the Spalart–Shur correction term*,
  J. Turbomach. 131, 041010 (2009)
- P. R. Spalart, *Strategies for turbulence modelling and simulations*,
  Int. J. Heat Fluid Flow 21, 252–263 (2000)
- P. R. Spalart, C. L. Rumsey, *Effective inflow conditions for turbulence
  models in aerodynamic calculations*, AIAA J. 45, 2544–2553 (2007)
- P. R. Spalart, M. Shur, *On the sensitization of turbulence models to
  rotation and curvature*, Aerosp. Sci. Technol. 1, 297–302 (1997)
- D. B. Spalding, *A single formula for the "law of the wall"*, J. Appl.
  Mech. 28, 455–458 (1961)
- S. Wallin, A. V. Johansson, *An explicit algebraic Reynolds stress model
  for incompressible and compressible turbulent flows*, J. Fluid Mech. 403,
  89–132 (2000)
- F. M. White, *Viscous fluid flow*, 3rd ed., McGraw-Hill (2006)
- D. C. Wilcox, *Turbulence modeling for CFD*, 3rd ed., DCW Industries
  (2006) (Q4 only)

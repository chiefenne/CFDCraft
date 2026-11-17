# Step 20 – Data-driven models: assessment for Foehncraft (end of 2026)

Requirements step 20 (2026-11-17), asked for by the maintainer: "what could
we do with the latest data-driven models at the end of 2026". This page is an
**assessment with proposals only**. Nothing on it is decided; the proposed
items (D165–D169, L37–L38, P13–P15) become decisions only when the maintainer
says so, and no code is written before the maintainer's go (project rule).

It extends, and does not repeat, the physics-based learning assessment of
2026-11-16 that gave D104 (model-correction fields), D105 (stepping API),
D106 (training-data export) and L26 (linear stability), and it builds on the
stage-6 specification ([`stage6-turbulence.md`](stage6-turbulence.md),
D129–D139) and the stage-7 specification
([`stage7-energy-buoyancy-cht.md`](stage7-energy-buoyancy-cht.md),
D150–D164). Clean-room rule (D21): methods from the cited papers only; no
code of other projects is read or copied for this. Other codes are not named;
papers are cited by authors.

Contents

1. [Short answer](#1-short-answer)
2. [State of the field at the end of 2026](#2-state-of-the-field-at-the-end-of-2026)
3. [Real or hype, judged for industrial accuracy](#3-real-or-hype-judged-for-industrial-accuracy)
4. [What Foehncraft already has](#4-what-foehncraft-already-has)
5. [Reproducibility, the reproducible mode (D93) and GPU kernels](#5-reproducibility-the-reproducible-mode-d93-and-gpu-kernels)
6. [Acceptance of learned components](#6-acceptance-of-learned-components)
7. [Licences and data](#7-licences-and-data)
8. [Proposed items and their stages](#8-proposed-items-and-their-stages)
9. [Questions for the maintainer](#9-questions-for-the-maintainer)
10. [References](#10-references)

## 1. Short answer

Prepare now, build after milestone 1, and keep every learned model outside
the solver: it comes in as a field (D104) or, new here, as a **closed-form
expression** found by symbolic regression and compiled like any other
expression (D53). That route is the one with real, published gains on
separated internal flows; it keeps bitwise reproducibility and the GPU kernels
untouched, and it needs no neural-network inference inside a run, which the
maintainer excluded on 2026-11-16.

Ranking after the literature check (changed from the first pass in the
project chat):

| # | Area | Verdict for Foehncraft | When |
|---|---|---|---|
| 1 | **Learned closure corrections on SST and ζ-f**, by conditioned field inversion and symbolic regression | Real gains on separation, and the only learned item that can meet industrial acceptance today. Plan it | Expression features in stage 6 or 8 (D166); campaign after v1.0 (L37) |
| 2 | **Neural-operator surrogates for design sweeps** in one product family | Real for ranking variants inside the training family, not for final numbers. Foehncraft is the data generator and the referee | Data side via D106 and D168; surrogate outside the solver (P14) |
| 3 | **Learned solver components** (AMG smoother, warm start) | Real but small: published wall-clock gains of 4–37 % for a learned smoother, and warm starts help little where the last step or iterate is already a good guess. Watch, do not build | P13 with a measured threshold |
| 4 | **Pretrained flow foundation models** | Promising for few-shot transfer, far from engineering accuracy zero-shot. Hype for industrial use in 2026 | Watch; only as a consumer of our datasets |
| 5 | **Physics-informed networks as solvers** | Not competitive with finite volumes for industrial RANS. Unchanged verdict | Not planned |
| 6 | Learned wall models and SGS models for LES | Active research with first a posteriori successes; belongs to LES (L5) | P15, with L5 |

The learned preconditioner moved down from rank 2: the published gains are
modest, every learned smoother would have to be frozen input data to keep
D31 and D93, and our pressure solve already starts from the previous time
step or outer iteration.

## 2. State of the field at the end of 2026

### 2.1 Learned turbulence closures (RANS)

- **Field inversion and machine learning** (Parish & Duraisamy 2016; Singh,
  Medida & Duraisamy 2017): a per-cell multiplier β on a production term is
  found by inverting reference data with an adjoint; a model β(features) is
  then trained. This is the D104 route. The main weakness in the literature
  is generalisation: corrections learned on separated flows spoil attached
  boundary layers (Duraisamy 2021; Rumsey et al. 2022; Spalart 2023).
- **Conditioned field inversion with symbolic regression** (Wu, Zhang & Zhang,
  AIAA Journal 2025): β is multiplied by a shielding function that is zero
  inside attached boundary layers (the delayed-DES f_d form), the inverted
  fields from two separated cases (NASA hump, curved backward-facing step)
  are turned into closed-form expressions by symbolic regression, and the
  model is tested on periodic hills, a two-element airfoil, a notchback car,
  the Ahmed body, the flat plate and a NACA 0012. The conditioned model
  improves the separated cases over SST 2003 **and keeps the flat-plate skin
  friction and the log layer**, which the unconditioned model loses. This is
  the most useful single result for us: SST 2003 is our baseline (D130), the
  output is an expression (D53), and the flat plate is in our suite (D139).
- **Sparse regression of algebraic stress corrections** (Weatheritt &
  Sandberg 2016; Schmelzer, Dwight & Cinnella 2020): explicit tensor
  polynomials in invariant bases, again closed form.
- **Ensemble (derivative-free) inversion** (Zhang, Xiao, Luo & He, JFM 2022):
  ensemble Kalman methods learn the model from indirect observations (velocity
  profiles, wall quantities) without an adjoint. Useful before our adjoint
  (L17) exists.
- **Tensor-basis and equivariant networks** (Ling, Kurzawski & Templeton
  2016; Miller, Kommalapati, Moser & Koumoutsakos 2025): equivariant networks
  for the rapid pressure–strain with exact linear tensor constraints reach
  errors below 0.1 % in rapid distortion, but only in that limit. Injecting a
  learned Reynolds stress directly is ill-conditioned (Wu, Xiao, Sun & Wang
  2019): small stress errors give large velocity errors, so the learned part
  must enter implicitly (eddy-viscosity part) or as a correction to a
  transported quantity.

### 2.2 Learned linear-solver components

- Learned multigrid prolongation and smoothers (Greenfeld et al. 2019; Luz
  et al. 2020) and, for an unstructured-mesh pressure AMG, a graph-network
  smoother that predicts polynomial coefficients of a sparse approximate
  inverse (Chillón, Lidtke, Doan & Font, arXiv 2026): the solver stays
  linear, V-cycles drop, wall-clock gains are **4 % to 37 %**, and the model
  transfers to meshes 128 times larger than the training meshes.
- Neural-operator warm starts for Krylov solvers (Eshaghi et al., arXiv 2025):
  the classical solver, discretisation and preconditioner are kept, so the
  final accuracy is the solver's; iteration counts fall 2–10× **against a
  cold start**. In a transient PISO or a SIMPLE loop the initial guess is the
  last solution, so most of that gain is already in hand.

### 2.3 Surrogates and foundation models

- Mesh-based surrogates (Pfaff et al. 2021), neural operators (Li et al.
  2021), transformers on point clouds (Wu et al. 2024; Alkin et al. 2025) now
  train on tens of millions of cells per sample. On an automotive set of about
  2000 cases with 45 M volume cells, a transformer surrogate trains in about
  14 GPU hours and predicts a full case in about 34 s (Alkin et al. 2025).
  Open automotive data sets exist (Ashton et al. 2024; Elrefaie et al. 2024).
- Cross-domain pretraining (Zhou, Barati Farimani, Ho & Morel, arXiv 2026):
  a model pretrained on six steady CFD data sets needs about 8× fewer samples
  for a new family, but its **zero-shot** normalised RMSE on a held-out body
  is 0.31 (0.09 after eight samples), and on its own training family it is
  worse than a specialist (0.134 against 0.048).
- A benchmark of foundation models for flow around complex geometries
  (Rabeh et al., Communications Engineering 2025): pretrained models beat
  neural operators by about one order of magnitude in MSE when data are
  scarce, but **all** models fail on Reynolds numbers outside the training
  range. A newer benchmark asks for well-posed test problems for aerodynamic
  foundation models (AeroZoo, arXiv 2611.17152, November 2026; read at title
  and abstract level only).
- Weak baselines: many published speed-ups of machine learning for fluid PDEs
  vanish against a well-tuned classical solver at equal accuracy (McGreivy &
  Hammett, Nature Machine Intelligence 2024). Every claim we make must be
  measured against our own solver at equal accuracy.

### 2.4 Lesson from weather prediction

The most successful hybrid systems of 2024–2026 keep a conservative numerical
core and learn only the parametrised physics (Kochkov et al., Nature 2024).
That is the D104 and D166 design: the finite-volume core stays exact and
conservative; the closure is what is learned.

## 3. Real or hype, judged for industrial accuracy

Industrial accuracy here means: a validated quantity (pressure loss,
separation length, heat transfer coefficient, force) within the model band
of D139 plus U_num plus U_ref, on a case not used for training, with grid
convergence shown.

| Claim | Real? | Why |
|---|---|---|
| Learned correction improves separated RANS flows | **Real**, within the flow class trained on | Conditioned FI + SR: better separation, flat plate kept (Wu, Zhang & Zhang 2025). Still RANS: no gain where the eddy-viscosity form itself is wrong (strong secondary flows of the second kind need anisotropy) |
| One learned model fixes RANS for all flows | Hype | No published model passes a broad suite; the "do no harm" cases are where unconditioned models fail |
| Surrogate replaces CFD in design | Partly real | For ranking variants inside the trained family, yes; absolute values and anything out of family, no (Rabeh et al. 2025) |
| Foundation model gives CFD answers zero-shot | Hype in 2026 | nRMSE ≈ 0.3 zero-shot (Zhou et al. 2026); fine for previews, not for decisions |
| Learned preconditioner makes the solver much faster | Overstated | 4–37 % wall-clock (Chillón et al. 2026); large factors only against weak baselines (McGreivy & Hammett 2024) |
| PINNs replace finite-volume solvers | Hype | Not competitive in accuracy or cost for 100 M cell RANS; the research groups that started them moved to hybrids |

## 4. What Foehncraft already has

- **Fields in, nothing linked in:** D104 correction fields (β on production,
  additive source) on every transport equation of D116, with the hook already
  named in the stage-6 framework (D129).
- **Expressions compiled to kernels:** D53 compiles TOML or Python expression
  objects to CPU and GPU kernels; no user code is linked. A symbolic-regression
  model is exactly such an expression.
- **Differentiable kernels:** D59 (pure kernels, dual numbers, smooth variants
  of min/max of D129, Aᵀ) and the steady adjoint L17 give dJ/dβ per cell.
- **Data side:** D105 stepping API, D106 dataset export with invariant
  features, D62 studies, D67 level sweeps (fidelity levels), D85 provenance,
  D31/D93 bitwise reproducibility.

Gaps found:

1. The expression language of D53 has space, time and solution fields, but no
   **named invariant features** (strain and rotation invariants, wall distance,
   the shielding function f_d, turbulence time scale ratios). A learned
   expression needs exactly these.
2. No **acceptance rule** for learned components: nothing says when a learned
   model may be used in a validation result or shipped.
3. No rule for **licences of training data and models**.
4. Inversion needs either the adjoint (L17, later) or an **ensemble method**
   that only needs many runs; the latter is not yet a study type.
5. Learned SGS or time-stepping corrections need a **transient** adjoint; L17
   names it only as "later".

## 5. Reproducibility, the reproducible mode (D93) and GPU kernels

| Entry route | D31 / D93 | GPU | Note |
|---|---|---|---|
| Field (D104), read from a file | Holds: input data | Unchanged | File hash in the provenance (D85) |
| Expression (D53, D166) | Holds: compiled like any expression, fixed evaluation order per cell | Same kernel path as other expressions | Smooth form of any min/max for the adjoint (D129) |
| Surrogate as initial state (D52) | Holds for the run after start; the start state is input data | Unchanged | Final result is the solver's, the surrogate only shortens the path |
| Learned AMG smoother (P13) | Only if its coefficients are computed once on the host and frozen as input data | Sparse polynomial: ordinary SpMV | Training not part of the run |
| Neural inference inside a run | Excluded (owner 2026-11-16) | — | Symbolic regression gives the benefit without it |

Training itself happens outside the solver and need not be bitwise
reproducible; what must be reproducible is the **model artefact** (expression
text or field file, with its hash) and therefore every run that uses it.
Stochastic training (ensemble inversion, random initialisation) takes its
seeds from the counter-based generator of D47 and records them.

## 6. Acceptance of learned components

Proposed as D165. A learned component (correction field, learned expression,
learned solver component, surrogate used for a result) is accepted for a
validation result or for shipping only if all of the following hold, with
bands fixed **before** the runs, as D139 does:

1. **Do no harm.** The verification and "attached-flow" validation cases of
   stage 6 (flat plate, channel at Re_τ 180/395/590, D139) stay inside their
   bands with the model switched on. For a learned solver component: the
   converged solution is bitwise equal to the one without it (only the path
   changes), checked in D93 mode.
2. **Improvement beyond uncertainty.** On the target case the comparison
   error E = S − D must shrink by more than the validation uncertainty
   u_val = √(u_num² + u_input² + u_D²) (ASME V&V 20-2009), with u_num from the
   GCI of D67 (Celik et al. 2008).
3. **Held-out case.** At least one case of the same flow class that was not
   used in training or inversion, and one of a different class.
4. **Grid convergence with the model on.** The model must not break
   convergence: three grids give a monotone or oscillatory-converged verdict
   and a GCI_fine below the D139 limit. Features that depend on the grid
   (cell size, y⁺ of the first cell) are not allowed as model inputs.
5. **Realisability and invariance.** Inputs are Galilean- and
   rotation-invariant and dimensionless; ν_t ≥ 0, k ≥ 0 and the Reynolds
   stress realisable after the correction; the model is bounded outside its
   training range (clip to the training feature ranges, reported as a health
   number D69).
6. **Model card.** Training data (with hashes and licences, §7), baseline model
   and version, feature definitions, training ranges, the acceptance results
   above, and the Foehncraft version they were obtained with.
7. **Surrogates** never fill a result summary; a surrogate value carries
   `source = "surrogate"` and an uncertainty (as proposed in the 2026-11-16
   assessment for the agent interface).

## 7. Licences and data

- **Reference data for inversion.** Channel DNS (Moser, Kim & Mansour 1999;
  Lee & Moser 2015), the curated data set for data-driven turbulence modelling
  (McConkey, Yee & Lien 2021), ERCOFTAC and NASA turbulence-modelling cases.
  Their terms differ and must be recorded per file.
- **Large surrogate data sets** are published under different licences, some
  share-alike, some non-commercial. A model trained on non-commercial data
  cannot be shipped with the MIT-licensed solver for unrestricted
  use. Proposed rule (D168): every imported data set and every model artefact
  carries an SPDX licence identifier and a source hash; models shipped with
  Foehncraft are trained only on data whose licence allows redistribution of
  derived models; licences are checked when the data are imported, not later.
- **Own data first.** Foehncraft can generate its own training data with D62
  and D106 under its own licence (fidelity levels from D67). That keeps the
  data clean and is the main thing an open solver can offer this field.
- **No proprietary data** in the repository. Users may train on their own
  data; the solver only reads the resulting fields or expressions.

## 8. Proposed items and their stages

Nothing below is decided. Proposals:

| ID | Proposed | Content | Stage |
|---|---|---|---|
| D165 | Acceptance of learned components | Section 6: do-no-harm on the D139 suite, improvement beyond u_val, held-out cases, grid convergence with the model on, realisability and invariance, model card, surrogate labelling | Recorded now, applies from the first learned model |
| D166 | Invariant features in expressions | D53 expressions may use named, dimensionless, invariant features: the five Pope invariants of the normalised strain and rotation tensors, wall distance and d⁺ based on k (not on the cell), ratio of turbulent to mean time scale, ν_t/ν, the shielding function f_d (delayed DES form), Q criterion, pressure-gradient and curvature indicators; the same definitions as the features of D106, so a model trained on exported data runs unchanged. Smooth forms for the adjoint (D129). With D104 this is the whole interface for learned closures | Stage 6 if cheap (the expression compiler and D106 feature code are shared), otherwise stage 8 |
| D167 | Ensemble inversion as a study type | `foehncraft.study` gets an ensemble Kalman inversion (Zhang et al. 2022) of D104 fields or of expression constants, using D105 and D62; seeds from D47; works without the adjoint and is replaced or complemented by adjoint field inversion when L17 exists | After D105 is complete (stage 6 or 8) |
| D168 | Data and model licences | Section 7: SPDX licence and hash on every imported data set and model artefact; shipped models only from redistributable data; recorded in the provenance (D85) and the model card (D165) | With D106 |
| D169 | Learned components never change results silently | A learned solver component (P13) or a surrogate start (D52) may change only the path, never the converged answer; the CI check runs the same case with and without it in D93 mode and compares bits | With P13 or the first surrogate start |
| L37 | Data-driven closure campaign | Conditioned field inversion and symbolic regression for SST 2003 and ζ-f: inversion first by D167, later by the adjoint (L17); training cases from open separated-flow data (hump, curved and backward-facing step, periodic hills); targets the asymmetric diffuser and the U-bend of D45 (curvature, P8); accepted by D165; result shipped as expressions with model cards, never as the default model | After v1.0 |
| L38 | Transient adjoint for learned unsteady corrections | Discrete adjoint of the time stepping with checkpointing (named "later" in L17), so that SGS and coarse-grid corrections can be trained through many steps (unrolled training, Franz et al. 2026) and used through D104 | After L17 and L5 |
| P13 | Learned AMG smoother or warm start | Only if (a) coefficients are computed on the host and frozen as input data, (b) D169 holds, and (c) it gives at least 30 % wall-clock on our own pressure benchmarks against our tuned aggregation AMG at equal tolerance | Watch item |
| P14 | Surrogate-assisted studies | Fit a surrogate (Gaussian process or an external neural surrogate) to D62 results for active learning and previews; the solver confirms every chosen point; agent "predict" tool with `source = "surrogate"` | After D62 and L18 |
| P15 | Learned wall and SGS models for LES | Through D104/D166 fields or expressions; only with L5 | With L5 |

Not proposed: neural-network inference inside the solver (excluded by the
maintainer on 2026-11-16; D166 covers the need); PINN solvers; any learned
model as a default.

Order if adopted: D165 and D168 now (text only); D166 with the stage-6
expression and D106 feature work, or in stage 8; D167 when D105 is complete;
L37 after v1.0 as the first learned-model campaign; P13–P15 stay watch items.

## 9. Questions for the maintainer

1. Adopt D165 (acceptance rules) and D168 (licences) now as decisions?
2. D166 (invariant features in expressions): stage 6 (with the D104 hook) or
   stage 8?
3. D167 (ensemble inversion): adopt, so that inversion works before the
   adjoint?
4. L37 (closure campaign after v1.0): first target the asymmetric diffuser or
   the U-bend, and SST 2003 first or ζ-f first?
5. P13: keep the learned smoother as a watch item with the 30 % threshold, or
   drop it?

## 10. References

- Alkin, B., Kurle, R., Serrano, L., Just, D. & Brandstetter, J. (2025). AB-UPT
  for automotive and aerospace applications. arXiv:2510.15808. (Model:
  arXiv:2502.09692.)
- ASME (2009). V&V 20-2009, Standard for verification and validation in
  computational fluid dynamics and heat transfer.
- Ashton, N. et al. (2024). DrivAerML: high-fidelity computational fluid
  dynamics dataset for road-car external aerodynamics. arXiv:2408.11969.
- Celik, I. B. et al. (2008). Procedure for estimation and reporting of
  uncertainty due to discretization in CFD applications. J. Fluids Eng. 130.
- Chillón, E., Lidtke, A. K., Doan, N. A. K. & Font, B. (2026). Acceleration of
  an algebraic multigrid pressure solver using graph neural networks.
  arXiv:2606.19251.
- Duraisamy, K. (2021). Perspectives on machine learning-augmented
  Reynolds-averaged and large eddy simulation models of turbulence. Phys. Rev.
  Fluids 6, 050504.
- Elrefaie, M. et al. (2024). DrivAerNet++: a large-scale multimodal car
  dataset. NeurIPS Datasets and Benchmarks.
- Eshaghi, M. S., Anitescu, C., Valizadeh, N., Wang, Y., Zhuang, X. &
  Rabczuk, T. (2025). NOWS: neural operator warm starts for accelerating
  iterative solvers. arXiv:2511.02481.
- Franz, Wei, Guastoni & Thuerey (2026). Differentiable GPU PISO solver and
  learned SGS model. J. Comput. Phys. 544 (cited in the 2026-11-16 assessment).
- Greenfeld, D., Galun, M., Basri, R., Yavneh, I. & Kimmel, R. (2019). Learning
  to optimize multigrid PDE solvers. ICML.
- Kochkov, D. et al. (2024). Neural general circulation models for weather and
  climate. Nature 632.
- Lee, M. & Moser, R. D. (2015). Direct numerical simulation of turbulent
  channel flow up to Re_τ ≈ 5200. J. Fluid Mech. 774.
- Li, Z. et al. (2021). Fourier neural operator for parametric PDEs. ICLR.
- Ling, J., Kurzawski, A. & Templeton, J. (2016). Reynolds averaged turbulence
  modelling using deep neural networks with embedded invariance. J. Fluid
  Mech. 807.
- Luz, I., Galun, M., Maron, H., Basri, R. & Yavneh, I. (2020). Learning
  algebraic multigrid using graph neural networks. ICML.
- McConkey, R., Yee, E. & Lien, F.-S. (2021). A curated dataset for data-driven
  turbulence modelling. Sci. Data 8, 255.
- McGreivy, N. & Hammett, A. (2024). Weak baselines and reporting biases lead
  to overoptimism in machine learning for fluid-related partial differential
  equations. Nat. Mach. Intell. 6.
- Miller, A., Kommalapati, S., Moser, R. & Koumoutsakos, P. (2025). Structure
  tensor RANS turbulence models with equivariant neural networks.
  arXiv:2511.09769.
- Moser, R. D., Kim, J. & Mansour, N. N. (1999). DNS of turbulent channel flow
  up to Re_τ = 590. Phys. Fluids 11.
- Parish, E. J. & Duraisamy, K. (2016). A paradigm for data-driven predictive
  modeling using field inversion and machine learning. J. Comput. Phys. 305.
- Pfaff, T., Fortunato, M., Sanchez-Gonzalez, A. & Battaglia, P. (2021).
  Learning mesh-based simulation with graph networks. ICLR.
- Rabeh, A. et al. (2025). Benchmarking scientific machine-learning approaches
  for flow prediction around complex geometries. Commun. Eng. 4, 182.
- Rumsey, C. L., Coleman, G. N. & Wang, L. (2022). In search of data-driven
  improvements to RANS models applied to separated flows. AIAA SciTech.
- Schmelzer, M., Dwight, R. P. & Cinnella, P. (2020). Discovery of algebraic
  Reynolds-stress models using sparse symbolic regression. Flow Turbul.
  Combust. 104.
- Singh, A. P., Medida, S. & Duraisamy, K. (2017). Machine-learning-augmented
  predictive modeling of turbulent separated flows over airfoils. AIAA J. 55.
- Spalart, P. R. (2023). An old-fashioned framework for machine learning in
  turbulence modeling. arXiv preprint.
- AeroZoo (2026). A benchmark of well-posed flow problems for foundation models
  of aerodynamics. arXiv:2611.17152.
- Weatheritt, J. & Sandberg, R. (2016). A novel evolutionary algorithm applied
  to algebraic modifications of the RANS stress–strain relationship. J.
  Comput. Phys. 325.
- Wu, C., Zhang, S. & Zhang, Y. (2025). Development of a generalizable
  data-driven turbulence model: conditioned field inversion and symbolic
  regression. AIAA J. doi:10.2514/1.J064416 (arXiv:2402.16355).
- Wu, H., Luo, H., Wang, H., Wang, J. & Long, M. (2024). Transolver: a fast
  transformer solver for PDEs on general geometries. ICML.
- Wu, J., Xiao, H., Sun, R. & Wang, Q. (2019). Reynolds-averaged Navier–Stokes
  equations with explicit data-driven Reynolds stress closure can be
  ill-conditioned. J. Fluid Mech. 869.
- Zhang, X.-L., Xiao, H., Luo, X. & He, G. (2022). Ensemble Kalman method for
  learning turbulence models from indirect observation data. J. Fluid Mech.
  949.
- Zhou, A., Barati Farimani, A., Ho, S. & Morel, R. (2026). Cross-domain
  pretraining for steady-state neural CFD surrogates. arXiv:2610.10398.

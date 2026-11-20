# Changelog

## Unreleased

Nothing yet.

## Summary

| Version | Stage | Highlights |
|---|---|---|
| 0.6 (2026-11-20) | 6 RANS turbulence, heat transfer, rotating machinery | Six RANS models: k-ω SST 2003 (with the curvature and rotation correction and QCR2000/QCR2013), standard and realizable k-ε, k-ε-ζ-f, k-ε-ζ-f₀ and the elliptic-blending k-ε; wall treatments `low-re`, `automatic`, `scalable`, `standard`, `compound`, `compound-2023` and rough walls; passive scalars with turbulent diffusivity and thermal wall functions (also for high Prandtl numbers); periodic boundaries; rotating machinery with the moving reference frame, sliding meshes (conformal stepping and imprinted interfaces of any shape of revolution, several rotors), torque, shaft power and machine performance; the coupled pressure–velocity algorithm with pseudo-time steps and multigrid for momentum and turbulence; SU2, OpenFOAM and CGNS mesh readers; derived turbulence fields, transient statistics, the grid convergence study of the models; the new log, case format 0.6 and `[[schedule]]`. **TMR flat plate C_f +0.37 % (SST) and within −2.6 to +2.4 % for first-cell y⁺ 0.24 to 85; channel at Re_τ 395 U_b⁺ −0.7 % (SST) and EB C_f 97.8–99.1 % of DNS; backward-facing step x_r/H +0.5 % (k-ε against the published k-ε) and +0.6 % (k-ε-ζ-f₀ against the experiment); impinging jet Nu₀ −0.2 % (SST); heated channel Nu within −2.6 to +1.0 % of DNS; GPU bitwise equal to the CPU for 14 turbulence cases on a Tesla T4; the coupled algorithm settles the step in about 70 instead of 600 iterations.** Gallery 26–129, recipes 36–129 |
| 0.5 (2026-11-16) | 5 Transient flow, run control, Python | Transient incompressible flow: PISO with variable-step BDF2 and Crank–Nicolson, outer correctors for large time steps, a time step from the integration error, output at exact times, restart files and control files; the whole flow loop on NVIDIA GPUs; adjoint-ready kernels (dual numbers, transposed matrices, checkpointing); the run contract for scripts and agents (`--json`, `--set`, `result.json`, run record, exit codes and signals, case schema, diagnostic codes, estimate, limits, provenance); the conservation ledger and the divergence report; the bounded central scheme and the boundary conditions of D29 (slip, wedge, tables, backflow, convective outlet); derived fields; expressions, the material database and regions; the Python package `foehncraft` with studies; the grid convergence index of every sweep and a generated verification manual; all scripts in Python. **Taylor–Green time order 2.0–2.4 and space order 2.1; Schäfer & Turek 2D-1 inside all reference ranges, C_D extrapolated to 5.5797 ± 0.0004 from three meshes; GPU flow loop bitwise equal to the CPU on a Tesla T4 (722 ns per cell and iteration on 970 225 polyhedra); a case built in Python is bitwise its TOML file.** Commands `foehn materials`, `foehn post derive`, `foehn schema`; gallery 14–25, recipes 21–35 |
| 0.4 (2026-11-15) | 4 Laminar steady flow | First flow solver: steady incompressible laminar flow with SIMPLEC (and SIMPLE) on polyhedral meshes, relaxation-independent Rhie–Chow interpolation (Majumdar 1988), CG + AMG for the pressure correction, boundary conditions of D29 (walls also moving, velocity and mass-flow inlets, static and total pressure inlets and outlets, symmetry) with a check of their combination and guidance (D30), TOML case file with "did you mean" (D39), log, CSV, wall and line output (D43). **Couette exact to the iteration tolerance on all mesh types; Poiseuille observed order 1.80–2.00; lid-driven cavity within 0.4–1.4 % of Ghia et al. at Re 100 and 1000 in 2D and within 2.5 % on 3D tetrahedra and polyhedra.** New crate `cfdcraft-flow`, command `cfdcraft flow`, thirteen case files, gallery 12–13, recipes 18–20. CPU only; the GPU flow loop follows with stage 5 |
| 0.3 (2026-11-14) | 3 Pressure solver | Pressure linear solver: smoothed-aggregation algebraic multigrid as the preconditioner of CG (Vaněk, Mandel, Brezina 1996), symmetric V-cycle with Chebyshev or Jacobi smoothing, float32 hierarchy with float64 Krylov iteration (mixed precision, D37), one kernel source for the CPU and the NVIDIA GPU, optional hypre 3.2.1 reference (BoomerAMG-PCG, built from source). **10–13 iterations from 512 to 274 625 cells on all mesh types; same solution as hypre to 1.6e-9; on a Tesla T4 6.5 ns per cell and iteration against 10.5 for hypre at 912 673 cells, all GPU results bitwise equal to the CPU.** New command `cfdcraft poisson`, four Poisson case files, gallery 10–11, recipes 15–17 |
| 0.2 (2026-11-14) | 2 Operators on a scalar equation | Finite-volume operators on a scalar transport equation, one source for the CPU and the NVIDIA GPU: least-squares gradients (exact for linear fields, error 6e-15 on all mesh types), linear upwind with Venkatakrishnan and Venkatakrishnan–Wang limiters (no new extrema at K = 0), over-relaxed non-orthogonal diffusion, reproducible gather-based assembly into sliced ELL, BiCGStab. Observed order 1.87–2.08 on hex, tet, poly and mixed meshes in 2D and 3D. **GPU and CPU results bitwise equal on a Tesla T4** (912 673 polyhedra, 311 ns per cell and iteration). New command `cfdcraft scalar`, seven case files, gallery 8–9, recipes 12–14 |
| 0.1 (2026-11-14) | 1 Mesh | Polyhedral face-based mesh with exact geometry (every cell closed to 1e-15, volumes exact), Gmsh 4.1 reader, built-in test mesher (hex, tet, poly, mixed; 2D and 3D; graded, refined, perturbed), quality report, VTKHDF output for ParaView, METIS and RCB partitioning with a solver-ready halo layout and multi-field halo exchange (tested over MPI on 2–4 ranks), particle point location on all mesh types, golden-value test matrix, `cfdcraft mesh` command. G1 rerun on polyhedral and shuffled meshes: bitwise equal; quiet-machine verdict open (O15) |
| 0.0 (2026-11-14) | 0 Foundations and gate G1 | Cargo workspace, CI on Linux, macOS and Windows with MPI runs on 4 ranks, the G1 kernel benchmark. **G1 GPU verdict passed on a Tesla T4: single-source Rust kernels within 1.0 % of hand-written CUDA C, bitwise equal results**; **CPU verdict passed on an i9-9900K: single-source bodies within 1.8 % of hand-written loops, in two runs** |

Release 0.6 closes stage 6: all six turbulence models, the wall
treatments, passive scalars with thermal wall functions and rotating
machinery are verified on the cases of the stage plan (TMR flat plate,
channel against DNS, heated channel and pipe, backward-facing step of
Driver & Seegmiller, impinging jet of the ERCOFTAC database,
Taylor–Couette flow, stirred vessels), with the findings and the open
points listed under "Known limits" of the 0.6 section and in the
[stage 6 verification notes](docs/verification/stage6.md). The generated
[verification manual](docs/verification/manual.md) lists every case with
its pass bands. The Cargo version of the crates is now 0.6.0; the release
is the Git tag `v0.6` (D46).

Release 0.5 closes stage 5: the acceptance checks of the stage plan pass
(Taylor–Green time and space order 2, also at Courant number 9 with outer
correctors; the cylinder benchmark of Schäfer & Turek, at Re 20 inside all
reference ranges, at Re 100 with C_L,max 0.5 % below its range on 24 576
cells; restart bit-identical; expressions from TOML and from Python
give the same run; the energy budget of the cavity closes; the GPU flow
loop gives bitwise the CPU results). The generated
[verification manual](docs/verification/manual.md) lists every case with
its pass bands. The Cargo version of the crates is now 0.5.0; the release
is the Git tag `v0.5` (D46).

Release 0.4 closes stage 4: all acceptance tests of the stage plan pass
(Couette at round-off, that is at the iteration tolerance; Poiseuille order
1.8–2.2 on all mesh types; lid-driven cavity in 2D and 3D at Re 100 and 1000
against Ghia et al.). The Cargo version of the crates is now 0.4.0; the
release is the Git tag `v0.4` (D46).

Release 0.3 closes stage 3: all acceptance tests of the stage plan pass
(iteration count nearly independent of the mesh size, same converged result
as hypre, GPU time per cell compared with hypre and reported). The Cargo
version of the crates is now 0.3.0; the release is the Git tag `v0.3` (D46).

Release 0.2 closes stage 2: all acceptance tests of the stage plan pass
(exact gradients, no new extrema, observed order 1.8–2.2 on all mesh types,
bitwise reproducible results, GPU and CPU bitwise equal on a Tesla T4). The
Cargo version of the crates is now 0.2.0; the release is the Git tag `v0.2`
(D46).

Release 0.1 closes stage 1: all acceptance tests of the stage plan pass
(closed cells, exact volumes, Gmsh files from Gmsh load, VTKHDF files load
with VTK's reader, point location on all mesh types). The Cargo version of
the crates is now 0.1.0; the release is the Git tag `v0.1` (D46).

Release 0.0 closes stage 0: both verdicts of gate G1 are in (GPU: Tesla T4,
CPU: i9-9900K, 2026-11-14, both PASS; O13 closed). The Cargo version of the
crates stays 0.0.0; the release is the Git tag `v0.0` (D46).

---

## 0.6 (2026-11-20) – stage 6: RANS turbulence, heat transfer and rotating machinery

**Highlights.**

- **RANS turbulence, six models.** k-ω SST 2003 (Menter, Kuntz & Langtry)
  with production limiters, sustaining terms and rough walls; standard and
  realizable k-ε; k-ε-ζ-f (Hanjalić, Popovac & Hadžiabdić 2004);
  k-ε-ζ-f₀ (Popovac 2023); the elliptic-blending k-ε (Billard & Laurence
  2012). One kernel source for the CPU and the GPU, the exact or the
  Poisson wall distance, inlet turbulence from intensity and length scale
  or viscosity ratio.
- **Wall treatments.** `low-re`, the SST `automatic` treatment that holds
  C_f within −2.6 to +2.4 % for first-cell y⁺ from 0.24 to 85,
  `scalable`, `standard`, `compound` and `compound-2023` wall functions,
  rough walls (Aupoix for SST, Cebeci & Bradshaw for the k-ε family).
- **SST options.** The curvature and rotation correction (Smirnov & Menter
  2009), a whole domain in a turning frame (`flow.frame_rotation`), and the
  quadratic constitutive relation QCR2000/QCR2013 for the secondary flow of
  non-circular ducts.
- **Heat and mass transfer in turbulent flow.** Passive scalars with
  ν_t/Pr_t diffusion, the Kader and Jayatilleke thermal wall laws, the
  Kays–Crawford turbulent Prandtl number and a high-Prandtl variant with
  variable properties (P-function); periodic boundaries with a heat-flux
  jump for developed heated channels and pipes.
- **Rotating machinery.** The moving reference frame (steady and
  transient), the sliding mesh by conformal stepping and by imprinting the
  interface at every step, several rotors at their own speeds, interfaces
  of any shape of revolution (cylinders, discs, cones), every turbulence
  model on all of them, speeds as `omega` or `rpm`; torque, shaft power,
  the angular-momentum balance and machine performance (head,
  efficiencies, Euler work); relative and blade velocities as derived
  fields.
- **Faster steady runs.** The coupled pressure–velocity algorithm
  (matrix-free, flexible GMRES with a SIMPLE-type preconditioner) with
  pseudo-time steps and a reduction multigrid for momentum and the
  turbulence equations: the backward-facing step settles in about 70
  outer iterations instead of about 600 with SIMPLEC.
- **Meshes from other tools.** SU2, OpenFOAM polyMesh and CGNS (HDF5)
  readers behind one mesh-format interface, besides Gmsh; removed boxes in
  the built-in mesher (a backward-facing step) and annuli for rotors.
- **Results and verification.** Derived turbulence fields (time and length
  scales, Reynolds stresses, wall fields with the viscous weight), time
  means, RMS and covariances of transient runs, the grid convergence study
  of all models on the plate and the channel, the CPU–GPU parity run of the
  turbulence models (14 cases bitwise equal on a Tesla T4), and the
  published sources of all reference data
  (`docs/verification/DATA_SOURCES.md`).
- **Run experience.** The rewritten log with framed blocks and a
  configurable table, case format 0.6 with `foehn case upgrade`, solver
  settings per equation, settings that change during a run
  (`[[schedule]]`), automatic remediation of skewed tetrahedra (D50).
- **Requirements for what comes next.** Stage 7 (energy, buoyancy,
  conjugate heat transfer, thin walls and the heat of the wall film),
  data-driven models (D165–D169), LES and high-Prandtl flows, the
  elliptic-blending Reynolds-stress model (L47) and an applications track.

**Validation in numbers.** TMR flat plate C_f(0.97) +0.37 % from the TMR
value with SST (GCI 1.87 %), −2.6 to −0.3 % from White with scalable k-ε
for first-cell y⁺ 1 to 100, +0.5 % with EB; channel at Re_τ 395 U_b⁺
−0.7 % from the DNS with SST, EB C_f 97.8 / 99.0 / 99.1 % of the DNS at
Re_τ 180 / 395 / 590; backward-facing step x_r/H 5.372 with k-ε (+0.5 %
from the published k-ε, GCI 1.94 %), 6.296 with k-ε-ζ-f₀ (+0.6 % from the
experiment), 6.747 with SST (+7.8 %); impinging jet Nu₀ −0.2 % (SST) and
+7.2 % (k-ε-ζ-f₀) from the experiment; heated channel Nu within −2.6 to
+1.0 % of the DNS at Pr 0.71 for Re_τ 180 to 1020; heated pipe Nu 182.6
against Gnielinski's 180.2; Taylor–Couette flow in a moving frame at
observed order 2.3–2.4.

**Known limits.**

- **k-ε-ζ-f₀ near the wall.** Below first-cell y⁺ 1 the channel is not
  grid-convergent: U_b⁺ drifts up with refinement (17.79, 17.84, 18.11,
  18.39 on levels 0 to 3 against the DNS 17.55), because the k sink of the
  model is linear in k at the wall. A consistent sink is being worked out
  and will come in a later release; meshes with first-cell y⁺ of about 1
  give the values above.
- **k-ε-ζ-f** does not converge on the finest channel mesh of the GCI study
  in 40 000 iterations, and its plate boundary layer stays laminar from the
  TMR free stream.
- **ζ-f and EB on the step and the jet.** The reattachment length of ζ-f
  (7.74) and EB (7.56) is outside its band on the middle mesh (their
  finest mesh is still to run); Nu₀ of the jet is +22.9 % (ζ-f) and +76 %
  (EB).
- **SST plate order.** The observed order of C_f on the TMR plate is 0.43,
  below the band of the specification (the GCI is 1.87 %).
- **Realizable k-ε** gives 3 to 7 % low friction; the compound treatment
  overpredicts the wall shear below y⁺ 30.
- **Open checks that need the maintainer or data:** the GPU run of the two
  SST options (curvature correction, QCR) on the Colab T4; the GCI levels 1
  to 3 of the corrected study on the maintainer's machine; two DNS values
  (the rotating channel of Kristoffersen & Andersson 1993, the square duct
  of Pinelli et al. 2010); the NASA wall-mounted hump data (host not
  reachable from the build machine); manufactured solutions per
  turbulence model.

The parts in detail, in the order of their pull requests (the first
entries, before the fragments of `changelog.d/`, are at the end):

- **Pull requests no longer conflict on bookkeeping.** Changelog entries
  are fragments in `changelog.d/` that the release collects
  (`tools/release/changelog.py`); new gallery and recipe entries take the
  number of their pull request; the number tables of `docs/gallery.md` and
  `docs/recipes.md` are regenerated after each merge by the workflow
  `docs-index.yml` (`tools/docs/indexes.py`); topic pages merge with git's
  union driver. CI checks the lock file of the Python crate first. Rules in
  CONTRIBUTING.md, "Working in parallel".
- **Stage 6 specified (RANS turbulence and wall treatments).**
  `docs/requirements/stage6-turbulence.md` with decisions D129–D139:
  k-ω SST 2003 with Kato–Launder and Durbin production options, standard
  and realizable k-ε, k-ε-ζ-f, all-y⁺ momentum wall treatments, the
  thermal law of the wall (Kader, Jayatilleke) with constant or
  Kays–Crawford turbulent Prandtl number, turbulence numerics, checks,
  case keys, and a verification suite with GCI pass bands (flat plate,
  channel at Re_τ 180/395/590, backward-facing step, impinging jet). The
  elliptic-blending k-ε is open (O17).
- **Stage 7 specified: energy equation, buoyancy, conjugate heat transfer.**
  Decisions D150–D164 with the full text in
  `docs/requirements/stage7-energy-buoyancy-cht.md`: enthalpy conserved and
  the temperature correction solved, low-Mach variable density, SGDH and
  GGDH, Boussinesq and variable-density buoyancy with a balanced-force
  interpolation, fluid–solid interfaces in one temperature matrix, 25
  verification and validation cases with pass bands. The maintainer
  confirmed the new stage order (O18): energy as stage 7, v1.0 after
  stage 8.
- **Data-driven models assessed (step 20, open step O19).**
  `docs/requirements/step20-data-driven-models.md` reviews the literature of
  2025–2026 on learned RANS closure corrections, learned linear-solver
  components, neural surrogates and flow foundation models, judges them for
  industrial accuracy and proposes, not yet decided, D165–D169 (acceptance
  rules for learned components, invariant features in expressions, ensemble
  inversion, data and model licences, learned components that change only the
  path), L37–L38 and P13–P15. Learned models stay outside the solver and come
  in only as fields or compiled expressions.
- **Passive scalars (stage 6 part 6a).** `[scalars.<name>]` adds a scalar
  carried by the flow (D51), solved with the generic scalar transport
  equation of stage 6 (D116) in its bounded form, steady and transient,
  with single-source kernels for the CPU and the GPU. Boundary conditions per scalar in
  `[boundary.<name>.scalar.<scalar>]` (`value` or `flux`, `SCAL001`),
  sources, an `exact` solution with errors and observed orders, `bounds`
  with the bounded BDF2 (D28 (1)). Scalars appear in `residuals.csv`,
  `history.csv`, the console table (`residual.scalars`), the fields, the
  lines and the restart files. The laminar heated channel gives
  Nu = 8.2398 against 140/17 = 8.2353 and second order (gallery 63,
  recipe 63).
- **Exact wall distance.** A bounding volume hierarchy over the wall faces
  finds the nearest wall of every cell (D48, D108 part 1), deterministic
  for any number of threads; the derived field `wall-distance` writes it.
- **Smooth limits.** Hard and C¹ min, max and abs in
  `foehncraft_core::math` for the turbulence limiters (spec 2.6).
- **Data-driven models decided (O19 closed).** Adopted D165 (acceptance
  rules for learned components: do no harm on the stage-6 suite, improvement
  beyond the validation uncertainty, held-out cases, grid convergence with
  the model on, model card), D166 (invariant features in expressions, stage 6
  pull request 6b), D167 (ensemble inversion as a study type), D168 (data and
  model licences), L37 (closure campaign for SST and ζ-f on the diffuser and
  the U-bend, after v1.0) and P13 (learned smoother as a watch item).
- **Stage 6 specification extended.** The curvature and rotation
  correction of SST (Smirnov & Menter 2009) and the quadratic constitutive
  relation QCR2000/QCR2013 as options (D170, pull request 6h), the explicit
  algebraic Reynolds-stress model of Wallin & Johansson with the k-ω of
  Hellsten as a later item (D171), and the k-ω of Wilcox (2006) as an open
  question.
- **Learned components may change only the path (D169).** A learned solver
  component or a surrogate start reaches the same converged answer with or
  without it; CI checks this in reproducible mode.
- **Periodic boundaries (stage 6 part 6p).** Two boundaries of
  `type = "periodic"` that name each other as `partner` form a
  translational pair (D29): the faces are matched by their centroids, each
  side keeps its own faces, and every face loop, gradient and limiter sees
  the cell behind the partner face shifted by the period. The flow through
  the pair is driven by its `mass_flow` (a uniform body force adjusted
  after every pressure correction) or by a fixed `pressure_gradient`; the
  driving force is in the log, the report and the restart file. Checks
  BC030 (pairs that do not match or name each other) and BC031 (a pair that
  nothing drives). Generated boxes keep their perturbed periodic sides
  translated copies. Results of cases without periodic boundaries are
  unchanged bit for bit.
- **Heat-flux jump of the passive scalars (D109 (1)).**
  `periodic = "heat-flux"` in a scalar table of a periodic boundary carries
  the scalar across the pair with the jump Δ = Q/ṁ of the wall fluxes and
  sources, `bulk` fixes its level; refused with fixed wall or inlet values
  (D109 (4)). The periodic heated channel gives Nu = 8.2460 against
  140/17 = 8.2353 and the short periodic heated pipe, one diameter long,
  Nu = 4.3691 against 48/11 = 4.3636 (D111), both at second order; the
  periodic Poiseuille flow equals the open channel (gallery 67, recipe 67).
- **k-ω SST 2003 (stage 6 part 6b, D130).** `[turbulence] model =
  "k-omega-sst"` adds the SST model of Menter, Kuntz & Langtry (2003) on
  the scalar path of 6a: point-implicit (Newton) destruction, the
  production limiter (`menter`, `kato-launder`, `durbin`, `none`), the
  sustaining terms (`ambient`), the ν_t cap, constant overrides by name
  (TURB007), smooth limits for the adjoint. Wall treatments `low-re`
  (k = 0 and ω = 60ν/(β₁ d₁²) on the wall face) and `automatic` (blended
  viscous and log-law u_τ, ω and production in the first cell, fixed with
  a fractional weight in the row assembly); rough walls after Aupoix
  (`roughness`, BC027). Inlet turbulence from an intensity with a length
  scale or a viscosity ratio, or k and ω directly (BC028, BC029, TURB009).
  The wall distance is exact (6a) or the Poisson distance
  (`turbulence.wall_distance = "pde"`). Log and report: the first-cell y⁺
  per wall, the largest ν_t/ν, floor and cap counts (TURB001, TURB002,
  TURB008); `run.jsonl` rows `turbulence` and `turbulence_floors`.
- **Derived field `turbulence-features` (D166).** The invariant features
  of data-driven corrections (Pope's five invariants of the normalised
  strain and rotation, wall Reynolds number, time-scale and viscosity
  ratios, F₁, Q, the streamwise pressure gradient, the curvature of the
  strain) as datasets of the field files.
- **Verification of SST (gallery 68, recipe 68).** TMR flat plate:
  C_f(0.97) = 2.7007e-3 on 256 × 192 cells, +0.37 % from the TMR value,
  C_D +0.41 %; observed order 0.43 (GCI 1.87 %), below the band of the
  specification, recorded as a finding. Channel at Re_τ = 395:
  U_b⁺ = 17.43, −0.7 % from the DNS of Moser, Kim & Mansour. The
  automatic treatment holds C_f within 2.3 % for first-cell y⁺ of 1.7, 27
  and 92 and overpredicts it by up to 16 % in the buffer layer
  (TURB001 warns there). Decay of k and ω at second order in time; restart
  and 1 against 4 threads bit-identical. Curvature correction and QCR
  (D170) follow in 6h.
- **LES and high-Prandtl assessment (step 21).** Review of the
  scale-resolving and high-Pr prompt against the code and the literature
  (`docs/requirements/step21-les-high-prandtl.md`). Adopted D172 (turbulent
  Prandtl number as a closure field), D173 (variable properties, semi-local
  scaling and an integrated thermal wall function, stage 7), D174 (per-face
  convection blending), D175 (equation-free turbulence-model slot), L39–L43
  (LES numerics, SGS models, wall-modelled LES, divergence-free synthetic
  inlet, LES statistics and a campaign test tier, after v1.0), P16–P17 and
  O20. PANS is dropped altogether (X2, was P7).
- **High-Prandtl wall treatment from the TU Graz P-function paper.** D172
  gains the `kays-crawford-hp` turbulent Prandtl number (modified
  Kays–Crawford model of Irrenfried & Steiner with a Pr-dependent wall
  value), and D173's `"integrated"` thermal wall function is specified as
  the revised P-function (full van Driest mixing length with that Pr_t),
  with reference values checked against the paper's figures.
- **Standard and realizable k-ε (stage 6 part 6c, D131).** `[turbulence]
  model = "k-epsilon"` (Launder & Spalding 1974) or
  `"realizable-k-epsilon"` (Shih et al. 1995, variable C_μ and the
  destruction C₂ε²/(k + √(νε))) in one kernel. The standard model is
  bounded by Durbin's realizability limit by default (`production_limiter`
  also `kato-launder`, `menter` as min(P_k, 10ε), `none`); the ε
  destruction is linearised in the Picard form by default. ε at inlets and
  in `[initial]`, or ω converted with ε = C_μ k ω; sustaining terms with
  `ambient = { k, epsilon }`; a restart file of SST continues with k-ε and
  back. Log, residual and output columns name ε.
- **Wall functions (D134).** `wall_treatment = "scalable"` (default of the
  k-ε family, any y⁺, accurate from y⁺ 30), `"standard"` (30 ≤ y⁺ ≤ 300)
  and `"compound"` (Popovac & Hanjalić blend; with the high-Re k-ε claimed
  from y⁺ 30, owner decision: resolved walls go to ζ-f and EB). Rough walls
  shift the log law by ΔB of Cebeci & Bradshaw (`roughness`,
  `roughness_constant`, default 0.5). TURB001 bands per treatment.
- **Verification (gallery 71, recipe 71).** TMR flat plate with the
  scalable treatment: C_f(0.97) within −2.6 to −0.3 % of White for every
  first-cell y⁺ from 1 to 100. Channel at Re_τ = 395: U_b⁺ +4.3 % from the
  DNS at y⁺ 32. Rough pipe and the Moody study
  (`examples/studies/moody_study.py`): smooth pipes −2 to −5 % from
  Colebrook, fully rough +8 to +11 % with C_s = 0.5 (+1 to +3 % with 0.4).
  Findings: the realizable model gives 3 to 7 % low friction; the compound
  treatment overpredicts the wall shear below y⁺ 30.
- **High-Prandtl wall treatment with variable properties.** D172's
  `kays-crawford-hp` and D173's `"integrated"` thermal wall function now
  default to the variable-property P-function model of Steiner & Irrenfried
  (2019). Our 1D pipe check is within 3.1 % of the published DNS Nusselt
  numbers. The CC BY 4.0 thesis of Irrenfried (2020) supplies the reference
  data for the high-Prandtl checks, which answers O20.
- **k-ε-ζ-f model (stage 6 part 6d, D132).** `[turbulence] model =
  "k-epsilon-zeta-f"` (Hanjalić, Popovac & Hadžiabdić 2004): ζ = v²/k and
  the elliptic relaxation f added to k and ε, ν_t = C_μ ζ k T with the
  realizability and Kolmogorov bounds of T and L. Wall treatments
  `compound` (default, any y⁺: k and ζ by zero gradient, ε, P_k and f_w
  blended between sublayer and log law) and `low-re` (k = ζ = 0,
  ε_w = 2νk₁/d₁², f_w = −2νζ₁/d₁²). `zeta` at inlets (BC028 with another
  model) and in `[initial]`; ζ, f, T and v² in the output, the restart
  file and the log columns; a k-ε restart file continues with ζ = 2/3,
  f = 0. A solve lowers k, ε or ζ by at most a factor 10 per outer
  iteration (positivity safeguard).
- **Multigrid of the turbulence equations (D139).** Hierarchy slots on
  both backends (pressure, f, k, ε or ω, ζ). The elliptic f equation runs
  CG with the multigrid by default (`[solver.elliptic]`); the transport
  equations take it with `solver.turbulence.preconditioner = "amg"`,
  built from the symmetric part ½(A + Aᵀ). All hierarchies rebuild on the
  schedule of the pressure (`solver.amg_rebuild`), and their matrices go
  into the restart file (bit-identical continuation). The run-time table
  shows the turbulence model and the passive scalars.
- **Verification (gallery 73, recipe 73).** Channel at Re_τ = 395 resolved
  to y⁺ 0.27: U_b⁺ +1.9 % (`compound`) and +1.4 % (`low-re`) from the DNS.
  TMR flat plate at y⁺ 1: C_f(0.97) +5.1 % from White (finding, outside
  ±3 %). Findings: from the TMR free stream the boundary layer stays
  laminar (the case starts turbulent); near the leading edge the floor of
  k acts (TURB002).
- **k-ε-ζ-f₀ model (stage 6 part 6d2, D132).** `[turbulence] model =
  "k-epsilon-zeta-f0"` after Popovac (2023): f₀ = f + (N_ζ − 1)ζ/T takes
  the place of f and is zero at the wall, ν_t = (C_μ/ζ_eq) ζ k T with the
  constants of the paper. No code of the author is available; the default
  deviates from the printed model where the validation asked for it
  (specification 3.3.1): the sign of the (N_ζ − 1) term of the f₀
  equation, ε/k in the (N_ζ − 1) terms, the k sink k/T_K with the
  Kolmogorov-bounded k/ε, and the realizability coefficient a = 0.6 √2
  of T (with 0.6 the bound halves ν_t in separated shear layers). New
  constants `n_zeta` (ζ-f₀), `c_epsilon1_zeta` and `a` (both ζ-f
  models).
- **Compound wall treatment of 2023.** `wall_treatment = "compound-2023"`
  (default of ζ-f₀, option of ζ-f): the blend
  u⁺ = [y*^{−10/3} + u⁺_log^{−10/3}]^{−3/10} with the non-equilibrium
  factor ψ of the log law from convection and pressure gradient, ε and
  P_k of the first cell from the slope of the blend, f₀,w = 0.
- **Removed boxes in the built-in mesher.** `mesh.generate.remove =
  [{ lower, upper, name }]` cuts boxes out of the hex, tet and mixed
  meshes; the exposed faces form a new boundary (a backward-facing step).
- **Verification (gallery 74, recipe 74).** Channel at Re_τ = 395:
  U_b⁺ +2.7 / +1.4 / −4.2 / −1.5 % at first-cell y⁺ 0.3 / 1 / 10 / 30.
  TMR flat plate: C_f(0.97) −2.5 / +0.9 / −2.0 % from White's correlation
  at y⁺ 1 / 10 / 30, every field converged without a floor, and the
  boundary layer turns turbulent from the TMR free stream (ζ-f stays
  laminar). Backward-facing step of Driver & Seegmiller (new case
  `backstep_zeta_f0`): x_r/H = 6.38, +2 % from the experiment (ζ-f +25 %,
  SST +9 %).
- **Sayings**: two new entries in `docs/sayings.md`, *Deemed prudent* and *Residuals got stuck*.
- **Backward-facing step: three models, two field pictures (gallery 74).**
  `examples/studies/backstep/{keps,zf0,sst}.toml` run the step of Driver &
  Seegmiller on 32 800 cells with k-ε and scalable wall functions,
  k-ε-ζ-f₀ and k-ω SST: x_r/H = 5.36, 6.35 and 6.80 against the
  experiment's 6.26. The gallery shows the velocity with streamlines and
  the turbulent kinetic energy, one panel per model.
- **Figure scripts for 2D fields.** `tools/plot/backstep_paraview.py` makes
  the two pictures in ParaView's Python (`pvpython`; run in ParaView 5.11
  and 6.2), `tools/plot/backstep_chart.py` with matplotlib and h5py
  (recipe 74).
- **Backward-facing step: settled results and half the iterations.** The
  step cases start with the momentum under-relaxation 0.5 and switch to
  0.8 at iteration 100 (`[[schedule]]`), which halves the iterations to a
  settled reattachment length: about 600 on 8 200 cells and 1 750 on
  32 800 cells for k-ε, k-ε-ζ-f₀ and SST. The fine study cases run 3 000
  iterations, `backstep_zeta_f0` 1 500 instead of 4 000. Run to the end,
  SST reattaches at x_r/H = 6.80 (+9 %), not 5.79; k-ε at 5.36 and
  ζ-f₀ at 6.35 on the fine mesh. The figures of gallery 74 are redrawn
  from the settled runs; recipe 74 shows how to check when a result
  has settled.
- **Coupled pressure–velocity algorithm (D176, L11).** `[solver] algorithm
  = "coupled"` solves the momentum equations and the continuity equation
  of a steady outer iteration as one linear system for the increments of
  velocity and pressure: matrix-free with the kernels of the segregated
  solver (no block matrix), flexible GMRES with a SIMPLE-type
  preconditioner, `[solver.coupled]`. The discrete equations and the
  converged result are those of SIMPLEC, and the face fluxes satisfy
  continuity after every outer iteration; the kinetic-energy ledger of
  the cavity at Re = 1000 is the same with both algorithms. CPU and GPU
  from one source (two new kernels).
- **Pseudo-time steps for steady runs (D177, D178).** With the coupled
  algorithm the momentum equations take one pseudo-time step for all
  cells (automatic: the volume of the domain over its boundary area, over
  a velocity scale; `[solver.pseudo_time]`) and the relaxation (default
  0.98) only caps the local Courant number; the turbulence equations take
  the step min(Δτ, 0.5 T) with the time scale T of the model. The
  backward-facing step settles in 65 (k-ε-ζ-f₀), 75 (k-ε) and 80 (SST)
  iterations instead of 525 to 675, on 32 800 cells in 175 instead
  of 1 600; the channel with k-ε-ζ-f₀ converges in 323 instead of about
  2 100, the cavity at Re = 1000 in 452 instead of 1 397 (counts good to
  about 10 %: they react to the last digit of the step). In run time the
  step gains a factor of 1.4 to 1.5 so far, because a coupled iteration
  costs five to seven SIMPLEC iterations; small laminar cases do not gain
  yet (L44, a multigrid for the momentum sweep, is the next step). The
  turbulence step is available with SIMPLEC too
  (`solver.pseudo_time.turbulence = 0.5`, off by default: 230 iterations
  on the step). New cases
  `backstep_zeta_f0_coupled.toml` and `cavity_re1000_coupled.toml`,
  gallery 78, recipe 78, requirements step 22.
- **Log.** The column `K` (`linear.coupled`) shows the Krylov iterations
  of the coupled solve; `residuals.csv` and the run record carry them as
  `krylov` for coupled runs, and the summary of a run shows the
  pseudo-time step of its last iteration (`pseudo_step` in the run
  record). The keys of `[solver.coupled]` and
  `[solver.pseudo_time]` may change during a run.
- **Chart tool.** `tools.plot.settling_chart` draws a monitored value
  against the outer iteration for several runs, with the ±1 % band and
  the iteration from which each run stays inside it.
- **Intermediate results and balance tables no longer change a steady
  run (D31).** A `.result_N` request and the balance table of an
  iteration (`[balance] every`) updated the boundary values and gradients
  for their output and left them in place, so the next iteration started
  from values that had seen one more update: the path of a steady run
  depended on its output schedule (a duct with 24 × 8 cells: 7 % in the
  momentum residual of the iteration after a result; the step case 10⁻⁶
  with SIMPLEC and 10⁻⁴ with the coupled algorithm), though not its
  converged result. Both now put the boundary values, gradients and
  limiters of the iteration back. A run with result files or balance
  tables repeats the plain run bit by bit, with SIMPLEC and with the
  coupled algorithm (test
  `intermediate_results_and_balances_do_not_change_the_run`). Transient
  runs were not affected. The figure of gallery 78 is redrawn from
  histories recorded this way: the SST step settles in 80 coupled
  iterations (85 before). The step on 32 800 cells, run to the end
  with both algorithms (x_r/H = 6.352), settles in 175 coupled iterations
  against 1 600 with SIMPLEC, in 235 s against 327 s.
- **Multigrid for the momentum sweep of the coupled algorithm (D179, the
  momentum part of L44).** The preconditioner of the coupled algorithm
  now applies one V-cycle of a reduction multigrid to each velocity
  component instead of BiCGStab: classical strength of connection, the
  first pass of the Ruge–Stüben splitting into coarse and fine points, a
  restriction that approximates the ideal one by a two-term Neumann
  series (Manteuffel, Ruge & Southworth 2018), direct interpolation, and
  Jacobi relaxation on the fine points only (`amg::reduction`). The
  cycle runs on the level kernels of the pressure multigrid (CPU and GPU
  from one source, float32 hierarchy); the setup is deterministic. On
  the backward-facing step the Krylov iterations per outer iteration
  halve (10 → 5.3 on 8 200 cells, 14 → 6.0 on 32 800), the fine mesh
  settles in 150 instead of 175 iterations, and the run time to the
  settled result on one thread is 7.0 s (SIMPLEC 18.7 s) and 74 s
  (SIMPLEC 243 s, the BiCGStab sweep 131 s); the cavity at Re = 1000 on
  128² converges in 12.1 s (BiCGStab sweep 18.8 s, SIMPLEC 28.3 s).
  Converged results do not change. New keys `[solver.coupled]
  momentum_sweep` (`"multigrid"`, the default, or `"bicgstab"`),
  `momentum_cycles`, `momentum_reuse` and the table
  `[solver.coupled.multigrid]`. Gallery 80, recipe 80, requirements
  step 22 sections 3.1 and 7.7, verification stage 6 part 6f.
- **Kept hierarchies.** A hierarchy of the momentum sweeps is kept while
  the momentum matrix has changed by at most 2 % since its setup
  (`momentum_reuse`) and while its cycle still halves the residual;
  components with the same matrix share one, a component without a
  residual (the third one of a planar case) needs none. The restart file
  holds the matrices of the kept hierarchies (`momentum_amg_use`,
  `momentum_amg_matrix_<n>`), so a continued run stays bit-identical;
  older restart files still load. If a setup fails, BiCGStab does the
  sweeps of that outer iteration and the log carries a warning. The
  report of a run counts the setups of the momentum sweeps and gives
  the time of the turbulence model and of the scalars among its
  seconds per phase.
- **One thread without a pool.** On a pool of one thread (meshes below
  20 000 cells, or `--threads 1`) every kernel call and every level
  operation of a multigrid cycle was handed to the worker thread and
  waited for, about 2 000 times per outer iteration of a cavity run.
  Such loops, and loops of a single task on any number of threads (the
  coarse levels of a multigrid cycle), now run on the calling thread.
  The values are the same bit by bit. Small and medium cases gain with
  both algorithms: the step on 8 200 cells takes 18.7 s instead of 31 s
  with SIMPLEC (525 iterations) and 8.3 s instead of 21 s with the
  coupled algorithm and the BiCGStab sweep (65 iterations); the step on
  32 800 cells with SIMPLEC 243 s instead of 327 s.
- **Chart tool.** `tools.plot.krylov_chart` draws the Krylov iterations
  of the coupled solve against the outer iteration for several runs.
- **Multigrid for the turbulence equations and larger steps with it
  (D180, the turbulence part of L44).** With the coupled algorithm the
  transport equations of the k-ε family (k-ε, realizable k-ε, k-ε-ζ-f,
  k-ε-ζ-f₀) take one cycle of the reduction multigrid of D179 as the
  preconditioner of BiCGStab, with kept hierarchies (`[solver.turbulence]
  preconditioner = "reduction"`, the default there, and `reuse`), and
  the relaxations of the momentum and turbulence equations default to
  0.995 instead of 0.98. The relaxation acts as a local step that
  shrinks with the cells: at 0.98 it made the outer iterations grow with
  the mesh. The backward-facing step settles in 40 iterations on 8 200
  cells and 65 on 32 800 (65 and 150 before), k-ε in 45 (75), the
  channel converges in 231 (323). One thread: the fine step takes 56 s
  instead of 74 s; the coarse steps about the same time as before (7.5 s
  against 7.0 s, k-ε 5.7 s against 4.5 s), since an iteration costs
  more. SST keeps Jacobi and 0.98 (its steady state on the step is not
  held with accurate solves of its equations), laminar cases keep 0.98.
  Converged results do not change. The restart file holds the matrices
  of the kept hierarchies; a failed setup falls back to Jacobi with a
  warning; the report counts the setups. Gallery 81, recipe 81,
  requirements step 22 sections 4.2 and 7.8, verification stage 6 part
  6g.
- **Coupled algorithm: looser pressure solves in the preconditioner, and
  the Schur complement measured.** `[solver.coupled] pressure_tolerance`
  defaults to 0.3 instead of 0.1: the Krylov iteration corrects what the
  pressure solve of the preconditioner leaves, so the Krylov count stays
  and the CG iterations fall by a fifth to more than a quarter. One thread, to the settled
  result: the coarse step takes 7.4 s instead of 8.0 s, k-ε 5.9 s instead
  of 6.4 s; the fine step is unchanged (57.1 s against 56.9 s), since
  its pressure solves stop at their iteration limit either way; the
  results do not change. Two published alternatives to the SIMPLE-type Schur complement
  were measured and do not pay (SIMPLER: a second pressure solve for a
  few Krylov iterations; the least-squares commutator: fine on the
  cavity, at the Krylov limit on the step). Requirements step 22
  section 7.9; the next step is a pressure hierarchy refreshed in every
  outer iteration.
- **Coupled algorithm: the pressure hierarchy of the preconditioner
  follows its matrix (D181, stage 6 part 6h).** The multigrid hierarchy
  of the pressure solve inside the coupled preconditioner is built anew
  when its matrix has changed by more than
  `[solver.coupled] pressure_reuse` (5 %) since its setup, besides every
  `amg_rebuild` iterations as before (`pressure_reuse = 1` restores the
  periodic rebuild alone). A hierarchy up to nine iterations old needed
  9 to 20 CG iterations per solve where a fresh one needs two. One
  thread, to the settled result: the fine backward-facing step (32 800
  cells) takes 39.8 s instead of 57.1 s, 6.1 times faster than SIMPLEC;
  the coarse step 6.2 s instead of 7.4 s, k-ε 4.7 s instead of 5.9 s.
  The settle counts and the converged results do not change; laminar
  cases do not trigger the rule. Requirements step 22 section 7.10,
  gallery and recipe 83.
- **Multigrid setup: flat row builders (L44, no change of results).** The
  reduction multigrid of the coupled algorithm (momentum and turbulence
  equations) builds its strong graph, interpolation, restriction and
  Galerkin products into flat arrays, blocks of rows in parallel joined in
  row order, with one dense accumulator per worker instead of a vector
  per row. The matrices are bit-identical on any number of threads, and
  runs give the same residual history bit for bit. One thread: the
  multigrid setups of the fine backward-facing step take 4.9 s instead of
  6.9 s, the run 33.3 s instead of 39.8 s (7.3 times faster than
  SIMPLEC); the coarse step 5.7 s instead of 6.2 s. Requirements step 22
  section 7.11.
- **Passive scalars of a turbulent flow: turbulent diffusivity and thermal
  wall function (D135, stage 6 part 6f).** A scalar now diffuses with
  ν/Pr + ν_t/Pr_t on every face, and a wall face carries the conductance
  h_w = ρ u*/T⁺ of the wall law `scalars.<name>.wall_function` (`kader`,
  the default, valid for all y⁺ and any Pr, `jayatilleke`, or `none`).
  `turbulent_prandtl` (0.85) and `turbulent_prandtl_model` (`constant`
  or `kays-crawford`) set Pr_t. Fixed-flux walls report the wall value
  from T⁺. The wall file gains `<name>_wall`, `<name>_flux`, `<name>_h`
  and `<name>_t_plus`. The heated pipe at Re = 10⁵, Pr = 0.71 gives Nu
  182.6 (Kays–Crawford 181.3) against Gnielinski's 180.2; the flat
  plate St is +1.7 % to +7.5 % from the correlation for 0.25 ≤ x ≤ 1.
  Cases without a turbulence model do not change. Requirements step 6
  section 5, gallery and recipe 85.
- **Heated channel against DNS (stage 6 part 6f, specification 9.4).**
  `channel_heated_keps.toml` is the Re_τ = 395 channel with a passive
  temperature (Pr 0.71, uniform wall heat flux, periodic heat-flux jump).
  Its Nusselt number on the width 2h is 36.01 with the constant Pr_t and
  35.65 with Kays–Crawford, against 36.39 of the DNS of Abe, Kawamura &
  Matsuo (−1.0 % and −2.0 %). `examples/studies/channel_heat_dns.py`
  fetches the DNS tables of the JAXA database into a git-ignored folder
  and recomputes their Nusselt numbers for Re_τ 180 to 1020 at Pr 0.71
  and 0.025; no DNS data enters the repository (JAXA reserves the
  rights). Gallery and recipe 86.
- **Heated channel at every DNS point (stage 6 part 6f).** The channel of
  `channel_heated_keps.toml` at Re_τ 180, 395, 640 and 1020 and Pr 0.71
  and 0.025 against the DNS of Abe, Kawamura & Matsuo, on uniform meshes
  with first-cell y⁺ ≈ 30: Nu within −2.6 … +1.0 % at Pr 0.71 with both
  Pr_t models; at Pr 0.025 +5.9 … +7.2 % with Kays–Crawford and
  +12 … +22 % with the constant Pr_t. New case `channel_heated_pr0025.toml`;
  `examples/studies/channel_heat_dns.py --run` runs all sixteen points.
  The thermal law now takes T⁺ at the scalable treatment's floor
  ỹ⁺ = max(y⁺, 11.06), as the momentum law does (channel Nu on a y⁺ ≈ 8
  mesh +10 % instead of +28 %). CG stops when the preconditioned residual
  vanishes below the float32 range of the multigrid (an exactly
  one-dimensional flow), as BiCGStab did. The PTX module now has the
  entry points of the two passive-scalar kernels of pull request 85
  (`scalar_turb_faces`, `scalar_wall_faces`): without them a GPU run did
  not load the module; a test checks that the module has every entry the
  CUDA backend loads. Gallery and recipe 87.
- **Elliptic-blending k-ε (stage 6 part 6e, D133).** New model
  `elliptic-blending-k-epsilon`, the BL-v²/k model of Billard & Laurence
  (2012) with the equations and constants of the paper: k, ε, φ = v²/k
  and the blending parameter α (α − L²∇²α = 1, α = 0 at walls), the E
  term in the k equation and the variable C*_ε2, ν/2 in the molecular
  diffusion, ε_w = νk/y²; wall treatments `compound` (default) and
  `low-re`. It runs on the fields, solvers and multigrid of k-ε-ζ-f (φ in
  the place of ζ, α in the place of f); new element-wise kernels
  `squares_add` and `scale_vector` (also in the PTX module). Channel C_f
  at 97.8 / 99.0 / 99.1 % of the DNS at Re_τ 180 / 395 / 590 (the paper:
  100 / 98 / 100 %); TMR plate C_f(0.97) +0.5 % from White; backward-facing
  step x_r/H = 7.59, +21 % (finding: the variable C*_ε2;
  `constants = { c_epsilon4 = 1.83 }` gives 6.24 and a plate 4.9 % high).
  Output fields `phi` and `alpha`, log columns `ResPhi` and `ResAlpha`;
  the key `zeta` sets φ. New cases `channel_eb.toml`, `tmr_plate_eb.toml`,
  `backstep_eb.toml`; specification 3.4 corrected after the paper.
  Gallery and recipe 88.
- **Backward-facing step against the experiment (stage 6 part 6g,
  specification 9.5).** New case `backstep_validation.toml`: the step of
  Driver & Seegmiller (1985) in the layout of NASA's Turbulence Modeling
  Resource (walls from x = −110 H, outlet at 50 H, no-slip top wall at
  9 H). New study `examples/studies/backstep_validation.py`: six models
  on three meshes, the reattachment length with its GCI, C_p and C_f at
  the experiment's stations and the profiles at x/H = 1, 4, 6, 10; the
  experimental data of the TMR (CC0) in `examples/studies/backstep/data`.
  New output `[[output.wall_crossings]]`: the first point along a
  direction where the wall shear changes sign (reattachment or
  separation), in the report and with its observed order and GCI in a
  level sweep; `skip` passes the corner vortices behind a step.
  Results on the finest mesh: SST x_r/H 6.747 (+7.8 %), k-ε 5.372
  (+0.5 % from the published k-ε, GCI 1.94 %), realizable k-ε 6.407,
  k-ε-ζ-f₀ 6.296 (+0.6 %); ζ-f and EB 7.74 and 7.56 on the middle mesh,
  outside the band (their finest mesh follows). The run banner names
  the turbulence model and wall treatment instead of "laminar". Gallery
  and recipe 89.
- **Impinging jet at H/D = 2 (stage 6 part 6g, specification 9.6).** New
  cases `jet_pipe_profile.toml` (the developed pipe flow at Re = 23 750
  on a periodic wedge) and `impinging_jet.toml` (the round jet from that
  pipe on a heated plate, axisymmetric wedge, a passive temperature with
  a uniform wall flux). New study `examples/studies/impinging_jet.py`:
  six models on three meshes, Nu₀ with its GCI, Nu(r/D) with its
  secondary maximum, the stagnation-point anomaly with and without the
  production limiters, and the wall-jet velocity profiles; it compares
  with the experiment once the data are in
  `examples/studies/impinging_jet/data`. On 41 600 cells (with the
  fix of #91) Nu₀ at the axis: ζ-f 173.8 (GCI 0.12 %), SST 141.2
  (2.9 %), EB 249.2 (0.19 %); Nu at r/D = 0.1, the second probe `r01`:
  ζ-f 173.9, EB 205.0, SST 141.1 (GCIs 0.03 to 0.40 %), ζ-f₀ 151.7
  (GCI 5.3 %). The k-ε without a limiter has the anomaly (Nu₀ 211), its
  Durbin default removes 71 % of the excess k; SST has none (its ν_t
  bound holds on the axis). Gallery and recipe 90.
- **Inlet tables with the turbulence of the profile (specification 7.1).**
  The header of a `profile = "table"` file names its columns: `u_n` or
  `u_x`, `u_y`, `u_z`, and `k` with `epsilon` or `omega`, and `zeta`.
  The turbulence columns set the inlet's k, ε or ω and ζ face by face; a
  boundary that gives turbulence values as well is refused (BC029), an
  unknown column too (BC019). `python3 -m foehncraft.profiles LINE.csv
  TABLE.csv` makes the table from a line of a periodic pipe or channel
  run (new module `foehncraft.profiles`).
- **Wall probes report the passive scalars.** Each `[[output.wall_probes]]`
  gives the wall value and the diffusive flux of every scalar there, in
  the log and the report, and as GCI quantities of a level sweep
  (`wall_probes.<name>.<scalar>_wall`, `_flux`).
- **Fix: the thermal wall law of the k-ε family and ζ-f/EB takes
  u* = C_μ^¼ k_P^½ (specification 5.2), not u_τ.** u_τ vanishes with the
  wall shear, so the heat transfer of a stagnation point or a
  reattachment line went to zero (the k-ε Nu₀ of the jet was 37 instead
  of 156). SST keeps its blended u_τ. The heated channel moves by −0.1 to
  −1.2 points (still within 2.9 % of the DNS at Pr 0.71), the heated pipe
  by +0.3; the flat plate is +8 to +9 % for 0.4 ≤ x ≤ 1, and +14 % at
  x = 0.29 just after the k-ε start of the turbulent layer, where k
  overshoots (the test window moves from 0.25 to 0.4).
- **Fix: the full viscous stress on symmetry planes, slip walls and wedge
  planes (stage 6 part 6g).** With a turbulence model the momentum rows
  take the transposed part of the stress, ∇·[μ_eff (∇u)ᵀ], on the
  interior faces; on symmetry faces it was missing. There u_n = 0, so
  that part is μ_eff |S| ∂u_n/∂n along the normal, and the normal stress
  on the face is 2μ_eff ∂u_n/∂n (implicit, as before with the factor 1).
  On the two planes of a wedge it is half of the hoop stress
  −2μ_eff v/r²: the cells next to the axis of the impinging jet had v/r
  too high and Nu too low, more so on finer meshes (the open item of
  #90). Laminar runs do not change. New test
  `tests/symmetry_stress.rs`: a linear flow with constant stress balances
  in every cell, on a wedge and on a symmetry plane (the residual was
  4 % and 65 % of a face flux).
- **Validation: the impinging jet against the ERCOFTAC data (stage 6 part
  6g, part 4).** The six models of the jet study are compared with
  the experiment (ERCOFTAC classic database, case 025, H/D = 2: Nu/Re^0.7
  of Baughn et al. scaled to Re 23 750, and the mean velocity of Cooper
  et al.). Nu₀ of SST (−0.2 %) and ζ-f₀ (+7.2 %) is inside its band, that
  of ζ-f (+22.9 %) and EB (+76 %) outside; the RMS of Nu(r/D) is met only
  by SST (13.7 %); standard k-ε without a limiter is +49 % (the band asks
  for more than +50 %). The data are not in the repository (the wiki of the
  database is licensed CC BY-NC-SA 4.0): `examples/studies/impinging_jet/ercoftac_to_csv.py` makes the
  CSV files from the database files, `.gitignore` keeps them out, and the
  study's new `--no-experiment` makes the gallery charts without the
  points. Verification page (finding 8), specification 9.6 and gallery 90
  updated.
- **Requirements: sliding mesh in both variants (D182).** The sliding mesh
  is planned together with the moving reference frame in stage 6 (6g parts 5
  to 7): conformal stepping (variant A, D56 moves from stage 8) and general
  sliding with the interface imprinted at every step (variant B, was L4).
  Specification, L4, D56, D58 and the stage rows updated; documentation only.
- **Moving reference frame (D55) with laminar Taylor–Couette flow (stage 6
  part 6g, part 5).** `[[mrf]]` tables turn a cell zone steadily about an
  axis: the solved velocity is the absolute one, momentum, turbulence
  variables and scalars are convected with the flux relative to the frame,
  the rotating cells get the Coriolis source, continuity and the pressure
  correction stay on the absolute flux. New: `[mesh.annulus]` (the gap
  between two cylinders with an optional `rotor` zone), `[analytic]
  type = "taylor-couette"`, the diagnostics `MRF001` to `MRF003`, the case
  `taylor_couette_mrf.toml` (L2 error 0.69 % on 1 024 cells, observed order
  2.3 and 2.4, the same in every frame), spec 9.7, verification page, gallery
  and recipe 94. Steady only; the sliding mesh follows (D182).
- **Sliding mesh with conformal stepping (D56, D182, stage 6 part 6g, part
  6).** `[[sliding]]` turns a cell zone of a ring mesh by whole sectors in
  every time step (Δt = shift · 2π/(N |Ω|), refused with the right value
  otherwise, `SLD003`), so interface cells meet new neighbours without
  interpolation and the mesh stays conformal. The solved velocity is the
  absolute one, momentum and scalars are convected with the flux relative to
  the moving mesh (the frame flux of D55, without the Coriolis source), and
  only the time term sees the motion: the old time levels of a rotor cell are
  those of the cell that stood at its place. Laminar flow and passive
  scalars; turbulence models, pressure boundaries on the rotor and steady
  cases are refused (`SLD001` to `SLD004`). New: `[initial] from_analytic`
  for transient runs, the case `taylor_couette_sliding.toml` (L2 error
  1.35 % on 512 cells, observed order 2.2 to 2.3 from 128 to 8 192 cells),
  spec 9.7, verification page, gallery and recipe 95. The variant with a free
  pitch ratio (imprinting) follows.
- **Sliding mesh by imprinting (D182 variant B, stage 6 part 6g, part 7).**
  A rotor and a stator with different numbers of sectors meet on an
  interface that is imprinted at the rotor angle of every time step, so the
  time step is free. `[mesh.annulus] n_theta_stator` makes the two rings
  with the boundaries `rotor_interface` and `stator_interface`;
  `[[sliding]] rotor_interface / stator_interface` names them. Every
  intersection of a rotor face with a stator face is an ordinary internal
  face (`foehncraft_mesh::imprint`: polygons clipped in the unrolled
  surface of revolution), the cells close, mass and momentum flux are
  conserved to round-off, and the solver is built on that mesh in every step
  with the state, the old time levels and the face fluxes carried over. The
  scheme is fully implicit in the step (`extrapolation = false`, at least
  two outer iterations, `SLD003`); `SLD005` and `SLD006` report an interface
  that cannot be imprinted. Taylor–Couette flow from rest with 32 rotor and
  48 stator sectors: L2 error 6.74, 1.29 and 0.25 % on 160, 640 and 2 560
  cells (order 2.4); equal pitch agrees with conformal stepping (1.350 % and
  1.354 %); solid-body rotation is kept; a passive scalar ends at 0.5 %.
  Laminar flow and passive scalars, one rotor, a fixed time step, no
  restart, CPU. Also: the VTKHDF writer writes cells whose faces do not
  close on the nodes of their shape as polyhedra, and probes in a moving
  rotor find their cell again in every step. New case
  `taylor_couette_imprint.toml`, gallery and recipe 96, spec 9.7.
- **Faster rebuild of the sliding mesh by imprinting (D182 variant B).**
  The solver and its pressure multigrid are still built for every step (the
  hierarchy of the step before was measured and does not fit: the rotor
  turns by about one cell per step, and the conjugate gradients needed 5 to
  11 times the iterations), but each build costs less: the dense inverse of
  the coarsest multigrid level takes a third of the work over contiguous
  rows (5 ms to 0.9 ms for 175 unknowns, in every setup of every run),
  multigrid setups below 20 000 rows run on one helper thread instead of
  handing each short parallel loop to the pool, the mesh geometry and the
  discretisation of a step are computed once instead of twice, and the
  state is moved instead of copied. The case file `taylor_couette_imprint.toml`
  runs in 5.6 s instead of 8.4 s, 2 560 cells in 45 s instead of 51 s, with
  the same results (L2 errors unchanged, histories bit-identical or to
  round-off). Verification page stage 6, finding 3 of the sliding mesh by
  imprinting.
- **Walls of a turning zone push the fluid (fix).** A wall in a zone of
  `[[mrf]]` or `[[sliding]]` kept only the tangential part of its `velocity`,
  the rule for walls that stand (D53, `BC010`), so a blade, which moves
  normal to itself, slid along the fluid instead of pushing it. Such a wall
  now moves normal to itself with the zone (Ω × r) and along itself with the
  given velocity: absolute flux ρ **S**·**V**, relative flux zero, as D55
  states. Cylinders that turn about their axis (gallery 94 to 96) are not
  affected. New test: a box turning with its frame carries its fluid (L2
  1.8·10⁻² before, 1.3·10⁻³ after).
- **A stirred vessel with blades and baffles, and field figures of rotating
  machinery (gallery 98, recipe 98).** `examples/meshes/mixer.py` meshes a
  2D baffled vessel with a four-blade rotor in Gmsh (conformal for the
  frame, split 64/96 for the sliding mesh); `examples/studies/mixer/` has the
  steady frame case (Re 100) and the sliding-mesh start-up over three
  revolutions (Re 400); `examples/studies/mixer.py` runs both and draws the
  figures. `tools/plot/rotor_fields.py` draws any 2D rotor run: speed with
  streamlines, the speed relative to the blades in the rotor, vorticity,
  pressure, time series and zooms with the cell edges. Gallery 94, 95 and 96
  have field figures now too.
- **The flow a blade sees: relative and frame velocities of rotating
  machinery (D183, gallery 99, recipe 99).** With an `[[mrf]]` or a
  `[[sliding]]` table, `output.derived` takes five new fields in the
  notation of turbomachinery (absolute **c**, blade speed **u** = Ω × r,
  relative **w** = **c** − **u**): `relative-velocity` (**w** in each turning
  zone, **c** elsewhere, so a pump shows the relative flow in the impeller
  and the absolute flow in the volute), `frame-velocity`, `turbomachinery`
  (c_r, c_u, c_a, c_m, w_u about the axis and the flow angles α and β),
  `relative-total-pressure` and `rothalpy`. `output.frame = "<zone>"` shows
  every cell in the frame of one zone. Sample lines get the columns of **w**
  and of the cylindrical components, point probes of a transient run those
  of **w**; on the sliding mesh at the current position of the cell.
  `tools/plot/rotor_fields.py` draws any dataset of the file, also in the
  rotor only (`NAME@rotor`).
- **Fix:** line samples and the reference comparison of the sliding mesh by
  imprinting located their points on the mesh of the start of the run, whose
  rotor cells had moved on; they use the mesh of the end of the run now.
- **Turbulence on rotating machinery: the moving reference frame and both
  sliding meshes (gallery 100, recipe 100).** Every turbulence model now runs
  with `[[sliding]]` (it was refused, `SLD004`) and is checked with
  `[[mrf]]`. With conformal stepping the old levels of k, ε, ω, ζ and f move
  with the rotor and each step starts from the moved fields; with imprinting
  the turbulence state is carried through the rebuild of every step. The wall
  distance of the turning rotor is kept per rigid part (one bounding volume
  hierarchy each) and only the cells nearer to the edge of their part than to
  its walls search the other part in each step, with a capped query that
  starts from the last nearest face. Several independent `[[mrf]]` zones are
  tested (two zones turning in opposite senses).
- **Fix:** the pressure multigrid built a hierarchy of NaN when a coarse row
  had a non-positive lumped diagonal in the prolongator smoothing (`GEN001` in
  the first iteration of the turbulent stirred vessel with the frame); such a
  row keeps its own diagonal now.
- **Torque, shaft power and machine performance of rotating machinery
  (D184, gallery 101, recipe 101).** For `[[mrf]]` and `[[sliding]]`: the
  torque of the fluid on every wall about the axis of its rotor (the walls of
  a rotor are its `walls`, by default the walls inside its zone), in the
  boundary table, `boundaries.csv`, `forces.csv` and the report; the shaft
  power −Ω·T of each rotor, booked in the step and the final energy balance
  (`balance.csv` column `shaft_power`, the rows of the ledger tables) and in
  the compact balance line; the angular momentum balance about the axis of
  each rotor. `[[output.machines]]` names a rotor, an inlet and an outlet and
  gives the volume flow, the total pressure rise, the head, the hydraulic and
  total-to-static efficiency and the Euler work ṁ Ω Δ(r c_u). New log
  columns `rotor.<zone>.torque|power|angle` (the angle of every sliding rotor,
  in degrees, each with its own speed) and `machine.<name>.flow|pressure_rise|head|efficiency`
  (all appear by default). The convected kinetic energy of the ledger uses the
  flux relative to a turning zone. The torque of a moving wall includes the
  transposed viscous stress: the Laplacian form of the momentum rows misses
  the −u_θ/r of the shear stress (a turning cylinder: 35 % too small before,
  0.5 % against the exact torque now). `tools/plot/rotor_movie.py` writes the
  frames of an animation of a field series in ParaView.
- **Several sliding rotors and the transient moving reference frame (D185,
  gallery 102, recipe 102).** One `[[sliding]]` table per rotor, each with
  its own zone, speed, axis, origin and interface to the stator: all
  interfaces are imprinted together in every step at the angles of their
  rotors (variant B of D182; conformal stepping stays one rotor, an
  interface between two rotors is refused). The exact wall distance of
  turbulence models follows every rotor. The log headings of rotor columns
  name the rotor. `[mesh.annulus] outer_rotor_layers` and `n_theta_outer`
  make a second rotor ring at the outer cylinder; `examples/meshes/twin_mixer.py`
  makes a vessel with two stirrers. `[[mrf]]` now runs in transient cases
  too (the transient frozen rotor), also next to sliding rotors in other
  zones; `MRF001` refuses a zone that a frame and a sliding mesh both turn
  (it refused every transient frame before). New cases
  `examples/cases/flow/couette_two_rotors.toml`,
  `examples/studies/twin_mixer/sliding.toml`,
  `examples/studies/mixer/mrf_transient.toml`.
- **Sliding interfaces of any rotationally symmetric shape (D186, gallery
  103, recipe 103).** One interface per rotor, of any surface of revolution
  about its axis whose meridian is made of straight pieces: a cylinder, a
  disc (an axial gap between blade rows), a cone, or a chain of them such as
  a pot of disc, cylinder and disc around an impeller; several separate
  pieces (discs before and after a rotor row) and discs through the axis
  too. The meridian is found from the mesh; the faces are clipped in the
  unrolled (angle, meridian length) coordinates with a polar chart at the
  axis, and the corners of the pieces are inserted into the neighbouring
  face edges, so every cell of an imprinted mesh closes to round-off on any
  shape. The sliding log line names the shape. New built-in mesh
  `[mesh.revolved]`: a block grid in (r, z) turned about the axis, with
  rotor blocks (own pitch), solid blocks (walls such as an impeller disc),
  cones and wedges at the axis. New cases
  `examples/cases/flow/axial_stage_uniform.toml`,
  `examples/cases/flow/pot_solid_body.toml`,
  `examples/studies/revolved/disc_pot.toml`,
  `examples/studies/revolved/swirl_stage.toml`; new plot tool
  `tools/plot/meridional_fields.py`.
- **README up to date, highlight images, gallery links that land.** The
  README lists what works on main ("What works today": turbulence models
  and wall treatments, heat transfer with passive scalars and thermal wall
  functions, the coupled algorithm, rotating machinery, …) under a banner
  and six highlight tiles linked to their gallery entries
  (`assets/readme/`). Gallery and README images carry their width and
  height (`tools/docs/images.py`, checked in CI), so a link to an entry far
  down a topic page no longer lands thousands of pixels too high while the
  images above it load. Recipe 33's anchor is `python-api` (it clashed with
  the page heading). New rule in CONTRIBUTING.md and the pull request
  template: every pull request keeps the README and the docs current.
- **Sayings**: new entry *Better to know what's not happening* in `docs/sayings.md`.
- **PISO corrector of SIMPLEC type with outer iterations: stable on thin
  wedge cells at an axis at high viscosity (gallery 103).** With
  `outer_iterations` > 1 the pressure corrections of PISO use the SIMPLEC
  coefficient D^C = V/(a_P/α − Σ|a_N|) with the consistent neighbour
  correction Σ|a_N| (u'_N − u'_P)/a_P between two correctors, instead of
  Issa's D = V/a_P with H(u')/a_P. Where a_P is almost all coupling to the
  neighbours (the wedge cells about an axis at a fine pitch and high
  viscosity, a_P about 50 times the time term), Issa's form let the
  pressure grow in every outer iteration: the disc impeller in a pot
  diverged at ν = 1e-4 m²/s. It now runs at that viscosity
  (`examples/studies/revolved/disc_pot.toml` is back at ν = 1e-4). Plain
  PISO keeps Issa's form, which damps more there (bitwise unchanged); new
  key `solver.corrector` (`auto`, `issa`, `simplec`). A converged step is
  unchanged. New case `examples/cases/flow/pot_solid_body_viscous.toml`:
  solid-body rotation at ν = 1e-2 m²/s, L2 error 0.11 % (2.0 before).
- **Rotor images with five blades (gallery 96, 98 to 100, 102, README).** The
  flow fields of the four-blade stirred vessel drew a four-armed pinwheel
  that can be misread as a symbol. Every flow image of a four-blade rotor
  is replaced by the same case on a vessel with five blades and three
  baffles, and the twin vessel (three and four blades) by three and five
  blades. The four-blade meshes, case files, tests and the numbers and
  torque plots from them stay. `examples/meshes/mixer.py --blades
  --baffles`, `examples/meshes/twin_mixer.py --right-blades`,
  `examples/studies/mixer.py --blades --rotors` (all gallery and README
  rotor images in one script). New rule in CONTRIBUTING.md: no images that
  could be read as a symbol; rotors in published images have three or five
  blades.
- **Derived fields of the turbulence models (D72, gallery and recipe 108).**
  `[output] derived` and `foehn post derive` take `turbulence-scales`
  (ε, ω, the time scale of the model, k^{3/2}/ε, ν_t/ν for every model),
  `reynolds-stress` (the Boussinesq estimate ⅔ k δᵢⱼ − 2 ν_t Sᵢⱼ, six
  components), `wall-normal` (∇α/|∇α| of the elliptic-blending k-ε, for the
  GGDH of stage 7) and `turbulence-wall` (u_τ of the wall treatment and y*
  on the walls). New plot `tools/plot/turbulence_fields.py`; tests
  `tests/turbulence_derived.rs`.
- **GPU parity check of the turbulence models (specification 6.6).**
  `tools/bench/stage6_cloud.py` and `tools/notebooks/stage6_colab.ipynb`
  run every turbulence model on the TMR plate and the channel, the heated
  channel and pipe and the coupled step on the CPU and the GPU and compare
  each pair bitwise; `--cpu-only` compares one thread with four.
  `tools/check/compare_runs.py` no longer reports the hash of
  `case.resolved.toml` in `result.json` as a difference when only the
  device keys of `[run]` differ (`--threads`, `--backend`).
- **Statistics of transient runs (D71, gallery and recipe 110).**
  `[output.statistics]` takes time means, RMS and covariances of the fields
  (`statistics.vtkhdf`: velocity, pressure, k, scalars, the resolved
  Reynolds shear stresses and scalar fluxes) and of the probes and force
  coefficients (`statistics.csv`), weighted with the time step, from
  `start` to the end. Each series gets the 95 % confidence interval of its
  mean from batch means and a verdict from the lag-1 correlation of the
  batch means; a window that is too short warns with STAT001.
  `phase_period` writes phase averages (`statistics_phase.csv`). Checked on
  the decaying Taylor–Green vortex against its exact time moments
  (`tests/statistics.rs`).
- **Grid convergence of the turbulence models (spec 9.3, 9.4).**
  `examples/studies/plate_channel_gci.py` runs all six models on three
  meshes of the TMR flat plate and the Re_τ 395 channel and reports the
  observed order, the extrapolated value and GCI_fine (Celik 2008) of
  C_f, C_D and U_b⁺ against White, the TMR and the DNS. Findings on the
  stage-6 verification page: no model is in the asymptotic range on the
  plate meshes (SST within 0.14 % of the TMR value, GCI 0.47 %), ζ-f does
  not converge on the finer plates, the realizable k-ε with scalable wall
  functions moves with every refinement of the channel.
- **SST `automatic` wall treatment holds C_f in the buffer layer (6g).**
  The near-wall ω is blended with the exponent 1.2 and imposed, and the
  first cell takes the equilibrium production of Reichardt's law, both with
  the weight tanh((y⁺/5)²). On the TMR plate C_f now lies within −2.6 to
  +2.4 % of the reference for first-cell y⁺ 0.24 to 85 (before: up to
  +15.6 % at y⁺ 4.6). TURB001 warns for `automatic` only above y⁺ 300.
  The sweep is `examples/studies/sst_wall_sweep.py`; gallery 68 shows it.
- **SST `low-re`: the wall value of ω checked.** A sweep of the factor in
  ω_w = F · 6ν/(β₁ d₁²) from 1 to 100 on the plate and the channel keeps
  F = 10, which has the smallest mesh dependence (verification manual).
- **Curvature and rotation correction of SST (6h, D170, gallery 113,
  recipe 113).** `turbulence.curvature_correction = true` multiplies the
  production of k and ω by the rotation function f_r1 of Spalart & Shur
  (1997) in the form of Smirnov & Menter (2009), from r* = S/Ω and the
  material derivative of the strain rate (constants `c_r1`, `c_r2`, `c_r3`).
  It sees the rotation of `[[mrf]]` zones; not with `[[sliding]]` yet.
  In pure shear the result is unchanged to the last digit.
- **A whole domain in a turning frame.** `flow.frame_rotation = [Ωx, Ωy, Ωz]`
  solves for the velocity relative to a frame that turns at Ω, with the
  Coriolis force (the centrifugal force stays in the pressure), for flows
  that are periodic in the rotating frame. The new case
  `channel_rotating_sst.toml` is the channel at Re_b 5 800 in spanwise
  rotation (Ro 0.15). There the wall shear stress of the unstable side is
  1.63 times that of the stable side with the correction (1.00 without).
- The derived field `curvature-correction` writes f_r1. TURB010 notes the
  correction together with Kato–Launder production.
- **Requirements step 24: elliptic-blending Reynolds-stress model as a later
  item (L47).** Review of a supplied EB-RSM specification (Manceau 2015) on
  the new page `docs/requirements/step24-eb-rsm.md`, with the rotation term
  per frame formulation, tensor boundary conditions, realizability by
  anisotropy scaling and the checks; the algebraic heat-flux companion is
  P18, open points O21–O23. No code.
- **Quadratic constitutive relation of SST (6h, D170, gallery 115,
  recipe 115).** `turbulence.qcr = "qcr2000"` adds the non-linear stress
  of Spalart (2000) to the momentum equations, and `"qcr2013"` adds the
  isotropic term of Mani, Babcock, Winkler & Spalart (2013) (constants
  `c_cr1`, `c_cr2`). The stress goes onto the momentum faces as an
  explicit stress, the path the later algebraic stress model reuses. The
  new case `square_duct_qcr.toml` (Re_b 3 500) shows the eight corner
  vortices with a largest secondary velocity of 1.56 % of U_b. Linear
  SST gives none.
- **Requirements: O22 answered.** The rotating channel of the EB-RSM (L47)
  and the EARSM (D171) runs on the whole-domain relative frame
  `flow.frame_rotation` of PR #113; the step 24 page points to it, and the
  frame term C_ij in the stress equations stays part of L47.
- **Validation reference data: sources and attribution.** New page
  `docs/verification/DATA_SOURCES.md` lists the published sources of the
  reference values of the validation cases (channel, square duct, rotating
  channel, EB-RSM benchmark) with their citations and the quantities used. The
  repository keeps the quoted numbers only, not the data sets; no fetch script
  ships.
- **Validation data sources: citations verified.** The DOIs of Manceau 2015 (IJHFF 51), Pinelli et al. 2010 (JFM 644) and Gavrilakis 1992 (JFM 244) are confirmed in `docs/verification/DATA_SOURCES.md`; no entry is left as "to verify".
- **Validation data sources: channel DNS citation.** `docs/verification/DATA_SOURCES.md` now cites Lee & Moser (2015, JFM 774) for the Re_τ = 180–5200 channel set, with its DOI. The Re_τ 395 and 590 entries keep the Moser, Kim & Mansour (1999) citation, since the 2015 set has no case at those values.
- **Validation data sources: MKM 1999 archive.** `docs/verification/DATA_SOURCES.md` now names the Moser, Kim & Mansour (1999) archive page and its Re_τ = 395 profile listing as the sources of the Re_τ 395 and 590 reference values. The Re_τ = 590 listing is still marked "to verify".
- **ζ-f on the finest meshes.** The multigrid setup of an equation with a
  wide range of coefficients (the f equation of k-ε-ζ-f on a fine mesh)
  refused a healthy coarsest level as singular, because its last pivot was
  compared with the largest diagonal. For the turbulence equations each pivot is
  now compared with its own diagonal (the pressure keeps the strict test). A transport or elliptic equation whose multigrid setup
  fails no longer stops the run: its solve takes the Jacobi preconditioner
  and the log notes it. The GCI study runs the level-3 ζ-f plate and channel with
  the relaxations 0.5, since both diverge with the defaults.
- **GCI levels 1 to 3 on the maintainer's machine.** The plate and channel
  results of the finest three meshes of the turbulence GCI study are in the
  stage 6 verification notes, with findings: SST plate GCI 1.87 %
  (+0.37 % from the TMR value), k-ε plate 2.39 %.
- **Level-3 channels rerun in the cloud.** SST, ζ-f₀ and EB reproduce the
  maintainer's values exactly; ζ-f does not converge on level 3 in 40 000
  iterations. Recorded in the stage 6 verification notes (the cause of the
  drop follows in #128).
- **Industrial applications track proposed (step 25).** A proposal page for an "Applications" gallery with its cases and recipes next to the validation and feature gallery: layout, numbering, an entry template, checks without reference data (balances, convergence, mesh sensitivity, wall resolution, plausibility) with three trust levels, 14 candidate cases that run with the features on main, and the CI plan. Nothing is implemented yet ([step 25 page](docs/requirements/step25-industrial-applications.md)).
- **SU2 meshes and one interface for every mesh format (D187, gallery 127,
  recipe 127).** `mesh.file` reads SU2 native files (`.su2`, ASCII) besides
  Gmsh: tetrahedra, hexahedra, prisms and pyramids in 3D; a 2D file
  (triangles, quadrilaterals) is extruded one cell layer in z with the planes
  `zmin` and `zmax` (`mesh.thickness`, default the median edge length).
  Every `MARKER_TAG` becomes a boundary; `[mesh.boundaries]` maps it to the
  name of a boundary condition. The readers sit behind one format layer
  (`foehncraft_mesh::format`): the format comes from the new key
  `mesh.format` (`gmsh`, `su2`), else from the extension, else from the
  content, and every reader hands cells and named faces to one shared mesh
  builder, so a further format is one parser. `foehn mesh` takes an `.su2`
  file directly. The Gmsh reader moved behind the layer unchanged (a test
  compares fingerprints of all Gmsh fixtures). Gallery 127: a NACA 0012
  C-mesh from PyAero, read from SU2, with a laminar flow at Re 1000.
- **GCI study: the level-3 channel drop resolved.** The drop of U_b⁺ on the finest channel meshes (SST 17.346 → 17.009, EB 17.866 → 16.720) came from runs that stopped at the scaled residual 10⁻⁶ while U_b⁺ was still moving (up to 2.9 % from the converged value). `examples/studies/plate_channel_gci.py` now converges the channel to 10⁻⁹ and refines it in y only (same y nodes, two cells in x); converged, SST (GCI 0.12 %, p 1.0) and EB (below 0.01 %, p 2.1) meet the 0.5 % band on levels 1 to 3. Finding: ζ-f₀ is not grid-convergent below first-cell y⁺ 1 (its wall-cell k does not vanish like y²); the plate results are unchanged.
- **OpenFOAM and CGNS meshes (D188, gallery 129, recipe 129).** `mesh.file`
  reads an OpenFOAM case, `constant` or `polyMesh` directory (ASCII or
  binary, polyhedral cells, patches as boundaries, `cellZones` as cell zones)
  and CGNS HDF5 files (hexahedra, tetrahedra, prisms, pyramids, NGON/NFACE
  polyhedra; BCs and Gmsh physical names as boundaries). `mesh.format` takes
  `cgns` and `openfoam`. Decomposed OpenFOAM cases, `.gz` files and ADF
  files are refused with the fix.
- **The viscous weight of the wall treatment as a wall field (spec 2.3,
  D72).** `turbulence-wall` of `[output] derived` and `foehn post derive`
  gains the boundary dataset `blending`: 1 where the wall is resolved
  (`low-re`, the sublayer of `standard`), 0 where the log law holds (the
  log layer, `scalable`), in between 1 − tanh((y⁺/5)²) of `automatic`,
  e^{−Γ} of `compound` and Γ = dU⁺/dy⁺ of `compound-2023`. It shows which
  walls of a mesh are resolved. The restart file stores it as
  `turb_wall_blending`; from an older restart file `turbulence-wall` is
  refused with CASE005 (continue the run for one iteration).
- **Requirements: thin walls (D189) and the heat exchange of the wall film
  (D190), stage 7.** Thin walls on wall boundaries, baffles and region
  interfaces with layers, heat capacity, sources and in-plane shell
  conduction, solved in the one energy matrix; the wall film exchanges heat
  with plain walls, thin walls and solids (evaporation cools the wall,
  condensation heats it), with ledger pairs and per-face slots built in
  stage 7. Checks V7.26–V7.32, PR 7t, open step O24. No code.
- **Requirements: O24 closed.** The maintainer accepted the thin-wall and
  film points TW1–TW5 as recommended (baffles from face zones, thin walls in
  the energy matrix, shell edges on solids refused first, radiation slot,
  film evaporation needs the vapour species of L2); TW6 (a film without
  sprays) waits until the wall film exists. No code.
- **CI: shorter waits with parallel work.** Pull requests run on Linux
  only (about 5 minutes); macOS and Windows join when a pull request touches
  build files, the toolchain, the workflow, the command line, MPI, VTKHDF,
  restart, run control or the console. Pushes to main run only what the
  gate did not (macOS and Windows tests and quick cases, MPI on macOS, the
  long cases, debug tests), skip docs-only merges, and a newer merge
  cancels the older run. The full tier runs nightly, by hand and with the
  label `full-ci`. HDF5 file locking is off in CI (flaky on macOS).
- **Banner font.** The name of the code on top of the log is set in the
  "Alpha" figlet font by Lennert Stock ("Foehn" over "craft", 133 columns
  with the box), as in the owner's example; the glyphs are in
  `banner_art.rs` with the credit. Only text changes.
- **Outer iterations in the log.** Every logged time step shows the
  residuals and the linear iterations of each outer iteration under its
  box. New `docs/sayings.md` (*Residuen-gaffen*), linked from the README.
- **Log look (D125).** The log opens with "Foehncraft" in block letters
  (see the banner font entry above). Every block is a box of asterisks as wide as its
  content and centred on the page, which has no fixed width any more
  (the banner is 133 columns). Summaries have blank rows, a heading and
  dashed rules; a logged time step starts with a box (step, time, dt,
  Courant number, outer iterations). Messages read
  `%FOEHN-I-CASE, text`: severity and the capitalised name of the source
  file that writes them. Only text changes; results, CSV and JSON are
  bitwise the same. Gallery 27 updated.
- **Cell-quality remediation (D50).** Steady runs on tetrahedra with skewed
  faces (Gmsh meshes without `Mesh.OptimizeNetgen = 1`) diverged with the
  least-squares pressure gradient. The solver now measures the face
  skewness at start-up; when faces exceed `solver.skewness_limit` (0.9),
  `solver.remediation = "auto"` (the default) switches to the Gauss
  pressure gradient (new key `solver.pressure_gradient`) and the
  Venkatakrishnan limiter, and reports it (`MESH001`, `MESH002`,
  `flow.log`, `flow.md`, `flow.json`). Meshes without such faces run
  bitwise as before (all catalogue cases). The raw Gmsh T-junction that
  diverged in iteration 87 converges; gallery 26, recipe 36, case
  `flow/pipe_gmsh_raw.toml` with its mesh in the repository.
- **Divergence report.** A pressure multigrid that breaks down in the
  middle of a run (the pressure-correction matrix no longer positive
  definite) is reported as a divergence (`RUN014`, `divergence.md`)
  instead of "the pressure level is not fixed"; `DIV003` names skewed
  faces anywhere in the mesh when the remediation is off.
- The mesh quality report counts the faces with a skewness above 0.9.
- `examples/meshes/tjunction_netgen.geo`: the T-junction with the Netgen
  optimisation; `tjunction.geo` loses a refinement field that Gmsh 4.15
  rejects (it never applied, the mesh is unchanged).
- `tools/plot/residuals_chart.py`: residual histories of several runs in
  one chart.

**Log output (D125–D127).** The log of `foehn flow` is rewritten so that a log passed on alone shows the run: a banner with version and build, framed blocks (input and output files, run control, mesh information with quality, fluid properties, solution settings as a matrix per equation, boundary conditions), tables with repeating header rows for iterations and time steps, summaries (convergence, verification, results per boundary, ledger, run time, warnings) and a closing frame with the exit code. Frames use asterisk rules; `[log] style = "box"` draws them with box characters. Console and `flow.log` are identical. Gallery 27, recipe 37.

**Configurable log table (D125, part 2).** New `[log] columns` (a list) chooses the columns of the iteration and time-step table: the residuals, the linear iterations, values at probes (`probe.<name>.<quantity>`) or in one cell (`cell.<index>.<quantity>`), the largest and smallest value of a quantity with the cell where it occurs (`max.velocity_magnitude`, `min.pressure`), the mass flow of a boundary and the time quantities; `[log] digits` sets the significant digits. Without `[log] columns` the table is the one of before, and no result changes with the columns (CPU and GPU), so `[log]` is not part of the case hash. `run.jsonl` and the machine formats go to 1.1 (additive): a `columns` row describes the table, a `table_row` row carries the values of every table line, and a `section` row carries every framed block of the log with its label, value and unit, so a monitor draws what the terminal shows. Gallery 28, recipe 38, case `flow/cavity_re100_columns.toml`.

**Case-file format 0.6 (D126).** New key `case_format` (absent = 0.5); old keys are mapped on reading and listed in the log; `foehn case upgrade [--write|--check]` rewrites files and keeps their comments.

**Renamed: "patch" is now "boundary" (D126).** Keys `[boundary.<boundary>]` stay as they were; the names that change are the file `patches.csv` → `boundaries.csv`, the field `patches` → `boundaries` in `result.json`, `flow.json` and the ledger, the quantity names `patches.*`, the VTKHDF array `patch` → `boundary`, and the words in the log and the reports. Scripts that read the old names must follow.

**Solver settings per equation (D127).** `[solver.momentum]` and `[solver.pressure]` replace the flat keys: `relaxation_velocity` → `momentum.relaxation`, `relaxation_pressure` → `pressure.relaxation`, `momentum_tolerance` → `momentum.tolerance`, `momentum_max_iter` → `momentum.max_iter`, `pressure_tolerance` → `pressure.tolerance`, `pressure_max_iter` → `pressure.max_iter`, `pressure_preconditioner` → `pressure.preconditioner`, `amg` → `pressure.amg`, and `scheme`, `limiter`, `limiter_k`, `central_blending`, `central_bounding` → `momentum.*`. All 73 example case files are upgraded. Old files still run (format 0.5).

**Settings that change during a run (D128).** New table `[[schedule]]`: entries with `at_iteration` (steady), `at_step` or `at_time` (transient) and `set = { key = value }` change settings from a prepared point of the run, so a whole run is set up before it starts (under-relaxation, tolerances, iteration limits, time step, output and log cadences). Schedulable are the keys a re-read may change (D40) without `run.limits.*` and `run.reread_every`; everything else (mesh, boundaries, schemes, fluid) is refused when the case is read, with the reason (`SCHED001`). Every change is a line in `flow.log`, a row in the "Run control" frame and an `event` row `RUN015` in `run.jsonl` (`key`, `old`, `new`); entries that have fired stay in force over a re-read, a restart applies the entries before its start. Case format stays 0.6 (a new optional table). Without `[[schedule]]` results, files and logs are bitwise unchanged. Gallery 29, recipe 39, case `flow/cavity_re100_schedule.toml`.

CPU results are bitwise unchanged against 0.5 (checked with `tools/check/compare_runs.py --across-versions` on the steady and the transient cavity; only the renamed names differ).

## 0.5 (2026-11-16) – stage 5: transient laminar flow, run control and Python

**Highlights.**

- **Transient flow.** PISO with variable-step BDF2 and Crank–Nicolson,
  outer correctors (Courant number 9 at time order 2), an adaptive time
  step from the integration error with exact step rejection, output at
  exact simulated times, potential-flow start, restart files that
  continue a run bit for bit, control files and re-reading the case file
  during a run, probes and force histories. Verified on the Taylor–Green
  vortex (time order 2.0–2.4, space order 2.1) and the cylinder benchmark
  of Schäfer & Turek (Re 20 and Re 100).
- **GPU.** The whole flow loop runs on NVIDIA GPUs from the same kernel
  source as the CPU, bitwise equal to the CPU on a Tesla T4.
- **Adjoint-ready kernels.** Dual-number derivatives of every kernel,
  transposed matrices as gathers, binomial checkpointing (D59).
- **Run contract and health.** `--json`, `--set`, `result.json`, the run
  record, exit codes and signals, the case schema, diagnostic codes, the
  estimate before a run, limits, provenance; the mass and kinetic-energy
  ledger and a divergence report with ranked causes and fixes.
- **Numerics and boundary conditions.** The bounded central scheme, slip
  walls, wedges, inlet profiles from tables, backflow treatment and
  convective outlets; derived fields (vorticity, Q, λ₂, wall shear stress,
  y⁺) with the output and from restart files.
- **Expressions and materials.** Boundary values, body forces and the
  viscosity as expressions of x, y, z and t; a material database with
  sources; regions.
- **Python.** The package `foehncraft`: cases, expression objects, runs in
  process with callbacks, results, materials and parametric studies; all
  scripts of the repository are Python now.
- **Verification.** The grid convergence index of every sweep of three or
  more meshes or time steps, and the generated verification manual.

The parts in detail, newest first:

**Scripts in Python (stage 5 part 6d, D61 (2)).** The shell scripts are
ported to Python and `tools/` is split into groups that run as modules
from the repository root: `tools/check/` (`ci_cases`, `compare_runs`,
`backend`, `c_free`, `kernels`, `vtkhdf`), `tools/bench/` (`build_hypre`
and the cloud GPU drivers with their shared steps in `cloud`),
`tools/mesh/`, `tools/plot/` (the charts), `tools/release/`
(`verification_manual`) and `tools/notebooks/` (the Colab notebooks).
`tools/ci_cases.sh SET` is now `python3 -m tools.check.ci_cases SET`,
`gpu/ptx-kernels/build.sh` is `python3 gpu/ptx-kernels/build.py`; CI,
docs, examples and error messages are updated, and CI lints the scripts
(ruff, syntax and names). See `tools/README.md`.

**Discretisation uncertainty and the verification manual (stage 5 part
6c, D67, D68).** A sweep of three or more meshes or time steps reports
for every monitored quantity (patch mass flows, mean pressures and forces,
force coefficients, benchmark values, viscous and numerical dissipation)
the observed order by the fixed-point iteration of Celik et al. (2008),
the extrapolated value, the band φ1 ± U, the grid convergence index and a
verdict (monotonic, oscillatory, divergent, unchanged): a table in
`flow.md` and `uncertainty` in `flow.json` and `result.json` (an
addition; `Result.uncertainty` in Python). Level sweeps no longer need an
exact solution or a benchmark; `expected_order` without `[analytic]` is
now refused (CASE005) instead of being ignored. New case
`cylinder_re20_gci.toml` (the benchmark 2D-1 on 6 144 to 98 304 cells,
`long` set): C_D extrapolates to 5.5797 ± 0.0004 inside the reference
range. New script `tools/release/verification_manual.py` writes
`docs/verification/manual.md`, `results.json` and a chart of the observed
orders from the case files and a catalogue run. Gallery 25, recipes 34
and 35 on new Verification pages; `tools/gci_chart.py`. Runs are bitwise
unchanged.

**The Python package (stage 5 part 6b, D61, D62).** `python/` holds the
package `foehncraft` with a native module (`crates/foehncraft-py`, PyO3,
built with maturin, abi3 for Python 3.10+): `Case` (read, change with
dotted keys, write; the TOML round trip holds for every example case;
relative paths kept when saved elsewhere), expression objects with Python
operators that render the strings of the case file (`foehncraft.expr`,
evaluated with the solver's own compiler), `run` in process or through
`foehn` with callbacks on every iteration or time step that can stop the
run, `Result` (summary, CSV tables, VTKHDF fields), the material database
and material builders, and `foehncraft.study` (grid, list, Latin
hypercube, Sobol; one case file per variant, `study.csv`, resume). The
engine gains `engine::Observer` (sees every row of the run record) and
`RunOptions::observer`; runs are bitwise unchanged. New CI job "Python
package" (ruff, pytest). New scripts `examples/python/poiseuille_inflow.py`
and `examples/studies/couette_viscosity_study.py`; gallery 24, recipe 33
on the new Python pages, case catalogue page "Python scripts".

**Expressions, materials and regions (stage 5 part 6a, D53, D102, D112,
D121, D124 W3/W4).** Boundary values (velocity, normal velocity, pressure,
total pressure, mass flow), the body force and the viscosity may be
strings of a small expression language in x, y, z and t, with
`[constants]` and piecewise linear `[tables]` (inline or from a file); the
strings are checked when the case is read (EXPR001–EXPR005, with the
place and a "did you mean") and compiled to a stack-machine program that
the solver evaluates with new in-house elementary functions
(`foehncraft_core::math`, D93 (6)), the same on every platform. Fluid
properties from a built-in material database (air, water, engine oil,
ethylene glycol, glycerol, each with its source) with `[fluid] material`,
`temperature` and `operating_pressure`; `[materials.<name>]` with `base`
and the forms constant, polynomial, Sutherland, power law, table (linear
or monotone cubic) and expression; the new command `foehn materials`
lists and evaluates them (also `--json`). `[regions.<name>]` assigns cell
zones to materials; non-fluid regions and unknown zones are refused
(REG001, REG002); different viscosities per region give a viscosity per
face. Walls take a `mass_flux` (blowing and suction), checked for closed
volumes (BC026), and expose per-face exchange coefficients
(`derived::wall_exchange`). The resolved property program is written into
`case.resolved.toml`; `properties_hash` and `materials_database` join the
provenance. The random-number key gains a stream word (`rng::keyed`;
meshes unchanged). `foehn flow --check` also checks the regions against
the mesh. New exact solutions `couette-viscosity` and
`couette-transpiration`; new cases `flow/couette_viscosity.toml`,
`couette_transpiration.toml`, `poiseuille_water.toml`,
`poiseuille_expression.toml`, `taylor_green_expression.toml` and
`channel_ramp.toml`; gallery 23, recipes 31 and 32 on the new
"Expressions and materials" pages; `tools/property_chart.py`. Cases
without the new keys run bitwise as before.

**Derived fields of stage 5 (part 5c, D72).** `[output] derived` adds
fields computed from the solution to every `fields*.vtkhdf`: `vorticity`
(components and magnitude), `q-criterion`, `lambda2`, `enstrophy`,
`strain-rate`, `pressure-coefficient` (reference state in
`[output.reference]`: velocity, pressure, density), `total-pressure`, and
on the patch blocks of the file `wall-shear-stress` and `y-plus`. The new
command `foehn post derive CASE` writes the same file afterwards from a
restart file (`--restart`, default `latest`; `--fields`, `--file`), equal
byte for byte to the file of the run. The VTKHDF writer gains boundary
face fields on the patch blocks; the `units` attribute lists the derived
datasets; the estimate counts them; `output.derived` and
`output.reference` may change during a run. `tools/render_mesh.py` gains
`--diverging`. New case `flow/couette_derived.toml` (exact on any mesh);
gallery 22, recipe 30 on the new Post-processing pages. Without `derived`
every output is as before.

**Boundary conditions of stage 5 (part 5b, D29).** New condition types
`slip-wall`, `wedge` and `convective-outlet` (∂u/∂t + U_c ∂u/∂n = 0 in
transient runs, the outlet flow balanced with the inflow, U_c by default
the mean outflow speed); inlet profiles from a table (`profile = "table"`,
`table`, `table_coordinate` x, y, z or the radius about an axis, the
table's sha256 in the case hash); backflow treatment at static and total
pressure outlets (`backflow = "total-pressure"`, the new default,
`"zero-gradient"`, the treatment before, or `"suppress"`, which closes the
faces where the flow would enter). Axisymmetric wedge meshes from any 2D
mesh of the box generator (`mesh.generate.wedge`, the axis without a
patch), and the exact solutions `uniform` and `pipe-poiseuille`. Seven new
diagnostic codes BC019–BC025. The restart file of a run with a
convective outlet gains the datasets `ub_old_x|y|z`. New cases
`flow/channel_slip.toml`, `poiseuille_table.toml` (with
`poiseuille_profile.csv`), `pipe_wedge.toml`, `cylinder_re100_short.toml`
and `cylinder_re100_convective.toml`; gallery 21, recipe 29. Runs whose
outlets see no backflow are bitwise as before; the Poiseuille cases on
perturbed prisms, polygons and mixed cells, where a few outlet faces see
inflow in the first iterations, change by 2e-6 relative
(`backflow = "zero-gradient"` gives the earlier results bitwise).

**Wall film requirements (D123, D124).** Requirements for the wall film of the Lagrangian sprays (L12): an Eulerian thin film on the wall faces with parcels for splash, rebound, separation and stripping; impingement regimes (Bai–Gosman, Kuhnke), capillarity and dry patches in the first version, CHT coupling of the heated film, ledger and reproducibility rules, staging F1 to F5 and validation cases; Milestone 1 keeps a surface mesh and surface halo of the wall faces, per-face wall boundary values including a wall mass flux, and gravity as a case vector. Documentation only.

**Numerics of stage 5 (part 5a).** A bounded central convection scheme
(`solver.scheme = "central"`, D25) with the bound β and a blend with linear
upwind; the time step from the error estimate of BDF2 and Crank–Nicolson
(`time.error_tolerance`, D65), with steps above the tolerance taken back
exactly and repeated shorter (`solves[].transient.rejected` in
`result.json`, an `error` value in the step rows of `run.jsonl`); output in
simulated time that the steps land on (`output.fields_interval`,
`fields_times`, `monitor_interval`, `run.restart_interval`,
`restart_times`, D66); and a potential-flow start
(`initial.potential_flow`, D52). The restart file gains the attribute
`dt_prev2` and, with the error estimate, the datasets `u_old2_x|y|z`;
older files still read. New cases `flow/taylor_green_error.toml`,
`flow/cylinder_potential.toml` and the variant `bounded-central` of
`scalar/oblique_step_2d.toml`; gallery 20, recipe 28,
`tools/time_step_chart.py`. Runs without the new keys are bitwise as
before.

**Material property requirements (D94–D103, L23–L25).** Requirements for a
material database and the handling of properties: database, case materials
and a resolved, hashed property program; property forms (constant,
polynomial, NASA-7 and NASA-9, Sutherland, tables with linear or monotone
cubic interpolation, polynomial c_p expressions, mixing rules, kinetic
theory, equations of state); checks; GPU evaluation; SI only; no property
changes at run time; a built-in starter set with IAPWS water and steam;
importers, table generators and a generic interface for external property
sources. Implementation starts in stage 5 PR 6 (μ and ρ) and continues with L1.

**Renamed from CFDCraft (and the short-lived Ennecraft) to Foehncraft**, after
the foehn wind of Innsbruck (see the README), so that the name cannot be
confused with other projects. Crates are now `foehncraft-*`, the command is `foehn`, the
environment variables are `FOEHNCRAFT_*`, and the repository is
`chiefenne/Foehncraft`. Entries below keep the names of their time.

**Faster CI.** Pull requests get a fast gate (about 5 minutes instead of
about 25): a `ci` build profile (release optimisation without LTO, results
bitwise equal), the case catalogue split over parallel jobs
(`tools/ci_cases.sh`), and a VTKHDF check that no longer computes unused
volumes of field files (11.5 minutes to seconds). Push to main, a nightly run
at 02:07 Vienna time and pull requests with the label `full-ci` run the full
tier: `release` profile, debug tests, every case including the cylinder at
Re 100, and the slow tests. One summary check, "CI OK". Main warms the `ci`-profile cache that pull
requests read (without it every PR compiled its dependencies from scratch).

**Run contract for scripts and agents (stage 5, part 4a).** `foehn flow`
can be driven by a program: `--set key=value` overrides (TOML values),
global `--json` (one JSON object on standard output, human text on
standard error) and `--quiet`, exit codes 0 finished, 1 error, 2 failed
checks, 3 stopped; SIGUSR1 writes a restart file, SIGTERM, SIGINT and
SIGUSR2 stop cleanly, `--restart latest` continues bitwise; budgets in
`[run.limits]` (wall time, steps or iterations, output volume). Every run
writes `result.json` (status, overrides, provenance with commit and case
hash, diagnostics, results, SHA-256 of the outputs; byte-identical between
identical runs) and the run record `run.jsonl`. 41 stable diagnostic
codes with key paths, locations and `--set` fixes; `foehn schema` prints
the JSON Schema of the case file, `docs/keys.md` is generated from it;
`--check` estimates memory, output and wall time; `--read-root` and
`--write-root` confine the files; CSV, VTKHDF and restart files describe
themselves; HDF5 files have no timestamps. The CLI is a thin layer on
`foehncraft_flow::engine`; `foehncraft-flow` without default features is
checked C-free in CI. Formats and codes: `docs/interfaces.md`; recipe 26,
gallery 18, case `cavity_re100_budget.toml`. No result changes: the
catalogue gives bitwise the fields of main.

**Health checks: conservation ledger and divergence report (stage 5,
part 4b).** Every flow run books mass and kinetic energy (D88): storage,
the terms of each patch (viscous work, convected energy, pressure work,
mass flow), the body force, the viscous dissipation of the discrete
operator, and the remainder, the numerical dissipation, which the second
law wants positive (D90 part 4), split into its pressure part (the
momentum interpolation) and the part of the convection and time
schemes. Steady runs give the
table every `[balance] every` iterations and at the end; transient runs a
compact line with every progress line and a time-integrated ledger of the
run (trapezoidal rule). `balance.csv`, `balance.svg`, `balance` rows in
`run.jsonl`, `solves[].balance` in `result.json`, tolerances with the
codes `BAL001`–`BAL003` (the new key `[balance] tolerance_energy` for
the negative numerical dissipation) and an optional stop. Measured: on the lid-driven
cavity the numerical dissipation is positive, about 1.5 % of the lid
power and falling with the mesh; on the Taylor–Green vortex it falls with
order 2.5, its pressure part with order 3. A run that diverges (D63) now writes its last valid state and
restart file, a VTKHDF snapshot of the failing cells and their
neighbours, and a report (`divergence.md`, `divergence.json`) with the
likely causes ranked by their evidence and a `--set` fix for each, and
ends with `RUN014` (exit code 1). A linear-solver breakdown counts as a
divergence, and a NaN in the continuity residual no longer slips past the
check. The restart file can hold point coordinates and cell volumes
(D91). Recipe 27, gallery 19, case `cavity_re100_energy.toml`,
`tools/balance_chart.py`, `tools/energy_budget_chart.py`. No result
changes: the catalogue gives bitwise the fields of part 4a.

### Adjoint-ready kernels (part 3)

Preparation for a discrete adjoint solver (D59; the solver itself is L17).
No result changes: all flow cases give bitwise the files of part 2.

- **Dual numbers:** `foehncraft_core::Dual<N>` implements `Real`, so every
  kernel body also runs on dual numbers and returns the exact derivatives
  of its results, from the same source. Switches (upwind direction,
  limiter branch, clipping) follow the values, i.e. they are frozen at the
  linearisation state; the limiter array ψ can be held constant.
  **Jacobians of 11 face, cell and boundary kernels agree with central
  differences in every input array (320 580 derivatives, largest relative
  difference 5e-7).**
- **Shape derivatives:** the geometry of `foehncraft-mesh` (face areas and
  centroids, cell volumes and centroids) is generic over the coordinate
  type, so derivatives with respect to the node coordinates come from the
  same code (165 024 derivatives against central differences, 2.8e-7).
- **Transposed matrices:** a transpose map of the matrix storage gives Aᵀ
  in the storage of A and the product Aᵀx without storing it, both gathers
  without atomics (`foehncraft_kernels::adjoint`; flow-device calls
  `transpose`, `matrix_product`, `transposed_product`; PTX entries
  `flow_transpose`, `flow_matrix_product`, `flow_transposed_product`).
  ⟨Av, w⟩ = ⟨v, Aᵀw⟩ to 1e-14 on the matrices of a SIMPLEC solve; CG with
  the multigrid and BiCGStab solve the transposed systems.
- **Binomial checkpointing** (`foehncraft_flow::checkpoint`, Griewank &
  Walther 2000): an optimal schedule (checked against a dynamic program),
  a driver for any primal run and the in-memory snapshot of the flow
  solver; the reverse sweep of a BDF2 run meets every level bitwise.
- Gallery 17, recipe 25, `tools/checkpoint_chart.py`, report
  [stage 5, part 3](docs/verification/stage5.md#adjoint-ready-kernels-part-3).

### GPU flow loop (part 2)

- **The flow loop on the GPU:** `run.backend = "cuda"` or
  `foehn flow --backend cuda` runs every flow case, steady and transient,
  on an NVIDIA GPU (compute capability 7.0 or newer), with all fields,
  matrices and old levels in GPU memory; results bitwise equal to the CPU
  (D31). **Tesla T4 (Colab): all 8 cases bitwise equal to the CPU; 722 ns
  per cell and SIMPLEC iteration on 970 225 polyhedra (7.1 times the two
  CPU threads of the VM), 1 049 ns per cell and PISO step on 98 304 cells
  (9.2 times).** Gallery 16.
- **One source down to the argument lists:** the 27 operations of the flow
  loop are entries of `foehncraft_kernels::entry` (a `#[repr(C)]` argument
  block and a body per index). The CPU loops and the PTX entries `flow_*`
  run the same body with the same block; the flow device is written once
  over an executor (`device::kernel`).
- **CUDA in one crate (D60):** new crate `foehncraft-cuda` holds all CUDA
  code (flow executor, the scalar and multigrid backend, the G1 backends
  and the CUDA C reference, moved from `foehncraft-fv` and
  `foehncraft-bench`); the solver crates find it through registered
  factories. CI lint `tools/check_backend.py`. The feature is now
  `foehncraft-cuda/cuda` (the CLI's `cuda` feature is unchanged).
- New tools `tools/compare_runs.py` (bitwise comparison of two runs),
  `tools/stage5_cloud.py` and `tools/stage5_colab.ipynb` (CPU against GPU
  on eight cases up to one million polyhedra); new cases
  `cavity_poly3d_gpu.toml`, `cylinder_re100_gpu.toml`; recipe 24.
- The CPU results of all flow cases are bitwise those of part 1. A missing
  NVIDIA driver now gives a message instead of a panic.

### Transient laminar flow (part 1)

Stage 5 comes in several pull requests and is released as 0.5 when it is
complete. The first part, the transient core and its benchmarks
([stage 5](docs/verification/stage5.md)):

- **PISO (Issa 1986)** with variable-step BDF2, Crank–Nicolson (first step
  BDF1, the start of Rannacher 1984) and BDF1; two correctors with the
  neighbour correction between them; optional outer iterations per step
  (PIMPLE-type) with a residual tolerance; the convecting flux and the
  deferred corrections extrapolated to the new time, which keeps plain
  PISO second order. A case is transient when it has a `[time]` table.
- **Adaptive time step** from the largest Courant number, with a growth
  limit inside the zero-stability range of BDF2.
- **Rhie–Chow time-step term (D27)** whose steady state is exactly the
  steady solution, for any Δt (cavity: equal to SIMPLEC to 2e-12).
- **Restart files (D42, HDF5):** a continued run is bit-identical to the
  run without the interruption. **Control files and re-reading (D40):**
  `.stop`, `.result_N`, `.read_toml_N`, `run.reread_every`. **Monitors
  (D43):** `history.csv`, `forces.csv`, `probes.csv`, field series,
  force coefficients.
- **Built-in cylinder mesh** `[mesh.cylinder]` (O-grid and wake block) and
  the parabolic inlet profile; the Schäfer & Turek benchmarks 2D-1 and
  2D-2 with their evaluation (C_D, C_L, St, Δp).
- **Results:** Taylor–Green time order 2.21 (BDF2) and 2.00
  (Crank–Nicolson) on the finest pair, space order 2.09, order 2 at
  Courant numbers up to 9 with outer iterations; cylinder Re 20 inside all
  reference ranges on 98 304 cells; Re 100 on 6 144 cells with C_D,max,
  St and Δp inside the ranges, on 24 576 cells C_L,max 0.5 % below its
  range and the others inside.
- `cfdcraft flow --restart FILE --end-time T`; seven new case files,
  gallery 14–15, recipes 21–23. The Release workflow now refreshes the
  cached README badges (recipe 6).

## 0.4 (2026-11-15) – stage 4: laminar steady flow

New crate `cfdcraft-flow`, kernel bodies `cfdcraft-kernels::flow` (one
source for the CPU and the GPU) and the command `cfdcraft flow`:

- **SIMPLEC and SIMPLE (D19).** Segregated outer iterations on the
  collocated, cell-centred polyhedral grid: momentum with implicit upwind
  and deferred linear-upwind convection, orthogonal diffusion implicit and
  the non-orthogonal and skewness corrections deferred, implicit
  under-relaxation; pressure correction with the SIMPLEC factor
  D^C = V/(a_P/α − Σ|a_N|), solved with CG and the stage-3 multigrid
  (hierarchy rebuilt every 10 outer iterations); reference cell when no
  boundary fixes the pressure.
- **Rhie–Chow interpolation (D27)** with the relaxation term of Majumdar
  (1988): the converged result does not depend on the under-relaxation
  (lid force equal to 3e-11 for α_u = 0.7 and 0.9).
- **Boundary conditions (D29):** wall (also moving), velocity inlet (vector,
  normal speed or the exact profile), mass-flow inlet, static- and
  total-pressure inlet and outlet (freely combined), symmetry (implicit per
  component). Periodic boundaries follow with stage 6.
- **Boundary-condition check (D30)** before every run and with
  `cfdcraft flow --check`: missing and unknown patches with "did you mean",
  net mass flow without a pressure boundary, over-specified flow paths,
  wall normal velocity, unused turbulence values, non-planar symmetry,
  reference-cell note; the guidance table of sensible combinations in the
  case catalogue.
- **Case file (D39):** `[fluid]`, `[flow]`, `[initial]`, `[pressure]`,
  `[boundary.<patch>]`, `[solver]`, `[analytic]`, `[verification]`, `[run]`,
  `[output]`; unknown keys and values with "did you mean"; the resolved case
  with all defaults in `case.resolved.toml`.
- **Output (D43):** `flow.log` with a header of what is written and how
  large, `residuals.csv`, `patches.csv` (mass flow, mean pressure, force),
  `wall_<patch>.csv` (pressure and wall shear stress per face, in global
  face order), `line_<name>.csv`, `fields.vtkhdf`, `flow.md`, `flow.json`
  and SVG charts (residuals, convergence, centre lines against the
  reference).
- **Verification:** exact solutions (Couette, Poiseuille) and the
  centre-line data of Ghia et al. (1982) for Re 100 and 1000; thirteen case
  files (CI runs all of them on Linux and the quick ones on macOS and
  Windows); tests for Couette, Poiseuille order, the cavity, relaxation
  independence, quick convergence with symmetry planes, the check, "did you
  mean" and bitwise equality on 1 and 4 threads.
- `Backend::set_matrix` (new matrix values on the same pattern) on the CPU
  and CUDA backends; `tools/render_mesh.py --no-edges`.
- Documentation: [stage 4 verification](docs/verification/stage4.md),
  gallery 12–13, recipes 18–20, flow cases, keys, boundary conditions and
  guidance in the case catalogue.

---

## 0.3 (2026-11-14) – stage 3: pressure solver

New module `cfdcraft-fv::amg`, kernel bodies `cfdcraft-kernels::amg` (one
source for the CPU and the GPU, float32 and float64) and the command
`cfdcraft poisson`:

- **Smoothed-aggregation multigrid (D35).** Strength of connection,
  three-phase greedy aggregation in the cell order, filtered and
  Jacobi-smoothed prolongator, Galerkin coarse operators, explicit
  restriction (gather, no atomics), dense coarsest level; deterministic host
  setup in float64. Plain aggregation as an option.
- **Symmetric V-cycle** with a Chebyshev (default) or damped Jacobi smoother
  in D⁻¹A, λ_max from a fixed power iteration.
- **Mixed precision (D37).** The hierarchy (including a float32 copy of the
  fine matrix) in float32, CG in float64; `precision = "double"` for
  comparison. R10 allows f32 only in the multigrid kernels and rejects
  contracted, approximate and flush-to-zero f32 arithmetic.
- **Flexible CG** (Notay 2000). `[solver] linear_solver = "cg"` and
  `preconditioner = "amg"` in `cfdcraft scalar` (also AMG with BiCGStab);
  the report shows the linear iterations and the multigrid setup time.
- **GPU.** PTX entries `amg_*_f32`, `amg_*_f64`, `blas_axpy`, `blas_xpby`;
  the CUDA backend runs the V-cycle with all vectors resident.
- **hypre reference (D20).** Cargo feature `hypre`: a C shim over hypre's
  public interface (BoomerAMG with PCG), `tools/build_hypre.sh` builds hypre
  from source (CPU or CUDA); a CI job builds it and compares.
- **`cfdcraft poisson`.** Mesh sweeps with several solvers and backends;
  report of iterations, setup and solve times, ns per cell and iteration,
  hierarchy, differences to hypre, bitwise comparison across backends;
  checks of iteration growth and agreement.
- **Verification.** Stage-3 tests (iterations nearly independent of the
  mesh size on all mesh types, mixed precision equals float64,
  SPD V-cycle, bitwise reproducibility, agreement with hypre), five Poisson
  case files and one scalar case, gallery 10–11, recipes 15–17, the
  explanation page `docs/verification/stage3.md`, `tools/stage3_colab.ipynb`.

## 0.2 (2026-11-14) – stage 2: operators on a scalar equation

New crate `cfdcraft-fv` (features `hdf5`, default on, and `cuda`) and new
kernel bodies in `cfdcraft-kernels` (`fv`, `blas`), one source for the CPU
and the GPU:

- **Least-squares gradient (D24).** Weighted least squares over the face
  neighbours; weights per cell-face entry precomputed on the host;
  fixed-value, fixed-gradient and symmetry boundaries keep linear fields
  exact on every mesh type.
- **Convection (D25).** First-order upwind and linear upwind as a deferred
  correction; Venkatakrishnan and Venkatakrishnan–Wang limiters.
- **Diffusion (D26).** Over-relaxed non-orthogonal correction, optionally
  limited.
- **Assembly (D31, D37).** Face kernel for the flux coefficients, cell
  kernel that gathers them into a sliced-ELL row (C = 32, σ = 1), the
  right-hand side and the Jacobi diagonal; implicit under-relaxation. No
  atomics; results independent of the thread count.
- **BiCGStab (D36)** with Jacobi preconditioning; dot products in fixed
  blocks of 256, the same order on the CPU and the GPU.
- **Backends.** CPU (rayon) and NVIDIA GPU (`cuda`: the bodies as PTX entry
  points `fv_*`, `blas_*` in `gpu/ptx-kernels`, all fields resident on the
  GPU).
- **`cfdcraft scalar`.** Case file with `[problem]` (manufactured solutions,
  boundary profiles), `[solver]`, `[run]` (backends compared with each
  other), `[verification]` (refinement sweeps, expected order, bounds),
  `[[variants]]` and `[output]` (VTKHDF fields, sampling line). Reports in
  Markdown and JSON with time per phase and ns per cell and iteration (R6),
  convergence and profile plots as SVG.
- **Verification.** Tests for exact gradients, no new extrema, order of
  accuracy 1.8–2.2 on hex, tet, poly and mixed meshes in 2D and 3D
  (release builds, R7), bitwise reproducibility; seven scalar case files,
  gallery examples 8 and 9, recipes 12–14, the explanation page
  `docs/verification/stage2.md`.
- **CI.** Every scalar case except the GPU one runs as a test (R5);
  `tools/check_kernels.sh` refuses atomics, f32 and fused or approximate f64
  arithmetic in the kernels and the PTX (R10).
- **GPU check.** `tools/stage2_colab.ipynb` and `tools/stage2_cloud.sh`
  (CPU against two GPU runs on 0.9 million polyhedra). Result on a Tesla
  T4 (Colab): both GPU runs, the four oblique-step schemes and the sweep
  are bitwise equal to the CPU; 7.1 s against 158.8 s on the 2-thread
  Colab VM.

---

## 0.1 (2026-11-14) – stage 1: mesh

- Release workflow (`.github/workflows/release.yml`): creates the tag
  `v<version>` on `main` and a GitHub release with the changelog section of
  that version as notes; dry-run by default (recipe 6).

New crate `cfdcraft-mesh` (features `hdf5` and `metis`, both on by default;
`mpi` for the MPI halo test):

- **Mesh structure (D23).** Faces as node lists, `owner` and `neighbour`
  per face, the faces of each cell with their orientation; patches and cell
  zones. Order: inner faces, processor faces, boundary faces by patch;
  cells inner, owned cells touching the halo, halo cells by source rank,
  each group in reverse Cuthill-McKee order (R1). u64 global ids for every
  cell, face and node in input order (R2). `validate()` checks index ranges,
  ordering and topological closure (every edge used twice).
- **Geometry.** Face area vectors and centroids from a triangle fan about
  the node average, cell volumes and centroids from tetrahedra on the same
  fan: Σ S = 0 to round-off for every cell, volumes add up to the domain
  volume.
- **Quality report.** Closure, total against enclosed volume, patch areas,
  non-orthogonality, skewness, face warp, aspect ratio, volume ratio,
  inverted faces, non-positive volumes; Markdown and JSON.
- **Gmsh reader (D33).** MSH 4.1 ASCII and binary (either endianness),
  tetrahedra, hexahedra, prisms, pyramids; physical surfaces as patches,
  physical volumes as cell zones, periodic links; refuses higher-order
  elements and other versions with a clear message.
- **Test mesher (D34).** Box meshes of hexahedra, tetrahedra (Kuhn split),
  polyhedra (dual of the tetrahedra), mixed cells; 3D or 2D one cell thick
  (D3); exponential grading per direction, one- or two-sided; nested
  refinement levels; node perturbation by the counter-based generator,
  boundary nodes kept on their faces. Conversion of any tetrahedral mesh into
  polyhedra with the exact volume (also on curved walls).
- **Renumbering.** Reverse Cuthill-McKee and shuffle, with the ids kept.
- **Partitioning (D32).** METIS (vendored) and recursive coordinate
  bisection; imbalance, cut faces and neighbours in the report.
- **Decomposition and halo (D38, R3).** Local meshes with halo depth as a
  parameter; per-neighbour send lists and contiguous receive ranges;
  several fields per message; host-staged transport now, GPU-aware selected
  at run time later on the same plan. Halo geometry filled by exchange.
- **VTKHDF writer (D41).** VTKHDF 2.4 `PartitionedDataSetCollection`:
  block `internal` with cell fields and one block per patch; standard VTK
  cells and `VTK_POLYHEDRON`; global-id order. HDF5 built from source and
  linked statically (needs CMake).
- **Point location and tracking (D47).** Tetrahedral-decomposition walk
  from cell to cell through the faces, with a bin search where the walk
  stops at a non-convex boundary; reports where a track leaves the mesh or
  the rank.

Tests: closure, volume and patch areas of 16 generated variants; Gmsh ASCII
against binary, zones, prisms, pyramids, periodic links, refusals; VTKHDF of
all types; point location on all types; decomposition with METIS and RCB on
1, 2, 3 and 5 parts and 1 and 2 halo layers; **golden values (R4)** for 32
cases of mesh type × dimension × perturbation × refinement × 1, 2 and 4
ranks, checking that owned volumes sum to the total and that every halo
cell has exactly one owner; the MPI halo exchange on 2 and 4 ranks in CI.

CLI: `cfdcraft mesh CASE.toml | FILE.msh` builds or reads, renumbers,
checks, partitions and writes the report and `mesh.vtkhdf`. G1:
`mesh.type = "poly" | "tet" | "hex" | "mixed"` runs the benchmark on
stage-1 meshes.

Tools: `tools/check_vtkhdf.py` (checks files with VTK's reader),
`tools/render_mesh.py` (gallery pictures), `tools/gmsh_test_meshes.py`
(the Gmsh test files), `examples/meshes/tjunction.geo`.

Documentation: gallery examples 3–7, mesh cases and keys in the case
catalogue, recipes 7–11, the G1 stage-1 rerun, implementation notes in
D23, D32, D38, D41 and D47, new open item O15.

Findings:

- Polyhedral boundary faces merged across a curved wall lose volume
  (0.34 % on the T-junction); keeping one planar piece per triangle there
  makes the dual exactly as large as the tetrahedral mesh.
- VTK computes the volume of a polyhedron as that of the convex hull of its
  nodes, so the VTKHDF check computes it from VTK's faces instead.
- **G1 on polyhedra (shared cloud machine, 2 million cells):** bitwise
  equal; the single-source bodies are 25–42 % faster on cell gather, equal
  on SpMV-CSR, 4–30 % slower on face flux, as on a hexahedral box of the
  same size. Shuffled order costs a factor 3 to 4. Quiet-machine and GPU
  reruns are open (O15).

---

## 0.0 (2026-11-14) – stage 0: foundations and gate G1

Workspace (`Cargo.toml`, Rust 2024 edition, minimum Rust 1.89):

- `cfdcraft-core` (`no_std`): the generic floating-point type `Real`
  (f64 reference path, f32 possible, D7), 32-bit local and 64-bit global
  indices, pairwise summation in a fixed order (D31), counter-based random
  numbers keyed by identity, step and event (D47).
- `cfdcraft-kernels` (`no_std`): the kernel bodies, one source for CPU and
  GPU: face flux, cell gather, SpMV in CSR and SELL-C-σ. Gather-based, no
  atomics. Feature `checked-indexing` for bounds-checked loads.
- `cfdcraft-parallel`: `Communicator` trait with `SerialComm` and `MpiComm`
  (feature `mpi`, system MPI library, D9); reproducible all-rank sums;
  multi-threaded sums bitwise equal to the sequential sum for any number of
  threads.
- `cfdcraft-bench`: gate G1. Synthetic periodic face-based mesh (natural or
  shuffled numbering), CSR and SELL-C-σ matrices, backends `cpu-serial`,
  `cpu-threads`, `cpu-native`, and with feature `cuda` `cuda-reference`
  (hand-written CUDA C, compiled by NVRTC at run time) and `cuda-rust` (the
  kernel bodies compiled to PTX). TOML case files, JSON, Markdown and SVG
  reports, G1 verdict per candidate/reference pair.
- `cfdcraft-cli`: the `cfdcraft` command with `info`, `bench g1` and
  `bench backends`.
- `gpu/ptx-kernels`: the kernel bodies as PTX entry points, built with the
  nightly toolchain (`build.sh`).

Continuous integration (`.github/workflows/ci.yml`, D22): format, clippy and
docs; tests in debug and release on Linux, macOS and Windows; minimum Rust
version; MPI smoke test with 2 and 4 ranks on Linux and macOS (O12); PTX
build of the single-source kernels; G1 smoke run on every platform.

Cloud GPU runs: `tools/g1_cloud.sh` (installs Rust, builds PTX and CLI for
the GPU found, runs G1 twice) and `tools/g1_colab.ipynb` (the same on a free
Google Colab T4).

CLI: `cfdcraft info` reports the communicator and the rank count.

Documentation: [gallery](docs/gallery.md), [case catalogue](docs/cases.md),
[recipes](docs/recipes.md), [G1 procedure and results](docs/benchmarks/g1.md).

Findings:

- **G1 GPU verdict, Tesla T4 (Google Colab), 10 million cells: PASS.** The
  single-source Rust kernels, compiled to PTX, run within 1.0 % of the
  hand-written CUDA C kernels on all four kernels (220–255 GB/s of 320 GB/s),
  in two runs, with results bitwise equal to the CPU. D8: one kernel source;
  D11 stands. The single-source GPU route needs nightly Rust (O14).
- **G1 CPU verdict, Intel i9-9900K (macOS), 10 million cells: PASS** in two
  runs: largest time ratio single-source against hand-written Rust loops
  1.018, results bitwise equal, 20–26 GB/s.
- The rustc NVPTX target supports compute capability 7.0 (Volta) and newer
  only; G1 needs such a GPU for the single-source route (O13). The
  maintainer's GTX 960 (compute capability 5.2) runs the CPU verdict and the
  CUDA reference (`g1_gpu_reference_4m.toml`, recipe 2); the GPU verdict is
  taken on a free Colab T4 (recipe 3).
- On the CPU the single-source bodies match hand-written loops once the
  parallel loop uses large enough work packets (64 Ki outputs): with 4 Ki the
  face-flux kernel was about 30 % slower on 10 million cells.
- A branching max/min in the upwind flux cost bandwidth on random mass-flux
  signs; `Real::max_of`/`min_of` are now branch-free.

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
  SST −7 %).

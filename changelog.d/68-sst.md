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

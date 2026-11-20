- **ζ-f and ζ-f₀ converge with the mesh at the wall (stage 6, maintainer
  decision 2026-11-20).** The sink of the k equation is now
  max(k/T_K, 2νk/d²), with d the wall distance and T_K the time scale
  with its Kolmogorov bound. It equals ε at the wall, so k falls like y²,
  and it is k/T_K away from walls, which keeps k positive at a laminar
  leading edge. The sink k/T_K of 6d2 let k fall like y only, and the
  wall value of ε then grew with every refinement: the converged channel
  of ζ-f₀ went from U_b⁺ 17.79 to 18.39 over four meshes. Now it gives
  17.34, 17.54, 17.58, 17.58 (+0.2 % from the DNS on 768 cells across).
  ζ-f₀ takes C_ε1 = 1.4 (1 + 0.009/ζ) in place of the printed
  1.4 (1 + 0.012/√ζ), which gives U_b⁺ −11 % with a consistent sink. The
  flat plate (C_f +1.8, +0.2, −0.2 % from White on three meshes, GCI 0.2 %) and the
  step (x_r/H 6.42) stay in their bands. The impinging jet converges with
  the mesh, but at Nu₀ = 196, +38 % from the experiment (the
  stagnation-point excess of k; finding). The ζ-f plate now converges
  (C_f +0.05 % from White on the middle mesh; before, it stopped at the
  iteration limit at +24 %).
- **New keys:** `turbulence.k_sink` (`"wall"`, the default; `"epsilon"`,
  the papers' −ε; `"kolmogorov"`, the old ζ-f₀ default) and
  `[turbulence.constants] c_epsilon1_power` (1, or ½ for the printed
  C_ε1″/√ζ of Popovac 2023).
- **Fix:** the ambient (sustaining) source of k in ζ-f and ζ-f₀ is the
  sink at the ambient state, so the free stream holds its k for every sink
  form.
- **The channel runs of the GCI study use the coupled algorithm**, at
  most 3 000 iterations each: ζ-f₀ converges to 10⁻⁹ in 436 to 883
  iterations instead of up to 35 000 with SIMPLEC, the k-ε models in
  under 700, ζ-f in under 2 600, with the same U_b⁺. SST on the finest
  mesh and EB on the two finest stop at the limit with a residual near
  10⁻⁷, within 0.3 % of their converged values (finding; the docs keep
  the SIMPLEC values of #128). The plate runs of ζ-f and ζ-f₀ use it too
  (178 to 496 iterations).
- **The channel case files of ζ-f and ζ-f₀ use the coupled algorithm**
  (about 230 and 370 iterations instead of 2 100 and 3 700).
- **The ζ-f plate converges** (`require_convergence = true`): C_f +2.0 %
  from White on 64 × 48 cells, +0.03 % on 128 × 96 (before +5.1 % with
  the k floor of the leading edge).
- **The inlet table of the impinging jet** (`impinging_jet_profile.csv`)
  is made again from the pipe flow with the new sink.
- **Gallery 141 and recipe 141**: the sink of k with two field images
  (the ζ-f plate near its leading edge, the wall layer of the jet).
- **Finding:** in the test of a turbulence patch carried by the rotor
  (`tests/sliding.rs`), the imprinted run of ζ-f₀ at a Courant number of
  1.7 drives ζ to its cap and diverges with any sink of k that is ε at the
  wall (the new default and the papers' −ε); k/T_K held there. The
  imprinted run of the test now takes half the time step (Courant 0.85).

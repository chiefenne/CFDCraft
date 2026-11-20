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
  flat plate (C_f +1.8, +0.15, −0.3 % from White on three meshes) and the
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
- **The channel runs of the GCI study use the coupled algorithm.** They
  converge to 10⁻⁹ in a few hundred iterations instead of up to 91 000
  with SIMPLE, with the same U_b⁺ to five digits.

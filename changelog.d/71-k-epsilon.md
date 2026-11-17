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

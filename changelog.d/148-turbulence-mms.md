- **Manufactured solutions of every turbulence model (stage 6
  specification 9.2).** `turbulence::mms` solves the turbulence equations
  of SST, the standard and the realizable k-ε, k-ε-ζ-f and the
  elliptic-blending k-ε alone in a frozen cellular flow with smooth
  manufactured k, ω/ε, ζ/φ, f/α and wall distance and every source term
  active; the L2 errors converge at order 1.8–2.3 on hexahedra, polyhedra
  and triangular prisms. New: `Flow::turbulence_iteration`,
  `Turbulence::set_added_zeta_sources` (with an added right-hand side of
  the f/α equation in the kernel) and `Turbulence::set_zeta_fields`.
- **Fix: realizable k-ε C_μ with the deviatoric strain.** W̃ of Shih et
  al. now comes from the trace-free strain; with the discrete trace it
  jumped across its range where the strain vanishes.

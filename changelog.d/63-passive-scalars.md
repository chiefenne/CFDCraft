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

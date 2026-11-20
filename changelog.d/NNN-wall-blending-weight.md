- **The viscous weight of the wall treatment as a wall field (spec 2.3,
  D72).** `turbulence-wall` of `[output] derived` and `foehn post derive`
  gains the boundary dataset `blending`: 1 where the wall is resolved
  (`low-re`, the sublayer of `standard`), 0 where the log law holds (the
  log layer, `scalable`), in between 1 − tanh((y⁺/5)²) of `automatic`,
  e^{−Γ} of `compound` and Γ = dU⁺/dy⁺ of `compound-2023`. It shows which
  walls of a mesh are resolved. The restart file stores it as
  `turb_wall_blending`; from an older restart file `turbulence-wall` is
  refused with CASE005 (continue the run for one iteration).

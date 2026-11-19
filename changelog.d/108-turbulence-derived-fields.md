- **Derived fields of the turbulence models (D72, gallery and recipe 108).**
  `[output] derived` and `foehn post derive` take `turbulence-scales`
  (ε, ω, the time scale of the model, k^{3/2}/ε, ν_t/ν for every model),
  `reynolds-stress` (the Boussinesq estimate ⅔ k δᵢⱼ − 2 ν_t Sᵢⱼ, six
  components), `wall-normal` (∇α/|∇α| of the elliptic-blending k-ε, for the
  GGDH of stage 7) and `turbulence-wall` (u_τ of the wall treatment and y*
  on the walls). New plot `tools/plot/turbulence_fields.py`; tests
  `tests/turbulence_derived.rs`.

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

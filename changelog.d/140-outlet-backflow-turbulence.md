- **Turbulence values of the backflow at pressure outlets (BC028, gallery
  and recipe 140).** A `static-pressure-outlet` or `total-pressure-outlet`
  that gives `k` with `omega`/`epsilon`, `turbulence_intensity` with
  `turbulence_length_scale` or `viscosity_ratio`, or `zeta` now uses them on
  the faces where the flow enters (fixed values), decided per face in every
  outer iteration; the faces where it leaves stay zero gradient. Before,
  the entering fluid always carried the values of the cell and BC028
  warned that the values were unused. New case
  `backstep_outlet_backflow.toml`, plot `tools/plot/outlet_backflow.py`.
- **Fix: restarts of turbulent runs with pressure boundaries.** The
  modified pressure p̂ = p + ⅔ρk prescribed at the pressure boundaries came
  from the start fields until the first turbulence solve after a restart,
  so a continued run was not bit-identical to the uninterrupted one. It is
  now set from the loaded k.

- **PISO corrector of SIMPLEC type: stable on thin wedge cells at an axis at
  high viscosity (gallery 103).** The pressure corrections of PISO use the
  SIMPLEC coefficient D^C = V/(a_P/α − Σ|a_N|) with the consistent
  neighbour correction Σ|a_N| (u'_N − u'_P)/a_P between two correctors,
  instead of Issa's D = V/a_P with H(u')/a_P. Where a_P is almost all
  coupling to the neighbours (the wedge cells about an axis at a fine
  pitch and high viscosity, a_P about 50 times the time term), Issa's form
  let the pressure grow in every outer iteration: the disc impeller in a pot
  diverged at ν = 1e-4 m²/s. It now runs at that viscosity
  (`examples/studies/revolved/disc_pot.toml` is back at ν = 1e-4). The
  converged step is unchanged. New case
  `examples/cases/flow/pot_solid_body_viscous.toml`: solid-body rotation at
  ν = 1e-2 m²/s, L2 error 0.11 % (2.0 before).

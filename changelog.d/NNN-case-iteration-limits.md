- **Iteration limits of the steady cases.** No case file allows 10 000
  iterations any more. `cylinder_re20_gci.toml` runs the coupled algorithm
  (51, 90 and 225 outer iterations instead of 960, 3 007 and 9 876 with
  SIMPLEC at the momentum relaxation 0.5; about a minute instead of 12 on
  four cores) with `max_iterations = 1000`; C_D and Δp agree to 10⁻⁵, the
  lift moves by up to 0.8 % with the pseudo-time step (gallery 25).
  `cavity_re1000.toml` and `cavity_re1000_poly3d.toml` keep SIMPLEC
  (1 397 and 1 240 iterations) with `max_iterations = 2000`;
  `tmr_plate_sst.toml` keeps SIMPLEC on level 1 (635) with
  `max_iterations = 2000`, and its four-level sweep runs with
  `solver.algorithm = "coupled"` (622 and 411 iterations on levels 2 and 3
  instead of 2 292 and more than 4 400).

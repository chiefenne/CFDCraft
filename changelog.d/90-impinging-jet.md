- **Impinging jet at H/D = 2 (stage 6 part 6g, specification 9.6).** New
  cases `jet_pipe_profile.toml` (the developed pipe flow at Re = 23 750
  on a periodic wedge) and `impinging_jet.toml` (the round jet from that
  pipe on a heated plate, axisymmetric wedge, a passive temperature with
  a uniform wall flux). New study `examples/studies/impinging_jet.py`:
  six models on three meshes, Nu₀ with its GCI, Nu(r/D) with its
  secondary maximum, the stagnation-point anomaly with and without the
  production limiters, and the wall-jet velocity profiles; it compares
  with the experiment once the data are in
  `examples/studies/impinging_jet/data`. On 41 600 cells (with the
  fix of #91) Nu₀ at the axis: ζ-f 173.8 (GCI 0.12 %), SST 141.2
  (2.9 %), EB 249.2 (0.19 %); Nu at r/D = 0.1, the second probe `r01`:
  ζ-f 173.9, EB 205.0, SST 141.1 (GCIs 0.03 to 0.40 %), ζ-f₀ 151.7
  (GCI 5.3 %). The k-ε without a limiter has the anomaly (Nu₀ 211), its
  Durbin default removes 71 % of the excess k; SST has none (its ν_t
  bound holds on the axis). Gallery and recipe 90.
- **Inlet tables with the turbulence of the profile (specification 7.1).**
  The header of a `profile = "table"` file names its columns: `u_n` or
  `u_x`, `u_y`, `u_z`, and `k` with `epsilon` or `omega`, and `zeta`.
  The turbulence columns set the inlet's k, ε or ω and ζ face by face; a
  boundary that gives turbulence values as well is refused (BC029), an
  unknown column too (BC019). `python3 -m foehncraft.profiles LINE.csv
  TABLE.csv` makes the table from a line of a periodic pipe or channel
  run (new module `foehncraft.profiles`).
- **Wall probes report the passive scalars.** Each `[[output.wall_probes]]`
  gives the wall value and the diffusive flux of every scalar there, in
  the log and the report, and as GCI quantities of a level sweep
  (`wall_probes.<name>.<scalar>_wall`, `_flux`).
- **Fix: the thermal wall law of the k-ε family and ζ-f/EB takes
  u* = C_μ^¼ k_P^½ (specification 5.2), not u_τ.** u_τ vanishes with the
  wall shear, so the heat transfer of a stagnation point or a
  reattachment line went to zero (the k-ε Nu₀ of the jet was 37 instead
  of 156). SST keeps its blended u_τ. The heated channel moves by −0.1 to
  −1.2 points (still within 2.9 % of the DNS at Pr 0.71), the heated pipe
  by +0.3; the flat plate is +8 to +9 % for 0.4 ≤ x ≤ 1, and +14 % at
  x = 0.29 just after the k-ε start of the turbulent layer, where k
  overshoots (the test window moves from 0.25 to 0.4).

- **Backward-facing step: three models, two field pictures (gallery 74).**
  `examples/studies/backstep/{keps,zf0,sst}.toml` run the step of Driver &
  Seegmiller on 32 800 cells with k-ε and scalable wall functions,
  k-ε-ζ-f₀ and k-ω SST: x_r/H = 5.34, 6.38 and 5.79 against the
  experiment's 6.26. The gallery shows the velocity with streamlines and
  the turbulent kinetic energy, one panel per model.
- **Figure scripts for 2D fields.** `tools/plot/backstep_paraview.py` makes
  the two pictures in ParaView's Python (`pvpython`; run in ParaView 5.11
  and 6.2), `tools/plot/backstep_chart.py` with matplotlib and h5py
  (recipe 74).

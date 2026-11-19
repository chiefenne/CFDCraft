- **SST `automatic` wall treatment holds C_f in the buffer layer (6g).**
  The near-wall ω is blended with the exponent 1.2 and imposed, and the
  first cell takes the equilibrium production of Reichardt's law, both with
  the weight tanh((y⁺/5)²). On the TMR plate C_f now lies within −2.6 to
  +2.4 % of the reference for first-cell y⁺ 0.24 to 85 (before: up to
  +15.6 % at y⁺ 4.6). TURB001 warns for `automatic` only above y⁺ 300.
  The sweep is `examples/studies/sst_wall_sweep.py`; gallery 68 shows it.
- **SST `low-re`: the wall value of ω checked.** A sweep of the factor in
  ω_w = F · 6ν/(β₁ d₁²) from 1 to 100 on the plate and the channel keeps
  F = 10, which has the smallest mesh dependence (verification manual).

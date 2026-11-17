- **Backward-facing step: settled results and half the iterations.** The
  step cases start with the momentum under-relaxation 0.5 and switch to
  0.8 at iteration 100 (`[[schedule]]`), which halves the iterations to a
  settled reattachment length: about 600 on 8 200 cells and 1 750 on
  32 800 cells for k-ε, k-ε-ζ-f₀ and SST. The fine study cases run 3 000
  iterations, `backstep_zeta_f0` 1 500 instead of 4 000. Run to the end,
  SST reattaches at x_r/H = 6.80 (+9 %), not 5.79; k-ε at 5.36 and
  ζ-f₀ at 6.35 on the fine mesh. The figures of gallery 74 are redrawn
  from the settled runs; recipe 74 shows how to check when a result
  has settled.

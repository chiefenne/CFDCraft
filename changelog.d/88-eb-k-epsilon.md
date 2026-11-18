- **Elliptic-blending k-ε (stage 6 part 6e, D133).** New model
  `elliptic-blending-k-epsilon`, the BL-v²/k model of Billard & Laurence
  (2012) with the equations and constants of the paper: k, ε, φ = v²/k
  and the blending parameter α (α − L²∇²α = 1, α = 0 at walls), the E
  term in the k equation and the variable C*_ε2, ν/2 in the molecular
  diffusion, ε_w = νk/y²; wall treatments `compound` (default) and
  `low-re`. It runs on the fields, solvers and multigrid of k-ε-ζ-f (φ in
  the place of ζ, α in the place of f); new element-wise kernels
  `squares_add` and `scale_vector` (also in the PTX module). Channel C_f
  at 97.8 / 99.0 / 99.1 % of the DNS at Re_τ 180 / 395 / 590 (the paper:
  100 / 98 / 100 %); TMR plate C_f(0.97) +0.5 % from White; backward-facing
  step x_r/H = 7.59, +21 % (finding: the variable C*_ε2;
  `constants = { c_epsilon4 = 1.83 }` gives 6.24 and a plate 4.9 % high).
  Output fields `phi` and `alpha`, log columns `ResPhi` and `ResAlpha`;
  the key `zeta` sets φ. New cases `channel_eb.toml`, `tmr_plate_eb.toml`,
  `backstep_eb.toml`; specification 3.4 corrected after the paper.
  Gallery and recipe 88.

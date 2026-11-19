- **Quadratic constitutive relation of SST (6h, D170, gallery 115,
  recipe 115).** `turbulence.qcr = "qcr2000"` adds the non-linear stress
  of Spalart (2000) to the momentum equations, and `"qcr2013"` adds the
  isotropic term of Mani, Babcock, Winkler & Spalart (2013) (constants
  `c_cr1`, `c_cr2`). The stress goes onto the momentum faces as an
  explicit stress, the path the later algebraic stress model reuses. The
  new case `square_duct_qcr.toml` (Re_b 3 500) shows the eight corner
  vortices with a largest secondary velocity of 1.56 % of U_b. Linear
  SST gives none.

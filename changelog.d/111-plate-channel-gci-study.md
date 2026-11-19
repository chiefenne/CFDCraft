- **Grid convergence of the turbulence models (spec 9.3, 9.4).**
  `examples/studies/plate_channel_gci.py` runs all six models on three
  meshes of the TMR flat plate and the Re_τ 395 channel and reports the
  observed order, the extrapolated value and GCI_fine (Celik 2008) of
  C_f, C_D and U_b⁺ against White, the TMR and the DNS. Findings on the
  stage-6 verification page: no model is in the asymptotic range on the
  plate meshes (SST within 0.14 % of the TMR value, GCI 0.47 %), ζ-f does
  not converge on the finer plates, the realizable k-ε with scalable wall
  functions moves with every refinement of the channel.

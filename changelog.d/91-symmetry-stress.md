- **Fix: the full viscous stress on symmetry planes, slip walls and wedge
  planes (stage 6 part 6g).** With a turbulence model the momentum rows
  take the transposed part of the stress, ∇·[μ_eff (∇u)ᵀ], on the
  interior faces; on symmetry faces it was missing. There u_n = 0, so
  that part is μ_eff |S| ∂u_n/∂n along the normal, and the normal stress
  on the face is 2μ_eff ∂u_n/∂n (implicit, as before with the factor 1).
  On the two planes of a wedge it is half of the hoop stress
  −2μ_eff v/r²: the cells next to the axis of the impinging jet had v/r
  too high and Nu too low, more so on finer meshes (the open item of
  #90). Laminar runs do not change. New test
  `tests/symmetry_stress.rs`: a linear flow with constant stress balances
  in every cell, on a wedge and on a symmetry plane (the residual was
  4 % and 65 % of a face flux).

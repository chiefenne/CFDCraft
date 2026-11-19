- **Turbulence on rotating machinery: the moving reference frame and both
  sliding meshes (gallery 100, recipe 100).** Every turbulence model now runs
  with `[[sliding]]` (it was refused, `SLD004`) and is checked with
  `[[mrf]]`. With conformal stepping the old levels of k, ε, ω, ζ and f move
  with the rotor and each step starts from the moved fields; with imprinting
  the turbulence state is carried through the rebuild of every step. The wall
  distance of the turning rotor is kept per rigid part (one bounding volume
  hierarchy each) and only the cells nearer to the edge of their part than to
  its walls search the other part in each step, with a capped query that
  starts from the last nearest face. Several independent `[[mrf]]` zones are
  tested (two zones turning in opposite senses).
- **Fix:** the pressure multigrid built a hierarchy of NaN when a coarse row
  had a non-positive lumped diagonal in the prolongator smoothing (`GEN001` in
  the first iteration of the turbulent stirred vessel with the frame); such a
  row keeps its own diagonal now.

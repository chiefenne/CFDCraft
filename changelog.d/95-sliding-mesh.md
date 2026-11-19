- **Sliding mesh with conformal stepping (D56, D182, stage 6 part 6g, part
  6).** `[[sliding]]` turns a cell zone of a ring mesh by whole sectors in
  every time step (Δt = shift · 2π/(N |Ω|), refused with the right value
  otherwise, `SLD003`), so interface cells meet new neighbours without
  interpolation and the mesh stays conformal. The solved velocity is the
  absolute one, momentum and scalars are convected with the flux relative to
  the moving mesh (the frame flux of D55, without the Coriolis source), and
  only the time term sees the motion: the old time levels of a rotor cell are
  those of the cell that stood at its place. Laminar flow and passive
  scalars; turbulence models, pressure boundaries on the rotor and steady
  cases are refused (`SLD001` to `SLD004`). New: `[initial] from_analytic`
  for transient runs, the case `taylor_couette_sliding.toml` (L2 error
  1.35 % on 512 cells, observed order 2.2 to 2.3 from 128 to 8 192 cells),
  spec 9.7, verification page, gallery and recipe 95. The variant with a free
  pitch ratio (imprinting) follows.

- **Moving reference frame (D55) with laminar Taylor–Couette flow (stage 6
  part 6g, part 5).** `[[mrf]]` tables turn a cell zone steadily about an
  axis: the solved velocity is the absolute one, momentum, turbulence
  variables and scalars are convected with the flux relative to the frame,
  the rotating cells get the Coriolis source, continuity and the pressure
  correction stay on the absolute flux. New: `[mesh.annulus]` (the gap
  between two cylinders with an optional `rotor` zone), `[analytic]
  type = "taylor-couette"`, the diagnostics `MRF001` to `MRF003`, the case
  `taylor_couette_mrf.toml` (L2 error 0.69 % on 1 024 cells, observed order
  2.3 and 2.4, the same in every frame), spec 9.7, verification page, gallery
  and recipe 94. Steady only; the sliding mesh follows (D182).

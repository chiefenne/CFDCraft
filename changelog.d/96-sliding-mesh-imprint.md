- **Sliding mesh by imprinting (D182 variant B, stage 6 part 6g, part 7).**
  A rotor and a stator with different numbers of sectors meet on an
  interface that is imprinted at the rotor angle of every time step, so the
  time step is free. `[mesh.annulus] n_theta_stator` makes the two rings
  with the boundaries `rotor_interface` and `stator_interface`;
  `[[sliding]] rotor_interface / stator_interface` names them. Every
  intersection of a rotor face with a stator face is an ordinary internal
  face (`foehncraft_mesh::imprint`: polygons clipped in the unrolled
  surface of revolution), the cells close, mass and momentum flux are
  conserved to round-off, and the solver is built on that mesh in every step
  with the state, the old time levels and the face fluxes carried over. The
  scheme is fully implicit in the step (`extrapolation = false`, at least
  two outer iterations, `SLD003`); `SLD005` and `SLD006` report an interface
  that cannot be imprinted. Taylor–Couette flow from rest with 32 rotor and
  48 stator sectors: L2 error 6.74, 1.29 and 0.25 % on 160, 640 and 2 560
  cells (order 2.4); equal pitch agrees with conformal stepping (1.350 % and
  1.354 %); solid-body rotation is kept; a passive scalar ends at 0.5 %.
  Laminar flow and passive scalars, one rotor, a fixed time step, no
  restart, CPU. Also: the VTKHDF writer writes cells whose faces do not
  close on the nodes of their shape as polyhedra, and probes in a moving
  rotor find their cell again in every step. New case
  `taylor_couette_imprint.toml`, gallery and recipe 96, spec 9.7.

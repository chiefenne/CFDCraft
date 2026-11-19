- **Several sliding rotors and the transient moving reference frame (D185,
  gallery 102, recipe 102).** One `[[sliding]]` table per rotor, each with
  its own zone, speed, axis, origin and interface to the stator: all
  interfaces are imprinted together in every step at the angles of their
  rotors (variant B of D182; conformal stepping stays one rotor, an
  interface between two rotors is refused). The exact wall distance of
  turbulence models follows every rotor. The log headings of rotor columns
  name the rotor. `[mesh.annulus] outer_rotor_layers` and `n_theta_outer`
  make a second rotor ring at the outer cylinder; `examples/meshes/twin_mixer.py`
  makes a vessel with two stirrers. `[[mrf]]` now runs in transient cases
  too (the transient frozen rotor), also next to sliding rotors in other
  zones; `MRF001` refuses a zone that a frame and a sliding mesh both turn
  (it refused every transient frame before). New cases
  `examples/cases/flow/couette_two_rotors.toml`,
  `examples/studies/twin_mixer/sliding.toml`,
  `examples/studies/mixer/mrf_transient.toml`.

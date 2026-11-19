- **Walls of a turning zone push the fluid (fix).** A wall in a zone of
  `[[mrf]]` or `[[sliding]]` kept only the tangential part of its `velocity`,
  the rule for walls that stand (D53, `BC010`), so a blade, which moves
  normal to itself, slid along the fluid instead of pushing it. Such a wall
  now moves normal to itself with the zone (Ω × r) and along itself with the
  given velocity: absolute flux ρ **S**·**V**, relative flux zero, as D55
  states. Cylinders that turn about their axis (gallery 94 to 96) are not
  affected. New test: a box turning with its frame carries its fluid (L2
  1.8·10⁻² before, 1.3·10⁻³ after).
- **A stirred vessel with blades and baffles, and field figures of rotating
  machinery (gallery 98, recipe 98).** `examples/meshes/mixer.py` meshes a
  2D baffled vessel with a four-blade rotor in Gmsh (conformal for the
  frame, split 64/96 for the sliding mesh); `examples/studies/mixer/` has the
  steady frame case (Re 100) and the sliding-mesh start-up over three
  revolutions (Re 400); `examples/studies/mixer.py` runs both and draws the
  figures. `tools/plot/rotor_fields.py` draws any 2D rotor run: speed with
  streamlines, the speed relative to the blades in the rotor, vorticity,
  pressure, time series and zooms with the cell edges. Gallery 94, 95 and 96
  have field figures now too.

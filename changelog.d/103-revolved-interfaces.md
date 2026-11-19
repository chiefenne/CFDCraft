- **Sliding interfaces of any rotationally symmetric shape (D186, gallery
  103, recipe 103).** One interface per rotor, of any surface of revolution
  about its axis whose meridian is made of straight pieces: a cylinder, a
  disc (an axial gap between blade rows), a cone, or a chain of them such as
  a pot of disc, cylinder and disc around an impeller; several separate
  pieces (discs before and after a rotor row) and discs through the axis
  too. The meridian is found from the mesh; the faces are clipped in the
  unrolled (angle, meridian length) coordinates with a polar chart at the
  axis, and the corners of the pieces are inserted into the neighbouring
  face edges, so every cell of an imprinted mesh closes to round-off on any
  shape. The sliding log line names the shape. New built-in mesh
  `[mesh.revolved]`: a block grid in (r, z) turned about the axis, with
  rotor blocks (own pitch), solid blocks (walls such as an impeller disc),
  cones and wedges at the axis. New cases
  `examples/cases/flow/axial_stage_uniform.toml`,
  `examples/cases/flow/pot_solid_body.toml`,
  `examples/studies/revolved/disc_pot.toml`,
  `examples/studies/revolved/swirl_stage.toml`; new plot tool
  `tools/plot/meridional_fields.py`.

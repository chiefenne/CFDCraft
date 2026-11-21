- **Solid regions and conjugate heat transfer (stage 7d1).** A region of
  kind `solid` (`[regions.<name>] kind = "solid"`, its cell zones and a
  solid material with density, specific heat and conductivity) now joins
  the energy equation: the flow runs on the fluid cells, where the faces
  against a solid become the no-slip wall `interface-<region>`, and the
  temperature runs on all cells as one matrix, with the series conductance
  of the two half cells at every face between materials. The report gives
  the heat flow and mean temperature of each interface; the heat balance
  of the whole domain closes to round-off. The box generator takes
  `zones = [{ name, lower, upper }]` to mark cell zones by centroid.
  Steady laminar cases only for now: transient runs, turbulence, restart
  and passive scalars together with solids are refused (ENE006), as are a
  solid without a complete material (MAT006), solids without `[energy]`
  (REG001) and unknown zones (REG002). Check V7.4 (a composite wall with
  conductivity ratios 1 : 100 : 10⁴, exact to round-off) and a heated wall
  under a channel flow (gallery and recipe 156).

- **Energy equation of the fluid (stage 7a1).** `[energy]` solves the
  temperature of the fluid, conservative in the sensible enthalpy, with
  c_p(T), k(T) and μ(T) from the material or from `[fluid] specific_heat`
  and `conductivity` (μ(T) reaches the momentum equation;
  `temperature_dependent = false` keeps the reference values). Viscous
  heating, volumetric sources and exact solutions. Walls take
  `temperature`, `heat_flux`, `heat_flow` (W) or
  `heat_transfer_coefficient` with `ambient_temperature` in
  `[boundary.<name>.energy]`; inlets take `temperature`. The report has a
  heat balance and the heat flow of every boundary. New diagnostics
  ENE001 to ENE003. Checks V7.1 and V7.2 at second order (gallery and
  recipe 154).

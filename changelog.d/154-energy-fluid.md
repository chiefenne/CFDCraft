- **Energy equation of the fluid (stage 7a1).** `[energy]` solves the
  temperature of the fluid, conservative in the sensible enthalpy, with
  c_p(T), k(T) and μ(T) from the material or from `[fluid] specific_heat`
  and `conductivity` (μ(T) reaches the momentum equation;
  `temperature_dependent = false` keeps the reference values). Viscous
  heating, volumetric sources and exact solutions. Walls take
  `temperature`, `heat_flux`, `heat_flow` (W) or
  `heat_transfer_coefficient` with `ambient_temperature` in
  `[boundary.<name>.energy]`; inlets take `temperature`. The report has a
  heat balance and the heat flow of every boundary. The condition follows
  from the keys (no `type` key). The thermal set-up check gives THM001 to THM005 and THM009 to
  THM011 (undetermined temperature level, heat with no way out,
  conditions on the wrong boundary kind, values out of range, adiabatic
  walls listed, tables without `[energy]`, malformed tables, inlets
  without a temperature); ENE006 for settings this version does not
  have. Checks V7.1 and V7.2 at second order (gallery and
  recipe 154).

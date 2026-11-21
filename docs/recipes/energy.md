# Energy equation and heat transfer

Part of the [Recipes](../recipes.md): the energy equation of the fluid (stage 7, D150), its properties, thermal wall conditions and heat balance. The keys are on the [key page](../keys.md).

| # | Recipe |
|---|---|
| 154 | [Solve the temperature of a flow with temperature-dependent properties](#energy-fluid) |

---

<a id="energy-fluid"></a>
## 154. Solve the temperature of a flow with temperature-dependent properties

**The equation.** An `[energy]` table adds the temperature T of the fluid,
solved after the pressure correction of every outer iteration:

ρ c_p ∂T/∂t + ∇·(ṁ c_p T) − T ∇·(ṁ c_p) − ∇·(k ∇T) = q''' + Φ.

The matrix is in T, and a correction on the right-hand side makes the
converged result conservative in the sensible enthalpy h_s = ∫c_p dT
(zero at 298.15 K), so a c_p that changes with T does not lose heat. The
equation uses the convection schemes, limiters, turbulent Prandtl number
and thermal wall functions of the passive scalars (recipe 63).

```toml
[fluid]
material = "air"               # c_p(T), k(T) and μ(T) from the database
# or the values themselves (numbers or expressions of T):
# density = 1.0
# viscosity = 0.01
# specific_heat = 1000.0       # J/(kg K)
# conductivity = "0.02 + 7e-5 * (T - 300)"   # W/(m K)

[energy]
initial = 300.0                # K, a number or an expression of x, y, z
source = 0.0                   # q''' (W/m³), expression of x, y, z (and t)
viscous_heating = false        # Φ = μ (∇u + ∇uᵀ):∇u
temperature_dependent = true   # refresh c_p, k, μ from T every outer iteration
exact = "..."                  # optional: L2 error and observed order of a sweep

[solver.energy]
relaxation = 1.0               # default 0.9 steady
tolerance = 1e-8
max_iter = 200
```

`[fluid] specific_heat` and `conductivity` replace those of the material.
With `temperature_dependent = true` (the default) μ(T) of the material
also reaches the momentum equation; the density stays constant until the
variable-density step of stage 7 (7c).

**Thermal boundary conditions.** Walls are adiabatic unless a
`[boundary.<name>.energy]` table gives exactly one condition:

```toml
[boundary.hot.energy]
temperature = 320.0            # K, a number or an expression

[boundary.heated.energy]
heat_flux = 500.0              # W/m² into the fluid

[boundary.heater.energy]
heat_flow = 0.5                # W over the whole boundary

[boundary.outer.energy]
heat_transfer_coefficient = 10.0   # W/(m² K)
ambient_temperature = 293.15       # K: q = h (T∞ − T_w)

[boundary.inlet.energy]
temperature = 300.0            # inlets need one (THM011)
```

The condition follows from the keys, so there is no `type` key. Walls and
slip walls take any of them; walls without a table are adiabatic and the
check lists them once (THM005). Outlets keep a zero gradient, symmetry
planes and wedges pass no heat, and periodic pairs carry T across: an
energy table there is refused (THM003), as a flux at an inlet is.

The check before the run (`foehn flow --check`) also refuses:
- two conditions or none in one table, or a table of an unknown boundary
  (THM010), and an inlet without a temperature (THM011);
- a temperature, ambient temperature or start temperature ≤ 0 K, and a
  heat-transfer coefficient ≤ 0 (THM004);
- a steady case with no inflow and no temperature or convection boundary.
  That is THM001 when the temperature level is left open, or THM002 when a
  heat input has no way out;
- a boundary table without `[energy]` (THM009);
- a passive scalar named `temperature`, `[regions]` together with
  `[energy]` (until the conjugate step 7d1), or an orthotropic fluid
  conductivity (ENE006), and a fluid without c_p or k (MAT007).

**The heat balance.** The run report (`flow.md`, `flow.json`,
`result.json`) and the log give, for the final state, the conduction
through the boundaries, the enthalpy carried in and out by the flow, the
sources (q''' and Φ) and the imbalance, and per boundary the mean
temperature, the heat flux (W/m²) and the heat flow (W). A Nusselt number
follows as Nu = q L/(k ΔT) from the wall flux; in the heated cavity of
gallery 154 it is 2.05.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cavity_heated.toml --out output/cavity_heated
target/release/foehn flow examples/cases/flow/couette_viscous_heating.toml --out output/couette
target/release/foehn flow examples/cases/flow/conduction_mms.toml --out output/conduction
```

**Limits of this step.** The property refresh and the enthalpy correction
run on the host: on a GPU every outer iteration copies the temperature
down and the properties up. With a turbulence model the wall functions and
the turbulence sources still use the reference viscosity. The
streamwise-periodic heat-flux mode of the passive scalars is not yet
available for the energy equation.

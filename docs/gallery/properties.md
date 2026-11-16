# Expressions and materials

Part of the [Example gallery](../gallery.md): boundary values, body forces and viscosities given as expressions (D53), fluid properties from the material database and regions of cell zones (D102, D112), and walls with a mass flux (D124 W3).

| # | Example | Stage | Feature |
|---|---|---|---|
| 23 | [Expressions and materials: variable viscosity, blowing and suction, water from the database](#expressions-materials) | 5 | `[constants]`, `[tables]`, expressions, `[materials]`, `[regions]`, `foehn materials`, wall `mass_flux` |

---

<a id="expressions-materials"></a>
## 23. Expressions and materials: variable viscosity, blowing and suction, water from the database

![Left: plane Couette flow with the viscosity μ0 (1 + 4 y) from an expression. Middle: Couette flow with blowing through the lower wall and suction through the upper, Re = 4. Right: the viscosity of water in the material database, with the point of poiseuille_water.toml](../figures/gallery/expressions_materials.svg)

**A viscosity that varies in space.** `couette_viscosity.toml` gives the
viscosity as the string `"mu0 * (1 + a * y)"` with the constants
`mu0 = 0.1` and `a = 4`. The shear stress of plane Couette flow is the
same on every layer, so the velocity follows the viscosity:
u = U ln(1 + a y)/ln(1 + a), steeper near the lower wall where the fluid
is thinner (left; the dashed line is the constant-viscosity profile). The
viscosity is evaluated on every face; the observed order of the velocity
error is 1.83–1.96 on three meshes.

**Walls with a mass flux.** `couette_transpiration.toml` blows
0.2 kg/(m² s) through the lower wall and sucks the same through the upper
(`mass_flux = 0.2` and `-0.2`, D124 W3). The cross flow V = 0.2 m/s
pushes the shear layer against the upper wall:
u = U (e^(Re η) − 1)/(e^Re − 1), Re = V H/ν = 4 (middle). Observed order
1.86–1.96. The same key carries the evaporating mass of a wall film and the
species flux of a catalytic wall in later stages. A closed volume whose
walls add mass is refused before the run (BC026).

**Water from the database.** `poiseuille_water.toml` names
`material = "water"` at `temperature = 313.15` K: the density
992.2 kg/m³ and the viscosity 0.6527 mPa s come from the tables of the
built-in database (IAPWS-95 and the IAPWS 2008 viscosity formulation,
monotone cubic in T; right). `foehn materials` lists the database and
evaluates an entry at any temperature. The run writes the resolved
property program, with its sources, into `case.resolved.toml` and its
hash into every output file.

```
target/release/foehn flow examples/cases/flow/couette_viscosity.toml
target/release/foehn flow examples/cases/flow/couette_transpiration.toml
target/release/foehn flow examples/cases/flow/poiseuille_water.toml
target/release/foehn materials
target/release/foehn materials water --temperature 313.15
# the chart: the middle meshes with a line through the cell centres
for c in couette_viscosity couette_transpiration; do
  target/release/foehn flow examples/cases/flow/$c.toml --set 'verification.levels=[1]' \
      --set 'output.lines=[{name="mid",from=[1.0,0.03125,0.05],to=[1.0,0.96875,0.05],points=16}]'
done
python3 tools/property_chart.py docs/figures/gallery/expressions_materials.svg
```

**Expressions against the built-in profiles.** Three more cases check
the expressions against the runs they replace: `poiseuille_expression.toml`
writes the inlet parabola as `"6 * u_mean * (y / h) * (1 - y / h)"` and
drives half of the flow by a body force expression; the velocity errors
equal those of `poiseuille_hex.toml` on all three meshes.
`taylor_green_expression.toml` prescribes the decaying vortex on the four
sides as expressions of x, y and t, evaluated at each time step (error
5.5e-4 at t = 0.4 s, as with `profile = "analytic"`).
`channel_ramp.toml` ramps the inflow with a table `ramp(t)`.

| Case | What the expression does | Check |
|---|---|---|
| `couette_viscosity.toml` | μ(y) | order 1.83–1.96 |
| `couette_transpiration.toml` | wall mass flux, blowing and suction | order 1.86–1.96 |
| `poiseuille_water.toml` | water at 313.15 K from the database, a region of the zone `fluid` | error 4.0e-4 |
| `poiseuille_expression.toml` | inlet profile and body force | errors of `poiseuille_hex.toml`, order 1.89–1.98 |
| `taylor_green_expression.toml` | velocity on the sides, in x, y and t | error 5.5e-4 |
| `channel_ramp.toml` | inflow ramped by a table in t | error 1.6e-3 |

How to use them: [recipe 31](../recipes/properties.md#expressions),
[recipe 32](../recipes/properties.md#materials).

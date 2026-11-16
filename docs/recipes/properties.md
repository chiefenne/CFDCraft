# Expressions and materials

Part of the [Recipes](../recipes.md): boundary values, body forces and viscosities as expressions, and fluid properties from the material database. The keys are on the [key page](../cases/flow-keys.md#expressions-materials-and-regions-d53-d102-d112).

| # | Recipe |
|---|---|
| 31 | [Write boundary values, body forces and the viscosity as expressions](#expressions) |
| 32 | [Take the fluid properties from the material database](#materials) |

---

<a id="expressions"></a>
## 31. Write boundary values, body forces and the viscosity as expressions

**A value as a string.** Any velocity, pressure, mass flux, mass flow,
body force or viscosity may be a string instead of a number (D53). Name
the numbers once in `[constants]`:

```toml
[constants]
u_max = 1.5
h = 0.41
omega = "2 * pi * 4.0"        # constants may use other constants

[boundary.inlet]
type = "velocity-inlet"
velocity = ["4 * u_max * y * (h - y) / h^2", 0.0, 0.0]

[boundary.outlet]
type = "static-pressure-outlet"
pressure = "100 * sin(omega * t)"   # transient cases only

[flow]
body_force = [0.0, "-1.2 * 9.81", 0.0]

[fluid]
viscosity = "1e-3 * (1 + 4 * y)"    # x, y, z; not t
```

**Tables.** A `[tables.<name>]` is a function of one argument, piecewise
linear and constant beyond its ends; inline or from a two-column file:

```toml
[tables.ramp]
x = [0.0, 0.5]   # t (s)
y = [0.0, 1.0]

[tables.lift]
file = "lift_curve.csv"   # relative to the case file

[boundary.inlet]
type = "velocity-inlet"
velocity = ["6 * y * (1 - y) * ramp(t)", 0.0, 0.0]
```

**Check before you run.** `foehn flow --check` parses and compiles every string
and reports the mistakes with their place and the nearest name:

```
target/release/foehn flow --check examples/cases/flow/channel_ramp.toml
target/release/foehn flow examples/cases/flow/channel_ramp.toml
target/release/foehn flow examples/cases/flow/taylor_green_expression.toml
target/release/foehn flow examples/cases/flow/couette_viscosity.toml
```

| Code | Meaning |
|---|---|
| EXPR001 | syntax error (with a caret under the place) |
| EXPR002 | unknown or reserved name, cycle between constants (with "did you mean") |
| EXPR003 | a variable the key does not take: `t` in a steady case or in the viscosity, x, y, z in a mass flow |
| EXPR004 | a table that is missing or not increasing |
| EXPR005 | a value that is not finite where it is evaluated |

**Walls with a mass flux.** `mass_flux` on a wall (kg/(m² s), into the
domain) blows or sucks through it; the wall's tangential velocity stays
as given. The flow needs a pressure boundary that takes up the net mass
(else BC026):

```toml
[boundary.porous_wall]
type = "wall"
mass_flux = "0.2 * step(x, 0.5)"   # blowing downstream of x = 0.5 m
```

The reserved names: the variables `x y z t`, the constants `pi e`, the
functions, `if clamp and or not`, and the names of solution fields that
later stages bring into expressions (`u v w p T rho mu k omega epsilon
nut Y`).

<a id="materials"></a>
## 32. Take the fluid properties from the material database

**List and evaluate.** The database has dry air (ideal gas, Sutherland),
water (IAPWS tables), engine oil, ethylene glycol and glycerol, each with
the source of its data and its temperature range:

```
target/release/foehn materials
target/release/foehn materials water --temperature 313.15
target/release/foehn materials air --temperature 350 --pressure 2e5
target/release/foehn --json materials        # for scripts
```

**Use one.** Name it in `[fluid]`; until the energy equation (L1) the
properties are evaluated once, at `temperature` and `operating_pressure`:

```toml
[fluid]
material = "water"
temperature = 313.15          # K, default 293.15
operating_pressure = 101325   # Pa, the default
```

**Change a property.** A `[materials.<name>]` table with `base` takes a
database entry and overrides what it names; a new material needs no base:

```toml
[materials.warm-oil]
base = "engine-oil"
viscosity = { form = "power-law", value_ref = 0.25, t_ref = 320.0, exponent = -6.0, source = "own measurement" }

[materials.test-fluid]
kind = "liquid"
density = 1050.0
viscosity = { form = "table", temperature = [280.0, 300.0, 320.0], values = [4e-3, 2.5e-3, 1.7e-3], interpolation = "monotone-cubic" }
```

**Regions.** `[regions.<name>]` gives cell zones their material (D112):
fluid regions only until the energy equation (REG001); every zone in one
region (REG002). A region with another viscosity gives a viscosity per
face:

```toml
[regions.core]
kind = "fluid"
zones = ["fluid"]
material = "water"
```

```
target/release/foehn flow examples/cases/flow/poiseuille_water.toml
grep -A12 "resolved property program" output/poiseuille_water/case.resolved.toml
```

`case.resolved.toml` lists the resolved properties with their form,
temperature and source; `properties_hash` and `materials_database` in the
provenance of every output file identify the property program and the
database it came from.

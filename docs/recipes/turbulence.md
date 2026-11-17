# Turbulence models

Part of the [Recipes](../recipes.md): run a case with a RANS turbulence model, choose the wall treatment, set inlet turbulence and roughness, and check the near-wall mesh. Every key is on the [key page](../cases/flow-keys.md#turbulence-d129d139).

| # | Recipe |
|---|---|
| 68 | [Run a case with the k-ω SST model](#sst) |

---

<a id="sst"></a>
## 68. Run a case with the k-ω SST model

Switch the model on and pick the wall treatment by the mesh: `low-re`
when the first cell lies at y⁺ ≲ 1, `automatic` (the default) for any
other mesh, best with the first cell at y⁺ ≳ 25 (in the buffer layer,
3 ≲ y⁺ ≲ 20, it overpredicts the wall shear by up to 16 %, and TURB001
says so):

```toml
[turbulence]
model = "k-omega-sst"
wall_treatment = "low-re"        # or "automatic"
# production_limiter = "menter"  # "kato-launder", "durbin", "none"
# ambient = { k = 2.25e-7, omega = 125.0 }   # sustaining terms, external flows
```

Give the inlets their turbulence: an intensity with a length scale or a
viscosity ratio, or k and ω directly. Without them the inlet takes
I = 5 % and ν_t/ν = 10 (TURB009 names it). Rough walls take the
sand-grain roughness:

```toml
[boundary.inlet]
type = "velocity-inlet"
velocity = [10.0, 0.0, 0.0]
turbulence_intensity = 0.02
turbulence_length_scale = 0.01   # or: viscosity_ratio = 10.0, or: k = ..., omega = ...
[boundary.wall]
type = "wall"
roughness = 1e-4                 # k_s (m), optional
```

The start values come from the first inlet unless `[initial]` gives
`k` and `omega` (numbers or expressions) or `turbulence_intensity` and
`viscosity_ratio`. Check the case, run it and read the first-cell y⁺ per
wall, the largest ν_t/ν and the floor counts at the end of the log:

```
target/release/foehn flow examples/cases/flow/tmr_plate_sst.toml --check
target/release/foehn flow examples/cases/flow/tmr_plate_sst.toml
grep -A3 "^turbulence:" output/tmr_plate_sst/flow.log
```

The report has the same under `turbulence` in `result.json`; the wall
shear stress and y⁺ per wall face are in `wall_<boundary>.csv`
(`output.wall_csv = true`). For a mesh study of C_f, run the case on four
meshes and read the discretisation uncertainty in `flow.md`:

```
target/release/foehn flow examples/cases/flow/tmr_plate_sst.toml --out output/plate \
    --set 'verification.levels=[0, 1, 2, 3]'
```

A fully developed channel needs one periodic module (recipe 67) and the
`low-re` treatment with y⁺ ≲ 1:

```
target/release/foehn flow examples/cases/flow/channel_sst.toml
grep "periodic xmin" output/channel_sst/flow.log
```

u_τ = √(β h/ρ) from the driving force β and the half height h, so
Re_τ = u_τ h/ν and U_b⁺ = U_b/u_τ. `output.derived = ["turbulence-features"]`
writes the invariant features of D166 (the inputs of data-driven
corrections) with the fields. The wall distance is exact by default;
`turbulence.wall_distance = "pde"` takes the Poisson distance instead
(cheaper on very large meshes, exact on plane walls).

[↑ back to top](#turbulence-models)

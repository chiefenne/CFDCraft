# Turbulence models

Part of the [Recipes](../recipes.md): run a case with a RANS turbulence model, choose the wall treatment, set inlet turbulence and roughness, and check the near-wall mesh. Every key is on the [key page](../cases/flow-keys.md#turbulence-d129d139).

| # | Recipe |
|---|---|
| 68 | [Run a case with the k-ω SST model](#sst) |
| 71 | [Run a case with the k-ε family and wall functions](#k-epsilon) |
| 73 | [Run a case with k-ε-ζ-f and resolved walls](#zeta-f) |
| 74 | [Run a case with k-ε-ζ-f₀ at any first-cell y⁺](#zeta-f0) |

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

---

<a id="k-epsilon"></a>
## 71. Run a case with the k-ε family and wall functions

The standard k-ε (`k-epsilon`) or the realizable k-ε
(`realizable-k-epsilon`) with wall functions: put the first cell in the
log layer, y⁺ ≳ 30. `scalable` (the default) takes any y⁺ without
failing and stays close below 30; `standard` needs 30 ≤ y⁺ ≤ 300;
`compound` is offered, but with these high-Re models it is claimed from
y⁺ 30 only. Walls resolved to y⁺ ≈ 1 need SST `low-re` now, ζ-f or EB
later (TURB001 warns):

```toml
[turbulence]
model = "k-epsilon"              # or "realizable-k-epsilon"
wall_treatment = "scalable"      # "standard", "compound"
# production_limiter = "durbin"  # the default of the standard model; "kato-launder", "menter" (min(P, 10ε)), "none"
# ambient = { k = 2.25e-7, epsilon = 2.53125e-6 }   # sustaining terms, external flows
```

Inlets take ε directly or ω (converted with ε = C_μ k ω), or the
intensity with a length scale or viscosity ratio, as for SST:

```toml
[boundary.inlet]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
k = 2.25e-7
epsilon = 2.53125e-6             # or: omega = 125.0
```

With free-stream values that decay fast (the TMR plate: ν_t/ν = 0.009)
give the sustaining terms the inlet values, or the boundary layer starts
from nothing. Run the plate, the channel and the rough pipe of the
gallery:

```
target/release/foehn flow examples/cases/flow/tmr_plate_keps.toml
target/release/foehn flow examples/cases/flow/channel_keps.toml
target/release/foehn flow examples/cases/flow/pipe_rough_keps.toml
grep -E "wall probe|periodic xmin|first-cell" output/*_keps/flow.log
```

C_f of the plate is in the wall probe line; in the channel
U_b⁺ = U_b/√(F h/ρ) with the driving force F; in the pipe the Darcy
friction factor is f = 2 F D/(ρ U²). A run that comes from SST starts
from its restart file: `[run] restart = 'output/plate/restart.h5'` with
`model = "k-epsilon"` converts ω to ε.

Rough walls: `roughness` is the sand-grain k_s; C_s of Cebeci & Bradshaw
is `roughness_constant` (0.5). For pipes with Colebrook-type commercial
roughness, 0.4 puts the fully rough friction factor within 3 % of
Colebrook's (0.5: 8 to 11 % above):

```toml
[boundary.wall]
type = "wall"
roughness = 1e-3
roughness_constant = 0.4
```

The Moody chart as a parametric study (16 runs, about 8 minutes with four
jobs; needs the Python package or `FOEHN=target/release/foehn`):

```
python3 examples/studies/moody_study.py output/moody --jobs 4 --chart output/moody.svg
python3 examples/studies/moody_study.py output/moody_cs04 --cs 0.4 --jobs 4
```

The y⁺ study of the plate and the figure of gallery 71:

```
for r in "1 17193.2" "2 7983.56" "5 2859.14" "10 1298.87" "30 361.051" "60 156.835" "100 83.3061"; do set -- $r
  for t in scalable compound; do
    target/release/foehn flow examples/cases/flow/tmr_plate_keps.toml --out output/plate_keps/${t}_$1 \
        --set "turbulence.wall_treatment=\"$t\"" --set "mesh.generate.grading.y.ratio=$2"
  done
done
python3 -m tools.plot.keps_chart docs/figures/gallery/keps.svg output/plate_keps output/moody/study.csv
```

[↑ back to top](#turbulence-models)

---

<a id="zeta-f"></a>
## 73. Run a case with k-ε-ζ-f and resolved walls

k-ε-ζ-f (`k-epsilon-zeta-f`) for walls resolved to y⁺ ≈ 1, where the
k-ε family needs wall functions. `compound` (the default) takes any y⁺;
`low-re` needs y⁺ ≲ 2 (TURB001):

```toml
[turbulence]
model = "k-epsilon-zeta-f"
wall_treatment = "compound"      # or "low-re"
# ambient = { k = 2.25e-7, epsilon = 2.53125e-6 }   # sustaining terms, external flows
```

Inlets take k with ε (or ω, converted) as for k-ε, and ζ = v²/k
(default 2/3, isotropic turbulence); `[initial]` takes `zeta` too:

```toml
[boundary.inlet]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
k = 1.0e-4
epsilon = 4.5e-5
zeta = 0.667
```

A weak free stream (ν_t/ν well below 1) leaves the boundary layer
laminar with this model: start from a turbulent field (`[initial] k`,
`epsilon` with ν_t/ν ≈ 10), as the plate of gallery 73 does, or give the
inlet more turbulence.

The f equation runs CG with its own multigrid hierarchy; the transport
equations take the multigrid too, which pays on large meshes:

```toml
[solver.elliptic]
preconditioner = "amg"           # default; "jacobi"
tolerance = 1e-2                 # default steady (transient 1e-8)

[solver.turbulence]
preconditioner = "amg"           # default "jacobi"
```

Run the channel and the plate of the gallery:

```
target/release/foehn flow examples/cases/flow/channel_zeta_f.toml
target/release/foehn flow examples/cases/flow/tmr_plate_zeta_f.toml
grep -E "wall probe|periodic xmin|first-cell|floors" output/*_zeta_f/flow.log
```

U_b⁺ of the channel is 17.55/√F with the driving force F; C_f of the plate
is in the wall probe line. The profile line of the channel
(`line_profile.csv`) has ζ and f next to k, ε and ν_t; the fields add
`zeta`, `f`, `turbulent_time_scale` and v² = ζk. A k-ε restart file
continues with ζ = 2/3 and f = 0.

[↑ back to top](#turbulence-models)

---

<a id="zeta-f0"></a>
## 74. Run a case with k-ε-ζ-f₀ at any first-cell y⁺

k-ε-ζ-f₀ (`k-epsilon-zeta-f0`, Popovac 2023) for meshes whose first-cell
y⁺ is not controlled, from resolved walls to the log layer.
`compound-2023` is its default treatment; `compound` and `low-re` (y⁺ ≲ 2)
are also accepted, and ζ-f takes `compound-2023` too:

```toml
[turbulence]
model = "k-epsilon-zeta-f0"
wall_treatment = "compound-2023"   # default; "compound", "low-re"
# constants = { n_zeta = 2.0, c_epsilon1_zeta = 0.012 }   # the defaults
```

Inlets and `[initial]` take k, ε (or ω) and ζ as for ζ-f (recipe 73).
Unlike ζ-f, a weak free stream is fine: the TMR plate starts from
ν_t/ν = 0.009 and turns turbulent on its own.

A step or another geometry with a cut-out corner can come from the
built-in mesher: `remove` cuts boxes out of the mesh and names the new
faces:

```toml
[mesh.generate]
type = "hex"
dimension = 2
lower = [-85.0, 0.0, 0.0]
upper = [40.0, 9.0, 0.1]
cells = [250, 160, 1]
remove = [ { lower = [-85.0, 0.0, 0.0], upper = [0.0, 1.0, 0.1], name = "step" } ]

[boundary.step]
type = "wall"
```

Run the channel, the plate and the step of the gallery:

```
target/release/foehn flow examples/cases/flow/channel_zeta_f0.toml
target/release/foehn flow examples/cases/flow/tmr_plate_zeta_f0.toml
target/release/foehn flow examples/cases/flow/backstep_zeta_f0.toml
grep -E "wall probe|periodic xmin|first-cell|floors" output/*_zeta_f0/flow.log
```

The reattachment length of the step is where the wall shear of `ymin`
in `wall_ymin.csv` (behind the step, x > 0) changes sign from negative
to positive.

[↑ back to top](#turbulence-models)

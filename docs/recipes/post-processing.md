# Post-processing

Part of the [Recipes](../recipes.md): derived fields with the output or afterwards from a restart file. The keys are on the [key page](../cases/flow-keys.md#derived-fields-d72).

| # | Recipe |
|---|---|
| 30 | [Write vorticity, Q, λ₂, wall shear stress and other derived fields](#derived-fields) |
| 108 | [Write the Reynolds stress, the scales and the wall values of a turbulence model](#turbulence-derived-fields) |

---

<a id="derived-fields"></a>
## 30. Write vorticity, Q, λ₂, wall shear stress and other derived fields

**With the output.** List the fields in `[output] derived`; every
`fields*.vtkhdf` of the run then carries them (D72):

```toml
[output]
derived = ["vorticity", "q-criterion", "lambda2", "wall-shear-stress", "y-plus"]
```

The names: `vorticity`, `q-criterion`, `lambda2`, `enstrophy`,
`strain-rate`, `pressure-coefficient`, `total-pressure`,
`wall-shear-stress`, `y-plus`, and with a turning zone `relative-velocity`,
`frame-velocity`, `turbomachinery`, `relative-total-pressure`, `rothalpy`
([recipe 99](flow.md#rotor-fields)). The order of the list does not matter.
Cell fields go to the `internal` block; the wall shear stress and y⁺ go to
the boundary blocks (`boundary/<boundary>` in ParaView's block selector), zero on
boundaries that are not walls. The `units` attribute of the file lists every
dataset with its unit. The list may change during a run (re-read, D40).

**The pressure coefficient** needs a reference state; the density defaults
to `[fluid] density`, the velocity to that of `[output.coefficients]`:

```toml
[output.reference]
velocity = 1.0   # U_ref, m/s
pressure = 0.0   # p_ref, Pa
```

**After the run, from a restart file.** `foehn post derive` reads the
restart file and the case and writes the fields with the derived fields,
the file the run would have written at that iteration or step:

```
target/release/foehn flow examples/cases/flow/couette_derived.toml
target/release/foehn post derive examples/cases/flow/couette_derived.toml
# output/couette_derived/derived.vtkhdf (from restart.h5)

target/release/foehn post derive CASE.toml --restart output/CASE/restart_000400.h5 \
    --fields vorticity,q-criterion --file vortices.vtkhdf
target/release/foehn --json post derive CASE.toml   # summary as JSON
```

Without `--fields` it writes `[output] derived` of the case, or every
field when that list is empty (the pressure coefficient only with a
reference velocity). `--restart` defaults to `latest`, the newest restart
file of the output directory; the file is `derived.vtkhdf` for
`restart.h5` and `derived_<N>.vtkhdf` for `restart_<N>.h5`. A run keeps
its restart file with `[run] restart_final = true` (or
`--set run.restart_final=true`).

**What they mean.** Q > 0 and λ₂ < 0 mark vortex cores (Hunt, Wray and
Moin 1988; Jeong and Hussain 1995); both are zero in pure shear, and Q is
also positive in boundary layers on convex walls. The wall values are those
of `wall_<boundary>.csv`: τ_w the tangential viscous stress of the fluid on
the wall, y⁺ = d √(τ_w/ρ)/ν with d the wall distance of the cell centroid.
Pictures: [gallery 22](../gallery/post-processing.md#derived-fields).

---

<a id="turbulence-derived-fields"></a>
## 108. Write the Reynolds stress, the scales and the wall values of a turbulence model

With a turbulence model (`[turbulence] model`), four more names go into
`[output] derived` or `foehn post derive --fields` (D72, specification
2.3):

```toml
[output]
derived = ["turbulence-scales", "reynolds-stress", "turbulence-wall"]
# the elliptic-blending k-ε also: "wall-normal"
```

- `turbulence-scales`: ε, ω, the time scale T, the length scale k^{3/2}/ε
  and ν_t/ν for every model. SST solves ω, so ε is β* k ω; the k-ε family
  solves ε, so ω is ε/(C_μ k). T is the one the model uses: 1/(β* ω) for
  SST, k/ε for the k-ε family, the bounded T of ζ-f and EB (also written
  as `turbulent_time_scale` with the fields).
- `reynolds-stress`: the six components of ⅔ k δᵢⱼ − 2 ν_t Sᵢⱼ, the
  stress the momentum equation sees (the eddy-viscosity models carry no
  other); compare `reynolds_stress_xy` with a measured −u'v'.
- `wall-normal`: n = ∇α/|∇α| of the elliptic-blending k-ε, the wall
  direction of its blending (0 where ∇α vanishes); refused with the other
  models.
- `turbulence-wall`: u_τ of the wall treatment and y* = d C_μ^¼ k_P^½/ν on
  the boundary blocks (zero off walls). y* is the wall coordinate of the
  k-ε wall functions and of the compound treatments; y⁺ is `y-plus`. The
  heat-transfer coefficient and T⁺ of each scalar are the columns
  `<name>_h` and `<name>_t_plus` of `wall_<boundary>.csv`.

A laminar case refuses all four with CASE005. Pictures:
[gallery 108](../gallery/post-processing.md#turbulence-derived-fields).

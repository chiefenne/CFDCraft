# Post-processing

Part of the [Recipes](../recipes.md): derived fields with the output or afterwards from a restart file. The keys are on the [key page](../cases/flow-keys.md#derived-fields-d72).

| # | Recipe |
|---|---|
| 30 | [Write vorticity, Q, λ₂, wall shear stress and other derived fields](#derived-fields) |

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
`wall-shear-stress`, `y-plus`. The order of the list does not matter.
Cell fields go to the `internal` block; the wall shear stress and y⁺ go to
the patch blocks (`boundary/<patch>` in ParaView's block selector), zero on
patches that are not walls. The `units` attribute of the file lists every
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
of `wall_<patch>.csv`: τ_w the tangential viscous stress of the fluid on
the wall, y⁺ = d √(τ_w/ρ)/ν with d the wall distance of the cell centroid.
Pictures: [gallery 22](../gallery/post-processing.md#derived-fields).

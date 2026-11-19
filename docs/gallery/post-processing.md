# Post-processing: derived fields

Part of the [Example gallery](../gallery.md): fields derived from the solution, written with the fields or afterwards from a restart file (D72).

| # | Example | Stage | Feature |
|---|---|---|---|
| 22 | [Derived fields: vorticity, Q, λ₂, wall shear stress and y⁺](#derived-fields) | 5 | `[output] derived`, `foehn post derive` |
| 108 | [Derived fields of the turbulence models: Reynolds stress, scales, wall values](#turbulence-derived-fields) | 6 | `turbulence-scales`, `reynolds-stress`, `wall-normal`, `turbulence-wall` |

---

<a id="derived-fields"></a>
## 22. Derived fields: vorticity, Q, λ₂, wall shear stress and y⁺

<img src="../figures/gallery/cylinder_vorticity.png" alt="The vortex street behind the cylinder at Re = 100, t = 8 s: vorticity ω_z" width="800" height="509">

<img src="../figures/gallery/cylinder_q_criterion.png" alt="The same instant: the Q criterion, Q &gt; 0 in the vortex cores" width="800" height="509">

**Vortices behind a cylinder.** The short channel of
[gallery 21](boundary-conditions.md#boundary-conditions) at Re = 100,
t = 8 s, 4 928 cells. The run wrote only its restart file;
`foehn post derive` read it afterwards and wrote the fields with the
derived fields, the same numbers the run writes with `[output] derived`
(the two files are equal byte for byte). The vorticity shows the shear
layers leaving the cylinder and rolling up into the street of alternating
vortices; the Q criterion (the rotation rate beats the strain rate) marks
their cores. Q is also positive in the boundary layer on the curved
surface, where it is (u/r) ∂u/∂r: a property of the criterion, not a
vortex. On the cylinder wall, the boundary block of the file carries the wall
shear stress (largest 0.46 Pa) and y⁺ (at most 0.87).

```
target/release/foehn flow examples/cases/flow/cylinder_re100_short.toml --set run.restart_final=true
target/release/foehn post derive examples/cases/flow/cylinder_re100_short.toml \
    --set output.reference.velocity=1.0
# output/cylinder_re100_short/derived.vtkhdf: open it in ParaView
python3 -m tools.mesh.render_mesh output/cylinder_re100_short/derived.vtkhdf vorticity.png \
    --field vorticity_z --top --no-edges --diverging --range -100 100 --size 1100x700 --zoom 1.3
python3 -m tools.mesh.render_mesh output/cylinder_re100_short/derived.vtkhdf q.png \
    --field q_criterion --top --no-edges --range 0 1000 --size 1100x700 --zoom 1.3
```

**Exact on every mesh.** In plane Couette flow the velocity gradient is
constant, and the least-squares gradient reproduces it on any mesh, so
every derived field is exact: ω_z = −1/s, Q = λ₂ = 0 (shear is as much
strain as rotation), shear rate 1/s, wall shear stress μU/H = 0.1 Pa and
y⁺ = d √(τ_w/ρ)/ν from the wall distance d of each cell, to 1e-8 on
hexahedra, tetrahedra, polyhedra and the mixed mesh. On the Taylor–Green
vortex the vorticity, Q and λ₂ converge to the values of the exact
velocity gradient at second order (observed 2.3–2.7 on 16² to 64² cells).

```
target/release/foehn flow examples/cases/flow/couette_derived.toml
target/release/foehn post derive examples/cases/flow/couette_derived.toml
cmp output/couette_derived/fields.vtkhdf output/couette_derived/derived.vtkhdf && echo same
```

| Name in `derived` | Datasets | Unit |
|---|---|---|
| `vorticity` | `vorticity_x`, `vorticity_y`, `vorticity_z`, `vorticity_magnitude` | 1/s |
| `q-criterion` | `q_criterion`: ½ (‖Ω‖² − ‖S‖²) | 1/s² |
| `lambda2` | `lambda2`: the middle eigenvalue of S² + Ω² | 1/s² |
| `enstrophy` | `enstrophy`: ½ ω·ω | 1/s² |
| `strain-rate` | `strain_rate`: √(2 S:S) | 1/s |
| `pressure-coefficient` | `pressure_coefficient`: (p − p_ref)/(½ ρ_ref U_ref²) | - |
| `total-pressure` | `total_pressure`: p + ½ ρ u·u | Pa |
| `wall-shear-stress` | on the boundary blocks: `wall_shear_stress_x`, `_y`, `_z`, `wall_shear_stress_magnitude` | Pa |
| `y-plus` | on the boundary blocks: `y_plus` | - |

How to use them: [recipe 30](../recipes/post-processing.md#derived-fields).

---

<a id="turbulence-derived-fields"></a>
## 108. Derived fields of the turbulence models: Reynolds stress, scales, wall values

<img src="../figures/gallery/turbulence_derived.png" alt="The backward-facing step with k-ε-ζ-f₀ on 8 200 cells, x/H from −3 to 14: the Reynolds shear stress u'v'/U_ref² of the Boussinesq estimate (down to −0.011 in the shear layer at x/H 4 to 7), the eddy-viscosity ratio ν_t/ν (up to 490 near the wall downstream of reattachment) and the time scale T of the model" width="800" height="576">

**What the model carries.** The step of
[gallery 74](turbulence.md#zeta-f0) (Driver & Seegmiller, Re_H = 36 000,
k-ε-ζ-f₀, 8 200 cells, 1 500 iterations, 66 s) kept its restart
file; `foehn post derive` added the quantities of the turbulence model.
The Reynolds shear stress −2 ν_t S_xy grows along the shear layer that
leaves the step corner and is largest, −0.011 U_ref², at x/H 4 to 7, just
upstream of the reattachment at 6.4 H; the eddy viscosity is largest near
the lower wall downstream of the reattachment, where the new boundary
layer grows under the turbulent remains of the shear layer (ν_t/ν up to
490); the time scale of the model is short in the shear layer and long in
the slowly decaying free stream and the recirculation bubble.

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0.toml --set run.restart_final=true
target/release/foehn post derive examples/cases/flow/backstep_zeta_f0.toml \
    --fields turbulence-scales,reynolds-stress,turbulence-wall
python3 -m tools.plot.turbulence_fields turbulence_derived.png \
    output/backstep_zeta_f0/derived.vtkhdf
```

| Name in `derived` | Datasets | Unit |
|---|---|---|
| `turbulence-scales` | `dissipation_rate` ε (β* k ω for SST), `specific_dissipation_rate` ω (ε/(C_μ k) for the ε models), `time_scale` T of the model (1/(β* ω), k/ε, or the bounded T of ζ-f and EB), `length_scale` k^{3/2}/ε, `viscosity_ratio` ν_t/ν | m²/s³, 1/s, s, m, - |
| `reynolds-stress` | `reynolds_stress_xx`, `_yy`, `_zz`, `_xy`, `_yz`, `_xz`: ⅔ k δᵢⱼ − 2 ν_t Sᵢⱼ | m²/s² |
| `wall-normal` | `wall_normal_x`, `_y`, `_z`: ∇α/\|∇α\| of the elliptic-blending k-ε | - |
| `turbulence-wall` | on the boundary blocks: `friction_velocity` u_τ of the wall treatment, `y_star` = d C_μ^¼ k_P^½/ν | m/s, - |

**Checked on the channel** of every model family (SST, k-ε, EB;
`tests/turbulence_derived.rs`): ε = C_μ k ω, T and L from their
definitions, the trace of the Reynolds stress 2k, u'v' negative below the
centre line and positive above, u_τ that of the wall treatment, and the
EB wall normal a unit vector pointing away from the nearer wall.

How to use them: [recipe 108](../recipes/post-processing.md#turbulence-derived-fields).

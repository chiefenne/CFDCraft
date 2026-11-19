# Python scripts (`examples/python/`, `examples/studies/`)

Part of the [Case catalogue](../cases.md): scripts that build, run and study cases with the `foehncraft` Python package (D61, D62). They start from the case files of `examples/cases/` and write the case files they run, so every result can be repeated with `foehn flow`.

| Script | Base case | Shows | Where |
|---|---|---|---|
| `python/poiseuille_inflow.py` | `flow/poiseuille_hex.toml`, coarse mesh, 512 cells | the inflow as a Python expression, a monitor of the outer iterations and a stop criterion of its own (stopped after 145 iterations) | [gallery 24](../gallery/python.md#python-study), [recipe 33](../recipes/python.md#python) |
| `studies/couette_viscosity_study.py` | `flow/couette_viscosity.toml`, 2 048 cells | a study over six viscosity gradients: one case file and output per variant, `study.csv`, resume; the wall shear stress against μ0 U a / ln(1 + a) | [gallery 24](../gallery/python.md#python-study), [recipe 33](../recipes/python.md#python) |
| `studies/moody_study.py` | `flow/pipe_rough_keps.toml`, 12 to 160 cells | the Moody chart: 16 rough-pipe runs over Re 3·10⁴ to 10⁷ and k_s/D 0 to 10⁻² with the k-ε model, f from the driving force against Colebrook, an SVG chart; `--cs` sets C_s | [gallery 71](../gallery/turbulence.md#k-epsilon), [recipe 71](../recipes/turbulence.md#k-epsilon) |
| `studies/mixer.py` | `studies/mixer/mrf.toml` and `sliding.toml`, 16 326 cells | the stirred vessel of gallery 98: builds the Gmsh meshes, runs the frame and the sliding mesh (a finished run is kept) and draws the five figures with `tools/plot/rotor_fields.py` (speed, relative speed in the rotor, vorticity, the start-up, the imprinted interface, and the relative flow from the fields of D183: \|w\|, β, rothalpy) | [gallery 98](../gallery/flow.md#stirred-vessel), [gallery 99](../gallery/flow.md#rotor-fields), [recipe 98](../recipes/flow.md#stirred-vessel) |

The tests of the package (`python/tests/`) run the example cases in
process and through the command line: a case built in Python and the
same case written by hand give the same fields bit for bit.

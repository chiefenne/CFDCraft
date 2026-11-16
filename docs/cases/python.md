# Python scripts (`examples/python/`, `examples/studies/`)

Part of the [Case catalogue](../cases.md): scripts that build, run and study cases with the `foehncraft` Python package (D61, D62). They start from the case files of `examples/cases/` and write the case files they run, so every result can be repeated with `foehn flow`.

| Script | Base case | Shows | Where |
|---|---|---|---|
| `python/poiseuille_inflow.py` | `flow/poiseuille_hex.toml`, coarse mesh, 512 cells | the inflow as a Python expression, a monitor of the outer iterations and a stop criterion of its own (stopped after 145 iterations) | [gallery 24](../gallery/python.md#python-study), [recipe 33](../recipes/python.md#python) |
| `studies/couette_viscosity_study.py` | `flow/couette_viscosity.toml`, 2 048 cells | a study over six viscosity gradients: one case file and output per variant, `study.csv`, resume; the wall shear stress against μ0 U a / ln(1 + a) | [gallery 24](../gallery/python.md#python-study), [recipe 33](../recipes/python.md#python) |

The tests of the package (`python/tests/`) run the example cases in
process and through the command line: a case built in Python and the
same case written by hand give the same fields bit for bit.

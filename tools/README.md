# Tools

The scripts around the solver, all in Python (D61). Run them from the
repository root as modules, `python3 -m tools.<group>.<name>`; each one
prints its usage with `--help` or in its docstring.

| Group | Scripts | What for |
|---|---|---|
| `tools/check/` | `ci_cases`, `compare_runs`, `backend`, `c_free`, `kernels`, `vtkhdf` | the case catalogue as CI runs it, bitwise comparison of runs (D31), the backend lint (D60), the C-free core (D87), the kernel rules (R10), VTKHDF files in VTK's reader |
| `tools/bench/` | `build_hypre`, `cloud`, `g1_cloud`, `stage2_cloud`, `stage3_cloud`, `stage5_cloud`, `stage6_cloud` | hypre from source (D20); the drivers for a cloud GPU (Colab, Kaggle, any Linux machine with an NVIDIA GPU), with their shared steps in `cloud` |
| `tools/mesh/` | `gmsh_test_meshes`, `render_mesh` | the Gmsh test meshes of the mesh tests; headless pictures of meshes and fields (L21) |
| `tools/plot/` | `*_chart`, `backstep_paraview` | the charts of the gallery and the verification pages; `backstep_paraview` runs in ParaView's Python (`pvpython tools/plot/backstep_paraview.py`) |
| `tools/release/` | `verification_manual` | the generated verification manual (D68) |
| `tools/notebooks/` | `*_colab.ipynb` | Google Colab notebooks that run the cloud drivers on a free T4 |

`tools/requirements-ci.txt` lists the Python packages of the CI checks
(VTK, h5py, NumPy). The PTX kernels are built by
`python3 gpu/ptx-kernels/build.py`. The package `foehncraft` itself lives in
[`python/`](../python/README.md).

# Example gallery

What Foehncraft computes, one picture per case, with the commands that
reproduce it. Command sequences for everyday tasks are in
[Recipes](recipes.md); all case files are listed in the
[Case catalogue](cases.md).

All commands run from the repository root. Case files are in
`examples/cases/`.

The gallery is split by topic; each page has its own table of contents. The numbers are global and stay the same, so "gallery 12" always means the same entry.

| Page | What it covers | Examples |
|---|---|---|
| [Kernel benchmarks and GPU runs](gallery/benchmarks.md) | gate G1 on CPUs and GPUs, on box, polyhedral and tetrahedral meshes, and the flow loop on a GPU | 1, 2, 7, 16 |
| [Meshes](gallery/meshes.md) | the built-in test mesher in 2D and 3D, partitioning, the Gmsh reader and flow on poor Gmsh meshes | 3, 4, 5, 6, 26 |
| [Scalar transport and the pressure solver](gallery/scalar-pressure.md) | convection–diffusion order of accuracy, convection schemes and limiters, the pressure multigrid | 8, 9, 10, 11 |
| [Incompressible flow](gallery/flow.md) | steady and transient laminar verification and benchmarks: cavity, Poiseuille, Taylor–Green, cylinder, time-step control | 12, 13, 14, 15, 20 |
| [Boundary conditions](gallery/boundary-conditions.md) | slip walls, axisymmetric wedges, inlet profiles from tables, backflow at pressure outlets, convective outlets | 21 |
| [Expressions and materials](gallery/properties.md) | boundary values, body forces and viscosities as expressions, the material database, regions, walls with a mass flux | 23 |
| [Python](gallery/python.md) | cases built, run and studied from Python: expressions with Python operators, callbacks, study tables | 24 |
| [Verification](gallery/verification.md) | the discretisation uncertainty of mesh sweeps (grid convergence index) and the generated verification manual | 25 |
| [Post-processing](gallery/post-processing.md) | derived fields: vorticity, Q, λ₂, wall shear stress, y⁺, from the run or a restart file | 22 |
| [Run control, health and adjoint readiness](gallery/runs.md) | adjoint-ready kernels, runs that a program drives, balances and the divergence report | 17, 18, 19 |

## All examples by number

| # | Example | Page |
|---|---|---|
| <a id="g1-cpu"></a>1 | [Kernel benchmark G1 on the CPU](gallery/benchmarks.md#g1-cpu) | [Kernel benchmarks and GPU runs](gallery/benchmarks.md) |
| <a id="g1-gpu"></a>2 | [Kernel benchmark G1 on a Tesla T4](gallery/benchmarks.md#g1-gpu) | [Kernel benchmarks and GPU runs](gallery/benchmarks.md) |
| <a id="mesh-types-3d"></a>3 | [The four test-mesh types in 3D](gallery/meshes.md#mesh-types-3d) | [Meshes](gallery/meshes.md) |
| <a id="poly-partition"></a>4 | [A polyhedral mesh in four parts](gallery/meshes.md#poly-partition) | [Meshes](gallery/meshes.md) |
| <a id="mesh-types-2d"></a>5 | [The four test-mesh types in 2D](gallery/meshes.md#mesh-types-2d) | [Meshes](gallery/meshes.md) |
| <a id="tjunction"></a>6 | [T-junction from Gmsh](gallery/meshes.md#tjunction) | [Meshes](gallery/meshes.md) |
| <a id="g1-stage1"></a>7 | [G1 on polyhedral and tetrahedral meshes](gallery/benchmarks.md#g1-stage1) | [Kernel benchmarks and GPU runs](gallery/benchmarks.md) |
| <a id="scalar-mms"></a>8 | [Convection–diffusion, manufactured solution](gallery/scalar-pressure.md#scalar-mms) | [Scalar transport and the pressure solver](gallery/scalar-pressure.md) |
| <a id="oblique-step"></a>9 | [Oblique step: five convection schemes](gallery/scalar-pressure.md#oblique-step) | [Scalar transport and the pressure solver](gallery/scalar-pressure.md) |
| <a id="poisson-iterations"></a>10 | [Pressure solver: iterations over the mesh size](gallery/scalar-pressure.md#poisson-iterations) | [Scalar transport and the pressure solver](gallery/scalar-pressure.md) |
| <a id="poisson-amg-mms"></a>11 | [Poisson manufactured solution with the multigrid](gallery/scalar-pressure.md#poisson-amg-mms) | [Scalar transport and the pressure solver](gallery/scalar-pressure.md) |
| <a id="cavity"></a>12 | [Lid-driven cavity against Ghia et al.](gallery/flow.md#cavity) | [Incompressible flow](gallery/flow.md) |
| <a id="poiseuille"></a>13 | [Poiseuille flow: order of accuracy of the flow solver](gallery/flow.md#poiseuille) | [Incompressible flow](gallery/flow.md) |
| <a id="taylor-green"></a>14 | [Taylor–Green vortex: second order in time and space](gallery/flow.md#taylor-green) | [Incompressible flow](gallery/flow.md) |
| <a id="cylinder"></a>15 | [Cylinder in a channel: the Schäfer & Turek benchmark](gallery/flow.md#cylinder) | [Incompressible flow](gallery/flow.md) |
| <a id="gpu-flow"></a>16 | [The flow loop on a GPU](gallery/benchmarks.md#gpu-flow) | [Kernel benchmarks and GPU runs](gallery/benchmarks.md) |
| <a id="adjoint-ready"></a>17 | [Adjoint-ready kernels: derivatives and checkpointing](gallery/runs.md#adjoint-ready) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="run-contract"></a>18 | [A run that a program can drive](gallery/runs.md#run-contract) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="health"></a>19 | [Balances and the divergence report](gallery/runs.md#health) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="time-control"></a>20 | [Time step from the error, output at exact times, potential-flow start](gallery/flow.md#time-control) | [Incompressible flow](gallery/flow.md) |
| <a id="boundary-conditions"></a>21 | [Boundary conditions: slip walls, wedges, tables, backflow, convective outlets](gallery/boundary-conditions.md#boundary-conditions) | [Boundary conditions](gallery/boundary-conditions.md) |
| <a id="derived-fields"></a>22 | [Derived fields: vorticity, Q, λ₂, wall shear stress and y⁺](gallery/post-processing.md#derived-fields) | [Post-processing](gallery/post-processing.md) |
| <a id="expressions-materials"></a>23 | [Expressions and materials: variable viscosity, blowing and suction, water from the database](gallery/properties.md#expressions-materials) | [Expressions and materials](gallery/properties.md) |
| <a id="python-study"></a>24 | [A parametric study from Python: Couette flow with a viscosity gradient](gallery/python.md#python-study) | [Python](gallery/python.md) |
| <a id="gci"></a>25 | [Grid convergence index: Couette flow and the cylinder at Re = 20](gallery/verification.md#gci) | [Verification](gallery/verification.md) |
| <a id="skewed-tets"></a>26 | [Flow on raw Gmsh tetrahedra: the cell-quality remediation](gallery/meshes.md#skewed-tets) | [Meshes](gallery/meshes.md) |

## Adding a feature

Every feature ships with a gallery example, case files and a recipe (D18,
as in PyAeroFlow); see [CONTRIBUTING.md](../CONTRIBUTING.md). Put the example on the topic page that fits, or start a new page when a new physics or topic begins (energy, turbulence, reacting flows, sprays, ...). Give it the next free number, an `<a id>` anchor, a row in that page's table of contents and a row in both tables above.

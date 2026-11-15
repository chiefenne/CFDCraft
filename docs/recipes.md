# Recipes

Short command sequences for everyday tasks. Each recipe names the case-file
keys it uses; the keys of the benchmark cases are listed in the
[Case catalogue](cases.md#case-file-keys). Pictures of complete cases are in
the [Example gallery](gallery.md).

1. [Build, test and run the G1 benchmark on your CPU](#g1-cpu)
2. [Run G1 on Windows: CPU verdict and the CUDA reference](#g1-windows)
3. [Run gate G1 on an NVIDIA GPU](#g1-gpu) (also: free Colab T4)
4. [Run the MPI test with several ranks on one machine](#mpi)
5. [Check a change the way CI does](#ci-local)
6. [Make a release (tag and GitHub release)](#release)
7. [Generate a test mesh and check its quality](#mesh-generate)
8. [Read a Gmsh mesh and convert it into polyhedra](#mesh-gmsh)
9. [Partition a mesh and test the halo exchange](#mesh-partition)
10. [Rerun G1 on polyhedral and tetrahedral meshes](#g1-stage1)
11. [Update the golden values of the mesher](#golden)
12. [Check the order of accuracy of the scalar operators](#scalar-mms)
13. [Compare convection schemes and limiters](#scalar-schemes)
14. [Run the scalar operators on an NVIDIA GPU](#scalar-gpu) (also: free Colab T4)
15. [Study the pressure solver on a Poisson problem](#poisson-study)
16. [Run the pressure solver on an NVIDIA GPU against hypre](#poisson-gpu) (also: free Colab T4)
17. [Build hypre and compare with it](#hypre)
18. [Verify the flow solver on Couette and Poiseuille flow](#flow-verification)
19. [Run the lid-driven cavity and compare with Ghia et al.](#flow-cavity)
20. [Set up the boundary conditions of a flow case](#flow-boundaries)

---

<a id="g1-cpu"></a>
## 1. Build, test and run the G1 benchmark on your CPU

Needs Rust 1.89 or newer (`rustup update stable`), a C compiler and CMake:
the HDF5 library behind the VTKHDF writer (D41) and METIS (D32) are built
from source on the first build, which takes a few minutes (macOS:
`brew install cmake`; Ubuntu: `sudo apt install cmake`; Windows: Visual
Studio Build Tools and `winget install Kitware.CMake`). Without them,
`cargo build --release -p cfdcraft-cli --no-default-features` builds
everything except the VTKHDF output and METIS (RCB partitioning remains).

```
git clone https://github.com/chiefenne/CFDCraft.git && cd CFDCraft
cargo test --workspace
cargo build --release -p cfdcraft-cli
target/release/cfdcraft info
target/release/cfdcraft bench g1 examples/cases/bench/g1_ci_small.toml   # seconds
target/release/cfdcraft bench g1 examples/cases/bench/g1_cpu_10m.toml    # about a minute, 4 GB
```

Options override the case file: `--cells 2000000`, `--ordering shuffled`,
`--threads 8`, `--repeat 20`, `--backends cpu-threads,cpu-native`,
`--out DIR` (keys `mesh.cells`, `mesh.ordering`, `run.threads`,
`run.repeat`, `run.backends`, `output.directory`). The report is
`DIR/g1.md`, `g1.json` and `g1.svg`.

[↑ back to top](#recipes)

---

<a id="g1-windows"></a>
## 2. Run G1 on Windows: CPU verdict and the CUDA reference

For a Windows PC with any NVIDIA GPU (the maintainer's machine: i7-4790,
GTX 960). Install Rust with `rustup-init.exe` from https://rustup.rs (it
offers the Visual Studio C++ build tools if they are missing) and Git. In
PowerShell:

```
git clone https://github.com/chiefenne/CFDCraft.git
cd CFDCraft
cargo test --workspace --release
cargo build --release -p cfdcraft-cli
target\release\cfdcraft info
target\release\cfdcraft bench g1 examples\cases\bench\g1_cpu_10m.toml
```

Close other programs first; run the last command twice. It prints
`G1 cpu-threads against cpu-native: PASS` or `FAIL` and writes
`output\g1_cpu_10m\g1.md`, `g1.json` and `g1.svg`.

The CUDA reference (needs an up-to-date NVIDIA driver and the CUDA 12
toolkit from https://developer.nvidia.com/cuda-downloads, for its NVRTC
library; on the GTX 960 use CUDA 12, because CUDA 13 dropped Maxwell):

```
cargo build --release -p cfdcraft-cli --features cuda
target\release\cfdcraft bench backends        # lists cuda-reference
target\release\cfdcraft bench g1 examples\cases\bench\g1_gpu_reference_4m.toml
```

The check column must say `bitwise` for every kernel. Send the two `g1.md`
files; they go into [G1 results](benchmarks/g1.md#maintainers-machine).
Keys: `run.backends`, `run.threads`, `mesh.cells`.

[↑ back to top](#recipes)

---

<a id="g1-gpu"></a>
## 3. Run gate G1 on an NVIDIA GPU

Needs the NVIDIA driver, a CUDA 12 toolkit (only its NVRTC library is used,
for the hand-written reference) and, for the single-source route, a GPU of
compute capability 7.0 (Volta) or newer. `nvidia-smi` shows the GPU.

```
gpu/ptx-kernels/build.sh                       # installs nightly + the nvptx target on first use
cargo build --release -p cfdcraft-cli --features cuda
target/release/cfdcraft bench backends         # must list cuda-reference and cuda-rust
target/release/cfdcraft bench g1 examples/cases/bench/g1_gpu_10m.toml
```

The last lines print `G1 cuda-rust against cuda-reference: PASS` or `FAIL`.
Send `output/g1_gpu_10m/g1.md` and `g1.json`; they go into
[G1 results](benchmarks/g1.md#results). For the exact architecture:
`gpu/ptx-kernels/build.sh sm_86` (your GPU's compute capability). Keys:
`run.backends`, `run.device`, `run.rust_ptx`, `run.tolerance`.

On Windows, run the build script from Git Bash or WSL, or the two `cargo`
commands in it by hand.

**Without a GPU of your own: Google Colab (free T4).** Open
[`tools/g1_colab.ipynb` in Colab](https://colab.research.google.com/github/chiefenne/CFDCraft/blob/main/tools/g1_colab.ipynb),
choose *Runtime → Change runtime type → T4 GPU*, then *Runtime → Run all*.
About 10 minutes; the last cell shows both reports. On any other Linux
machine with a GPU (Kaggle with "GPU T4 x2", Lightning.ai):

```
git clone --depth 1 https://github.com/chiefenne/CFDCraft.git && cd CFDCraft
bash tools/g1_cloud.sh            # installs Rust, builds, runs g1_gpu_10m.toml twice
```

[↑ back to top](#recipes)

---

<a id="mpi"></a>
## 4. Run the MPI test with several ranks on one machine

Needs an MPI library and libclang (Ubuntu:
`sudo apt install libopenmpi-dev openmpi-bin libclang-dev`; macOS:
`brew install open-mpi`).

```
cargo test -p cfdcraft-parallel --features mpi
cargo build --release -p cfdcraft-parallel --features mpi --example mpi_smoke
mpirun --oversubscribe -n 4 target/release/examples/mpi_smoke
```

Prints `MPI smoke test passed` after a ring exchange (the step of a halo
exchange) and a reproducible sum that gives the same bits on every rank.

[↑ back to top](#recipes)

---

<a id="ci-local"></a>
## 5. Check a change the way CI does

```
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo clippy --workspace --all-targets --features cfdcraft-bench/cuda,cfdcraft-fv/cuda -- -D warnings
cargo test --workspace --release
cargo test -p cfdcraft-mesh --no-default-features      # without HDF5 and METIS
cargo test -p cfdcraft-fv --no-default-features
gpu/ptx-kernels/build.sh && tools/check_kernels.sh     # PTX, and R10: no atomics, no f32, no FMA
for f in examples/cases/mesh/*_box.toml examples/cases/mesh/*_2d.toml; do
  target/release/cfdcraft mesh "$f" || break
done
for f in examples/cases/scalar/*.toml; do
  case "$f" in *_gpu.toml) continue;; esac
  target/release/cfdcraft scalar "$f" || break
done
pip install vtk numpy && python3 tools/check_vtkhdf.py output/*/mesh.vtkhdf output/*/fields.vtkhdf output/*/*/fields.vtkhdf
```

The workflow is `.github/workflows/ci.yml`; it runs on every pull request on
Linux, macOS and Windows.

[↑ back to top](#recipes)

---

<a id="release"></a>
## 6. Make a release (tag and GitHub release)

A release is a Git tag `v<stage>` (D46) plus a GitHub release. Steps:

1. A pull request moves the `## Unreleased` section of `CHANGELOG.md` to
   `## <version> (YYYY-MM-DD) – ...` and adds a row to the summary table.
2. After the merge: GitHub, Actions, **Release**, *Run workflow*, version
   `0.1`, first with *dry run* checked (prints the notes), then without.

The workflow refuses a version without a changelog section, and creates
only the release if the tag already exists. The assistant starts it through
the GitHub API in the same way.

---

<a id="mesh-generate"></a>
## 7. Generate a test mesh and check its quality

The built-in mesher (D34) makes a box of hexahedra, tetrahedra, polyhedra
or mixed cells, in 3D or 2D (one cell thick), graded, refined and
perturbed. Keys `[mesh.generate]` ([Case catalogue](cases.md#mesh-case-file-keys)).

```
cargo build --release -p cfdcraft-cli
target/release/cfdcraft mesh examples/cases/mesh/poly_box.toml
cat output/poly_box/mesh_report.md
```

The report gives the verdict (valid or not), the cell types, the closure
|Σ S| / Σ |S| of the worst cell, the total volume against the volume
enclosed by the boundary, the patch areas, non-orthogonality, skewness,
face warp, aspect and volume ratio with the cell or face where the maximum
sits, and the point-location test. `mesh.vtkhdf` holds the mesh with these
fields per cell; open it in ParaView (5.13 or newer reads VTKHDF 2.4 with
polyhedra) or render a picture:

```
pip install vtk
python3 tools/render_mesh.py output/poly_box/mesh.vtkhdf poly.png --cut --field non_orthogonality
```

A grid-convergence family: the same case with `refine = 0, 1, 2` gives
nested meshes (every interval halved), the grading kept.

[↑ back to top](#recipes)

---

<a id="mesh-gmsh"></a>
## 8. Read a Gmsh mesh and convert it into polyhedra

Needs Gmsh 4.x (`pip install gmsh` or gmsh.info). The reader takes MSH 4.1,
ASCII or binary, first-order elements (tetrahedra, hexahedra, prisms,
pyramids), physical groups and periodic links (D33).

```
gmsh -3 examples/meshes/tjunction.geo -format msh41 -o examples/meshes/tjunction.msh
target/release/cfdcraft mesh examples/meshes/tjunction.msh             # defaults, no case file
target/release/cfdcraft mesh examples/cases/mesh/tjunction_tet.toml    # RCM, 4 parts
target/release/cfdcraft mesh examples/cases/mesh/tjunction_poly.toml   # converted into polyhedra
```

In your own case file: `mesh.file` (relative to the case file),
`mesh.scale = 0.001` for a mesh in millimetres, `mesh.polyhedral = true`
for the conversion of a tetrahedral mesh into polyhedra (one cell per node,
same volume), `mesh.renumber = "rcm"` for a band-limited cell order, and
`[mesh.patches]` to rename physical groups. Give every boundary surface a
physical group in Gmsh; surfaces without one become `surface_<tag>`.

[↑ back to top](#recipes)

---

<a id="mesh-partition"></a>
## 9. Partition a mesh and test the halo exchange

```
target/release/cfdcraft mesh examples/cases/mesh/poly_partition.toml
target/release/cfdcraft mesh examples/cases/mesh/poly_partition.toml --parts 8 --method rcb
```

Keys `partition.parts`, `partition.method` (`metis` or `rcb`),
`partition.halo_layers`. The report lists cells, halo cells, processor
faces and neighbours per rank, the imbalance and the cut faces;
`mesh.vtkhdf` has the field `rank`.

The halo exchange over MPI (needs MPI as in recipe 4): every rank builds
its part with a two-layer halo and fills it from its neighbours.

```
cargo build --release -p cfdcraft-mesh --features mpi --example mpi_halo
mpirun --oversubscribe -n 4 target/release/examples/mpi_halo
```

Prints `MPI halo exchange: PASS` when the halo values, volumes and centroids
equal the global mesh and the volume sum gives the same bits on every rank.

[↑ back to top](#recipes)

---

<a id="g1-stage1"></a>
## 10. Rerun G1 on polyhedral and tetrahedral meshes

D8 asks for G1 to be rerun on polyhedral and shuffled meshes in stage 1.
The cloud runs (gallery 7) are on a shared machine; the verdict needs a
quiet one. On the iMac (about 3 GB of memory each, a few minutes in all):

```
git pull && cargo build --release -p cfdcraft-cli
for c in g1_cpu_poly_2m g1_cpu_poly_2m_shuffled g1_cpu_tet_2m; do
  target/release/cfdcraft bench g1 examples/cases/bench/$c.toml
  target/release/cfdcraft bench g1 examples/cases/bench/$c.toml   # second run
done
```

Send `output/g1_cpu_*/g1.md`. The key `mesh.type` (`poly`, `tet`, `hex`,
`mixed`, `periodic-box`) chooses the mesh, `mesh.cells` its size. On a GPU
of compute capability 7.0 or newer (or Colab, recipe 3):
`cfdcraft bench g1 examples/cases/bench/g1_gpu_poly_4m.toml`.

[↑ back to top](#recipes)

---

<a id="golden"></a>
## 11. Update the golden values of the mesher

`crates/cfdcraft-mesh/tests/golden/mesher.json` holds the cell, face and
node counts, volumes, quality metrics and partition statistics of every
combination of mesh type, dimension, perturbation, refinement and 1, 2 and
4 ranks (R4). A difference fails the test. Only for an intended change:

```
CFDCRAFT_BLESS=1 cargo test -p cfdcraft-mesh --test golden
git diff crates/cfdcraft-mesh/tests/golden/mesher.json
```

and say in the commit message why the values changed.

[↑ back to top](#recipes)

---

<a id="scalar-mms"></a>
## 12. Check the order of accuracy of the scalar operators

The manufactured-solution sweeps of stage 2 ([explanation](verification/stage2.md)):
each case solves on three or four nested refinements, prints the observed
order and fails if it is outside 1.8–2.2.

```
cargo build --release -p cfdcraft-cli
for m in hex tet poly mixed; do
  target/release/cfdcraft scalar examples/cases/scalar/mms_${m}_box.toml
done
target/release/cfdcraft scalar examples/cases/scalar/mms_poly_2d.toml
cat output/mms_poly_box/scalar.md          # errors, orders, time per phase
```

Keys: `verification.levels` (refinement levels of `[mesh.generate]`,
override with `--levels 1,2,3,4`), `verification.expected_order`,
`problem.manufactured`, `problem.boundary.<patch>` with `profile =
"manufactured"`, `mesh.generate.perturb`. The same checks as tests (one
level finer with `--ignored`):

```
cargo test --release -p cfdcraft-fv
cargo test --release -p cfdcraft-fv -- --ignored --nocapture   # about 1.6 M tetrahedra, minutes
cargo run --release -p cfdcraft-fv --example sweep_table         # the table of the explanation page
```

[↑ back to top](#recipes)

---

<a id="scalar-schemes"></a>
## 13. Compare convection schemes and limiters

One case, several solver settings: `[[variants]]` runs the case once per
variant, writes the fields of each into its own directory and the profiles
along `output.profile` into one plot.

```
target/release/cfdcraft scalar examples/cases/scalar/oblique_step_2d.toml
cat output/oblique_step_2d/scalar.md       # φ range per scheme
# profile plot: output/oblique_step_2d/profile.svg
# fields: output/oblique_step_2d/<variant>/fields.vtkhdf (open in ParaView, field "phi" and "limiter")
```

Keys: `solver.scheme` (`upwind`, `linear-upwind`), `solver.limiter`
(`none`, `venkatakrishnan`, `venkatakrishnan-wang`), `solver.limiter_k`
(0: strictly bounded; larger: smoother, small overshoots),
`[[variants]].bounds`, `output.profile`. To try another K, add a variant:

```toml
[[variants]]
name = "venkatakrishnan-k03"
solver = { scheme = "linear-upwind", limiter = "venkatakrishnan", limiter_k = 0.3, max_outer = 400 }
```

[↑ back to top](#recipes)

---

<a id="scalar-gpu"></a>
## 14. Run the scalar operators on an NVIDIA GPU

Without a GPU of your own: open
[`tools/stage2_colab.ipynb`](https://colab.research.google.com/github/chiefenne/CFDCraft/blob/main/tools/stage2_colab.ipynb)
in Google Colab, choose a T4 runtime, set `BRANCH` if the stage-2 branch is
not merged yet, and run all cells (about 15 minutes). On any Linux machine
with an NVIDIA GPU of compute capability 7.0 or newer:

```
bash tools/stage2_cloud.sh
```

By hand:

```
gpu/ptx-kernels/build.sh sm_75                          # your GPU's sm_XY (nvidia-smi: compute_cap)
cargo build --release -p cfdcraft-cli --features cuda
target/release/cfdcraft scalar examples/cases/scalar/mms_poly_gpu.toml
target/release/cfdcraft scalar examples/cases/scalar/oblique_step_2d.toml --backends cpu,cuda
target/release/cfdcraft scalar examples/cases/scalar/mms_poly_box.toml --backends cuda
```

`mms_poly_gpu.toml` solves on the CPU and twice on the GPU (`run.backends =
["cpu", "cuda", "cuda"]`); the column "vs first backend" of the report says
"bitwise equal" when a GPU run reproduces the CPU run bit for bit, and the
run fails if any difference exceeds 1e-10. Keys: `run.backends`,
`run.rust_ptx` (or `--rust-ptx`).

[↑ back to top](#recipes)

---

<a id="poisson-study"></a>
## 15. Study the pressure solver on a Poisson problem

The mesh sweeps of stage 3 ([explanation](verification/stage3.md)): each
case solves the Poisson problem on three or four refinement levels with CG
and the multigrid (mixed precision and float64) and with Jacobi-CG, on 4
threads and on 1, prints the iteration ratio between the finest and the
coarsest level and fails if it exceeds 1.5.

```
cargo build --release -p cfdcraft-cli
for m in hex tet poly mixed; do
  target/release/cfdcraft poisson examples/cases/poisson/poisson_${m}_box.toml
done
cat output/poisson_poly_box/poisson.md      # iterations, times, hierarchy, differences
# plot: output/poisson_poly_box/iterations.svg
```

Try other multigrid settings from the command line or in `[linear.amg]`:

```
target/release/cfdcraft poisson examples/cases/poisson/poisson_poly_box.toml --solvers amg,amg-plain,jacobi
target/release/cfdcraft poisson examples/cases/poisson/poisson_poly_box.toml --levels 4 --threads 1 --repeat 3
```

```toml
[linear.amg]
smoother = "jacobi"        # or "chebyshev"
degree = 3                 # smoother steps before and after the coarse correction
strength = 0.02            # more aggressive coarsening
precision = "double"       # float64 hierarchy
```

The same solver in `cfdcraft scalar` (gallery 11):

```toml
[solver]
linear_solver = "cg"       # symmetric matrices only (no convection)
preconditioner = "amg"     # also with linear_solver = "bicgstab"
amg = { smoother = "chebyshev", degree = 2, precision = "mixed" }
```

The checks as tests (the mesh sweep in release builds):

```
cargo test --release -p cfdcraft-fv --test stage3
```

[↑ back to top](#recipes)

---

<a id="poisson-gpu"></a>
## 16. Run the pressure solver on an NVIDIA GPU against hypre

Without a GPU of your own: open
[`tools/stage3_colab.ipynb`](https://colab.research.google.com/github/chiefenne/CFDCraft/blob/main/tools/stage3_colab.ipynb)
in Google Colab, choose a T4 runtime, set `BRANCH` if the stage-3 branch is
not merged yet, and run all cells (about 30 minutes, most of it compiling
hypre with CUDA). On any Linux machine with an NVIDIA GPU of compute
capability 7.0 or newer and the CUDA toolkit:

```
bash tools/stage3_cloud.sh
```

By hand (sm_75 for a T4; `nvidia-smi --query-gpu=compute_cap --format=csv`):

```
tools/build_hypre.sh $PWD/target/hypre-cuda75 75
gpu/ptx-kernels/build.sh sm_75
HYPRE_DIR=$PWD/target/hypre-cuda75 cargo build --release -p cfdcraft-cli --features cuda,hypre
target/release/cfdcraft poisson examples/cases/poisson/poisson_poly_gpu.toml
target/release/cfdcraft poisson examples/cases/poisson/poisson_poly_box.toml --backends cuda --solvers amg,amg-double,hypre
```

`poisson_poly_gpu.toml` solves 912 673 polyhedra with the in-house solvers
on the CPU and on the GPU and with hypre on the GPU; the report gives the
setup and solve times and ns per cell and iteration of each, and the
column "Bitwise vs first backend" says whether the GPU solve equals the CPU
solve bit for bit. Keys: `study.backends`, `study.solvers`,
`study.rust_ptx` (or `--rust-ptx`), `study.repeat`.

[↑ back to top](#recipes)

---

<a id="hypre"></a>
## 17. Build hypre and compare with it

hypre (MIT or Apache-2.0) is the optional reference solver (D20). It is
built from source, sequential (no MPI), as a shared library; CFDCraft
links it only with the cargo feature `hypre` and finds it through
`HYPRE_DIR`. Needs git, CMake 3.21 or newer and a C compiler (on macOS:
`brew install cmake`).

```
tools/build_hypre.sh $PWD/target/hypre                 # about 1–3 minutes
export HYPRE_DIR=$PWD/target/hypre
cargo build --release -p cfdcraft-cli --features hypre
target/release/cfdcraft info                            # features: ..., hypre
for m in hex tet poly mixed; do
  target/release/cfdcraft poisson examples/cases/poisson/poisson_${m}_box.toml \
      --solvers amg,amg-double,amg-plain,jacobi,hypre --backends cpu,cpu-serial \
      --out output/poisson_${m}_hypre
done
cargo test --release -p cfdcraft-fv --features hypre --test stage3
```

With hypre in `study.solvers` the differences in the report are measured
against hypre's solution, and the run fails above `study.max_difference`
(1e-6). hypre runs on one CPU thread; compare its times with the
`cpu-serial` rows. BoomerAMG settings: `[hypre]` in the case file
([keys](cases.md#poisson-case-file-keys)).

[↑ back to top](#recipes)

---

<a id="flow-verification"></a>
## 18. Verify the flow solver on Couette and Poiseuille flow

The exact solutions of stage 4 ([explanation](verification/stage4.md)):
Couette flow must come out at the iteration error on every mesh type,
Poiseuille flow with an observed order of 1.8–2.2. Each case checks its own
result and fails otherwise.

```
cargo build --release -p cfdcraft-cli
for m in hex tet poly mixed; do
  target/release/cfdcraft flow examples/cases/flow/couette_${m}.toml
done
for m in hex tet poly mixed; do
  target/release/cfdcraft flow examples/cases/flow/poiseuille_${m}.toml
done
cat output/poiseuille_poly/flow.md          # errors per level, observed order, patches
# plots: output/poiseuille_poly/convergence.svg, residuals.svg
```

Keys used: `[analytic]` (`type = "couette"` or `"poiseuille"`),
`profile = "analytic"` at the velocity inlet, `[verification]` `levels`,
`expected_order`, `max_error`. More levels or other settings from the
command line:

```
target/release/cfdcraft flow examples/cases/flow/poiseuille_hex.toml --levels 0,1,2,3
target/release/cfdcraft flow examples/cases/flow/poiseuille_hex.toml --threads 1 --max-iterations 8000
```

The checks as tests (the solving ones in release builds):

```
cargo test --release -p cfdcraft-flow --test stage4
```

[↑ back to top](#recipes)

---

<a id="flow-cavity"></a>
## 19. Run the lid-driven cavity and compare with Ghia et al.

```
cargo build --release -p cfdcraft-cli
target/release/cfdcraft flow examples/cases/flow/cavity_re100.toml       # about 5 s
target/release/cfdcraft flow examples/cases/flow/cavity_re1000.toml      # about 1 min
target/release/cfdcraft flow examples/cases/flow/cavity_re100_poly3d.toml
cat output/cavity_re1000/flow.md
# centre lines against Ghia et al.: output/cavity_re1000/centreline_u.svg, centreline_v.svg
# velocity and pressure along the lines: line_vertical.csv, line_horizontal.csv
# wall pressure and shear stress: wall_ymax.csv (lid) and the other walls
pip install vtk
python3 tools/render_mesh.py output/cavity_re1000/fields.vtkhdf cavity.png --field speed --top --no-edges
```

Keys used: `[verification]` `reference = "ghia-re100"` or `"ghia-re1000"`,
`reference_tolerance` (lid speeds), `lid_velocity`; `[[output.lines]]`
for the sampling lines. The Reynolds number is ρ U L / μ with the unit
cavity; the case is refused if it does not match the reference.

Convergence: the outer iterations grow with the mesh and the Reynolds
number. A larger `relaxation_velocity` (0.9 to 0.95) needs fewer; lower it
(0.7, or `scheme = "upwind"` for the first iterations) if a run diverges:

```toml
[solver]
relaxation_velocity = 0.9
tolerance = 1e-8
max_iterations = 10000
```

[↑ back to top](#recipes)

---

<a id="flow-boundaries"></a>
## 20. Set up the boundary conditions of a flow case

Every patch of the mesh needs a `[boundary.<patch>]` table. Check a case
before running it: the check (D30) names each boundary that is missing,
unknown (with "did you mean"), over-specified or unsuitable, and says what
to change; `--check` stops after it.

```
target/release/cfdcraft mesh examples/cases/mesh/hex_2d.toml     # patch names: output/hex_2d/mesh_report.md
target/release/cfdcraft flow my_case.toml --check
```

A channel with a known flow rate:

```toml
[boundary.inlet]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]        # or normal_velocity = 1.0, or mass-flow-inlet
[boundary.outlet]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.walls]
type = "wall"
```

A pressure-driven flow (the flow rate is a result):

```toml
[boundary.inlet]
type = "total-pressure-inlet"
total_pressure = 10.0
[boundary.outlet]
type = "static-pressure-outlet"
pressure = 0.0
```

An enclosed flow (no pressure boundary: the reference cell holds the level):

```toml
[boundary.lid]
type = "wall"
velocity = [1.0, 0.0, 0.0]
[pressure]
reference_point = [0.5, 0.5, 0.0]
reference_value = 0.0
```

The table of sensible combinations and the cases each one suits is in the
[case catalogue](cases.md#sensible-combinations-d30); every condition and
key is listed in [boundary conditions](cases.md#boundary-conditions-d29).

[↑ back to top](#recipes)

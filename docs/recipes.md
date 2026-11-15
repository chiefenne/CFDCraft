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
21. [Run a transient flow case and check its time order](#transient-flow)
22. [Restart a run, stop it and steer it while it runs](#restart)
23. [Run the cylinder benchmark of Schäfer & Turek](#cylinder)
24. [Run a flow case on an NVIDIA GPU and compare it with the CPU](#flow-gpu) (also: free Colab T4)

---

<a id="g1-cpu"></a>
## 1. Build, test and run the G1 benchmark on your CPU

Needs Rust 1.89 or newer (`rustup update stable`), a C compiler and CMake:
the HDF5 library behind the VTKHDF writer (D41) and METIS (D32) are built
from source on the first build, which takes a few minutes (macOS:
`brew install cmake`; Ubuntu: `sudo apt install cmake`; Windows: Visual
Studio Build Tools and `winget install Kitware.CMake`). Without them,
`cargo build --release -p foehncraft-cli --no-default-features` builds
everything except the VTKHDF output and METIS (RCB partitioning remains).

```
git clone https://github.com/chiefenne/Foehncraft.git && cd Foehncraft
cargo test --workspace
cargo build --release -p foehncraft-cli
target/release/foehn info
target/release/foehn bench g1 examples/cases/bench/g1_ci_small.toml   # seconds
target/release/foehn bench g1 examples/cases/bench/g1_cpu_10m.toml    # about a minute, 4 GB
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
git clone https://github.com/chiefenne/Foehncraft.git
cd Foehncraft
cargo test --workspace --release
cargo build --release -p foehncraft-cli
target\release\foehn info
target\release\foehn bench g1 examples\cases\bench\g1_cpu_10m.toml
```

Close other programs first; run the last command twice. It prints
`G1 cpu-threads against cpu-native: PASS` or `FAIL` and writes
`output\g1_cpu_10m\g1.md`, `g1.json` and `g1.svg`.

The CUDA reference (needs an up-to-date NVIDIA driver and the CUDA 12
toolkit from https://developer.nvidia.com/cuda-downloads, for its NVRTC
library; on the GTX 960 use CUDA 12, because CUDA 13 dropped Maxwell):

```
cargo build --release -p foehncraft-cli --features cuda
target\release\foehn bench backends        # lists cuda-reference
target\release\foehn bench g1 examples\cases\bench\g1_gpu_reference_4m.toml
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
cargo build --release -p foehncraft-cli --features cuda
target/release/foehn bench backends         # must list cuda-reference and cuda-rust
target/release/foehn bench g1 examples/cases/bench/g1_gpu_10m.toml
```

The last lines print `G1 cuda-rust against cuda-reference: PASS` or `FAIL`.
Send `output/g1_gpu_10m/g1.md` and `g1.json`; they go into
[G1 results](benchmarks/g1.md#results). For the exact architecture:
`gpu/ptx-kernels/build.sh sm_86` (your GPU's compute capability). Keys:
`run.backends`, `run.device`, `run.rust_ptx`, `run.tolerance`.

On Windows, run the build script from Git Bash or WSL, or the two `cargo`
commands in it by hand.

**Without a GPU of your own: Google Colab (free T4).** Open
[`tools/g1_colab.ipynb` in Colab](https://colab.research.google.com/github/chiefenne/Foehncraft/blob/main/tools/g1_colab.ipynb),
choose *Runtime → Change runtime type → T4 GPU*, then *Runtime → Run all*.
About 10 minutes; the last cell shows both reports. On any other Linux
machine with a GPU (Kaggle with "GPU T4 x2", Lightning.ai):

```
git clone --depth 1 https://github.com/chiefenne/Foehncraft.git && cd Foehncraft
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
cargo test -p foehncraft-parallel --features mpi
cargo build --release -p foehncraft-parallel --features mpi --example mpi_smoke
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
tools/ci_cases.sh --check                               # every flow case is in a CI set
cargo clippy --workspace --all-targets -- -D warnings
cargo clippy --workspace --all-targets --features foehncraft-cuda/cuda,foehncraft-cli/cuda -- -D warnings
python3 tools/check_backend.py                          # D60: CUDA only in crates/foehncraft-cuda and gpu/
cargo test --workspace --profile ci                     # what a pull request runs (--release on main)
cargo test -p foehncraft-mesh --no-default-features      # without HDF5 and METIS
cargo test -p foehncraft-fv --no-default-features
gpu/ptx-kernels/build.sh && tools/check_kernels.sh     # PTX, and R10: no atomics, no f32, no FMA
cargo build --profile ci -p foehncraft-cli
for s in linux-1 linux-2 linux-3; do                    # add "long" for the cylinder at Re 100
  tools/ci_cases.sh "$s" target/ci/foehn || break
done
pip install -r tools/requirements-ci.txt
python3 tools/check_vtkhdf.py output/*/mesh.vtkhdf output/*/fields.vtkhdf output/*/*/fields.vtkhdf
```

The workflow is `.github/workflows/ci.yml`. It has two tiers:

- **Pull requests** (fast gate, about 5 minutes): the `ci` profile (release
  optimisation without LTO, results bitwise equal to `release`), all tests
  including the bitwise-reproducibility tests on Linux, macOS and Windows,
  the case catalogue in three parallel Linux jobs, the quick cases on macOS
  and Windows, hypre, MPI and PTX. Code jobs are skipped when only `docs/`,
  Markdown files, `assets/` or `LICENSE` change.
- **Full tier**: every push to main, every night at 02:07 Vienna time (01:07
  UTC), a pull request with the label `full-ci`, or by hand (Actions, CI,
  *Run workflow*): the `release` profile, debug tests, the cylinder at
  Re 100, and at night also the tests marked slow (`--ignored`).

The check "CI OK" sums up all jobs.

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

The release badge of the README updates by itself: the last step of the
workflow waits until img.shields.io reports the new tag and then purges
GitHub's cached copies of the README images (GitHub serves them through its
image proxy, which otherwise keeps the old badge for hours or days). If the
badge still shows the old version, run that purge by hand:

```
curl -s https://github.com/chiefenne/Foehncraft | grep -oE 'https://camo\.githubusercontent\.com/[0-9a-f]+/[0-9a-f]+' | sort -u | xargs -n1 curl -s -X PURGE
```

---

<a id="mesh-generate"></a>
## 7. Generate a test mesh and check its quality

The built-in mesher (D34) makes a box of hexahedra, tetrahedra, polyhedra
or mixed cells, in 3D or 2D (one cell thick), graded, refined and
perturbed. Keys `[mesh.generate]` ([Case catalogue](cases.md#mesh-case-file-keys)).

```
cargo build --release -p foehncraft-cli
target/release/foehn mesh examples/cases/mesh/poly_box.toml
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
target/release/foehn mesh examples/meshes/tjunction.msh             # defaults, no case file
target/release/foehn mesh examples/cases/mesh/tjunction_tet.toml    # RCM, 4 parts
target/release/foehn mesh examples/cases/mesh/tjunction_poly.toml   # converted into polyhedra
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
target/release/foehn mesh examples/cases/mesh/poly_partition.toml
target/release/foehn mesh examples/cases/mesh/poly_partition.toml --parts 8 --method rcb
```

Keys `partition.parts`, `partition.method` (`metis` or `rcb`),
`partition.halo_layers`. The report lists cells, halo cells, processor
faces and neighbours per rank, the imbalance and the cut faces;
`mesh.vtkhdf` has the field `rank`.

The halo exchange over MPI (needs MPI as in recipe 4): every rank builds
its part with a two-layer halo and fills it from its neighbours.

```
cargo build --release -p foehncraft-mesh --features mpi --example mpi_halo
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
git pull && cargo build --release -p foehncraft-cli
for c in g1_cpu_poly_2m g1_cpu_poly_2m_shuffled g1_cpu_tet_2m; do
  target/release/foehn bench g1 examples/cases/bench/$c.toml
  target/release/foehn bench g1 examples/cases/bench/$c.toml   # second run
done
```

Send `output/g1_cpu_*/g1.md`. The key `mesh.type` (`poly`, `tet`, `hex`,
`mixed`, `periodic-box`) chooses the mesh, `mesh.cells` its size. On a GPU
of compute capability 7.0 or newer (or Colab, recipe 3):
`foehn bench g1 examples/cases/bench/g1_gpu_poly_4m.toml`.

[↑ back to top](#recipes)

---

<a id="golden"></a>
## 11. Update the golden values of the mesher

`crates/foehncraft-mesh/tests/golden/mesher.json` holds the cell, face and
node counts, volumes, quality metrics and partition statistics of every
combination of mesh type, dimension, perturbation, refinement and 1, 2 and
4 ranks (R4). A difference fails the test. Only for an intended change:

```
FOEHNCRAFT_BLESS=1 cargo test -p foehncraft-mesh --test golden
git diff crates/foehncraft-mesh/tests/golden/mesher.json
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
cargo build --release -p foehncraft-cli
for m in hex tet poly mixed; do
  target/release/foehn scalar examples/cases/scalar/mms_${m}_box.toml
done
target/release/foehn scalar examples/cases/scalar/mms_poly_2d.toml
cat output/mms_poly_box/scalar.md          # errors, orders, time per phase
```

Keys: `verification.levels` (refinement levels of `[mesh.generate]`,
override with `--levels 1,2,3,4`), `verification.expected_order`,
`problem.manufactured`, `problem.boundary.<patch>` with `profile =
"manufactured"`, `mesh.generate.perturb`. The same checks as tests (one
level finer with `--ignored`):

```
cargo test --release -p foehncraft-fv
cargo test --release -p foehncraft-fv -- --ignored --nocapture   # about 1.6 M tetrahedra, minutes
cargo run --release -p foehncraft-fv --example sweep_table         # the table of the explanation page
```

[↑ back to top](#recipes)

---

<a id="scalar-schemes"></a>
## 13. Compare convection schemes and limiters

One case, several solver settings: `[[variants]]` runs the case once per
variant, writes the fields of each into its own directory and the profiles
along `output.profile` into one plot.

```
target/release/foehn scalar examples/cases/scalar/oblique_step_2d.toml
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
[`tools/stage2_colab.ipynb`](https://colab.research.google.com/github/chiefenne/Foehncraft/blob/main/tools/stage2_colab.ipynb)
in Google Colab, choose a T4 runtime, set `BRANCH` if the stage-2 branch is
not merged yet, and run all cells (about 15 minutes). On any Linux machine
with an NVIDIA GPU of compute capability 7.0 or newer:

```
bash tools/stage2_cloud.sh
```

By hand:

```
gpu/ptx-kernels/build.sh sm_75                          # your GPU's sm_XY (nvidia-smi: compute_cap)
cargo build --release -p foehncraft-cli --features cuda
target/release/foehn scalar examples/cases/scalar/mms_poly_gpu.toml
target/release/foehn scalar examples/cases/scalar/oblique_step_2d.toml --backends cpu,cuda
target/release/foehn scalar examples/cases/scalar/mms_poly_box.toml --backends cuda
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
cargo build --release -p foehncraft-cli
for m in hex tet poly mixed; do
  target/release/foehn poisson examples/cases/poisson/poisson_${m}_box.toml
done
cat output/poisson_poly_box/poisson.md      # iterations, times, hierarchy, differences
# plot: output/poisson_poly_box/iterations.svg
```

Try other multigrid settings from the command line or in `[linear.amg]`:

```
target/release/foehn poisson examples/cases/poisson/poisson_poly_box.toml --solvers amg,amg-plain,jacobi
target/release/foehn poisson examples/cases/poisson/poisson_poly_box.toml --levels 4 --threads 1 --repeat 3
```

```toml
[linear.amg]
smoother = "jacobi"        # or "chebyshev"
degree = 3                 # smoother steps before and after the coarse correction
strength = 0.02            # more aggressive coarsening
precision = "double"       # float64 hierarchy
```

The same solver in `foehn scalar` (gallery 11):

```toml
[solver]
linear_solver = "cg"       # symmetric matrices only (no convection)
preconditioner = "amg"     # also with linear_solver = "bicgstab"
amg = { smoother = "chebyshev", degree = 2, precision = "mixed" }
```

The checks as tests (the mesh sweep in release builds):

```
cargo test --release -p foehncraft-fv --test stage3
```

[↑ back to top](#recipes)

---

<a id="poisson-gpu"></a>
## 16. Run the pressure solver on an NVIDIA GPU against hypre

Without a GPU of your own: open
[`tools/stage3_colab.ipynb`](https://colab.research.google.com/github/chiefenne/Foehncraft/blob/main/tools/stage3_colab.ipynb)
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
HYPRE_DIR=$PWD/target/hypre-cuda75 cargo build --release -p foehncraft-cli --features cuda,hypre
target/release/foehn poisson examples/cases/poisson/poisson_poly_gpu.toml
target/release/foehn poisson examples/cases/poisson/poisson_poly_box.toml --backends cuda --solvers amg,amg-double,hypre
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
built from source, sequential (no MPI), as a shared library; Foehncraft
links it only with the cargo feature `hypre` and finds it through
`HYPRE_DIR`. Needs git, CMake 3.21 or newer and a C compiler (on macOS:
`brew install cmake`).

```
tools/build_hypre.sh $PWD/target/hypre                 # about 1–3 minutes
export HYPRE_DIR=$PWD/target/hypre
cargo build --release -p foehncraft-cli --features hypre
target/release/foehn info                            # features: ..., hypre
for m in hex tet poly mixed; do
  target/release/foehn poisson examples/cases/poisson/poisson_${m}_box.toml \
      --solvers amg,amg-double,amg-plain,jacobi,hypre --backends cpu,cpu-serial \
      --out output/poisson_${m}_hypre
done
cargo test --release -p foehncraft-fv --features hypre --test stage3
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
cargo build --release -p foehncraft-cli
for m in hex tet poly mixed; do
  target/release/foehn flow examples/cases/flow/couette_${m}.toml
done
for m in hex tet poly mixed; do
  target/release/foehn flow examples/cases/flow/poiseuille_${m}.toml
done
cat output/poiseuille_poly/flow.md          # errors per level, observed order, patches
# plots: output/poiseuille_poly/convergence.svg, residuals.svg
```

Keys used: `[analytic]` (`type = "couette"` or `"poiseuille"`),
`profile = "analytic"` at the velocity inlet, `[verification]` `levels`,
`expected_order`, `max_error`. More levels or other settings from the
command line:

```
target/release/foehn flow examples/cases/flow/poiseuille_hex.toml --levels 0,1,2,3
target/release/foehn flow examples/cases/flow/poiseuille_hex.toml --threads 1 --max-iterations 8000
```

The checks as tests (the solving ones in release builds):

```
cargo test --release -p foehncraft-flow --test stage4
```

[↑ back to top](#recipes)

---

<a id="flow-cavity"></a>
## 19. Run the lid-driven cavity and compare with Ghia et al.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cavity_re100.toml       # about 5 s
target/release/foehn flow examples/cases/flow/cavity_re1000.toml      # about 1 min
target/release/foehn flow examples/cases/flow/cavity_re100_poly3d.toml
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
target/release/foehn mesh examples/cases/mesh/hex_2d.toml     # patch names: output/hex_2d/mesh_report.md
target/release/foehn flow my_case.toml --check
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

---

<a id="transient-flow"></a>
## 21. Run a transient flow case and check its time order

A case is transient when it has a `[time]` table; it then runs PISO with
BDF2 (or Crank–Nicolson) instead of the steady SIMPLEC iterations.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/taylor_green_bdf2.toml    # about 40 s
target/release/foehn flow examples/cases/flow/taylor_green_cn.toml
target/release/foehn flow examples/cases/flow/taylor_green_space.toml   # space order, about 30 s
target/release/foehn flow examples/cases/flow/taylor_green_cfl10.toml   # Courant 1 to 9, outer iterations
cat output/taylor_green_bdf2/flow.md
# plot: output/taylor_green_bdf2/time_convergence.svg, per step: history.csv
```

The smallest transient case file adds this to a steady one:

```toml
[time]
scheme = "bdf2"          # or "crank-nicolson", "bdf1"
time_step = 0.01         # s
end_time = 2.0           # or steps = 200
adaptive = true          # optional: Δt from the largest Courant number
max_cfl = 0.8
```

Plain PISO (one outer iteration, two correctors) is second order in time
up to Courant numbers of about 1. For larger steps add outer iterations,
which end when the residuals of a step are below `outer_tolerance`:

```toml
[solver]
outer_iterations = 50
outer_tolerance = 1e-8
relaxation_velocity = 0.9   # below 1 only with outer_iterations > 1
```

A time-order sweep runs the case once per step of
`verification.time_steps` (largest first) on the same mesh and reports
the order from the differences of successive solutions; with
`expected_time_order` it fails outside that range:

```toml
[verification]
time_steps = [0.01, 0.005, 0.0025, 0.00125]
expected_time_order = [1.8, 2.3]
```

Keys: [transient cases](cases.md#flow-case-file-keys); the method and the
results: [stage 5](verification/stage5.md).

[↑ back to top](#recipes)

---

<a id="restart"></a>
## 22. Restart a run, stop it and steer it while it runs

Every run writes `restart.h5` at its end (`run.restart_final`), and on a
schedule if asked. A run continued from it gives bitwise the result of
the run without the interruption (same mesh, same build):

```toml
[run]
restart_every = 500          # restart.h5 every 500 steps
restart_at = [2000, 4000]    # also restart_002000.h5, restart_004000.h5
```

```
target/release/foehn flow my_case.toml
# continue to a later end time, into another directory:
target/release/foehn flow my_case.toml --restart output/my_case/restart.h5 --end-time 20 --out output/my_case_2
```

or, in the case file, `[run] restart = "output/my_case/restart.h5"` (a
path in the case file is relative to the case file, on the command line
to the current directory).

While a case runs, files in its output directory steer it (the solver
never deletes them; remove them before the next run):

```
touch output/my_case/.stop            # write restart.h5 and the results, then end
touch output/my_case/.result_1500     # write fields_001500.vtkhdf at step 1500
touch output/my_case/.read_toml_1200  # re-read my_case.toml at step 1200
```

A re-read applies the output and restart schedule, the stop criteria, the
relaxation, the linear-solver tolerances and the time-step settings; the
log lists each change (`key: old -> new`) and each change it cannot apply
during a run. `run.reread_every = N` re-reads every N steps without a
control file.

[↑ back to top](#recipes)

---

<a id="cylinder"></a>
## 23. Run the cylinder benchmark of Schäfer & Turek

```
target/release/foehn flow examples/cases/flow/cylinder_re20.toml     # steady, Re 20, about 1 min
target/release/foehn flow examples/cases/flow/cylinder_re100.toml    # periodic, Re 100, about 8 min
cat output/cylinder_re100/flow.md
# C_D, C_L per step: output/cylinder_re100/forces.csv, chart: coefficients.svg
# probes in front of and behind the cylinder: probes.csv
pip install vtk
python3 tools/render_mesh.py output/cylinder_re20/fields.vtkhdf cylinder.png --field speed --top --no-edges --size 1800x420 --zoom 2.9
```

The mesh is built in (`[mesh.cylinder]`); `refine = 1` and `2` give 24 576
and 98 304 cells. Re 100 with `refine = 1` takes about 40 minutes on two
cores and brings C_L,max to 0.5 % of its range; with a finer mesh lower the
first `time_step` (the adaptive step starts from it). On 98 304 cells the steady case needs
`relaxation_velocity = 0.5` and about 12 000 outer iterations:

```toml
[mesh.cylinder]
refine = 2
[solver]
relaxation_velocity = 0.5
max_iterations = 12000
```

Keys: `[verification]` `benchmark = "schafer-turek-2d1"` or
`"schafer-turek-2d2"`, `benchmark_tolerance` (relative to the midpoints of
the reference ranges), `evaluation_start` (2D-2). Force coefficients of
any case: `[output.coefficients]` with `patches`, `velocity`, `area`;
point histories: `[[output.probes]]` with `name` and `point`.

[↑ back to top](#recipes)

---

<a id="flow-gpu"></a>
## 24. Run a flow case on an NVIDIA GPU and compare it with the CPU

Every flow case runs on the GPU as it is: the whole flow loop (gradients,
limiters, momentum and pressure-correction assembly, Rhie–Chow fluxes,
corrections, the BiCGStab and CG + multigrid solves) stays on the device,
and the results are bitwise those of the CPU (D31). Needs an NVIDIA GPU of
compute capability 7.0 or newer with its driver; the build needs no CUDA
toolkit, the kernels need the nightly toolchain (installed by
`rust-toolchain.toml` on first use). From the repository root:

```
gpu/ptx-kernels/build.sh sm_75                       # PTX for your GPU (sm_75: T4, sm_86: RTX 30xx, ...)
cargo build --release -p foehncraft-cli --features cuda
target/release/foehn info                             # features: cuda ...
target/release/foehn flow examples/cases/flow/cavity_re100.toml --backend cuda
```

To check the GPU against the CPU, run the case twice and compare the output
directories (`pip install h5py numpy` once):

```
target/release/foehn flow examples/cases/flow/cylinder_re100_gpu.toml --backend cpu  --out output/gpu/cpu
target/release/foehn flow examples/cases/flow/cylinder_re100_gpu.toml --backend cuda --out output/gpu/cuda
python3 tools/compare_runs.py output/gpu/cpu output/gpu/cuda    # "bitwise equal", exit code 0
```

`python3 tools/stage5_cloud.py` does this for eight cases, from Couette
flow to the cavity on one million polyhedra, and writes
`output/stage5_gpu/report.md` with the times per cell and iteration.
**Without a GPU of your own:** open
[`tools/stage5_colab.ipynb`](https://colab.research.google.com/github/chiefenne/Foehncraft/blob/main/tools/stage5_colab.ipynb)
in Google Colab, choose a T4 runtime, set `BRANCH` if the stage-5 branch
is not merged yet, and run all cells (about 15 minutes).

Keys: `run.backend` (`cpu`, `cuda`), `run.ptx` (the PTX module); the
command-line options `--backend` and `--ptx` override them. A build
without the `cuda` feature, or a machine without the NVIDIA driver, stops
with a message instead of running on the CPU.

[↑ back to top](#recipes)

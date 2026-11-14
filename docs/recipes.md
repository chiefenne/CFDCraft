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

---

<a id="g1-cpu"></a>
## 1. Build, test and run the G1 benchmark on your CPU

Needs Rust 1.89 or newer (`rustup update stable`).

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
cargo clippy --workspace --all-targets --features cfdcraft-bench/cuda -- -D warnings
cargo test --workspace --release
gpu/ptx-kernels/build.sh
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

[↑ back to top](#recipes)

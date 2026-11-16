# Kernel benchmark G1

Part of the [Recipes](../recipes.md): build and run gate G1 on a CPU, on Windows and on an NVIDIA GPU, and rerun it on unstructured meshes.

| # | Recipe |
|---|---|
| 1 | [Build, test and run the G1 benchmark on your CPU](#g1-cpu) |
| 2 | [Run G1 on Windows: CPU verdict and the CUDA reference](#g1-windows) |
| 3 | [Run gate G1 on an NVIDIA GPU](#g1-gpu) (also: free Colab T4) |
| 10 | [Rerun G1 on polyhedral and tetrahedral meshes](#g1-stage1) |

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

[↑ back to top](#kernel-benchmark-g1)

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
files; they go into [G1 results](../benchmarks/g1.md#maintainers-machine).
Keys: `run.backends`, `run.threads`, `mesh.cells`.

[↑ back to top](#kernel-benchmark-g1)

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
[G1 results](../benchmarks/g1.md#results). For the exact architecture:
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

[↑ back to top](#kernel-benchmark-g1)

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

[↑ back to top](#kernel-benchmark-g1)

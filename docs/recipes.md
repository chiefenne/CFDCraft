# Recipes

Short command sequences for everyday tasks. Each recipe names the case-file
keys it uses; the keys of the benchmark cases are listed in the
[Case catalogue](cases.md#case-file-keys). Pictures of complete cases are in
the [Example gallery](gallery.md).

1. [Build, test and run the G1 benchmark on your CPU](#g1-cpu)
2. [Run gate G1 on an NVIDIA GPU](#g1-gpu)
3. [Run the MPI test with several ranks on one machine](#mpi)
4. [Check a change the way CI does](#ci-local)

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

<a id="g1-gpu"></a>
## 2. Run gate G1 on an NVIDIA GPU

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

[↑ back to top](#recipes)

---

<a id="mpi"></a>
## 3. Run the MPI test with several ranks on one machine

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
## 4. Check a change the way CI does

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

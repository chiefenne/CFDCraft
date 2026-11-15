#!/usr/bin/env bash
# Gate G1 on a cloud GPU (Google Colab, Kaggle, Lightning.ai, any Linux
# machine with an NVIDIA GPU of compute capability 7.0 or newer and the CUDA
# 12 libraries). Run from the repository root:
#
#   bash tools/g1_cloud.sh [CASE.toml]     # default: examples/cases/bench/g1_gpu_10m.toml
#
# Installs Rust if needed, builds the single-source PTX kernels for the GPU
# found, builds the CLI with CUDA, runs the case twice and prints both
# reports. Takes about 10 minutes on a free Colab T4, most of it compiling.
set -euo pipefail
case_file=${1:-examples/cases/bench/g1_gpu_10m.toml}

echo "== GPU"
nvidia-smi --query-gpu=name,compute_cap,memory.total,driver_version --format=csv,noheader
cc=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d ' .')
if [ "$cc" -lt 70 ]; then
    echo "compute capability ${cc:0:1}.${cc:1} is older than 7.0: the single-source Rust route cannot run here (O13)." >&2
    echo "On Kaggle choose the accelerator 'GPU T4 x2', not 'P100'." >&2
    exit 1
fi

echo "== Rust"
if ! command -v cargo >/dev/null 2>&1; then
    [ -f "$HOME/.cargo/env" ] || curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
    . "$HOME/.cargo/env"
fi
cargo --version

# NVRTC (for the CUDA C reference) from the CUDA toolkit of the image
for d in /usr/local/cuda/lib64 /usr/local/cuda/targets/x86_64-linux/lib; do
    [ -d "$d" ] && export LD_LIBRARY_PATH="$d:${LD_LIBRARY_PATH:-}"
done

echo "== Single-source kernels -> PTX for sm_$cc"
gpu/ptx-kernels/build.sh "sm_$cc"

echo "== CLI with CUDA"
cargo build --release --locked -p inncraft-cli --features cuda
target/release/inncraft info

for run in 1 2; do
    echo "== G1 run $run of 2: $case_file"
    target/release/inncraft bench g1 "$case_file" --out "output/g1_cloud_run$run"
done

echo
echo "== Reports (copy everything below into the project thread)"
for run in 1 2; do
    echo
    echo "### Run $run"
    cat "output/g1_cloud_run$run/g1.md"
done

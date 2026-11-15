#!/usr/bin/env bash
# Stage 2 on a cloud GPU (Google Colab, Kaggle, any Linux machine with an
# NVIDIA GPU of compute capability 7.0 or newer). Run from the repository
# root:
#
#   bash tools/stage2_cloud.sh
#
# Installs Rust if needed, builds the single-source PTX kernels for the GPU
# found, builds the CLI with CUDA and runs the stage-2 GPU checks: the
# polyhedral manufactured-solution case on about 0.9 million cells (CPU and
# twice on the GPU: GPU against CPU and GPU run to run), the oblique step
# with all four schemes on CPU and GPU, and the polyhedral order-of-accuracy
# sweep on the GPU alone. Takes about 15 minutes on a free Colab T4, most of
# it compiling.
set -euo pipefail

echo "== GPU"
nvidia-smi --query-gpu=name,compute_cap,memory.total,driver_version --format=csv,noheader
cc=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d ' .')
if [ "$cc" -lt 70 ]; then
    echo "compute capability ${cc:0:1}.${cc:1} is older than 7.0: the single-source Rust route cannot run here (O13)." >&2
    exit 1
fi

echo "== Rust"
if ! command -v cargo >/dev/null 2>&1; then
    [ -f "$HOME/.cargo/env" ] || curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
    . "$HOME/.cargo/env"
fi
cargo --version
for d in /usr/local/cuda/lib64 /usr/local/cuda/targets/x86_64-linux/lib; do
    [ -d "$d" ] && export LD_LIBRARY_PATH="$d:${LD_LIBRARY_PATH:-}"
done

echo "== Single-source kernels -> PTX for sm_$cc"
gpu/ptx-kernels/build.sh "sm_$cc"

echo "== CLI with CUDA"
cargo build --release --locked -p ennecraft-cli --features cuda
target/release/ennecraft info

status=0
run() {
    echo "== ennecraft scalar $*"
    target/release/ennecraft scalar "$@" || status=1
}
run examples/cases/scalar/mms_poly_gpu.toml --out output/stage2_gpu_mms
run examples/cases/scalar/oblique_step_2d.toml --backends cpu,cuda --out output/stage2_gpu_oblique
run examples/cases/scalar/mms_poly_box.toml --backends cuda --out output/stage2_gpu_sweep

echo
echo "== Reports (copy everything below into the project thread)"
for d in stage2_gpu_mms stage2_gpu_oblique stage2_gpu_sweep; do
    echo
    cat "output/$d/scalar.md"
done
exit $status

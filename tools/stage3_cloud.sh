#!/usr/bin/env bash
# Stage 3 on a cloud GPU (Google Colab, Kaggle, any Linux machine with an
# NVIDIA GPU of compute capability 7.0 or newer and the CUDA toolkit). Run
# from the repository root:
#
#   bash tools/stage3_cloud.sh
#
# Installs Rust if needed, builds hypre with CUDA for the GPU found (the
# reference, D20), builds the single-source PTX kernels and the CLI with
# CUDA and hypre, and runs the stage-3 GPU checks:
#
# 1. the pressure solver on about 0.9 million polyhedra: CG with the
#    smoothed-aggregation multigrid (mixed precision and float64) and
#    Jacobi-CG on the CPU and on the GPU, hypre's BoomerAMG-PCG on the GPU;
#    times per cell and iteration, GPU against CPU bitwise;
# 2. the polyhedral mesh sweep (729 to 274 625 cells) on the GPU with
#    hypre: iterations nearly independent of the mesh size;
# 3. the stage-3 tests on the CPU against the GPU build of hypre.
#
# Takes about 30 minutes on a free Colab T4, most of it compiling hypre.
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
[ -d /usr/local/cuda/bin ] && export PATH="/usr/local/cuda/bin:$PATH"

echo "== hypre with CUDA for sm_$cc (about 20 minutes)"
hypre=$PWD/target/hypre-cuda$cc
if [ ! -f "$hypre/include/HYPRE.h" ]; then
    tools/build_hypre.sh "$hypre" "$cc"
fi
export HYPRE_DIR=$hypre

echo "== Single-source kernels -> PTX for sm_$cc"
gpu/ptx-kernels/build.sh "sm_$cc"

echo "== CLI with CUDA and hypre"
cargo build --release --locked -p cfdcraft-cli --features cuda,hypre
target/release/cfdcraft info

status=0
run() {
    echo "== cfdcraft poisson $*"
    target/release/cfdcraft poisson "$@" || status=1
}
run examples/cases/poisson/poisson_poly_gpu.toml --out output/stage3_gpu_poly
run examples/cases/poisson/poisson_poly_box.toml --backends cuda --solvers amg,amg-double,hypre \
    --out output/stage3_gpu_sweep

echo "== Stage-3 tests against the GPU build of hypre"
cargo test --release --locked -p cfdcraft-fv --features hypre --test stage3 || status=1

echo
echo "== Reports (copy everything below into the project thread)"
for d in stage3_gpu_poly stage3_gpu_sweep; do
    echo
    cat "output/$d/poisson.md"
done
exit $status

#!/usr/bin/env sh
# Builds the single-source kernels (G1, stage 2 and stage 3) as a PTX module for
# NVIDIA GPUs.
#
#   gpu/ptx-kernels/build.sh            # PTX for sm_70 (Volta), JIT-compiled for newer GPUs
#   gpu/ptx-kernels/build.sh sm_86      # PTX for a given architecture (rebuilds `core`)
#
# Needs rustup; the nightly toolchain, the nvptx64-nvidia-cuda target and the
# llvm-bitcode-linker component are installed on first use
# (rust-toolchain.toml). The rustc NVPTX target supports sm_70 and newer.
# Output: target/nvptx64-nvidia-cuda/release/inncraft_ptx_kernels.ptx in the
# repository root (the default `run.rust_ptx` of the G1 case files).
set -eu
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
cd "$here"
if [ "${1:-}" != "" ] && [ "$1" != "sm_70" ]; then
    rustup component add rust-src --toolchain nightly >/dev/null
    RUSTFLAGS="-C target-cpu=$1" cargo build --release --target nvptx64-nvidia-cuda \
        -Zbuild-std=core --target-dir "$root/target" ${CARGO_FEATURES:-}
else
    cargo build --release --target nvptx64-nvidia-cuda --target-dir "$root/target" ${CARGO_FEATURES:-}
fi
ptx="$root/target/nvptx64-nvidia-cuda/release/inncraft_ptx_kernels.ptx"
for k in g1_face_flux g1_cell_gather g1_spmv_csr g1_spmv_sell \
    fv_gradient fv_limiter fv_faces fv_assemble fv_residual \
    blas_block_dot blas_block_min_max blas_update_p blas_sub_scaled blas_update_x blas_scale \
    blas_axpy blas_xpby \
    amg_smooth_first_f32 amg_smooth_step_f32 amg_add_f32 amg_residual_f32 amg_spmv_f32 \
    amg_prolong_f32 amg_dense_f32 amg_to_f32 amg_from_f32 \
    amg_smooth_first_f64 amg_smooth_step_f64 amg_add_f64 amg_residual_f64 amg_spmv_f64 \
    amg_prolong_f64 amg_dense_f64 amg_to_f64 amg_from_f64; do
    grep -q ".entry $k" "$ptx" || { echo "missing entry point $k in $ptx" >&2; exit 1; }
done
echo "wrote $ptx ($(grep -m1 '^\.target' "$ptx"))"

#!/usr/bin/env sh
# R10: the kernels stay reproducible (D31) and on the float64 path (D7).
#
#   tools/check_kernels.sh [PTX]   # default: the module of gpu/ptx-kernels/build.sh
#
# Fails when
# - the kernel sources use atomics (no floating-point atomics, no atomic
#   scatter: every output is gathered by one index);
# - the kernel sources use f32 (the reference path is float64);
# - the PTX module contains atomic or reduction instructions (atom., red.),
#   fused multiply-adds on f64 (fma.rn.f64: results would differ from the
#   CPU) or approximate f64 division.
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
ptx=${1:-$root/target/nvptx64-nvidia-cuda/release/cfdcraft_ptx_kernels.ptx}
fail=0
src="$root/crates/cfdcraft-kernels/src $root/gpu/ptx-kernels/src"
if grep -rnE 'Atomic|atomic_|core::intrinsics::atomic' $src; then
    echo "R10: atomics in the kernel sources (above)" >&2; fail=1
fi
if grep -rnE '\bf32\b' $src | grep -vE '^\S+:\s*//|impl Real for f32|//!|///'; then
    echo "R10: f32 in the kernel sources (above)" >&2; fail=1
fi
if [ -f "$ptx" ]; then
    if grep -nE '^\s*(atom|red)\.' "$ptx"; then
        echo "R10: atomic instructions in $ptx (above)" >&2; fail=1
    fi
    if grep -nE 'fma\.rn\.f64|div\.(approx|full)\.f64' "$ptx"; then
        echo "R10: contracted or approximate f64 arithmetic in $ptx (above)" >&2; fail=1
    fi
    echo "PTX checked: $ptx"
else
    echo "no PTX module at $ptx (build it with gpu/ptx-kernels/build.sh); sources checked only"
fi
[ $fail -eq 0 ] && echo "R10 kernel checks passed"
exit $fail

#!/usr/bin/env sh
# R10: the kernels stay reproducible (D31) and on the float64 path (D7).
#
#   tools/check_kernels.sh [PTX]   # default: the module of gpu/ptx-kernels/build.sh
#
# Fails when
# - the kernel sources use atomics (no floating-point atomics, no atomic
#   scatter: every output is gathered by one index);
# - the kernel sources use f32 outside the multigrid modules (`amg.rs`: the
#   float32 hierarchy of D37; everything else stays on the float64 path);
# - the PTX module contains atomic or reduction instructions (atom., red.),
#   fused multiply-adds (fma.rn.f64, fma.rn.f32: results would differ from
#   the CPU), approximate division or reciprocals, or flush-to-zero f32
#   arithmetic (.ftz).
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
ptx=${1:-$root/target/nvptx64-nvidia-cuda/release/ennecraft_ptx_kernels.ptx}
fail=0
src="$root/crates/ennecraft-kernels/src $root/gpu/ptx-kernels/src"
if grep -rnE 'Atomic|atomic_|core::intrinsics::atomic' $src; then
    echo "R10: atomics in the kernel sources (above)" >&2; fail=1
fi
if grep -rnE '\bf32\b' $src | grep -vE '/amg\.rs:' | grep -vE '^\S+:\s*//|impl Real for f32|//!|///'; then
    echo "R10: f32 outside the multigrid kernels (above)" >&2; fail=1
fi
if [ -f "$ptx" ]; then
    if grep -nE '^\s*(atom|red)\.' "$ptx"; then
        echo "R10: atomic instructions in $ptx (above)" >&2; fail=1
    fi
    if grep -nE 'fma\.rn\.f(64|32)|div\.(approx|full)\.f(64|32)|rcp\.approx|\.ftz\.' "$ptx"; then
        echo "R10: contracted, approximate or flush-to-zero arithmetic in $ptx (above)" >&2; fail=1
    fi
    echo "PTX checked: $ptx"
else
    echo "no PTX module at $ptx (build it with gpu/ptx-kernels/build.sh); sources checked only"
fi
[ $fail -eq 0 ] && echo "R10 kernel checks passed"
exit $fail

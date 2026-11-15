#!/usr/bin/env bash
# Builds hypre (https://github.com/hypre-space/hypre, MIT or Apache-2.0),
# the optional reference solver of `ennecraft poisson` (D20), from source:
# sequential (no MPI), shared library, release.
#
#   tools/build_hypre.sh                 # CPU build into target/hypre
#   tools/build_hypre.sh PREFIX          # CPU build into PREFIX
#   tools/build_hypre.sh PREFIX 75       # CUDA build for sm_75 (Tesla T4)
#
# Needs git, CMake 3.21 or newer and a C compiler (and the CUDA toolkit for
# a GPU build). Then build Ennecraft with
#
#   HYPRE_DIR=PREFIX cargo build --release -p ennecraft-cli --features hypre
#
# HYPRE_VERSION selects the release tag (default v3.2.1).
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
prefix=${1:-$root/target/hypre}
arch=${2:-}
version=${HYPRE_VERSION:-v3.2.1}
src=$root/target/hypre-src-$version
# A source tree restored from a cache can be incomplete: re-clone unless the
# CMake project is really there.
if [ ! -f "$src/src/CMakeLists.txt" ]; then
    rm -rf "$src"
    git -c advice.detachedHead=false clone --quiet --depth 1 -b "$version" https://github.com/hypre-space/hypre "$src"
fi
build=$src/build${arch:+-cuda$arch}
opts=(-DHYPRE_ENABLE_MPI=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release
      -DCMAKE_INSTALL_PREFIX="$prefix")
if [ -n "$arch" ]; then
    # Umpire (a GPU memory pool) is optional; without it hypre allocates with
    # cudaMalloc, which is what the Colab build needs (no Umpire on the image).
    opts+=(-DHYPRE_ENABLE_CUDA=ON -DHYPRE_ENABLE_UMPIRE=OFF -DCMAKE_CUDA_ARCHITECTURES="$arch")
fi
cmake -S "$src/src" -B "$build" "${opts[@]}" >/dev/null
cmake --build "$build" -j "$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)" >/dev/null
cmake --install "$build" >/dev/null
echo "hypre $version installed in $prefix${arch:+ (CUDA, sm_$arch)}"

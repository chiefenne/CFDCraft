"""Builds hypre (https://github.com/hypre-space/hypre, MIT or Apache-2.0),
the optional reference solver of `foehn poisson` (D20), from source:
sequential (no MPI), shared library, release.

    python -m tools.bench.build_hypre                 # CPU build into target/hypre
    python -m tools.bench.build_hypre PREFIX          # CPU build into PREFIX
    python -m tools.bench.build_hypre PREFIX 75       # CUDA build for sm_75 (Tesla T4)

Needs git, CMake 3.21 or newer and a C compiler (and the CUDA toolkit for
a GPU build). Then build Foehncraft with

    HYPRE_DIR=PREFIX cargo build --release -p foehncraft-cli --features hypre

HYPRE_VERSION selects the release tag (default v3.2.1).
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run(cmd, **kw):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, **kw)


def main(argv) -> int:
    prefix = Path(argv[0]).resolve() if argv else ROOT / "target" / "hypre"
    arch = argv[1] if len(argv) > 1 else ""
    version = os.environ.get("HYPRE_VERSION", "v3.2.1")
    src = ROOT / "target" / f"hypre-src-{version}"
    # a source tree restored from a cache can be incomplete: re-clone unless
    # the CMake project is really there
    if not (src / "src" / "CMakeLists.txt").is_file():
        shutil.rmtree(src, ignore_errors=True)
        subprocess.run(["git", "-c", "advice.detachedHead=false", "clone", "--quiet", "--depth", "1",
                        "-b", version, "https://github.com/hypre-space/hypre", str(src)], check=True)
    build = src / (f"build-cuda{arch}" if arch else "build")
    opts = ["-DHYPRE_ENABLE_MPI=OFF", "-DBUILD_SHARED_LIBS=ON", "-DCMAKE_BUILD_TYPE=Release",
            f"-DCMAKE_INSTALL_PREFIX={prefix}"]
    if arch:
        # Umpire (a GPU memory pool) is optional; without it hypre allocates
        # with cudaMalloc, which is what the Colab build needs
        opts += ["-DHYPRE_ENABLE_CUDA=ON", "-DHYPRE_ENABLE_UMPIRE=OFF", f"-DCMAKE_CUDA_ARCHITECTURES={arch}"]
    run(["cmake", "-S", str(src / "src"), "-B", str(build), *opts])
    run(["cmake", "--build", str(build), "-j", str(os.cpu_count() or 2)])
    run(["cmake", "--install", str(build)])
    print(f"hypre {version} installed in {prefix}" + (f" (CUDA, sm_{arch})" if arch else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

"""Stage 3 on a cloud GPU (Google Colab, Kaggle, any Linux machine with an
NVIDIA GPU of compute capability 7.0 or newer and the CUDA toolkit). Run
from the repository root:

    python3 -m tools.bench.stage3_cloud

Installs Rust if needed, builds hypre with CUDA for the GPU found (the
reference, D20), builds the single-source PTX kernels and the CLI with
CUDA and hypre, and runs the stage-3 GPU checks:

1. the pressure solver on about 0.9 million polyhedra: CG with the
   smoothed-aggregation multigrid (mixed precision and float64) and
   Jacobi-CG on the CPU and on the GPU, hypre's BoomerAMG-PCG on the GPU;
   times per cell and iteration, GPU against CPU bitwise;
2. the polyhedral mesh sweep (729 to 274 625 cells) on the GPU with
   hypre: iterations nearly independent of the mesh size;
3. the stage-3 tests on the CPU against the GPU build of hypre.

Takes about 30 minutes on a free Colab T4, most of it compiling hypre.
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.bench import build_hypre, cloud  # noqa: E402

RUNS = [
    (["examples/cases/poisson/poisson_poly_gpu.toml"], "stage3_gpu_poly"),
    (["examples/cases/poisson/poisson_poly_box.toml", "--backends", "cuda",
      "--solvers", "amg,amg-double,hypre"], "stage3_gpu_sweep"),
]


def main() -> int:
    cc = cloud.gpu()
    cloud.rust()
    cloud.cuda_paths(bin_dir=True)
    print(f"== hypre with CUDA for sm_{cc} (about 20 minutes)", flush=True)
    hypre = cloud.ROOT / "target" / f"hypre-cuda{cc}"
    if not (hypre / "include" / "HYPRE.h").is_file():
        build_hypre.main([str(hypre), cc])
    os.environ["HYPRE_DIR"] = str(hypre)
    cloud.kernels(cc)
    cloud.cli("cuda,hypre")
    status = 0
    for args, out in RUNS:
        print(f"== foehn poisson {' '.join(args)} --out output/{out}", flush=True)
        if cloud.sh([str(cloud.FOEHN), "poisson", *args, "--out", f"output/{out}"], check=False):
            status = 1
    print("== Stage-3 tests against the GPU build of hypre", flush=True)
    if cloud.sh(["cargo", "test", "--release", "--locked", "-p", "foehncraft-fv", "--features", "hypre",
                 "--test", "stage3"], check=False):
        status = 1
    cloud.reports([cloud.ROOT / f"output/{out}/poisson.md" for _, out in RUNS])
    return status


if __name__ == "__main__":
    sys.exit(main())

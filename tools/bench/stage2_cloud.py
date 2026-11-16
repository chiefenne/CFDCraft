"""Stage 2 on a cloud GPU (Google Colab, Kaggle, any Linux machine with an
NVIDIA GPU of compute capability 7.0 or newer). Run from the repository
root:

    python3 -m tools.bench.stage2_cloud

Installs Rust if needed, builds the single-source PTX kernels for the GPU
found, builds the CLI with CUDA and runs the stage-2 GPU checks: the
polyhedral manufactured-solution case on about 0.9 million cells (CPU and
twice on the GPU: GPU against CPU and GPU run to run), the oblique step
with all its schemes on CPU and GPU, and the polyhedral order-of-accuracy
sweep on the GPU alone. Takes about 15 minutes on a free Colab T4, most of
it compiling.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.bench import cloud  # noqa: E402

RUNS = [
    (["examples/cases/scalar/mms_poly_gpu.toml"], "stage2_gpu_mms"),
    (["examples/cases/scalar/oblique_step_2d.toml", "--backends", "cpu,cuda"], "stage2_gpu_oblique"),
    (["examples/cases/scalar/mms_poly_box.toml", "--backends", "cuda"], "stage2_gpu_sweep"),
]


def main() -> int:
    cc = cloud.gpu()
    cloud.rust()
    cloud.cuda_paths()
    cloud.kernels(cc)
    cloud.cli("cuda")
    status = 0
    for args, out in RUNS:
        print(f"== foehn scalar {' '.join(args)} --out output/{out}", flush=True)
        if cloud.sh([str(cloud.FOEHN), "scalar", *args, "--out", f"output/{out}"], check=False):
            status = 1
    cloud.reports([cloud.ROOT / f"output/{out}/scalar.md" for _, out in RUNS])
    return status


if __name__ == "__main__":
    sys.exit(main())

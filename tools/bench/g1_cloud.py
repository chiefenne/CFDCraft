"""Gate G1 on a cloud GPU (Google Colab, Kaggle, Lightning.ai, any Linux
machine with an NVIDIA GPU of compute capability 7.0 or newer and the CUDA
12 libraries). Run from the repository root:

    python3 -m tools.bench.g1_cloud [CASE.toml]     # default: examples/cases/bench/g1_gpu_10m.toml

Installs Rust if needed, builds the single-source PTX kernels for the GPU
found, builds the CLI with CUDA, runs the case twice and prints both
reports. Takes about 10 minutes on a free Colab T4, most of it compiling.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.bench import cloud  # noqa: E402


def main(argv) -> int:
    case = argv[0] if argv else "examples/cases/bench/g1_gpu_10m.toml"
    cc = cloud.gpu(kaggle_hint=True)
    cloud.rust()
    cloud.cuda_paths()
    cloud.kernels(cc)
    cloud.cli("cuda")
    for run in (1, 2):
        print(f"== G1 run {run} of 2: {case}", flush=True)
        cloud.sh([str(cloud.FOEHN), "bench", "g1", case, "--out", f"output/g1_cloud_run{run}"])
    print("\n== Reports (copy everything below into the project thread)")
    for run in (1, 2):
        print(f"\n### Run {run}")
        print((cloud.ROOT / f"output/g1_cloud_run{run}/g1.md").read_text(encoding="utf-8"), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

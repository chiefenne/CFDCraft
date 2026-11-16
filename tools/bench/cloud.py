"""Shared steps of the cloud GPU drivers (Google Colab, Kaggle, any Linux
machine with an NVIDIA GPU of compute capability 7.0 or newer): check the
GPU, install Rust if needed, find the CUDA libraries, build the PTX kernels
and the CLI with CUDA."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FOEHN = ROOT / "target" / "release" / "foehn"


def sh(cmd: list[str], check: bool = True, **kw) -> int:
    """Runs a command from the repository root, output on the terminal."""
    r = subprocess.run(cmd, cwd=ROOT, **kw)
    if check and r.returncode != 0:
        raise SystemExit(f"failed ({r.returncode}): {' '.join(map(str, cmd))}")
    return r.returncode


def gpu(kaggle_hint: bool = False) -> str:
    """Prints the GPU and returns its compute capability without the dot
    (`75`); exits when it is older than 7.0."""
    print("== GPU", flush=True)
    q = ["nvidia-smi", "--format=csv,noheader"]
    sh(q + ["--query-gpu=name,compute_cap,memory.total,driver_version"])
    cc = subprocess.run(q + ["--query-gpu=compute_cap"], capture_output=True, text=True, check=True)
    cc = cc.stdout.splitlines()[0].strip().replace(".", "")
    if int(cc) < 70:
        print(f"compute capability {cc[0]}.{cc[1:]} is older than 7.0: the single-source Rust route "
              "cannot run here (O13).", file=sys.stderr)
        if kaggle_hint:
            print("On Kaggle choose the accelerator 'GPU T4 x2', not 'P100'.", file=sys.stderr)
        raise SystemExit(1)
    return cc


def rust() -> None:
    print("== Rust", flush=True)
    cargo_bin = Path.home() / ".cargo" / "bin"
    if shutil.which("cargo") is None:
        if not (cargo_bin / "cargo").exists():
            sh(["sh", "-c", "curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal"])
        os.environ["PATH"] = f"{cargo_bin}{os.pathsep}{os.environ['PATH']}"
    sh(["cargo", "--version"])


def cuda_paths(bin_dir: bool = False) -> None:
    """NVRTC and the CUDA libraries of the image (and nvcc for hypre)."""
    for d in ("/usr/local/cuda/lib64", "/usr/local/cuda/targets/x86_64-linux/lib"):
        if os.path.isdir(d):
            os.environ["LD_LIBRARY_PATH"] = d + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")
    if bin_dir and os.path.isdir("/usr/local/cuda/bin"):
        os.environ["PATH"] = "/usr/local/cuda/bin" + os.pathsep + os.environ["PATH"]


def kernels(cc: str) -> None:
    print(f"== Single-source kernels -> PTX for sm_{cc}", flush=True)
    sh([sys.executable, "gpu/ptx-kernels/build.py", f"sm_{cc}"])


def cli(features: str = "cuda") -> None:
    print("== CLI with " + " and ".join(features.split(",")), flush=True)
    sh(["cargo", "build", "--release", "--locked", "-p", "foehncraft-cli", "--features", features])
    sh([str(FOEHN), "info"])


def reports(paths: list[Path]) -> None:
    print("\n== Reports (copy everything below into the project thread)")
    for p in paths:
        print()
        print(Path(p).read_text(encoding="utf-8"), end="")

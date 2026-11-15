#!/usr/bin/env python3
"""Stage 5 on a cloud GPU (Google Colab, Kaggle, any Linux machine with an
NVIDIA GPU of compute capability 7.0 or newer and its driver). Run from the
repository root:

    python3 tools/stage5_cloud.py [--quick]

Installs Rust if needed, builds the single-source PTX kernels for the GPU
found and the CLI with the `cuda` feature, then runs flow cases twice, with
the flow loop on the CPU (`--backend cpu`) and on the GPU (`--backend
cuda`), and compares the two runs bitwise with tools/compare_runs.py (D31):
fields, residual histories, forces, probes, reports and restart files.

Cases: Couette and Poiseuille flows on mixed and polyhedral meshes, the
lid-driven cavity steady and transient, Taylor-Green with Crank-Nicolson,
the start of the cylinder flow at Re = 100 on 98 304 cells, and the cavity
on about one million polyhedra (30 SIMPLEC iterations; the time per cell
and iteration of CPU and GPU). `--quick` skips the two large cases.

Writes output/stage5_gpu/report.md (also printed) and exits with 1 when a
run fails or a GPU run differs from its CPU run. Takes about 15 minutes on
a free Colab T4.
"""
import json
import os
import pathlib
import shutil
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "output" / "stage5_gpu"
F = "examples/cases/flow/"
CASES = [
    F + "couette_poly.toml",
    F + "poiseuille_mixed.toml",
    F + "cavity_re100.toml",
    F + "cavity_re100_poly3d.toml",
    F + "cavity_re100_transient.toml",
    F + "taylor_green_cn.toml",
]
LARGE = [F + "cylinder_re100_gpu.toml", F + "cavity_poly3d_gpu.toml"]


def sh(cmd, **kw):
    print("$ " + " ".join(cmd), flush=True)
    return subprocess.run(cmd, cwd=ROOT, **kw)


def gpu():
    q = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,compute_cap,memory.total,driver_version", "--format=csv,noheader"],
        capture_output=True, text=True, check=True,
    ).stdout.strip().splitlines()[0]
    name, cc = q.split(",")[0].strip(), q.split(",")[1].strip()
    print("== GPU: " + q, flush=True)
    sm = cc.replace(".", "")
    if int(sm) < 70:
        sys.exit(f"compute capability {cc} is older than 7.0: the single-source Rust route cannot run here (O13)")
    return name, sm


def rust():
    cargo = shutil.which("cargo") or str(pathlib.Path.home() / ".cargo/bin/cargo")
    if not pathlib.Path(cargo).exists():
        sh(["sh", "-c", "curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal"], check=True)
    os.environ["PATH"] = str(pathlib.Path.home() / ".cargo/bin") + os.pathsep + os.environ["PATH"]
    sh(["cargo", "--version"], check=True)


def stats(d):
    j = json.loads((d / "flow.json").read_text())
    s = j["solves"][-1]
    t = s.get("transient") or {}
    per = t.get("ns_per_cell_step") or s.get("ns_per_cell_iteration")
    unit = "step" if t else "iteration"
    return s["cells"], j["seconds"], per, unit


def main():
    quick = "--quick" in sys.argv
    name, sm = gpu()
    rust()
    print(f"== Single-source kernels -> PTX for sm_{sm}", flush=True)
    sh(["gpu/ptx-kernels/build.sh", f"sm_{sm}"], check=True)
    print("== CLI with CUDA", flush=True)
    sh(["cargo", "build", "--release", "--locked", "-p", "foehncraft-cli", "--features", "cuda"], check=True)
    foehn = str(ROOT / "target/release/foehn")
    sh([foehn, "info"], check=True)
    shutil.rmtree(OUT, ignore_errors=True)
    rows, status = [], 0
    for case in CASES + ([] if quick else LARGE):
        n = pathlib.Path(case).stem
        times = {}
        ok = True
        for be in ("cpu", "cuda"):
            d = OUT / n / be
            print(f"== {n} on {be}", flush=True)
            t0 = time.time()
            r = sh([foehn, "flow", case, "--backend", be, "--out", str(d)])
            times[be] = time.time() - t0
            if r.returncode != 0:
                print(f"!! {n} on {be} failed (exit code {r.returncode})", flush=True)
                ok = False
        if not ok:
            rows.append(f"| {n} | | | | | | failed |")
            status = 1
            continue
        c = sh([sys.executable, "tools/compare_runs.py", str(OUT / n / "cpu"), str(OUT / n / "cuda")])
        equal = c.returncode == 0
        status |= 0 if equal else 1
        cells, _, per_cpu, unit = stats(OUT / n / "cpu")
        _, _, per_gpu, _ = stats(OUT / n / "cuda")
        rows.append(
            f"| {n} | {cells} | {times['cpu']:.1f} | {times['cuda']:.1f} | "
            f"{per_cpu:.0f} | {per_gpu:.0f} ns/cell/{unit} | {'yes' if equal else '**no**'} |"
        )
    head = (
        f"# Stage 5: flow loop on the GPU\n\nGPU: {name} (sm_{sm}); CPU: {os.cpu_count()} threads. "
        "Wall times include mesh generation and output.\n\n"
        "| case | cells | CPU s | GPU s | CPU ns | GPU ns | bitwise equal |\n|---|---:|---:|---:|---:|---:|---|\n"
    )
    rep = head + "\n".join(rows) + "\n"
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.md").write_text(rep)
    print("\n" + rep)
    print("all GPU runs bitwise equal to the CPU runs" if status == 0 else "FAILED: see above")
    return status


if __name__ == "__main__":
    sys.exit(main())

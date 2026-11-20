#!/usr/bin/env python3
"""Stage 6 on a cloud GPU (Google Colab, Kaggle, any Linux machine with an
NVIDIA GPU of compute capability 7.0 or newer and its driver): the parity
check of specification 6.6 for the turbulence models. Run from the
repository root:

    python3 -m tools.bench.stage6_cloud [--quick] [--iterations N]
    python3 -m tools.bench.stage6_cloud --cpu-only     # no GPU: 1 against 4 threads

Installs Rust if needed, builds the single-source PTX kernels for the GPU
found and the CLI with the `cuda` feature, then runs every turbulence model
on the TMR flat plate and the channel at Re_tau 395 twice, with the flow
loop on the CPU (`--backend cpu`) and on the GPU (`--backend cuda`), and
compares the two runs bitwise with tools/check/compare_runs.py (D31):
fields, residual histories, forces, wall values, reports and restart files
(the turbulence fields, their old levels and the multigrid matrices).
Further cases: the realizable k-epsilon and SST on the plate, the heated
channel and pipe (passive scalar with the thermal wall function), the
coupled algorithm on the backward-facing step. Each run stops after N outer
iterations (default 300; the check is on the bits, not on convergence).

`--quick` runs the channels only. `--cpu-only` builds without CUDA and
compares one thread against four on the CPU (the thread parity of 6.6,
also a test of this script without a GPU).

Writes output/stage6_gpu/report.md (also printed) and exits with 1 when a
run fails or a pair of runs differs. Takes about 20 minutes on a free
Colab T4 (the build about 8).
"""
import json
import os
import pathlib
import shutil
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "stage6_gpu"
F = "examples/cases/flow/"
# (case file, extra --set values); the channels first (--quick)
CHANNELS = [
    (F + "channel_sst.toml", []),
    (F + "channel_keps.toml", []),
    (F + "channel_zeta_f.toml", []),
    (F + "channel_zeta_f0.toml", []),
    (F + "channel_eb.toml", []),
]
OTHERS = [
    (F + "tmr_plate_sst.toml", []),
    (F + "tmr_plate_keps.toml", []),
    (F + "tmr_plate_keps.toml", ['turbulence.model="realizable-k-epsilon"']),
    (F + "tmr_plate_zeta_f.toml", []),
    (F + "tmr_plate_zeta_f0.toml", []),
    (F + "tmr_plate_eb.toml", []),
    (F + "channel_heated_keps.toml", []),
    (F + "pipe_heated_keps.toml", []),
    (F + "backstep_zeta_f0_coupled.toml", []),
    # the SST options of D170: curvature correction in a turning frame, QCR
    (F + "channel_rotating_sst.toml", []),
    (F + "square_duct_qcr.toml", []),
]


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
    return s["cells"], s.get("ns_per_cell_iteration")


def label(case, sets):
    n = pathlib.Path(case).stem
    for s in sets:
        n += "_" + s.split("=")[1].strip('"').replace("-k-epsilon", "")
    return n


def main():
    quick = "--quick" in sys.argv
    cpu_only = "--cpu-only" in sys.argv
    iterations = 300
    if "--iterations" in sys.argv:
        iterations = int(sys.argv[sys.argv.index("--iterations") + 1])
    rust()
    if cpu_only:
        name, sm = "none", None
        pair = (("cpu", ["--threads", "1"]), ("cpu4", ["--threads", "4"]))
        sh(["cargo", "build", "--release", "--locked", "-p", "foehncraft-cli"], check=True)
    else:
        name, sm = gpu()
        pair = (("cpu", []), ("cuda", []))
        print(f"== Single-source kernels -> PTX for sm_{sm}", flush=True)
        sh([sys.executable, "gpu/ptx-kernels/build.py", f"sm_{sm}"], check=True)
        print("== CLI with CUDA", flush=True)
        sh(["cargo", "build", "--release", "--locked", "-p", "foehncraft-cli", "--features", "cuda"], check=True)
    foehn = str(ROOT / "target/release/foehn")
    sh([foehn, "info"], check=True)
    shutil.rmtree(OUT, ignore_errors=True)
    common = [
        f"solver.max_iterations={iterations}",
        "verification.require_convergence=false",
        "output.vtkhdf=true",
        "run.restart_final=true",
    ]
    rows, status = [], 0
    for case, sets in CHANNELS + ([] if quick else OTHERS):
        n = label(case, sets)
        times, ok = {}, True
        for tag, extra in pair:
            be = "cuda" if tag == "cuda" else "cpu"
            d = OUT / n / tag
            print(f"== {n} on {tag}", flush=True)
            cmd = [foehn, "flow", case, "--backend", be, "--out", str(d)]
            for s in common + sets:
                cmd += ["--set", s]
            t0 = time.time()
            r = sh(cmd + extra)
            times[tag] = time.time() - t0
            if r.returncode != 0:
                print(f"!! {n} on {tag} failed (exit code {r.returncode})", flush=True)
                ok = False
        a, b = pair[0][0], pair[1][0]
        if not ok:
            rows.append(f"| {n} | | | | | | failed |")
            status = 1
            continue
        c = sh([sys.executable, "tools/check/compare_runs.py", str(OUT / n / a), str(OUT / n / b)])
        equal = c.returncode == 0
        status |= 0 if equal else 1
        cells, per_a = stats(OUT / n / a)
        _, per_b = stats(OUT / n / b)
        fmt = lambda v: f"{v:.0f}" if v else ""  # noqa: E731
        rows.append(
            f"| {n} | {cells} | {times[a]:.1f} | {times[b]:.1f} | {fmt(per_a)} | {fmt(per_b)} | "
            f"{'yes' if equal else '**no**'} |"
        )
    what = "one thread against four on the CPU" if cpu_only else f"CPU against GPU {name} (sm_{sm})"
    head = (
        f"# Stage 6: turbulence models, {what}\n\n"
        f"CPU: {os.cpu_count()} threads; {iterations} outer iterations per run. "
        "Wall times include mesh generation and output; ns per cell and iteration.\n\n"
        f"| case | cells | {a} s | {b} s | {a} ns | {b} ns | bitwise equal |\n"
        "|---|---:|---:|---:|---:|---:|---|\n"
    )
    rep = head + "\n".join(rows) + "\n"
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.md").write_text(rep)
    print("\n" + rep)
    print("all pairs bitwise equal" if status == 0 else "FAILED: see above")
    return status


if __name__ == "__main__":
    sys.exit(main())

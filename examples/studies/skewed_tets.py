"""The skewed-tetrahedra test tank of the cell-quality remediation (D50):
the same steady MRF run with the bounded face velocity of the mass flux
(the remediation as it is) and with the interpolated one (the rest of the
remediation kept), 300 outer iterations each, and the comparison.

    pip install gmsh
    cargo build --release -p foehncraft-cli
    python3 examples/studies/skewed_tets.py [OUT] [--iterations 300] [--threads N]

Builds examples/studies/skewed_tets/stirred_tank_raw.msh with its mesh
script when it is missing (about 36k tetrahedra, faces with a skewness up
to 1.3 in the blade-hub and sparger corners), runs
examples/studies/skewed_tets/case.toml into OUT/bounded and
OUT/interpolated (default OUT = output/skewed_tets; a few minutes each on
four cores) and prints the largest residual and the impeller torque every
50 iterations. With the interpolated face velocity the torque runs away
(from about -1.1 N m to -100 N m and more by iteration 300); with the
bounded one the residuals keep falling.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "examples" / "studies" / "skewed_tets"
KEEP = [
    "--set", "solver.remediation=off",
    "--set", "solver.pressure_gradient=gauss",
    "--set", "solver.momentum.limiter=venkatakrishnan",
]


def history(out: Path):
    """(iteration, largest residual, torque) of the logged iterations of
    run.jsonl (every `solver.log_every`)."""
    rows = []
    with open(out / "run.jsonl") as f:
        for line in f:
            r = json.loads(line)
            if r.get("type") != "table_row" or r.get("table") != "iteration":
                continue
            v = r["values"]
            res = [x for k, x in v.items() if k.startswith("residual.") and x is not None]
            tq = [x for k, x in v.items() if "torque" in k and x is not None]
            rows.append((r["index"], max(res), tq[0] if tq else float("nan")))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", nargs="?", default="output/skewed_tets")
    ap.add_argument("--iterations", type=int, default=300)
    ap.add_argument("--threads", type=int, default=None)
    a = ap.parse_args()
    foehn = ROOT / "target" / "release" / ("foehn.exe" if os.name == "nt" else "foehn")
    if not foehn.is_file():
        raise SystemExit(f"{foehn} is missing: cargo build --release -p foehncraft-cli")
    if not (HERE / "stirred_tank_raw.msh").is_file():
        subprocess.run([sys.executable, str(HERE / "mesh.py")], check=True)
    out = Path(a.out)
    runs = {"bounded": [], "interpolated": KEEP + ["--set", "solver.face_velocity=interpolated"]}
    for name, extra in runs.items():
        cmd = [str(foehn), "flow", str(HERE / "case.toml"), "--out", str(out / name),
               "--set", f"solver.max_iterations={a.iterations}", "--set", "solver.log_every=50"] + extra
        if a.threads:
            cmd += ["--threads", str(a.threads)]
        print(" ".join(cmd), flush=True)
        # a divergence ends the run with an error; the history is still written
        subprocess.run(cmd, cwd=ROOT, stdout=subprocess.DEVNULL)
    print(f"\n{'iteration':>9} | {'bounded: residual':>17} {'torque':>10} | {'interpolated: residual':>22} {'torque':>10}")
    hb, hi = (dict((i, (r, t)) for i, r, t in history(out / n)) for n in runs)
    for it in range(50, a.iterations + 1, 50):
        b, i = hb.get(it), hi.get(it)
        fb = f"{b[0]:17.2e} {b[1]:10.3g}" if b else f"{'-':>17} {'-':>10}"
        fi = f"{i[0]:22.2e} {i[1]:10.3g}" if i else f"{'diverged':>22} {'-':>10}"
        print(f"{it:9d} | {fb} | {fi}")


if __name__ == "__main__":
    main()

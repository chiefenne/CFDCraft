"""The torque histories of the rotors of transient runs (gallery 101, 102).

    pip install matplotlib numpy
    python3 -m tools.plot.rotor_torque OUT.png RUN[:LABEL][@PANEL] [...] [options]

RUN is the output directory of a transient `foehn flow` run with a turning
zone (`[[mrf]]` with `[time]`, `[[sliding]]`); its `forces.csv` has a column
`<rotor>_torque` per rotor (N m, about the axis along the rotor's own Ω, so
a rotor that drives the fluid has a negative torque). Every rotor of every
run is one line, labelled `LABEL rotor` (LABEL defaults to the directory
name). With `--panels` each run gets its own panel, or the panel PANEL
(0, 1, ...) when given, so that two runs can share one.

Options: --rotors rotor_left,rotor_right (only these), --tmin (s, leaves out
the impulsive start of a run from rest), --panels,
--titles (panel titles, comma separated), --dpi, --width (inches).
"""
import argparse
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read(run):
    """Time and the torque columns of the forces.csv of `run`."""
    with open(os.path.join(run, "forces.csv")) as f:
        rows = list(csv.reader(line for line in f if not line.startswith("#")))
    head, data = rows[0], np.array(rows[1:], dtype=float)
    t = data[:, head.index("time")]
    walls = {h[: -len("_torque")] for h in head if h.endswith("_torque")}
    # the rotors: a torque column with a power column next to it
    rotors = [w for w in walls if f"{w}_power" in head]
    return t, {r: data[:, head.index(f"{r}_torque")] for r in sorted(rotors)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--rotors", default=None)
    ap.add_argument("--tmin", type=float, default=0.0)
    ap.add_argument("--panels", action="store_true")
    ap.add_argument("--titles", default=None)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--width", type=float, default=None)
    a = ap.parse_args()
    keep = a.rotors.split(",") if a.rotors else None
    runs = []
    for k, spec in enumerate(a.runs):
        spec, _, panel = spec.partition("@")
        path, _, label = spec.partition(":")
        runs.append((path, label or os.path.basename(os.path.normpath(path)),
                     int(panel) if panel else k))
    n = max(r[2] for r in runs) + 1 if a.panels else 1
    width = a.width or 5.5 * n
    fig, axes = plt.subplots(1, n, figsize=(width, 3.6), squeeze=False)
    titles = a.titles.split(",") if a.titles else []
    shared = {p for p in (r[2] for r in runs) if sum(q[2] == p for q in runs) > 1}
    for path, label, k in runs:
        ax = axes[0][k if a.panels else 0]
        t, torques = read(path)
        m = t >= a.tmin
        for r, q in torques.items():
            if keep and r not in keep:
                continue
            ax.plot(t[m], q[m] * 1e3, lw=1.2, label=f"{label} {r}" if not a.panels or k in shared else r)
        if a.panels and k < len(titles):
            ax.set_title(titles[k], fontsize=10)
    for ax in axes[0]:
        ax.axhline(0.0, color="0.5", lw=0.6)
        ax.set_xlabel("t (s)")
        ax.set_ylabel("torque of the fluid on the rotor (mN m)")
        ax.legend(fontsize=8, frameon=False)
        ax.grid(alpha=0.3)
    if not a.panels and titles:
        axes[0][0].set_title(titles[0], fontsize=10)
    fig.tight_layout()
    fig.savefig(a.out, dpi=a.dpi)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()

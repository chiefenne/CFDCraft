"""The channel in spanwise rotation (`channel_rotating_sst.toml`, Ro = 0.15,
gallery 113): the velocity and the turbulence kinetic energy across the
channel with and without the curvature correction of SST, from the
`profile` line of two runs.

    pip install matplotlib numpy
    python3 -m tools.plot.rotating_channel_chart OUT.svg ON_DIR OFF_DIR

e.g.

    target/release/foehn flow examples/cases/flow/channel_rotating_sst.toml --out output/rot_on
    target/release/foehn flow examples/cases/flow/channel_rotating_sst.toml --out output/rot_off \\
        --set turbulence.curvature_correction=false
    python3 -m tools.plot.rotating_channel_chart docs/figures/gallery/rotating_channel.svg \\
        output/rot_on output/rot_off
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

U_B = 14.95


def profile(d):
    with open(os.path.join(d, "line_profile.csv")) as f:
        rows = [ln for ln in f if not ln.startswith("#")]
    head = rows[0].strip().split(",")
    a = np.array([[float(v) for v in ln.split(",")] for ln in rows[1:]])
    col = {n: a[:, i] for i, n in enumerate(head)}
    return col["y"], col["u"], col["k"]


def walls(d):
    s = json.load(open(os.path.join(d, "flow.json")))["solves"][-1]
    f = {b["name"]: b["force"][0] for b in s["boundaries"]}
    return f["ymin"], f["ymax"]


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    out, on, off = sys.argv[1:]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    for d, label, colour, ls in ((off, "SST", "0.35", "--"), (on, "SST with the correction", "#c0392b", "-")):
        y, u, k = profile(d)
        lo, hi = walls(d)
        a1.plot(u / U_B, y, ls, color=colour, label=f"{label}: τ_w lower/upper {lo / hi:.2f}")
        a2.plot(k / U_B**2 * 1e3, y, ls, color=colour, label=label)
    for a in (a1, a2):
        a.axhline(1.0, color="0.85", lw=0.8, zorder=0)
        a.set_ylim(0, 2)
        a.grid(alpha=0.3)
    a1.set_xlabel("u / U_b")
    a1.set_ylabel("y / h   (lower wall: unstable side)")
    a2.set_xlabel("k / U_b² · 10³")
    a1.legend(fontsize=8, loc="center left")
    a2.legend(fontsize=8, loc="center right")
    fig.suptitle("Channel in spanwise rotation, Re_b 5 800, Ro = 2Ωh/U_b = 0.15", fontsize=11)
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()

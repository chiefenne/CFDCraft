"""Passive scalars (D51, D116) and the bounded BDF2 (D28 (1)): the
temperature profile of the laminar heated channel, its Nusselt number on
three meshes, and the head of a bounded front with and without the
bounds. Prints the numbers and writes the gallery chart.

    python3 -m tools.plot.scalar_chart docs/figures/gallery/passive_scalars.svg SWEEP MID FRONT FREE

SWEEP is the output directory of channel_heated.toml (three meshes), MID
the same case on its middle mesh with a line `mid` at x = 1 m through the
cell centres, FRONT and FREE the output directories of
scalar_front_bdf2.toml with and without the bounds:

    target/release/foehn flow examples/cases/flow/channel_heated.toml --out SWEEP
    target/release/foehn flow examples/cases/flow/channel_heated.toml --out MID \\
        --set 'verification.levels=[1]' \\
        --set 'output.lines=[{name="mid",from=[1.0,0.015625,0.05],to=[1.0,0.984375,0.05],points=32}]'
    target/release/foehn flow examples/cases/flow/scalar_front_bdf2.toml --out FRONT
    target/release/foehn flow examples/cases/flow/scalar_front_bdf2.toml --out FREE \\
        --set 'scalars.c.bounds=[-1e9, 1e9]'
"""
import csv
import json
import os
import sys

from tools.plot.property_chart import panel

Q, GAMMA, U, H, L = 1e-3, 1e-3, 1.0, 1.0, 2.0
G = 2 * Q / (U * H)
NU = 140 / 17


def theta(x, y):
    e = y / H
    return G * x + H * (2 * e**3 - e**4 - e)


def line(d, name, col):
    with open(os.path.join(d, f"line_{name}.csv")) as f:
        rows = [r for r in f if not r.startswith("#")]
    return [(float(r["x"]), float(r["y"]), float(r[col])) for r in csv.DictReader(rows)]


def nusselt(d):
    with open(os.path.join(d, "result.json")) as f:
        r = json.load(f)
    out = []
    for s in r["solves"]:
        sc = s["scalars"][0]
        b = {x["boundary"]: x for x in sc["boundaries"]}
        tw = 0.5 * (b["ymin"]["mean"] + b["ymax"]["mean"])
        tb = b["xmin"]["mixed_mean"] + G * L / 2
        out.append((s["cells"], 2 * Q / (GAMMA * (tw - tb)), sc["l2"]))
    return out


def main():
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    dst, sweep, mid, front, free = sys.argv[1:6]
    nu = nusselt(sweep)
    for cells, n, l2 in nu:
        print(f"{cells:6d} cells: Nu {n:.4f} ({100 * (n / NU - 1):+.3f} %), L2 error of theta {l2:.3e} K")
    prof = line(mid, "mid", "theta")
    e = max(abs(t - theta(x, y)) for x, y, t in prof)
    print(f"largest difference on the line x = 1 m from the exact profile: {e:.3e} K")
    fb = [(x, 1e4 * (c - 1)) for x, _, c in line(front, "axis", "c") if x < 0.8]
    ff = [(x, 1e4 * (c - 1)) for x, _, c in line(free, "axis", "c") if x < 0.8]
    print(f"largest c - 1 with bounds {max(v for _, v in fb) * 1e-4:.3e}, without {max(v for _, v in ff) * 1e-4:.3e}")
    # the panel shows the head: the front itself falls out of it (clipped)
    def head(pts):
        n = next((i for i, (_, v) in enumerate(pts) if v < -5.0), len(pts) - 1)
        return [(x, max(v, -5.0)) for x, v in pts[: n + 1]]

    fb, ff = head(fb), head(ff)
    w, h = 1020, 380
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
           'font-family="Helvetica, Arial, sans-serif" font-size="12">',
           f'<rect width="{w}" height="{h}" fill="white"/>']
    ys = [i / 100 for i in range(101)]
    panel(out, 70, 50, 250, 260, "heated channel, x = 1 m", "θ (K)", "y (m)", (-0.32, 0.02), (0, 1),
          [("line", [(theta(1.0, y), y) for y in ys], "#1f6fb4", "exact"),
           ("dot", [(t, y) for _, y, t in prof], "#d9822b", "64 x 32 cells")],
          ([-0.3, -0.2, -0.1, 0], [0, 0.25, 0.5, 0.75, 1]), legend=12)
    pts = [(i, n) for i, (_, n, _) in enumerate(nu)]
    panel(out, 410, 50, 250, 260, "Nusselt number", "mesh (cells)", "Nu", (-0.4, len(nu) - 0.6), (8.2, 8.34),
          [("dash", [(-0.4, NU), (len(nu) - 0.6, NU)], "#1f6fb4", "exact 140/17"),
           ("line", pts, "#d9822b", "Foehncraft"),
           ("dot", pts, "#d9822b", "")],
          ([], [8.2, 8.24, 8.28, 8.32]), legend=110)
    sx = lambda v: 410 + (v + 0.4) / (len(nu) - 0.2) * 250
    for i, (cells, _, _) in enumerate(nu):
        out.append(f'<text x="{sx(i):.1f}" y="329" text-anchor="middle">{cells}</text>')
    panel(out, 750, 50, 250, 260, "head of a bounded front", "x (m)", "(c − 1) x 10⁴", (0, 0.8), (-5, 4),
          [("line", ff, "#c0392b", "BDF2"),
           ("line", fb, "#1f6fb4", "bounded BDF2"),
           ("dash", [(0, 0), (0.8, 0)], "#888888", "upper bound")],
          ([0, 0.2, 0.4, 0.6, 0.8], [-4, -2, 0, 2, 4]), legend=12)
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()

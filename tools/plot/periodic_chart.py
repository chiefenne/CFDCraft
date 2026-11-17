"""Periodic boundaries (D29, D109 (1), D111): the temperature across the
short periodic heated pipe against the exact profile, the Nusselt numbers
of the periodic heated channel and pipe on three meshes, and the driving
force that holds the mass flow of the periodic Poiseuille flow and of the
pipe. Prints the numbers and writes the gallery chart.

    python3 -m tools.plot.periodic_chart docs/figures/gallery/periodic.svg POIS CHAN PIPE MID

POIS, CHAN and PIPE are the output directories of poiseuille_periodic.toml,
channel_heated_periodic.toml and pipe_heated_periodic.toml (three meshes
each), MID the pipe on its middle mesh with a line `radius` at x = 0.5 m:

    target/release/foehn flow examples/cases/flow/poiseuille_periodic.toml --out POIS
    target/release/foehn flow examples/cases/flow/channel_heated_periodic.toml --out CHAN
    target/release/foehn flow examples/cases/flow/pipe_heated_periodic.toml --out PIPE
    target/release/foehn flow examples/cases/flow/pipe_heated_periodic.toml --out MID \\
        --set 'verification.levels=[1]' \\
        --set 'output.lines=[{name="radius",from=[0.5,0.0078125,0.0],to=[0.5,0.4921875,0.0],points=32}]'
"""
import csv
import json
import os
import sys

from tools.plot.property_chart import panel

Q, GAMMA, R = 1e-3, 1e-3, 0.5
NU_CHANNEL, NU_PIPE = 140 / 17, 48 / 11
BETA = {"pois": 1.2, "pipe": 3.2}


def solves(d):
    with open(os.path.join(d, "result.json")) as f:
        return json.load(f)["solves"]


def sides(s):
    return {x["boundary"]: x for x in s["scalars"][0]["boundaries"]}


def nu_channel(s):
    b = sides(s)
    # mean wall temperature at x = 0: G L / 2 below the area mean
    tw = 0.5 * (b["ymin"]["mean"] + b["ymax"]["mean"]) - 2e-3 * 2.0 / 2
    return 2 * Q / (GAMMA * (tw - b["xmin"]["mixed_mean"]))


def nu_pipe(s):
    b = sides(s)
    tb = 0.5 * (b["xmin"]["mixed_mean"] + b["xmax"]["mixed_mean"])
    return Q * 2 * R / (GAMMA * (b["ymax"]["mean"] - tb))


def theta_pipe(r):
    return Q / (GAMMA * R) * (r * r - r**4 / (4 * R * R)) + 2 * Q / R * 0.5


def line(d, name, col):
    with open(os.path.join(d, f"line_{name}.csv")) as f:
        rows = [r for r in f if not r.startswith("#")]
    return [(float(r["y"]), float(r[col])) for r in csv.DictReader(rows)]


def main():
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    dst, pois, chan, pipe, mid = sys.argv[1:6]
    sp, sc, sq = solves(pois), solves(chan), solves(pipe)
    for name, ss, f, exact in (("channel", sc, nu_channel, NU_CHANNEL), ("pipe", sq, nu_pipe, NU_PIPE)):
        for s in ss:
            n = f(s)
            print(f"{name} {s['cells']:6d} cells: Nu {n:.4f} ({100 * (n / exact - 1):+.3f} %), "
                  f"L2 error of theta {s['scalars'][0]['l2']:.3e} K")
    for name, ss, b in (("Poiseuille", sp, BETA["pois"]), ("pipe", sq, BETA["pipe"])):
        for s in ss:
            p = s["periodic"][0]
            print(f"{name} {s['cells']:6d} cells: driving force {p['driving_force']:.6f} N/m³ "
                  f"({100 * (p['driving_force'] / b - 1):+.3f} %), mass flow {p['mass_flow']:.7e} kg/s")
    prof = line(mid, "radius", "theta")
    e = max(abs(t - theta_pipe(r)) for r, t in prof)
    print(f"largest difference on the line x = 0.5 m (the values of the cells the points lie in) from the exact profile: {e:.3e} K")
    w, h = 1020, 380
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
           'font-family="Helvetica, Arial, sans-serif" font-size="12">',
           f'<rect width="{w}" height="{h}" fill="white"/>']
    rs = [i / 200 for i in range(101)]
    panel(out, 70, 50, 250, 260, "periodic pipe, x = 0.5 m", "θ (K)", "r (m)", (0.0, 0.4), (0, 0.5),
          [("line", [(theta_pipe(r), r) for r in rs], "#1f6fb4", "exact"),
           ("dot", [(t, r) for r, t in prof], "#d9822b", "16 x 16 cells")],
          ([0, 0.1, 0.2, 0.3, 0.4], [0, 0.1, 0.2, 0.3, 0.4, 0.5]), legend=12)
    n = len(sc)
    ch = [(i, 100 * (nu_channel(s) / NU_CHANNEL - 1)) for i, s in enumerate(sc)]
    pp = [(i, 100 * (nu_pipe(s) / NU_PIPE - 1)) for i, s in enumerate(sq)]
    panel(out, 410, 50, 250, 260, "Nusselt number", "mesh", "Nu / Nu_exact − 1 (%)", (-0.4, n - 0.6), (-0.5, 3.5),
          [("dash", [(-0.4, 0), (n - 0.6, 0)], "#888888", "exact"),
           ("line", ch, "#d9822b", "channel, 140/17"),
           ("line", pp, "#1f6fb4", "pipe, 48/11"),
           ("dot", ch, "#d9822b", ""),
           ("dot", pp, "#1f6fb4", "")],
          ([], [0, 1, 2, 3]), legend=110)
    sx = lambda v: 410 + (v + 0.4) / (n - 0.2) * 250
    for i in range(n):
        out.append(f'<text x="{sx(i):.1f}" y="329" text-anchor="middle">{i}</text>')
    bp = [(i, 100 * (s["periodic"][0]["driving_force"] / BETA["pois"] - 1)) for i, s in enumerate(sp)]
    bq = [(i, 100 * (s["periodic"][0]["driving_force"] / BETA["pipe"] - 1)) for i, s in enumerate(sq)]
    panel(out, 750, 50, 250, 260, "driving force at fixed mass flow", "mesh", "β / β_exact − 1 (%)",
          (-0.4, n - 0.6), (-1.2, 0.8),
          [("dash", [(-0.4, 0), (n - 0.6, 0)], "#888888", "exact"),
           ("line", bp, "#d9822b", "channel, 12 μU/h²"),
           ("line", bq, "#1f6fb4", "pipe, 8 μU/R²"),
           ("dot", bp, "#d9822b", ""),
           ("dot", bq, "#1f6fb4", "")],
          ([], [-1.2, -0.8, -0.4, 0, 0.4, 0.8]), legend=110)
    sx = lambda v: 750 + (v + 0.4) / (n - 0.2) * 250
    for i in range(n):
        out.append(f'<text x="{sx(i):.1f}" y="329" text-anchor="middle">{i}</text>')
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()

"""The kinetic-energy form of the second-law check (D90 part 4) over the
mesh size: the numerical dissipation of the steady lid-driven cavity at
Re 100 (relative to the lid power) and of the Taylor–Green vortex
(absolute, over 0.4 s at a fixed Courant number, with its pressure and
scheme parts), each on several meshes. Prints one line per run and writes
the gallery chart.

    python3 tools/energy_budget_chart.py docs/figures/gallery/energy_budget.svg [FOEHN] [OUTDIR]

FOEHN defaults to target/release/foehn, OUTDIR to output/energy_budget.
"""
import json
import math
import os
import subprocess
import sys

CAVITY = [16, 32, 64, 128]
TAYLOR_GREEN = [16, 32, 64]


def run(foehn, case, out, sets):
    args = [foehn, "--quiet", "flow", f"examples/cases/flow/{case}.toml", "--out", out]
    for s in sets:
        args += ["--set", s]
    subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(os.path.join(out, "result.json")) as f:
        return json.load(f)


def ledger(res):
    return res["solves"][-1]["balance"]


def chart(cav, tg, tgp, tgs, dst):
    w, h, left, right, top, bottom = 640, 440, 80, 190, 50, 60
    xs = [1.0 / n for n, _ in cav + tg]
    ys = [v for _, v in cav + tg + tgp + tgs]
    x0, x1 = math.floor(math.log10(min(xs))), math.ceil(math.log10(max(xs)))
    y0, y1 = math.floor(math.log10(min(ys))), math.ceil(math.log10(max(ys)))
    sx = lambda v: left + (math.log10(v) - x0) / (x1 - x0) * (w - left - right)
    sy = lambda v: top + (1 - (math.log10(v) - y0) / (y1 - y0)) * (h - top - bottom)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
           'font-family="Helvetica, Arial, sans-serif" font-size="12">',
           f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
           f'<text x="{(left + w - right) / 2}" y="26" text-anchor="middle" font-size="15">'
           'Numerical dissipation of kinetic energy over the mesh size</text>',
           f'<rect x="{left}" y="{top}" width="{w - left - right}" height="{h - top - bottom}" fill="none" stroke="#444"/>']
    for k in range(x0, x1 + 1):
        out.append(f'<line x1="{sx(10.0 ** k):.1f}" y1="{top}" x2="{sx(10.0 ** k):.1f}" y2="{h - bottom}" stroke="#e5e5e5"/>')
        out.append(f'<text x="{sx(10.0 ** k):.1f}" y="{h - bottom + 18}" text-anchor="middle">1e{k}</text>')
    for k in range(y0, y1 + 1):
        out.append(f'<line x1="{left}" y1="{sy(10.0 ** k):.1f}" x2="{w - right}" y2="{sy(10.0 ** k):.1f}" stroke="#e5e5e5"/>')
        out.append(f'<text x="{left - 8}" y="{sy(10.0 ** k) + 4:.1f}" text-anchor="end">1e{k}</text>')
    for pts, colour, dash in ((cav, "#1f6fb4", ""), (tg, "#d9822b", ""),
                              (tgp, "#d9822b", ' stroke-dasharray="6 3"'),
                              (tgs, "#d9822b", ' stroke-dasharray="1 3"')):
        line = " ".join(f"{sx(1.0 / n):.1f},{sy(v):.1f}" for n, v in pts)
        out.append(f'<polyline points="{line}" fill="none" stroke="{colour}" stroke-width="1.8"{dash}/>')
        for n, v in pts:
            r = 4 if not dash else 2.5
            out.append(f'<circle cx="{sx(1.0 / n):.1f}" cy="{sy(v):.1f}" r="{r}" fill="{colour}"/>')
    # slope 2 through the coarsest Taylor–Green point
    n0, v0 = tg[0]
    n1 = tg[-1][0]
    out.append(f'<line x1="{sx(1.0 / n0):.1f}" y1="{sy(v0):.1f}" x2="{sx(1.0 / n1):.1f}" '
               f'y2="{sy(v0 * (n0 / n1) ** 2):.1f}" stroke="#888" stroke-dasharray="4 4"/>')
    out.append(f'<text x="{(left + w - right) / 2}" y="{h - 16}" text-anchor="middle">mesh size h (m)</text>')
    out.append(f'<text x="18" y="{(top + h - bottom) / 2}" text-anchor="middle" '
               f'transform="rotate(-90 18 {(top + h - bottom) / 2})">numerical dissipation</text>')
    legend = (("#1f6fb4", "", "cavity Re 100:", "share of lid power (-)"),
              ("#d9822b", "", "Taylor–Green:", "over the run (J)"),
              ("#d9822b", "6 3", "  pressure part", ""),
              ("#d9822b", "1 3", "  scheme part", ""),
              ("#888", "4 4", "slope 2", ""))
    yy = top + 20
    for c, dash, a, b in legend:
        if dash:
            out.append(f'<line x1="{w - right + 10}" y1="{yy - 4}" x2="{w - right + 30}" y2="{yy - 4}" '
                       f'stroke="{c}" stroke-width="1.8" stroke-dasharray="{dash}"/>')
        else:
            out.append(f'<circle cx="{w - right + 20}" cy="{yy - 4}" r="4" fill="{c}"/>')
        out.append(f'<text x="{w - right + 32}" y="{yy}">{a}</text>')
        if b:
            out.append(f'<text x="{w - right + 32}" y="{yy + 15}">{b}</text>')
        yy += 36 if b else 22
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    dst = sys.argv[1]
    foehn = sys.argv[2] if len(sys.argv) > 2 else "target/release/foehn"
    outdir = sys.argv[3] if len(sys.argv) > 3 else "output/energy_budget"
    cav, tg, tgp, tgs = [], [], [], []
    for n in CAVITY:
        r = run(foehn, "cavity_re100_energy", os.path.join(outdir, f"cavity_{n}"),
                [f"mesh.generate.cells=[{n},{n},1]", "solver.max_iterations=5000",
                 "verification.reference_tolerance=1"])
        k = ledger(r)["last"]["kinetic_energy"]
        cav.append((n, k["relative"]))
        print(f"cavity {n:4d}^2  lid power {k['boundary']:.4e} W  viscous {k['viscous_dissipation']:.4e} W  "
              f"numerical {k['numerical_dissipation']:.4e} W ({100 * k['relative']:.2f} %): "
              f"pressure {k['pressure_dissipation']:+.4e} W, schemes {k['scheme_dissipation']:+.4e} W")
    for n in TAYLOR_GREEN:
        r = run(foehn, "taylor_green_bdf2", os.path.join(outdir, f"tg_{n}"),
                [f"mesh.generate.cells=[{n},{n},1]", "verification.time_steps=[]",
                 f"time.time_step={0.32 / n}", "time.end_time=0.4"])
        b = ledger(r)
        k = b["run"]["kinetic_energy"]
        tg.append((n, k["numerical_dissipation"]))
        tgp.append((n, k["pressure_dissipation"]))
        tgs.append((n, k["scheme_dissipation"]))
        print(f"taylor-green {n:4d}^2  ΔK {k['storage']:.4e} J  viscous {k['viscous_dissipation']:.4e} J  "
              f"numerical {k['numerical_dissipation']:.4e} J ({100 * k['relative']:.2f} %): "
              f"pressure {k['pressure_dissipation']:+.4e} J, schemes {k['scheme_dissipation']:+.4e} J; "
              f"{b['negative_steps']} of {b['steps']} steps negative")
    chart(cav, tg, tgp, tgs, dst)

    def orders(pts):
        return [round(math.log(pts[i][1] / pts[i + 1][1]) / math.log(2.0), 2) for i in range(len(pts) - 1)]
    print(f"Taylor–Green orders: total {orders(tg)}, pressure {orders(tgp)}, schemes {orders(tgs)}; chart {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

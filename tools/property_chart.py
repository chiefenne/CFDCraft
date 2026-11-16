"""Expressions and materials (D53, D102): the velocity profiles of the two
Couette cases with a non-trivial exact solution, and the viscosity of
water from the material database. Prints the largest profile error and
writes the gallery chart.

    python3 tools/property_chart.py docs/figures/gallery/expressions_materials.svg [DIR DIR]

The DIRs are the output directories of couette_viscosity.toml and
couette_transpiration.toml run on their middle mesh (32 x 16 cells) with
a line `mid` at x = 1 m through the nominal cell centres (default
output/<case>); the line takes the value of the cell it passes:

    target/release/foehn flow examples/cases/flow/couette_viscosity.toml \\
        --set 'verification.levels=[1]' \\
        --set 'output.lines=[{name="mid",from=[1.0,0.03125,0.05],to=[1.0,0.96875,0.05],points=16}]'

(the same for couette_transpiration.toml).
"""
import csv
import math
import os
import sys
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "crates", "foehncraft-flow", "data", "materials.toml")


def couette_viscosity(y, a=4.0):
    # μ = μ0 (1 + a y), constant shear stress: u = ln(1 + a y) / ln(1 + a)
    return math.log1p(a * y) / math.log1p(a)


def couette_transpiration(y, re=0.2 / 0.05):
    # blowing V through the lower wall, suction at the upper: u = (e^{Re y} - 1)/(e^{Re} - 1)
    return math.expm1(re * y) / math.expm1(re)


def line(d):
    with open(os.path.join(d, "line_mid.csv")) as f:
        rows = [r for r in f if not r.startswith("#")]
    return [(float(r["y"]), float(r["u"])) for r in csv.DictReader(rows)]


def water():
    with open(DATA, "rb") as f:
        m = tomllib.load(f)["materials"]["water"]["viscosity"]
    return list(zip(m["temperature"], m["values"]))


def panel(out, x0, y0, w, h, title, xlab, ylab, xr, yr, curves, ticks, legend=12):
    sx = lambda v: x0 + (v - xr[0]) / (xr[1] - xr[0]) * w
    sy = lambda v: y0 + h - (v - yr[0]) / (yr[1] - yr[0]) * h
    out.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="none" stroke="#444"/>')
    out.append(f'<text x="{x0 + w / 2}" y="{y0 - 12}" text-anchor="middle" font-weight="bold">{title}</text>')
    out.append(f'<text x="{x0 + w / 2}" y="{y0 + h + 38}" text-anchor="middle">{xlab}</text>')
    out.append(f'<text transform="translate({x0 - 46},{y0 + h / 2}) rotate(-90)" text-anchor="middle">{ylab}</text>')
    for v in ticks[0]:
        out.append(f'<line x1="{sx(v):.1f}" y1="{y0 + h}" x2="{sx(v):.1f}" y2="{y0 + h + 5}" stroke="#444"/>')
        out.append(f'<text x="{sx(v):.1f}" y="{y0 + h + 19}" text-anchor="middle">{v:g}</text>')
    for v in ticks[1]:
        out.append(f'<line x1="{x0 - 5}" y1="{sy(v):.1f}" x2="{x0}" y2="{sy(v):.1f}" stroke="#444"/>')
        out.append(f'<text x="{x0 - 8}" y="{sy(v) + 4:.1f}" text-anchor="end">{v:g}</text>')
    ly = y0 + 18
    for kind, pts, colour, label in curves:
        if kind == "line":
            d = " ".join(f"{sx(a):.1f},{sy(b):.1f}" for a, b in pts)
            out.append(f'<polyline points="{d}" fill="none" stroke="{colour}" stroke-width="2"/>')
        elif kind == "dash":
            d = " ".join(f"{sx(a):.1f},{sy(b):.1f}" for a, b in pts)
            out.append(f'<polyline points="{d}" fill="none" stroke="{colour}" stroke-width="1.5" stroke-dasharray="5,4"/>')
        else:
            for a, b in pts:
                out.append(f'<circle cx="{sx(a):.1f}" cy="{sy(b):.1f}" r="3.2" fill="none" stroke="{colour}" stroke-width="1.5"/>')
        out.append(f'<text x="{x0 + legend}" y="{ly}" fill="{colour}">{label}</text>')
        ly += 16


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    dst = sys.argv[1]
    dirs = sys.argv[2:4] or ["output/couette_viscosity", "output/couette_transpiration"]
    vis, tra = line(dirs[0]), line(dirs[1])
    ev = max(abs(u - couette_viscosity(y)) for y, u in vis)
    et = max(abs(u - couette_transpiration(y)) for y, u in tra)
    print(f"largest difference of the cell values on the line from the exact u: variable viscosity {ev:.3e} m/s, transpiration {et:.3e} m/s")
    ys = [i / 100 for i in range(101)]
    wt = water()
    mu_313 = [mu for t, mu in wt if abs(t - 313.15) < 1e-9]
    w, h = 1020, 380
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
           'font-family="Helvetica, Arial, sans-serif" font-size="12">',
           f'<rect width="{w}" height="{h}" fill="white"/>']
    yt = [0, 0.25, 0.5, 0.75, 1]
    panel(out, 70, 50, 250, 260, "viscosity μ0 (1 + 4 y)", "y (m)", "u (m/s)", (0, 1), (0, 1),
          [("line", [(y, couette_viscosity(y)) for y in ys], "#1f6fb4", "exact"),
           ("dot", vis, "#d9822b", "32 x 16 cells"),
           ("dash", [(0, 0), (1, 1)], "#888888", "constant μ")], (yt, yt))
    panel(out, 410, 50, 250, 260, "blowing and suction, Re = 4", "y (m)", "u (m/s)", (0, 1), (0, 1),
          [("line", [(y, couette_transpiration(y)) for y in ys], "#1f6fb4", "exact"),
           ("dot", tra, "#d9822b", "32 x 16 cells"),
           ("dash", [(0, 0), (1, 1)], "#888888", "no transpiration")], (yt, yt))
    pts = [(t, 1e3 * mu) for t, mu in wt]
    curves = [("line", pts, "#1f6fb4", "water, IAPWS 2008"),
              ("dot", pts, "#1f6fb4", "table points")]
    if mu_313:
        curves.append(("dot", [(313.15, 1e3 * mu_313[0])], "#d9822b", "poiseuille_water, 313.15 K"))
    panel(out, 750, 50, 250, 260, "material database", "T (K)", "μ (mPa s)", (270, 375), (0, 1.9),
          curves, ([280, 300, 320, 340, 360], [0, 0.5, 1.0, 1.5]), legend=80)
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()

"""The outlet of the cylinder benchmark at Re = 100 (D29): the lift
coefficient over the last second of the long channel (outlet 20 D behind
the cylinder) and of the short channel (5.5 D) with a static pressure
outlet (backflow as total pressure) and with a convective outlet. Prints
the benchmark values of each run and writes the gallery chart.

    python3 -m tools.plot.outlet_chart docs/figures/gallery/cylinder_outlets.svg [DIR ...]

The DIRs are the output directories of cylinder_re100.toml,
cylinder_re100_short.toml and cylinder_re100_convective.toml (default
output/<case>); run the cases first.
"""
import csv
import json
import os
import sys

RUNS = [  # label, colour, default directory
    ("long channel, pressure outlet 20 D behind", "#888888", "output/cylinder_re100"),
    ("short channel, pressure outlet 5.5 D behind", "#1f6fb4", "output/cylinder_re100_short"),
    ("short channel, convective outlet 5.5 D behind", "#d9822b", "output/cylinder_re100_convective"),
]
# C_L = 2 F_y / (ρ Ū² D depth) with ρ = 1, Ū = 1 m/s, D = 0.1 m, depth 0.01 m
SCALE = 2.0 / (1.0 * 1.0 * 0.1 * 0.01)
T0, T1 = 7.0, 8.0


def lift(d):
    with open(os.path.join(d, "forces.csv")) as f:
        rows = [r for r in f if not r.startswith("#")]
    return [(float(r["time"]), SCALE * float(r["cylinder_fy"])) for r in csv.DictReader(rows)
            if T0 <= float(r["time"]) <= T1]


def benchmark(d):
    with open(os.path.join(d, "flow.json")) as f:
        rep = json.load(f)
    s = rep["solves"][-1]
    return s.get("benchmark") or {}


def chart(series, dst):
    w, h, left, right, top, bottom = 760, 420, 70, 290, 50, 60
    cl = [c for _, _, s in series for _, c in s]
    y0, y1 = -1.2, 1.2
    assert min(cl) > y0 and max(cl) < y1
    sx = lambda v: left + (v - T0) / (T1 - T0) * (w - left - right)
    sy = lambda v: top + (1 - (v - y0) / (y1 - y0)) * (h - top - bottom)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
           'font-family="Helvetica, Arial, sans-serif" font-size="12">',
           f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
           f'<text x="{(left + w - right) / 2}" y="26" text-anchor="middle" font-size="15">'
           'Lift coefficient at Re = 100: long and short channel</text>',
           f'<rect x="{left}" y="{top}" width="{w - left - right}" height="{h - top - bottom}" fill="none" stroke="#444"/>']
    for k in range(0, 6):
        x = T0 + 0.2 * k
        out.append(f'<line x1="{sx(x):.1f}" y1="{top}" x2="{sx(x):.1f}" y2="{h - bottom}" stroke="#e5e5e5"/>')
        out.append(f'<text x="{sx(x):.1f}" y="{h - bottom + 18}" text-anchor="middle">{x:.1f}</text>')
    for k in range(-2, 3):
        y = 0.5 * k
        out.append(f'<line x1="{left}" y1="{sy(y):.1f}" x2="{w - right}" y2="{sy(y):.1f}" stroke="#e5e5e5"/>')
        out.append(f'<text x="{left - 8}" y="{sy(y) + 4:.1f}" text-anchor="end">{y:.1f}</text>')
    for label, colour, s in series:
        line = " ".join(f"{sx(t):.1f},{sy(c):.1f}" for t, c in s)
        out.append(f'<polyline points="{line}" fill="none" stroke="{colour}" stroke-width="1.6"/>')
    out.append(f'<text x="{(left + w - right) / 2}" y="{h - 16}" text-anchor="middle">time t (s)</text>')
    out.append(f'<text x="18" y="{(top + h - bottom) / 2}" text-anchor="middle" '
               f'transform="rotate(-90 18 {(top + h - bottom) / 2})">lift coefficient C_L</text>')
    yy = top + 20
    for label, colour, s in series:
        out.append(f'<line x1="{w - right + 10}" y1="{yy - 4}" x2="{w - right + 30}" y2="{yy - 4}" '
                   f'stroke="{colour}" stroke-width="1.8"/>')
        out.append(f'<text x="{w - right + 36}" y="{yy}">{label.split(", ")[0]}</text>')
        out.append(f'<text x="{w - right + 36}" y="{yy + 15}">{label.split(", ")[1]}</text>')
        yy += 42
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    dst = sys.argv[1]
    dirs = sys.argv[2:] or [d for _, _, d in RUNS]
    series = []
    for (label, colour, _), d in zip(RUNS, dirs):
        series.append((label, colour, lift(d)))
        b = benchmark(d)
        print(f"{label}: " + ", ".join(f"{v['range']['name']} {v['value']:.4f}" for v in b.get("values", [])))
    chart(series, dst)
    print(f"chart {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

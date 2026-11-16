"""Time histories of the conservation ledger of a transient run (D89):
the kinetic-energy terms of balance.csv (boundary power, viscous and
numerical dissipation, storage) and the relative mass imbalance, drawn as
an SVG chart. The bars of the end of the run are in balance.svg, which the
run writes itself; PNG output follows with the Python package (D61, L21).

    python3 tools/balance_chart.py OUTDIR [CHART]

OUTDIR is the output directory of a transient run; CHART defaults to
OUTDIR/balance_history.svg. No dependencies beyond the standard library.
"""
import csv
import math
import os
import sys

TERMS = [
    ("boundary_power", "boundary power", "#2f6fd0"),
    ("viscous_dissipation", "viscous dissipation", "#d0582f"),
    ("numerical_dissipation", "numerical dissipation", "#8a4fc8"),
    ("storage", "storage dK/dt", "#2e9a5a"),
]


def read(path):
    with open(path) as f:
        rows = [r for r in csv.reader(line for line in f if not line.startswith("#"))]
    head, data = rows[0], rows[1:]
    if "time" not in head:
        sys.exit(f"{path} is the ledger of a steady run (no time column)")
    cols = {h: [float(r[i]) for r in data] for i, h in enumerate(head)}
    return cols


def panel(out, cols, keys, x0, y0, w, h, title, log):
    t = cols["time"]
    vals = [v for k, _, _ in keys for v in cols[k] if not log or v != 0.0]
    if log:
        vals = [abs(v) for v in vals]
        lo = math.floor(math.log10(max(min(vals), 1e-300))) if vals else -20
        hi = math.ceil(math.log10(max(vals))) if vals else 0
        hi = max(hi, lo + 1)
        fy = lambda v: y0 + h - (math.log10(max(abs(v), 10.0 ** lo)) - lo) / (hi - lo) * h
        ticks = [(10.0 ** k, f"1e{k}") for k in range(lo, hi + 1)]
    else:
        lo, hi = (min(vals), max(vals)) if vals else (0.0, 1.0)
        pad = 0.05 * max(hi - lo, 1e-300)
        lo, hi = lo - pad, hi + pad
        fy = lambda v: y0 + h - (v - lo) / (hi - lo) * h
        ticks = [(lo + (hi - lo) * i / 4, f"{lo + (hi - lo) * i / 4:.2e}") for i in range(5)]
    t0, t1 = t[0], t[-1] if t[-1] > t[0] else t[0] + 1
    fx = lambda v: x0 + (v - t0) / (t1 - t0) * w
    out.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="none" stroke="#444"/>')
    out.append(f'<text x="{x0 + w / 2}" y="{y0 - 8}" text-anchor="middle" font-size="13">{title}</text>')
    for v, label in ticks:
        y = fy(v)
        out.append(f'<line x1="{x0}" x2="{x0 + w}" y1="{y:.1f}" y2="{y:.1f}" stroke="#eee"/>')
        out.append(f'<text x="{x0 - 6}" y="{y + 4:.1f}" text-anchor="end" font-size="11">{label}</text>')
    for i, (k, label, colour) in enumerate(keys):
        pts = " ".join(f"{fx(a):.1f},{fy(b):.1f}" for a, b in zip(t, cols[k]))
        out.append(f'<polyline points="{pts}" fill="none" stroke="{colour}" stroke-width="1.6"/>')
        ly = y0 + 16 + 16 * i
        out.append(f'<line x1="{x0 + w + 12}" x2="{x0 + w + 30}" y1="{ly}" y2="{ly}" stroke="{colour}" stroke-width="2"/>')
        out.append(f'<text x="{x0 + w + 36}" y="{ly + 4}" font-size="12">{label}</text>')
    out.append(f'<text x="{x0}" y="{y0 + h + 16}" font-size="11">{t0:.3g} s</text>')
    out.append(f'<text x="{x0 + w}" y="{y0 + h + 16}" text-anchor="end" font-size="11">{t1:.3g} s</text>')


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    d = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else os.path.join(d, "balance_history.svg")
    cols = read(os.path.join(d, "balance.csv"))
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="760" height="620" viewBox="0 0 760 620" '
           'font-family="sans-serif" font-size="12"><rect width="760" height="620" fill="#fff"/>']
    panel(out, cols, TERMS, 90, 40, 460, 300, "kinetic energy (W)", False)
    panel(out, cols, [("mass_relative", "relative mass imbalance", "#555555")], 90, 400, 460, 170,
          "mass (|relative imbalance|)", True)
    out.append('<text x="320" y="610" text-anchor="middle">time</text></svg>')
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"{len(cols['time'])} steps; chart {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

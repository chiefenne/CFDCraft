"""Residual histories of several steady runs in one SVG chart, for example
a case with and without the cell-quality remediation (D50, gallery 26):
the largest scaled momentum residual and the continuity residual against
the outer iteration, on a logarithmic scale.

    python3 -m tools.plot.residuals_chart CHART OUTDIR[=LABEL] [OUTDIR[=LABEL] ...]

Each OUTDIR is the output directory of a steady flow run. The residuals
come from its residuals.csv; a run that diverged writes no residuals.csv,
then they come from the iteration lines of its flow.log (run it with
`--set solver.log_every=1` for every iteration). No dependencies beyond
the standard library.
"""
import csv
import math
import os
import re
import sys

COLOURS = ["#d0582f", "#2f6fd0", "#2e9a5a", "#8a4fc8"]
LINE = re.compile(r"^iter\s+(\d+)\s+momentum\s+(\S+)\s+(\S+)\s+(\S+)\s+continuity\s+(\S+)")


def read(d):
    """(iterations, largest momentum residual, continuity, diverged)."""
    path = os.path.join(d, "residuals.csv")
    if os.path.isfile(path):
        with open(path) as f:
            rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
        it = [int(r["iteration"]) for r in rows]
        mom = [max(float(r[k]) for k in ("momentum_x", "momentum_y", "momentum_z")) for r in rows]
        con = [float(r["continuity"]) for r in rows]
        return it, mom, con, False
    it, mom, con = [], [], []
    with open(os.path.join(d, "flow.log")) as f:
        for line in f:
            m = LINE.match(line.strip())
            if m:
                it.append(int(m.group(1)))
                mom.append(max(float(m.group(k)) for k in (2, 3, 4)))
                con.append(float(m.group(5)))
    if not it:
        sys.exit(f"{d}: neither residuals.csv nor iteration lines in flow.log")
    return it, mom, con, os.path.isfile(os.path.join(d, "divergence.json"))


def main(argv):
    if len(argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    dst, runs = argv[0], []
    for a in argv[1:]:
        d, _, label = a.partition("=")
        runs.append((label or os.path.basename(os.path.normpath(d)),) + read(d))
    vals = [v for r in runs for v in r[2] + r[3] if v > 0]
    lo = math.floor(math.log10(min(vals)))
    hi = max(math.ceil(math.log10(max(vals))), lo + 1)
    n = max(max(r[1]) for r in runs)
    x0, y0, w, h = 80, 30, 520, 330
    fx = lambda i: x0 + (i - 1) / max(n - 1, 1) * w
    fy = lambda v: y0 + h - (math.log10(max(v, 10.0 ** lo)) - lo) / (hi - lo) * h
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="900" height="410" viewBox="0 0 900 410" '
           'font-family="sans-serif" font-size="12"><rect width="900" height="410" fill="#fff"/>',
           f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="none" stroke="#444"/>']
    for k in range(lo, hi + 1):
        y = fy(10.0 ** k)
        out.append(f'<line x1="{x0}" x2="{x0 + w}" y1="{y:.1f}" y2="{y:.1f}" stroke="#eee"/>')
        out.append(f'<text x="{x0 - 6}" y="{y + 4:.1f}" text-anchor="end" font-size="11">1e{k}</text>')
    step = next((s for s in (1, 2, 5, 10, 20, 25, 50, 100, 200, 500, 1000, 2000, 5000) if n / s <= 8), 10000)
    for i in range(step, n + 1, step):
        out.append(f'<text x="{fx(i):.1f}" y="{y0 + h + 16}" text-anchor="middle" font-size="11">{i}</text>')
    legend = 0
    for j, (label, it, mom, con, diverged) in enumerate(runs):
        c = COLOURS[j % len(COLOURS)]
        for vals, dash, what in ((con, "", "continuity"), (mom, ' stroke-dasharray="5,3"', "momentum (largest)")):
            pts = " ".join(f"{fx(a):.1f},{fy(b):.1f}" for a, b in zip(it, vals))
            out.append(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="1.5"{dash}/>')
            ly = y0 + 14 + 18 * legend
            out.append(f'<line x1="{x0 + w + 14}" x2="{x0 + w + 38}" y1="{ly}" y2="{ly}" stroke="{c}" stroke-width="2"{dash}/>')
            out.append(f'<text x="{x0 + w + 44}" y="{ly + 4}" font-size="12">{label}: {what}</text>')
            legend += 1
        if diverged:
            x, y = fx(it[-1]), fy(con[-1])
            out.append(f'<text x="{x:.1f}" y="{y - 10:.1f}" text-anchor="middle" font-size="16" fill="{c}">×</text>')
            out.append(f'<text x="{x:.1f}" y="{y - 26:.1f}" text-anchor="middle" font-size="11" fill="{c}">diverged</text>')
    out.append(f'<text x="{x0 + w / 2}" y="{y0 + h + 36}" text-anchor="middle">outer iteration</text>')
    out.append(f'<text x="18" y="{y0 + h / 2}" text-anchor="middle" transform="rotate(-90 18 {y0 + h / 2})">scaled residual</text>')
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"{len(runs)} runs, up to {n} iterations; chart {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

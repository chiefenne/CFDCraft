"""The cost of the coupled solve along a steady run: the Krylov iterations
of each outer iteration for several runs, one panel per group, as an SVG
chart (gallery 80: the backward-facing step with the BiCGStab and with
the multigrid momentum sweep).

    python3 -m tools.plot.krylov_chart CHART PANEL:LABEL=OUTDIR [...]

Each OUTDIR is the output directory of a steady run with
`algorithm = "coupled"`; the counts come from the column `krylov` of its
residuals.csv. Runs with the same PANEL share a panel; LABEL names the
line, and the legend adds the mean count per outer iteration. No
dependencies beyond the standard library.
"""
import csv
import math
import os
import sys

COLOURS = ["#d0582f", "#2f6fd0", "#2e9a5a", "#8a4fc8"]


def read(d):
    """[(iteration, Krylov iterations)] of a run."""
    path = os.path.join(d, "residuals.csv")
    with open(path) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    if not rows or "krylov" not in rows[0]:
        sys.exit(f"{path}: no column krylov (not a coupled run)")
    return [(int(r["iteration"]), int(r["krylov"])) for r in rows]


def ticks(hi, about):
    """A step of 1, 2 or 5 times a power of ten with about `about` ticks."""
    raw = hi / about
    p = 10 ** math.floor(math.log10(raw))
    return next(m * p for m in (1, 2, 5, 10) if m * p >= raw)


def main(argv):
    if len(argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    dst = argv[0]
    panels = {}
    for a in argv[1:]:
        head, _, path = a.partition("=")
        panel, _, label = head.partition(":")
        panels.setdefault(panel, []).append((label, read(path)))
    pw, ph, gap, x0, y0 = 330, 220, 56, 56, 40
    width = x0 + len(panels) * (pw + gap)
    height = y0 + ph + 100
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'font-family="sans-serif" font-size="12">',
        f'<rect width="{width}" height="{height}" fill="white"/>',
    ]
    vhi = 1.05 * max(v for runs in panels.values() for _, rows in runs for _, v in rows)
    for k, (panel, runs) in enumerate(panels.items()):
        px = x0 + k * (pw + gap)
        ihi = max(i for _, rows in runs for i, _ in rows)
        fx = lambda i: px + i / ihi * pw
        fy = lambda v: y0 + ph - v / vhi * ph
        out.append(f'<rect x="{px}" y="{y0}" width="{pw}" height="{ph}" fill="none" stroke="#444"/>')
        out.append(
            f'<text x="{px + pw / 2}" y="{y0 - 12}" text-anchor="middle" font-size="13">{panel}</text>'
        )
        step = ticks(ihi, 6)
        t = 0
        while t <= ihi:
            x = fx(t)
            out.append(f'<line x1="{x:.1f}" y1="{y0 + ph}" x2="{x:.1f}" y2="{y0 + ph + 4}" stroke="#444"/>')
            out.append(f'<text x="{x:.1f}" y="{y0 + ph + 17}" text-anchor="middle">{t:g}</text>')
            t += step
        step = ticks(vhi, 6)
        v = 0
        while v <= vhi:
            y = fy(v)
            out.append(f'<line x1="{px}" y1="{y:.1f}" x2="{px + pw}" y2="{y:.1f}" stroke="#e4e4e4"/>')
            out.append(f'<text x="{px - 7}" y="{y + 4:.1f}" text-anchor="end">{v:g}</text>')
            v += step
        for j, (label, rows) in enumerate(runs):
            c = COLOURS[j % len(COLOURS)]
            pts = " ".join(f"{fx(i):.1f},{fy(v):.1f}" for i, v in rows)
            out.append(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="1.4"/>')
            mean = sum(v for _, v in rows) / len(rows)
            ly = y0 + ph + 62 + 17 * j
            out.append(f'<line x1="{px}" y1="{ly - 4}" x2="{px + 22}" y2="{ly - 4}" stroke="{c}" stroke-width="1.8"/>')
            out.append(f'<text x="{px + 28}" y="{ly}">{label}: {mean:.1f} per outer iteration</text>')
        out.append(
            f'<text x="{px + pw / 2}" y="{y0 + ph + 36}" text-anchor="middle" fill="#444">outer iteration</text>'
        )
    out.append(
        f'<text x="16" y="{y0 + ph / 2}" text-anchor="middle" '
        f'transform="rotate(-90 16 {y0 + ph / 2})">Krylov iterations of the coupled solve</text>'
    )
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

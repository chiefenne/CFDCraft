"""How fast a monitored value settles: the value against the outer
iteration for several runs, one panel per group, as an SVG chart (gallery
78: the reattachment length of the backward-facing step with SIMPLEC and
with the coupled algorithm).

    python3 -m tools.plot.settling_chart CHART YLABEL PANEL:LABEL=FILE.csv [...]

Each FILE.csv has two columns, the iteration and the value (a header line
and rows without a finite value are skipped). Runs with the same PANEL
share a panel; LABEL names the line.
The band of ±1 % around the last value of the longest run of a panel is
shaded, and the iteration from which a run stays inside it is marked. The
iteration axis is logarithmic. No dependencies beyond the standard
library.
"""
import math
import sys

COLOURS = ["#2f6fd0", "#d0582f", "#2e9a5a", "#8a4fc8"]


def read(path):
    """[(iteration, value)] of a two-column CSV file."""
    rows = []
    with open(path) as f:
        for line in f:
            p = line.replace(";", ",").split(",")
            try:
                i, v = float(p[0]), float(p[1])
            except (ValueError, IndexError):
                continue
            if math.isfinite(v):
                rows.append((i, v))
    if not rows:
        sys.exit(f"{path}: no rows")
    return rows


def settled(rows, final, band=0.01):
    """The first iteration from which the value stays within the band."""
    it = None
    for i, v in rows:
        if abs(v - final) <= band * abs(final):
            it = i if it is None else it
        else:
            it = None
    return it


def main(argv):
    if len(argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    dst, ylabel = argv[0], argv[1]
    panels = {}
    for a in argv[2:]:
        head, _, path = a.partition("=")
        panel, _, label = head.partition(":")
        panels.setdefault(panel, []).append((label, read(path)))
    pw, ph, gap, x0, y0 = 270, 220, 56, 62, 40
    width = x0 + len(panels) * (pw + gap)
    height = y0 + ph + 100
    its = [i for runs in panels.values() for _, rows in runs for i, _ in rows]
    lo, hi = math.log10(max(min(its), 1.0)), math.log10(max(its))
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'font-family="sans-serif" font-size="12">',
        f'<rect width="{width}" height="{height}" fill="white"/>',
    ]
    for k, (panel, runs) in enumerate(panels.items()):
        px = x0 + k * (pw + gap)
        final = max(runs, key=lambda r: r[1][-1][0])[1][-1][1]
        vals = [v for _, rows in runs for _, v in rows]
        vlo, vhi = min(vals + [0.98 * final]), max(vals + [1.02 * final])
        pad = 0.05 * (vhi - vlo)
        vlo, vhi = vlo - pad, vhi + pad
        fx = lambda i: px + (math.log10(max(i, 1.0)) - lo) / (hi - lo) * pw
        fy = lambda v: y0 + ph - (v - vlo) / (vhi - vlo) * ph
        out.append(
            f'<rect x="{px}" y="{fy(1.01 * final):.1f}" width="{pw}" '
            f'height="{fy(0.99 * final) - fy(1.01 * final):.1f}" fill="#e9e9e9"/>'
        )
        out.append(f'<rect x="{px}" y="{y0}" width="{pw}" height="{ph}" fill="none" stroke="#444"/>')
        out.append(
            f'<text x="{px + pw / 2}" y="{y0 - 12}" text-anchor="middle" font-size="13">{panel}</text>'
        )
        d = 1
        while d <= 10 ** hi:
            for m in (1, 2, 5):
                t = m * d
                if 10 ** lo <= t <= 10 ** hi:
                    x = fx(t)
                    out.append(f'<line x1="{x:.1f}" y1="{y0 + ph}" x2="{x:.1f}" y2="{y0 + ph + 4}" stroke="#444"/>')
                    out.append(f'<text x="{x:.1f}" y="{y0 + ph + 17}" text-anchor="middle">{t}</text>')
            d *= 10
        step = 10 ** math.floor(math.log10(vhi - vlo))
        if (vhi - vlo) / step < 3:
            step /= 2
        v = math.ceil(vlo / step) * step
        while v <= vhi:
            y = fy(v)
            out.append(f'<line x1="{px - 4}" y1="{y:.1f}" x2="{px}" y2="{y:.1f}" stroke="#444"/>')
            out.append(f'<text x="{px - 7}" y="{y + 4:.1f}" text-anchor="end">{v:g}</text>')
            v += step
        for j, (label, rows) in enumerate(runs):
            c = COLOURS[j % len(COLOURS)]
            pts = " ".join(f"{fx(i):.1f},{fy(v):.1f}" for i, v in rows)
            out.append(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="1.8"/>')
            s = settled(rows, final)
            note = ""
            if s is not None:
                out.append(
                    f'<line x1="{fx(s):.1f}" y1="{y0}" x2="{fx(s):.1f}" y2="{y0 + ph}" '
                    f'stroke="{c}" stroke-dasharray="4 3"/>'
                )
                note = f": within 1 % from {s:g}"
            ly = y0 + ph + 62 + 17 * j
            out.append(f'<line x1="{px}" y1="{ly - 4}" x2="{px + 22}" y2="{ly - 4}" stroke="{c}" stroke-width="1.8"/>')
            out.append(f'<text x="{px + 28}" y="{ly}">{label}{note}</text>')
        out.append(
            f'<text x="{px + pw / 2}" y="{y0 + ph + 36}" text-anchor="middle" fill="#444">outer iteration</text>'
        )
    out.append(
        f'<text x="16" y="{y0 + ph / 2}" text-anchor="middle" '
        f'transform="rotate(-90 16 {y0 + ph / 2})">{ylabel}</text>'
    )
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

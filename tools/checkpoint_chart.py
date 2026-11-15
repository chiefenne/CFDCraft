"""Draws the gallery chart of binomial checkpointing (D59): primal steps
per time step of a reverse sweep over n steps against the number of
snapshots, from the closed form t n - beta(s + 2, t - 1) of the schedule in
crates/foehncraft-flow/src/checkpoint.rs (its tests check the closed form
against a dynamic program).

    python3 tools/checkpoint_chart.py docs/figures/gallery/checkpoint_steps.svg
"""
import math
import sys


def beta(s, t):
    return math.comb(s + t, t) if t >= 0 else 0


def steps(n, s):
    """Primal steps of the reverse sweep of n steps with s snapshots."""
    if n <= 1:
        return 0
    if s == 0:
        return n * (n - 1) // 2
    t = 0
    while beta(s + 1, t) < n:
        t += 1
    return t * n - beta(s + 2, t - 1)


def main(dst):
    w, h, left, right, top, bottom = 680, 420, 70, 150, 50, 60
    smax = 40
    ymax = 12.0
    x = lambda s: left + (s - 1) / (smax - 1) * (w - left - right)
    y = lambda v: top + (1 - v / ymax) * (h - top - bottom)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
           'font-family="Helvetica, Arial, sans-serif" font-size="12">',
           f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
           f'<text x="{(left + w - right) / 2}" y="26" text-anchor="middle" font-size="15">'
           'Reverse sweep: primal steps per time step</text>',
           f'<rect x="{left}" y="{top}" width="{w - left - right}" height="{h - top - bottom}" fill="none" stroke="#444"/>']
    for v in range(0, 13, 2):
        out.append(f'<line x1="{left}" x2="{w - right}" y1="{y(v):.1f}" y2="{y(v):.1f}" stroke="#ddd"/>'
                   f'<text x="{left - 6}" y="{y(v) + 4:.1f}" text-anchor="end">{v}</text>')
    for s in (1, 5, 10, 15, 20, 25, 30, 35, 40):
        out.append(f'<text x="{x(s):.1f}" y="{h - bottom + 18}" text-anchor="middle">{s}</text>')
    out.append(f'<text x="{(left + w - right) / 2}" y="{h - 16}" text-anchor="middle">snapshots in memory</text>')
    colors = ["#2a6fdb", "#d9822b", "#2f9e44"]
    for k, n in enumerate((1_000, 10_000, 100_000)):
        pts = []
        for s in range(1, smax + 1):
            v = min(steps(n, s) / n, ymax)
            pts.append(f"{x(s):.1f},{y(v):.1f}")
        out.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{colors[k]}" stroke-width="2"/>')
        ly = top + 20 + 22 * k
        out.append(f'<line x1="{w - right + 14}" x2="{w - right + 38}" y1="{ly}" y2="{ly}" stroke="{colors[k]}" stroke-width="2"/>'
                   f'<text x="{w - right + 44}" y="{ly + 4}">{n:,} steps</text>')
    out.append("</svg>")
    open(dst, "w").write("".join(out))
    for n in (1_000, 10_000, 100_000):
        print(n, [(s, round(steps(n, s) / n, 2)) for s in (5, 10, 20, 40)])


main(sys.argv[1])

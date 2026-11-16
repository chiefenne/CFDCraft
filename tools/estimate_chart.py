"""Checks the estimate of `foehn flow --check` (D81) against runs: the peak
resident memory and the sizes of fields.vtkhdf and restart.h5 of short
runs of catalogue cases (40 outer iterations, or 20 time steps), drawn
against the estimate. Prints one line per case and writes the gallery
chart. Unix only (peak memory from wait4).

    python3 tools/estimate_chart.py docs/figures/gallery/estimate.svg [FOEHN] [OUTDIR]

FOEHN defaults to target/release/foehn, OUTDIR to output/estimate_check.
"""
import json
import math
import os
import subprocess
import sys

CASES = [
    "couette_hex", "poiseuille_tet", "poiseuille_mixed", "cavity_re100",
    "cavity_re100_tet3d", "cavity_re100_poly3d", "cavity_re1000_poly3d",
    "cylinder_re20", "cylinder_re100", "cavity_re100_transient", "taylor_green_bdf2",
]
TRANSIENT = ("transient", "taylor", "cylinder_re100")


def measure(foehn, case, outdir):
    path = f"examples/cases/flow/{case}.toml"
    out = os.path.join(outdir, case)
    short = (["--set", "time.steps=20", "--set", "time.end_time=1e9", "--set", "verification.levels=[]"]
             if any(k in case for k in TRANSIENT) else ["--max-iterations", "40"])
    args = [foehn, "--quiet", "flow", path, "--out", out, "--threads", "2"] + short
    p = subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    _, _, ru = os.wait4(p.pid, 0)
    with open(os.path.join(out, "result.json")) as f:
        sizes = {x["name"]: x["bytes"] for x in json.load(f)["files"]}
    check = subprocess.run([foehn, "--json", "flow", path, "--check", "--out", out] + short,
                           capture_output=True, check=True).stdout
    e = json.loads(check)["estimate"]
    est = {o["name"]: o["bytes"] for o in e["outputs"]}
    rss = ru.ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    rows = [("memory", rss, e["host_memory"])]
    for name in ("fields.vtkhdf", "restart.h5"):
        if name in sizes and name in est:
            rows.append((name, sizes[name], est[name]))
    return e["cells"], e["faces"], rows


def chart(points, dst):
    w, h, left, right, top, bottom = 640, 460, 80, 170, 50, 60
    lo, hi = 1e4, 1e9
    sx = lambda v: left + (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * (w - left - right)
    sy = lambda v: top + (1 - (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))) * (h - top - bottom)
    colour = {"memory": "#1f6fb4", "fields.vtkhdf": "#d9822b", "restart.h5": "#2e9e5b"}
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
           'font-family="Helvetica, Arial, sans-serif" font-size="12">',
           f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
           f'<text x="{(left + w - right) / 2}" y="26" text-anchor="middle" font-size="15">'
           'Estimate of foehn flow --check against the run</text>',
           f'<rect x="{left}" y="{top}" width="{w - left - right}" height="{h - top - bottom}" fill="none" stroke="#444"/>']
    for k in range(4, 10):
        v = 10.0 ** k
        out.append(f'<line x1="{sx(v):.1f}" y1="{top}" x2="{sx(v):.1f}" y2="{h - bottom}" stroke="#e5e5e5"/>')
        out.append(f'<line x1="{left}" y1="{sy(v):.1f}" x2="{w - right}" y2="{sy(v):.1f}" stroke="#e5e5e5"/>')
        label = ["10 kB", "100 kB", "1 MB", "10 MB", "100 MB", "1 GB"][k - 4]
        out.append(f'<text x="{sx(v):.1f}" y="{h - bottom + 18}" text-anchor="middle">{label}</text>')
        out.append(f'<text x="{left - 8}" y="{sy(v) + 4:.1f}" text-anchor="end">{label}</text>')
    # the diagonal and a band of +-20 %
    for f, dash in ((1.0, ""), (1.2, ' stroke-dasharray="4 4"'), (1 / 1.2, ' stroke-dasharray="4 4"')):
        x0, x1 = max(lo, lo / f), min(hi, hi / f)  # y = f x inside the frame
        out.append(f'<line x1="{sx(x0):.1f}" y1="{sy(x0 * f):.1f}" x2="{sx(x1):.1f}" y2="{sy(x1 * f):.1f}" '
                   f'stroke="#888"{dash}/>')
    for kind, measured, estimated in points:
        out.append(f'<circle cx="{sx(measured):.1f}" cy="{sy(estimated):.1f}" r="4" fill="{colour[kind]}" '
                   'fill-opacity="0.85"/>')
    out.append(f'<text x="{(left + w - right) / 2}" y="{h - 16}" text-anchor="middle">measured</text>')
    out.append(f'<text x="18" y="{(top + h - bottom) / 2}" text-anchor="middle" '
               f'transform="rotate(-90 18 {(top + h - bottom) / 2})">estimated</text>')
    for i, (kind, label) in enumerate((("memory", "peak memory"), ("fields.vtkhdf", "fields.vtkhdf"),
                                       ("restart.h5", "restart.h5"))):
        yy = top + 20 + 22 * i
        out.append(f'<circle cx="{w - right + 20}" cy="{yy - 4}" r="4" fill="{colour[kind]}"/>')
        out.append(f'<text x="{w - right + 30}" y="{yy}">{label}</text>')
    yy = top + 20 + 22 * 3
    out.append(f'<line x1="{w - right + 12}" y1="{yy - 4}" x2="{w - right + 28}" y2="{yy - 4}" stroke="#888" '
               'stroke-dasharray="4 4"/>')
    out.append(f'<text x="{w - right + 30}" y="{yy}">± 20 %</text>')
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    dst = sys.argv[1]
    foehn = sys.argv[2] if len(sys.argv) > 2 else "target/release/foehn"
    outdir = sys.argv[3] if len(sys.argv) > 3 else "output/estimate_check"
    points, worst = [], 0.0
    for c in CASES:
        cells, faces, rows = measure(foehn, c, outdir)
        parts = []
        for kind, m, e in rows:
            points.append((kind, m, e))
            worst = max(worst, abs(e / m - 1))
            parts.append(f"{kind} {m / 1e6:.2f} MB est {e / 1e6:.2f} ({100 * (e / m - 1):+.0f} %)")
        print(f"{c:24s} {cells:7d} cells {faces:7d} faces  " + "  ".join(parts))
    chart(points, dst)
    print(f"largest deviation {100 * worst:.0f} %; chart {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

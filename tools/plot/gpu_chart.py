"""Draws the gallery chart of the GPU flow loop (time per cell and
iteration or step, CPU against GPU, log axis) from the report table of
tools/bench/stage5_cloud.py.

    python3 -m tools.plot.gpu_chart docs/benchmarks/stage5_gpu_t4.md docs/figures/gallery/gpu_flow_times.svg
"""
import math
import re
import sys


def rows(path):
    out = []
    for line in open(path):
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        if len(c) == 7 and c[1].isdigit():
            gpu = int(re.match(r"\d+", c[5]).group())
            out.append((c[0], int(c[1]), int(c[4]), gpu))
    return out


def main(src, dst):
    r = rows(src)
    w, h, left, right, top = 760, 60 + 46 * len(r) + 60, 190, 40, 56
    lo, hi = math.log10(300), math.log10(40000)
    x = lambda v: left + (math.log10(v) - lo) / (hi - lo) * (w - left - right)
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
         'font-family="Helvetica, Arial, sans-serif" font-size="12">',
         f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
         f'<text x="{w/2}" y="26" text-anchor="middle" font-size="15">'
         'Flow loop, time per cell and iteration (step): CPU 2 threads, Tesla T4</text>']
    for e in range(3, 5):
        for m in (1, 2, 5):
            v = m * 10**e if e > 2 else m * 100
            if v < 300 or v > 40000:
                continue
            s.append(f'<line x1="{x(v):.1f}" x2="{x(v):.1f}" y1="{top-8}" y2="{h-50}" stroke="#ddd"/>'
                     f'<text x="{x(v):.1f}" y="{h-34}" text-anchor="middle">{v}</text>')
    s.append(f'<text x="{(left+w-right)/2}" y="{h-12}" text-anchor="middle">ns per cell and iteration (or step)</text>')
    for k, (name, n, cpu, gpu) in enumerate(r):
        y = top + 46 * k
        s.append(f'<text x="{left-8}" y="{y+17}" text-anchor="end">{name}</text>'
                 f'<text x="{left-8}" y="{y+32}" text-anchor="end" fill="#666">{n:,} cells, ×{cpu/gpu:.1f}</text>')
        for j, (v, col, lab) in enumerate(((cpu, "#8a8a8a", "CPU"), (gpu, "#2a6fdb", "GPU"))):
            yy = y + 4 + 15 * j
            s.append(f'<rect x="{left}" y="{yy}" width="{x(v)-left:.1f}" height="13" fill="{col}"/>'
                     f'<text x="{x(v)+4:.1f}" y="{yy+11}">{v}</text>')
    s.append(f'<rect x="{w-right-110}" y="{top-34}" width="12" height="12" fill="#8a8a8a"/><text x="{w-right-94}" y="{top-24}">CPU</text>'
             f'<rect x="{w-right-56}" y="{top-34}" width="12" height="12" fill="#2a6fdb"/><text x="{w-right-40}" y="{top-24}">GPU</text>')
    s.append("</svg>")
    open(dst, "w").write("".join(s))


main(*sys.argv[1:3])

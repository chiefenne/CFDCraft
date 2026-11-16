"""Discretisation uncertainty of two level sweeps (D67): the wall shear
force of Couette flow with μ = μ0 (1 + a y) against its exact value, and
the drag coefficient of the cylinder at Re = 20 against the reference range
of Schäfer & Turek (1996). Each panel shows the three solutions over the
relative cell size h/h1, the curve φ_ext + (φ1 − φ_ext) (h/h1)^p through
them, the extrapolated value and the band φ1 ± U of the grid convergence
index. Prints the numbers and writes the gallery chart.

    python3 -m tools.plot.gci_chart docs/figures/gallery/gci.svg [COUETTE CYLINDER]

COUETTE and CYLINDER are the output directories of couette_viscosity.toml
and cylinder_re20_gci.toml (default output/<case>).
"""
import json
import math
import os
import sys

MU0, U, A, AREA = 0.1, 1.0, 4.0, 0.2  # couette_viscosity.toml, lower wall 2 m x 0.1 m


def estimate(result, quantity):
    for e in result["uncertainty"]:
        if e["quantity"] == quantity:
            return e
    raise SystemExit(f"no estimate of {quantity} in the result (a sweep of three meshes?)")


def panel(out, x0, y0, w, h, e, title, ylabel, exact=None, band=None, scale=1.0):
    """One panel: the three solutions at h/h1 = 1, r21, r21 r32."""
    r21, r32 = e["ratios"]
    hs = [1.0, r21, r21 * r32]
    phi = [v * scale for v in e["values"]]
    ext = e["extrapolated"] * scale
    u = e["band"] * scale
    p = e["order"]
    lo = min(phi + [ext, phi[0] - u] + ([exact] if exact else []) + (band or []))
    hi = max(phi + [ext, phi[0] + u] + ([exact] if exact else []) + (band or []))
    pad = 0.08 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    xmax = hs[-1] * 1.1
    sx = lambda v: x0 + v / xmax * w  # noqa: E731
    sy = lambda v: y0 + (1 - (v - lo) / (hi - lo)) * h  # noqa: E731
    out.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="none" stroke="#444"/>')
    out.append(f'<text x="{x0 + w / 2}" y="{y0 - 10}" text-anchor="middle" font-weight="bold">{title}</text>')
    out.append(f'<text x="{x0 + w / 2}" y="{y0 + h + 36}" text-anchor="middle">relative cell size h / h1</text>')
    out.append(f'<text transform="translate({x0 - 58},{y0 + h / 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>')
    for v in (0, 1, 2, 3, 4):
        if v <= xmax:
            out.append(f'<text x="{sx(v):.1f}" y="{y0 + h + 16}" text-anchor="middle">{v}</text>')
    for i in range(5):
        v = lo + (hi - lo) * (i + 0.5) / 5
        out.append(f'<text x="{x0 - 6}" y="{sy(v) + 4:.1f}" text-anchor="end">{v:.4g}</text>')
    if band:
        out.append(f'<rect x="{x0}" y="{sy(band[1]):.1f}" width="{w}" height="{sy(band[0]) - sy(band[1]):.1f}" '
                   'fill="#2a9d8f" fill-opacity="0.18"/>')
        out.append(f'<text x="{x0 + w - 6}" y="{sy(band[1]) - 4:.1f}" text-anchor="end" fill="#2a9d8f">'
                   'reference range</text>')
    if exact is not None:
        out.append(f'<line x1="{x0}" y1="{sy(exact):.1f}" x2="{x0 + w}" y2="{sy(exact):.1f}" '
                   'stroke="#2a9d8f" stroke-dasharray="5 3"/>')
        out.append(f'<text x="{x0 + w - 6}" y="{sy(exact) - 5:.1f}" text-anchor="end" fill="#2a9d8f">exact</text>')
    pts = [(t, ext + (phi[0] - ext) * t**p) for t in [xmax * i / 200 for i in range(201)]]
    pts = [(t, v) for t, v in pts if lo <= v <= hi]
    d = " ".join(f"{sx(t):.1f},{sy(v):.1f}" for t, v in pts)
    out.append(f'<polyline points="{d}" fill="none" stroke="#1f6fb4" stroke-width="1.5"/>')
    out.append(f'<line x1="{sx(1):.1f}" y1="{sy(phi[0] - u):.1f}" x2="{sx(1):.1f}" y2="{sy(phi[0] + u):.1f}" '
               'stroke="#d9822b" stroke-width="2"/>')
    for t in (sy(phi[0] - u), sy(phi[0] + u)):
        out.append(f'<line x1="{sx(1) - 5:.1f}" y1="{t:.1f}" x2="{sx(1) + 5:.1f}" y2="{t:.1f}" stroke="#d9822b" stroke-width="2"/>')
    for t, v in zip(hs, phi):
        out.append(f'<circle cx="{sx(t):.1f}" cy="{sy(v):.1f}" r="4.5" fill="#1f6fb4"/>')
    out.append(f'<rect x="{sx(0) - 4:.1f}" y="{sy(ext) - 4:.1f}" width="8" height="8" fill="#d9822b"/>')
    out.append(f'<text x="{x0 + 10}" y="{y0 + h - 28}">p = {p:.2f}, φ_ext = {ext:.6g}</text>')
    gci = e.get("gci")
    out.append(f'<text x="{x0 + 10}" y="{y0 + h - 12}">φ1 = {phi[0]:.6g} ± {u:.2g}'
               + (f" (GCI {100 * gci:.2f} %)" if gci is not None else "") + "</text>")


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    dst = argv[0]
    couette = argv[1] if len(argv) > 1 else os.path.join("output", "couette_viscosity")
    cyl = argv[2] if len(argv) > 2 else os.path.join("output", "cylinder_re20_gci")
    rc = json.load(open(os.path.join(couette, "result.json"), encoding="utf-8"))
    ry = json.load(open(os.path.join(cyl, "result.json"), encoding="utf-8"))
    ec = estimate(rc, "patches.ymin.force_x")
    ed = estimate(ry, "benchmark.C_D")
    exact = MU0 * U * A / math.log1p(A) * AREA
    rng = next(v["range"] for v in ry["solves"][-1]["benchmark"]["values"] if v["range"]["name"] == "C_D")
    print(f"Couette wall force: finest {ec['values'][0]:.6e}, extrapolated {ec['extrapolated']:.6e}, "
          f"exact {exact:.6e}; order {ec['order']:.3f}, band {ec['band']:.2e}")
    print(f"cylinder C_D: finest {ed['values'][0]:.5f}, extrapolated {ed['extrapolated']:.5f}, "
          f"reference {rng['lower']}-{rng['upper']}; order {ed['order']:.3f}, band {ed['band']:.2e}")
    w, h = 900, 400
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
        'font-family="Helvetica, Arial, sans-serif" font-size="12">',
        f'<rect width="{w}" height="{h}" fill="white"/>',
    ]
    panel(out, 90, 50, 320, 280, ec, "Couette flow, μ = μ0 (1 + 4 y): lower-wall force",
          "F_x (N)", exact=exact)
    panel(out, 560, 50, 320, 280, ed, "Cylinder at Re = 20: drag coefficient", "C_D",
          band=[rng["lower"], rng["upper"]])
    out.append(f'<text x="{w / 2}" y="{h - 8}" text-anchor="middle" fill="#555">circles: the three meshes; '
               'square: extrapolated to h = 0; bar: φ1 ± U (Celik et al. 2008)</text>')
    out.append("</svg>")
    with open(dst, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main(sys.argv[1:])

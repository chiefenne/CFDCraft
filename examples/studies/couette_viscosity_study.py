"""A parametric study from Python (D61, D62): plane Couette flow with the
viscosity μ0 (1 + a y), for six viscosity gradients a.

    python examples/studies/couette_viscosity_study.py [OUT] [--chart FILE.svg]

The base case is examples/cases/flow/couette_viscosity.toml on its middle
mesh (32 x 16 cells). For each a the study writes a case file of its own
(OUT/v0000/v0000.toml, ...), runs it, and collects the numbers of the run
into OUT/study.csv. The shear stress is the same on every layer,
τ = μ0 U a / ln(1 + a); the study compares the wall shear stress of the
runs (the force on the lower wall over its area) with it. Run it again and
the finished variants are read back instead of run (resume).

Needs the package (pip install ./python) or the `foehn` executable
(FOEHN=target/release/foehn, backend "cli").
"""

import math
import sys
from pathlib import Path

from foehncraft import Case
from foehncraft.study import Study, listing

ROOT = Path(__file__).resolve().parents[2]
MU0, U = 0.1, 1.0
GRADIENTS = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]


def exact_shear(a: float) -> float:
    return MU0 * U * a / math.log1p(a)


def chart(rows, dst: Path) -> None:
    """Wall shear stress against a: exact curve and the runs (SVG)."""
    w, h, left, right, top, bottom = 640, 380, 70, 30, 40, 55
    a0, a1, t0, t1 = 0.0, 17.0, 0.0, 0.7
    sx = lambda v: left + (v - a0) / (a1 - a0) * (w - left - right)  # noqa: E731
    sy = lambda v: top + (1 - (v - t0) / (t1 - t0)) * (h - top - bottom)  # noqa: E731
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" '
        f'height="{h}" font-family="Helvetica, Arial, sans-serif" font-size="12">',
        f'<rect width="{w}" height="{h}" fill="white"/>',
        f'<rect x="{left}" y="{top}" width="{w - left - right}" height="{h - top - bottom}" '
        'fill="none" stroke="#444"/>',
        f'<text x="{w / 2}" y="22" text-anchor="middle" font-weight="bold">'
        "Couette flow with μ = μ0 (1 + a y): wall shear stress</text>",
        f'<text x="{w / 2}" y="{h - 12}" text-anchor="middle">viscosity gradient a (1/m)</text>',
        f'<text transform="translate(20,{h / 2}) rotate(-90)" text-anchor="middle">τ_w (Pa)</text>',
    ]
    for v in (0, 4, 8, 12, 16):
        out.append(f'<text x="{sx(v)}" y="{h - bottom + 18}" text-anchor="middle">{v}</text>')
    for v in (0, 0.2, 0.4, 0.6):
        out.append(f'<text x="{left - 8}" y="{sy(v) + 4}" text-anchor="end">{v:g}</text>')
    pts = [(a, exact_shear(a)) for a in [0.05 * i for i in range(1, 341)]]
    d = " ".join(f"{sx(a):.1f},{sy(v):.1f}" for a, v in pts)
    out.append(f'<polyline points="{d}" fill="none" stroke="#1f6fb4" stroke-width="2"/>')
    for r in rows:
        a, tau = r["a"], r["tau"]
        out.append(
            f'<circle cx="{sx(a):.1f}" cy="{sy(tau):.1f}" r="4" fill="none" '
            'stroke="#d9822b" stroke-width="1.8"/>'
        )
    out.append(f'<text x="{left + 14}" y="{top + 20}" fill="#1f6fb4">exact μ0 U a / ln(1 + a)</text>')
    out.append(
        f'<text x="{left + 14}" y="{top + 36}" fill="#d9822b">study: 6 runs, 32 x 16 cells</text>'
    )
    out.append("</svg>")
    dst.write_text("\n".join(out) + "\n")


def main(argv):
    out = Path(argv[0]) if argv and not argv[0].startswith("--") else Path("output/couette_study")
    base = Case.load(ROOT / "examples/cases/flow/couette_viscosity.toml")
    base["verification"] = {"levels": [1], "require_convergence": True}
    variants = listing({"constants.a": a, "analytic.gradient": a} for a in GRADIENTS)
    rows = Study(base, variants, out).run()
    results = []
    for a, r in zip(GRADIENTS, rows):
        tau = r["boundaries.ymin.force.0"] / r["boundaries.ymin.area"]
        results.append({"a": a, "tau": tau})
        print(
            f"a = {a:5.1f}  {r['status']:<9}  τ_w {tau:.6f} Pa, exact {exact_shear(a):.6f} "
            f"({tau / exact_shear(a) - 1:+.2e}), velocity error L2 {r['errors.l2']:.2e}"
        )
    print(f"table: {out / 'study.csv'}")
    if "--chart" in argv:
        dst = Path(argv[argv.index("--chart") + 1])
        chart(results, dst)
        print(f"wrote {dst}")


if __name__ == "__main__":
    main(sys.argv[1:])

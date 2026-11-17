"""The Moody chart from a parametric study (D62, stage 6 part 6c): the
friction factor of fully developed rough pipe flow with the standard k-ε
model, scalable wall functions and the roughness shift of Cebeci &
Bradshaw, against Colebrook's equation.

    python examples/studies/moody_study.py [OUT] [--chart FILE.svg] [--cs 0.5] [--jobs 4]

The base case is examples/cases/flow/pipe_rough_keps.toml (wedge pipe,
D = 1 m, U = 1 m/s, periodic module). The study runs Re = 3·10⁴, 10⁵,
10⁶ and 10⁷ (ν = 1/Re) with k_s/D = 0, 10⁻⁴, 10⁻³ and 10⁻², 16 runs. The
mesh has N uniform cells across the radius (6, 12, 40, 80 by Re, so that
the first-cell y⁺ lies between 60 and 1000), fewer where the roughness
would reach beyond the first-cell distance (BC027). The friction factor
is f = 2 F D/(ρ U²) with the driving force F of the periodic pair.
`--cs` sets the roughness constant C_s (default 0.5).

Findings of 6c (C_s = 0.5): smooth pipes −2 … −5 % against Colebrook,
fully rough +8 … +11 % (C_s = 0.4: +1 … +3 %); in the transitional range
the shift follows the inflectional curve of Nikuradse's sand grains and
lies up to 17 % below Colebrook's monotonic interpolation of commercial
pipes.

Needs the package (pip install ./python) or the `foehn` executable
(FOEHN=target/release/foehn, backend "cli"). About 25 minutes on one
core (`--jobs 4`: 8 minutes).
"""

import math
import sys
from pathlib import Path

from foehncraft import Case
from foehncraft.study import Study, listing

ROOT = Path(__file__).resolve().parents[2]
CELLS = {3e4: 6, 1e5: 12, 1e6: 40, 1e7: 80}
ROUGHNESS = [0.0, 1e-4, 1e-3, 1e-2]


def colebrook(re: float, ks: float) -> float:
    """Darcy friction factor of Colebrook's equation (k_s relative to D)."""
    f = 0.02
    for _ in range(100):
        f = (-2.0 * math.log10(ks / 3.7 + 2.51 / (re * math.sqrt(f)))) ** -2
    return f


def cells(re: float, ks: float) -> int:
    """Cells across the radius: the first-cell distance R/(2N) above 1.2 k_s."""
    n = CELLS[re]
    if ks > 0.0 and 0.25 / n <= 1.2 * ks:
        n = int(0.25 / (1.2 * ks))
    return n


def chart(rows, dst: Path) -> None:
    """f over Re on log axes: Colebrook's curves and the runs (SVG)."""
    w, h, left, right, top, bottom = 680, 440, 70, 150, 40, 55
    lx0, lx1, ly0, ly1 = math.log10(1e4), math.log10(3e7), math.log10(0.006), math.log10(0.06)

    def sx(re):
        return left + (math.log10(re) - lx0) / (lx1 - lx0) * (w - left - right)

    def sy(f):
        return top + (1 - (math.log10(f) - ly0) / (ly1 - ly0)) * (h - top - bottom)

    colours = ["#1f6fb4", "#2a9d5b", "#d9822b", "#b03a48"]
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" '
        f'height="{h}" font-family="Helvetica, Arial, sans-serif" font-size="12">',
        f'<rect width="{w}" height="{h}" fill="white"/>',
        f'<rect x="{left}" y="{top}" width="{w - left - right}" height="{h - top - bottom}" '
        'fill="none" stroke="#444"/>',
        f'<text x="{(w - right + left) / 2}" y="22" text-anchor="middle" font-weight="bold">'
        "Moody chart: k-ε, scalable wall functions, Cebeci–Bradshaw roughness</text>",
        f'<text x="{(w - right + left) / 2}" y="{h - 12}" text-anchor="middle">Re = U D / ν</text>',
        f'<text transform="translate(20,{h / 2}) rotate(-90)" text-anchor="middle">'
        "Darcy friction factor f</text>",
    ]
    for e in (4, 5, 6, 7):
        out.append(f'<text x="{sx(10**e)}" y="{h - bottom + 18}" text-anchor="middle">1e{e}</text>')
    for f in (0.006, 0.01, 0.02, 0.04, 0.06):
        out.append(f'<text x="{left - 8}" y="{sy(f) + 4}" text-anchor="end">{f:g}</text>')
    for ks, col in zip(ROUGHNESS, colours):
        pts = [10 ** (lx0 + i * (lx1 - lx0) / 200) for i in range(201)]
        d = " ".join(f"{sx(r):.1f},{sy(colebrook(r, ks)):.1f}" for r in pts)
        out.append(f'<polyline points="{d}" fill="none" stroke="{col}" stroke-width="1.6"/>')
        label = "smooth" if ks == 0 else f"k_s/D = {ks:g}"
        out.append(f'<text x="{w - right + 10}" y="{sy(colebrook(3e7, ks)) + 4}" fill="{col}">{label}</text>')
    for r in rows:
        col = colours[ROUGHNESS.index(r["ks"])]
        if r["f"] is not None:
            out.append(
                f'<circle cx="{sx(r["re"]):.1f}" cy="{sy(r["f"]):.1f}" r="4.5" fill="none" '
                f'stroke="{col}" stroke-width="1.8"/>'
            )
    out.append(f'<text x="{left + 12}" y="{h - bottom - 12}" fill="#444">lines: Colebrook; circles: runs</text>')
    out.append("</svg>")
    dst.write_text("\n".join(out) + "\n")


def main(argv):
    out = Path(argv[0]) if argv and not argv[0].startswith("--") else Path("output/moody_study")
    cs = float(argv[argv.index("--cs") + 1]) if "--cs" in argv else 0.5
    jobs = int(argv[argv.index("--jobs") + 1]) if "--jobs" in argv else 1
    base = Case.load(ROOT / "examples/cases/flow/pipe_rough_keps.toml")
    base["boundary.ymax.roughness_constant"] = cs
    cases = [(re, ks) for re in CELLS for ks in ROUGHNESS]
    variants = listing(
        {
            "fluid.viscosity": 1.0 / re,
            "boundary.ymax.roughness": ks,
            "mesh.generate.cells": [2, cells(re, ks), 1],
            "solver.max_iterations": 6000,
        }
        for re, ks in cases
    )
    rows = Study(base, variants, out).run(jobs=jobs)
    results = []
    for (re, ks), r in zip(cases, rows):
        force = r.get("periodic.0.driving_force")
        f = 2.0 * force if force is not None else None
        ref = colebrook(re, ks)
        results.append({"re": re, "ks": ks, "f": f})
        dev = f"{f / ref - 1:+.1%}" if f is not None else "-"
        print(f"Re {re:8.0e}  k_s/D {ks:6.0e}  {r['status']:<9}  f {f or float('nan'):.5f}  "
              f"Colebrook {ref:.5f}  ({dev})")
    print(f"table: {out / 'study.csv'}")
    if "--chart" in argv:
        dst = Path(argv[argv.index("--chart") + 1])
        chart(results, dst)
        print(f"wrote {dst}")


if __name__ == "__main__":
    main(sys.argv[1:])

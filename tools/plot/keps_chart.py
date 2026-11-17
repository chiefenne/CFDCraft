"""Gallery 71: the k-ε family on the flat plate and in rough pipes (SVG).

    python3 -m tools.plot.keps_chart OUT.svg PLATE_DIR STUDY_CSV

PLATE_DIR holds the runs of the y⁺ study of tmr_plate_keps.toml, one
output directory per run named `<treatment>_<y+>` (scalable_1,
compound_30, ...; recipe 71 has the loop). STUDY_CSV is the table of
examples/studies/moody_study.py. Left: C_f(0.97) against the first-cell
y⁺ with White's correlation and its ±5 % band. Right: the Moody chart,
Colebrook's curves and the runs.
"""

import csv
import json
import math
import sys
from pathlib import Path

WHITE = 0.455 / math.log(0.06 * 0.970084 / 2.0e-7) ** 2
ROUGH = [0.0, 1e-4, 1e-3, 1e-2]
COLOURS = ["#1f6fb4", "#2a9d5b", "#d9822b", "#b03a48"]


def colebrook(re, ks):
    f = 0.02
    for _ in range(100):
        f = (-2.0 * math.log10(ks / 3.7 + 2.51 / (re * math.sqrt(f)))) ** -2
    return f


def plate(root: Path):
    runs = {}
    for d in sorted(root.iterdir()):
        r = d / "result.json"
        if not r.exists() or "_" not in d.name:
            continue
        treatment = d.name.split("_")[0]
        p = json.loads(r.read_text())["solves"][0]["wall_probes"][0]
        runs.setdefault(treatment, []).append((p["y_plus"], p["skin_friction"]))
    return {k: sorted(v) for k, v in runs.items()}


def moody(path: Path):
    rows = []
    with path.open() as fh:
        for r in csv.DictReader(fh):
            force = r.get("periodic.0.driving_force")
            if r["status"] != "finished" or not force:
                continue
            rows.append((1.0 / float(r["fluid.viscosity"]), float(r["boundary.ymax.roughness"]), 2.0 * float(force)))
    return rows


def main(argv):
    dst, plate_dir, study = Path(argv[0]), Path(argv[1]), Path(argv[2])
    runs, pipes = plate(plate_dir), moody(study)
    w, h = 1080, 430
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
        'font-family="Helvetica, Arial, sans-serif" font-size="12">',
        f'<rect width="{w}" height="{h}" fill="white"/>',
    ]
    # left: C_f over y+
    l, t, pw, ph = 70, 40, 400, 320
    x0, x1, y0, y1 = math.log10(0.7), math.log10(150), 0.0022, 0.0048

    def sx(v):
        return l + (math.log10(v) - x0) / (x1 - x0) * pw

    def sy(v):
        return t + (1 - (v - y0) / (y1 - y0)) * ph

    out += [
        f'<rect x="{l}" y="{t}" width="{pw}" height="{ph}" fill="none" stroke="#444"/>',
        f'<text x="{l + pw / 2}" y="24" text-anchor="middle" font-weight="bold">'
        "TMR flat plate, standard k-ε: C_f at x = 0.97</text>",
        f'<text x="{l + pw / 2}" y="{t + ph + 40}" text-anchor="middle">first-cell y⁺</text>',
        f'<text transform="translate(18,{t + ph / 2}) rotate(-90)" text-anchor="middle">C_f</text>',
        f'<rect x="{l}" y="{sy(1.05 * WHITE):.1f}" width="{pw}" '
        f'height="{sy(0.95 * WHITE) - sy(1.05 * WHITE):.1f}" fill="#888" opacity="0.15"/>',
        f'<line x1="{l}" x2="{l + pw}" y1="{sy(WHITE):.1f}" y2="{sy(WHITE):.1f}" stroke="#444" stroke-dasharray="5,4"/>',
        f'<text x="{l + pw - 6}" y="{sy(1.05 * WHITE) - 5:.1f}" text-anchor="end" fill="#444">White ±5 %</text>',
    ]
    for v in (1, 3, 10, 30, 100):
        out.append(f'<text x="{sx(v):.1f}" y="{t + ph + 18}" text-anchor="middle">{v}</text>')
    for v in (0.0025, 0.003, 0.0035, 0.004, 0.0045):
        out.append(f'<text x="{l - 6}" y="{sy(v) + 4:.1f}" text-anchor="end">{v:g}</text>')
    for i, (name, pts) in enumerate(sorted(runs.items(), reverse=True)):
        col = COLOURS[i * 3 % 4]
        d = " ".join(f"{sx(a):.1f},{sy(min(b, y1)):.1f}" for a, b in pts)
        out.append(f'<polyline points="{d}" fill="none" stroke="{col}" stroke-width="1.6"/>')
        for a, b in pts:
            out.append(f'<circle cx="{sx(a):.1f}" cy="{sy(min(b, y1)):.1f}" r="4" fill="white" stroke="{col}" stroke-width="1.8"/>')
        out.append(f'<text x="{l + 12}" y="{t + ph - 30 + 16 * i}" fill="{col}">{name}</text>')
    # right: Moody chart
    l = 590
    lx0, lx1, ly0, ly1 = math.log10(1e4), math.log10(3e7), math.log10(0.006), math.log10(0.06)

    def mx(re):
        return l + (math.log10(re) - lx0) / (lx1 - lx0) * pw

    def my(f):
        return t + (1 - (math.log10(f) - ly0) / (ly1 - ly0)) * ph

    out += [
        f'<rect x="{l}" y="{t}" width="{pw}" height="{ph}" fill="none" stroke="#444"/>',
        f'<text x="{l + pw / 2}" y="24" text-anchor="middle" font-weight="bold">'
        "Rough pipes: Darcy f against Colebrook</text>",
        f'<text x="{l + pw / 2}" y="{t + ph + 40}" text-anchor="middle">Re = U D/ν</text>',
        f'<text transform="translate({l - 52},{t + ph / 2}) rotate(-90)" text-anchor="middle">f</text>',
    ]
    for e in (4, 5, 6, 7):
        out.append(f'<text x="{mx(10**e):.1f}" y="{t + ph + 18}" text-anchor="middle">1e{e}</text>')
    for f in (0.006, 0.01, 0.02, 0.04):
        out.append(f'<text x="{l - 6}" y="{my(f) + 4:.1f}" text-anchor="end">{f:g}</text>')
    for ks, col in zip(ROUGH, COLOURS):
        pts = [10 ** (lx0 + i * (lx1 - lx0) / 200) for i in range(201)]
        d = " ".join(f"{mx(r):.1f},{my(colebrook(r, ks)):.1f}" for r in pts)
        out.append(f'<polyline points="{d}" fill="none" stroke="{col}" stroke-width="1.6"/>')
        label = "smooth" if ks == 0 else f"k_s/D {ks:g}"
        out.append(f'<text x="{l + pw + 6}" y="{my(colebrook(3e7, ks)) + 4:.1f}" fill="{col}">{label}</text>')
    for re, ks, f in pipes:
        col = COLOURS[ROUGH.index(ks)] if ks in ROUGH else "#444"
        out.append(f'<circle cx="{mx(re):.1f}" cy="{my(f):.1f}" r="4.5" fill="white" stroke="{col}" stroke-width="1.8"/>')
    out.append(f'<text x="{l + 10}" y="{t + ph - 10}" fill="#444">lines: Colebrook; circles: k-ε, C_s = 0.5</text>')
    out.append("</svg>")
    dst.write_text("\n".join(out) + "\n")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main(sys.argv[1:])

"""Validation of the turbulence models on the backward-facing step of
Driver & Seegmiller (1985), stage 6 part 6g, specification 9.5: the
reattachment length with its observed order and GCI on three meshes, the
pressure and skin-friction coefficients along the lower wall and the
velocity and shear-stress profiles at x/H = 1, 4, 6 and 10 against the
experiment.

    python3 examples/studies/backstep_validation.py [OUT] [--models sst,keps,...]
        [--levels 0,1,2] [--jobs 3] [--report-only] [--chart DIR]

The base case is examples/cases/flow/backstep_validation.toml (the layout
of NASA's Turbulence Modeling Resource, k-ε-ζ-f₀, coupled algorithm); the
study runs each model with `verification.levels` (default [0, 1, 2]:
11 040, 44 160 and 176 640 cells), so that the solver reports the GCI of
the reattachment length (the wall crossing `xr`, D67), and evaluates the
finest mesh against the experiment. Each model runs in OUT/<model>
(default output/backstep_validation); a finished run with an unchanged
case is not run again. `--report-only` only evaluates what is there.
`--chart DIR` writes backstep_validation_wall.svg and
backstep_validation_profiles.svg.

Models (key, label, wall treatment, algorithm):

    sst   k-ω SST                automatic   coupled
    keps  k-ε                    scalable    coupled
    rke   realizable k-ε         scalable    coupled
    zf    k-ε-ζ-f (2004)         compound    SIMPLEC (the coupled run does not converge)
    zf0   k-ε-ζ-f₀ (2023)        compound-2023 coupled
    eb    elliptic-blending k-ε  compound    SIMPLEC (the coupled run does not converge)

The coefficients use the reference speed of the experiment, the velocity
at x/H = −4 in the middle of the channel (y/H = 5); C_p is shifted to
zero at x/H = 37.5, the last station of the shifted data of the TMR (the
data are shifted the same way). The shear stress of the models is
−u'v' = ν_t ∂u/∂y on the lines (∂v/∂x left out), with the mean ν_t and
the least-squares slope of u over ±0.05 H, plotted from y/H = 0.05.

Experimental data: examples/studies/backstep/data (from the TMR, CC0;
see its README). Needs the package (pip install ./python) or the `foehn`
executable (FOEHN=target/release/foehn).
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
from itertools import pairwise
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "examples/studies/backstep/data"
BASE = ROOT / "examples/cases/flow/backstep_validation.toml"
X_EXP, X_EXP_ERR = 6.26, 0.10          # Driver & Seegmiller (1985), x_r/H
U_STATION = (-4.0, 5.0)                # reference speed: x/H, y/H
CP_ZERO = 37.5                         # C_p shifted to 0 here (the data too)
STATIONS = (1, 4, 6, 10)

COUPLED = {"solver.max_iterations": 800}
SIMPLEC = {
    "solver.algorithm": "simplec",
    "solver.max_iterations": 40000,
    "solver.momentum.relaxation": 0.5,
    # the pressure CG at its default limit of 100 iterations lets the
    # finer meshes diverge in the first iterations (RUN014, DIV006)
    "solver.pressure.max_iter": 200,
    # 0.5 while the bubble forms, 0.8 from iteration 100 (6d2)
    "schedule": [{"at_iteration": 100, "set": {"solver": {"momentum": {"relaxation": 0.8}}}}],
}
MODELS = {
    "sst": ("k-ω SST", {
        "turbulence.model": "k-omega-sst", "turbulence.wall_treatment": "automatic",
        "turbulence.ambient": {"k": 1.5e-4, "omega": 0.54}, **COUPLED}),
    "keps": ("k-ε", {
        "turbulence.model": "k-epsilon", "turbulence.wall_treatment": "scalable", **COUPLED}),
    "rke": ("realizable k-ε", {
        "turbulence.model": "realizable-k-epsilon", "turbulence.wall_treatment": "scalable",
        **COUPLED}),
    "zf": ("k-ε-ζ-f", {
        "turbulence.model": "k-epsilon-zeta-f", "turbulence.wall_treatment": "compound",
        # the positivity floor holds k in 44 cells, and the residual of k
        # stays at 1.4e-6 on the middle mesh (every other one below 1e-6)
        **SIMPLEC, "solver.tolerance": 2e-6}),
    "zf0": ("k-ε-ζ-f₀", dict(COUPLED)),
    "eb": ("elliptic-blending k-ε", {
        "turbulence.model": "elliptic-blending-k-epsilon", "turbulence.wall_treatment": "compound",
        **SIMPLEC}),
}
COLOURS = {"sst": "#1f6fb4", "keps": "#8c564b", "rke": "#d9822b", "zf": "#9467bd",
           "zf0": "#2a9d5b", "eb": "#0fa3b1"}

# Published x_r/H of the standard k-ε with wall functions on this case
# (selected in 6g, specification 9.5): Mizukami, Georgiadis & Cannon
# (1993), NASA Lewis, table 1: 5.47 and 5.22 from two solvers; the band
# ±5 % is taken about their mean.
KEPS_PUBLISHED = 0.5 * (5.47 + 5.22)


# -- reading ---------------------------------------------------------------


def _numbers(path: Path):
    """Tecplot-like data: zones of number rows; returns {zone: rows}."""
    zones: dict[str, list[list[float]]] = {}
    name = ""
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.lower().startswith("variables"):
            continue
        if s.lower().startswith("zone"):
            name = s.split("=", 1)[1].strip().strip('"')
            continue
        zones.setdefault(name, []).append([float(v) for v in s.split()])
    return zones


def experiment():
    """C_p (lower and upper wall, shifted), C_f (lower wall) and the
    profiles (y/H, u/U_ref, 1000 u'v'/U_ref²) by station."""
    cp = _numbers(DATA / "cp.expnew.dat")
    cf = _numbers(DATA / "cf.exp.dat")[""]
    prof = _numbers(DATA / "profiles.exp.dat")
    stations = {}
    for z, rows in prof.items():
        x = float(z.split("=")[1])
        stations[x] = [(r[1], r[2], r[6]) for r in rows]
    return {
        "cp_lower": [(r[0], r[1]) for r in cp["bottom wall"]],
        "cp_upper": [(r[0], r[1]) for r in cp["top wall"]],
        "cf": [(r[0], r[1], r[2]) for r in cf],
        "profiles": stations,
    }


def _csv(path: Path) -> list[dict[str, float]]:
    with open(path) as f:
        rows = [r for r in f if not r.startswith("#")]
    return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(rows)]


def interp(x: float, pts: list[tuple[float, float]]) -> float:
    """Linear interpolation in points sorted by their first coordinate."""
    if x <= pts[0][0]:
        return pts[0][1]
    for (x0, y0), (x1, y1) in pairwise(pts):
        if x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0) if x1 > x0 else y1
    return pts[-1][1]


def shear(rows: list[dict[str, float]], i: int, half: float = 0.05) -> float:
    """ν_t ∂u/∂y at row i of a line: the mean ν_t and the least-squares
    slope of u over the samples within ±half of its y (the line samples
    cell values, so a plain difference of neighbours is noisy)."""
    y0 = rows[i]["y"]
    w = [r for r in rows if abs(r["y"] - y0) <= half + 1e-12]
    if len(w) < 2:
        return 0.0
    my = sum(r["y"] for r in w) / len(w)
    mu = sum(r["u"] for r in w) / len(w)
    syy = sum((r["y"] - my) ** 2 for r in w)
    dudy = sum((r["y"] - my) * (r["u"] - mu) for r in w) / syy if syy > 0 else 0.0
    return sum(r["nu_t"] for r in w) / len(w) * dudy


def evaluate(d: Path, exp) -> dict | None:
    """The numbers of one model from its run directory (finest mesh)."""
    if not (d / "flow.json").exists():
        return None
    rep = json.loads((d / "flow.json").read_text())
    solves = rep.get("solves", [])
    xr = [next((c["position"] for c in s.get("wall_crossings", []) if c["name"] == "xr"), None)
          for s in solves]
    unc = next((e for e in rep.get("uncertainty", [])
                if e["quantity"] == "wall_crossings.xr.position"), None)
    top = _csv(d / "line_x-4.csv")
    uref = interp(U_STATION[1], [(r["y"], r["u"]) for r in top])
    q = 0.5 * uref * uref
    step = [r for r in _csv(d / "wall_step.csv") if abs(r["y"] - 1.0) < 1e-9]
    lower = sorted(step + _csv(d / "wall_ymin.csv"), key=lambda r: r["x"])
    upper = sorted(_csv(d / "wall_ymax.csv"), key=lambda r: r["x"])
    p0 = interp(CP_ZERO, [(r["x"], r["pressure"]) for r in lower])
    cp = [(r["x"], (r["pressure"] - p0) / q) for r in lower]
    cpu = [(r["x"], (r["pressure"] - p0) / q) for r in upper]
    cf = [(r["x"], r["tau_x"] / q) for r in lower]
    dcp = [interp(x, cp) - c for x, c in exp["cp_lower"]]
    dcpu = [interp(x, cpu) - c for x, c in exp["cp_upper"]]
    dcf = [interp(x, cf) - c for x, c, _ in exp["cf"]]
    profiles, du = {}, {}
    for s in STATIONS:
        rows = _csv(d / f"line_x{s}.csv")
        pts = [(r["y"], r["u"] / uref) for r in rows]
        uv = [(r["y"], 1000.0 * shear(rows, i) / (uref * uref)) for i, r in enumerate(rows)]
        profiles[s] = (pts, uv)
        ex = exp["profiles"][float(s)]
        du[s] = math.sqrt(sum((interp(y, pts) - u) ** 2 for y, u, _ in ex) / len(ex))
    return {
        "cells": [s["cells"] for s in solves],
        "iterations": [s["iterations"] for s in solves],
        "converged": [s.get("converged") for s in solves],
        "xr": xr,
        "uncertainty": unc,
        "uref": uref,
        "cf_m4": interp(-3.956, cf),
        "cp_rms": math.sqrt(sum(v * v for v in dcp) / len(dcp)),
        "cp_rms_upper": math.sqrt(sum(v * v for v in dcpu) / len(dcpu)),
        "cf_rms": math.sqrt(sum(v * v for v in dcf) / len(dcf)),
        "u_rms": du,
        "cp": cp, "cp_upper": cpu, "cf": cf, "profiles": profiles,
    }


# -- report ----------------------------------------------------------------


def _f(v, fmt, unit=""):
    """A number with a true minus sign, or – when there is none."""
    return "–" if v is None else format(v, fmt).replace("-", "−") + unit


# The bands of specification 9.5: x_r/H ±10 % of the experiment for SST,
# ζ-f (both forms) and EB (GCI_fine ≤ 2 %, p in [1, 3], C_p RMS ≤ 0.02),
# ±5 % of the published result for the standard k-ε (GCI_fine ≤ 2 %);
# the realizable k-ε is reported.
BANDS = {"sst": (X_EXP, 0.10, True), "zf": (X_EXP, 0.10, True), "zf0": (X_EXP, 0.10, True),
         "eb": (X_EXP, 0.10, True), "keps": (KEPS_PUBLISHED, 0.05, False)}


def verdict(key: str, r: dict) -> str:
    """The checks of specification 9.5 for one model, as text."""
    e = r["uncertainty"] or {}
    conv = (f"convergence {e['verdict']}"
            if e.get("verdict") and e["verdict"] != "monotonic" else "")
    if key not in BANDS:
        return ", ".join(v for v in ("reported", conv) if v)
    ref, band, strict = BANDS[key]
    fin = r["xr"][-1] if r["xr"] else None
    gci, p = e.get("gci"), e.get("order")
    out = []
    if fin is not None:
        out.append(f"x_r {'in' if abs(fin / ref - 1) <= band else 'outside'} ±{100 * band:.0f} %")
    if conv:
        out.append(conv)
    if gci is not None:
        out.append(f"GCI {'≤' if gci <= 0.02 else '>'} 2 %")
    if strict and p is not None:
        out.append(f"p {'in' if 1.0 <= p <= 3.0 else 'outside'} [1, 3]")
    if strict:
        out.append(f"C_p {'≤' if r['cp_rms'] <= 0.02 else '>'} 0.02")
    return ", ".join(out)


def _gci(g):
    """The GCI in per cent; below 0.01 % as a bound."""
    if g is None:
        return "–"
    return "< 0.01 %" if 100 * g < 0.01 else f"{100 * g:.2f} %"


def table(res: dict) -> str:
    """The Markdown tables of the results."""
    out = [
        ("| Model | x_r/H on the meshes | x_r/H finest | vs 6.26 | order p | GCI_fine | C_p RMS "
         "| C_f RMS ×10³ | C_f(−4) ×10³ | u/U_ref RMS at x/H 1, 4, 6, 10 | Specification 9.5 |"),
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for key, (label, _) in MODELS.items():
        r = res.get(key)
        if r is None:
            continue
        e = r["uncertainty"] or {}
        xr = [v for v in r["xr"] if v is not None]
        fin = r["xr"][-1] if r["xr"] else None
        out.append(
            f"| {label} | {', '.join(f'{v:.3f}' for v in xr)} | {_f(fin, '.3f')} | "
            f"{_f(fin and 100 * (fin / X_EXP - 1), '+.1f', ' %')} | {_f(e.get('order'), '.2f')} | "
            f"{_gci(e.get('gci'))} | {r['cp_rms']:.4f} | "
            f"{1e3 * r['cf_rms']:.3f} | {_f(1e3 * r['cf_m4'], '.3f')} | "
            + ", ".join(f"{r['u_rms'][s]:.3f}" for s in STATIONS)
            + f" | {verdict(key, r)} |"
        )
    out += [
        "",
        "| Model | Cells | Iterations | Converged |",
        "|---|---|---|---|",
    ]
    for key, (label, _) in MODELS.items():
        r = res.get(key)
        if r is None:
            continue
        conv = ["yes" if c else "no" for c in r["converged"]]
        out.append(f"| {label} | {', '.join(f'{c:,}'.replace(',', ' ') for c in r['cells'])} | "
                   f"{', '.join(str(i) for i in r['iterations'])} | {', '.join(conv)} |")
    return "\n".join(out) + "\n"


# -- charts ----------------------------------------------------------------


def _axes(out, x0, y0, w, h, title, xlab, ylab, xr, yr, ticks):
    def sx(v):
        return x0 + (v - xr[0]) / (xr[1] - xr[0]) * w

    def sy(v):
        return y0 + h - (v - yr[0]) / (yr[1] - yr[0]) * h

    out.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="none" stroke="#444"/>')
    out.append(f'<text x="{x0 + w / 2}" y="{y0 - 10}" text-anchor="middle" font-weight="bold">{title}</text>')
    out.append(f'<text x="{x0 + w / 2}" y="{y0 + h + 34}" text-anchor="middle">{xlab}</text>')
    out.append(f'<text transform="translate({x0 - 44},{y0 + h / 2}) rotate(-90)" '
               f'text-anchor="middle">{ylab}</text>')
    for v in ticks[0]:
        out.append(f'<line x1="{sx(v):.1f}" y1="{y0 + h}" x2="{sx(v):.1f}" y2="{y0 + h + 5}" stroke="#444"/>')
        out.append(f'<text x="{sx(v):.1f}" y="{y0 + h + 18}" text-anchor="middle">{v:g}</text>')
    for v in ticks[1]:
        out.append(f'<line x1="{x0 - 5}" y1="{sy(v):.1f}" x2="{x0}" y2="{sy(v):.1f}" stroke="#444"/>')
        out.append(f'<text x="{x0 - 8}" y="{sy(v) + 4:.1f}" text-anchor="end">{v:g}</text>')
    clip = f"c{x0}_{y0}"
    out.append(f'<clipPath id="{clip}"><rect x="{x0}" y="{y0}" width="{w}" height="{h}"/></clipPath>')

    def line(pts, colour, dash=""):
        d = " ".join(f"{sx(a):.1f},{sy(b):.1f}" for a, b in pts)
        extra = f' stroke-dasharray="{dash}"' if dash else ""
        out.append(f'<polyline points="{d}" fill="none" stroke="{colour}" stroke-width="1.6"'
                   f'{extra} clip-path="url(#{clip})"/>')

    def dots(pts, colour="#000"):
        for a, b in pts:
            if xr[0] <= a <= xr[1] and yr[0] <= b <= yr[1]:
                out.append(f'<circle cx="{sx(a):.1f}" cy="{sy(b):.1f}" r="3" fill="none" '
                           f'stroke="{colour}" stroke-width="1.3"/>')

    return line, dots, sx, sy


def _thin(pts, xr, step=0.0):
    """Points inside the x range, at least `step` apart."""
    out, last = [], None
    for a, b in pts:
        if xr[0] - 1 <= a <= xr[1] + 1 and (last is None or a - last >= step):
            out.append((a, b))
            last = a
    return out


def _legend(out, x, y, keys):
    for i, k in enumerate(keys):
        out.append(f'<line x1="{x}" y1="{y + 16 * i - 4}" x2="{x + 22}" y2="{y + 16 * i - 4}" '
                   f'stroke="{COLOURS[k]}" stroke-width="2"/>')
        out.append(f'<text x="{x + 28}" y="{y + 16 * i}">{MODELS[k][0]}</text>')
    i = len(keys)
    out.append(f'<circle cx="{x + 11}" cy="{y + 16 * i - 4}" r="3" fill="none" stroke="#000"/>')
    out.append(f'<text x="{x + 28}" y="{y + 16 * i}">Driver &amp; Seegmiller (1985)</text>')


def _svg(w, h):
    return [(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" '
             f'height="{h}" font-family="Helvetica, Arial, sans-serif" font-size="12">'),
            f'<rect width="{w}" height="{h}" fill="white"/>']


def charts(res: dict, exp, dst: Path) -> list[Path]:
    keys = [k for k in MODELS if res.get(k)]
    dst.mkdir(parents=True, exist_ok=True)
    w, h = 1220, 360
    out = _svg(w, h)
    xr = (-6.0, 36.0)
    line, dots, _, _ = _axes(out, 70, 40, 420, 270, "pressure coefficient, lower wall", "x/H",
                             'C<tspan baseline-shift="sub" font-size="9">p</tspan>', xr, (-0.25, 0.08), ([0, 10, 20, 30], [-0.2, -0.1, 0, 0.05]))
    for k in keys:
        line(_thin(res[k]["cp"], xr, 0.05), COLOURS[k])
    dots(exp["cp_lower"])
    line, dots, sx, sy = _axes(out, 580, 40, 420, 270, "skin friction, lower wall", "x/H",
                               'C<tspan baseline-shift="sub" font-size="9">f</tspan> × 10³', xr, (-1.5, 3.5), ([0, 10, 20, 30], [-1, 0, 1, 2, 3]))
    out.append(f'<line x1="{sx(xr[0])}" y1="{sy(0):.1f}" x2="{sx(xr[1])}" y2="{sy(0):.1f}" '
               'stroke="#999" stroke-dasharray="3,3"/>')
    for k in keys:
        line([(a, 1e3 * b) for a, b in _thin(res[k]["cf"], xr, 0.05)], COLOURS[k])
    dots([(x, 1e3 * c) for x, c, _ in exp["cf"]])
    _legend(out, 1020, 70, keys)
    out.append("</svg>")
    p1 = dst / "backstep_validation_wall.svg"
    p1.write_text("\n".join(out) + "\n")

    w, h = 1220, 700
    out = _svg(w, h)
    for i, s in enumerate(STATIONS):
        x0 = 70 + 240 * i
        line, dots, _, _ = _axes(out, x0, 40, 180, 260, f"x/H = {s}", 'u/U<tspan baseline-shift="sub" font-size="9">ref</tspan>',
                                 "y/H" if i == 0 else "", (-0.3, 1.1), (0, 3.0),
                                 ([0, 0.5, 1], [0, 1, 2, 3]))
        for k in keys:
            line([(u, y) for y, u in res[k]["profiles"][s][0] if y <= 3.05], COLOURS[k])
        dots([(u, y) for y, u, _ in exp["profiles"][float(s)]])
        line, dots, _, _ = _axes(out, x0, 380, 180, 260, f"x/H = {s}", "−1000 u'v'/U" '<tspan baseline-shift="sub" font-size="9">ref</tspan>²',
                                 "y/H" if i == 0 else "", (-2, 14), (0, 3.0),
                                 ([0, 5, 10], [0, 1, 2, 3]))
        for k in keys:
            line([(uv, y) for y, uv in res[k]["profiles"][s][1] if 0.05 <= y <= 3.05], COLOURS[k])
        dots([(-uv, y) for y, _, uv in exp["profiles"][float(s)]])
    _legend(out, 1010, 70, keys)
    out.append("</svg>")
    p2 = dst / "backstep_validation_profiles.svg"
    p2.write_text("\n".join(out) + "\n")
    return [p1, p2]


# -- main ------------------------------------------------------------------


def _arg(argv, name, default):
    return argv[argv.index(name) + 1] if name in argv else default


def main(argv):
    out = Path(argv[0]) if argv and not argv[0].startswith("--") else Path("output/backstep_validation")
    models = _arg(argv, "--models", ",".join(MODELS)).split(",")
    levels = [int(v) for v in _arg(argv, "--levels", "0,1,2").split(",")]
    jobs = int(_arg(argv, "--jobs", "1"))
    unknown = [m for m in models if m not in MODELS]
    if unknown:
        sys.exit(f"unknown models {unknown}; known: {', '.join(MODELS)}")
    if "--report-only" not in argv:
        from foehncraft import Case
        from foehncraft.study import Study, listing

        base = Case.load(BASE)
        threads = max(1, (os.cpu_count() or 1) // max(jobs, 1))
        variants = listing({**MODELS[m][1], "verification.levels": levels, "run.threads": threads}
                           for m in models)
        Study(base, variants, out, name=lambda i, v: models[i]).run(jobs=jobs)
    exp = experiment()
    res = {m: evaluate(out / m, exp) for m in MODELS}
    text = table(res)
    print(text)
    (out / "summary.md").write_text(text)
    keep = {m: {k: v for k, v in r.items() if k not in ("cp", "cp_upper", "cf", "profiles")}
            for m, r in res.items() if r}
    (out / "summary.json").write_text(json.dumps(keep, indent=1) + "\n")
    if "--chart" in argv:
        for p in charts(res, exp, Path(_arg(argv, "--chart", "."))):
            print(f"wrote {p}")


if __name__ == "__main__":
    main(sys.argv[1:])

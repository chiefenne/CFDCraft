"""Validation of the turbulence models on the round jet from a pipe that
impinges on a heated plate at H/D = 2, Re = 23 750, stage 6 part 6g,
specification 9.6: the stagnation Nusselt number Nu₀ with its observed
order and GCI on three meshes, Nu(r/D) out to r/D = 4 with its secondary
maximum, the stagnation-point anomaly with and without the production
limiters, and the wall-jet velocity profiles at r/D = 0.5, 1, 2, 2.5
and 3.

    python3 examples/studies/impinging_jet.py [OUT] [--models zf0,zf,...]
        [--levels 0,1,2] [--jobs 2] [--no-limiters] [--report-only]
        [--chart DIR]

Each model first runs the developed pipe flow of
examples/cases/flow/jet_pipe_profile.toml (OUT/pipe_<model>), whose line
`profile` becomes the inlet table of the jet (OUT/table_<model>.csv,
`foehncraft.profiles`); the jet, examples/cases/flow/impinging_jet.toml,
then runs with `verification.levels` (default [0, 1, 2]: 2 600, 10 400
and 41 600 cells), so that the solver reports the GCI of θ_w at the wall
probe `stag` (the face next to the axis; the GCI of Nu₀ = q D/(Γ θ_w) is
the same to first order). The finest mesh is evaluated. The limiter cases
run on the mesh of the case file (10 400 cells) in OUT/limiter: the
standard k-ε with `none`, `kato-launder` and its default `durbin`, SST
with `none`, `menter` (its default), `kato-launder` and `durbin`, and
k-ε-ζ-f₀ as the reference of the excess k. A finished run with an
unchanged case is not run again; `--report-only` only evaluates what is
there. `--chart DIR` writes impinging_jet_nusselt.svg and
impinging_jet_profiles.svg.

Models (key, label, wall treatment, algorithm):

    zf0   k-ε-ζ-f₀ (2023)          compound-2023  coupled
    zf    k-ε-ζ-f (2004)           compound       SIMPLEC (the coupled run does not converge)
    eb    elliptic-blending k-ε    compound       SIMPLEC (the coupled run does not converge)
    sst   k-ω SST, Kato–Launder    automatic      coupled
    keps  k-ε                      scalable       coupled
    rke   realizable k-ε           scalable       coupled (tolerance 1e-4: the residuals stall)

The excess k of the anomaly is the largest k on the axis within
x/D ≤ 0.25 of the plate above that of k-ε-ζ-f₀ on the same mesh; a
limiter passes when it removes at least 70 % of the excess of the run
without it.

Experimental data: examples/studies/impinging_jet/data (see its README).
Without them the study reports the models only. Needs the package
(pip install ./python) or the `foehn` executable
(FOEHN=target/release/foehn).
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
PIPE = ROOT / "examples/cases/flow/jet_pipe_profile.toml"
JET = ROOT / "examples/cases/flow/impinging_jet.toml"
DATA = ROOT / "examples/studies/impinging_jet/data"
NU = 4.2105263157894736e-5     # ν, Re = U D/ν = 23 750 (U = 1 m/s, D = 1 m)
GAMMA = NU / 0.71              # Γ = ν/Pr
Q = 1e-3                       # wall flux q/(ρ c_p), K m/s
STATIONS = ("0.5", "1", "2", "2.5", "3")
AXIS = 0.25                    # x/D of the stagnation region (excess k)

SIMPLEC = {
    "solver.algorithm": "simplec",
    "solver.max_iterations": 20000,
    "solver.momentum.relaxation": 0.5,
    "solver.pressure.max_iter": 200,
    # 0.5 while the jet turns, 0.8 from iteration 100 (as on the step)
    "schedule": [{"at_iteration": 100, "set": {"solver": {"momentum": {"relaxation": 0.8}}}}],
}
# (label, turbulence keys of both runs, solver keys of the jet)
MODELS = {
    "zf0": ("k-ε-ζ-f₀", {"turbulence.model": "k-epsilon-zeta-f0",
                         "turbulence.wall_treatment": "compound-2023"}, {}),
    "zf": ("k-ε-ζ-f", {"turbulence.model": "k-epsilon-zeta-f",
                       "turbulence.wall_treatment": "compound"}, dict(SIMPLEC)),
    # the floor of k holds a few cells, its residual stays at 4.5e-6 on
    # the middle mesh (every other one below 1e-7)
    "eb": ("elliptic-blending k-ε", {"turbulence.model": "elliptic-blending-k-epsilon",
                                     "turbulence.wall_treatment": "compound"},
           {**SIMPLEC, "solver.tolerance": 1e-5}),
    "sst": ("k-ω SST, Kato–Launder", {"turbulence.model": "k-omega-sst",
                                      "turbulence.wall_treatment": "automatic",
                                      "turbulence.production_limiter": "kato-launder"}, {}),
    "keps": ("k-ε", {"turbulence.model": "k-epsilon",
                     "turbulence.wall_treatment": "scalable"}, {}),
    # the residuals stall at 5e-5 (coupled) or 1e-4 (SIMPLEC) with a
    # steady field
    "rke": ("realizable k-ε", {"turbulence.model": "realizable-k-epsilon",
                               "turbulence.wall_treatment": "scalable"},
            {"solver.tolerance": 1e-4}),
}
COLOURS = {"sst": "#1f6fb4", "keps": "#8c564b", "rke": "#d9822b", "zf": "#9467bd",
           "zf0": "#2a9d5b", "eb": "#0fa3b1"}
# (model, production limiter); None: the model's default
LIMITERS = [("zf0", None), ("keps", "none"), ("keps", "kato-launder"), ("keps", "durbin"),
            ("sst", "none"), ("sst", "menter"), ("sst", "kato-launder"), ("sst", "durbin")]


# -- reading ---------------------------------------------------------------


def _csv(path: Path) -> list[dict[str, float]]:
    with open(path) as f:
        rows = [r for r in f if not r.startswith("#")]
    return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(rows)]


def experiment() -> dict:
    """The data of examples/studies/impinging_jet/data, where present:
    `nusselt` [(r/D, Nu)] and `velocity` {station: [(x/D, u_r/U)]}."""
    out: dict = {"nusselt": None, "velocity": {}}
    p = DATA / "nusselt.csv"
    if p.exists():
        out["nusselt"] = sorted((r["r_D"], r["Nu"]) for r in _csv(p))
    for s in STATIONS:
        p = DATA / f"velocity_r{s}.csv"
        if p.exists():
            out["velocity"][s] = sorted((r["x_D"], r["u_U"]) for r in _csv(p))
    return out


def interp(x: float, pts: list[tuple[float, float]]) -> float:
    """Linear interpolation in points sorted by their first coordinate."""
    if x <= pts[0][0]:
        return pts[0][1]
    for (x0, y0), (x1, y1) in pairwise(pts):
        if x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0) if x1 > x0 else y1
    return pts[-1][1]


def nusselt(d: Path) -> list[tuple[float, float]]:
    """Nu(r/D) = q D/(Γ θ_w) along the plate of the run in `d`."""
    return sorted((r["y"], Q / (GAMMA * r["theta_wall"])) for r in _csv(d / "wall_xmin.csv"))


def secondary(nu: list[tuple[float, float]]) -> tuple | None:
    """The secondary maximum (r/D, Nu) for 0.8 ≤ r/D ≤ 4: the largest
    local maximum at least 1 % above the lowest Nu before it (from
    r/D = 0.8 on), or None."""
    lo = math.inf
    best = None
    for (_, a), (r, b), (_, c) in zip(nu, nu[1:], nu[2:]):
        if r < 0.8 or r > 4.0:
            continue
        lo = min(lo, b)
        if b > a and b >= c and b >= 1.01 * lo and (best is None or b > best[1]):
            best = (r, b)
    return best


def axis_k(d: Path) -> float:
    """The largest k on the axis within x/D ≤ AXIS of the plate."""
    return max(r["k"] for r in _csv(d / "line_axis.csv") if r["s"] <= AXIS)


def wall_jet(d: Path, s: str) -> tuple[list, float, float]:
    """The radial velocity over the distance from the plate at r/D = s,
    its maximum and the distance of the maximum."""
    rows = _csv(d / f"line_r{s}.csv")
    pts = [(r["s"], r["v"]) for r in rows]
    um, xm = max((v, x) for x, v in pts)
    return pts, um, xm


def probe(solve: dict) -> float | None:
    """Nu₀ from θ_w of the wall probe `stag`."""
    for p in solve.get("wall_probes", []):
        if p["name"] == "stag":
            for s in p.get("scalars", []):
                if s["name"] == "theta":
                    return Q / (GAMMA * s["wall"])
    return None


def evaluate(d: Path, exp: dict) -> dict | None:
    """The numbers of one model from its run directory (finest mesh)."""
    if not (d / "flow.json").exists():
        return None
    rep = json.loads((d / "flow.json").read_text())
    solves = rep.get("solves", [])
    unc = next((e for e in rep.get("uncertainty", [])
                if e["quantity"] == "wall_probes.stag.theta_wall"), None)
    nu = nusselt(d)
    out = {
        "cells": [s["cells"] for s in solves],
        "iterations": [s["iterations"] for s in solves],
        "converged": [s.get("converged") for s in solves],
        "nu0": [probe(s) for s in solves],
        "uncertainty": unc,
        "secondary": secondary(nu),
        "nu_at": {r: interp(r, nu) for r in (0.5, 1.0, 2.0, 3.0, 4.0)},
        "k_axis": axis_k(d),
        "nu": nu,
        "jet": {s: wall_jet(d, s) for s in STATIONS},
        "nu_rms": None,
    }
    if exp["nusselt"]:
        pts = [(r, n) for r, n in exp["nusselt"] if r <= 4.0]
        out["nu_rms"] = math.sqrt(sum((interp(r, nu) / n - 1.0) ** 2 for r, n in pts) / len(pts))
    return out


def evaluate_limiters(root: Path) -> dict:
    """Nu₀ and the excess k of the limiter cases."""
    out = {}
    for m, lim in LIMITERS:
        d = root / _lim_name(m, lim)
        if (d / "flow.json").exists():
            rep = json.loads((d / "flow.json").read_text())
            out[(m, lim)] = {"nu0": probe(rep["solves"][-1]), "k": axis_k(d),
                             "converged": rep["solves"][-1].get("converged")}
    ref = out.get(("zf0", None))
    for m, lims in (("keps", ("kato-launder", "durbin")),
                    ("sst", ("menter", "kato-launder", "durbin"))):
        none = out.get((m, "none"))
        for lim in lims:
            r = out.get((m, lim))
            if r and none and ref and none["k"] > ref["k"]:
                r["removed"] = (none["k"] - r["k"]) / (none["k"] - ref["k"])
    return out


def _lim_name(m: str, lim: str | None) -> str:
    return m if lim is None else f"{m}_{lim}"


# -- report ----------------------------------------------------------------


def _f(v, fmt, unit=""):
    """A number with a true minus sign, or – when there is none."""
    return "–" if v is None else format(v, fmt).replace("-", "−") + unit


def _gci(g):
    if g is None:
        return "–"
    return "< 0.01 %" if 100 * g < 0.01 else f"{100 * g:.2f} %"


# The bands of specification 9.6 (Nu₀ band, secondary maximum checked):
# ζ-f (both forms) and EB ±15 %, SST with a limiter ±20 %, GCI_fine ≤ 3 %;
# Nu(r/D) RMS ≤ 15 % for these; the k-ε family is reported.
BANDS = {"zf0": (0.15, True), "zf": (0.15, True), "eb": (0.15, True), "sst": (0.20, False)}


def verdict(key: str, r: dict, exp: dict) -> str:
    e = r["uncertainty"] or {}
    out = []
    if e.get("verdict") and e["verdict"] != "monotonic":
        out.append(f"convergence {e['verdict']}")
    if key not in BANDS:
        return ", ".join(out + ["reported"])
    band, second = BANDS[key]
    if e.get("gci") is not None:
        out.append(f"GCI {'≤' if e['gci'] <= 0.03 else '>'} 3 %")
    if exp["nusselt"]:
        n0 = exp["nusselt"][0][1]
        fin = r["nu0"][-1]
        if fin is not None:
            out.append(f"Nu₀ {'in' if abs(fin / n0 - 1) <= band else 'outside'} "
                       f"±{100 * band:.0f} %")
        if second:
            es = secondary(exp["nusselt"])
            s = r["secondary"]
            if s is None:
                out.append("no secondary maximum")
            elif es is not None:
                out.append(f"secondary maximum {'within' if abs(s[0] - es[0]) <= 0.5 else 'off'}"
                           " ±0.5 D")
        if r["nu_rms"] is not None:
            out.append(f"RMS {'≤' if r['nu_rms'] <= 0.15 else '>'} 15 %")
    else:
        out.append("no experimental data")
    return ", ".join(out)


def table(res: dict, lim: dict, exp: dict) -> str:
    """The Markdown tables of the results."""
    n0 = exp["nusselt"][0][1] if exp["nusselt"] else None
    out = [
        ("| Model | Nu₀ on the meshes | Nu₀ finest | vs experiment | order p | GCI_fine "
         "| secondary maximum r/D, Nu | Nu at r/D 1, 2, 3, 4 | Nu RMS 0–4 | Specification 9.6 |"),
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for key, (label, _, _) in MODELS.items():
        r = res.get(key)
        if r is None:
            continue
        e = r["uncertainty"] or {}
        fin = r["nu0"][-1]
        s = r["secondary"]
        out.append(
            f"| {label} | {', '.join(_f(v, '.1f') for v in r['nu0'])} | {_f(fin, '.1f')} | "
            f"{_f(n0 and fin and 100 * (fin / n0 - 1), '+.1f', ' %')} | "
            f"{_f(e.get('order'), '.2f')} | {_gci(e.get('gci'))} | "
            f"{'none' if s is None else f'{s[0]:.2f}, {s[1]:.1f}'} | "
            + ", ".join(f"{r['nu_at'][x]:.1f}" for x in (1.0, 2.0, 3.0, 4.0))
            + f" | {_f(r['nu_rms'] and 100 * r['nu_rms'], '.1f', ' %')} | {verdict(key, r, exp)} |"
        )
    out += ["", "| Model | Cells | Iterations | Converged |", "|---|---|---|---|"]
    for key, (label, _, _) in MODELS.items():
        r = res.get(key)
        if r is None:
            continue
        conv = ["yes" if c else "no" for c in r["converged"]]
        out.append(f"| {label} | {', '.join(f'{c:,}'.replace(',', ' ') for c in r['cells'])} | "
                   f"{', '.join(str(i) for i in r['iterations'])} | {', '.join(conv)} |")
    out += ["", "| Wall jet: u_r,max/U at x/D of the maximum | "
            + " | ".join(f"r/D = {s}" for s in STATIONS) + " |",
            "|---|" + "---|" * len(STATIONS)]
    for key, (label, _, _) in MODELS.items():
        r = res.get(key)
        if r is None:
            continue
        out.append(f"| {label} | " + " | ".join(
            f"{r['jet'][s][1]:.3f} at {r['jet'][s][2]:.3f}" for s in STATIONS) + " |")
    if lim:
        out += ["", "| Limiter case (10 400 cells) | Nu₀ | vs experiment | k_max, x/D ≤ 0.25 "
                "| excess k removed | Specification 9.6 |", "|---|---|---|---|---|---|"]
        for (m, lm), r in lim.items():
            label = MODELS[m][0].split(",")[0] + (f", `{lm}`" if lm else " (reference)")
            dev = n0 and r["nu0"] and r["nu0"] / n0 - 1
            check = ""
            if lm == "none" and dev is not None:
                check = f"anomaly {'reproduced' if dev > 0.5 else 'not reproduced'} (> +50 %)"
            elif "removed" in r:
                check = f"{'≥' if r['removed'] >= 0.7 else '<'} 70 % removed"
            out.append(f"| {label} | {_f(r['nu0'], '.1f')} | {_f(dev and 100 * dev, '+.1f', ' %')}"
                       f" | {r['k']:.4f} | {_f(r.get('removed') and 100 * r['removed'], '.0f', ' %')}"
                       f" | {check} |")
    return "\n".join(out) + "\n"


# -- charts ----------------------------------------------------------------


def _axes(out, x0, y0, w, h, title, xlab, ylab, xr, yr, ticks):
    def sx(v):
        return x0 + (v - xr[0]) / (xr[1] - xr[0]) * w

    def sy(v):
        return y0 + h - (v - yr[0]) / (yr[1] - yr[0]) * h

    out.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="none" stroke="#444"/>')
    out.append(f'<text x="{x0 + w / 2}" y="{y0 - 10}" text-anchor="middle" font-weight="bold">'
               f'{title}</text>')
    out.append(f'<text x="{x0 + w / 2}" y="{y0 + h + 34}" text-anchor="middle">{xlab}</text>')
    out.append(f'<text transform="translate({x0 - 44},{y0 + h / 2}) rotate(-90)" '
               f'text-anchor="middle">{ylab}</text>')
    for v in ticks[0]:
        out.append(f'<line x1="{sx(v):.1f}" y1="{y0 + h}" x2="{sx(v):.1f}" y2="{y0 + h + 5}" '
                   'stroke="#444"/>')
        out.append(f'<text x="{sx(v):.1f}" y="{y0 + h + 18}" text-anchor="middle">{v:g}</text>')
    for v in ticks[1]:
        out.append(f'<line x1="{x0 - 5}" y1="{sy(v):.1f}" x2="{x0}" y2="{sy(v):.1f}" '
                   'stroke="#444"/>')
        out.append(f'<text x="{x0 - 8}" y="{sy(v) + 4:.1f}" text-anchor="end">{v:g}</text>')
    clip = f"c{x0}_{y0}"
    out.append(f'<clipPath id="{clip}"><rect x="{x0}" y="{y0}" width="{w}" height="{h}"/>'
               '</clipPath>')

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

    return line, dots


def _svg(w, h):
    return [(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" '
             f'height="{h}" font-family="Helvetica, Arial, sans-serif" font-size="12">'),
            f'<rect width="{w}" height="{h}" fill="white"/>']


def _legend(out, x, y, entries, data):
    for i, (label, colour, dash) in enumerate(entries):
        extra = f' stroke-dasharray="{dash}"' if dash else ""
        out.append(f'<line x1="{x}" y1="{y + 16 * i - 4}" x2="{x + 22}" y2="{y + 16 * i - 4}" '
                   f'stroke="{colour}" stroke-width="2"{extra}/>')
        out.append(f'<text x="{x + 28}" y="{y + 16 * i}">{label}</text>')
    if data:
        i = len(entries)
        out.append(f'<circle cx="{x + 11}" cy="{y + 16 * i - 4}" r="3" fill="none" '
                   'stroke="#000"/>')
        out.append(f'<text x="{x + 28}" y="{y + 16 * i}">{data}</text>')


def charts(res: dict, lim_root: Path, exp: dict, dst: Path) -> list[Path]:
    keys = [k for k in MODELS if res.get(k)]
    dst.mkdir(parents=True, exist_ok=True)
    w, h = 1220, 380
    out = _svg(w, h)
    xr, yr, ticks = (0.0, 4.0), (0.0, 260.0), ([0, 1, 2, 3, 4], [0, 50, 100, 150, 200, 250])
    line, dots = _axes(out, 70, 40, 420, 270, "Nusselt number, models", "r/D", "Nu", xr, yr,
                       ticks)
    for k in keys:
        line([(r, n) for r, n in res[k]["nu"] if r <= 4.05], COLOURS[k])
    if exp["nusselt"]:
        dots(exp["nusselt"])
    line, dots = _axes(out, 580, 40, 420, 270, "stagnation-point anomaly (10 400 cells)", "r/D",
                       "Nu", xr, yr, ticks)
    styles = {"none": "", "menter": "2,3", "kato-launder": "7,4", "durbin": "1,3"}
    entries = []
    for m, lm in LIMITERS:
        d = lim_root / _lim_name(m, lm)
        if lm is None or not (d / "wall_xmin.csv").exists():
            continue
        line([(r, n) for r, n in nusselt(d) if r <= 4.05], COLOURS[m], styles[lm])
        entries.append((f"{MODELS[m][0].split(',')[0]}, {lm}", COLOURS[m], styles[lm]))
    if exp["nusselt"]:
        dots(exp["nusselt"])
    data = "Baughn et al. (1991)" if exp["nusselt"] else ""
    _legend(out, 1020, 60, [(MODELS[k][0], COLOURS[k], "") for k in keys], data)
    _legend(out, 1020, 80 + 16 * (len(keys) + 1), entries, "")
    out.append("</svg>")
    p1 = dst / "impinging_jet_nusselt.svg"
    p1.write_text("\n".join(out) + "\n")

    w, h = 1220, 380
    out = _svg(w, h)
    for i, s in enumerate(STATIONS):
        x0 = 70 + 190 * i
        line, dots = _axes(out, x0, 40, 150, 270, f"r/D = {s}", "u_r/U", "x/D" if i == 0 else "",
                           (-0.1, 1.2), (0.0, 0.3), ([0, 0.5, 1], [0, 0.1, 0.2, 0.3]))
        for k in keys:
            line([(v, x) for x, v in res[k]["jet"][s][0] if x <= 0.31], COLOURS[k])
        if s in exp["velocity"]:
            dots([(u, x) for x, u in exp["velocity"][s]])
    data = "Cooper et al. (1993)" if exp["velocity"] else ""
    _legend(out, 1030, 60, [(MODELS[k][0], COLOURS[k], "") for k in keys], data)
    out.append("</svg>")
    p2 = dst / "impinging_jet_profiles.svg"
    p2.write_text("\n".join(out) + "\n")
    return [p1, p2]


# -- main ------------------------------------------------------------------


def _arg(argv, name, default):
    return argv[argv.index(name) + 1] if name in argv else default


def main(argv):
    out = (Path(argv[0]) if argv and not argv[0].startswith("--")
           else Path("output/impinging_jet")).resolve()
    models = _arg(argv, "--models", ",".join(MODELS)).split(",")
    levels = [int(v) for v in _arg(argv, "--levels", "0,1,2").split(",")]
    jobs = int(_arg(argv, "--jobs", "1"))
    unknown = [m for m in models if m not in MODELS]
    if unknown:
        sys.exit(f"unknown models {unknown}; known: {', '.join(MODELS)}")
    lim_root = out / "limiter"
    if "--report-only" not in argv:
        from foehncraft import Case
        from foehncraft.profiles import inlet_table
        from foehncraft.study import Study, listing

        threads = max(1, (os.cpu_count() or 1) // max(jobs, 1))
        need = sorted(set(models) | ({m for m, _ in LIMITERS}
                                     if "--no-limiters" not in argv else set()))
        # the developed pipe flow of each model, then its inlet table
        Study(Case.load(PIPE), listing({**MODELS[m][1], "run.threads": threads} for m in need),
              out, name=lambda i, v: f"pipe_{need[i]}").run(jobs=jobs)
        for m in need:
            inlet_table(out / f"pipe_{m}" / "line_profile.csv", out / f"table_{m}.csv",
                        header=f"inlet table of the developed pipe flow, {MODELS[m][0]}")

        def jet(m, extra):
            return {**MODELS[m][1], **MODELS[m][2], "boundary.xmax.table":
                    str(out / f"table_{m}.csv"), "run.threads": threads, **extra}

        base = Case.load(JET)
        Study(base, listing(jet(m, {"verification.levels": levels}) for m in models), out,
              name=lambda i, v: models[i]).run(jobs=jobs)
        if "--no-limiters" not in argv:
            lim = [jet(m, {} if lm is None else {"turbulence.production_limiter": lm})
                   for m, lm in LIMITERS]
            Study(base, listing(lim), lim_root,
                  name=lambda i, v: _lim_name(*LIMITERS[i])).run(jobs=jobs)
    exp = experiment()
    res = {m: evaluate(out / m, exp) for m in MODELS}
    lim = evaluate_limiters(lim_root)
    text = table(res, lim, exp)
    print(text)
    (out / "summary.md").write_text(text)
    keep = {m: {k: v for k, v in r.items() if k not in ("nu", "jet")}
            for m, r in res.items() if r}
    keep["limiters"] = {_lim_name(*k): v for k, v in lim.items()}
    (out / "summary.json").write_text(json.dumps(keep, indent=1) + "\n")
    if "--chart" in argv:
        for p in charts(res, lim_root, exp, Path(_arg(argv, "--chart", "."))):
            print(f"wrote {p}")


if __name__ == "__main__":
    main(sys.argv[1:])

"""The time step from the time-integration error (D65) on the Taylor–Green
vortex: Δt over t for the error tolerances 1e-4 to 1e-6 and for the Courant number
alone, with the rejected steps marked, and the velocity error at the end
of each run. Prints one line per run and writes the gallery chart.

    python3 -m tools.plot.time_step_chart docs/figures/gallery/time_step_error.svg [FOEHN] [OUTDIR]

FOEHN defaults to target/release/foehn, OUTDIR to output/time_step_error.
The runs are those of examples/cases/flow/taylor_green_error.toml with the
monitor and field schedules moved to the end time, so that only the error
and the Courant number choose the steps.
"""
import csv
import json
import math
import os
import re
import subprocess
import sys

CASE = "examples/cases/flow/taylor_green_error.toml"
RUNS = [  # label, colour, extra --set
    ("tolerance 1e-4", "#1f6fb4", ["time.error_tolerance=1e-4"]),
    ("tolerance 1e-5", "#d9822b", ["time.error_tolerance=1e-5"]),
    ("tolerance 1e-6", "#2e9e5b", ["time.error_tolerance=1e-6"]),
]
REJECTED = re.compile(r"step (\d+) rejected: error ([0-9.]+) of the tolerance at Δt ([0-9.e+-]+) s")


def run(foehn, out, sets, case=CASE):
    args = [foehn, "--quiet", "flow", case, "--out", out, "--no-vtkhdf",
            "--set", "output.monitor_interval=0.4", "--set", "output.fields_interval=0.4"]
    for s in sets:
        args += ["--set", s]
    subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    with open(os.path.join(out, "result.json")) as f:
        res = json.load(f)
    steps = []
    with open(os.path.join(out, "history.csv")) as f:
        rows = [r for r in f if not r.startswith("#")]
    for r in csv.DictReader(rows):
        steps.append((float(r["time"]), float(r["dt"])))
    rejected = []
    t_before = {0: 0.0}
    for i, (t, dt) in enumerate(steps):
        t_before[i + 1] = t
    with open(os.path.join(out, "flow.log")) as f:
        for line in f:
            m = REJECTED.search(line)
            if m:
                n, dt = int(m.group(1)), float(m.group(3))
                rejected.append((t_before[n - 1] + dt, dt))
    return res, steps, rejected


def chart(series, dst):
    w, h, left, right, top, bottom = 680, 440, 80, 210, 50, 60
    ts = [t for _, _, s, _, _ in series for t, _ in s]
    ds = [d for _, _, s, r, _ in series for _, d in s + r]
    x0, x1 = 0.0, max(ts)
    y0, y1 = math.floor(math.log10(min(ds))), math.ceil(math.log10(max(ds)))
    sx = lambda v: left + (v - x0) / (x1 - x0) * (w - left - right)
    sy = lambda v: top + (1 - (math.log10(v) - y0) / (y1 - y0)) * (h - top - bottom)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
           'font-family="Helvetica, Arial, sans-serif" font-size="12">',
           f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
           f'<text x="{(left + w - right) / 2}" y="26" text-anchor="middle" font-size="15">'
           'Time step from the BDF2 error estimate, Taylor–Green vortex</text>',
           f'<rect x="{left}" y="{top}" width="{w - left - right}" height="{h - top - bottom}" fill="none" stroke="#444"/>']
    for k in range(0, 5):
        x = 0.1 * k
        out.append(f'<line x1="{sx(x):.1f}" y1="{top}" x2="{sx(x):.1f}" y2="{h - bottom}" stroke="#e5e5e5"/>')
        out.append(f'<text x="{sx(x):.1f}" y="{h - bottom + 18}" text-anchor="middle">{x:.1f}</text>')
    for k in range(y0, y1 + 1):
        out.append(f'<line x1="{left}" y1="{sy(10.0 ** k):.1f}" x2="{w - right}" y2="{sy(10.0 ** k):.1f}" stroke="#e5e5e5"/>')
        out.append(f'<text x="{left - 8}" y="{sy(10.0 ** k) + 4:.1f}" text-anchor="end">1e{k}</text>')
    for label, colour, steps, rejected, _ in series:
        line = " ".join(f"{sx(t):.1f},{sy(d):.1f}" for t, d in steps)
        out.append(f'<polyline points="{line}" fill="none" stroke="{colour}" stroke-width="1.8"/>')
        for t, d in rejected:
            x, y = sx(t), sy(d)
            out.append(f'<path d="M{x - 4:.1f},{y - 4:.1f} L{x + 4:.1f},{y + 4:.1f} M{x - 4:.1f},{y + 4:.1f} '
                       f'L{x + 4:.1f},{y - 4:.1f}" stroke="{colour}" stroke-width="1.6"/>')
    out.append(f'<text x="{(left + w - right) / 2}" y="{h - 16}" text-anchor="middle">time t (s)</text>')
    out.append(f'<text x="18" y="{(top + h - bottom) / 2}" text-anchor="middle" '
               f'transform="rotate(-90 18 {(top + h - bottom) / 2})">time step Δt (s)</text>')
    yy = top + 20
    for label, colour, steps, rejected, err in series:
        out.append(f'<line x1="{w - right + 10}" y1="{yy - 4}" x2="{w - right + 30}" y2="{yy - 4}" '
                   f'stroke="{colour}" stroke-width="1.8"/>')
        out.append(f'<text x="{w - right + 36}" y="{yy}">{label}</text>')
        out.append(f'<text x="{w - right + 36}" y="{yy + 15}">{len(steps)} steps, L2 error {err:.2e}</text>')
        yy += 40
    out.append(f'<path d="M{w - right + 16},{yy - 8} L{w - right + 24},{yy} M{w - right + 16},{yy} '
               f'L{w - right + 24},{yy - 8}" stroke="#444" stroke-width="1.6"/>')
    out.append(f'<text x="{w - right + 36}" y="{yy}">rejected step</text>')
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    dst = sys.argv[1]
    foehn = sys.argv[2] if len(sys.argv) > 2 else "target/release/foehn"
    outdir = sys.argv[3] if len(sys.argv) > 3 else "output/time_step_error"
    series = []
    for i, (label, colour, sets) in enumerate(RUNS):
        res, steps, rejected = run(foehn, os.path.join(outdir, f"run_{i}"), sets)
        s = res["solves"][-1]
        err = s["errors"]["l2"]
        series.append((label, colour, steps, rejected, err))
        print(f"{label}: {len(steps)} steps, {len(rejected)} rejected, Δt {s['transient']['dt_min']:.3e}–"
              f"{s['transient']['dt_max']:.3e} s, largest Courant number {s['transient']['cfl_max']:.2f}, "
              f"velocity error L2 {err:.3e}")
    # the case without the error tolerance: the Courant number and the
    # growth limit alone
    os.makedirs(outdir, exist_ok=True)
    plain = os.path.join(outdir, "courant.toml")
    with open(CASE) as f, open(plain, "w") as g:
        g.writelines(l for l in f if not l.startswith("error_tolerance"))
    res, steps, rejected = run(foehn, os.path.join(outdir, "courant"), [], plain)
    err = res["solves"][-1]["errors"]["l2"]
    series.append(("Courant number only", "#888888", steps, rejected, err))
    print(f"Courant number only: {len(steps)} steps, velocity error L2 {err:.3e}")
    chart(series, dst)
    print(f"chart {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

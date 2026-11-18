"""Channel heat-transfer DNS data against our heated channel (stage 6
part 6f, D135, specification 9.4): fetches the DNS tables of the JAXA DNS
database, recomputes their bulk velocity and Nusselt number and, with
`--run`, runs `examples/cases/flow/channel_heated_keps.toml` and compares.

    python examples/studies/channel_heat_dns.py [OUT] [--run]

The data are the work of H. Abe, H. Kawamura and Y. Matsuo (Re_τ 180 to
1020, Pr 0.71 and 0.025, uniform wall heat flux on both walls) and the
property of JAXA: https://jaxa-dns-database.jaxa.jp/channelflow.html says
"Rights for all information including the data, pictures and charts in
this database are reserved by JAXA" and asks for citation of Abe et al.
(2001, 2004, 2009, 2019). The files are therefore downloaded into OUT
(default `output/dns_jaxa`, which git ignores) and never committed; the
repository holds only the published Nusselt number that the test
`heated_channel_meets_the_dns_nusselt_number` compares with. Read the
terms on that page before using the data for anything else.

Normalisation of the tables: y⁺ = y u_τ/ν, U⁺, Θ⁺ = (T_w − T)/T_τ with
T_τ = q_w/(ρ c_p u_τ). From them U_b⁺ = (1/Re_τ) ∫ U⁺ dy⁺ over the half
channel, Θ_b⁺ = ∫ U⁺ Θ⁺ dy⁺ / ∫ U⁺ dy⁺ and Nu_2h = 2 Re_τ Pr/Θ_b⁺
(Nu_4h = 2 Nu_2h on the hydraulic diameter 4h). The first line of each
printed row compares this with the Nu of the file's header.

Findings (2026-11-18): the integrated Nu agrees with the header's to the
precision of the tables for the cases listed; the header of
ch1020pr0025.dat gives Re_bulk = 41441 and Nu = 278.71, the Re_bulk of
the Pr = 0.71 case and a Nu that no profile supports (the script prints
the integrated value instead).
"""

import re
import subprocess
import sys
import urllib.request
from pathlib import Path

BASE = "https://jaxa-dns-database.jaxa.jp/channelflow/"
ROOT = Path(__file__).resolve().parents[2]
CASES = [
    (180, "071"),
    (395, "071"),
    (640, "071"),
    (1020, "071"),
    (180, "0025"),
    (395, "0025"),
    (640, "0025"),
    (1020, "0025"),
]


def fetch(out: Path, name: str) -> Path:
    path = out / name
    if not path.exists():
        with urllib.request.urlopen(BASE + name, timeout=60) as r:
            path.write_bytes(r.read())
    return path


def first_table(path: Path):
    """The rows (floats) of the first table of the file."""
    rows, on = [], False
    for line in path.read_text(errors="replace").splitlines():
        if not on:
            on = re.match(r"\s*j\s+y\+", line) is not None
            continue
        p = line.split()
        if not p:
            break
        try:
            rows.append([float(x) for x in p])
        except ValueError:
            break
    return rows


def header_value(path: Path, key: str) -> float:
    m = re.search(
        rf"^\s+{key}\s*=\s*([0-9.eE+-]+)\s*$",
        path.read_text(errors="replace"),
        re.MULTILINE,
    )
    return float(m.group(1))


def trapz(x, f):
    return sum(0.5 * (f[i] + f[i + 1]) * (x[i + 1] - x[i]) for i in range(len(x) - 1))


def integrals(re_tau, vel, thermal):
    y = [0.0] + [r[1] for r in vel]
    u = [0.0] + [r[2] for r in vel]
    th = [0.0] + [r[2] for r in thermal]
    n = min(len(y), len(th))
    y, u, th = y[:n], u[:n], th[:n]
    if y[-1] < re_tau:  # the table stops just short of the centre line
        y.append(re_tau)
        u.append(u[-1])
        th.append(th[-1])
    ub = trapz(y, u) / re_tau
    thb = trapz(y, [a * b for a, b in zip(u, th)]) / trapz(y, u)
    return ub, thb


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    out = Path(args[0]) if args else ROOT / "output" / "dns_jaxa"
    out.mkdir(parents=True, exist_ok=True)
    print(
        f"{'Re_tau':>6} {'Pr':>6} {'U_b+':>7} {'Re_b':>7} {'Re_b file':>9} "
        f"{'Th_b+':>7} {'Nu_2h':>7} {'Nu file':>8} {'Nu_4h':>7}"
    )
    nu_by = {}
    for re_tau, pr in CASES:
        v = fetch(out, f"ch{re_tau}.dat")
        t = fetch(out, f"ch{re_tau}pr{pr}.dat")
        fetch(out, f"nut_at_Prt_{re_tau}pr{pr}.dat")
        prn = header_value(t, "Pr")
        ub, thb = integrals(re_tau, first_table(v), first_table(t))
        nu = 2 * re_tau * prn / thb
        # the header's Nu where it agrees with the profiles (all but one case)
        nu_file = header_value(t, "Nu")
        nu_by[(re_tau, pr)] = nu_file if abs(nu_file / nu - 1) < 0.05 else nu
        print(
            f"{re_tau:6d} {prn:6.3f} {ub:7.3f} {2 * re_tau * ub:7.0f} "
            f"{header_value(t, 'Re_bulk'):9.0f} {thb:7.3f} {nu:7.2f} "
            f"{header_value(t, 'Nu'):8.2f} {2 * nu:7.2f}"
        )
    if "--run" in sys.argv:
        exe = ROOT / "target" / "release" / "foehn"
        base = (
            ROOT / "examples" / "cases" / "flow" / "channel_heated_keps.toml"
        ).read_text()
        print(
            f"\n{'Re_tau':>6} {'Pr':>6} {'cells':>5} {'Pr_t':>14} {'y+_P':>6} "
            f"{'Re_tau run':>10} {'Nu_2h':>7} {'DNS':>7} {'diff':>7}"
        )
        for re_tau, pr in CASES:
            t = out / f"ch{re_tau}pr{pr}.dat"
            re_b = header_value(t, "Re_bulk")
            ref = nu_by[(re_tau, pr)]
            n = CELLS[re_tau]
            for model in ("constant", "kays-crawford"):
                text = variant(base, re_tau, header_value(t, "Pr"), re_b, n, model)
                name = f"ch{re_tau}pr{pr}_{model}"
                case = out / f"{name}.toml"
                case.write_text(text)
                run_dir = out / "run" / name
                subprocess.run(
                    [str(exe), "flow", str(case), "--out", str(run_dir)],
                    check=True,
                    capture_output=True,
                )
                nu, tau = nusselt(run_dir, re_tau, header_value(t, "Pr"), re_b)
                print(
                    f"{re_tau:6d} {header_value(t, 'Pr'):6.3f} {n:5d} {model:>14} "
                    f"{tau / n:6.1f} {tau:10.1f} {nu:7.2f} {ref:7.2f} "
                    f"{100 * (nu / ref - 1):+6.1f}%"
                )


# Uniform meshes across the channel (2h = 2 m) with the first cell centre at
# y+ = Re_tau/N of about 30, the claimed band of the scalable wall functions.
CELLS = {180: 6, 395: 12, 640: 20, 1020: 32}


def variant(base: str, re_tau: int, pr: float, re_b: float, n: int, model: str) -> str:
    """The case file of `channel_heated_keps.toml` at another Re_tau and Pr:
    nu = 1/Re_tau (u_tau = 1 m/s, h = 1 m), Gamma = nu/Pr, the mass flow of
    the DNS bulk Reynolds number Re_b = 2 h U_b/nu on the 0.2 m^2 section and
    N uniform cells across."""
    ub = re_b / (2 * re_tau)
    rep = {
        "cells = [2, 12, 1]": f"cells = [2, {n}, 1]",
        "viscosity = 0.0025316455696202532   # 1/395": f"viscosity = {1 / re_tau!r}",
        "diffusivity = 0.0035661218424962853   # ν/Pr, Pr = 0.71": f"diffusivity = {1 / (re_tau * pr)!r}",
        "velocity = [15.0, 0.0, 0.0]": f"velocity = [{ub:.4f}, 0.0, 0.0]",
        "mass_flow = 3.51": f"mass_flow = {0.2 * ub!r}",
        "turbulent_prandtl = 0.85\n": f'turbulent_prandtl = 0.85\nturbulent_prandtl_model = "{model}"\n',
    }
    for a, b in rep.items():
        assert a in base, a
        base = base.replace(a, b)
    return base


def nusselt(run_dir: Path, re_tau: int, pr: float, re_b: float):
    """Nu_2h = 2 h q/(Gamma (T_w - T_b)) from the result file, with T_b = 0 at
    xmin and T_w the wall mean there (minus Delta/2, Delta = Q/m), and the
    Re_tau of the run from the driving force F = tau_w/h."""
    import json

    res = json.loads((run_dir / "result.json").read_text())
    sol = res["solves"][-1] if "solves" in res else res
    b = {x["boundary"]: x for x in sol["scalars"][0]["boundaries"]}
    q, mdot = 1e-3, 0.2 * re_b / (2 * re_tau)
    delta = q * 0.02 / mdot
    tw = 0.5 * (b["ymin"]["mean"] + b["ymax"]["mean"]) - delta / 2
    nu = 2 * q / (1 / (re_tau * pr) * (tw - b["xmin"]["mixed_mean"]))
    force = sol["periodic"][0]["driving_force"]  # F = tau_w/h, rho = 1, h = 1
    tau = re_tau * abs(force) ** 0.5
    return nu, tau


if __name__ == "__main__":
    main()

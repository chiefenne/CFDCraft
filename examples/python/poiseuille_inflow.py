"""A case built in Python and run in process (D61): Poiseuille flow with
the inflow profile written as a Python expression, a monitor that prints
the residuals, and a stop criterion of our own.

    python examples/python/poiseuille_inflow.py [OUT]

The base case is examples/cases/flow/poiseuille_hex.toml on its coarse
mesh. Python replaces `profile = "analytic"` by the parabola built from
expression objects, writes the case file into the output directory (the
run can be repeated from it with `foehn flow OUT/poiseuille_py.toml`) and
runs it. The callback sees every outer iteration and stops the run once
the continuity residual is below 1e-8, before the solver tolerance of the
case file (1e-9): the run ends as "stopped" (exit code 3) with its restart
file and the final output, as with a `.stop` file.
"""

import sys
from pathlib import Path

from foehncraft import Case, run
from foehncraft.expr import constant, y

ROOT = Path(__file__).resolve().parents[2]


def main(argv):
    out = Path(argv[0]) if argv else Path("output/poiseuille_py")
    case = Case.load(ROOT / "examples/cases/flow/poiseuille_hex.toml")
    case.name = "poiseuille_py"
    case["verification"] = {"max_error": 2e-3}
    case["constants"] = {"u_mean": 1.0, "h": 1.0}
    um, h = constant("u_mean"), constant("h")
    profile = 6 * um * (y / h) * (1 - y / h)
    case["boundary.xmin.velocity"] = [profile, 0.0, 0.0]
    del case["boundary.xmin.profile"]
    print(f"inflow: {profile}")

    def monitor(row):
        if row["iteration"] % 50 == 0:
            print(f"  iteration {row['iteration']:4d}  continuity {row['continuity']:.3e}")
        return row["continuity"] < 1e-8  # True stops the run

    res = run(case, out=out, on_iteration=monitor)
    print(f"{res.status} (exit code {res.exit_code}) after {res.value('iterations')} iterations")
    print(f"velocity error L2 {res.value('errors.l2'):.4e}")
    print(f"case file: {out / 'poiseuille_py.toml'}")
    return 0 if res.exit_code in (0, 3) else res.exit_code


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

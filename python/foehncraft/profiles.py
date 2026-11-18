"""Inlet tables from sampled profiles (stage 6 specification 7.1, D29).

A developed profile of a periodic pipe or channel run (D111) feeds a jet
or a duct as a table profile with turbulence columns:

    python3 -m foehncraft.profiles output/jet_pipe/line_profile.csv \\
        examples/cases/flow/impinging_jet_profile.csv

`inlet_table` reads a line of `[[output.lines]]` that runs from the axis
(or the centre plane) to the wall. The solver samples the turbulence
variables as cell values, so the samples of one cell form a run of equal
values: the table takes one row per run, at the mean distance of its
samples (the centre of the cell's chord), with the velocity there, and
ends with a wall row at the last sample (u = 0, k = 0, ζ = 0; ε and ω
of the wall cell).
The columns are the coordinate `r` (the distance along the line), `u_n`
(the speed along the line's flow direction), `k`, `epsilon` or `omega`,
and `zeta` (ζ of k-ε-ζ-f, φ of the elliptic-blending k-ε) when the model
has it.
"""

from __future__ import annotations

import sys
from pathlib import Path

from .results import read_csv

#: The turbulence columns of a line and their names in an inlet table.
TURBULENCE = (
    ("k", "k"),
    ("epsilon", "epsilon"),
    ("omega", "omega"),
    ("zeta", "zeta"),
    ("phi", "zeta"),
)


def inlet_table(
    line: dict[str, list] | str | Path,
    out: str | Path | None = None,
    velocity: str = "u",
    header: str = "",
) -> dict[str, list]:
    """The inlet table of the line `line` (a file or the columns of
    `read_csv`), written to `out` when given. `velocity` is the column of
    the speed (`u` for a pipe along x); its sign is kept."""
    cols = read_csv(line) if isinstance(line, (str, Path)) else line
    s = [float(v) for v in cols["s"]]
    names = [(src, dst) for src, dst in TURBULENCE if src in cols]
    if not any(src == "k" for src, _ in names):
        raise ValueError("the line has no column `k` (a laminar run?)")
    keys = [src for src, _ in names]

    def same(a, b) -> bool:
        return abs(float(a) - float(b)) <= 1e-12 * max(abs(float(a)), abs(float(b)))

    # the last sample lies on the wall: the wall row stands for it
    n = len(s) - 1
    rows: list[list[float]] = []
    start = 0
    for i in range(1, n + 1):
        if i == n or not all(same(cols[c][i], cols[c][start]) for c in keys):
            run = range(start, i)
            centre = sum(s[j] for j in run) / len(run)
            # the velocity at the sample nearest the centre (it is
            # reconstructed inside the cell, unlike the turbulence)
            j = min(run, key=lambda j: abs(s[j] - centre))
            rows.append([centre, float(cols[velocity][j])] + [float(cols[c][start]) for c in keys])
            start = i
    wall = [s[-1], 0.0]
    for src, dst in names:
        wall.append(0.0 if dst in ("k", "zeta") else rows[-1][2 + keys.index(src)])
    rows.append(wall)
    table = {"r": [r[0] for r in rows], "u_n": [r[1] for r in rows]}
    for n, (_, dst) in enumerate(names):
        table[dst] = [r[2 + n] for r in rows]
    if out is not None:
        write(table, out, header)
    return table


def write(table: dict[str, list], out: str | Path, header: str = "") -> None:
    """Writes the columns `table` as an inlet table with comment lines
    `header`."""
    names = list(table)
    lines = [f"# {h}" if h else "#" for h in header.splitlines()] if header else []
    lines.append(", ".join(names))
    for i in range(len(table[names[0]])):
        lines.append(", ".join(repr(float(table[n][i])) for n in names))
    Path(out).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print("usage: python3 -m foehncraft.profiles LINE.csv TABLE.csv", file=sys.stderr)
        return 2
    src = Path(args[0])
    first = src.read_text(encoding="utf-8").splitlines()[0].lstrip("# ")
    t = inlet_table(src, args[1], header=f"inlet table from {src.name} ({first})")
    print(f"{args[1]}: {len(t['r'])} rows, columns {', '.join(t)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

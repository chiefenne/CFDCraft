"""Parametric studies (D62): variants of a base case, run and tabulated.

    from foehncraft import Case
    from foehncraft.study import Study, grid

    base = Case.load("examples/cases/flow/couette_hex.toml")
    study = Study(base, grid({"fluid.viscosity": [0.05, 0.1, 0.2],
                              "boundary.ymax.velocity": [[1.0, 0, 0], [2.0, 0, 0]]}),
                  root="output/couette_study")
    table = study.run(jobs=2)            # also writes output/couette_study/study.csv
    table[0]["errors.l2"]

The variants come from a full grid, a list, a Latin hypercube or a Sobol
sequence (the last two with SciPy). Each variant is a case file of its
own (`<root>/<id>/<id>.toml`) with its output next to it, so every row
of the table can be repeated from its case file. A study that is run
again skips the variants whose `result.json` belongs to the same case
file (resume), and a failed variant is a row with its status, never the
end of the study.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from ..case import Case
from ..results import Result

Variant = dict[str, Any]

__all__ = ["Study", "grid", "listing", "latin_hypercube", "sobol", "write_table"]


def grid(params: dict[str, Sequence]) -> list[Variant]:
    """Every combination of the values (the last key varies fastest)."""
    keys = list(params)
    return [dict(zip(keys, combo)) for combo in itertools.product(*(params[k] for k in keys))]


def listing(rows: Iterable[Variant]) -> list[Variant]:
    """The variants as given (a list of {key: value})."""
    return [dict(r) for r in rows]


def _scale(
    unit: Sequence[Sequence[float]], bounds: dict[str, tuple[float, float]]
) -> list[Variant]:
    keys = list(bounds)
    out = []
    for u in unit:
        out.append(
            {k: bounds[k][0] + float(ui) * (bounds[k][1] - bounds[k][0]) for k, ui in zip(keys, u)}
        )
    return out


def latin_hypercube(n: int, bounds: dict[str, tuple[float, float]], seed: int = 0) -> list[Variant]:
    """`n` variants of a Latin hypercube in the box `bounds` (SciPy)."""
    from scipy.stats import qmc

    return _scale(qmc.LatinHypercube(d=len(bounds), seed=seed).random(n), bounds)


def sobol(n: int, bounds: dict[str, tuple[float, float]], seed: int | None = None) -> list[Variant]:
    """`n` points of a Sobol sequence in the box `bounds` (SciPy;
    scrambled with `seed`, plain without; `n` a power of two keeps the
    balance properties)."""
    from scipy.stats import qmc

    s = qmc.Sobol(d=len(bounds), scramble=seed is not None, seed=seed)
    return _scale(s.random(n), bounds)


def _flat_value(v: Any) -> Any:
    if isinstance(v, (list, tuple)):
        return json.dumps([_flat_value(x) for x in v])
    if hasattr(v, "__toml__"):
        return v.__toml__()
    return v


class Study:
    """Variants of `base` (dicts of dotted case keys) under `root`."""

    def __init__(
        self,
        base: Case,
        variants: Sequence[Variant],
        root,
        *,
        name: Callable[[int, Variant], str] | None = None,
        columns: Sequence[str] | None = None,
    ):
        self.base = base
        self.variants = list(variants)
        self.root = Path(root)
        self.name = name or (lambda i, v: f"v{i:04d}")
        self.columns = columns

    def cases(self) -> list[tuple[str, Case]]:
        """The case of every variant with its id."""
        out = []
        for i, v in enumerate(self.variants):
            vid = self.name(i, v)
            c = self.base.with_values(v)
            c.name = vid
            out.append((vid, c))
        return out

    def _one(self, vid: str, case: Case, run_kw: dict) -> dict:
        from ..run import FoehncraftError, run

        d = self.root / vid
        path = case.save(d / f"{vid}.toml")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        stamp = d / ".case_sha256"
        res = None
        if (d / "result.json").exists() and stamp.exists() and stamp.read_text() == digest:
            res = Result.load(d)
        else:
            try:
                res = run(path, out=d, **run_kw)
                stamp.write_text(digest)
            except (FoehncraftError, OSError) as e:
                return {"id": vid, "status": "error", "exit_code": 1, "message": str(e)}
        row: dict[str, Any] = {"id": vid, "status": res.status, "exit_code": res.exit_code}
        flat = res.flat()
        row.update(flat)
        if res.orders:
            row["orders"] = json.dumps(res.orders)
        return row

    def run(self, jobs: int = 1, table: str = "study.csv", **run_kw) -> list[dict]:
        """Runs every variant (`jobs` at once; with the command line
        backend the runs are separate processes, in process they share
        the threads of the solver) and writes the table `<root>/<table>`.
        Returns the rows: the parameters, the status and the numbers of
        the last solve with dotted keys."""
        self.root.mkdir(parents=True, exist_ok=True)
        cases = self.cases()
        if jobs > 1:
            run_kw.setdefault("backend", "cli")
            with ThreadPoolExecutor(jobs) as ex:
                rows = list(ex.map(lambda vc: self._one(vc[0], vc[1], run_kw), cases))
        else:
            rows = [self._one(vid, c, run_kw) for vid, c in cases]
        for row, v in zip(rows, self.variants):
            for k, x in v.items():
                row[k] = _flat_value(x)
        write_table(
            self.root / table, rows, list(self.variants[0]) if self.variants else [], self.columns
        )
        return rows


def write_table(
    path, rows: list[dict], params: Sequence[str] = (), columns: Sequence[str] | None = None
) -> Path:
    """Writes the study table as CSV: id, the parameters, status, exit
    code, then `columns` (default every number of the rows)."""
    head = ["id", *params, "status", "exit_code"]
    if columns is None:
        rest: list[str] = []
        for r in rows:
            for k in r:
                if k not in head and k not in rest:
                    rest.append(k)
    else:
        rest = list(columns)
    p = Path(path)
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=head + rest, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in head + rest})
    return p

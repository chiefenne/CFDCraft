"""Running cases from Python (D61): in process or through the command line.

    from foehncraft import Case, run

    case = Case.load("examples/cases/flow/cavity_re100.toml")
    res = run(case, out="output/cavity_py",
              on_iteration=lambda row: row["continuity"] < 1e-8)  # True stops
    res.ok, res.status, res.solves[0]["iterations"]

Two backends with the same results and files:

- `"native"`: the solver runs in this process (the native module of the
  package). The callbacks see every row of the run record (`run.jsonl`)
  as it is written; a true return value stops the run as a stop request
  does (restart file and final output, exit code 3).
- `"cli"`: the `foehn` command line (from `FOEHN`, `PATH` or the
  repository's `target/release`). The callbacks see the rows of
  `run.jsonl` while it grows; a true return value writes the `.stop`
  control file (D40).

`backend="auto"` takes the native module when it is there. Calls go to
Python between outer iterations or time steps, never per cell, so the
speed and the reproducibility of the solver are kept (D53, D31).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import threading
import time
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from .case import Case
from .results import Result

Row = dict[str, Any]
Callback = Callable[[Row], Any]


class FoehncraftError(RuntimeError):
    """An error of the solver, with its diagnostic code (D78)."""

    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


def native_available() -> bool:
    """True when the native module is built."""
    try:
        from . import _native  # noqa: F401
    except ImportError:
        return False
    return True


def find_foehn(foehn: str | None = None) -> str:
    """The `foehn` executable: the argument, `FOEHN`, `PATH`, or the
    release build of the repository this package lives in."""
    for c in (foehn, os.environ.get("FOEHN")):
        if c:
            return c
    on_path = shutil.which("foehn")
    if on_path:
        return on_path
    root = Path(__file__).resolve().parents[2]
    for profile in ("release", "ci"):
        p = root / "target" / profile / ("foehn.exe" if os.name == "nt" else "foehn")
        if p.exists():
            return str(p)
    raise FileNotFoundError(
        "no `foehn` executable: build it (cargo build --release -p foehncraft-cli), "
        "put it on PATH or set FOEHN"
    )


def _backend(backend: str) -> str:
    if backend == "auto":
        return "native" if native_available() else "cli"
    if backend not in ("native", "cli"):
        raise ValueError(f"backend must be 'auto', 'native' or 'cli', not {backend!r}")
    if backend == "native" and not native_available():
        raise ImportError(
            "the native module is not built: pip install ./python (maturin) or use backend='cli'"
        )
    return backend


def _output_dir(case: Case, out) -> Path:
    if out is not None:
        return Path(out)
    d = case.get("output.directory")
    if d:
        return case.base_dir / d
    return Path("output") / case.name


def case_file(case, out=None) -> tuple[Path, Path]:
    """The case file to run and the output directory. A `Case` is written
    into the output directory as `<name>.toml` (its relative paths
    rewritten), so every run can be repeated from its case file."""
    if isinstance(case, (str, os.PathLike)):
        path = Path(case)
        c = Case.load(path)
        return path, Path(out) if out is not None else _output_dir(c, None)
    if not isinstance(case, Case):
        raise TypeError(f"run() takes a Case or a path, not {type(case).__name__}")
    d = _output_dir(case, out)
    d.mkdir(parents=True, exist_ok=True)
    return case.save(d / f"{case.name}.toml"), d


def _dispatch(on_row, on_iteration, on_step) -> Callable[[Row], bool]:
    def f(row: Row) -> bool:
        stop = False
        if on_row is not None:
            stop = bool(on_row(row)) or stop
        kind = row.get("type")
        if kind == "iteration" and on_iteration is not None:
            stop = bool(on_iteration(row)) or stop
        if kind == "step" and on_step is not None:
            stop = bool(on_step(row)) or stop
        return stop

    return f


def run(
    case,
    out=None,
    *,
    sets: Sequence[str] | dict = (),
    on_row: Callback | None = None,
    on_iteration: Callback | None = None,
    on_step: Callback | None = None,
    log: Callable[[str], Any] | None = None,
    backend: str = "auto",
    foehn: str | None = None,
) -> Result:
    """Runs `case` (a `Case` or a case file) and returns its `Result`.

    `sets` are `--set` overrides (`["fluid.viscosity=1e-3"]` or a dict).
    `on_row` sees every row of the run record, `on_iteration` the outer
    iterations of a steady run, `on_step` the time steps of a transient
    run; any of them returning a true value stops the run. `log` gets the
    progress lines.
    """
    sets = _sets(sets)
    path, out_dir = case_file(case, out)
    callback = (
        _dispatch(on_row, on_iteration, on_step) if (on_row or on_iteration or on_step) else None
    )
    if _backend(backend) == "native":
        from . import _native

        observer = None
        if callback is not None:

            def observer(text: str) -> bool:
                return callback(json.loads(text))

        try:
            summary = _native.run(str(path), list(sets), str(out_dir), observer, log)
        except _native.FoehncraftError as e:
            code, message = e.args
            raise FoehncraftError(code, message) from None
        return Result(out_dir, json.loads(summary))
    return _run_cli(path, out_dir, sets, callback, log, foehn)


def _sets(sets) -> list[str]:
    from . import _toml

    if isinstance(sets, dict):
        return [f"{k}={_toml.value(v)}" for k, v in sets.items()]
    return [str(s) for s in sets]


def _run_cli(path: Path, out_dir: Path, sets, callback, log, foehn) -> Result:
    exe = find_foehn(foehn)
    cmd = [exe, "--json", "flow", str(path), "--out", str(out_dir)]
    for s in sets:
        cmd += ["--set", s]
    out_dir.mkdir(parents=True, exist_ok=True)
    record = out_dir / "run.jsonl"
    started = time.time()
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8"
    )
    done = threading.Event()

    def follow():
        # the run record grows while the run goes on
        pos = 0
        while True:
            finished = done.is_set()
            if record.exists() and record.stat().st_mtime >= started - 1:
                with open(record, encoding="utf-8") as f:
                    f.seek(pos)
                    while True:
                        line = f.readline()
                        if not line.endswith("\n"):
                            break
                        pos = f.tell()
                        if callback(json.loads(line)):
                            (out_dir / ".stop").touch()
            if finished:
                return
            time.sleep(0.05)

    follower = threading.Thread(target=follow, daemon=True) if callback else None
    if follower:
        follower.start()

    def lines():
        for line in proc.stderr:
            if log is not None:
                log(line.rstrip("\n"))

    logger = threading.Thread(target=lines, daemon=True)
    logger.start()
    stdout = proc.stdout.read()
    proc.wait()
    logger.join()
    done.set()
    if follower:
        follower.join()
    try:
        summary = json.loads(stdout)
    except json.JSONDecodeError:
        raise FoehncraftError(
            "RUN001", f"{exe} gave no JSON (exit code {proc.returncode})"
        ) from None
    if summary.get("format") != "foehncraft.result":
        diags = summary.get("diagnostics") or [{}]
        d = diags[0]
        raise FoehncraftError(d.get("code", "RUN001"), d.get("message", json.dumps(summary)))
    stop = out_dir / ".stop"
    if stop.exists():
        stop.unlink()
    return Result(out_dir, summary)


def check(case, out=None, *, sets=(), backend: str = "auto", foehn: str | None = None) -> dict:
    """The check report of a case (`foehncraft.check`): validation, the
    boundary-condition check on the first mesh and the estimate."""
    sets = _sets(sets)
    if isinstance(case, Case):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            # the relative paths are rewritten to the temporary directory
            path = case.save(Path(tmp) / f"{case.name}.toml")
            return _check(path, out, sets, backend, foehn)
    return _check(Path(case), out, sets, backend, foehn)


def _check(path: Path, out, sets, backend, foehn) -> dict:
    if _backend(backend) == "native":
        from . import _native

        return json.loads(_native.check(str(path), list(sets), None if out is None else str(out)))
    cmd = [find_foehn(foehn), "--json", "flow", str(path), "--check"]
    if out is not None:
        cmd += ["--out", str(out)]
    for s in sets:
        cmd += ["--set", s]
    p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8")
    return json.loads(p.stdout)

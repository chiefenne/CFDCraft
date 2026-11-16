"""Case files in Python (D61): read, change, build and write them.

A `Case` is the nested table of a TOML case file. Keys are reached with
dots, the way `--set` names them::

    from foehncraft import Case
    from foehncraft.expr import y, constant

    case = Case.load("examples/cases/flow/poiseuille_hex.toml")
    case["fluid.viscosity"] = 0.05
    case["boundary.xmin.velocity"] = [6 * y * (1 - y), 0.0, 0.0]
    del case["boundary.xmin.profile"]
    case.save("output/study/variant_1.toml")

The TOML a case writes reads back to the same table (round trip), and a
run of the written file is the run of the case: what Python sets, the
case file holds. Relative paths of the case (mesh file, tables) are
rewritten when the case is saved to another directory, so they keep
pointing at the same files.
"""

from __future__ import annotations

import copy
import fnmatch
import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from . import _toml

#: Keys whose values are paths relative to the case file (glob patterns
#: over the dotted key).
PATH_KEYS = (
    "mesh.file",
    "tables.*.file",
    "boundary.*.table",
    "materials.*.density.file",
    "materials.*.viscosity.file",
    "run.restart",
    "output.directory",
)


def _split(key) -> list[str]:
    if isinstance(key, (tuple, list)):
        return [str(k) for k in key]
    return str(key).split(".")


class Case:
    """A flow case: the table of a case file and the directory its relative
    paths start from."""

    def __init__(self, data: dict | None = None, base_dir=".", name: str = "case"):
        self.data: dict[str, Any] = data if data is not None else {}
        self.base_dir = Path(base_dir)
        self.name = name

    # -- reading and writing
    @classmethod
    def load(cls, path) -> Case:
        """Reads a case file."""
        p = Path(path)
        return cls(_toml.load(p), p.parent, p.stem)

    @classmethod
    def from_toml(cls, text: str, base_dir=".", name: str = "case") -> Case:
        """A case from TOML text; relative paths start at `base_dir`."""
        return cls(_toml.loads(text), base_dir, name)

    def to_toml(self, header: str = "") -> str:
        """The case as TOML text."""
        return _toml.dumps(self.data, header)

    def save(self, path, header: str = "written by the foehncraft Python package") -> Path:
        """Writes the case file `path`; relative paths of the case are
        rewritten to stay valid from the new directory. Returns the path."""
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        moved = self.relocated(p.parent)
        p.write_text(moved.to_toml(header), encoding="utf-8")
        return p

    def relocated(self, new_dir) -> Case:
        """A copy whose relative paths start at `new_dir`."""
        new_dir = Path(new_dir)
        c = self.copy()
        c.base_dir = new_dir
        for key, value in self.paths():
            if not isinstance(value, str) or value == "latest" or os.path.isabs(value):
                continue
            target = (self.base_dir / value).resolve()
            try:
                rel = os.path.relpath(target, new_dir.resolve())
            except ValueError:  # another drive on Windows
                rel = str(target)
            c[key] = Path(rel).as_posix()
        return c

    def paths(self) -> Iterator[tuple[str, Any]]:
        """The keys that hold paths, with their values."""
        for key, value in self.items():
            if any(fnmatch.fnmatchcase(key, pat) for pat in PATH_KEYS):
                yield key, value

    def copy(self) -> Case:
        """A deep copy."""
        return Case(copy.deepcopy(self.data), self.base_dir, self.name)

    # -- dotted keys
    def __getitem__(self, key):
        d: Any = self.data
        for k in _split(key):
            if not isinstance(d, dict) or k not in d:
                raise KeyError(key)
            d = d[k]
        return d

    def __setitem__(self, key, value) -> None:
        parts = _split(key)
        d = self.data
        for k in parts[:-1]:
            d = d.setdefault(k, {})
            if not isinstance(d, dict):
                raise KeyError(f"{key}: `{k}` is not a table")
        d[parts[-1]] = value

    def __delitem__(self, key) -> None:
        parts = _split(key)
        d = self.data
        for k in parts[:-1]:
            d = d[k]
        del d[parts[-1]]

    def __contains__(self, key) -> bool:
        try:
            self[key]
        except KeyError:
            return False
        return True

    def get(self, key, default=None):
        """The value of `key`, or `default`."""
        try:
            return self[key]
        except KeyError:
            return default

    def set(self, values: dict | None = None, **kw) -> Case:
        """Sets several dotted keys (also as keywords with `__` for dots);
        returns the case."""
        for k, v in {**(values or {}), **{k.replace("__", "."): v for k, v in kw.items()}}.items():
            self[k] = v
        return self

    def with_values(self, values: dict) -> Case:
        """A copy with the dotted keys of `values` set."""
        return self.copy().set(values)

    def items(self, prefix: str = "") -> Iterator[tuple[str, Any]]:
        """Every leaf value with its dotted key (arrays of tables by index)."""

        def walk(d, path):
            for k, v in d.items():
                key = f"{path}.{k}" if path else k
                if isinstance(v, dict):
                    yield from walk(v, key)
                else:
                    yield key, v

        for k, v in walk(self.data, ""):
            if k.startswith(prefix):
                yield k, v

    def __eq__(self, other) -> bool:
        return isinstance(other, Case) and _plain(self.data) == _plain(other.data)

    def __repr__(self) -> str:
        return f"Case({self.name!r}, base_dir={str(self.base_dir)!r})"

    # -- the solver
    def check(self, out=None, sets=()) -> dict:
        """The check report of the solver (`foehncraft.check`)."""
        from .run import check

        return check(self, out=out, sets=sets)

    def run(self, out=None, **kw):
        """Runs the case; see `foehncraft.run.run`."""
        from .run import run

        return run(self, out=out, **kw)


def _plain(v):
    """Values with expression objects as their strings (for comparison)."""
    if isinstance(v, dict):
        return {k: _plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    if hasattr(v, "__toml__"):
        return v.__toml__()
    return v

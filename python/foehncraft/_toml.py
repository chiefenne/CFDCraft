"""Reading and writing TOML for case files.

Reading uses `tomllib` (Python 3.11+) or `tomli`. Writing is a small
writer for the subset of TOML that case files use: tables, arrays of
tables, inline tables inside arrays, strings, booleans, integers and
floats. Floats are written with `repr`, the shortest text that reads back
to the same bits, so a value set in Python reaches the solver unchanged.
"""

from __future__ import annotations

import json
import math
import re
from typing import Any

try:  # Python 3.11+
    import tomllib as _tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    import tomli as _tomllib

_BARE = re.compile(r"^[A-Za-z0-9_-]+$")


def loads(text: str) -> dict[str, Any]:
    """Parses TOML text."""
    return _tomllib.loads(text)


def load(path) -> dict[str, Any]:
    """Reads a TOML file."""
    with open(path, "rb") as f:
        return _tomllib.load(f)


def _key(k: str) -> str:
    return k if _BARE.match(k) else json.dumps(k, ensure_ascii=False)


def _is_table(v: Any) -> bool:
    return isinstance(v, dict)


def _is_table_array(v: Any) -> bool:
    return isinstance(v, list) and len(v) > 0 and all(isinstance(x, dict) for x in v)


def value(v: Any) -> str:
    """One value as TOML text."""
    # expression objects and other objects with a TOML form
    to_toml = getattr(v, "__toml__", None)
    if to_toml is not None:
        return value(to_toml())
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if math.isnan(v):
            return "nan"
        if math.isinf(v):
            return "inf" if v > 0 else "-inf"
        return repr(v)
    if isinstance(v, str):
        return json.dumps(v, ensure_ascii=False)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(value(x) for x in v) + "]"
    if isinstance(v, dict):
        return "{ " + ", ".join(f"{_key(k)} = {value(x)}" for k, x in v.items()) + " }"
    # numpy scalars and similar
    if hasattr(v, "item"):
        return value(v.item())
    raise TypeError(f"cannot write {type(v).__name__} to TOML: {v!r}")


def _table(out: list[str], path: list[str], d: dict[str, Any]) -> None:
    plain = [(k, v) for k, v in d.items() if not (_is_table(v) or _is_table_array(v))]
    tables = [(k, v) for k, v in d.items() if _is_table(v)]
    arrays = [(k, v) for k, v in d.items() if _is_table_array(v)]
    if path and (plain or not (tables or arrays)):
        out.append("")
        out.append("[" + ".".join(_key(p) for p in path) + "]")
    for k, v in plain:
        if v is None:
            continue
        out.append(f"{_key(k)} = {value(v)}")
    for k, v in tables:
        _table(out, [*path, k], v)
    for k, v in arrays:
        for item in v:
            out.append("")
            out.append("[[" + ".".join(_key(p) for p in [*path, k]) + "]]")
            for kk, vv in item.items():
                if vv is None:
                    continue
                out.append(f"{_key(kk)} = {value(vv)}")


def dumps(d: dict[str, Any], header: str = "") -> str:
    """Writes a nested dict as TOML text; `header` becomes comment lines."""
    out: list[str] = [f"# {line}".rstrip() for line in header.splitlines()]
    _table(out, [], d)
    text = "\n".join(out).strip("\n") + "\n"
    return text

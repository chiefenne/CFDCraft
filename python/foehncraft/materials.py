"""The material database and materials of a case (D100, D102).

    from foehncraft import materials

    materials.database(temperature=313.15)       # every entry, evaluated
    w = materials.get("water", temperature=313.15)
    w["density"]["value"], w["viscosity"]["value"]

    # a material of the case, built on a database entry
    case["materials.warm-oil"] = materials.material(
        base="engine-oil",
        viscosity=materials.table([300.0, 320.0, 340.0], [0.49, 0.14, 0.053],
                                  interpolation="monotone-cubic", source="own data"),
    )
    case["fluid.material"] = "warm-oil"

The values come from the solver itself (the native module or
`foehn --json materials`), so Python shows what a run uses.
"""

from __future__ import annotations

import json
import subprocess
from collections.abc import Sequence
from typing import Any


def database(
    temperature: float = 293.15, pressure: float = 101325.0, *, backend: str = "auto"
) -> list[dict]:
    """Every database entry evaluated at `temperature` (K) and `pressure` (Pa)."""
    return _query(None, temperature, pressure, backend)["materials"]


def get(
    name: str, temperature: float = 293.15, pressure: float = 101325.0, *, backend: str = "auto"
) -> dict:
    """One database entry evaluated at `temperature` and `pressure`."""
    return _query(name, temperature, pressure, backend)["materials"][0]


def database_sha256(*, backend: str = "auto") -> str:
    """The SHA-256 of the database the solver was built with."""
    return _query(None, 293.15, 101325.0, backend)["database_sha256"]


def _query(name, temperature, pressure, backend) -> dict:
    from .run import FoehncraftError, _backend, find_foehn

    if _backend(backend) == "native":
        from . import _native

        try:
            return json.loads(_native.material_database(name, temperature, pressure))
        except _native.FoehncraftError as e:
            raise FoehncraftError(*e.args) from None
    cmd = [find_foehn(), "--json", "materials"]
    if name:
        cmd.append(name)
    cmd += ["--temperature", repr(float(temperature)), "--pressure", repr(float(pressure))]
    p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8")
    return json.loads(p.stdout)


def _drop_none(d: dict) -> dict:
    return {k: v for k, v in d.items() if v is not None}


def material(
    *,
    base: str | None = None,
    kind: str | None = None,
    description: str | None = None,
    density: Any = None,
    viscosity: Any = None,
    molar_mass: float | None = None,
    temperature_range: Sequence[float] | None = None,
) -> dict:
    """A `[materials.<name>]` table. Each property is a number, an
    expression (an `Expr` or a string) or a form from the functions below."""
    return _drop_none(
        {
            "base": base,
            "kind": kind,
            "description": description,
            "molar_mass": molar_mass,
            "temperature_range": list(temperature_range) if temperature_range else None,
            "density": density,
            "viscosity": viscosity,
        }
    )


def constant(value: float, source: str | None = None) -> dict:
    """A constant property with its source."""
    return _drop_none({"form": "constant", "value": float(value), "source": source})


def polynomial(coefficients: Sequence[float], *, range=None, source=None) -> dict:  # noqa: A002
    """Σ aᵢ Tⁱ, T in K."""
    return _drop_none(
        {
            "form": "polynomial",
            "coefficients": [float(c) for c in coefficients],
            "range": list(range) if range else None,
            "source": source,
        }
    )


def sutherland(mu_ref: float, t_ref: float, s: float, *, range=None, source=None) -> dict:  # noqa: A002
    """μ = μ_ref (T/T_ref)^1.5 (T_ref + S)/(T + S)."""
    return _drop_none(
        {
            "form": "sutherland",
            "mu_ref": float(mu_ref),
            "t_ref": float(t_ref),
            "s": float(s),
            "range": list(range) if range else None,
            "source": source,
        }
    )


def power_law(value_ref: float, t_ref: float, exponent: float, *, range=None, source=None) -> dict:  # noqa: A002
    """value_ref (T/T_ref)^exponent."""
    return _drop_none(
        {
            "form": "power-law",
            "value_ref": float(value_ref),
            "t_ref": float(t_ref),
            "exponent": float(exponent),
            "range": list(range) if range else None,
            "source": source,
        }
    )


def table(
    temperature: Sequence[float] | None = None,
    values: Sequence[float] | None = None,
    *,
    file: str | None = None,
    interpolation: str = "linear",
    source: str | None = None,
) -> dict:
    """A table in T (K), `linear` or `monotone-cubic`, clamped outside its
    range with a warning; inline or from a two-column file."""
    return _drop_none(
        {
            "form": "table",
            "temperature": [float(v) for v in temperature] if temperature is not None else None,
            "values": [float(v) for v in values] if values is not None else None,
            "file": file,
            "interpolation": interpolation,
            "source": source,
        }
    )


def ideal_gas(source: str | None = None) -> dict:
    """ρ = p M / (R T) (needs `molar_mass`)."""
    return _drop_none({"form": "ideal-gas", "source": source})

"""Case files: round trip with TOML, dotted keys, relocated paths."""

from pathlib import Path

import pytest
from conftest import CASES, ROOT

from foehncraft import Case, _toml
from foehncraft.expr import y

ALL = sorted(
    p
    for d in ("flow", "scalar", "poisson", "mesh")
    for p in (ROOT / "examples/cases" / d).glob("*.toml")
)


@pytest.mark.parametrize("path", ALL, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_every_example_case_round_trips(path):
    c = Case.load(path)
    text = c.to_toml()
    assert _toml.loads(text) == c.data
    # and once more: writing is stable
    assert Case.from_toml(text).to_toml() == text


def test_dotted_keys():
    c = Case.load(CASES / "poiseuille_hex.toml")
    assert c["fluid.viscosity"] == 0.1
    c["fluid.viscosity"] = 0.05
    c["boundary.xmin.velocity"] = [6 * y * (1 - y), 0.0, 0.0]
    del c["boundary.xmin.profile"]
    assert "boundary.xmin.profile" not in c
    assert c.get("no.such.key", 7) == 7
    c.set(solver__tolerance=1e-8)
    assert c["solver.tolerance"] == 1e-8
    data = _toml.loads(c.to_toml())
    assert data["boundary"]["xmin"]["velocity"] == ["6 * y * (1 - y)", 0.0, 0.0]
    assert data["fluid"]["viscosity"] == 0.05


def test_floats_keep_their_bits():
    c = Case()
    vals = [0.1, 1 / 3, 2.0**-1074, 1.7976931348623157e308, -0.0, 1e-5, 123456789.123456789]
    c["constants"] = {f"c{i}": v for i, v in enumerate(vals)}
    back = _toml.loads(c.to_toml())["constants"]
    for i, v in enumerate(vals):
        assert back[f"c{i}"].hex() == v.hex()


def test_saved_elsewhere_keeps_the_relative_paths(tmp_path):
    c = Case.load(CASES / "poiseuille_table.toml")
    p = c.save(tmp_path / "deep" / "dir" / "case.toml")
    back = Case.load(p)
    table = back["boundary.xmin.table"]
    assert (p.parent / table).resolve() == (CASES / "poiseuille_profile.csv").resolve()
    assert Path(c["boundary.xmin.table"]) == Path("poiseuille_profile.csv")

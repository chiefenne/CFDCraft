"""The material database from Python and the study layer (D62)."""

import csv

import pytest
from conftest import CASES, needs_native

from foehncraft import Case, materials
from foehncraft.study import Study, grid, latin_hypercube, listing, sobol


@needs_native
def test_database():
    db = {m["name"]: m for m in materials.database(300.0)}
    assert {"air", "water", "engine-oil", "ethylene-glycol", "glycerol"} <= db.keys()
    assert db["water"]["density"]["value"] == pytest.approx(996.56, abs=0.01)
    assert len(materials.database_sha256()) == 64
    hot = materials.get("water", 400.0)
    assert hot["warnings"], "outside the range of the table: clamped with a warning"


@needs_native
def test_a_material_built_in_python(tmp_path):
    c = Case.load(CASES / "poiseuille_water.toml")
    c["materials.warm-water"] = materials.material(
        base="water",
        viscosity=materials.table([300.0, 320.0], [0.85e-3, 0.58e-3], interpolation="linear"),
    )
    c["fluid.material"] = "warm-water"
    c["regions.channel.material"] = "warm-water"
    rep = c.check()
    assert rep["ok"], rep["diagnostics"]


def test_variants():
    g = grid({"a": [1, 2], "b": ["x", "y", "z"]})
    assert len(g) == 6 and g[1] == {"a": 1, "b": "y"}
    assert listing([{"a": 1}]) == [{"a": 1}]
    pytest.importorskip("scipy")
    lhs = latin_hypercube(8, {"a": (0.0, 1.0), "b": (10.0, 20.0)}, seed=3)
    assert len(lhs) == 8 and all(0 <= v["a"] <= 1 and 10 <= v["b"] <= 20 for v in lhs)
    # one sample in every eighth of each axis
    assert sorted(int(v["a"] * 8) for v in lhs) == list(range(8))
    s = sobol(4, {"a": (0.0, 1.0)})
    assert sorted(v["a"] for v in s) == [0.0, 0.25, 0.5, 0.75]


@needs_native
def test_study_runs_tabulates_and_resumes(tmp_path):
    base = Case.load(CASES / "couette_hex.toml")
    st = Study(base, grid({"fluid.viscosity": [0.05, 0.1]}), tmp_path / "s")
    rows = st.run()
    assert [r["status"] for r in rows] == ["finished", "finished"]
    with open(tmp_path / "s" / "study.csv") as f:
        table = list(csv.DictReader(f))
    assert [t["fluid.viscosity"] for t in table] == ["0.05", "0.1"]
    assert float(table[0]["errors.l2"]) < 1e-8
    # run again: nothing is run, the results are read back
    stamp = (tmp_path / "s" / "v0000" / "result.json").stat().st_mtime_ns
    st.run()
    assert (tmp_path / "s" / "v0000" / "result.json").stat().st_mtime_ns == stamp

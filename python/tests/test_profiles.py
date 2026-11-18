"""Inlet tables from sampled profiles (stage 6 specification 7.1): one row
per run of equal turbulence samples (a cell), at the mean distance of the
run, and a wall row."""

import math

from foehncraft.profiles import inlet_table, write


def line():
    """Five samples of two cells and the wall sample: the turbulence
    columns are cell values, u is reconstructed."""
    return {
        "s": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5],
        "u": [1.2, 1.19, 1.15, 1.0, 0.8, 0.0],
        "k": [0.01, 0.01, 0.01, 0.02, 0.02, 0.02],
        "epsilon": [0.1, 0.1, 0.1, 0.4, 0.4, 0.4],
        "phi": [0.5, 0.5, 0.5, 0.3, 0.3, 0.3],
    }


def test_one_row_per_cell_and_a_wall_row(tmp_path):
    t = inlet_table(line())
    assert list(t) == ["r", "u_n", "k", "epsilon", "zeta"]
    # the cells: samples 0-2 (mean 0.1) and 3-4 (mean 0.35, u at the
    # first of the two samples nearest it); the wall row at the last sample
    assert all(math.isclose(a, b) for a, b in zip(t["r"], [0.1, 0.35, 0.5], strict=True))
    assert t["u_n"] == [1.19, 1.0, 0.0]
    assert t["k"] == [0.01, 0.02, 0.0]
    assert t["epsilon"] == [0.1, 0.4, 0.4]
    assert t["zeta"] == [0.5, 0.3, 0.0]
    write(t, tmp_path / "t.csv", header="a test\ntable")
    text = (tmp_path / "t.csv").read_text()
    assert text.startswith("# a test\n# table\nr, u_n, k, epsilon, zeta\n")
    rows = [[float(v) for v in r.split(",")] for r in text.splitlines()[3:]]
    assert [r[2] for r in rows] == t["k"]


def test_a_laminar_line_is_refused():
    cols = {k: v for k, v in line().items() if k in ("s", "u")}
    try:
        inlet_table(cols)
    except ValueError as e:
        assert "`k`" in str(e)
    else:
        raise AssertionError("a line without k was accepted")

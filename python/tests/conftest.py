"""Shared fixtures: the repository, the example cases and the backends."""

import os
from pathlib import Path

import pytest

from foehncraft import native_available

ROOT = Path(__file__).resolve().parents[2]
CASES = ROOT / "examples" / "cases" / "flow"


def foehn_path():
    """The `foehn` executable for the command-line tests, or None."""
    from foehncraft import find_foehn

    try:
        return find_foehn()
    except FileNotFoundError:
        return None


needs_native = pytest.mark.skipif(not native_available(), reason="native module not built")
needs_foehn = pytest.mark.skipif(foehn_path() is None, reason="no foehn executable (set FOEHN)")


@pytest.fixture
def cases():
    return CASES


@pytest.fixture(autouse=True)
def _in_tmp(tmp_path, monkeypatch):
    # runs write output/<name> below the working directory
    monkeypatch.chdir(tmp_path)
    os.environ.setdefault("RAYON_NUM_THREADS", "2")

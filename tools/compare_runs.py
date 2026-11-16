#!/usr/bin/env python3
"""Bitwise comparison of two runs of the same case (D31), for example the
CPU and the GPU flow loop, or two thread counts.

    python3 tools/compare_runs.py [--across-versions] DIR_A DIR_B

DIR_A and DIR_B are output directories of `foehn` (or parents of several,
compared file by file under the same relative paths). Every number must be
bitwise equal:
- HDF5 files (fields.vtkhdf, restart.h5): every dataset, compared as bits;
- CSV files: byte for byte, except the comment lines (`#`, D75), which
  name the version;
- JSON reports: every value, except timings and the device description;
- case.resolved.toml: every line except `backend`, `ptx` and `threads`
  of `[run]`.
Logs, Markdown reports, SVG plots and the run record (run.jsonl) carry
timings and are skipped.

--across-versions compares the runs of two versions of Foehncraft (a
branch against main): the files of DIR_A only, and in JSON reports only
the keys both have, so that new files and new keys are not differences.

Exit code 0 when all files agree, 1 otherwise. Needs h5py and numpy.
"""
import json
import pathlib
import sys

import h5py
import numpy as np

SKIP_SUFFIX = (".log", ".md", ".svg", ".txt", ".jsonl")
SKIP_KEYS = ("seconds", "ns_per_", "timings", "device", "backend", "threads", "gpu")
COMMON = False
RUN_KEYS = ("backend", "ptx", "threads")


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if not any(s in k for s in SKIP_KEYS)}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def common(a, b):
    """a and b reduced to the keys both have (--across-versions)."""
    if isinstance(a, dict) and isinstance(b, dict):
        ks = [k for k in a if k in b]
        pairs = [common(a[k], b[k]) for k in ks]
        return dict(zip(ks, (p[0] for p in pairs))), dict(zip(ks, (p[1] for p in pairs)))
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        pairs = [common(x, y) for x, y in zip(a, b)]
        return [p[0] for p in pairs], [p[1] for p in pairs]
    return a, b


def data_lines(p):
    return [l for l in p.read_bytes().split(b"\n") if not l.startswith(b"#")]


def toml_lines(p):
    out, section = [], ""
    for line in p.read_text().splitlines():
        s = line.strip()
        if s.startswith("["):
            section = s
        elif section == "[run]" and s.split("=")[0].strip() in RUN_KEYS:
            continue
        out.append(line)
    return out


def bits(x):
    x = np.ascontiguousarray(x)
    return x.view(np.uint8) if x.dtype.kind in "fc" else x


def h5_diff(a, b, rel):
    bad = []

    def visit(name, obj):
        if not isinstance(obj, h5py.Dataset):
            return
        if name not in b:
            bad.append(f"{rel}: dataset {name} missing")
            return
        x, y = obj[()], b[name][()]
        if np.shape(x) != np.shape(y):
            bad.append(f"{rel}: {name} shapes {np.shape(x)} and {np.shape(y)}")
        elif isinstance(x, np.ndarray):
            if not np.array_equal(bits(x), bits(y)):
                n = int(np.count_nonzero(bits(x).reshape(len(x), -1) != bits(y).reshape(len(y), -1)) if x.ndim else 1)
                d = float(np.max(np.abs(x - y))) if x.dtype.kind == "f" else "n/a"
                bad.append(f"{rel}: {name} differs ({n} entries, max |a - b| = {d})")
        elif x != y:
            bad.append(f"{rel}: {name} differs")

    a.visititems(visit)
    return bad


def compare(da, db):
    bad, n = [], 0
    for pa in sorted(p for p in da.rglob("*") if p.is_file()):
        rel = pa.relative_to(da)
        pb = db / rel
        if pa.suffix in SKIP_SUFFIX:
            continue
        if not pb.exists():
            bad.append(f"{rel}: missing in {db}")
            continue
        n += 1
        if pa.suffix == ".json":
            x, y = strip(json.loads(pa.read_text())), strip(json.loads(pb.read_text()))
            if COMMON:
                x, y = common(x, y)
            if x != y:
                bad.append(f"{rel}: values differ")
        elif pa.suffix in (".vtkhdf", ".h5", ".hdf5"):
            with h5py.File(pa, "r") as A, h5py.File(pb, "r") as B:
                bad += h5_diff(A, B, rel)
        elif pa.name == "case.resolved.toml":
            if toml_lines(pa) != toml_lines(pb):
                bad.append(f"{rel}: resolved case differs")
        elif pa.suffix == ".csv":
            if data_lines(pa) != data_lines(pb):
                bad.append(f"{rel}: bytes differ")
        elif pa.read_bytes() != pb.read_bytes():
            bad.append(f"{rel}: bytes differ")
    return n, bad


def main():
    global COMMON
    args = sys.argv[1:]
    if args and args[0] == "--across-versions":
        COMMON, args = True, args[1:]
    if len(args) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    da, db = pathlib.Path(args[0]), pathlib.Path(args[1])
    n, bad = compare(da, db)
    for b in bad:
        print("  " + b)
    print(f"{da} vs {db}: {n} files, " + ("bitwise equal" if not bad else f"{len(bad)} differences"))
    return 1 if bad or n == 0 else 0


if __name__ == "__main__":
    sys.exit(main())

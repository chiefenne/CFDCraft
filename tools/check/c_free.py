"""D87: the flow core without default features depends on no C library.

    python -m tools.check.c_free

Fails if the normal and build dependency tree of foehncraft-flow (no
default features, so no HDF5) holds a `-sys` crate or a crate that builds
C code (`cc`, `cmake`, `pkg-config`).
"""

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
C_CRATE = re.compile(r"^([A-Za-z0-9_-]+-sys|cc|cmake|pkg-config) ")


def main() -> int:
    out = subprocess.run(
        ["cargo", "tree", "-p", "foehncraft-flow", "--no-default-features", "-e", "normal,build",
         "--prefix", "none", "--format", "{p}"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout
    tree = sorted({line for line in out.splitlines() if line.strip()})
    bad = [line for line in tree if C_CRATE.match(line)]
    if bad:
        print("foehncraft-flow (no default features) depends on C code:")
        print("\n".join(bad))
        return 1
    print(f"foehncraft-flow (no default features) is C-free: {len(tree)} crates")
    return 0


if __name__ == "__main__":
    sys.exit(main())

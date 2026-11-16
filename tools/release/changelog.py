"""Moves the changelog fragments of changelog.d/ into CHANGELOG.md.

    python -m tools.release.changelog           write them to the top of "## Unreleased", delete them
    python -m tools.release.changelog --print   print "## Unreleased" with the fragments, change nothing

Fragments are named <PR number>-<short-name>.md and are taken in the order
of their numbers (see changelog.d/README.md).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FRAGMENTS = ROOT / "changelog.d"
CHANGELOG = ROOT / "CHANGELOG.md"
HEAD = "## Unreleased"


def fragments() -> list[Path]:
    files = [p for p in FRAGMENTS.glob("*.md") if p.name != "README.md"]
    for p in files:
        if not re.match(r"^\d+-[\w.-]+\.md$", p.name):
            raise SystemExit(f"{p.name}: expected <PR number>-<short-name>.md")
    return sorted(files, key=lambda p: (int(p.name.split("-", 1)[0]), p.name))


def merged() -> str:
    text = CHANGELOG.read_text(encoding="utf-8")
    if HEAD not in text:
        raise SystemExit(f"CHANGELOG.md has no '{HEAD}' section")
    parts = [p.read_text(encoding="utf-8").strip() for p in fragments()]
    parts = [p for p in parts if p]
    if not parts:
        return text
    return text.replace(HEAD + "\n\n", HEAD + "\n\n" + "\n".join(parts) + "\n", 1)


def main(argv: list[str]) -> int:
    new = merged()
    if "--print" in argv:
        start = new.index(HEAD)
        nxt = new.find("\n## ", start + len(HEAD))
        print(new[start : nxt if nxt > 0 else None].rstrip())
        return 0
    files = fragments()
    CHANGELOG.write_text(new, encoding="utf-8")
    for p in files:
        p.unlink()
    print(f"{len(files)} fragments moved into CHANGELOG.md")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

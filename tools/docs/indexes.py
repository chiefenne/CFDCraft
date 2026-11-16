"""Regenerates the number tables of the gallery and recipe index pages
from their topic pages.

    python -m tools.docs.indexes           rewrite docs/gallery.md and docs/recipes.md
    python -m tools.docs.indexes --check   exit 1 when they are not current

An entry on a topic page (docs/gallery/*.md, docs/recipes/*.md) is

    <a id="short-name"></a>
    ## 54. Title

and its number is the number of the pull request that adds it (a second
entry of the same pull request is 54b). From these the script writes two
things on each index page:

- the "All ... by number" table (one row per entry, sorted by number,
  with the anchor of the entry on its row);
- the last column of the page table (the numbers on each topic page).

Pull requests leave both tables alone: the workflow docs-index.yml runs
this script after every merge to main and commits the result, so that
parallel pull requests never conflict on the index pages.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"

ENTRY = re.compile(r'<a id="([^"]+)"></a>\s*\n##\s+(\d+)([a-z]?)\.\s+(.+?)\s*$', re.M)

# index page, topic folder, header of the by-number table
INDEXES = [
    ("gallery.md", "gallery", "| # | Example | Page |"),
    ("recipes.md", "recipes", "| # | Recipe | Page |"),
]


def entries(folder: str) -> tuple[list[tuple[int, str, str, str, str]], dict[str, list[str]]]:
    """(number, suffix, anchor, title, page file) of every entry, and the
    numbers per page in the order they stand on the page."""
    rows = []
    per_page: dict[str, list[str]] = {}
    for page in sorted((DOCS / folder).glob("*.md")):
        text = page.read_text(encoding="utf-8")
        rel = f"{folder}/{page.name}"
        found = ENTRY.findall(text)
        per_page[rel] = [f"{n}{s}" for _, n, s, _ in found]
        for anchor, n, s, title in found:
            rows.append((int(n), s, anchor, title, rel))
    rows.sort(key=lambda r: (r[0], r[1]))
    seen: dict[str, str] = {}
    for n, s, anchor, _, rel in rows:
        key = f"{n}{s}"
        if key in seen:
            raise SystemExit(f"number {key} used twice: {seen[key]} and {rel}")
        seen[key] = rel
    return rows, per_page


def render(index: str, folder: str, header: str) -> str:
    path = DOCS / index
    text = path.read_text(encoding="utf-8")
    rows, per_page = entries(folder)

    # the by-number table: header, rule, rows up to the first blank line
    lines = text.split("\n")
    try:
        h = lines.index(header)
    except ValueError:
        raise SystemExit(f"{index}: table header {header!r} not found") from None
    end = h + 2
    while end < len(lines) and lines[end].startswith("|"):
        end += 1
    # the page table above it: rows "| [Name](folder/page.md) | ... | numbers |";
    # its link texts name the pages in the by-number table too
    link = re.compile(r"^\| \[([^\]]+)\]\((" + re.escape(folder) + r"/[^)#]+\.md)\) \|")
    names: dict[str, str] = {}
    for i, line in enumerate(lines[:h]):
        m = link.match(line)
        if m and m.group(2) in per_page:
            names[m.group(2)] = m.group(1)
            cells = line.split(" | ")
            cells[-1] = ", ".join(per_page[m.group(2)]) + " |"
            lines[i] = " | ".join(cells)
    missing = sorted(set(per_page) - set(names))
    if missing:
        raise SystemExit(f"{index}: no row in the page table for {', '.join(missing)}")

    body = [
        f'| <a id="{anchor}"></a>{n}{s} | [{title}]({rel}#{anchor}) | [{names[rel]}]({rel}) |'
        for n, s, anchor, title, rel in rows
    ]
    lines[h + 2 : end] = body
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    check = "--check" in argv
    stale = []
    for index, folder, header in INDEXES:
        path = DOCS / index
        new = render(index, folder, header)
        if new != path.read_text(encoding="utf-8"):
            stale.append(index)
            if not check:
                path.write_text(new, encoding="utf-8")
    if check and stale:
        print("not current: " + ", ".join(stale) + " (run python -m tools.docs.indexes)")
        return 1
    print(("updated: " + ", ".join(stale)) if stale and not check else "indexes current")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

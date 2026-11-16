# Changelog fragments

Each pull request adds its changelog entry here as a file of its own,
instead of editing `CHANGELOG.md`, so that parallel pull requests never
conflict on the changelog.

- Name: `<PR number>-<short-name>.md`, for example `56-ci-tiers.md`.
- Content: the bullets for the "Unreleased" section, written as they will
  stand in `CHANGELOG.md` (bold lead, plain sentences). No heading.
- The release pull request runs `python -m tools.release.changelog`, which
  puts all fragments, sorted by number, at the top of "## Unreleased" in
  `CHANGELOG.md` and deletes them; then the section is renamed to the
  version as before (recipe 6).

`python -m tools.release.changelog --print` shows the "Unreleased" section
with the fragments as it would read now.

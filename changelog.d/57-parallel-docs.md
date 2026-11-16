- **Pull requests no longer conflict on bookkeeping.** Changelog entries
  are fragments in `changelog.d/` that the release collects
  (`tools/release/changelog.py`); new gallery and recipe entries take the
  number of their pull request; the number tables of `docs/gallery.md` and
  `docs/recipes.md` are regenerated after each merge by the workflow
  `docs-index.yml` (`tools/docs/indexes.py`); topic pages merge with git's
  union driver. CI checks the lock file of the Python crate first. Rules in
  CONTRIBUTING.md, "Working in parallel".

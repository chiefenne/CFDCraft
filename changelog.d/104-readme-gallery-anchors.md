- **README up to date, highlight images, gallery links that land.** The
  README lists what works on main ("What works today": turbulence models
  and wall treatments, heat transfer with passive scalars and thermal wall
  functions, the coupled algorithm, rotating machinery, …) under a banner
  and six highlight tiles linked to their gallery entries
  (`assets/readme/`). Gallery and README images carry their width and
  height (`tools/docs/images.py`, checked in CI), so a link to an entry far
  down a topic page no longer lands thousands of pixels too high while the
  images above it load. Recipe 33's anchor is `python-api` (it clashed with
  the page heading). New rule in CONTRIBUTING.md and the pull request
  template: every pull request keeps the README and the docs current.

# Contributing

## Documentation for every feature (D18)

A feature is not done without three pieces of documentation, in the same
form as in PyAeroFlow:

1. **A gallery example**: one picture with the commands that reproduce it,
   in [docs/gallery.md](docs/gallery.md).
2. **Case files**: one or more commented case files in `examples/cases/`,
   listed in the [case catalogue](docs/cases.md) with the mesh they need
   and where they are discussed.
3. **A recipe**: a short command sequence for the task, naming the
   case-file keys it uses, in [docs/recipes.md](docs/recipes.md).

### Where the pages go

Each of the three is an index page with a one-line summary per topic page;
the entries themselves live on topic pages, each with its own table of
contents:

| Index | Topic pages | Split by |
|---|---|---|
| [docs/gallery.md](docs/gallery.md) | `docs/gallery/*.md` | feature (benchmarks, meshes, scalar and pressure, flow, boundary conditions, post-processing, run control) |
| [docs/recipes.md](docs/recipes.md) | `docs/recipes/*.md` | task area (benchmarks, development, meshes, scalar and pressure, GPU, flow, boundary conditions, post-processing, run control) |
| [docs/cases.md](docs/cases.md) | `docs/cases/*.md` | physics, following the folders of `examples/cases/` |

When you add an entry:

- put it on the topic page that fits; start a new topic page when a new
  physics or topic begins (energy, turbulence, reacting flows, sprays,
  multiphase, ...) and add a row for it to the index page;
- gallery examples and recipes get the next free number; numbers are
  global and never reused, so "gallery 12" or "recipe 17" keeps its
  meaning in other pages, the changelog and commit messages;
- give each entry an `<a id="short-name"></a>` anchor right above its
  heading, a row in its page's table of contents and a row in the
  "by number" table of the index page; the anchor names are unique per
  section and stay valid when an entry moves to another page;
- link to an entry on its topic page (`gallery/flow.md#cavity`), not to the
  index; links to an anchor of the index page still work and land on the
  entry's row there;
- in the case catalogue, write the "Where" column as links
  (`[gallery 12](../gallery/flow.md#cavity)`), and add new keys to the
  page's key table; the complete key list `docs/keys.md` is generated from
  the code (`foehn schema`).

Keep a topic page below a few hundred lines; when one grows past that,
split it into two topic pages the same way.

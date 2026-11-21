# ArtCraft development method: comparison and challenges

Read on 2026-11-20. The maintainer asked that Foehncraft follow the way the
[ArtCraft](https://getartcraft.com/) Crafting Apps are built as closely as
possible, and that their practice be used to challenge the maintainer's own
decisions. This note records what they do, where Foehncraft matches, where it
differs, and the questions that follow. It changes no code and no decision;
each challenge waits for the maintainer.

**Framing.** ArtCraft rebuilds the features and the look and feel of
established desktop products. Foehncraft is measured against the industry
standards of general purpose CFD, but does not take over their flaws, in the
solver or in the later desktop app (L7): the maintainer decides what goes in.
So the recommendation is to adopt ArtCraft's **development method**, not to
copy the behaviour of industry-standard solvers.

**Sources** (public, read only): the Craft Anything page and the news page on
getartcraft.com; the READMEs of `storytold/craft-anything`,
`storytold/craft-launcher` and `storytold/artcraft`; in `craft-anything` the
files `AGENTS.md`, `CONTRIBUTING.md`, `docs/STATUS.md`, `docs/ROADMAP.md`,
`docs/tool-status.md`, `docs/extending.md`, `docs/pitfalls.md`; in
`storytold/craft-contrib` the `README.md` and the standards `app-agents.md`,
`clean-room.md`, `testing.md`, `pull-requests.md`, `claims.md`,
`licences.md`, `naming.md`. On 2026-11-21 also craft.wiki: the Craft
Anything guides (overview, getting started, install and update, workflows,
troubleshooting), the "Connect an AI assistant" guide, the automation and AI
agents pages of two sibling apps, and the about page (section 6). Craft
Anything was at version 0.6.1, the launcher at 0.19.2.

## 1. What ArtCraft does

| # | Practice | Where it is written |
|---|---|---|
| A1 | **One engine, thin front ends.** A pure-Rust engine behind the desktop app (egui), a headless CLI for scripts and agents, an MCP server and a WebAssembly build. "Both modes drive the same engine." | README, "Under the hood" |
| A2 | **Everything is a tool.** Each user-visible behaviour is a `ToolSpec` in a registry; the UI, the CLI, the MCP server and the AI agent dispatch by id. `<app>-cli commands` lists every tool, `<app>-cli mcp` serves them. 320 tools. | `AGENTS.md`, README "Built for agents" |
| A3 | **Headless first, flows are tools.** Anything a person can start from the UI is reachable from `Engine::call` and the CLI; orchestration, persistence and process control never live in the UI crate or the apps. Panels change state only through registered commands. | `ROADMAP.md` architecture rules 3; `app-agents.md` 14 |
| A4 | **Seams, not matches.** Traits plus a registry instead of enums; catalogues as typed data tables (TOML through `include_str!`); every file written and read back has a version, a JSON schema and a migration path; decisions read typed fields, never English text. | `ROADMAP.md` architecture rules 1, 2, 4, 5 |
| A5 | **Enforced layering.** Crates in layers L0 to L6, lowest first; only the top layer uses the GUI toolkit; `cargo xtask layers` fails the build on a wrong dependency. | `AGENTS.md` layout |
| A6 | **Confined agent entry point.** `<app>-cli mcp` must refuse to start without `--automation-read-root` and `--automation-write-root` and refuses file access outside them. | `app-agents.md` 18 |
| A7 | **Approval gates.** Before running a command of a kind not seen in the project, the app asks; "always allow" is per kind. "Review requests for actions outside the rules you set before they run." | README "Status", product page |
| A8 | **Evidence records.** Results come back with confidence, coverage and limitations; large results come back as a summary first and are paged. | README "Built for agents" |
| A9 | **Honest status.** `STATUS.md` says what works and *how it was checked* (machine, test counts, what did not run); `tool-status.md` grades every tool Verified, Partial, Fixed, Unverified or New. Rule: never present partial or placeholder behaviour as success, in the GUI, the CLI or automation results. | `STATUS.md`, `tool-status.md`, `app-agents.md` 15 |
| A10 | **Roadmap with ids.** Items `R<milestone>.<n>` with an owner, a status and small pending changes; `cargo xtask pending` lists open work; ids are never reused; design points marked *decision* wait for a maintainer. | `ROADMAP.md` |
| A11 | **Pitfalls file.** Past mistakes and their fixes in `docs/pitfalls.md` (about 100 kB), read by people and agents before work. | `AGENTS.md`, `CONTRIBUTING.md` |
| A12 | **One rule file for agents.** Every repository carries the same `AGENTS.md` (RFC 2119 wording, copied from a shared standard with a checksum), followed by a "This repository" section with layout, conventions and the exact checks to run. `CLAUDE.md` only points to it. | `AGENTS.md`, `craft-contrib` |
| A13 | **Clean room as a written standard.** Learn only from specifications, public documentation, papers and behaviour specifications; never read GPL, LGPL or AGPL source (copyleft tools may run as separate processes for tests); never copy code, wording, presets or pixels. Observing a proprietary product needs two separate roles: an observer who writes only a behaviour specification, and an implementer who never runs the product. Numbers measured from a proprietary product's output need their own decision before use. Every asset and fixture has recorded provenance. | `clean-room.md` |
| A14 | **No other company's names.** Feature names must not use another company's trademarked product or feature names. | `app-agents.md` 22, `naming.md` |
| A15 | **Testing standard.** A test stays only if one can name its invariant, a plausible change that breaks it, and that no other test catches it. Assert relations, not chosen constants or timings. A bug fix adds a test that fails without it; a crash fix adds its input to the crate's hostile-input sweep. Never delete, weaken or skip a test to get green; cuts go in their own PR listing every removed test. | `testing.md`, `app-agents.md` 27-32 |
| A16 | **Code rules.** `unsafe_code = "forbid"` workspace-wide (or only in named modules, each block with a `// SAFETY:` comment); non-test code does not panic; lints are fixed, not silenced; every external program goes through one process module; `cargo xtask secrets` before commit. | `AGENTS.md`, `app-agents.md` 13 |
| A17 | **Pull requests.** One issue per PR; moving code and changing behaviour in separate PRs; CI changes in their own PR; the description says what the user gains, the cause, the test, and what ran and what did not ("a check that did not run must not be reported as passed"). Changelog as one fragment per change. | `pull-requests.md`, `app-agents.md` 17, 19, 33 |
| A18 | **Operator approval.** An agent asks its operator before any step a new commit cannot undo: publishing, releasing, merging, deleting a branch or tag. | `app-agents.md` 4 |
| A19 | **Releases.** Builds for Windows, macOS, Linux (AppImage, Flatpak, deb, rpm, tarball, Arch), FreeBSD and the web; x86_64, ARM64, RISC-V; code-signed and notarized; `SHA256SUMS.txt` on every release; the CLI shipped standalone and inside the app. One launcher installs and updates every app and checks every download against the checksums. `<app>-cli doctor` checks the install and its external engines. | READMEs, launcher README |
| A20 | **Decisions by RFC.** Shared rules change only through an RFC in `craft-contrib`, with a 3 or 7 day comment period and approval by a maintainer who is not the author. | `craft-contrib` README |
| A21 | **For the later GUI.** All user-facing text through a translation system with English as the key; keyboard use without traps; accessibility tree with roles and names; contrast 4.5:1 for text; colours, sizes and fonts only from theme tokens. | `app-agents.md` 14, 16, 21 |

## 2. Where Foehncraft matches

- **Rust, one engine, CLI first** (D4, D11, D86): the CLI, the Python package
  (D61) and the future MCP server (L18) and desktop app (L7) are planned on one
  engine API. Same shape as A1.
- **Machine interface** (D77-D85): `foehn schema` (JSON Schema of the case
  keys, D77), `--json` on every subcommand with one JSON object on standard
  output and human text on standard error (D79), diagnostic codes with
  suggested fixes (D78), versioned machine formats (D80), exit codes (D64),
  read and write roots (D82), run limits (D83), byte-identical files with
  `SOURCE_DATE_EPOCH` (D84), provenance (D85). This is A4 (versioned files,
  structured results) and most of A6 already.
- **C-free validation crate** (D87) for a web build (P3): the same idea as
  ArtCraft's WebAssembly parts of the engine.
- **Evidence and honest results**: the cases check themselves against pass
  bands, the verification manual and `results.json` are generated (D68), and
  the README keeps "What works today" current. Close to A8 and A9.
- **Changelog fragments** (`changelog.d/`): the same as A17.
- **No other companies' product names** in the repository: the same as A14.
- **Clean room** in principle: methods from the published literature, no
  copied code (D10, `docs/verification/DATA_SOURCES.md` for reference data).
- **Registries for multiphysics** (D112-D122): regions, quantity and term
  registries as data, close to A4.

## 3. Where Foehncraft differs

| # | ArtCraft | Foehncraft now | Assessment |
|---|---|---|---|
| F1 | Every behaviour is a registered tool with one id across CLI, MCP and GUI; `commands` lists them (A2) | Hand-written `clap` subcommands; the engine API of D86 is a set of functions; no catalogue | Biggest structural gap for the agent-first rule. Cheapest to fix before the MCP server and the GUI exist. |
| F2 | The MCP entry point refuses to start without read and write roots (A6) | Roots are optional: confinement applies only when a root is given (D82) | Fine for a person at a terminal, too open for an agent. |
| F3 | One `AGENTS.md` in the repository with all rules and the exact checks (A12) | `CONTRIBUTING.md` holds the documentation rules; many working rules (merging, progress notes, iteration budgets, image rules, naming) live in the maintainer's project notes outside the repository | An agent or contributor who clones the repository does not see them. |
| F4 | Pitfalls in the repository (A11) | Lessons (coupled-solver experiments not to repeat, skewed-tet instability, observers that perturb a run, iteration budgets) are in project notes outside the repository | Same as F3. |
| F5 | Plans with ids and pending changes in the repository (A10) | Decisions D/L/P/X/O in `docs/requirements.md` (good), but stage plans such as the stage 7 plan live outside the repository | Partly the same as F3. |
| F6 | Status per tool: Verified, Partial, Unverified, and how each was checked (A9) | Verification manual per case; no per-feature grade; partial features are described in prose | Example: GPU paths verified on one card type only; k-ε-ζ-f₀ near-wall k not grid-convergent. These are known and written down, but not as a grade that a result summary can carry. |
| F7 | Never read copyleft source; observer and implementer separated; measured output of proprietary products needs its own decision (A13) | Open-source code may be read where the licence permits, never copied | Much open CFD source is under copyleft licences. "Licence permitting" is weaker than A13 and invites shaped-like-it code. |
| F8 | `unsafe` forbidden, or only in named modules with `// SAFETY:` comments; no panics in non-test code (A16) | `unsafe` in the kernels, CUDA, finite-volume and bench crates; no workspace `forbid`; about 300 `unwrap`/`panic!` in library sources | `unsafe` is needed for the GPU and SIMD paths; the rule can still name the modules. Panics matter once an agent drives the solver through MCP: it should get a diagnostic code, not a dead process. |
| F9 | Rust only, `cargo xtask` for repository tasks (A5, A16) | Python for documentation indexes, image checks, the changelog, the verification manual and plots (`tools/`); C libraries allowed where they cost no speed (D9: MPI, HDF5, METIS) | The C libraries are a justified difference (D9). The Python development tools are not required by any decision. |
| F10 | Layering checked by a tool (A5) | Crate structure exists, not checked | Low cost, protects D86 and D87. |
| F11 | Tests judged by the keep rule; hostile-input sweep for every reader (A15) | Verification against analytic solutions and published data with pass bands; unit tests | For physics, exact reference values *are* the requirement, so the relational-assertion rule does not transfer to verification cases. The keep rule and the hostile-input sweep do transfer, to mesh and case readers (`.msh`, VTK, multi-zone, case TOML). |
| F12 | Agent asks before merging or releasing (A18) | Maintainer rule: the agent merges itself when CI is green, and bumps and tags versions | A deliberate maintainer decision for a one-maintainer project; recorded here as a difference, not a challenge. |
| F13 | Signed binaries for every platform, checksums, launcher, `doctor` (A19) | Tags and GitHub releases with notes; no binaries; build from source | Fine before v1.0; an agent today has to build from source. |
| F14 | Approval per kind of action (A7) | The solver starts no processes (D82), so there is nothing to approve inside it; but a run can cost hours of GPU time | The relevant gate for a CFD solver is cost, not commands. |
| F15 | Large results as summary first, then paged (A8) | Result summary (D62) plus files on disk | Needed for L18: a field of 100 million cells cannot go through MCP. |

## 4. Challenges for the maintainer

Each item is a question for a decision, ordered by value for the agent-first
rule. None is started.

1. **Tool registry before the GUI (F1).** Make every CLI command a registered
   engine tool with an id, typed input (from the D77 schema) and typed output
   (D80), and add `foehn commands` (JSON catalogue). Then the MCP server (L18),
   the Python package and the desktop app dispatch the same ids, and nothing
   can exist only in the GUI. Decide now whether D86 becomes "one registry of
   tools" instead of "one engine API".
2. **Mandatory roots for agents (F2).** `foehn mcp` (L18) refuses to start
   without `--read-root` and `--write-root`; the plain CLI keeps today's
   behaviour.
3. **`AGENTS.md` in the repository (F3, F4, F5).** One file with the
   standing rules written as MUST/SHOULD, the layout, and the exact checks to
   run; a `docs/pitfalls.md` with the lessons; stage plans moved into
   `docs/requirements/`. Today the rules work only inside this project's own
   tooling.
4. **A cost gate instead of a command gate (F14).** An agent-started run whose
   D81 estimate exceeds a limit the maintainer sets (cell-hours, GPU-hours,
   output volume) needs approval; `[run.limits]` (D83) is mandatory for runs
   started through MCP. This also writes down the existing rule to look at a
   steady run after 2000 to 3000 iterations.
5. **Verification grades in results (F6).** A grade per feature (verified,
   partial, unverified, and how it was checked) in a generated status page,
   and in every result summary a `limitations` list (for example: model
   outside its validated range, mesh study not grid-converged, backend not
   verified on this GPU). The rule "never present partial as success"
   applies to `--json` output too.
6. **A stricter clean-room rule (F7).** Do not read copyleft CFD source at
   all; learn from papers, documented formats and public reference data. For
   the GUI later: an observer writes a behaviour specification of an
   industry-standard workflow (layout, steps), the implementer works only from
   it, and each specification has a column "flaw, not adopted" that the
   maintainer fills in. Comparisons with numbers produced by industry-standard
   solvers need a decision of their own before they are used, as ArtCraft
   requires.
7. **Panics and `unsafe` (F8).** Forbid `unsafe` in all crates except the
   named GPU and SIMD modules, each block with a `// SAFETY:` comment; deny
   `unwrap` and `panic!` in the case, mesh-reader and engine crates step by
   step, so that agents get diagnostic codes.
8. **`cargo xtask` for repository tasks (F9, F10).** Move the Python
   development tools (indexes, images, changelog, verification manual) to a
   Rust `xtask` and add a layering check. Python stays as a user front end
   (D61) and for plots in studies.
9. **Hostile-input sweep for readers (F11).** Every mesh and case reader gets
   a corpus of broken inputs that must end in a diagnostic, not a panic, and
   adopt the keep rule for unit tests; verification cases keep their exact
   reference bands.
10. **PR description "what ran and what did not" (A17).** Add this line to
    the pull request template; it matches the maintainer's rule to report
    test outcomes faithfully.
11. **Binaries and `foehn doctor` (F13).** From v0.7 or v1.0: CLI binaries for
    Linux, macOS and Windows with `SHA256SUMS.txt`, and `foehn doctor` that
    checks the build features, the GPU, MPI and HDF5. An agent can then pin a
    version without building.
12. **Paged results for MCP (F15).** Design the L18 result types as summary
    first, with paging for monitors, field statistics and probes, never raw
    fields.
13. **For the desktop app (L7), later.** Adopt A21 from the first GUI commit:
    translation keys in English, keyboard use, accessibility tree, theme
    tokens, and "panels change state only through registered tools".

## 5. From craft.wiki (2026-11-21)

What the user guides add to sections 1 to 4:

| # | Practice | Bearing on Foehncraft |
|---|---|---|
| W1 | **Approvals per kind of action.** Auto mode is on by default, but the builder never deletes outside its project, changes system settings or publishes. Per project, each kind of action is set to Default, Allow, Ask or Refuse: reading the web, other commands, files outside the project, pushing and publishing, deleting outside the project, system changes, other tools. "Allow once" or "always for this project". | Sharpens challenge 4: a short list of kinds for a solver driven by an agent (start a run, exceed the cost estimate, write outside the write root, overwrite or delete earlier results, publish), each with the same four settings. |
| W2 | **Spend shown and bounded.** Each project shows its run time, tokens and an estimated cost; paid extras have a budget per project, raised on request. | A GPU-hour or cell-hour budget per study, shown in the run record next to the D81 estimate and spent against `[run.limits]` (D83). |
| W3 | **"Good enough" and refining rounds** are settings, not hidden defaults. | Same idea as our iteration budgets and the rule to diagnose a steady run at 2000 to 3000 iterations: make the budget a visible setting of a study. |
| W4 | **Action journal and replay.** The app lists every command run in a session, by the person or the assistant; `<app>-cli run --cmd 'id={json}'` or `--script` replays commands without a window. | New challenge 14. |
| W5 | **Live control.** `<app> --control <port>` listens only on the local computer and needs an access token; `<app>-cli mcp --connect` lets an assistant work in the open app, `send <port>` sends one request. Notes typed during a build are read at the next step. | Matches L19 (control at step boundaries); add the token. |
| W6 | **Every app has the same guide set:** overview, getting started, install and update, workflows, AI agents (or automation), troubleshooting with known issues linked by issue number; plus one step-by-step "Connect an AI assistant" guide for people who are not programmers. The wiki is a public repository with reviewed edits. | New challenge 15. |
| W7 | **Safety advice to users:** give the assistant a folder of its own, keep copies, read what it asks before approving, and treat text inside files and web pages as information, not orders. | Belongs in our agent guide; case files, meshes and reference data can carry text too. |
| W8 | **`check` subcommand** recomputes stored files and compares them with the values saved in them. | Matches our self-checking verification cases. |

## 6. More challenges from the wiki

14. **Command journal (W4).** Every engine call of a session, by a person or
    an agent, goes into the run record (D69) as typed tool calls, and
    `foehn run --script` replays them. Then a study made by an agent can be
    rerun by a person and the other way round. Builds on challenge 1.
15. **One guide set (W6, W7).** Alongside gallery, cases and recipes (D18):
    `docs/agents.md` (connect an assistant to `foehn mcp`, step by step, with
    the safety advice of W7), a troubleshooting page with known issues
    linked to their issues, and an install and update page.
16. **Approval kinds and budgets (W1, W2, W3).** Turn challenge 4 into a
    table of action kinds with Default, Allow, Ask or Refuse, a budget per
    study, and the iteration budget as a visible setting.
17. **Token on the control channel (W5).** L19 listens on the local computer
    only and requires an access token readable only by the user.

## 7. Watching ArtCraft

A weekly check (Monday morning, Vienna time) reads their news page, the
Craft Anything page and wiki, the craft-contrib standards and the releases of the
repositories above, and reports only changes that affect a decision here
(new shared helper tools, changed agent rules, new standards). This note is
updated when such a change is adopted.

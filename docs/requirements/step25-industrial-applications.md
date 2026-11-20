# Step 25 – Industrial applications: gallery, cases and recipes (proposal)

Requirements step 25 (2026-11-20). The maintainer asked for an
industrial-level gallery with its recipes and cases, **in parallel** to
the validation and feature gallery, and for a structure first. This page
is that structure. Nothing here is implemented yet: the pages, tools and
case files below are written only on the maintainer's go (project rule).

The decisions to take are collected in [section 10](#10-decisions-for-the-maintainer).

Contents

1. [Short answer](#1-short-answer)
2. [Why a second track](#2-why-a-second-track)
3. [Where the entries live](#3-where-the-entries-live)
4. [Numbering](#4-numbering)
5. [The entry template](#5-the-entry-template)
6. [Trust without reference data: the checks](#6-trust-without-reference-data-the-checks)
7. [Candidate cases on main today](#7-candidate-cases-on-main-today)
8. [Cases that wait for later physics](#8-cases-that-wait-for-later-physics)
9. [CI, run time and repository size](#9-ci-run-time-and-repository-size)
10. [Decisions for the maintainer](#10-decisions-for-the-maintainer)

## 1. Short answer

- A new track called **Applications**, with its own index page
  `docs/applications.md` and topic pages by industry sector in
  `docs/applications/`, next to (not inside) the feature gallery.
- **One self-contained entry per application**: problem, geometry and mesh
  source, setup, run cost, 1–3 field images, engineering results, the
  checks that make the result trustworthy, and the limits. The run commands
  sit in the entry; generic workflows (CAD to mesh to run to report) go to
  one new recipe topic page.
- Case files in a new folder `examples/applications/<slug>/` (case file,
  mesh script or `.geo`, post-processing script), listed on a new case
  catalogue page `docs/cases/applications.md`.
- **Numbers are the pull-request numbers**, as for the gallery, written
  "application 131": the noun keeps them apart from "gallery 131".
- **No reference data needed**, but every entry passes a fixed set of
  checks (balances, iterative convergence, mesh sensitivity with the GCI
  of D67, wall resolution for the chosen wall treatment, plausibility
  against an engineering correlation where one exists) and shows a
  **check card** with the results.
- **First batch of six** (section 7): orifice plate, ribbed cooling channel,
  thermal mixing tee, stirred tank in 3D, double pipe bend, offset-strip-fin
  module. All run with features on main.
- **CI stays fast**: pull requests run nothing new except a docs check; a
  coarse smoke variant of every application runs on main, nightly and with
  the `full-ci` label; production runs happen on the maintainer's iMac or
  in a cloud session, never in CI.

## 2. Why a second track

The feature gallery answers *does it work, and how accurately*: small
canonical cases, reference data, pass bands, one feature per entry. An
industrial user asks a different question: *can I solve my problem with
it, how long does it take, and how far can I trust the number I get*.
That needs

- realistic 3D geometry from a CAD-like source (Gmsh `.geo`/OpenCASCADE or
  a Python mesh script), not a unit box;
- several features at once (Gmsh tetrahedra or polyhedra, a RANS model
  with a wall treatment, MRF, passive scalars, periodic modules);
- engineering outputs (loss coefficient, power number, head, mixing index,
  Nusselt number, uniformity index) instead of profiles against a table;
- run cost and mesh size stated up front;
- an honest statement of what is trusted and what is not.

Mixing both into one gallery would blur the meaning of an entry: a
validation entry is a claim backed by data, an application entry is a
demonstration backed by checks. Keeping two tracks keeps both claims
clean, and an application may link to the validation entries its models
rest on (for example the ribbed channel to the channel against DNS,
gallery 68, and the heated pipe, gallery 85).

## 3. Where the entries live

Proposed layout, following the conventions of CONTRIBUTING.md and PR #42:

| Place | Content |
|---|---|
| `docs/applications.md` | index page: what the track is, how it differs from the gallery, the trust legend (section 6), the topic-page table, the "by number" table (generated) |
| `docs/applications/process-piping.md` | pipes, bends, tees, orifices, valves, manifolds, static mixers |
| `docs/applications/rotating-machinery.md` | stirred tanks, fans, pumps |
| `docs/applications/thermal-management.md` | cooling channels, heat-exchanger modules, heat sinks, impingement cooling |
| `docs/applications/hvac-buildings.md` | room ventilation, ducts, air distribution |
| later pages | energy conversion (reacting flows), sprays, multiphase, external aerodynamics, when the physics exists |
| `docs/cases/applications.md` | the case catalogue page of `examples/applications/`: files, mesh source, cells, run cost, where discussed |
| `docs/recipes/applications.md` | generic workflows: geometry to Gmsh mesh with boundary layers, mesh checks before a long run, a production run on the iMac with restarts and `.result_N`, the check card (section 6), images for the docs |
| `examples/applications/<slug>/` | `case.toml` (production), `case_ci.toml` (coarse smoke), `mesh.geo` or `mesh.py`, `post.py` (engineering numbers, check card, images), `README.md` (three lines) |
| `docs/figures/applications/` | the images of the entries (section 9 for sizes) |

The README gets one line and one picture row once the first entries exist;
`docs/gallery.md` gets one sentence pointing to the new track. The
`tools.docs.indexes` script gets a third index (`applications.md`,
`applications/`), so the by-number table is regenerated after a merge as
for the gallery and recipes, and pull requests never touch it.

**Alternative considered:** a tag or a section "Industrial" on the existing
topic pages. Rejected: the existing pages are split by feature, an
application uses many features, so it would have no natural home, and
the gallery index would mix two kinds of claims.

## 4. Numbering

- An application entry takes the number of the pull request that adds it,
  `<a id="slug"></a>` + `## 131. Title`, exactly the form the index script
  already parses; a second entry of the same pull request is 131b.
- In text it is "application 131" (short "app 131"), as "gallery 54" and
  "recipe 54" already share number 54. The numbers are global and never
  reused.
- Its recipe on `docs/recipes/applications.md` and its case rows use the
  same number, so application 131 has recipe 131 and its case files.
- A later improvement of an entry (finer mesh, a new model) edits the entry
  and notes the pull request in its history line; it keeps its number.

## 5. The entry template

Every entry has the same headings in the same order, so a reader can
compare entries at a glance. Sketch for one entry:

```markdown
<a id="orifice-plate"></a>
## 131. Orifice plate in a water pipe: discharge coefficient and pressure loss

<img src="../figures/applications/orifice_speed.jpg" alt="..." width="800" height="300">

**Problem.** What an engineer wants to know, in two to four sentences
(sizing a flow meter: Δp at the design flow, permanent pressure loss,
recovery length downstream).

| | |
|---|---|
| Sector | process and piping |
| Geometry | D = 100 mm, β = d/D = 0.6, plate 3 mm with a 45° bevel; source `examples/applications/orifice_plate/mesh.py` |
| Mesh | axisymmetric wedge, 42 000 cells (fine), y⁺ ≈ 1 at the plate; coarse/medium/fine 10 k/21 k/42 k |
| Fluid | water at 20 °C from the material database |
| Boundary conditions | mass-flow inlet (fully developed profile from a periodic run), static-pressure outlet 30 D downstream, walls |
| Models | k-ω SST, `low-re` walls; steady, coupled solver |
| Run cost | 4 min on the iMac (8 cores), 220 iterations; smoke variant 20 s |
| Confidence | **checked** (all checks of section 6 pass; plausibility against the ISO 5167-2 discharge coefficient within 2 %) |
| Case files | `applications/orifice_plate/case.toml`, `case_ci.toml` |

### Results
The engineering numbers in a table, with the GCI band from the three meshes.
1–3 field images (speed with streamlines, pressure, a scalar or vorticity).

### Check card
The table of section 6, with values and pass/fail.

### What to trust, what not
Two to five bullets: which numbers are robust, which depend on the model
(reattachment length, peak wall heat flux), what the model cannot see.

### Run it
The commands from the repository root (mesh, run, post), production and
smoke variant.

### Limits and next steps
What changes with later physics (CHT, buoyancy, LES), and the validation
entries this case rests on.
```

The field images follow the colourful-fluid-dynamics rule (velocity,
pressure, vorticity, scalar) and the rule against symbol-like images
(CONTRIBUTING.md): no four-fold patterns, rotors with three or five
blades, every image looked at in full size and as a thumbnail.

## 6. Trust without reference data: the checks

A validation entry is trusted because it matches data. An application
entry is trusted because it passes checks that need no data. The checks
use tools already on main.

| Id | Check | Tool on main | Pass rule (proposal) |
|---|---|---|---|
| C1 | Mass balance | ledger (D88), `result.json` | net imbalance below 1e-6 of the throughput |
| C2 | Energy balance | kinetic-energy ledger (D90 part 4), scalar ledger | boundary power = dissipation + numerical dissipation within `[balance] tolerance_energy`; numerical dissipation not negative; scalar balance below 1e-5 of the flux |
| C3 | Iterative convergence | monitors in the log, D69 health numbers | engineering quantities change by less than 0.1 % over the last 10 % of the iterations; the residual alone is not used (see the steady iteration budget) |
| C4 | Mesh sensitivity | `foehncraft_flow::uncertainty` (D67, gallery 25) | three meshes with a ratio of at least 1.3; GCI of each reported quantity below 5 % (or the band is reported and the entry says "sensitive") |
| C5 | Wall resolution | `y-plus`, `turbulence-wall` derived fields | y⁺ in the range of the chosen wall treatment on at least 95 % of the wall area; the rest shown |
| C6 | Mesh quality | quality report and D50 remediation | no cell left unremediated; skewness and non-orthogonality histograms in the entry |
| C7 | Plausibility | `post.py` | against an engineering correlation or handbook value where one exists (discharge coefficient, Dean-flow loss, power number, Nusselt correlation, j and f factors); labelled "plausibility", **not** validation; the deviation is reported, not banded |
| C8 | Independence of the domain | one extra run | outlet moved 50 % further: quantity changes by less than 1 % |
| C9 | Health | divergence report, health numbers | no warning left open; every remaining warning explained in the entry |

Each entry shows one of three trust levels in its table and on the index:

- **checked**: C1–C6 and C9 pass, C7 within the uncertainty of the
  correlation where one exists;
- **checked, model-sensitive**: the checks pass, but the result depends on
  the turbulence model by more than the GCI (shown with a second model);
- **demonstration**: runs and looks right, checks incomplete (for example a
  single mesh because the fine mesh is too large for the iMac); allowed for
  at most the first version of an entry.

`post.py` writes the check card as Markdown and as `checks.json`, so the
card in the docs is generated, not typed. A small `tools/check/app_card.py`
(written with the first entry) holds the common part.

## 7. Candidate cases on main today

Only features on main are used: Gmsh meshes (tetrahedra, polyhedra,
prism layers from Gmsh), the built-in mesher with annulus and revolved
meshes, wedges, periodic pairs with heat flux (D109), the RANS models of
stage 6 with their wall treatments, rough walls, QCR and curvature
correction, MRF and the sliding mesh, passive scalars with sources,
transient statistics, forces, torque and machine performance, derived
fields. Temperature is a **passive scalar** until stage 7, so heat transfer
is valid only for small temperature differences (constant properties, no
buoyancy); each such entry says so.

Cells and run costs are first estimates from the cases run so far (for
example 16 000 cells of the turbulent stirred vessel in 10 minutes on four
cores); the first entry measures the cost per cell and iteration on the
iMac and this table is corrected.

### Priority 1: the first batch

| # | Case | Sector | Features used | Engineering result | Plausibility (C7) | Cells (fine) | Cost (iMac, estimate) |
|---|---|---|---|---|---|---|---|
| A1 | **Orifice plate** in a water pipe, β = 0.6 | process, metering | wedge, SST `low-re`, coupled, material database | discharge coefficient, permanent pressure loss, recovery length | ISO 5167-2 discharge coefficient | 40–60 k | minutes |
| A2 | **Ribbed cooling channel**, periodic module, ribs at 45° | gas-turbine blade cooling | periodic pair with mass flow and heat-flux scalar (D109), ζ-f₀ and SST, Kays–Crawford Pr_t | friction factor and Nusselt enhancement over a smooth duct, wall heat-transfer map | Dittus–Boelter / Gnielinski ratio ranges of rib-roughened channels | 300–600 k | 2–5 h |
| A3 | **Thermal mixing tee**, hot branch into cold main pipe | power plants, process | Gmsh prisms + tets, SST, passive scalar, transient statistics (mean, RMS) | mixing index along the pipe, wall temperature RMS (thermal-fatigue screening) | jet momentum ratio regimes (wall jet, deflecting, impinging) | 0.5–1 M | steady 2–4 h; transient overnight |
| A4 | **Stirred tank in 3D**, three-blade hydrofoil or five-blade pitched turbine, three baffles | chemical and food process | built-in or Gmsh rotor zone, MRF steady, then sliding, SST, torque and power (D184), scalar for blending | power number, pumping number, blend time from a tracer | power-number ranges of the impeller type | 300–800 k | MRF 2–4 h; sliding overnight |
| A5 | **Double pipe bend** out of plane (two 90° bends, R/D = 1.5) | piping, flow meters downstream | Gmsh hex-dominant or polyhedra, SST and QCR, derived fields | loss coefficient per bend, swirl angle and profile distortion downstream (flow-meter installation) | handbook bend loss coefficients | 400–800 k | 2–4 h |
| A6 | **Offset-strip-fin module** of a plate-fin heat exchanger | heat exchangers, HVAC | built-in or Gmsh periodic module, laminar and transitional Re, heat-flux scalar | Colburn j and Fanning f over Re (5 points) | Manglik–Bergles correlation | 50–200 k per point | minutes to 1 h per point |

A1 and A6 are cheap enough to make the first entry and to measure the
cost per cell. A2 and A3 show the near-wall models and the scalar, A4
rotating machinery, A5 secondary flow with QCR.

### Priority 2

| # | Case | Features used | Engineering result | Cells | Cost |
|---|---|---|---|---|---|
| A7 | Pin-fin heat sink in a duct (the 1.3 M pin-fin mesh of the README picture, turbulent) | Gmsh, SST or ζ-f₀, heat-flux scalar | pressure drop, Nusselt number per row | 1.3 M | overnight |
| A8 | Ventilated office: ceiling diffuser and exhaust, isothermal | Gmsh, SST, scalar with unit source = **mean age of air** | age of air, air-change effectiveness, draught risk from speed and turbulence | 0.5–1 M | 3–6 h |
| A9 | Intake manifold or distribution header with five runners | Gmsh, SST, mass-flow inlet, pressure outlets | flow split and uniformity index, loss per runner | 0.5–1 M | 3–6 h |
| A10 | Axial fan in a duct, five or seven blades | Gmsh rotor zone, MRF, `[[output.machines]]` | fan curve at three operating points: pressure rise, efficiency | 0.5–1.5 M | 3 × 3–6 h |
| A11 | Impingement cooling array (row of jets on a heated plate) | builds on gallery 90, SST and ζ-f₀, scalar | Nusselt map and average, cross-flow effect | 0.5–1 M | 3–6 h |

### Priority 3 (heavier, after the first batch has measured the cost)

| # | Case | Note |
|---|---|---|
| A12 | Centrifugal pump with volute, five- or seven-blade impeller | MRF frozen rotor at three flows, head and efficiency (D184); transient sliding for the blade-passing pressure pulsation later; 1–3 M cells |
| A13 | Butterfly valve at three openings | loss and torque coefficients; 0.5–1 M cells per opening |
| A14 | Static mixer (helical elements) in a pipe | scalar mixing (coefficient of variation) per element; laminar and turbulent |

## 8. Cases that wait for later physics

Listed so that the sectors exist on the index from the start, each marked
"after stage N":

| Case | Needs |
|---|---|
| Heat sink with conduction in the fins, cold plate of a battery pack, electronics box | conjugate heat transfer (stage 7, D157–D158) |
| Data-centre aisle, displacement ventilation, heated room | buoyancy (stage 7, D153–D155) |
| Shell-and-tube heat exchanger, filters, catalysts, packed beds | porous regions (REG001 lifts with stage 7) |
| Cyclone separator | strong swirl needs a Reynolds-stress model (L47); with eddy-viscosity models it would be published only as a counter-example |
| Burners, furnaces | reacting flows (L2) |
| Spray cooling, fuel injection, scrubbers | sprays (L12) |
| Bubble columns, sloshing tanks | multiphase (L3) |
| Vehicles, wings, buildings in wind | external aerodynamics (after milestone 1) |
| Valve noise, fan tonal noise | LES (L5) |

## 9. CI, run time and repository size

**Pull requests stay as fast as today.** Nothing new runs in a pull
request except the existing docs checks (`tools.docs.indexes --check`,
`tools.docs.images --check`), which cover the new pages, and a cheap
parse of every application case file (`foehn case upgrade --check`,
seconds).

**Smoke variants on main, nightly and `full-ci`.** Every application has
a `case_ci.toml` on a coarse mesh (target below 60 s on the Linux runner)
that runs a few hundred iterations and asserts only robustness: no
divergence, balances C1–C2, engineering quantity inside a broad band
(± 20 % of the production value). They form a new set `applications` in
`tools/check/ci_cases.py`, run like `long` (main, nightly, `full-ci`),
with a budget of 5 minutes for the whole set. Gmsh is installed (`pip
install gmsh`) only in that job; meshes are generated, never stored.

**Production runs never in CI.** They run on the maintainer's iMac (long
cases, overnight) or in a cloud session for runs up to about an hour. The
run writes `result.json` and `checks.json`; the pull request of the entry
commits `checks.json` (a few kB) and the images, not fields or meshes. A
docs check compares the card on the page with `checks.json`, so the page
cannot drift from the run that made it.

**Repository size.** `docs/figures/gallery` is 11 MB today. Proposal for
the new track: JPEG for field images, at most 3 per entry, each below
300 kB at 1600 × (≤ 900) pixels, so 14 entries add about 10 MB; no
`.msh`, no `.vtkhdf`, no restart files in the repository.

**Run-time statement.** Every entry states cells, iterations or steps,
wall time and machine (iMac i9-9900K, 8 cores, or the cloud machine), so
readers can scale it to their hardware.

## 10. Decisions for the maintainer

| Id | Question | Proposal |
|---|---|---|
| Q1 | Name of the track | "Applications" (alternatives: "Industrial cases", "Engineering gallery") |
| Q2 | Separate track vs a tag in the existing gallery | separate track, section 3 |
| Q3 | Numbering | the pull-request number, read as "application N", section 4 |
| Q4 | Trust levels and checks C1–C9 | as section 6; the pass rules are proposals |
| Q5 | First batch | A1–A6, starting with A1 (orifice plate) to measure the cost, then A6, A2, A3, A4, A5 |
| Q6 | Stirred tank impeller | three-blade hydrofoil or five-blade pitched turbine with three baffles, following the image rule (four baffles are common in industry, but their flow images are four-fold) |
| Q7 | CI | smoke variants in a new `applications` set on main, nightly and `full-ci`; production runs off CI |
| Q8 | Where production runs happen | the iMac for runs over an hour, the cloud for shorter ones; the maintainer sends `result.json` and the fields for the images |

On the maintainer's go, the first pull request adds the empty structure
(index page, topic pages, case and recipe pages, the third index in
`tools.docs.indexes`, the `applications` CI set and `app_card.py`) together
with application A1, so the structure is tested by a real entry.

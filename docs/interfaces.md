# Interfaces

What scripts, studies and agents can rely on when they drive `foehn`: the
command line, the exit codes and signals, the files a run writes, the
diagnostic codes and the engine API behind the command line. Everything
here is versioned: a file or a JSON object carries `format` and
`format_version`, and a change that breaks a reader raises the major
version. The case-file keys are in [keys.md](keys.md) (generated from the
code by `foehn schema --markdown`).

## Command line (D62, D80)

| Option | Effect |
|---|---|
| `--json` | One JSON object on standard output at the end (the result summary, the check report, the command's report or an error object). Log lines go to standard error. |
| `--quiet` | No log lines on the terminal (the files are written as usual). |
| `--set KEY=VALUE` | Changes one case key before the run, repeatable. The value is read as TOML (`--set fluid.viscosity=1e-3`, `--set 'time.adaptive=false'`); if it is not valid TOML it is taken as a string (`--set run.backend=cuda`). Applied in order, after the case file and before the validation. The overrides are recorded in the result and in the run record and survive a re-read of the case file (D40). |
| `--restart latest` | Continues from the newest restart file in the output directory: the largest step or iteration in the file wins, `restart.h5` wins a tie. No file is `RUN011`. |
| `--read-root DIR`, `--write-root DIR` | Confines the files of a run (D82): every file read (case, mesh, restart) must lie below a read root and the output directory below a write root. A path outside is `CASE007`. Without roots nothing is confined. |
| `--check` | Validates the case, checks the boundary conditions and estimates the run (D81) without running it. Exit code 0 when the case can run, 1 when it has an error. |
| `foehn schema` | The JSON Schema of the case file (draft 2020-12, with `x-unit` and `x-runtime`). `--markdown` prints the key reference instead. |

## Exit codes (D64)

| Code | Status | Meaning |
|---|---|---|
| 0 | `finished` | The run ended as planned and every check of the case passed. |
| 1 | `error` | The case could not be read or run, or the solver failed (diverged, a file could not be written). |
| 2 | `failed-checks` | The run ended, but a check failed: not converged (`RUN006`), an order, error, reference or benchmark band, a mass balance (`RUN007`). |
| 3 | `stopped` | Stopped on request before its end: the `.stop` file, a signal or a limit of `[run.limits]`. The restart file and the output of the last iteration or step are written; `--restart latest` continues bitwise as if the run had not stopped. |

## Signals (D64, unix)

| Signal | Effect |
|---|---|
| `SIGUSR1` | Writes `restart.h5` after the current iteration or step; the run goes on. |
| `SIGTERM`, `SIGINT`, `SIGUSR2` | Stops after the current iteration or step with the restart file and the final output (exit 3, `RUN002`). |
| second `SIGINT` | Ends at once (exit 3) without writing anything more. |

A batch system that sends SIGTERM before it kills a job (or SIGUSR2 a set
time before the end) thus gets a clean stop that the next job continues
with `--restart latest`.

## Limits (D83)

`[run.limits]` stops a run with exit code 3 when it reaches a budget:

| Key | Unit | Code |
|---|---|---|
| `max_wall_time` | s | `RUN003` |
| `max_steps` (transient) | | `RUN004` |
| `max_iterations` (steady, outer iterations) | | `RUN004` |
| `max_output` (bytes in the output directory) | B | `RUN005` |

The limits may change during a run (a re-read of the case file, D40).

## Files of a run

| File | Format | Contents |
|---|---|---|
| `result.json` | `foehncraft.result` | The result summary (D62): status, exit code, stop, overrides, provenance, diagnostics, one entry per solve (mesh hash, iterations, residuals, mass imbalance, patches, coefficients, errors, transient summary), observed orders, and the output files with sizes and SHA-256 hashes. No timings: two identical runs give the same bytes. |
| `run.jsonl` | `foehncraft.run_record` | The run record (D69), one JSON object per line, flushed at every log line, so it can be followed while the run goes on (`tail -f`). |
| `flow.log`, `flow.md`, `flow.json` | | The log and the report for people, with timings. |
| `fields.vtkhdf`, `fields_N.vtkhdf` | VTKHDF 2.x | Fields, with the root attributes `foehncraft_format`, `foehncraft_format_version`, `foehncraft_version`, `case`, `units` and `step`/`time` or `iteration`. |
| `restart.h5` | Foehncraft restart | With `foehncraft_version` and `units`. |
| `*.csv` | `foehncraft-<kind>` | Residuals, history, forces, probes, patches, wall and line data, each with a comment preamble (below). |
| `case.read_toml_N.toml` | TOML | The settings in force after a re-read (D40). |

### Rows of `run.jsonl`

| `type` | Fields |
|---|---|
| `header` | `format`, `format_version`, `case`, `date`, `provenance`, `threads`, `overrides`, `case_file`, `settings` (the resolved case) |
| `iteration` | `iteration`, `momentum` (3), `continuity`, `bicgstab` (3), `cg` |
| `step` | `step`, `time`, `dt`, `cfl`, `momentum` (3), `continuity`, `outer`, `cg` |
| `event` | `code`, `severity`, `message`, `at` (iteration or step) |
| `timings` | `seconds`, `cells`, `steps`, `ns_per_cell_step`, `solves` (per solve: `mesh`, `device`, `timings`, `ns_per_cell_iteration`, `ns_per_cell_step`) |
| `result` | `result` (the object of `result.json`) |
| `end` | `status`, `exit_code` |

A reader skips a `type` it does not know; new types and new fields do not
raise the major version.

### Self-describing files (D75)

Every CSV starts with comment lines that a reader may skip (`#`):

```
# format foehncraft-history 1; Foehncraft 0.5.0; case cylinder_re100; step 400; time 40 s
# units: time s, Fx N, Fy N, Fz N
time,Fx,Fy,Fz
```

### Byte-identical output (D84)

Two runs of the same case with the same build give byte-identical files,
on any thread count and on CPU and GPU (D31), except the measured files
`flow.log`, `flow.md`, `flow.json` and `run.jsonl`, which hold timings.
`result.json` holds nothing measured and is identical too; it lists every
other output file with its SHA-256 hash. HDF5 files are written
without object timestamps. Dates (in reports and in `run.jsonl`) come from
`SOURCE_DATE_EPOCH` when it is set.

## Provenance (D85)

The `provenance` object in `result.json`, in the header of `run.jsonl` and
in `foehn info --json`: `version`, `commit` (`-dirty` with uncommitted
changes), `features`, `backends`, `backend`, `platform`, `cpu`, `gpu`,
`ranks`, `gpus`, `reproducible_mode` and `case_hash` (the SHA-256 of the
resolved case without `run.threads`, `run.backend`, `run.ptx` and
`output.directory`, which do not change the results). Two results with the
same `case_hash` and `commit` are comparable bit for bit.

## JSON objects of the command line

| `format` | From |
|---|---|
| `foehncraft.result` | `foehn flow --json` (also `result.json`) |
| `foehncraft.check` | `foehn flow --check --json`: `case_file`, `ok`, `diagnostics`, `estimate` (cells, faces, steps, host and device memory, output bytes and files, wall time and its basis) |
| `foehncraft.error` | any command with `--json` that fails before a run: `command`, `exit_code`, `diagnostics` |
| `foehncraft.info` | `foehn info --json` |
| `foehncraft.command` | the other commands with `--json` (`mesh`, `scalar`, `poisson`, ...): `command` and the report |

## Diagnostics (D78)

A diagnostic is an object:

```json
{
  "code": "CASE002",
  "severity": "error",
  "key": "mesh.generate.cels",
  "location": [4, 1],
  "message": "unknown key `cels` in [mesh.generate]; did you mean `cells`?",
  "fix": null
}
```

`key` is the case key (`section.key`), `boundary` the patch of a
boundary-condition finding, `location` line and column in the case file
(for syntax errors and unknown keys), and `fix` a list of `--set`
arguments that resolve it, where one exists (`{"set": ["time.adaptive=false"]}`).
A code keeps its meaning; a retired code stays in the table.

| Code | Severity | Title |
|---|---|---|
| `CASE001` | error | The case file is not valid TOML. |
| `CASE002` | error | A key that the case file does not have (with the nearest key as a hint). |
| `CASE003` | error | A value of the wrong type or an unknown choice. |
| `CASE004` | error | A value outside its range. |
| `CASE005` | error | Settings that do not fit together. |
| `CASE006` | error | A file that cannot be read: the case, the mesh or the restart file. |
| `CASE007` | error | A path outside the read roots or the write roots (D82). |
| `CASE008` | error | A `--set key=value` that cannot be applied. |
| `CASE009` | error | A key that the case needs is missing. |
| `BC001` | error | A mesh patch without a boundary condition. |
| `BC002` | error | A boundary condition for a patch the mesh does not have. |
| `BC003` | error | A condition that a later stage brings. |
| `BC004` | error | A velocity inlet with none or more than one of `velocity`, `normal_velocity` and `profile`. |
| `BC005` | error | `profile = "analytic"` without an `[analytic]` table. |
| `BC006` | error | `profile = "parabolic"` without `max_velocity` or `profile_axes`. |
| `BC007` | warning | `max_velocity` or `profile_axes` without a parabolic profile. |
| `BC008` | warning | A mass-flow inlet with a negative mass flow. |
| `BC009` | warning | Turbulence values in a laminar case. |
| `BC010` | warning | A wall velocity with a component normal to the wall. |
| `BC011` | note | The exact solution leaves the domain through a velocity inlet. |
| `BC012` | warning | A velocity inlet whose velocity points out of the domain. |
| `BC013` | warning | A symmetry patch that is not planar. |
| `BC014` | error | A net mass flow through the prescribed boundaries and no pressure boundary to balance it. |
| `BC015` | note | No pressure boundary: the pressure level is fixed in the reference cell. |
| `BC016` | error | A flow path whose flow rate is over-specified. |
| `BC017` | warning | Nothing drives a flow: the result is the fluid at rest. |
| `BC018` | note | No prescribed inflow: the flow rate follows from the pressure boundaries. |
| `RUN001` | note | Stopped by the `.stop` control file. |
| `RUN002` | note | Stopped by a signal (SIGTERM, SIGINT or SIGUSR2). |
| `RUN003` | warning | Stopped at the wall-time limit (`[run.limits] max_wall_time`). |
| `RUN004` | warning | Stopped at the step or iteration limit (`[run.limits] max_steps`, `max_iterations`). |
| `RUN005` | warning | Stopped at the output-volume limit (`[run.limits] max_output`). |
| `RUN006` | error | A steady solve that did not converge to the tolerance. |
| `RUN007` | error | A verification check of the case failed (order, error, reference or benchmark band, mass balance). |
| `RUN008` | warning | A change of the re-read case file that cannot be applied during a run. |
| `RUN009` | warning | A re-read case file that does not validate; the run goes on with the settings it had. |
| `RUN010` | warning | A restart file in another cell or face order: the continuation is not bit-identical. |
| `RUN011` | error | `--restart latest` found no restart file in the output directory. |
| `RUN012` | note | No benchmark of this hardware: the estimate has no wall time. |
| `RUN013` | error | An output file could not be written. |
| `GEN001` | error | An error without a more specific code. |

## Engine API (D86)

The command line is a thin layer on `foehncraft_flow::engine`; a Rust
program (and later the Python bindings, D61) drives a run the
same way:

```rust
use foehncraft_flow::engine::{self, RunOptions, Roots};

let case = engine::load(path, &["fluid.viscosity=1e-3".into()], &Roots::default())?;
let opts = RunOptions::default();               // interrupts, build info, roots
let stop = opts.interrupts.clone();              // stop.request_stop(StopSignal::Api) from another thread
let out = engine::run(&case, "cavity", dir, &opts, |line| eprintln!("{line}"))?;
std::process::exit(out.exit_code());             // out.summary is the result.json object
```

| Item | Purpose |
|---|---|
| `load(path, sets, roots)` | Reads a case with overrides, inside the roots. |
| `check_file(path, sets, roots, dir, build)` | The check report (never fails: every problem is a diagnostic). |
| `estimate(case, dir, build)` | The estimate (D81). |
| `run(case, name, dir, opts, log)` | Runs; writes `result.json` and `run.jsonl`; returns the outcome (report and summary). |
| `Interrupts` | `request_stop`, `request_restart`, and the flags a signal handler sets. |
| `latest_restart(dir)` | The newest restart file. |
| `read_result(dir)` | `result.json` of an earlier run. |
| `diag::of_error(&e)` | The diagnostic of an error. |

## C-free core (D87)

`foehncraft-flow` without default features (no HDF5) depends on no C
library: `tools/check_c_free.sh` checks the dependency tree in CI. HDF5
(output and restart) and METIS (partitioning) are optional.

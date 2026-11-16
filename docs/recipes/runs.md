# Run control and health

Part of the [Recipes](../recipes.md): restart, stop and steer runs, drive runs from a script or an agent, balances and divergence reports.

| # | Recipe |
|---|---|
| 22 | [Restart a run, stop it and steer it while it runs](#restart) |
| 26 | [Drive runs from a script or an agent](#run-contract) |
| 27 | [Check the balances and read a divergence report](#health) |

---

<a id="restart"></a>
## 22. Restart a run, stop it and steer it while it runs

Every run writes `restart.h5` at its end (`run.restart_final`), and on a
schedule if asked. A run continued from it gives bitwise the result of
the run without the interruption (same mesh, same build):

```toml
[run]
restart_every = 500          # restart.h5 every 500 steps
restart_at = [2000, 4000]    # also restart_002000.h5, restart_004000.h5
```

```
target/release/foehn flow my_case.toml
# continue to a later end time, into another directory:
target/release/foehn flow my_case.toml --restart output/my_case/restart.h5 --end-time 20 --out output/my_case_2
```

or, in the case file, `[run] restart = "output/my_case/restart.h5"` (a
path in the case file is relative to the case file, on the command line
to the current directory).

While a case runs, files in its output directory steer it (the solver
never deletes them; remove them before the next run):

```
touch output/my_case/.stop            # write restart.h5 and the results, then end
touch output/my_case/.result_1500     # write fields_001500.vtkhdf at step 1500
touch output/my_case/.read_toml_1200  # re-read my_case.toml at step 1200
```

A re-read applies the output and restart schedule, the stop criteria, the
relaxation, the linear-solver tolerances and the time-step settings; the
log lists each change (`key: old -> new`) and each change it cannot apply
during a run. `run.reread_every = N` re-reads every N steps without a
control file.

[↑ back to top](#run-control-and-health)

---

<a id="run-contract"></a>
## 26. Drive runs from a script or an agent

`foehn` keeps a contract for programs that run it (D62, D64, D78; the
details are in [interfaces.md](../interfaces.md)): one JSON object on
standard output with `--json`, an exit code that says how the run ended,
codes for every finding, and a result file without timings. From the
repository root, after `cargo build --release`:

```
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml > result.json; echo "exit $?"
python3 -c "import json; r = json.load(open('result.json')); print(r['status'], r['solves'][0]['iterations'])"
```

| Exit code | Status | Next step of a script |
|---|---|---|
| 0 | `finished` | use the result |
| 1 | `error` | read `diagnostics`: code, key, message and a fix where one exists |
| 2 | `failed-checks` | the run ended, a check failed (`RUN006` not converged, `RUN007` a check of the case) |
| 3 | `stopped` | a signal, `.stop` or a limit: continue with `--restart latest` |

**Change keys without editing the file**: `--set KEY=VALUE`, repeatable,
the value in TOML (a bare word is a string). The overrides are in the
result and survive a re-read of the case file.

```
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml --set fluid.viscosity=0.005 --set solver.relaxation_velocity=0.8 --out output/re200
```

**Check before running** (D81): the findings with codes and an estimate
of memory, output and, where this hardware has a benchmark, wall time:

```
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml --check
```

A finding with a fix gives the `--set` arguments that resolve it:

```
target/release/foehn --json flow examples/cases/flow/taylor_green_cfl10.toml --check --set time.max_growth=3
```

**Budgets** (D83): `[run.limits]` (`max_wall_time` in s,
`max_iterations`, `max_steps`, `max_output` in bytes) stops a run with its
restart file and results and exit code 3. Continue it with
`--restart latest`, bitwise as if it had not stopped:

```
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml --set run.limits.max_iterations=300 --out output/budget; echo "exit $?"
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml --restart latest --out output/budget; echo "exit $?"
```

**Signals** (unix): SIGUSR1 writes `restart.h5` and the run goes on;
SIGTERM, SIGINT and SIGUSR2 stop it cleanly (exit 3); a second Ctrl-C
ends at once. A batch job that is sent SIGTERM before its end is
continued by the next job with `--restart latest`:

```
target/release/foehn --quiet flow examples/cases/flow/cavity_re100_transient.toml --out output/job & pid=$!
sleep 5; kill -USR1 $pid      # a restart file now
sleep 5; kill -TERM $pid; wait $pid; echo "exit $?"
target/release/foehn --quiet flow examples/cases/flow/cavity_re100_transient.toml --restart latest --out output/job
```

**Follow a run**: `run.jsonl` gets one JSON object per iteration or step,
flushed at every log line (`tail -f output/job/run.jsonl`).

**Confine the files** (D82): `--read-root DIR` and `--write-root DIR`
refuse every path outside (`CASE007`), for agents that must not touch the
rest of the machine:

```
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml --read-root examples --write-root output
```

**Compare results**: two identical runs give byte-identical outputs
(D84); `result.json` lists their SHA-256 hashes, and
`provenance.case_hash` with `provenance.commit` says whether two results
are comparable. The keys and their units: `foehn schema` (JSON Schema)
and [keys.md](../keys.md).

[↑ back to top](#run-control-and-health)

---

<a id="health"></a>
## 27. Check the balances and read a divergence report

Every flow run books mass and kinetic energy in a ledger (D88, D89;
details in [interfaces.md](../interfaces.md#conservation-ledger-d88-d89-d90-part-4)).
From the repository root, after `cargo build --release`:

```
target/release/foehn flow examples/cases/flow/cavity_re100_energy.toml --out output/energy
```

The log shows the full table every 100 iterations (`[balance] every`) and
at the end; for the steady cavity the lid power goes into the viscous
dissipation and, about 1.5 %, into the numerical dissipation, which the
second law wants positive (D90 part 4); its two parts, the pressure part
of the momentum interpolation and the part of the convection scheme,
follow in the table. The same
table is in `flow.md`, the history in `balance.csv`, the bars in
`balance.svg`, the values under `solves[].balance` in `result.json`.

A transient run logs one compact balance line with every progress line
and a ledger of the whole run at the end; draw the time histories with

```
target/release/foehn flow examples/cases/flow/taylor_green_bdf2.toml --set "verification.time_steps=[]" --out output/tg
python3 -m tools.plot.balance_chart output/tg
```

**Tolerances**: `[balance] tolerance_step` (1e-8, the relative mass
imbalance of a step or of the final steady state), `tolerance_run` (1e-4,
the mass over the run) and `tolerance_energy` (1e-2, the negative
numerical dissipation relative to the largest term) give the warnings
`BAL001`–`BAL003`; `action = "stop"` ends a transient run at the
step with exit code 3 instead (here with a tolerance no run meets: the
Taylor–Green vortex has inflow and outflow through every side, and its
mass imbalance is at round-off, about 1e-17):

```
target/release/foehn flow examples/cases/flow/taylor_green_bdf2.toml --set "verification.time_steps=[]" --set balance.tolerance_step=1e-30 --set 'balance.action="stop"' --out output/stop; echo "exit $?"
```

**The energy budget over the mesh** (the gallery chart):

```
python3 -m tools.plot.energy_budget_chart output/energy_budget.svg
```

**A divergence** (D63) ends the run with exit code 1 and `RUN014`, after
it wrote the last valid state (`fields_last_valid.vtkhdf`, `restart.h5`),
the failing cells with three layers of neighbours (`divergence.vtkhdf`,
threshold on `divergence_layer` in ParaView) and the report
`divergence.md` / `divergence.json` with the likely causes ranked and what
to change. Provoke one: SIMPLE with no relaxation

```
target/release/foehn flow examples/cases/flow/cavity_re1000.toml --set 'solver.algorithm="simple"' --set solver.relaxation_velocity=0.99 --set solver.relaxation_pressure=1.0 --out output/diverge; echo "exit $?"
cat output/diverge/divergence.md
```

The first cause is `DIV002` (relaxation) with the fix
`--set solver.relaxation_velocity=0.69 --set solver.relaxation_pressure=0.31`,
which runs. A time step far too large for a transient case points to
`DIV001` with a smaller time step:

```
target/release/foehn flow examples/cases/flow/cavity_re100_transient.toml --set time.adaptive=false --set time.time_step=5.0 --set solver.correctors=1 --set 'time.scheme="crank-nicolson"' --out output/diverge_t; echo "exit $?"
```

With `--json`, the `RUN014` diagnostic in `result.json` carries the fix
of the most likely cause as `--set` arguments, so a script can retry.

[↑ back to top](#run-control-and-health)

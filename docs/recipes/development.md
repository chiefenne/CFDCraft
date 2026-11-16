# Development: tests, CI and releases

Part of the [Recipes](../recipes.md): MPI tests, checking a change the way CI does, releases, derivative checks of the kernels.

| # | Recipe |
|---|---|
| 4 | [Run the MPI test with several ranks on one machine](#mpi) |
| 5 | [Check a change the way CI does](#ci-local) |
| 6 | [Make a release (tag and GitHub release)](#release) |
| 25 | [Check the derivatives of the kernels (adjoint readiness)](#adjoint-ready) |

---

<a id="mpi"></a>
## 4. Run the MPI test with several ranks on one machine

Needs an MPI library and libclang (Ubuntu:
`sudo apt install libopenmpi-dev openmpi-bin libclang-dev`; macOS:
`brew install open-mpi`).

```
cargo test -p foehncraft-parallel --features mpi
cargo build --release -p foehncraft-parallel --features mpi --example mpi_smoke
mpirun --oversubscribe -n 4 target/release/examples/mpi_smoke
```

Prints `MPI smoke test passed` after a ring exchange (the step of a halo
exchange) and a reproducible sum that gives the same bits on every rank.

[↑ back to top](#development-tests-ci-and-releases)

---

<a id="ci-local"></a>
## 5. Check a change the way CI does

```
cargo fmt --all --check
tools/ci_cases.sh --check                               # every flow case is in a CI set
cargo clippy --workspace --all-targets -- -D warnings
cargo clippy --workspace --all-targets --features foehncraft-cuda/cuda,foehncraft-cli/cuda -- -D warnings
python3 tools/check_backend.py                          # D60: CUDA only in crates/foehncraft-cuda and gpu/
cargo test --workspace --profile ci                     # what a pull request runs (--release on main)
cargo test -p foehncraft-mesh --no-default-features      # without HDF5 and METIS
cargo test -p foehncraft-fv --no-default-features
gpu/ptx-kernels/build.sh && tools/check_kernels.sh     # PTX, and R10: no atomics, no f32, no FMA
cargo build --profile ci -p foehncraft-cli
for s in linux-1 linux-2 linux-3; do                    # add "long" for the cylinder at Re 100
  tools/ci_cases.sh "$s" target/ci/foehn || break
done
pip install -r tools/requirements-ci.txt
python3 tools/check_vtkhdf.py output/*/mesh.vtkhdf output/*/fields.vtkhdf output/*/*/fields.vtkhdf
```

The workflow is `.github/workflows/ci.yml`. It has two tiers:

- **Pull requests** (fast gate, about 5 minutes): the `ci` profile (release
  optimisation without LTO, results bitwise equal to `release`), all tests
  including the bitwise-reproducibility tests on Linux, macOS and Windows,
  the case catalogue in three parallel Linux jobs, the quick cases on macOS
  and Windows, hypre, MPI and PTX. Code jobs are skipped when only `docs/`,
  Markdown files, `assets/` or `LICENSE` change.
- **Full tier**: every push to main, every night at 02:07 Vienna time (01:07
  UTC), a pull request with the label `full-ci`, or by hand (Actions, CI,
  *Run workflow*): the `release` profile, debug tests, the cylinder at
  Re 100, and in nightly and manual runs also the tests marked slow
  (`--ignored`). Main also builds the dependencies in the `ci` profile
  once, so that pull requests start from a warm cache.

The check "CI OK" sums up all jobs.

[↑ back to top](#development-tests-ci-and-releases)

---

<a id="release"></a>
## 6. Make a release (tag and GitHub release)

A release is a Git tag `v<stage>` (D46) plus a GitHub release. Steps:

1. A pull request moves the `## Unreleased` section of `CHANGELOG.md` to
   `## <version> (YYYY-MM-DD) – ...` and adds a row to the summary table.
2. After the merge: GitHub, Actions, **Release**, *Run workflow*, version
   `0.1`, first with *dry run* checked (prints the notes), then without.

The workflow refuses a version without a changelog section, and creates
only the release if the tag already exists. The assistant starts it through
the GitHub API in the same way.

The release badge of the README updates by itself: the last step of the
workflow waits until img.shields.io reports the new tag and then purges
GitHub's cached copies of the README images (GitHub serves them through its
image proxy, which otherwise keeps the old badge for hours or days). If the
badge still shows the old version, run that purge by hand:

```
curl -s https://github.com/chiefenne/Foehncraft | grep -oE 'https://camo\.githubusercontent\.com/[0-9a-f]+/[0-9a-f]+' | sort -u | xargs -n1 curl -s -X PURGE
```

---

<a id="adjoint-ready"></a>
## 25. Check the derivatives of the kernels (adjoint readiness)

The kernels are prepared for a discrete adjoint solver (D59): they run on
dual numbers, matrices apply transposed, and transient runs reverse with
binomial checkpointing. The checks are tests; from the repository root:

```
cargo test --release -p foehncraft-flow --test adjoint -- --nocapture
cargo test --release -p foehncraft-flow --lib -- adjoint_tests checkpoint
cargo test --release -p foehncraft-mesh shape_derivatives -- --nocapture
```

The first prints one line per kernel, for example
`momentum_face: 101700 derivatives, largest relative difference 5.28e-8`.

**Derivatives of a kernel in your own code.** Instantiate the body on
`foehncraft_core::Dual<N>` instead of `f64`: give the inputs you
differentiate by the derivative 1 in one of the N directions
(`Dual::variable(x, k)`), all others as constants (`Dual::constant(x)`),
and read the N derivatives of each result from its `.d` (the value is
`.v`, the same bits as the `f64` result). Leave the limiter array ψ
constant for the frozen-limiter linearisation; switches such as the upwind
direction follow the values.

**Transposed matrices** on a flow device: `transpose(val, out)` writes Aᵀ
in the storage of A (for the solvers and the multigrid),
`transposed_product(val, x, out)` computes Aᵀx without storing it.

**Reversing a transient run**: implement `foehncraft_flow::checkpoint::Primal`
(or use `FlowRun` for the flow solver with a fixed time step) and call
`checkpoint::reverse(run, steps, snapshots, adjoint)`; the adjoint step i
is called for i = n − 1, …, 0 with level i current. `schedule(n, s)` lists
the actions, `python3 tools/checkpoint_chart.py out.svg` draws the cost
against the snapshots.

[↑ back to top](#development-tests-ci-and-releases)

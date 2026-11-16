# Scalar transport and the pressure solver

Part of the [Recipes](../recipes.md): order of accuracy of the scalar operators, convection schemes, pressure-solver studies, hypre.

| # | Recipe |
|---|---|
| 12 | [Check the order of accuracy of the scalar operators](#scalar-mms) |
| 13 | [Compare convection schemes and limiters](#scalar-schemes) |
| 15 | [Study the pressure solver on a Poisson problem](#poisson-study) |
| 17 | [Build hypre and compare with it](#hypre) |

---

<a id="scalar-mms"></a>
## 12. Check the order of accuracy of the scalar operators

The manufactured-solution sweeps of stage 2 ([explanation](../verification/stage2.md)):
each case solves on three or four nested refinements, prints the observed
order and fails if it is outside 1.8–2.2.

```
cargo build --release -p foehncraft-cli
for m in hex tet poly mixed; do
  target/release/foehn scalar examples/cases/scalar/mms_${m}_box.toml
done
target/release/foehn scalar examples/cases/scalar/mms_poly_2d.toml
cat output/mms_poly_box/scalar.md          # errors, orders, time per phase
```

Keys: `verification.levels` (refinement levels of `[mesh.generate]`,
override with `--levels 1,2,3,4`), `verification.expected_order`,
`problem.manufactured`, `problem.boundary.<patch>` with `profile =
"manufactured"`, `mesh.generate.perturb`. The same checks as tests (one
level finer with `--ignored`):

```
cargo test --release -p foehncraft-fv
cargo test --release -p foehncraft-fv -- --ignored --nocapture   # about 1.6 M tetrahedra, minutes
cargo run --release -p foehncraft-fv --example sweep_table         # the table of the explanation page
```

[↑ back to top](#scalar-transport-and-the-pressure-solver)

---

<a id="scalar-schemes"></a>
## 13. Compare convection schemes and limiters

One case, several solver settings: `[[variants]]` runs the case once per
variant, writes the fields of each into its own directory and the profiles
along `output.profile` into one plot.

```
target/release/foehn scalar examples/cases/scalar/oblique_step_2d.toml
cat output/oblique_step_2d/scalar.md       # φ range per scheme
# profile plot: output/oblique_step_2d/profile.svg
# fields: output/oblique_step_2d/<variant>/fields.vtkhdf (open in ParaView, field "phi" and "limiter")
```

Keys: `solver.scheme` (`upwind`, `linear-upwind`, `central`),
`solver.limiter` (`none`, `venkatakrishnan`, `venkatakrishnan-wang`),
`solver.limiter_k` (0: strictly bounded; larger: smoother, small
overshoots), `[[variants]].bounds`, `output.profile`. The bounded central
scheme (stage 5, D25) takes `solver.central_bounding` β (0.1; the variant
uses 0.5 for this pure convection) and `solver.central_blending` (1; the
share of the central correction, the rest linear upwind); on these
distorted polygons it needs the limiter as well, which scales its
correction in each upwind cell. To try another K, add a variant:

```toml
[[variants]]
name = "venkatakrishnan-k03"
solver = { scheme = "linear-upwind", limiter = "venkatakrishnan", limiter_k = 0.3, max_outer = 400 }
```

[↑ back to top](#scalar-transport-and-the-pressure-solver)

---

<a id="poisson-study"></a>
## 15. Study the pressure solver on a Poisson problem

The mesh sweeps of stage 3 ([explanation](../verification/stage3.md)): each
case solves the Poisson problem on three or four refinement levels with CG
and the multigrid (mixed precision and float64) and with Jacobi-CG, on 4
threads and on 1, prints the iteration ratio between the finest and the
coarsest level and fails if it exceeds 1.5.

```
cargo build --release -p foehncraft-cli
for m in hex tet poly mixed; do
  target/release/foehn poisson examples/cases/poisson/poisson_${m}_box.toml
done
cat output/poisson_poly_box/poisson.md      # iterations, times, hierarchy, differences
# plot: output/poisson_poly_box/iterations.svg
```

Try other multigrid settings from the command line or in `[linear.amg]`:

```
target/release/foehn poisson examples/cases/poisson/poisson_poly_box.toml --solvers amg,amg-plain,jacobi
target/release/foehn poisson examples/cases/poisson/poisson_poly_box.toml --levels 4 --threads 1 --repeat 3
```

```toml
[linear.amg]
smoother = "jacobi"        # or "chebyshev"
degree = 3                 # smoother steps before and after the coarse correction
strength = 0.02            # more aggressive coarsening
precision = "double"       # float64 hierarchy
```

The same solver in `foehn scalar` (gallery 11):

```toml
[solver]
linear_solver = "cg"       # symmetric matrices only (no convection)
preconditioner = "amg"     # also with linear_solver = "bicgstab"
amg = { smoother = "chebyshev", degree = 2, precision = "mixed" }
```

The checks as tests (the mesh sweep in release builds):

```
cargo test --release -p foehncraft-fv --test stage3
```

[↑ back to top](#scalar-transport-and-the-pressure-solver)

---

<a id="hypre"></a>
## 17. Build hypre and compare with it

hypre (MIT or Apache-2.0) is the optional reference solver (D20). It is
built from source, sequential (no MPI), as a shared library; Foehncraft
links it only with the cargo feature `hypre` and finds it through
`HYPRE_DIR`. Needs git, CMake 3.21 or newer and a C compiler (on macOS:
`brew install cmake`).

```
python3 -m tools.bench.build_hypre $PWD/target/hypre                 # about 1–3 minutes
export HYPRE_DIR=$PWD/target/hypre
cargo build --release -p foehncraft-cli --features hypre
target/release/foehn info                            # features: ..., hypre
for m in hex tet poly mixed; do
  target/release/foehn poisson examples/cases/poisson/poisson_${m}_box.toml \
      --solvers amg,amg-double,amg-plain,jacobi,hypre --backends cpu,cpu-serial \
      --out output/poisson_${m}_hypre
done
cargo test --release -p foehncraft-fv --features hypre --test stage3
```

With hypre in `study.solvers` the differences in the report are measured
against hypre's solution, and the run fails above `study.max_difference`
(1e-6). hypre runs on one CPU thread; compare its times with the
`cpu-serial` rows. BoomerAMG settings: `[hypre]` in the case file
([keys](../cases/poisson.md#poisson-case-file-keys)).

[↑ back to top](#scalar-transport-and-the-pressure-solver)

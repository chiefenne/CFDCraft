- **Multigrid for the momentum sweep of the coupled algorithm (D179, the
  momentum part of L44).** The preconditioner of the coupled algorithm
  now applies one V-cycle of a reduction multigrid to each velocity
  component instead of BiCGStab: classical strength of connection, the
  first pass of the Ruge–Stüben splitting into coarse and fine points, a
  restriction that approximates the ideal one by a two-term Neumann
  series (Manteuffel, Ruge & Southworth 2018), direct interpolation, and
  Jacobi relaxation on the fine points only (`amg::reduction`). The
  cycle runs on the level kernels of the pressure multigrid (CPU and GPU
  from one source, float32 hierarchy); the setup is deterministic. On
  the backward-facing step the Krylov iterations per outer iteration
  halve (10 → 5.3 on 8 200 cells, 14 → 6.0 on 32 800), the fine mesh
  settles in 150 instead of 175 iterations, and the run time to the
  settled result on one thread is 7.0 s (SIMPLEC 18.7 s) and 74 s
  (SIMPLEC 243 s, the BiCGStab sweep 131 s); the cavity at Re = 1000 on
  128² converges in 12.1 s (BiCGStab sweep 18.8 s, SIMPLEC 28.3 s).
  Converged results do not change. New keys `[solver.coupled]
  momentum_sweep` (`"multigrid"`, the default, or `"bicgstab"`),
  `momentum_cycles`, `momentum_reuse` and the table
  `[solver.coupled.multigrid]`. Gallery 80, recipe 80, requirements
  step 22 sections 3.1 and 7.7, verification stage 6 part 6f.
- **Kept hierarchies.** A hierarchy of the momentum sweeps is kept while
  the momentum matrix has changed by at most 2 % since its setup
  (`momentum_reuse`) and while its cycle still halves the residual;
  components with the same matrix share one, a component without a
  residual (the third one of a planar case) needs none. The restart file
  holds the matrices of the kept hierarchies (`momentum_amg_use`,
  `momentum_amg_matrix_<n>`), so a continued run stays bit-identical;
  older restart files still load. If a setup fails, BiCGStab does the
  sweeps of that outer iteration and the log carries a warning. The
  report of a run counts the setups of the momentum sweeps.
- **One thread without a pool.** On a pool of one thread (meshes below
  20 000 cells, or `--threads 1`) every kernel call and every level
  operation of a multigrid cycle was handed to the worker thread and
  waited for, about 2 000 times per outer iteration of a cavity run.
  Such loops, and loops of a single task on any number of threads (the
  coarse levels of a multigrid cycle), now run on the calling thread.
  The values are the same bit by bit. Small and medium cases gain with
  both algorithms: the step on 8 200 cells takes 18.7 s instead of 31 s
  with SIMPLEC (525 iterations) and 8.3 s instead of 21 s with the
  coupled algorithm and the BiCGStab sweep (65 iterations); the step on
  32 800 cells with SIMPLEC 243 s instead of 327 s.
- **Chart tool.** `tools.plot.krylov_chart` draws the Krylov iterations
  of the coupled solve against the outer iteration for several runs.

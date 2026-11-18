- **Coupled pressure–velocity algorithm (D176, L11).** `[solver] algorithm
  = "coupled"` solves the momentum equations and the continuity equation
  of a steady outer iteration as one linear system for the increments of
  velocity and pressure: matrix-free with the kernels of the segregated
  solver (no block matrix), flexible GMRES with a SIMPLE-type
  preconditioner, `[solver.coupled]`. The discrete equations and the
  converged result are those of SIMPLEC, and the face fluxes satisfy
  continuity after every outer iteration; the kinetic-energy ledger of
  the cavity at Re = 1000 is the same with both algorithms. CPU and GPU
  from one source (two new kernels).
- **Pseudo-time steps for steady runs (D177, D178).** With the coupled
  algorithm the momentum equations take one pseudo-time step for all
  cells (automatic: the volume of the domain over its boundary area, over
  a velocity scale; `[solver.pseudo_time]`) and the relaxation (default
  0.98) only caps the local Courant number; the turbulence equations take
  the step min(Δτ, 0.5 T) with the time scale T of the model. The
  backward-facing step settles in 65 (k-ε-ζ-f₀), 75 (k-ε) and 80 (SST)
  iterations instead of 525 to 675, on 32 800 cells in 175 instead
  of 1 600; the channel with k-ε-ζ-f₀ converges in 323 instead of about
  2 100, the cavity at Re = 1000 in 452 instead of 1 397 (counts good to
  about 10 %: they react to the last digit of the step). In run time the
  step gains a factor of 1.4 to 1.5 so far, because a coupled iteration
  costs five to seven SIMPLEC iterations; small laminar cases do not gain
  yet (L44, a multigrid for the momentum sweep, is the next step). The
  turbulence step is available with SIMPLEC too
  (`solver.pseudo_time.turbulence = 0.5`, off by default: 230 iterations
  on the step). New cases
  `backstep_zeta_f0_coupled.toml` and `cavity_re1000_coupled.toml`,
  gallery 78, recipe 78, requirements step 22.
- **Log.** The column `K` (`linear.coupled`) shows the Krylov iterations
  of the coupled solve; `residuals.csv` and the run record carry them as
  `krylov` for coupled runs, and the summary of a run shows the
  pseudo-time step of its last iteration (`pseudo_step` in the run
  record). The keys of `[solver.coupled]` and
  `[solver.pseudo_time]` may change during a run.
- **Chart tool.** `tools.plot.settling_chart` draws a monitored value
  against the outer iteration for several runs, with the ±1 % band and
  the iteration from which each run stays inside it.

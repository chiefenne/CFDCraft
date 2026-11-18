- **Coupled algorithm: looser pressure solves in the preconditioner, and
  the Schur complement measured.** `[solver.coupled] pressure_tolerance`
  defaults to 0.3 instead of 0.1: the Krylov iteration corrects what the
  pressure solve of the preconditioner leaves, so the Krylov count stays
  and the CG iterations fall by a fifth to more than a quarter. One thread, to the settled
  result: the coarse step takes 7.4 s instead of 8.0 s, k-ε 5.9 s instead
  of 6.4 s; the fine step is unchanged (57.1 s against 56.9 s), since
  its pressure solves stop at their iteration limit either way; the
  results do not change. Two published alternatives to the SIMPLE-type Schur complement
  were measured and do not pay (SIMPLER: a second pressure solve for a
  few Krylov iterations; the least-squares commutator: fine on the
  cavity, at the Krylov limit on the step). Requirements step 22
  section 7.9; the next step is a pressure hierarchy refreshed in every
  outer iteration.

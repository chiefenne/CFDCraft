- **Coupled algorithm: the pressure hierarchy of the preconditioner
  follows its matrix (D181, stage 6 part 6h).** The multigrid hierarchy
  of the pressure solve inside the coupled preconditioner is built anew
  when its matrix has changed by more than
  `[solver.coupled] pressure_reuse` (5 %) since its setup, besides every
  `amg_rebuild` iterations as before (`pressure_reuse = 1` restores the
  periodic rebuild alone). A hierarchy up to nine iterations old needed
  9 to 20 CG iterations per solve where a fresh one needs two. One
  thread, to the settled result: the fine backward-facing step (32 800
  cells) takes 39.8 s instead of 57.1 s, 6.1 times faster than SIMPLEC;
  the coarse step 6.2 s instead of 7.4 s, k-ε 4.7 s instead of 5.9 s.
  The settle counts and the converged results do not change; laminar
  cases do not trigger the rule. Requirements step 22 section 7.10,
  gallery and recipe 83.

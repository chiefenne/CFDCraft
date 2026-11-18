- **Multigrid for the turbulence equations and larger steps with it
  (D180, the turbulence part of L44).** With the coupled algorithm the
  transport equations of the k-ε family (k-ε, realizable k-ε, k-ε-ζ-f,
  k-ε-ζ-f₀) take one cycle of the reduction multigrid of D179 as the
  preconditioner of BiCGStab, with kept hierarchies (`[solver.turbulence]
  preconditioner = "reduction"`, the default there, and `reuse`), and
  the relaxations of the momentum and turbulence equations default to
  0.995 instead of 0.98. The relaxation acts as a local step that
  shrinks with the cells: at 0.98 it made the outer iterations grow with
  the mesh. The backward-facing step settles in 40 iterations on 8 200
  cells and 65 on 32 800 (65 and 150 before), k-ε in 45 (75), the
  channel converges in 231 (323). One thread: the fine step takes 56 s
  instead of 74 s; the coarse steps about the same time as before (7.5 s
  against 7.0 s, k-ε 5.7 s against 4.5 s), since an iteration costs
  more. SST keeps Jacobi and 0.98 (its steady state on the step is not
  held with accurate solves of its equations), laminar cases keep 0.98.
  Converged results do not change. The restart file holds the matrices
  of the kept hierarchies; a failed setup falls back to Jacobi with a
  warning; the report counts the setups. Gallery 81, recipe 81,
  requirements step 22 sections 4.2 and 7.8, verification stage 6 part
  6g.

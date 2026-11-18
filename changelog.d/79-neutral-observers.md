- **Intermediate results and balance tables no longer change a steady
  run (D31).** A `.result_N` request and the balance table of an
  iteration (`[balance] every`) updated the boundary values and gradients
  for their output and left them in place, so the next iteration started
  from values that had seen one more update: the path of a steady run
  depended on its output schedule (a duct with 24 × 8 cells: 7 % in the
  momentum residual of the iteration after a result; the step case 10⁻⁶
  with SIMPLEC and 10⁻⁴ with the coupled algorithm), though not its
  converged result. Both now put the boundary values, gradients and
  limiters of the iteration back. A run with result files or balance
  tables repeats the plain run bit by bit, with SIMPLEC and with the
  coupled algorithm (test
  `intermediate_results_and_balances_do_not_change_the_run`). Transient
  runs were not affected. The figure of gallery 78 is redrawn from
  histories recorded this way: the SST step settles in 80 coupled
  iterations (85 before). The step on 32 800 cells, run to the end
  with both algorithms (x_r/H = 6.352), settles in 175 coupled iterations
  against 1 600 with SIMPLEC, in 235 s against 327 s.

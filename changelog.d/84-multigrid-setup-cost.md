- **Multigrid setup: flat row builders (L44, no change of results).** The
  reduction multigrid of the coupled algorithm (momentum and turbulence
  equations) builds its strong graph, interpolation, restriction and
  Galerkin products into flat arrays, blocks of rows in parallel joined in
  row order, with one dense accumulator per worker instead of a vector
  per row. The matrices are bit-identical on any number of threads, and
  runs give the same residual history bit for bit. One thread: the
  multigrid setups of the fine backward-facing step take 4.9 s instead of
  6.9 s, the run 33.3 s instead of 39.8 s (7.3 times faster than
  SIMPLEC); the coarse step 5.7 s instead of 6.2 s. Requirements step 22
  section 7.11.

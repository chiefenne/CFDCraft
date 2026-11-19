- **Faster rebuild of the sliding mesh by imprinting (D182 variant B).**
  The solver and its pressure multigrid are still built for every step (the
  hierarchy of the step before was measured and does not fit: the rotor
  turns by about one cell per step, and the conjugate gradients needed 5 to
  11 times the iterations), but each build costs less: the dense inverse of
  the coarsest multigrid level takes a third of the work over contiguous
  rows (5 ms to 0.9 ms for 175 unknowns, in every setup of every run),
  multigrid setups below 20 000 rows run on one helper thread instead of
  handing each short parallel loop to the pool, the mesh geometry and the
  discretisation of a step are computed once instead of twice, and the
  state is moved instead of copied. The case file `taylor_couette_imprint.toml`
  runs in 5.6 s instead of 8.4 s, 2 560 cells in 45 s instead of 51 s, with
  the same results (L2 errors unchanged, histories bit-identical or to
  round-off). Verification page stage 6, finding 3 of the sliding mesh by
  imprinting.

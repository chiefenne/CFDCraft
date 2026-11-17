- **Periodic boundaries (stage 6 part 6p).** Two boundaries of
  `type = "periodic"` that name each other as `partner` form a
  translational pair (D29): the faces are matched by their centroids, each
  side keeps its own faces, and every face loop, gradient and limiter sees
  the cell behind the partner face shifted by the period. The flow through
  the pair is driven by its `mass_flow` (a uniform body force adjusted
  after every pressure correction) or by a fixed `pressure_gradient`; the
  driving force is in the log, the report and the restart file. Checks
  BC030 (pairs that do not match or name each other) and BC031 (a pair that
  nothing drives). Generated boxes keep their perturbed periodic sides
  translated copies. Results of cases without periodic boundaries are
  unchanged bit for bit.
- **Heat-flux jump of the passive scalars (D109 (1)).**
  `periodic = "heat-flux"` in a scalar table of a periodic boundary carries
  the scalar across the pair with the jump Δ = Q/ṁ of the wall fluxes and
  sources, `bulk` fixes its level; refused with fixed wall or inlet values
  (D109 (4)). The periodic heated channel gives Nu = 8.2460 against
  140/17 = 8.2353 and the short periodic heated pipe, one diameter long,
  Nu = 4.3691 against 48/11 = 4.3636 (D111), both at second order; the
  periodic Poiseuille flow equals the open channel (gallery 67, recipe 67).

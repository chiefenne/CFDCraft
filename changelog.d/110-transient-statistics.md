- **Statistics of transient runs (D71, gallery and recipe 110).**
  `[output.statistics]` takes time means, RMS and covariances of the fields
  (`statistics.vtkhdf`: velocity, pressure, k, scalars, the resolved
  Reynolds shear stresses and scalar fluxes) and of the probes and force
  coefficients (`statistics.csv`), weighted with the time step, from
  `start` to the end. Each series gets the 95 % confidence interval of its
  mean from batch means and a verdict from the lag-1 correlation of the
  batch means; a window that is too short warns with STAT001.
  `phase_period` writes phase averages (`statistics_phase.csv`). Checked on
  the decaying Taylor–Green vortex against its exact time moments
  (`tests/statistics.rs`).

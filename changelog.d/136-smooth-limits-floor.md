- **Fix: smooth limits keep floors (spec 2.6).** With
  `turbulence.smooth_limits = true` the SST flat plate did not converge:
  the smooth width followed the larger operand, so the floor of CD_kω
  became a large, non-monotone value. The width of `smin`/`smax` now
  follows the smaller operand, m = |a||b|/(|a| + |b|), in a form free of
  cancellation. Runs without smooth limits are unchanged. New test: the
  flat-plate C_f with smooth limits is within 0.090 % of the hard form
  (band 0.1 %).

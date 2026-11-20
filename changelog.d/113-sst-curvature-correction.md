- **Curvature and rotation correction of SST (6h, D170, gallery 113,
  recipe 113).** `turbulence.curvature_correction = true` multiplies the
  production of k and ω by the rotation function f_r1 of Spalart & Shur
  (1997) in the form of Smirnov & Menter (2009), from r* = S/Ω and the
  material derivative of the strain rate (constants `c_r1`, `c_r2`, `c_r3`).
  It sees the rotation of `[[mrf]]` zones; not with `[[sliding]]` yet.
  In pure shear the result is unchanged to the last digit.
- **A whole domain in a turning frame.** `flow.frame_rotation = [Ωx, Ωy, Ωz]`
  solves for the velocity relative to a frame that turns at Ω, with the
  Coriolis force (the centrifugal force stays in the pressure), for flows
  that are periodic in the rotating frame. The new case
  `channel_rotating_sst.toml` is the channel at Re_b 5 800 in spanwise
  rotation (Ro 0.15). There the wall shear stress of the unstable side is
  1.63 times that of the stable side with the correction (1.00 without).
- The derived field `curvature-correction` writes f_r1. TURB010 notes the
  correction together with Kato–Launder production.

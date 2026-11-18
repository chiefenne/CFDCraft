- **Passive scalars of a turbulent flow: turbulent diffusivity and thermal
  wall function (D135, stage 6 part 6f).** A scalar now diffuses with
  ν/Pr + ν_t/Pr_t on every face, and a wall face carries the conductance
  h_w = ρ u*/T⁺ of the wall law `scalars.<name>.wall_function` (`kader`,
  the default, valid for all y⁺ and any Pr, `jayatilleke`, or `none`).
  `turbulent_prandtl` (0.85) and `turbulent_prandtl_model` (`constant`
  or `kays-crawford`) set Pr_t. Fixed-flux walls report the wall value
  from T⁺. The wall file gains `<name>_wall`, `<name>_flux`, `<name>_h`
  and `<name>_t_plus`. The heated pipe at Re = 10⁵, Pr = 0.71 gives Nu
  182.6 (Kays–Crawford 181.3) against Gnielinski's 180.2; the flat
  plate St is +1.7 % to +7.5 % from the correlation for 0.25 ≤ x ≤ 1.
  Cases without a turbulence model do not change. Requirements step 6
  section 5, gallery and recipe 85.

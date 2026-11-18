# Experimental data of the impinging jet (specification 9.6)

`examples/studies/impinging_jet.py` compares its runs with the files of
this folder when they are here; without them it reports the models only.
Each file is a CSV with `#` comment lines (the source, the figure or
table, and how the numbers were taken) and a header:

| File | Columns | Data |
|---|---|---|
| `nusselt.csv` | `r_D`, `Nu` | Nu(r/D) on the plate, round pipe jet at H/D = 2, Re = 23 750, uniform wall heat flux: Baughn & Shimizu (1989), Baughn, Hechanova & Yan (1991) |
| `velocity_r0.5.csv`, `velocity_r1.csv`, `velocity_r2.csv`, `velocity_r2.5.csv`, `velocity_r3.csv` | `x_D`, `u_U` | the mean radial velocity over the distance from the plate at r/D = 0.5, 1, 2, 2.5 and 3, H/D = 2, Re = 23 000: Cooper, Jackson, Launder & Liao (1993) |

The first row of `nusselt.csv` is the stagnation point (r/D = 0); its Nu
is the reference of the Nu₀ bands. Only data whose licence allows a copy
in this repository go here (the licence of the source is named in the
file's comments).

References:

- J. W. Baughn and S. Shimizu, Heat transfer measurements from a surface
  with uniform heat flux and an impinging jet, J. Heat Transfer 111
  (1989) 1096–1098.
- J. W. Baughn, A. E. Hechanova and X. Yan, An experimental study of
  entrainment effects on the heat transfer from a flat surface to a
  heated circular impinging jet, J. Heat Transfer 113 (1991) 1023–1025.
- D. Cooper, D. C. Jackson, B. E. Launder and G. X. Liao, Impinging jet
  studies for turbulence model assessment — I. Flow-field experiments,
  Int. J. Heat Mass Transfer 36 (1993) 2675–2684.

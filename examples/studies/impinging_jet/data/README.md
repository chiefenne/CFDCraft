# Experimental data of the impinging jet (specification 9.6)

`examples/studies/impinging_jet.py` compares its runs with the files of
this folder when they are here; without them (or with `--no-experiment`)
it reports the models only. The CSV files are **not in the repository**:
the wiki of their source is licensed CC BY-NC-SA 4.0 (non-commercial,
share-alike), which does not fit the MIT licence of this repository, so
`.gitignore` keeps them out.

**Where they come from.** The ERCOFTAC classic database, case 025
(impinging jet; http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case025;
the terms of that site apply). Download the data files of the case, put
the `ij2lr-*.dat` files (H/D = 2, Re = 23 000) in a folder DIR and run

```
python3 examples/studies/impinging_jet/ercoftac_to_csv.py DIR
```

which writes this folder's files, each a CSV with `#` comment lines (the
source file and how the numbers were made) and a header:

| File | Columns | Data |
|---|---|---|
| `nusselt.csv` | `r_D`, `Nu` | Nu(r/D) on the plate, round pipe jet at H/D = 2, uniform wall heat flux: `ij2lr-nuss.dat`, Nu/Re^0.7 of Baughn et al. at Re = 23 000 (labelled 1992 in the file), multiplied by 23 750^0.7 for the Re of the study (`--re` changes it) |
| `velocity_r0.5.csv`, `velocity_r1.csv`, `velocity_r2.csv`, `velocity_r2.5.csv`, `velocity_r3.csv` | `x_D`, `u_U` | the single-wire mean velocity U/U_bulk over the distance from the plate at r/D = 0.5, 1, 2, 2.5 and 3, H/D = 2, Re = 23 000: `ij2lr-<nn>-sw-mu.dat` of Cooper et al. (labelled 1992 in the file) |

The first row of `nusselt.csv` is the stagnation point (r/D = 0); its Nu
is the reference of the Nu₀ bands. The same database has H/D = 6 and
Re = 70 000 (`ij2hr`, `ij6lr`, `ij6hr`) and the velocity variances and
shear stress (`-uu`, `-vv`, `-uv` files); the study does not use them yet.
The charts of the gallery are made with `--no-experiment`, so that they
do not carry the points.

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

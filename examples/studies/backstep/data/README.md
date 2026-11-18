# Backward-facing step of Driver & Seegmiller (1985): experimental data

The three files are the experimental data of the 2D backward-facing step
validation case of NASA's Turbulence Modeling Resource (TMR), unchanged:

| File | Content |
|---|---|
| `cf.exp.dat` | skin friction C_f along the lower wall with its uncertainty, x/H from −3.96 to 36 |
| `cp.expnew.dat` | pressure coefficient along the lower and the upper wall; `cp` shifted to 0 near x/H = 40 (the column used here), `cp_orig` referenced to x/H = −4 |
| `profiles.exp.dat` | mean velocity u/U_ref, v/U_ref and Reynolds stresses (1000 u'u'/U_ref² …, triple moments) at x/H = −4, 1, 4, 6 and 10 |

U_ref is the speed in the middle of the channel near x/H = −4. The
reattachment length of the experiment is x/H = 6.26 ± 0.10.

**Source.** NASA Turbulence Modeling Resource, case "2D Backward Facing
Step" (2DBFS), https://tmbwg.github.io/turbmodels/backstep_val.html; the
files are in the folder `Backstep_validation` of the repository
https://github.com/TMBWG/turbmodels (commit
20a39f549dbff988a6accc421ac84dafc4487695, retrieved 2026-11-18;
SHA-256 `cd9b6434…` cf.exp.dat, `2b3ec9f8…` cp.expnew.dat, `23e8b44d…`
profiles.exp.dat).

**Reference.** D. M. Driver and H. L. Seegmiller, "Features of
Reattaching Turbulent Shear Layer in Divergent Channel Flow," AIAA
Journal 23 (2), 1985, pp. 163–171, https://doi.org/10.2514/3.8890. The
C_p shift follows Eça et al., AIAA-2005-4728, AIAA-2007-4089,
AIAA-2009-3647.

**Licence.** The TMR repository is published under the Creative Commons
CC0 1.0 Universal public-domain dedication; the experiment was carried
out at NASA Ames Research Center. The files may be redistributed with
this repository (MIT); cite the paper above when you use them.

Used by `examples/studies/backstep_validation.py` (stage 6 part 6g,
specification 9.5).

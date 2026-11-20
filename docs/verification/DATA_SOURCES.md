# Validation reference data: sources and attribution

The validation cases of Foehncraft compare against reference numbers from
published papers and published data sets. This page lists each source, the
quantities used, and how it is cited. **The repository holds only the numbers
that the tests and the validation pages quote, with their citation. It does
not hold the data sets themselves**, and it holds no figures, tables or prose
from the papers.

Status of each entry: **verified** means the citation (DOI or publisher page)
was found in a bibliographic search; **to verify** means the DOI, page range
or archive link still has to be checked against the paper or the archive
before the entry is used in a test.

The source code and the validation scripts are under the MIT licence (root
`LICENSE`). The reference values are factual scientific data points. If you
use them in your own work, cite the original authors and papers below, not
Foehncraft alone.

## 1. Plane channel flow, Re_τ = 395 and 590

### 1.1 DNS reference (Moser, Kim & Mansour 1999)

- Citation: R. D. Moser, J. Kim, N. N. Mansour, *Direct numerical simulation
  of turbulent channel flow up to Re_τ = 590*, Physics of Fluids 11(4),
  943–945 (1999).
- DOI: 10.1063/1.869966 (**verified**)
- Archive: the DNS data page of the authors' group,
  https://turbulence.oden.utexas.edu/MKM_1999.html (**verified**: the page
  names the friction Reynolds numbers 180, 395 and 590, cites Phys. Fluids
  11(4), page 943, and offers the statistics as a gzipped tar file and a
  browsable tree)
- Re_τ = 395 profiles: https://turbulence.oden.utexas.edu/data/MKM/chan395/profiles/
  (**verified**: the listing has chan395.flat, chan395.means,
  chan395.reystress, chan395.skew, chan395.velp and chan395.vortvar)
- Re_τ = 590 profiles: the listing is **to verify**; the 2015 channel set
  (section 1.3) has no 395 or 590 case, so these values come from this archive
- Quantities used: mean streamwise velocity U⁺(y⁺); Reynolds stress
  components R₁₁⁺, R₂₂⁺, R₃₃⁺, R₁₂⁺; turbulent kinetic energy k⁺.

### 1.2 Model benchmark (Manceau 2015)

- Citation: R. Manceau, *Recent progress in the development of the Elliptic
  Blending Reynolds-stress model*, International Journal of Heat and Fluid
  Flow 51, 195–220 (2015).
- DOI: 10.1016/j.ijheatfluidflow.2014.09.002 (**verified**: the DOI resolves
  at doi.org, and library records give the same volume and pages; open
  copy: HAL hal-01092931)
- Quantities used: the model's own channel results at Re_τ = 395 and 590
  (mean velocity and stress profiles), read from the paper's tables or
  figures. Digitised values are recorded with the figure they come from.

### 1.3 Newer channel DNS set (Lee & Moser 2015)

- Citation: M. Lee, R. D. Moser, *Direct numerical simulation of turbulent
  channel flow up to Re_τ = 5200*, Journal of Fluid Mechanics 774, 395–415
  (2015).
- DOI: 10.1017/S0022112015002682 (**verified**: doi.org resolves to the
  publisher page for this article)
- Archive: the channel data page of the authors' group,
  https://turbulence.oden.utexas.edu/channel2015/ (the README asks that the
  data source be cited as above). The set has Re_τ = 180, 550, 1000, 2000 and
  5200. **It has no case at Re_τ = 395 or 590**, so the entries in section 1
  above keep their citation to Moser, Kim & Mansour (1999).
- Raw velocity fields are provided on request by the authors; statistics are
  published as text files on the archive. Nothing from the archive is copied
  into the repository. Only numbers quoted by a test or validation page, with
  this citation, may enter the repository.

## 2. Square duct (secondary flow, corner vortices)

### 2.1 Pinelli, Uhlmann, Sekimoto & Kawahara (2010)

- Citation: A. Pinelli, M. Uhlmann, A. Sekimoto, G. Kawahara, *Reynolds
  number dependence of mean flow structure in square duct turbulence*,
  Journal of Fluid Mechanics 644, 107–122 (2010).
- DOI: 10.1017/S0022112009992242 (**verified**: the publisher page carries
  this DOI, volume 644 and pages 107–122)
- Quantities used: mean streamwise velocity and the secondary velocity along
  the wall bisector and the corner diagonal, at the Reynolds numbers that
  the paper tabulates; the Re_b values to be taken from the paper (the
  abstract gives none).
- The data archive named in the draft (a KIT repository) has not been
  confirmed; it is removed until it is checked.

### 2.2 Gavrilakis (1992)

- Citation: S. Gavrilakis, *Numerical simulation of low-Reynolds-number
  turbulent flow through a straight square duct*, Journal of Fluid Mechanics
  244, 101–129 (1992).
- DOI: 10.1017/S0022112092002982 (**verified**: the publisher page carries
  this DOI, volume 244 and pages 101–129)
- Quantities used: secondary flow profiles and the wall shear distribution
  near the corners, as printed in the paper. The abstract states a bulk
  Reynolds number of 4410 (bulk velocity and hydraulic diameter); Re_τ is to
  be taken from the paper, not from this page.

## 3. Rotating channel

### 3.1 Kristoffersen & Andersson (1993)

- Citation: R. Kristoffersen, H. I. Andersson, *Direct simulations of
  low-Reynolds-number turbulent flow in a rotating channel*, Journal of Fluid
  Mechanics 256, 163–197 (1993).
- DOI: 10.1017/S0022112093002757 (**verified**, the publisher page carries
  this DOI)
- Quantities used: the asymmetric mean velocity profiles, the friction
  velocity ratio between the suction and the pressure side, and the
  Reynolds stresses across the channel, as printed in the paper for the
  rotation numbers it reports.
- The ERCOFTAC Classic Collection case is **not** used: its data carries
  CC BY-NC-SA 4.0, and it stays out of this repository. Values come from the
  paper.

## 4. Rules for adding an entry

1. A value goes into the repository only when a test or a validation page
   quotes it, with its citation next to it.
2. Digitised values from a figure are marked as such and name the figure.
3. Full data sets, archives and fetch scripts are not added to the
   repository. Data needed for a validation run is downloaded in the
   project for that run only, and stays out of the repository.
4. A new entry needs a verified DOI or publisher page, or it is marked
   **to verify** until it is.

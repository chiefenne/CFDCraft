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
  https://turbulence.oden.utexas.edu/MKM_1999.html (**verified** as the page
  of the data; the table layout to be checked when used)
- Quantities used: mean streamwise velocity U⁺(y⁺); Reynolds stress
  components R₁₁⁺, R₂₂⁺, R₃₃⁺, R₁₂⁺; turbulent kinetic energy k⁺.

### 1.2 Model benchmark (Manceau 2015)

- Citation: R. Manceau, *Recent progress in the development of the Elliptic
  Blending Reynolds-stress model*, International Journal of Heat and Fluid
  Flow 51, 195–220 (2015).
- DOI: to verify (the DOI in the draft could not be confirmed in a search)
- Quantities used: the model's own channel results at Re_τ = 395 and 590
  (mean velocity and stress profiles), read from the paper's tables or
  figures. Digitised values are recorded with the figure they come from.

## 2. Square duct (secondary flow, corner vortices)

### 2.1 Pinelli, Uhlmann, Sekimoto & Kawahara (2010)

- Citation: A. Pinelli, M. Uhlmann, A. Sekimoto, G. Kawahara, *Reynolds
  number dependence of mean flow structure in square duct turbulence*,
  Journal of Fluid Mechanics 644, 107–122 (2010).
- DOI: to verify (not confirmed in the search)
- Quantities used: mean streamwise velocity and the secondary velocity along
  the wall bisector and the corner diagonal, at the Reynolds numbers that
  the paper tabulates; the Re_b values to be taken from the paper.
- The data archive named in the draft (a KIT repository) has not been
  confirmed; it is removed until it is checked.

### 2.2 Gavrilakis (1992)

- Citation: S. Gavrilakis, *Numerical simulation of low-Reynolds-number
  turbulent flow through a straight square duct*, Journal of Fluid Mechanics
  244, 101–129 (1992).
- DOI: to verify (publisher page found, DOI not read from it)
- Quantities used: secondary flow profiles and the wall shear distribution
  near the corners, as printed in the paper. The Re_b and Re_τ of the case are
  to be taken from the paper, not from this page.

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

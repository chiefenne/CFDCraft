"""Make the CSV files of the impinging-jet study from the ERCOFTAC data.

    python3 examples/studies/impinging_jet/ercoftac_to_csv.py DIR [--re 23750]

DIR holds the files `ij2lr-*.dat` of the ERCOFTAC classic database, case 025
(impinging jet, H/D = 2, Re = 23 000). The script writes nusselt.csv and
velocity_r<s>.csv in the folder `data` next to it, in the format examples/studies/impinging_jet.py
reads (see the README of this folder). Nu/Re^0.7 of the source is multiplied by
RE^0.7 of the study (default 23 750), the normalisation of the source itself.
The files are not part of the repository: the licence of the source is open.
"""

import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent / "data"
RE_FILES = 23000.0
# (r/D, the two digits of the file name)
STATIONS = (("0.5", "05"), ("1", "10"), ("2", "20"), ("2.5", "25"), ("3", "30"))


def read(path: Path) -> list[tuple[float, float]]:
    """The two columns of a data file, without its `#` lines."""
    rows = []
    for line in path.read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            a, b = line.split()[:2]
            rows.append((float(a), float(b)))
    return rows


def main(argv: list[str]) -> None:
    if not argv or argv[0].startswith("--"):
        sys.exit(__doc__)
    src = Path(argv[0])
    re_run = float(argv[argv.index("--re") + 1]) if "--re" in argv else 23750.0
    scale = re_run**0.7
    with open(OUT / "nusselt.csv", "w") as f:
        f.write("# ERCOFTAC classic database, case 025 (impinging jet), file ij2lr-nuss.dat:\n")
        f.write(f"# H/D = 2, Re = {RE_FILES:.0f}, Nu/Re^0.7 of Baughn et al. (labelled 1992 in the file).\n")
        f.write(f"# Nu = (Nu/Re^0.7) * {re_run:.0f}^0.7 = (Nu/Re^0.7) * {scale:.2f}, the Re of the study.\n")
        f.write("r_D,Nu\n")
        for r, v in read(src / "ij2lr-nuss.dat"):
            f.write(f"{r:.3f},{v * scale:.3f}\n")
    for s, n in STATIONS:
        p = src / f"ij2lr-{n}-sw-mu.dat"
        with open(OUT / f"velocity_r{s}.csv", "w") as f:
            f.write(f"# ERCOFTAC classic database, case 025 (impinging jet), file {p.name}:\n")
            f.write(f"# H/D = 2, Re = {RE_FILES:.0f}, single-wire mean velocity U/U_bulk at r/D = {s}\n")
            f.write("# (Cooper et al., labelled 1992 in the file); x_D is the distance from the plate over D.\n")
            f.write("x_D,u_U\n")
            for x, u in read(p):
                f.write(f"{x:.6g},{u:.6g}\n")
    print(f"wrote nusselt.csv and {len(STATIONS)} velocity files to {OUT}")


if __name__ == "__main__":
    main(sys.argv[1:])

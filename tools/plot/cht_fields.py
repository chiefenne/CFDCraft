"""The first conjugate heat transfer case of stage 7d1 (gallery entry of
the PR): a heated wall cooled by a laminar channel flow. Left the
temperature in the wall and the fluid with the streamlines; right the
temperature of the wall's top cells along the channel for three
conductivity ratios k_s/k_f.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.cht_fields OUT_DIR RUN_PREFIX

with RUN_PREFIX the output directories `<prefix>_0.1`, `<prefix>_10` and
`<prefix>_100`, e.g.

    for r in 0.1 10 100; do
      k=$(python3 -c "print($r * 14.084507042253522)")
      target/release/foehn flow examples/cases/flow/cht_heated_wall.toml \\
          --out output/cht_$r --set "materials.wall-solid.conductivity=$k"
    done
    python3 -m tools.plot.cht_fields docs/figures/gallery output/cht

The figure is cht_heated_wall.png.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tools.plot.energy_fields import field, grid, line
from tools.plot.turbulence_fields import read


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    out, prefix = sys.argv[1:3]
    os.makedirs(out, exist_ok=True)
    fig, axes = plt.subplots(
        1, 2, figsize=(14.5, 4.4), constrained_layout=True, gridspec_kw={"width_ratios": [1.6, 1]}
    )
    quads, _, d = read(os.path.join(f"{prefix}_10", "fields.vtkhdf"))
    ax = axes[0]
    field(ax, quads, d["temperature"], "inferno", "temperature (K)", fig)
    xs, ys, ux = grid(quads, d["velocity_x"])
    _, _, uy = grid(quads, d["velocity_y"])
    fluid = ys > 0.2
    ax.streamplot(
        xs, ys[fluid], ux[fluid], uy[fluid], color="w", linewidth=0.5, density=0.9, arrowsize=0.6
    )
    ax.axhline(0.2, color="c", lw=1.0, ls="--")
    ax.text(3.95, 0.1, "solid wall, k_s/k_f = 10", color="w", ha="right", va="center", fontsize=9)
    ax.text(0.05, 1.12, "fluid, Re = 100, Pr = 0.71, 300 K in", color="w", va="top", fontsize=9)
    ax.set_title("Heated wall cooled by a channel flow: 1000 W/m² into the wall from below")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax = axes[1]
    for r, c in (("0.1", "#d95f02"), ("10", "#1b9e77"), ("100", "#7570b3")):
        p = os.path.join(f"{prefix}_{r}", "line_wall_top.csv")
        if not os.path.exists(p):
            continue
        v = line(p)
        ax.plot(v["x"], v["temperature"], "-", color=c, lw=1.6, label=f"k_s/k_f = {r}")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("temperature of the wall's top cells (K)")
    ax.set_title("The wall spreads the heat towards the inlet")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out, "cht_heated_wall.png"), dpi=120)


if __name__ == "__main__":
    main()

"""Turbulence values of the backflow at a pressure outlet (BC028, gallery
140): k behind the backward-facing step whose outlet sits inside the
recirculation, with the outlet's backflow values and without them (the
entering fluid then carries the k of the cell next to the outlet). The
black line is u = 0, the edge of the recirculation; the outlet is the
right edge of each panel.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.outlet_backflow OUT.png WITH.vtkhdf WITHOUT.vtkhdf

e.g.

    target/release/foehn flow examples/cases/flow/backstep_outlet_backflow.toml
    # the same case without the lines k and omega of [boundary.xmax]:
    target/release/foehn flow backstep_outlet_backflow_default.toml
    python3 -m tools.plot.outlet_backflow docs/figures/gallery/outlet_backflow.png \\
        output/backstep_outlet_backflow/fields.vtkhdf \\
        output/backstep_outlet_backflow_default/fields.vtkhdf
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

from tools.plot.turbulence_fields import read  # noqa: E402

WINDOW = (-1.0, 4.0, 0.0, 2.0)


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    out, with_values, without = sys.argv[1:]
    runs = [("outlet gives k and ω of the entering fluid", with_values),
            ("default: the entering fluid carries k and ω of the cell", without)]
    loaded = [read(p) for _, p in runs]
    quads = loaded[0][0]
    cen = quads.mean(axis=1)
    keep = (quads[:, :, 0].max(axis=1) > WINDOW[0]) & (quads[:, :, 1].min(axis=1) < WINDOW[3])
    hi = max(np.percentile(d["k"][keep], 99.5) for _, _, d in loaded)
    fig, axes = plt.subplots(2, 1, figsize=(6.2, 5.0), sharex=True, constrained_layout=True)
    for ax, (title, _), (q, _, d) in zip(axes, runs, loaded):
        pc = PolyCollection(q[keep], array=d["k"][keep], cmap="inferno", edgecolors="face",
                            linewidths=0.2)
        pc.set_clim(0.0, hi)
        ax.add_collection(pc)
        sel = keep & (cen[:, 0] > 0.0)
        ax.tricontour(cen[sel, 0], cen[sel, 1], d["velocity_x"][sel], levels=[0.0], colors="w", linewidths=1.2)
        ax.fill([WINDOW[0], 0, 0, WINDOW[0]], [0, 0, 1, 1], color="0.55", zorder=3)
        ax.set_xlim(WINDOW[0], WINDOW[1])
        ax.set_ylim(WINDOW[2], WINDOW[3])
        ax.set_aspect("equal")
        ax.set_ylabel("y/H")
        ax.set_title(title, fontsize=10)
        fig.colorbar(pc, ax=ax, label="k / U²", shrink=0.9, pad=0.01)
    axes[-1].set_xlabel("x/H  (outlet at x/H = 4)")
    fig.suptitle("Backflow at a pressure outlet, SST: k behind the step, white line u = 0",
                 fontsize=11)
    fig.savefig(out, dpi=110)
    print(out)


if __name__ == "__main__":
    main()

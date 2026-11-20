"""Writes examples/meshes/taylor_couette_zones.su2: the gap between two coaxial
cylinders as a two-zone SU2 file for a sliding-mesh run.

Zone 1 is the rotor (r 0.05 .. 0.075), zone 2 the stator (r 0.075 .. 0.10); each
has its own number of sectors (rotor 32, stator 48, so the interface does not
match), one layer of hexahedra in z (depth 0.01) and 8 radial layers.
The two zones do not share nodes: the interface is the marker rotor_interface
of zone 1 and the marker stator_interface of zone 2, as `[[sliding]]` wants it.
Markers: inner, outer, zmin, zmax (in both zones for the last two).

    python3 -m tools.mesh.taylor_couette_zones examples/meshes/taylor_couette_zones.su2
"""

import math
import sys

N_R = 8
DEPTH = 0.01


def zone(r0, r1, inner, outer, n_theta):
    N_THETA = n_theta
    nodes = []
    for k in (0.0, DEPTH):
        for j in range(N_R + 1):
            r = r0 + (r1 - r0) * j / N_R
            for i in range(N_THETA):
                t = 2.0 * math.pi * i / N_THETA
                nodes.append((r * math.cos(t), r * math.sin(t), k))

    def n(i, j, k):
        return (i % N_THETA) + N_THETA * (j + (N_R + 1) * k)

    cells, zmin, zmax = [], [], []
    for j in range(N_R):
        for i in range(N_THETA):
            q = [n(i, j, 0), n(i, j + 1, 0), n(i + 1, j + 1, 0), n(i + 1, j, 0)]
            cells.append(q + [v + N_THETA * (N_R + 1) for v in q])
            zmin.append(q)
            zmax.append([v + N_THETA * (N_R + 1) for v in q])
    wall_in = [[n(i, 0, 0), n(i + 1, 0, 0), n(i + 1, 0, 1), n(i, 0, 1)] for i in range(N_THETA)]
    wall_out = [[n(i, N_R, 0), n(i + 1, N_R, 0), n(i + 1, N_R, 1), n(i, N_R, 1)] for i in range(N_THETA)]
    return nodes, cells, {inner: wall_in, outer: wall_out, "zmin": zmin, "zmax": zmax}


def main(path):
    zones = [
        zone(0.05, 0.075, "inner", "rotor_interface", 32),
        zone(0.075, 0.10, "stator_interface", "outer", 48),
    ]
    out = [f"NZONE= {len(zones)}"]
    for z, (nodes, cells, markers) in enumerate(zones):
        out += [f"IZONE= {z + 1}", "NDIME= 3", f"NELEM= {len(cells)}"]
        out += [f"12 {' '.join(map(str, c))} {i}" for i, c in enumerate(cells)]
        out += [f"NPOIN= {len(nodes)}"]
        out += [f"{x:.12g} {y:.12g} {w:.12g} {i}" for i, (x, y, w) in enumerate(nodes)]
        out += [f"NMARK= {len(markers)}"]
        for name, faces in markers.items():
            out += [f"MARKER_TAG= {name}", f"MARKER_ELEMS= {len(faces)}"]
            out += [f"9 {' '.join(map(str, f))}" for f in faces]
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])

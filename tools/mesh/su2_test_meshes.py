"""Writes the small SU2 test meshes of crates/foehncraft-mesh/tests/data.

Plain Python, no dependencies. Run from the repository root:
    python3 -m tools.mesh.su2_test_meshes

- cubes_mixed_3d.su2: four unit cubes in a row along x: a hexahedron, two
  prisms, six pyramids about the cube centre, and one pyramid plus ten
  tetrahedra. Markers xmin, xmax and walls. Element and node indices at the
  ends of the lines, comments, and a trailing FFD section as SU2 files carry
  them.
- channel_2d.su2: the 2D rectangle [0, 2] x [0, 1], quadrilaterals on the
  left half and triangles on the right, nodes before elements (the order
  some writers use). Markers inlet, outlet, lower, upper.
"""
import os

OUT = "crates/foehncraft-mesh/tests/data"


def boundary_faces(cells, faces_of):
    """Faces used by one cell only, with the nodes in the cell's order."""
    seen = {}
    for c in cells:
        for f in faces_of(c):
            k = frozenset(f)
            seen.setdefault(k, []).append(f)
    return [v[0] for v in seen.values() if len(v) == 1]


def cubes_mixed_3d():
    pts = []
    index = {}

    def node(x, y, z):
        key = (x, y, z)
        if key not in index:
            index[key] = len(pts)
            pts.append(key)
        return index[key]

    cells = []  # (vtk type, nodes)
    for k in range(4):
        c = {(i, j, l): node(k + i, j, l) for i in (0, 1) for j in (0, 1) for l in (0, 1)}
        # the six faces of the cube as loops
        quads = [
            [c[0, 0, 0], c[0, 1, 0], c[0, 1, 1], c[0, 0, 1]],
            [c[1, 0, 0], c[1, 0, 1], c[1, 1, 1], c[1, 1, 0]],
            [c[0, 0, 0], c[0, 0, 1], c[1, 0, 1], c[1, 0, 0]],
            [c[0, 1, 0], c[1, 1, 0], c[1, 1, 1], c[0, 1, 1]],
            [c[0, 0, 0], c[1, 0, 0], c[1, 1, 0], c[0, 1, 0]],
            [c[0, 0, 1], c[0, 1, 1], c[1, 1, 1], c[1, 0, 1]],
        ]
        if k == 0:
            cells.append((12, [c[0, 0, 0], c[1, 0, 0], c[1, 1, 0], c[0, 1, 0],
                               c[0, 0, 1], c[1, 0, 1], c[1, 1, 1], c[0, 1, 1]]))
        elif k == 1:
            for t in ([(0, 0), (1, 0), (1, 1)], [(0, 0), (1, 1), (0, 1)]):
                cells.append((13, [c[i, j, 0] for i, j in t] + [c[i, j, 1] for i, j in t]))
        else:
            centre = node(k + 0.5, 0.5, 0.5)
            if k == 2:
                for q in quads:
                    cells.append((14, q + [centre]))
            else:
                cells.append((14, quads[0] + [centre]))
                for q in quads[1:]:
                    cells.append((10, [q[0], q[1], q[2], centre]))
                    cells.append((10, [q[0], q[2], q[3], centre]))

    def faces_of(cell):
        t, v = cell
        if t == 10:
            return [[v[0], v[2], v[1]], [v[0], v[1], v[3]], [v[0], v[3], v[2]], [v[1], v[2], v[3]]]
        if t == 12:
            return [[v[0], v[3], v[2], v[1]], [v[4], v[5], v[6], v[7]], [v[0], v[1], v[5], v[4]],
                    [v[1], v[2], v[6], v[5]], [v[2], v[3], v[7], v[6]], [v[3], v[0], v[4], v[7]]]
        if t == 13:
            return [[v[0], v[2], v[1]], [v[3], v[4], v[5]], [v[0], v[1], v[4], v[3]],
                    [v[1], v[2], v[5], v[4]], [v[2], v[0], v[3], v[5]]]
        return [[v[0], v[3], v[2], v[1]], [v[0], v[1], v[4]], [v[1], v[2], v[4]],
                [v[2], v[3], v[4]], [v[3], v[0], v[4]]]

    markers = {"xmin": [], "xmax": [], "walls": []}
    for f in boundary_faces(cells, faces_of):
        xs = [pts[n][0] for n in f]
        name = "xmin" if max(xs) == 0 else "xmax" if min(xs) == 4 else "walls"
        markers[name].append(f)

    with open(os.path.join(OUT, "cubes_mixed_3d.su2"), "w") as o:
        o.write("% Four unit cubes: hexahedron, prisms, pyramids, tetrahedra\n")
        o.write("% (tools/mesh/su2_test_meshes.py)\n%\nNDIME= 3\n%\n")
        o.write(f"NELEM= {len(cells)}\n")
        for i, (t, v) in enumerate(cells):
            o.write(f"{t}\t" + "\t".join(str(n) for n in v) + f"\t{i}\n")
        o.write(f"NPOIN= {len(pts)}\n")
        for i, p in enumerate(pts):
            o.write(f"  {p[0]:.15e}  {p[1]:.15e}  {p[2]:.15e}  {i}\n")
        o.write(f"NMARK= {len(markers)}\n")
        for name, fs in markers.items():
            o.write(f"MARKER_TAG= {name}\nMARKER_ELEMS= {len(fs)}\n")
            for f in fs:
                o.write(f"{5 if len(f) == 3 else 9} " + " ".join(str(n) for n in f) + "\n")
        o.write("%\n% FFD box data follow the mesh in some files\nFFD_NBOX= 0\nFFD_NLEVEL= 0\n")


def channel_2d():
    nx, ny = 4, 2
    pts = [(2 * i / nx, j / ny) for j in range(ny + 1) for i in range(nx + 1)]

    def n(i, j):
        return j * (nx + 1) + i

    cells = []
    for j in range(ny):
        for i in range(nx):
            a, b, c, d = n(i, j), n(i + 1, j), n(i + 1, j + 1), n(i, j + 1)
            if i < nx // 2:
                cells.append((9, [a, b, c, d]))
            else:
                cells += [(5, [a, b, c]), (5, [a, c, d])]
    markers = {
        "inlet": [[n(0, j + 1), n(0, j)] for j in range(ny)],
        "outlet": [[n(nx, j), n(nx, j + 1)] for j in range(ny)],
        "lower": [[n(i, 0), n(i + 1, 0)] for i in range(nx)],
        "upper": [[n(i + 1, ny), n(i, ny)] for i in range(nx)],
    }
    with open(os.path.join(OUT, "channel_2d.su2"), "w") as o:
        o.write("%\n% 2D channel, quadrilaterals and triangles\n%\nNDIME= 2\n")
        o.write(f"NPOIN= {len(pts)}\n")
        for i, p in enumerate(pts):
            o.write(f"{p[0]: .8e} {p[1]: .8e} {i}\n")
        o.write(f"NELEM= {len(cells)}\n")
        for i, (t, v) in enumerate(cells):
            o.write(f"{t} " + " ".join(str(x) for x in v) + f" {i}\n")
        o.write(f"NMARK= {len(markers)}\n")
        for name, es in markers.items():
            o.write(f"MARKER_TAG= {name}\nMARKER_ELEMS= {len(es)}\n")
            for e in es:
                o.write(f"3 {e[0]} {e[1]}\n")


if __name__ == "__main__":
    cubes_mixed_3d()
    channel_2d()

"""Writes the multi-zone test meshes of crates/foehncraft-mesh/tests/data/multizone:
a rotor block (x 0..1) and a stator block (x 1..2), each a unit-section hexahedral
zone, as SU2 files (NZONE= 2) and as a two-zone CGNS file.

Needs numpy and h5py (CGNS only). Run from the repository root:
    python3 -m tools.mesh.multizone_fixtures crates/foehncraft-mesh/tests/data/multizone

Files:
    rotor_stator.su2        interface 2 x 2 faces on the rotor side, 3 x 3 on the
                            stator side (non-matching, as in a sliding run);
                            markers inlet, rotor_interface, stator_interface,
                            outlet, wall (the wall marker is in both zones)
    rotor_stator_conf.su2   both sides 2 x 2 (conforming), markers iface_a and iface_b
    rotor_stator.cgns       conforming, zones rotor and stator, BCs inlet,
                            rotor_interface, stator_interface, outlet, wall
"""

import sys
from pathlib import Path

import numpy as np


def block(x0, x1, nx, ny, nz):
    """Hexahedra of the box [x0, x1] x [0, 1] x [0, 1]; returns the nodes, the
    cells and the faces of the sides xmin, xmax and the four walls."""
    pts = [
        (x0 + (x1 - x0) * i / nx, j / ny, k / nz)
        for k in range(nz + 1)
        for j in range(ny + 1)
        for i in range(nx + 1)
    ]

    def n(i, j, k):
        return i + (nx + 1) * (j + (ny + 1) * k)

    cells = [
        [
            n(i, j, k),
            n(i + 1, j, k),
            n(i + 1, j + 1, k),
            n(i, j + 1, k),
            n(i, j, k + 1),
            n(i + 1, j, k + 1),
            n(i + 1, j + 1, k + 1),
            n(i, j + 1, k + 1),
        ]
        for k in range(nz)
        for j in range(ny)
        for i in range(nx)
    ]
    side = {"xmin": [], "xmax": [], "wall": []}
    for k in range(nz):
        for j in range(ny):
            side["xmin"].append([n(0, j, k), n(0, j + 1, k), n(0, j + 1, k + 1), n(0, j, k + 1)])
            side["xmax"].append([n(nx, j, k), n(nx, j + 1, k), n(nx, j + 1, k + 1), n(nx, j, k + 1)])
    for k in range(nz):
        for i in range(nx):
            side["wall"].append([n(i, 0, k), n(i + 1, 0, k), n(i + 1, 0, k + 1), n(i, 0, k + 1)])
            side["wall"].append([n(i, ny, k), n(i + 1, ny, k), n(i + 1, ny, k + 1), n(i, ny, k + 1)])
    for j in range(ny):
        for i in range(nx):
            side["wall"].append([n(i, j, 0), n(i + 1, j, 0), n(i + 1, j + 1, 0), n(i, j + 1, 0)])
            side["wall"].append([n(i, j, nz), n(i + 1, j, nz), n(i + 1, j + 1, nz), n(i, j + 1, nz)])
    return pts, cells, side


def write_su2(path, zones):
    """zones: list of (nodes, cells, {marker: faces})."""
    out = [f"NZONE= {len(zones)}"]
    for z, (pts, cells, markers) in enumerate(zones):
        out += [f"IZONE= {z + 1}", "NDIME= 3", f"NELEM= {len(cells)}"]
        out += [f"12 {' '.join(map(str, c))} {i}" for i, c in enumerate(cells)]
        out += [f"NPOIN= {len(pts)}"]
        out += [f"{x} {y} {zz} {i}" for i, (x, y, zz) in enumerate(pts)]
        out += [f"NMARK= {len(markers)}"]
        for name, faces in markers.items():
            out += [f"MARKER_TAG= {name}", f"MARKER_ELEMS= {len(faces)}"]
            out += [f"9 {' '.join(map(str, f))}" for f in faces]
    Path(path).write_text("\n".join(out) + "\n")


def write_cgns(path, zones):
    """zones: list of (name, nodes, cells, {BC: faces})."""
    import h5py

    from tools.mesh.cgns_ngon_fixture import node, text

    f = h5py.File(path, "w")
    f.create_dataset(" format", data=text("IEEE_LITTLE_32"))
    f.create_dataset(" hdf5version", data=text("HDF5 Version 1.14.0"))
    node(f, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4", np.array([4.0], dtype=np.float32))
    base = node(f, "Base", "CGNSBase_t", "I4", np.array([3, 3], dtype=np.int32))
    for name, pts, cells, bcs in zones:
        z = node(base, name, "Zone_t", "I4", np.array([[len(pts)], [len(cells)], [0]], dtype=np.int32))
        node(z, "ZoneType", "ZoneType_t", "C1", text("Unstructured"))
        gc = node(z, "GridCoordinates", "GridCoordinates_t", "MT")
        a = np.array(pts, dtype=float)
        for i, c in enumerate("XYZ"):
            node(gc, f"Coordinate{c}", "DataArray_t", "R8", a[:, i].copy())

        def section(sname, kind, first, items):
            s = node(z, sname, "Elements_t", "I4", np.array([kind, 0], dtype=np.int32))
            node(s, "ElementRange", "IndexRange_t", "I4", np.array([first, first + len(items) - 1], dtype=np.int32))
            conn = (np.array(items) + 1).astype(np.int32).ravel()
            node(s, "ElementConnectivity", "DataArray_t", "I4", conn)

        section("Hexa", 17, 1, cells)
        zbc = node(z, "ZoneBC", "ZoneBC_t", "MT")
        first = len(cells) + 1
        faces = [x for fs in bcs.values() for x in fs]
        section("Quad", 7, first, faces)
        for bname, fs in bcs.items():
            b = node(zbc, bname, "BC_t", "C1", text("BCWall"))
            node(b, "PointRange", "IndexRange_t", "I4", np.array([[first], [first + len(fs) - 1]], dtype=np.int32))
            node(b, "GridLocation", "GridLocation_t", "C1", text("FaceCenter"))
            first += len(fs)
    f.close()


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    # non-matching interface: rotor 2 x 2, stator 3 x 3
    rp, rc, rs = block(0.0, 1.0, 2, 2, 2)
    sp, sc, ss = block(1.0, 2.0, 3, 3, 3)
    write_su2(
        out / "rotor_stator.su2",
        [
            (rp, rc, {"inlet": rs["xmin"], "rotor_interface": rs["xmax"], "wall": rs["wall"]}),
            (sp, sc, {"stator_interface": ss["xmin"], "outlet": ss["xmax"], "wall": ss["wall"]}),
        ],
    )
    # conforming
    rp, rc, rs = block(0.0, 1.0, 2, 2, 2)
    sp, sc, ss = block(1.0, 2.0, 2, 2, 2)
    write_su2(
        out / "rotor_stator_conf.su2",
        [
            (rp, rc, {"inlet": rs["xmin"], "iface_a": rs["xmax"], "wall": rs["wall"]}),
            (sp, sc, {"iface_b": ss["xmin"], "outlet": ss["xmax"], "wall": ss["wall"]}),
        ],
    )
    write_cgns(
        out / "rotor_stator.cgns",
        [
            ("rotor", rp, rc, {"inlet": rs["xmin"], "rotor_interface": rs["xmax"], "wall": rs["wall"]}),
            ("stator", sp, sc, {"stator_interface": ss["xmin"], "outlet": ss["xmax"], "wall": ss["wall"]}),
        ],
    )


if __name__ == "__main__":
    main(sys.argv[1])

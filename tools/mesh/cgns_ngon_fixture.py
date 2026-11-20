"""Writes a polyhedral CGNS test file (NGON_n and NFACE_n sections) from an
ASCII OpenFOAM polyMesh directory.

Needs numpy and h5py. Run from the repository root:
    python3 -m tools.mesh.cgns_ngon_fixture \
        crates/foehncraft-mesh/tests/data/openfoam/dual_polyhedra/constant/polyMesh \
        crates/foehncraft-mesh/tests/data/cgns/dual_ngon.cgns

CGNS 4 layout: the faces are one NGON_n section with ElementStartOffset, the
cells one NFACE_n section whose connectivity holds signed face numbers (positive:
the face normal points out of the cell). Each patch becomes one BC of the same
name, given as a face range.
"""
import re
import sys

import h5py
import numpy as np


def body(path):
    text = open(path).read()
    text = re.sub(r"//.*", "", text)
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = text[text.index("}", text.index("FoamFile")) + 1:]
    return text


def lists(path):
    """Numbers of an ASCII list file: N ( ... ), faces as n(a b c)."""
    t = body(path)
    return t


def foam_numbers(path):
    t = body(path)
    n, rest = re.match(r"\s*(\d+)\s*\((.*)\)\s*$", t, re.S).groups()
    return int(n), rest


def node(parent, name, label, kind, data=None):
    g = parent.create_group(name)
    g.attrs["label"] = np.bytes_(label)
    g.attrs["name"] = np.bytes_(name)
    g.attrs["type"] = np.bytes_(kind)
    g.attrs["flags"] = np.array([1], dtype=np.int32)
    if data is not None:
        g.create_dataset(" data", data=data)
    return g


def text(s):
    return np.frombuffer(s.encode(), dtype=np.int8)


def main(src, dst):
    n, rest = foam_numbers(f"{src}/points")
    pts = np.array(re.findall(r"\(([^()]+)\)", rest), dtype=object)
    pts = np.array([[float(v) for v in p.split()] for p in pts])
    n, rest = foam_numbers(f"{src}/faces")
    faces = [[int(v) + 1 for v in f.split()] for f in re.findall(r"\d+\s*\(([^()]+)\)", rest)]
    assert len(faces) == n
    owner = np.array(foam_numbers(f"{src}/owner")[1].split(), dtype=int)
    neigh = np.array(foam_numbers(f"{src}/neighbour")[1].split(), dtype=int)
    n_cells = max(owner.max(), neigh.max()) + 1
    cell_faces = [[] for _ in range(n_cells)]
    for f, c in enumerate(owner):
        cell_faces[c].append(f + 1)
    for f, c in enumerate(neigh):
        cell_faces[c].append(-(f + 1))
    patches = re.findall(
        r"(\w+)\s*\{[^}]*?nFaces\s+(\d+);\s*startFace\s+(\d+);", foam_numbers(f"{src}/boundary")[1] if False else body(f"{src}/boundary")
    )

    f = h5py.File(dst, "w")
    f.create_dataset(" format", data=text("IEEE_LITTLE_32"))
    f.create_dataset(" hdf5version", data=text("HDF5 Version 1.14.0"))
    node(f, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4", np.array([4.0], dtype=np.float32))
    base = node(f, "Base", "CGNSBase_t", "I4", np.array([3, 3], dtype=np.int32))
    zone = node(base, "Zone", "Zone_t", "I4", np.array([[len(pts)], [n_cells], [0]], dtype=np.int32))
    node(zone, "ZoneType", "ZoneType_t", "C1", text("Unstructured"))
    gc = node(zone, "GridCoordinates", "GridCoordinates_t", "MT")
    for i, c in enumerate("XYZ"):
        node(gc, f"Coordinate{c}", "DataArray_t", "R8", pts[:, i].copy())

    def section(name, kind, first, items):
        s = node(zone, name, "Elements_t", "I4", np.array([kind, 0], dtype=np.int32))
        node(s, "ElementRange", "IndexRange_t", "I4", np.array([first, first + len(items) - 1], dtype=np.int32))
        off = np.concatenate([[0], np.cumsum([len(x) for x in items])]).astype(np.int32)
        node(s, "ElementStartOffset", "DataArray_t", "I4", off)
        node(s, "ElementConnectivity", "DataArray_t", "I4", np.concatenate(items).astype(np.int32))

    section("NGON", 22, 1, faces)
    section("NFACE", 23, len(faces) + 1, cell_faces)
    bc = node(zone, "ZoneBC", "ZoneBC_t", "MT")
    for name, nf, start in patches:
        nf, start = int(nf), int(start)
        b = node(bc, name, "BC_t", "C1", text("BCWall"))
        node(b, "PointRange", "IndexRange_t", "I4", np.array([[start + 1], [start + nf]], dtype=np.int32))
        node(b, "GridLocation", "GridLocation_t", "C1", text("FaceCenter"))
    f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

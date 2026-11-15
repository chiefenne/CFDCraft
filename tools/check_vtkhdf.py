"""Reads Inncraft VTKHDF files with VTK's own reader (the reader ParaView
uses) and checks them:

- the file loads, with the block `internal` and one block per patch;
- the cell volumes agree with the field `volume` written by Inncraft. For
  standard cells VTK's own cell volumes are used (this checks the node order
  of tetrahedra, hexahedra, prisms and pyramids); for polyhedra the volume is
  computed here from the faces as VTK reads them, by the divergence theorem
  (this checks that every face is there and points outward; VTK's own
  polyhedron volume is the volume of the convex hull of the nodes, which
  differs for warped faces).

    pip install vtk
    python3 tools/check_vtkhdf.py FILE.vtkhdf [...]
"""
import sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy as v2n


def polyhedron_volume(cell):
    vol = 0.0
    for i in range(cell.GetNumberOfFaces()):
        f = cell.GetFace(i)
        p = np.array([f.GetPoints().GetPoint(k) for k in range(f.GetNumberOfPoints())])
        m = p.mean(axis=0)
        for k in range(len(p)):
            q = p[(k + 1) % len(p)]
            s = 0.5 * np.cross(p[k] - m, q - m)
            vol += np.dot((p[k] + q + m) / 3.0, s) / 3.0
    return vol


def check(path):
    r = vtk.vtkHDFReader()
    r.SetFileName(path)
    r.Update()
    out = r.GetOutputDataObject(0)
    assert out.IsA("vtkPartitionedDataSetCollection"), out.GetClassName()
    names = [out.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME())
             for i in range(out.GetNumberOfPartitionedDataSets())]
    grid = out.GetPartitionedDataSet(names.index("internal")).GetPartition(0)
    n = grid.GetNumberOfCells()
    sizes = vtk.vtkCellSizeFilter()
    sizes.SetInputData(grid)
    sizes.Update()
    vol = v2n(sizes.GetOutput().GetCellData().GetArray("Volume")).copy()
    types = np.array([grid.GetCellType(i) for i in range(n)])
    for i in np.nonzero(types == vtk.VTK_POLYHEDRON)[0]:
        vol[i] = polyhedron_volume(grid.GetCell(i))
    msg = f"{path}: {n} cells, {len(names) - 1} patches"
    ours = grid.GetCellData().GetArray("volume")
    if ours is not None:
        ours = v2n(ours)
        err = abs(vol - ours).max() / abs(ours).max()
        msg += f", max relative volume difference {err:.1e}"
        assert err < 1e-8, msg
    msg += f", VTK cell types {sorted(set(types.tolist()))}"
    print(msg)


for p in sys.argv[1:]:
    check(p)

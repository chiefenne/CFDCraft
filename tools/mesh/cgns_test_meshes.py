"""Writes the Gmsh CGNS test meshes of crates/foehncraft-mesh/tests/data/cgns.

Needs the gmsh Python module (pip install gmsh). Run in the target directory:
    python3 /path/to/tools/mesh/cgns_test_meshes.py
Writes box_hex.cgns, box_tet.cgns and box_prism.cgns (unit cube, physical
groups xmin..zmax and fluid).
"""

import gmsh

gmsh.initialize(["-v", "0"])


def box(name, kind):
    gmsh.clear()
    gmsh.model.add(name)
    occ = gmsh.model.occ
    if kind == "prism":
        pts = [occ.addPoint(x, y, 0) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
        ls = [occ.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)]
        s = occ.addPlaneSurface([occ.addCurveLoop(ls)])
        out = occ.extrude([(2, s)], 0, 0, 1, numElements=[3], recombine=True)
        occ.synchronize()
        vol = next(t for d, t in out if d == 3)
    else:
        vol = occ.addBox(0, 0, 0, 1, 1, 1)
        occ.synchronize()
    names = {}
    for dim, tag in gmsh.model.getBoundary([(3, vol)], oriented=False):
        x, y, z = occ.getCenterOfMass(dim, tag)
        for c, lo, hi in (
            (x, "xmin", "xmax"),
            (y, "ymin", "ymax"),
            (z, "zmin", "zmax"),
        ):
            if abs(c) < 1e-9:
                names[tag] = lo
            if abs(c - 1) < 1e-9:
                names[tag] = hi
    for tag, n in sorted(names.items()):
        gmsh.model.addPhysicalGroup(2, [tag], name=n)
    gmsh.model.addPhysicalGroup(3, [vol], name="fluid")
    if kind == "hex":
        for d, t in gmsh.model.getEntities(1):
            gmsh.model.mesh.setTransfiniteCurve(t, 3)
        for d, t in gmsh.model.getEntities(2):
            gmsh.model.mesh.setTransfiniteSurface(t)
            gmsh.model.mesh.setRecombine(2, t)
        gmsh.model.mesh.setTransfiniteVolume(vol)
    else:
        gmsh.option.setNumber("Mesh.MeshSizeFactor", 2.5)
    gmsh.model.mesh.generate(3)


for kind in ["hex", "tet", "prism"]:
    box("box_" + kind, kind)
    gmsh.write(f"box_{kind}.cgns")

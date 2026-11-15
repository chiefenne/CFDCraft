"""Writes the small Gmsh test meshes of crates/inncraft-mesh/tests/data.

Needs the Gmsh Python module (pip install gmsh). Run from the repository
root:  python3 tools/gmsh_test_meshes.py
"""
import gmsh, os

OUT = "crates/inncraft-mesh/tests/data"


def save(name, binary):
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 1 if binary else 0)
    gmsh.write(os.path.join(OUT, name))


def box_tet():
    """Unit box, tetrahedra, six named patches, one fluid zone."""
    gmsh.model.add("box_tet")
    v = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()
    names = {}
    for dim, tag in gmsh.model.getBoundary([(3, v)], oriented=False):
        x, y, z = gmsh.model.occ.getCenterOfMass(dim, tag)
        for c, lo, hi in ((x, "xmin", "xmax"), (y, "ymin", "ymax"), (z, "zmin", "zmax")):
            if abs(c) < 1e-9: names[tag] = lo
            if abs(c - 1) < 1e-9: names[tag] = hi
    for tag, n in sorted(names.items()):
        gmsh.model.addPhysicalGroup(2, [tag], name=n)
    gmsh.model.addPhysicalGroup(3, [v], name="fluid")
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.25)
    gmsh.model.mesh.generate(3)
    save("box_tet_ascii.msh", False)
    save("box_tet_binary.msh", True)


def channel_mixed():
    """Extruded and recombined channel: hexahedra and prisms in two volumes
    (two zones), with a periodic pair of z planes."""
    gmsh.model.add("channel_mixed")
    occ = gmsh.model.geo
    p = [occ.addPoint(x, y, 0) for x, y in ((0, 0), (1, 0), (1, 1), (0, 1), (2, 0), (2, 1))]
    l = [occ.addLine(p[0], p[1]), occ.addLine(p[1], p[2]), occ.addLine(p[2], p[3]), occ.addLine(p[3], p[0]),
         occ.addLine(p[1], p[4]), occ.addLine(p[4], p[5]), occ.addLine(p[5], p[2])]
    s1 = occ.addPlaneSurface([occ.addCurveLoop([l[0], l[1], l[2], l[3]])])
    s2 = occ.addPlaneSurface([occ.addCurveLoop([l[4], l[5], l[6], -l[1]])])
    for c in l:
        occ.mesh.setTransfiniteCurve(c, 5)
    occ.mesh.setTransfiniteSurface(s1)
    occ.mesh.setRecombine(2, s1)  # quads -> hexahedra; s2 stays triangles -> prisms
    e1 = occ.extrude([(2, s1)], 0, 0, 0.5, numElements=[2], recombine=True)
    e2 = occ.extrude([(2, s2)], 0, 0, 0.5, numElements=[2], recombine=True)
    occ.synchronize()
    gmsh.model.addPhysicalGroup(3, [e1[1][1]], name="left")
    gmsh.model.addPhysicalGroup(3, [e2[1][1]], name="right")
    gmsh.model.addPhysicalGroup(2, [s1, s2], name="back")
    gmsh.model.addPhysicalGroup(2, [e1[0][1], e2[0][1]], name="front")
    gmsh.model.addPhysicalGroup(2, [e1[5][1]], name="inlet")
    gmsh.model.addPhysicalGroup(2, [e2[3][1]], name="outlet")
    gmsh.model.addPhysicalGroup(2, [e1[2][1], e1[4][1], e2[2][1], e2[4][1]], name="walls")
    gmsh.model.mesh.setPeriodic(2, [e1[0][1]], [s1], [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0.5, 0, 0, 0, 1])
    gmsh.model.mesh.generate(3)
    save("channel_mixed_ascii.msh", False)
    save("channel_mixed_binary.msh", True)


def pyramid_block():
    """A hexahedral block next to a tetrahedral block, joined by pyramids
    (Gmsh's transition), all patches unnamed (no physical surfaces)."""
    gmsh.model.add("pyramids")
    a = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    b = gmsh.model.occ.addBox(1, 0, 0, 1, 1, 1)
    gmsh.model.occ.fragment([(3, a)], [(3, b)])
    gmsh.model.occ.synchronize()
    for c in gmsh.model.getEntities(1):
        gmsh.model.mesh.setTransfiniteCurve(c[1], 4)
    for s in gmsh.model.getBoundary([(3, a)], oriented=False):
        gmsh.model.mesh.setTransfiniteSurface(s[1])
        gmsh.model.mesh.setRecombine(2, s[1])
    gmsh.model.mesh.setTransfiniteVolume(a)
    gmsh.model.mesh.generate(3)
    save("pyramids_ascii.msh", False)


gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
os.makedirs(OUT, exist_ok=True)
box_tet()
channel_mixed()
pyramid_block()
gmsh.finalize()

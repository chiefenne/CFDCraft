"""A vessel with two stirrers in 2D (one cell layer): two rotors with their
own blades, each inside its own interface circle, in a stator that is the
rest of the vessel. The mesh of the sliding mesh with several rotors by
imprinting (D182 variant B, D185): each rotor turns at its own speed about
its own axis.

    pip install gmsh
    python3 examples/meshes/twin_mixer.py              # default size
    python3 examples/meshes/twin_mixer.py --size 1.2e-3
    python3 examples/meshes/twin_mixer.py --right-blades 5   # twin_mixer5.msh

Writes `twin_mixer.msh` next to this script: the vessel is a stadium (two
half circles of radius 0.10 m joined by straight walls, 0.36 m long), the
rotors sit at x = -0.08 m (zone `rotor_left`, three blades) and x = +0.08 m
(zone `rotor_right`, four blades), and each meets the stator on its circle
of radius 0.06 m with its own nodes (64 rotor and 96 stator intervals).
`--right-blades 5` writes `twin_mixer5.msh` with five blades on the right
rotor: the mesh of the gallery images (CONTRIBUTING.md, "Images": no
four-fold rotor patterns in published images).

Boundaries: `blades_left`, `blades_right` (blades and hubs, turning with
their rotors), `vessel` (the vessel wall), `left_rotor_interface`,
`left_stator_interface`, `right_rotor_interface`, `right_stator_interface`,
`zmin`, `zmax`. Lengths in m.
"""
import argparse
import math
import os

import gmsh

R_VESSEL = 0.10      # radius of the half circles of the vessel
X_ROTOR = 0.08       # distance of the axes from the middle
R_IFACE = 0.06       # interface radius
R_HUB = 0.010        # hub radius
R_TIP = 0.045        # blade tip radius
T_BLADE = 0.005      # blade thickness
DEPTH = 0.005        # thickness of the cell layer
N_ROTOR = 64         # intervals of an interface circle, rotor side
N_STATOR = 96        # stator side
ROTORS = [("left", -X_ROTOR, 3), ("right", X_ROTOR, 4)]


def rotor(x0, blades):
    occ = gmsh.model.occ
    disk = (2, occ.addDisk(x0, 0, 0, R_IFACE, R_IFACE))
    tools = [(2, occ.addDisk(x0, 0, 0, R_HUB, R_HUB))]
    for k in range(blades):
        s = occ.addRectangle(0.5 * R_HUB, -0.5 * T_BLADE, 0.0, R_TIP - 0.5 * R_HUB, T_BLADE)
        occ.rotate([(2, s)], 0, 0, 0, 0, 0, 1, 2.0 * math.pi * k / blades)
        occ.translate([(2, s)], x0, 0, 0)
        tools.append((2, s))
    out, _ = occ.cut([disk], tools)
    return out


def stator():
    occ = gmsh.model.occ
    a = (2, occ.addDisk(-X_ROTOR, 0, 0, R_VESSEL, R_VESSEL))
    b = (2, occ.addDisk(X_ROTOR, 0, 0, R_VESSEL, R_VESSEL))
    mid = (2, occ.addRectangle(-X_ROTOR, -R_VESSEL, 0, 2 * X_ROTOR, 2 * R_VESSEL))
    vessel, _ = occ.fuse([a], [b, mid])
    holes = [(2, occ.addDisk(x0, 0, 0, R_IFACE, R_IFACE)) for _, x0, _ in ROTORS]
    out, _ = occ.cut(vessel, holes)
    return out


def on_circle(dim, tag, x0, r):
    """True for a curve or surface that lies on the circle (cylinder) r
    about (x0, 0)."""
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim, tag)
    return (abs(xmax - x0 - r) < 1e-6 and abs(xmin - x0 + r) < 1e-6
            and abs(ymax - r) < 1e-6 and abs(ymin + r) < 1e-6)


def which_rotor(dim, tag):
    """The rotor whose interface disk holds the entity, or None."""
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim, tag)
    for name, x0, _ in ROTORS:
        if xmin > x0 - R_IFACE - 1e-6 and xmax < x0 + R_IFACE + 1e-6 and ymax < R_IFACE + 1e-6:
            return name
    return None


def build(size, path):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("twin_mixer")
    occ = gmsh.model.occ
    rotor_surfs = {name: rotor(x0, n) for name, x0, n in ROTORS}
    stator_surfs = stator()
    occ.synchronize()
    # interface curves: a fixed number of intervals on each side
    for name, x0, _ in ROTORS:
        for s, n in [(rotor_surfs[name], N_ROTOR), (stator_surfs, N_STATOR)]:
            for c in gmsh.model.getBoundary(s, oriented=False):
                if on_circle(1, c[1], x0, R_IFACE):
                    gmsh.model.mesh.setTransfiniteCurve(c[1], n + 1)
    surfaces = gmsh.model.getEntities(2)
    gmsh.option.setNumber("Mesh.MeshSizeMax", size)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.25 * size)
    # finer near the blades
    blades = [c[1] for name in rotor_surfs
              for c in gmsh.model.getBoundary(rotor_surfs[name], oriented=False)
              if not any(on_circle(1, c[1], x0, R_IFACE) for _, x0, _ in ROTORS)]
    f = gmsh.model.mesh.field
    d = f.add("Distance")
    f.setNumbers(d, "CurvesList", blades)
    f.setNumber(d, "Sampling", 200)
    t = f.add("Threshold")
    f.setNumber(t, "InField", d)
    f.setNumber(t, "SizeMin", 0.3 * size)
    f.setNumber(t, "SizeMax", size)
    f.setNumber(t, "DistMin", 0.002)
    f.setNumber(t, "DistMax", 0.02)
    f.setAsBackgroundMesh(t)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 8)          # frontal-Delaunay for quads
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
    # one layer of cells along z
    occ.extrude(surfaces, 0, 0, DEPTH, numElements=[1], recombine=True)
    occ.synchronize()
    vols = {"rotor_left": [], "rotor_right": [], "stator": []}
    for e in gmsh.model.getEntities(3):
        r = which_rotor(*e)
        vols["stator" if r is None else "rotor_" + r].append(e[1])
    groups = {"blades_left": [], "blades_right": [], "vessel": [], "zmin": [], "zmax": [],
              "left_rotor_interface": [], "left_stator_interface": [],
              "right_rotor_interface": [], "right_stator_interface": []}
    seen = {}
    for name, vs in vols.items():
        for v in vs:
            for b in gmsh.model.getBoundary([(3, v)], oriented=False):
                seen.setdefault(b[1], []).append(name)
    for tag, owners in seen.items():
        bb = gmsh.model.getBoundingBox(2, tag)
        owner = owners[0]
        if bb[5] < 1e-6:
            groups["zmin"].append(tag)
        elif bb[2] > DEPTH - 1e-6:
            groups["zmax"].append(tag)
        else:
            iface = [n for n, x0, _ in ROTORS if on_circle(2, tag, x0, R_IFACE)]
            if iface:
                side = "stator" if owner == "stator" else "rotor"
                groups[f"{iface[0]}_{side}_interface"].append(tag)
            elif owner == "stator":
                groups["vessel"].append(tag)
            else:
                groups["blades_" + owner.split("_")[1]].append(tag)
    for name, tags in groups.items():
        if tags:
            gmsh.model.addPhysicalGroup(2, tags, name=name)
    for name, vs in vols.items():
        gmsh.model.addPhysicalGroup(3, vs, name=name)
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(path)
    n = sum(len(gmsh.model.mesh.getElements(3, v)[1][0]) for vs in vols.values() for v in vs)
    gmsh.finalize()
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--size", type=float, default=2.5e-3, help="largest cell size (m)")
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--right-blades", type=int, default=4,
                    help="number of blades of the right rotor")
    a = ap.parse_args()
    global ROTORS
    ROTORS = [("left", -X_ROTOR, 3), ("right", X_ROTOR, a.right_blades)]
    name = "twin_mixer.msh" if a.right_blades == 4 else f"twin_mixer{a.right_blades}.msh"
    path = os.path.join(a.out, name)
    print(f"{path}: {build(a.size, path)} cells")


if __name__ == "__main__":
    main()

"""A baffled stirred vessel in 2D (one cell layer): a rotor with four blades
(`--blades`) on a hub inside a circular interface, and a stator ring with
four baffles (`--baffles`) on the vessel wall. Meshes for the moving
reference frame and for the sliding mesh by imprinting (D55, D182).

    pip install gmsh
    python3 examples/meshes/mixer.py              # both meshes, default size
    python3 examples/meshes/mixer.py --size 1.2e-3
    python3 examples/meshes/mixer.py --blades 5 --baffles 3   # mixer5_*.msh

The five-blade vessel with three baffles is the one of the gallery images
(CONTRIBUTING.md, "Images": no four-fold rotor patterns in published
images); the four-blade meshes are the ones of the torque studies.

Writes, next to this script (`mixer5_...` for five blades),

- `mixer_conformal.msh`: rotor and stator share the nodes of the interface
  circle (interior faces between the cell zones `rotor` and `stator`), for
  `[[mrf]]`;
- `mixer_split.msh`: rotor and stator are two meshes that meet on the
  interface circle with their own nodes (64 rotor and 96 stator intervals
  around it), boundaries `rotor_interface` and `stator_interface`, for
  `[[sliding]] rotor_interface / stator_interface`.

Boundaries: `blades` (blades and hub, rotating with the rotor), `vessel`
(vessel wall and baffles), `zmin`, `zmax`. Lengths in m.
"""
import argparse
import math
import os

import gmsh

R_VESSEL = 0.10      # vessel radius
R_IFACE = 0.06       # interface radius
R_HUB = 0.012        # hub radius
R_TIP = 0.045        # blade tip radius
T_BLADE = 0.005      # blade thickness
L_BAFFLE = 0.012     # baffle length (radial)
T_BAFFLE = 0.005     # baffle thickness
DEPTH = 0.005        # thickness of the cell layer
N_ROTOR = 64         # intervals of the interface circle, rotor side
N_STATOR = 96        # stator side (split mesh)


def rect(x0, x1, half, angle):
    """A rectangle [x0, x1] x [-half, half] turned by `angle` about z."""
    s = gmsh.model.occ.addRectangle(x0, -half, 0.0, x1 - x0, 2.0 * half)
    gmsh.model.occ.rotate([(2, s)], 0, 0, 0, 0, 0, 1, angle)
    return (2, s)


def rotor(blades):
    occ = gmsh.model.occ
    disk = (2, occ.addDisk(0, 0, 0, R_IFACE, R_IFACE))
    tools = [(2, occ.addDisk(0, 0, 0, R_HUB, R_HUB))]
    for k in range(blades):
        tools.append(rect(0.5 * R_HUB, R_TIP, 0.5 * T_BLADE, 2.0 * math.pi * k / blades))
    out, _ = occ.cut([disk], tools)
    return out


def stator(hole, baffles):
    occ = gmsh.model.occ
    ring = (2, occ.addDisk(0, 0, 0, R_VESSEL, R_VESSEL))
    tools = [(2, occ.addDisk(0, 0, 0, hole, hole))]
    for k in range(baffles):
        # baffles through the wall (between the blades' angles for four)
        tools.append(rect(R_VESSEL - L_BAFFLE, R_VESSEL + 0.01, 0.5 * T_BAFFLE,
                          2.0 * math.pi * (k + 0.5) / baffles))
    out, _ = occ.cut([ring], tools)
    return out


def on_circle(dim, tag, r):
    """True for a curve or surface that lies on the circle (cylinder) r."""
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim, tag)
    return abs(xmax - r) < 1e-6 and abs(xmin + r) < 1e-6 and abs(ymax - r) < 1e-6


def build(split, size, path, blades=4, baffles=4):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("mixer")
    occ = gmsh.model.occ
    rot = rotor(blades)
    sta = stator(R_IFACE, baffles)
    if not split:
        occ.fragment(rot, sta)
    occ.synchronize()
    surfaces = gmsh.model.getEntities(2)
    rotor_surfs = [s for s in surfaces if gmsh.model.getBoundingBox(*s)[3] < R_IFACE + 1e-6]
    stator_surfs = [s for s in surfaces if s not in rotor_surfs]
    # interface curves: a fixed number of intervals on each side
    for s, n in [(rotor_surfs, N_ROTOR), (stator_surfs, N_STATOR if split else N_ROTOR)]:
        for c in gmsh.model.getBoundary(s, oriented=False):
            if on_circle(1, c[1], R_IFACE):
                gmsh.model.mesh.setTransfiniteCurve(c[1], n + 1)
    gmsh.option.setNumber("Mesh.MeshSizeMax", size)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.25 * size)
    # finer near the blades and baffles
    walls = [c[1] for c in gmsh.model.getBoundary(surfaces, oriented=False)
             if not on_circle(1, c[1], R_IFACE)]
    f = gmsh.model.mesh.field
    d = f.add("Distance")
    f.setNumbers(d, "CurvesList", walls)
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
    # (one call for all surfaces, so that the conformal mesh shares the
    # lateral faces of the interface)
    occ.extrude(rotor_surfs + stator_surfs, 0, 0, DEPTH, numElements=[1], recombine=True)
    if not split:
        occ.removeAllDuplicates()
    occ.synchronize()
    vols = {"rotor": [], "stator": []}
    for e in gmsh.model.getEntities(3):
        inside = gmsh.model.getBoundingBox(*e)[3] < R_IFACE + 1e-6
        vols["rotor" if inside else "stator"].append(e[1])
    groups = {"blades": [], "vessel": [], "zmin": [], "zmax": [],
              "rotor_interface": [], "stator_interface": [], }
    seen = {}
    for name, vs in vols.items():
        for v in vs:
            for b in gmsh.model.getBoundary([(3, v)], oriented=False):
                tag = b[1]
                seen.setdefault(tag, []).append(name)
    for tag, owners in seen.items():
        bb = gmsh.model.getBoundingBox(2, tag)
        if bb[5] < 1e-6:
            groups["zmin"].append(tag)
        elif bb[2] > DEPTH - 1e-6:
            groups["zmax"].append(tag)
        elif len(owners) == 2:
            continue                                    # conformal interface
        elif on_circle(2, tag, R_IFACE):
            groups[owners[0] + "_interface"].append(tag)
        elif owners[0] == "rotor":
            groups["blades"].append(tag)
        else:
            groups["vessel"].append(tag)
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
    ap.add_argument("--blades", type=int, default=4, help="number of blades")
    ap.add_argument("--baffles", type=int, default=4, help="number of baffles")
    a = ap.parse_args()
    stem = "mixer" if a.blades == 4 else f"mixer{a.blades}"
    for split, name in [(False, f"{stem}_conformal.msh"), (True, f"{stem}_split.msh")]:
        path = os.path.join(a.out, name)
        n = build(split, a.size, path, a.blades, a.baffles)
        print(f"{path}: {n} cells")


if __name__ == "__main__":
    main()

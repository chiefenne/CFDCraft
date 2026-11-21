"""An industrial baffled stirred tank in 3D: a dished (2:1 elliptical) bottom,
a centre shaft with a three-blade hydrofoil impeller (twisted, tapered
blades pumping down), three baffles with a gap to the wall, a sparger ring
on three legs below the impeller and a feed pipe dipping in from the top.
The rotor zone is a cylinder around the impeller whose faces the mesh
shares with the stator (conformal interface), for the moving reference
frame (`[[mrf]]`, D55).

    pip install gmsh
    python3 examples/applications/stirred_tank/mesh.py                 # level 1 (coarse)
    python3 examples/applications/stirred_tank/mesh.py --level 2       # medium
    python3 examples/applications/stirred_tank/mesh.py --level 3       # fine (iMac)
    python3 examples/applications/stirred_tank/mesh.py --level 0       # CI smoke mesh

Each level divides every cell size by 1.3 (about 2.2 times the cells), the
ratio the grid convergence index needs (D67). Writes `tank_L<level>.msh`
next to this script. Tetrahedra, optimised with Netgen; no prism layers
(the case uses wall functions).

Boundaries: `impeller` (hub, blades and the shaft inside the rotor zone),
`shaft` (the shaft above it), `vessel` (shell and dished bottom), `baffles`,
`sparger` (ring and legs), `feed_pipe`, `surface` (the liquid level).
Cell zones: `rotor`, `stator`. Lengths in m.

The tank (T = 1 m) is a pilot-to-production size: liquid height H = T,
impeller diameter D = 0.4 T, off-bottom clearance C = T/3 at the impeller
mid-plane. Three baffles and three blades: no four-fold patterns in the
images (CONTRIBUTING.md).
"""
import argparse
import math
import os

import gmsh

T = 1.0                  # tank diameter
R = 0.5 * T              # tank radius
H = 1.0 * T              # liquid height above the bottom apex
HEAD = 0.25 * T          # depth of the 2:1 elliptical head
D = 0.4 * T              # impeller diameter
Z_IMP = T / 3.0          # impeller mid-plane above the bottom apex
R_HUB, L_HUB = 0.045, 0.12  # the hub reaches well past the blade roots (edges away from the rims)
R_FILLET = 0.010         # rounded hub rims
R_SHAFT = 0.025
N_BLADES = 3
N_SECTION = 48           # points of the spline of one blade section
H_JUNCTION = 0.35        # cell size at the blade-hub corners, times the blade size
R_ZONE, Z0_ZONE, Z1_ZONE = 0.25, 0.245, 0.425   # rotor zone cylinder
R_COLUMN = 0.06          # the rotor zone goes on around the shaft up to the liquid level
N_BAFFLES = 3
W_BAFFLE = T / 10.0      # baffle width
GAP_BAFFLE = 0.015       # gap between baffle and wall
T_BAFFLE = 0.010         # baffle thickness
Z_BAFFLE = 0.27          # baffles start just above the head
R_RING, R_TUBE, Z_RING = 0.13, 0.012, 0.17      # sparger torus
R_LEG = 0.008
R_FEED, RAD_FEED, Z_FEED = 0.02, 0.36, 0.55     # feed pipe: radius, position, lower end
A_FEED = 2.0 * math.pi / N_BAFFLES               # between the first two baffles

# blade sections: radius, chord, thickness, pitch angle from the horizontal (deg)
SECTIONS = [
    (0.035, 0.085, 0.014, 46.0),
    (0.080, 0.110, 0.012, 38.0),
    (0.130, 0.110, 0.010, 30.0),
    (0.175, 0.095, 0.009, 25.0),
    (0.200, 0.070, 0.008, 22.0),
]

LEVEL_SIZES = {  # bulk, rotor zone, blades, baffles and sparger, feed and shaft
    0: (0.10, 0.035, 0.016, 0.030, 0.04),
}


def sizes(level):
    """Cell sizes of a level: level 1 is the coarse mesh, each level / 1.3."""
    if level in LEVEL_SIZES:
        return LEVEL_SIZES[level]
    base = (0.060, 0.020, 0.009, 0.016, 0.020)
    f = 1.3 ** (level - 1)
    return tuple(s / f for s in base)


def blade(angle):
    """One twisted, tapered blade along +x through round-edged sections, turned
    by `angle` about z. Leading edge up and forward (+y): pumps down for a
    positive rotation about +z."""
    occ = gmsh.model.occ
    wires = []
    for r, chord, thick, pitch in SECTIONS:
        # a flat blade of constant thickness with round edges (radius t/2),
        # as one closed spline per section: a section of lines and arcs
        # leaves seam curves on the blade that meet the hub in short edges,
        # and a thin ellipse meets it in a cusp; both give skewed cells there
        p = math.radians(pitch)
        c = (0.0, math.cos(p), math.sin(p))               # chord direction
        n = (0.0, -math.sin(p), math.cos(p))              # thickness direction
        a, h = 0.5 * chord - 0.5 * thick, 0.5 * thick
        pts = []
        for i in range(N_SECTION):
            # arc length parameter around the stadium, from the trailing edge
            s_ = i / N_SECTION * (4.0 * a + 2.0 * math.pi * h)
            if s_ < 2.0 * a:                               # lower flat side
                u, v = -a + s_, -h
            elif s_ < 2.0 * a + math.pi * h:               # leading edge
                phi = (s_ - 2.0 * a) / h - 0.5 * math.pi
                u, v = a + h * math.cos(phi), h * math.sin(phi)
            elif s_ < 4.0 * a + math.pi * h:               # upper flat side
                u, v = a - (s_ - 2.0 * a - math.pi * h), h
            else:                                          # trailing edge
                phi = (s_ - 4.0 * a - math.pi * h) / h + 0.5 * math.pi
                u, v = -a + h * math.cos(phi), h * math.sin(phi)
            pts.append(occ.addPoint(r, u * c[1] + v * n[1], Z_IMP + u * c[2] + v * n[2]))
        curve = occ.addSpline(pts + [pts[0]])
        wires.append(occ.addWire([curve]))
    out = occ.addThruSections(wires, makeSolid=True, makeRuled=False)
    vols = [e for e in out if e[0] == 3]
    occ.rotate(vols, 0, 0, 0, 0, 0, 1, angle)
    return vols


def geometry():
    occ = gmsh.model.occ
    # tank: cylinder above the head, ellipsoid for the head
    shell = occ.addCylinder(0, 0, HEAD, 0, 0, H - HEAD, R)
    head = occ.addSphere(0, 0, HEAD, R)
    occ.dilate([(3, head)], 0, 0, HEAD, 1.0, 1.0, HEAD / R)
    cut_box = occ.addBox(-R, -R, HEAD, 2 * R, 2 * R, R)
    head_low, _ = occ.cut([(3, head)], [(3, cut_box)])
    tank, _ = occ.fuse([(3, shell)], head_low)

    # impeller: hub, blades, shaft up through the liquid level
    hub = occ.addCylinder(0, 0, Z_IMP - 0.5 * L_HUB, 0, 0, L_HUB, R_HUB)
    # rounded rims: sharp hub edges give sliver cells (skewness above 1)
    occ.synchronize()
    faces = gmsh.model.getBoundary([(3, hub)], oriented=False)
    edges = sorted(set(abs(c) for f in faces
                       for _, c in gmsh.model.getBoundary([f], oriented=False)
                       if gmsh.model.getType(1, abs(c)) == "Circle"))
    hub = occ.fillet([hub], edges, [R_FILLET])[0][1]
    shaft = occ.addCylinder(0, 0, Z_IMP, 0, 0, H - Z_IMP + 0.05, R_SHAFT)
    blades = []
    for k in range(N_BLADES):
        blades += blade(2.0 * math.pi * k / N_BLADES)
    impeller, _ = occ.fuse([(3, hub)], blades + [(3, shaft)])

    # baffles, offset by half a pitch from the blades' start angle
    obstacles = []
    for k in range(N_BAFFLES):
        a = 2.0 * math.pi * (k + 0.5) / N_BAFFLES
        b = occ.addBox(R - GAP_BAFFLE - W_BAFFLE, -0.5 * T_BAFFLE, Z_BAFFLE,
                       W_BAFFLE, T_BAFFLE, H - Z_BAFFLE + 0.05)
        occ.rotate([(3, b)], 0, 0, 0, 0, 0, 1, a)
        obstacles.append((3, b))
    # sparger ring on three legs, legs between the baffles
    ring = occ.addTorus(0, 0, Z_RING, R_RING, R_TUBE)
    legs = []
    for k in range(3):
        a = 2.0 * math.pi * k / 3 + math.pi / 6
        legs.append((3, occ.addCylinder(R_RING * math.cos(a), R_RING * math.sin(a), -0.05,
                                        0, 0, Z_RING + 0.05, R_LEG)))
    sparger, _ = occ.fuse([(3, ring)], legs)
    obstacles += sparger
    # feed pipe from above, between two baffles
    feed = occ.addCylinder(RAD_FEED * math.cos(A_FEED), RAD_FEED * math.sin(A_FEED), Z_FEED,
                           0, 0, H - Z_FEED + 0.05, R_FEED)
    obstacles.append((3, feed))

    fluid, _ = occ.cut(tank, impeller + obstacles)
    # the rotor zone: a cylinder made of three 120° sectors (cut planes
    # between the blades), so that each lateral face has four sides and
    # takes a structured mesh of vertical quadrilaterals (see interface())
    # The zone goes on as a column around the shaft up to the liquid level,
    # so that the whole shaft turns with the frame (a wall at rest in it):
    # a turning curved wall in the fixed frame destabilises the steady run.
    # The column has no blades, its flow is axisymmetric and the frozen
    # rotor is exact there.
    corners = [(0, Z0_ZONE), (R_ZONE, Z0_ZONE), (R_ZONE, Z1_ZONE), (R_COLUMN, Z1_ZONE),
               (R_COLUMN, H), (0, H)]
    p = [occ.addPoint(r, 0, z) for r, z in corners]
    loop = occ.addCurveLoop([occ.addLine(p[i], p[(i + 1) % len(p)]) for i in range(len(p))])
    face = occ.addPlaneSurface([loop])
    occ.rotate([(2, face)], 0, 0, 0, 0, 0, 1, math.pi / N_BLADES)
    sectors = []
    for k in range(N_BLADES):
        f = occ.copy([(2, face)])
        occ.rotate(f, 0, 0, 0, 0, 0, 1, 2.0 * math.pi * k / N_BLADES)
        sectors += [e for e in occ.revolve(f, 0, 0, 0, 0, 0, 1, 2.0 * math.pi / N_BLADES)
                    if e[0] == 3]
    occ.remove([(2, face)], recursive=True)
    out, children = occ.fragment(fluid, sectors)
    occ.synchronize()
    fluid_children = set(t for _, t in children[0])
    zone_children = set(t for c in children[1:] for _, t in c)
    rotor = sorted(fluid_children & zone_children)
    stator = sorted(fluid_children - zone_children)
    junk = sorted(zone_children - fluid_children)  # impeller inside the zone
    if junk:
        occ.remove([(3, t) for t in junk], recursive=True)
        occ.synchronize()
    return rotor, stator


def classify(rotor, stator):
    """Physical groups of the boundary surfaces, by owner and position."""
    owners, count = {}, {}
    for name, vs in (("rotor", rotor), ("stator", stator)):
        for v in vs:
            for _, s in gmsh.model.getBoundary([(3, v)], oriented=False):
                owners.setdefault(s, set()).add(name)
                count[s] = count.get(s, 0) + 1
    groups = {n: [] for n in ("impeller", "shaft", "vessel", "baffles", "sparger",
                              "feed_pipe", "surface")}
    for s, own in owners.items():
        if count[s] == 2:
            continue                  # conformal interface, or between two rotor sectors
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, s)
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        rc, span = math.hypot(cx, cy), max(x1 - x0, y1 - y0)
        if z0 > H - 1e-6:
            groups["surface"].append(s)
        elif max(abs(x0), abs(x1), abs(y0), abs(y1)) < 1.01 * R_SHAFT:
            groups["shaft"].append(s)                 # in pieces, cut by the sector planes
        elif own == {"rotor"}:
            groups["impeller"].append(s)
        elif span > 0.9 * T:
            groups["vessel"].append(s)
        elif z1 < Z_RING + 2 * R_TUBE:
            groups["sparger"].append(s)
        elif (abs(rc - RAD_FEED) < 2 * R_FEED and span < 3 * R_FEED
              and abs(math.atan2(cy, cx) - A_FEED) < 0.2):
            groups["feed_pipe"].append(s)
        elif rc > R - GAP_BAFFLE - W_BAFFLE - 0.02:
            groups["baffles"].append(s)
        else:
            raise SystemExit(f"surface {s} not classified: bbox {x0, y0, z0, x1, y1, z1}")
    for name, tags in groups.items():
        if not tags:
            raise SystemExit(f"no surfaces for {name}")
        gmsh.model.addPhysicalGroup(2, sorted(tags), name=name)
    gmsh.model.addPhysicalGroup(3, rotor, name="rotor")
    gmsh.model.addPhysicalGroup(3, stator, name="stator")
    return groups


def interface(rotor, stator, level):
    """A structured mesh on the lateral interface of the rotor zone. The
    frame flux of a face is ρ (Ω × x_f) · S_f at its centroid; on the
    interface it should vanish, or the relative flux of the cells next to it
    does not close (MRF003). A vertical rectangle on a chord has zero flux,
    a triangle on the cylinder does not. Gmsh cannot put pyramids on quads
    of a surface between two volumes, so each rectangle is split into two
    triangles: their fluxes are equal and opposite (the rectangle's is
    zero), and the closure error of a cell is about Δθ/6 of its flux for
    an interval Δθ around the axis, which falls with the mesh (2.6 % on
    level 0, under 1 % on level 3)."""
    h_zone, h_feed = sizes(level)[1], sizes(level)[4]
    rot = set(s for v in rotor for _, s in gmsh.model.getBoundary([(3, v)], oriented=False))
    sta = set(s for v in stator for _, s in gmsh.model.getBoundary([(3, v)], oriented=False))
    lateral = []
    for s in sorted(rot & sta):
        z0, z1 = gmsh.model.getBoundingBox(2, s)[2], gmsh.model.getBoundingBox(2, s)[5]
        if z1 - z0 < 1e-6:
            continue                                  # top or bottom disc: any faces will do
        column = z0 > Z1_ZONE - 1e-6                  # the lateral face of the shaft column
        radius, h = (R_COLUMN, h_feed) if column else (R_ZONE, h_zone)
        n_theta = max(4, math.ceil(2.0 * math.pi * radius / N_BLADES / h))
        n_z = max(2, math.ceil((z1 - z0) / h))
        for _, c in gmsh.model.getBoundary([(2, s)], oriented=False):
            bb = gmsh.model.getBoundingBox(1, c)
            vertical = bb[5] - bb[2] > 1e-6
            gmsh.model.mesh.setTransfiniteCurve(c, (n_z if vertical else n_theta) + 1)
        gmsh.model.mesh.setTransfiniteSurface(s, "Alternate")
        lateral.append(s)
    if len(lateral) != 2 * N_BLADES:
        raise SystemExit(f"expected {2 * N_BLADES} lateral interface faces, found {lateral}")


def junction_curves(impeller):
    """The curves where a blade surface meets the hub cylinder."""
    hub = {s for s in impeller if gmsh.model.getType(2, s) == "Cylinder"}
    curves = set()
    for s in impeller:
        if s in hub:
            continue
        for _, c in gmsh.model.getBoundary([(2, s)], oriented=False):
            up, _ = gmsh.model.getAdjacencies(1, abs(c))
            if any(u in hub for u in up):
                curves.add(abs(c))
    if not curves:
        raise SystemExit("no curves found where the blades meet the hub")
    return sorted(curves)


def size_fields(groups, level):
    h_bulk, h_zone, h_blade, h_obst, h_feed = sizes(level)
    f = gmsh.model.mesh.field
    fields = []

    def near(surfs, h, dist_min, dist_max):
        d = f.add("Distance")
        f.setNumbers(d, "SurfacesList", surfs)
        f.setNumber(d, "Sampling", 40)
        t = f.add("Threshold")
        f.setNumber(t, "InField", d)
        f.setNumber(t, "SizeMin", h)
        f.setNumber(t, "SizeMax", h_bulk)
        f.setNumber(t, "DistMin", dist_min)
        f.setNumber(t, "DistMax", dist_max)
        fields.append(t)

    near(groups["impeller"], h_blade, 0.01, 0.12)
    # the corners where the blades meet the hub: cells of the blade size
    # there have faces skewed above 1 (D50), where the pressure blew up
    d = f.add("Distance")
    f.setNumbers(d, "CurvesList", junction_curves(groups["impeller"]))
    f.setNumber(d, "Sampling", 60)
    t = f.add("Threshold")
    f.setNumber(t, "InField", d)
    f.setNumber(t, "SizeMin", H_JUNCTION * h_blade)
    f.setNumber(t, "SizeMax", h_bulk)
    f.setNumber(t, "DistMin", 0.004)
    f.setNumber(t, "DistMax", 0.02)
    fields.append(t)
    near(groups["baffles"] + groups["sparger"], h_obst, 0.01, 0.12)
    near(groups["feed_pipe"] + groups["shaft"], h_feed, 0.01, 0.10)
    near(groups["vessel"], 1.4 * h_obst, 0.01, 0.10)
    c = f.add("Cylinder")
    f.setNumber(c, "Radius", R_ZONE + 0.03)
    f.setNumber(c, "XCenter", 0.0)
    f.setNumber(c, "YCenter", 0.0)
    f.setNumber(c, "ZCenter", 0.5 * (Z0_ZONE + Z1_ZONE))
    f.setNumber(c, "XAxis", 0.0)
    f.setNumber(c, "YAxis", 0.0)
    f.setNumber(c, "ZAxis", 0.5 * (Z1_ZONE - Z0_ZONE) + 0.03)
    f.setNumber(c, "VIn", h_zone)
    f.setNumber(c, "VOut", h_bulk)
    fields.append(c)
    m = f.add("Min")
    f.setNumbers(m, "FieldsList", fields)
    f.setAsBackgroundMesh(m)
    # sizes graded from the walls: the blade-hub and sparger corners keep cells
    # with skewness below 1.1 (with 0 they reach 1.4 and the solver blows up)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 10)   # small radii: legs, pipes, edges
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_bulk)
    gmsh.option.setNumber("Mesh.MeshSizeMin", min(0.5, H_JUNCTION) * h_blade)


def build(level, path, threads, verbose=False):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.Verbosity", 4 if verbose else 2)
    gmsh.option.setNumber("General.NumThreads", threads)
    gmsh.model.add("stirred_tank")
    rotor, stator = geometry()
    groups = classify(rotor, stator)
    size_fields(groups, level)
    interface(rotor, stator, level)
    for s in groups["impeller"]:                      # lofted blade surfaces
        gmsh.model.mesh.setAlgorithm(2, s, 1)         # MeshAdapt: robust on B-splines
    gmsh.option.setNumber("Mesh.Algorithm", 6)        # Frontal-Delaunay surfaces
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)      # Delaunay: robust boundary recovery
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(path)
    n = 0
    for v in rotor + stator:
        types, tags, _ = gmsh.model.mesh.getElements(3, v)
        n += sum(len(t) for t in tags)
    gmsh.finalize()
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--level", type=int, default=1, help="0 smoke, 1 coarse, 2 medium, 3 fine")
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args()
    path = os.path.join(a.out, f"tank_L{a.level}.msh")
    n = build(a.level, path, a.threads)
    print(f"{path}: {n} cells")


if __name__ == "__main__":
    main()

# Meshes

Part of the [Recipes](../recipes.md): generate, read (Gmsh, SU2, CGNS, OpenFOAM), convert and partition meshes, update the mesher's golden values, run flows on poor Gmsh meshes.

| # | Recipe |
|---|---|
| 7 | [Generate a test mesh and check its quality](#mesh-generate) |
| 8 | [Read a Gmsh mesh and convert it into polyhedra](#mesh-gmsh) |
| 9 | [Partition a mesh and test the halo exchange](#mesh-partition) |
| 11 | [Update the golden values of the mesher](#golden) |
| 36 | [Run a flow on raw Gmsh tetrahedra (skewed faces)](#skewed-tets) |
| 127 | [Read an SU2 mesh and map its markers to boundary conditions](#mesh-su2) |
| 129 | [Read a CGNS or OpenFOAM mesh](#mesh-cgns-openfoam) |

---

<a id="mesh-generate"></a>
## 7. Generate a test mesh and check its quality

The built-in mesher (D34) makes a box of hexahedra, tetrahedra, polyhedra
or mixed cells, in 3D or 2D (one cell thick), graded, refined and
perturbed. Keys `[mesh.generate]` ([Case catalogue](../cases/meshes.md#mesh-case-file-keys)).

```
cargo build --release -p foehncraft-cli
target/release/foehn mesh examples/cases/mesh/poly_box.toml
cat output/poly_box/mesh_report.md
```

The report gives the verdict (valid or not), the cell types, the closure
|Σ S| / Σ |S| of the worst cell, the total volume against the volume
enclosed by the boundary, the boundary areas, non-orthogonality, skewness,
face warp, aspect and volume ratio with the cell or face where the maximum
sits, and the point-location test. `mesh.vtkhdf` holds the mesh with these
fields per cell; open it in ParaView (5.13 or newer reads VTKHDF 2.4 with
polyhedra) or render a picture:

```
pip install vtk
python3 -m tools.mesh.render_mesh output/poly_box/mesh.vtkhdf poly.png --cut --field non_orthogonality
```

A grid-convergence family: the same case with `refine = 0, 1, 2` gives
nested meshes (every interval halved), the grading kept.

[↑ back to top](#meshes)

---

<a id="mesh-gmsh"></a>
## 8. Read a Gmsh mesh and convert it into polyhedra

Needs Gmsh 4.x (`pip install gmsh` or gmsh.info). The reader takes MSH 4.1,
ASCII or binary, first-order elements (tetrahedra, hexahedra, prisms,
pyramids), physical groups and periodic links (D33).

```
gmsh -3 examples/meshes/tjunction.geo -format msh41 -o examples/meshes/tjunction.msh
target/release/foehn mesh examples/meshes/tjunction.msh             # defaults, no case file
target/release/foehn mesh examples/cases/mesh/tjunction_tet.toml    # RCM, 4 parts
target/release/foehn mesh examples/cases/mesh/tjunction_poly.toml   # converted into polyhedra
```

In your own case file: `mesh.file` (relative to the case file),
`mesh.scale = 0.001` for a mesh in millimetres, `mesh.polyhedral = true`
for the conversion of a tetrahedral mesh into polyhedra (one cell per node,
same volume), `mesh.renumber = "rcm"` for a band-limited cell order, and
`[mesh.boundaries]` to rename physical groups. Give every boundary surface a
physical group in Gmsh; surfaces without one become `surface_<tag>`.

[↑ back to top](#meshes)

---

<a id="mesh-partition"></a>
## 9. Partition a mesh and test the halo exchange

```
target/release/foehn mesh examples/cases/mesh/poly_partition.toml
target/release/foehn mesh examples/cases/mesh/poly_partition.toml --parts 8 --method rcb
```

Keys `partition.parts`, `partition.method` (`metis` or `rcb`),
`partition.halo_layers`. The report lists cells, halo cells, processor
faces and neighbours per rank, the imbalance and the cut faces;
`mesh.vtkhdf` has the field `rank`.

The halo exchange over MPI (needs MPI as in recipe 4): every rank builds
its part with a two-layer halo and fills it from its neighbours.

```
cargo build --release -p foehncraft-mesh --features mpi --example mpi_halo
mpirun --oversubscribe -n 4 target/release/examples/mpi_halo
```

Prints `MPI halo exchange: PASS` when the halo values, volumes and centroids
equal the global mesh and the volume sum gives the same bits on every rank.

[↑ back to top](#meshes)

---

<a id="golden"></a>
## 11. Update the golden values of the mesher

`crates/foehncraft-mesh/tests/golden/mesher.json` holds the cell, face and
node counts, volumes, quality metrics and partition statistics of every
combination of mesh type, dimension, perturbation, refinement and 1, 2 and
4 ranks (R4). A difference fails the test. Only for an intended change:

```
FOEHNCRAFT_BLESS=1 cargo test -p foehncraft-mesh --test golden
git diff crates/foehncraft-mesh/tests/golden/mesher.json
```

and say in the commit message why the values changed.

[↑ back to top](#meshes)

---

<a id="skewed-tets"></a>
## 36. Run a flow on raw Gmsh tetrahedra (skewed faces)

Gmsh's tetrahedra without its optimisation pass have faces of a skewness
above 0.9, mostly at curved walls; the steady solver diverges on them with
the least-squares pressure gradient ([gallery 26](../gallery/meshes.md#skewed-tets)).
First choice: optimise the mesh. Add `Mesh.OptimizeNetgen = 1;` to the .geo
file (or `-optimize_netgen` on the command line) and check the count in the
quality report:

```
gmsh -3 examples/meshes/tjunction_netgen.geo -format msh41 -o examples/meshes/tjunction_netgen.msh
target/release/foehn mesh examples/meshes/tjunction_netgen.msh
grep "skewness above 0.9" output/tjunction_netgen/mesh_report.md
```

When the mesh has such faces, the flow solver treats them itself
(`solver.remediation = "auto"`, the default): it switches to the Gauss
pressure gradient and the Venkatakrishnan limiter for the whole run and
says so at the start of `flow.log` and in `flow.md`, `flow.json`
(`remediation` of each solve) and the diagnostics (`MESH001`, `MESH002`):

```
target/release/foehn flow examples/cases/flow/pipe_gmsh_raw.toml
grep -E "mesh quality|remediation" output/pipe_gmsh_raw/flow.log
```

Keys (`[solver]`): `remediation = "auto"` or `"off"` (keep the schemes as
set; the run on `pipe_gmsh_raw` then diverges and `divergence.md` names the
skewed faces), `skewness_limit = 0.9` (a smaller value treats milder
faces), `pressure_gradient = "least-squares"` or `"gauss"` (set by hand,
also on a good mesh). A limiter set in the case stays, and first-order
`upwind` needs none. Least squares is more accurate on good meshes (the
cylinder lift at Re = 20 is 3.5 % lower with the Gauss gradient), so treat
the remediation as a remedy and fix the mesh for production runs.

```
target/release/foehn flow examples/cases/flow/pipe_gmsh_raw.toml --set solver.remediation=off \
    --out output/pipe_gmsh_raw_off
cat output/pipe_gmsh_raw_off/divergence.md
```

[↑ back to top](#meshes)

---

<a id="mesh-su2"></a>
## 127. Read an SU2 mesh and map its markers to boundary conditions

The SU2 reader (D187) takes native `.su2` files, ASCII, one zone, 2D or
3D, first-order elements: tetrahedra, hexahedra, prisms and pyramids in 3D
with triangle and quadrilateral markers; triangles and quadrilaterals in 2D
with line markers. A 2D file is extruded one cell layer in z (triangles to
prisms, quadrilaterals to hexahedra), and the planes become the boundaries
`zmin` and `zmax`, as for the built-in 2D meshes. Every `MARKER_TAG` is a
boundary of that name.

```
target/release/foehn mesh examples/meshes/naca0012.su2               # defaults, no case file
target/release/foehn mesh examples/cases/mesh/naca0012_su2.toml      # thickness, RCM, 4 parts
target/release/foehn flow examples/cases/flow/naca0012_su2_re1000.toml
```

(`examples/meshes/naca0012.su2` is a PyAero mesh in the repository, see
[gallery 127](../gallery/meshes.md#su2-airfoil).) In your own case file:

```toml
[mesh]
file = "wing.su2"          # the format follows from the extension, or:
format = "su2"             # gmsh | su2 | cgns | openfoam, for a file with another extension
thickness = 0.01           # 2D files: layer thickness, in the file's units
scale = 0.001              # millimetres to metres (after the extrusion)

[mesh.boundaries]          # marker names to the names of the conditions
farfield = "inlet"
wing = "wall"

[boundary.wall]
type = "wall"
```

Without `thickness` the layer is as thick as the median edge of the 2D
cells. Boundary faces that no marker names go to `default` with a warning.
A marker called `zmin` or `zmax` in a 2D file merges with the extruded
plane of that name: rename it with `[mesh.boundaries]`. Multi-zone files
(`NZONE` > 1) are refused: write one file per zone. The format layer
picks the reader from `mesh.format`, else the extension (`.msh`, `.su2`),
else the first lines of the file; a new format is one reader in
`crates/foehncraft-mesh/src/format/` (module docs of `format`).

[↑ back to top](#meshes)

---

<a id="mesh-cgns-openfoam"></a>
## 129. Read a CGNS or OpenFOAM mesh

Two more readers sit behind the same interface (D188).

- **OpenFOAM.** Give the case directory, its `constant` directory or the
  `polyMesh` directory itself: `file = "channel"`. The files `points`,
  `faces`, `owner`, `neighbour` and, if present, `cellZones` are read, ASCII
  or binary (32 or 64 bit labels, either byte order), with polyhedral cells
  kept as they are. Patches become boundaries by name, `cellZones` become
  cell zones. `empty` patches stay ordinary boundaries (the reader says so;
  use a `symmetry` boundary for 2D work). Decomposed cases (`processor`
  patches) and `.gz` files are refused: reconstruct and `gunzip` first.
- **CGNS.** HDF5 files (CGNS 3 or 4 as written by Gmsh and others), one
  base, one unstructured zone, 3D cells TETRA_4, PYRA_5, PENTA_6, HEXA_8
  and polyhedra as NGON_n/NFACE_n. The BCs of `ZoneBC`, given on faces,
  become boundaries; for Gmsh files the physical-group name is taken from
  the family. Old ADF files are refused: convert them with `adf2hdf`.

```
target/release/foehn mesh crates/foehncraft-mesh/tests/data/openfoam/dual_polyhedra
target/release/foehn mesh crates/foehncraft-mesh/tests/data/cgns/box_tet.cgns
```

```toml
[mesh]
file = "channel"            # an OpenFOAM case, or "wing.cgns"

[mesh.boundaries]           # patch or BC names to the names of the conditions
walls = "wall"
```

Faces that no patch or BC names go to `default` with a warning. Multi-zone
and 2D CGNS files and structured zones are refused.

[↑ back to top](#meshes)

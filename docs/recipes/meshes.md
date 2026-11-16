# Meshes

Part of the [Recipes](../recipes.md): generate, read, convert and partition meshes, update the mesher's golden values.

| # | Recipe |
|---|---|
| 7 | [Generate a test mesh and check its quality](#mesh-generate) |
| 8 | [Read a Gmsh mesh and convert it into polyhedra](#mesh-gmsh) |
| 9 | [Partition a mesh and test the halo exchange](#mesh-partition) |
| 11 | [Update the golden values of the mesher](#golden) |

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

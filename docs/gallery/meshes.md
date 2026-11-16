# Meshes

Part of the [Example gallery](../gallery.md): the built-in test mesher in 2D and 3D, partitioning and the Gmsh reader.

| # | Example | Stage | Feature |
|---|---|---|---|
| 3 | [The four test-mesh types in 3D](#mesh-types-3d) | 1 | built-in mesher, grading, perturbation, quality report |
| 4 | [A polyhedral mesh in four parts](#poly-partition) | 1 | partitioning, halo layout, point location |
| 5 | [The four test-mesh types in 2D](#mesh-types-2d) | 1 | one-cell-thick 2D meshes, boundary renaming |
| 6 | [T-junction from Gmsh](#tjunction) | 1 | Gmsh reader, conversion to polyhedra, renumbering |

---

<a id="mesh-types-3d"></a>
## 3. The four test-mesh types in 3D

![Hexahedral, tetrahedral, polyhedral and mixed test meshes of a box, cut open](../figures/gallery/mesh_types_3d.png)

The built-in test mesher (D34) makes four mesh types of the same box
2 × 1 × 1: hexahedra, tetrahedra (six per hexahedron), polyhedra (the dual
of the tetrahedral mesh, one cell per node) and a mixed mesh of hexahedra,
prisms and tetrahedra, with polyhedra where a hexahedron meets the
triangles of its neighbours. All four are clustered towards the walls
y = 0 and y = 1 (middle-to-wall ratio 4) and have perturbed nodes (15 % of
the local interval), so the geometry is never trivially orthogonal. Every
run writes the quality report (`mesh_report.md` and `.json`): closure of
every cell, total volume against the enclosed volume, non-orthogonality,
skewness, face warp, aspect and volume ratio, inverted faces and negative
volumes. All four meshes have the exact volume 2 and closure below 1e-15.

| | |
|---|---|
| Meshes | hex 576, tet 3 456, poly 819, mixed 1 728 cells |
| Case files | `mesh/hex_box.toml`, `mesh/tet_box.toml`, `mesh/poly_box.toml`, `mesh/mixed_box.toml` |

```
cargo build --release -p foehncraft-cli
for t in hex tet poly mixed; do
  target/release/foehn mesh examples/cases/mesh/${t}_box.toml
done
# output/<case>/mesh_report.md, mesh_report.json, mesh.vtkhdf (open in ParaView)
python3 -m tools.mesh.render_mesh output/poly_box/mesh.vtkhdf poly_box.png --cut
```

[↑ back to top](#meshes)

---

<a id="poly-partition"></a>
## 4. A polyhedral mesh in four parts

![Polyhedral unit cube in four METIS parts](../figures/gallery/poly_partition.png)

A polyhedral mesh of the unit cube (2 197 cells, nodes perturbed by 20 %)
split into four parts with METIS (D32), one colour per rank. Each part is
laid out as solver ranks will see it: inner cells, then cells touching the
halo, then the halo cells by source rank, each group in reverse
Cuthill-McKee order; faces inner, processor, then boundary by boundary. The
report lists the cells, halo cells, processor faces and neighbours of each
rank (imbalance 1.001, 1 230 cut faces). The point-location test of D47
throws 5 000 random points and finds each one by walking from a random cell
through the faces; all 5 000 are found in the right cell.

| | |
|---|---|
| Mesh | poly, 12 × 12 × 12 base grid, 2 197 cells |
| Case file | `mesh/poly_partition.toml` |

```
target/release/foehn mesh examples/cases/mesh/poly_partition.toml
target/release/foehn mesh examples/cases/mesh/poly_partition.toml --method rcb --parts 6
python3 -m tools.mesh.render_mesh output/poly_partition/mesh.vtkhdf part.png --cut --field rank --categorical
```

[↑ back to top](#meshes)

---

<a id="mesh-types-2d"></a>
## 5. The four test-mesh types in 2D

![Hexahedral, tetrahedral, polyhedral and mixed 2D meshes](../figures/gallery/mesh_types_2d.png)

2D meshes are one cell thick in z (D3): quadrilaterals become hexahedra,
triangles become prisms, polygons become prismatic polyhedra. The mixed
mesh combines both. Clustering towards x = 0 and towards the walls, nodes
perturbed in x and y. The case files rename the generator's boundaries
(`xmin` to `inlet` and so on) with a `[mesh.boundaries]` table; `zmin` and
`zmax` take the symmetry condition in a 2D case.

| | |
|---|---|
| Meshes | hex 160, tet 320, poly 187, mixed 240 cells |
| Case files | `mesh/hex_2d.toml`, `mesh/tet_2d.toml`, `mesh/poly_2d.toml`, `mesh/mixed_2d.toml` |

```
for t in hex tet poly mixed; do
  target/release/foehn mesh examples/cases/mesh/${t}_2d.toml
done
python3 -m tools.mesh.render_mesh output/poly_2d/mesh.vtkhdf poly_2d.png --top
```

[↑ back to top](#meshes)

---

<a id="tjunction"></a>
## 6. T-junction from Gmsh

![T-junction, tetrahedral mesh from Gmsh, cut at the mid plane, four ranks](../figures/gallery/tjunction_tet.png)

![The same T-junction converted into polyhedra](../figures/gallery/tjunction_poly.png)

The T-junction mixer of stage 7 (D45) meshed with Gmsh
(`examples/meshes/tjunction.geo`, main pipe 100 mm diameter, branch 50 mm)
and read by the Gmsh v4 reader (D33). Physical groups become boundaries
(`inlet_main`, `inlet_branch`, `outlet`, `wall`) and the cell zone `fluid`.
The second case converts the 77 776 tetrahedra into 15 927 polyhedra, one
per node; on the curved wall the boundary faces follow the triangles of
the tetrahedral mesh, so both meshes have the same volume to the last
digit (8.310580753620e-3 m³). Both are renumbered by reverse
Cuthill-McKee and split into four parts. The point-location walk stops at
the curved, non-convex wall for some points (176 of 1 195); the bin search
finds those, none lands in a wrong cell.

| | |
|---|---|
| Meshes | tet 77 776 cells, poly 15 927 cells (7 to 23 faces, mean 15) |
| Case files | `mesh/tjunction_tet.toml`, `mesh/tjunction_poly.toml` |

```
gmsh -3 examples/meshes/tjunction.geo -format msh41 -o examples/meshes/tjunction.msh
target/release/foehn mesh examples/cases/mesh/tjunction_tet.toml
target/release/foehn mesh examples/cases/mesh/tjunction_poly.toml
python3 -m tools.mesh.render_mesh output/tjunction_poly/mesh.vtkhdf tj.png \
    --field rank --categorical --half --top --crop 0.2 0.8 -0.06 0.31
```

[↑ back to top](#meshes)

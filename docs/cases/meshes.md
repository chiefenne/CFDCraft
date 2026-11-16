# Meshes (`examples/cases/mesh/`)

Part of the [Case catalogue](../cases.md): build or read a mesh, check, renumber and partition it; mesh keys.

On this page:

- [Mesh case-file keys](#mesh-case-file-keys)

Stage 1: build or read a mesh, check it, renumber and partition it, write
the quality report (`mesh_report.md`, `mesh_report.json`) and `mesh.vtkhdf`
for ParaView. Run with `foehn mesh FILE`; a Gmsh file can be given
directly instead of a case file (`foehn mesh pipe.msh`), command-line
options override the file (`foehn mesh --help`).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `mesh/hex_box.toml` | generated, 576 hexahedra, graded, perturbed | hexahedral test mesh, quality report | [gallery 3](../gallery/meshes.md#mesh-types-3d) |
| `mesh/tet_box.toml` | generated, 3 456 tetrahedra | tetrahedral test mesh | [gallery 3](../gallery/meshes.md#mesh-types-3d) |
| `mesh/poly_box.toml` | generated, 819 polyhedra | polyhedral test mesh (dual of tetrahedra) | [gallery 3](../gallery/meshes.md#mesh-types-3d) |
| `mesh/mixed_box.toml` | generated, 1 728 mixed cells | hexahedra, prisms, tetrahedra and transitional polyhedra in one mesh | [gallery 3](../gallery/meshes.md#mesh-types-3d) |
| `mesh/hex_2d.toml`, `tet_2d.toml`, `poly_2d.toml`, `mixed_2d.toml` | generated, 2D (one cell thick) | the four types in 2D, patch renaming | [gallery 5](../gallery/meshes.md#mesh-types-2d) |
| `mesh/poly_partition.toml` | generated, 2 197 polyhedra | METIS partition into 4 parts, halo layout, point location | [gallery 4](../gallery/meshes.md#poly-partition), [recipe 9](../recipes/meshes.md#mesh-partition) |
| `mesh/tjunction_tet.toml` | Gmsh, `meshes/tjunction.geo`, 77 776 tetrahedra | Gmsh reader, physical groups as patches and zones, RCM renumbering, 4 parts | [gallery 6](../gallery/meshes.md#tjunction), [recipe 8](../recipes/meshes.md#mesh-gmsh) |
| `mesh/tjunction_poly.toml` | Gmsh, the same, converted into 15 927 polyhedra | conversion of tetrahedra into polyhedra with the exact volume | [gallery 6](../gallery/meshes.md#tjunction), [recipe 8](../recipes/meshes.md#mesh-gmsh) |

## Mesh case-file keys

Every key is optional; a case needs either `mesh.file` or a
`[mesh.generate]` table. Paths are relative to the case file.

| Key | Default | Meaning |
|---|---|---|
| `mesh.file` | – | Gmsh MSH 4.1 file, ASCII or binary (D33) |
| `mesh.scale` | 1.0 | factor on the coordinates of the file (0.001: millimetres to metres) |
| `mesh.polyhedral` | false | convert a tetrahedral mesh into polyhedra (dual mesh, one cell per node) |
| `mesh.renumber` | `"none"` | cell order: `none` (as read), `rcm` (reverse Cuthill-McKee), `shuffle` |
| `mesh.seed` | 2026 | seed of `shuffle` |
| `[mesh.patches]` | – | patch renaming, `old = "new"` |
| `mesh.generate.type` | `"hex"` | `hex`, `tet`, `poly` or `mixed` |
| `mesh.generate.dimension` | 3 | 2: one cell layer thick in z (D3) |
| `mesh.generate.lower`, `upper` | `[0,0,0]`, `[1,1,1]` | corners of the box |
| `mesh.generate.cells` | `[8, 8, 8]` | intervals of the base hexahedral grid per direction (z ignored in 2D) |
| `mesh.generate.grading.x` (`.y`, `.z`) | `{ ratio = 1 }` | `ratio`: last to first interval; with `symmetric = true`: middle to end intervals |
| `mesh.generate.refine` | 0 | halve every interval this many times (nested meshes for grid-convergence tests) |
| `mesh.generate.perturb` | 0 | random node displacement as a fraction of the local interval, 0 to 0.3; boundary nodes stay on their faces |
| `mesh.generate.wedge` | none | revolve the 2D mesh (`dimension = 2`, `lower[1] ≥ 0`) into an axisymmetric wedge of this angle (degrees, up to 90) about the x axis; `zmin` and `zmax` become the wedge planes, `ymin` disappears on the axis |
| `mesh.generate.seed` | 2026 | seed of the perturbation (counter-based, the same on every machine) |
| `partition.parts` | 1 | number of parts (1: no decomposition) |
| `partition.method` | `"metis"` | `metis` or `rcb` (recursive coordinate bisection) |
| `partition.halo_layers` | 1 | halo depth |
| `check.locate_points` | 1000 | random points of the point-location test (0: none) |
| `check.seed` | 2026 | seed of the points |
| `output.directory` | `output/<case name>` | report directory |
| `output.vtkhdf` | true | write `mesh.vtkhdf` (fields: volume, non-orthogonality, skewness, aspect ratio, volume ratio, closure, rank) |

The generator names the six sides `xmin`, `xmax`, `ymin`, `ymax`, `zmin`,
`zmax`. A Gmsh file gives its surface physical groups as patches (or
`surface_<tag>` for surfaces without one) and its volume physical groups as
cell zones.

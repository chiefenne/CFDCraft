# Tutorial: your first mesh

This takes about ten minutes. You build Foehncraft, generate a small
polyhedral mesh, check it, partition it and look at it in ParaView. It is the
shortest path through the code; the [example gallery](gallery.md) shows the
results and the [recipes](recipes.md) give the other tasks.

## 1. Build

You need Rust 1.89 or newer, a C compiler and CMake (see
[recipe 1](recipes/benchmarks.md#g1-cpu)).

```
git clone https://github.com/chiefenne/Foehncraft.git && cd Foehncraft
cargo build --release -p foehncraft-cli
target/release/foehn info
```

## 2. Generate, check and write a mesh

```
target/release/foehn mesh examples/cases/mesh/poly_partition.toml
```

This makes a polyhedral mesh of the unit cube (2 197 cells), splits it into
four parts with METIS and writes the report to
`output/poly_partition/mesh_report.md`. Open the report: the verdict should
say **valid**, and the total volume should be 1.

## 3. Check it independently

```
pip install vtk numpy
python3 tools/check_vtkhdf.py output/poly_partition/mesh.vtkhdf
```

The check reads the file with VTK, the library ParaView uses, and compares
the cell volumes with those written by Foehncraft. A relative difference near
1e-15 means the file is right.

## 4. Look at it in ParaView

Open `output/poly_partition/mesh.vtkhdf` in ParaView 5.13 or newer. Colour
the cells by `rank` to see the four parts, or by `volume`. Use the `internal`
block for the cells and the patch blocks for the walls.

Next steps:

- change the mesh in `examples/cases/mesh/poly_partition.toml`, for example
  `perturb` or `cells`, and run again;
- try the other mesh types with [`poly_box.toml`](cases/meshes.md#meshes-examplescasesmesh),
  or read a Gmsh file as in [recipe 8](recipes/meshes.md#mesh-gmsh);
- see what the mesh is good for in the [gallery](gallery.md).

[↑ back to top](#tutorial-your-first-mesh)

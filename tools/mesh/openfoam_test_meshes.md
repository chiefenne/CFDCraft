# OpenFOAM test meshes

The meshes of `crates/foehncraft-mesh/tests/data/openfoam` were made with
OpenFOAM's own utilities (not needed to run the tests):

- `channel_ascii`: `blockMesh` of the box [0, 2] x [0, 1] x [0, 1] with
  (4, 2, 2) cells; patches `inlet`, `outlet`, `walls` (wall) and `sides`
  (empty); cell zones `left` and `right` from `topoSet` and `setsToZones`.
- `channel_binary`: the same mesh after `foamFormatConvert` with
  `writeFormat binary;` in `controlDict` (the faces become `faceCompactList`).
- `dual_polyhedra`: a Gmsh tetrahedral unit cube (`.msh` 2.2 from
  `cgns_test_meshes.py`), `gmshToFoam`, then `polyDualMesh 80 -overwrite`
  (needs `system/fvSchemes` and `fvSolution`): 81 polyhedral cells.

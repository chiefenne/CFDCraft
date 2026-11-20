- **SU2 meshes and one interface for every mesh format (D187, gallery 127,
  recipe 127).** `mesh.file` reads SU2 native files (`.su2`, ASCII) besides
  Gmsh: tetrahedra, hexahedra, prisms and pyramids in 3D; a 2D file
  (triangles, quadrilaterals) is extruded one cell layer in z with the planes
  `zmin` and `zmax` (`mesh.thickness`, default the median edge length).
  Every `MARKER_TAG` becomes a boundary; `[mesh.boundaries]` maps it to the
  name of a boundary condition. The readers sit behind one format layer
  (`foehncraft_mesh::format`): the format comes from the new key
  `mesh.format` (`gmsh`, `su2`), else from the extension, else from the
  content, and every reader hands cells and named faces to one shared mesh
  builder, so a further format is one parser. `foehn mesh` takes an `.su2`
  file directly. The Gmsh reader moved behind the layer unchanged (a test
  compares fingerprints of all Gmsh fixtures). Gallery 127: a NACA 0012
  C-mesh from PyAero, read from SU2, with a laminar flow at Re 1000.

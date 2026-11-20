- **Industry `.msh` and VTK meshes, VTU output (D191, D192, gallery 134,
  recipe 134).** `mesh.file` reads the ASCII face-based industry `.msh`
  (polyhedra, zones; it shares `.msh` with Gmsh, the content decides) and VTK
  grids, legacy `.vtk` and XML `.vtu` in every encoding including zlib.
  `mesh.format` takes `industry_msh` and `vtk`. `output.vtu` and
  `foehn mesh --vtu` write the mesh as a `.vtu` with polyhedra and the names
  of zones and boundaries.

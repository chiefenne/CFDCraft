- **OpenFOAM and CGNS meshes (D188, gallery 129, recipe 129).** `mesh.file`
  reads an OpenFOAM case, `constant` or `polyMesh` directory (ASCII or
  binary, polyhedral cells, patches as boundaries, `cellZones` as cell zones)
  and CGNS HDF5 files (hexahedra, tetrahedra, prisms, pyramids, NGON/NFACE
  polyhedra; BCs and Gmsh physical names as boundaries). `mesh.format` takes
  `cgns` and `openfoam`. Decomposed OpenFOAM cases, `.gz` files and ADF
  files are refused with the fix.

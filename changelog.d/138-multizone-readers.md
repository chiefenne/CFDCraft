- **Multi-zone SU2 and CGNS meshes for rotor/stator runs (D197, gallery 138,
  recipe 138).** `mesh.file` reads all zones of an SU2 file (`NZONE= n`) or a
  CGNS base into one mesh, each zone a cell zone, the nodes of different zones
  apart, so the interface is two boundaries for `[[sliding]]`. New key
  `mesh.merge_nodes` joins coincident nodes for touching stationary zones.

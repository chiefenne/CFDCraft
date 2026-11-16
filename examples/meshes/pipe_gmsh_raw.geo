// Pipe of diameter 0.1 m and length 1 m in raw Gmsh tetrahedra (no
// optimisation), for the cell-quality remediation (D50): the case
// examples/cases/flow/pipe_gmsh_raw.toml and crates/foehncraft-flow/tests/
// quality.rs. Some faces have a skewness above 0.9, on which the steady
// solver diverges with the least-squares pressure gradient. The mesh file
// is in the repository (made with Gmsh 4.15.2; another version may give
// another mesh):
//
//   gmsh -3 examples/meshes/pipe_gmsh_raw.geo -format msh41 -o examples/meshes/pipe_gmsh_raw.msh

SetFactory("OpenCASCADE");
Mesh.MeshSizeMax = 0.03;
Mesh.MeshSizeMin = 0.002;
Mesh.MshFileVersion = 4.1;
Cylinder(1) = {0, 0, 0, 1.0, 0, 0, 0.05};
eps = 1e-6;
in_main() = Surface In BoundingBox{-eps, -1, -1, eps, 1, 1};
outlet() = Surface In BoundingBox{1 - eps, -1, -1, 1 + eps, 1, 1};
all() = Boundary{ Volume{1}; };
wall() = all();
wall() -= in_main();
wall() -= outlet();
Physical Surface("inlet_main") = in_main();
Physical Surface("outlet") = outlet();
Physical Surface("wall") = wall();
Physical Volume("fluid") = {1};

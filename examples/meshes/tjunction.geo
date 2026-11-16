// T-junction mixer: a main pipe (diameter 0.1 m, length 1 m) with a branch
// (diameter 0.05 m, length 0.3 m) joining at x = 0.3 m; tetrahedra.
// The internal-flow case of D45 (stage 7) will use this geometry.
//
//   gmsh -3 examples/meshes/tjunction.geo -format msh41 -o examples/meshes/tjunction.msh
//
// Physical groups name the boundaries: inlet_main, inlet_branch, outlet,
// wall; the volume is the zone fluid.

SetFactory("OpenCASCADE");
Mesh.MeshSizeMax = 0.01;     // halve for a finer mesh
Mesh.MeshSizeMin = 0.002;
Mesh.MshFileVersion = 4.1;
Mesh.ElementOrder = 1;

Cylinder(1) = {0, 0, 0, 1.0, 0, 0, 0.05};
Cylinder(2) = {0.3, 0, 0, 0, 0.3, 0, 0.025};
BooleanUnion(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

eps = 1e-6;
in_main()   = Surface In BoundingBox{-eps, -1, -1, eps, 1, 1};
outlet()    = Surface In BoundingBox{1 - eps, -1, -1, 1 + eps, 1, 1};
in_branch() = Surface In BoundingBox{0, 0.3 - eps, -1, 1, 0.3 + eps, 1};
all()       = Boundary{ Volume{3}; };
wall()      = all();
wall() -= in_main();
wall() -= outlet();
wall() -= in_branch();

Physical Surface("inlet_main") = in_main();
Physical Surface("inlet_branch") = in_branch();
Physical Surface("outlet") = outlet();
Physical Surface("wall") = wall();
Physical Volume("fluid") = {3};

// finer near the junction
Field[1] = Distance;
Field[1].PointsList = {};
Field[1].CurvesList = {};
Field[1].SurfacesList = in_branch();
Field[2] = Ball;
Field[2].Radius = 0.08;
Field[2].XCenter = 0.3; Field[2].YCenter = 0.05; Field[2].ZCenter = 0;
Field[2].VIn = 0.005; Field[2].VOut = 0.01;
Background Field = 2;

// The T-junction of tjunction.geo with the tetrahedra optimised by Netgen:
// no faces with a skewness above 0.9, so the flow solver runs its schemes
// as set (no cell-quality remediation, D50). Recommended for flow runs.
//
//   gmsh -3 examples/meshes/tjunction_netgen.geo -format msh41 -o examples/meshes/tjunction_netgen.msh

Include "tjunction.geo";
Mesh.OptimizeNetgen = 1;

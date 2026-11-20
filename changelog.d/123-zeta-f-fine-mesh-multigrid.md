- **ζ-f on the finest meshes.** The multigrid setup of an equation with a
  wide range of coefficients (the f equation of k-ε-ζ-f on a fine mesh)
  refused a healthy coarsest level as singular, because its last pivot was
  compared with the largest diagonal. For the turbulence equations each pivot is
  now compared with its own diagonal (the pressure keeps the strict test). A transport or elliptic equation whose multigrid setup
  fails no longer stops the run: its solve takes the Jacobi preconditioner
  and the log notes it. The GCI study runs the level-3 ζ-f plate and channel with
  the relaxations 0.5, since both diverge with the defaults.

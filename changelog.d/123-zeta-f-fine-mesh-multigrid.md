- **ζ-f on the finest meshes.** The multigrid setup of an equation with a
  wide range of coefficients (the f equation of k-ε-ζ-f on a fine mesh)
  refused a healthy coarsest level as singular, because its last pivot was
  compared with the largest diagonal. Each pivot is now compared with its
  own diagonal. A transport or elliptic equation whose multigrid setup
  fails no longer stops the run: its solve takes the Jacobi preconditioner
  and the log notes it. The GCI study runs the level-3 ζ-f plate with the
  relaxations 0.5, since it diverges with the defaults.

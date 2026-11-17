- **k-ε-ζ-f model (stage 6 part 6d, D132).** `[turbulence] model =
  "k-epsilon-zeta-f"` (Hanjalić, Popovac & Hadžiabdić 2004): ζ = v²/k and
  the elliptic relaxation f added to k and ε, ν_t = C_μ ζ k T with the
  realizability and Kolmogorov bounds of T and L. Wall treatments
  `compound` (default, any y⁺: k and ζ by zero gradient, ε, P_k and f_w
  blended between sublayer and log law) and `low-re` (k = ζ = 0,
  ε_w = 2νk₁/d₁², f_w = −2νζ₁/d₁²). `zeta` at inlets (BC028 with another
  model) and in `[initial]`; ζ, f, T and v² in the output, the restart
  file and the log columns; a k-ε restart file continues with ζ = 2/3,
  f = 0. A solve lowers k, ε or ζ by at most a factor 10 per outer
  iteration (positivity safeguard).
- **Multigrid of the turbulence equations (D139).** Hierarchy slots on
  both backends (pressure, f, k, ε or ω, ζ). The elliptic f equation runs
  CG with the multigrid by default (`[solver.elliptic]`); the transport
  equations take it with `solver.turbulence.preconditioner = "amg"`,
  built from the symmetric part ½(A + Aᵀ). All hierarchies rebuild on the
  schedule of the pressure (`solver.amg_rebuild`), and their matrices go
  into the restart file (bit-identical continuation). The run-time table
  shows the turbulence model and the passive scalars.
- **Verification (gallery 73, recipe 73).** Channel at Re_τ = 395 resolved
  to y⁺ 0.27: U_b⁺ +1.9 % (`compound`) and +1.4 % (`low-re`) from the DNS.
  TMR flat plate at y⁺ 1: C_f(0.97) +5.1 % from White (finding, outside
  ±3 %). Findings: from the TMR free stream the boundary layer stays
  laminar (the case starts turbulent); near the leading edge the floor of
  k acts (TURB002).

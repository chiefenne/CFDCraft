- **Heated channel at every DNS point (stage 6 part 6f).** The channel of
  `channel_heated_keps.toml` at Re_τ 180, 395, 640 and 1020 and Pr 0.71
  and 0.025 against the DNS of Abe, Kawamura & Matsuo, on uniform meshes
  with first-cell y⁺ ≈ 30: Nu within −2.6 … +1.0 % at Pr 0.71 with both
  Pr_t models; at Pr 0.025 +5.9 … +7.2 % with Kays–Crawford and
  +12 … +22 % with the constant Pr_t. New case `channel_heated_pr0025.toml`;
  `examples/studies/channel_heat_dns.py --run` runs all sixteen points.
  The thermal law now takes T⁺ at the scalable treatment's floor
  ỹ⁺ = max(y⁺, 11.06), as the momentum law does (channel Nu on a y⁺ ≈ 8
  mesh +10 % instead of +28 %). CG stops when the preconditioned residual
  vanishes below the float32 range of the multigrid (an exactly
  one-dimensional flow), as BiCGStab did. The PTX module now has the
  entry points of the two passive-scalar kernels of pull request 85
  (`scalar_turb_faces`, `scalar_wall_faces`): without them a GPU run did
  not load the module; a test checks that the module has every entry the
  CUDA backend loads. Gallery and recipe 87.

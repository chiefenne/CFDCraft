- **Heated channel against DNS (stage 6 part 6f, specification 9.4).**
  `channel_heated_keps.toml` is the Re_τ = 395 channel with a passive
  temperature (Pr 0.71, uniform wall heat flux, periodic heat-flux jump).
  Its Nusselt number on the width 2h is 36.01 with the constant Pr_t and
  35.65 with Kays–Crawford, against 36.39 of the DNS of Abe, Kawamura &
  Matsuo (−1.0 % and −2.0 %). `examples/studies/channel_heat_dns.py`
  fetches the DNS tables of the JAXA database into a git-ignored folder
  and recomputes their Nusselt numbers for Re_τ 180 to 1020 at Pr 0.71
  and 0.025; no DNS data enters the repository (JAXA reserves the
  rights). Gallery and recipe 86.

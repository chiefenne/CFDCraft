- **GPU parity check of the turbulence models (specification 6.6).**
  `tools/bench/stage6_cloud.py` and `tools/notebooks/stage6_colab.ipynb`
  run every turbulence model on the TMR plate and the channel, the heated
  channel and pipe and the coupled step on the CPU and the GPU and compare
  each pair bitwise; `--cpu-only` compares one thread with four.
  `tools/check/compare_runs.py` no longer reports the hash of
  `case.resolved.toml` in `result.json` as a difference when only the
  device keys of `[run]` differ (`--threads`, `--backend`).

# Benchmarks (`examples/cases/bench/`)

Part of the [Case catalogue](../cases.md): gate G1, the kernel benchmark; benchmark keys.

On this page:

- [Case-file keys](#case-file-keys)

Gate G1, the kernel benchmark of stage 0 ([G1](../benchmarks/g1.md)). Run with
`foehn bench g1 FILE`; every key is optional, command-line options
override the file (`foehn bench g1 --help`).

| Case file | Backends | Mesh | Shows | Where |
|---|---|---|---|---|
| `bench/g1_ci_small.toml` | CPU | synthetic, 0.2 M cells, shuffled | smoke test of the harness: all backends bitwise equal; run by CI on every platform | CI, [recipe 1](../recipes/benchmarks.md#g1-cpu) |
| `bench/g1_cpu_10m.toml` | CPU | synthetic, 10 M cells, natural | the CPU half of G1: single-source bodies against hand-written loops | [gallery 1](../gallery/benchmarks.md#g1-cpu), G1 results |
| `bench/g1_cpu_10m_shuffled.toml` | CPU | synthetic, 10 M cells, shuffled | worst-case memory access, the cost of a badly numbered mesh | [gallery 1](../gallery/benchmarks.md#g1-cpu) |
| `bench/g1_gpu_10m.toml` | GPU and CPU | synthetic, 10 M cells, natural | **the G1 verdict**: single-source Rust PTX against hand-written CUDA C | [recipe 3](../recipes/benchmarks.md#g1-gpu), G1 procedure |
| `bench/g1_gpu_reference_4m.toml` | GPU (reference only) and CPU | synthetic, 4 M cells, natural | the CUDA path on a small or pre-Volta GPU: GPU results equal the CPU bitwise; no verdict | [recipe 2](../recipes/benchmarks.md#g1-windows), G1 maintainer's machine |
| `bench/g1_gpu_50m.toml` | GPU | synthetic, 50 M cells, natural | G1 at the upper end of the stage-plan range (17 GB GPU memory) | [recipe 3](../recipes/benchmarks.md#g1-gpu) |
| `bench/g1_cpu_poly_2m.toml` | CPU | generated, 2 M polyhedra, natural | G1 rerun on polyhedra (stage 1): uneven rows, 14 faces per cell | [gallery 7](../gallery/benchmarks.md#g1-stage1), [recipe 10](../recipes/benchmarks.md#g1-stage1) |
| `bench/g1_cpu_poly_2m_shuffled.toml` | CPU | generated, 2 M polyhedra, shuffled | the same, worst-case ordering | [gallery 7](../gallery/benchmarks.md#g1-stage1), [recipe 10](../recipes/benchmarks.md#g1-stage1) |
| `bench/g1_cpu_tet_2m.toml` | CPU | generated, 2 M tetrahedra, natural | G1 rerun on tetrahedra | [gallery 7](../gallery/benchmarks.md#g1-stage1), [recipe 10](../recipes/benchmarks.md#g1-stage1) |
| `bench/g1_gpu_poly_4m.toml` | GPU and CPU | generated, 4 M polyhedra, natural | the GPU verdict rerun on polyhedra (fits a Colab T4) | [recipe 10](../recipes/benchmarks.md#g1-stage1) |

## Case-file keys

| Key | Default | Meaning |
|---|---|---|
| `mesh.type` | `"periodic-box"` | `periodic-box` (stage 0), or a generated mesh of the unit cube: `hex`, `tet`, `poly`, `mixed` (inner faces only, no boundary) |
| `mesh.cells` | 10 000 000 | approximate cell count (of the near-cubic box, or of the generated mesh) |
| `mesh.nx`, `mesh.ny`, `mesh.nz` | – | exact box dimensions (all three, at least 3 each) |
| `mesh.ordering` | `"natural"` | `"natural"` or `"shuffled"` cell numbering |
| `mesh.seed` | 2026 | seed of the fields and of the shuffle |
| `run.backends` | the three CPU backends | see `foehn bench backends` |
| `run.kernels` | all four | `face-flux`, `cell-gather`, `spmv-csr`, `spmv-sell` |
| `run.repeat` | 10 | timed runs per kernel (median and minimum reported) |
| `run.warmup` | 2 | untimed runs first |
| `run.threads` | all cores | CPU threads |
| `run.device` | 0 | CUDA device ordinal |
| `run.rust_ptx` | `target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx` | PTX of the single-source GPU kernels |
| `run.tolerance` | 1e-12 | allowed relative difference to `cpu-serial` |
| `sell.chunk` | 32 | SELL slice height C |
| `sell.sigma` | 1 | SELL sorting window σ |
| `output.directory` | `"output/g1"` | report directory |
| `output.name` | `"g1"` | base name of `NAME.json`, `NAME.md`, `NAME.svg` |

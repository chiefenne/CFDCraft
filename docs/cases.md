# Case catalogue

The case files in `examples/cases/` with what each one is, the mesh it needs
and where it is discussed. Pictures are in the [Example gallery](gallery.md),
command sequences in the [Recipes](recipes.md).

How to read the **Mesh** column:

* **synthetic** – the mesh is generated in memory from the `[mesh]` section
  of the case file; no mesh file is needed.

## Benchmarks (`examples/cases/bench/`)

Gate G1, the kernel benchmark of stage 0 ([G1](benchmarks/g1.md)). Run with
`cfdcraft bench g1 FILE`; every key is optional, command-line options
override the file (`cfdcraft bench g1 --help`).

| Case file | Backends | Mesh | Shows | Where |
|---|---|---|---|---|
| `bench/g1_ci_small.toml` | CPU | synthetic, 0.2 M cells, shuffled | smoke test of the harness: all backends bitwise equal; run by CI on every platform | CI, recipe 1 |
| `bench/g1_cpu_10m.toml` | CPU | synthetic, 10 M cells, natural | the CPU half of G1: single-source bodies against hand-written loops | gallery 1, G1 results |
| `bench/g1_cpu_10m_shuffled.toml` | CPU | synthetic, 10 M cells, shuffled | worst-case memory access, the cost of a badly numbered mesh | gallery 1 |
| `bench/g1_gpu_10m.toml` | GPU and CPU | synthetic, 10 M cells, natural | **the G1 verdict**: single-source Rust PTX against hand-written CUDA C | recipe 3, G1 procedure |
| `bench/g1_gpu_reference_4m.toml` | GPU (reference only) and CPU | synthetic, 4 M cells, natural | the CUDA path on a small or pre-Volta GPU: GPU results equal the CPU bitwise; no verdict | recipe 2, G1 maintainer's machine |
| `bench/g1_gpu_50m.toml` | GPU | synthetic, 50 M cells, natural | G1 at the upper end of the stage-plan range (17 GB GPU memory) | recipe 3 |

### Case-file keys

| Key | Default | Meaning |
|---|---|---|
| `mesh.cells` | 10 000 000 | approximate cell count of a near-cubic periodic box |
| `mesh.nx`, `mesh.ny`, `mesh.nz` | – | exact box dimensions (all three, at least 3 each) |
| `mesh.ordering` | `"natural"` | `"natural"` or `"shuffled"` cell numbering |
| `mesh.seed` | 2026 | seed of the fields and of the shuffle |
| `run.backends` | the three CPU backends | see `cfdcraft bench backends` |
| `run.kernels` | all four | `face-flux`, `cell-gather`, `spmv-csr`, `spmv-sell` |
| `run.repeat` | 10 | timed runs per kernel (median and minimum reported) |
| `run.warmup` | 2 | untimed runs first |
| `run.threads` | all cores | CPU threads |
| `run.device` | 0 | CUDA device ordinal |
| `run.rust_ptx` | `target/nvptx64-nvidia-cuda/release/cfdcraft_ptx_kernels.ptx` | PTX of the single-source GPU kernels |
| `run.tolerance` | 1e-12 | allowed relative difference to `cpu-serial` |
| `sell.chunk` | 32 | SELL slice height C |
| `sell.sigma` | 1 | SELL sorting window σ |
| `output.directory` | `"output/g1"` | report directory |
| `output.name` | `"g1"` | base name of `NAME.json`, `NAME.md`, `NAME.svg` |

# Example gallery

What CFDCraft computes, one picture per case, with the commands that
reproduce it. Command sequences for everyday tasks are in
[Recipes](recipes.md); all case files are listed in the
[Case catalogue](cases.md).

All commands run from the repository root. Case files are in
`examples/cases/`.

| # | Example | Stage | Feature |
|---|---|---|---|
| 1 | [Kernel benchmark G1 on the CPU](#g1-cpu) | 0 | gate G1, reproducible kernels |
| 2 | [Kernel benchmark G1 on a Tesla T4](#g1-gpu) | 0 | gate G1, single-source Rust on the GPU |

---

<a id="g1-cpu"></a>
## 1. Kernel benchmark G1 on the CPU

![G1 effective bandwidth of the four kernels on 10 million cells](figures/gallery/g1_cpu_10m.svg)

The four kernels that dominate a finite-volume solver — face flux, cell
gather and the sparse matrix-vector product in CSR and SELL-C-σ form — on a
periodic box of 10 million hexahedra in face-based form. The single-source
kernel bodies (`cpu-threads`) run as fast as hand-written Rust loops
(`cpu-native`) and give bitwise the same results as the single-thread run
(`cpu-serial`), for any number of threads. The bars are effective bandwidth:
bytes that must move over the median time. Procedure, numbers and reading:
[G1](benchmarks/g1.md).

| | |
|---|---|
| Mesh | periodic box 215 × 215 × 216 = 9 984 600 cells, 29 953 800 faces, natural numbering |
| Machine | Intel Xeon @ 2.10 GHz, 4 virtual cores (shared cloud machine) |
| Case file | `bench/g1_cpu_10m.toml` |

```
cargo build --release -p cfdcraft-cli
target/release/cfdcraft bench g1 examples/cases/bench/g1_cpu_10m.toml
# picture: output/g1_cpu_10m/g1.svg, numbers: g1.md and g1.json
```

The same run with a randomly numbered mesh (`bench/g1_cpu_10m_shuffled.toml`)
shows what a badly ordered mesh costs: three to five times less bandwidth.

![G1 with shuffled cell numbering](figures/gallery/g1_cpu_10m_shuffled.svg)

[↑ back to top](#example-gallery)

---

<a id="g1-gpu"></a>
## 2. Kernel benchmark G1 on a Tesla T4

![G1 effective bandwidth on a Tesla T4, 10 million cells](figures/gallery/g1_gpu_t4_10m.svg)

The same four kernels on a Tesla T4 (Google Colab, free tier). Blue: the
hand-written CUDA C kernels; orange: the single-source Rust kernels compiled
to PTX. Both reach 220–255 GB/s of the T4's 320 GB/s, the Rust kernels are
within 1.0 % of the CUDA C kernels (gate limit: 10 %), and all results are
bitwise equal to the CPU. The two short bars per group are the CPU backends
of the Colab machine (2 virtual cores). Reading and limits:
[G1](benchmarks/g1.md#gpu-verdict-tesla-t4-google-colab-2026-11-14).

| | |
|---|---|
| Mesh | periodic box 215 × 215 × 216 = 9 984 600 cells, natural numbering |
| Machine | Tesla T4 (sm_75, 16 GB), Intel Xeon @ 2.00 GHz, 2 threads |
| Case file | `bench/g1_gpu_10m.toml` |

```
bash tools/g1_cloud.sh      # on any Linux machine with an NVIDIA GPU of compute capability 7.0+
# or open tools/g1_colab.ipynb in Google Colab (recipe 3)
```

[↑ back to top](#example-gallery)

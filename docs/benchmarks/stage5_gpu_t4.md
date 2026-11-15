# Stage 5: flow loop on the GPU
GPU: Tesla T4 (sm_75); CPU: 2 threads. Wall times include mesh generation and output.
| case | cells | CPU s | GPU s | CPU ns | GPU ns | bitwise equal |
|---|---:|---:|---:|---:|---:|---|
| couette_poly | 153 | 0.2 | 1.5 | 12358 | 20594 ns/cell/iteration | yes |
| poiseuille_mixed | 12288 | 43.1 | 15.3 | 2925 | 855 ns/cell/iteration | yes |
| cavity_re100 | 4096 | 7.3 | 2.9 | 4011 | 1389 ns/cell/iteration | yes |
| cavity_re100_poly3d | 5445 | 6.4 | 2.3 | 5848 | 1559 ns/cell/iteration | yes |
| cavity_re100_transient | 4096 | 28.9 | 13.0 | 14085 | 6140 ns/cell/step | yes |
| taylor_green_cn | 4096 | 46.4 | 15.1 | 8678 | 2621 ns/cell/step | yes |
| cylinder_re100_gpu | 98304 | 191.8 | 22.2 | 9692 | 1049 ns/cell/step | yes |
| cavity_poly3d_gpu | 970225 | 200.1 | 136.4 | 5116 | 722 ns/cell/iteration | yes |
all GPU runs bitwise equal to the CPU runs

Measured 2026-11-15 on Google Colab (Tesla T4, 2 CPU threads), branch `stage5-gpu-flow-loop`, `tools/stage5_colab.ipynb`.

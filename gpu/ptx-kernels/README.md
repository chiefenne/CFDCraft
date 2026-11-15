# Single-source GPU kernels

The kernel bodies of [`foehncraft-kernels`](../../crates/foehncraft-kernels),
wrapped as PTX entry points for NVIDIA GPUs. The same bodies run in the CPU
backends. Entry points: `g1_*` (gate G1, stage 0; `src/lib.rs`) and `fv_*`,
`blas_*` (the scalar operators and Krylov vector operations of stage 2,
`src/fv.rs`, launched by `crates/foehncraft-fv/src/backend/cuda.rs`).

For G1, this crate is the "single-source Rust" half of
gate G1 (requirements D8, D11); the other half is the hand-written CUDA C
reference in [`crates/foehncraft-bench/kernels/g1_reference.cu`](../../crates/foehncraft-bench/kernels/g1_reference.cu).

```
gpu/ptx-kernels/build.sh          # sm_70 PTX, runs on Volta and every newer GPU
gpu/ptx-kernels/build.sh sm_86    # tuned for one architecture
```

Needs the nightly toolchain (`rust-toolchain.toml` installs it): the
`ptx-kernel` ABI and the NVPTX intrinsics are not stable Rust yet. The rustc
NVPTX target supports compute capability 7.0 (Volta) and newer; older GPUs
can run the CUDA reference but not this module.

Not a member of the main Cargo workspace, because it is built for another
target with another toolchain. CI builds it on every pull request and checks
that all entry points are present; `tools/check_kernels.sh` then checks
that the PTX has no atomic instructions and no fused or approximate f64
arithmetic (R10).

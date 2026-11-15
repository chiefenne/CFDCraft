// Gate G1: hand-written CUDA C reference kernels.
//
// The way an experienced CUDA programmer writes these kernels: one thread
// per output, __restrict__ and read-only loads (__ldg), 32-bit indices.
// They compute exactly what the single-source Rust bodies in
// crates/inncraft-kernels compute, in the same order, so the results are
// bitwise comparable when both are compiled without contraction of a*b+c
// into fused multiply-add (NVRTC option --fmad=false, set by the harness).
//
// The argument lists are identical to the Rust PTX entry points
// (gpu/ptx-kernels); the size arguments that the Rust side needs to build
// its slices are unused here.
//
// Compiled at run time by NVRTC (crates/inncraft-bench/src/backend/cuda.rs).

#define FLAG_BIT 0x80000000u

extern "C" __global__ void g1_face_flux(
    unsigned n_faces, unsigned n_cells,
    const unsigned* __restrict__ owner, const unsigned* __restrict__ neighbour,
    const double* __restrict__ phi, const double* __restrict__ gamma,
    const double* __restrict__ mass_flux, double* __restrict__ flux)
{
    unsigned f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= n_faces) return;
    unsigned p = __ldg(owner + f);
    unsigned n = __ldg(neighbour + f);
    double phi_p = __ldg(phi + p);
    double phi_n = __ldg(phi + n);
    double m = __ldg(mass_flux + f);
    double g = __ldg(gamma + f);
    double mp = m > 0.0 ? m : 0.0;
    double mn = m < 0.0 ? m : 0.0;
    flux[f] = mp * phi_p + mn * phi_n - g * (phi_n - phi_p);
}

extern "C" __global__ void g1_cell_gather(
    unsigned n_cells, unsigned n_entries, unsigned n_faces,
    const unsigned* __restrict__ cell_ptr, const unsigned* __restrict__ cell_faces,
    const double* __restrict__ flux, double* __restrict__ out)
{
    unsigned c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c >= n_cells) return;
    unsigned lo = __ldg(cell_ptr + c);
    unsigned hi = __ldg(cell_ptr + c + 1);
    double s = 0.0;
    for (unsigned k = lo; k < hi; ++k) {
        unsigned e = __ldg(cell_faces + k);
        double v = __ldg(flux + (e & ~FLAG_BIT));
        if (e & FLAG_BIT) s -= v; else s += v;
    }
    out[c] = s;
}

extern "C" __global__ void g1_spmv_csr(
    unsigned n, unsigned nnz,
    const unsigned* __restrict__ row_ptr, const unsigned* __restrict__ col,
    const double* __restrict__ val, const double* __restrict__ x,
    double* __restrict__ y)
{
    unsigned r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= n) return;
    unsigned lo = __ldg(row_ptr + r);
    unsigned hi = __ldg(row_ptr + r + 1);
    double s = 0.0;
    for (unsigned k = lo; k < hi; ++k)
        s += __ldg(val + k) * __ldg(x + __ldg(col + k));
    y[r] = s;
}

extern "C" __global__ void g1_spmv_sell(
    unsigned n, unsigned chunk, unsigned n_slices, unsigned stored,
    const unsigned* __restrict__ slice_ptr, const unsigned* __restrict__ slice_len,
    const unsigned* __restrict__ col, const double* __restrict__ val,
    const double* __restrict__ x, double* __restrict__ y_sorted)
{
    unsigned r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= n) return;
    unsigned s_idx = r / chunk;
    unsigned lane = r - s_idx * chunk;
    unsigned base = __ldg(slice_ptr + s_idx) + lane;
    unsigned len = __ldg(slice_len + s_idx);
    double s = 0.0;
    for (unsigned j = 0; j < len; ++j) {
        unsigned k = base + j * chunk;
        s += __ldg(val + k) * __ldg(x + __ldg(col + k));
    }
    y_sorted[r] = s;
}

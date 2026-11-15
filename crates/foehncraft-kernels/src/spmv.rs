//! Sparse matrix-vector product y = A x, the core of every Krylov and
//! multigrid iteration (D35, D36), in two storage formats (D37):
//!
//! - CSR (compressed sparse row), one thread per row;
//! - SELL-C-σ (Kreutzer, Hager, Wellein, Fehske, Bishop, *A unified sparse
//!   matrix data format for efficient general sparse matrix-vector
//!   multiplication on modern processors with wide SIMD units*, SIAM J. Sci.
//!   Comput. 36, C401–C423, 2014): rows sorted by length inside windows of σ
//!   rows, grouped into slices of C rows, each slice stored column-major and
//!   padded to its longest row. With C equal to the GPU warp width, the loads
//!   of one warp are coalesced.

use crate::ld;
use foehncraft_core::Real;

/// Row `r` of y = A x for a CSR matrix.
///
/// # Safety
/// `r + 1 < row_ptr.len()`, `row_ptr` non-decreasing with entries
/// `<= col.len() == val.len()`, every column index `< x.len()`.
#[inline(always)]
pub unsafe fn csr_row<T: Real>(r: usize, row_ptr: &[u32], col: &[u32], val: &[T], x: &[T]) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(row_ptr, r) as usize;
        let hi = ld(row_ptr, r + 1) as usize;
        let mut s = T::ZERO;
        let mut k = lo;
        while k < hi {
            s += ld(val, k) * ld(x, ld(col, k) as usize);
            k += 1;
        }
        s
    }
}

/// Sorted row `r` of y = A x for a SELL-C-σ matrix with slice height
/// `chunk`. The result belongs to the original row `perm[r]`; the caller
/// writes it there. Padding entries have value zero and a valid column.
///
/// # Safety
/// `r / chunk < slice_ptr.len() == slice_len.len()`, and for every
/// `j < slice_len[s]` the entry `slice_ptr[s] + j * chunk + r % chunk` is
/// `< col.len() == val.len()`; every column index `< x.len()`.
#[inline(always)]
pub unsafe fn sell_row<T: Real>(
    r: usize,
    chunk: usize,
    slice_ptr: &[u32],
    slice_len: &[u32],
    col: &[u32],
    val: &[T],
    x: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let s_idx = r / chunk;
        let lane = r - s_idx * chunk;
        let base = ld(slice_ptr, s_idx) as usize + lane;
        let len = ld(slice_len, s_idx) as usize;
        let mut s = T::ZERO;
        let mut j = 0;
        while j < len {
            let k = base + j * chunk;
            s += ld(val, k) * ld(x, ld(col, k) as usize);
            j += 1;
        }
        s
    }
}

//! Stage 3: the multigrid V-cycle bodies of `ennecraft_kernels::amg` as PTX
//! entry points, instantiated for float32 (`*_f32`, the mixed-precision
//! hierarchy of D37) and float64 (`*_f64`). The host side is
//! `crates/ennecraft-fv/src/backend/cuda.rs`; the argument lists there must
//! match these signatures exactly. This is the only module with f32
//! arithmetic (R10).
//!
//! Matrices are sliced-ELL (σ = 1): `chunk`, `n_slices`, `stored`, then the
//! slice pointers, slice lengths, columns and values. Every output is written
//! by one thread; the in-place updates read only their own entry.

use crate::tid;
use ennecraft_kernels::amg::{self, Sell};
use core::slice::from_raw_parts as sl;

/// Matrix from raw parts.
///
/// # Safety
/// Device arrays of the given lengths.
#[inline(always)]
unsafe fn sell<'a, T>(
    chunk: u32,
    n_slices: u32,
    stored: u32,
    slice_ptr: *const u32,
    slice_len: *const u32,
    col: *const u32,
    val: *const T,
) -> Sell<'a, T> {
    let (ns, st) = (n_slices as usize, stored as usize);
    // SAFETY: guaranteed by the caller.
    unsafe {
        Sell {
            chunk: chunk as usize,
            slice_ptr: sl(slice_ptr, ns),
            slice_len: sl(slice_len, ns),
            col: sl(col, st),
            val: sl(val, st),
        }
    }
}

macro_rules! amg_entries {
    ($t:ty, $first:ident, $step:ident, $add:ident, $res:ident, $spmv:ident, $prolong:ident,
     $dense:ident, $to:ident, $from:ident $(,)?) => {
        /// First smoother step from zero: d = x = c₂ D⁻¹ b.
        #[unsafe(no_mangle)]
        pub unsafe extern "ptx-kernel" fn $first(
            n: u32,
            c2: $t,
            dinv: *const $t,
            b: *const $t,
            d: *mut $t,
            x: *mut $t,
        ) {
            let i = tid();
            let n = n as usize;
            if i >= n {
                return;
            }
            // SAFETY: see the module docs.
            unsafe {
                let v = amg::smooth_first(i, sl(dinv, n), sl(b, n), c2);
                *d.add(i) = v;
                *x.add(i) = v;
            }
        }

        /// Smoother step: d = c₁ d + c₂ D⁻¹ (b − A x), in place.
        #[unsafe(no_mangle)]
        pub unsafe extern "ptx-kernel" fn $step(
            n: u32,
            chunk: u32,
            n_slices: u32,
            stored: u32,
            slice_ptr: *const u32,
            slice_len: *const u32,
            col: *const u32,
            val: *const $t,
            c1: $t,
            c2: $t,
            dinv: *const $t,
            b: *const $t,
            x: *const $t,
            d: *mut $t,
        ) {
            let i = tid();
            let n = n as usize;
            if i >= n {
                return;
            }
            // SAFETY: see the module docs.
            unsafe {
                let a = sell(chunk, n_slices, stored, slice_ptr, slice_len, col, val);
                let v = amg::smooth_step(i, &a, sl(dinv, n), sl(b, n), sl(x, n), *d.add(i), c1, c2);
                *d.add(i) = v;
            }
        }

        /// x = x + d, in place.
        #[unsafe(no_mangle)]
        pub unsafe extern "ptx-kernel" fn $add(n: u32, d: *const $t, x: *mut $t) {
            let i = tid();
            if i >= n as usize {
                return;
            }
            // SAFETY: see the module docs.
            unsafe {
                *x.add(i) = amg::add(*x.add(i), *d.add(i));
            }
        }

        /// w = b − A x.
        #[unsafe(no_mangle)]
        pub unsafe extern "ptx-kernel" fn $res(
            n: u32,
            chunk: u32,
            n_slices: u32,
            stored: u32,
            slice_ptr: *const u32,
            slice_len: *const u32,
            col: *const u32,
            val: *const $t,
            b: *const $t,
            x: *const $t,
            w: *mut $t,
        ) {
            let i = tid();
            let n = n as usize;
            if i >= n {
                return;
            }
            // SAFETY: see the module docs.
            unsafe {
                let a = sell(chunk, n_slices, stored, slice_ptr, slice_len, col, val);
                *w.add(i) = amg::residual(i, &a, sl(b, n), sl(x, n));
            }
        }

        /// y = M x for a rectangular matrix M (the restriction).
        #[unsafe(no_mangle)]
        pub unsafe extern "ptx-kernel" fn $spmv(
            n_rows: u32,
            n_cols: u32,
            chunk: u32,
            n_slices: u32,
            stored: u32,
            slice_ptr: *const u32,
            slice_len: *const u32,
            col: *const u32,
            val: *const $t,
            x: *const $t,
            y: *mut $t,
        ) {
            let i = tid();
            if i >= n_rows as usize {
                return;
            }
            // SAFETY: see the module docs.
            unsafe {
                let m = sell(chunk, n_slices, stored, slice_ptr, slice_len, col, val);
                *y.add(i) = m.row(i, sl(x, n_cols as usize));
            }
        }

        /// x = x + P x_c, in place.
        #[unsafe(no_mangle)]
        pub unsafe extern "ptx-kernel" fn $prolong(
            n_rows: u32,
            n_cols: u32,
            chunk: u32,
            n_slices: u32,
            stored: u32,
            slice_ptr: *const u32,
            slice_len: *const u32,
            col: *const u32,
            val: *const $t,
            xc: *const $t,
            x: *mut $t,
        ) {
            let i = tid();
            if i >= n_rows as usize {
                return;
            }
            // SAFETY: see the module docs.
            unsafe {
                let p = sell(chunk, n_slices, stored, slice_ptr, slice_len, col, val);
                *x.add(i) = amg::prolong_add(i, &p, sl(xc, n_cols as usize), *x.add(i));
            }
        }

        /// x = M b with a dense n × n matrix M (the coarsest level).
        #[unsafe(no_mangle)]
        pub unsafe extern "ptx-kernel" fn $dense(n: u32, m: *const $t, b: *const $t, x: *mut $t) {
            let i = tid();
            let n = n as usize;
            if i >= n {
                return;
            }
            // SAFETY: see the module docs.
            unsafe {
                *x.add(i) = amg::dense_row(i, n, sl(m, n * n), sl(b, n));
            }
        }

        /// dst = src rounded from float64 to the hierarchy precision.
        #[unsafe(no_mangle)]
        pub unsafe extern "ptx-kernel" fn $to(n: u32, src: *const f64, dst: *mut $t) {
            let i = tid();
            let n = n as usize;
            if i >= n {
                return;
            }
            // SAFETY: see the module docs.
            unsafe {
                *dst.add(i) = amg::convert::<f64, $t>(i, sl(src, n));
            }
        }

        /// dst = src widened to float64.
        #[unsafe(no_mangle)]
        pub unsafe extern "ptx-kernel" fn $from(n: u32, src: *const $t, dst: *mut f64) {
            let i = tid();
            let n = n as usize;
            if i >= n {
                return;
            }
            // SAFETY: see the module docs.
            unsafe {
                *dst.add(i) = amg::convert::<$t, f64>(i, sl(src, n));
            }
        }
    };
}

amg_entries!(
    f32,
    amg_smooth_first_f32,
    amg_smooth_step_f32,
    amg_add_f32,
    amg_residual_f32,
    amg_spmv_f32,
    amg_prolong_f32,
    amg_dense_f32,
    amg_to_f32,
    amg_from_f32,
);

amg_entries!(
    f64,
    amg_smooth_first_f64,
    amg_smooth_step_f64,
    amg_add_f64,
    amg_residual_f64,
    amg_spmv_f64,
    amg_prolong_f64,
    amg_dense_f64,
    amg_to_f64,
    amg_from_f64,
);

//! Vector operations of the Krylov solvers (D36), one output per index.
//!
//! The fused updates of BiCGStab (van der Vorst 1992) are written out so that
//! each step is one pass over memory. Dot products are summed in blocks of
//! [`cfdcraft_core::repro::BLOCK`] values in sequential order
//! ([`block_dot`]); the block sums are then combined in the fixed tree of
//! [`cfdcraft_core::repro::combine_block_sums`]. The CPU and the GPU use the
//! same blocks, so a dot product is bitwise the same on both (D31).

use crate::ld;
use cfdcraft_core::Real;
use cfdcraft_core::repro::BLOCK;

/// Sum of `a[i] b[i]` over block `blk` (indices `blk*BLOCK..` up to the end).
///
/// # Safety
/// `blk * BLOCK < a.len() == b.len()`.
#[inline(always)]
pub unsafe fn block_dot<T: Real>(blk: usize, a: &[T], b: &[T]) -> T {
    let lo = blk * BLOCK;
    let hi = if lo + BLOCK < a.len() {
        lo + BLOCK
    } else {
        a.len()
    };
    let mut s = T::ZERO;
    let mut i = lo;
    while i < hi {
        // SAFETY: i < a.len() == b.len().
        unsafe {
            s += ld(a, i) * ld(b, i);
        }
        i += 1;
    }
    s
}

/// Minimum and maximum of block `blk` of `a` (order-independent).
///
/// # Safety
/// `blk * BLOCK < a.len()`.
#[inline(always)]
pub unsafe fn block_min_max<T: Real>(blk: usize, a: &[T]) -> (T, T) {
    let lo = blk * BLOCK;
    let hi = if lo + BLOCK < a.len() {
        lo + BLOCK
    } else {
        a.len()
    };
    // SAFETY: lo < a.len().
    let first = unsafe { ld(a, lo) };
    let (mut mn, mut mx) = (first, first);
    let mut i = lo + 1;
    while i < hi {
        // SAFETY: i < a.len().
        let v = unsafe { ld(a, i) };
        mn = mn.min_of(v);
        mx = mx.max_of(v);
        i += 1;
    }
    (mn, mx)
}

/// BiCGStab direction update: p_i = r_i + β (p_i − ω v_i).
///
/// # Safety
/// `i` is inside every slice.
#[inline(always)]
pub unsafe fn update_p<T: Real>(i: usize, r: &[T], p: &[T], v: &[T], beta: T, omega: T) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { ld(r, i) + beta * (ld(p, i) - omega * ld(v, i)) }
}

/// y_i = x_i − α z_i.
///
/// # Safety
/// `i` is inside every slice.
#[inline(always)]
pub unsafe fn sub_scaled<T: Real>(i: usize, x: &[T], z: &[T], alpha: T) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { ld(x, i) - alpha * ld(z, i) }
}

/// BiCGStab solution update: x_i + α p̂_i + ω ŝ_i.
///
/// # Safety
/// `i` is inside every slice.
#[inline(always)]
pub unsafe fn update_x<T: Real>(i: usize, x: &[T], ph: &[T], sh: &[T], alpha: T, omega: T) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { ld(x, i) + alpha * ld(ph, i) + omega * ld(sh, i) }
}

/// Jacobi preconditioner: d_i x_i with d = 1/diag.
///
/// # Safety
/// `i` is inside both slices.
#[inline(always)]
pub unsafe fn scale<T: Real>(i: usize, d: &[T], x: &[T]) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { ld(d, i) * ld(x, i) }
}

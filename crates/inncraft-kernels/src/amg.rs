//! Kernel bodies of the algebraic multigrid V-cycle (stage 3, D35, D37).
//!
//! Every matrix of the hierarchy (the level operators A_l, the
//! prolongations P_l and the restrictions R_l = P_lᵀ, stored explicitly) is
//! a sliced-ELL matrix, so every matrix operation is a row gather by
//! [`crate::spmv::sell_row`]. The bodies are generic over the floating-point
//! type: the hierarchy runs in float32 in the default mixed-precision mode
//! and in float64 for comparison (D37). This is the only kernel module that
//! is instantiated for float32 (R10).
//!
//! The smoothers are polynomial in D⁻¹A (damped Jacobi, Chebyshev), written
//! as a sequence of steps
//!
//! d ← c₁ d + c₂ D⁻¹ (b − A x),   x ← x + d,
//!
//! with coefficients (c₁, c₂) per step computed on the host. The first step
//! from a zero initial guess needs no product with A ([`smooth_first`]).

use crate::ld;
use crate::spmv::sell_row;
use inncraft_core::Real;

/// A sliced-ELL matrix (SELL-C-σ with σ = 1) as borrowed slices.
#[derive(Clone, Copy)]
pub struct Sell<'a, T> {
    /// Slice height C.
    pub chunk: usize,
    /// Start of each slice in `col` and `val`.
    pub slice_ptr: &'a [u32],
    /// Slots per row of each slice.
    pub slice_len: &'a [u32],
    /// Column of every slot (padding: a valid column with value zero).
    pub col: &'a [u32],
    /// Value of every slot.
    pub val: &'a [T],
}

impl<T: Real> Sell<'_, T> {
    /// Row `r` of the product with `x`.
    ///
    /// # Safety
    /// As [`crate::spmv::sell_row`].
    #[inline(always)]
    pub unsafe fn row(&self, r: usize, x: &[T]) -> T {
        // SAFETY: guaranteed by the caller.
        unsafe {
            sell_row(
                r,
                self.chunk,
                self.slice_ptr,
                self.slice_len,
                self.col,
                self.val,
                x,
            )
        }
    }
}

/// Smoother step, row `i`: c₁ dᵢ + c₂ D⁻¹ᵢ (bᵢ − (A x)ᵢ), with `d_i` the
/// current dᵢ. The result replaces dᵢ (it reads `d` at row `i` only, so in
/// place); it is added to `x` afterwards with [`add`], because the step
/// reads `x` at other rows.
///
/// # Safety
/// `i` is a row of `a`; `dinv` and `b` have an entry for `i`; the columns of
/// `a` are inside `x`.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub unsafe fn smooth_step<T: Real>(
    i: usize,
    a: &Sell<T>,
    dinv: &[T],
    b: &[T],
    x: &[T],
    d_i: T,
    c1: T,
    c2: T,
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let r = ld(b, i) - a.row(i, x);
        c1 * d_i + c2 * (ld(dinv, i) * r)
    }
}

/// First smoother step from x = 0, row `i`: c₂ D⁻¹ᵢ bᵢ (both d and x).
///
/// # Safety
/// `i` is inside both slices.
#[inline(always)]
pub unsafe fn smooth_first<T: Real>(i: usize, dinv: &[T], b: &[T], c2: T) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { c2 * (ld(dinv, i) * ld(b, i)) }
}

/// xᵢ + dᵢ (the update after a smoother step).
#[inline(always)]
pub fn add<T: Real>(x_i: T, d_i: T) -> T {
    x_i + d_i
}

/// Residual row: bᵢ − (A x)ᵢ.
///
/// # Safety
/// As [`smooth_step`].
#[inline(always)]
pub unsafe fn residual<T: Real>(i: usize, a: &Sell<T>, b: &[T], x: &[T]) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { ld(b, i) - a.row(i, x) }
}

/// Prolongation and correction, fine row `i`: xᵢ + (P x_c)ᵢ with `x_i` the
/// current xᵢ (written in place).
///
/// # Safety
/// `i` is a row of `p`; its columns are inside `xc`.
#[inline(always)]
pub unsafe fn prolong_add<T: Real>(i: usize, p: &Sell<T>, xc: &[T], x_i: T) -> T {
    // SAFETY: guaranteed by the caller.
    x_i + unsafe { p.row(i, xc) }
}

/// Row `i` of a dense n × n matrix times `b` (the inverse of the coarsest
/// operator), summed in column order.
///
/// # Safety
/// `i < n`, `m.len() >= n * n`, `b.len() >= n`.
#[inline(always)]
pub unsafe fn dense_row<T: Real>(i: usize, n: usize, m: &[T], b: &[T]) -> T {
    let mut s = T::ZERO;
    let base = i * n;
    let mut j = 0;
    while j < n {
        // SAFETY: guaranteed by the caller.
        unsafe {
            s += ld(m, base + j) * ld(b, j);
        }
        j += 1;
    }
    s
}

/// Conversion between precisions, round to nearest (the float64 → float32
/// step at the entry of the V-cycle and back at its exit).
///
/// # Safety
/// `i < src.len()`.
#[inline(always)]
pub unsafe fn convert<S: Real, T: Real>(i: usize, src: &[S]) -> T {
    // SAFETY: guaranteed by the caller.
    T::from_f64(unsafe { ld(src, i) }.to_f64())
}

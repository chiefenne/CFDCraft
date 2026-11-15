//! Transposed matrices for the adjoint (requirement D59, item 4): the
//! matrix of a finite-volume operator applied transposed, as a gather
//! without atomics (D31).
//!
//! The matrix structure of a cell-centred discretisation is symmetric:
//! interior face f couples owner P and neighbour N through two slots, a_PN
//! in row P and a_NP in row N (see [`crate::fv::SellLayout`]). The
//! *transpose map* gives, for every slot, the slot of the mirrored entry
//! (the diagonal and padding slots map to themselves). With it
//!
//! - [`transposed_value`] builds Aᵀ in the same storage as A (one gather per
//!   slot), so every solver and the multigrid work on Aᵀ unchanged;
//! - [`transposed_row_product`] computes (Aᵀx)_c without storing Aᵀ (one
//!   gather per row, the faces of the cell in stored order).
//!
//! Each result is written by one index and summed in a fixed order, so CPU
//! and GPU agree bitwise and the transposed product is reproducible like
//! the forward one.

use crate::fv::{SellLayout, Topology};
use crate::{ld, unpack_face};
use foehncraft_core::Real;

/// Slot `k` of Aᵀ in the storage of A: `val[tmap[k]]`.
///
/// # Safety
/// `k < tmap.len()`, every entry of `tmap` is `< val.len()`.
#[inline(always)]
pub unsafe fn transposed_value<T: Real>(k: usize, tmap: &[u32], val: &[T]) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { ld(val, ld(tmap, k) as usize) }
}

/// Row `c` of A x for the matrix `val` in `layout`: the diagonal times
/// x_c, then the off-diagonal slots (the interior faces of the cell in
/// stored order).
///
/// # Safety
/// `c` is a local cell of `topo`; `val` is the storage of `layout`; `x`
/// has an entry per cell.
#[inline(always)]
pub unsafe fn row_product<T: Real>(
    c: usize,
    topo: &Topology,
    layout: &SellLayout,
    val: &[T],
    x: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { row(c, topo, layout, &[], val, x, false) }
}

/// Row `c` of Aᵀ x, without storing Aᵀ: a_cc x_c + Σ_N a_Nc x_N, a_Nc read
/// from the slot of row N that the transpose map gives.
///
/// # Safety
/// As [`row_product`]; `tmap` is the transpose map of `layout`.
#[inline(always)]
pub unsafe fn transposed_row_product<T: Real>(
    c: usize,
    topo: &Topology,
    layout: &SellLayout,
    tmap: &[u32],
    val: &[T],
    x: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { row(c, topo, layout, tmap, val, x, true) }
}

#[inline(always)]
unsafe fn row<T: Real>(
    c: usize,
    topo: &Topology,
    layout: &SellLayout,
    tmap: &[u32],
    val: &[T],
    x: &[T],
    transposed: bool,
) -> T {
    // SAFETY: guaranteed by the callers.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut s = ld(val, layout.pos(c, 0)) * ld(x, c);
        let mut j = 1;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            if f < topo.n_interior {
                let o = if is_nb {
                    ld(topo.owner, f)
                } else {
                    ld(topo.neighbour, f)
                };
                let slot = layout.pos(c, j);
                let a = if transposed {
                    ld(val, ld(tmap, slot) as usize)
                } else {
                    ld(val, slot)
                };
                s += a * ld(x, o as usize);
                j += 1;
            }
            k += 1;
        }
        s
    }
}

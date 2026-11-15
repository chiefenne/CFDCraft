//! Face-based flux loop: the core pattern of a cell-centred finite-volume
//! assembly (D23).
//!
//! Pass 1, one thread per face: the face flux of a scalar φ with first-order
//! upwind convection and central diffusion,
//!
//! F_f = max(ṁ_f, 0) φ_P + min(ṁ_f, 0) φ_N − Γ_f (φ_N − φ_P),
//!
//! where P is the owner and N the neighbour cell, ṁ_f the mass flux (positive
//! from P to N) and Γ_f the diffusion coefficient times face area over
//! distance. Pass 2, one thread per cell: the net outflow, a gather over the
//! faces of the cell in a fixed order. Higher-order schemes (stage 2) add
//! gradient loads but keep the same access pattern.

use crate::{ld, unpack_face};
use inncraft_core::Real;

/// Flux through face `f`.
///
/// # Safety
/// `f < owner.len() == neighbour.len() == gamma.len() == mass_flux.len()`,
/// and every owner and neighbour index is `< phi.len()`.
#[inline(always)]
pub unsafe fn face_flux<T: Real>(
    f: usize,
    owner: &[u32],
    neighbour: &[u32],
    phi: &[T],
    gamma: &[T],
    mass_flux: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let p = ld(owner, f) as usize;
        let n = ld(neighbour, f) as usize;
        let phi_p = ld(phi, p);
        let phi_n = ld(phi, n);
        let m = ld(mass_flux, f);
        let g = ld(gamma, f);
        m.max_of(T::ZERO) * phi_p + m.min_of(T::ZERO) * phi_n - g * (phi_n - phi_p)
    }
}

/// Net outflow of cell `c`: the sum of the fluxes of its faces, with the sign
/// of the face orientation, in the stored order.
///
/// # Safety
/// `c + 1 < cell_ptr.len()`, the entries of `cell_ptr` are non-decreasing and
/// `<= cell_faces.len()`, and every unpacked face index is `< flux.len()`.
#[inline(always)]
pub unsafe fn cell_net_flux<T: Real>(
    c: usize,
    cell_ptr: &[u32],
    cell_faces: &[u32],
    flux: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(cell_ptr, c) as usize;
        let hi = ld(cell_ptr, c + 1) as usize;
        let mut s = T::ZERO;
        let mut k = lo;
        while k < hi {
            let (f, cell_is_neighbour) = unpack_face(ld(cell_faces, k));
            let v = ld(flux, f);
            if cell_is_neighbour {
                s -= v;
            } else {
                s += v;
            }
            k += 1;
        }
        s
    }
}

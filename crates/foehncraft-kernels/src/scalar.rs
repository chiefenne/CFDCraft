//! Kernel bodies of the generic scalar transport equation (D116), stage 6:
//!
//! ρ ∂φ/∂t + ∇·(ṁ φ) − φ ∇·ṁ − ∇·(ρΓ ∇φ) = ρ (S_u + S_p φ),  S_p ≤ 0,
//!
//! for passive scalars (D51) and, from PR 6b on, the turbulence variables.
//! The face coefficients are those of [`crate::fv::face_coefficients`]
//! (upwind with the linear-upwind or bounded-central deferred correction,
//! over-relaxed non-orthogonal diffusion); this module adds the row: the
//! bounded form, the point-implicit source and the time derivative.
//!
//! - **Bounded form** (stage 6 specification, 6.1): the term −φ_P Σ_f ṁ_f
//!   (the net mass outflow of the cell from the current fluxes) is
//!   implicit, so that the upwind convection operator has zero row sums
//!   and stays an M-matrix while continuity is not yet converged.
//! - **Sources** (6.2): S_u explicit, S_p ≤ 0 implicit in the diagonal
//!   (Patankar 1980); a positive S_p is taken as zero, so the diagonal can
//!   only grow.
//! - **Time** (D28): backward differences
//!   ∂φ/∂t ≈ (a₀ φⁿ⁺¹ + a₁ φⁿ + a₂ φⁿ⁻¹)/Δt. **Bounded BDF2** (D28 (1)):
//!   the backward difference is backward Euler from the old value
//!   φ* = −(a₁ φⁿ + a₂ φⁿ⁻¹)/a₀ with the step Δt/a₀, so it keeps the bounds
//!   of the scalar when φ* lies inside them. Where φ* leaves the bounds
//!   [lo, hi], the cell falls back to first order (a₀, a₁, a₂) = (1, −1, 0)
//!   for this step: a select per cell, no data-dependent branch beyond it.
//!
//! As everywhere in the crate the bodies are gathers in the stored face
//! order, so CPU loops and GPU threads give bitwise the same rows (D31).

use crate::fv::{SellLayout, Topology};
use crate::{ld, unpack_face};
use foehncraft_core::Real;

/// Inputs of [`scalar_row`].
#[derive(Clone, Copy)]
pub struct ScalarRowInputs<'a, T> {
    /// Connectivity.
    pub topo: Topology<'a>,
    /// Face coefficients c_P of [`crate::fv::face_coefficients`].
    pub cp: &'a [T],
    /// Face coefficients c_N.
    pub cn: &'a [T],
    /// Deferred parts b_f.
    pub bf: &'a [T],
    /// The convecting face mass flux (kg/s) of the bounded form.
    pub mass_flux: &'a [T],
    /// Cell volumes.
    pub volume: &'a [T],
    /// Current cell values (relaxation, old values of the source split).
    pub phi: &'a [T],
    /// Cell density ρ, one entry per cell.
    pub rho: &'a [T],
    /// Explicit source S_u per cell (the unit of φ per second); empty: none.
    pub su: &'a [T],
    /// Implicit source coefficient S_p ≤ 0 per cell (1/s); empty: none.
    pub sp: &'a [T],
    /// Under-relaxation factor α.
    pub alpha: T,
    /// Matrix layout.
    pub layout: SellLayout<'a>,
    /// True in a time step: the time derivative enters.
    pub transient: bool,
    /// 1/Δt.
    pub inv_dt: T,
    /// Coefficients a₀, a₁, a₂ of the backward difference.
    pub a: [T; 3],
    /// φⁿ.
    pub old: &'a [T],
    /// φⁿ⁻¹ (read when a₂ ≠ 0).
    pub old2: &'a [T],
    /// Lower bound of the bounded BDF2 (−∞: none).
    pub lo: T,
    /// Upper bound (+∞: none).
    pub hi: T,
    /// Weight m ∈ [0, 1] of the fixed value per cell (the ω of the
    /// automatic wall treatment of SST in the near-wall cells, stage 6
    /// specification 4.4): the row becomes (1 − m)·(transport row) +
    /// m a_P (φ_P − value), so m = 1 fixes the value and m = 0 leaves the
    /// row; empty: none.
    pub fix_mask: &'a [T],
    /// The fixed values (read where `fix_mask` is above 0).
    pub fix_val: &'a [T],
    /// A positive variable (k, ω): the sum D of the deferred corrections
    /// enters as S_u = max(D, 0) and S_p = min(D, 0)/φ* (stage 6
    /// specification 6.2), so that it cannot make the right-hand side
    /// negative.
    pub positive: bool,
    /// Pseudo-time coefficient ρ V/Δτ per cell (kg/s) of a steady
    /// iteration ([`pseudo_coefficient`]): the term ρ V (φ − φ^prev)/Δτ,
    /// zero at convergence; empty: none.
    pub pseudo: &'a [T],
}

/// Outputs of [`scalar_row`] for one cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScalarRow<T> {
    /// Relaxed diagonal a_P/α (also written into slot 0 of the row).
    pub diag: T,
    /// Right-hand side.
    pub rhs: T,
    /// 1 where the bounded BDF2 fell back to first order, else 0.
    pub first_order: T,
}

/// Row `c` of the scalar equation: the off-diagonals into `val` (SELL
/// slots 1.., as [`crate::fv::assemble_row`]), slot 0 the relaxed diagonal;
/// returns the diagonal, the right-hand side and the order flag.
///
/// # Safety
/// As [`crate::fv::assemble_row`]; `mass_flux` has one entry per face,
/// `volume`, `phi` (and `su`, `sp`, `old`, `old2` when used) one per cell.
#[inline(always)]
pub unsafe fn scalar_row<T: Real>(c: usize, inp: &ScalarRowInputs<T>, val: *mut T) -> ScalarRow<T> {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let topo = &inp.topo;
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut diag = T::ZERO;
        let mut b = T::ZERO;
        let mut dc = T::ZERO;
        let mut net = T::ZERO;
        let mut j = 1;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            let m = ld(inp.mass_flux, f);
            if is_nb {
                diag -= ld(inp.cn, f);
                dc += ld(inp.bf, f);
                net -= m;
                *val.add(inp.layout.pos(c, j)) = -ld(inp.cp, f);
                j += 1;
            } else {
                diag += ld(inp.cp, f);
                dc -= ld(inp.bf, f);
                net += m;
                if f < topo.n_interior {
                    *val.add(inp.layout.pos(c, j)) = ld(inp.cn, f);
                    j += 1;
                }
            }
            k += 1;
        }
        // bounded form: −φ_P Σ ±ṁ
        diag -= net;
        if inp.positive {
            let phi = ld(inp.phi, c).max_of(T::from_f64(1e-300));
            b += dc.max_of(T::ZERO);
            diag -= dc.min_of(T::ZERO) / phi;
        } else {
            b += dc;
        }
        let v = ld(inp.volume, c);
        let rv = ld(inp.rho, c) * v;
        if !inp.su.is_empty() {
            b += rv * ld(inp.su, c);
        }
        if !inp.sp.is_empty() {
            diag -= rv * ld(inp.sp, c).min_of(T::ZERO);
        }
        let mut first_order = T::ZERO;
        if inp.transient {
            let [a0, a1, a2] = inp.a;
            let p1 = ld(inp.old, c);
            let p2 = if a2 == T::ZERO {
                T::ZERO
            } else {
                ld(inp.old2, c)
            };
            let star = -(a1 * p1 + a2 * p2) / a0;
            let inside = star >= inp.lo && star <= inp.hi;
            let (a0, rest) = if inside {
                (a0, -(a1 * p1 + a2 * p2))
            } else {
                first_order = T::ONE;
                (T::ONE, p1)
            };
            let tau = rv * inp.inv_dt;
            diag += tau * a0;
            b += tau * rest;
        }
        if !inp.pseudo.is_empty() {
            // pseudo-time term: the old value is that of the last iteration
            let tau = ld(inp.pseudo, c);
            diag += tau;
            b += tau * ld(inp.phi, c);
        }
        let dr = diag / inp.alpha;
        b += (dr - diag) * ld(inp.phi, c);
        // a (partly) fixed value: (1 − m)·row + m a_P (φ_P − value), the
        // off-diagonals scaled by 1 − m (the row stays diagonally dominant;
        // m = 1: φ_P = value); the rows of the neighbours keep their coupling
        if !inp.fix_mask.is_empty() {
            let m = ld(inp.fix_mask, c);
            let keep = T::ONE - m;
            let mut i = 1;
            while i < j {
                let at = val.add(inp.layout.pos(c, i));
                *at *= keep;
                i += 1;
            }
            b = keep * b + m * dr * ld(inp.fix_val, c);
        }
        *val.add(inp.layout.pos(c, 0)) = dr;
        ScalarRow {
            diag: dr,
            rhs: b,
            first_order,
        }
    }
}

/// Parts of the residual scale of a scalar row, as those of the momentum
/// ([`crate::flow::residual_scale`]): (a_P φ_P, b − (a_P/α − a_P) φ_P) with
/// the unrelaxed a_P and b, so that the scale is finite also for φ = 0.
/// `pseudo` is the pseudo-time coefficient of the row (0: none), which is
/// taken out of both parts: the scale does not depend on the pseudo-time
/// step.
#[inline(always)]
pub fn residual_parts<T: Real>(phi: T, rhs: T, diag_relaxed: T, alpha: T, pseudo: T) -> (T, T) {
    let a = diag_relaxed * alpha;
    (
        (a - pseudo) * phi,
        rhs - (diag_relaxed - a) * phi - pseudo * phi,
    )
}

/// Pseudo-time coefficient ρ V/Δτ of a cell with the step
/// Δτ = min(Δτ_max, T/f): `rho_v` = ρ V, `inv_dt` = 1/Δτ_max and the rate
/// f/T given as `factor · num / den` (T = k/ε: num = ε, den = k; T = 1/(β*ω):
/// num = ω, den = 1; a stored T: num = 1, den = T). A cell with `den` ≤ 0
/// takes Δτ_max.
#[inline(always)]
pub fn pseudo_coefficient<T: Real>(rho_v: T, inv_dt: T, factor: T, num: T, den: T) -> T {
    if den > T::ZERO {
        rho_v * inv_dt.max_of(factor * num / den)
    } else {
        rho_v * inv_dt
    }
}

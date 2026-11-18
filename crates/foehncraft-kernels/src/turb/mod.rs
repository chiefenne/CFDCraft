//! Kernel bodies of the turbulence models (stage 6, D129–D139), generic
//! over [`foehncraft_core::Real`] so that the dual-number Jacobian (D59)
//! comes from the same source as the residual.
//!
//! - [`invariants`]: strain-rate and rotation tensors of the velocity
//!   gradient, their invariants, the Pope invariants and the other
//!   dimensionless features of D166;
//! - [`sst`]: the k-ω SST model, 2003 version (D130), per cell: eddy
//!   viscosity, blending, the point-implicit source split of the k and ω
//!   equations (section 6.2 of the stage 6 specification);
//! - [`wall`]: the low-Reynolds and automatic wall treatments of SST
//!   (D134) per wall face, and the rough-wall values of Aupoix (2015);
//!   the wall functions of the k-ε family and the compound treatment;
//! - [`thermal`]: the thermal wall functions of Kader and Jayatilleke and the
//!   turbulent Prandtl number of Kays & Crawford (D135);
//! - [`keps`]: the standard and the realizable k-ε model (D131);
//! - [`zeta_f`]: the k-ε-ζ-f model (D132), its time and length scales and
//!   the right-hand side of the elliptic f equation.
//!
//! Switches (min, max) go through [`smin`] and [`smax`]: the published
//! hard form for a width w = 0, the C¹ form of section 2.6 otherwise.

pub mod invariants;
pub mod keps;
pub mod sst;
pub mod thermal;
pub mod wall;
pub mod zeta_f;

use foehncraft_core::Real;

/// min(a, b): hard (w = 0, a select) or the smooth form
/// ½(a + b − √((a − b)² + (w m)²)), m = ½(|a| + |b|) (stage 6
/// specification 2.6; `foehncraft_core::math::smin`).
#[inline(always)]
pub fn smin<T: Real>(a: T, b: T, w: T) -> T {
    if w > T::ZERO {
        let half = T::from_f64(0.5);
        let m = half * (a.abs() + b.abs());
        let d = a - b;
        let e = w * m;
        half * (a + b - (d * d + e * e).sqrt_of())
    } else {
        a.min_of(b)
    }
}

/// max(a, b), hard or smooth (see [`smin`]).
#[inline(always)]
pub fn smax<T: Real>(a: T, b: T, w: T) -> T {
    -smin(-a, -b, w)
}

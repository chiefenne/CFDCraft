//! The k-ω SST model, 2003 version (D130; stage 6 specification 3.1, 6.2):
//! F. R. Menter, M. Kuntz, R. Langtry, *Ten years of industrial experience
//! with the SST turbulence model*, Turbulence, Heat and Mass Transfer 4,
//! 625–632 (2003), with the constants of F. R. Menter, *Two-equation
//! eddy-viscosity turbulence models for engineering applications*, AIAA J.
//! 32, 1598–1605 (1994) where the 2003 paper keeps them.
//!
//! Per unit mass, constant density:
//!
//!   Dk/Dt = P̃_k − β* k ω + ∇·[(ν + σ_k ν_t)∇k]
//!
//!   Dω/Dt = (γ/ν_t) P̃_k − β ω² + ∇·[(ν + σ_ω ν_t)∇ω]
//!           + 2(1 − F₁) σ_ω2 (1/ω) ∇k·∇ω
//!
//! ν_t = a₁ k / max(a₁ ω, S F₂), S = (2 S_ij S_ij)^½; P_k = ν_t S²,
//! P̃_k = min(P_k, 10 β* k ω) in both equations (the 2003 limiter; the ω
//! equation takes the limited P̃_k, the correction of the 2003 paper
//! recorded on the NASA Turbulence Modeling Resource page "The Menter
//! Shear Stress Transport Turbulence Model", section SST-2003, retrieved
//! 2026-11-17).
//!
//! Constants checked on 2026-11-17 against the TMR SST page (which
//! reproduces Menter 1994, table of constants, and the SST-2003 changes):
//! β* = 0.09, a₁ = 0.31, κ = 0.41; set 1 σ_k1 = 0.85, σ_ω1 = 0.5,
//! β₁ = 0.075, γ₁ = 5/9 (2003: set directly); set 2 σ_k2 = 1.0,
//! σ_ω2 = 0.856, β₂ = 0.0828, γ₂ = 0.44 (2003: set directly); the floor
//! of CD_kω 10⁻¹⁰ (2003; 10⁻²⁰ in 1994). The TMR writes CD_kω with ρ; the
//! per-unit-mass form here has the floor in 1/s² (ρ = 1 kg/m³ there for
//! the comparison; the floor only acts in the free stream).
//!
//! The source split (Patankar 1980; specification 6.2) keeps S_u ≥ 0 and
//! S_p ≤ 0 in every cell, with select operations only, so the matrices
//! stay M-matrices and k, ω stay positive.

use super::invariants::{Tensor, s2_omega2};
use super::{smax, smin};
use foehncraft_core::Real;

/// Production form: the 2003 limiter P̃_k = min(ν_t S², 10 β* k ω).
pub const LIMITER_MENTER: u32 = 0;
/// Kato & Launder (1993): P_k = ν_t S Ω, then the 2003 limiter.
pub const LIMITER_KATO_LAUNDER: u32 = 1;
/// Durbin (1996): ν_t ≤ C_lim k/(√6 |S|) (realizability), then the 2003
/// limiter.
pub const LIMITER_DURBIN: u32 = 2;
/// No limiter: P_k = ν_t S² (the stagnation-point anomaly is reproduced).
pub const LIMITER_NONE: u32 = 3;

/// The constants of the model (the paper's by default; a case may change
/// them, TURB007).
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SstConstants {
    /// β* = 0.09 (Menter 1994).
    pub beta_star: f64,
    /// a₁ = 0.31 (Menter 1994, Bradshaw's structural constant).
    pub a1: f64,
    /// κ = 0.41.
    pub kappa: f64,
    /// γ₁ = 5/9 (Menter, Kuntz & Langtry 2003).
    pub gamma1: f64,
    /// γ₂ = 0.44 (2003).
    pub gamma2: f64,
    /// β₁ = 0.075.
    pub beta1: f64,
    /// β₂ = 0.0828.
    pub beta2: f64,
    /// σ_k1 = 0.85.
    pub sigma_k1: f64,
    /// σ_k2 = 1.0.
    pub sigma_k2: f64,
    /// σ_ω1 = 0.5.
    pub sigma_w1: f64,
    /// σ_ω2 = 0.856.
    pub sigma_w2: f64,
    /// Factor of the production limiter, 10 (2003; 20 in 1994).
    pub limiter_factor: f64,
    /// Floor of CD_kω, 10⁻¹⁰ 1/s² (2003).
    pub cd_floor: f64,
}

impl Default for SstConstants {
    fn default() -> Self {
        Self {
            beta_star: 0.09,
            a1: 0.31,
            kappa: 0.41,
            gamma1: 5.0 / 9.0,
            gamma2: 0.44,
            beta1: 0.075,
            beta2: 0.0828,
            sigma_k1: 0.85,
            sigma_k2: 1.0,
            sigma_w1: 0.5,
            sigma_w2: 0.856,
            limiter_factor: 10.0,
            cd_floor: 1e-10,
        }
    }
}

/// Run-time options of the model.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SstOptions {
    /// Production form ([`LIMITER_MENTER`], …).
    pub limiter: u32,
    /// 1: Newton linearisation of −β ω² (S_p = −2βω*), 0: Picard (−βω*).
    pub newton: u32,
    /// C_lim of the Durbin bound.
    pub c_lim: f64,
    /// Relative width of the smooth switches (0: the published hard form).
    pub smooth: f64,
    /// Ambient k of the sustaining terms (Spalart & Rumsey 2007); 0: off.
    pub k_amb: f64,
    /// Ambient ω of the sustaining terms; 0: off.
    pub w_amb: f64,
    /// Largest ν_t/ν (`turbulence.max_viscosity_ratio`).
    pub max_ratio: f64,
}

impl Default for SstOptions {
    fn default() -> Self {
        Self {
            limiter: LIMITER_MENTER,
            newton: 1,
            c_lim: 1.0,
            smooth: 0.0,
            k_amb: 0.0,
            w_amb: 0.0,
            max_ratio: 1e5,
        }
    }
}

/// Inputs of [`sst_cell`] for one cell.
#[derive(Clone, Copy, Debug)]
pub struct SstCellInputs<T> {
    /// k (m²/s²).
    pub k: T,
    /// ω (1/s).
    pub w: T,
    /// ∇k.
    pub grad_k: [T; 3],
    /// ∇ω.
    pub grad_w: [T; 3],
    /// Velocity gradient g_ij = ∂u_i/∂x_j.
    pub grad_u: Tensor<T>,
    /// Wall distance d (m; +∞ without walls).
    pub d: T,
    /// ν (m²/s).
    pub nu: T,
}

/// Results of [`sst_cell`] for one cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SstCell<T> {
    /// ν_t (m²/s).
    pub nu_t: T,
    /// F₁.
    pub f1: T,
    /// The production entering the k equation (P̃_k, m²/s³).
    pub pk: T,
    /// S_u of the k equation (m²/s³).
    pub su_k: T,
    /// S_p ≤ 0 of the k equation (1/s).
    pub sp_k: T,
    /// S_u of the ω equation (1/s²).
    pub su_w: T,
    /// S_p ≤ 0 of the ω equation (1/s).
    pub sp_w: T,
    /// σ_k ν_t (diffusivity of k without ν).
    pub dk: T,
    /// σ_ω ν_t.
    pub dw: T,
    /// 1 where the cap `max_ratio` of ν_t/ν acted, else 0.
    pub capped: T,
}

/// ν_t, F₁ and the source split of the k and ω equations in one cell,
/// with the values k*, ω* of the previous outer iteration.
#[inline(always)]
pub fn sst_cell<T: Real>(i: &SstCellInputs<T>, c: &SstConstants, o: &SstOptions) -> SstCell<T> {
    let f = T::from_f64;
    let w_s = f(o.smooth);
    let (k, w, d, nu) = (i.k, i.w, i.d, i.nu);
    let bs = f(c.beta_star);
    let a1 = f(c.a1);
    let (s2, o2) = s2_omega2(&i.grad_u);
    let s = s2.sqrt_of();
    let sk = k.max_of(T::ZERO).sqrt_of();
    let dkdw = i.grad_k[0] * i.grad_w[0] + i.grad_k[1] * i.grad_w[1] + i.grad_k[2] * i.grad_w[2];
    let two = f(2.0);
    let sw2 = f(c.sigma_w2);
    // CD_kω = max(2 σ_ω2 (1/ω) ∇k·∇ω, 10⁻¹⁰)
    let cross = two * sw2 * dkdw / w;
    let cd = smax(cross, f(c.cd_floor), w_s);
    let d2 = d * d;
    // arg₁ = min[max(√k/(β* ω d), 500 ν/(d² ω)), 4 σ_ω2 k/(CD_kω d²)]
    let t1 = sk / (bs * w * d);
    let t2 = f(500.0) * nu / (d2 * w);
    let arg1 = smin(smax(t1, t2, w_s), f(4.0) * sw2 * k / (cd * d2), w_s);
    let a4 = arg1 * arg1 * arg1 * arg1;
    let f1 = a4.tanh_of();
    // arg₂ = max(2√k/(β* ω d), 500 ν/(d² ω))
    let arg2 = smax(two * t1, t2, w_s);
    let f2 = (arg2 * arg2).tanh_of();
    // ν_t = k / r with r = max(a₁ω, S F₂)/a₁ (and the Durbin bound and
    // the cap of ν_t/ν as further maxima of r), so that no division by k
    // is needed: P̃_k/ν_t = min(S², 10 β* ω r)
    let mut r = smax(a1 * w, s * f2, w_s) / a1;
    if o.limiter == LIMITER_DURBIN {
        // ν_t ≤ C_lim k/(√6 |S|), |S| = S/√2: k/ν_t ≥ √3 S / C_lim
        r = smax(r, f(foehncraft_core::math::sqrt(3.0) / o.c_lim) * s, w_s);
    }
    let r_cap = k / (f(o.max_ratio) * nu);
    let capped = if r_cap > r { T::ONE } else { T::ZERO };
    r = r.max_of(r_cap);
    let nu_t = k.max_of(T::ZERO) / r;
    // production per ν_t
    let p_over = if o.limiter == LIMITER_KATO_LAUNDER {
        s * o2.sqrt_of()
    } else {
        s2
    };
    let p_over = if o.limiter == LIMITER_NONE {
        p_over
    } else {
        smin(p_over, f(c.limiter_factor) * bs * w * r, w_s)
    };
    let pk = nu_t * p_over;
    let one_f1 = T::ONE - f1;
    let blend = |a: f64, b: f64| f1 * f(a) + one_f1 * f(b);
    let gamma = blend(c.gamma1, c.gamma2);
    let beta = blend(c.beta1, c.beta2);
    let sig_k = blend(c.sigma_k1, c.sigma_k2);
    let sig_w = blend(c.sigma_w1, c.sigma_w2);
    // k: S_u = P̃_k (+ β* ω_amb k_amb), S_p = −β* ω*
    let amb_k = f(c.beta_star * o.w_amb * o.k_amb);
    let su_k = pk + amb_k;
    let sp_k = -bs * w;
    // ω: (γ/ν_t) P̃_k = γ P̃_k/ν_t; −β ω² (Newton or Picard); the cross
    // diffusion split by its sign; the sustaining term β ω_amb²
    let cdt = one_f1 * cross;
    let amb_w = beta * f(o.w_amb * o.w_amb);
    let (su_d, sp_d) = if o.newton != 0 {
        (beta * w * w, -two * beta * w)
    } else {
        (T::ZERO, -beta * w)
    };
    let su_w = gamma * p_over + su_d + cdt.max_of(T::ZERO) + amb_w;
    let sp_w = sp_d + cdt.min_of(T::ZERO) / w;
    SstCell {
        nu_t,
        f1,
        pk,
        su_k,
        sp_k,
        su_w,
        sp_w,
        dk: sig_k * nu_t,
        dw: sig_w * nu_t,
        capped,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use foehncraft_core::Dual;

    fn base() -> SstCellInputs<f64> {
        let mut g = [[0.0; 3]; 3];
        g[0][1] = 50.0;
        SstCellInputs {
            k: 0.01,
            w: 100.0,
            grad_k: [0.0, -0.3, 0.0],
            grad_w: [0.0, -40.0, 0.0],
            grad_u: g,
            d: 0.01,
            nu: 1.5e-5,
        }
    }

    #[test]
    fn eddy_viscosity_and_production_in_equilibrium_shear() {
        // far from walls (F₁ = F₂ = 0): ν_t = k/ω where a₁ω > S F₂
        let mut i = base();
        i.d = f64::INFINITY;
        i.grad_k = [0.0; 3];
        let c = SstConstants::default();
        let r = sst_cell(&i, &c, &SstOptions::default());
        assert_eq!(r.f1, 0.0);
        assert!((r.nu_t - 1e-4).abs() < 1e-18);
        // P = ν_t S² = 1e-4 · 2500 = 0.25 > 10 β* k ω = 0.9? no: 0.25 < 0.9
        assert!((r.pk - 0.25).abs() < 1e-15);
        assert!((r.sp_k + 0.09 * 100.0).abs() < 1e-12);
        // ω: γ₂ S² + β₂ ω² (Newton), S_p = −2β₂ω
        assert!((r.su_w - (0.44 * 2500.0 + 0.0828 * 1e4)).abs() < 1e-9);
        assert!((r.sp_w + 2.0 * 0.0828 * 100.0).abs() < 1e-12);
        assert!((r.dk - 1e-4).abs() < 1e-18 && (r.dw - 0.856e-4).abs() < 1e-18);
    }

    #[test]
    fn limiter_caps_the_production() {
        let mut i = base();
        i.d = f64::INFINITY;
        i.grad_u[0][1] = 1000.0;
        let c = SstConstants::default();
        let r = sst_cell(&i, &c, &SstOptions::default());
        // ν_t = a₁k/max(a₁ω, S F₂) = k/ω (F₂ = 0); P = 1e-4·1e6 = 100 >
        // 10 β* k ω = 0.9
        assert!((r.pk - 0.9).abs() < 1e-12, "{}", r.pk);
        let none = SstOptions {
            limiter: LIMITER_NONE,
            ..Default::default()
        };
        assert!((sst_cell(&i, &c, &none).pk - 100.0).abs() < 1e-9);
    }

    #[test]
    fn kato_launder_removes_production_without_vorticity() {
        // plane strain (stagnation flow): u = a x, v = −a y
        let mut i = base();
        i.d = f64::INFINITY;
        i.grad_u = [[10.0, 0.0, 0.0], [0.0, -10.0, 0.0], [0.0; 3]];
        let c = SstConstants::default();
        let kl = SstOptions {
            limiter: LIMITER_KATO_LAUNDER,
            ..Default::default()
        };
        assert_eq!(sst_cell(&i, &c, &kl).pk, 0.0);
        assert!(sst_cell(&i, &c, &SstOptions::default()).pk > 0.0);
    }

    #[test]
    fn durbin_bounds_the_eddy_viscosity() {
        let mut i = base();
        i.d = f64::INFINITY;
        i.w = 1.0;
        i.grad_u = [[100.0, 0.0, 0.0], [0.0, -100.0, 0.0], [0.0; 3]];
        let c = SstConstants::default();
        let o = SstOptions {
            limiter: LIMITER_DURBIN,
            ..Default::default()
        };
        let r = sst_cell(&i, &c, &o);
        // S = √(2·(100² + 100²)) = 200, |S| = S/√2
        let bound = i.k / (6f64.sqrt() * 200.0 / 2f64.sqrt());
        assert!((r.nu_t - bound).abs() < 1e-15 * bound.max(1.0));
        assert!(r.nu_t < i.k / i.w);
    }

    #[test]
    fn near_wall_blending_is_one() {
        let mut i = base();
        i.d = 1e-3;
        let c = SstConstants::default();
        let r = sst_cell(&i, &c, &SstOptions::default());
        assert!(r.f1 > 0.99, "{}", r.f1);
        assert!(r.su_k >= 0.0 && r.sp_k <= 0.0 && r.su_w >= 0.0 && r.sp_w <= 0.0);
    }

    #[test]
    fn sustaining_terms_balance_the_ambient_decay() {
        // at k = k_amb, ω = ω_amb and no shear: S_u + S_p φ = 0 in both
        let mut i = base();
        i.d = f64::INFINITY;
        i.grad_u = [[0.0; 3]; 3];
        i.grad_k = [0.0; 3];
        i.k = 1e-6;
        i.w = 5.0;
        let c = SstConstants::default();
        let o = SstOptions {
            k_amb: 1e-6,
            w_amb: 5.0,
            ..Default::default()
        };
        let r = sst_cell(&i, &c, &o);
        assert!((r.su_k + r.sp_k * i.k).abs() < 1e-20);
        assert!((r.su_w + r.sp_w * i.w).abs() < 1e-12);
    }

    #[test]
    fn smooth_switches_stay_close_to_the_hard_ones() {
        let i = base();
        let c = SstConstants::default();
        let hard = sst_cell(&i, &c, &SstOptions::default());
        let smooth = sst_cell(
            &i,
            &c,
            &SstOptions {
                smooth: 1e-3,
                ..Default::default()
            },
        );
        assert!((hard.nu_t - smooth.nu_t).abs() < 1e-3 * hard.nu_t);
        assert!((hard.f1 - smooth.f1).abs() < 1e-3);
    }

    #[test]
    fn dual_numbers_give_the_derivative_of_the_sources() {
        // d(su_w)/dω by dual numbers against central differences, smooth
        // switches (the linearisation of D59)
        let i = base();
        let c = SstConstants::default();
        let o = SstOptions {
            smooth: 1e-3,
            ..Default::default()
        };
        let lift = |x: f64| Dual::<1>::constant(x);
        let di = SstCellInputs {
            k: lift(i.k),
            w: Dual::<1>::variable(i.w, 0),
            grad_k: i.grad_k.map(lift),
            grad_w: i.grad_w.map(lift),
            grad_u: i.grad_u.map(|r| r.map(lift)),
            d: lift(i.d),
            nu: lift(i.nu),
        };
        let r = sst_cell(&di, &c, &o);
        let h = 1e-4;
        let at = |w: f64| {
            let mut j = i;
            j.w = w;
            sst_cell(&j, &c, &o)
        };
        let (p, m) = (at(i.w + h), at(i.w - h));
        for (dv, fd) in [
            (r.su_w.d[0], (p.su_w - m.su_w) / (2.0 * h)),
            (r.nu_t.d[0], (p.nu_t - m.nu_t) / (2.0 * h)),
            (r.f1.d[0], (p.f1 - m.f1) / (2.0 * h)),
        ] {
            assert!((dv - fd).abs() <= 1e-6 * fd.abs().max(1e-8), "{dv} {fd}");
        }
    }
}

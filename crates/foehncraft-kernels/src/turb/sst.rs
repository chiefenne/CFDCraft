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
    /// c_r1 = 1.0 of the curvature correction (Smirnov & Menter 2009).
    pub c_r1: f64,
    /// c_r2 = 2.0.
    pub c_r2: f64,
    /// c_r3 = 1.0.
    pub c_r3: f64,
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
            c_r1: 1.0,
            c_r2: 2.0,
            c_r3: 1.0,
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
    /// The factor f_r1 of the production of k and ω (the curvature
    /// correction, [`rotation_function`]; 1 without it).
    pub fr1: T,
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
    // the curvature correction multiplies the production of both
    // equations (Smirnov & Menter 2009)
    let p_over = p_over * i.fr1;
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

/// ε_ijk.
#[inline(always)]
fn levi(i: usize, j: usize, k: usize) -> f64 {
    if i == j || j == k || i == k {
        0.0
    } else if (i + 1) % 3 == j {
        1.0
    } else {
        -1.0
    }
}

/// The rotation function f_r1 of the curvature and rotation correction of
/// Spalart & Shur (1997) in the SST form of Smirnov & Menter (2009)
/// (specification 3.6, D170):
///
///   f_r1 = max[0, min(f_rotation, 1.25)],
///   f_rotation = (1 + c_r1) 2r*/(1 + r*) [1 − c_r3 atan(c_r2 r̃)] − c_r1,
///   r* = S/Ω,
///   r̃ = 2 Ω_ik S_jk [DS_ij/Dt + (ε_imn S_jn + ε_jmn S_in) Ω_m^rot] / (Ω D³),
///   D² = max(S², 0.09 ω²),
///
/// S = (2 S_ij S_ij)^½, Ω = (2 Ω_ij Ω_ij)^½. Ω_ij is the rotation tensor of
/// the **absolute** velocity: the gradient `g` of the absolute velocity
/// gives it directly (a moving reference frame, D55), a velocity relative
/// to a frame turning at `spin` gives Ω_ij = W_ij − ε_ijk spin_k (the
/// absolute vorticity is the relative one plus 2 spin). With that sign the
/// side of a channel in spanwise rotation whose mean vorticity opposes the
/// rotation gets f_r1 > 1 from both r* and r̃ (the unstable side).
/// `dsdt` holds DS_ij/Dt in the order xx, yy, zz, xy, yz, zx: in the frame
/// that turns at `omega_rot`, the convective derivative with the velocity
/// relative to that frame plus the time derivative. `smooth` > 0 takes the
/// smooth limits of section 2.6 for the clip.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub fn rotation_function<T: Real>(
    g: &Tensor<T>,
    dsdt: &[T; 6],
    omega_rot: [T; 3],
    spin: [T; 3],
    w: T,
    c: &SstConstants,
    smooth: f64,
) -> T {
    let f = T::from_f64;
    let (s, mut wt) = super::invariants::strain_rotation(g);
    let mut i = 0;
    while i < 3 {
        let mut j = 0;
        while j < 3 {
            let mut k = 0;
            while k < 3 {
                let e = levi(i, j, k);
                if e != 0.0 {
                    wt[i][j] = wt[i][j] - f(e) * spin[k];
                }
                k += 1;
            }
            j += 1;
        }
        i += 1;
    }
    let two = f(2.0);
    let s2 = two * super::invariants::contract(&s, &s);
    let o2 = two * super::invariants::contract(&wt, &wt);
    let tiny = f(1e-300);
    let sm = s2.sqrt_of();
    let om = o2.sqrt_of().max_of(tiny);
    let r_star = sm / om;
    let d2 = s2.max_of(f(0.09) * w * w);
    let d3 = (d2 * d2.sqrt_of()).max_of(tiny);
    let ds = [
        [dsdt[0], dsdt[3], dsdt[5]],
        [dsdt[3], dsdt[1], dsdt[4]],
        [dsdt[5], dsdt[4], dsdt[2]],
    ];
    // Σ_ij A_ij T_ij with A_ij = Ω_ik S_jk and T_ij = DS_ij/Dt +
    // (ε_imn S_jn + ε_jmn S_in) Ω_m^rot
    let mut num = T::ZERO;
    let mut i = 0;
    while i < 3 {
        let mut j = 0;
        while j < 3 {
            let mut a = T::ZERO;
            let mut k = 0;
            while k < 3 {
                a = a + wt[i][k] * s[j][k];
                k += 1;
            }
            let mut t = ds[i][j];
            let mut m = 0;
            while m < 3 {
                let mut n = 0;
                while n < 3 {
                    let e1 = levi(i, m, n);
                    let e2 = levi(j, m, n);
                    if e1 != 0.0 {
                        t = t + f(e1) * s[j][n] * omega_rot[m];
                    }
                    if e2 != 0.0 {
                        t = t + f(e2) * s[i][n] * omega_rot[m];
                    }
                    n += 1;
                }
                m += 1;
            }
            num = num + a * t;
            j += 1;
        }
        i += 1;
    }
    let r_t = two * num / (om * d3);
    let fr = (T::ONE + f(c.c_r1)) * two * r_star / (T::ONE + r_star)
        * (T::ONE - f(c.c_r3) * (f(c.c_r2) * r_t).atan_of())
        - f(c.c_r1);
    let w_s = f(smooth);
    smax(smin(fr, f(1.25), w_s), T::ZERO, w_s)
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
            fr1: 1.0,
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
            fr1: lift(1.0),
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

    fn shear(gamma: f64) -> Tensor<f64> {
        let mut g = [[0.0; 3]; 3];
        g[0][1] = gamma;
        g
    }

    #[test]
    fn rotation_function_is_one_in_pure_shear() {
        // r* = 1 and r̃ = 0: SST exactly (the first check of 9.8)
        let c = SstConstants::default();
        let f = rotation_function(&shear(50.0), &[0.0; 6], [0.0; 3], [0.0; 3], 100.0, &c, 0.0);
        assert!((f - 1.0).abs() < 1e-14, "{f}");
        let f = rotation_function(&shear(-7.0), &[0.0; 6], [0.0; 3], [0.0; 3], 3.0, &c, 0.0);
        assert!((f - 1.0).abs() < 1e-14, "{f}");
    }

    #[test]
    fn spanwise_rotation_raises_the_unstable_side() {
        // channel in spanwise rotation Ω_z > 0, relative velocity u(y):
        // the lower wall (dU/dy > 0, mean vorticity −dU/dy opposes Ω) is
        // the unstable side, the upper wall the stable one
        let c = SstConstants::default();
        let rot = [0.0, 0.0, 2.0];
        let lower = rotation_function(&shear(50.0), &[0.0; 6], rot, rot, 100.0, &c, 0.0);
        let upper = rotation_function(&shear(-50.0), &[0.0; 6], rot, rot, 100.0, &c, 0.0);
        assert!(lower > 1.0 && lower <= 1.25, "{lower}");
        assert!(upper < 1.0 && upper >= 0.0, "{upper}");
        // r* by hand: S = 50, Ω = |50 − 4| = 46 on the lower side
        let rs: f64 = 50.0 / 46.0;
        // r̃ = −2 (γ/2 − Ω_z) γ² Ω_z/(Ω D³), D = max(S, 0.3 ω) = 50
        let rt = -2.0 * (25.0 - 2.0) * 2500.0 * 2.0 / (46.0 * 50.0f64.powi(3));
        let want =
            2.0 * 2.0 * rs / (1.0 + rs) * (1.0 - foehncraft_core::math::atan(2.0 * rt)) - 1.0;
        assert!((lower - want.min(1.25)).abs() < 1e-12, "{lower} {want}");
    }

    #[test]
    fn absolute_and_relative_frames_agree() {
        // the same flow seen in the inertial frame (absolute gradient with
        // the solid-body part, spin 0) and in the turning frame (relative
        // gradient, spin Ω) gives the same rotation tensor and f_r1
        let c = SstConstants::default();
        let om = 3.0;
        let rel = shear(40.0);
        let mut abs = rel;
        // u_abs = u_rel + Ω × r: ∂u_x/∂y −= Ω, ∂u_y/∂x += Ω
        abs[0][1] -= om;
        abs[1][0] += om;
        let ds = [0.3, -0.2, 0.0, 1.5, 0.0, 0.0];
        let rot = [0.0, 0.0, om];
        let a = rotation_function(&abs, &ds, rot, [0.0; 3], 80.0, &c, 0.0);
        let b = rotation_function(&rel, &ds, rot, rot, 80.0, &c, 0.0);
        assert!((a - b).abs() < 1e-14, "{a} {b}");
    }

    #[test]
    fn rotation_function_scales_the_production() {
        let i = base();
        let c = SstConstants::default();
        let o = SstOptions::default();
        let r1 = sst_cell(&i, &c, &o);
        let mut j = i;
        j.fr1 = 0.5;
        let r2 = sst_cell(&j, &c, &o);
        assert!((r2.pk - 0.5 * r1.pk).abs() < 1e-15 * r1.pk.abs().max(1.0));
        assert_eq!(r2.nu_t, r1.nu_t);
        assert!(r2.su_w < r1.su_w);
    }
}

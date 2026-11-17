//! The standard and the realizable k-ε model (D131; stage 6 specification
//! 3.2, 6.2): one implementation with two closures (D12).
//!
//! Per unit mass, constant density:
//!
//!   Dk/Dt = P_k − ε + ∇·[(ν + ν_t/σ_k)∇k]
//!
//! **Standard** (B. E. Launder, D. B. Spalding, *The numerical computation
//! of turbulent flows*, Comput. Methods Appl. Mech. Eng. 3, 269–289
//! (1974)): ν_t = C_μ k²/ε,
//!
//!   Dε/Dt = (ε/k)(C_ε1 P_k − C_ε2 ε) + ∇·[(ν + ν_t/σ_ε)∇ε],
//!
//! C_μ = 0.09, C_ε1 = 1.44, C_ε2 = 1.92, σ_k = 1.0, σ_ε = 1.3 (the
//! constants of the 1974 paper, its table; checked on 2026-11-17 against
//! their reproduction in S. B. Pope, *Turbulent Flows*, CUP 2000, eq.
//! 10.53, which cites the paper). Production options: Kato & Launder
//! (1993) P_k = ν_t S Ω; the realizability bound of Durbin (1996)
//! ν_t ≤ C_lim k/(√6 |S|); and the clip P̃_k = min(P_k, 10 ε), the
//! analogue of the SST limiter (10 β* k ω = 10 ε).
//!
//! **Realizable** (T.-H. Shih, W. W. Liou, A. Shabbir, Z. Yang, J. Zhu,
//! *A new k-ε eddy viscosity model for high Reynolds number turbulent
//! flows*, Computers & Fluids 24, 227–238 (1995); report version NASA
//! TM-106721, ICOMP-94-21, CMOTT-94-6 (1994)):
//!
//!   Dε/Dt = C₁ S ε − C₂ ε²/(k + √(νε)) + ∇·[(ν + ν_t/σ_ε)∇ε],
//!
//! C₁ = max(0.43, η/(η + 5)), η = S k/ε, S = (2 S_ij S_ij)^½, C₂ = 1.9,
//! σ_k = 1.0, σ_ε = 1.2; ν_t = C_μ k²/ε with
//! C_μ = 1/(A₀ + A_s U* k/ε), U* = (S_ij S_ij + W_ij W_ij)^½,
//! A_s = √6 cos φ, φ = ⅓ arccos(√6 W̃), W̃ = S_ij S_jk S_ki/S̃³,
//! S̃ = (S_ij S_ij)^½, A₀ = 4.0. W̃ is clipped to [−1/√6, 1/√6] (its
//! range for a trace-free S_ij; round-off can leave it). Constants checked
//! on 2026-11-17 against NASA TM-106721 (table 2: A₀ = 4.0, C₂ = 1.9,
//! σ_k = 1.0, σ_ε = 1.2; eqs. 19–21 and 25). In a rotating zone the paper
//! takes W_ij relative to the frame; the solved velocity is the absolute
//! one (D55), whose rotation tensor is the paper's Ω̃_ij.
//!
//! The source split (Patankar 1980; specification 6.2) keeps S_u ≥ 0 and
//! S_p ≤ 0 in every cell, with select operations only.

use super::invariants::{Tensor, s2_omega2, strain_rotation};
use super::{smax, smin};
use foehncraft_core::Real;

/// The standard model (Launder & Spalding 1974).
pub const STANDARD: u32 = 0;
/// The realizable model (Shih et al. 1995).
pub const REALIZABLE: u32 = 1;

/// Production form: P_k = ν_t S² (the default of the k-ε family).
pub const LIMITER_NONE: u32 = 0;
/// Kato & Launder (1993): P_k = ν_t S Ω.
pub const LIMITER_KATO_LAUNDER: u32 = 1;
/// Durbin (1996): ν_t ≤ C_lim k/(√6 |S|).
pub const LIMITER_DURBIN: u32 = 2;
/// The clip P̃_k = min(ν_t S², 10 ε).
pub const LIMITER_CLIP: u32 = 3;

/// The constants of the model (the paper's by default; a case may change
/// them, TURB007).
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KepsConstants {
    /// C_μ = 0.09 (standard; the realizable model computes its own and
    /// takes this one in the wall functions and at inlets).
    pub c_mu: f64,
    /// C_ε1 = 1.44 (standard).
    pub c_e1: f64,
    /// C_ε2 = 1.92 (standard), C₂ = 1.9 (realizable).
    pub c_e2: f64,
    /// σ_k = 1.0.
    pub sigma_k: f64,
    /// σ_ε = 1.3 (standard), 1.2 (realizable).
    pub sigma_e: f64,
    /// A₀ = 4.0 (realizable).
    pub a0: f64,
    /// Lower bound of C₁, 0.43 (realizable).
    pub c1_min: f64,
    /// κ = 0.41 (wall functions).
    pub kappa: f64,
}

impl KepsConstants {
    /// The constants of the standard model.
    pub const fn standard() -> Self {
        Self {
            c_mu: 0.09,
            c_e1: 1.44,
            c_e2: 1.92,
            sigma_k: 1.0,
            sigma_e: 1.3,
            a0: 4.0,
            c1_min: 0.43,
            kappa: 0.41,
        }
    }

    /// The constants of the realizable model.
    pub const fn realizable() -> Self {
        Self {
            c_e2: 1.9,
            sigma_e: 1.2,
            ..Self::standard()
        }
    }
}

impl Default for KepsConstants {
    fn default() -> Self {
        Self::standard()
    }
}

/// Run-time options of the model.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KepsOptions {
    /// [`STANDARD`] or [`REALIZABLE`].
    pub variant: u32,
    /// Production form ([`LIMITER_NONE`], …; standard only).
    pub limiter: u32,
    /// 1: Newton linearisation of the destruction of ε
    /// (S_p = −2 C ε*/k*), 0: Picard (−C ε*/k*).
    pub newton: u32,
    /// C_lim of the Durbin bound.
    pub c_lim: f64,
    /// Relative width of the smooth switches (0: the published hard form).
    pub smooth: f64,
    /// Ambient k of the sustaining terms (Spalart & Rumsey 2007); 0: off.
    pub k_amb: f64,
    /// Ambient ε of the sustaining terms; 0: off.
    pub e_amb: f64,
    /// Largest ν_t/ν (`turbulence.max_viscosity_ratio`).
    pub max_ratio: f64,
}

impl Default for KepsOptions {
    fn default() -> Self {
        Self {
            variant: STANDARD,
            limiter: LIMITER_NONE,
            newton: 1,
            c_lim: 1.0,
            smooth: 0.0,
            k_amb: 0.0,
            e_amb: 0.0,
            max_ratio: 1e5,
        }
    }
}

/// Inputs of [`keps_cell`] for one cell.
#[derive(Clone, Copy, Debug)]
pub struct KepsCellInputs<T> {
    /// k (m²/s²).
    pub k: T,
    /// ε (m²/s³).
    pub e: T,
    /// Velocity gradient g_ij = ∂u_i/∂x_j.
    pub grad_u: Tensor<T>,
    /// ν (m²/s).
    pub nu: T,
}

/// Results of [`keps_cell`] for one cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KepsCell<T> {
    /// ν_t (m²/s).
    pub nu_t: T,
    /// C_μ (constant for the standard model).
    pub c_mu: T,
    /// The production entering the k equation (m²/s³).
    pub pk: T,
    /// S_u of the k equation (m²/s³).
    pub su_k: T,
    /// S_p ≤ 0 of the k equation (1/s).
    pub sp_k: T,
    /// S_u of the ε equation (m²/s⁴).
    pub su_e: T,
    /// S_p ≤ 0 of the ε equation (1/s).
    pub sp_e: T,
    /// ν_t/σ_k (diffusivity of k without ν).
    pub dk: T,
    /// ν_t/σ_ε.
    pub de: T,
    /// 1 where the cap `max_ratio` of ν_t/ν acted, else 0.
    pub capped: T,
}

/// The variable C_μ of the realizable model from the velocity gradient
/// and the time scale k/ε (Shih et al. 1995, eqs. 19–21 of the report).
#[inline(always)]
pub fn realizable_c_mu<T: Real>(g: &Tensor<T>, k_over_e: T, a0: f64, w: T) -> T {
    let f = T::from_f64;
    let (s, r) = strain_rotation(g);
    let mut ss = T::ZERO;
    let mut rr = T::ZERO;
    let mut sss = T::ZERO;
    let mut i = 0;
    while i < 3 {
        let mut j = 0;
        while j < 3 {
            ss += s[i][j] * s[i][j];
            rr += r[i][j] * r[i][j];
            let mut l = 0;
            while l < 3 {
                sss += s[i][j] * s[j][l] * s[l][i];
                l += 1;
            }
            j += 1;
        }
        i += 1;
    }
    let u_star = (ss + rr).sqrt_of();
    // W̃ = S_ij S_jk S_ki / S̃³ (0 without strain), clipped to ±1/√6
    let st = ss.sqrt_of();
    let st3 = st * st * st;
    let tiny = f(1e-300);
    let wt = if st3 > tiny { sss / st3 } else { T::ZERO };
    let lim = f(1.0 / foehncraft_core::math::sqrt(6.0));
    let wt = smax(smin(wt, lim, w), -lim, w);
    let s6 = f(foehncraft_core::math::sqrt(6.0));
    let x = (s6 * wt).max_of(f(-1.0)).min_of(T::ONE);
    let a_s = s6 * x.cos_acos3_of();
    T::ONE / (f(a0) + a_s * u_star * k_over_e)
}

/// ν_t, P_k and the source split of the k and ε equations in one cell,
/// with the values k*, ε* of the previous outer iteration.
#[inline(always)]
pub fn keps_cell<T: Real>(
    i: &KepsCellInputs<T>,
    c: &KepsConstants,
    o: &KepsOptions,
) -> KepsCell<T> {
    let f = T::from_f64;
    let w_s = f(o.smooth);
    let tiny = f(1e-300);
    let (k, e, nu) = (i.k.max_of(T::ZERO), i.e.max_of(tiny), i.nu);
    let (s2, o2) = s2_omega2(&i.grad_u);
    let s = s2.sqrt_of();
    let k_e = k / e;
    let realizable = o.variant == REALIZABLE;
    let c_mu = if realizable {
        realizable_c_mu(&i.grad_u, k_e, c.a0, w_s)
    } else {
        f(c.c_mu)
    };
    // ν_t = C_μ k²/ε, the Durbin bound (standard) and the cap of ν_t/ν
    let mut nu_t = c_mu * k * k_e;
    if !realizable && o.limiter == LIMITER_DURBIN {
        // ν_t ≤ C_lim k/(√6 |S|), |S| = S/√2: C_lim k/(√3 S)
        let s3 = f(foehncraft_core::math::sqrt(3.0)) * s;
        if s3 > tiny {
            nu_t = smin(nu_t, f(o.c_lim) * k / s3, w_s);
        }
    }
    let cap = f(o.max_ratio) * nu;
    let capped = if nu_t > cap { T::ONE } else { T::ZERO };
    nu_t = nu_t.min_of(cap);
    // production
    let p_over = if !realizable && o.limiter == LIMITER_KATO_LAUNDER {
        s * o2.sqrt_of()
    } else {
        s2
    };
    let mut pk = nu_t * p_over;
    if !realizable && o.limiter == LIMITER_CLIP {
        pk = smin(pk, f(10.0) * e, w_s);
    }
    // k: S_u = P_k (+ ε_amb), S_p = −ε*/k*
    let k_div = k.max_of(tiny);
    let inv_t = e / k_div;
    let su_k = pk + f(o.e_amb);
    let sp_k = -inv_t;
    let (su_e, sp_e) = if realizable {
        // C₁ S ε − C₂ ε²/(k + √(νε)); C₁ = max(0.43, η/(η + 5))
        let eta = s * k_e;
        let c1 = smax(f(c.c1_min), eta / (eta + f(5.0)), w_s);
        let den = k + (nu * e).sqrt_of();
        let c2 = f(c.c_e2);
        let amb = if o.k_amb > 0.0 {
            let ea = o.e_amb;
            f(c.c_e2 * ea * ea / (o.k_amb + foehncraft_core::math::sqrt(i.nu.to_f64() * ea)))
        } else {
            T::ZERO
        };
        let r = c2 * e / den;
        if o.newton != 0 {
            (c1 * s * e + r * e + amb, -f(2.0) * r)
        } else {
            (c1 * s * e + amb, -r)
        }
    } else {
        // (C_ε1 P_k − C_ε2 ε)/T, T = k/ε
        let c2 = f(c.c_e2);
        let amb = if o.k_amb > 0.0 {
            f(c.c_e2 * o.e_amb * o.e_amb / o.k_amb)
        } else {
            T::ZERO
        };
        let r = c2 * inv_t;
        if o.newton != 0 {
            (f(c.c_e1) * pk * inv_t + r * e + amb, -f(2.0) * r)
        } else {
            (f(c.c_e1) * pk * inv_t + amb, -r)
        }
    };
    KepsCell {
        nu_t,
        c_mu,
        pk,
        su_k,
        sp_k,
        su_e,
        sp_e,
        dk: nu_t / f(c.sigma_k),
        de: nu_t / f(c.sigma_e),
        capped,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use foehncraft_core::Dual;

    fn shear(k: f64, e: f64, dudy: f64) -> KepsCellInputs<f64> {
        let mut g = [[0.0; 3]; 3];
        g[0][1] = dudy;
        KepsCellInputs {
            k,
            e,
            grad_u: g,
            nu: 1.5e-5,
        }
    }

    #[test]
    fn standard_model_in_simple_shear() {
        let i = shear(0.01, 0.1, 50.0);
        let c = KepsConstants::standard();
        let r = keps_cell(&i, &c, &KepsOptions::default());
        let nu_t = 0.09 * 1e-4 / 0.1;
        assert!((r.nu_t - nu_t).abs() < 1e-18);
        assert!((r.pk - nu_t * 2500.0).abs() < 1e-15);
        assert!((r.sp_k + 10.0).abs() < 1e-12);
        // Newton: S_u = C_ε1 P ε/k + C_ε2 ε²/k, S_p = −2 C_ε2 ε/k
        assert!((r.su_e - (1.44 * r.pk * 10.0 + 1.92 * 0.1 * 10.0)).abs() < 1e-12);
        assert!((r.sp_e + 2.0 * 1.92 * 10.0).abs() < 1e-12);
        // the implicit form at the current values: S_u + S_p ε = rhs
        let rhs = 10.0 * (1.44 * r.pk - 1.92 * 0.1);
        assert!((r.su_e + r.sp_e * 0.1 - rhs).abs() < 1e-12);
        assert!((r.dk - nu_t).abs() < 1e-18 && (r.de - nu_t / 1.3).abs() < 1e-18);
    }

    #[test]
    fn realizable_c_mu_in_the_log_layer_and_without_strain() {
        // equilibrium simple shear: S k/ε = 1/√C_μ ≈ 3.33 gives C_μ near
        // 0.09 (the paper's calibration point)
        let c = KepsConstants::realizable();
        let (k, e) = (1.0, 1.0);
        let dudy = 3.33;
        let o = KepsOptions {
            variant: REALIZABLE,
            ..Default::default()
        };
        let r = keps_cell(&shear(k, e, dudy), &c, &o);
        // simple shear: W̃ = 0, A_s = √6 cos(π/6) = 3/√2·… = √6·√3/2
        let a_s = 6f64.sqrt() * (core::f64::consts::PI / 6.0).cos();
        let u_star = dudy / 2f64.sqrt() * 2f64.sqrt(); // (S_ijS_ij + W_ijW_ij)^½ = |dudy|/√2·√2
        let exact = 1.0 / (4.0 + a_s * u_star * k / e);
        assert!((r.c_mu - exact).abs() < 1e-12, "{} {exact}", r.c_mu);
        assert!(r.c_mu > 0.08 && r.c_mu < 0.1, "{}", r.c_mu);
        // no strain: C_μ = 1/A₀
        let r0 = keps_cell(&shear(k, e, 0.0), &c, &o);
        assert!((r0.c_mu - 0.25).abs() < 1e-15);
        assert_eq!(r0.pk, 0.0);
    }

    #[test]
    fn realizable_c_mu_stays_realizable_under_strong_strain() {
        // the normal stresses ⅔k − 2ν_t λ stay ≥ 0: ν_t ≤ k/(3 λ_max) for
        // plane and axisymmetric strain with the largest eigenvalue a
        let c = KepsConstants::realizable();
        let o = KepsOptions {
            variant: REALIZABLE,
            ..Default::default()
        };
        for a in [1.0, 10.0, 1e3, 1e5] {
            let mut i = shear(1.0, 1.0, 0.0);
            i.grad_u = [[a, 0.0, 0.0], [0.0, -a, 0.0], [0.0; 3]];
            let r = keps_cell(&i, &c, &o);
            let bound = 1.0 / (3.0 * a);
            assert!(r.nu_t <= bound * (1.0 + 1e-12), "{a}: {} {bound}", r.nu_t);
            // axisymmetric strain reaches W̃ = 1/√6
            i.grad_u = [[a, 0.0, 0.0], [0.0, -0.5 * a, 0.0], [0.0, 0.0, -0.5 * a]];
            let r = keps_cell(&i, &c, &o);
            assert!(r.nu_t <= bound * (1.0 + 1e-12), "{a}: {} {bound}", r.nu_t);
            if a >= 1e3 {
                assert!(r.nu_t > 0.9 * bound);
            }
        }
    }

    #[test]
    fn kato_launder_and_durbin_act_at_a_stagnation_point() {
        let mut i = shear(1.0, 0.1, 0.0);
        i.grad_u = [[10.0, 0.0, 0.0], [0.0, -10.0, 0.0], [0.0; 3]];
        let c = KepsConstants::standard();
        let plain = keps_cell(&i, &c, &KepsOptions::default());
        let kl = KepsOptions {
            limiter: LIMITER_KATO_LAUNDER,
            ..Default::default()
        };
        assert_eq!(keps_cell(&i, &c, &kl).pk, 0.0);
        let du = KepsOptions {
            limiter: LIMITER_DURBIN,
            ..Default::default()
        };
        let r = keps_cell(&i, &c, &du);
        let s = (2.0f64 * 200.0).sqrt();
        assert!((r.nu_t - 1.0 / (3f64.sqrt() * s)).abs() < 1e-15);
        assert!(r.nu_t < plain.nu_t);
        let cl = KepsOptions {
            limiter: LIMITER_CLIP,
            ..Default::default()
        };
        assert!((keps_cell(&i, &c, &cl).pk - 1.0).abs() < 1e-15);
    }

    #[test]
    fn sustaining_terms_balance_the_ambient_decay() {
        for (variant, c) in [
            (STANDARD, KepsConstants::standard()),
            (REALIZABLE, KepsConstants::realizable()),
        ] {
            let i = shear(1e-4, 2e-3, 0.0);
            let o = KepsOptions {
                variant,
                k_amb: 1e-4,
                e_amb: 2e-3,
                ..Default::default()
            };
            let r = keps_cell(&i, &c, &o);
            assert!((r.su_k + r.sp_k * i.k).abs() < 1e-18);
            assert!((r.su_e + r.sp_e * i.e).abs() < 1e-15, "{variant}");
        }
    }

    #[test]
    fn sources_keep_their_signs() {
        for variant in [STANDARD, REALIZABLE] {
            for newton in [0, 1] {
                let o = KepsOptions {
                    variant,
                    newton,
                    ..Default::default()
                };
                let r = keps_cell(&shear(1e-3, 5e-2, 300.0), &KepsConstants::default(), &o);
                assert!(r.su_k >= 0.0 && r.sp_k <= 0.0 && r.su_e >= 0.0 && r.sp_e <= 0.0);
            }
        }
    }

    #[test]
    fn dual_numbers_give_the_derivative_of_the_sources() {
        let i = {
            let mut i = shear(0.02, 0.3, 40.0);
            i.grad_u[1][0] = 5.0;
            i.grad_u[0][0] = 3.0;
            i.grad_u[1][1] = -2.0;
            i.grad_u[2][2] = -1.0;
            i
        };
        let c = KepsConstants::realizable();
        let o = KepsOptions {
            variant: REALIZABLE,
            smooth: 1e-3,
            ..Default::default()
        };
        let lift = |x: f64| Dual::<1>::constant(x);
        for which in 0..2 {
            let di = KepsCellInputs {
                k: if which == 0 {
                    Dual::<1>::variable(i.k, 0)
                } else {
                    lift(i.k)
                },
                e: if which == 1 {
                    Dual::<1>::variable(i.e, 0)
                } else {
                    lift(i.e)
                },
                grad_u: i.grad_u.map(|r| r.map(lift)),
                nu: lift(i.nu),
            };
            let r = keps_cell(&di, &c, &o);
            let h = 1e-6 * if which == 0 { i.k } else { i.e };
            let at = |x: f64| {
                let mut j = i;
                if which == 0 {
                    j.k = x;
                } else {
                    j.e = x;
                }
                keps_cell(&j, &c, &o)
            };
            let x0 = if which == 0 { i.k } else { i.e };
            let (p, m) = (at(x0 + h), at(x0 - h));
            for (dv, fd) in [
                (r.su_e.d[0], (p.su_e - m.su_e) / (2.0 * h)),
                (r.nu_t.d[0], (p.nu_t - m.nu_t) / (2.0 * h)),
                (r.su_k.d[0], (p.su_k - m.su_k) / (2.0 * h)),
            ] {
                assert!(
                    (dv - fd).abs() <= 1e-5 * fd.abs().max(1e-8),
                    "{which}: {dv} {fd}"
                );
            }
        }
    }
}

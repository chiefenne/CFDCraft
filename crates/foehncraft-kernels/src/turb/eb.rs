//! The elliptic-blending k-ε (D133; stage 6 specification 3.4) per cell
//! (`form` = [`FORM_EB`] of the ζ-f constants; F. Billard,
//! D. Laurence, *A robust k-ε-v²/k elliptic blending turbulence model
//! applied to near-wall, separated and buoyant flows*, Int. J. Heat Fluid
//! Flow 33, 45–58 (2012), the "BL-v²/k" model) shares the fields and the
//! solvers: φ = v²/k in the slot of ζ and the blending parameter α in the
//! slot of f:
//!
//!   ν_t = C_μ φ k min(T, T_lim)
//!
//!   α − L²∇²α = 1, α_w = 0
//!
//!   Dk/Dt = P_k − ε − E + ∇·[(ν/2 + ν_t/σ_k)∇k],
//!   E = C_ε3 (1 − α)³ (k/ε) 2ν ν_t (∂_l ∂_j U_i)(∂_l ∂_j U_i)
//!
//!   Dε/Dt = (C_ε1 P_k − C*_ε2 ε)/T + ∇·[(ν/2 + ν_t/σ_ε)∇ε],
//!   C*_ε2 = C_ε2 + α³ (C_ε4 − C_ε2) tanh(|∇·(ν_t/σ_k ∇k)/ε|^{3/2})
//!
//!   Dφ/Dt = (1 − α³) f_w + α³ f_h − (φ/k) P_k + (2/k)(ν_t/σ_k) ∇φ·∇k
//!           + ∇·[(ν/2 + ν_t/σ_φ)∇φ],
//!   f_w = −(ε/2)(φ/k), f_h = −(1/T)(C₁ − 1 + C₂ P_k/ε)(φ − 2/3)
//!
//! with T = (k²/ε² + C_T² ν/ε)^½, L = C_L (k³/ε² + C_η² ν^{3/2}/ε^{1/2})^½,
//! T_lim = 0.6/(√6 C_μ φ (S_ij S_ij)^½) and the constants C_μ = 0.22,
//! C_ε1 = 1.44, C_ε2 = 1.83, C_ε3 = 2.3, C_ε4 = 0.4, C₁ = 1.7, C₂ = 0.9,
//! σ_k = 1, σ_ε = 1.5, σ_φ = 1, C_T = 4, C_L = 0.164, C_η = 75 (the
//! paper's Appendix A, eqs. (A.1)–(A.7), and Table A.10). The transported
//! ε is the homogeneous part of the dissipation; to keep the near-wall
//! balances, the paper halves the viscous diffusion of k, ε and φ (hence
//! also the ½ in f_w) and takes ε_w = ν k_P/y_P² (its eq. (18); k_w =
//! φ_w = α_w = 0). Checked against the paper on 2026-11-18 (the ◇ of the
//! specification 3.4); the draft of the specification, written from
//! memory, differed in the ν/2, the place of E, C*_ε2 and the smooth T
//! and L. The paper's channel result (its Table 6: C_f 100, 98 and 100 %
//! of the DNS at Re_τ = 180, 395 and 590) is the acceptance target.

use super::invariants::s2_omega2;
use super::smin;
#[cfg(doc)]
use super::zeta_f::FORM_EB;
use super::zeta_f::{ZETA_MAX, ZETA_MIN, ZetaFCell, ZetaFCellInputs, ZetaFConstants, ZetaFOptions};
use foehncraft_core::Real;

/// T and L of the elliptic-blending k-ε: T = (k²/ε² + C_T² ν/ε)^½,
/// L = C_L (k³/ε² + C_η² ν^{3/2}/ε^{1/2})^½.
#[inline(always)]
pub fn eb_scales<T: Real>(k: T, e: T, nu: T, c: &ZetaFConstants) -> (T, T) {
    let f = T::from_f64;
    let t_k = k / e;
    let t = (t_k * t_k + f(c.c_t * c.c_t) * nu / e).sqrt_of();
    let l = f(c.c_l)
        * (t_k * t_k * k + f(c.c_eta * c.c_eta) * nu * nu.sqrt_of() / e.sqrt_of()).sqrt_of();
    (t, l)
}

/// [`super::zeta_f::zeta_f_cell`] of the elliptic-blending k-ε ([`FORM_EB`]): `zeta` is
/// φ = v²/k, `f` is α (the solution of the α equation), and the inputs
/// `grad_k`, `grad_z`, `g2` and `dkt` enter. The source split keeps
/// S_u ≥ 0 and S_p ≤ 0: E/k, the near-wall f_w/φ, the homogeneous
/// α³ f_h and P_k/k are sinks; the cross diffusion goes to S_u or, divided
/// by φ*, to S_p by its sign.
#[inline(always)]
pub fn eb_cell<T: Real>(
    i: &ZetaFCellInputs<T>,
    c: &ZetaFConstants,
    o: &ZetaFOptions,
) -> ZetaFCell<T> {
    let f = T::from_f64;
    let w = f(o.smooth);
    let tiny = f(1e-300);
    let (k, e, nu) = (i.k.max_of(T::ZERO), i.e.max_of(tiny), i.nu);
    let phi = i.zeta.max_of(T::ZERO).min_of(f(ZETA_MAX));
    let al = i.f.max_of(T::ZERO).min_of(T::ONE);
    let (s2, _) = s2_omega2(&i.grad_u);
    let (t, l) = eb_scales(k, e, nu, c);
    // α − L²∇²α = 1 with L² outside the divergence, as the paper writes
    // it: ∇²α + (1 − α)/L² = 0 (diffusivity 1, S_u = 1/L², S_p = −1/L²)
    let inv_l2 = T::ONE / (l * l).max_of(tiny);
    // T_lim = a/(√6 C_μ φ (S_ij S_ij)^½), S_ij S_ij = S²/2
    let den = f(foehncraft_core::math::sqrt(6.0) * c.c_mu) * phi * (f(0.5) * s2).sqrt_of();
    let t_nu = if den > tiny {
        smin(t, f(c.a) / den, w)
    } else {
        t
    };
    let raw = f(c.c_mu) * phi * k * t_nu;
    let cap = f(o.max_ratio) * nu;
    let over = raw > cap;
    let capped = if over { T::ONE } else { T::ZERO };
    let nu_t = raw.min_of(cap);
    let k_div = k.max_of(tiny);
    // ν_t/k, finite as k → 0
    let nt_k = if over {
        cap / k_div
    } else {
        f(c.c_mu) * phi * t_nu
    };
    let m = i.wall_m;
    let pk = (T::ONE - m) * nu_t * s2 + i.wall_pk;
    let a3 = al * al * al;
    let b = T::ONE - al;
    let b3 = b * b * b;
    let inv_t = T::ONE / t;
    // k: P_k − ε − E; E/k = C_ε3 (1 − α)³ 2ν ν_t G/ε
    let su_k = pk + f(o.e_amb);
    let sp_k = -(e / k_div) - f(c.c_e3) * b3 * f(2.0) * nu * nu_t * i.g2.max_of(T::ZERO) / e;
    // ε: (C_ε1 P_k − C*_ε2 ε)/T
    let x = (i.dkt / e).abs();
    let ce2 = f(c.c_e2) + a3 * f(c.c_e4 - c.c_e2) * (x * x.sqrt_of()).tanh_of();
    let amb_e = if o.k_amb > 0.0 && o.e_amb > 0.0 {
        let (ka, ea, nu0) = (o.k_amb, o.e_amb, i.nu.to_f64());
        let ta = foehncraft_core::math::sqrt((ka / ea) * (ka / ea) + c.c_t * c.c_t * nu0 / ea);
        f(c.c_e2 * ea / ta)
    } else {
        T::ZERO
    };
    let r = ce2 * inv_t;
    let ce1p = f(c.c_e1) * pk * inv_t;
    let (su_e, sp_e) = if o.newton != 0 {
        (ce1p + r * e + amb_e, -f(2.0) * r)
    } else {
        (ce1p + amb_e, -r)
    };
    // φ: (1 − α³) f_w + α³ f_h − (φ/k) P_k + (2/k)(ν_t/σ_k) ∇φ·∇k
    let ch = (f(c.c1 - 1.0) + f(c.c2) * pk / e) * inv_t;
    let gk = i.grad_k;
    let gz = i.grad_z;
    let cross = f(2.0) * nt_k / f(c.sigma_k) * (gz[0] * gk[0] + gz[1] * gk[1] + gz[2] * gk[2]);
    let su_z = a3 * ch * f(2.0 / 3.0) + cross.max_of(T::ZERO);
    let sp_z = -(T::ONE - a3) * f(0.5) * e / k_div - a3 * ch - pk / k_div
        + cross.min_of(T::ZERO) / phi.max_of(f(ZETA_MIN));
    ZetaFCell {
        nu_t,
        t,
        l2: T::ONE,
        pk,
        su_k,
        sp_k,
        su_e,
        sp_e,
        su_z,
        sp_z,
        rhs_f: inv_l2,
        sp_f: -inv_l2,
        dk: nu_t / f(c.sigma_k),
        de: nu_t / f(c.sigma_e),
        dz: nu_t / f(c.sigma_z),
        capped,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::turb::zeta_f::zeta_f_cell;
    use foehncraft_core::Dual;

    #[allow(clippy::too_many_arguments)]
    fn cell(
        k: f64,
        e: f64,
        phi: f64,
        alpha: f64,
        dudy: f64,
        g2: f64,
        dkt: f64,
        cross: f64,
    ) -> ZetaFCellInputs<f64> {
        let mut g = [[0.0; 3]; 3];
        g[0][1] = dudy;
        ZetaFCellInputs {
            k,
            e,
            zeta: phi,
            f: alpha,
            grad_u: g,
            nu: 1.5e-5,
            wall_m: 0.0,
            wall_pk: 0.0,
            grad_k: [0.0, 2.0, 0.0],
            grad_z: [0.0, 0.5 * cross, 0.0],
            g2,
            dkt,
        }
    }

    /// The net sources are the terms of the paper's Appendix A, evaluated
    /// here from its formulas, and the split keeps S_u ≥ 0, S_p ≤ 0.
    #[test]
    fn sources_are_the_terms_of_the_paper() {
        let c = ZetaFConstants::billard_laurence_2012();
        let o = ZetaFOptions::default();
        for &(k, e, phi, al, s, g2, dkt, cross) in &[
            (0.01, 0.1, 0.4, 0.9, 30.0, 0.0, 0.0, 0.0),
            (0.02, 0.3, 0.1, 0.3, 200.0, 4e6, 0.05, 3.0),
            (1e-4, 0.02, 0.02, 0.05, 900.0, 1e8, -0.01, -5.0),
            (0.5, 0.2, 0.6, 1.0, 1.0, 1.0, 0.4, 0.2),
        ] {
            let i = cell(k, e, phi, al, s, g2, dkt, cross);
            let r = zeta_f_cell(&i, &c, &o);
            assert_eq!(r, eb_cell(&i, &c, &o));
            let nu = i.nu;
            let t = ((k / e).powi(2) + c.c_t * c.c_t * nu / e).sqrt();
            let l =
                c.c_l * (k.powi(3) / (e * e) + c.c_eta * c.c_eta * nu.powf(1.5) / e.sqrt()).sqrt();
            // S_ij S_ij = s²/2 in simple shear
            let t_lim = 0.6 / (6f64.sqrt() * c.c_mu * phi * (0.5 * s * s).sqrt());
            let nu_t = c.c_mu * phi * k * t.min(t_lim);
            let p = nu_t * s * s;
            let big_e = c.c_e3 * (1.0 - al).powi(3) * (k / e) * 2.0 * nu * nu_t * g2;
            let ce2 = c.c_e2 + al.powi(3) * (c.c_e4 - c.c_e2) * (dkt / e).abs().powf(1.5).tanh();
            let fw = -0.5 * e * phi / k;
            let fh = -(c.c1 - 1.0 + c.c2 * p / e) * (phi - 2.0 / 3.0) / t;
            let gdot = 2.0 * 0.5 * cross;
            let net_phi = (1.0 - al.powi(3)) * fw + al.powi(3) * fh - phi * p / k
                + 2.0 / k * nu_t / c.sigma_k * gdot;
            let close = |a: f64, b: f64| (a - b).abs() <= 1e-12 * b.abs().max(1e-30);
            assert!(close(r.nu_t, nu_t), "ν_t {} {nu_t}", r.nu_t);
            assert!(close(r.t, t) && close(r.l2, 1.0), "T, L²");
            assert!(close(r.su_k + r.sp_k * k, p - e - big_e), "k");
            assert!(close(r.su_e + r.sp_e * e, (c.c_e1 * p - ce2 * e) / t), "ε");
            assert!(close(r.su_z + r.sp_z * phi, net_phi), "φ");
            assert!(close(r.rhs_f + r.sp_f * al, (1.0 - al) / (l * l)), "α");
            assert!(r.su_k >= 0.0 && r.su_e >= 0.0 && r.su_z >= 0.0 && r.rhs_f > 0.0);
            assert!(r.sp_k <= 0.0 && r.sp_e <= 0.0 && r.sp_z <= 0.0 && r.sp_f < 0.0);
            assert!(close(r.dz, nu_t / c.sigma_z) && close(r.de, nu_t / c.sigma_e));
        }
    }

    /// Towards a wall (k → 0) T and L take their Kolmogorov bounds
    /// C_T (ν/ε)^½ and C_L C_η (ν³/ε)^¼.
    #[test]
    fn scales_take_the_kolmogorov_bounds_at_a_wall() {
        let c = ZetaFConstants::billard_laurence_2012();
        let (nu, e): (f64, f64) = (1.5e-5, 0.4);
        let (t, l) = eb_scales(1e-12, e, nu, &c);
        assert!((t / (c.c_t * (nu / e).sqrt()) - 1.0).abs() < 1e-9, "{t}");
        assert!(
            (l / (c.c_l * c.c_eta * (nu * nu * nu / e).powf(0.25)) - 1.0).abs() < 1e-9,
            "{l}"
        );
    }

    #[test]
    fn dual_numbers_give_the_derivative_of_the_sources() {
        let mut i = cell(0.02, 0.3, 0.35, 0.6, 40.0, 2e5, 0.08, 1.5);
        i.grad_u[1][0] = 5.0;
        i.grad_u[0][0] = 3.0;
        i.grad_u[1][1] = -2.0;
        i.grad_u[2][2] = -1.0;
        let c = ZetaFConstants::billard_laurence_2012();
        let o = ZetaFOptions {
            smooth: 1e-3,
            ..Default::default()
        };
        let lift = |x: f64| Dual::<1>::constant(x);
        for which in 0..3 {
            let pick = |j: usize, x: f64| {
                if j == which {
                    Dual::<1>::variable(x, 0)
                } else {
                    lift(x)
                }
            };
            let di = ZetaFCellInputs {
                k: pick(0, i.k),
                e: pick(1, i.e),
                zeta: pick(2, i.zeta),
                f: lift(i.f),
                grad_u: i.grad_u.map(|r| r.map(lift)),
                nu: lift(i.nu),
                wall_m: lift(0.0),
                wall_pk: lift(0.0),
                grad_k: i.grad_k.map(lift),
                grad_z: i.grad_z.map(lift),
                g2: lift(i.g2),
                dkt: lift(i.dkt),
            };
            let r = eb_cell(&di, &c, &o);
            let x0 = [i.k, i.e, i.zeta][which];
            let h = 1e-5 * x0;
            let at = |x: f64| {
                let mut j = i;
                match which {
                    0 => j.k = x,
                    1 => j.e = x,
                    _ => j.zeta = x,
                }
                eb_cell(&j, &c, &o)
            };
            let (p, m) = (at(x0 + h), at(x0 - h));
            for (dv, fd) in [
                (r.su_e.d[0], (p.su_e - m.su_e) / (2.0 * h)),
                (r.sp_k.d[0], (p.sp_k - m.sp_k) / (2.0 * h)),
                (r.nu_t.d[0], (p.nu_t - m.nu_t) / (2.0 * h)),
                (r.rhs_f.d[0], (p.rhs_f - m.rhs_f) / (2.0 * h)),
                (r.sp_z.d[0], (p.sp_z - m.sp_z) / (2.0 * h)),
            ] {
                assert!(
                    (dv - fd).abs() <= 1e-5 * fd.abs().max(1e-8),
                    "{which}: {dv} {fd}"
                );
            }
        }
    }
}

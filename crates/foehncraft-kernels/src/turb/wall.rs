//! Wall treatments of SST per wall face (D134; stage 6 specification 4.1,
//! 4.3, 4.5).
//!
//! - **`low-re`** (integration to the wall): τ_w = μ U_t/y_P; the first
//!   cell takes ω_P = 6ν/(β₁ y_P²), the near-wall solution of the ω
//!   equation (Menter 1994), and k_w = 0.
//! - **`automatic`** (Menter, Kuntz & Langtry 2003): u_τ = (u_τ,vis⁴ +
//!   u_τ,log⁴)^¼ with u_τ,vis = (ν U_t/y_P)^½ and u_τ,log from the log law
//!   U_t = (u_τ/κ) ln(E y_P u_τ/ν) (fixed number of Newton steps), and
//!   ω_P = (ω_vis² + ω_log²)^½, ω_vis = 6ν/(β₁ y_P²),
//!   ω_log = u_τ/(√β* κ y_P); the wall shear τ_w = ρ u_τ² enters the
//!   momentum equation as the effective wall viscosity μ_w = τ_w y_P/U_t,
//!   k has zero flux through the wall.
//! - **Rough walls** (Aupoix 2015, the set fitted to Colebrook's
//!   correlation): k_w = u_τ² k_w⁺ and ω_w = (u_τ²/ν) ω_w⁺ as wall values
//!   with
//!   k_w⁺ = max(0, (1/√β*) tanh[(ln(k_s⁺/30)/ln 10 + 1 − tanh(k_s⁺/125))
//!   tanh(k_s⁺/125)]),
//!   ω_w⁺ = 300/k_s⁺² [tanh(15/(4 k_s⁺))]⁻¹ + 191/k_s⁺ [1 − exp(−k_s⁺/250)],
//!   k_s⁺ = u_τ k_s/ν (B. Aupoix, *Roughness corrections for the k–ω shear
//!   stress transport model: status and proposals*, J. Fluids Eng. 137,
//!   021202 (2015); formulas checked on 2026-11-17 against their
//!   reproduction as eq. (11) of EUCASS2017-372, which names the source).
//!   With k_w > 0 the eddy viscosity at the wall is k_w/ω_w, so the wall
//!   shear is (μ + ρ ν_t,w) U_t/y_P.
//!
//! The k-ε family (stage 6 PR 6c; specification 4.3) has its own kernel,
//! [`keps_wall_face`], with the velocity scale u* = C_μ^¼ k_P^½ of the
//! near-wall cell in all three treatments (the equilibrium scale of
//! Launder & Spalding 1974): y* = u* y_P/ν, and the wall shear stress
//! τ_w/ρ = u* U_t/u⁺(y*), i.e. the effective wall viscosity
//! ν_w = ν y*/u⁺(y*). The viscous limit u⁺ = y* gives ν_w = ν exactly,
//! for any k; the log limit is the wall function of Launder & Spalding.
//! No iteration is needed, and u* stays finite where U_t vanishes
//! (separation and reattachment), unlike u_τ.
//!
//! - **`scalable`** (G. Grotjans, F. R. Menter, *Wall functions for general
//!   application CFD codes*, ECCOMAS 98, 1112–1117 (1998)): the log law at
//!   ỹ* = max(y*, 11.06), ε_P = u*³/(κ ỹ), P_k,P = (τ_w/ρ) u*/(κ ỹ) with
//!   ỹ = ỹ* ν/u*;
//! - **`standard`** (Launder & Spalding 1974): the log law for
//!   y* ≥ 11.06 and the viscous sublayer u⁺ = y* below (with P_k,P = 0
//!   there), ε_P = u*³/(κ y_P);
//! - **`compound`** (M. Popovac, K. Hanjalić, *Compound wall treatment for
//!   RANS computation of complex turbulent flows and heat transfer*, Flow
//!   Turbul. Combust. 78, 177–202 (2007)): u⁺ = e^{−Γ} y* + e^{−1/Γ} u⁺_log,
//!   Γ = 0.01 y*⁴/(1 + 5 y*) (the blending of Kader 1981 with Pr = 1),
//!   ε_P = e^{−Γ} 2νk_P/y_P² + e^{−1/Γ} u*³/(κ y_P), and the production
//!   of the cell e^{−Γ} P_k + e^{−1/Γ} (τ_w/ρ) u*/(κ y_P). Γ is evaluated
//!   with y* (the specification writes y⁺; the two agree in the log layer,
//!   and in the sublayer, where y* < y⁺, the viscous branch is exact).
//! - **Rough walls** (all three): u⁺_log = ln(E y*)/κ − ΔB(k_s⁺),
//!   k_s⁺ = u* k_s/ν, with ΔB of T. Cebeci, P. Bradshaw, *Momentum
//!   transfer in boundary layers*, Hemisphere (1977): 0 for k_s⁺ < 2.25,
//!   (1/κ) ln[(k_s⁺ − 2.25)/87.75 + C_s k_s⁺] sin[0.4258 (ln k_s⁺ − 0.811)]
//!   up to 90, (1/κ) ln(1 + C_s k_s⁺) above (the pieces join
//!   continuously). The log-law u⁺ is kept at least 1 (a roughness height
//!   above the first-cell centre would make it vanish; BC027 warns).

use super::sst::SstConstants;
use foehncraft_core::Real;

/// E of the log law u⁺ = (1/κ) ln(E y⁺) (B = 5.2 with κ = 0.41:
/// E = e^{κB} = 8.43; E = 9.793 is the value of the specification 4.2
/// with B ≈ 5.56). The specification fixes E = 9.793.
pub const LOG_LAW_E: f64 = 9.793;

/// Newton steps of the log law (fixed: no data-dependent loop on the GPU).
pub const NEWTON_STEPS: usize = 12;

/// Wall treatment: integration to the wall.
pub const LOW_RE: u32 = 0;
/// Wall treatment: the automatic blending of Menter et al. (2003).
pub const AUTOMATIC: u32 = 1;
/// Wall treatment (k-ε family): scalable wall functions (Grotjans & Menter
/// 1998).
pub const SCALABLE: u32 = 2;
/// Wall treatment (k-ε family): standard wall functions (Launder &
/// Spalding 1974).
pub const STANDARD: u32 = 3;
/// Wall treatment (k-ε family): the compound treatment (Popovac &
/// Hanjalić 2007).
pub const COMPOUND: u32 = 4;

/// The lower end of the log layer, y⁺ = 11.06 (where u⁺ = y⁺ meets the log
/// law with κ = 0.41, E = 9.793; the limit of the scalable wall functions
/// of Grotjans & Menter 1998).
pub const Y_PLUS_LOG: f64 = 11.06;

/// u_τ of the log law U = (u_τ/κ) ln(E max(y u_τ/ν, 11.06)), by Newton's
/// method from u₀ = κ U.
///
/// The floor of y⁺ at the lower end of the log layer makes u_τ,log vanish
/// in the viscous sublayer (u_τ,log = κU/ln(11.06 E) ≈ 0.088 U there,
/// against u_τ,vis = (νU/y)^½ ≫ U for y⁺ < 1), so that the quartic blend
/// of the automatic treatment tends to the viscous value; without it the
/// log law has a spurious root at E y⁺ ≈ 1 whenever the cell lies in the
/// sublayer. The function g(u) = (u/κ) ln(E max(y u/ν, 11.06)) − U is
/// convex and increasing, and g(κU) ≥ U (ln(11.06 E) − 1) > 0, so the
/// iterates fall monotonically to the root.
#[inline(always)]
pub fn log_law_u_tau<T: Real>(u: T, y: T, nu: T, kappa: T) -> T {
    let e = T::from_f64(LOG_LAW_E);
    let y_nu = y / nu;
    let floor = T::from_f64(Y_PLUS_LOG);
    let mut ut = kappa * u;
    let mut i = 0;
    while i < NEWTON_STEPS {
        let yp = y_nu * ut;
        let lin = yp < floor;
        let l = (e * yp.max_of(floor)).ln_of();
        let g = ut * l / kappa - u;
        // d/du [u ln(E y u/ν)] = ln + 1 in the log layer, ln(E 11.06) below
        let dg = if lin { l / kappa } else { (l + T::ONE) / kappa };
        ut -= g / dg;
        i += 1;
    }
    ut
}

/// Inputs of [`wall_face`] for one wall face.
#[derive(Clone, Copy, Debug)]
pub struct WallInputs<T> {
    /// Tangential speed of the near-wall cell relative to the wall (m/s).
    pub ut: T,
    /// Normal distance y_P of the cell centroid from the face (m).
    pub y: T,
    /// ν (m²/s).
    pub nu: T,
    /// Equivalent sand-grain roughness k_s (m), 0 for a smooth wall.
    pub ks: T,
    /// u_τ of the last outer iteration (rough walls: k_s⁺).
    pub u_tau_old: T,
}

/// Results of [`wall_face`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallFace<T> {
    /// u_τ (m/s).
    pub u_tau: T,
    /// y⁺ = u_τ y_P/ν.
    pub y_plus: T,
    /// ω of the near-wall cell (automatic treatment: weighted by
    /// 1 − `blending` in its row) (1/s).
    pub omega: T,
    /// ω at the wall face (Dirichlet value): smooth walls
    /// 10·6ν/(β₁ d₁²) with d₁ = 2 y_P, the wall-cell height (Menter 1994);
    /// rough walls the value of Aupoix (1/s).
    pub omega_wall: T,
    /// k at the wall (rough walls; 0 smooth).
    pub k_wall: T,
    /// Effective kinematic wall viscosity ν_w = τ_w y_P/(ρ U_t) (m²/s).
    pub nu_wall: T,
    /// The weight of the viscous part, u_τ,vis⁴/u_τ⁴ (1: viscous sublayer).
    pub blending: T,
    /// The production of k of the log layer at the near-wall cell,
    /// u_τ³/(κ y_P) (m²/s³), which the automatic treatment blends in with
    /// the weight 1 − `blending` (the cell's velocity gradient U_t/y_P is
    /// far too large once y_P is in the log layer).
    pub pk_wall: T,
}

/// k_w⁺ of Aupoix (2015), Colebrook set.
#[inline(always)]
pub fn aupoix_k_plus<T: Real>(ks_plus: T, beta_star: T) -> T {
    let f = T::from_f64;
    let t = (ks_plus / f(125.0)).tanh_of();
    let a = (ks_plus / f(30.0)).ln_of() / f(core::f64::consts::LN_10) + T::ONE - t;
    ((a * t).tanh_of() / beta_star.sqrt_of()).max_of(T::ZERO)
}

/// ω_w⁺ of Aupoix (2015), Colebrook set.
#[inline(always)]
pub fn aupoix_omega_plus<T: Real>(ks_plus: T) -> T {
    let f = T::from_f64;
    let ks2 = ks_plus * ks_plus;
    f(300.0) / ks2 / (f(15.0) / (f(4.0) * ks_plus)).tanh_of()
        + f(191.0) / ks_plus * (T::ONE - (-ks_plus / f(250.0)).exp_of())
}

/// The wall treatment `treatment` ([`LOW_RE`] or [`AUTOMATIC`]) for one
/// wall face.
#[inline(always)]
pub fn wall_face<T: Real>(i: &WallInputs<T>, c: &SstConstants, treatment: u32) -> WallFace<T> {
    let f = T::from_f64;
    let (ut, y, nu) = (i.ut, i.y, i.nu);
    let kappa = f(c.kappa);
    let w_vis = f(6.0) * nu / (f(c.beta1) * y * y);
    // 10·6ν/(β₁ (2 y_P)²)
    let omega_wall = f(2.5) * w_vis;
    let tiny = f(1e-300);
    if i.ks > T::ZERO {
        // rough wall: wall values of Aupoix with u_τ of the last iteration
        // (at least that of the viscous wall shear)
        let u_vis = (nu * ut / y).sqrt_of();
        let u_tau = i.u_tau_old.max_of(u_vis);
        let ksp = (u_tau * i.ks / nu).max_of(tiny);
        let u2 = u_tau * u_tau;
        let k_wall = u2 * aupoix_k_plus(ksp, f(c.beta_star));
        let omega = u2 / nu * aupoix_omega_plus(ksp);
        return WallFace {
            u_tau,
            y_plus: u_tau * y / nu,
            omega,
            omega_wall: omega,
            k_wall,
            nu_wall: nu + k_wall / omega,
            blending: T::ONE,
            pk_wall: T::ZERO,
        };
    }
    if treatment == LOW_RE {
        let u_tau = (nu * ut / y).sqrt_of();
        return WallFace {
            u_tau,
            y_plus: u_tau * y / nu,
            omega: w_vis,
            omega_wall,
            k_wall: T::ZERO,
            nu_wall: nu,
            blending: T::ONE,
            pk_wall: T::ZERO,
        };
    }
    let u_vis2 = nu * ut / y;
    let u_log = log_law_u_tau(ut, y, nu, kappa);
    let v4 = u_vis2 * u_vis2;
    let l2 = u_log * u_log;
    let l4 = l2 * l2;
    let u4 = v4 + l4;
    let u_tau = u4.sqrt_of().sqrt_of();
    let w_log = u_tau / (f(c.beta_star).sqrt_of() * kappa * y);
    let omega = (w_vis * w_vis + w_log * w_log).sqrt_of();
    // ν_w = u_τ² y/U_t (U_t > 0; ν in a fluid at rest)
    let nu_wall = if ut > T::ZERO {
        u4.sqrt_of() * y / ut
    } else {
        nu
    };
    WallFace {
        u_tau,
        y_plus: u_tau * y / nu,
        omega,
        omega_wall,
        k_wall: T::ZERO,
        nu_wall,
        blending: v4 / u4.max_of(tiny),
        pk_wall: u_tau * u_tau * u_tau / (kappa * y),
    }
}

/// ΔB(k_s⁺) of Cebeci & Bradshaw (1977) with the roughness constant
/// C_s (0 below k_s⁺ = 2.25, transitional up to 90, fully rough above).
#[inline(always)]
pub fn roughness_shift<T: Real>(ks_plus: T, cs: T, kappa: T) -> T {
    let f = T::from_f64;
    let lo = f(2.25);
    let hi = f(90.0);
    let k = ks_plus.max_of(lo);
    // transitional
    let tr = ((k - lo) / f(87.75) + cs * k).ln_of() * (f(0.4258) * (k.ln_of() - f(0.811))).sin_of();
    // fully rough
    let fr = (T::ONE + cs * k).ln_of();
    let v = if ks_plus < lo {
        T::ZERO
    } else if ks_plus < hi {
        tr
    } else {
        fr
    };
    v / kappa
}

/// The wall values of the k-ε-ζ-f model (Hanjalić, Popovac & Hadžiabdić
/// 2004) at a wall face from its near-wall cell: ε_w = 2ν k_P/y_P² and
/// f_w = −2ν ζ_P/y_P² (the limits for y → 0 with the first-cell values,
/// k ~ y² and ζ ~ y² near the wall). With the compound treatment
/// (Popovac & Hanjalić 2007) f_w is blended towards 0 with the viscous
/// weight e^{−Γ} (`viscous`, 1 for the integration to the wall). Returns
/// (ε_w, f_w).
#[inline(always)]
pub fn zeta_f_wall_values<T: Real>(k: T, zeta: T, y: T, nu: T, viscous: T) -> (T, T) {
    let two_nu_y2 = T::from_f64(2.0) * nu / (y * y);
    (
        two_nu_y2 * k.max_of(T::ZERO),
        -(viscous * two_nu_y2 * zeta.max_of(T::ZERO)),
    )
}

/// Γ = 0.01 y⁴/(1 + 5 y) of the compound blending (Kader 1981, Pr = 1).
#[inline(always)]
pub fn compound_gamma<T: Real>(y: T) -> T {
    let y2 = y * y;
    T::from_f64(0.01) * y2 * y2 / (T::ONE + T::from_f64(5.0) * y)
}

/// Inputs of [`keps_wall_face`] for one wall face.
#[derive(Clone, Copy, Debug)]
pub struct KepsWallInputs<T> {
    /// Tangential speed of the near-wall cell relative to the wall (m/s).
    pub ut: T,
    /// Normal distance y_P of the cell centroid from the face (m).
    pub y: T,
    /// ν (m²/s).
    pub nu: T,
    /// k of the near-wall cell (m²/s²).
    pub k: T,
    /// Equivalent sand-grain roughness k_s (m), 0 for a smooth wall.
    pub ks: T,
    /// C_s of Cebeci & Bradshaw.
    pub cs: T,
}

/// Results of [`keps_wall_face`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KepsWallFace<T> {
    /// u_τ = (τ_w/ρ)^½ (m/s).
    pub u_tau: T,
    /// y⁺ = u_τ y_P/ν.
    pub y_plus: T,
    /// y* = C_μ^¼ k_P^½ y_P/ν.
    pub y_star: T,
    /// ε of the near-wall cell (fixed value) (m²/s³).
    pub epsilon: T,
    /// Effective kinematic wall viscosity ν_w = τ_w y_P/(ρ U_t) (m²/s).
    pub nu_wall: T,
    /// Weight m_p of the wall production: the cell takes
    /// (1 − m_p) P_k + `pk_wall`.
    pub prod_weight: T,
    /// The wall part of the production (already weighted) (m²/s³).
    pub pk_wall: T,
    /// The viscous weight e^{−Γ} of the compound blending (1 in the
    /// sublayer of the standard wall function, else 0).
    pub viscous: T,
}

/// The wall treatment `treatment` ([`SCALABLE`], [`STANDARD`] or
/// [`COMPOUND`]) of the k-ε family for one wall face.
#[inline(always)]
pub fn keps_wall_face<T: Real>(
    i: &KepsWallInputs<T>,
    c: &super::keps::KepsConstants,
    treatment: u32,
) -> KepsWallFace<T> {
    let f = T::from_f64;
    let (ut, y, nu) = (i.ut, i.y, i.nu);
    let kappa = f(c.kappa);
    let tiny = f(1e-300);
    let k = i.k.max_of(T::ZERO);
    let u_star = f(foehncraft_core::math::sqrt(foehncraft_core::math::sqrt(
        c.c_mu,
    ))) * k.sqrt_of();
    let ys = u_star * y / nu;
    let shift = if i.ks > T::ZERO {
        roughness_shift(u_star * i.ks / nu, i.cs, kappa)
    } else {
        T::ZERO
    };
    let e = f(LOG_LAW_E);
    let floor = f(Y_PLUS_LOG);
    let log_plus = |yy: T| ((e * yy).ln_of() / kappa - shift).max_of(T::ONE);
    let u3 = u_star * u_star * u_star;
    // ratio r = y*/u⁺ (ν_w = ν r), the length of the log-layer values and
    // the weights
    let (r, y_len, viscous, w_log, eps_vis) = if treatment == COMPOUND {
        let g = compound_gamma(ys).max_of(tiny);
        let a = (-g).exp_of();
        let b = (-(T::ONE / g)).exp_of();
        let yl = ys.max_of(T::ONE);
        // u⁺/y* = a + b u⁺_log/y* (finite as y* → 0, where b = 0)
        let r = T::ONE / (a + b * log_plus(yl) / ys.max_of(tiny));
        (r, y, a, b, a * f(2.0) * nu * k / (y * y))
    } else if treatment == STANDARD {
        let log = ys >= floor;
        if log {
            (ys / log_plus(ys), y, T::ZERO, T::ONE, T::ZERO)
        } else {
            (T::ONE, y, T::ONE, T::ZERO, T::ZERO)
        }
    } else {
        // scalable: the log law at ỹ* = max(y*, 11.06), ỹ = ỹ* ν/u*
        let yt = ys.max_of(floor);
        let yl = y.max_of(floor * nu / u_star.max_of(tiny));
        (ys / log_plus(yt), yl, T::ZERO, T::ONE, T::ZERO)
    };
    let nu_wall = nu * r;
    // τ_w/ρ = ν_w U_t/y_P
    let tau = nu_wall * ut / y;
    let u_tau = tau.sqrt_of();
    let eps = eps_vis + w_log * u3 / (kappa * y_len);
    let pk_log = tau * u_star / (kappa * y_len);
    // standard: no production in the sublayer (w_log = 0 there)
    let (prod_weight, pk_wall) = if treatment == COMPOUND {
        (T::ONE - viscous, w_log * pk_log)
    } else {
        (T::ONE, w_log * pk_log)
    };
    KepsWallFace {
        u_tau,
        y_plus: u_tau * y / nu,
        y_star: ys,
        epsilon: eps.max_of(tiny),
        nu_wall,
        prod_weight,
        pk_wall,
        viscous,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_law_newton_converges_over_the_y_plus_range() {
        let (nu, kappa) = (1e-5, 0.41);
        for &yp in &[11.06, 11.5, 30.0, 100.0, 1e3, 1e4, 1e5] {
            // build U from a chosen u_τ = 1 and y = y⁺ ν
            let y = yp * nu;
            let u = (1.0 / kappa) * foehncraft_core::math::ln(LOG_LAW_E * yp);
            let ut = log_law_u_tau(u, y, nu, kappa);
            assert!((ut - 1.0).abs() < 1e-12, "y+ {yp}: {ut}");
        }
    }

    #[test]
    fn automatic_treatment_tends_to_its_limits() {
        let c = SstConstants::default();
        let nu = 1e-5;
        // viscous sublayer: y⁺ = 0.1 with u_τ = 1: U = y⁺ u_τ
        let y = 0.1 * nu;
        let i = WallInputs {
            ut: 0.1,
            y,
            nu,
            ks: 0.0,
            u_tau_old: 0.0,
        };
        let r = wall_face(&i, &c, AUTOMATIC);
        assert!((r.u_tau - 1.0).abs() < 1e-6, "{}", r.u_tau);
        assert!((r.nu_wall - nu).abs() < 1e-5 * nu);
        let w_vis = 6.0 * nu / (0.075 * y * y);
        assert!((r.omega - w_vis).abs() < 1e-3 * w_vis);
        assert!(r.blending > 0.999);
        // log layer: y⁺ = 300 on the log law
        let y = 300.0 * nu;
        let u = foehncraft_core::math::ln(LOG_LAW_E * 300.0) / 0.41;
        let i = WallInputs {
            ut: u,
            y,
            nu,
            ks: 0.0,
            u_tau_old: 0.0,
        };
        let r = wall_face(&i, &c, AUTOMATIC);
        // the viscous part adds (νU/y)²/4 relative: U/300 = 0.048, +0.06 %
        let v = u / 300.0;
        let blended = foehncraft_core::math::sqrt(foehncraft_core::math::sqrt(1.0 + v * v));
        assert!((r.u_tau - blended).abs() < 1e-12, "{}", r.u_tau);
        let w_log = r.u_tau / (0.3 * 0.41 * y);
        assert!((r.omega - w_log).abs() < 1e-3 * w_log);
        assert!((r.y_plus - 300.0 * r.u_tau).abs() < 1e-9);
        // low-re: the viscous values
        let r = wall_face(&i, &c, LOW_RE);
        assert_eq!(r.nu_wall, nu);
    }

    #[test]
    fn aupoix_values_of_fully_rough_and_smooth_limits() {
        // ω_w⁺ ~ 300/k_s⁺² · 4k_s⁺/15 ... for small k_s⁺ the smooth limit
        // grows like k_s⁺⁻²; fully rough: ω_w⁺ → 191/k_s⁺ + 300/k_s⁺²·
        let w = aupoix_omega_plus(1000.0_f64);
        let exact = 300.0 / 1e6 / foehncraft_core::math::tanh(15.0 / 4000.0)
            + 0.191 * (1.0 - foehncraft_core::math::exp(-4.0));
        assert!((w - exact).abs() < 1e-14);
        // k_w⁺ = 0 below the onset of roughness effects (k_s⁺ small)
        assert_eq!(aupoix_k_plus(1.0_f64, 0.09), 0.0);
        // fully rough: tanh(k_s⁺/125) = 1, k_w⁺ = tanh(log₁₀(k_s⁺/30))/√β*
        let k = aupoix_k_plus(1e4_f64, 0.09);
        use foehncraft_core::math as m;
        let exact = m::tanh(m::log10(1e4 / 30.0)) / 0.3;
        assert!((k - exact).abs() < 1e-12, "{k} {exact}");
    }

    #[test]
    fn keps_wall_functions_in_the_log_layer() {
        // equilibrium log layer: k = u_τ²/√C_μ, U = u_τ u⁺_log(y⁺): every
        // treatment returns u_τ, ε = u_τ³/(κ y), P = u_τ³/(κ y)
        let c = crate::turb::keps::KepsConstants::standard();
        let nu = 1e-5;
        for &yp in &[30.0, 100.0, 300.0] {
            let y = yp * nu;
            let u = foehncraft_core::math::ln(LOG_LAW_E * yp) / 0.41;
            let i = KepsWallInputs {
                ut: u,
                y,
                nu,
                k: 1.0 / 0.3,
                ks: 0.0,
                cs: 0.5,
            };
            for t in [SCALABLE, STANDARD, COMPOUND] {
                let r = keps_wall_face(&i, &c, t);
                // the compound blend reaches the log law slowly:
                // e^{−1/Γ} = 0.98 at y⁺ = 30
                let tol = if t == COMPOUND { 2.5e-2 } else { 1e-12 };
                assert!((r.u_tau - 1.0).abs() < tol, "{t} {yp}: {}", r.u_tau);
                assert!((r.y_star - yp).abs() < 1e-9 * yp);
                let e = 1.0 / (0.41 * y);
                assert!((r.epsilon - e).abs() < tol * e, "{t}: {} {e}", r.epsilon);
                assert!((r.pk_wall - e).abs() < 2.0 * tol * e, "{t}: {}", r.pk_wall);
                assert!((r.prod_weight - 1.0).abs() < tol);
            }
        }
    }

    #[test]
    fn keps_wall_treatments_in_the_viscous_sublayer() {
        let c = crate::turb::keps::KepsConstants::standard();
        let nu = 1e-5;
        // y⁺ = 0.5 with u_τ = 1, k⁺ = 0.1 y⁺²
        let y = 0.5 * nu;
        let k: f64 = 0.1 * 0.25;
        let i = KepsWallInputs {
            ut: 0.5,
            y,
            nu,
            k,
            ks: 0.0,
            cs: 0.5,
        };
        for t in [STANDARD, COMPOUND] {
            let r = keps_wall_face(&i, &c, t);
            // the viscous wall shear exactly
            assert!((r.nu_wall - nu).abs() < 1e-6 * nu, "{t}: {}", r.nu_wall);
            assert!((r.u_tau - 1.0).abs() < 1e-6);
        }
        let r = keps_wall_face(&i, &c, COMPOUND);
        assert!((r.epsilon - 2.0 * nu * k / (y * y)).abs() < 1e-6 * r.epsilon);
        assert!(r.prod_weight.abs() < 1e-6 && r.pk_wall == 0.0);
        // scalable takes the log law at ỹ* = max(y*, 11.06) (it moves the
        // wall; accurate for y⁺ ≥ 30 only)
        let r = keps_wall_face(&i, &c, SCALABLE);
        let m = foehncraft_core::math::sqrt;
        let ys = m(m(c.c_mu)) * m(k) * y / nu;
        let up = foehncraft_core::math::ln(LOG_LAW_E * Y_PLUS_LOG) / 0.41;
        assert!((r.nu_wall - nu * ys / up).abs() < 1e-12 * nu);
        // no k: still the viscous shear (compound, standard) and finite
        let i0 = KepsWallInputs::<f64> { k: 0.0, ..i };
        for t in [SCALABLE, STANDARD, COMPOUND] {
            let r = keps_wall_face(&i0, &c, t);
            assert!(r.nu_wall.is_finite() && r.epsilon.is_finite() && r.pk_wall.is_finite());
        }
        assert_eq!(keps_wall_face(&i0, &c, COMPOUND).nu_wall, nu);
    }

    #[test]
    fn cebeci_bradshaw_shift_is_continuous() {
        let (cs, kappa) = (0.5, 0.41);
        assert_eq!(roughness_shift(1.0_f64, cs, kappa), 0.0);
        let at = |x: f64| roughness_shift(x, cs, kappa);
        // (to the four digits of the constants 0.811 = ln 2.25)
        assert!(at(2.25 + 1e-9).abs() < 1e-4);
        assert!((at(90.0 - 1e-9) - at(90.0)).abs() < 1e-6);
        // fully rough: u⁺ = (1/κ) ln(E y⁺/(1 + C_s k_s⁺))
        let v = at(1000.0);
        assert!((v - foehncraft_core::math::ln(501.0) / 0.41).abs() < 1e-12);
    }

    #[test]
    fn keps_wall_face_derivatives_by_dual_numbers() {
        use foehncraft_core::Dual;
        let c = crate::turb::keps::KepsConstants::standard();
        let base = KepsWallInputs {
            ut: 5.0,
            y: 2e-4,
            nu: 1.5e-5,
            k: 0.08,
            ks: 1e-4,
            cs: 0.5,
        };
        let lift = |x: f64| Dual::<1>::constant(x);
        for t in [SCALABLE, STANDARD, COMPOUND] {
            let di = KepsWallInputs {
                ut: lift(base.ut),
                y: lift(base.y),
                nu: lift(base.nu),
                k: Dual::<1>::variable(base.k, 0),
                ks: lift(base.ks),
                cs: lift(base.cs),
            };
            let r = keps_wall_face(&di, &c, t);
            let h = 1e-7;
            let at = |k: f64| keps_wall_face(&KepsWallInputs { k, ..base }, &c, t);
            let (p, m) = (at(base.k + h), at(base.k - h));
            for (dv, fd) in [
                (r.epsilon.d[0], (p.epsilon - m.epsilon) / (2.0 * h)),
                (r.nu_wall.d[0], (p.nu_wall - m.nu_wall) / (2.0 * h)),
                (r.pk_wall.d[0], (p.pk_wall - m.pk_wall) / (2.0 * h)),
            ] {
                assert!(
                    (dv - fd).abs() <= 1e-5 * fd.abs().max(1e-10),
                    "{t}: {dv} {fd}"
                );
            }
        }
    }
}

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
}

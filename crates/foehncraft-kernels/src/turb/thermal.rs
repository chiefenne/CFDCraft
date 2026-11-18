//! Thermal wall functions and the turbulent Prandtl number (D135; stage 6
//! specification section 5), generic over [`Real`] like the other turbulence
//! kernels.
//!
//! - [`kader`]: the temperature (or scalar) profile of B. A. Kader, *Temperature
//!   and concentration profiles in fully turbulent boundary layers*, Int. J.
//!   Heat Mass Transfer 24, 1541–1544 (1981), valid for all y⁺ and any Pr:
//!   T⁺ = Pr y⁺ e^{−Γ} + [2.12 ln(1 + y⁺) + β(Pr)] e^{−1/Γ},
//!   β = (3.85 Pr^⅓ − 1.3)² + 2.12 ln Pr, Γ = 0.01 (Pr y⁺)⁴/(1 + 5 Pr³ y⁺);
//! - [`jayatilleke`]: the log-region form of C. L. V. Jayatilleke, *The
//!   influence of Prandtl number and surface roughness on the resistance of
//!   the laminar sublayer to momentum and heat transfer*, Prog. Heat Mass
//!   Transfer 1, 193–329 (1969): T⁺ = Pr_t [u⁺ + P(Pr/Pr_t)],
//!   P = 9.24 [(Pr/Pr_t)^¾ − 1][1 + 0.28 e^{−0.007 Pr/Pr_t}], blended with
//!   the conduction limit T⁺ = Pr y⁺ by the Γ of the Kader form;
//! - [`kays_crawford`]: W. M. Kays, *Turbulent Prandtl number — where are
//!   we?*, J. Heat Transfer 116, 284–295 (1994), the correlation of Kays &
//!   Crawford: 1/Pr_t = 1/(2 Pr_t∞) + C Pe_t/√Pr_t∞
//!   − (C Pe_t)² [1 − exp(−1/(C Pe_t √Pr_t∞))], Pe_t = (ν_t/ν) Pr.
//!
//! The heat-transfer coefficient of a wall face is h_w = ρ c_p u*/T⁺ (for a
//! passive scalar ρ u*/φ⁺); [`wall_conductance`] returns the face
//! diffusivity Γ_w = h_w y_P that replaces the molecular one, so the face
//! conductance of the scalar equation is h_w.

use super::wall::LOG_LAW_E;
use foehncraft_core::Real;

/// Constants of the Kader profile (the 2.12 is 1/κ_T of the log region).
const KADER_SLOPE: f64 = 2.12;
/// The thermal wall law of a scalar.
pub const NONE: u32 = 0;
/// See [`NONE`]: the profile of Kader (all y⁺, any Pr).
pub const KADER: u32 = 1;
/// See [`NONE`]: the log-region form of Jayatilleke.
pub const JAYATILLEKE: u32 = 2;

/// The turbulent Prandtl number: a constant.
pub const PRT_CONSTANT: u32 = 0;
/// The turbulent Prandtl number: Kays & Crawford.
pub const PRT_KAYS_CRAWFORD: u32 = 1;

/// Pr_t∞ and C of Kays & Crawford.
pub const KC_PRT_INF: f64 = 0.85;
/// See [`KC_PRT_INF`].
pub const KC_C: f64 = 0.3;

/// Kármán constant of the log region of Jayatilleke.
const KAPPA: f64 = 0.41;

#[inline(always)]
fn cbrt<T: Real>(x: T) -> T {
    (x.ln_of() / T::from_f64(3.0)).exp_of()
}

/// The blending variable Γ of Kader (1981): zero at the wall.
#[inline(always)]
fn gamma_t<T: Real>(yp: T, pr: T) -> T {
    let a = pr * yp;
    T::from_f64(0.01) * a * a * a * a / (T::ONE + T::from_f64(5.0) * pr * pr * pr * yp)
}

/// e^{−1/Γ} without the overflow at Γ → 0.
#[inline(always)]
fn exp_inv<T: Real>(g: T) -> T {
    let g = g.max_of(T::from_f64(1e-30));
    (-(T::ONE / g)).max_of(T::from_f64(-700.0)).exp_of()
}

/// T⁺ of Kader (1981) for y⁺ ≥ 0 and Pr > 0.
#[inline(always)]
pub fn kader<T: Real>(yp: T, pr: T) -> T {
    let s = T::from_f64(3.85) * cbrt(pr) - T::from_f64(1.3);
    let beta = s * s + T::from_f64(KADER_SLOPE) * pr.ln_of();
    let g = gamma_t(yp, pr);
    pr * yp * (-g).exp_of() + (T::from_f64(KADER_SLOPE) * (T::ONE + yp).ln_of() + beta) * exp_inv(g)
}

/// T⁺ of the log region of Jayatilleke (1969), blended with Pr y⁺ at the
/// wall; `prt` is the turbulent Prandtl number of the log layer.
#[inline(always)]
pub fn jayatilleke<T: Real>(yp: T, pr: T, prt: T) -> T {
    let r = pr / prt;
    let p = T::from_f64(9.24)
        * ((T::from_f64(0.75) * r.ln_of()).exp_of() - T::ONE)
        * (T::ONE + T::from_f64(0.28) * (-(T::from_f64(0.007) * r)).exp_of());
    let u = (T::from_f64(LOG_LAW_E) * yp.max_of(T::ONE / T::from_f64(LOG_LAW_E))).ln_of()
        / T::from_f64(KAPPA);
    let g = gamma_t(yp, pr);
    pr * yp * (-g).exp_of() + prt * (u + p) * exp_inv(g)
}

/// The wall function selected by `law` ([`KADER`], [`JAYATILLEKE`]); with
/// [`NONE`] the conduction limit Pr y⁺.
#[inline(always)]
pub fn t_plus<T: Real>(law: u32, yp: T, pr: T, prt: T) -> T {
    if law == KADER {
        kader(yp, pr)
    } else if law == JAYATILLEKE {
        jayatilleke(yp, pr, prt)
    } else {
        pr * yp
    }
}

/// The turbulent Prandtl number of Kays & Crawford for the ratio
/// ν_t/ν and the molecular Prandtl number `pr`.
#[inline(always)]
pub fn kays_crawford<T: Real>(nut_over_nu: T, pr: T) -> T {
    let pe = (nut_over_nu * pr).max_of(T::ZERO);
    if pe.to_f64() < 1e-12 {
        return T::from_f64(2.0 * KC_PRT_INF);
    }
    let c = T::from_f64(KC_C) * pe;
    let root = T::from_f64(KC_PRT_INF).sqrt_of();
    let inv = T::ONE / T::from_f64(2.0 * KC_PRT_INF) + c / root
        - c * c
            * (T::ONE
                - (-(T::ONE / (c * root)))
                    .max_of(T::from_f64(-700.0))
                    .exp_of());
    T::ONE / inv
}

/// The turbulent Prandtl number of `model` ([`PRT_CONSTANT`]: `prt`).
#[inline(always)]
pub fn prandtl_t<T: Real>(model: u32, prt: T, nut_over_nu: T, pr: T) -> T {
    if model == PRT_KAYS_CRAWFORD {
        kays_crawford(nut_over_nu, pr)
    } else {
        prt
    }
}

/// The wall diffusivity Γ_w = μ y⁺/T⁺ (= ρ u* y_P/T⁺ = h_w y_P) of a wall
/// face; μ/Pr (the molecular value) where y⁺ vanishes.
#[inline(always)]
pub fn wall_conductance<T: Real>(law: u32, mu: T, yp: T, pr: T, prt: T) -> T {
    wall_conductance_at(law, mu, yp, yp, pr, prt)
}

/// Γ_w = μ y⁺/T⁺(ỹ⁺): the conductance h_w = ρ u*/T⁺ with the thermal law
/// taken at ỹ⁺ ≥ y⁺ (the scalable treatment's virtual wall distance).
#[inline(always)]
pub fn wall_conductance_at<T: Real>(law: u32, mu: T, yp: T, yt: T, pr: T, prt: T) -> T {
    if yt.to_f64() < 1e-10 {
        return mu / pr;
    }
    mu * yp / t_plus(law, yt, pr, prt)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() <= 1e-12 * b.abs(), "{a} against {b}");
    }

    #[test]
    fn the_conductance_takes_the_thermal_law_at_the_floor() {
        let (mu, pr, prt) = (1e-3, 0.71, 0.85);
        for law in [KADER, JAYATILLEKE] {
            for yp in [0.5, 5.0, 30.0] {
                close(
                    wall_conductance_at(law, mu, yp, yp, pr, prt),
                    wall_conductance(law, mu, yp, pr, prt),
                );
                let yt = f64::max(yp, 11.06);
                close(
                    wall_conductance_at(law, mu, yp, yt, pr, prt),
                    mu * yp / t_plus(law, yt, pr, prt),
                );
            }
        }
        // the floor lowers the conductance below y⁺ = 11.06 only
        assert!(
            wall_conductance_at(KADER, mu, 5.0, 11.06, pr, prt)
                < wall_conductance(KADER, mu, 5.0, pr, prt)
        );
    }

    #[test]
    fn kader_matches_the_formula_evaluated_apart() {
        // values of the formula in extended precision (Python, 2026-11-18)
        let table: [(f64, [f64; 3]); 4] = [
            (
                0.01,
                [
                    0.009999999999000005,
                    0.09999990000504975,
                    0.9900547815368125,
                ],
            ),
            (
                0.71,
                [0.7093535141584629, 6.087968473786469, 13.604990779252585],
            ),
            (
                7.0,
                [6.902739125663779, 42.819908941411086, 50.68738208125302],
            ),
            (
                100.0,
                [83.79879006578386, 287.97181173595175, 294.1143234202137],
            ),
        ];
        for (pr, v) in table {
            for (yp, want) in [1.0, 10.0, 100.0].into_iter().zip(v) {
                close(kader(yp, pr), want);
            }
        }
    }

    #[test]
    fn limits_of_the_wall_functions() {
        // the conduction limit and the log slope 2.12 per ln y⁺
        close(kader(1e-3, 0.71), 0.71e-3);
        let slope = (kader(2e5, 0.71) - kader(1e5, 0.71)) / 2f64.ln();
        assert!((slope - 2.12).abs() < 1e-3, "{slope}");
        // wall diffusivity: molecular where the wall law is none
        close(wall_conductance(NONE, 2.0, 30.0, 0.71, 0.85), 2.0 / 0.71);
        close(wall_conductance(KADER, 2.0, 0.0, 0.71, 0.85), 2.0 / 0.71);
    }

    #[test]
    fn jayatilleke_matches_the_formula_evaluated_apart() {
        for (yp, want) in [
            (1.0, 0.7093535141584629),
            (10.0, 5.765681315425503),
            (100.0, 13.000639520594774),
        ] {
            close(jayatilleke(yp, 0.71, 0.85), want);
        }
    }

    #[test]
    fn kays_crawford_values_and_limits() {
        for (r, want) in [
            (1.0, 1.2916959734455915),
            (10.0, 0.9189982190150322),
            (100.0, 0.857183326364018),
        ] {
            close(kays_crawford(r, 0.71), want);
        }
        close(kays_crawford(0.0, 0.71), 1.7);
        assert_eq!(prandtl_t(PRT_CONSTANT, 0.85, 50.0, 0.71), 0.85);
    }
}

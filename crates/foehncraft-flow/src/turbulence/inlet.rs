//! Turbulence values at inlets and at the start (stage 6 specification
//! 7.1, 7.3; D137).
//!
//! An inlet takes one of
//!
//! - `turbulence_intensity` I with `turbulence_length_scale` l:
//!   k = (3/2)(I |u|)², ω = k^½/(C_μ^¼ l) with C_μ = β*;
//! - `turbulence_intensity` with `viscosity_ratio` r = ν_t/ν:
//!   k = (3/2)(I |u|)², ω = k/(r ν);
//! - `k` with `omega` (numbers or expressions of x, y, z and, in a
//!   transient case, t; D53), as the TMR cases give them;
//! - nothing: I = 0.05 and r = 10 (TURB009).
//!
//! |u| is the speed of the boundary face velocity of the current iteration
//! (the prescribed one at velocity and mass-flow inlets, the one the flow
//! takes at pressure inlets).

use crate::bc::{BoundaryCondition, InletTurbulence};
use crate::expr::{self, Env, Value, Vars};
use anyhow::Result;

/// Default intensity of an inlet without turbulence values.
pub const DEFAULT_INTENSITY: f64 = 0.05;
/// Default viscosity ratio of an inlet without turbulence values.
pub const DEFAULT_RATIO: f64 = 10.0;

/// How the turbulence of an inlet face is set.
#[derive(Debug, Clone)]
pub(crate) enum InletSpec {
    /// k and ω given.
    Direct { k: Value, omega: Value },
    /// Intensity and length scale.
    Length { i: f64, l: f64 },
    /// Intensity and viscosity ratio.
    Ratio { i: f64, r: f64 },
}

impl InletSpec {
    /// The spec of an inlet's values (validated before).
    pub fn of(t: &InletTurbulence, env: &Env, transient: bool, key: &str) -> Result<Self> {
        let time = if transient { Vars::T } else { Vars::NONE };
        Ok(match (&t.k, &t.omega) {
            (Some(k), Some(w)) => Self::Direct {
                k: expr::value(k, env, Vars::SPACE.with(time), &format!("{key}.k"))?,
                omega: expr::value(w, env, Vars::SPACE.with(time), &format!("{key}.omega"))?,
            },
            _ => {
                let i = t.turbulence_intensity.unwrap_or(DEFAULT_INTENSITY);
                match (t.turbulence_length_scale, t.viscosity_ratio) {
                    (Some(l), _) => Self::Length { i, l },
                    (None, r) => Self::Ratio {
                        i,
                        r: r.unwrap_or(DEFAULT_RATIO),
                    },
                }
            }
        })
    }

    /// True when the values depend on t.
    pub fn time_dependent(&self) -> bool {
        match self {
            Self::Direct { k, omega } => k.uses().time() || omega.uses().time(),
            _ => false,
        }
    }

    /// True when the values depend on the face speed.
    pub fn speed_dependent(&self) -> bool {
        !matches!(self, Self::Direct { .. })
    }

    /// k and ω at a face with the speed `speed`, the centroid `x` and the
    /// time `t`; ν the kinematic viscosity, `beta_star` = C_μ.
    pub fn at(&self, speed: f64, x: [f64; 3], t: f64, nu: f64, beta_star: f64) -> (f64, f64) {
        match self {
            Self::Direct { k, omega } => (k.at(x, t), omega.at(x, t)),
            Self::Length { i, l } => {
                let k = 1.5 * (i * speed) * (i * speed);
                (k, k.sqrt() / (beta_star.sqrt().sqrt() * l))
            }
            Self::Ratio { i, r } => {
                let k = 1.5 * (i * speed) * (i * speed);
                (k, k / (r * nu))
            }
        }
    }
}

/// True for the conditions that take turbulence values at their faces
/// (inlets); outlets take their backflow from the cell (zero gradient).
pub(crate) fn is_inlet(bc: &BoundaryCondition) -> bool {
    use BoundaryCondition as B;
    matches!(
        bc,
        B::VelocityInlet { .. }
            | B::MassFlowInlet { .. }
            | B::StaticPressureInlet { .. }
            | B::TotalPressureInlet { .. }
    )
}

/// The checks of the inlet turbulence of one boundary (BC029, TURB003,
/// TURB009): (level is error, code, message) per finding.
pub(crate) fn check(bc: &BoundaryCondition) -> Vec<(bool, &'static str, String)> {
    let t = bc.turbulence();
    let mut out = Vec::new();
    let direct = t.k.is_some() || t.omega.is_some();
    let scales = t.turbulence_intensity.is_some()
        || t.turbulence_length_scale.is_some()
        || t.viscosity_ratio.is_some();
    if direct && scales {
        out.push((
            true,
            crate::diag::TURB_INLET_CONFLICT,
            format!(
                "{} together: give `k` with `omega`, or `turbulence_intensity` with \
                 `turbulence_length_scale` or `viscosity_ratio`",
                t.given().join(", ")
            ),
        ));
    } else if t.k.is_some() != t.omega.is_some() {
        out.push((
            true,
            crate::diag::TURB_INLET_CONFLICT,
            "`k` and `omega` go together".into(),
        ));
    } else if t.turbulence_length_scale.is_some() && t.viscosity_ratio.is_some() {
        out.push((
            true,
            crate::diag::TURB_INLET_CONFLICT,
            "`turbulence_length_scale` and `viscosity_ratio` together: give one".into(),
        ));
    }
    if let Some(l) = t.turbulence_length_scale
        && !(l > 0.0 && l.is_finite())
    {
        out.push((
            true,
            crate::diag::CASE_RANGE,
            format!("turbulence_length_scale must be > 0, is {l}"),
        ));
    }
    if let Some(i) = t.turbulence_intensity {
        if !(i >= 0.0 && i.is_finite()) {
            out.push((
                true,
                crate::diag::CASE_RANGE,
                format!("turbulence_intensity must be >= 0, is {i}"),
            ));
        } else if i > 0.5 {
            out.push((
                false,
                crate::diag::TURB_INLET_RANGE,
                format!("turbulence_intensity = {i} is above 0.5"),
            ));
        }
    }
    if let Some(r) = t.viscosity_ratio {
        if !(r > 0.0 && r.is_finite()) {
            out.push((
                true,
                crate::diag::CASE_RANGE,
                format!("viscosity_ratio must be > 0, is {r}"),
            ));
        } else if !(1e-2..=1e5).contains(&r) {
            out.push((
                false,
                crate::diag::TURB_INLET_RANGE,
                format!("viscosity_ratio = {r} is outside 1e-2 … 1e5"),
            ));
        }
    }
    if is_inlet(bc) && !direct && !scales {
        out.push((
            false,
            crate::diag::TURB_INLET_DEFAULT,
            format!(
                "no turbulence values: turbulence_intensity = {DEFAULT_INTENSITY} and \
                 viscosity_ratio = {DEFAULT_RATIO}"
            ),
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intensity_with_length_scale_or_viscosity_ratio() {
        // I = 5 %, |u| = 10 m/s: k = 1.5 (0.5)² = 0.375 m²/s²
        let (k, w) = InletSpec::Length { i: 0.05, l: 0.1 }.at(10.0, [0.0; 3], 0.0, 1e-5, 0.09);
        assert!((k - 0.375).abs() < 1e-15);
        // ω = √k / (C_μ^¼ l), C_μ^¼ = 0.09^¼ = 0.5477
        assert!((w - 0.375f64.sqrt() / (0.09f64.powf(0.25) * 0.1)).abs() < 1e-12);
        let (k, w) = InletSpec::Ratio { i: 0.05, r: 10.0 }.at(10.0, [0.0; 3], 0.0, 1e-5, 0.09);
        assert!((k - 0.375).abs() < 1e-15);
        // ν_t = k/ω = r ν
        assert!((k / w - 1e-4).abs() < 1e-18);
    }
}

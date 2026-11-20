//! Turbulence values at inlets and at the start (stage 6 specification
//! 7.1, 7.3; D137).
//!
//! An inlet takes one of
//!
//! - `turbulence_intensity` I with `turbulence_length_scale` l:
//!   k = (3/2)(I |u|)², ω = k^½/(C_μ^¼ l) with C_μ = β*,
//!   ε = C_μ^¾ k^{3/2}/l;
//! - `turbulence_intensity` with `viscosity_ratio` r = ν_t/ν:
//!   k = (3/2)(I |u|)², ω = k/(r ν), ε = C_μ k²/(r ν);
//! - `k` with `omega` or `epsilon` (numbers or expressions of x, y, z
//!   and, in a transient case, t; D53), as the TMR cases give them;
//!   converted with ε = C_μ k ω when the model solves the other one;
//! - the columns `k` with `epsilon` or `omega` (and `zeta`) of the table
//!   of a table profile, interpolated at the face centroid like the
//!   velocity: a profile from a periodic pipe or channel run (D111) feeds
//!   a jet or a duct;
//! - nothing: I = 0.05 and r = 10 (TURB009).
//!
//! |u| is the speed of the boundary face velocity of the current iteration
//! (the prescribed one at velocity and mass-flow inlets, the one the flow
//! takes at pressure inlets).

use crate::bc::{BoundaryCondition, InletTurbulence, ProfileColumns, TableProfile};
use crate::expr::{self, Env, Value, Vars};
use anyhow::Result;

/// Default intensity of an inlet without turbulence values.
pub const DEFAULT_INTENSITY: f64 = 0.05;
/// Default viscosity ratio of an inlet without turbulence values.
pub const DEFAULT_RATIO: f64 = 10.0;

/// How the turbulence of an inlet face is set.
#[derive(Debug, Clone)]
pub(crate) enum InletSpec {
    /// k and ω (`eps` false) or ε (`eps` true) given.
    Direct { k: Value, second: Value, eps: bool },
    /// Intensity and length scale.
    Length { i: f64, l: f64 },
    /// Intensity and viscosity ratio.
    Ratio { i: f64, r: f64 },
    /// The columns of the table of a table profile.
    Table {
        profile: TableProfile,
        k: usize,
        second: usize,
        eps: bool,
        zeta: Option<usize>,
    },
}

impl InletSpec {
    /// The spec of the inlet `bc` (validated before): the turbulence
    /// columns of its table profile, if it has them, else
    /// [`Self::of`].
    pub fn of_inlet(bc: &BoundaryCondition, env: &Env, transient: bool, key: &str) -> Result<Self> {
        if let Some(p) = bc.table_profile()
            && let Ok(c) = p.table.columns()
            && let Some(k) = c.k
        {
            let (second, eps) = match (c.epsilon, c.omega) {
                (Some(e), _) => (e, true),
                (None, Some(w)) => (w, false),
                (None, None) => unreachable!("validated with the table"),
            };
            return Ok(Self::Table {
                profile: p,
                k,
                second,
                eps,
                zeta: c.zeta,
            });
        }
        Self::of(&bc.turbulence(), env, transient, key)
    }

    /// The spec of an inlet's values (validated before).
    pub fn of(t: &InletTurbulence, env: &Env, transient: bool, key: &str) -> Result<Self> {
        let time = if transient { Vars::T } else { Vars::NONE };
        let second = match (&t.omega, &t.epsilon) {
            (Some(w), _) => Some((w, "omega", false)),
            (None, Some(e)) => Some((e, "epsilon", true)),
            _ => None,
        };
        Ok(match (&t.k, second) {
            (Some(k), Some((w, name, eps))) => Self::Direct {
                k: expr::value(k, env, Vars::SPACE.with(time), &format!("{key}.k"))?,
                second: expr::value(w, env, Vars::SPACE.with(time), &format!("{key}.{name}"))?,
                eps,
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
            Self::Direct { k, second, .. } => k.uses().time() || second.uses().time(),
            _ => false,
        }
    }

    /// True when the values depend on the face speed.
    pub fn speed_dependent(&self) -> bool {
        matches!(self, Self::Length { .. } | Self::Ratio { .. })
    }

    /// ζ at the face centroid `x` from the table, if it has the column.
    pub fn zeta_at(&self, x: [f64; 3]) -> Option<f64> {
        match self {
            Self::Table {
                profile,
                zeta: Some(z),
                ..
            } => profile.at(x).map(|v| v[*z]),
            _ => None,
        }
    }

    /// k and ω (`eps` false) or ε (`eps` true) at a face with the speed
    /// `speed`, the centroid `x` and the time `t`; ν the kinematic
    /// viscosity, `c_mu` = C_μ (β* of SST).
    pub fn at(&self, speed: f64, x: [f64; 3], t: f64, nu: f64, c_mu: f64, eps: bool) -> (f64, f64) {
        let (k, w) = match self {
            Self::Table {
                profile,
                k,
                second,
                eps: e,
                ..
            } => {
                let row = profile.at(x).unwrap_or_default();
                let (k, v) = (
                    row.get(*k).copied().unwrap_or(0.0),
                    row.get(*second).copied().unwrap_or(0.0),
                );
                if *e == eps {
                    return (k, v);
                }
                if *e {
                    return (k, v / (c_mu * k).max(1e-300));
                }
                (k, v)
            }
            Self::Direct { k, second, eps: e } => {
                let (k, v) = (k.at(x, t), second.at(x, t));
                if *e == eps {
                    return (k, v);
                }
                if *e {
                    // ε given, ω wanted
                    return (k, v / (c_mu * k).max(1e-300));
                }
                (k, v)
            }
            Self::Length { i, l } => {
                let k = 1.5 * (i * speed) * (i * speed);
                if eps {
                    return (k, c_mu.sqrt().sqrt().powi(3) * k * k.sqrt() / l);
                }
                (k, k.sqrt() / (c_mu.sqrt().sqrt() * l))
            }
            Self::Ratio { i, r } => {
                let k = 1.5 * (i * speed) * (i * speed);
                if eps {
                    return (k, c_mu * k * k / (r * nu));
                }
                (k, k / (r * nu))
            }
        };
        // ω given, ε wanted
        if eps { (k, c_mu * k * w) } else { (k, w) }
    }
}

/// True for the conditions that take turbulence values at their faces
/// (inlets); outlets take the values of their backflow where the flow
/// enters, if given ([`has_backflow_values`]), else the cell's (zero
/// gradient).
/// True for a pressure outlet that gives turbulence values for its
/// backflow (k with ω or ε, or the scales).
pub(crate) fn has_backflow_values(bc: &BoundaryCondition) -> bool {
    let t = bc.turbulence();
    bc.is_pressure_outlet()
        && (t.k.is_some()
            || t.omega.is_some()
            || t.epsilon.is_some()
            || t.turbulence_intensity.is_some()
            || t.turbulence_length_scale.is_some()
            || t.viscosity_ratio.is_some())
}

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
    let direct = t.k.is_some() || t.omega.is_some() || t.epsilon.is_some();
    let scales = t.turbulence_intensity.is_some()
        || t.turbulence_length_scale.is_some()
        || t.viscosity_ratio.is_some();
    let columns: Option<ProfileColumns> = bc
        .table_profile()
        .and_then(|p| p.table.columns().ok())
        .filter(ProfileColumns::has_turbulence);
    if let Some(c) = columns {
        if direct || scales {
            out.push((
                true,
                crate::diag::TURB_INLET_CONFLICT,
                format!(
                    "the table of the profile has the turbulence columns, and the boundary \
                     gives {} as well: give one",
                    t.given().join(", ")
                ),
            ));
        }
        if c.zeta.is_some() && t.zeta.is_some() {
            out.push((
                true,
                crate::diag::TURB_INLET_CONFLICT,
                "the table of the profile has the column `zeta`, and the boundary gives `zeta` \
                 as well: give one"
                    .into(),
            ));
        }
        return out;
    }
    if direct && scales {
        out.push((
            true,
            crate::diag::TURB_INLET_CONFLICT,
            format!(
                "{} together: give `k` with `omega` or `epsilon`, or `turbulence_intensity` \
                 with `turbulence_length_scale` or `viscosity_ratio`",
                t.given().join(", ")
            ),
        ));
    } else if t.omega.is_some() && t.epsilon.is_some() {
        out.push((
            true,
            crate::diag::TURB_INLET_CONFLICT,
            "`omega` and `epsilon` together: give one (with `k`)".into(),
        ));
    } else if t.k.is_some() != (t.omega.is_some() || t.epsilon.is_some()) {
        out.push((
            true,
            crate::diag::TURB_INLET_CONFLICT,
            "`k` goes with `omega` or `epsilon`".into(),
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
    use crate::bc::{ProfileTable, TableCoordinate, VelocityColumns};
    use std::sync::Arc;

    fn columns(text: &str) -> std::result::Result<ProfileColumns, String> {
        ProfileTable::parse(text)?.columns()
    }

    #[test]
    fn a_header_names_the_columns_of_a_profile() {
        // the legacy forms: the header "y, u_n" and four unnamed columns
        let c = columns("y, u_n\n0, 0\n1, 1\n").unwrap();
        assert_eq!(c.velocity, VelocityColumns::Normal(0));
        assert!(!c.has_turbulence());
        let c = columns("0, 1, 0, 0\n1, 2, 0, 0\n").unwrap();
        assert_eq!(c.velocity, VelocityColumns::Components([0, 1, 2]));
        // any order, k with ε or ω, ζ
        let c = columns("r zeta epsilon k u_n\n0 0.5 1e-3 1e-2 1\n1 0.3 3e-3 3e-2 2\n").unwrap();
        assert_eq!(c.velocity, VelocityColumns::Normal(3));
        assert_eq!(
            (c.k, c.epsilon, c.omega, c.zeta),
            (Some(2), Some(1), None, Some(0))
        );
        let c = columns("y,u_x,u_y,u_z,k,omega\n0,1,0,0,1,2\n1,1,0,0,1,2\n").unwrap();
        assert_eq!(c.velocity, VelocityColumns::Components([0, 1, 2]));
        assert_eq!((c.k, c.omega), (Some(3), Some(4)));
        // refusals
        for (text, part) in [
            ("r, u_n, k, eps\n0,1,1,1\n1,1,1,1\n", "`eps`"),
            (
                "r, u_n, u_x, u_y, u_z\n0,1,1,0,0\n1,1,1,0,0\n",
                "one of the two",
            ),
            ("r, k, epsilon\n0,1,1\n1,1,1\n", "one of the two"),
            ("r, u_n, k\n0,1,1\n1,1,1\n", "goes with"),
            ("r, u_n, zeta\n0,1,1\n1,1,1\n", "goes with"),
            (
                "r, u_n, k, epsilon, omega\n0,1,1,1,1\n1,1,1,1,1\n",
                "together",
            ),
            ("r, u_n, k, k\n0,1,1,1\n1,1,1,1\n", "twice"),
            ("r, a, b\n0,1,1\n1,1,1\n", "without a header"),
        ] {
            let e = columns(text).unwrap_err();
            assert!(e.contains(part), "{text:?}: {e}");
        }
        // without a header only 2 or 4 columns
        assert!(
            ProfileTable::parse("0,1,1\n1,1,1\n")
                .unwrap_err()
                .contains("columns")
        );
    }

    #[test]
    fn the_turbulence_columns_are_interpolated_at_the_face() {
        let table = ProfileTable::parse(
            "r, u_n, k, epsilon, zeta\n0, 1, 0.01, 0.001, 0.5\n1, 2, 0.03, 0.003, 0.3\n",
        )
        .unwrap();
        let c = table.columns().unwrap();
        let spec = InletSpec::Table {
            profile: TableProfile {
                table: Arc::new(table),
                coordinate: TableCoordinate::Radius,
                origin: [0.0; 3],
                direction: Some([1.0, 0.0, 0.0]),
            },
            k: c.k.unwrap(),
            second: c.epsilon.unwrap(),
            eps: true,
            zeta: c.zeta,
        };
        assert!(!spec.speed_dependent() && !spec.time_dependent());
        // r = 0.25 at x = (7, 0.15, 0.2)
        let x = [7.0, 0.15, 0.2];
        let (k, e) = spec.at(99.0, x, 0.0, 1e-5, 0.09, true);
        assert!(
            (k - 0.015).abs() < 1e-15 && (e - 0.0015).abs() < 1e-15,
            "{k} {e}"
        );
        let (k, w) = spec.at(99.0, x, 0.0, 1e-5, 0.09, false);
        assert!((k - 0.015).abs() < 1e-15);
        assert!((w - 0.0015 / (0.09 * 0.015)).abs() < 1e-12, "{w}");
        assert!((spec.zeta_at(x).unwrap() - 0.45).abs() < 1e-15);
        // beyond the rows: the end values
        let (k, _) = spec.at(0.0, [0.0, 3.0, 0.0], 0.0, 1e-5, 0.09, true);
        assert_eq!(k, 0.03);
    }

    #[test]
    fn intensity_with_length_scale_or_viscosity_ratio() {
        // I = 5 %, |u| = 10 m/s: k = 1.5 (0.5)² = 0.375 m²/s²
        let (k, w) =
            InletSpec::Length { i: 0.05, l: 0.1 }.at(10.0, [0.0; 3], 0.0, 1e-5, 0.09, false);
        assert!((k - 0.375).abs() < 1e-15);
        // ω = √k / (C_μ^¼ l), C_μ^¼ = 0.09^¼ = 0.5477
        assert!((w - 0.375f64.sqrt() / (0.09f64.powf(0.25) * 0.1)).abs() < 1e-12);
        let (k, w) =
            InletSpec::Ratio { i: 0.05, r: 10.0 }.at(10.0, [0.0; 3], 0.0, 1e-5, 0.09, false);
        assert!((k - 0.375).abs() < 1e-15);
        // ν_t = k/ω = r ν
        assert!((k / w - 1e-4).abs() < 1e-18);
        // the ε forms: ε = C_μ k ω in both
        for spec in [
            InletSpec::Length { i: 0.05, l: 0.1 },
            InletSpec::Ratio { i: 0.05, r: 10.0 },
        ] {
            let (k, w) = spec.at(10.0, [0.0; 3], 0.0, 1e-5, 0.09, false);
            let (_, e) = spec.at(10.0, [0.0; 3], 0.0, 1e-5, 0.09, true);
            assert!((e - 0.09 * k * w).abs() < 1e-12 * e);
        }
    }
}

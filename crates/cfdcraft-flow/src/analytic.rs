//! Exact solutions and reference data of the stage-4 verification and
//! validation cases (D45):
//!
//! - plane Couette flow between two walls, the upper one moving in x;
//! - plane Poiseuille flow between two walls at rest, driven by a pressure
//!   drop in x;
//! - the lid-driven cavity centre-line velocities of Ghia, Ghia & Shin
//!   (1982), tables I and II, Re = 100 and 1000.

use serde::{Deserialize, Serialize};

/// An exact solution, `[analytic]`. The walls are the planes y = `lower`
/// and y = `upper`; the flow is in x.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "kebab-case", deny_unknown_fields)]
pub enum Analytic {
    /// Plane Couette flow: u = U (y − y₀)/(y₁ − y₀), p constant.
    Couette {
        /// y of the wall at rest.
        #[serde(default)]
        lower: f64,
        /// y of the moving wall.
        #[serde(default = "one")]
        upper: f64,
        /// Speed of the upper wall U (m/s).
        #[serde(default = "one")]
        wall_velocity: f64,
        /// The constant pressure (Pa).
        #[serde(default)]
        pressure: f64,
    },
    /// Plane Poiseuille flow: u = 6 U_m η (1 − η), η = (y − y₀)/H,
    /// p = p_out + G (x_out − x) with G = 12 μ U_m / H².
    Poiseuille {
        /// y of the lower wall.
        #[serde(default)]
        lower: f64,
        /// y of the upper wall.
        #[serde(default = "one")]
        upper: f64,
        /// Mean velocity U_m (m/s).
        #[serde(default = "one")]
        mean_velocity: f64,
        /// Dynamic viscosity μ (Pa s), for the pressure gradient (set by the
        /// case from `[fluid]`).
        #[serde(default)]
        viscosity: f64,
        /// x of the outlet.
        #[serde(default = "one")]
        x_outlet: f64,
        /// Pressure at the outlet (Pa).
        #[serde(default)]
        pressure_outlet: f64,
    },
}

fn one() -> f64 {
    1.0
}

impl Analytic {
    /// Exact velocity at `x`.
    pub fn velocity(&self, x: [f64; 3]) -> [f64; 3] {
        match *self {
            Self::Couette {
                lower,
                upper,
                wall_velocity,
                ..
            } => [wall_velocity * (x[1] - lower) / (upper - lower), 0.0, 0.0],
            Self::Poiseuille {
                lower,
                upper,
                mean_velocity,
                ..
            } => {
                let eta = (x[1] - lower) / (upper - lower);
                [6.0 * mean_velocity * eta * (1.0 - eta), 0.0, 0.0]
            }
        }
    }

    /// Exact pressure at `x`.
    pub fn pressure(&self, x: [f64; 3]) -> f64 {
        match *self {
            Self::Couette { pressure, .. } => pressure,
            Self::Poiseuille {
                lower,
                upper,
                mean_velocity,
                viscosity,
                x_outlet,
                pressure_outlet,
            } => {
                let h = upper - lower;
                pressure_outlet + 12.0 * viscosity * mean_velocity / (h * h) * (x_outlet - x[0])
            }
        }
    }

    /// Velocity scale for the error norms (the wall or the largest speed).
    pub fn velocity_scale(&self) -> f64 {
        match *self {
            Self::Couette { wall_velocity, .. } => wall_velocity.abs(),
            Self::Poiseuille { mean_velocity, .. } => 1.5 * mean_velocity.abs(),
        }
    }

    /// Sets the viscosity of the pressure gradient.
    pub fn with_viscosity(mut self, mu: f64) -> Self {
        if let Self::Poiseuille { viscosity, .. } = &mut self {
            *viscosity = mu;
        }
        self
    }
}

/// Ghia et al. (1982), table I: y on the vertical centre line x = 0.5.
pub const GHIA_Y: [f64; 17] = [
    1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719,
    0.1016, 0.0703, 0.0625, 0.0547, 0.0000,
];
/// Table I: u / U_lid at `GHIA_Y`, Re = 100.
pub const GHIA_U_RE100: [f64; 17] = [
    1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090,
    -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.00000,
];
/// Table I: u / U_lid at `GHIA_Y`, Re = 1000.
pub const GHIA_U_RE1000: [f64; 17] = [
    1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.06080, -0.10648,
    -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0.00000,
];
/// Ghia et al. (1982), table II: x on the horizontal centre line y = 0.5.
pub const GHIA_X: [f64; 17] = [
    1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563,
    0.0938, 0.0781, 0.0703, 0.0625, 0.0000,
];
/// Table II: v / U_lid at `GHIA_X`, Re = 100.
pub const GHIA_V_RE100: [f64; 17] = [
    0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454,
    0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000,
];
/// Table II: v / U_lid at `GHIA_X`, Re = 1000.
pub const GHIA_V_RE1000: [f64; 17] = [
    0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550, -0.42665, -0.31966, 0.02526,
    0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.00000,
];

/// Reference data of a case, `[verification] reference`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Reference {
    /// Ghia et al. 1982, Re = 100.
    GhiaRe100,
    /// Ghia et al. 1982, Re = 1000.
    GhiaRe1000,
}

impl Reference {
    /// (y, u) on the vertical and (x, v) on the horizontal centre line.
    pub fn centre_lines(self) -> ([f64; 17], [f64; 17], [f64; 17], [f64; 17]) {
        match self {
            Self::GhiaRe100 => (GHIA_Y, GHIA_U_RE100, GHIA_X, GHIA_V_RE100),
            Self::GhiaRe1000 => (GHIA_Y, GHIA_U_RE1000, GHIA_X, GHIA_V_RE1000),
        }
    }

    /// Reynolds number of the data.
    pub fn reynolds(self) -> f64 {
        match self {
            Self::GhiaRe100 => 100.0,
            Self::GhiaRe1000 => 1000.0,
        }
    }

    /// Label for reports.
    pub fn label(self) -> &'static str {
        match self {
            Self::GhiaRe100 => "Ghia et al. (1982), Re = 100",
            Self::GhiaRe1000 => "Ghia et al. (1982), Re = 1000",
        }
    }
}

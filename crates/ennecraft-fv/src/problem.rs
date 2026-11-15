//! The steady scalar transport problem of stage 2:
//!
//! ∇·(ρ u φ) − ∇·(Γ ∇φ) = s,
//!
//! with a prescribed, constant velocity u (so the discrete mass fluxes are
//! conservative to round-off on every mesh), constant ρ and Γ, a source s and
//! a boundary condition per patch. A manufactured solution (Roache 2002)
//! supplies s and the boundary values of the order-of-accuracy tests.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Kind of a boundary condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum BoundaryKind {
    /// Fixed value (Dirichlet).
    FixedValue,
    /// Fixed normal gradient ∂φ/∂n along the outward normal (Neumann).
    FixedGradient,
    /// Zero normal gradient (also the symmetry condition of a scalar, and the
    /// usual outflow condition).
    ZeroGradient,
}

/// Values of a boundary condition over the patch.
#[derive(Debug, Clone, PartialEq)]
pub enum Profile {
    /// One value.
    Constant {
        /// The value.
        value: f64,
    },
    /// `below` where `x[axis] < at`, else `above`.
    Step {
        /// 0, 1 or 2 for x, y, z.
        axis: usize,
        /// Position of the step.
        at: f64,
        /// Value below.
        below: f64,
        /// Value above.
        above: f64,
    },
    /// A linear field φ₀ + g·x (its normal gradient g·n for a fixed
    /// gradient).
    Linear {
        /// φ₀.
        value: f64,
        /// g.
        gradient: [f64; 3],
    },
    /// The manufactured solution of the problem (its normal gradient for a
    /// fixed gradient).
    Manufactured,
}

/// Boundary condition of one patch. In the case file a flat table:
/// `type`, `profile` (`constant`, `step`, `linear`, `manufactured`) and the
/// keys of the profile (`value`; `axis`, `at`, `below`, `above`;
/// `value`, `gradient`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "BoundaryToml", into = "BoundaryToml")]
pub struct Boundary {
    /// Kind.
    pub kind: BoundaryKind,
    /// Values (fixed value or fixed gradient); ignored for zero gradient.
    pub profile: Profile,
    /// Copy of the manufactured solution, filled in by the problem.
    pub exact: Option<Manufactured>,
}

/// Profile names of the case file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "kebab-case")]
pub enum ProfileName {
    /// [`Profile::Constant`].
    #[default]
    Constant,
    /// [`Profile::Step`].
    Step,
    /// [`Profile::Linear`].
    Linear,
    /// [`Profile::Manufactured`].
    Manufactured,
}

/// The case-file form of [`Boundary`].
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BoundaryToml {
    /// Kind.
    #[serde(rename = "type")]
    pub kind: Option<BoundaryKind>,
    /// Profile.
    #[serde(default)]
    pub profile: ProfileName,
    /// Value (constant, linear).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub value: Option<f64>,
    /// Gradient (linear).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gradient: Option<[f64; 3]>,
    /// Axis (step).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub axis: Option<usize>,
    /// Position (step).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub at: Option<f64>,
    /// Value below (step).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub below: Option<f64>,
    /// Value above (step).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub above: Option<f64>,
}

impl TryFrom<BoundaryToml> for Boundary {
    type Error = String;
    fn try_from(t: BoundaryToml) -> Result<Self, String> {
        let kind = t
            .kind
            .ok_or("`type` is required (fixed-value, fixed-gradient, zero-gradient)")?;
        let need = |v: Option<f64>, k: &str| v.ok_or(format!("profile needs `{k}`"));
        let profile = match t.profile {
            ProfileName::Constant => Profile::Constant {
                value: if kind == BoundaryKind::ZeroGradient {
                    0.0
                } else {
                    need(t.value, "value")?
                },
            },
            ProfileName::Step => {
                let axis = t.axis.ok_or("profile `step` needs `axis` (0, 1 or 2)")?;
                if axis > 2 {
                    return Err("`axis` must be 0, 1 or 2".into());
                }
                Profile::Step {
                    axis,
                    at: need(t.at, "at")?,
                    below: need(t.below, "below")?,
                    above: need(t.above, "above")?,
                }
            }
            ProfileName::Linear => Profile::Linear {
                value: t.value.unwrap_or(0.0),
                gradient: t.gradient.ok_or("profile `linear` needs `gradient`")?,
            },
            ProfileName::Manufactured => Profile::Manufactured,
        };
        Ok(Self {
            kind,
            profile,
            exact: None,
        })
    }
}

impl From<Boundary> for BoundaryToml {
    fn from(b: Boundary) -> Self {
        let mut t = BoundaryToml {
            kind: Some(b.kind),
            ..Default::default()
        };
        match b.profile {
            Profile::Constant { value } => t.value = Some(value),
            Profile::Step {
                axis,
                at,
                below,
                above,
            } => {
                t.profile = ProfileName::Step;
                (t.axis, t.at, t.below, t.above) = (Some(axis), Some(at), Some(below), Some(above));
            }
            Profile::Linear { value, gradient } => {
                t.profile = ProfileName::Linear;
                (t.value, t.gradient) = (Some(value), Some(gradient));
            }
            Profile::Manufactured => t.profile = ProfileName::Manufactured,
        }
        t
    }
}

fn zero_profile() -> Profile {
    Profile::Constant { value: 0.0 }
}

impl Boundary {
    /// Fixed value.
    pub fn value(v: f64) -> Self {
        Self {
            kind: BoundaryKind::FixedValue,
            profile: Profile::Constant { value: v },
            exact: None,
        }
    }

    /// Zero gradient.
    pub fn zero_gradient() -> Self {
        Self {
            kind: BoundaryKind::ZeroGradient,
            profile: zero_profile(),
            exact: None,
        }
    }

    /// Value of the condition at face centroid `x` with outward unit
    /// normal `n`: the value, or the normal gradient.
    pub fn value_at(&self, x: [f64; 3], n: [f64; 3]) -> f64 {
        let gradient = self.kind == BoundaryKind::FixedGradient;
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        if self.kind == BoundaryKind::ZeroGradient {
            return 0.0;
        }
        match &self.profile {
            Profile::Constant { value } => *value,
            Profile::Step {
                axis,
                at,
                below,
                above,
            } => {
                if x[*axis] < *at {
                    *below
                } else {
                    *above
                }
            }
            Profile::Linear { value, gradient: g } => {
                if gradient {
                    dot(*g, n)
                } else {
                    value + dot(*g, x)
                }
            }
            Profile::Manufactured => self.exact.as_ref().map_or(0.0, |m| {
                if gradient {
                    dot(m.grad(x), n)
                } else {
                    m.phi(x)
                }
            }),
        }
    }
}

/// Manufactured solution φ = φ₀ + A sin(k_x x + p_x) cos(k_y y + p_y)
/// cos(k_z z + p_z); k_z = 0 gives a 2D field.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Manufactured {
    /// Mean value φ₀.
    pub offset: f64,
    /// Amplitude A.
    pub amplitude: f64,
    /// Wave numbers (1/m).
    pub wavenumber: [f64; 3],
    /// Phases (rad).
    pub phase: [f64; 3],
}

impl Default for Manufactured {
    fn default() -> Self {
        let pi = std::f64::consts::PI;
        Self {
            offset: 1.0,
            amplitude: 1.0,
            wavenumber: [pi, pi, 0.5 * pi],
            phase: [0.5, 0.3, 0.2],
        }
    }
}

impl Manufactured {
    /// The same solution without z dependence (for 2D meshes).
    pub fn planar(mut self) -> Self {
        self.wavenumber[2] = 0.0;
        self.phase[2] = 0.0;
        self
    }

    /// φ(x).
    pub fn phi(&self, x: [f64; 3]) -> f64 {
        let [kx, ky, kz] = self.wavenumber;
        let [px, py, pz] = self.phase;
        self.offset
            + self.amplitude
                * (kx * x[0] + px).sin()
                * (ky * x[1] + py).cos()
                * (kz * x[2] + pz).cos()
    }

    /// ∇φ(x).
    pub fn grad(&self, x: [f64; 3]) -> [f64; 3] {
        let [kx, ky, kz] = self.wavenumber;
        let [px, py, pz] = self.phase;
        let (sx, cx) = (kx * x[0] + px).sin_cos();
        let (sy, cy) = (ky * x[1] + py).sin_cos();
        let (sz, cz) = (kz * x[2] + pz).sin_cos();
        let a = self.amplitude;
        [
            a * kx * cx * cy * cz,
            -a * ky * sx * sy * cz,
            -a * kz * sx * cy * sz,
        ]
    }

    /// ∇²φ(x).
    pub fn laplacian(&self, x: [f64; 3]) -> f64 {
        let [kx, ky, kz] = self.wavenumber;
        -(kx * kx + ky * ky + kz * kz) * (self.phi(x) - self.offset)
    }
}

/// The scalar problem.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Problem {
    /// Density ρ (kg/m³).
    pub density: f64,
    /// Diffusion coefficient Γ (kg/(m s)).
    pub diffusivity: f64,
    /// Velocity u (m/s), constant.
    pub velocity: [f64; 3],
    /// Constant volumetric source (units of φ times kg/(m³ s)), added to the
    /// manufactured source.
    pub source: f64,
    /// Manufactured solution: sets the source and the `manufactured`
    /// boundary values, and the error norms of the report.
    pub manufactured: Option<Manufactured>,
    /// Boundary condition per patch name.
    pub boundary: BTreeMap<String, Boundary>,
}

impl Default for Problem {
    fn default() -> Self {
        Self {
            density: 1.0,
            diffusivity: 0.05,
            velocity: [1.0, 0.6, 0.3],
            source: 0.0,
            manufactured: None,
            boundary: BTreeMap::new(),
        }
    }
}

impl Problem {
    /// Copies the manufactured solution into the boundary conditions that
    /// use it; call after reading a case.
    pub fn link(&mut self) {
        for b in self.boundary.values_mut() {
            b.exact = self.manufactured;
        }
    }

    /// A manufactured-solution problem with fixed values on the given
    /// patches and zero gradient (symmetry) on the others; with `planar` the
    /// solution has no z dependence.
    pub fn manufactured(fixed: &[&str], symmetric: &[&str], planar: bool) -> Self {
        let m = if planar {
            Manufactured::default().planar()
        } else {
            Manufactured::default()
        };
        let mut p = Self {
            manufactured: Some(m),
            ..Default::default()
        };
        if planar {
            p.velocity[2] = 0.0;
        }
        for name in fixed {
            p.boundary.insert(
                name.to_string(),
                Boundary {
                    kind: BoundaryKind::FixedValue,
                    profile: Profile::Manufactured,
                    exact: None,
                },
            );
        }
        for name in symmetric {
            p.boundary
                .insert(name.to_string(), Boundary::zero_gradient());
        }
        p.link();
        p
    }

    /// Velocity at `x`.
    pub fn velocity(&self, _x: [f64; 3]) -> [f64; 3] {
        self.velocity
    }

    /// Source density s(x).
    pub fn source(&self, x: [f64; 3]) -> f64 {
        let mut s = self.source;
        if let Some(m) = &self.manufactured {
            let g = m.grad(x);
            let u = self.velocity(x);
            s += self.density * (u[0] * g[0] + u[1] * g[1] + u[2] * g[2])
                - self.diffusivity * m.laplacian(x);
        }
        s
    }
}

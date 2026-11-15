//! Boundary conditions of the incompressible flow solver (D29) and the check
//! of their combination before a run (D30).
//!
//! Every patch of the mesh takes one condition, `[boundary.<patch>]` with a
//! `type`:
//!
//! | type | velocity | pressure | mass flux |
//! |---|---|---|---|
//! | `wall` | `velocity` (default 0; a moving wall keeps only its tangential part) | extrapolated | 0 |
//! | `velocity-inlet` | `velocity`, `normal_velocity`, `profile = "analytic"` or `profile = "parabolic"` (with `max_velocity`, `profile_axes`) | extrapolated | ρ u_b·S |
//! | `mass-flow-inlet` | normal, from `mass_flow` (kg/s) spread by area | extrapolated | prescribed |
//! | `static-pressure-inlet` | normal to the face, from the flux | `pressure` | Rhie–Chow |
//! | `total-pressure-inlet` | normal to the face, from the flux | p₀ − ½ρu_n² from `total_pressure` | Rhie–Chow |
//! | `static-pressure-outlet` | extrapolated | `pressure` | Rhie–Chow |
//! | `total-pressure-outlet` | extrapolated | p₀ − ½ρu_n² from `total_pressure` | Rhie–Chow |
//! | `symmetry` | no shear stress, zero normal component (implicit per component) | extrapolated | 0 |
//! | `periodic` | not available before stage 6 (needs the periodic face pairing) | | |
//!
//! Static and total pressure conditions combine freely, for example a static
//! pressure inlet with a total pressure outlet. Turbulence values at inlets
//! (`turbulence_intensity`, `turbulence_length_scale`, `viscosity_ratio`) are
//! accepted for the turbulence models of stage 6 and reported as unused in
//! laminar flow.

use crate::analytic::Analytic;
use anyhow::Result;
use inncraft_kernels::flow::{
    FLUX_FIXED, FLUX_PRESSURE, MOM_FIXED, MOM_OUTFLOW, MOM_SYMMETRY, PB_EXTRAPOLATE, PB_FIXED,
    PB_TOTAL, UB_EXTRAPOLATE, UB_FIXED, UB_FROM_FLUX, UB_SYMMETRY,
};
use inncraft_mesh::geometry::{dot, norm, scale, sub};
use inncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;

/// Inlet profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum InletProfile {
    /// The velocity of the `[analytic]` solution at the face centroid (at
    /// the time of the step in a transient case).
    Analytic,
    /// A parabola normal to the face, u_n = U_max Π_a 4 η_a (1 − η_a) with
    /// η_a = (x_a − x_a,min)/(x_a,max − x_a,min) over the patch extent along
    /// each axis a of `profile_axes` (developed laminar inflow; two axes for
    /// a duct).
    Parabolic,
}

/// A coordinate axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Axis {
    /// x.
    X,
    /// y.
    Y,
    /// z.
    Z,
}

/// Bounding box of the nodes of a patch.
pub type Extent = ([f64; 3], [f64; 3]);

/// The bounding box of the nodes of patch `p`.
pub fn patch_extent(mesh: &Mesh, p: &inncraft_mesh::Patch) -> Extent {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for f in p.range() {
        for &n in mesh.face(f) {
            let x = mesh.points[n as usize];
            for i in 0..3 {
                lo[i] = lo[i].min(x[i]);
                hi[i] = hi[i].max(x[i]);
            }
        }
    }
    (lo, hi)
}

/// Turbulence values at an inlet (used from stage 6 on).
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct InletTurbulence {
    /// Turbulence intensity (fraction, e.g. 0.05).
    pub turbulence_intensity: Option<f64>,
    /// Turbulence length scale (m).
    pub turbulence_length_scale: Option<f64>,
    /// Ratio of the eddy viscosity to the molecular one.
    pub viscosity_ratio: Option<f64>,
}

impl InletTurbulence {
    fn given(&self) -> Vec<&'static str> {
        let mut v = Vec::new();
        if self.turbulence_intensity.is_some() {
            v.push("turbulence_intensity");
        }
        if self.turbulence_length_scale.is_some() {
            v.push("turbulence_length_scale");
        }
        if self.viscosity_ratio.is_some() {
            v.push("viscosity_ratio");
        }
        v
    }
}

/// One boundary condition, `[boundary.<patch>]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "kebab-case", deny_unknown_fields)]
pub enum BoundaryCondition {
    /// No-slip wall, at rest or moving tangentially.
    Wall {
        /// Wall velocity (m/s).
        #[serde(default)]
        velocity: [f64; 3],
    },
    /// Prescribed velocity.
    VelocityInlet {
        /// Velocity vector (m/s).
        #[serde(default)]
        velocity: Option<[f64; 3]>,
        /// Speed normal to the face, into the domain (m/s).
        #[serde(default)]
        normal_velocity: Option<f64>,
        /// Profile instead of a constant.
        #[serde(default)]
        profile: Option<InletProfile>,
        /// Largest speed of the parabolic profile (m/s).
        #[serde(default)]
        max_velocity: Option<f64>,
        /// Axes across which the parabolic profile varies.
        #[serde(default)]
        profile_axes: Vec<Axis>,
        /// Turbulence values (stage 6).
        #[serde(default)]
        turbulence_intensity: Option<f64>,
        /// Turbulence length scale (stage 6).
        #[serde(default)]
        turbulence_length_scale: Option<f64>,
        /// Viscosity ratio (stage 6).
        #[serde(default)]
        viscosity_ratio: Option<f64>,
    },
    /// Prescribed mass flow (kg/s) into the domain, normal to the faces,
    /// spread by face area.
    MassFlowInlet {
        /// Mass flow into the domain (kg/s).
        mass_flow: f64,
        /// Turbulence values (stage 6).
        #[serde(default)]
        turbulence_intensity: Option<f64>,
        /// Turbulence length scale (stage 6).
        #[serde(default)]
        turbulence_length_scale: Option<f64>,
        /// Viscosity ratio (stage 6).
        #[serde(default)]
        viscosity_ratio: Option<f64>,
    },
    /// Prescribed static pressure, inflow normal to the face.
    StaticPressureInlet {
        /// Static pressure (Pa).
        pressure: f64,
        /// Turbulence values (stage 6).
        #[serde(default)]
        turbulence_intensity: Option<f64>,
        /// Turbulence length scale (stage 6).
        #[serde(default)]
        turbulence_length_scale: Option<f64>,
        /// Viscosity ratio (stage 6).
        #[serde(default)]
        viscosity_ratio: Option<f64>,
    },
    /// Prescribed total pressure, inflow normal to the face.
    TotalPressureInlet {
        /// Total pressure p₀ = p + ½ρ|u|² (Pa).
        total_pressure: f64,
        /// Turbulence values (stage 6).
        #[serde(default)]
        turbulence_intensity: Option<f64>,
        /// Turbulence length scale (stage 6).
        #[serde(default)]
        turbulence_length_scale: Option<f64>,
        /// Viscosity ratio (stage 6).
        #[serde(default)]
        viscosity_ratio: Option<f64>,
    },
    /// Prescribed static pressure, velocity extrapolated.
    StaticPressureOutlet {
        /// Static pressure (Pa).
        #[serde(default)]
        pressure: f64,
    },
    /// Prescribed total pressure, velocity extrapolated.
    TotalPressureOutlet {
        /// Total pressure (Pa).
        total_pressure: f64,
    },
    /// Symmetry plane (also the two side planes of a 2D case).
    Symmetry,
    /// Periodic pair (stage 6).
    Periodic {
        /// The partner patch.
        #[serde(default)]
        partner: Option<String>,
    },
}

impl BoundaryCondition {
    /// The `type` name.
    pub fn type_name(&self) -> &'static str {
        match self {
            Self::Wall { .. } => "wall",
            Self::VelocityInlet { .. } => "velocity-inlet",
            Self::MassFlowInlet { .. } => "mass-flow-inlet",
            Self::StaticPressureInlet { .. } => "static-pressure-inlet",
            Self::TotalPressureInlet { .. } => "total-pressure-inlet",
            Self::StaticPressureOutlet { .. } => "static-pressure-outlet",
            Self::TotalPressureOutlet { .. } => "total-pressure-outlet",
            Self::Symmetry => "symmetry",
            Self::Periodic { .. } => "periodic",
        }
    }

    /// True for the conditions that fix the pressure.
    pub fn is_pressure(&self) -> bool {
        matches!(
            self,
            Self::StaticPressureInlet { .. }
                | Self::TotalPressureInlet { .. }
                | Self::StaticPressureOutlet { .. }
                | Self::TotalPressureOutlet { .. }
        )
    }

    /// True for the pressure conditions meant for outflow.
    pub fn is_pressure_outlet(&self) -> bool {
        matches!(
            self,
            Self::StaticPressureOutlet { .. } | Self::TotalPressureOutlet { .. }
        )
    }

    fn turbulence(&self) -> InletTurbulence {
        match *self {
            Self::VelocityInlet {
                turbulence_intensity,
                turbulence_length_scale,
                viscosity_ratio,
                ..
            }
            | Self::MassFlowInlet {
                turbulence_intensity,
                turbulence_length_scale,
                viscosity_ratio,
                ..
            }
            | Self::StaticPressureInlet {
                turbulence_intensity,
                turbulence_length_scale,
                viscosity_ratio,
                ..
            }
            | Self::TotalPressureInlet {
                turbulence_intensity,
                turbulence_length_scale,
                viscosity_ratio,
                ..
            } => InletTurbulence {
                turbulence_intensity,
                turbulence_length_scale,
                viscosity_ratio,
            },
            _ => InletTurbulence::default(),
        }
    }
}

/// Severity of a finding of [`check`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Level {
    /// Information (for example how the pressure level is fixed).
    Note,
    /// Allowed but probably not intended.
    Warning,
    /// The run cannot start.
    Error,
}

/// One finding of [`check`]: what is wrong, where, and what to change.
#[derive(Debug, Clone, Serialize)]
pub struct Finding {
    /// Severity.
    pub level: Level,
    /// The boundary it concerns, if one.
    pub patch: Option<String>,
    /// What is wrong and what to change.
    pub message: String,
}

impl fmt::Display for Finding {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let l = match self.level {
            Level::Note => "note",
            Level::Warning => "warning",
            Level::Error => "error",
        };
        match &self.patch {
            Some(p) => write!(f, "{l}: boundary `{p}`: {}", self.message),
            None => write!(f, "{l}: {}", self.message),
        }
    }
}

/// Levenshtein distance, for the "did you mean" hints.
pub fn edit_distance(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut cur = vec![i + 1; b.len() + 1];
        for (j, &cb) in b.iter().enumerate() {
            cur[j + 1] = (prev[j] + usize::from(ca != cb))
                .min(prev[j + 1] + 1)
                .min(cur[j] + 1);
        }
        prev = cur;
    }
    prev[b.len()]
}

/// The candidate closest to `word`, if it is close enough to be a typo.
pub fn did_you_mean<'a>(
    word: &str,
    candidates: impl IntoIterator<Item = &'a str>,
) -> Option<&'a str> {
    candidates
        .into_iter()
        .map(|c| (edit_distance(word, c), c))
        .filter(|&(d, c)| d <= 2.max(c.len() / 3))
        .min_by_key(|&(d, _)| d)
        .map(|(_, c)| c)
}

fn hint(word: &str, candidates: &[&str]) -> String {
    match did_you_mean(word, candidates.iter().copied()) {
        Some(c) => format!(" (did you mean `{c}`?)"),
        None => String::new(),
    }
}

/// Mass flow (kg/s) out of the domain through the faces of every patch whose
/// condition prescribes the flux, for the initial velocity of the
/// conditions (walls and symmetry planes give 0).
fn prescribed_outflow(
    mesh: &Mesh,
    geo: &Geometry,
    bcs: &BTreeMap<String, BoundaryCondition>,
    rho: f64,
    analytic: Option<&Analytic>,
) -> (BTreeMap<String, f64>, f64) {
    let mut out = BTreeMap::new();
    let mut gross = 0.0;
    for p in &mesh.patches {
        let Some(bc) = bcs.get(&p.name) else { continue };
        let area: f64 = p.range().map(|f| norm(geo.face_area[f])).sum();
        let ext = patch_extent(mesh, p);
        let q = match bc {
            BoundaryCondition::VelocityInlet { .. } => p
                .range()
                .map(|f| {
                    let q = rho
                        * dot(
                            inlet_velocity(bc, geo, f, analytic, &ext).unwrap_or([0.0; 3]),
                            geo.face_area[f],
                        );
                    gross += q.abs();
                    q
                })
                .sum(),
            BoundaryCondition::MassFlowInlet { mass_flow, .. } => {
                gross += mass_flow.abs();
                if area > 0.0 { -mass_flow } else { 0.0 }
            }
            _ => continue,
        };
        out.insert(p.name.clone(), q);
    }
    (out, gross)
}

/// The prescribed velocity of face `f` of a velocity inlet whose patch has
/// the extent `ext` (None for other conditions).
pub fn inlet_velocity(
    bc: &BoundaryCondition,
    geo: &Geometry,
    f: usize,
    analytic: Option<&Analytic>,
    ext: &Extent,
) -> Option<[f64; 3]> {
    let BoundaryCondition::VelocityInlet {
        velocity,
        normal_velocity,
        profile,
        max_velocity,
        profile_axes,
        ..
    } = bc
    else {
        return None;
    };
    if let Some(InletProfile::Analytic) = profile {
        return analytic.map(|a| a.velocity(geo.face_centroid[f]));
    }
    if let Some(InletProfile::Parabolic) = profile {
        let x = geo.face_centroid[f];
        let mut un = (*max_velocity)?;
        for a in profile_axes {
            let i = *a as usize;
            let (lo, hi) = (ext.0[i], ext.1[i]);
            let eta = if hi > lo {
                (x[i] - lo) / (hi - lo)
            } else {
                0.5
            };
            un *= 4.0 * eta * (1.0 - eta);
        }
        let s = geo.face_area[f];
        return Some(scale(s, -un / norm(s)));
    }
    if let Some(v) = velocity {
        return Some(*v);
    }
    let s = geo.face_area[f];
    let n = scale(s, 1.0 / norm(s));
    normal_velocity.map(|un| scale(n, -un))
}

/// Checks the combination of boundary conditions on the whole domain (D30)
/// before the run. Errors stop the run; warnings and notes go to the log.
/// Every finding names the boundary and says what to change.
pub fn check(
    mesh: &Mesh,
    geo: &Geometry,
    bcs: &BTreeMap<String, BoundaryCondition>,
    rho: f64,
    analytic: Option<&Analytic>,
    laminar: bool,
) -> Vec<Finding> {
    let mut out = Vec::new();
    let mut push = |level, patch: Option<&str>, message: String| {
        out.push(Finding {
            level,
            patch: patch.map(str::to_string),
            message,
        })
    };
    let names: Vec<&str> = mesh.patches.iter().map(|p| p.name.as_str()).collect();
    // every patch has a condition and every condition a patch
    for p in &mesh.patches {
        if !bcs.contains_key(&p.name) {
            push(
                Level::Error,
                Some(&p.name),
                format!(
                    "no boundary condition: add a table [boundary.{}] with a `type` \
                     (wall, velocity-inlet, mass-flow-inlet, static-pressure-inlet, \
                     total-pressure-inlet, static-pressure-outlet, total-pressure-outlet, symmetry)",
                    p.name
                ),
            );
        }
    }
    for name in bcs.keys() {
        if mesh.patch(name).is_none() {
            push(
                Level::Error,
                Some(name),
                format!(
                    "the mesh has no such patch{}; its patches are {}",
                    hint(name, &names),
                    names.join(", ")
                ),
            );
        }
    }
    let present: Vec<(&str, &BoundaryCondition)> = mesh
        .patches
        .iter()
        .filter_map(|p| bcs.get(&p.name).map(|b| (p.name.as_str(), b)))
        .collect();
    // conditions that are not available yet or need more input
    for &(name, bc) in &present {
        match bc {
            BoundaryCondition::Periodic { .. } => push(
                Level::Error,
                Some(name),
                "periodic boundaries come with stage 6 (they need the periodic face \
                 pairing of the mesh); use a velocity-inlet with `profile = \"analytic\"` and a \
                 static-pressure-outlet for a developed flow until then"
                    .into(),
            ),
            BoundaryCondition::VelocityInlet {
                velocity,
                normal_velocity,
                profile,
                max_velocity,
                profile_axes,
                ..
            } => {
                let given = usize::from(velocity.is_some())
                    + usize::from(normal_velocity.is_some())
                    + usize::from(profile.is_some());
                if given != 1 {
                    push(
                        Level::Error,
                        Some(name),
                        format!(
                            "a velocity-inlet takes exactly one of `velocity`, `normal_velocity` \
                             and `profile` ({given} given)"
                        ),
                    );
                }
                if *profile == Some(InletProfile::Analytic) && analytic.is_none() {
                    push(
                        Level::Error,
                        Some(name),
                        "`profile = \"analytic\"` needs an [analytic] table".into(),
                    );
                }
                if *profile == Some(InletProfile::Parabolic)
                    && (max_velocity.is_none() || profile_axes.is_empty())
                {
                    push(
                        Level::Error,
                        Some(name),
                        "`profile = \"parabolic\"` needs `max_velocity` (m/s) and \
                         `profile_axes` (for example [\"y\"])"
                            .into(),
                    );
                }
                if *profile != Some(InletProfile::Parabolic)
                    && (max_velocity.is_some() || !profile_axes.is_empty())
                {
                    push(
                        Level::Warning,
                        Some(name),
                        "`max_velocity` and `profile_axes` are only used with \
                         `profile = \"parabolic\"`"
                            .into(),
                    );
                }
            }
            BoundaryCondition::MassFlowInlet { mass_flow, .. } if *mass_flow < 0.0 => push(
                Level::Warning,
                Some(name),
                format!(
                    "mass_flow = {mass_flow} kg/s is negative, so the flow leaves the domain \
                     here; for an outlet use static-pressure-outlet"
                ),
            ),
            _ => {}
        }
        let t = bc.turbulence().given();
        if laminar && !t.is_empty() {
            push(
                Level::Warning,
                Some(name),
                format!(
                    "{} unused: the flow is laminar (turbulence models come with stage 6)",
                    t.join(", ")
                ),
            );
        }
    }
    // wall velocity normal to the wall, velocity inlets with outflow,
    // non-planar symmetry planes
    for p in &mesh.patches {
        let Some(bc) = bcs.get(&p.name) else { continue };
        match bc {
            BoundaryCondition::Wall { velocity } if norm(*velocity) > 0.0 => {
                let worst = p
                    .range()
                    .map(|f| {
                        let s = geo.face_area[f];
                        dot(*velocity, s).abs() / norm(s)
                    })
                    .fold(0.0f64, f64::max);
                if worst > 1e-9 * norm(*velocity) {
                    push(
                        Level::Warning,
                        Some(&p.name),
                        format!(
                            "the wall velocity has a component of up to {worst:.3e} m/s normal to \
                             the wall; only the tangential part is used (a wall lets no mass through)"
                        ),
                    );
                }
            }
            BoundaryCondition::VelocityInlet { profile, .. } => {
                let ext = patch_extent(mesh, p);
                let out_faces = p
                    .range()
                    .filter(|&f| {
                        inlet_velocity(bc, geo, f, analytic, &ext)
                            .is_some_and(|u| dot(u, geo.face_area[f]) > 0.0)
                    })
                    .count();
                if out_faces > 0 && *profile == Some(InletProfile::Analytic) {
                    push(
                        Level::Note,
                        Some(&p.name),
                        format!(
                            "the exact solution leaves the domain on {out_faces} of {} faces \
                             (the velocity is prescribed there as well)",
                            p.len
                        ),
                    );
                } else if out_faces > 0 {
                    push(
                        Level::Warning,
                        Some(&p.name),
                        format!(
                            "the velocity points out of the domain on {out_faces} of {} faces; \
                             for an outlet use static-pressure-outlet, or reverse the velocity",
                            p.len
                        ),
                    );
                }
            }
            BoundaryCondition::Symmetry if p.len > 0 => {
                let n0 = {
                    let s = geo.face_area[p.start];
                    scale(s, 1.0 / norm(s))
                };
                let worst = p
                    .range()
                    .map(|f| {
                        let s = geo.face_area[f];
                        let n = scale(s, 1.0 / norm(s));
                        norm(sub(n, n0))
                    })
                    .fold(0.0f64, f64::max);
                if worst > 1e-6 {
                    push(
                        Level::Warning,
                        Some(&p.name),
                        format!(
                            "the symmetry patch is not planar (normals differ by up to {:.2}°); \
                             the condition is applied face by face",
                            2.0 * (worst / 2.0).asin().to_degrees()
                        ),
                    );
                }
            }
            _ => {}
        }
    }
    // the domain as a whole: pressure level, mass balance, flow paths
    let pressure: Vec<&str> = present
        .iter()
        .filter(|(_, b)| b.is_pressure())
        .map(|(n, _)| *n)
        .collect();
    let outlets: Vec<&str> = present
        .iter()
        .filter(|(_, b)| b.is_pressure_outlet())
        .map(|(n, _)| *n)
        .collect();
    // the scale of the balance is the flow through the faces, so that a
    // boundary that the flow enters and leaves (an exact solution) does not
    // make round-off look like an imbalance
    let (q, scale_q) = prescribed_outflow(mesh, geo, bcs, rho, analytic);
    let net: f64 = q.values().sum();
    let inflow: Vec<&str> = q
        .iter()
        .filter(|&(_, &v)| v < 0.0)
        .map(|(n, _)| n.as_str())
        .collect();
    if pressure.is_empty() {
        if scale_q > 0.0 && net.abs() > 1e-9 * scale_q {
            push(
                Level::Error,
                None,
                format!(
                    "the prescribed boundaries let a net {:.6e} kg/s {} the domain and there is no \
                     pressure boundary to balance it: make one boundary a static-pressure-outlet \
                     (or a total-pressure-outlet), or balance the inflow and outflow",
                    net.abs(),
                    if net < 0.0 { "into" } else { "out of" }
                ),
            );
        } else {
            push(
                Level::Note,
                None,
                "no pressure boundary: the pressure level is fixed in the reference cell \
                 ([pressure] reference_point, reference_value)"
                    .into(),
            );
        }
    } else if outlets.is_empty() && !inflow.is_empty() {
        push(
            Level::Error,
            None,
            format!(
                "the flow entering through {} can only leave through the pressure inlet(s) {}: \
                 the flow rate is over-specified on this flow path; change one of {} to \
                 static-pressure-outlet",
                inflow
                    .iter()
                    .map(|n| format!("`{n}`"))
                    .collect::<Vec<_>>()
                    .join(", "),
                pressure
                    .iter()
                    .map(|n| format!("`{n}`"))
                    .collect::<Vec<_>>()
                    .join(", "),
                pressure
                    .iter()
                    .map(|n| format!("`{n}`"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        );
    } else if inflow.is_empty() {
        let inlets: Vec<&str> = present
            .iter()
            .filter(|(_, b)| b.is_pressure() && !b.is_pressure_outlet())
            .map(|(n, _)| *n)
            .collect();
        let moving = present.iter().any(
            |(_, b)| matches!(b, BoundaryCondition::Wall { velocity } if velocity.iter().any(|v| *v != 0.0)),
        );
        if inlets.is_empty() && outlets.len() == 1 && !moving {
            push(
                Level::Warning,
                Some(outlets[0]),
                "the only open boundary is this outlet and no boundary drives a flow; unless a \
                 body force acts, the result is the fluid at rest (add an inlet or a moving wall)"
                    .into(),
            );
        } else {
            push(
                Level::Note,
                None,
                "no prescribed inflow: the flow rate follows from the pressure boundaries".into(),
            );
        }
    }
    out.sort_by_key(|f| std::cmp::Reverse(f.level));
    out
}

/// Per boundary face arrays of the conditions, in boundary-face order
/// (face index − number of interior faces).
#[derive(Debug, Clone, Default)]
pub struct BoundaryArrays {
    /// Momentum kind (`MOM_*`).
    pub mom_kind: Vec<u32>,
    /// Flux kind (`FLUX_*`).
    pub flux_kind: Vec<u32>,
    /// Update rule of the face velocity (`UB_*`).
    pub ub_rule: Vec<u32>,
    /// Update rule of the face pressure (`PB_*`).
    pub pb_rule: Vec<u32>,
    /// Prescribed face velocity (and the start value of the others).
    pub ub: [Vec<f64>; 3],
    /// Prescribed pressure or total pressure.
    pub pset: Vec<f64>,
    /// Patch of each boundary face.
    pub patch: Vec<u32>,
    /// 1 where the velocity has a zero normal gradient in the least-squares
    /// gradient (outlets), else 0 (the face value enters).
    pub u_zero_gradient: Vec<f64>,
    /// 1 where the pressure has a zero normal gradient in the least-squares
    /// gradient (walls, symmetry planes), else 0.
    pub p_zero_gradient: Vec<f64>,
    /// Stencil point of the velocity gradient per patch: `ZeroGradient` for
    /// the normal projection of the face centroid (outlets, symmetry
    /// planes), `FixedValue` for the face centroid.
    pub u_points: BTreeMap<String, inncraft_fv::problem::BoundaryKind>,
    /// Stencil point of the pressure gradient per patch (`ZeroGradient` at
    /// walls and symmetry planes).
    pub p_points: BTreeMap<String, inncraft_fv::problem::BoundaryKind>,
}

impl BoundaryArrays {
    /// True when some face fixes the pressure.
    pub fn has_pressure(&self) -> bool {
        self.flux_kind.contains(&FLUX_PRESSURE)
    }
}

/// The per-face arrays of the conditions (after [`check`] found no error).
pub fn arrays(
    mesh: &Mesh,
    geo: &Geometry,
    bcs: &BTreeMap<String, BoundaryCondition>,
    rho: f64,
    analytic: Option<&Analytic>,
) -> Result<BoundaryArrays> {
    let ni = mesh.n_interior_faces();
    let nb = mesh.n_boundary_faces();
    let mut a = BoundaryArrays {
        mom_kind: vec![MOM_FIXED; nb],
        flux_kind: vec![FLUX_FIXED; nb],
        ub_rule: vec![UB_FIXED; nb],
        pb_rule: vec![PB_EXTRAPOLATE; nb],
        ub: [vec![0.0; nb], vec![0.0; nb], vec![0.0; nb]],
        pset: vec![0.0; nb],
        patch: vec![0; nb],
        u_zero_gradient: vec![0.0; nb],
        p_zero_gradient: vec![0.0; nb],
        u_points: BTreeMap::new(),
        p_points: BTreeMap::new(),
    };
    use inncraft_fv::problem::BoundaryKind as K;
    for (pi, p) in mesh.patches.iter().enumerate() {
        let bc = bcs
            .get(&p.name)
            .ok_or_else(|| anyhow::anyhow!("no boundary condition for `{}`", p.name))?;
        let area: f64 = p.range().map(|f| norm(geo.face_area[f])).sum();
        let ext = patch_extent(mesh, p);
        let outlet = bc.is_pressure_outlet();
        let wall_like = matches!(
            bc,
            BoundaryCondition::Wall { .. } | BoundaryCondition::Symmetry
        );
        let symmetry = matches!(bc, BoundaryCondition::Symmetry);
        a.u_points.insert(
            p.name.clone(),
            if outlet || symmetry {
                K::ZeroGradient
            } else {
                K::FixedValue
            },
        );
        a.p_points.insert(
            p.name.clone(),
            if wall_like {
                K::ZeroGradient
            } else {
                K::FixedValue
            },
        );
        for f in p.range() {
            let fb = f - ni;
            a.patch[fb] = pi as u32;
            a.u_zero_gradient[fb] = if outlet { 1.0 } else { 0.0 };
            a.p_zero_gradient[fb] = if wall_like { 1.0 } else { 0.0 };
            let s = geo.face_area[f];
            let n = scale(s, 1.0 / norm(s));
            let mut set_ub = |v: [f64; 3]| {
                for (ub, vi) in a.ub.iter_mut().zip(v) {
                    ub[fb] = vi;
                }
            };
            match bc {
                BoundaryCondition::Wall { velocity } => {
                    set_ub(sub(*velocity, scale(n, dot(*velocity, n))));
                }
                BoundaryCondition::VelocityInlet { .. } => {
                    set_ub(inlet_velocity(bc, geo, f, analytic, &ext).unwrap_or([0.0; 3]));
                }
                BoundaryCondition::MassFlowInlet { mass_flow, .. } => {
                    set_ub(scale(n, -mass_flow / (rho * area)));
                }
                BoundaryCondition::StaticPressureInlet { pressure, .. } => {
                    a.ub_rule[fb] = UB_FROM_FLUX;
                    a.flux_kind[fb] = FLUX_PRESSURE;
                    a.pb_rule[fb] = PB_FIXED;
                    a.pset[fb] = *pressure;
                }
                BoundaryCondition::TotalPressureInlet { total_pressure, .. } => {
                    a.ub_rule[fb] = UB_FROM_FLUX;
                    a.flux_kind[fb] = FLUX_PRESSURE;
                    a.pb_rule[fb] = PB_TOTAL;
                    a.pset[fb] = *total_pressure;
                }
                BoundaryCondition::StaticPressureOutlet { pressure } => {
                    a.mom_kind[fb] = MOM_OUTFLOW;
                    a.ub_rule[fb] = UB_EXTRAPOLATE;
                    a.flux_kind[fb] = FLUX_PRESSURE;
                    a.pb_rule[fb] = PB_FIXED;
                    a.pset[fb] = *pressure;
                }
                BoundaryCondition::TotalPressureOutlet { total_pressure } => {
                    a.mom_kind[fb] = MOM_OUTFLOW;
                    a.ub_rule[fb] = UB_EXTRAPOLATE;
                    a.flux_kind[fb] = FLUX_PRESSURE;
                    a.pb_rule[fb] = PB_TOTAL;
                    a.pset[fb] = *total_pressure;
                }
                BoundaryCondition::Symmetry => {
                    a.mom_kind[fb] = MOM_SYMMETRY;
                    a.ub_rule[fb] = UB_SYMMETRY;
                }
                BoundaryCondition::Periodic { .. } => {
                    anyhow::bail!("periodic boundaries come with stage 6 (`{}`)", p.name)
                }
            }
        }
    }
    Ok(a)
}

/// Guidance: the sensible combinations and the cases each one suits (also
/// in `docs/cases.md`).
pub const GUIDANCE: &str = "\
| Combination | Suits | Pressure level |
|---|---|---|
| velocity-inlet + static-pressure-outlet (+ walls, symmetry) | internal flows with a known flow rate: channels, pipes, ducts, the Poiseuille case | outlet pressure |
| mass-flow-inlet + static-pressure-outlet | as above when the mass flow is known rather than a velocity profile | outlet pressure |
| total-pressure-inlet + static-pressure-outlet | pressure-driven flows: the flow rate is a result (fans, valves, plenum to plenum) | both |
| static-pressure-inlet + static-pressure-outlet | pressure-driven flows with a known static pressure drop | both |
| static-pressure-inlet + total-pressure-outlet | flows that exhaust into a large volume at a known total pressure | both |
| walls only, one wall moving (+ symmetry) | enclosed flows: the lid-driven cavity | reference cell |
| velocity-inlet + velocity-inlet + static-pressure-outlet | mixing of two streams: the T-junction | outlet pressure |
";

//! Passive scalars (D51) on the generic transport equation of D116
//! (stage 6, PR 6a):
//!
//! ∂(ρφ)/∂t + ∇·(ṁ φ) − φ ∇·ṁ − ∇·(ρΓ ∇φ) = ρ S,
//!
//! carried by the face mass fluxes of the flow solver and solved after its
//! pressure correction: once per outer iteration of a steady run, once per
//! outer iteration of a time step (with BDF1 or BDF2; a Crank–Nicolson run
//! takes BDF2 for the scalars, as the stage 6 specification sets for every
//! bounded variable). A temperature with constant properties is such a
//! scalar, with Γ = ν/Pr = λ/(ρ c_p).
//!
//! Case file:
//!
//! ```toml
//! [scalars.temperature]
//! diffusivity = 2.1e-5          # Γ (m²/s): number or expression of x, y, z
//! initial = 300.0               # start value (number or expression)
//! source = 0.0                  # S ([φ]/s): number or expression of x, y, z, t
//! bounds = [0.0, 1.0]           # optional: bounded BDF2 (D28 (1))
//!
//! [boundary.wall.scalar.temperature]   # per-physics sub-table (D116)
//! flux = 1.0e-3                 # diffusive flux into the domain ([φ] m/s)
//! # value = 310.0               # or a fixed value
//!
//! [solver.scalar]               # one block for every scalar
//! scheme = "linear-upwind"
//! limiter = "venkatakrishnan"
//! ```
//!
//! Boundaries without a sub-table: walls, slip walls, symmetry planes and
//! wedges have zero flux, inlets the value 0, outlets a zero normal
//! gradient (the backflow carries the owner value). An inlet takes `value`,
//! a wall-like boundary `value` or `flux`, an outlet neither.
//!
//! The device work lives in [`Scalars`]: one set of buffers per scalar
//! (value, gradient, limiter, the boundary coefficients φ_f = a φ_P + b of
//! `foehncraft_kernels::fv`, the least-squares weights of its own boundary
//! kinds, the face diffusivity, the source) and one set of work buffers
//! shared by all of them.

use crate::bc::BoundaryCondition;
use crate::device::{Buf, Buf3, Eps2, FlowDevice, ScalarFaceArgs, ScalarRowArgs, ScalarTime};
use crate::diag::{self, Diagnostic};
use crate::expr::{self, Env, Scalar, Value, Vars};
use anyhow::Result;
use foehncraft_fv::disc::Discretisation;
use foehncraft_fv::problem::BoundaryKind;
use foehncraft_fv::solver::{LimiterKind, Preconditioner, Scheme};
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// `[scalars.<name>]`: a passive scalar.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ScalarSection {
    /// Molecular diffusivity Γ (m²/s): ν/Pr for a temperature, ν/Sc for a
    /// concentration; a number or an expression of x, y, z (D53), evaluated
    /// at the face centroids.
    pub diffusivity: Scalar,
    /// Start value: a number or an expression of x, y, z (D53).
    #[serde(default)]
    pub initial: Scalar,
    /// Source per unit mass S (the unit of the scalar per second): a number
    /// or an expression of x, y, z and, in a transient case, t (D53);
    /// explicit in the equation.
    #[serde(default)]
    pub source: Scalar,
    /// Lower and upper bound of the scalar (a mixture fraction: [0, 1]).
    /// The second-order time integration falls back to first order in the
    /// cells where it would leave them (bounded BDF2, D28 (1)), and the log
    /// reports the extremes against them.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bounds: Option<[f64; 2]>,
    /// The exact solution, for the verification (an expression of x, y, z
    /// and, in a transient case, t): the report gives the L2 and largest
    /// error at the cell centroids, a mesh sweep the observed order
    /// (`[verification] scalar_order`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exact: Option<Scalar>,
}

/// `[boundary.<name>.scalar.<scalar>]`: the condition of a scalar on a
/// boundary (D116). One of the two keys.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ScalarBoundary {
    /// Fixed value of the scalar on the boundary: a number or an expression
    /// of x, y, z and, in a transient case, t (D53). Inlets and walls.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub value: Option<Scalar>,
    /// Diffusive flux into the domain per unit area, Γ ∂φ/∂n (the unit of
    /// the scalar times m/s; a temperature: q_w/(ρ c_p) in K m/s): a number
    /// or an expression (D53). Walls, slip walls, symmetry planes, wedges.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub flux: Option<Scalar>,
}

/// `[solver.scalar]`: the transport equations of the passive scalars.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct ScalarSettings {
    /// Convection scheme of the scalars.
    pub scheme: Scheme,
    /// `scheme = "central"`: fraction of the bounded central scheme.
    pub central_blending: f64,
    /// `scheme = "central"`: β of the normalised-variable bounding.
    pub central_bounding: f64,
    /// Limiter of the scalar gradients.
    pub limiter: LimiterKind,
    /// Constant K of the limiter.
    pub limiter_k: f64,
    /// Under-relaxation of the scalars (default 0.9 steady, 1 transient).
    pub relaxation: Option<f64>,
    /// Relative tolerance of each scalar solve (default 1e-3 steady, 1e-8
    /// transient).
    pub tolerance: Option<f64>,
    /// Largest number of BiCGStab iterations per scalar solve.
    pub max_iter: usize,
    /// Solves of each scalar per outer iteration (each with new deferred
    /// corrections).
    pub sweeps: usize,
}

impl Default for ScalarSettings {
    fn default() -> Self {
        Self {
            scheme: Scheme::LinearUpwind,
            central_blending: 1.0,
            central_bounding: 0.1,
            limiter: LimiterKind::Venkatakrishnan,
            limiter_k: 5.0,
            relaxation: None,
            tolerance: None,
            max_iter: 50,
            sweeps: 1,
        }
    }
}

impl ScalarSettings {
    /// α of the scalars.
    pub fn alpha(&self, transient: bool) -> f64 {
        self.relaxation.unwrap_or(if transient { 1.0 } else { 0.9 })
    }

    /// Relative tolerance of a solve.
    pub fn tol(&self, transient: bool) -> f64 {
        self.tolerance
            .unwrap_or(if transient { 1e-8 } else { 1e-3 })
    }

    /// Checks the ranges.
    pub fn validate(&self, transient: bool) -> Result<()> {
        let a = self.alpha(transient);
        crate::need!(
            a > 0.0 && a <= 1.0,
            CASE_RANGE,
            "[solver.scalar] relaxation must be in (0, 1], is {a}"
        );
        crate::need!(
            self.tol(transient) > 0.0,
            CASE_RANGE,
            "[solver.scalar] tolerance must be > 0"
        );
        crate::need!(
            self.max_iter >= 1 && self.sweeps >= 1,
            CASE_RANGE,
            "[solver.scalar] max_iter and sweeps must be >= 1"
        );
        crate::need!(
            self.limiter_k > 0.0,
            CASE_RANGE,
            "[solver.scalar] limiter_k must be > 0"
        );
        crate::need!(
            (0.0..=1.0).contains(&self.central_blending)
                && (0.0..1.0).contains(&self.central_bounding),
            CASE_RANGE,
            "[solver.scalar] central_blending must be in [0, 1] and central_bounding in [0, 1)"
        );
        Ok(())
    }
}

/// Names a scalar cannot take: the fields of the output and the variables
/// of the expressions.
const RESERVED: &[&str] = &[
    "velocity",
    "velocity_x",
    "velocity_y",
    "velocity_z",
    "speed",
    "pressure",
    "velocity_error",
    "wall_distance",
    "x",
    "y",
    "z",
    "t",
];

/// Moves the per-physics sub-tables `scalar` out of the `[boundary.<name>]`
/// tables of a parsed case (the flow conditions keep only their own keys)
/// and returns them per boundary.
pub fn extract_boundary_tables(
    v: &mut toml::Table,
) -> Result<BTreeMap<String, BTreeMap<String, ScalarBoundary>>> {
    let mut out = BTreeMap::new();
    let Some(toml::Value::Table(bs)) = v.get_mut("boundary") else {
        return Ok(out);
    };
    for (name, b) in bs.iter_mut() {
        let toml::Value::Table(t) = b else { continue };
        let Some(s) = t.remove("scalar") else {
            continue;
        };
        let m: BTreeMap<String, ScalarBoundary> = s.try_into().map_err(|e: toml::de::Error| {
            diag::err(
                Diagnostic::new(
                    diag::CASE_BAD_VALUE,
                    format!("[boundary.{name}.scalar]: {}", e.message()),
                )
                .key(format!("boundary.{name}.scalar")),
            )
        })?;
        out.insert(name.clone(), m);
    }
    Ok(out)
}

/// What a boundary of the flow is for a scalar.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    /// Walls, slip walls, symmetry planes, wedges: value or flux, default
    /// zero flux.
    WallLike,
    /// Inlets: value, default 0.
    Inlet,
    /// Outlets: zero normal gradient.
    Outlet,
}

fn side(b: &BoundaryCondition) -> Option<Side> {
    use BoundaryCondition as B;
    Some(match b {
        B::Wall { .. } | B::SlipWall | B::Symmetry | B::Wedge => Side::WallLike,
        B::VelocityInlet { .. }
        | B::MassFlowInlet { .. }
        | B::StaticPressureInlet { .. }
        | B::TotalPressureInlet { .. } => Side::Inlet,
        B::StaticPressureOutlet { .. }
        | B::TotalPressureOutlet { .. }
        | B::ConvectiveOutlet { .. } => Side::Outlet,
        B::Periodic { .. } => return None,
    })
}

/// Checks the `[scalars]` tables and the boundary sub-tables of a case.
pub fn validate(
    scalars: &BTreeMap<String, ScalarSection>,
    tables: &BTreeMap<String, BTreeMap<String, ScalarBoundary>>,
    boundary: &BTreeMap<String, BoundaryCondition>,
    settings: &ScalarSettings,
    env: &Env,
    transient: bool,
) -> Result<()> {
    let time = if transient { Vars::T } else { Vars::NONE };
    for (name, s) in scalars {
        let ok = name
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_alphabetic() || c == '_')
            && name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_');
        crate::need!(
            ok && !RESERVED.contains(&name.as_str())
                && !expr::reserved(name)
                && !crate::derived::Derived::ALL
                    .iter()
                    .flat_map(|d| d.datasets())
                    .any(|(n, _)| n == name)
                && !name.ends_with("_error"),
            CASE_BAD_VALUE,
            "[scalars.{name}]: a scalar name is letters, digits and _ (not starting with a digit) \
             and not one of the output fields (also the derived ones and names ending in \
             `_error`) or expression names ({})",
            RESERVED.join(", ")
        );
        let d = expr::value(
            &s.diffusivity,
            env,
            Vars::SPACE,
            &format!("scalars.{name}.diffusivity"),
        )?;
        if let Some(v) = d.constant() {
            crate::need!(
                v > 0.0 && v.is_finite(),
                CASE_RANGE,
                "[scalars.{name}] diffusivity must be > 0, is {v}"
            );
        }
        expr::value(
            &s.initial,
            env,
            Vars::SPACE,
            &format!("scalars.{name}.initial"),
        )?;
        if let Some(e) = &s.exact {
            expr::value(
                e,
                env,
                Vars::SPACE.with(time),
                &format!("scalars.{name}.exact"),
            )?;
        }
        expr::value(
            &s.source,
            env,
            Vars::SPACE.with(time),
            &format!("scalars.{name}.source"),
        )?;
        if let Some([lo, hi]) = s.bounds {
            crate::need!(
                lo < hi,
                CASE_RANGE,
                "[scalars.{name}] bounds must be [lower, upper] with lower < upper"
            );
        }
    }
    for (b, m) in tables {
        let Some(bc) = boundary.get(b) else {
            return Err(diag::err(
                Diagnostic::new(
                    diag::SCALAR_BOUNDARY,
                    format!(
                        "[boundary.{b}.scalar] belongs to no boundary condition `[boundary.{b}]`"
                    ),
                )
                .key(format!("boundary.{b}.scalar")),
            ));
        };
        for (n, t) in m {
            let key = format!("boundary.{b}.scalar.{n}");
            if !scalars.contains_key(n) {
                let names: Vec<&str> = scalars.keys().map(String::as_str).collect();
                return Err(diag::err(
                    Diagnostic::new(
                        diag::SCALAR_BOUNDARY,
                        format!(
                            "[{key}]: no scalar `{n}` in [scalars] (the case has {})",
                            if names.is_empty() {
                                "none".into()
                            } else {
                                names.join(", ")
                            }
                        ),
                    )
                    .key(key),
                ));
            }
            crate::need!(
                t.value.is_some() != t.flux.is_some(),
                SCALAR_BOUNDARY,
                "[{key}] needs one of `value` and `flux`"
            );
            let sd = side(bc);
            let allowed = match sd {
                Some(Side::WallLike) => true,
                Some(Side::Inlet) => t.value.is_some(),
                _ => false,
            };
            crate::need!(
                allowed,
                SCALAR_BOUNDARY,
                "[{key}]: a {} takes {} (outlets keep a zero normal gradient)",
                bc.type_name(),
                match sd {
                    Some(Side::Inlet) => "a `value`, not a flux",
                    Some(Side::Outlet) => "no scalar condition",
                    _ => "no scalar condition here",
                }
            );
            let vars = Vars::SPACE.with(time);
            if let Some(v) = &t.value {
                expr::value(v, env, vars, &format!("{key}.value"))?;
            }
            if let Some(v) = &t.flux {
                expr::value(v, env, vars, &format!("{key}.flux"))?;
            }
        }
    }
    if !scalars.is_empty() {
        settings.validate(transient)?;
        for (b, bc) in boundary {
            crate::need!(
                side(bc).is_some(),
                BC_NOT_AVAILABLE,
                "[boundary.{b}]: passive scalars on periodic boundaries come with periodic flow (stage 6)"
            );
        }
    }
    Ok(())
}

/// How the boundary faces of one scalar are set.
#[derive(Debug, Clone)]
enum FaceRule {
    /// Zero normal gradient (a = 1, b = 0).
    ZeroGradient,
    /// Fixed value.
    Value(Value),
    /// Fixed diffusive flux into the domain; with d_n/Γ_f of the face.
    Flux(Value, f64),
}

/// One scalar on the device.
#[derive(Debug, Clone)]
struct Field {
    name: String,
    phi: Buf,
    grad: Buf3,
    psi: Buf,
    bc_a: Buf,
    bc_b: Buf,
    lsq: Buf3,
    gamma: Buf,
    su: Option<Buf>,
    eps2: Option<Buf>,
    /// φⁿ, φⁿ⁻¹ and the φⁿ⁻¹ of the step before (taking a step back).
    old: Option<(Buf, Buf, Buf)>,
    /// The last valid value (D63).
    valid: Option<Buf>,
    bounds: (f64, f64),
    rules: Vec<FaceRule>,
    /// The source S, when it is not zero.
    source: Option<Value>,
    time_dependent: bool,
}

/// Device copy of one scalar on the host (restart, output).
#[derive(Debug, Clone, PartialEq)]
pub struct ScalarState {
    /// Name.
    pub name: String,
    /// Cell values.
    pub phi: Vec<f64>,
    /// Transient: the cell values of the level before (φⁿ⁻¹ of the next
    /// step).
    pub old: Option<Vec<f64>>,
}

/// Results of the scalar solves of one outer iteration.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ScalarResiduals {
    /// Scaled residual per scalar (before the solve, as the momentum's).
    pub residual: Vec<f64>,
    /// BiCGStab iterations per scalar (all sweeps).
    pub iterations: Vec<usize>,
    /// Cells where the bounded BDF2 fell back to first order, per scalar.
    pub first_order: Vec<usize>,
}

impl ScalarResiduals {
    /// No scalars.
    pub fn is_empty(&self) -> bool {
        self.residual.is_empty()
    }
}

/// Time terms of the scalars in a step.
#[derive(Debug, Clone, Copy)]
pub struct StepTime {
    /// Time of the new level (s).
    pub time: f64,
    /// Δt (s).
    pub dt: f64,
    /// a₀, a₁, a₂.
    pub a: [f64; 3],
}

/// The passive scalars of a flow on its device.
pub struct Scalars {
    fields: Vec<Field>,
    /// Settings.
    pub settings: ScalarSettings,
    nc: usize,
    cp: Buf,
    cn: Buf,
    bf: Buf,
    val: Buf,
    rhs: Buf,
    diag: Buf,
    order: Buf,
    /// Boundary face centroids (for the expressions).
    xf: Vec<[f64; 3]>,
    /// The last boundary time the time-dependent values were set for.
    t_set: f64,
}

impl Scalars {
    /// The scalars of a case on the device `dev` of a flow with the
    /// discretisation `disc`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        dev: &mut dyn FlowDevice,
        mesh: &Mesh,
        geo: &Geometry,
        disc: &Discretisation,
        rho: f64,
        scalars: &BTreeMap<String, ScalarSection>,
        tables: &BTreeMap<String, BTreeMap<String, ScalarBoundary>>,
        boundary: &BTreeMap<String, BoundaryCondition>,
        settings: &ScalarSettings,
        env: &Env,
        transient: bool,
    ) -> Result<Self> {
        let (nc, nf, ni) = (disc.n_cells, disc.n_faces, disc.n_interior);
        let nb = nf - ni;
        let stored = disc.sell.stored();
        let time = if transient { Vars::T } else { Vars::NONE };
        // boundary of each boundary face
        let mut fb_boundary = vec![0usize; nb];
        for (bi, p) in mesh.boundaries.iter().enumerate() {
            for f in p.range() {
                fb_boundary[f - ni] = bi;
            }
        }
        let xf: Vec<[f64; 3]> = (ni..nf).map(|f| geo.face_centroid[f]).collect();
        // d_n = (x_f − x_P)·n̂ per boundary face
        let dn: Vec<f64> = (ni..nf)
            .map(|f| {
                let (s, _) = geo.normal_area(f);
                let o = mesh.owner[f] as usize;
                let d = [
                    geo.face_centroid[f][0] - geo.cell_centroid[o][0],
                    geo.face_centroid[f][1] - geo.cell_centroid[o][1],
                    geo.face_centroid[f][2] - geo.cell_centroid[o][2],
                ];
                d[0] * s[0] + d[1] * s[1] + d[2] * s[2]
            })
            .collect();
        let mut fields = Vec::new();
        for (name, s) in scalars {
            let dif = expr::value(
                &s.diffusivity,
                env,
                Vars::SPACE,
                &format!("scalars.{name}.diffusivity"),
            )?;
            let gam_face: Vec<f64> = (0..nf).map(|f| dif.at(geo.face_centroid[f], 0.0)).collect();
            if let Some(f) = gam_face.iter().position(|g| !(*g > 0.0 && g.is_finite())) {
                return Err(diag::err(
                    Diagnostic::new(
                        diag::EXPR_VALUE,
                        format!(
                            "[scalars.{name}] diffusivity is {} at the face centroid {:?}; it must be > 0",
                            gam_face[f], geo.face_centroid[f]
                        ),
                    )
                    .key(format!("scalars.{name}.diffusivity")),
                ));
            }
            // the face rules and the stencil points of the gradient
            let mut points = BTreeMap::new();
            let mut rules = Vec::with_capacity(nb);
            let mut time_dependent = false;
            let mut per_boundary: Vec<FaceRule> = Vec::new();
            for p in &mesh.boundaries {
                let bc = boundary.get(&p.name);
                let t = tables.get(&p.name).and_then(|m| m.get(name));
                let key = format!("boundary.{}.scalar.{name}", p.name);
                let rule = match (bc.and_then(side), t) {
                    (_, Some(ScalarBoundary { value: Some(v), .. })) => FaceRule::Value(
                        expr::value(v, env, Vars::SPACE.with(time), &format!("{key}.value"))?,
                    ),
                    (_, Some(ScalarBoundary { flux: Some(v), .. })) => FaceRule::Flux(
                        expr::value(v, env, Vars::SPACE.with(time), &format!("{key}.flux"))?,
                        0.0,
                    ),
                    (Some(Side::Inlet), _) => FaceRule::Value(Value::Number(0.0)),
                    _ => FaceRule::ZeroGradient,
                };
                let kind = match rule {
                    FaceRule::Value(_) => BoundaryKind::FixedValue,
                    _ => BoundaryKind::ZeroGradient,
                };
                match &rule {
                    FaceRule::Value(v) | FaceRule::Flux(v, _) => {
                        time_dependent |= v.uses().time();
                    }
                    FaceRule::ZeroGradient => {}
                }
                points.insert(p.name.clone(), kind);
                per_boundary.push(rule);
            }
            for fb in 0..nb {
                let r = per_boundary[fb_boundary[fb]].clone();
                rules.push(match r {
                    FaceRule::Flux(v, _) => {
                        // b = g d_n with g = q/Γ_f
                        FaceRule::Flux(v, dn[fb] / gam_face[ni + fb])
                    }
                    r => r,
                });
            }
            let lsq = crate::simplec::disc_with(mesh, geo, &points)?.lsq;
            let mut put = |v: &[f64]| -> Result<Buf> {
                let b = dev.alloc(v.len())?;
                dev.upload(b, v)?;
                Ok(b)
            };
            let lsq_b = Buf3 {
                x: put(&lsq.x)?,
                y: put(&lsq.y)?,
                z: put(&lsq.z)?,
            };
            let gamma = put(&gam_face.iter().map(|g| rho * g).collect::<Vec<_>>())?;
            let a: Vec<f64> = rules
                .iter()
                .map(|r| {
                    if matches!(r, FaceRule::Value(_)) {
                        0.0
                    } else {
                        1.0
                    }
                })
                .collect();
            let bc_a = put(&a)?;
            let ini = expr::value(
                &s.initial,
                env,
                Vars::SPACE,
                &format!("scalars.{name}.initial"),
            )?;
            let phi0: Vec<f64> = (0..nc).map(|c| ini.at(geo.cell_centroid[c], 0.0)).collect();
            let phi = put(&phi0)?;
            let bc_b = dev.alloc(nb)?;
            let src = expr::value(
                &s.source,
                env,
                Vars::SPACE.with(time),
                &format!("scalars.{name}.source"),
            )?;
            let source = (src.constant() != Some(0.0)).then_some(src);
            time_dependent |= source.as_ref().is_some_and(|v| v.uses().time());
            let su = match &source {
                Some(_) => Some(dev.alloc(nc)?),
                None => None,
            };
            let eps2 = if settings.limiter == LimiterKind::Venkatakrishnan {
                let e = disc.venkat_eps2(settings.limiter_k);
                let b = dev.alloc(nc)?;
                dev.upload(b, &e)?;
                Some(b)
            } else {
                None
            };
            let psi = dev.alloc(nc)?;
            dev.fill(psi, 1.0)?;
            let old = if transient {
                Some((dev.alloc(nc)?, dev.alloc(nc)?, dev.alloc(nc)?))
            } else {
                None
            };
            let bounds = s
                .bounds
                .map_or((f64::NEG_INFINITY, f64::INFINITY), |[l, h]| (l, h));
            fields.push(Field {
                name: name.clone(),
                phi,
                grad: Buf3 {
                    x: dev.alloc(nc)?,
                    y: dev.alloc(nc)?,
                    z: dev.alloc(nc)?,
                },
                psi,
                bc_a,
                bc_b,
                lsq: lsq_b,
                gamma,
                su,
                eps2,
                old,
                valid: None,
                bounds,
                rules,
                source,
                time_dependent,
            });
        }
        let mut s = Self {
            fields,
            settings: settings.clone(),
            nc,
            cp: dev.alloc(nf)?,
            cn: dev.alloc(nf)?,
            bf: dev.alloc(nf)?,
            val: dev.alloc(stored)?,
            rhs: dev.alloc(nc)?,
            diag: dev.alloc(nc)?,
            order: dev.alloc(nc)?,
            xf,
            t_set: f64::NAN,
        };
        let cells: Vec<[f64; 3]> = geo.cell_centroid[..nc].to_vec();
        s.set_values(dev, &cells, 0.0, true)?;
        Ok(s)
    }

    /// Number of scalars.
    pub fn len(&self) -> usize {
        self.fields.len()
    }

    /// No scalars.
    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    /// The names, in the order of the residuals.
    pub fn names(&self) -> Vec<String> {
        self.fields.iter().map(|f| f.name.clone()).collect()
    }

    /// Boundary values and sources at time `t` (only those that depend on
    /// t, unless `all`).
    fn set_values(
        &mut self,
        dev: &mut dyn FlowDevice,
        cells: &[[f64; 3]],
        t: f64,
        all: bool,
    ) -> Result<()> {
        for f in &self.fields {
            if !all && !f.time_dependent {
                continue;
            }
            let b: Vec<f64> = f
                .rules
                .iter()
                .zip(&self.xf)
                .map(|(r, &x)| match r {
                    FaceRule::ZeroGradient => 0.0,
                    FaceRule::Value(v) => v.at(x, t),
                    FaceRule::Flux(v, dn_g) => v.at(x, t) * dn_g,
                })
                .collect();
            dev.upload(f.bc_b, &b)?;
            if let (Some(src), Some(su)) = (&f.source, f.su) {
                let v: Vec<f64> = cells.iter().map(|&x| src.at(x, t)).collect();
                dev.upload(su, &v)?;
            }
        }
        self.t_set = t;
        Ok(())
    }

    /// The scalar solves of one outer iteration with the fluxes
    /// `mass_flux`; `step` carries the time terms of a time step (and the
    /// new time, for time-dependent values; `cells` are the cell centroids).
    pub fn solve(
        &mut self,
        dev: &mut dyn FlowDevice,
        mass_flux: Buf,
        step: Option<StepTime>,
        cells: &[[f64; 3]],
    ) -> Result<ScalarResiduals> {
        if let Some(st) = step
            && st.time != self.t_set
        {
            self.set_values(dev, cells, st.time, false)?;
        }
        let s = self.settings.clone();
        let tr = step.is_some();
        let alpha = s.alpha(tr);
        let mut out = ScalarResiduals::default();
        for f in &self.fields {
            let mut res = f64::NAN;
            let mut its = 0;
            let mut first = 0;
            for sweep in 0..s.sweeps {
                dev.gradient_with(f.lsq, f.phi, f.bc_a, f.bc_b, f.grad)?;
                if s.scheme != Scheme::Upwind && s.limiter != LimiterKind::None {
                    let eps2 = match (s.limiter, f.eps2) {
                        (LimiterKind::Venkatakrishnan, Some(b)) => Eps2::Cells(b),
                        _ => {
                            let (mn, mx) = dev.min_max(f.phi)?;
                            let e = s.limiter_k * (mx - mn);
                            Eps2::Single(e * e)
                        }
                    };
                    dev.limiter_with(f.phi, f.bc_a, f.bc_b, f.grad, eps2, f.psi)?;
                }
                dev.scalar_faces(&ScalarFaceArgs {
                    bc_a: f.bc_a,
                    bc_b: f.bc_b,
                    mass_flux,
                    gamma: f.gamma,
                    phi: f.phi,
                    grad: f.grad,
                    psi: f.psi,
                    blend: if s.scheme == Scheme::Upwind { 0.0 } else { 1.0 },
                    central: if s.scheme == Scheme::Central {
                        s.central_blending
                    } else {
                        0.0
                    },
                    beta: s.central_bounding,
                    nonorth: 1.0,
                    cp: self.cp,
                    cn: self.cn,
                    bf: self.bf,
                })?;
                let time = match (step, f.old) {
                    (Some(st), Some((o1, o2, _))) => Some(ScalarTime {
                        inv_dt: 1.0 / st.dt,
                        a: st.a,
                        old: o1,
                        old2: o2,
                        bounds: f.bounds,
                    }),
                    _ => None,
                };
                dev.scalar_rows(&ScalarRowArgs {
                    cp: self.cp,
                    cn: self.cn,
                    bf: self.bf,
                    mass_flux,
                    phi: f.phi,
                    su: f.su,
                    sp: None,
                    alpha,
                    time,
                    val: self.val,
                    rhs: self.rhs,
                    diag: self.diag,
                    order: self.order,
                })?;
                let scale = dev.scalar_scale(f.phi, self.rhs, self.diag, alpha)?;
                let (it, r0) = dev.solve(
                    self.val,
                    self.rhs,
                    f.phi,
                    Preconditioner::Jacobi,
                    false,
                    s.tol(tr),
                    s.max_iter,
                )?;
                if sweep == 0 {
                    res = r0 / scale;
                    if time.is_some() {
                        first = dev.sum_sq(self.order)?.round() as usize;
                    }
                }
                its += it;
            }
            out.residual.push(res);
            out.iterations.push(its);
            out.first_order.push(first);
        }
        Ok(out)
    }

    /// Start of a time step: φⁿ⁻¹ ← φⁿ (when an old level exists), φⁿ ← φ;
    /// the φⁿ⁻¹ before is kept for [`step_back`](Self::step_back).
    ///
    /// With `extrapolate = Some(ω)` (BDF2 from the second step, ω =
    /// Δt/Δt_prev) the start value of the step is the linear extrapolation
    /// φ ← φⁿ + ω (φⁿ − φⁿ⁻¹): the deferred correction of the convection
    /// is evaluated from φ, so with one outer iteration per step (PISO) it
    /// would otherwise lag by a whole step and make the scheme first order
    /// in time; extrapolated, its lag error is O(Δt²).
    pub fn step_start(
        &mut self,
        dev: &mut dyn FlowDevice,
        has_old: bool,
        extrapolate: Option<f64>,
    ) -> Result<()> {
        for f in &self.fields {
            let Some((o1, o2, keep)) = f.old else {
                continue;
            };
            dev.copy(o2, keep)?;
            if has_old {
                dev.copy(o1, o2)?;
            }
            dev.copy(f.phi, o1)?;
            if has_old && let Some(w) = extrapolate {
                dev.lincomb(f.phi, 1.0 + w, o1, -w, o2)?;
            }
        }
        Ok(())
    }

    /// Takes the last step back exactly: φ ← φⁿ, φⁿ ← φⁿ⁻¹, φⁿ⁻¹ ← the
    /// one before.
    pub fn step_back(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for f in &self.fields {
            let Some((o1, o2, keep)) = f.old else {
                continue;
            };
            dev.copy(o1, f.phi)?;
            dev.copy(o2, o1)?;
            dev.copy(keep, o2)?;
        }
        Ok(())
    }

    /// Keeps a device copy of the values as the last valid state (D63).
    pub fn keep_valid(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for f in &mut self.fields {
            let v = match f.valid {
                Some(v) => v,
                None => dev.alloc(self.nc)?,
            };
            dev.copy(f.phi, v)?;
            f.valid = Some(v);
        }
        Ok(())
    }

    /// Puts the last valid values back.
    pub fn rewind(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for f in &self.fields {
            if let Some(v) = f.valid {
                dev.copy(v, f.phi)?;
            }
        }
        Ok(())
    }

    /// After a rewind of a time step: the old level of n−1 back into φⁿ.
    pub fn rewind_levels(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for f in &self.fields {
            if let Some((o1, o2, _)) = f.old {
                dev.copy(o2, o1)?;
            }
        }
        Ok(())
    }

    /// Host copies (restart, output).
    pub fn states(&self, dev: &mut dyn FlowDevice) -> Result<Vec<ScalarState>> {
        self.fields
            .iter()
            .map(|f| {
                Ok(ScalarState {
                    name: f.name.clone(),
                    phi: dev.download(f.phi)?,
                    old: match f.old {
                        Some((o1, _, _)) => Some(dev.download(o1)?),
                        None => None,
                    },
                })
            })
            .collect()
    }

    /// Loads host copies (a restart): by name; a scalar of the case that
    /// the file does not hold keeps its start value.
    pub fn load(
        &mut self,
        dev: &mut dyn FlowDevice,
        states: &[ScalarState],
    ) -> Result<Vec<String>> {
        let mut missing = Vec::new();
        for f in &self.fields {
            match states.iter().find(|s| s.name == f.name) {
                Some(s) => {
                    dev.upload(f.phi, &s.phi)?;
                    if let (Some((o1, _, _)), Some(o)) = (f.old, &s.old) {
                        dev.upload(o1, o)?;
                    }
                }
                None => missing.push(f.name.clone()),
            }
        }
        Ok(missing)
    }

    /// Per scalar: the cell values, the boundary coefficients a and b and
    /// the face diffusivity ρΓ.
    #[allow(clippy::type_complexity)]
    fn host_data(
        &self,
        dev: &mut dyn FlowDevice,
    ) -> Result<Vec<(String, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)>> {
        self.fields
            .iter()
            .map(|f| {
                Ok((
                    f.name.clone(),
                    dev.download(f.phi)?,
                    dev.download(f.bc_a)?,
                    dev.download(f.bc_b)?,
                    dev.download(f.gamma)?,
                ))
            })
            .collect()
    }

    /// The minimum and maximum of each scalar.
    pub fn extremes(&self, dev: &mut dyn FlowDevice) -> Result<Vec<(f64, f64)>> {
        self.fields.iter().map(|f| dev.min_max(f.phi)).collect()
    }
}

/// BDF coefficients of the scalars for a step: BDF2 with ω = Δt/Δt_prev
/// when `second`, else BDF1.
pub fn bdf(second: bool, omega: f64) -> [f64; 3] {
    if second {
        [
            (1.0 + 2.0 * omega) / (1.0 + omega),
            -(1.0 + omega),
            omega * omega / (1.0 + omega),
        ]
    } else {
        [1.0, -1.0, 0.0]
    }
}

/// Adds the passive scalars of `case` to the flow solver `fl` (no-op
/// without `[scalars]`); call before loading a restart.
pub fn attach(
    fl: &mut crate::simplec::Flow,
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
) -> Result<()> {
    if case.scalars.is_empty() {
        return Ok(());
    }
    let env = case.env()?;
    let transient = case.time.is_some();
    let s = Scalars::new(
        fl.dev.as_mut(),
        mesh,
        geo,
        &fl.disc,
        fl.phys.rho,
        &case.scalars,
        &case.boundary_scalar,
        &case.boundary,
        &case.solver.scalar,
        &env,
        transient,
    )?;
    fl.add_scalars(s);
    Ok(())
}

/// A scalar on one boundary, for the report.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct ScalarBoundarySummary {
    /// Boundary name.
    pub boundary: String,
    /// Area-weighted mean of the face values.
    pub mean: f64,
    /// Mean weighted with the mass flux through the faces (the mixed-mean
    /// or bulk value), where mass crosses the boundary.
    pub mixed_mean: Option<f64>,
    /// Area-weighted mean diffusive flux into the domain, Γ ∂φ/∂n (the unit
    /// of the scalar times m/s).
    pub flux: f64,
}

/// A scalar in the report of a solve.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct ScalarSummary {
    /// Name.
    pub name: String,
    /// Smallest cell value.
    pub min: f64,
    /// Largest cell value.
    pub max: f64,
    /// `bounds` of the scalar, if set.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bounds: Option<[f64; 2]>,
    /// Mesh size h (as the velocity errors).
    pub h: f64,
    /// Volume-weighted L2 error against `exact` (absolute).
    pub l2: Option<f64>,
    /// Largest error at a cell centroid against `exact` (absolute).
    pub linf: Option<f64>,
    /// Per boundary.
    pub boundaries: Vec<ScalarBoundarySummary>,
}

/// The scalars of `fl` for the report, at time `t` (for the exact
/// solutions); `dim` is 2 for a mesh one cell layer thick.
pub fn summaries(
    fl: &mut crate::simplec::Flow,
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    dim: u8,
    t: f64,
) -> Result<Vec<ScalarSummary>> {
    let Some(sc) = fl.scalars.as_ref() else {
        return Ok(Vec::new());
    };
    let data = sc.host_data(fl.dev.as_mut())?;
    let m = fl.dev.download(fl.f.mass_flux)?;
    let env = case.env()?;
    let transient = case.time.is_some();
    let time = if transient { Vars::T } else { Vars::NONE };
    let nc = mesh.n_cells;
    let ni = mesh.n_interior_faces();
    let rho = fl.phys.rho;
    let vol: f64 = geo.cell_volume[..nc].iter().sum();
    let h = if dim == 2 {
        let (lo, hi) = mesh.bounding_box();
        (vol / (hi[2] - lo[2]) / nc as f64).sqrt()
    } else {
        (vol / nc as f64).cbrt()
    };
    let mut out = Vec::new();
    for (name, phi, a, b, gamma) in data {
        let (mut mn, mut mx) = (f64::INFINITY, f64::NEG_INFINITY);
        for &v in &phi {
            mn = mn.min(v);
            mx = mx.max(v);
        }
        let (mut l2, mut linf) = (None, None);
        if let Some(e) = case.scalars.get(&name).and_then(|s| s.exact.as_ref()) {
            let ex = expr::value(
                e,
                &env,
                Vars::SPACE.with(time),
                &format!("scalars.{name}.exact"),
            )?;
            let (mut e2, mut ei) = (0.0, 0.0f64);
            for ((&v, &x), &vc) in phi[..nc]
                .iter()
                .zip(&geo.cell_centroid[..nc])
                .zip(&geo.cell_volume[..nc])
            {
                let d = v - ex.at(x, t);
                e2 += d * d * vc;
                ei = ei.max(d.abs());
            }
            l2 = Some((e2 / vol).sqrt());
            linf = Some(ei);
        }
        let mut bs = Vec::new();
        for p in &mesh.boundaries {
            let (mut area, mut sum, mut flux, mut mf, mut msum) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for f in p.range() {
                let fb = f - ni;
                let o = mesh.owner[f] as usize;
                let (_, ar) = geo.normal_area(f);
                let pf = a[fb] * phi[o] + b[fb];
                area += ar;
                sum += ar * pf;
                // the diffusive flux into the domain, orthogonal part (the
                // whole flux at a flux face): ρΓ |S|²/(d·S) (φ_f − φ_P)
                let dgeo = fl.disc.dgeo[f];
                flux += gamma[f] / rho * dgeo * (pf - phi[o]);
                mf += m[f].abs();
                msum += m[f].abs() * pf;
            }
            if area > 0.0 {
                bs.push(ScalarBoundarySummary {
                    boundary: p.name.clone(),
                    mean: sum / area,
                    mixed_mean: (mf > 0.0).then(|| msum / mf),
                    flux: flux / area,
                });
            }
        }
        out.push(ScalarSummary {
            bounds: case.scalars.get(&name).and_then(|s| s.bounds),
            name,
            min: mn,
            max: mx,
            h,
            l2,
            linf,
            boundaries: bs,
        });
    }
    Ok(out)
}

/// The error fields `<name>_error` = φ − φ_exact of the scalars with an
/// exact solution, for the field output; `at` as in the fields file (the
/// time of a transient step).
pub fn field_errors(
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    fl: &crate::simplec::Solution,
    at: Option<(usize, Option<f64>)>,
) -> Result<Vec<(String, Vec<f64>)>> {
    let mut out = Vec::new();
    if fl.state.scalars.is_empty() {
        return Ok(out);
    }
    let env = case.env()?;
    let time = if case.time.is_some() {
        Vars::T
    } else {
        Vars::NONE
    };
    let t = at.and_then(|(_, t)| t).unwrap_or(0.0);
    for st in &fl.state.scalars {
        if let Some(e) = case.scalars.get(&st.name).and_then(|s| s.exact.as_ref()) {
            let ex = expr::value(
                e,
                &env,
                Vars::SPACE.with(time),
                &format!("scalars.{}.exact", st.name),
            )?;
            let err = (0..mesh.n_cells)
                .map(|c| st.phi[c] - ex.at(geo.cell_centroid[c], t))
                .collect();
            out.push((format!("{}_error", st.name), err));
        }
    }
    Ok(out)
}

/// One log line of a scalar summary.
pub fn summary_line(s: &ScalarSummary) -> String {
    let mut l = format!("scalar {}: min {:.6e}, max {:.6e}", s.name, s.min, s.max);
    if let Some([lo, hi]) = s.bounds {
        let inside = s.min >= lo && s.max <= hi;
        l.push_str(&format!(
            " ({} the bounds [{lo}, {hi}])",
            if inside { "inside" } else { "OUTSIDE" }
        ));
    }
    if let (Some(a), Some(b)) = (s.l2, s.linf) {
        l.push_str(&format!(", error L2 {a:.4e}, max {b:.4e}"));
    }
    l
}

/// The checks of the scalars over the solves of a report: the largest
/// error and, for a mesh sweep, the observed orders.
pub fn check_sweep(
    rep: &mut crate::case::FlowReport,
    v: &crate::case::VerificationSection,
    levels: usize,
    log: &mut dyn FnMut(&str),
) {
    let Some(first) = rep.solves.first() else {
        return;
    };
    let names: Vec<String> = first
        .scalars
        .iter()
        .filter(|s| s.l2.is_some())
        .map(|s| s.name.clone())
        .collect();
    for name in names {
        let e: Vec<(f64, f64)> = rep
            .solves
            .iter()
            .filter_map(|s| s.scalars.iter().find(|x| x.name == name))
            .filter_map(|x| x.l2.map(|l2| (x.h, l2)))
            .collect();
        if let Some(tol) = v.scalar_max_error
            && let Some(&(_, l2)) = e.last()
            // a NaN error fails too
            && !matches!(
                l2.partial_cmp(&tol),
                Some(std::cmp::Ordering::Less | std::cmp::Ordering::Equal)
            )
        {
            rep.failures.push(format!(
                "scalar {name}: L2 error {l2:.4e} on the finest mesh above {tol}"
            ));
        }
        if levels >= 2 && e.len() == rep.solves.len() {
            let orders: Vec<f64> = e
                .windows(2)
                .map(|w| (w[0].1 / w[1].1).ln() / (w[0].0 / w[1].0).ln())
                .collect();
            log(&format!("scalar {name}: observed orders {orders:.3?}"));
            if let (Some([lo, hi]), Some(&o)) = (v.scalar_order, orders.last())
                && !(lo..=hi).contains(&o)
            {
                rep.failures.push(format!(
                    "scalar {name}: observed order {o:.3} on the finest pair, outside [{lo}, {hi}]"
                ));
            }
            rep.scalar_orders.push((name, orders));
        }
    }
}

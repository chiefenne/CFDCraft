//! Derived fields (D72): quantities computed from the solution when the
//! fields are written (`[output] derived = [...]`), or after the run from a
//! restart file (`foehn post derive`), with the same code, so both give the
//! same numbers.
//!
//! Cell fields, from the cell velocity gradient G with Gᵢⱼ = ∂uᵢ/∂xⱼ, its
//! symmetric part S = ½ (G + Gᵀ) (rate of strain) and its antisymmetric
//! part Ω = ½ (G − Gᵀ) (rate of rotation):
//!
//! | name | fields | definition | unit |
//! |---|---|---|---|
//! | `vorticity` | `vorticity_x`, `_y`, `_z`, `vorticity_magnitude` | ω = ∇ × u | 1/s |
//! | `q-criterion` | `q_criterion` | Q = ½ (‖Ω‖² − ‖S‖²); vortices where Q > 0 (Hunt, Wray and Moin 1988) | 1/s² |
//! | `lambda2` | `lambda2` | the middle eigenvalue of S² + Ω²; vortices where λ₂ < 0 (Jeong and Hussain 1995) | 1/s² |
//! | `enstrophy` | `enstrophy` | ½ ω·ω | 1/s² |
//! | `strain-rate` | `strain_rate` | √(2 S:S), the shear rate | 1/s |
//! | `pressure-coefficient` | `pressure_coefficient` | (p − p_ref) / (½ ρ_ref U_ref²) | - |
//! | `total-pressure` | `total_pressure` | p + ½ ρ u·u | Pa |
//! | `wall-distance` | `wall_distance` | the exact distance of the cell centroid to the nearest no-slip wall face (D48, a BVH closest-point search, D108) | m |
//! | `turbulence-features` | `feature_lambda1` … `feature_lambda5`, `feature_wall_reynolds`, `feature_time_ratio`, `feature_viscosity_ratio`, `feature_shielding`, `feature_q`, `feature_pressure_gradient`, `feature_curvature` | the invariant features of D166 (below); needs a turbulence model | - |
//!
//! Boundary face fields, on the boundary blocks of the file (zero on boundaries
//! that are not walls), from the same wall values as `wall_<boundary>.csv`:
//!
//! | name | fields | definition | unit |
//! |---|---|---|---|
//! | `wall-shear-stress` | `wall_shear_stress_x`, `_y`, `_z`, `wall_shear_stress_magnitude` | the tangential viscous stress of the fluid on the wall | Pa |
//! | `y-plus` | `y_plus` | y⁺ = d u_τ / ν with u_τ = √(τ_w / ρ) and d the wall distance of the cell centroid | - |
//!
//! The features of `turbulence-features` (D166; the definitions of the
//! export of D106), with the turbulence time scale τ = 1/(β* ω) = k/ε,
//! the turbulence length L = √k/(β* ω), the wall distance d of the model
//! and |S| = √(2 S:S):
//!
//! - λ₁ … λ₅: tr Ŝ², tr Ŵ², tr Ŝ³, tr Ŵ²Ŝ, tr Ŵ²Ŝ² of Ŝ = τS, Ŵ = τΩ
//!   (Pope 1975);
//! - wall Reynolds number d√k/ν (based on k, not on the cell);
//! - time-scale ratio τ |S| (turbulent over mean time scale);
//! - viscosity ratio ν_t/ν;
//! - shielding f_d = 1 − tanh((8 r_d)³), r_d = (ν_t + ν)/(‖G‖ κ² d²)
//!   (delayed DES, Spalart et al. 2006): 0 in attached boundary layers;
//! - Q = (‖Ω‖² − ‖S‖²)/(‖Ω‖² + ‖S‖²), in [−1, 1];
//! - pressure gradient along the flow u·∇p/(|u| |∇p|), in [−1, 1]
//!   (velocity in the frame of the walls);
//! - streamline curvature |u × (G u)| L/|u|³.
//!
//! Every quotient has its denominator floored at 10⁻³⁰⁰ (zero for a fluid at
//! rest). ‖·‖ is the Frobenius norm. The eigenvalues of λ₂ come from Jacobi
//! rotations, exact to round-off also where two of them are equal (as in
//! every 2D vortex). Turbulence quantities follow with the turbulence
//! models.

use crate::bc::BoundaryCondition;
use crate::case::FlowCase;
use crate::simplec::Solution;
use anyhow::Result;
use foehncraft_mesh::geometry::{dot, norm, sub};
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};

/// A derived field of `[output] derived`.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Serialize,
    Deserialize,
    schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Derived {
    /// ω = ∇ × u, components and magnitude (1/s).
    Vorticity,
    /// Q = ½ (‖Ω‖² − ‖S‖²) (1/s²).
    QCriterion,
    /// Middle eigenvalue of S² + Ω² (1/s²).
    Lambda2,
    /// ½ ω·ω (1/s²).
    Enstrophy,
    /// √(2 S:S) (1/s).
    StrainRate,
    /// (p − p_ref) / (½ ρ_ref U_ref²) (-).
    PressureCoefficient,
    /// p + ½ ρ u·u (Pa).
    TotalPressure,
    /// Exact distance of the cell centroid to the nearest wall (m).
    WallDistance,
    /// Wall shear stress on the walls, components and magnitude (Pa).
    WallShearStress,
    /// y⁺ of the wall-adjacent cells, on the walls (-).
    YPlus,
    /// The invariant features of D166 (-).
    TurbulenceFeatures,
}

impl Derived {
    /// Every derived field, in the order they are written.
    pub const ALL: [Derived; 11] = [
        Derived::Vorticity,
        Derived::QCriterion,
        Derived::Lambda2,
        Derived::Enstrophy,
        Derived::StrainRate,
        Derived::PressureCoefficient,
        Derived::TotalPressure,
        Derived::WallDistance,
        Derived::WallShearStress,
        Derived::YPlus,
        Derived::TurbulenceFeatures,
    ];

    /// The name in the case file.
    pub fn name(self) -> &'static str {
        match self {
            Derived::Vorticity => "vorticity",
            Derived::QCriterion => "q-criterion",
            Derived::Lambda2 => "lambda2",
            Derived::Enstrophy => "enstrophy",
            Derived::StrainRate => "strain-rate",
            Derived::PressureCoefficient => "pressure-coefficient",
            Derived::TotalPressure => "total-pressure",
            Derived::WallDistance => "wall-distance",
            Derived::WallShearStress => "wall-shear-stress",
            Derived::YPlus => "y-plus",
            Derived::TurbulenceFeatures => "turbulence-features",
        }
    }

    /// The field from its name in the case file.
    pub fn parse(s: &str) -> Option<Derived> {
        Self::ALL.into_iter().find(|d| d.name() == s)
    }

    /// The datasets it writes, with their units.
    pub fn datasets(self) -> &'static [(&'static str, &'static str)] {
        match self {
            Derived::Vorticity => &[
                ("vorticity_x", "1/s"),
                ("vorticity_y", "1/s"),
                ("vorticity_z", "1/s"),
                ("vorticity_magnitude", "1/s"),
            ],
            Derived::QCriterion => &[("q_criterion", "1/s²")],
            Derived::Lambda2 => &[("lambda2", "1/s²")],
            Derived::Enstrophy => &[("enstrophy", "1/s²")],
            Derived::StrainRate => &[("strain_rate", "1/s")],
            Derived::PressureCoefficient => &[("pressure_coefficient", "-")],
            Derived::TotalPressure => &[("total_pressure", "Pa")],
            Derived::WallDistance => &[("wall_distance", "m")],
            Derived::WallShearStress => &[
                ("wall_shear_stress_x", "Pa"),
                ("wall_shear_stress_y", "Pa"),
                ("wall_shear_stress_z", "Pa"),
                ("wall_shear_stress_magnitude", "Pa"),
            ],
            Derived::YPlus => &[("y_plus", "-")],
            Derived::TurbulenceFeatures => &[
                ("feature_lambda1", "-"),
                ("feature_lambda2", "-"),
                ("feature_lambda3", "-"),
                ("feature_lambda4", "-"),
                ("feature_lambda5", "-"),
                ("feature_wall_reynolds", "-"),
                ("feature_time_ratio", "-"),
                ("feature_viscosity_ratio", "-"),
                ("feature_shielding", "-"),
                ("feature_q", "-"),
                ("feature_pressure_gradient", "-"),
                ("feature_curvature", "-"),
            ],
        }
    }

    /// True for the fields on the boundary faces (the boundary blocks).
    pub fn on_faces(self) -> bool {
        matches!(self, Derived::WallShearStress | Derived::YPlus)
    }
}

/// `[output.reference]`: the reference state of the pressure coefficient.
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct ReferenceValues {
    /// Reference velocity U_ref (m/s); default `[output.coefficients]
    /// velocity`.
    pub velocity: Option<f64>,
    /// Reference pressure p_ref (Pa).
    pub pressure: f64,
    /// Reference density ρ_ref (kg/m³); default `[fluid] density`.
    pub density: Option<f64>,
}

/// The reference dynamic pressure ½ ρ_ref U_ref² of `case` and p_ref.
fn reference(case: &FlowCase) -> Option<(f64, f64)> {
    let r = &case.output.reference;
    let u = r
        .velocity
        .or(case.output.coefficients.as_ref().map(|c| c.velocity))?;
    let rho = r.density.unwrap_or(case.rho());
    Some((0.5 * rho * u * u, r.pressure))
}

/// Checks `[output] derived` and `[output.reference]`.
pub fn validate(case: &FlowCase) -> Result<()> {
    let r = &case.output.reference;
    if let Some(u) = r.velocity {
        crate::need!(
            u > 0.0 && u.is_finite(),
            CASE_RANGE,
            "[output.reference] velocity must be > 0 (m/s)"
        );
    }
    if let Some(d) = r.density {
        crate::need!(
            d > 0.0 && d.is_finite(),
            CASE_RANGE,
            "[output.reference] density must be > 0 (kg/m³)"
        );
    }
    crate::need!(
        r.pressure.is_finite(),
        CASE_RANGE,
        "[output.reference] pressure must be finite (Pa)"
    );
    validate_list(case, &case.output.derived)
}

/// Checks a list of derived fields against `case` (also for
/// `foehn post derive --fields`).
pub fn validate_list(case: &FlowCase, list: &[Derived]) -> Result<()> {
    if list.contains(&Derived::WallDistance) {
        crate::need!(
            case.boundary
                .values()
                .any(|b| matches!(b, BoundaryCondition::Wall { .. })),
            CASE_INCONSISTENT,
            "derived field `wall-distance` needs at least one wall boundary"
        );
    }
    if list.contains(&Derived::TurbulenceFeatures) {
        crate::need!(
            !case.turbulence.is_laminar(),
            CASE_INCONSISTENT,
            "derived field `turbulence-features` needs a turbulence model ([turbulence] model)"
        );
    }
    if list.contains(&Derived::PressureCoefficient) {
        crate::need!(
            reference(case).is_some(),
            CASE_INCONSISTENT,
            "derived field `pressure-coefficient` needs a reference velocity: set \
             [output.reference] velocity (or [output.coefficients] velocity)"
        );
    }
    Ok(())
}

/// Computed derived fields: cell fields (one value per cell) and boundary
/// face fields (one value per boundary face).
#[derive(Debug, Clone, Default)]
pub struct Fields {
    /// Cell fields, name and values.
    pub cells: Vec<(&'static str, Vec<f64>)>,
    /// Boundary face fields, name and values.
    pub faces: Vec<(&'static str, Vec<f64>)>,
}

impl Fields {
    /// The values of the field `name`.
    pub fn get(&self, name: &str) -> Option<&[f64]> {
        self.cells
            .iter()
            .chain(&self.faces)
            .find(|(n, _)| *n == name)
            .map(|(_, v)| v.as_slice())
    }

    /// The units of the datasets, `name unit; …`.
    pub fn units(list: &[Derived]) -> String {
        let mut s = String::new();
        for d in sorted(list) {
            for (n, u) in d.datasets() {
                s.push_str(&format!("; {n} {u}"));
            }
        }
        s
    }
}

/// `list` without repetitions, in the order of [`Derived::ALL`] (the file
/// does not depend on the order of the list).
pub fn sorted(list: &[Derived]) -> Vec<Derived> {
    Derived::ALL
        .into_iter()
        .filter(|d| list.contains(d))
        .collect()
}

/// The velocity gradient of cell `c`, G[i][j] = ∂uᵢ/∂xⱼ.
fn gradient(fl: &Solution, c: usize) -> [[f64; 3]; 3] {
    [
        fl.grad_u[0].get(c),
        fl.grad_u[1].get(c),
        fl.grad_u[2].get(c),
    ]
}

/// The cell quantities of one velocity gradient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kinematics {
    /// Vorticity ω (1/s).
    pub vorticity: [f64; 3],
    /// Q criterion (1/s²).
    pub q: f64,
    /// λ₂ (1/s²).
    pub lambda2: f64,
    /// Enstrophy ½ ω·ω (1/s²).
    pub enstrophy: f64,
    /// Shear rate √(2 S:S) (1/s).
    pub strain_rate: f64,
}

/// The derived quantities of the velocity gradient `g`, gᵢⱼ = ∂uᵢ/∂xⱼ.
pub fn kinematics(g: [[f64; 3]; 3]) -> Kinematics {
    let mut s = [[0.0; 3]; 3];
    let mut w = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            s[i][j] = 0.5 * (g[i][j] + g[j][i]);
            w[i][j] = 0.5 * (g[i][j] - g[j][i]);
        }
    }
    let (mut ss, mut ww) = (0.0, 0.0);
    for i in 0..3 {
        for j in 0..3 {
            ss += s[i][j] * s[i][j];
            ww += w[i][j] * w[i][j];
        }
    }
    let vorticity = [g[2][1] - g[1][2], g[0][2] - g[2][0], g[1][0] - g[0][1]];
    // S² + Ω², symmetric
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut v = 0.0;
            for k in 0..3 {
                v += s[i][k] * s[k][j] + w[i][k] * w[k][j];
            }
            m[i][j] = v;
        }
    }
    Kinematics {
        vorticity,
        q: 0.5 * (ww - ss),
        lambda2: symmetric_eigenvalues(m)[1],
        enstrophy: 0.5 * dot(vorticity, vorticity),
        strain_rate: (2.0 * ss).sqrt(),
    }
}

/// The features of D166 of one cell (the order of the datasets of
/// [`Derived::TurbulenceFeatures`]): velocity gradient `g`, velocity `u`,
/// pressure gradient `gp`, k, ω, ν_t, wall distance `d`, ν and β*.
#[allow(clippy::too_many_arguments)]
pub fn features(
    g: [[f64; 3]; 3],
    u: [f64; 3],
    gp: [f64; 3],
    k: f64,
    omega: f64,
    nu_t: f64,
    d: f64,
    nu: f64,
    beta_star: f64,
) -> [f64; 12] {
    use foehncraft_kernels::turb::invariants as inv;
    let tiny = 1e-300;
    let tau = 1.0 / (beta_star * omega).max(tiny);
    let l = inv::pope(&g, tau);
    let (s2, w2) = inv::s2_omega2(&g);
    let k = k.max(0.0);
    let sk = k.sqrt();
    let speed = norm(u);
    let a = [
        g[0][0] * u[0] + g[0][1] * u[1] + g[0][2] * u[2],
        g[1][0] * u[0] + g[1][1] * u[1] + g[1][2] * u[2],
        g[2][0] * u[0] + g[2][1] * u[1] + g[2][2] * u[2],
    ];
    let cross = [
        u[1] * a[2] - u[2] * a[1],
        u[2] * a[0] - u[0] * a[2],
        u[0] * a[1] - u[1] * a[0],
    ];
    let length = sk * tau;
    [
        l[0],
        l[1],
        l[2],
        l[3],
        l[4],
        if d.is_finite() { d * sk / nu } else { 1e300 },
        tau * s2.sqrt(),
        nu_t / nu,
        if d.is_finite() {
            inv::shielding(&g, nu_t, nu, d)
        } else {
            1.0
        },
        (w2 - s2) / (w2 + s2).max(tiny),
        dot(u, gp) / (speed * norm(gp)).max(tiny),
        norm(cross) * length / (speed * speed * speed).max(tiny),
    ]
}

/// The eigenvalues of the symmetric matrix `a`, ascending: cyclic Jacobi
/// rotations (Golub and Van Loan, Matrix Computations, section 8.5),
/// accurate to round-off also for repeated eigenvalues, where the closed
/// form of the cubic loses half the digits.
pub fn symmetric_eigenvalues(a: [[f64; 3]; 3]) -> [f64; 3] {
    let mut m = a;
    let scale = (0..3)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .fold(0.0f64, |s, (i, j)| s + a[i][j] * a[i][j])
        .sqrt();
    for _ in 0..32 {
        let off = m[0][1] * m[0][1] + m[0][2] * m[0][2] + m[1][2] * m[1][2];
        if off.sqrt() <= f64::EPSILON * 1e-3 * scale || off == 0.0 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if m[p][q] == 0.0 {
                continue;
            }
            // the rotation that zeroes m[p][q]
            let theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            let mut r = m;
            for k in 0..3 {
                r[k][p] = c * m[k][p] - s * m[k][q];
                r[k][q] = s * m[k][p] + c * m[k][q];
            }
            let mut n = r;
            for k in 0..3 {
                n[p][k] = c * r[p][k] - s * r[q][k];
                n[q][k] = s * r[p][k] + c * r[q][k];
            }
            n[p][q] = 0.0;
            n[q][p] = 0.0;
            m = n;
        }
    }
    let mut e = [m[0][0], m[1][1], m[2][2]];
    e.sort_by(f64::total_cmp);
    e
}

/// Computes the derived fields `list` of the solution `fl` (cell order and
/// boundary face order of `mesh`).
pub fn compute(
    case: &FlowCase,
    list: &[Derived],
    mesh: &Mesh,
    geo: &Geometry,
    fl: &Solution,
) -> Result<Fields> {
    validate_list(case, list)?;
    let list = sorted(list);
    let nc = mesh.n_cells;
    let mut out = Fields::default();
    let cellwise = list.iter().any(|d| {
        matches!(
            d,
            Derived::Vorticity
                | Derived::QCriterion
                | Derived::Lambda2
                | Derived::Enstrophy
                | Derived::StrainRate
        )
    });
    let kin: Vec<Kinematics> = if cellwise {
        (0..nc).map(|c| kinematics(gradient(fl, c))).collect()
    } else {
        Vec::new()
    };
    let rho = fl.phys.rho;
    let nb = mesh.n_boundary_faces();
    let walls = (list.contains(&Derived::WallShearStress) || list.contains(&Derived::YPlus))
        .then(|| wall_values(case, mesh, geo, fl));
    for d in list {
        match d {
            Derived::Vorticity => {
                for (i, n) in ["vorticity_x", "vorticity_y", "vorticity_z"]
                    .into_iter()
                    .enumerate()
                {
                    out.cells
                        .push((n, kin.iter().map(|k| k.vorticity[i]).collect()));
                }
                out.cells.push((
                    "vorticity_magnitude",
                    kin.iter().map(|k| norm(k.vorticity)).collect(),
                ));
            }
            Derived::QCriterion => out
                .cells
                .push(("q_criterion", kin.iter().map(|k| k.q).collect())),
            Derived::Lambda2 => out
                .cells
                .push(("lambda2", kin.iter().map(|k| k.lambda2).collect())),
            Derived::Enstrophy => out
                .cells
                .push(("enstrophy", kin.iter().map(|k| k.enstrophy).collect())),
            Derived::StrainRate => out
                .cells
                .push(("strain_rate", kin.iter().map(|k| k.strain_rate).collect())),
            Derived::PressureCoefficient => {
                let (q, p0) = reference(case).expect("checked by validate_list");
                out.cells.push((
                    "pressure_coefficient",
                    fl.p.iter().map(|&p| (p - p0) / q).collect(),
                ));
            }
            Derived::TotalPressure => out.cells.push((
                "total_pressure",
                (0..nc)
                    .map(|c| {
                        let u = fl.u.get(c);
                        fl.p[c] + 0.5 * rho * dot(u, u)
                    })
                    .collect(),
            )),
            Derived::WallDistance => out
                .cells
                .push(("wall_distance", wall_distance(case, mesh, geo)?)),
            Derived::WallShearStress => {
                let w = walls.as_ref().expect("computed above");
                for (i, n) in [
                    "wall_shear_stress_x",
                    "wall_shear_stress_y",
                    "wall_shear_stress_z",
                ]
                .into_iter()
                .enumerate()
                {
                    out.faces.push((n, w.iter().map(|v| v.0[i]).collect()));
                }
                out.faces.push((
                    "wall_shear_stress_magnitude",
                    w.iter().map(|v| norm(v.0)).collect(),
                ));
            }
            Derived::YPlus => {
                let w = walls.as_ref().expect("computed above");
                out.faces.push(("y_plus", w.iter().map(|v| v.1).collect()));
            }
            Derived::TurbulenceFeatures => {
                let t = fl.turbulence.as_ref().ok_or_else(|| {
                    crate::diag::err(crate::diag::Diagnostic::new(
                        crate::diag::CASE_INCONSISTENT,
                        "derived field `turbulence-features` needs the turbulence fields",
                    ))
                })?;
                // the k-ε family: ω = ε/(C_μ k), so that τ = 1/(β* ω) = k/ε
                let beta_star = case.turbulence.c_mu()?;
                let nu = fl.phys.mu / fl.phys.rho;
                let rows: Vec<[f64; 12]> = (0..nc)
                    .map(|c| {
                        let gp = fl.grad_p.get(c);
                        features(
                            gradient(fl, c),
                            fl.u.get(c),
                            gp,
                            t.k[c],
                            if t.epsilon {
                                t.omega[c] / (beta_star * t.k[c]).max(1e-300)
                            } else {
                                t.omega[c]
                            },
                            t.nu_t[c],
                            t.wall_distance.get(c).copied().unwrap_or(f64::INFINITY),
                            nu,
                            beta_star,
                        )
                    })
                    .collect();
                for (i, (n, _)) in d.datasets().iter().enumerate() {
                    out.cells.push((n, rows.iter().map(|r| r[i]).collect()));
                }
            }
        }
    }
    debug_assert!(out.faces.iter().all(|(_, v)| v.len() == nb));
    Ok(out)
}

/// The exact wall distance of every cell centroid (D48): the closest point
/// on the faces of the wall boundaries of `case`, by the BVH search of
/// D108.
pub fn wall_distance(case: &FlowCase, mesh: &Mesh, geo: &Geometry) -> Result<Vec<f64>> {
    let walls: Vec<usize> = mesh
        .boundaries
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            matches!(
                case.boundary.get(&p.name),
                Some(BoundaryCondition::Wall { .. })
            )
        })
        .map(|(i, _)| i)
        .collect();
    let d = foehncraft_mesh::query::wall_distance(mesh, geo, &walls).ok_or_else(|| {
        crate::diag::err(crate::diag::Diagnostic::new(
            crate::diag::CASE_INCONSISTENT,
            "derived field `wall-distance` needs at least one wall face in the mesh",
        ))
    })?;
    Ok(d.into_iter().map(|c| c.distance).collect())
}

/// What [`from_restart`] wrote.
#[derive(Debug, Clone, Serialize)]
pub struct PostReport {
    /// The file written.
    pub file: std::path::PathBuf,
    /// The restart file read.
    pub restart: std::path::PathBuf,
    /// Cells.
    pub cells: usize,
    /// Outer iterations of the solution.
    pub iteration: usize,
    /// Time step of a transient solution.
    pub step: Option<usize>,
    /// Simulated time of a transient solution (s).
    pub time: Option<f64>,
    /// The derived fields written.
    pub derived: Vec<Derived>,
    /// Their datasets.
    pub datasets: Vec<&'static str>,
}

/// `foehn post derive` (D72): reads the solution of `restart`, a restart
/// file of `case`, and writes the fields with the derived fields `list` to
/// `file`, as the run would have written them (the same file as
/// `fields.vtkhdf` of a run with `[output] derived = list` at the same
/// iteration or step).
pub fn from_restart(
    case: &FlowCase,
    restart: &std::path::Path,
    list: &[Derived],
    file: &std::path::Path,
) -> Result<PostReport> {
    validate_list(case, list)?;
    with_restart(case, restart, |l| {
        let mut c = case.clone();
        c.output.vtkhdf = true;
        c.output.derived = sorted(list);
        let at = match l.clock {
            Some(k) => Some((k.step, Some(k.time))),
            None => Some((l.iteration, None)),
        };
        crate::case::write_fields(&c, file, l.mesh, l.geo, l.sol, l.analytic.as_ref(), at)?;
        Ok(PostReport {
            file: file.to_path_buf(),
            restart: restart.to_path_buf(),
            cells: l.mesh.n_cells,
            iteration: l.iteration,
            step: l.clock.map(|k| k.step),
            time: l.clock.map(|k| k.time),
            datasets: c
                .output
                .derived
                .iter()
                .flat_map(|d| d.datasets().iter().map(|(n, _)| *n))
                .collect(),
            derived: c.output.derived,
        })
    })
}

/// The derived fields `list` of the solution in `restart`, a restart file
/// of `case`, with the mesh and geometry they belong to (for checks and
/// scripts; [`from_restart`] writes them).
pub fn evaluate(
    case: &FlowCase,
    restart: &std::path::Path,
    list: &[Derived],
) -> Result<(Mesh, Geometry, Fields)> {
    validate_list(case, list)?;
    with_restart(case, restart, |l| {
        let f = compute(case, list, l.mesh, l.geo, l.sol)?;
        Ok((l.mesh.clone(), l.geo.clone(), f))
    })
}

/// The [`WallExchange`] of the solution in `restart`, a restart file of
/// `case` (D124 W4; for checks and for the wall models of later stages).
pub fn wall_exchange_from_restart(
    case: &FlowCase,
    restart: &std::path::Path,
) -> Result<(Mesh, Geometry, WallExchange)> {
    with_restart(case, restart, |l| {
        Ok((
            l.mesh.clone(),
            l.geo.clone(),
            wall_exchange(case, l.mesh, l.geo, l.sol),
        ))
    })
}

/// A solution read back from a restart file.
struct Loaded<'a> {
    mesh: &'a Mesh,
    geo: &'a Geometry,
    sol: &'a Solution<'a>,
    analytic: Option<crate::analytic::Analytic>,
    iteration: usize,
    clock: Option<crate::transient::Clock>,
}

/// Reads `restart`, builds the mesh and the discretisation of `case` and
/// hands the solution to `f`.
fn with_restart<R>(
    case: &FlowCase,
    restart: &std::path::Path,
    f: impl FnOnce(&Loaded) -> Result<R>,
) -> Result<R> {
    use anyhow::{bail, ensure};
    ensure!(
        crate::restart::AVAILABLE,
        "post derive reads restart files and writes VTKHDF: this build has no HDF5 (feature `hdf5`)"
    );
    let level = match case.verification.levels.as_slice() {
        [] => None,
        [l] => Some(*l),
        _ => bail!(
            "a restart file holds a single run; this case is a sweep over levels \
             (choose one with --set verification.levels=[N])"
        ),
    };
    let (mesh, geo, _) = crate::case::mesh_for(case, level)?;
    let nc = mesh.n_cells;
    let r = crate::restart::read(restart)?;
    ensure!(
        r.geometry.is_none(),
        "{}: restart files of moving meshes are not read by post derive yet",
        restart.display()
    );
    let (r, _) = r.to_local(&mesh.cell_id[..nc], &mesh.face_id, mesh.n_interior_faces())?;
    let clock = r.transient.as_ref().map(|t| t.1);
    match (clock.is_some(), case.time.is_some()) {
        (true, false) => bail!(
            "{} is the restart file of a transient run; the case is steady (no [time] table)",
            restart.display()
        ),
        (false, true) => bail!(
            "{} is the restart file of a steady run; the case is transient",
            restart.display()
        ),
        _ => {}
    }
    let analytic = match clock {
        Some(c) => case.analytic_at(c.time),
        None => case.analytic(),
    };
    let phys = crate::case::physics(case)?;
    let phys = crate::case::physics_on(case, &phys, &mesh, &geo)?;
    let analytic0 = match &case.time {
        Some(t) => case.analytic_at(t.start_time),
        None => analytic,
    };
    let values = case.compiled_boundaries()?;
    let t_now = clock.map_or(case.time.map_or(0.0, |t| t.start_time), |c| c.time);
    let arrays = crate::bc::arrays(
        &mesh,
        &geo,
        &case.boundary,
        phys.rho,
        crate::bc::At {
            analytic: analytic0.as_ref(),
            values: &values,
            t: t_now,
            turning: &case.turning(),
        },
    )?;
    let (disc, _) = crate::simplec::discretisation(&mesh, &geo, &arrays)?;
    // the driving of the periodic pairs (D29) is part of the body force
    let mut phys = phys;
    for (d, b) in arrays.drives.iter().zip(&r.state.periodic_force) {
        for (f, t) in phys.body_force.iter_mut().zip(d.direction) {
            *f += b * t;
        }
    }
    let sol = Solution {
        state: r.state,
        disc: &disc,
        phys,
    };
    f(&Loaded {
        mesh: &mesh,
        geo: &geo,
        sol: &sol,
        analytic,
        iteration: r.iteration,
        clock,
    })
}

/// What the wall treatment gives a module coupled to the walls (D124 W4,
/// the wall film of D123 and later CHT): per wall face the gas-side shear
/// stress and the transfer coefficients of the wall treatment.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct WallExchange {
    /// The wall faces (mesh face indices), boundary by boundary.
    pub faces: Vec<usize>,
    /// Shear stress of the fluid on the wall, tangential (Pa).
    pub shear: Vec<[f64; 3]>,
    /// The momentum transfer coefficient μ_f/Δn of the wall treatment
    /// (Pa s/m): τ = coefficient · (u_P − u_wall) tangential for a
    /// laminar wall; the turbulence models of stage 6 replace it by their
    /// wall functions.
    pub momentum: Vec<f64>,
    /// The heat transfer coefficient (W/(m² K)); with the energy equation
    /// (L1).
    pub heat: Option<Vec<f64>>,
    /// The mass transfer coefficient (kg/(m² s)); with species (L2).
    pub mass: Option<Vec<f64>>,
}

/// The [`WallExchange`] of a solution (D124 W4).
pub fn wall_exchange(case: &FlowCase, mesh: &Mesh, geo: &Geometry, fl: &Solution) -> WallExchange {
    let ni = mesh.n_interior_faces();
    let mut w = WallExchange::default();
    for p in &mesh.boundaries {
        if !matches!(
            case.boundary.get(&p.name),
            Some(BoundaryCondition::Wall { .. })
        ) {
            continue;
        }
        for f in p.range() {
            let (_, _, tau) = crate::case::wall_face_values(fl, geo, f, ni);
            let s = geo.face_area[f];
            let c = fl.disc.owner[f] as usize;
            let dn = dot(sub(geo.face_centroid[f], geo.cell_centroid[c]), s).abs() / norm(s);
            w.faces.push(f);
            w.shear.push(tau);
            w.momentum.push(fl.mu_eff_at(f) / dn);
        }
    }
    w
}

/// Wall shear stress and y⁺ per boundary face (zero off the walls).
fn wall_values(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    fl: &Solution,
) -> Vec<([f64; 3], f64)> {
    let ni = mesh.n_interior_faces();
    let mut v = vec![([0.0; 3], 0.0); mesh.n_boundary_faces()];
    for p in &mesh.boundaries {
        if !matches!(
            case.boundary.get(&p.name),
            Some(BoundaryCondition::Wall { .. })
        ) {
            continue;
        }
        for f in p.range() {
            let (_, _, tau) = crate::case::wall_face_values(fl, geo, f, ni);
            let s = geo.face_area[f];
            let c = fl.disc.owner[f] as usize;
            let y = dot(sub(geo.face_centroid[f], geo.cell_centroid[c]), s).abs() / norm(s);
            let u_tau = (norm(tau) / fl.phys.rho).sqrt();
            let nu = fl.phys.mu_at(f) / fl.phys.rho;
            // with a turbulence model: the y⁺ of its wall treatment
            let yp = match &fl.state.turbulence {
                Some(t) => t.y_plus[f - ni],
                None => y * u_tau / nu,
            };
            v[f - ni] = (tau, yp);
        }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * (1.0 + b.abs())
    }

    #[test]
    fn solid_body_rotation() {
        // u = −Ω y, v = Ω x
        let om = 3.0;
        let k = kinematics([[0.0, -om, 0.0], [om, 0.0, 0.0], [0.0; 3]]);
        assert_eq!(k.vorticity, [0.0, 0.0, 2.0 * om]);
        assert!(close(k.q, om * om));
        assert!(close(k.lambda2, -om * om));
        assert!(close(k.enstrophy, 2.0 * om * om));
        assert!(close(k.strain_rate, 0.0));
    }

    #[test]
    fn simple_shear_is_no_vortex() {
        // u = γ y: as much strain as rotation
        let g = 2.5;
        let k = kinematics([[0.0, g, 0.0], [0.0; 3], [0.0; 3]]);
        assert_eq!(k.vorticity, [0.0, 0.0, -g]);
        assert!(close(k.q, 0.0));
        assert!(close(k.lambda2, 0.0));
        assert!(close(k.strain_rate, g));
        assert!(close(k.enstrophy, 0.5 * g * g));
    }

    #[test]
    fn pure_strain() {
        // u = a x, v = −a y: Q = −a², λ₂ = 0 (eigenvalues a², a², 0)
        let a = 1.5;
        let k = kinematics([[a, 0.0, 0.0], [0.0, -a, 0.0], [0.0; 3]]);
        assert_eq!(k.vorticity, [0.0; 3]);
        assert!(close(k.q, -a * a));
        assert!(close(k.lambda2, a * a));
        assert!(close(k.strain_rate, 2.0 * a));
    }

    #[test]
    fn eigenvalues_of_symmetric_matrices() {
        // A = R diag(e) Rᵀ with a rotation R
        let (c, s) = (0.3f64.cos(), 0.3f64.sin());
        let (c2, s2) = (1.1f64.cos(), 1.1f64.sin());
        let rz = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        let rx = [[1.0, 0.0, 0.0], [0.0, c2, -s2], [0.0, s2, c2]];
        let mul = |a: [[f64; 3]; 3], b: [[f64; 3]; 3]| {
            let mut m = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        m[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
            m
        };
        let r = mul(rz, rx);
        let mut rt = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                rt[i][j] = r[j][i];
            }
        }
        for e in [
            [-2.0, 0.5, 4.0],
            [1.0, 1.0, 3.0],
            [-1e-3, -1e-3, -1e-3],
            [0.0, 0.0, 7.0],
        ] {
            let d = [[e[0], 0.0, 0.0], [0.0, e[1], 0.0], [0.0, 0.0, e[2]]];
            let a = mul(mul(r, d), rt);
            let got = symmetric_eigenvalues(a);
            for i in 0..3 {
                assert!((got[i] - e[i]).abs() < 1e-14 * 8.0, "{got:?} {e:?}");
            }
        }
    }

    /// The wall distance in a channel between walls at y = 0 and y = 1
    /// (the ends are an inlet and an outlet, which do not count): min(y,
    /// 1 − y) at every centroid, exactly on hex cells; the polyhedral mesh
    /// gives the distance to its (plane) wall faces, the same function.
    #[test]
    fn wall_distance_of_a_channel() {
        for mesh in ["hex", "poly", "tet"] {
            let case = FlowCase::from_toml(
                &format!(
                    r#"
[mesh.generate]
type = "{mesh}"
dimension = 2
upper = [2.0, 1.0, 0.1]
cells = [8, 6, 1]
[fluid]
viscosity = 0.1
[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.xmax]
type = "static-pressure-outlet"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[output]
derived = ["wall-distance"]
"#
                ),
                std::path::Path::new("."),
            )
            .unwrap();
            validate(&case).unwrap();
            let (m, g, _) = crate::case::mesh_for(&case, None).unwrap();
            let d = wall_distance(&case, &m, &g).unwrap();
            assert_eq!(d.len(), m.n_cells);
            for (c, &dc) in d.iter().enumerate() {
                let y = g.cell_centroid[c][1];
                let e = y.min(1.0 - y);
                assert!((dc - e).abs() < 1e-12, "{mesh}: cell {c}, {dc} vs {e}");
            }
        }
    }

    #[test]
    fn names_round_trip() {
        for d in Derived::ALL {
            assert_eq!(Derived::parse(d.name()), Some(d));
            let t: toml::Value = toml::Value::try_from(d).unwrap();
            assert_eq!(t.as_str(), Some(d.name()));
        }
    }
}

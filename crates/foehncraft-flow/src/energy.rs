//! The energy equation of the fluid (stage 7, PR 7a1; D150, D151, D164):
//! the temperature T in the form weighted with ρ c_p,
//!
//! ρ c_p ∂T/∂t + ∇·(Ċ T) − T ∇·Ċ − ∇·(k ∇T) = q''' + Φ,   Ċ_f = ṁ_f c_p,f,
//!
//! solved on the generic scalar transport of D116 ([`crate::scalar`])
//! after the pressure correction of every outer iteration. The convection
//! term is that of the bounded form with the heat capacity flux Ċ instead
//! of the mass flux, so that it equals ρ c_p u·∇T also where c_p varies;
//! the diffusion is the conduction with k at the faces (the series
//! conductance between materials in 7d1 is then a face value as well).
//! c_p(T), k(T) and μ(T) of the material are refreshed per cell and face
//! from the current temperature before each energy solve (μ: the face
//! viscosity of the next momentum solve); ρ stays the one number of the
//! case in 7a1 (a density that follows T is 7c).
//!
//! Case file:
//!
//! ```toml
//! [fluid]
//! material = "air"              # needs specific_heat and conductivity
//!
//! [energy]
//! initial = 300.0               # K: number or expression of x, y, z
//! source = 0.0                  # q''' (W/m³): number or expression of x, y, z (t)
//! viscous_heating = true        # Φ = τ:∇u (default false)
//!
//! [boundary.inlet.energy]
//! temperature = 300.0           # K
//! [boundary.wall.energy]
//! heat_flux = 500.0             # W/m² into the fluid
//! # heat_flow = 2.0             # W into the fluid through the whole boundary
//! # heat_transfer_coefficient = 10.0   # W/(m² K), with
//! # ambient_temperature = 293.15       # K: q = h (T_amb − T_wall)
//! # temperature = 350.0
//!
//! [solver.energy]               # as [solver.scalar]
//! scheme = "linear-upwind"
//! ```
//!
//! Boundaries without an `energy` table: walls, slip walls, symmetry
//! planes and wedges are adiabatic, outlets have a zero normal gradient
//! (the backflow carries the owner value), periodic pairs carry T across;
//! inlets need a `temperature` (ENE002).

use crate::bc::BoundaryCondition;
use crate::diag::{self, Diagnostic};
use crate::expr::{self, Env, Scalar, Vars};
use crate::materials::{Conditions, Conductivity, Curve, MaterialSection, Property, curve, lookup};
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// The name of the temperature field (output, restart, log columns).
pub const FIELD: &str = "temperature";

/// `[energy]`: the energy equation of the fluid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EnergySection {
    /// Start temperature (K): a number or an expression of x, y, z (D53).
    pub initial: Scalar,
    /// Volumetric heat source q''' (W/m³): a number or an expression of
    /// x, y, z and, in a transient case, t (D53).
    #[serde(default)]
    pub source: Scalar,
    /// The viscous heating Φ = τ:∇u = μ_eff (∇u + ∇uᵀ):∇u as a source
    /// (default off; it matters where the Brinkman number μU²/(kΔT) is not
    /// small: oils, lubrication gaps, high-speed gas flows).
    #[serde(default)]
    pub viscous_heating: bool,
    /// Refresh c_p, k and μ from the temperature every outer iteration
    /// (default true). Off: the properties stay those at `[fluid]
    /// temperature`.
    #[serde(default = "yes")]
    pub temperature_dependent: bool,
    /// The exact temperature, for the verification (an expression of x,
    /// y, z and, in a transient case, t), as `exact` of a scalar.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exact: Option<Scalar>,
    /// Turbulent Prandtl number Pr_t of the `constant` model (D135); the
    /// turbulent conductivity is c_p μ_t/Pr_t.
    #[serde(default = "default_prt")]
    pub turbulent_prandtl: f64,
    /// Model of Pr_t: `constant` or `kays-crawford`.
    #[serde(default)]
    pub turbulent_prandtl_model: crate::scalar::PrandtlModel,
    /// Thermal wall function at the walls of a turbulent flow (D135).
    #[serde(default)]
    pub wall_function: crate::scalar::ThermalWallFunction,
}

fn yes() -> bool {
    true
}

fn default_prt() -> f64 {
    0.85
}

/// `[boundary.<name>.energy]`: the thermal condition of a boundary. One
/// of `temperature`, `heat_flux`, `heat_flow` and the pair
/// `heat_transfer_coefficient` + `ambient_temperature`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EnergyBoundary {
    /// Fixed temperature (K): a number or an expression of x, y, z and, in
    /// a transient case, t. Inlets and walls.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temperature: Option<Scalar>,
    /// Heat flux into the fluid (W/m²): a number or an expression. Walls,
    /// slip walls, symmetry planes, wedges.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heat_flux: Option<Scalar>,
    /// Heat flow into the fluid through the whole boundary (W), spread as
    /// a uniform flux over its area.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heat_flow: Option<f64>,
    /// Heat transfer coefficient h (W/(m² K)) to an ambient at
    /// `ambient_temperature`: the flux into the fluid is
    /// h (T_amb − T_wall).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heat_transfer_coefficient: Option<Scalar>,
    /// The ambient temperature of `heat_transfer_coefficient` (K).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ambient_temperature: Option<Scalar>,
}

/// Moves the per-physics sub-tables `energy` out of the
/// `[boundary.<name>]` tables of a parsed case and returns them.
pub fn extract_boundary_tables(v: &mut toml::Table) -> Result<BTreeMap<String, EnergyBoundary>> {
    let mut out = BTreeMap::new();
    let Some(toml::Value::Table(bs)) = v.get_mut("boundary") else {
        return Ok(out);
    };
    for (name, b) in bs.iter_mut() {
        let toml::Value::Table(t) = b else { continue };
        let Some(s) = t.remove("energy") else {
            continue;
        };
        let e: EnergyBoundary = s.try_into().map_err(|e: toml::de::Error| {
            diag::err(
                Diagnostic::new(
                    diag::ENE_BOUNDARY,
                    format!("[boundary.{name}.energy]: {}", e.message()),
                )
                .key(format!("boundary.{name}.energy")),
            )
        })?;
        out.insert(name.clone(), e);
    }
    Ok(out)
}

/// What a boundary of the flow is for the energy equation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Side {
    WallLike,
    Inlet,
    Outlet,
    Periodic,
}

pub(crate) fn side(b: &BoundaryCondition) -> Side {
    use BoundaryCondition as B;
    match b {
        B::Wall { .. } | B::SlipWall | B::Symmetry | B::Wedge => Side::WallLike,
        B::VelocityInlet { .. }
        | B::MassFlowInlet { .. }
        | B::StaticPressureInlet { .. }
        | B::TotalPressureInlet { .. } => Side::Inlet,
        B::StaticPressureOutlet { .. }
        | B::TotalPressureOutlet { .. }
        | B::ConvectiveOutlet { .. } => Side::Outlet,
        B::Periodic { .. } => Side::Periodic,
    }
}

/// The condition of one boundary, compiled.
#[derive(Debug, Clone)]
pub(crate) enum Condition {
    /// Zero normal gradient (adiabatic walls, outlets).
    ZeroGradient,
    /// Fixed temperature.
    Temperature(expr::Value),
    /// Heat flux into the fluid (W/m²).
    HeatFlux(expr::Value),
    /// h and T_amb.
    Convection(expr::Value, expr::Value),
}

impl EnergyBoundary {
    fn count(&self) -> usize {
        usize::from(self.temperature.is_some())
            + usize::from(self.heat_flux.is_some())
            + usize::from(self.heat_flow.is_some())
            + usize::from(self.heat_transfer_coefficient.is_some())
    }
}

/// The compiled condition of every boundary of the case (by name), with
/// `area` the area of each boundary (for `heat_flow`).
pub(crate) fn conditions(
    case: &crate::case::FlowCase,
    env: &Env,
    area: &BTreeMap<String, f64>,
) -> Result<BTreeMap<String, Condition>> {
    let vars = Vars::SPACE.with(if case.time.is_some() {
        Vars::T
    } else {
        Vars::NONE
    });
    let mut out = BTreeMap::new();
    for (name, bc) in &case.boundary {
        let key = format!("boundary.{name}.energy");
        let c = match case.boundary_energy.get(name) {
            None => Condition::ZeroGradient,
            Some(t) => {
                if let Some(v) = &t.temperature {
                    Condition::Temperature(expr::value(
                        v,
                        env,
                        vars,
                        &format!("{key}.temperature"),
                    )?)
                } else if let Some(v) = &t.heat_flux {
                    Condition::HeatFlux(expr::value(v, env, vars, &format!("{key}.heat_flux"))?)
                } else if let Some(q) = t.heat_flow {
                    let a = area.get(name).copied().unwrap_or(0.0);
                    crate::need!(
                        a > 0.0,
                        ENE_BOUNDARY,
                        "[{key}] heat_flow: the boundary has no area"
                    );
                    Condition::HeatFlux(expr::Value::Number(q / a))
                } else if let (Some(h), Some(ta)) =
                    (&t.heat_transfer_coefficient, &t.ambient_temperature)
                {
                    Condition::Convection(
                        expr::value(h, env, vars, &format!("{key}.heat_transfer_coefficient"))?,
                        expr::value(ta, env, vars, &format!("{key}.ambient_temperature"))?,
                    )
                } else {
                    Condition::ZeroGradient
                }
            }
        };
        let _ = bc;
        out.insert(name.clone(), c);
    }
    Ok(out)
}

/// Checks `[energy]` and the `energy` sub-tables of the boundaries.
pub fn validate(case: &crate::case::FlowCase, env: &Env) -> Result<()> {
    let Some(e) = &case.energy else {
        crate::need!(
            case.boundary_energy.is_empty(),
            ENE_BOUNDARY,
            "[boundary.{}.energy] without [energy]: the energy equation is solved only with an \
             [energy] table",
            case.boundary_energy
                .keys()
                .next()
                .map_or("", String::as_str)
        );
        return Ok(());
    };
    let transient = case.time.is_some();
    let time = if transient { Vars::T } else { Vars::NONE };
    expr::value(&e.initial, env, Vars::SPACE, "energy.initial")?;
    expr::value(&e.source, env, Vars::SPACE.with(time), "energy.source")?;
    if let Some(x) = &e.exact {
        expr::value(x, env, Vars::SPACE.with(time), "energy.exact")?;
    }
    crate::need!(
        e.turbulent_prandtl > 0.0 && e.turbulent_prandtl.is_finite(),
        CASE_RANGE,
        "[energy] turbulent_prandtl must be > 0, is {}",
        e.turbulent_prandtl
    );
    crate::need!(
        !case.scalars.contains_key(FIELD),
        ENE_UNSUPPORTED,
        "[scalars.{FIELD}] and [energy]: the energy equation solves the field `{FIELD}`; give \
         the passive scalar another name"
    );
    crate::need!(
        case.regions.is_empty(),
        ENE_UNSUPPORTED,
        "[energy] with [regions]: the energy equation takes the one `[fluid]` material in this \
         version (several fluid materials and solid regions come with 7d1)"
    );
    case.solver.energy.validate(transient)?;
    for (b, t) in &case.boundary_energy {
        let key = format!("boundary.{b}.energy");
        let Some(bc) = case.boundary.get(b) else {
            return Err(diag::err(
                Diagnostic::new(
                    diag::ENE_BOUNDARY,
                    format!("[{key}] belongs to no boundary condition `[boundary.{b}]`"),
                )
                .key(key),
            ));
        };
        crate::need!(
            t.count() == 1,
            ENE_BOUNDARY,
            "[{key}] needs exactly one of `temperature`, `heat_flux`, `heat_flow` and \
             `heat_transfer_coefficient` (with `ambient_temperature`)"
        );
        crate::need!(
            t.heat_transfer_coefficient.is_some() == t.ambient_temperature.is_some(),
            ENE_BOUNDARY,
            "[{key}]: `heat_transfer_coefficient` and `ambient_temperature` go together"
        );
        let sd = side(bc);
        let ok = match sd {
            Side::WallLike => true,
            Side::Inlet => t.temperature.is_some(),
            Side::Outlet | Side::Periodic => false,
        };
        crate::need!(
            ok,
            ENE_BOUNDARY,
            "[{key}]: a {} takes {}",
            bc.type_name(),
            match sd {
                Side::Inlet => "a `temperature`, not a flux",
                Side::Outlet => "no energy condition (outlets keep a zero normal gradient)",
                _ =>
                    "no energy condition (T crosses the periodic pair; the periodic heat-flux \
                      mode of the energy equation is not in this version)",
            }
        );
        let vars = Vars::SPACE.with(time);
        for (k, v) in [
            ("temperature", &t.temperature),
            ("heat_flux", &t.heat_flux),
            ("heat_transfer_coefficient", &t.heat_transfer_coefficient),
            ("ambient_temperature", &t.ambient_temperature),
        ] {
            if let Some(v) = v {
                let x = expr::value(v, env, vars, &format!("{key}.{k}"))?;
                if k == "heat_transfer_coefficient"
                    && let Some(h) = x.constant()
                {
                    crate::need!(
                        h > 0.0 && h.is_finite(),
                        CASE_RANGE,
                        "[{key}] heat_transfer_coefficient must be > 0, is {h}"
                    );
                }
            }
        }
        if let Some(q) = t.heat_flow {
            crate::need!(
                q.is_finite(),
                CASE_RANGE,
                "[{key}] heat_flow must be a finite number"
            );
        }
    }
    for (b, bc) in &case.boundary {
        crate::need!(
            side(bc) != Side::Inlet
                || case
                    .boundary_energy
                    .get(b)
                    .is_some_and(|t| t.temperature.is_some()),
            ENE_INLET,
            "[boundary.{b}] is an inlet: with [energy] it needs `[boundary.{b}.energy] \
             temperature`"
        );
    }
    fluid(case, env)?;
    Ok(())
}

/// The thermal properties of the fluid as functions of T.
#[derive(Debug, Clone)]
pub struct Fluid {
    /// The material (`None`: the `[fluid]` shorthand).
    pub material: Option<String>,
    /// ρ (kg/m³), constant in 7a1.
    pub rho: f64,
    /// c_p(T) (J/(kg K)).
    pub cp: Curve,
    /// k(T) (W/(m K)).
    pub k: Curve,
    /// μ(T) (Pa s), when the viscosity is a form of T (else the flow keeps
    /// its viscosity).
    pub mu: Option<Curve>,
}

/// The thermal properties of the `[fluid]` (its material, with the
/// `[fluid]` numbers replacing the material's): MAT007 without c_p or k.
pub fn fluid(case: &crate::case::FlowCase, env: &Env) -> Result<Fluid> {
    let cond: Conditions = case.fluid.conditions();
    let (m, owner) = match &case.fluid.material {
        Some(n) => (lookup(n, &case.materials)?, format!("materials.{n}")),
        None => (MaterialSection::default(), "fluid".to_string()),
    };
    let pick = |given: &Option<Scalar>, own: &Option<Property>| -> Option<(Property, bool)> {
        match (given, own) {
            (Some(s), _) => Some((Property::from(s.clone()), true)),
            (None, Some(p)) => Some((p.clone(), false)),
            (None, None) => None,
        }
    };
    let own_k = match &m.conductivity {
        Some(Conductivity::Isotropic(p)) => Some(p.clone()),
        Some(Conductivity::Orthotropic(_)) => {
            return Err(diag::err(
                Diagnostic::new(
                    diag::ENE_UNSUPPORTED,
                    format!(
                        "{owner}: an orthotropic conductivity belongs to a solid; a fluid takes \
                         one value"
                    ),
                )
                .key(format!("{owner}.conductivity")),
            ));
        }
        None => None,
    };
    let cp = pick(&case.fluid.specific_heat, &m.specific_heat);
    let k = pick(&case.fluid.conductivity, &own_k);
    let missing: Vec<&str> = [
        ("specific_heat", cp.is_none()),
        ("conductivity", k.is_none()),
    ]
    .iter()
    .filter(|x| x.1)
    .map(|x| x.0)
    .collect();
    if !missing.is_empty() {
        let what = match &case.fluid.material {
            Some(n) => format!(
                "material `{n}` has no {}; the energy equation needs c_p and k of the fluid \
                 (give them in [materials.{n}] with `base = \"{n}\"`, or in [fluid])",
                missing.join(" and no ")
            ),
            None => format!(
                "[fluid] has no {}; the energy equation needs c_p and k of the fluid (give \
                 `specific_heat` and `conductivity` in [fluid], or a `material`)",
                missing.join(" and no ")
            ),
        };
        return Err(diag::err(
            Diagnostic::new(diag::MAT_THERMAL, what).key("fluid"),
        ));
    }
    let (cp, cp_given) = cp.unwrap_or((Property::Number(f64::NAN), false));
    let (k, k_given) = k.unwrap_or((Property::Number(f64::NAN), false));
    let label = |given: bool| {
        if given {
            "fluid".to_string()
        } else {
            owner.clone()
        }
    };
    let cp = curve(&cp, "specific_heat", &label(cp_given), &m, env, cond)?;
    let k = curve(&k, "conductivity", &label(k_given), &m, env, cond)?;
    // μ(T): only a form of T of the material's own viscosity (a [fluid]
    // viscosity, a number or an expression of x, y, z stays as it is)
    let mu = match (&case.fluid.viscosity, &m.viscosity) {
        (None, Some(p @ Property::Form(f)))
            if !matches!(f, crate::materials::Form::Expression { .. }) =>
        {
            let c = curve(p, "viscosity", &owner, &m, env, cond)?;
            (!c.is_constant()).then_some(c)
        }
        _ => None,
    };
    Ok(Fluid {
        material: case.fluid.material.clone(),
        rho: case.rho(),
        cp,
        k,
        mu,
    })
}

/// The temperature of zero sensible enthalpy (K).
pub const H_REF_T: f64 = 298.15;

/// The sensible enthalpy h_s(T) = ∫ c_p dT from [`H_REF_T`] of a c_p
/// curve: exact for a constant c_p, else the trapezoidal integral on a
/// grid of [`Enthalpy::STEP`] K, linear between the grid points (and with
/// the end values of c_p outside the grid). The energy equation is
/// conservative in this h_s (D150): the convection and storage terms of the
/// matrix in T are corrected with the difference of the h-form at the
/// current temperature, so the residual at convergence is that of h_s.
#[derive(Debug, Clone, PartialEq)]
pub struct Enthalpy {
    cp: Option<f64>,
    t0: f64,
    h: Vec<f64>,
    cp_lo: f64,
    cp_hi: f64,
}

impl Enthalpy {
    /// The grid step (K).
    pub const STEP: f64 = 0.25;

    /// h_s of a constant c_p.
    pub fn constant(cp: f64) -> Self {
        Self {
            cp: Some(cp),
            t0: 0.0,
            h: Vec::new(),
            cp_lo: cp,
            cp_hi: cp,
        }
    }

    /// h_s of the curve `cp`.
    pub fn of(cp: &Curve) -> Self {
        if cp.is_constant() {
            return Self::constant(cp.at(H_REF_T));
        }
        let [lo, hi] = Curve::SAMPLE_RANGE;
        let n = ((hi - lo) / Self::STEP).round() as usize + 1;
        let mut h = Vec::with_capacity(n);
        let mut acc = 0.0;
        let mut prev = cp.at(lo);
        h.push(0.0);
        for i in 1..n {
            let c = cp.at(lo + i as f64 * Self::STEP);
            acc += 0.5 * (prev + c) * Self::STEP;
            h.push(acc);
            prev = c;
        }
        let mut e = Self {
            cp: None,
            t0: lo,
            h,
            cp_lo: cp.at(lo),
            cp_hi: prev,
        };
        let shift = e.at(H_REF_T);
        for v in &mut e.h {
            *v -= shift;
        }
        e
    }

    /// True for a constant c_p (no correction needed).
    pub fn is_constant(&self) -> bool {
        self.cp.is_some()
    }

    /// h_s (J/kg) at `t` (K).
    pub fn at(&self, t: f64) -> f64 {
        if let Some(cp) = self.cp {
            return cp * (t - H_REF_T);
        }
        let n = self.h.len();
        let x = (t - self.t0) / Self::STEP;
        if x <= 0.0 {
            return self.h[0] + self.cp_lo * (t - self.t0);
        }
        let last = (n - 1) as f64;
        if x >= last {
            return self.h[n - 1] + self.cp_hi * (t - (self.t0 + last * Self::STEP));
        }
        let i = (x.floor() as usize).min(n - 2);
        let w = x - i as f64;
        self.h[i] + w * (self.h[i + 1] - self.h[i])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::Env;
    use crate::materials::{CaseMaterials, Conditions, Property, curve};

    #[test]
    fn enthalpy_integrates_the_specific_heat() {
        let env = Env::default();
        let m = MaterialSection::default();
        let cond = Conditions::default();
        // c_p = 1000 + 0.5 T: h = 1000 (T − T_r) + 0.25 (T² − T_r²)
        let cp = curve(
            &Property::Expr("1000 + 0.5 * T".into()),
            "specific_heat",
            "test",
            &m,
            &env,
            cond,
        )
        .unwrap();
        let h = Enthalpy::of(&cp);
        assert!(!h.is_constant());
        for t in [250.0, 298.15, 300.0, 412.3, 1500.0] {
            let exact = 1000.0 * (t - H_REF_T) + 0.25 * (t * t - H_REF_T * H_REF_T);
            // linear interpolation of the 0.25 K table: |err| ≤ c_p' Δ²/8 = 0.004 J/kg
            assert!((h.at(t) - exact).abs() < 5e-3, "{t}: {} {exact}", h.at(t));
        }
        let c = Enthalpy::constant(4184.0);
        assert_eq!(c.at(H_REF_T + 1.0), 4184.0);
        let _ = CaseMaterials::default();
    }
}

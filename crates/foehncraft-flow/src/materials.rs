//! Material properties (D94–D102), stage-5 part: density and viscosity.
//!
//! Three layers (D94):
//!
//! 1. the **database**: named materials with their data, validity ranges
//!    and sources, built into the binary from `data/materials.toml`
//!    (air, water, an engine oil, ethylene glycol, glycerol), read-only;
//! 2. the **case materials** `[materials.<name>]`: new materials, or
//!    changes of a database material with `base = "<name>"`, and the
//!    assignment of materials to the fluid (`[fluid] material`) and to
//!    regions of cell zones (`[regions.<name>]`, D112);
//! 3. the **resolved property program** ([`PropertySet`]): every property
//!    of every region as one value or one compiled expression (D53), with
//!    its source, written into the resolved case and hashed into the run
//!    record (D85). The solver sees only this; the database is not read
//!    again during the run, and properties do not change during a run
//!    (D100).
//!
//! Forms (D95) evaluated in this version: a number, an expression (of x,
//! y, z for the viscosity; `T` is the case temperature), `constant`,
//! `polynomial` in T, `sutherland`, `power-law`, `table` in T (linear or
//! the monotone cubic of Steffen 1990, clamped with a warning outside the
//! table, D96) and `ideal-gas` for the density. The flow is isothermal
//! and incompressible until the energy equation (L1): every form is
//! evaluated once at `[fluid] temperature` (default 293.15 K) and the
//! operating pressure `[fluid] operating_pressure` (default 101 325 Pa),
//! the density must be one number for the whole fluid (a varying density
//! comes with L1), and only the viscosity may vary in space. The thermal
//! data (c_p, k, NASA-7/9, IAPWS-95) follow with L1.
//!
//! `[fluid] density` and `viscosity` stay valid as the shorthand of one
//! constant material on the whole mesh (D100), so every older case keeps
//! running with the same numbers.

use crate::diag::{self, Diagnostic};
use crate::expr::{self, Env, Scalar, Value, Vars};
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::PathBuf;
use std::sync::{Arc, OnceLock};

/// Universal gas constant (J/(mol K), CODATA 2018, exact).
pub const R_GAS: f64 = 8.314_462_618;
/// The default case temperature (K).
pub const DEFAULT_TEMPERATURE: f64 = 293.15;
/// The default operating pressure (Pa).
pub const DEFAULT_PRESSURE: f64 = 101_325.0;

/// Interpolation of a property table (D96).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Interpolation {
    /// Piecewise linear.
    #[default]
    Linear,
    /// Monotone piecewise cubic (Steffen 1990): no overshoot between the
    /// points, continuous first derivative.
    MonotoneCubic,
}

/// A property given by a form (D95).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(tag = "form", rename_all = "kebab-case", deny_unknown_fields)]
pub enum Form {
    /// A constant.
    Constant {
        /// The value (SI).
        value: f64,
        /// Where the value comes from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// Σ aᵢ Tⁱ (Horner form), T in K.
    Polynomial {
        /// a₀, a₁, … (SI).
        coefficients: Vec<f64>,
        /// Validity range in T (K).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        range: Option<[f64; 2]>,
        /// Where the coefficients come from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// μ = μ_ref (T/T_ref)^1.5 (T_ref + S)/(T + S) (Sutherland 1893), for
    /// the viscosity of a gas.
    Sutherland {
        /// μ_ref (Pa s).
        mu_ref: f64,
        /// T_ref (K).
        t_ref: f64,
        /// S (K).
        s: f64,
        /// Validity range in T (K).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        range: Option<[f64; 2]>,
        /// Where the coefficients come from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// value_ref (T/T_ref)^exponent.
    PowerLaw {
        /// The value at T_ref (SI).
        value_ref: f64,
        /// T_ref (K).
        t_ref: f64,
        /// The exponent.
        exponent: f64,
        /// Validity range in T (K).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        range: Option<[f64; 2]>,
        /// Where the coefficients come from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// A table in T (K), clamped with a warning outside (D96).
    Table {
        /// Temperatures (K), strictly increasing.
        #[serde(default)]
        temperature: Vec<f64>,
        /// Values (SI).
        #[serde(default)]
        values: Vec<f64>,
        /// A file of two columns (T, value) instead of the arrays.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        file: Option<PathBuf>,
        /// Interpolation.
        #[serde(default)]
        interpolation: Interpolation,
        /// Where the data come from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// An expression (D53); `T` is the case temperature.
    Expression {
        /// The expression.
        expression: String,
        /// Where it comes from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// ρ = p M / (R T) at the operating pressure (for the density).
    IdealGas {
        /// Molar mass (kg/mol); default the material's `molar_mass`.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        molar_mass: Option<f64>,
        /// Where it comes from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// NASA 7-coefficient polynomials (McBride, Zehe & Gordon 2002); with
    /// the energy equation (L1).
    Nasa7 {
        /// The coefficients.
        #[serde(default)]
        coefficients: Vec<Vec<f64>>,
    },
    /// NASA 9-coefficient polynomials; with the energy equation (L1).
    Nasa9 {
        /// The coefficients.
        #[serde(default)]
        coefficients: Vec<Vec<f64>>,
    },
}

impl Form {
    /// The form's name.
    pub fn name(&self) -> &'static str {
        match self {
            Form::Constant { .. } => "constant",
            Form::Polynomial { .. } => "polynomial",
            Form::Sutherland { .. } => "sutherland",
            Form::PowerLaw { .. } => "power-law",
            Form::Table { .. } => "table",
            Form::Expression { .. } => "expression",
            Form::IdealGas { .. } => "ideal-gas",
            Form::Nasa7 { .. } => "nasa-7",
            Form::Nasa9 { .. } => "nasa-9",
        }
    }

    fn source(&self) -> Option<&str> {
        match self {
            Form::Constant { source, .. }
            | Form::Polynomial { source, .. }
            | Form::Sutherland { source, .. }
            | Form::PowerLaw { source, .. }
            | Form::Table { source, .. }
            | Form::Expression { source, .. }
            | Form::IdealGas { source, .. } => source.as_deref(),
            _ => None,
        }
    }

    fn range(&self) -> Option<[f64; 2]> {
        match self {
            Form::Polynomial { range, .. }
            | Form::Sutherland { range, .. }
            | Form::PowerLaw { range, .. } => *range,
            Form::Table { temperature, .. } if temperature.len() >= 2 => {
                Some([temperature[0], temperature[temperature.len() - 1]])
            }
            _ => None,
        }
    }
}

/// A property: a number, an expression string, or a form.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum Property {
    /// A number (SI).
    Number(f64),
    /// An expression (D53).
    Expr(String),
    /// A form of D95.
    Form(Form),
}

impl From<Scalar> for Property {
    fn from(s: Scalar) -> Self {
        match s {
            Scalar::Number(v) => Property::Number(v),
            Scalar::Expr(e) => Property::Expr(e),
        }
    }
}

/// The kind of a material.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum MaterialKind {
    /// A gas.
    Gas,
    /// A liquid.
    Liquid,
    /// A solid (with the energy equation, L1).
    Solid,
}

/// A material, `[materials.<name>]` of the case or of the database.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct MaterialSection {
    /// A material (of the case or the database) this one starts from;
    /// the keys given here replace its keys.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub base: Option<String>,
    /// One line on what it is.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    /// Gas, liquid or solid.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub kind: Option<MaterialKind>,
    /// Temperature range of the data (K).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperature_range: Option<[f64; 2]>,
    /// Molar mass (kg/mol).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub molar_mass: Option<f64>,
    /// Where the molar mass comes from.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub molar_mass_source: Option<String>,
    /// Density ρ (kg/m³).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub density: Option<Property>,
    /// Dynamic viscosity μ (Pa s).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub viscosity: Option<Property>,
}

/// The kind of a region (D112).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum RegionKind {
    /// A fluid region: the flow is solved on it.
    #[default]
    Fluid,
    /// A solid (conjugate heat transfer, with L1).
    Solid,
    /// A porous zone (with L1 and D55).
    Porous,
}

/// `[regions.<name>]`: cell zones grouped into a region with a kind and a
/// material (D112).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct RegionSection {
    /// Fluid, solid or porous; this version solves fluids only.
    pub kind: RegionKind,
    /// The cell zones of the region (mesh zone names).
    pub zones: Vec<String>,
    /// The material; default the material of `[fluid]`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub material: Option<String>,
}

/// The database file.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct DatabaseFile {
    #[allow(dead_code)]
    format: String,
    #[allow(dead_code)]
    format_version: String,
    materials: BTreeMap<String, MaterialSection>,
}

const DATABASE_TOML: &str = include_str!("../data/materials.toml");

/// The built-in material database.
pub fn database() -> &'static BTreeMap<String, MaterialSection> {
    static DB: OnceLock<BTreeMap<String, MaterialSection>> = OnceLock::new();
    DB.get_or_init(|| {
        let f: DatabaseFile =
            toml::from_str(DATABASE_TOML).expect("the built-in material database parses");
        f.materials
    })
}

/// SHA-256 of the built-in database (the run record names it, D85).
pub fn database_hash() -> String {
    crate::engine::sha256(DATABASE_TOML.as_bytes())
}

/// The inputs of the property resolution that the case gives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Conditions {
    /// The isothermal temperature (K).
    pub temperature: f64,
    /// The operating pressure (Pa).
    pub pressure: f64,
}

impl Default for Conditions {
    fn default() -> Self {
        Self {
            temperature: DEFAULT_TEMPERATURE,
            pressure: DEFAULT_PRESSURE,
        }
    }
}

fn mat_err(code: &'static str, key: &str, m: String) -> anyhow::Error {
    diag::err(Diagnostic::new(code, m).key(key))
}

/// The material `name` with its `base` chain merged (case materials first,
/// then the database).
pub fn lookup(name: &str, case: &BTreeMap<String, MaterialSection>) -> Result<MaterialSection> {
    let mut seen: Vec<(String, bool)> = Vec::new();
    let mut layers: Vec<MaterialSection> = Vec::new();
    let mut cur = name.to_string();
    // a database material's base is a database material; a case material
    // may take the name of a database one and change it, its `base` of
    // the same name then means the database entry
    let mut db_only = false;
    loop {
        if seen.contains(&(cur.clone(), db_only)) || seen.len() > 32 {
            let names: Vec<&str> = seen.iter().map(|s| s.0.as_str()).collect();
            return Err(mat_err(
                diag::MAT_BASE,
                &format!("materials.{name}.base"),
                format!(
                    "the materials {} → {cur} refer to each other through `base` in a cycle",
                    names.join(" → ")
                ),
            ));
        }
        seen.push((cur.clone(), db_only));
        let (m, in_case) = match (db_only, case.get(&cur), database().get(&cur)) {
            (false, Some(m), _) => (m, true),
            (_, _, Some(m)) => (m, false),
            _ => {
                let names: Vec<&str> = case
                    .keys()
                    .chain(database().keys())
                    .map(String::as_str)
                    .collect();
                let hint = match crate::bc::did_you_mean(&cur, names.iter().copied()) {
                    Some(c) => format!(" (did you mean `{c}`?)"),
                    None => String::new(),
                };
                return Err(mat_err(
                    diag::MAT_UNKNOWN,
                    "fluid.material",
                    format!("unknown material `{cur}`{hint}; `foehn materials` lists the database"),
                ));
            }
        };
        layers.push(m.clone());
        match &m.base {
            Some(b) => {
                db_only = !in_case || *b == cur;
                cur = b.clone();
            }
            None => break,
        }
    }
    // merge from the root up
    let mut out = MaterialSection::default();
    for m in layers.into_iter().rev() {
        let m2 = m.clone();
        out.description = m2.description.or(out.description);
        out.kind = m2.kind.or(out.kind);
        out.temperature_range = m2.temperature_range.or(out.temperature_range);
        out.molar_mass = m2.molar_mass.or(out.molar_mass);
        out.molar_mass_source = m2.molar_mass_source.or(out.molar_mass_source);
        out.density = m2.density.or(out.density);
        out.viscosity = m2.viscosity.or(out.viscosity);
    }
    out.base = None;
    Ok(out)
}

/// One resolved property.
#[derive(Debug, Clone, PartialEq)]
pub struct Resolved {
    /// The value (a number, or a program in x, y, z).
    pub value: Value,
    /// How it was obtained (form, T, source), for the resolved case.
    pub origin: String,
}

/// The properties of one material at the case conditions.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialProperties {
    /// The material (`None`: the `[fluid]` shorthand without a material).
    pub material: Option<String>,
    /// Density.
    pub density: Resolved,
    /// Viscosity.
    pub viscosity: Resolved,
}

/// Monotone cubic interpolation of Steffen (1990), *A simple method for
/// monotonic interpolation in one dimension*, A&A 239, 443–450.
pub fn steffen(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[n - 1] {
        return ys[n - 1];
    }
    let i = xs.partition_point(|&v| v <= x) - 1;
    if n == 2 {
        return ys[0] + (x - xs[0]) / (xs[1] - xs[0]) * (ys[1] - ys[0]);
    }
    let h = |j: usize| xs[j + 1] - xs[j];
    let s = |j: usize| (ys[j + 1] - ys[j]) / h(j);
    let sign = |v: f64| {
        if v > 0.0 {
            1.0
        } else if v < 0.0 {
            -1.0
        } else {
            0.0
        }
    };
    let slope = |j: usize| -> f64 {
        if j == 0 {
            let p = s(0) * (1.0 + h(0) / (h(0) + h(1))) - s(1) * h(0) / (h(0) + h(1));
            if p * s(0) <= 0.0 {
                0.0
            } else if p.abs() > 2.0 * s(0).abs() {
                2.0 * s(0)
            } else {
                p
            }
        } else if j == n - 1 {
            let (a, b) = (n - 2, n - 3);
            let p = s(a) * (1.0 + h(a) / (h(a) + h(b))) - s(b) * h(a) / (h(a) + h(b));
            if p * s(a) <= 0.0 {
                0.0
            } else if p.abs() > 2.0 * s(a).abs() {
                2.0 * s(a)
            } else {
                p
            }
        } else {
            let (s0, s1) = (s(j - 1), s(j));
            let p = (s0 * h(j) + s1 * h(j - 1)) / (h(j - 1) + h(j));
            (sign(s0) + sign(s1)) * s0.abs().min(s1.abs()).min(0.5 * p.abs())
        }
    };
    let (y0, y1) = (slope(i), slope(i + 1));
    let hi = h(i);
    let si = s(i);
    let a = (y0 + y1 - 2.0 * si) / (hi * hi);
    let b = (3.0 * si - 2.0 * y0 - y1) / hi;
    let d = x - xs[i];
    ((a * d + b) * d + y0) * d + ys[i]
}

/// Evaluates a property of `material` (name for messages) at the
/// conditions. `space` allows x, y, z (the viscosity).
#[allow(clippy::too_many_arguments)]
fn evaluate(
    p: &Property,
    what: &str,
    owner: &str,
    m: &MaterialSection,
    env: &Env,
    cond: Conditions,
    space: bool,
    warnings: &mut Vec<Diagnostic>,
) -> Result<Resolved> {
    let key = format!("{owner}.{what}");
    let t = cond.temperature;
    let mut tenv = env.clone();
    tenv.constants.insert("T".into(), t);
    let allowed = if space { Vars::SPACE } else { Vars::NONE };
    let warn_range = |range: Option<[f64; 2]>, warnings: &mut Vec<Diagnostic>| -> f64 {
        if let Some([lo, hi]) = range
            && (t < lo || t > hi)
        {
            warnings.push(
                Diagnostic::new(
                    diag::MAT_RANGE,
                    format!(
                        "{key}: T = {t} K is outside the data range [{lo}, {hi}] K; \
                         the property is taken at {} K",
                        t.clamp(lo, hi)
                    ),
                )
                .key(key.clone()),
            );
            return t.clamp(lo, hi);
        }
        t
    };
    let num = |v: f64, origin: String| Resolved {
        value: Value::Number(v),
        origin,
    };
    let r = match p {
        Property::Number(v) => num(*v, format!("{v:?}")),
        Property::Expr(s) => {
            let v = expr::value(&Scalar::Expr(s.clone()), &tenv, allowed, &key).map_err(|e| {
                if !space && format!("{e:#}").contains("but this key takes no variable") {
                    mat_err(
                        diag::MAT_FORM,
                        &key,
                        format!(
                            "{key}: `{s}` varies in space; the density is one number for the \
                             whole fluid until the energy equation (L1) brings a variable density"
                        ),
                    )
                } else {
                    e
                }
            })?;
            let origin = format!("expression \"{s}\" at T = {t} K");
            match v.constant() {
                Some(c) => num(c, origin),
                None => Resolved { value: v, origin },
            }
        }
        Property::Form(f) => {
            let src = f.source().map(|s| format!("; {s}")).unwrap_or_default();
            let tt = warn_range(f.range(), warnings);
            let origin = |extra: &str| format!("{}{extra} at T = {tt} K{src}", f.name());
            match f {
                Form::Constant { value, .. } => num(*value, origin("")),
                Form::Polynomial { coefficients, .. } => {
                    let mut v = 0.0;
                    for c in coefficients.iter().rev() {
                        v = v * tt + c;
                    }
                    num(v, origin(""))
                }
                Form::Sutherland {
                    mu_ref, t_ref, s, ..
                } => {
                    let r = tt / t_ref;
                    num(
                        mu_ref * r * foehncraft_core::math::sqrt(r) * (t_ref + s) / (tt + s),
                        origin(""),
                    )
                }
                Form::PowerLaw {
                    value_ref,
                    t_ref,
                    exponent,
                    ..
                } => num(
                    value_ref * foehncraft_core::math::pow(tt / t_ref, *exponent),
                    origin(""),
                ),
                Form::Table {
                    temperature,
                    values,
                    file,
                    interpolation,
                    ..
                } => {
                    if file.is_some() {
                        return Err(mat_err(
                            diag::MAT_FORM,
                            &key,
                            format!(
                                "{key}: property tables from files come with the importers of D103; give `temperature` and `values`"
                            ),
                        ));
                    }
                    check_table(temperature, values, &key)?;
                    let v = match interpolation {
                        Interpolation::Linear => linear(temperature, values, tt),
                        Interpolation::MonotoneCubic => steffen(temperature, values, tt),
                    };
                    let name = match interpolation {
                        Interpolation::Linear => " (linear)",
                        Interpolation::MonotoneCubic => " (monotone-cubic)",
                    };
                    num(v, origin(name))
                }
                Form::Expression { expression, .. } => {
                    let inner = Property::Expr(expression.clone());
                    let mut r = evaluate(&inner, what, owner, m, env, cond, space, warnings)?;
                    r.origin.push_str(&src);
                    r
                }
                Form::IdealGas { molar_mass, .. } => {
                    let mm = molar_mass.or(m.molar_mass).ok_or_else(|| {
                        mat_err(
                            diag::MAT_FORM,
                            &key,
                            format!("{key}: the ideal-gas density needs `molar_mass`"),
                        )
                    })?;
                    num(
                        cond.pressure * mm / (R_GAS * t),
                        format!(
                            "ideal-gas, p = {} Pa, M = {mm} kg/mol, T = {t} K{src}",
                            cond.pressure
                        ),
                    )
                }
                Form::Nasa7 { .. } | Form::Nasa9 { .. } => {
                    return Err(mat_err(
                        diag::MAT_FORM,
                        &key,
                        format!(
                            "{key}: the {} form gives c_p, h and s and comes with the energy \
                             equation (L1)",
                            f.name()
                        ),
                    ));
                }
            }
        }
    };
    if let Some(v) = r.value.constant()
        && !(v > 0.0 && v.is_finite())
    {
        return Err(mat_err(
            diag::MAT_VALUE,
            &key,
            format!("{key} = {v} at T = {t} K must be positive and finite"),
        ));
    }
    Ok(r)
}

fn linear(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[n - 1] {
        return ys[n - 1];
    }
    let i = xs.partition_point(|&v| v <= x) - 1;
    ys[i] + (x - xs[i]) / (xs[i + 1] - xs[i]) * (ys[i + 1] - ys[i])
}

fn check_table(t: &[f64], v: &[f64], key: &str) -> Result<()> {
    if t.len() < 2 || t.len() != v.len() {
        return Err(mat_err(
            diag::MAT_FORM,
            key,
            format!(
                "{key}: a table needs `temperature` and `values` of the same length, at least 2 \
                 (has {} and {})",
                t.len(),
                v.len()
            ),
        ));
    }
    if t.iter().chain(v).any(|x| !x.is_finite()) || t.windows(2).any(|w| w[1] <= w[0]) {
        return Err(mat_err(
            diag::MAT_FORM,
            key,
            format!("{key}: the temperatures must be finite and increase strictly"),
        ));
    }
    Ok(())
}

/// The properties of a material (`None`: the shorthand of `[fluid]`), with
/// the `[fluid]` density and viscosity replacing the material's.
pub fn resolve_material(
    name: Option<&str>,
    case_materials: &BTreeMap<String, MaterialSection>,
    fluid_density: Option<&Scalar>,
    fluid_viscosity: Option<&Scalar>,
    env: &Env,
    cond: Conditions,
    warnings: &mut Vec<Diagnostic>,
) -> Result<MaterialProperties> {
    let (m, owner) = match name {
        Some(n) => (
            lookup(n, case_materials)?,
            if case_materials.contains_key(n) {
                format!("materials.{n}")
            } else {
                format!("material {n}")
            },
        ),
        None => (MaterialSection::default(), "fluid".to_string()),
    };
    if m.kind == Some(MaterialKind::Solid) {
        return Err(mat_err(
            diag::REG_NOT_FLUID,
            "fluid.material",
            format!(
                "`{}` is a solid; solids come with the energy equation (L1)",
                name.unwrap_or("")
            ),
        ));
    }
    let pick = |given: Option<&Scalar>, own: &Option<Property>, default: f64| -> (Property, bool) {
        match (given, own) {
            (Some(s), _) => (Property::from(s.clone()), true),
            (None, Some(p)) => (p.clone(), false),
            (None, None) => (Property::Number(default), false),
        }
    };
    let (rho_p, rho_given) = pick(fluid_density, &m.density, 1.0);
    let (mu_p, mu_given) = pick(fluid_viscosity, &m.viscosity, 1e-3);
    if name.is_some() {
        for (what, own, given) in [
            ("density", &m.density, rho_given),
            ("viscosity", &m.viscosity, mu_given),
        ] {
            if own.is_none() && !given {
                return Err(mat_err(
                    diag::MAT_FORM,
                    &format!("{owner}.{what}"),
                    format!("{owner} has no {what}; give it there or in [fluid] {what}"),
                ));
            }
        }
    }
    let rho_owner = if rho_given {
        "fluid".to_string()
    } else {
        owner.clone()
    };
    let mu_owner = if mu_given {
        "fluid".to_string()
    } else {
        owner.clone()
    };
    let density = evaluate(
        &rho_p, "density", &rho_owner, &m, env, cond, false, warnings,
    )?;
    let viscosity = evaluate(&mu_p, "viscosity", &mu_owner, &m, env, cond, true, warnings)?;
    Ok(MaterialProperties {
        material: name.map(str::to_string),
        density,
        viscosity,
    })
}

/// The resolved property program of a case (D94 (3)): the default fluid
/// and one entry per region.
#[derive(Debug, Clone, PartialEq)]
pub struct PropertySet {
    /// Conditions of the evaluation.
    pub conditions: Conditions,
    /// The `[fluid]` material (on every zone outside the regions).
    pub fluid: MaterialProperties,
    /// `(region name, zones, properties)` in name order.
    pub regions: Vec<(String, Vec<String>, MaterialProperties)>,
    /// Range warnings (MAT003).
    pub warnings: Vec<Diagnostic>,
}

impl PropertySet {
    /// The density (one number: the same in every region, checked).
    pub fn density(&self) -> f64 {
        self.fluid.density.value.constant().unwrap_or(f64::NAN)
    }

    /// The viscosity of the `[fluid]` material, as the reference value of
    /// Reynolds numbers and exact solutions; for a viscosity that varies in
    /// space, its value at the origin.
    pub fn reference_viscosity(&self) -> f64 {
        self.fluid.viscosity.value.at([0.0; 3], 0.0)
    }

    /// True when the viscosity is not one number on the whole mesh.
    pub fn viscosity_varies(&self) -> bool {
        let mu0 = self.fluid.viscosity.value.constant();
        mu0.is_none()
            || self
                .regions
                .iter()
                .any(|r| r.2.viscosity.value.constant().map(f64::to_bits) != mu0.map(f64::to_bits))
    }

    /// The program as text: what the resolved case shows and the case hash
    /// includes.
    pub fn text(&self) -> String {
        let mut s = String::new();
        let c = self.conditions;
        let _ = writeln!(
            s,
            "# resolved property program (D94): T = {} K, p = {} Pa, database {}",
            c.temperature,
            c.pressure,
            &database_hash()[..16]
        );
        let mut one = |label: &str, m: &MaterialProperties| {
            let mat = m.material.as_deref().unwrap_or("(fluid shorthand)");
            for (what, r) in [("density", &m.density), ("viscosity", &m.viscosity)] {
                let v = match &r.value {
                    Value::Number(v) => format!("{v:?}"),
                    Value::Program(p) => format!("\"{}\"", p.source),
                };
                let _ = writeln!(s, "#   {label} [{mat}] {what} = {v}  ({})", r.origin);
            }
        };
        one("fluid", &self.fluid);
        for (name, zones, m) in &self.regions {
            one(&format!("region {name} {zones:?}"), m);
        }
        s
    }

    /// The viscosity per cell of a mesh (zone by zone; expressions at the
    /// cell centroids).
    pub fn cell_viscosity(
        &self,
        mesh: &foehncraft_mesh::Mesh,
        geo: &foehncraft_mesh::Geometry,
    ) -> Vec<f64> {
        let zone_mu = self.zone_values(mesh);
        (0..mesh.n_cells)
            .map(|c| zone_mu[mesh.cell_zone[c] as usize].at(geo.cell_centroid[c], 0.0))
            .collect()
    }

    /// The viscosity of each cell zone of `mesh`.
    fn zone_values<'a>(&'a self, mesh: &foehncraft_mesh::Mesh) -> Vec<&'a Value> {
        mesh.zones
            .iter()
            .map(|z| {
                self.regions
                    .iter()
                    .find(|r| r.1.iter().any(|n| n == z))
                    .map_or(&self.fluid.viscosity.value, |r| &r.2.viscosity.value)
            })
            .collect()
    }

    /// The viscosity per face (`None` when it is one number everywhere,
    /// which keeps the solver's constant path bit for bit). Inside a zone
    /// an expression is evaluated at the face centroid; between zones the
    /// harmonic mean of the two cell values (D99).
    pub fn face_viscosity(
        &self,
        mesh: &foehncraft_mesh::Mesh,
        geo: &foehncraft_mesh::Geometry,
    ) -> Option<Arc<Vec<f64>>> {
        if !self.viscosity_varies() {
            return None;
        }
        let zone_mu = self.zone_values(mesh);
        let ni = mesh.n_interior_faces();
        let mut mu = Vec::with_capacity(mesh.n_faces());
        let cell = |c: usize| zone_mu[mesh.cell_zone[c] as usize].at(geo.cell_centroid[c], 0.0);
        for f in 0..mesh.n_faces() {
            let o = mesh.owner[f] as usize;
            let zo = mesh.cell_zone[o];
            if f < ni {
                let n = mesh.neighbour[f] as usize;
                if mesh.cell_zone[n] != zo {
                    let (a, b) = (cell(o), cell(n));
                    mu.push(2.0 * a * b / (a + b));
                    continue;
                }
            }
            mu.push(zone_mu[zo as usize].at(geo.face_centroid[f], 0.0));
        }
        Some(Arc::new(mu))
    }
}

/// The JSON description of a database material at the conditions (for
/// `foehn materials`).
pub fn describe(name: &str, cond: Conditions) -> Result<serde_json::Value> {
    let m = lookup(name, &BTreeMap::new())?;
    let mut warnings = Vec::new();
    let p = resolve_material(
        Some(name),
        &BTreeMap::new(),
        None,
        None,
        &Env::default(),
        cond,
        &mut warnings,
    )?;
    let prop = |r: &Resolved| {
        serde_json::json!({
            "value": r.value.constant(),
            "origin": r.origin,
        })
    };
    Ok(serde_json::json!({
        "name": name,
        "description": m.description,
        "kind": m.kind,
        "temperature_range": m.temperature_range,
        "molar_mass": m.molar_mass,
        "temperature": cond.temperature,
        "operating_pressure": cond.pressure,
        "density": prop(&p.density),
        "viscosity": prop(&p.viscosity),
        "forms": {
            "density": m.density.as_ref().map(|d| match d { Property::Form(f) => f.name(), Property::Number(_) => "number", Property::Expr(_) => "expression" }),
            "viscosity": m.viscosity.as_ref().map(|d| match d { Property::Form(f) => f.name(), Property::Number(_) => "number", Property::Expr(_) => "expression" }),
        },
        "warnings": warnings.iter().map(|w| w.message.clone()).collect::<Vec<_>>(),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn at(name: &str, t: f64) -> MaterialProperties {
        let mut w = Vec::new();
        resolve_material(
            Some(name),
            &BTreeMap::new(),
            None,
            None,
            &Env::default(),
            Conditions {
                temperature: t,
                pressure: DEFAULT_PRESSURE,
            },
            &mut w,
        )
        .unwrap()
    }

    #[test]
    fn database_parses_and_resolves() {
        for name in database().keys() {
            let m = at(name, 300.0);
            let (r, mu) = (
                m.density.value.constant().unwrap(),
                m.viscosity.value.constant().unwrap(),
            );
            assert!(r > 0.0 && mu > 0.0, "{name}: {r} {mu}");
        }
        // air at 300 K: ρ = 1.1766, μ = 1.846e-5 (White, Table A.2: 1.846e-5)
        let a = at("air", 300.0);
        let rho = a.density.value.constant().unwrap();
        let mu = a.viscosity.value.constant().unwrap();
        assert!((rho - 1.1766).abs() < 1e-3, "{rho}");
        assert!((mu - 1.846e-5).abs() < 0.005e-5, "{mu}");
        // water at the table points is exact
        let w = at("water", 293.15);
        assert_eq!(w.density.value.constant(), Some(998.21));
        assert_eq!(w.viscosity.value.constant(), Some(1.0016e-3));
    }

    #[test]
    fn steffen_is_monotone_and_exact_at_points() {
        let x = [0.0, 1.0, 2.0, 3.0, 4.0];
        let y = [0.0, 0.1, 5.0, 5.1, 5.2];
        for (xi, yi) in x.iter().zip(&y) {
            assert_eq!(steffen(&x, &y, *xi), *yi);
        }
        let mut prev = -1.0;
        for i in 0..=400 {
            let v = steffen(&x, &y, 4.0 * i as f64 / 400.0);
            assert!(v >= prev - 1e-15, "not monotone at {i}: {v} < {prev}");
            assert!((0.0..=5.2).contains(&v));
            prev = v;
        }
    }

    #[test]
    fn range_warning_and_base() {
        let mut w = Vec::new();
        let mut case = BTreeMap::new();
        case.insert(
            "thin-oil".to_string(),
            MaterialSection {
                base: Some("engine-oil".into()),
                viscosity: Some(Property::Number(0.05)),
                ..Default::default()
            },
        );
        let p = resolve_material(
            Some("thin-oil"),
            &case,
            None,
            None,
            &Env::default(),
            Conditions {
                temperature: 260.0,
                pressure: DEFAULT_PRESSURE,
            },
            &mut w,
        )
        .unwrap();
        assert_eq!(p.viscosity.value.constant(), Some(0.05));
        assert_eq!(p.density.value.constant(), Some(899.1));
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].code, diag::MAT_RANGE);
        // a cycle
        let mut case = BTreeMap::new();
        case.insert(
            "a".to_string(),
            MaterialSection {
                base: Some("b".into()),
                ..Default::default()
            },
        );
        case.insert(
            "b".to_string(),
            MaterialSection {
                base: Some("a".into()),
                ..Default::default()
            },
        );
        let e = resolve_material(
            Some("a"),
            &case,
            None,
            None,
            &Env::default(),
            Conditions::default(),
            &mut w,
        )
        .unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::MAT_BASE);
        // an unknown name with a hint
        let e = resolve_material(
            Some("watr"),
            &BTreeMap::new(),
            None,
            None,
            &Env::default(),
            Conditions::default(),
            &mut w,
        )
        .unwrap_err();
        assert!(format!("{e:#}").contains("did you mean `water`"));
    }

    #[test]
    fn a_varying_density_is_refused() {
        let mut w = Vec::new();
        let e = resolve_material(
            None,
            &BTreeMap::new(),
            Some(&Scalar::Expr("1000 + x".into())),
            None,
            &Env::default(),
            Conditions::default(),
            &mut w,
        )
        .unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::MAT_FORM);
        // T in an expression is the case temperature
        let p = resolve_material(
            None,
            &BTreeMap::new(),
            None,
            Some(&Scalar::Expr("1e-3 * T / 293.15".into())),
            &Env::default(),
            Conditions::default(),
            &mut w,
        )
        .unwrap();
        assert_eq!(p.viscosity.value.constant(), Some(1e-3));
    }
}

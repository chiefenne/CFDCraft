//! The JSON Schema of the flow case file (D77), generated from the Rust
//! structs that parse it (`schemars`, from the same types, serde
//! attributes and doc comments), so that it cannot drift from the reader.
//!
//! Every key has its type, its default, its range where the schema can say
//! it, and its description; two annotations are added:
//!
//! - `x-unit`: the SI unit (D44), taken from the description where it ends
//!   in a unit in parentheses (`Time step Δt (s)`);
//! - `x-runtime`: `true` for the keys a re-read of the case file may change
//!   during a run (D40, [`crate::control::RUNTIME_KEYS`]).
//!
//! [`markdown`] renders the key reference of `docs/keys.md` from it (a
//! test keeps the file up to date).

use crate::case::FlowCase;
use crate::control::RUNTIME_KEYS;
use serde_json::{Value, json};

/// Version of the schema envelope (D80).
pub const FORMAT_VERSION: &str = "1.0";

/// Units recognised at the end of a description, in parentheses.
const UNITS: &[&str] = &[
    "m", "m/s", "m²", "m³", "s", "Pa", "Pa s", "kg/m³", "kg/s", "N/m³", "N", "1/s", "bytes",
];

/// The JSON Schema (draft 2020-12) of the flow case file.
pub fn json_schema() -> Value {
    let mut v = schemars::schema_for!(FlowCase).to_value();
    if let Some(m) = v.as_object_mut() {
        m.insert(
            "$id".into(),
            json!(format!(
                "https://github.com/chiefenne/Foehncraft/schema/flow-case-{}.json",
                foehncraft_core::VERSION
            )),
        );
        m.insert("title".into(), json!("Foehncraft flow case"));
    }
    let defs = v.get("$defs").cloned().unwrap_or(Value::Null);
    let defs = annotate_defs(&v, &defs);
    if let Some(props) = v.get_mut("properties").and_then(Value::as_object_mut) {
        for (k, p) in props.iter_mut() {
            mark(p, k);
        }
    }
    if let Some(m) = v.as_object_mut()
        && !defs.is_null()
    {
        m.insert("$defs".into(), defs);
    }
    v
}

/// The schema with the format envelope of `foehn schema --json`.
pub fn envelope(schema: Value) -> Value {
    json!({
        "format": "foehncraft.schema",
        "format_version": FORMAT_VERSION,
        "foehncraft_version": foehncraft_core::VERSION,
        "schema": schema,
    })
}

/// The unit at the end of a description: the first parenthesis whose
/// content is a known unit.
pub fn unit_of(description: &str) -> Option<&'static str> {
    let mut rest = description;
    while let Some(i) = rest.find('(') {
        let after = &rest[i + 1..];
        let Some(j) = after.find(')') else { break };
        let inner = after[..j].trim();
        let inner = inner.split(';').next().unwrap_or(inner).trim();
        if let Some(u) = UNITS.iter().find(|u| **u == inner) {
            return Some(u);
        }
        rest = &after[j + 1..];
    }
    None
}

fn resolve<'a>(s: &'a Value, defs: &'a Value) -> &'a Value {
    if let Some(r) = s.get("$ref").and_then(Value::as_str)
        && let Some(name) = r.strip_prefix("#/$defs/")
        && let Some(d) = defs.get(name)
    {
        return d;
    }
    s
}

/// The schema of a struct-like value inside `s` (through `$ref`, an
/// optional `anyOf` with null, or `allOf`).
fn object_of<'a>(s: &'a Value, defs: &'a Value) -> Option<&'a Value> {
    let s = resolve(s, defs);
    let tagged = s
        .get("oneOf")
        .and_then(Value::as_array)
        .is_some_and(|a| a.iter().any(|x| x.get("properties").is_some()));
    if s.get("properties").is_some() || tagged {
        return Some(s);
    }
    for k in ["anyOf", "allOf"] {
        if let Some(a) = s.get(k).and_then(Value::as_array) {
            for x in a {
                if let Some(o) = object_of(x, defs) {
                    return Some(o);
                }
            }
        }
    }
    None
}

/// Adds `x-unit` and `x-runtime` to a property of the root.
fn mark(p: &mut Value, path: &str) {
    let Some(m) = p.as_object_mut() else { return };
    let desc = m
        .get("description")
        .and_then(Value::as_str)
        .unwrap_or("")
        .to_string();
    if let Some(u) = unit_of(&desc) {
        m.insert("x-unit".into(), json!(u));
    }
    if RUNTIME_KEYS.contains(&path) {
        m.insert("x-runtime".into(), json!(true));
    }
}

/// (definition, property, unit, runtime) for every property reached from
/// the root.
fn collect(
    s: &Value,
    defs: &Value,
    path: &str,
    def: Option<&str>,
    out: &mut Vec<(String, String, Option<&'static str>, bool)>,
) {
    let target = resolve(s, defs);
    let def_name = s
        .get("$ref")
        .and_then(Value::as_str)
        .and_then(|r| r.strip_prefix("#/$defs/"))
        .or(def);
    if let Some(props) = target.get("properties").and_then(Value::as_object) {
        for (k, p) in props {
            let sub = if path.is_empty() {
                k.clone()
            } else {
                format!("{path}.{k}")
            };
            if let Some(d) = def_name {
                let desc = p.get("description").and_then(Value::as_str).unwrap_or("");
                out.push((
                    d.to_string(),
                    k.clone(),
                    unit_of(desc),
                    RUNTIME_KEYS.contains(&sub.as_str()),
                ));
            }
            for inner in refs_in(p) {
                collect(inner, defs, &sub, None, out);
            }
        }
    }
    for k in ["oneOf", "anyOf", "allOf"] {
        if let Some(a) = target.get(k).and_then(Value::as_array) {
            for x in a {
                collect(x, defs, path, def_name, out);
            }
        }
    }
    if let Some(ap) = target.get("additionalProperties")
        && ap.is_object()
    {
        for inner in refs_in(ap) {
            collect(inner, defs, &format!("{path}.<name>"), None, out);
        }
    }
}

/// The `$ref` schemas directly inside a property (itself, `anyOf`, array
/// items, map values).
fn refs_in(p: &Value) -> Vec<&Value> {
    let mut v = Vec::new();
    if p.get("$ref").is_some() {
        v.push(p);
    }
    for k in ["anyOf", "oneOf", "allOf"] {
        if let Some(a) = p.get(k).and_then(Value::as_array) {
            for x in a {
                v.extend(refs_in(x));
            }
        }
    }
    if let Some(i) = p.get("items") {
        v.extend(refs_in(i));
    }
    if let Some(ap) = p.get("additionalProperties") {
        v.extend(refs_in(ap));
    }
    v
}

fn annotate_defs(root: &Value, defs: &Value) -> Value {
    let mut marks = Vec::new();
    collect(root, defs, "", None, &mut marks);
    let mut defs = defs.clone();
    for (d, k, unit, runtime) in marks {
        let Some(p) = defs
            .get_mut(&d)
            .and_then(|x| x.get_mut("properties"))
            .and_then(|x| x.get_mut(&k))
            .and_then(Value::as_object_mut)
        else {
            continue;
        };
        if let Some(u) = unit {
            p.insert("x-unit".into(), json!(u));
        }
        if runtime {
            p.insert("x-runtime".into(), json!(true));
        }
    }
    // tagged variants (oneOf entries) carry their properties inline
    if let Some(m) = defs.as_object_mut() {
        for (_, d) in m.iter_mut() {
            if let Some(a) = d.get_mut("oneOf").and_then(Value::as_array_mut) {
                for x in a {
                    if let Some(props) = x.get_mut("properties").and_then(Value::as_object_mut) {
                        for (_, p) in props.iter_mut() {
                            let desc = p.get("description").and_then(Value::as_str).unwrap_or("");
                            if let Some(u) = unit_of(desc)
                                && let Some(pm) = p.as_object_mut()
                            {
                                pm.insert("x-unit".into(), json!(u));
                            }
                        }
                    }
                }
            }
        }
    }
    defs
}

/// One row of the key reference.
#[derive(Debug, Clone)]
pub struct KeyRow {
    /// Dotted path (`<name>` for a map key, `[]` for an array of tables).
    pub path: String,
    /// Type.
    pub ty: String,
    /// Default (JSON).
    pub default: String,
    /// Unit.
    pub unit: String,
    /// May change during a run.
    pub runtime: bool,
    /// Description.
    pub description: String,
}

/// Every key of the case file, in schema order.
pub fn keys() -> Vec<KeyRow> {
    let s = json_schema();
    let defs = s.get("$defs").cloned().unwrap_or(Value::Null);
    let mut out = Vec::new();
    walk(&s, &defs, "", &mut out, 0);
    out
}

fn type_of(p: &Value, defs: &Value) -> String {
    let p = resolve(p, defs);
    if let Some(e) = p.get("enum").and_then(Value::as_array) {
        return e
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(" \\| ");
    }
    if let Some(c) = p.get("const") {
        return c.to_string();
    }
    if let Some(a) = p
        .get("oneOf")
        .or_else(|| p.get("anyOf"))
        .and_then(Value::as_array)
    {
        let parts: Vec<String> = a
            .iter()
            .map(|x| type_of(x, defs))
            .filter(|t| t != "null")
            .collect();
        return parts.join(" \\| ");
    }
    let is_array = match p.get("type") {
        Some(Value::String(t)) => t == "array",
        Some(Value::Array(a)) => a.iter().any(|t| t == "array"),
        _ => false,
    };
    match p.get("type") {
        _ if is_array => {
            let item = p.get("items").map_or("value".into(), |i| type_of(i, defs));
            match (p.get("minItems"), p.get("maxItems")) {
                (Some(a), Some(b)) if a == b => format!("[{item}; {a}]"),
                _ => format!("[{item}]"),
            }
        }
        Some(Value::String(t)) => t.clone(),
        Some(Value::Array(a)) => a
            .iter()
            .filter_map(Value::as_str)
            .filter(|t| *t != "null")
            .collect::<Vec<_>>()
            .join(" \\| "),
        _ => "table".into(),
    }
}

fn first_sentence(d: &str) -> String {
    let d = d.replace('\n', " ");
    d.trim().to_string()
}

fn walk(s: &Value, defs: &Value, path: &str, out: &mut Vec<KeyRow>, depth: usize) {
    if depth > 12 {
        return;
    }
    let t = resolve(s, defs);
    let required: Vec<&str> = t
        .get("required")
        .and_then(Value::as_array)
        .map(|a| a.iter().filter_map(Value::as_str).collect())
        .unwrap_or_default();
    if let Some(props) = t.get("properties").and_then(Value::as_object) {
        for (k, p) in props {
            let sub = if path.is_empty() {
                k.clone()
            } else {
                format!("{path}.{k}")
            };
            let desc = p
                .get("description")
                .or_else(|| resolve(p, defs).get("description"))
                .and_then(Value::as_str)
                .unwrap_or("");
            // a table: recurse
            let is_array = resolve(p, defs).get("type") == Some(&json!("array"))
                || p.get("anyOf").is_some_and(|a| {
                    a.as_array()
                        .is_some_and(|a| a.iter().any(|x| x.get("type") == Some(&json!("array"))))
                });
            let items = resolve(p, defs).get("items");
            if let Some(o) = object_of(p, defs)
                && !is_array
            {
                walk(o, defs, &sub, out, depth + 1);
                continue;
            }
            if is_array
                && let Some(i) = items
                && let Some(o) = object_of(i, defs)
            {
                walk(o, defs, &format!("{sub}[]"), out, depth + 1);
                continue;
            }
            if let Some(ap) = resolve(p, defs).get("additionalProperties")
                && let Some(o) = object_of(ap, defs)
            {
                walk(o, defs, &format!("{sub}.<name>"), out, depth + 1);
                continue;
            }
            let default = p.get("default").map(|d| d.to_string()).unwrap_or_else(|| {
                if required.contains(&k.as_str()) {
                    "required".into()
                } else {
                    String::new()
                }
            });
            out.push(KeyRow {
                path: sub.clone(),
                ty: type_of(p, defs),
                default,
                unit: p
                    .get("x-unit")
                    .and_then(Value::as_str)
                    .map(str::to_string)
                    .or_else(|| unit_of(desc).map(str::to_string))
                    .unwrap_or_default(),
                runtime: RUNTIME_KEYS.contains(&sub.as_str()),
                description: first_sentence(desc),
            });
        }
    }
    // tagged enums: one block per variant
    if let Some(a) = t.get("oneOf").and_then(Value::as_array) {
        for x in a {
            let tag = x
                .get("properties")
                .and_then(|p| p.get("type"))
                .and_then(|t| t.get("const"))
                .and_then(Value::as_str);
            match tag {
                Some(tag) => {
                    let mut x = x.clone();
                    if let Some(p) = x.get_mut("properties").and_then(Value::as_object_mut) {
                        p.remove("type");
                    }
                    let before = out.len();
                    walk(
                        &x,
                        defs,
                        &format!("{path} (type = \"{tag}\")"),
                        out,
                        depth + 1,
                    );
                    if out.len() == before {
                        out.push(KeyRow {
                            path: format!("{path} (type = \"{tag}\")"),
                            ty: "table".into(),
                            default: String::new(),
                            unit: String::new(),
                            runtime: false,
                            description: x
                                .get("description")
                                .and_then(Value::as_str)
                                .map(first_sentence)
                                .unwrap_or_default(),
                        });
                    }
                }
                None => walk(x, defs, path, out, depth + 1),
            }
        }
    }
}

/// The key reference (`docs/keys.md`).
pub fn markdown() -> String {
    let mut s = String::from(
        "# Flow case keys\n\n\
         Generated by `foehn schema --markdown` from the structs that read the case file \
         (D77); do not edit by hand. Units are SI (D44). *Run time*: the key may change \
         during a run, by re-reading the case file (D40). `<name>` is a patch name, `[]` an \
         array of tables. The JSON Schema itself: `foehn schema`.\n\n\
         | Key | Type | Default | Unit | Run time | Description |\n\
         |---|---|---|---|---|---|\n",
    );
    for k in keys() {
        let esc = |x: &str| x.replace('|', "\\|");
        s.push_str(&format!(
            "| `{}` | {} | {} | {} | {} | {} |\n",
            k.path,
            k.ty,
            if k.default.is_empty() {
                String::new()
            } else {
                format!("`{}`", k.default)
            },
            k.unit,
            if k.runtime { "yes" } else { "" },
            esc(&k.description)
        ));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn find<'a>(rows: &'a [KeyRow], path: &str) -> &'a KeyRow {
        rows.iter()
            .find(|r| r.path == path)
            .unwrap_or_else(|| panic!("no key {path}"))
    }

    #[test]
    fn schema_has_types_defaults_units_and_runtime_flags() {
        let rows = keys();
        let d = find(&rows, "fluid.density");
        assert_eq!(d.unit, "kg/m³");
        // no default of its own: 1, or the material's (D100)
        assert_eq!(d.default, "");
        assert!(!d.runtime);
        let t = find(&rows, "solver.tolerance");
        assert!(t.runtime);
        assert_eq!(t.default, "1e-8");
        assert_eq!(find(&rows, "time.time_step").unit, "s");
        assert_eq!(find(&rows, "run.limits.max_wall_time").unit, "s");
        assert!(find(&rows, "run.limits.max_steps").runtime);
        assert!(find(&rows, "mesh.renumber").ty.contains("\"rcm\""));
        // every run-time key exists in the schema
        for k in RUNTIME_KEYS {
            assert!(
                rows.iter().any(|r| r.path == *k),
                "run-time key {k} not in the schema"
            );
        }
        // boundary conditions per variant
        assert!(
            rows.iter()
                .any(|r| r.path.starts_with("boundary.<name> (type = \"wall\")"))
        );
        // the annotations are in the schema itself
        let s = json_schema();
        let text = s.to_string();
        assert!(text.contains("\"x-runtime\":true"));
        assert!(text.contains("\"x-unit\":\"kg/m³\""));
    }

    #[test]
    fn units_from_descriptions() {
        assert_eq!(
            unit_of("Time step Δt (s); with `adaptive` the first one."),
            Some("s")
        );
        assert_eq!(unit_of("Density ρ (kg/m³)."), Some("kg/m³"));
        assert_eq!(unit_of("Number of points (one per line)."), None);
    }

    #[test]
    fn key_reference_is_up_to_date() {
        // a Windows checkout may have CRLF line ends
        let doc = include_str!("../../../docs/keys.md").replace("\r\n", "\n");
        assert!(
            doc == markdown(),
            "docs/keys.md is out of date: regenerate it with\n  \
             cargo run --release -p foehncraft-cli -- schema --markdown > docs/keys.md"
        );
    }
}

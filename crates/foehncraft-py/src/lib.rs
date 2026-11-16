//! The native module `foehncraft._native` of the Python package (D61).
//!
//! A thin layer over the engine API of `foehncraft-flow` (D86): check and
//! run a case file in process, see every row of the run record while the
//! run goes on (monitors and stop criteria at iteration or step level),
//! evaluate expressions (D53) with the solver's own compiler and
//! functions, list the material database (D102) and return the case-file
//! schema (D77). Everything crosses the boundary as JSON text, the same
//! versioned formats the command line writes (D80), so the Python side
//! and `foehn --json` see the same data.
//!
//! The run releases the GIL; the callbacks take it back between outer
//! iterations or time steps, never per cell, so GPU speed and
//! reproducibility are kept (D53, D31).

use foehncraft_flow::case::output_dir;
use foehncraft_flow::engine::{
    self, BuildInfo, Interrupts, Observer, Roots, RunOptions, StopSignal,
};
use foehncraft_flow::expr::{self, Env, Scalar, TableSection, Vars};
use foehncraft_flow::{diag, materials, schema};
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

create_exception!(
    _native,
    FoehncraftError,
    PyException,
    "An error of the solver with its diagnostic code (D78): args = (code, message)."
);

/// The error of the solver as a Python exception with its code.
fn py_err(e: anyhow::Error) -> PyErr {
    let d = diag::of_error(&e);
    FoehncraftError::new_err((d.code.to_string(), format!("{e:#}")))
}

/// Build and host for the provenance (D85).
fn build_info() -> BuildInfo {
    let mut features = vec!["python".to_string()];
    if cfg!(feature = "hdf5") {
        features.push("hdf5".into());
    }
    BuildInfo {
        commit: option_env!("FOEHNCRAFT_GIT_COMMIT")
            .unwrap_or("unknown")
            .into(),
        features,
        platform: format!("{}-{}", std::env::consts::OS, std::env::consts::ARCH),
        cpu: String::new(),
        gpu: None,
    }
}

/// The version of the solver the module was built from.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// Checks a case file: validation, boundary conditions on the first mesh
/// and the estimate. Returns the check report (`foehncraft.check`) as JSON.
#[pyfunction]
#[pyo3(signature = (path, sets = Vec::new(), out = None))]
fn check(py: Python<'_>, path: PathBuf, sets: Vec<String>, out: Option<PathBuf>) -> String {
    py.allow_threads(|| {
        let rep = engine::check_file(
            &path,
            &sets,
            &Roots::default(),
            out.as_deref(),
            &build_info(),
        );
        serde_json::to_string(&rep).unwrap_or_default()
    })
}

/// Runs a case file in this process and returns the result summary
/// (`foehncraft.result`, also written to `result.json`) as JSON.
///
/// `observer(row)` gets every row of the run record as JSON text (the rows
/// of `run.jsonl`: header, iteration, step, balance, event, end); a true
/// return value stops the run as a stop request does (restart file, final
/// output, exit code 3). `log(line)` gets the progress lines. An exception
/// in a callback stops the run and is raised after it. Ctrl-C stops the
/// run the same way.
#[pyfunction]
#[pyo3(signature = (path, sets = Vec::new(), out = None, observer = None, log = None))]
fn run(
    py: Python<'_>,
    path: PathBuf,
    sets: Vec<String>,
    out: Option<PathBuf>,
    observer: Option<PyObject>,
    log: Option<PyObject>,
) -> PyResult<String> {
    let roots = Roots::default();
    let case = engine::load(&path, &sets, &roots).map_err(py_err)?;
    let (name, dir) = output_dir(&case, &path, out);
    let interrupts = Interrupts::new();
    let failed: Arc<Mutex<Option<PyErr>>> = Arc::new(Mutex::new(None));
    let observer = observer.map(|cb| {
        let stop = interrupts.clone();
        let failed = failed.clone();
        Observer::new(move |row: &serde_json::Value| {
            Python::with_gil(|py| {
                let text = row.to_string();
                match cb.call1(py, (text,)).and_then(|r| r.is_truthy(py)) {
                    Ok(true) => stop.request_stop(StopSignal::Api),
                    Ok(false) => {}
                    Err(e) => {
                        failed.lock().unwrap().get_or_insert(e);
                        stop.request_stop(StopSignal::Api);
                    }
                }
                if py.check_signals().is_err() {
                    stop.request_stop(StopSignal::Int);
                }
            })
        })
    });
    let opts = RunOptions {
        interrupts: interrupts.clone(),
        build: build_info(),
        roots,
        observer,
    };
    let lg_failed = failed.clone();
    let lg_stop = interrupts.clone();
    let res = py.allow_threads(move || {
        engine::run(&case, &name, &dir, &opts, |line: &str| {
            Python::with_gil(|py| {
                if let Some(cb) = &log
                    && let Err(e) = cb.call1(py, (line,))
                {
                    lg_failed.lock().unwrap().get_or_insert(e);
                    lg_stop.request_stop(StopSignal::Api);
                }
                if py.check_signals().is_err() {
                    lg_stop.request_stop(StopSignal::Int);
                }
            })
        })
    });
    if let Some(e) = failed.lock().unwrap().take() {
        return Err(e);
    }
    let o = res.map_err(py_err)?;
    serde_json::to_string(&o.summary).map_err(|e| py_err(e.into()))
}

/// The constants and tables an expression may refer to, as in a case file.
#[derive(serde::Deserialize, Default)]
#[serde(deny_unknown_fields)]
struct EnvText {
    #[serde(default)]
    constants: BTreeMap<String, Scalar>,
    #[serde(default)]
    tables: BTreeMap<String, TableSection>,
}

fn variables(v: &str) -> PyResult<Vars> {
    let mut vars = Vars::NONE;
    for c in v.chars() {
        vars = vars.with(match c {
            'x' => Vars::X,
            'y' => Vars::Y,
            'z' => Vars::Z,
            't' => Vars::T,
            _ => {
                return Err(FoehncraftError::new_err((
                    "EXPR003".to_string(),
                    format!("unknown variable `{c}` (x, y, z, t)"),
                )));
            }
        });
    }
    Ok(vars)
}

/// Compiles the expression `source` with the `[constants]` and `[tables]`
/// of the TOML text `env` (relative table files from `base`) and the
/// variables it may use (a string of `x`, `y`, `z`, `t`), and evaluates it
/// at the `points` (x, y, z, t). Returns the values and the variables the
/// expression uses. Errors carry the codes EXPR001–EXPR005.
#[pyfunction]
#[pyo3(signature = (source, env = "", allowed = "xyzt", points = Vec::new(), base = None))]
fn expression(
    source: &str,
    env: &str,
    allowed: &str,
    points: Vec<[f64; 4]>,
    base: Option<PathBuf>,
) -> PyResult<(Vec<f64>, String)> {
    let mut e: EnvText = toml::from_str(env).map_err(|e| {
        FoehncraftError::new_err(("CASE002".to_string(), format!("constants and tables: {e}")))
    })?;
    let base = base.unwrap_or_else(|| PathBuf::from("."));
    for (name, t) in &mut e.tables {
        if let Some(f) = &t.file
            && f.is_relative()
        {
            t.file = Some(base.join(f));
        }
        t.load(name).map_err(py_err)?;
    }
    let env = Env::new(&e.constants, &e.tables).map_err(py_err)?;
    let p = expr::compile(source, &env, variables(allowed)?, "expression").map_err(py_err)?;
    let values = points
        .iter()
        .map(|q| expr::Value::Program(Arc::new(p.clone())).at([q[0], q[1], q[2]], q[3]))
        .collect();
    let uses = [
        ('x', Vars::X),
        ('y', Vars::Y),
        ('z', Vars::Z),
        ('t', Vars::T),
    ]
    .iter()
    .filter(|(_, v)| p.uses.contains(*v))
    .map(|(c, _)| *c)
    .collect();
    Ok((values, uses))
}

/// The material database (`foehncraft.materials`) as JSON, every entry
/// evaluated at `temperature` (K) and `pressure` (Pa); `name` picks one.
#[pyfunction]
#[pyo3(signature = (name = None, temperature = 293.15, pressure = 101325.0))]
fn material_database(name: Option<String>, temperature: f64, pressure: f64) -> PyResult<String> {
    if !(temperature > 0.0 && pressure > 0.0) {
        return Err(FoehncraftError::new_err((
            "MAT005".to_string(),
            "temperature and pressure must be > 0".to_string(),
        )));
    }
    let cond = materials::Conditions {
        temperature,
        pressure,
    };
    let names: Vec<String> = match name {
        Some(n) => vec![n],
        None => materials::database().keys().cloned().collect(),
    };
    let mut list = Vec::new();
    for n in &names {
        list.push(materials::describe(n, cond).map_err(py_err)?);
    }
    Ok(serde_json::json!({
        "format": "foehncraft.materials",
        "format_version": engine::FORMAT_VERSION,
        "database_sha256": materials::database_hash(),
        "materials": list,
    })
    .to_string())
}

/// The JSON Schema of the flow case file (D77), as `foehn schema` writes it.
#[pyfunction]
fn case_schema() -> String {
    schema::envelope(schema::json_schema()).to_string()
}

/// The diagnostic codes (D78) with their severity and meaning, as JSON.
#[pyfunction]
fn codes() -> String {
    let list: Vec<serde_json::Value> = diag::CODES
        .iter()
        .map(|c| serde_json::json!({"code": c.code, "severity": c.severity, "title": c.title}))
        .collect();
    serde_json::Value::Array(list).to_string()
}

/// Reads a case file the way a run does and returns it, with the defaults
/// filled in and the overrides applied, as TOML text (the
/// `case.resolved.toml` of a run).
#[pyfunction]
#[pyo3(signature = (path, sets = Vec::new()))]
fn resolved(path: PathBuf, sets: Vec<String>) -> PyResult<String> {
    let case = engine::load(Path::new(&path), &sets, &Roots::default()).map_err(py_err)?;
    case.resolved_toml().map_err(py_err)
}

/// The native module.
#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("FoehncraftError", m.py().get_type::<FoehncraftError>())?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(check, m)?)?;
    m.add_function(wrap_pyfunction!(run, m)?)?;
    m.add_function(wrap_pyfunction!(expression, m)?)?;
    m.add_function(wrap_pyfunction!(material_database, m)?)?;
    m.add_function(wrap_pyfunction!(case_schema, m)?)?;
    m.add_function(wrap_pyfunction!(codes, m)?)?;
    m.add_function(wrap_pyfunction!(resolved, m)?)?;
    Ok(())
}

//! The run contract of `foehn flow` (D62, D64, D78, D81, D84): JSON on
//! standard output, `result.json`, exit codes, `--check --json`, limits,
//! `--restart latest` and, on unix, a clean stop on SIGTERM.

use serde_json::Value;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const FOEHN: &str = env!("CARGO_BIN_EXE_foehn");

/// A fresh directory for one test.
fn dir(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("run_contract")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

const MESH_AND_WALLS: &str = r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.0625]
cells = [16, 16, 1]

[fluid]
density = 1.0
viscosity = 0.01

[boundary.ymax]
type = "wall"
velocity = [1.0, 0.0, 0.0]
[boundary.xmin]
type = "wall"
[boundary.xmax]
type = "wall"
[boundary.ymin]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
"#;

/// A steady cavity on 16 x 16 cells.
fn steady(d: &Path) -> PathBuf {
    let p = d.join("cavity.toml");
    let text = format!(
        "{MESH_AND_WALLS}\n[solver]\nrelaxation_velocity = 0.9\ntolerance = 1e-8\nmax_iterations = 2000\n"
    );
    std::fs::write(&p, text).unwrap();
    p
}

/// A transient cavity on 16 x 16 cells.
#[cfg(unix)]
fn transient(d: &Path, steps: usize) -> PathBuf {
    let p = d.join("cavity_transient.toml");
    let text = format!(
        "{MESH_AND_WALLS}\n[solver]\nlog_every = 1000\n\n[time]\nscheme = \"bdf2\"\ntime_step = 0.01\nsteps = {steps}\n"
    );
    std::fs::write(&p, text).unwrap();
    p
}

fn foehn(args: &[&str]) -> Output {
    Command::new(FOEHN).args(args).output().expect("foehn runs")
}

fn json(o: &Output) -> Value {
    serde_json::from_slice(&o.stdout).unwrap_or_else(|e| {
        panic!(
            "stdout is not one JSON object ({e}):\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&o.stdout),
            String::from_utf8_lossy(&o.stderr)
        )
    })
}

fn bytes(p: &Path) -> Vec<u8> {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

#[test]
fn json_result_and_byte_identical_files() {
    let d = dir("identical");
    let case = steady(&d);
    let (a, b) = (d.join("a"), d.join("b"));
    for out in [&a, &b] {
        let o = foehn(&["--json", "--quiet", "flow", s(&case), "--out", s(out)]);
        assert_eq!(
            o.status.code(),
            Some(0),
            "{}",
            String::from_utf8_lossy(&o.stderr)
        );
        let v = json(&o);
        assert_eq!(v["format"], "foehncraft.result");
        assert_eq!(v["status"], "finished");
        assert_eq!(v["exit_code"], 0);
        assert_eq!(v["solves"][0]["cells"], 256);
        assert!(v["provenance"]["case_hash"].as_str().unwrap().len() == 64);
        // result.json is the object on standard output
        let file: Value = serde_json::from_slice(&bytes(&out.join("result.json"))).unwrap();
        assert_eq!(file, v);
        // the run record starts with the header and ends with the end row
        let rec = std::fs::read_to_string(out.join("run.jsonl")).unwrap();
        let rows: Vec<Value> = rec
            .lines()
            .map(|l| serde_json::from_str(l).unwrap())
            .collect();
        assert_eq!(rows[0]["type"], "header");
        assert_eq!(rows[0]["format"], "foehncraft.run_record");
        assert!(rows.iter().any(|r| r["type"] == "iteration"));
        assert_eq!(rows.last().unwrap()["type"], "end");
    }
    for f in ["result.json", "fields.vtkhdf", "residuals.csv"] {
        assert!(
            bytes(&a.join(f)) == bytes(&b.join(f)),
            "{f} differs between two runs"
        );
    }
    // result.json lists the hash of every output except the measured files
    let v: Value = serde_json::from_slice(&bytes(&a.join("result.json"))).unwrap();
    let names: Vec<&str> = v["files"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| f["name"].as_str().unwrap())
        .collect();
    assert!(names.contains(&"fields.vtkhdf"), "{names:?}");
    assert!(
        !names.contains(&"flow.log") && !names.contains(&"run.jsonl"),
        "{names:?}"
    );
}

#[test]
fn exit_codes() {
    let d = dir("exit_codes");
    let case = steady(&d);
    // 1: a value out of range, as a JSON error object with code and key
    let o = foehn(&[
        "--json",
        "--quiet",
        "flow",
        s(&case),
        "--set",
        "fluid.viscosity=-1",
        "--out",
        s(&d.join("e1")),
    ]);
    assert_eq!(o.status.code(), Some(1));
    let v = json(&o);
    assert_eq!(v["format"], "foehncraft.error");
    assert_eq!(v["diagnostics"][0]["code"], "CASE004");
    assert_eq!(v["diagnostics"][0]["key"], "fluid.viscosity");
    // 1: an unknown key, with its location and a hint
    let bad = d.join("typo.toml");
    std::fs::write(
        &bad,
        std::fs::read_to_string(&case)
            .unwrap()
            .replace("cells = [16", "cels = [16"),
    )
    .unwrap();
    let o = foehn(&[
        "--json",
        "--quiet",
        "flow",
        s(&bad),
        "--out",
        s(&d.join("e1b")),
    ]);
    assert_eq!(o.status.code(), Some(1));
    let v = json(&o);
    assert_eq!(v["diagnostics"][0]["code"], "CASE002");
    assert_eq!(v["diagnostics"][0]["key"], "mesh.generate.cels");
    // 2: the run ends, but it did not converge
    let o = foehn(&[
        "--json",
        "--quiet",
        "flow",
        s(&case),
        "--set",
        "solver.max_iterations=5",
        "--set",
        "verification.require_convergence=true",
        "--out",
        s(&d.join("e2")),
    ]);
    assert_eq!(
        o.status.code(),
        Some(2),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let v = json(&o);
    assert_eq!(v["status"], "failed-checks");
    assert!(
        v["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|d| d["code"] == "RUN006")
    );
    // 3: stopped at a limit
    let o = foehn(&[
        "--json",
        "--quiet",
        "flow",
        s(&case),
        "--set",
        "run.limits.max_iterations=10",
        "--out",
        s(&d.join("e3")),
    ]);
    assert_eq!(o.status.code(), Some(3));
    let v = json(&o);
    assert_eq!(v["status"], "stopped");
    assert_eq!(v["stop"]["code"], "RUN004");
    // 3: the .stop control file of an earlier run blocks the start
    std::fs::write(d.join("e3").join(".stop"), "").unwrap();
    let o = foehn(&[
        "--json",
        "--quiet",
        "flow",
        s(&case),
        "--out",
        s(&d.join("e3")),
    ]);
    assert_eq!(o.status.code(), Some(1));
}

#[test]
fn check_json_with_estimate() {
    let d = dir("check");
    let case = steady(&d);
    let o = foehn(&["--json", "flow", s(&case), "--check"]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let v = json(&o);
    assert_eq!(v["format"], "foehncraft.check");
    assert_eq!(v["ok"], true);
    // the pressure level is held in a reference cell: a note with its code
    assert!(
        v["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|d| d["code"] == "BC015")
    );
    let e = &v["estimate"];
    assert_eq!(e["cells"], 256);
    assert!(e["host_memory"].as_f64().unwrap() > 1e6);
    assert!(e["output_bytes"].as_f64().unwrap() > 0.0);
    // a broken case is a report with ok = false, not a crash
    let o = foehn(&[
        "--json",
        "flow",
        s(&case),
        "--check",
        "--set",
        "boundary.xmin.type=\"nozzle\"",
    ]);
    assert_ne!(o.status.code(), Some(0));
    let v = json(&o);
    assert_eq!(v["ok"], false);
}

#[test]
fn schema_is_json_schema() {
    let o = foehn(&["schema"]);
    assert_eq!(o.status.code(), Some(0));
    let v: Value = serde_json::from_slice(&o.stdout).unwrap();
    assert!(v["$schema"].as_str().unwrap().contains("2020-12"));
    assert!(v["properties"]["fluid"].is_object());
}

#[test]
fn steady_restart_latest_is_bitwise() {
    let d = dir("steady_restart");
    let case = steady(&d);
    let (straight, split) = (d.join("straight"), d.join("split"));
    let o = foehn(&["--quiet", "flow", s(&case), "--out", s(&straight)]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let o = foehn(&[
        "--quiet",
        "flow",
        s(&case),
        "--set",
        "run.limits.max_iterations=40",
        "--out",
        s(&split),
    ]);
    assert_eq!(o.status.code(), Some(3));
    assert!(split.join("restart.h5").exists());
    let o = foehn(&[
        "--quiet",
        "flow",
        s(&case),
        "--restart",
        "latest",
        "--out",
        s(&split),
    ]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    for f in ["fields.vtkhdf", "restart.h5"] {
        assert!(
            bytes(&straight.join(f)) == bytes(&split.join(f)),
            "{f}: stop and restart differs from the straight run"
        );
    }
}

#[cfg(unix)]
#[test]
fn sigterm_stops_cleanly_and_restart_latest_continues_bitwise() {
    let d = dir("sigterm");
    let case = transient(&d, 600);
    let (straight, split) = (d.join("straight"), d.join("split"));
    let o = foehn(&["--quiet", "flow", s(&case), "--out", s(&straight)]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );

    let mut child = Command::new(FOEHN)
        .args(["--quiet", "flow", s(&case), "--out", s(&split)])
        .spawn()
        .unwrap();
    // wait for the first step rows of the run record
    let rec = split.join("run.jsonl");
    let t0 = std::time::Instant::now();
    loop {
        let n = std::fs::read_to_string(&rec).map_or(0, |t| t.matches("\"step\"").count());
        if n > 0 {
            break;
        }
        assert!(t0.elapsed().as_secs() < 120, "no step rows in run.jsonl");
        std::thread::sleep(std::time::Duration::from_millis(20));
    }
    let k = Command::new("kill")
        .args(["-TERM", &child.id().to_string()])
        .status()
        .unwrap();
    assert!(k.success());
    let st = child.wait().unwrap();
    assert_eq!(st.code(), Some(3));
    let v: Value = serde_json::from_slice(&bytes(&split.join("result.json"))).unwrap();
    assert_eq!(v["stop"]["code"], "RUN002");
    assert_eq!(v["stop"]["signal"], "SIGTERM");

    let o = foehn(&[
        "--quiet",
        "flow",
        s(&case),
        "--restart",
        "latest",
        "--out",
        s(&split),
    ]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    for f in ["fields.vtkhdf", "restart.h5"] {
        assert!(
            bytes(&straight.join(f)) == bytes(&split.join(f)),
            "{f}: stop and restart differs from the straight run"
        );
    }
}

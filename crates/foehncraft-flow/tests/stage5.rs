//! Stage 5 acceptance: transient flow with PISO (docs/verification/stage5.md).
//!
//! The full Taylor–Green sweeps and the cylinder benchmark also run from the
//! case files in CI; the tests here use small meshes.
#![cfg(feature = "hdf5")]

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::restart;
use foehncraft_fv::disc::Soa3;
use std::path::{Path, PathBuf};

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_stage5_{name}_{}", std::process::id()))
}

fn run_in(case: &FlowCase, name: &str, dir: &Path) -> FlowReport {
    run(case, name, dir, |_| {}).unwrap()
}

/// The lid-driven cavity on an `n` x `n` hex mesh with a `[time]` table.
fn cavity(n: usize, time: &str, extra: &str) -> FlowCase {
    let text = format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [{n}, {n}, 1]
[fluid]
viscosity = 0.01
[boundary.xmin]
type = "wall"
[boundary.xmax]
type = "wall"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
velocity = [1.0, 0.0, 0.0]
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[output]
vtkhdf = false
wall_csv = false
[time]
{time}
{extra}
"#
    );
    FlowCase::from_toml(&text, Path::new(".")).unwrap()
}

/// The Taylor–Green vortex on the unit square with the exact velocity on
/// the four sides.
fn taylor_green(n: usize, time: &str, extra: &str) -> FlowCase {
    let text = format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [{n}, {n}, 1]
[fluid]
viscosity = 0.01
[flow]
convection = true
[boundary.xmin]
type = "velocity-inlet"
profile = "analytic"
[boundary.xmax]
type = "velocity-inlet"
profile = "analytic"
[boundary.ymin]
type = "velocity-inlet"
profile = "analytic"
[boundary.ymax]
type = "velocity-inlet"
profile = "analytic"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[analytic]
type = "taylor-green"
[output]
vtkhdf = false
wall_csv = false
[time]
{time}
{extra}
"#
    );
    FlowCase::from_toml(&text, Path::new(".")).unwrap()
}

fn same3(a: &Soa3, b: &Soa3) -> bool {
    a.x == b.x && a.y == b.y && a.z == b.z
}

fn same_restart(a: &restart::Restart, b: &restart::Restart) -> bool {
    let (s, t) = (&a.state, &b.state);
    let fields = same3(&s.u, &t.u)
        && s.p == t.p
        && s.mass_flux == t.mass_flux
        && same3(&s.ub, &t.ub)
        && s.pb == t.pb
        && (0..3).all(|i| same3(&s.grad_u[i], &t.grad_u[i]))
        && same3(&s.grad_p, &t.grad_p)
        && same3(&s.psi, &t.psi);
    let levels = match (&a.transient, &b.transient) {
        (Some((o, c, dt, _)), Some((p, d, du, _))) => {
            o.u_old == p.u_old
                && o.m_old == p.m_old
                && o.x_cur == p.x_cur
                && o.x_old == p.x_old
                && o.g_old == p.g_old
                && o.u_old2 == p.u_old2
                && c == d
                && dt == du
        }
        _ => false,
    };
    fields && levels && a.iteration == b.iteration && a.amg_matrix == b.amg_matrix
}

/// A run continued from a restart file ends bit-identically to the run
/// without the interruption (D42): BDF2 with the adaptive time step, so the
/// clock, the old levels, the extrapolation and the multigrid all carry
/// over.
#[test]
fn restart_continues_bit_identically() {
    let time = "scheme = \"bdf2\"\ntime_step = 0.01\nsteps = 24\nadaptive = true\nmax_cfl = 2.0";
    let full = cavity(
        16,
        time,
        "[run]\nrestart_at = [10]\n[solver]\namg_rebuild = 5",
    );
    let da = out("restart_full");
    let ra = run_in(&full, "full", &da);
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    let from = da.join("restart_000010.h5");
    let cont = cavity(
        16,
        time,
        &format!(
            "[run]\nrestart = '{}'\n[solver]\namg_rebuild = 5",
            from.display()
        ),
    );
    let db = out("restart_cont");
    let rb = run_in(&cont, "cont", &db);
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    let b = restart::read(&db.join("restart.h5")).unwrap();
    assert_eq!(
        a.transient.as_ref().unwrap().1.step,
        24,
        "the full run ends at step 24"
    );
    assert!(same_restart(&a, &b), "the continued run differs");
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// The time step from the error estimate (D65) with rejected steps: the
/// run steps back exactly and continues from a restart file bit-identically
/// (the level uⁿ⁻² of the predictor and the step before the last are part
/// of the restart).
#[test]
fn error_controlled_steps_restart_bit_identically() {
    let time = "scheme = \"bdf2\"\ntime_step = 0.02\nsteps = 30\nadaptive = true\nmax_cfl = 2.0\n\
                error_tolerance = 1e-3";
    let full = cavity(16, time, "[run]\nrestart_at = [12]");
    let da = out("restart_err_full");
    let ra = run_in(&full, "full", &da);
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    let log = std::fs::read_to_string(da.join("flow.log")).unwrap();
    assert!(
        log.contains(" rejected: error "),
        "no step was rejected:\n{log}"
    );
    assert!(
        log.contains("  error 0."),
        "no error estimate in the log:\n{log}"
    );
    let from = da.join("restart_000012.h5");
    let r = restart::read(&from).unwrap();
    let (o, c, ..) = r.transient.as_ref().unwrap();
    assert!(o.u_old2.is_some() && c.dt_prev2 > 0.0);
    let cont = cavity(16, time, &format!("[run]\nrestart = '{}'", from.display()));
    let db = out("restart_err_cont");
    let rb = run_in(&cont, "cont", &db);
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    let b = restart::read(&db.join("restart.h5")).unwrap();
    assert!(same_restart(&a, &b), "the continued run differs");
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// Crank–Nicolson with two outer iterations continues bit-identically too
/// (the explicit operator of the old level and the outer loop).
#[test]
fn restart_crank_nicolson_outer() {
    let time = "scheme = \"crank-nicolson\"\ntime_step = 0.01\nsteps = 12";
    let extra = "[solver]\nouter_iterations = 2";
    let full = taylor_green(16, time, &format!("[run]\nrestart_at = [5]\n{extra}"));
    let da = out("restart_cn_full");
    run_in(&full, "full", &da);
    let from = da.join("restart_000005.h5");
    let cont = taylor_green(
        16,
        time,
        &format!("[run]\nrestart = '{}'\n{extra}", from.display()),
    );
    let db = out("restart_cn_cont");
    run_in(&cont, "cont", &db);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    let b = restart::read(&db.join("restart.h5")).unwrap();
    assert!(same_restart(&a, &b), "the continued run differs");
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// `.stop` ends the run after the current step with the restart file and
/// the results (D40).
#[test]
fn stop_file_ends_the_run() {
    let case = cavity(
        12,
        "time_step = 0.01\nsteps = 50",
        "[solver]\nlog_every = 1",
    );
    let dir = out("stop");
    std::fs::create_dir_all(&dir).unwrap();
    let stop = dir.join(".stop");
    let rep = run(&case, "stop", &dir, |l: &str| {
        if l.starts_with("step") && l.split_whitespace().nth(1) == Some("5") {
            std::fs::write(&stop, "").unwrap();
        }
    })
    .unwrap();
    assert!(rep.stopped);
    let r = restart::read(&dir.join("restart.h5")).unwrap();
    assert_eq!(r.transient.unwrap().1.step, 5);
    // a start with an existing .stop fails at once
    let e = run(&case, "stop", &dir, |_| {}).unwrap_err();
    assert!(format!("{e:#}").contains(".stop"), "{e:#}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// `.result_N` writes the fields of step N; `.read_toml_N` re-reads the
/// case file at step N and applies only what may change during a run.
#[test]
fn result_and_reread_files() {
    let dir = out("reread");
    std::fs::create_dir_all(&dir).unwrap();
    let text = |correctors: usize, viscosity: f64, steps: usize| {
        format!(
            r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [12, 12, 1]
[fluid]
viscosity = {viscosity}
[boundary.xmin]
type = "wall"
[boundary.xmax]
type = "wall"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
velocity = [1.0, 0.0, 0.0]
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[output]
wall_csv = false
[solver]
log_every = 1
correctors = {correctors}
[time]
time_step = 0.01
steps = {steps}
"#
        )
    };
    let file = dir.join("case.toml");
    std::fs::write(&file, text(2, 0.01, 20)).unwrap();
    std::fs::write(dir.join(".result_3"), "").unwrap();
    std::fs::write(dir.join(".read_toml_4"), "").unwrap();
    let case = FlowCase::from_file(&file).unwrap();
    let rep = run(&case, "reread", &dir, |l: &str| {
        if l.starts_with("step") && l.split_whitespace().nth(1) == Some("2") {
            std::fs::write(&file, text(3, 0.02, 8)).unwrap();
        }
    })
    .unwrap();
    assert!(!rep.stopped);
    let log = std::fs::read_to_string(dir.join("flow.log")).unwrap();
    assert!(log.contains("changed solver.correctors: 2 -> 3"), "{log}");
    assert!(log.contains("changed time.steps: 20 -> 8"), "{log}");
    assert!(log.contains("not applied"), "{log}");
    assert!(log.contains("fluid.viscosity"), "{log}");
    assert!(dir.join("fields_000003.vtkhdf").exists());
    assert!(dir.join("case.read_toml_4.toml").exists());
    let r = restart::read(&dir.join("restart.h5")).unwrap();
    assert_eq!(r.transient.unwrap().1.step, 8);
    let _ = std::fs::remove_dir_all(&dir);
}

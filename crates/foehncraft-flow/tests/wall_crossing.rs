//! Stage 6, part 6g: wall crossings (`[[output.wall_crossings]]`), the
//! first point along a direction where the wall shear changes sign, and
//! their discretisation uncertainty in a level sweep (D67),
//! docs/verification/stage6.md. The validation case is
//! `backstep_validation.toml`; the tests here use a laminar step.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_crossing_{name}_{}", std::process::id()))
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

/// A laminar step: a channel of height 1 from x = −2 opens at x = 0 to
/// height 2 (step height 1), parabolic inflow of mean speed 2/3, ν = 0.01
/// (Re = 133 on the mean speed and the inlet height): one bubble behind
/// the step. `crossings` are the `[[output.wall_crossings]]` entries.
fn laminar_step(crossings: &str, extra: &str) -> String {
    format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
lower = [-2.0, 0.0, 0.0]
upper = [12.0, 2.0, 0.1]
cells = [56, 16, 1]
remove = [ {{ lower = [-2.0, 0.0, 0.0], upper = [0.0, 1.0, 0.1], name = "step" }} ]
[fluid]
viscosity = 0.01
[boundary.xmin]
type = "velocity-inlet"
profile = "parabolic"
max_velocity = 1.0
profile_axes = ["y"]
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.step]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver]
algorithm = "coupled"
tolerance = 1e-8
max_iterations = 200
[output]
vtkhdf = false
wall_csv = true
{crossings}
{extra}
"#
    )
}

/// Runs a case; returns the report, the log and the (x, τ_x) rows of
/// `wall_ymin.csv`, ordered by x.
fn run_step(text: &str, name: &str) -> (FlowReport, String, Vec<(f64, f64)>) {
    let dir = out(name);
    let mut log = String::new();
    let rep = run(&parse(text), name, &dir, |l| {
        log.push_str(l);
        log.push('\n');
    })
    .unwrap();
    let csv = std::fs::read_to_string(dir.join("wall_ymin.csv")).unwrap_or_default();
    let _ = std::fs::remove_dir_all(&dir);
    let mut rows: Vec<&str> = csv.lines().filter(|l| !l.starts_with('#')).collect();
    let mut wall = Vec::new();
    if !rows.is_empty() {
        let head: Vec<&str> = rows.remove(0).split(',').collect();
        let ix = head.iter().position(|h| *h == "x").unwrap();
        let it = head.iter().position(|h| *h == "tau_x").unwrap();
        for r in rows {
            let v: Vec<f64> = r.split(',').map(|s| s.parse().unwrap()).collect();
            wall.push((v[ix], v[it]));
        }
    }
    wall.sort_by(|a, b| a.0.total_cmp(&b.0));
    (rep, log, wall)
}

/// The reattachment point is the sign change of the wall shear that the
/// wall file shows, linear between the two face centroids; searched
/// upstream from the outlet (direction −x) it is the same point, at the
/// distance from the start; a kind that the shear does not show gives
/// no position.
#[test]
fn reattachment_is_the_sign_change_of_the_wall_shear() {
    let crossings = r#"wall_crossings = [
  { name = "xr", boundary = "ymin", from = [0.0, 0.0, 0.05] },
  { name = "back", boundary = "ymin", from = [12.0, 0.0, 0.05], direction = [-1.0, 0.0, 0.0] },
  { name = "none", boundary = "ymin", kind = "separation", from = [6.0, 0.0, 0.05] },
  { name = "row", boundary = "ymin", from = [0.0, 0.0, 0.05], width = 0.1 },
]"#;
    let (rep, log, wall) = run_step(&laminar_step(crossings, ""), "laminar");
    let s = rep.solves.last().unwrap();
    let by = |n: &str| s.wall_crossings.iter().find(|c| c.name == n).unwrap();
    let (x0, t0, x1, t1) = wall
        .windows(2)
        .find(|w| w[0].1 < 0.0 && w[1].1 >= 0.0)
        .map(|w| (w[0].0, w[0].1, w[1].0, w[1].1))
        .expect("a bubble behind the step");
    let want = x0 + (x1 - x0) * (-t0) / (t1 - t0);
    let xr = by("xr");
    assert_eq!(xr.kind, "reattachment");
    let got = xr.position.unwrap();
    assert!(
        (got - want).abs() < 1e-12,
        "x_r {got}, from the wall file {want}"
    );
    assert!((1.0..8.0).contains(&got), "x_r/H {got}");
    let at = xr.at.unwrap();
    assert!(
        (at[0] - want).abs() < 1e-12 && at[1].abs() < 1e-12,
        "at {at:?}"
    );
    let back = by("back").position.unwrap();
    assert!(
        (12.0 - back - want).abs() < 1e-12,
        "searched upstream: {back}"
    );
    assert_eq!(by("none").position, None);
    assert_eq!(by("none").at, None);
    assert!((by("row").position.unwrap() - want).abs() < 1e-12);
    // in the report data and in the log
    let js = serde_json::to_value(s).unwrap();
    assert_eq!(js["wall_crossings"][0]["name"], "xr");
    assert!(
        log.contains("wall crossing xr (reattachment) on ymin"),
        "{log}"
    );
    assert!(log.contains("wall crossing none (separation) on ymin: the shear does not"));
}

/// Three meshes: the reattachment length is a monitored quantity of the
/// sweep, with its observed order and GCI.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn level_sweep_gives_the_uncertainty_of_a_crossing() {
    let crossings =
        r#"wall_crossings = [ { name = "xr", boundary = "ymin", from = [0.0, 0.0, 0.05] } ]"#;
    let text = laminar_step(crossings, "[verification]\nlevels = [0, 1, 2]\n")
        .replace("max_iterations = 200", "max_iterations = 400");
    let (rep, _, _) = run_step(&text, "sweep");
    assert_eq!(rep.uncertainty_of, "meshes");
    let e = rep
        .uncertainty
        .iter()
        .find(|e| e.quantity == "wall_crossings.xr.position")
        .expect("the reattachment length is monitored");
    let xr: Vec<f64> = rep
        .solves
        .iter()
        .map(|s| s.wall_crossings[0].position.unwrap())
        .collect();
    assert_eq!(e.values, [xr[2], xr[1], xr[0]], "finest first");
    assert!(e.gci.is_some_and(|g| g < 0.05), "GCI {:?}", e.gci);
    assert!(rep.to_markdown().contains("wall_crossings.xr.position"));
}

/// A zero direction, a width that is not positive and an unknown kind are
/// refused before the run.
#[test]
fn bad_crossing_input_is_refused() {
    for (entry, code) in [
        (
            r#"{ name = "a", boundary = "ymin", from = [0.0, 0.0, 0.0], direction = [0.0, 0.0, 0.0] }"#,
            "CASE004",
        ),
        (
            r#"{ name = "a", boundary = "ymin", from = [0.0, 0.0, 0.0], width = 0.0 }"#,
            "CASE004",
        ),
        (
            r#"{ name = "a", boundary = "ymin", from = [0.0, 0.0, 0.0], kind = "transition" }"#,
            "CASE",
        ),
    ] {
        let text = laminar_step(&format!("wall_crossings = [ {entry} ]"), "");
        let e = FlowCase::from_toml(&text, &cases()).expect_err(entry);
        let got = diag::of_error(&e).code.to_string();
        assert!(got.starts_with(code), "{entry}: {got} ({e:#})");
    }
}

//! Stage 6, part 6h: the curvature and rotation correction of SST
//! (Smirnov & Menter 2009, D170, specification 3.6 and 9.8) and the
//! rotating frame of the whole domain (`flow.frame_rotation`).

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_cc_{name}_{}", std::process::id()))
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

fn refused(text: &str) -> (String, String) {
    let e = match FlowCase::from_toml(text, &cases()) {
        Err(e) => e,
        Ok(c) => match run(&c, "refused", &out("refused"), |_| {}) {
            Err(e) => e,
            Ok(_) => panic!("the case was accepted"),
        },
    };
    (diag::of_error(&e).code.to_string(), format!("{e:#}"))
}

/// The rotating channel of the gallery, small: 24 cells across.
fn small_rotating(iterations: usize, correction: bool) -> String {
    text("channel_rotating_sst.toml")
        .replace("cells = [2, 96, 1]", "cells = [2, 24, 1]")
        .replace("ratio = 30.0", "ratio = 6.0")
        .replace(
            "max_iterations = 6000",
            &format!("max_iterations = {iterations}"),
        )
        .replace(
            "curvature_correction = true",
            &format!("curvature_correction = {correction}"),
        )
        + "\n[verification]\nrequire_convergence = false\n"
}

/// The wall shear forces (x) on the lower and the upper wall.
fn wall_forces(r: &FlowReport) -> (f64, f64) {
    let s = r.solves.last().unwrap();
    let f = |n: &str| {
        s.boundaries
            .iter()
            .find(|b| b.name == n)
            .map(|b| b.force[0])
            .unwrap()
    };
    (f("ymin"), f("ymax"))
}

/// Without rotation the channel is pure shear (r* = 1, r̃ = 0): the
/// correction leaves the result unchanged (band of 9.8: 1e-6).
#[test]
fn pure_shear_is_unchanged() {
    let base = small_rotating(60, false).replace("frame_rotation = [0.0, 0.0, 1.121]", "");
    let on = base.replace(
        "curvature_correction = false",
        "curvature_correction = true",
    );
    let a = run_case(&parse(&base), "shear_off");
    let b = run_case(&parse(&on), "shear_on");
    let (fa, fb) = (
        a.solves[0].periodic[0].driving_force,
        b.solves[0].periodic[0].driving_force,
    );
    assert!((fa - fb).abs() <= 1e-6 * fa.abs(), "{fa} {fb}");
}

/// Spanwise rotation: linear SST keeps the two walls equal, the
/// correction raises the wall shear stress on the unstable (lower) side
/// and lowers it on the stable side.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn spanwise_rotation_splits_the_walls() {
    let (lo, hi) = wall_forces(&run_case(&parse(&small_rotating(800, false)), "rot_off"));
    assert!((lo - hi).abs() <= 1e-6 * lo.abs(), "{lo} {hi}");
    let (lo, hi) = wall_forces(&run_case(&parse(&small_rotating(800, true)), "rot_on"));
    assert!(lo > 1.2 * hi && hi > 0.0, "{lo} {hi}");
}

/// The option belongs to SST; the rotating frame of the whole domain does
/// not combine with a moving reference frame.
#[test]
fn bad_curvature_input_is_refused() {
    let t = small_rotating(10, true);
    let (code, msg) = refused(
        &t.replace("model = \"k-omega-sst\"", "model = \"k-epsilon\"")
            .replace(
                "wall_treatment = \"low-re\"",
                "wall_treatment = \"scalable\"",
            ),
    );
    assert_eq!(code, "CASE005", "{msg}");
    assert!(msg.contains("curvature_correction"), "{msg}");
    let (code, msg) = refused(&format!(
        "{t}\n[[mrf]]\naxis = [0.0, 0.0, 1.0]\nomega = 1.0\n"
    ));
    assert!(msg.contains("frame_rotation"), "{code}: {msg}");
}

/// The square duct of the gallery, small: 12 × 12 cells on the quarter.
fn small_duct(iterations: usize, qcr: &str) -> String {
    text("square_duct_qcr.toml")
        .replace("cells = [2, 72, 72]", "cells = [2, 12, 12]")
        .replace(
            "ratio = 12.0 }, z = { ratio = 12.0",
            "ratio = 3.0 }, z = { ratio = 3.0",
        )
        .replace(
            "max_iterations = 3000",
            &format!("max_iterations = {iterations}"),
        )
        .replace("qcr = \"qcr2000\"", &format!("qcr = \"{qcr}\""))
        .replace("vtkhdf = true", "vtkhdf = false")
        + "\n[verification]\nrequire_convergence = false\n"
}

/// The largest secondary velocity (√(v² + w²)) on the bisector line.
fn secondary(case: &str, name: &str) -> f64 {
    let dir = out(name);
    run(&parse(case), name, &dir, |_| {}).unwrap();
    let t = std::fs::read_to_string(dir.join("line_bisector.csv")).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let mut rows = t.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = rows.next().unwrap().split(',').collect();
    let col = |n: &str| head.iter().position(|h| *h == n).unwrap();
    let (iv, iw) = (col("v"), col("w"));
    rows.map(|l| {
        let v: Vec<f64> = l.split(',').map(|x| x.parse().unwrap()).collect();
        v[iv].hypot(v[iw])
    })
    .fold(0.0, f64::max)
}

/// Linear SST gives no secondary flow in a square duct; the QCR (2000 and
/// 2013) drives the corner vortices (specification 9.8).
#[test]
fn qcr_drives_the_secondary_flow_of_a_duct() {
    let off = secondary(&small_duct(150, "none"), "duct_off");
    assert!(off < 1e-10, "{off}");
    for (mode, name) in [("qcr2000", "duct_2000"), ("qcr2013", "duct_2013")] {
        let on = secondary(&small_duct(150, mode), name);
        assert!(on > 1e-4, "{mode}: {on}");
    }
}

/// The QCR belongs to SST.
#[test]
fn qcr_is_refused_without_sst() {
    let (code, msg) = refused(
        &small_duct(10, "qcr2000")
            .replace("model = \"k-omega-sst\"", "model = \"k-epsilon\"")
            .replace(
                "wall_treatment = \"low-re\"",
                "wall_treatment = \"scalable\"",
            ),
    );
    assert_eq!(code, "CASE005", "{msg}");
    assert!(msg.contains("turbulence.qcr"), "{msg}");
}

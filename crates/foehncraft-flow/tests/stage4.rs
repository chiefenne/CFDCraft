//! Stage 4 acceptance: laminar steady flow with SIMPLEC (docs/verification/stage4.md).
//!
//! The full sweeps and the cavity cases also run from the case files in CI;
//! the heavier tests here run in release builds only.

use foehncraft_flow::case::{FlowCase, FlowReport, check_case, run};
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_stage4_{name}_{}", std::process::id()))
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

fn load(file: &str) -> FlowCase {
    FlowCase::from_file(&cases().join(file)).unwrap()
}

/// A small case: the 2D unit square with one cell layer and four boundaries.
fn square(n: usize, kind: &str, extra: &str) -> FlowCase {
    let text = format!(
        r#"
[mesh.generate]
type = "{kind}"
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
{extra}
"#
    );
    FlowCase::from_toml(&text, Path::new(".")).unwrap()
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn couette_is_exact_on_all_mesh_types() {
    for t in ["hex", "tet", "poly", "mixed"] {
        let rep = run_case(&load(&format!("couette_{t}.toml")), &format!("couette_{t}"));
        assert!(rep.failures.is_empty(), "{t}: {:?}", rep.failures);
        let e = rep.solves[0].errors.as_ref().unwrap();
        assert!(e.l2 < 1e-9, "{t}: L2 error {:e}", e.l2);
    }
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn poiseuille_is_second_order_on_hexahedra() {
    let mut case = load("poiseuille_hex.toml");
    case.verification.levels = vec![0, 1];
    let rep = run_case(&case, "poiseuille_hex");
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let o = rep.orders[0];
    assert!((1.8..=2.2).contains(&o), "observed order {o}");
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn cavity_re100_agrees_with_ghia() {
    let mut case = load("cavity_re100.toml");
    case.mesh.generate.as_mut().unwrap().cells = [32, 32, 1];
    case.mesh.generate.as_mut().unwrap().upper[2] = 1.0 / 32.0;
    case.output.lines.clear();
    case.output.vtkhdf = false;
    case.verification.reference_tolerance = Some(0.04);
    let rep = run_case(&case, "cavity_re100");
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let r = rep.solves[0].reference.as_ref().unwrap();
    assert!(r.u_max_dev < 0.04 && r.v_max_dev < 0.04, "{r:?}");
}

/// Symmetry planes are implicit per component: a regression test for the
/// lagged treatment that acted as a strong extra under-relaxation (about
/// 1700 outer iterations for this case instead of under 100).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn stokes_cavity_converges_quickly() {
    let case = square(
        16,
        "hex",
        "[flow]\nconvection = false\n[solver]\ntolerance = 1e-8\n[solver.momentum]\nrelaxation = 0.8\n",
    );
    let rep = run_case(&case, "stokes_cavity");
    let s = &rep.solves[0];
    assert!(
        s.converged && s.iterations < 200,
        "{} iterations",
        s.iterations
    );
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn results_are_bitwise_equal_on_1_and_4_threads() {
    // above the 20 000 cells below which one thread is used
    let mk = |threads: usize| {
        square(
            150,
            "hex",
            &format!(
                "[solver]\nmax_iterations = 15\n[run]\nthreads = {threads}\n[verification]\nrequire_convergence = false\n"
            ),
        )
    };
    let a = run_case(&mk(1), "threads1");
    let b = run_case(&mk(4), "threads4");
    let (ra, rb) = (&a.solves[0].residuals, &b.solves[0].residuals);
    assert_eq!(ra.momentum.map(f64::to_bits), rb.momentum.map(f64::to_bits));
    assert_eq!(ra.continuity.to_bits(), rb.continuity.to_bits());
    let (pa, pb) = (&a.solves[0].boundaries, &b.solves[0].boundaries);
    for (x, y) in pa.iter().zip(pb) {
        assert_eq!(
            x.force.map(f64::to_bits),
            y.force.map(f64::to_bits),
            "{}",
            x.name
        );
    }
}

#[test]
fn unknown_keys_and_values_get_a_hint() {
    let e = FlowCase::from_toml("[fluid]\nviscosty = 1.0\n", Path::new(".")).unwrap_err();
    assert!(e.to_string().contains("did you mean `viscosity`"), "{e}");
    let e = FlowCase::from_toml(
        "[boundary.xmin]\ntype = \"velocity-inlett\"\nvelocity = [1.0, 0.0, 0.0]\n",
        Path::new("."),
    )
    .unwrap_err();
    assert!(
        e.to_string().contains("did you mean `velocity-inlet`"),
        "{e}"
    );
}

#[test]
fn boundary_check_finds_unknown_and_missing_boundaries() {
    let mut case = square(4, "hex", "");
    let lid = case.boundary.remove("ymax").unwrap();
    case.boundary.insert("ymx".into(), lid);
    let mut lines = Vec::new();
    let r = check_case(&case, Path::new("case.toml"), |l| lines.push(l.to_string()));
    assert!(r.is_err());
    let all = lines.join("\n");
    assert!(all.contains("did you mean `ymax`"), "{all}");
}

#[test]
fn boundary_check_rejects_an_unbalanced_closed_domain() {
    let mut case = square(4, "hex", "");
    case.boundary.insert(
        "xmin".into(),
        toml::from_str("type = \"velocity-inlet\"\nvelocity = [1.0, 0.0, 0.0]").unwrap(),
    );
    let mut lines = Vec::new();
    assert!(check_case(&case, Path::new("case.toml"), |l| lines.push(l.to_string())).is_err());
    assert!(
        lines.iter().any(|l| l.contains("no pressure boundary")),
        "{lines:?}"
    );
}

#[test]
fn boundary_check_accepts_the_cavity_with_a_note() {
    let case = square(4, "hex", "");
    let mut lines = Vec::new();
    check_case(&case, Path::new("case.toml"), |l| lines.push(l.to_string())).unwrap();
    assert!(
        lines.iter().any(|l| l.contains("reference cell")),
        "{lines:?}"
    );
}

/// D27: the converged result does not depend on the under-relaxation (the
/// Rhie–Chow flux carries the relaxation term of Majumdar 1988).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn converged_result_does_not_depend_on_relaxation() {
    let lid_force = |alpha: f64| {
        let case = square(
            16,
            "hex",
            &format!(
                "[solver]\ntolerance = 1e-11\nmax_iterations = 5000\n[solver.momentum]\nrelaxation = {alpha}\n"
            ),
        );
        let rep = run_case(&case, &format!("relax{alpha}"));
        assert!(rep.solves[0].converged);
        rep.solves[0]
            .boundaries
            .iter()
            .find(|p| p.name == "ymax")
            .unwrap()
            .force
    };
    let (a, b) = (lid_force(0.6), lid_force(0.9));
    for i in 0..2 {
        assert!((a[i] - b[i]).abs() <= 1e-8 * a[0].abs(), "{a:?} {b:?}");
    }
}

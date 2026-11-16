//! Stage 5, part 6c: the discretisation uncertainty of sweeps (D67),
//! docs/verification/stage5.md.

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::diag;
use foehncraft_flow::uncertainty::Verdict;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

/// A case file with `\n` line ends (a Windows checkout has `\r\n`).
fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

/// Couette flow with μ = μ0 (1 + a y) on three meshes: the wall shear
/// force converges at second order, and its extrapolated value lies
/// closer to the exact μ0 U a / ln(1 + a) times the wall area than the
/// finest solution, inside the band of the GCI.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn wall_shear_of_a_level_sweep() {
    let case = FlowCase::from_toml(&text("couette_viscosity.toml"), &cases()).unwrap();
    let dir = std::env::temp_dir().join(format!("foehncraft_gci_{}", std::process::id()));
    let rep = run(&case, "couette_viscosity", &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    assert_eq!(rep.uncertainty_of, "meshes");
    let e = rep
        .uncertainty
        .iter()
        .find(|e| e.quantity == "boundaries.ymin.force_x")
        .expect("the lower-wall force is monitored");
    assert_eq!(e.verdict, Verdict::Monotonic);
    let p = e.order.unwrap();
    assert!((1.5..2.5).contains(&p), "order {p}");
    let (mu0, u, a, area) = (0.1, 1.0, 4.0, 2.0 * 0.1);
    let exact = mu0 * u * a / f64::ln_1p(a) * area;
    let ext = e.extrapolated.unwrap();
    let band = e.band.unwrap();
    assert!(
        (ext - exact).abs() < 0.1 * (e.values[0] - exact).abs(),
        "extrapolated {ext}, finest {}, exact {exact}",
        e.values[0]
    );
    assert!(
        (e.values[0] - exact).abs() <= band,
        "exact {exact} outside φ1 ± {band}"
    );
    // the error norms have their own orders and are not repeated
    assert!(
        rep.uncertainty
            .iter()
            .all(|e| !e.quantity.starts_with("errors"))
    );
    assert!(rep.to_markdown().contains("Discretisation uncertainty"));
}

/// A sweep needs no exact solution (its quantities get their uncertainty),
/// but an expected order does.
#[test]
fn sweeps_without_an_exact_solution() {
    let t = text("cavity_re100.toml")
        .replace("[verification]\n", "[verification]\nlevels = [0, 1, 2]\n");
    assert!(FlowCase::from_toml(&t, &cases()).is_ok());
    let bad = t.replace(
        "levels = [0, 1, 2]\n",
        "levels = [0, 1, 2]\nexpected_order = [1.8, 2.2]\n",
    );
    let e = FlowCase::from_toml(&bad, &cases()).unwrap_err();
    assert_eq!(diag::of_error(&e).code, "CASE005", "{e:#}");
    assert!(format!("{e:#}").contains("expected_order needs an [analytic] solution"));
}

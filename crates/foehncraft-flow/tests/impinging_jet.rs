//! Stage 6 part 6g (specification 9.6): the round jet from a pipe that
//! impinges on a heated plate at H/D = 2 (`impinging_jet.toml`), on the
//! coarsest mesh of the study (2 600 cells). The validation against the
//! experiment is `examples/studies/impinging_jet.py`; here the shape of
//! Nu(r/D) and the thermal wall law at the stagnation point.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use std::path::{Path, PathBuf};

/// Γ = ν/Pr of the case (ν = 1/23 750 m²/s, Pr = 0.71).
const GAMMA: f64 = 1.0 / 23750.0 / 0.71;

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

/// The jet with `model` and its wall treatment on 2 600 cells: the report
/// and Nu(r/D) = q D/(Γ θ_w) along the plate.
fn jet(model: &str, treatment: &str) -> (FlowReport, Vec<(f64, f64)>) {
    let mut t = std::fs::read_to_string(cases().join("impinging_jet.toml"))
        .unwrap()
        .replace("\r\n", "\n");
    for (a, b) in [
        ("refine = 1 ", "refine = 0 ".to_string()),
        (
            "model = \"k-epsilon-zeta-f0\"",
            format!("model = \"{model}\""),
        ),
        (
            "wall_treatment = \"compound-2023\"",
            format!("wall_treatment = \"{treatment}\""),
        ),
    ] {
        assert!(t.contains(a), "{a}");
        t = t.replace(a, &b);
    }
    let case = FlowCase::from_toml(&t, &cases()).unwrap();
    let dir = std::env::temp_dir().join(format!(
        "foehncraft_impinging_jet_{model}_{}",
        std::process::id()
    ));
    let r = run(&case, "jet", &dir, |_| {}).unwrap();
    let csv = std::fs::read_to_string(dir.join("wall_xmin.csv")).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let mut lines = csv.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split(',').collect();
    let col = |n: &str| head.iter().position(|h| *h == n).unwrap();
    let (cr, cw) = (col("y"), col("theta_wall"));
    let mut nu: Vec<(f64, f64)> = lines
        .map(|l| {
            let v: Vec<f64> = l.split(',').map(|s| s.parse().unwrap()).collect();
            (v[cr], 1e-3 / (GAMMA * v[cw]))
        })
        .collect();
    nu.sort_by(|a, b| a.0.total_cmp(&b.0));
    (r, nu)
}

fn at(nu: &[(f64, f64)], r: f64) -> f64 {
    let i = nu.iter().position(|p| p.0 >= r).unwrap();
    let ((r0, a), (r1, b)) = (nu[i - 1], nu[i]);
    a + (b - a) * (r - r0) / (r1 - r0)
}

/// k-ε-ζ-f₀: converged, the wall probe `stag` at the axis gives the Nu₀ of
/// the face next to it and `r01` the Nu at r/D = 0.1, and Nu(r/D) falls to
/// a minimum and rises again to a secondary maximum near r/D = 2 (128 at
/// r/D 2.0 on this mesh).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn the_jet_has_its_secondary_maximum() {
    let (r, nu) = jet("k-epsilon-zeta-f0", "compound-2023");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let probes = &r.solves[0].wall_probes;
    assert_eq!(probes[0].name, "stag");
    assert_eq!(probes[1].name, "r01");
    let nu0 = 1e-3 / (GAMMA * probes[0].scalars[0].wall);
    assert!((nu0 / nu[0].1 - 1.0).abs() < 1e-9, "{nu0} {}", nu[0].1);
    // the probe interpolates θ_w, the wall table Nu: equal to 0.1 %
    let (nu01, near) = (1e-3 / (GAMMA * probes[1].scalars[0].wall), at(&nu, 0.1));
    assert!((nu01 / near - 1.0).abs() < 1e-3, "{nu01} {near}");
    assert!((150.0..185.0).contains(&nu0), "Nu₀ {nu0}");
    let (rmin, nmin) = nu
        .iter()
        .filter(|p| (0.8..1.8).contains(&p.0))
        .fold((0.0, f64::INFINITY), |m, p| if p.1 < m.1 { *p } else { m });
    let (rmax, nmax) = nu
        .iter()
        .filter(|p| p.0 > rmin && p.0 < 3.0)
        .fold((0.0, 0.0), |m, p| if p.1 > m.1 { *p } else { m });
    eprintln!("Nu₀ {nu0:.1}, minimum {nmin:.1} at {rmin:.2}, maximum {nmax:.1} at {rmax:.2}");
    assert!(
        (1.5..2.5).contains(&rmax),
        "secondary maximum at r/D {rmax}"
    );
    assert!(nmax > 1.03 * nmin, "{nmax} against the minimum {nmin}");
}

/// The thermal law of the k-ε family takes u* = C_μ^¼ k_P^½ (section 5.2
/// of the stage-6 specification), which stays finite at the stagnation
/// point: Nu₀ stays near Nu(r/D = 0.5). With u_τ, which vanishes with the
/// wall shear, Nu₀ was a quarter of it.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn the_k_epsilon_thermal_law_keeps_the_stagnation_point() {
    let (r, nu) = jet("k-epsilon", "scalable");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let (nu0, half) = (nu[0].1, at(&nu, 0.5));
    eprintln!("k-ε: Nu₀ {nu0:.1}, Nu(0.5) {half:.1}");
    assert!(nu0 > 0.9 * half, "Nu₀ {nu0} against Nu(0.5) {half}");
}

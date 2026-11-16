//! Stage 5, part 6a: expressions (D53), materials and regions (D102, D112)
//! and the wall mass flux (D124 W3), docs/verification/stage5.md.
//!
//! The case files of the gallery also run in CI; the tests here use small
//! meshes or the same case files.

use foehncraft_flow::case::{FlowCase, FlowReport, check_case, run};
use foehncraft_flow::diag;
use foehncraft_flow::materials::{self, Conditions};
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_expr_{name}_{}", std::process::id()))
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

/// A case file with `\n` line ends (a Windows checkout has `\r\n`).
fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

/// The code of the error a case text is refused with, when it loads or
/// when its property program is resolved.
fn refused(text: &str) -> (String, String) {
    let e = match FlowCase::from_toml(text, &cases()) {
        Err(e) => e,
        Ok(c) => c.properties().expect_err("the case was accepted"),
    };
    (diag::of_error(&e).code.to_string(), format!("{e:#}"))
}

/// The lines of the boundary-condition check of a case.
fn check_lines(case: &FlowCase) -> (bool, String) {
    let mut lines = Vec::new();
    let r = check_case(case, Path::new("case.toml"), |l| lines.push(l.to_string()));
    if let Err(e) = &r {
        lines.push(format!("{}: {e:#}", diag::of_error(e).code));
    }
    (r.is_ok(), lines.join("\n"))
}

/// The Poiseuille case on its coarse mesh only.
fn poiseuille_coarse(t: &str) -> String {
    t.replace(
        "levels = [0, 1, 2]\nexpected_order = [1.8, 2.2]\n",
        "max_error = 2e-3\n",
    )
}

/// The exact parabola written as an expression gives the error of the
/// analytic profile, and half of the pressure gradient replaced by a body
/// force expression gives the same velocity.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn expressions_match_the_analytic_profile_and_the_pressure_gradient() {
    let hex = poiseuille_coarse(&text("poiseuille_hex.toml"));
    let exact = run_case(&parse(&hex), "exact");
    let inlet = hex.replace(
        "profile = \"analytic\"     # the exact parabola",
        "velocity = [\"6 * (y / 1.0) * (1 - y / 1.0)\", 0.0, 0.0]",
    );
    assert_ne!(inlet, hex);
    let expr = run_case(&parse(&inlet), "inlet");
    let body = run_case(
        &parse(&poiseuille_coarse(&text("poiseuille_expression.toml"))),
        "body",
    );
    for r in [&exact, &expr, &body] {
        assert!(r.failures.is_empty(), "{}: {:?}", r.name, r.failures);
    }
    let l2 = |r: &FlowReport| r.solves[0].errors.as_ref().unwrap().l2;
    let (a, b, c) = (l2(&exact), l2(&expr), l2(&body));
    assert!((a - b).abs() < 1e-9 * a, "analytic {a:e}, expression {b:e}");
    assert!((a - c).abs() < 1e-6 * a, "analytic {a:e}, body force {c:e}");
}

/// Errors in expressions, constants and tables carry their codes.
#[test]
fn bad_expressions_are_refused() {
    let t = text("poiseuille_expression.toml");
    // a syntax error
    let (code, msg) = refused(&t.replace("\"0.5 * gradient\"", "\"0.5 * (gradient\""));
    assert_eq!(code, "EXPR001", "{msg}");
    // an unknown name, with the nearest name as a hint
    let (code, msg) = refused(&t.replace("\"0.5 * gradient\"", "\"0.5 * gradiant\""));
    assert_eq!(code, "EXPR002", "{msg}");
    assert!(msg.contains("gradient"), "{msg}");
    // a reserved name as a constant
    let (code, msg) = refused(&t.replace("u_mean = 1.0", "u_mean = 1.0\npi = 3.0"));
    assert_eq!(code, "EXPR002", "{msg}");
    // time in a steady case
    let (code, msg) = refused(&t.replace("\"0.5 * gradient\"", "\"0.5 * gradient * t\""));
    assert_eq!(code, "EXPR003", "{msg}");
    // constants that refer to each other in a cycle
    let (code, msg) = refused(&t.replace("h = 1.0", "h = \"2 * gradient\""));
    assert_eq!(code, "EXPR002", "{msg}");
    assert!(msg.contains("cycle"), "{msg}");
    // a table that is not there
    let (code, msg) = refused(&format!("{t}\n[tables.lift]\nfile = \"no_such.csv\"\n"));
    assert_eq!(code, "EXPR004", "{msg}");
}

/// The materials database, case materials and the regions of D112.
#[test]
fn materials_and_regions() {
    // the database: water at 300 K (IAPWS-95 996.5 kg/m³, IAPWS-2008
    // 0.8538 mPa s) and air (ideal gas, Sutherland)
    let cond = Conditions {
        temperature: 300.0,
        pressure: 101325.0,
    };
    let mut warn = Vec::new();
    let env = foehncraft_flow::expr::Env::default();
    let num = |r: &materials::Resolved| r.value.constant().unwrap();
    let water = materials::resolve_material(
        Some("water"),
        &Default::default(),
        None,
        None,
        &env,
        cond,
        &mut warn,
    )
    .unwrap();
    let (rho, mu) = (num(&water.density), num(&water.viscosity));
    assert!((rho - 996.5).abs() < 0.2, "{rho}");
    assert!((mu - 0.8538e-3).abs() < 2e-6, "{mu}");
    let air = materials::resolve_material(
        Some("air"),
        &Default::default(),
        None,
        None,
        &env,
        cond,
        &mut warn,
    )
    .unwrap();
    let (rho, mu) = (num(&air.density), num(&air.viscosity));
    assert!((rho - 1.1766).abs() < 1e-3, "{rho}");
    assert!((mu - 1.846e-5).abs() < 5e-8, "{mu}");
    assert!(warn.is_empty(), "{warn:?}");
    // outside the range of the data: clamped, with a warning
    let hot = Conditions {
        temperature: 400.0,
        ..cond
    };
    materials::resolve_material(
        Some("water"),
        &Default::default(),
        None,
        None,
        &env,
        hot,
        &mut warn,
    )
    .unwrap();
    assert!(warn.iter().any(|w| w.code == "MAT003"), "{warn:?}");

    let t = text("poiseuille_water.toml");
    // the case resolves water at its temperature
    let c = parse(&t);
    let p = c.properties().unwrap();
    assert!(p.density() > 990.0 && p.density() < 1000.0);
    assert!(
        c.resolved_toml()
            .unwrap()
            .contains("resolved property program")
    );
    // an unknown material, with the nearest name as a hint
    let (code, msg) = refused(&t.replace("material = \"water\"", "material = \"watr\""));
    assert_eq!(code, "MAT001", "{msg}");
    assert!(msg.contains("water"), "{msg}");
    // a solid region is refused until the energy equation (L1)
    let (code, msg) = refused(&t.replace("kind = \"fluid\"", "kind = \"solid\""));
    assert_eq!(code, "REG001", "{msg}");
    // a zone the mesh does not have
    let bad = parse(&t.replace("zones = [\"fluid\"]", "zones = [\"fluids\"]"));
    let (ok, all) = check_lines(&bad);
    assert!(!ok && all.contains("REG002"), "{all}");
}

/// A wall that adds mass needs a boundary that takes it out (BC026);
/// suction balanced by blowing is accepted.
#[test]
fn wall_mass_flux_in_a_closed_volume() {
    let t = text("couette_transpiration.toml");
    let (ok, all) = check_lines(&parse(&t));
    assert!(ok, "{all}");
    // the outlet closed and the suction turned into blowing: the walls add
    // mass that cannot leave
    let closed = t
        .replace(
            "type = \"static-pressure-outlet\"\npressure = 0.0",
            "type = \"wall\"",
        )
        .replace("mass_flux = -0.2", "mass_flux = 0.2")
        .replace("profile = \"analytic\"", "velocity = [0.0, 0.0, 0.0]");
    assert_ne!(closed, t);
    let (ok, all) = check_lines(&parse(&closed));
    assert!(!ok && all.contains("BC026"), "{all}");
}

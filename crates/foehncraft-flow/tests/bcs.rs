//! Stage 5, part 5b: the boundary conditions of D29 (slip wall, wedge,
//! inlet profile from a table, backflow treatment at pressure outlets,
//! convective outlet), docs/verification/stage5.md.
//!
//! The case files of the gallery also run in CI; the tests here use small
//! meshes or the same case files.

use foehncraft_flow::case::{FlowCase, FlowReport, check_case, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_bcs_{name}_{}", std::process::id()))
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

/// A case file with `\n` line ends (a Windows checkout has `\r\n`), so
/// the edits of the tests match.
fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

/// The lines of the boundary-condition check of a case.
fn check_lines(case: &FlowCase) -> (bool, String) {
    let mut lines = Vec::new();
    let ok = check_case(case, Path::new("case.toml"), |l| lines.push(l.to_string())).is_ok();
    (ok, lines.join("\n"))
}

/// The uniform-flow channel between slip walls with `outlet` in place of
/// the static pressure outlet and `ymax` in place of the upper slip wall.
fn slip_channel(kind: &str, outlet: &str, ymax: &str, analytic: &str, extra: &str) -> FlowCase {
    let mut t = text("channel_slip.toml")
        .replace("type = \"mixed\"", &format!("type = \"{kind}\""))
        .replace(
            "[boundary.xmax]\ntype = \"static-pressure-outlet\"\npressure = 0.0",
            &format!("[boundary.xmax]\n{outlet}"),
        )
        .replace(
            "[boundary.ymax]\ntype = \"slip-wall\"",
            &format!("[boundary.ymax]\n{ymax}"),
        );
    if !analytic.is_empty() {
        let i = t.find("[analytic]").unwrap();
        let j = t.find("[verification]").unwrap();
        t.replace_range(i..j, &format!("[analytic]\n{analytic}\n\n"));
    }
    t.push_str(extra);
    parse(&t)
}

/// A slip wall carries no normal flow and no shear: the uniform flow is the
/// exact solution on every mesh type, held to round-off.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn slip_walls_hold_the_uniform_flow_on_every_mesh_type() {
    let outlet = "type = \"static-pressure-outlet\"\npressure = 0.0";
    for kind in ["hex", "tet", "poly", "mixed"] {
        let case = slip_channel(kind, outlet, "type = \"slip-wall\"", "", "");
        let rep = run_case(&case, &format!("slip_{kind}"));
        assert!(rep.failures.is_empty(), "{kind}: {:?}", rep.failures);
        let e = rep.solves[0].errors.as_ref().unwrap();
        assert!(e.l2 < 1e-9 && e.pressure_linf < 1e-9, "{kind}: {e:?}");
    }
}

/// Hagen–Poiseuille flow on the wedge: second order on three nested
/// meshes (the case fails outside 1.8–2.2).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn wedge_pipe_is_second_order() {
    let rep = run_case(&parse(&text("pipe_wedge.toml")), "pipe_wedge");
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let o = *rep.orders.last().unwrap();
    assert!((1.8..=2.2).contains(&o), "observed order {o}");
}

/// The wedge on polyhedral cells: the same pipe on a perturbed polygonal
/// mesh (prisms of polygons; the cells next to the axis are polyhedra with
/// an edge on it) converges and stays close to the exact profile.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn wedge_pipe_on_polyhedra() {
    let t = text("pipe_wedge.toml")
        .replace("type = \"hex\"", "type = \"poly\"")
        .replace("levels = [0, 1, 2]", "levels = [1]");
    let rep = run_case(&parse(&t), "pipe_wedge_poly");
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let e = rep.solves[0].errors.as_ref().unwrap();
    assert!(e.l2 < 2e-3, "{e:?}");
}

/// A table whose rows sit at the inlet face centres gives the run with the
/// exact profile; the table enters the case hash.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn table_inlet_matches_the_exact_profile() {
    let t = text("poiseuille_table.toml");
    let table = run_case(&parse(&t), "table");
    let exact = run_case(
        &parse(&t.replace("profile = \"table\"", "profile = \"analytic\"")),
        "table_exact",
    );
    assert!(table.failures.is_empty(), "{:?}", table.failures);
    let (a, b) = (
        table.solves[0].errors.as_ref().unwrap().l2,
        exact.solves[0].errors.as_ref().unwrap().l2,
    );
    assert!((a - b).abs() < 1e-6 * b, "table {a:e}, exact profile {b:e}");
}

/// A table that is missing or malformed is refused when the case loads.
#[test]
fn bad_tables_are_refused() {
    let t = text("poiseuille_table.toml");
    let missing = t.replace("poiseuille_profile.csv", "no_such_table.csv");
    let e = FlowCase::from_toml(&missing, &cases()).unwrap_err();
    assert_eq!(diag::of_error(&e).code, "BC019", "{e:#}");
    let dir = out("bad_table");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("t.csv"), "y, u\n0.0, 1.0\n0.5, 2.0, 3.0\n").unwrap();
    let e = FlowCase::from_toml(&t.replace("poiseuille_profile.csv", "t.csv"), &dir).unwrap_err();
    assert!(format!("{e:#}").contains("columns"), "{e:#}");
    std::fs::write(dir.join("t.csv"), "0.5, 1.0\n0.2, 2.0\n").unwrap();
    let e = FlowCase::from_toml(&t.replace("poiseuille_profile.csv", "t.csv"), &dir).unwrap_err();
    assert!(format!("{e:#}").contains("does not increase"), "{e:#}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// The whole outlet sees inflow (the flow reversed): with the default
/// backflow treatment the outlet pressure is the total pressure of the
/// entering fluid, so the static pressure is p_set − ½ρU² and the uniform
/// flow is exact; with `zero-gradient` the static pressure is p_set.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn backflow_takes_the_outlet_pressure_as_total_pressure() {
    let reversed = |mode: &str, p: f64, extra: &str| {
        let mut t = text("channel_slip.toml")
            .replace(
                "type = \"velocity-inlet\"\nprofile = \"analytic\"     # the uniform flow",
                "type = \"velocity-inlet\"\nvelocity = [-1.0, 0.0, 0.0]",
            )
            .replace(
                "pressure = 0.0\n[boundary.ymin]",
                &format!("pressure = 0.0\nbackflow = \"{mode}\"\n[boundary.ymin]"),
            )
            .replace(
                "velocity = [1.0, 0.0, 0.0]\nlength = 1.0",
                &format!("velocity = [-1.0, 0.0, 0.0]\npressure = {p}\nlength = 1.0"),
            );
        t.push_str(extra);
        parse(&t)
    };
    let rep = run_case(&reversed("total-pressure", -0.5, ""), "backflow_total");
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let e = rep.solves[0].errors.as_ref().unwrap();
    assert!(e.l2 < 1e-9 && e.pressure_linf < 1e-9, "{e:?}");
    // zero-gradient converges slowly (the reason for the default); the
    // pressure level is that of the outlet
    let mut case = reversed("zero-gradient", 0.0, "");
    case.verification.require_convergence = Some(false);
    case.verification.max_error = None;
    case.solver.max_iterations = 1500;
    let rep = run_case(&case, "backflow_zero_gradient");
    let e = rep.solves[0].errors.as_ref().unwrap();
    assert!(e.pressure_linf < 0.02, "{e:?}");
}

/// A side outlet held above the pressure in the channel: with backflow
/// suppressed the closed faces act as a slip wall, so the uniform flow
/// stays exact and no mass crosses the outlet; without suppression fluid
/// enters through it.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn suppressed_backflow_closes_the_outlet() {
    let outlet = "type = \"static-pressure-outlet\"\npressure = 0.0";
    let side = |mode: &str| {
        format!("type = \"static-pressure-outlet\"\npressure = 0.1\nbackflow = \"{mode}\"")
    };
    let rep = run_case(
        &slip_channel("mixed", outlet, &side("suppress"), "", ""),
        "suppress",
    );
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let e = rep.solves[0].errors.as_ref().unwrap();
    assert!(e.l2 < 1e-9, "{e:?}");
    let ymax = rep.solves[0]
        .boundaries
        .iter()
        .find(|p| p.name == "ymax")
        .unwrap();
    assert!(ymax.mass_flow.abs() < 1e-12, "{}", ymax.mass_flow);
    let mut open = slip_channel("mixed", outlet, &side("total-pressure"), "", "");
    open.verification.max_error = None;
    let rep = run_case(&open, "suppress_open");
    let ymax = rep.solves[0]
        .boundaries
        .iter()
        .find(|p| p.name == "ymax")
        .unwrap();
    assert!(ymax.mass_flow < -1e-3, "{}", ymax.mass_flow);
}

/// The transient channel with a convective outlet and no pressure
/// boundary: the outlet flow equals the inflow to round-off in every step
/// and the start-up disturbance leaves through the outlet.
fn convective_channel(steps: usize, extra: &str) -> FlowCase {
    let mut t = text("channel_slip.toml")
        .replace(
            "[boundary.xmax]\ntype = \"static-pressure-outlet\"\npressure = 0.0",
            "[boundary.xmax]\ntype = \"convective-outlet\"",
        )
        .replace("max_iterations = 2000\n", "")
        .replace(
            "[analytic]",
            &format!("[time]\nscheme = \"bdf2\"\ntime_step = 0.05\nsteps = {steps}\n\n[analytic]"),
        );
    t.push_str(extra);
    parse(&t)
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn convective_outlet_balances_the_flow_and_lets_the_start_leave() {
    let mut case = convective_channel(80, "");
    case.verification.max_error = Some(1e-4);
    let rep = run_case(&case, "convective");
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let s = &rep.solves[0];
    assert!(s.mass_imbalance < 1e-12, "{}", s.mass_imbalance);
    let xmax = s.boundaries.iter().find(|p| p.name == "xmax").unwrap();
    let xmin = s.boundaries.iter().find(|p| p.name == "xmin").unwrap();
    assert!((xmax.mass_flow + xmin.mass_flow).abs() < 1e-12 * xmax.mass_flow);
}

/// The checks of the new conditions.
#[test]
fn checks_of_the_new_conditions() {
    // a convective outlet needs a transient case
    let steady = text("channel_slip.toml").replace(
        "type = \"static-pressure-outlet\"\npressure = 0.0",
        "type = \"convective-outlet\"",
    );
    let (ok, all) = check_lines(&parse(&steady));
    assert!(!ok && all.contains("BC021"), "{all}");
    // without inflow it needs its convection velocity
    let t = text("channel_slip.toml")
        .replace(
            "type = \"static-pressure-outlet\"\npressure = 0.0",
            "type = \"convective-outlet\"",
        )
        .replace(
            "type = \"velocity-inlet\"\nprofile = \"analytic\"     # the uniform flow",
            "type = \"convective-outlet\"",
        )
        .replace(
            "[analytic]",
            "[time]\ntime_step = 0.05\nsteps = 2\n\n[analytic]",
        );
    let (ok, all) = check_lines(&parse(&t));
    assert!(!ok && all.contains("BC022"), "{all}");
    // one wedge boundary is not a wedge
    let one = text("pipe_wedge.toml").replace(
        "[boundary.zmax]\ntype = \"wedge\"",
        "[boundary.zmax]\ntype = \"symmetry\"",
    );
    let (ok, all) = check_lines(&parse(&one));
    assert!(!ok && all.contains("BC023"), "{all}");
    // a wide wedge warns
    let wide = text("pipe_wedge.toml").replace("wedge = 0.5 ", "wedge = 20.0 ");
    let (ok, all) = check_lines(&parse(&wide));
    assert!(ok && all.contains("BC024"), "{all}");
}

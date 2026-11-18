//! Stage 6 part 6g (specification 7.1, D29): the turbulence of an inlet
//! from the columns of its table profile, `k` with `epsilon` or `omega`
//! and `zeta`, so that a profile from a periodic pipe or channel run
//! (D111) feeds a jet or a duct.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "foehncraft_inlet_table_{name}_{}",
        std::process::id()
    ))
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

/// The flat plate of ζ-f₀ on a small mesh, a few iterations.
fn small_plate() -> String {
    text("tmr_plate_zeta_f0.toml")
        .replace("cells = [64, 48, 1]", "cells = [16, 12, 1]")
        .replace("cells = 16, ratio = 0.1", "cells = 4, ratio = 0.1")
        .replace("cells = 48, ratio = 20.0", "cells = 12, ratio = 20.0")
        .replace("max_iterations = 4000", "max_iterations = 15")
        .replace("require_convergence = true", "require_convergence = false")
}

/// The inlet `xmin` of the plate as a table profile of the file `file`.
fn with_table(case: &str, file: &str) -> String {
    let old = "[boundary.xmin]\ntype = \"velocity-inlet\"\nvelocity = [1.0, 0.0, 0.0]\n\
               k = 2.25e-7\nepsilon = 2.53125e-6\n";
    assert!(case.contains(old), "the inlet of the plate case changed");
    case.replace(
        old,
        &format!(
            "[boundary.xmin]\ntype = \"velocity-inlet\"\nprofile = \"table\"\ntable = \"{file}\"\n\
             table_coordinate = \"y\"\n"
        ),
    )
}

/// A table with constant columns gives the run of the constant inlet,
/// bit for bit: the velocity, k, ε and ζ reach the faces as the numbers
/// of the boundary; ten times the k of the table gives another run.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn a_table_of_constant_columns_is_the_constant_inlet() {
    let dir = out("tables");
    std::fs::create_dir_all(&dir).unwrap();
    let plate = small_plate();
    let constant = FlowCase::from_toml(&plate, &cases()).unwrap();
    let a = run_case(&constant, "constant");
    let zeta = 2.0f64 / 3.0;
    std::fs::write(
        dir.join("eps.csv"),
        format!(
            "# the free stream of the plate\ny, zeta, u_x, u_y, u_z, epsilon, k\n\
             0.0, {zeta:?}, 1.0, 0.0, 0.0, 2.53125e-6, 2.25e-7\n\
             1.0, {zeta:?}, 1.0, 0.0, 0.0, 2.53125e-6, 2.25e-7\n"
        ),
    )
    .unwrap();
    let table = FlowCase::from_toml(&with_table(&plate, "eps.csv"), &dir).unwrap();
    let b = run_case(&table, "table");
    assert!(
        a.failures.is_empty() && b.failures.is_empty(),
        "{:?} {:?}",
        a.failures,
        b.failures
    );
    assert_eq!(a.solves[0].iterations, b.solves[0].iterations);
    assert_eq!(
        format!("{:?}", a.solves[0].residuals),
        format!("{:?}", b.solves[0].residuals)
    );
    let forces = |r: &FlowReport| {
        r.solves[0]
            .boundaries
            .iter()
            .map(|s| (s.name.clone(), s.force, s.mass_flow))
            .collect::<Vec<_>>()
    };
    assert_eq!(forces(&a), forces(&b));
    // the columns are used: ten times the k gives another run
    std::fs::write(
        dir.join("k10.csv"),
        format!(
            "y, zeta, u_x, u_y, u_z, epsilon, k\n\
             0.0, {zeta:?}, 1.0, 0.0, 0.0, 2.53125e-6, 2.25e-6\n\
             1.0, {zeta:?}, 1.0, 0.0, 0.0, 2.53125e-6, 2.25e-6\n"
        ),
    )
    .unwrap();
    let c = run_case(
        &FlowCase::from_toml(&with_table(&plate, "k10.csv"), &dir).unwrap(),
        "k10",
    );
    assert_ne!(
        format!("{:?}", a.solves[0].residuals),
        format!("{:?}", c.solves[0].residuals)
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// A table that also gets turbulence values from the boundary, or whose
/// header names a column that a profile does not have, is refused when
/// the case loads.
#[test]
fn conflicting_or_unknown_columns_are_refused() {
    let dir = out("refused");
    std::fs::create_dir_all(&dir).unwrap();
    let plate = small_plate();
    std::fs::write(
        dir.join("t.csv"),
        "y, u_n, k, epsilon\n0.0, 1.0, 1e-6, 1e-6\n1.0, 1.0, 1e-6, 1e-6\n",
    )
    .unwrap();
    let both = with_table(&plate, "t.csv").replace(
        "table_coordinate = \"y\"\n",
        "table_coordinate = \"y\"\nturbulence_intensity = 0.01\n",
    );
    let e = FlowCase::from_toml(&both, &dir)
        .and_then(|c| run(&c, "both", &out("both"), |_| {}).map(|_| ()))
        .unwrap_err();
    assert_eq!(diag::of_error(&e).code, "BC029", "{e:#}");
    assert!(format!("{e:#}").contains("turbulence columns"), "{e:#}");
    let zeta = with_table(&plate, "t.csv").replace(
        "table_coordinate = \"y\"\n",
        "table_coordinate = \"y\"\nzeta = 0.5\n",
    );
    std::fs::write(
        dir.join("t.csv"),
        "y, u_n, k, epsilon, zeta\n0.0, 1.0, 1e-6, 1e-6, 0.5\n1.0, 1.0, 1e-6, 1e-6, 0.5\n",
    )
    .unwrap();
    let e = FlowCase::from_toml(&zeta, &dir)
        .and_then(|c| run(&c, "zeta", &out("zeta"), |_| {}).map(|_| ()))
        .unwrap_err();
    assert_eq!(diag::of_error(&e).code, "BC029", "{e:#}");
    std::fs::write(
        dir.join("t.csv"),
        "y, u_n, k, eps\n0.0, 1.0, 1e-6, 1e-6\n1.0, 1.0, 1e-6, 1e-6\n",
    )
    .unwrap();
    let e = FlowCase::from_toml(&with_table(&plate, "t.csv"), &dir).unwrap_err();
    assert_eq!(diag::of_error(&e).code, "BC019", "{e:#}");
    assert!(format!("{e:#}").contains("`eps`"), "{e:#}");
    let _ = std::fs::remove_dir_all(&dir);
}

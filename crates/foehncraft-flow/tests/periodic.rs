//! Stage 6, part 6p: translational periodic boundaries (D29) for the flow
//! and the passive scalars, with the heat-flux jump of D109 (1),
//! docs/verification/stage6.md.
//!
//! The case files of the gallery (`poiseuille_periodic.toml`,
//! `channel_heated_periodic.toml`, `pipe_heated_periodic.toml` on three
//! meshes each) also run in CI; the tests here use fewer meshes.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_periodic_{name}_{}", std::process::id()))
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

/// A case file with `\n` line ends (a Windows checkout has `\r\n`).
fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

/// The code and message a case text is refused with.
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

/// Periodic plane Poiseuille flow equals the open channel (spec 6p
/// acceptance): driven by its mass flow, the velocity error is that of the
/// open channel within a factor of two on the same meshes and the driving
/// force converges to 12 μ U / h² = 1.2 Pa/m; driven by that pressure
/// gradient on uniform meshes, the mass flow converges to 0.1 kg/s at
/// second order.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn periodic_poiseuille_equals_the_open_channel() {
    let two = |t: String| t.replace("levels = [0, 1, 2]", "levels = [0, 1]");
    let per = run_case(&parse(&two(text("poiseuille_periodic.toml"))), "mass");
    let open = run_case(&parse(&two(text("poiseuille_hex.toml"))), "open");
    assert!(per.failures.is_empty(), "{:?}", per.failures);
    for (p, o) in per.solves.iter().zip(&open.solves) {
        let (ep, eo) = (p.errors.as_ref().unwrap().l2, o.errors.as_ref().unwrap().l2);
        assert!(ep < 2.0 * eo, "periodic {ep:.3e}, open {eo:.3e}");
        let d = &p.periodic[0];
        assert_eq!((d.from.as_str(), d.to.as_str()), ("xmin", "xmax"));
        assert!((d.mass_flow - 0.1).abs() < 1e-12, "{}", d.mass_flow);
    }
    let beta = per.solves[1].periodic[0].driving_force;
    assert!((beta / 1.2 - 1.0).abs() < 2e-3, "driving force {beta}");

    let pg = two(text("poiseuille_periodic.toml"))
        .replace("mass_flow = 0.1 ", "pressure_gradient = 1.2 ")
        .replace("perturb = 0.15", "perturb = 0.0")
        .replace("expected_order = [1.8, 2.2]\n", "");
    let r = run_case(&parse(&pg), "gradient");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let e: Vec<f64> = r
        .solves
        .iter()
        .map(|s| {
            assert_eq!(s.periodic[0].driving_force, 1.2);
            (s.periodic[0].mass_flow / 0.1 - 1.0).abs()
        })
        .collect();
    let order = (e[0] / e[1]).log2();
    assert!((1.8..=2.2).contains(&order), "mass-flow errors {e:?}");
}

/// The periodic heated channel (D111, spec 9.2): the heat-flux jump of
/// D109 (1) gives the temperature of the open channel; Nu approaches
/// 140/17 at second order, the jump is Q/ṁ = 4e-3 K and the bulk value at
/// the periodic plane is the given level.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn periodic_heated_channel_nusselt_number() {
    // two meshes: the first pair of the sweep converges slightly faster
    // than second order
    let t = text("channel_heated_periodic.toml")
        .replace("levels = [0, 1, 2]", "levels = [0, 1]")
        .replace("scalar_order = [1.8, 2.2]", "scalar_order = [1.8, 2.6]");
    let r = run_case(&parse(&t), "heated");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let (g, l, q, gamma) = (2e-3, 2.0, 1e-3, 1e-3);
    let exact = 140.0 / 17.0;
    let nu: Vec<f64> = r
        .solves
        .iter()
        .map(|s| {
            let sc = &s.scalars[0];
            let b = |n: &str| sc.boundaries.iter().find(|b| b.boundary == n).unwrap();
            let (ta, tb) = (b("xmin").mixed_mean.unwrap(), b("xmax").mixed_mean.unwrap());
            assert!((ta + 17.0 / 70.0).abs() < 1e-12, "bulk at xmin {ta}");
            assert!((tb - ta - g * l).abs() < 1e-12, "jump {}", tb - ta);
            assert!((b("ymin").flux - q).abs() < 1e-12, "{}", b("ymin").flux);
            let tw = 0.5 * (b("ymin").mean + b("ymax").mean) - g * l / 2.0;
            2.0 * q / (gamma * (tw - ta))
        })
        .collect();
    let order = ((nu[0] - exact) / (nu[1] - exact)).log2();
    assert!((1.6..=2.6).contains(&order), "Nu {nu:?}");
    assert!((nu[1] / exact - 1.0).abs() < 1e-2, "Nu {nu:?}");
    let o = &r.scalar_orders[0].1;
    assert!((1.8..=2.6).contains(&o[0]), "{o:?}");
}

/// Wrong periodic pairs and periodic scalar tables are refused with their
/// codes.
#[test]
fn bad_periodic_conditions_are_refused() {
    let t = text("channel_heated_periodic.toml");
    // a periodic boundary without a periodic partner
    let (code, msg) = refused(&t.replace(
        "[boundary.xmax]\ntype = \"periodic\"\n",
        "[boundary.xmax]\ntype = \"wall\"\n",
    ));
    assert_eq!(code, "BC030", "{msg}");
    // the driving on both sides
    let (code, msg) = refused(&t.replace(
        "[boundary.xmax]\ntype = \"periodic\"\n",
        "[boundary.xmax]\ntype = \"periodic\"\nmass_flow = 0.1\n",
    ));
    assert_eq!(code, "BC030", "{msg}");
    // a pair whose faces do not match (y periodic with x)
    let (code, msg) = refused(
        &t.replace("partner = \"xmax\"", "partner = \"ymax\"")
            .replace(
                "[boundary.xmax]\ntype = \"periodic\"\n",
                "[boundary.xmax]\ntype = \"periodic\"\npartner = \"ymin\"\n",
            )
            .replace(
                "[boundary.ymin]\ntype = \"wall\"\n[boundary.ymin.scalar.theta]\nflux = \"q\"\n",
                "[boundary.ymin]\ntype = \"periodic\"\n",
            )
            .replace(
                "[boundary.ymax]\ntype = \"wall\"\n[boundary.ymax.scalar.theta]\nflux = \"q\"\n",
                "[boundary.ymax]\ntype = \"periodic\"\n",
            ),
    );
    assert_eq!(code, "BC030", "{msg}");
    // a value on a periodic boundary
    let (code, msg) = refused(&t.replace("periodic = \"heat-flux\"", "value = 0.0"));
    assert_eq!(code, "SCAL001", "{msg}");
    // the mode on both sides of the pair
    let (code, msg) = refused(&t.replace(
        "[boundary.xmax]\ntype = \"periodic\"\n",
        "[boundary.xmax]\ntype = \"periodic\"\n[boundary.xmax.scalar.theta]\nperiodic = \"periodic\"\n",
    ));
    assert_eq!(code, "SCAL001", "{msg}");
    // the heat-flux mode with a fixed wall value (D109 (4))
    let (code, msg) = refused(&t.replace(
        "[boundary.ymax.scalar.theta]\nflux = \"q\"",
        "[boundary.ymax.scalar.theta]\nvalue = 0.0",
    ));
    assert_eq!(code, "SCAL001", "{msg}");
    assert!(msg.contains("heat-flux"), "{msg}");
    // `periodic` on a wall
    let (code, msg) = refused(&t.replace(
        "[boundary.ymax.scalar.theta]\nflux = \"q\"",
        "[boundary.ymax.scalar.theta]\nperiodic = \"periodic\"",
    ));
    assert_eq!(code, "SCAL001", "{msg}");
}

/// A periodic pair that nothing drives gives the warning BC031 (a
/// decaying flow), and the resolved case reads back to the same case.
#[test]
fn undriven_pair_warns_and_resolved_case_reads_back() {
    let t = text("poiseuille_periodic.toml");
    let c = parse(&t);
    let back = FlowCase::from_toml(&c.resolved_toml().unwrap(), &cases()).unwrap();
    assert_eq!(back.boundary, c.boundary);
    let quiet = t
        .replace("mass_flow = 0.1 ", "# ")
        .replace("levels = [0, 1, 2]", "levels = [0]")
        .replace("max_iterations = 5000", "max_iterations = 3")
        .replace("require_convergence = true", "require_convergence = false")
        .replace("expected_order = [1.8, 2.2]\n", "");
    let r = run_case(&parse(&quiet), "undriven");
    assert!(
        r.findings.iter().any(|w| w.code == "BC031"),
        "{:?}",
        r.findings
    );
}

/// A transient periodic run with a heat-flux scalar continued from a
/// restart file ends bit-identically to the run without the interruption
/// (D42): the driving force is part of the restart file.
#[cfg(feature = "hdf5")]
#[test]
fn periodic_restart_continues_bit_identically() {
    use foehncraft_flow::restart;
    let case = |extra: &str| {
        let t = text("channel_heated_periodic.toml")
            .replace("cells = [32, 16, 1]", "cells = [8, 6, 1]")
            .replace("levels = [0, 1, 2]", "levels = [0]")
            .replace("expected_order = [1.8, 2.2]\n", "")
            .replace("scalar_order = [1.8, 2.2]\n", "")
            .replace("require_convergence = true", "require_convergence = false")
            .replace("relaxation = 0.9", "relaxation = 1.0");
        parse(&format!(
            "{t}\n[time]\nscheme = \"bdf2\"\ntime_step = 0.05\nsteps = 12\n{extra}\n"
        ))
    };
    let da = out("restart_full");
    let ra = run(&case("[run]\nrestart_at = [6]"), "full", &da, |_| {}).unwrap();
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    let from = da.join("restart_000006.h5");
    let db = out("restart_cont");
    let rb = run(
        &case(&format!("[run]\nrestart = '{}'", from.display())),
        "cont",
        &db,
        |_| {},
    )
    .unwrap();
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    let b = restart::read(&db.join("restart.h5")).unwrap();
    assert_eq!(a.state.periodic_force.len(), 1);
    assert!(a.state.periodic_force[0] > 0.0);
    assert_eq!(a.state.periodic_force, b.state.periodic_force);
    assert_eq!(
        a.state.scalars, b.state.scalars,
        "the continued scalar differs"
    );
    assert_eq!(a.state.u.x, b.state.u.x);
    assert_eq!(a.state.p, b.state.p);
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

//! Stage 6, part 6d: the k-ε-ζ-f model (D132), its wall treatments
//! (D134), the elliptic f equation and the scalar multigrid of the
//! turbulence equations (D139), docs/verification/stage6.md.
//!
//! The channel and the flat plate of the gallery (`channel_zeta_f.toml`,
//! `tmr_plate_zeta_f.toml`) are the validation cases; the tests here use
//! small meshes.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_zeta_f_{name}_{}", std::process::id()))
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

/// The channel of the gallery on 24 cells with few iterations.
fn small_channel(extra: &str) -> String {
    text("channel_zeta_f.toml")
        .replace("cells = [2, 96, 1]", "cells = [2, 24, 1]")
        .replace("ratio = 60.0", "ratio = 8.0")
        .replace("algorithm = \"coupled\"\n", "")
        .replace("max_iterations = 1000", "max_iterations = 30")
        + "\n[verification]\nrequire_convergence = false\n"
        + extra
}

/// The box at rest of the decay tests, with the initial values of
/// `initial`.
#[cfg(feature = "hdf5")]
fn decay_box(initial: &str, dt: f64, n: usize) -> String {
    format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.5]
cells = [2, 2, 1]
[fluid]
viscosity = 1e-3
[turbulence]
model = "k-epsilon-zeta-f"
[initial]
velocity = [0.0, 0.0, 0.0]
{initial}
[boundary.xmin]
type = "symmetry"
[boundary.xmax]
type = "symmetry"
[boundary.ymin]
type = "symmetry"
[boundary.ymax]
type = "symmetry"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver]
outer_iterations = 4
[solver.turbulence]
tolerance = 1e-12
[solver.elliptic]
tolerance = 1e-12
[time]
scheme = "bdf2"
time_step = {dt}
steps = {n}
[output]
vtkhdf = false
wall_csv = false
"#
    )
}

/// Homogeneous decay in a box at rest (specification 9.2): with P = 0 and
/// T = k/ε, k and ε follow the k-ε solution with C_ε2 = 1.9, and ζ
/// relaxes to 2/3 through f = (C₁ − 1)(2/3 − ζ)/T:
/// 2/3 − ζ = (2/3 − ζ₀)(1 + (C_ε2 − 1) ε₀ t/k₀)^(−(C₁ − 1)/(C_ε2 − 1)).
/// BDF2 is second order in time for all three (the first step is BDF1).
#[cfg(feature = "hdf5")]
#[test]
fn zeta_f_homogeneous_decay_is_second_order_in_time() {
    use foehncraft_flow::restart;
    let (k0, e0, z0, t_end): (f64, f64, f64, f64) = (1.0, 1.0, 0.3, 0.5);
    let (c2, c1): (f64, f64) = (1.9, 1.4);
    let a = 1.0 + (c2 - 1.0) * e0 * t_end / k0;
    let exact_k = k0 * a.powf(-1.0 / (c2 - 1.0));
    let exact_e = e0 * a.powf(-c2 / (c2 - 1.0));
    let exact_z = 2.0 / 3.0 - (2.0 / 3.0 - z0) * a.powf(-(c1 - 1.0) / (c2 - 1.0));
    let mut errors = Vec::new();
    for n in [10usize, 20, 40] {
        let c = parse(&decay_box(
            &format!("k = {k0}\nepsilon = {e0}\nzeta = {z0}"),
            t_end / n as f64,
            n,
        ));
        let d = out(&format!("decay{n}"));
        let r = run(&c, "decay", &d, |_| {}).unwrap();
        assert!(r.failures.is_empty(), "{:?}", r.failures);
        let rs = restart::read(&d.join("restart.h5")).unwrap();
        let t = rs.state.turbulence.unwrap();
        let _ = std::fs::remove_dir_all(&d);
        let z = t.zeta.expect("the restart file holds ζ and f");
        assert!(z.zeta_old.is_some(), "ζ of the level before");
        assert!(
            t.k.iter().all(|&k| (k - t.k[0]).abs() < 1e-14 * t.k[0]),
            "{:?}",
            t.k
        );
        errors.push((
            (t.k[0] - exact_k).abs() / exact_k,
            (t.omega[0] - exact_e).abs() / exact_e,
            (z.zeta[0] - exact_z).abs() / (2.0 / 3.0 - exact_z),
        ));
    }
    for i in 1..errors.len() {
        let pk = (errors[i - 1].0 / errors[i].0).log2();
        let pe = (errors[i - 1].1 / errors[i].1).log2();
        let pz = (errors[i - 1].2 / errors[i].2).log2();
        eprintln!(
            "decay: errors {:?}, orders k {pk:.3}, ε {pe:.3}, ζ {pz:.3}",
            errors[i]
        );
        assert!((1.8..=2.2).contains(&pk), "order of k {pk}");
        assert!((1.8..=2.2).contains(&pe), "order of ε {pe}");
        assert!((1.8..=2.2).contains(&pz), "order of ζ {pz}");
    }
    assert!(
        errors[2].0 < 1e-3 && errors[2].1 < 1e-3 && errors[2].2 < 1e-3,
        "{errors:?}"
    );
}

/// The sustaining terms hold k and ε at their ambient values in a box at
/// rest, and ζ = 2/3 with f = 0 is a fixed point (specification 9.2).
#[cfg(feature = "hdf5")]
#[test]
fn zeta_f_sustaining_terms_hold_the_ambient_state() {
    use foehncraft_flow::restart;
    let t = decay_box("k = 0.01\nepsilon = 0.002", 0.05, 10).replace(
        "model = \"k-epsilon-zeta-f\"",
        "model = \"k-epsilon-zeta-f\"\nambient = { k = 0.01, epsilon = 0.002 }",
    );
    let d = out("ambient");
    let r = run(&parse(&t), "ambient", &d, |_| {}).unwrap();
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = restart::read(&d.join("restart.h5"))
        .unwrap()
        .state
        .turbulence
        .unwrap();
    let _ = std::fs::remove_dir_all(&d);
    let z = s.zeta.unwrap();
    for c in 0..s.k.len() {
        assert!((s.k[c] / 0.01 - 1.0).abs() < 1e-10, "k {}", s.k[c]);
        assert!((s.omega[c] / 0.002 - 1.0).abs() < 1e-10, "ε {}", s.omega[c]);
        assert!((z.zeta[c] - 2.0 / 3.0).abs() < 1e-12, "ζ {}", z.zeta[c]);
        assert!(z.f[c].abs() < 1e-12, "f {}", z.f[c]);
    }
}

/// A run continued from a restart file ends bit-identically to the run
/// without the interruption (D42): ζ, f and their old levels are part of
/// the restart file. The multigrid of f (default) and of the transport
/// equations (`preconditioner = "amg"`) rebuild on the schedule of the
/// pressure, so the continued run builds the same hierarchies.
#[cfg(feature = "hdf5")]
#[test]
fn zeta_f_restart_continues_bit_identically() {
    use foehncraft_flow::restart;
    let case = |extra: &str| {
        parse(&small_channel(&format!(
            "[solver.turbulence]\npreconditioner = \"amg\"\n\
             [time]\nscheme = \"bdf2\"\ntime_step = 0.002\nsteps = 10\n{extra}\n"
        )))
    };
    let da = out("restart_full");
    let ra = run(&case("[run]\nrestart_at = [5]"), "full", &da, |_| {}).unwrap();
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    let from = da.join("restart_000005.h5");
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
    let (ta, tb) = (a.state.turbulence.unwrap(), b.state.turbulence.unwrap());
    assert!(ta.zeta.as_ref().unwrap().zeta_old.is_some());
    assert_eq!(ta, tb, "the continued turbulence state differs");
    assert_eq!(a.state.u.x, b.state.u.x);
    assert_eq!(a.state.p, b.state.p);
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// A run of k-ε-ζ-f continues from the restart file of a k-ε run with
/// ζ = 2/3 and f = 0 (specification 6.4).
#[cfg(feature = "hdf5")]
#[test]
fn keps_restart_continues_with_zeta_f() {
    use foehncraft_flow::restart;
    let keps = small_channel("")
        .replace("model = \"k-epsilon-zeta-f\"", "model = \"k-epsilon\"")
        .replace(
            "wall_treatment = \"compound\"",
            "wall_treatment = \"scalable\"",
        );
    let da = out("restart_keps");
    let ra = run(&parse(&keps), "keps", &da, |_| {}).unwrap();
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    assert!(
        restart::read(&da.join("restart.h5"))
            .unwrap()
            .state
            .turbulence
            .unwrap()
            .zeta
            .is_none()
    );
    let db = out("restart_zf");
    let rb = run(
        &parse(&small_channel(&format!(
            "[run]\nrestart = '{}'\n",
            da.join("restart.h5").display()
        ))),
        "zf",
        &db,
        |_| {},
    )
    .unwrap();
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let t = restart::read(&db.join("restart.h5"))
        .unwrap()
        .state
        .turbulence
        .unwrap();
    let z = t.zeta.expect("ζ and f in the new file");
    assert!(
        z.zeta
            .iter()
            .all(|v| v.is_finite() && (0.0..=2.0).contains(v))
    );
    assert!(z.f.iter().all(|v| v.is_finite()));
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// One and four threads give the same bits (D31) with k-ε-ζ-f, the
/// compound treatment and the multigrid of every turbulence equation, on
/// a plate mesh above the 20 000 cells below which one thread is used.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn zeta_f_results_are_bitwise_equal_on_1_and_4_threads() {
    let mk = |threads: usize| {
        let mut c = parse(
            &(text("tmr_plate_zeta_f.toml")
                .replace("cells = [64, 48, 1]", "cells = [160, 128, 1]")
                .replace("cells = 16, ratio", "cells = 40, ratio")
                .replace("cells = 48, ratio", "cells = 120, ratio")
                .replace("max_iterations = 4000", "max_iterations = 8")
                + "\n[solver.turbulence]\npreconditioner = \"amg\"\n"
                + &format!("\n[run]\nthreads = {threads}\n")),
        );
        c.verification.require_convergence = Some(false);
        c
    };
    let a = run_case(&mk(1), "threads1");
    let b = run_case(&mk(4), "threads4");
    assert!(a.solves[0].cells >= 20_000);
    let (ra, rb) = (&a.solves[0].residuals, &b.solves[0].residuals);
    assert_eq!(ra.momentum.map(f64::to_bits), rb.momentum.map(f64::to_bits));
    assert_eq!(ra.continuity.to_bits(), rb.continuity.to_bits());
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(ra.turbulence.residual.len(), 4, "k, ε, f, ζ");
    assert_eq!(bits(&ra.turbulence.residual), bits(&rb.turbulence.residual));
    let (wa, wb) = (&a.solves[0].wall_probes[0], &b.solves[0].wall_probes[0]);
    assert_eq!(wa.magnitude.to_bits(), wb.magnitude.to_bits());
}

/// The channel of the gallery on its own mesh: U_b⁺ from the driving force
/// within the ±3 % band (specification 9.4); the low-Re treatment runs as
/// well.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn zeta_f_channel_meets_its_band() {
    let r = run_case(&parse(&text("channel_zeta_f.toml")), "channel");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let f = r.solves[0].periodic[0].driving_force;
    let ub_plus = 17.55 / f.sqrt();
    let dev = ub_plus / 17.55 - 1.0;
    eprintln!(
        "channel k-ε-ζ-f compound: U_b⁺ {ub_plus:.3} ({:+.2} %)",
        100.0 * dev
    );
    assert!(dev.abs() < 0.03, "U_b⁺ {ub_plus}");
    let c = parse(&small_channel("").replace(
        "wall_treatment = \"compound\"",
        "wall_treatment = \"low-re\"",
    ));
    let r = run_case(&c, "channel_low_re");
    let s = &r.solves[0];
    assert_eq!(s.residuals.turbulence.residual.len(), 4);
    assert!(
        s.residuals
            .turbulence
            .residual
            .iter()
            .all(|r| r.is_finite())
    );
}

/// Wrong combinations are refused with their codes.
#[test]
fn bad_zeta_f_input_is_refused() {
    let t = small_channel("");
    // a wall function of the k-ε family
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"compound\"",
        "wall_treatment = \"scalable\"",
    ));
    assert_eq!(code, "CASE005", "{msg}");
    // a production limiter
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"compound\"",
        "wall_treatment = \"compound\"\nproduction_limiter = \"kato-launder\"",
    ));
    assert_eq!(code, "CASE005", "{msg}");
    // ambient values of SST
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"compound\"",
        "wall_treatment = \"compound\"\nambient = { k = 1e-3, omega = 1.0 }",
    ));
    assert_eq!(code, "CASE005", "{msg}");
    // an unknown constant (one of the k-ε family that ζ-f does not have)
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"compound\"",
        "wall_treatment = \"compound\"\nconstants = { c_epsilon3 = 1.0 }",
    ));
    assert_eq!(code, "CASE002", "{msg}");
}

/// ζ at an inlet: in (0, 2] for k-ε-ζ-f; a warning (BC028) for another
/// model, which ignores it.
#[test]
fn inlet_zeta_checks() {
    let base = text("tmr_plate_zeta_f.toml")
        .replace("require_convergence = true", "require_convergence = false")
        .replace("max_iterations = 4000", "max_iterations = 1");
    let inlet = "k = 2.25e-7\nepsilon = 2.53125e-6\n[boundary.xmax]";
    let (code, msg) = refused(&base.replace(
        inlet,
        "k = 2.25e-7\nepsilon = 2.53125e-6\nzeta = 2.5\n[boundary.xmax]",
    ));
    // one of the findings of the boundary-condition check
    assert!(
        msg.contains("CASE004") && msg.contains("zeta = 2.5"),
        "{code}: {msg}"
    );
    let r = run_case(
        &parse(&base.replace(
            inlet,
            "k = 2.25e-7\nepsilon = 2.53125e-6\nzeta = 0.5\n[boundary.xmax]",
        )),
        "inlet_zeta",
    );
    assert!(r.failures.is_empty(), "{:?}", r.failures);
}

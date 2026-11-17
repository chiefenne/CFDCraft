//! Stage 6, part 6c: the standard and realizable k-ε models (D131), their
//! wall treatments (D134), ε at inlets and in restart files,
//! docs/verification/stage6.md.
//!
//! The flat plate, the channel and the rough pipe of the gallery
//! (`tmr_plate_keps.toml`, `channel_keps.toml`, `pipe_rough_keps.toml`)
//! are the validation cases; the tests here use small meshes.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_keps_{name}_{}", std::process::id()))
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

/// The channel of the gallery with few iterations.
fn small_channel(extra: &str) -> String {
    text("channel_keps.toml").replace("max_iterations = 3000", "max_iterations = 40")
        + "\n[verification]\nrequire_convergence = false\n"
        + extra
}

/// The box at rest of the decay tests, with `model` and the initial
/// values of `initial`.
#[cfg(feature = "hdf5")]
fn decay_box(model: &str, initial: &str, dt: f64, n: usize) -> String {
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
model = "{model}"
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
[solver.turbulence]
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

/// Homogeneous decay of the standard k-ε in a box at rest (specification
/// 9.2): k = k₀ (1 + (C_ε2 − 1) ε₀ t/k₀)^(−1/(C_ε2 − 1)) and
/// ε = ε₀ (1 + (C_ε2 − 1) ε₀ t/k₀)^(−C_ε2/(C_ε2 − 1)). BDF2 is second
/// order in time (the first step is BDF1).
#[cfg(feature = "hdf5")]
#[test]
fn keps_homogeneous_decay_is_second_order_in_time() {
    use foehncraft_flow::restart;
    let (k0, e0, t_end): (f64, f64, f64) = (1.0, 1.0, 0.5);
    let c2: f64 = 1.92;
    let a = 1.0 + (c2 - 1.0) * e0 * t_end / k0;
    let exact_k = k0 * a.powf(-1.0 / (c2 - 1.0));
    let exact_e = e0 * a.powf(-c2 / (c2 - 1.0));
    let mut errors = Vec::new();
    for n in [10usize, 20, 40] {
        let c = parse(&decay_box(
            "k-epsilon",
            &format!("k = {k0}\nepsilon = {e0}"),
            t_end / n as f64,
            n,
        ));
        let d = out(&format!("decay{n}"));
        let r = run(&c, "decay", &d, |_| {}).unwrap();
        assert!(r.failures.is_empty(), "{:?}", r.failures);
        let rs = restart::read(&d.join("restart.h5")).unwrap();
        let t = rs.state.turbulence.unwrap();
        let _ = std::fs::remove_dir_all(&d);
        assert!(t.epsilon, "the restart file holds ε");
        assert!(
            t.k.iter().all(|&k| (k - t.k[0]).abs() < 1e-14 * t.k[0]),
            "{:?}",
            t.k
        );
        errors.push((
            (t.k[0] - exact_k).abs() / exact_k,
            (t.omega[0] - exact_e).abs() / exact_e,
        ));
    }
    for i in 1..errors.len() {
        let pk = (errors[i - 1].0 / errors[i].0).log2();
        let pe = (errors[i - 1].1 / errors[i].1).log2();
        eprintln!("decay: errors {:?}, orders k {pk:.3}, ε {pe:.3}", errors[i]);
        assert!((1.8..=2.2).contains(&pk), "order of k {pk}");
        assert!((1.8..=2.2).contains(&pe), "order of ε {pe}");
    }
    assert!(errors[2].0 < 1e-3 && errors[2].1 < 1e-3, "{errors:?}");
}

/// The sustaining terms hold k and ε at their ambient values in a box at
/// rest (specification 9.2, both models).
#[cfg(feature = "hdf5")]
#[test]
fn keps_sustaining_terms_hold_the_ambient_state() {
    use foehncraft_flow::restart;
    for model in ["k-epsilon", "realizable-k-epsilon"] {
        let t = decay_box(model, "k = 0.01\nepsilon = 0.002", 0.05, 10).replace(
            &format!("model = \"{model}\""),
            &format!("model = \"{model}\"\nambient = {{ k = 0.01, epsilon = 0.002 }}"),
        );
        let d = out(&format!("ambient_{model}"));
        let r = run(&parse(&t), "ambient", &d, |_| {}).unwrap();
        assert!(r.failures.is_empty(), "{:?}", r.failures);
        let s = restart::read(&d.join("restart.h5"))
            .unwrap()
            .state
            .turbulence
            .unwrap();
        let _ = std::fs::remove_dir_all(&d);
        for (&k, &e) in s.k.iter().zip(&s.omega) {
            assert!((k / 0.01 - 1.0).abs() < 1e-10, "{model}: k {k}");
            assert!((e / 0.002 - 1.0).abs() < 1e-10, "{model}: ε {e}");
        }
    }
}

/// A run of the k-ε model continues from the restart file of an SST run:
/// ε = C_μ k ω on the way in, the new file holds ε.
#[cfg(feature = "hdf5")]
#[test]
fn sst_restart_continues_with_k_epsilon() {
    use foehncraft_flow::restart;
    let sst = text("channel_sst.toml")
        .replace("cells = [2, 96, 1]", "cells = [2, 12, 1]")
        .replace("grading = { y = { ratio = 60.0, symmetric = true } }", "")
        .replace(
            "wall_treatment = \"low-re\"",
            "wall_treatment = \"automatic\"",
        )
        .replace("max_iterations = 3000", "max_iterations = 20")
        + "\n[verification]\nrequire_convergence = false\n";
    let da = out("restart_sst");
    let ra = run(&parse(&sst), "sst", &da, |_| {}).unwrap();
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    assert!(!a.state.turbulence.as_ref().unwrap().epsilon);
    let db = out("restart_keps");
    let rb = run(
        &parse(&small_channel(&format!(
            "[run]\nrestart = '{}'\n",
            da.join("restart.h5").display()
        ))),
        "keps",
        &db,
        |_| {},
    )
    .unwrap();
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let b = restart::read(&db.join("restart.h5")).unwrap();
    let t = b.state.turbulence.unwrap();
    assert!(t.epsilon);
    assert!(
        t.k.iter()
            .chain(&t.omega)
            .all(|v| v.is_finite() && *v > 0.0)
    );
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// One and four threads give the same bits (D31) with the k-ε model and
/// the compound treatment, on a plate mesh above the 20 000 cells below
/// which one thread is used.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn keps_results_are_bitwise_equal_on_1_and_4_threads() {
    let mk = |threads: usize| {
        let mut c = parse(
            &(text("tmr_plate_keps.toml")
                .replace("cells = [64, 48, 1]", "cells = [160, 128, 1]")
                .replace("cells = 16, ratio", "cells = 40, ratio")
                .replace("cells = 48, ratio", "cells = 120, ratio")
                .replace(
                    "wall_treatment = \"scalable\"",
                    "wall_treatment = \"compound\"",
                )
                .replace("max_iterations = 3000", "max_iterations = 8")
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
    assert_eq!(bits(&ra.turbulence.residual), bits(&rb.turbulence.residual));
    let (wa, wb) = (&a.solves[0].wall_probes[0], &b.solves[0].wall_probes[0]);
    assert_eq!(wa.magnitude.to_bits(), wb.magnitude.to_bits());
}

/// The channel of the gallery on its own mesh: U_b⁺ from the driving force
/// within the ±5 % band of the scalable wall functions (specification
/// 9.4), realizable and standard treatment run as well.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn keps_channel_meets_its_band() {
    let r = run_case(&parse(&text("channel_keps.toml")), "channel");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let f = r.solves[0].periodic[0].driving_force;
    let ub_plus = 17.55 / f.sqrt();
    let dev = ub_plus / 17.55 - 1.0;
    eprintln!(
        "channel k-ε scalable: U_b⁺ {ub_plus:.3} ({:+.2} %)",
        100.0 * dev
    );
    assert!(dev.abs() < 0.05, "U_b⁺ {ub_plus}");
    for (model, treatment) in [
        ("realizable-k-epsilon", "scalable"),
        ("k-epsilon", "standard"),
    ] {
        let c = parse(
            &small_channel("")
                .replace("model = \"k-epsilon\"", &format!("model = \"{model}\""))
                .replace(
                    "wall_treatment = \"scalable\"",
                    &format!("wall_treatment = \"{treatment}\""),
                ),
        );
        let r = run_case(&c, &format!("channel_{model}_{treatment}"));
        let s = &r.solves[0];
        let t = s.turbulence.as_ref().expect("turbulence summary");
        assert_eq!(t.walls.len(), 2);
        assert!(
            s.residuals
                .turbulence
                .residual
                .iter()
                .all(|r| r.is_finite()),
            "{model} {treatment}"
        );
    }
}

/// Wrong combinations are refused with their codes.
#[test]
fn bad_k_epsilon_input_is_refused() {
    let t = small_channel("");
    // a treatment of SST
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"scalable\"",
        "wall_treatment = \"low-re\"",
    ));
    assert_eq!(code, "CASE005", "{msg}");
    assert!(msg.contains("k-ε family"), "{msg}");
    // a treatment of the k-ε family with SST
    let (code, msg) = refused(&t.replace(
        "model = \"k-epsilon\"\nwall_treatment = \"scalable\"",
        "model = \"k-omega-sst\"\nwall_treatment = \"compound\"",
    ));
    assert_eq!(code, "CASE005", "{msg}");
    // the realizable model with a production limiter
    let (code, msg) = refused(&t.replace(
        "model = \"k-epsilon\"",
        "model = \"realizable-k-epsilon\"\nproduction_limiter = \"kato-launder\"",
    ));
    assert_eq!(code, "CASE005", "{msg}");
    // ambient values of SST
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"scalable\"",
        "wall_treatment = \"scalable\"\nambient = { k = 1e-3, omega = 1.0 }",
    ));
    assert_eq!(code, "CASE005", "{msg}");
    // an unknown constant of the k-ε family
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"scalable\"",
        "wall_treatment = \"scalable\"\nconstants = { c_epsilon3 = 1.0 }",
    ));
    assert_eq!(code, "CASE002", "{msg}");
}

/// An inlet takes ε or ω (converted with C_μ), never both (BC029).
#[test]
fn inlet_epsilon_checks() {
    let base = text("tmr_plate_keps.toml")
        .replace("require_convergence = true", "require_convergence = false")
        .replace("max_iterations = 3000", "max_iterations = 1");
    let (code, msg) = refused(&base.replace(
        "k = 2.25e-7\nepsilon = 2.53125e-6\n[boundary.xmax]",
        "k = 2.25e-7\nepsilon = 2.53125e-6\nomega = 125.0\n[boundary.xmax]",
    ));
    assert_eq!(code, "BC029", "{msg}");
    let r = run_case(
        &parse(&base.replace(
            "k = 2.25e-7\nepsilon = 2.53125e-6\n[boundary.xmax]",
            "k = 2.25e-7\nomega = 125.0\n[boundary.xmax]",
        )),
        "inlet_omega",
    );
    assert!(r.failures.is_empty(), "{:?}", r.failures);
}

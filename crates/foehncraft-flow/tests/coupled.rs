//! The coupled pressure–velocity algorithm and the pseudo-time steps of a
//! steady run (D176 to D178, docs/verification/stage6.md, part 6e).
//!
//! The backward-facing step with the three turbulence models and the
//! cavity at Re = 1000 run from the case files in CI; the tests here use
//! small meshes.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_coupled_{name}_{}", std::process::id()))
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

fn cavity_text(n: usize, viscosity: f64, extra: &str) -> String {
    format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [{n}, {n}, 1]
[fluid]
viscosity = {viscosity}
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
    )
}

/// The lid-driven cavity on `n` x `n` hexahedra, one cell layer.
fn cavity(n: usize, viscosity: f64, extra: &str) -> FlowCase {
    FlowCase::from_toml(&cavity_text(n, viscosity, extra), Path::new(".")).unwrap()
}

/// The code and message a case text is refused with.
fn refused(text: &str) -> (String, String) {
    let e = match FlowCase::from_toml(text, Path::new(".")) {
        Err(e) => e,
        Ok(c) => match run(&c, "refused", &out("refused"), |_| {}) {
            Err(e) => e,
            Ok(_) => panic!("the case was accepted"),
        },
    };
    (diag::of_error(&e).code.to_string(), format!("{e:#}"))
}

fn lid_force(rep: &FlowReport) -> [f64; 3] {
    rep.solves[0]
        .boundaries
        .iter()
        .find(|p| p.name == "ymax")
        .unwrap()
        .force
}

/// D176: the coupled algorithm solves the discrete equations of SIMPLEC.
/// The converged result is the same, whatever the pseudo-time step
/// (automatic, given, none) and the cap α_u on the local step: the face
/// flux depends on neither (D27 with the current iterate as the old
/// level).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn coupled_result_is_the_simplec_result() {
    let force = |name: &str, solver: &str| {
        let case = cavity(
            16,
            0.01,
            &format!("[solver]\ntolerance = 1e-11\nmax_iterations = 5000\n{solver}"),
        );
        let rep = run_case(&case, name);
        assert!(rep.solves[0].converged, "{name} did not converge");
        (lid_force(&rep), rep.solves[0].iterations)
    };
    let (reference, n_simplec) = force("simplec", "[solver.momentum]\nrelaxation = 0.8\n");
    let variants = [
        ("automatic", "algorithm = \"coupled\"\n".to_string()),
        (
            "small_step",
            "algorithm = \"coupled\"\n[solver.pseudo_time]\nstep = 0.02\n".to_string(),
        ),
        (
            "large_step",
            "algorithm = \"coupled\"\n[solver.momentum]\nrelaxation = 0.9\n[solver.pseudo_time]\nfactor = 20.0\n"
                .to_string(),
        ),
        (
            "relaxation_only",
            "algorithm = \"coupled\"\n[solver.momentum]\nrelaxation = 0.9\n[solver.pseudo_time]\nmomentum = false\n"
                .to_string(),
        ),
    ];
    for (name, solver) in variants {
        let (f, n) = force(name, &solver);
        eprintln!("{name}: {n} iterations (SIMPLEC {n_simplec})");
        for i in 0..2 {
            assert!(
                (f[i] - reference[i]).abs() <= 1e-8 * reference[0].abs(),
                "{name}: lid force {f:?}, SIMPLEC {reference:?}"
            );
        }
    }
}

/// The pseudo-time step sets the number of outer iterations (D177): on a
/// laminar cavity at Re = 400 the automatic step (the volume over the
/// boundary area, 0.25 L/U) is a careful one, and a larger step needs
/// fewer iterations than SIMPLEC at its largest stable relaxation. (Each
/// iteration then costs more in the momentum sweeps of the
/// preconditioner: step 22 of the requirements, section 7.)
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn a_larger_step_needs_fewer_iterations() {
    let iterations = |name: &str, solver: &str| {
        let case = cavity(
            48,
            0.0025,
            &format!("[solver]\ntolerance = 1e-8\nmax_iterations = 5000\n{solver}"),
        );
        let rep = run_case(&case, name);
        assert!(rep.solves[0].converged, "{name} did not converge");
        rep.solves[0].iterations
    };
    let simplec = iterations("count_simplec", "[solver.momentum]\nrelaxation = 0.9\n");
    let coupled = |factor: f64| {
        iterations(
            &format!("count_coupled_{factor}"),
            &format!("algorithm = \"coupled\"\n[solver.pseudo_time]\nfactor = {factor}\n"),
        )
    };
    let (one, four, twenty) = (coupled(1.0), coupled(4.0), coupled(20.0));
    eprintln!(
        "cavity Re 400, 48 x 48: SIMPLEC {simplec}, coupled {one} (factor 1), {four} (4), {twenty} (20) iterations"
    );
    assert!(one < 400, "{one} iterations with the automatic step");
    assert!(
        10 * four <= 7 * simplec && twenty < four,
        "coupled {four} (factor 4), {twenty} (factor 20), SIMPLEC {simplec} iterations"
    );
}

/// D177: the automatic step is L/U with L the volume over the area of the
/// boundaries without symmetry planes (a quarter in the unit cavity) and U
/// the largest prescribed boundary speed (the lid). The symmetry planes of
/// a mesh of one cell layer carry the whole field and must not enter U:
/// the step is still 0.25 s after the flow has developed for some
/// iterations, and 0.125 s at twice the lid speed. The step is compared,
/// not the iteration count of a run: that reacts to the last digit of the
/// step.
#[test]
fn automatic_step_is_the_length_over_the_prescribed_speed() {
    let step = |name: &str, lid: f64, pseudo: &str| {
        let text = cavity_text(
            16,
            0.0025 * lid,
            &format!(
                "[solver]\nalgorithm = \"coupled\"\nmax_iterations = 10\n{pseudo}[verification]\nrequire_convergence = false\n"
            ),
        )
        .replace(
            "velocity = [1.0, 0.0, 0.0]",
            &format!("velocity = [{lid:?}, 0.0, 0.0]"),
        );
        let rep = run_case(&FlowCase::from_toml(&text, Path::new(".")).unwrap(), name);
        rep.solves[0].residuals.pseudo_step
    };
    for (lid, expected) in [(1.0, 0.25), (2.0, 0.125)] {
        let dt = step(&format!("step_automatic_{lid}"), lid, "").expect("a step");
        assert!(
            (dt - expected).abs() < 1e-12,
            "lid speed {lid}: automatic step {dt} s"
        );
    }
    assert_eq!(
        step("step_factor", 1.0, "[solver.pseudo_time]\nfactor = 4.0\n").map(|dt| dt.round()),
        Some(1.0)
    );
    assert_eq!(
        step("step_given", 1.0, "[solver.pseudo_time]\nstep = 0.3\n"),
        Some(0.3)
    );
    assert_eq!(
        step("step_none", 1.0, "[solver.pseudo_time]\nmomentum = false\n"),
        None
    );
}

/// D31: the result does not depend on the number of threads.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn coupled_results_are_bitwise_equal_on_1_and_4_threads() {
    // above the 20 000 cells below which one thread is used
    let mk = |threads: usize| {
        cavity(
            150,
            0.01,
            &format!(
                "[solver]\nalgorithm = \"coupled\"\nmax_iterations = 8\n[run]\nthreads = {threads}\n[verification]\nrequire_convergence = false\n"
            ),
        )
    };
    let a = run_case(&mk(1), "threads1");
    let b = run_case(&mk(4), "threads4");
    let (ra, rb) = (&a.solves[0].residuals, &b.solves[0].residuals);
    assert_eq!(ra.momentum.map(f64::to_bits), rb.momentum.map(f64::to_bits));
    assert_eq!(ra.continuity.to_bits(), rb.continuity.to_bits());
    assert_eq!(ra.coupled_iterations, rb.coupled_iterations);
    assert!(ra.coupled_iterations > 0);
    for (x, y) in a.solves[0].boundaries.iter().zip(&b.solves[0].boundaries) {
        assert_eq!(
            x.force.map(f64::to_bits),
            y.force.map(f64::to_bits),
            "{}",
            x.name
        );
    }
}

/// D42: a coupled run continued from a restart file ends bit-identically
/// to the run without the interruption (the algorithm keeps no state
/// between outer iterations but the fields and the multigrid).
#[cfg(feature = "hdf5")]
#[test]
fn coupled_restart_continues_bit_identically() {
    use foehncraft_flow::restart;
    let solver = "[solver]\nalgorithm = \"coupled\"\nmax_iterations = 24\namg_rebuild = 5\n[verification]\nrequire_convergence = false\n";
    let full = cavity(16, 0.01, &format!("[run]\nrestart_at = [10]\n{solver}"));
    let da = out("restart_full");
    let ra = run(&full, "full", &da, |_| {}).unwrap();
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    let from = da.join("restart_000010.h5");
    let cont = cavity(
        16,
        0.01,
        &format!("[run]\nrestart = '{}'\n{solver}", from.display()),
    );
    let db = out("restart_cont");
    let rb = run(&cont, "cont", &db, |_| {}).unwrap();
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    let b = restart::read(&db.join("restart.h5")).unwrap();
    assert_eq!(a.iteration, 24);
    assert_eq!(b.iteration, 24);
    let (s, t) = (&a.state, &b.state);
    assert!(
        s.u.x == t.u.x && s.u.y == t.u.y && s.u.z == t.u.z,
        "the continued velocity differs"
    );
    assert!(s.p == t.p, "the continued pressure differs");
    assert!(s.mass_flux == t.mass_flux, "the continued fluxes differ");
    assert!(a.amg_matrix == b.amg_matrix, "the multigrid matrix differs");
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// The turbulence equations take the step min(Δτ, c_T T) with the coupled
/// algorithm (D178): the channel of the gallery meets its band in a
/// fraction of the iterations of SIMPLEC (about 2 100).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn coupled_channel_meets_its_band() {
    let text = std::fs::read_to_string(cases().join("channel_zeta_f0.toml"))
        .unwrap()
        .replace("\r\n", "\n")
        .replace(
            "max_iterations = 5000",
            "max_iterations = 5000\nalgorithm = \"coupled\"",
        );
    let r = run_case(&FlowCase::from_toml(&text, &cases()).unwrap(), "channel");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let f = r.solves[0].periodic[0].driving_force;
    let ub_plus = 17.55 / f.sqrt();
    let dev = ub_plus / 17.55 - 1.0;
    eprintln!(
        "channel k-ε-ζ-f₀, coupled: U_b⁺ {ub_plus:.3} ({:+.2} %), {} iterations",
        100.0 * dev,
        r.solves[0].iterations
    );
    assert!(dev.abs() < 0.03, "U_b⁺ {ub_plus}");
    assert!(
        r.solves[0].iterations < 1000,
        "{} iterations",
        r.solves[0].iterations
    );
}

#[test]
fn coupled_settings_are_checked() {
    // steady only
    let (code, msg) = refused(&cavity_text(
        8,
        0.01,
        "[time]\ntime_step = 0.01\nsteps = 2\n[solver]\nalgorithm = \"coupled\"\n",
    ));
    assert_eq!(code, diag::CASE_INCONSISTENT, "{msg}");
    assert!(msg.contains("steady algorithm"), "{msg}");
    // the cap on the local step
    let (code, msg) = refused(&cavity_text(
        8,
        0.01,
        "[solver]\nalgorithm = \"coupled\"\n[solver.momentum]\nrelaxation = 1.0\n",
    ));
    assert_eq!(code, diag::CASE_RANGE, "{msg}");
    assert!(msg.contains("relaxation < 1"), "{msg}");
    for (extra, key) in [
        ("[solver.pseudo_time]\nfactor = 0.0\n", "factor"),
        ("[solver.pseudo_time]\nstep = -1.0\n", "step"),
        ("[solver.pseudo_time]\nturbulence = -0.5\n", "turbulence"),
        (
            "[solver]\nalgorithm = \"coupled\"\n[solver.coupled]\nmax_iter = 500\n",
            "max_iter",
        ),
        (
            "[solver]\nalgorithm = \"coupled\"\n[solver.coupled]\nmomentum_tolerance = 1.0\n",
            "momentum_tolerance",
        ),
    ] {
        let (code, msg) = refused(&cavity_text(8, 0.01, extra));
        assert_eq!(code, diag::CASE_RANGE, "{msg}");
        assert!(msg.contains(key), "{msg}");
    }
    // the defaults: the cap 0.98, c_T 0.5 with the coupled algorithm only
    let c = cavity(8, 0.01, "[solver]\nalgorithm = \"coupled\"\n");
    assert_eq!(c.solver.alpha_u(false), 0.98);
    assert_eq!(
        c.solver
            .pseudo_time
            .turbulence_factor(c.solver.algorithm(false)),
        0.5
    );
    let c = cavity(8, 0.01, "");
    assert_eq!(c.solver.alpha_u(false), 0.7);
    assert_eq!(
        c.solver
            .pseudo_time
            .turbulence_factor(c.solver.algorithm(false)),
        0.0
    );
}

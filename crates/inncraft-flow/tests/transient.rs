//! Stage 5 acceptance without restart files (docs/verification/stage5.md):
//! the steady state of the time stepping and the time order on small
//! meshes. The full sweeps run from the case files in CI.

use inncraft_flow::case::{FlowCase, FlowReport, run};
use std::path::Path;

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir =
        std::env::temp_dir().join(format!("inncraft_transient_{name}_{}", std::process::id()));
    let r = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    r
}

const CAVITY: &str = r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [24, 24, 1]
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
"#;

fn lid_force(r: &FlowReport) -> [f64; 3] {
    let s = r.solves.last().unwrap();
    s.patches.iter().find(|p| p.name == "ymax").unwrap().force
}

/// D27: the Rhie–Chow flux carries the time-step term in the form whose
/// steady state does not depend on Δt, so time stepping to the steady state
/// (BDF2 and Crank–Nicolson, Courant number about 5) ends at the SIMPLEC
/// solution.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn steady_state_of_the_time_stepping_is_the_simplec_solution() {
    let steady = FlowCase::from_toml(
        &format!("{CAVITY}\n[solver]\ntolerance = 1e-12\nmax_iterations = 5000\nmomentum_tolerance = 1e-10\npressure_tolerance = 1e-10\n"),
        Path::new("."),
    )
    .unwrap();
    let rs = run_case(&steady, "simplec");
    assert!(rs.solves[0].converged);
    let f0 = lid_force(&rs);
    for scheme in ["bdf2", "crank-nicolson"] {
        let case = FlowCase::from_toml(
            &format!("{CAVITY}\n[time]\nscheme = \"{scheme}\"\ntime_step = 0.2\nsteps = 400\n"),
            Path::new("."),
        )
        .unwrap();
        let rt = run_case(&case, scheme);
        let f = lid_force(&rt);
        let d = (0..2).map(|i| (f[i] - f0[i]).abs()).fold(0.0, f64::max);
        let scale = f0[0].hypot(f0[1]);
        assert!(d < 1e-9 * scale, "{scheme}: lid force {f:?} against {f0:?}");
    }
}

fn taylor_green(scheme: &str) -> FlowCase {
    let text = format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.015625]
cells = [64, 64, 1]
[fluid]
viscosity = 0.01
[boundary.xmin]
type = "velocity-inlet"
profile = "analytic"
[boundary.xmax]
type = "velocity-inlet"
profile = "analytic"
[boundary.ymin]
type = "velocity-inlet"
profile = "analytic"
[boundary.ymax]
type = "velocity-inlet"
profile = "analytic"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[analytic]
type = "taylor-green"
[time]
scheme = "{scheme}"
time_step = 0.02
end_time = 0.1
[verification]
time_steps = [0.02, 0.01, 0.005, 0.0025]
[output]
vtkhdf = false
wall_csv = false
"#
    );
    FlowCase::from_toml(&text, Path::new(".")).unwrap()
}

/// D28: plain PISO (one outer iteration, two correctors) is second order in
/// time with BDF2 and with Crank–Nicolson. BDF1 is first order; at these steps
/// its observed order is still 1.3–1.4, falling towards 1, so the test only
/// checks that it stays well below 2.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn taylor_green_time_order() {
    for (scheme, lo, hi) in [
        ("bdf2", 1.8, 2.6),
        ("crank-nicolson", 1.8, 2.6),
        ("bdf1", 0.9, 1.6),
    ] {
        let r = run_case(&taylor_green(scheme), scheme);
        let o = *r.time_orders.last().unwrap();
        assert!(
            o > lo && o < hi,
            "{scheme}: time orders {:?}",
            r.time_orders
        );
    }
}

/// The steady cylinder benchmark 2D-1 of Schäfer & Turek on the coarsest
/// built-in mesh (6 144 cells): C_D and Δp inside the reference ranges,
/// C_L (a small difference of large pressure forces) within 10 %; the
/// case file `cylinder_re20.toml` checks the finer mesh to 0.5 %.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn cylinder_re20_on_the_coarse_mesh() {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cases/flow/cylinder_re20.toml"
    ))
    .unwrap()
    .replace("refine = 1", "refine = 0")
    .replace("benchmark_tolerance = 0.005", "benchmark_tolerance = 0.1");
    let case = FlowCase::from_toml(&text, Path::new(".")).unwrap();
    let r = run_case(&case, "cylinder_re20");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let b = r.solves[0].benchmark.as_ref().unwrap();
    assert_eq!(b.values[0].deviation, 0.0, "C_D {}", b.values[0].value);
    assert_eq!(b.values[2].deviation, 0.0, "Δp {}", b.values[2].value);
}

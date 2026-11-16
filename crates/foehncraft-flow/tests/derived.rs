//! Derived fields (D72): exact values of the Couette flow on every mesh
//! type, the vortex criteria of the Taylor–Green vortex against the exact
//! velocity gradient, `foehn post derive` against the fields written by the
//! run, and the checks of the case keys (docs/verification/stage5.md).
#![cfg(feature = "hdf5")]

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::derived::{self, Derived, kinematics};
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_derived_{name}_{}", std::process::id()))
}

/// Every derived field, as a `--set`.
const ALL: &str = r#"output.derived=["vorticity","q-criterion","lambda2","enstrophy","strain-rate","pressure-coefficient","total-pressure","wall-shear-stress","y-plus"]"#;

fn load(file: &str, sets: &[&str]) -> FlowCase {
    let sets: Vec<String> = sets.iter().map(|s| s.to_string()).collect();
    FlowCase::from_file_with(&cases().join(file), &sets).unwrap()
}

/// Plane Couette flow u = U y / H: the velocity gradient is constant and
/// the least-squares gradient exact on any mesh, so every derived field is
/// exact: ω_z = −U/H, Q = λ₂ = 0, the shear rate U/H, the wall shear stress
/// μ U / H, and y⁺ = d √(τ_w/ρ) / ν with the wall distance d of the cells.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn couette_derived_fields_are_exact_on_every_mesh_type() {
    let sets = [
        ALL,
        "output.reference.velocity=1.0",
        "output.reference.pressure=0.5",
        "run.restart_final=true",
    ];
    for t in ["hex", "tet", "poly", "mixed"] {
        let name = format!("couette_{t}");
        let case = load(&format!("{name}.toml"), &sets);
        let dir = out(&name);
        let rep = run(&case, &name, &dir, |_| {}).unwrap();
        assert!(rep.failures.is_empty(), "{t}: {:?}", rep.failures);
        let (mesh, geo, f) =
            derived::evaluate(&case, &dir.join("restart.h5"), &case.output.derived).unwrap();
        let g = |n: &str| f.get(n).unwrap();
        let tol = 1e-8;
        for (n, exact) in [
            ("vorticity_x", 0.0),
            ("vorticity_y", 0.0),
            ("vorticity_z", -1.0),
            ("vorticity_magnitude", 1.0),
            ("q_criterion", 0.0),
            ("lambda2", 0.0),
            ("enstrophy", 0.5),
            ("strain_rate", 1.0),
        ] {
            let e = g(n).iter().fold(0.0f64, |m, v| m.max((v - exact).abs()));
            assert!(e < tol, "{t}: {n} off by {e:e}");
        }
        // total pressure and pressure coefficient from p and u
        let (cp, pt) = (g("pressure_coefficient"), g("total_pressure"));
        let ex = rep.solves[0].errors.as_ref().map(|e| e.l2);
        assert!(ex.is_some_and(|e| e < 1e-9), "{t}: {ex:?}");
        for c in 0..mesh.n_cells {
            let y = geo.cell_centroid[c][1];
            // p ≈ 0 (the exact pressure is constant, the outlet sets 0), u = y
            assert!(
                (pt[c] - 0.5 * y * y).abs() < 1e-8,
                "{t}: total pressure {}",
                pt[c]
            );
            assert!(
                (cp[c] - (0.0 - 0.5) / 0.5).abs() < 1e-8,
                "{t}: cp {}",
                cp[c]
            );
        }
        // walls: τ = μ U / H = 0.1 Pa along ±x, y⁺ from the wall distance
        let ni = mesh.n_interior_faces();
        let (tx, tm, yp) = (
            g("wall_shear_stress_x"),
            g("wall_shear_stress_magnitude"),
            g("y_plus"),
        );
        let mut walls = 0;
        for pa in &mesh.patches {
            for fc in pa.range() {
                let b = fc - ni;
                match pa.name.as_str() {
                    "ymin" | "ymax" => {
                        let sign = if pa.name == "ymin" { 1.0 } else { -1.0 };
                        assert!((tx[b] - sign * 0.1).abs() < 1e-9, "{t}: tau {}", tx[b]);
                        assert!((tm[b] - 0.1).abs() < 1e-9);
                        let c = mesh.owner[fc] as usize;
                        let d = if pa.name == "ymin" {
                            geo.cell_centroid[c][1]
                        } else {
                            1.0 - geo.cell_centroid[c][1]
                        };
                        let exact = d * (0.1f64).sqrt() / 0.1;
                        assert!((yp[b] - exact).abs() < 1e-8 * exact, "{t}: y+ {}", yp[b]);
                        walls += 1;
                    }
                    _ => assert!(tm[b] == 0.0 && yp[b] == 0.0, "{t}: not a wall"),
                }
            }
        }
        assert!(walls > 0);
        let _ = std::fs::remove_dir_all(&dir);
    }
}

/// The Taylor–Green vortex at t = 0.1 s on n × n cells: the vortex
/// criteria against the same formulas applied to the exact velocity
/// gradient at the cell centroids.
fn taylor_green_errors(n: usize) -> [f64; 3] {
    let text = format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [{n}, {n}, 1]
[fluid]
viscosity = 0.01
[flow]
convection = true
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
[output]
vtkhdf = false
wall_csv = false
derived = ["vorticity", "q-criterion", "lambda2"]
[run]
restart_final = true
[time]
scheme = "bdf2"
time_step = 0.0025
end_time = 0.1
"#
    );
    let case = FlowCase::from_toml(&text, Path::new(".")).unwrap();
    let name = format!("taylor_green_{n}");
    let dir = out(&name);
    let rep = run(&case, &name, &dir, |_| {}).unwrap();
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let (mesh, geo, f) =
        derived::evaluate(&case, &dir.join("restart.h5"), &case.output.derived).unwrap();
    let a = case.analytic_at(0.1).unwrap();
    let h = 1e-6;
    let (mut eq, mut el, mut ew, mut sq, mut sw) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let (q, l2, wz) = (
        f.get("q_criterion").unwrap(),
        f.get("lambda2").unwrap(),
        f.get("vorticity_z").unwrap(),
    );
    for c in 0..mesh.n_cells {
        let x = geo.cell_centroid[c];
        // exact gradient by central differences of the exact velocity
        let mut g = [[0.0; 3]; 3];
        for j in 0..2 {
            let (mut xp, mut xm) = (x, x);
            xp[j] += h;
            xm[j] -= h;
            let (up, um) = (a.velocity(xp), a.velocity(xm));
            for i in 0..3 {
                g[i][j] = (up[i] - um[i]) / (2.0 * h);
            }
        }
        let k = kinematics(g);
        eq += (q[c] - k.q).powi(2);
        el += (l2[c] - k.lambda2).powi(2);
        ew += (wz[c] - k.vorticity[2]).powi(2);
        sq += k.q * k.q;
        sw += k.vorticity[2] * k.vorticity[2];
    }
    let _ = std::fs::remove_dir_all(&dir);
    let nc = mesh.n_cells as f64;
    let qs = (sq / nc).sqrt();
    // RMS errors relative to the RMS of the exact Q (and of ω for ω)
    [
        (eq / nc).sqrt() / qs,
        (el / nc).sqrt() / qs,
        (ew / sw).sqrt(),
    ]
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn taylor_green_vortex_criteria_converge() {
    let e: Vec<[f64; 3]> = [16, 32, 64]
        .iter()
        .map(|&n| taylor_green_errors(n))
        .collect();
    for k in 0..3 {
        let o1 = (e[0][k] / e[1][k]).log2();
        let o2 = (e[1][k] / e[2][k]).log2();
        eprintln!(
            "field {k}: errors {:.3e} {:.3e} {:.3e}, orders {o1:.2} {o2:.2}",
            e[0][k], e[1][k], e[2][k]
        );
        assert!(o2 > 1.8, "field {k}: order {o2} ({e:?})");
        assert!(e[2][k] < 0.05, "field {k}: {e:?}");
    }
}

/// `foehn post derive` writes the file the run writes with `[output]
/// derived` at the same iteration or step: byte for byte (D84), steady and
/// transient.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn post_derive_gives_the_file_of_the_run() {
    for (file, extra) in [
        ("couette_poly.toml", vec![]),
        (
            "taylor_green_bdf2.toml",
            vec!["time.end_time=0.05", "verification.time_steps=[]"],
        ),
    ] {
        let mut sets = vec![
            r#"output.derived=["vorticity","q-criterion","lambda2","enstrophy","strain-rate","total-pressure","wall-shear-stress","y-plus"]"#,
            "run.restart_final=true",
            "output.vtkhdf=true",
        ];
        sets.extend(extra);
        let case = load(file, &sets);
        let name = file.trim_end_matches(".toml");
        let dir = out(&format!("post_{name}"));
        let rep = run(&case, name, &dir, |_| {}).unwrap();
        assert!(rep.failures.is_empty(), "{name}: {:?}", rep.failures);
        let post = dir.join("derived.vtkhdf");
        let r = derived::from_restart(&case, &dir.join("restart.h5"), &case.output.derived, &post)
            .unwrap();
        assert_eq!(r.derived.len(), 8);
        let a = std::fs::read(dir.join("fields.vtkhdf")).unwrap();
        let b = std::fs::read(&post).unwrap();
        assert!(a == b, "{name}: the files differ");
        let _ = std::fs::remove_dir_all(&dir);
    }
}

#[test]
fn checks_of_the_derived_keys() {
    let base = std::fs::read_to_string(cases().join("couette_hex.toml")).unwrap();
    let with = |extra: &str| FlowCase::from_toml(&format!("{base}\n{extra}"), &cases());
    // the names of the case file
    let c = with("[output]\nderived = [\"y-plus\", \"vorticity\", \"vorticity\"]").unwrap();
    assert_eq!(
        derived::sorted(&c.output.derived),
        vec![Derived::Vorticity, Derived::YPlus]
    );
    let e = with("[output]\nderived = [\"helicity\"]").unwrap_err();
    assert!(format!("{e:#}").contains("helicity"), "{e:#}");
    // the pressure coefficient needs a reference velocity
    let e = with("[output]\nderived = [\"pressure-coefficient\"]").unwrap_err();
    let d = foehncraft_flow::diag::of_error(&e);
    assert_eq!(d.code, "CASE005", "{e:#}");
    assert!(d.message.contains("reference velocity"), "{}", d.message);
    with("[output]\nderived = [\"pressure-coefficient\"]\n[output.reference]\nvelocity = 2.0")
        .unwrap();
    // or takes the one of the force coefficients
    with(
        "[output]\nderived = [\"pressure-coefficient\"]\n[output.coefficients]\npatches = [\"ymin\"]\nvelocity = 1.0\narea = 1.0",
    )
    .unwrap();
    let e = with("[output.reference]\nvelocity = -1.0").unwrap_err();
    assert_eq!(foehncraft_flow::diag::of_error(&e).code, "CASE004");
    let e = with("[output.reference]\ndensity = 0.0").unwrap_err();
    assert_eq!(foehncraft_flow::diag::of_error(&e).code, "CASE004");
}

//! Stage 6, part 6f: the thermal wall function and the turbulent Prandtl
//! number of a passive scalar (D135), docs/verification/stage6.md.
//!
//! The wall-function formulas are tested against their reference values in
//! `foehncraft-kernels` (`turb::thermal`); here the coupled solver: the
//! heated pipe of the gallery (`pipe_heated_keps.toml`) against the
//! correlation of Gnielinski (1976).

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir =
        std::env::temp_dir().join(format!("foehncraft_thermal_{name}_{}", std::process::id()));
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

/// Nu = q D/(Γ (T_w − T_b)) of the periodic pipe: T_b = 0 at xmin and T_w
/// there is the wall mean minus G L/2.
fn nusselt(r: &FlowReport) -> f64 {
    let sc = &r.solves.last().unwrap().scalars[0];
    let b = |n: &str| sc.boundaries.iter().find(|b| b.boundary == n).unwrap();
    let (q, d, l, u, radius, nu_f, pr) = (1e-3, 1.0, 0.1, 1.0, 0.5, 1e-5, 0.71);
    let g = 2.0 * q / (u * radius);
    let tw = b("ymax").mean - g * l / 2.0;
    let tb = b("xmin").mixed_mean.unwrap();
    q * d / (nu_f / pr * (tw - tb))
}

/// Gnielinski (1976), Pr = 0.71.
fn gnielinski(re: f64, pr: f64) -> f64 {
    let f = (0.79 * re.ln() - 1.64).powi(-2);
    f / 8.0 * (re - 1000.0) * pr / (1.0 + 12.7 * (f / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0))
}

/// The heated pipe at Re = 10⁵ and Pr = 0.71 against Gnielinski with the
/// constant Pr_t and with Kays & Crawford (±10 %), and the wall function
/// itself: without it the molecular conduction across the first cell
/// (y⁺ ≈ 32) gives a Nusselt number far below.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn heated_pipe_meets_gnielinski() {
    let base = text("pipe_heated_keps.toml");
    let want = gnielinski(1e5, 0.71);
    for (name, extra) in [
        ("constant", ""),
        ("kc", "turbulent_prandtl_model = \"kays-crawford\"\n"),
        ("jayatilleke", "wall_function = \"jayatilleke\"\n"),
    ] {
        let mut t = base.clone();
        if name == "jayatilleke" {
            t = t.replace("wall_function = \"kader\"\n", "");
        }
        t = t.replace(
            "turbulent_prandtl = 0.85\n",
            &format!("turbulent_prandtl = 0.85\n{extra}"),
        );
        let r = run_case(&FlowCase::from_toml(&t, &cases()).unwrap(), name);
        assert!(r.failures.is_empty(), "{name}: {:?}", r.failures);
        let nu = nusselt(&r);
        eprintln!(
            "heated pipe {name}: Nu {nu:.2} against {want:.2} ({:+.1} %)",
            100.0 * (nu / want - 1.0)
        );
        assert!(
            (nu / want - 1.0).abs() < 0.10,
            "{name}: Nu {nu} against {want}"
        );
    }
    let t = base.replace("wall_function = \"kader\"", "wall_function = \"none\"");
    let r = run_case(&FlowCase::from_toml(&t, &cases()).unwrap(), "none");
    let nu = nusselt(&r);
    assert!(nu < 0.8 * want, "without the wall function Nu is {nu}");
}

/// The keys are checked.
#[test]
fn bad_thermal_keys_are_refused() {
    let base = text("pipe_heated_keps.toml");
    for (from, to) in [
        ("turbulent_prandtl = 0.85", "turbulent_prandtl = 0.0"),
        ("wall_function = \"kader\"", "wall_function = \"log\""),
    ] {
        let t = base.replace(from, to);
        let bad = match FlowCase::from_toml(&t, &cases()) {
            Err(_) => true,
            Ok(c) => run(
                &c,
                "bad",
                &std::env::temp_dir().join("foehncraft_thermal_bad"),
                |_| {},
            )
            .is_err(),
        };
        assert!(bad, "{to} was accepted");
    }
}

/// The isothermal flat plate (`plate_heated_keps.toml`, Pr = 0.71): the
/// Stanton number of the wall file against St = 0.0287 Re_x^−0.2 Pr^−0.4
/// (Kays, Crawford & Weigand 2005) for 0.25 ≤ x ≤ 1 (Re_x from 1.25·10⁶ to
/// 5·10⁶), ±10 %.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn heated_plate_meets_the_stanton_correlation() {
    let case = FlowCase::from_toml(&text("plate_heated_keps.toml"), &cases()).unwrap();
    let dir = std::env::temp_dir().join(format!("foehncraft_thermal_plate_{}", std::process::id()));
    let r = run(&case, "plate", &dir, |_| {}).unwrap();
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let csv = std::fs::read_to_string(dir.join("wall_plate.csv")).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let mut lines = csv.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split(',').collect();
    let col = |n: &str| head.iter().position(|h| *h == n).unwrap();
    let (cx, cq, cw) = (col("x"), col("theta_flux"), col("theta_wall"));
    let mut n = 0;
    for l in lines {
        let v: Vec<f64> = l.split(',').map(|s| s.parse().unwrap()).collect();
        let x = v[cx];
        if !(0.25..=1.0).contains(&x) {
            continue;
        }
        // flux into the domain is positive for a hot wall; U = 1, θ_∞ = 0
        let st = v[cq] / v[cw];
        let want = 0.0287 * (5e6 * x).powf(-0.2) * 0.71f64.powf(-0.4);
        assert!(
            (st / want - 1.0).abs() < 0.10,
            "x {x}: St {st} against {want}"
        );
        n += 1;
    }
    assert!(n > 5, "{n} wall faces");
}

/// The heated channel of `channel_heated_keps.toml` at Re_τ = 180 … 1020,
/// ν = 1/Re_τ (u_τ = 1 m/s, h = 1 m), Γ = ν/Pr, the mass flow of the DNS
/// bulk Reynolds number and `cells` uniform cells across 2h (first-cell
/// y⁺ ≈ 30), with `model` for Pr_t; returns Nu_2h and the y⁺ of the run.
fn heated_channel(re_tau: f64, pr: f64, re_b: f64, cells: usize, model: &str) -> (f64, f64) {
    let ub = re_b / (2.0 * re_tau);
    let mut t = text("channel_heated_keps.toml");
    for (a, b) in [
        (
            "cells = [2, 12, 1]".to_string(),
            format!("cells = [2, {cells}, 1]"),
        ),
        (
            "viscosity = 0.0025316455696202532   # 1/395".to_string(),
            format!("viscosity = {:?}", 1.0 / re_tau),
        ),
        (
            "diffusivity = 0.0035661218424962853   # ν/Pr, Pr = 0.71".to_string(),
            format!("diffusivity = {:?}", 1.0 / (re_tau * pr)),
        ),
        (
            "velocity = [15.0, 0.0, 0.0]".to_string(),
            format!("velocity = [{ub:?}, 0.0, 0.0]"),
        ),
        (
            "mass_flow = 3.51".to_string(),
            format!("mass_flow = {:?}", 0.2 * ub),
        ),
        (
            "turbulent_prandtl = 0.85\n".to_string(),
            format!("turbulent_prandtl = 0.85\nturbulent_prandtl_model = \"{model}\"\n"),
        ),
    ] {
        assert!(t.contains(&a), "{a}");
        t = t.replace(&a, &b);
    }
    let r = run_case(
        &FlowCase::from_toml(&t, &cases()).unwrap(),
        &format!("channel_{re_tau}_{pr}_{cells}_{model}"),
    );
    assert!(
        r.failures.is_empty(),
        "{re_tau} {pr} {model}: {:?}",
        r.failures
    );
    let sol = r.solves.last().unwrap();
    let sc = &sol.scalars[0];
    let b = |n: &str| sc.boundaries.iter().find(|b| b.boundary == n).unwrap();
    // Δ = Q/ṁ with Q = q × wall area (2 × 0.1 m × 0.1 m)
    let q = 1e-3;
    let delta = q * 0.02 / (0.2 * ub);
    let tw = 0.5 * (b("ymin").mean + b("ymax").mean) - delta / 2.0;
    let tb = b("xmin").mixed_mean.unwrap();
    let nu = 2.0 * q / (1.0 / (re_tau * pr) * (tw - tb));
    // u_τ = √(F h/ρ) from the driving force
    let y_plus = re_tau * sol.periodic[0].driving_force.abs().sqrt() / cells as f64;
    (nu, y_plus)
}

/// The heated channel against the DNS of Abe, Kawamura & Matsuo (JAXA DNS
/// database: Re_τ 180, 395, 640, 1020 at Pr 0.71 and 0.025, uniform wall
/// heat flux on both walls), on the channel width 2h, with the bulk Reynolds
/// number of each DNS. The reference values are the published Nusselt
/// numbers of the database's file headers, except at Re_τ = 1020,
/// Pr = 0.025, whose header (Nu = 278.71) no profile supports: there the
/// Nu of the database's own profiles (`examples/studies/channel_heat_dns.py`
/// fetches them; the profiles are not part of the repository).
///
/// Pr 0.71: both Pr_t models within ±7 % (specification 9.4; measured
/// −2.6 … +1.0 %). Pr 0.025 (liquid metal): Kays & Crawford within ±8 %
/// (+5.9 … +7.2 %), the constant Pr_t = 0.85 too high by 12 … 22 %, since
/// it ignores the rise of Pr_t at small Péclet numbers.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn heated_channel_meets_the_dns_nusselt_number() {
    // (Re_τ, Pr, Re_b, Nu_2h of the DNS, cells across)
    let table = [
        (180.0, 0.71, 5546.0, 18.25, 6),
        (395.0, 0.71, 13849.0, 36.39, 12),
        (640.0, 0.71, 24095.0, 55.57, 20),
        (1020.0, 0.71, 41441.0, 83.00, 32),
        (180.0, 0.025, 5731.0, 5.48, 6),
        (395.0, 0.025, 14147.0, 6.60, 12),
        (640.0, 0.025, 24428.0, 7.89, 20),
        (1020.0, 0.025, 41441.0, 9.81, 32),
    ];
    for (re_tau, pr, re_b, nu_dns, cells) in table {
        let mut nus = Vec::new();
        for model in ["constant", "kays-crawford"] {
            let (nu, y_plus) = heated_channel(re_tau, pr, re_b, cells, model);
            let err = nu / nu_dns - 1.0;
            eprintln!(
                "heated channel Re_τ {re_tau} Pr {pr} {model}: y⁺ {y_plus:.1}, Nu {nu:.2} \
                 against {nu_dns} ({:+.1} %)",
                100.0 * err
            );
            assert!((25.0..40.0).contains(&y_plus), "y⁺ {y_plus}");
            let band = if pr > 0.5 {
                0.07
            } else if model == "constant" {
                f64::INFINITY
            } else {
                0.08
            };
            assert!(
                err.abs() < band,
                "Re_τ {re_tau} Pr {pr} {model}: Nu {nu} against {nu_dns}"
            );
            nus.push(nu);
        }
        if pr < 0.5 {
            // the constant Pr_t over-predicts the turbulent heat flux of a
            // liquid metal; Kays & Crawford corrects most of it
            assert!(nus[0] > 1.1 * nu_dns && nus[1] < nus[0], "Pr {pr}: {nus:?}");
        }
    }
}

/// With the scalable momentum treatment the thermal law takes T⁺ at
/// ỹ⁺ = max(y⁺, 11.06), as the momentum law does: on a fine mesh (first-cell
/// y⁺ ≈ 8, 48 cells) the Nusselt number of the channel at Re_τ = 395,
/// Pr = 0.71 is +10 % from the DNS (+28 % with T⁺ at the true y⁺). The
/// rest is the k-ε model without damping below y⁺ ≈ 30; such meshes take
/// ζ-f or EB with the compound treatment.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn the_thermal_law_follows_the_scalable_floor() {
    let (nu, y_plus) = heated_channel(395.0, 0.71, 13849.0, 48, "kays-crawford");
    eprintln!("heated channel, 48 cells: y⁺ {y_plus:.1}, Nu {nu:.2} against 36.39");
    assert!(y_plus < 11.06, "y⁺ {y_plus}");
    assert!((nu / 36.39 - 1.0).abs() < 0.13, "Nu {nu}");
}

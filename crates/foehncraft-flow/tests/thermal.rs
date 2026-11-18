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

/// The heated channel at Re_τ = 395, Pr = 0.71 (`channel_heated_keps.toml`)
/// against the DNS of Abe, Kawamura & Matsuo (Nu = 36.39 on the channel
/// width 2h, JAXA DNS database, Re_b = 13 849), with the constant Pr_t and
/// with Kays & Crawford (±7 %, specification 9.4). The reference value is a
/// published number; the DNS profiles are not part of the repository
/// (`examples/studies/channel_heat_dns.py` fetches them).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn heated_channel_meets_the_dns_nusselt_number() {
    let base = text("channel_heated_keps.toml");
    let nu_dns = 36.39;
    for (name, extra) in [
        ("constant", ""),
        ("kc", "turbulent_prandtl_model = \"kays-crawford\"\n"),
    ] {
        let t = base.replace(
            "turbulent_prandtl = 0.85\n",
            &format!("turbulent_prandtl = 0.85\n{extra}"),
        );
        let r = run_case(
            &FlowCase::from_toml(&t, &cases()).unwrap(),
            &format!("channel_{name}"),
        );
        assert!(r.failures.is_empty(), "{name}: {:?}", r.failures);
        let sc = &r.solves.last().unwrap().scalars[0];
        let b = |n: &str| sc.boundaries.iter().find(|b| b.boundary == n).unwrap();
        // Δ = Q/ṁ with Q = q × wall area (2 × 0.1 m × 0.1 m), ṁ = 3.51 kg/s
        let (q, delta) = (1e-3, 1e-3 * 0.02 / 3.51);
        let tw = 0.5 * (b("ymin").mean + b("ymax").mean) - delta / 2.0;
        let tb = b("xmin").mixed_mean.unwrap();
        let nu = 2.0 * q / (1.0 / 395.0 / 0.71 * (tw - tb));
        eprintln!(
            "heated channel {name}: Nu {nu:.2} against {nu_dns} ({:+.1} %)",
            100.0 * (nu / nu_dns - 1.0)
        );
        assert!(
            (nu / nu_dns - 1.0).abs() < 0.07,
            "{name}: Nu {nu} against {nu_dns}"
        );
    }
}

//! Statistics of transient runs (D71): the time mean, RMS and covariance
//! of the fields and probes against the exact values of the decaying
//! Taylor–Green vortex, the batch-means verdict of a trend, the checks of
//! the keys.

#![cfg(feature = "hdf5")]

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_stat_{name}_{}", std::process::id()))
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

const NU: f64 = 0.01;
const T0: f64 = 0.1;
const T1: f64 = 0.4;

/// The time moments of f(t) = e^{−a t} over [T0, T1], a = 2ν k²: mean and
/// variance.
fn decay_moments() -> (f64, f64) {
    let k = 2.0 * std::f64::consts::PI;
    let a = 2.0 * NU * k * k;
    let len = T1 - T0;
    let m = ((-a * T0).exp() - (-a * T1).exp()) / (a * len);
    let m2 = ((-2.0 * a * T0).exp() - (-2.0 * a * T1).exp()) / (2.0 * a * len);
    (m, m2 - m * m)
}

fn u0(x: [f64; 2]) -> [f64; 2] {
    let k = 2.0 * std::f64::consts::PI;
    [
        -(k * x[0]).cos() * (k * x[1]).sin(),
        (k * x[0]).sin() * (k * x[1]).cos(),
    ]
}

/// The decaying Taylor–Green vortex (32² cells, Δt 0.005, statistics from
/// t = 0.1 to 0.4): the velocity at every point is u₀(x) f(t), so the time
/// mean is u₀ ⟨f⟩, the RMS |u₀| σ_f and the covariance of u and v
/// u₀ v₀ σ_f². The probe statistics meet these to the spatial error of the
/// mesh, the fields too (in the L2 norm); the decay is a trend, so the batch
/// means are correlated and the window is `too short` (STAT001).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn taylor_green_statistics_meet_the_exact_time_moments() {
    let t = text("taylor_green_bdf2.toml")
        .replace("cells = [64, 64, 1]", "cells = [32, 32, 1]")
        .replace("time_step = 0.01", "time_step = 0.005")
        .replace(
            "[verification]\n",
            "[verification]\nrequire_convergence = false\n",
        )
        .replace(
            "time_steps = [0.01, 0.005, 0.0025, 0.00125, 0.000625]\n",
            "",
        )
        .replace("expected_time_order = [1.8, 2.3]\n", "")
        + "\n[output.statistics]\nstart = 0.1\nbatches = 10\n\
           \n[[output.probes]]\nname = \"a\"\npoint = [0.125, 0.375, 0.0078125]\n\
           \n[[output.probes]]\nname = \"b\"\npoint = [0.3, 0.1, 0.0078125]\n";
    let case = FlowCase::from_toml(&t, &cases()).unwrap();
    let dir = out("tg");
    let rep = run(&case, "tg", &dir, |_| {}).unwrap();
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let s = rep.solves[0].statistics.as_ref().expect("statistics");
    assert!(
        (s.window[0] - T0).abs() < 1e-12 && (s.window[1] - T1).abs() < 1e-9,
        "{:?}",
        s.window
    );
    assert_eq!(s.samples, 60);
    let (mf, vf) = decay_moments();
    for (p, x) in [("a", [0.125, 0.375]), ("b", [0.3, 0.1])] {
        let e = u0(x);
        for (c, ex) in ["u", "v"].iter().zip(e) {
            let q = s
                .series
                .iter()
                .find(|q| q.quantity == format!("probes.{p}.{c}"))
                .unwrap();
            assert!(
                (q.mean - ex * mf).abs() < 0.01 * e[0].hypot(e[1]),
                "{p}.{c}: mean {} against {}",
                q.mean,
                ex * mf
            );
            assert!(
                (q.rms - ex.abs() * vf.sqrt()).abs() < 0.01 * e[0].hypot(e[1]) * vf.sqrt() + 1e-12,
                "{p}.{c}: rms {} against {}",
                q.rms,
                ex.abs() * vf.sqrt()
            );
            println!(
                "{p}.{c}: mean {:.6e} (exact {:.6e}), rms {:.6e} (exact {:.6e}), ±{:.2e}, r1 {:.2}",
                q.mean,
                ex * mf,
                q.rms,
                ex.abs() * vf.sqrt(),
                q.half_width,
                q.batch_correlation
            );
            assert_eq!(q.verdict, "too short", "{p}.{c}: a decay is a trend");
        }
    }
    assert!(
        rep.warnings.iter().any(|w| w.code == "STAT001"),
        "{:?}",
        rep.warnings
    );
    // the fields
    let f = hdf5_metno::File::open(dir.join("statistics.vtkhdf")).unwrap();
    let g = f.group("VTKHDF/internal").unwrap();
    let read = |n: &str| {
        g.dataset(&format!("CellData/{n}"))
            .unwrap()
            .read_raw::<f64>()
            .unwrap()
    };
    let (mu, mv, ru, cuv) = (
        read("mean_velocity_x"),
        read("mean_velocity_y"),
        read("rms_velocity_x"),
        read("covariance_velocity_xy"),
    );
    // the cell centres: the mean of the eight nodes of each hexahedron
    let pts = g.dataset("Points").unwrap().read_2d::<f64>().unwrap();
    let conn = g
        .dataset("Connectivity")
        .unwrap()
        .read_raw::<i64>()
        .unwrap();
    let centre = |c: usize| {
        let mut x = [0.0; 2];
        for &n in &conn[8 * c..8 * c + 8] {
            x[0] += pts[[n as usize, 0]] / 8.0;
            x[1] += pts[[n as usize, 1]] / 8.0;
        }
        x
    };
    let (mut e2, mut n2, mut r2, mut rn, mut c2, mut cn) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for c in 0..mu.len() {
        let e = u0(centre(c));
        e2 += (mu[c] - e[0] * mf).powi(2) + (mv[c] - e[1] * mf).powi(2);
        n2 += (e[0] * mf).powi(2) + (e[1] * mf).powi(2);
        r2 += (ru[c] - e[0].abs() * vf.sqrt()).powi(2);
        rn += e[0] * e[0] * vf;
        c2 += (cuv[c] - e[0] * e[1] * vf).powi(2);
        cn += (e[0] * e[1] * vf).powi(2);
    }
    let (em, er, ec) = ((e2 / n2).sqrt(), (r2 / rn).sqrt(), (c2 / cn).sqrt());
    println!("L2 relative: mean {em:.3e}, rms {er:.3e}, covariance {ec:.3e}");
    assert!(em < 0.01, "L2 of the mean {em}");
    assert!(er < 0.02, "L2 of the RMS {er}");
    assert!(ec < 0.03, "L2 of the covariance {ec}");
    let csv = std::fs::read_to_string(dir.join("statistics.csv")).unwrap();
    assert!(csv.contains("probes.a.u,"), "{csv}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// The section needs a transient run and at least two batches.
#[test]
fn statistics_keys_are_checked() {
    let code = |t: &str| {
        let e = FlowCase::from_toml(t, &cases())
            .and_then(|c| c.validate().map(|_| c))
            .unwrap_err();
        diag::of_error(&e).code.to_string()
    };
    let steady = text("couette_poly.toml") + "\n[output.statistics]\n";
    assert_eq!(code(&steady), "CASE005");
    let tg = text("taylor_green_bdf2.toml");
    assert_eq!(
        code(&(tg.clone() + "\n[output.statistics]\nbatches = 1\n")),
        "CASE004"
    );
    assert_eq!(
        code(&(tg + "\n[output.statistics]\nphase_period = 0.0\n")),
        "CASE004"
    );
}

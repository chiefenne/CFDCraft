//! Stage 6, part 6g: the derived fields of the turbulence models (D72,
//! specification 2.3): `turbulence-scales`, `reynolds-stress`,
//! `wall-normal` and `turbulence-wall`, checked against the turbulence
//! fields of the restart file on small channels of every model family.

#![cfg(feature = "hdf5")]
#![allow(clippy::needless_range_loop)]

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::derived::{self, Derived};
use foehncraft_flow::{diag, restart};
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_tderived_{name}_{}", std::process::id()))
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

/// The channel `file` on 24 cells across, `n` iterations, with the final
/// restart file.
fn channel(file: &str, n: usize) -> FlowCase {
    let t = text(file)
        .replace("cells = [2, 96, 1]", "cells = [2, 24, 1]")
        .replace("ratio = 60.0", "ratio = 8.0")
        .replace("max_iterations = 3000", &format!("max_iterations = {n}"))
        .replace("max_iterations = 5000", &format!("max_iterations = {n}"))
        + "\n[verification]\nrequire_convergence = false\n[run]\nrestart_final = true\n";
    FlowCase::from_toml(&t, &cases()).unwrap()
}

fn get<'a>(f: &'a derived::Fields, name: &str) -> &'a [f64] {
    f.get(name).unwrap_or_else(|| panic!("no dataset {name}"))
}

/// Every model family: ε and ω satisfy ε = C_μ k ω, T = k/ε (SST, k-ε) or
/// the bounded T of the model (EB), L = k^{3/2}/ε, ν_t/ν; the Reynolds
/// stress has the trace 2k (the flow is free of divergence) and the shear
/// component −ν_t du/dy, negative below the centre line and positive above;
/// on the walls u_τ is that of the wall treatment and y* = y C_μ^¼ k_P^½/ν;
/// EB: n = ∇α/|∇α| is the unit normal pointing away from the nearer wall.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn turbulence_fields_of_every_model() {
    for (file, c_mu, eb) in [
        ("channel_sst.toml", 0.09, false),
        ("channel_keps.toml", 0.09, false),
        ("channel_eb.toml", 0.22, true),
    ] {
        let case = channel(file, 1500);
        let name = file.trim_end_matches(".toml");
        let dir = out(name);
        let rep = run(&case, name, &dir, |_| {}).unwrap();
        assert!(rep.failures.is_empty(), "{name}: {:?}", rep.failures);
        let mut list = vec![
            Derived::TurbulenceScales,
            Derived::ReynoldsStress,
            Derived::TurbulenceWall,
        ];
        if eb {
            list.push(Derived::WallNormal);
        }
        let file = dir.join("restart.h5");
        let (mesh, geo, f) = derived::evaluate(&case, &file, &list).unwrap();
        let nc = mesh.n_cells;
        let ni = mesh.n_interior_faces();
        let r = restart::read(&file).unwrap();
        let (r, _) = r.to_local(&mesh.cell_id[..nc], &mesh.face_id, ni).unwrap();
        let t = r.state.turbulence.unwrap();
        let c_mu = if eb {
            case.turbulence.c_mu().unwrap()
        } else {
            c_mu
        };
        let nu = case.mu() / case.rho();
        let rel = |a: f64, b: f64| (a - b).abs() <= 1e-10 * (a.abs() + b.abs()) + 1e-300;
        let (eps, om, ts, ls, vr) = (
            get(&f, "dissipation_rate"),
            get(&f, "specific_dissipation_rate"),
            get(&f, "time_scale"),
            get(&f, "length_scale"),
            get(&f, "viscosity_ratio"),
        );
        for c in 0..nc {
            let k = t.k[c];
            assert!(
                rel(eps[c], c_mu * k * om[c]),
                "{name} cell {c}: ε = C_μ k ω"
            );
            assert!(rel(ls[c], k.powf(1.5) / eps[c]), "{name} cell {c}: L");
            assert!(rel(vr[c], t.nu_t[c] / nu), "{name} cell {c}: ν_t/ν");
            match &t.zeta {
                Some(z) => assert_eq!(ts[c], z.t[c], "{name} cell {c}: T of the model"),
                None => assert!(rel(ts[c], k / eps[c]), "{name} cell {c}: T = k/ε"),
            }
        }
        let (xx, yy, zz, xy) = (
            get(&f, "reynolds_stress_xx"),
            get(&f, "reynolds_stress_yy"),
            get(&f, "reynolds_stress_zz"),
            get(&f, "reynolds_stress_xy"),
        );
        let (yz, xz) = (get(&f, "reynolds_stress_yz"), get(&f, "reynolds_stress_xz"));
        for c in 0..nc {
            let k = t.k[c];
            let tr = xx[c] + yy[c] + zz[c];
            assert!(
                (tr - 2.0 * k).abs() <= 1e-6 * k + 1e-12,
                "{name} cell {c}: trace {tr} against 2k {}",
                2.0 * k
            );
            let y = geo.cell_centroid[c][1];
            if (y - 1.0).abs() > 0.05 {
                assert!(
                    xy[c] * (y - 1.0) > 0.0,
                    "{name} cell {c}: sign of u'v' at y {y}"
                );
            }
            assert!(
                yz[c].abs() <= 1e-9 * k && xz[c].abs() <= 1e-9 * k,
                "{name} cell {c}"
            );
        }
        let (ut, ys) = (get(&f, "friction_velocity"), get(&f, "y_star"));
        let mut walls = 0;
        for p in &mesh.boundaries {
            let wall = p.name == "ymin" || p.name == "ymax";
            for face in p.range() {
                let fb = face - ni;
                if !wall {
                    assert_eq!((ut[fb], ys[fb]), (0.0, 0.0), "{name}: off the walls");
                    continue;
                }
                walls += 1;
                assert_eq!(ut[fb], t.u_tau[fb], "{name}: u_τ of the wall treatment");
                let c = mesh.owner[face] as usize;
                let y = (geo.cell_centroid[c][1] - geo.face_centroid[face][1]).abs();
                assert!(
                    rel(ys[fb], y * c_mu.sqrt().sqrt() * t.k[c].max(0.0).sqrt() / nu),
                    "{name}: y*"
                );
            }
        }
        assert!(walls > 0);
        if eb {
            let n = [
                get(&f, "wall_normal_x"),
                get(&f, "wall_normal_y"),
                get(&f, "wall_normal_z"),
            ];
            for c in 0..nc {
                let y = geo.cell_centroid[c][1];
                let m = (n[0][c].powi(2) + n[1][c].powi(2) + n[2][c].powi(2)).sqrt();
                assert!((m - 1.0).abs() < 1e-12, "{name} cell {c}: |n| = {m}");
                if (y - 1.0).abs() > 0.05 {
                    assert!(
                        n[1][c] * (1.0 - y).signum() > 0.999,
                        "{name} cell {c}: n at y {y}: {:?}",
                        [n[0][c], n[1][c], n[2][c]]
                    );
                }
            }
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}

/// The turbulence fields need a turbulence model, the wall normal the
/// elliptic-blending k-ε.
#[test]
fn turbulence_fields_are_refused_without_their_model() {
    let code = |text: &str| {
        let e = FlowCase::from_toml(text, &cases())
            .and_then(|c| derived::validate(&c).map(|_| c))
            .unwrap_err();
        (diag::of_error(&e).code.to_string(), format!("{e:#}"))
    };
    let sst = text("channel_sst.toml");
    let (c, m) = code(&sst.replace(
        "wall_csv = true",
        "wall_csv = true\nderived = [\"wall-normal\"]",
    ));
    assert_eq!(c, "CASE005", "{m}");
    assert!(m.contains("elliptic-blending"), "{m}");
    let laminar = text("couette_poly.toml");
    for d in ["turbulence-scales", "reynolds-stress", "turbulence-wall"] {
        let t = if laminar.contains("[output]") {
            laminar.replace("[output]", &format!("[output]\nderived = [\"{d}\"]"))
        } else {
            format!("{laminar}\n[output]\nderived = [\"{d}\"]\n")
        };
        let (c, m) = code(&t);
        assert_eq!(c, "CASE005", "{d}: {m}");
        assert!(m.contains("needs a turbulence model"), "{d}: {m}");
    }
}

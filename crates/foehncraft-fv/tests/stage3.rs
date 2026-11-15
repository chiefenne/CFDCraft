//! Acceptance checks of stage 3 (docs/requirements.md, stage plan row 3):
//! the pressure solver, CG with the smoothed-aggregation multigrid (D35) and
//! its mixed-precision hierarchy (D37).
//!
//! The mesh sweeps run in release builds (CI runs the tests in debug and in
//! release); the agreement with hypre (D20) needs the feature `hypre`.

use foehncraft_fv::amg::{AmgSettings, Precision, Smoother};
use foehncraft_fv::backend::{Backend, V, cpu::CpuBackend};
use foehncraft_fv::poisson::{assemble, solve_inhouse};
use foehncraft_fv::solver::{LinearSolver, Preconditioner, Settings, precondition, setup_amg};
use foehncraft_fv::verify::*;
use foehncraft_fv::{Discretisation, Problem};
use foehncraft_mesh::generate::Kind;
use std::path::Path;

const KINDS: [Kind; 4] = [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed];

fn poisson_problem() -> Problem {
    let mut p = mms_problem(3);
    p.velocity = [0.0; 3];
    p.diffusivity = 1.0;
    p
}

/// A backend with the Poisson system assembled.
fn system(kind: Kind, cells: usize, refine: u32, threads: Option<usize>) -> CpuBackend {
    let (m, g) = test_mesh(kind, 3, cells, refine, 0.15).unwrap();
    let d = Discretisation::new(&m, &g, &poisson_problem()).unwrap();
    let mut be = CpuBackend::new(d, threads).unwrap();
    assemble(&mut be).unwrap();
    be
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// The CG iterations of the multigrid stay nearly constant over a factor of
/// 64 in the number of cells on every mesh type (stage 3 acceptance), while
/// Jacobi-preconditioned CG grows with the mesh.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn amg_iterations_nearly_independent_of_mesh_size() {
    for kind in KINDS {
        let mut its = Vec::new();
        let mut jac = Vec::new();
        for refine in [1, 3] {
            let mut be = system(kind, 4, refine, None);
            let s = AmgSettings::default();
            let (it, ..) = solve_inhouse(&mut be, Preconditioner::Amg, &s, 1e-10, 500).unwrap();
            its.push(it);
            let (it, ..) = solve_inhouse(&mut be, Preconditioner::Jacobi, &s, 1e-10, 5000).unwrap();
            jac.push(it);
        }
        let ratio = its[1] as f64 / its[0] as f64;
        assert!(
            its[1] <= 16 && ratio <= 1.5,
            "{}: AMG iterations {its:?}",
            kind.name()
        );
        assert!(
            jac[1] as f64 > 2.5 * jac[0] as f64,
            "{}: Jacobi {jac:?}",
            kind.name()
        );
    }
}

/// Mixed precision (float32 hierarchy) converges to the same solution as the
/// float64 hierarchy, in the same number of iterations (±1), and the
/// preconditioner itself differs only at float32 rounding.
#[test]
fn mixed_precision_gives_the_float64_result() {
    for kind in KINDS {
        let mut be = system(kind, 4, 1, None);
        let mut sol = Vec::new();
        let mut its = Vec::new();
        let mut z = Vec::new();
        for precision in [Precision::Mixed, Precision::Double] {
            let s = AmgSettings {
                precision,
                ..AmgSettings::default()
            };
            let (it, ..) = solve_inhouse(&mut be, Preconditioner::Amg, &s, 1e-11, 200).unwrap();
            its.push(it as i64);
            sol.push(be.download(V::Phi).unwrap());
            be.copy(V::Rhs, V::R).unwrap();
            precondition(&mut be, Preconditioner::Amg, V::R, V::Z).unwrap();
            z.push(be.download(V::Z).unwrap());
        }
        assert!((its[0] - its[1]).abs() <= 1, "{}: {its:?}", kind.name());
        let scale = sol[1].iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let d = sol[0]
            .iter()
            .zip(&sol[1])
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        assert!(
            d <= 1e-8 * scale,
            "{}: solutions differ by {d:e}",
            kind.name()
        );
        let zs = z[1].iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let dz = z[0]
            .iter()
            .zip(&z[1])
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        assert!(
            dz > 0.0 && dz <= 1e-5 * zs,
            "{}: V-cycles differ by {dz:e} of {zs:e}",
            kind.name()
        );
    }
}

/// The float64 V-cycle is a symmetric positive definite operator, with both
/// smoothers: (M⁻¹x, y) = (x, M⁻¹y) and (M⁻¹x, x) > 0.
#[test]
fn vcycle_is_symmetric_positive_definite() {
    for smoother in [Smoother::Chebyshev, Smoother::Jacobi] {
        let mut be = system(Kind::Mixed, 4, 1, None);
        let s = AmgSettings {
            precision: Precision::Double,
            smoother,
            ..AmgSettings::default()
        };
        setup_amg(&mut be, &s).unwrap();
        let n = be.n_cells();
        let x: Vec<f64> = (0..n)
            .map(|i| foehncraft_core::rng::uniform([i as u64, 1, 2]) - 0.5)
            .collect();
        let y: Vec<f64> = (0..n)
            .map(|i| foehncraft_core::rng::uniform([i as u64, 3, 4]) - 0.5)
            .collect();
        let mut apply = |v: &[f64]| {
            be.upload(V::R, v).unwrap();
            precondition(&mut be, Preconditioner::Amg, V::R, V::Z).unwrap();
            be.download(V::Z).unwrap()
        };
        let (mx, my) = (apply(&x), apply(&y));
        let (a, b) = (dot(&mx, &y), dot(&x, &my));
        assert!(
            (a - b).abs() <= 1e-12 * a.abs().max(b.abs()),
            "{smoother:?}: {a} vs {b}"
        );
        assert!(dot(&mx, &x) > 0.0 && dot(&my, &y) > 0.0);
    }
}

/// The hierarchy setup and the solve give bitwise the same solution on 1
/// and 4 threads (D31).
#[test]
fn amg_solve_is_bitwise_reproducible() {
    for kind in KINDS {
        let mut reference: Option<Vec<u64>> = None;
        for threads in [1, 4] {
            let mut be = system(kind, 4, 1, Some(threads));
            solve_inhouse(
                &mut be,
                Preconditioner::Amg,
                &AmgSettings::default(),
                1e-10,
                200,
            )
            .unwrap();
            let bits: Vec<u64> = be
                .download(V::Phi)
                .unwrap()
                .iter()
                .map(|v| v.to_bits())
                .collect();
            match &reference {
                None => reference = Some(bits),
                Some(r) => assert!(r == &bits, "{}: {threads} threads differ", kind.name()),
            }
        }
    }
}

/// CG with the multigrid inside the outer deferred-correction iterations of
/// the scalar solver (with the non-orthogonal correction) reaches the same
/// solution as BiCGStab with Jacobi, so the discretisation error is the
/// same.
#[test]
fn cg_amg_in_the_scalar_solver_matches_bicgstab() {
    let (m, g) = test_mesh(Kind::Poly, 3, 4, 1, 0.15).unwrap();
    let p = poisson_problem();
    let mut errs = Vec::new();
    for (ls, pc) in [
        (LinearSolver::Bicgstab, Preconditioner::Jacobi),
        (LinearSolver::Cg, Preconditioner::Amg),
        (LinearSolver::Bicgstab, Preconditioner::Amg),
    ] {
        let s = Settings {
            linear_solver: ls,
            preconditioner: pc,
            linear_tolerance: 1e-6,
            tolerance: 1e-11,
            ..Settings::default()
        };
        let (r, _) = manufactured("poly", &m, &g, 3, &p, &s, "cpu", None, Path::new("")).unwrap();
        errs.push(r.l2);
    }
    for e in &errs[1..] {
        assert!((e - errs[0]).abs() <= 1e-9 * errs[0], "{errs:?}");
    }
}

/// The in-house solution agrees with hypre's BoomerAMG-PCG (D20) to the
/// solver tolerance on every mesh type.
#[cfg(feature = "hypre")]
#[test]
fn same_converged_result_as_hypre() {
    use foehncraft_fv::hypre;
    for kind in KINDS {
        let mut be = system(kind, 4, 2, None);
        let val = be.matrix_values().unwrap();
        let a = be.discretisation().unwrap().matrix_csr(&val);
        let b = be.download(V::Rhs).unwrap();
        let h = hypre::solve(&a, &b, 1e-11, 500, &Default::default()).unwrap();
        assert!(h.iterations > 0 && h.iterations < 100, "{}", h.iterations);
        solve_inhouse(
            &mut be,
            Preconditioner::Amg,
            &AmgSettings::default(),
            1e-11,
            200,
        )
        .unwrap();
        let x = be.download(V::Phi).unwrap();
        let scale = h.x.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let d = x
            .iter()
            .zip(&h.x)
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        assert!(
            d <= 1e-8 * scale,
            "{}: differs from hypre by {d:e}",
            kind.name()
        );
    }
}

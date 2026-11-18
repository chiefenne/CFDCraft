//! The reduction multigrid of transport matrices (D179,
//! docs/verification/stage6.md, part 6f) on a backend: the hierarchy of a
//! convection–diffusion matrix of the finite-volume discretisation runs on
//! the level operations of the aggregation multigrid, in both precisions.
//! (The setup itself is tested in `amg::reduction`.)

use foehncraft_fv::Discretisation;
use foehncraft_fv::amg::{Precision, ReductionSettings};
use foehncraft_fv::backend::{Backend, FaceSettings, V, cpu::CpuBackend};
use foehncraft_fv::solver::{cycles, setup_reduction};
use foehncraft_fv::verify::*;
use foehncraft_mesh::generate::Kind;

/// A backend with the upwind convection–diffusion system of the
/// manufactured problem assembled, at a cell Péclet number of about
/// `peclet` on the hexahedral mesh.
fn system(kind: Kind, refine: u32, peclet: f64) -> CpuBackend {
    let (m, g) = test_mesh(kind, 3, 4, refine, 0.15).unwrap();
    let mut p = mms_problem(3);
    let speed = p.velocity.iter().map(|v| v * v).sum::<f64>().sqrt();
    assert!(speed > 0.0, "the manufactured problem has a velocity");
    let h = 1.0 / (4 << refine) as f64;
    p.diffusivity = speed * h / peclet;
    let d = Discretisation::new(&m, &g, &p).unwrap();
    let mut be = CpuBackend::new(d, None).unwrap();
    let n = be.n_cells();
    be.upload(V::Phi, &vec![0.0; n]).unwrap();
    be.gradient().unwrap();
    be.faces(FaceSettings {
        blend: 0.0,
        central: 0.0,
        beta: 0.0,
        nonorth: false,
        nonorth_limit: -1.0,
    })
    .unwrap();
    be.assemble(1.0).unwrap();
    be
}

/// ‖b − A x‖ / ‖b‖ after `n` stationary cycles from x = 0.
fn reduction(be: &mut CpuBackend, n: usize) -> f64 {
    let cells = be.n_cells();
    let b: Vec<f64> = (0..cells)
        .map(|i| foehncraft_core::rng::uniform([i as u64, 17, 4]) - 0.5)
        .collect();
    be.upload(V::Rhs, &b).unwrap();
    cycles(be, n).unwrap();
    be.residual(V::Phi, V::Rhs, V::R).unwrap();
    let r = be.dot(V::R, V::R).unwrap().sqrt();
    let b = be.dot(V::Rhs, V::Rhs).unwrap().sqrt();
    r / b
}

/// The cycle converges on every mesh type, from diffusion-dominated to
/// convection-dominated matrices, and the rate holds on the finer mesh.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn cycles_converge_on_every_mesh_type() {
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        for (refine, peclet) in [(2, 0.5), (2, 20.0), (3, 20.0)] {
            let mut be = system(kind, refine, peclet);
            let info = setup_reduction(&mut be, &ReductionSettings::default()).unwrap();
            let red = reduction(&mut be, 8);
            eprintln!(
                "{} refine {refine} Pe {peclet}: {} cells, levels {:?}, operator complexity \
                 {:.2}, factor per cycle {:.2}",
                kind.name(),
                be.n_cells(),
                info.rows.len(),
                info.operator_complexity,
                red.powf(0.125)
            );
            assert!(
                red < 1e-3,
                "{} refine {refine} Pe {peclet}: residual {red:.2e} after 8 cycles",
                kind.name()
            );
            assert!(info.operator_complexity < 5.5, "{info:?}");
        }
    }
}

/// The hierarchy in float32 (the default) converges like the one in
/// float64 down to the rounding of float32; one cycle is a fixed linear
/// operator (the same bits for the same right-hand side); a second slot
/// leaves the first untouched.
#[test]
fn precision_and_slots() {
    let mut be = system(Kind::Poly, 1, 10.0);
    be.amg_slot(5).unwrap();
    setup_reduction(&mut be, &ReductionSettings::default()).unwrap();
    let mixed = reduction(&mut be, 6);
    let z1 = be.download(V::Phi).unwrap();
    let again = reduction(&mut be, 6);
    assert_eq!(mixed.to_bits(), again.to_bits());
    be.amg_slot(6).unwrap();
    let double = ReductionSettings {
        precision: Precision::Double,
        ..ReductionSettings::default()
    };
    setup_reduction(&mut be, &double).unwrap();
    let exact = reduction(&mut be, 6);
    eprintln!("after 6 cycles: float32 hierarchy {mixed:.2e}, float64 {exact:.2e}");
    assert!(exact < 1e-3 && mixed < 1e-3);
    assert!(
        mixed < 30.0 * exact + 1e-6,
        "float32 {mixed:.2e}, float64 {exact:.2e}"
    );
    be.amg_slot(5).unwrap();
    reduction(&mut be, 6);
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(&z1), bits(&be.download(V::Phi).unwrap()));
}

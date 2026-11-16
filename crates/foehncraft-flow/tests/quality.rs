//! The cell-quality remediation (D50): a pipe in raw Gmsh tetrahedra
//! (`examples/meshes/pipe_gmsh_raw.geo`, faces with a skewness up to 1.23)
//! diverges with the least-squares pressure gradient; with
//! `solver.remediation = "auto"` the solver switches to the Gauss pressure
//! gradient and the Venkatakrishnan limiter, reports it and converges.

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn out(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("foehncraft_quality_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

/// Developing laminar pipe flow at Re = 100 on the raw Gmsh mesh.
fn pipe(remediation: &str) -> FlowCase {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/meshes");
    FlowCase::from_toml(
        &format!(
            r#"
[mesh]
file = "pipe_gmsh_raw.msh"
renumber = "rcm"
[fluid]
viscosity = 1e-3
[boundary.inlet_main]
type = "velocity-inlet"
normal_velocity = 1.0
[boundary.outlet]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.wall]
type = "wall"
[solver]
tolerance = 1e-6
max_iterations = 400
remediation = "{remediation}"
[output]
vtkhdf = false
wall_csv = false
"#
        ),
        &dir,
    )
    .unwrap()
}

#[test]
fn raw_gmsh_tetrahedra_converge_with_the_remediation() {
    let d = out("auto");
    let r = run(&pipe("auto"), "pipe", &d, |_| {}).unwrap();
    let s = &r.solves[0];
    assert!(s.converged, "{} iterations", s.iterations);
    let q = s.remediation.as_ref().expect("a remediation report");
    assert!(q.skewed_faces > 0 && q.poor_cells > 0, "{q:?}");
    assert!(q.max_skewness > 0.9, "{q:?}");
    assert_eq!(
        q.changes,
        [
            "solver.pressure_gradient: \"least-squares\" -> \"gauss\"",
            "solver.limiter: \"none\" -> \"venkatakrishnan\""
        ]
    );
    let codes: Vec<&str> = r.warnings.iter().map(|w| w.code).collect();
    assert!(codes.contains(&diag::MESH_SKEWED), "{codes:?}");
    assert!(codes.contains(&diag::MESH_REMEDIATION), "{codes:?}");
    // the report says what was changed and how to avoid it
    let md = std::fs::read_to_string(d.join("flow.md")).unwrap();
    assert!(md.contains("Cell quality (D50)") && md.contains("Mesh.OptimizeNetgen"));
    // the mass balance of the converged state
    assert!(s.mass_imbalance < 1e-6, "{}", s.mass_imbalance);
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn raw_gmsh_tetrahedra_fail_without_the_remediation() {
    let d = out("off");
    let e = match run(&pipe("off"), "pipe", &d, |_| {}) {
        Ok(r) => panic!(
            "the least-squares pressure gradient ran through (converged: {})",
            r.solves.iter().all(|s| s.converged)
        ),
        Err(e) => e,
    };
    // a divergence (the pressure-correction matrix stops being positive
    // definite before a value is non-finite), not a missing pressure level
    assert_eq!(diag::of_error(&e).code, diag::RUN_DIVERGED, "{e:#}");
    // the report names the skewed faces first and the fix
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(d.join("divergence.json")).unwrap()).unwrap();
    let first = &j["causes"][0];
    assert_eq!(first["code"], "DIV003", "{j:#}");
    assert!(
        first["fix"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f == "solver.remediation=\"auto\""),
        "{first:#}"
    );
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn the_remediation_reports_poor_faces_but_changes_nothing_when_off() {
    let d = out("off_report");
    let mut c = pipe("off");
    c.solver.max_iterations = 1;
    let r = run(&c, "pipe", &d, |_| {}).unwrap();
    let q = r.solves[0].remediation.as_ref().unwrap();
    assert!(q.skewed_faces > 0 && q.changes.is_empty(), "{q:?}");
    let codes: Vec<&str> = r.warnings.iter().map(|w| w.code).collect();
    assert!(codes.contains(&diag::MESH_SKEWED) && !codes.contains(&diag::MESH_REMEDIATION));
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn the_gauss_weights_give_the_exact_gradient_of_a_linear_field_without_skewness() {
    use foehncraft_fv::disc::Discretisation;
    use foehncraft_fv::problem::{Boundary, Problem};
    use foehncraft_mesh::Geometry;
    use foehncraft_mesh::generate::{BoxSpec, Grading, Kind};
    let mut spec = BoxSpec {
        kind: Kind::Hex,
        upper: [2.0, 1.0, 0.5],
        cells: [7, 5, 3],
        ..Default::default()
    };
    spec.grading.x = Grading {
        ratio: 3.0,
        symmetric: false,
    };
    spec.grading.y = Grading {
        ratio: 4.0,
        symmetric: true,
    };
    let m = spec.generate().unwrap();
    let g = Geometry::compute(&m);
    let mut prob = Problem::default();
    for p in &m.patches {
        prob.boundary.insert(p.name.clone(), Boundary::value(0.0));
    }
    let d = Discretisation::new(&m, &g, &prob).unwrap();
    let w = foehncraft_flow::simplec::gauss_weights(&d);
    let a = [0.7, -1.3, 2.1];
    let phi = |x: [f64; 3]| 0.4 + a[0] * x[0] + a[1] * x[1] + a[2] * x[2];
    let ni = d.n_interior;
    for c in 0..d.n_cells {
        let pc = phi(g.cell_centroid[c]);
        let mut grad = [0.0; 3];
        for k in d.cell_ptr[c] as usize..d.cell_ptr[c + 1] as usize {
            let (f, is_nb) = foehncraft_kernels::unpack_face(d.cell_faces[k]);
            let other = if f >= ni {
                phi(g.face_centroid[f])
            } else if is_nb {
                phi(g.cell_centroid[d.owner[f] as usize])
            } else {
                phi(g.cell_centroid[d.neighbour[f] as usize])
            };
            for (i, gi) in grad.iter_mut().enumerate() {
                *gi += [w.x[k], w.y[k], w.z[k]][i] * (other - pc);
            }
        }
        for i in 0..3 {
            assert!((grad[i] - a[i]).abs() < 1e-12, "cell {c}: {grad:?}");
        }
    }
    // and no face of this mesh is skewed
    let skew = foehncraft_flow::simplec::face_skewness(&d);
    assert!(skew.iter().all(|&s| s < 1e-12));
}

#[test]
fn the_remediation_changes_only_what_the_case_does_not_set() {
    use foehncraft_flow::simplec::{PressureGradient, Remediation, Settings};
    use foehncraft_fv::solver::{LimiterKind, Scheme};
    let s = Settings::default();
    // a mesh without poor faces, or the remediation off: nothing changes
    assert!(s.remedied(0.9).1.is_empty());
    let off = Settings {
        remediation: Remediation::Off,
        ..Settings::default()
    };
    assert!(off.remedied(1.3).1.is_empty());
    // poor faces: the Gauss pressure gradient and the limiter
    let (r, ch) = s.remedied(1.3);
    assert_eq!(ch.len(), 2);
    assert_eq!(r.pressure_gradient, PressureGradient::Gauss);
    assert_eq!(r.limiter, LimiterKind::Venkatakrishnan);
    // first-order upwind needs no limiter; a limiter of the case stays
    let up = Settings {
        scheme: Scheme::Upwind,
        ..Settings::default()
    };
    assert_eq!(up.remedied(1.3).1.len(), 1);
    let wang = Settings {
        limiter: LimiterKind::VenkatakrishnanWang,
        ..Settings::default()
    };
    let (r, ch) = wang.remedied(1.3);
    assert_eq!(ch.len(), 1);
    assert_eq!(r.limiter, LimiterKind::VenkatakrishnanWang);
    // a stricter limit treats milder faces
    let strict = Settings {
        skewness_limit: 0.5,
        ..Settings::default()
    };
    assert_eq!(strict.remedied(0.6).1.len(), 2);
}

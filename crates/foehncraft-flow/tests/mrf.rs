//! The moving reference frame (D55, stage 6 part 5 of 6g): laminar
//! Taylor–Couette flow against the analytic solution, in any frame.
//!
//! The flow between a turning inner cylinder and a fixed outer one is
//! steady in the inertial frame and in every frame that turns about the
//! same axis, so a rotating zone (`[[mrf]]`) must give the same solution
//! whatever its angular velocity, up to the discretisation error. A fluid
//! in solid-body rotation is at rest in the frame that turns with it: the
//! relative flux is zero, the convective terms vanish, and the Coriolis
//! source balances the pressure gradient.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::frame::{self, MrfSection};
use foehncraft_fv::Discretisation;
use foehncraft_fv::problem::{Boundary, Problem};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::annulus::AnnulusSpec;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_mrf_{name}_{}", std::process::id()))
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

/// The case file with `\n` line ends, one level (refine 1, 1 024 cells).
fn text() -> String {
    std::fs::read_to_string(cases().join("taylor_couette_mrf.toml"))
        .unwrap()
        .replace("\r\n", "\n")
        .replace("levels = [0, 1, 2]\nexpected_order = [1.8, 2.6]\n", "")
}

fn parse(t: &str) -> FlowCase {
    FlowCase::from_toml(t, &cases()).unwrap()
}

fn l2(r: &FlowReport) -> f64 {
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = &r.solves[0];
    assert!(
        s.converged,
        "not converged after {} iterations",
        s.iterations
    );
    s.errors.as_ref().expect("an analytic solution").l2
}

fn discretisation(mesh: &foehncraft_mesh::Mesh, geo: &Geometry) -> Discretisation {
    let mut prob = Problem::default();
    for p in &mesh.boundaries {
        prob.boundary
            .insert(p.name.clone(), Boundary::zero_gradient());
    }
    Discretisation::new(mesh, geo, &prob).unwrap()
}

/// The frame flux closes in every cell when the interface and the
/// boundaries of the zone are surfaces of revolution about its axis, and
/// does not when the axis is moved off the centre of the mesh.
#[test]
fn the_relative_flux_closes_next_to_a_surface_of_revolution() {
    let spec = AnnulusSpec {
        n_theta: 24,
        n_r: 6,
        rotor_layers: 3,
        radial_ratio: 2.0,
        ..Default::default()
    };
    let mesh = spec.generate().unwrap();
    let geo = Geometry::compute(&mesh);
    let d = discretisation(&mesh, &geo);
    let on_axis = MrfSection {
        zone: Some("rotor".into()),
        omega: Some(10.0),
        ..Default::default()
    };
    let cells = frame::assign_cells(std::slice::from_ref(&on_axis), &mesh).unwrap();
    assert_eq!(cells.iter().filter(|c| c.is_some()).count(), 24 * 3);
    let f = frame::frame_fluxes(std::slice::from_ref(&on_axis), &cells, &d, 1.2).unwrap();
    assert!(f.worst_imbalance < 1e-12, "{}", f.worst_imbalance);
    assert_eq!(f.unbalanced_cells, 0);
    // the faces between rotor cells turn the fluid through them: some
    // frame flux is non-zero
    assert!(f.flux.iter().any(|m| m.abs() > 1e-6));
    // an axis through another point: the interface is a circle about the
    // centre, not about that point
    let off = MrfSection {
        origin: [0.02, 0.0, 0.0],
        ..on_axis
    };
    let f = frame::frame_fluxes(&[off], &cells, &d, 1.2).unwrap();
    assert!(f.worst_imbalance > 1e-3 && f.unbalanced_cells > 0, "{f:?}");
}

#[test]
fn zones_are_checked() {
    let mesh = AnnulusSpec {
        rotor_layers: 2,
        ..Default::default()
    }
    .generate()
    .unwrap();
    let rotor = |omega: Option<f64>, rpm: Option<f64>, zone: Option<&str>| MrfSection {
        zone: zone.map(String::from),
        omega,
        rpm,
        ..Default::default()
    };
    // a zone the mesh does not have, with a hint
    let e = frame::assign_cells(&[rotor(Some(1.0), None, Some("rotr"))], &mesh).unwrap_err();
    let m = format!("{e:#}");
    assert!(
        m.contains("no cell zone `rotr`") && m.contains("rotor"),
        "{m}"
    );
    // two entries on the same cells
    let e = frame::assign_cells(
        &[
            rotor(Some(1.0), None, Some("rotor")),
            rotor(Some(2.0), None, None),
        ],
        &mesh,
    )
    .unwrap_err();
    assert!(format!("{e:#}").contains("both rotate"), "{e:#}");
    // zones given: the cells of the rotor only
    let c = frame::assign_cells(&[rotor(Some(1.0), None, Some("rotor"))], &mesh).unwrap();
    assert_eq!(c.iter().filter(|c| c.is_some()).count(), 64 * 2);
}

#[test]
fn a_transient_case_is_refused() {
    let t = text() + "\n[time]\nscheme = \"bdf2\"\ntime_step = 0.01\nsteps = 10\n";
    let e = FlowCase::from_toml(&t, &cases())
        .and_then(|c| c.validate())
        .unwrap_err();
    assert_eq!(foehncraft_flow::diag::of_error(&e).code, "MRF001", "{e:#}");
}

/// The error against the analytic solution does not depend on the angular
/// velocity of the frame: in the inertial frame, in the frame of the
/// inner cylinder, in the opposite one and in a faster one.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn taylor_couette_does_not_depend_on_the_frame() {
    let base = parse(&text());
    // no frame, and a frame at rest: bitwise the same run
    let mut inertial = base.clone();
    inertial.mrf.clear();
    let e_inertial = l2(&run_case(&inertial, "inertial"));
    let mut at_rest = base.clone();
    at_rest.mrf[0].omega = Some(0.0);
    let e_rest = l2(&run_case(&at_rest, "at_rest"));
    assert_eq!(e_rest, e_inertial, "a frame at rest must change nothing");

    let mut errors = vec![("inertial", e_inertial)];
    for (name, omega) in [
        ("with_the_cylinder", 10.0),
        ("opposite", -10.0),
        ("faster", 25.0),
    ] {
        let mut c = base.clone();
        c.mrf[0].omega = Some(omega);
        errors.push((name, l2(&run_case(&c, name))));
    }
    let lo = errors.iter().map(|e| e.1).fold(f64::INFINITY, f64::min);
    let hi = errors.iter().map(|e| e.1).fold(0.0_f64, f64::max);
    eprintln!("L2 errors by frame: {errors:?}");
    assert!(hi < 0.012, "{errors:?}");
    assert!(hi < 2.0 * lo, "{errors:?}");
}

/// Fluid in solid-body rotation, both cylinders and the whole mesh turning
/// at ω: in the frame that turns with it nothing flows relative to the
/// frame; only the Coriolis source and the pressure gradient remain, and
/// the convective terms vanish. The velocity error is that of the pressure
/// gradient on a polygonal ring (second order, the same in both frames).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn solid_body_rotation_in_its_own_frame() {
    let t = text()
        .replace("zone = \"rotor\"\n", "")
        .replace("rotor_layers = 4\n", "")
        .replace(
            "[boundary.outer]\ntype = \"wall\"\n",
            "[boundary.outer]\ntype = \"wall\"\nvelocity = [\"-spin * y\", \"spin * x\", 0.0]\n",
        )
        .replace(
            "inner_omega = 10.0",
            "inner_omega = 10.0\nouter_omega = 10.0",
        );
    let with_frame = parse(&t);
    let mut inertial = with_frame.clone();
    inertial.mrf.clear();
    let e_frame = l2(&run_case(&with_frame, "solid_frame"));
    let e_inertial = l2(&run_case(&inertial, "solid_inertial"));
    eprintln!("solid body: frame {e_frame:e}, inertial {e_inertial:e}");
    assert!(
        e_frame < 0.02 && e_inertial < 0.02,
        "{e_frame} {e_inertial}"
    );
}

/// A passive scalar between the cylinders (θ = 1 inside, 0 outside, pure
/// conduction: θ = ln(r₂/r)/ln(r₂/r₁)) is convected with the relative flux
/// in the rotor and with the absolute one outside; its error against the
/// exact solution does not depend on the frame.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn a_passive_scalar_does_not_depend_on_the_frame() {
    let t = text()
        + "\n[scalars.theta]\ndiffusivity = 1e-3\ninitial = 0.0\n\
           exact = \"log(0.1 / sqrt(x^2 + y^2)) / log(2.0)\"\n\
           [boundary.inner.scalar.theta]\nvalue = 1.0\n\
           [boundary.outer.scalar.theta]\nvalue = 0.0\n";
    let base = parse(&t);
    let theta_l2 = |c: &FlowCase, name: &str| -> f64 {
        let r = run_case(c, name);
        assert!(r.failures.is_empty(), "{:?}", r.failures);
        r.solves[0].scalars[0].l2.expect("exact solution of theta")
    };
    let mut inertial = base.clone();
    inertial.mrf.clear();
    let e0 = theta_l2(&inertial, "theta_inertial");
    let e1 = theta_l2(&base, "theta_frame");
    let mut back = base.clone();
    back.mrf[0].omega = Some(-10.0);
    let e2 = theta_l2(&back, "theta_back");
    eprintln!("theta L2: inertial {e0:e}, frame {e1:e}, opposite {e2:e}");
    for e in [e0, e1, e2] {
        assert!(e < 5e-3, "{e0} {e1} {e2}");
    }
}

/// A square box that turns as a rigid body with its frame carries its fluid
/// with it: u = Ω × r everywhere, at rest in the frame. Its walls move
/// normal to themselves (a blade or a baffle of a rotor does the same), so
/// their whole velocity counts, not only the tangential part (before the
/// fix of PR #98 the walls slid along the fluid without pushing it: L2
/// 1.8e-2 on this mesh, 1.3e-3 after; what is left is the wall pressure,
/// extrapolated with zero normal gradient where the true one is ρΩ²r·n).
#[test]
fn a_turning_box_carries_its_fluid() {
    let t = r#"
case_format = "0.6"

[mesh.generate]
type = "hex"
dimension = 2
lower = [-0.05, -0.05, 0.0]
upper = [0.05, 0.05, 0.005]
cells = [16, 16, 1]

[fluid]
density = 1.0
viscosity = 1e-3

[constants]
spin = 10.0

[flow]
convection = true

[[mrf]]
axis = [0.0, 0.0, 1.0]
omega = 10.0

[boundary.xmin]
type = "wall"
velocity = ["-spin * y", "spin * x", 0.0]
[boundary.xmax]
type = "wall"
velocity = ["-spin * y", "spin * x", 0.0]
[boundary.ymin]
type = "wall"
velocity = ["-spin * y", "spin * x", 0.0]
[boundary.ymax]
type = "wall"
velocity = ["-spin * y", "spin * x", 0.0]
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"

[solver]
tolerance = 1e-10
max_iterations = 2000

[analytic]
type = "taylor-couette"
inner_radius = 0.01
outer_radius = 1.0
inner_omega = 10.0
outer_omega = 10.0
"#;
    let case = parse(t);
    let rep = run_case(&case, "turning_box");
    // no BC010 warning: the normal velocity of the walls is that of the zone
    assert!(
        !rep.findings.iter().any(|f| f.code == "BC010"),
        "{:?}",
        rep.findings
    );
    let e = l2(&rep);
    eprintln!("turning box: L2 {e:e}");
    assert!(e < 3e-3, "{e}");
}

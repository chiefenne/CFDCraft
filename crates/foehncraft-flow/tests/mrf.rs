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

/// A frame and a sliding mesh do not turn the same cells (MRF001); a
/// transient case with a frame is accepted (D185).
#[test]
fn a_transient_frame_is_accepted() {
    let t = text() + "\n[time]\nscheme = \"bdf2\"\ntime_step = 0.01\nsteps = 10\n";
    FlowCase::from_toml(&t, &cases())
        .and_then(|c| c.validate())
        .unwrap();
}

/// The text of the transient frozen rotor: the spin-up from rest to `steps`
/// steps of 0.02 s.
fn transient(steps: usize) -> String {
    text().replace(
        "\n[[mrf]]\nzone",
        &format!(
            "\n[time]\nscheme = \"bdf2\"\ntime_step = 0.02\nsteps = {steps}\nstart_time = 0.0\n\n[[mrf]]\nzone"
        ),
    )
    .replace("[solver.momentum]\nrelaxation = 0.8\n", "")
}

/// The spin-up of Taylor–Couette flow from rest (D185): the gap is a body
/// of revolution, so the spin-up in a rotating frame is the spin-up of the
/// inertial frame, at every moment. At t = 0.5 s (a tenth of the viscous
/// time of the gap) the flow is far from the steady state, and its
/// distance from it agrees in the inertial frame, the frame of the cylinder
/// and the opposite frame to the discretisation error; at t = 20 s the
/// transient frame has reached the steady solution of the steady frame.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn the_spin_up_does_not_depend_on_the_frame() {
    let base = parse(&transient(25));
    let mut inertial = base.clone();
    inertial.mrf.clear();
    let e_inertial = l2(&run_case(&inertial, "spin_inertial"));
    let mut e = vec![("inertial", e_inertial)];
    for (name, omega) in [("with_the_cylinder", 10.0), ("opposite", -10.0)] {
        let mut c = base.clone();
        c.mrf[0].omega = Some(omega);
        e.push((name, l2(&run_case(&c, &format!("spin_{name}")))));
    }
    eprintln!("spin-up at t = 0.5 s, L2 distance from the steady state: {e:?}");
    assert!(e_inertial > 0.1, "{e:?}");
    for (name, v) in &e {
        assert!((v / e_inertial - 1.0).abs() < 0.01, "{name}: {e:?}");
    }
    // the end: the steady solution
    let steady = l2(&run_case(&parse(&text()), "spin_steady"));
    let end = l2(&run_case(&parse(&transient(1000)), "spin_end"));
    eprintln!("steady frame {steady:e}, transient frame at t = 20 s {end:e}");
    assert!((end / steady - 1.0).abs() < 0.05, "{steady} {end}");
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

/// The fields of rotating machinery (D183) on Taylor–Couette flow with the
/// inner half of the gap in the frame (1 024 cells): the frame velocity is
/// Ω × x in the rotor and zero in the stator, so c_u − w_u = Ω r and 0;
/// c_u follows the exact u_θ = A r + B/r to the discretisation error, c_r
/// and the flow angle α vanish with it; ρ I = p + ½ρ|w|² − ½ρ|u|² ties the
/// relative total pressure to the rothalpy. With `output.frame = "rotor"`
/// the stator is seen in the frame of the rotor, too.
#[cfg(feature = "hdf5")]
#[test]
fn relative_and_frame_velocities_of_taylor_couette() {
    use foehncraft_flow::derived;
    let t = text().replace(
        "[output]\nvtkhdf = false\n",
        "[output]\nvtkhdf = false\nderived = [\"relative-velocity\", \"frame-velocity\", \
         \"turbomachinery\", \"relative-total-pressure\", \"rothalpy\"]\n",
    ) + "\n[run]\nrestart_final = true\n";
    assert!(t.contains("rothalpy"), "the case file changed");
    let (omega, rho, r_iface) = (10.0, 1.0, 0.075);
    let exact = |r: f64| -10.0 * 0.0025 / 0.0075 * r + 10.0 * 0.0025 * 0.01 / 0.0075 / r;
    for fixed in [false, true] {
        let mut case = parse(&t);
        if fixed {
            case.output.frame = Some("rotor".into());
        }
        let name = if fixed { "fields_fixed" } else { "fields" };
        let dir = out(name);
        let rep = run(&case, name, &dir, |_| {}).unwrap();
        assert!(rep.failures.is_empty(), "{:?}", rep.failures);
        let (mesh, geo, f) =
            derived::evaluate(&case, &dir.join("restart.h5"), &case.output.derived).unwrap();
        let _ = std::fs::remove_dir_all(&dir);
        let g = |n: &str| f.get(n).unwrap_or_else(|| panic!("{n}"));
        let (cu, wu, cr) = (
            g("velocity_circumferential"),
            g("relative_velocity_circumferential"),
            g("velocity_radial"),
        );
        let (fx, fy, fm) = (
            g("frame_velocity_x"),
            g("frame_velocity_y"),
            g("frame_velocity_magnitude"),
        );
        let (wm, ptr, rot, alpha) = (
            g("relative_velocity_magnitude"),
            g("relative_total_pressure"),
            g("rothalpy"),
            g("flow_angle"),
        );
        let mut worst = [0.0f64; 3];
        for c in 0..mesh.n_cells {
            let x = geo.cell_centroid[c];
            let r = x[0].hypot(x[1]);
            let turning = fixed || r < r_iface;
            let ub = if turning { omega * r } else { 0.0 };
            assert!((cu[c] - wu[c] - ub).abs() < 1e-12, "c_u − w_u at r = {r}");
            let (ex, ey) = if turning {
                (-omega * x[1], omega * x[0])
            } else {
                (0.0, 0.0)
            };
            assert!((fx[c] - ex).abs() < 1e-12 && (fy[c] - ey).abs() < 1e-12);
            assert!((fm[c] - ub).abs() < 1e-12);
            // |w| ≥ |w_u|, the rest is the small radial velocity
            assert!(wm[c] + 1e-12 >= wu[c].abs());
            // ρ I = p + ½ρ|w|² − ½ρ|u|²
            assert!((ptr[c] - rho * rot[c] - 0.5 * rho * ub * ub).abs() < 1e-9);
            worst[0] = worst[0].max((cu[c] - exact(r)).abs());
            worst[1] = worst[1].max(cr[c].abs());
            worst[2] = worst[2].max(alpha[c].abs());
        }
        eprintln!(
            "frame {fixed}: |c_u − u_θ| {:.2e}, |c_r| {:.2e}, |α| {:.2e}°",
            worst[0], worst[1], worst[2]
        );
        // to the discretisation error (L∞ 2.95 % of the wall speed 0.5 m/s)
        assert!(worst[0] < 0.035 * 0.5, "{worst:?}");
        assert!(worst[1] < 0.035 * 0.5, "{worst:?}");
        assert!(worst[2] < 5.0, "{worst:?}");
    }
    // a derived field of rotors without a turning zone, and an unknown frame
    let mut c = parse(&t);
    c.mrf.clear();
    let e = format!("{:#}", c.validate().unwrap_err());
    assert!(e.contains("needs a turning zone"), "{e}");
    let mut c = parse(&t);
    c.output.frame = Some("rotr".into());
    let e = format!("{:#}", c.validate().unwrap_err());
    assert!(e.contains("did you mean \"rotor\""), "{e}");
}

/// Runs a case and reads the rows of its line `line` (`line_<line>.csv`).
fn line_rows(
    case: &FlowCase,
    name: &str,
    line: &str,
) -> Vec<std::collections::HashMap<String, f64>> {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    assert!(rep.failures.is_empty(), "{name}: {:?}", rep.failures);
    assert!(rep.solves[0].converged, "{name}: not converged");
    let text = std::fs::read_to_string(dir.join(format!("line_{line}.csv"))).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let mut rows = text.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<String> = rows.next().unwrap().split(',').map(String::from).collect();
    rows.map(|l| {
        head.iter()
            .cloned()
            .zip(l.split(',').map(|v| v.parse::<f64>().unwrap()))
            .collect()
    })
    .collect()
}

/// Turbulent Taylor–Couette flow (k-ω SST, Re = Ω r₁ (r₂ − r₁)/ν = 12 500,
/// steady) with the inner half of the gap in the frame of the cylinder and
/// in the inertial frame: the turbulence is convected with the relative
/// flux in the rotor and produced by the strain, which does not depend on
/// the frame, so the profiles of u_θ, k and ν_t across the gap agree to the
/// discretisation error (0.2 %, 0.4 % and 0.3 % of the largest value here,
/// less with more sectors). The first cells are at y⁺ 0.8: with the first
/// cell in the buffer layer (y⁺ 5.6, TURB001) its k depends on the frame by
/// 13 %, and in the inertial frame alone by a factor of 3 between 32 and
/// 256 sectors, a cell that the automatic wall treatment does not resolve.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn turbulent_taylor_couette_does_not_depend_on_the_frame() {
    let t = text()
        .replace(
            "n_r = 8\nrotor_layers = 4\n",
            "n_r = 16\nrotor_layers = 8\nradial_ratio = 30.0\nsymmetric = true\n",
        )
        .replace(
            "density = 1.0\nviscosity = 5e-4",
            "density = 1000.0\nviscosity = 2e-3",
        )
        .replace(
            "[flow]\n",
            "[turbulence]\nmodel = \"k-omega-sst\"\nwall_treatment = \"automatic\"\n\n\
             [initial]\nturbulence_intensity = 0.05\nviscosity_ratio = 10.0\n\
             velocity = [0.0, 0.0, 0.0]\n\n[flow]\n",
        )
        .replace(
            "[solver]\ntolerance = 1e-10\nmax_iterations = 3000\n",
            "[solver]\ntolerance = 1e-7\nmax_iterations = 4000\n",
        )
        .replace(
            "[output]\n",
            "[output]\nlines = [{ name = \"gap\", from = [0.0505, 0.0, 0.005], \
             to = [0.0995, 0.0, 0.005], points = 40 }]\n",
        );
    let mut frame = parse(&t);
    frame.analytic = None;
    frame.verification = Default::default();
    assert!(!frame.turbulence.is_laminar() && frame.output.lines.len() == 1);
    let mut inertial = frame.clone();
    inertial.mrf.clear();
    let a = line_rows(&inertial, "turb_inertial", "gap");
    let b = line_rows(&frame, "turb_frame", "gap");
    let d: Vec<f64> = ["v", "k", "nu_t"]
        .iter()
        .map(|q| {
            let scale = a.iter().map(|r| r[*q].abs()).fold(0.0, f64::max);
            a.iter()
                .zip(&b)
                .map(|(x, y)| (x[*q] - y[*q]).abs())
                .fold(0.0, f64::max)
                / scale
        })
        .collect();
    eprintln!("turbulent Taylor–Couette, frame against inertial: u_θ, k, ν_t {d:?}");
    assert!(d[0] < 0.005 && d[1] < 0.015 && d[2] < 0.015, "{d:?}");
}

/// Two independent rotating zones: the inner half of the gap turns with
/// the inner cylinder (10 rad/s), the outer half with the outer one, which
/// turns the other way (−4 rad/s). Every frame sees its own relative flow,
/// and the solution is the inertial one, A r + B/r for the two cylinder
/// speeds, to the discretisation error.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn two_independent_frames() {
    let t = text()
        .replace(
            "[boundary.outer]\ntype = \"wall\"\n",
            "[boundary.outer]\ntype = \"wall\"\nvelocity = [\"4.0 * y\", \"-4.0 * x\", 0.0]\n",
        )
        .replace(
            "inner_omega = 10.0",
            "inner_omega = 10.0\nouter_omega = -4.0",
        )
        + "\n[[mrf]]\nzone = \"stator\"\naxis = [0.0, 0.0, 1.0]\nomega = -4.0\n";
    let two = parse(&t);
    assert_eq!(two.mrf.len(), 2);
    let mut inertial = two.clone();
    inertial.mrf.clear();
    let e_two = l2(&run_case(&two, "two_frames"));
    let e_inertial = l2(&run_case(&inertial, "two_inertial"));
    eprintln!("two frames: L2 {e_two:e}, inertial {e_inertial:e}");
    assert!(e_two < 0.012, "{e_two}");
    assert!(
        e_two < 2.0 * e_inertial && e_inertial < 2.0 * e_two,
        "{e_two} {e_inertial}"
    );
}

/// Torque and shaft power of the inner cylinder of Taylor–Couette flow
/// (D184): T = −4π μ ω r₁² r₂² L / (r₂² − r₁²), in any frame; the shaft
/// power −Ω·T equals the work of the wall on the fluid and the power the
/// kinetic-energy ledger books for it; the stator carries the opposite
/// torque.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn torque_and_shaft_power_of_taylor_couette_flow() {
    let (r1, r2, depth, mu, w) = (0.05_f64, 0.10_f64, 0.01_f64, 5e-4_f64, 10.0_f64);
    let t_exact =
        -4.0 * std::f64::consts::PI * mu * w * r1 * r1 * r2 * r2 * depth / (r2 * r2 - r1 * r1);
    let base = parse(&text());
    for (name, omega) in [
        ("t_inner", 10.0),
        ("t_inertial_zone", 0.0),
        ("t_faster", 25.0),
    ] {
        let mut c = base.clone();
        c.mrf[0].omega = Some(omega);
        let rep = run_case(&c, name);
        let s = &rep.solves[0];
        let v = s.rotors.as_ref().expect("a rotor");
        let r = v.rotor("rotor").expect("the rotor");
        assert_eq!(r.walls, ["inner"]);
        let e = (r.torque / t_exact - 1.0).abs();
        eprintln!(
            "{name}: torque {:e} (exact {t_exact:e}, {:.2} %), power {:e}, wall work {:e}, slip {:e}, angular momentum {:?}",
            r.torque,
            100.0 * e,
            r.power,
            r.wall_work,
            r.slip,
            r.angular_momentum
        );
        assert!(e < 0.03, "{name}: torque {} against {t_exact}", r.torque);
        // the torque does not depend on the frame; the power is Ω·T of the zone
        assert!(
            (r.power + omega.abs() * r.torque).abs() <= 1e-12 * r.power.abs() + 1e-300
                || omega == 0.0,
            "{name}"
        );
        // the outer wall carries the opposite torque
        let outer = s.boundaries.iter().find(|b| b.name == "outer").unwrap();
        let inner = s.boundaries.iter().find(|b| b.name == "inner").unwrap();
        assert!(
            (inner.torque.unwrap() + outer.torque.unwrap()).abs() < 0.02 * t_exact.abs(),
            "{name}: {:?} {:?}",
            inner.torque,
            outer.torque
        );
        // when the zone turns with the wall, the shaft power is its work on
        // the fluid and the ledger books it
        if omega == 10.0 {
            assert!(r.slip < 1e-12, "{}", r.slip);
            assert!(
                ((r.wall_work - r.power) / r.power).abs() < 1e-10,
                "{} {}",
                r.wall_work,
                r.power
            );
            let b = s.balance.as_ref().unwrap();
            let k = &b.last.kinetic_energy;
            let inner_power: f64 = b
                .last
                .boundaries
                .iter()
                .find(|p| p.boundary == "inner")
                .unwrap()
                .power();
            eprintln!(
                "ledger: shaft {:e}, inner wall power {:e}, dissipation {:e}",
                k.shaft, inner_power, k.viscous_dissipation
            );
            let tr = r.transposed_work;
            assert!(
                ((k.shaft - inner_power - tr) / k.shaft).abs() < 5e-3,
                "{} {} {tr}",
                k.shaft,
                inner_power
            );
            // the shaft power is the exact dissipation of the flow
            assert!(((k.shaft - w * t_exact.abs()) / k.shaft).abs() < 0.03);
        }
    }
}

/// The performance values of a machine (D184) on a flow with a known
/// pressure loss: plane Poiseuille flow through a zone that does not turn
/// (Ω = 0): Δp_t = −12 μ U L / h², Q = U h d, no shaft power, no Euler work.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn machine_values_of_poiseuille_flow() {
    let t = std::fs::read_to_string(cases().join("poiseuille_hex.toml"))
        .unwrap()
        .replace("\r\n", "\n")
        .replace("levels = [0, 1, 2]\nexpected_order = [1.8, 2.2]\n", "")
        + "\n[[mrf]]\naxis = [0.0, 0.0, 1.0]\nomega = 0.0\n\n[[output.machines]]\nname = \"duct\"\ninlet = [\"xmin\"]\noutlet = [\"xmax\"]\n";
    let rep = run_case(&parse(&t), "machine_poiseuille");
    let v = rep.solves[0].rotors.as_ref().expect("values");
    let m = v.machine("duct").expect("the machine");
    eprintln!("{m:?}");
    let dpt = -12.0 * 0.1 * 1.0 * 2.0 / 1.0;
    assert!(
        (m.total_pressure_rise / dpt - 1.0).abs() < 0.02,
        "{} {dpt}",
        m.total_pressure_rise
    );
    assert!(
        (m.volume_flow - 0.1).abs() < 5e-3 * 0.1,
        "{}",
        m.volume_flow
    );
    assert!((m.head - dpt / 9.80665).abs() < 0.02 * (dpt / 9.80665).abs());
    assert!(
        (m.hydraulic_power / (0.1 * dpt) - 1.0).abs() < 0.02,
        "{}",
        m.hydraulic_power
    );
    assert_eq!(m.shaft_power, 0.0);
    assert!(m.efficiency.is_nan());
    assert!(m.euler_work.abs() < 1e-9, "{}", m.euler_work);
}

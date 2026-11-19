//! The sliding mesh with conformal stepping (D56, D182, stage 6 part 6 of
//! 6g): the rotor of a ring mesh turns by whole sectors per time step.
//!
//! The checks: the shift of the rotor cells (a bijection, N of them are
//! the identity), the refusals of the case file, laminar Taylor–Couette
//! flow from rest against the exact solution on a series of meshes
//! (second order, Δt ∝ h), a fluid in solid-body rotation that the scheme
//! keeps, a rotor with a shift of two sectors, and a passive scalar.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use foehncraft_flow::frame::{self, MrfSection};
use foehncraft_flow::sliding::{SlidingSection, plan};
use foehncraft_fv::Discretisation;
use foehncraft_fv::problem::{Boundary, Problem};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::annulus::AnnulusSpec;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_sliding_{name}_{}", std::process::id()))
}

fn try_run(case: &FlowCase, name: &str) -> anyhow::Result<FlowReport> {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {});
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    try_run(case, name).unwrap()
}

fn text() -> String {
    std::fs::read_to_string(cases().join("taylor_couette_sliding.toml"))
        .unwrap()
        .replace("\r\n", "\n")
}

fn parse(t: &str) -> FlowCase {
    FlowCase::from_toml(t, &cases()).unwrap()
}

fn code(e: &anyhow::Error) -> String {
    diag::of_error(e).code.to_string()
}

fn l2(r: &FlowReport) -> f64 {
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = &r.solves[0];
    assert!(s.converged);
    s.errors.as_ref().expect("an analytic solution").l2
}

/// The case at refinement `level` (16 · 2^level sectors) with the time step
/// of one sector per step, to the time `end`.
fn at_level(base: &FlowCase, level: u32, shift: usize, end: f64) -> FlowCase {
    let mut c = base.clone();
    c.verification.max_error = None;
    c.mesh.annulus.as_mut().unwrap().refine = level;
    c.sliding[0].shift = shift;
    let sectors = 16usize << level;
    let dt = c.sliding[0].time_step(sectors).unwrap();
    let mut t = c.time.unwrap();
    t.time_step = dt;
    t.steps = Some((end / dt).round() as usize);
    t.end_time = None;
    c.time = Some(t);
    c
}

fn discretisation(mesh: &foehncraft_mesh::Mesh, geo: &Geometry) -> Discretisation {
    let mut prob = Problem::default();
    for p in &mesh.boundaries {
        prob.boundary
            .insert(p.name.clone(), Boundary::zero_gradient());
    }
    Discretisation::new(mesh, geo, &prob).unwrap()
}

/// The rotor of a ring of 24 sectors: found from the geometry, shifted by
/// whole sectors (a bijection that returns after N steps).
#[test]
fn the_rotor_cells_are_shifted_by_whole_sectors() {
    let mesh = AnnulusSpec {
        n_theta: 24,
        n_r: 6,
        rotor_layers: 3,
        radial_ratio: 2.0,
        ..Default::default()
    }
    .generate()
    .unwrap();
    let geo = Geometry::compute(&mesh);
    let d = discretisation(&mesh, &geo);
    let entry = SlidingSection {
        zone: Some("rotor".into()),
        omega: Some(10.0),
        ..Default::default()
    };
    let cell_zone = frame::assign_cells(&[entry.frame()], &mesh).unwrap();
    let rotor: Vec<bool> = cell_zone.iter().map(Option::is_some).collect();
    let p = plan(&entry, &rotor, &d).unwrap();
    assert_eq!((p.sectors, p.rotor_cells), (24, 24 * 3));
    assert!((p.angle - std::f64::consts::TAU / 24.0).abs() < 1e-14);
    // a permutation of the cells that fixes the stator and every ring
    let mut seen = vec![false; p.cell.len()];
    for (i, &c) in p.cell.iter().enumerate() {
        assert!(!seen[c], "cell {c} twice");
        seen[c] = true;
        assert_eq!(rotor[i], rotor[c], "the shift keeps the rotor a rotor");
        if !rotor[i] {
            assert_eq!(i, c);
        } else {
            let (a, b) = (d.cell_centroid.get(i), d.cell_centroid.get(c));
            let (ra, rb) = (a[0].hypot(a[1]), b[0].hypot(b[1]));
            assert!((ra - rb).abs() < 1e-12 && (a[2] - b[2]).abs() < 1e-12);
        }
    }
    // N steps are a whole turn: the identity
    let mut at: Vec<usize> = (0..p.cell.len()).collect();
    for _ in 0..24 {
        at = at.iter().map(|&i| p.cell[i]).collect();
    }
    assert!(at.iter().enumerate().all(|(i, &c)| i == c));
    // two sectors per step: the shift twice
    let two = SlidingSection {
        shift: 2,
        ..entry.clone()
    };
    let p2 = plan(&two, &rotor, &d).unwrap();
    assert!(
        p2.cell
            .iter()
            .enumerate()
            .all(|(i, &c)| c == p.cell[p.cell[i]])
    );
    assert!((p2.angle - 2.0 * p.angle).abs() < 1e-14);
    // a whole turn per step is no step
    let whole = SlidingSection {
        shift: 24,
        ..entry.clone()
    };
    assert_eq!(code(&plan(&whole, &rotor, &d).unwrap_err()), "SLD002");
    // an axis that does not go through the centre: no cell has an image
    let off = SlidingSection {
        origin: [0.013, 0.0, 0.0],
        ..entry.clone()
    };
    let e = plan(&off, &rotor, &d).unwrap_err();
    assert_eq!(code(&e), "SLD004", "{e:#}");
    // the angular velocity sets the sense only: the shift is the same set of
    // pairs turned the other way
    let back = SlidingSection {
        omega: Some(-10.0),
        ..entry
    };
    let pb = plan(&back, &rotor, &d).unwrap();
    assert!(
        pb.cell.iter().enumerate().all(|(i, &c)| p.cell[c] == i),
        "the opposite turn is the inverse shift"
    );
}

/// The checks of the tables and of the time stepping.
#[test]
fn the_case_is_checked() {
    let valid = |t: &str| -> Result<(), String> {
        FlowCase::from_toml(t, &cases())
            .and_then(|c| c.validate())
            .map_err(|e| code(&e))
    };
    let base = text();
    valid(&base).unwrap();
    // a steady case
    let no_time = {
        let i = base.find("\n[time]\n").unwrap();
        let j = base.find("\n[[sliding]]\n").unwrap();
        format!("{}{}", &base[..i], &base[j..])
    };
    assert_eq!(valid(&no_time).unwrap_err(), "SLD001");
    // the entry
    let two = base.clone() + "\n[[sliding]]\nomega = 1.0\n";
    assert_eq!(valid(&two).unwrap_err(), "SLD002");
    let both = base.replace(
        "omega = 10.0\n\n[boundary",
        "omega = 10.0\nrpm = 5.0\n\n[boundary",
    );
    assert_eq!(valid(&both).unwrap_err(), "SLD002");
    let none = base.replace("omega = 10.0\n\n[boundary", "\n[boundary");
    assert_eq!(valid(&none).unwrap_err(), "SLD002");
    let stopped = base.replace("omega = 10.0\n\n[boundary", "omega = 0.0\n\n[boundary");
    assert_eq!(valid(&stopped).unwrap_err(), "SLD002");
    let zero = base.replace(
        "omega = 10.0\n\n[boundary",
        "omega = 10.0\nshift = 0\n\n[boundary",
    );
    assert_eq!(valid(&zero).unwrap_err(), "SLD002");
    // the stepping
    let adaptive = base.replace("[time]\n", "[time]\nadaptive = true\n");
    assert_eq!(valid(&adaptive).unwrap_err(), "SLD003");
    let sweep = base.replace("[verification]\n", "[verification]\nlevels = [0, 1]\n");
    assert_eq!(valid(&sweep).unwrap_err(), "SLD003");
}

/// What the mesh decides: the time step, the end time, the zone, the
/// boundaries and the physics that the sliding mesh does not do yet.
#[test]
fn the_mesh_and_the_step_are_checked() {
    let base = text();
    let refused = |t: String| -> (String, String) {
        let c = parse(&t);
        let e = try_run(&c, "refused").unwrap_err();
        (code(&e), format!("{e:#}"))
    };
    // the time step that the sectors and the speed give is named
    let (c, m) = refused(base.replace("time_step = 0.019634954084936207", "time_step = 0.02"));
    assert_eq!(c, "SLD003", "{m}");
    assert!(m.contains("1.9634954085e-2"), "{m}");
    // an end time that is not a whole number of steps
    let (c, _) = refused(base.replace("steps = 510", "end_time = 10.0"));
    assert_eq!(c, "SLD003");
    // the zone
    let (c, m) = refused(base.replace("zone = \"rotor\"", "zone = \"rotr\""));
    assert_eq!(c, "SLD002", "{m}");
    assert!(m.contains("did you mean `rotor`"), "{m}");
    // a turbulence model
    let (c, _) = refused(base.replace(
        "[flow]\n",
        "[turbulence]\nmodel = \"k-omega-sst\"\n\n[flow]\n",
    ));
    assert_eq!(c, "SLD004");
    // a pressure boundary on the rotor
    let t = base.replace(
        "[boundary.inner]\ntype = \"wall\"\nvelocity = [\"-spin * y\", \"spin * x\", 0.0]\n",
        "[boundary.inner]\ntype = \"static-pressure-outlet\"\npressure = 0.0\n",
    );
    let (c, m) = refused(t);
    assert_eq!(c, "SLD004", "{m}");
    assert!(m.contains("boundary `inner`"), "{m}");
}

/// Taylor–Couette flow from rest: the error against the exact solution at
/// the end of the spin-up falls with the second power of the sector angle
/// and the time step (128, 512 and 2 048 cells).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn taylor_couette_converges_with_second_order() {
    let base = parse(&text());
    let e: Vec<f64> = (0..3)
        .map(|l| l2(&run_case(&at_level(&base, l, 1, 10.0), &format!("tc{l}"))))
        .collect();
    let order: Vec<f64> = e.windows(2).map(|w| (w[0] / w[1]).log2()).collect();
    eprintln!("sliding Taylor–Couette: L2 {e:?}, orders {order:?}");
    assert!(e[1] < 0.02 && e[2] < 0.004, "{e:?}");
    assert!(
        order.iter().all(|o| (1.9..2.7).contains(o)),
        "orders {order:?}"
    );
}

/// Two sectors per step: the same flow with a time step twice as long and
/// four times the error of the angle (second order in the step).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn two_sectors_per_step() {
    let base = parse(&text());
    let one = l2(&run_case(&at_level(&base, 1, 1, 10.0), "shift1"));
    let two = l2(&run_case(&at_level(&base, 1, 2, 10.0), "shift2"));
    eprintln!("L2: one sector per step {one:e}, two {two:e}");
    assert!(two < 0.08 && two > one, "{one} {two}");
}

/// A fluid in solid-body rotation with the whole mesh turning, started from
/// the exact solution: the scheme keeps it (the pressure gradient turns the
/// vectors of the cells that move), with an error that does not grow.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn solid_body_rotation_is_kept() {
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
        )
        .replace("[flow]\n", "[initial]\nfrom_analytic = true\n\n[flow]\n")
        .replace("max_error = 0.02", "max_error = 0.1");
    let base = parse(&t);
    // 0.6 s and 3 s: the error stays where it settled (the first two steps
    // of a run are backward differences of first order)
    let short = l2(&run_case(&at_level(&base, 1, 1, 0.6), "solid_short"));
    let long = l2(&run_case(&at_level(&base, 1, 1, 3.0), "solid_long"));
    eprintln!("solid body, L2: {short:e} at 0.6 s, {long:e} at 3 s");
    assert!(short < 0.03 && long < 0.03, "{short} {long}");
    assert!(long < 1.5 * short, "{short} {long}");
}

/// A passive scalar between the cylinders (θ = 1 inside, 0 outside, pure
/// conduction) with the rotor turning through it: the old level of a moved
/// cell is what that cell held, and the scalar ends in the exact solution.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn a_passive_scalar_ends_in_the_exact_solution() {
    let t = text()
        + "\n[scalars.theta]\ndiffusivity = 1e-3\ninitial = 0.0\n\
           exact = \"log(0.1 / sqrt(x^2 + y^2)) / log(2.0)\"\n\
           [boundary.inner.scalar.theta]\nvalue = 1.0\n\
           [boundary.outer.scalar.theta]\nvalue = 0.0\n";
    let base = parse(&t);
    let r = run_case(&at_level(&base, 1, 1, 12.0), "theta");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let e = r.solves[0].scalars[0].l2.expect("exact solution of theta");
    eprintln!("theta L2 {e:e}");
    assert!(e < 1.5e-2, "{e}");
}

/// The rotating zone as a moving reference frame of the same case is
/// steady only: the table of a sliding mesh needs a time, the table of the
/// frame must not have one (D55, D182).
#[test]
fn a_frame_and_a_sliding_mesh_do_not_mix() {
    let t = text() + "\n[[mrf]]\nzone = \"rotor\"\nomega = 10.0\n";
    let e = FlowCase::from_toml(&t, &cases())
        .and_then(|c| c.validate())
        .unwrap_err();
    assert_eq!(code(&e), "MRF001", "{e:#}");
    let _ = MrfSection::default();
}

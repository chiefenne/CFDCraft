//! The transient `foehn flow` case (a case file with a `[time]` table):
//! the time loop with PISO ([`crate::transient`]), the adaptive time step,
//! restart files (D42), control files and re-reading the case file (D40),
//! probes, force histories and field series (D43), and the transient
//! verification: order sweeps in space (levels) and time (time steps)
//! against the Taylor–Green vortex, and the cylinder benchmark of Schäfer &
//! Turek (1996).
//!
//! Output in addition to the steady one: `history.csv` (per step: time, Δt,
//! Courant number, residuals, outer and CG iterations), `forces.csv` (force
//! per wall and the coefficients), `probes.csv`, `fields_<step>.vtkhdf`,
//! `restart.h5`.

use crate::analytic::Benchmark;
use crate::bc::{self, BoundaryCondition, InletProfile, Level, inlet_velocity, patch_extent};
use crate::case::{
    ErrorNorms, FlowCase, FlowReport, Log, Solve, dimension, error_norms, locate, mesh_for,
    patch_summaries, write_fields, write_outputs,
};
use crate::control::Control;
use crate::engine::Steer;
use crate::monitor::{self, History, Probes};
use crate::plot::{self, Series};
use crate::restart::{self, Restart};
use crate::simplec::{self, Flow, Physics, Solution, Target};
use crate::transient::{StepReport, TimeScheme, Transient};
use anyhow::{Context, Result, bail};
use foehncraft_fv::disc::Soa3;
use foehncraft_mesh::geometry::norm;
use foehncraft_mesh::locate::Locator;
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;
use std::path::Path;
use std::time::Instant;

/// Summary of a transient run.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TransientSummary {
    /// Time integration.
    pub scheme: String,
    /// Time steps done in this run.
    pub steps: usize,
    /// Step number at the end (with a restart: counted from the start of
    /// the first run).
    pub last_step: usize,
    /// Time at the end (s).
    pub time: f64,
    /// Smallest and largest time step (s).
    pub dt_min: f64,
    /// See `dt_min`.
    pub dt_max: f64,
    /// Largest Courant number.
    pub cfl_max: f64,
    /// Mean outer iterations per step.
    pub outer_mean: f64,
    /// Mean CG iterations per step.
    pub cg_mean: f64,
    /// Restart file the run continued from.
    pub restart_from: Option<String>,
    /// Ended by `.stop`.
    pub stopped: bool,
    /// Nanoseconds per cell and time step.
    pub ns_per_cell_step: f64,
}

/// Time-dependent boundary faces: boundary index and the face.
struct AnalyticFaces {
    faces: Vec<(usize, usize)>,
    ub: [Vec<f64>; 3],
}

impl AnalyticFaces {
    fn new(case: &FlowCase, mesh: &foehncraft_mesh::Mesh, ub: &[Vec<f64>; 3]) -> Self {
        let ni = mesh.n_interior_faces();
        let mut faces = Vec::new();
        for p in &mesh.patches {
            if let Some(BoundaryCondition::VelocityInlet {
                profile: Some(InletProfile::Analytic),
                ..
            }) = case.boundary.get(&p.name)
            {
                faces.extend(p.range().map(|f| (f - ni, f)));
            }
        }
        Self {
            faces,
            ub: ub.clone(),
        }
    }

    /// Sets the analytic velocity of time `t` on the faces.
    fn update(
        &mut self,
        case: &FlowCase,
        mesh: &foehncraft_mesh::Mesh,
        geo: &foehncraft_mesh::Geometry,
        fl: &mut Flow,
        t: f64,
    ) -> Result<()> {
        let Some(a) = case.analytic_at(t) else {
            return Ok(());
        };
        if self.faces.is_empty() || !a.is_transient() {
            return Ok(());
        }
        for p in &mesh.patches {
            let Some(b) = case.boundary.get(&p.name) else {
                continue;
            };
            let ext = patch_extent(mesh, p);
            for f in p.range() {
                if let Some((fb, _)) = self.faces.iter().find(|x| x.1 == f) {
                    let v = inlet_velocity(b, geo, f, Some(&a), &ext).unwrap_or([0.0; 3]);
                    for (ub, vi) in self.ub.iter_mut().zip(v) {
                        ub[*fb] = vi;
                    }
                }
            }
        }
        fl.set_boundary_velocity(&self.ub)
    }
}

/// A run of the transient sweep: a refinement level and a time step.
#[derive(Debug, Clone, Copy)]
struct RunSpec {
    level: Option<u32>,
    dt: f64,
}

/// The result of one run that the sweeps need.
struct RunOut {
    solve: Solve,
    u: Soa3,
    volume: Vec<f64>,
    stopped: bool,
}

/// Runs a transient case; `log` receives progress lines, `steer` carries
/// the requests, limits and run record (crate::engine).
pub fn run(
    case: &FlowCase,
    name: &str,
    out: &Path,
    steer: &mut Steer,
    log: impl FnMut(&str),
) -> Result<FlowReport> {
    let t0 = Instant::now();
    std::fs::create_dir_all(out)?;
    let control = Control::new(out)?;
    let mut log = log;
    let mut lg = Log {
        out: &mut log,
        text: String::new(),
    };
    let time = case.time.context("a transient case needs a [time] table")?;
    let v = &case.verification;
    let specs: Vec<RunSpec> = if !v.time_steps.is_empty() {
        let level = v.levels.first().copied();
        v.time_steps
            .iter()
            .map(|&dt| RunSpec { level, dt })
            .collect()
    } else if v.levels.is_empty() {
        vec![RunSpec {
            level: None,
            dt: time.time_step,
        }]
    } else {
        v.levels
            .iter()
            .map(|&l| RunSpec {
                level: Some(l),
                dt: time.time_step,
            })
            .collect()
    };
    if specs.len() > 1 && case.run.restart.is_some() {
        bail!("[run] restart continues a single run, not a sweep over levels or time steps");
    }
    let analytic = case.analytic_at(time.start_time);
    let reynolds = case.reynolds();
    let mut rep = FlowReport {
        name: name.into(),
        reynolds,
        findings: Vec::new(),
        solves: Vec::new(),
        orders: Vec::new(),
        time_orders: Vec::new(),
        stopped: false,
        stop: None,
        warnings: Vec::new(),
        failures: Vec::new(),
        seconds: 0.0,
    };
    std::fs::write(
        out.join("case.resolved.toml"),
        toml::to_string_pretty(case)?,
    )?;
    let mut outs: Vec<RunOut> = Vec::new();
    let last = specs.len() - 1;
    for (ri, spec) in specs.iter().enumerate() {
        let (mesh, geo, desc) = mesh_for(case, spec.level)?;
        let findings = bc::check(
            &mesh,
            &geo,
            &case.boundary,
            case.fluid.density,
            analytic.as_ref(),
            true,
        );
        if ri == 0 {
            for f in &findings {
                lg.line(&f.to_string());
            }
            rep.findings = findings.clone();
        }
        if findings.iter().any(|f| f.level == Level::Error) {
            let msg: Vec<String> = findings
                .iter()
                .filter(|f| f.level == Level::Error)
                .map(|f| f.to_string())
                .collect();
            std::fs::write(out.join("flow.log"), &lg.text)?;
            bail!(
                "the boundary conditions cannot be run:\n  {}",
                msg.join("\n  ")
            );
        }
        if ri == 0 || ri == last {
            header(&mut lg, case, &mesh, &desc, ri == last, spec.dt);
        }
        let label = if specs.len() > 1 {
            format!("{desc}, Δt = {:e} s", spec.dt)
        } else {
            desc.clone()
        };
        lg.line(&format!("mesh  {label}"));
        let r = run_one(
            case,
            name,
            out,
            &mut lg,
            &control,
            steer,
            &mut rep.warnings,
            &mesh,
            &geo,
            &label,
            spec.dt,
            ri == last,
            specs.len() == 1,
        )?;
        rep.stopped |= r.stopped;
        outs.push(r);
        if rep.stopped {
            rep.stop = steer.stop.clone();
            break;
        }
    }
    // sweeps
    let v = &case.verification;
    if !v.time_steps.is_empty() && outs.len() == specs.len() {
        // differences of successive solutions on one mesh
        let diffs: Vec<f64> = outs
            .windows(2)
            .map(|w| {
                let (a, b) = (&w[0], &w[1]);
                let (mut e2, mut vol) = (0.0, 0.0);
                for c in 0..a.volume.len() {
                    let d = [
                        a.u.x[c] - b.u.x[c],
                        a.u.y[c] - b.u.y[c],
                        a.u.z[c] - b.u.z[c],
                    ];
                    e2 += (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) * a.volume[c];
                    vol += a.volume[c];
                }
                (e2 / vol).sqrt()
            })
            .collect();
        let dts = &v.time_steps;
        rep.time_orders = diffs
            .windows(2)
            .enumerate()
            .map(|(i, d)| (d[0] / d[1]).ln() / (dts[i] / dts[i + 1]).ln())
            .collect();
        lg.line(&format!(
            "time sweep: differences of successive solutions [{}], observed time orders {:.3?}",
            diffs
                .iter()
                .map(|d| format!("{d:.3e}"))
                .collect::<Vec<_>>()
                .join(", "),
            rep.time_orders
        ));
        if let (Some([lo, hi]), Some(&o)) = (v.expected_time_order, rep.time_orders.last())
            && !(lo..=hi).contains(&o)
        {
            rep.failures.push(format!(
                "observed time order {o:.3} on the smallest steps, outside [{lo}, {hi}]"
            ));
        }
        let pts: Vec<(f64, f64)> = dts[1..]
            .iter()
            .copied()
            .zip(diffs.iter().copied())
            .collect();
        let slope: Vec<(f64, f64)> = match (pts.first(), pts.last()) {
            (Some(&a), Some(&b)) => vec![(a.0, b.1 * (a.0 / b.0).powi(2)), b],
            _ => Vec::new(),
        };
        std::fs::write(
            out.join("time_convergence.svg"),
            plot::chart(
                &format!("{name}: change of the solution over the time step"),
                "Δt (s)",
                "‖u(Δt) − u(2Δt)‖ (m/s)",
                &[
                    Series {
                        label: "Foehncraft",
                        points: &pts,
                        markers: true,
                    },
                    Series {
                        label: "slope 2",
                        points: &slope,
                        markers: false,
                    },
                ],
                true,
                true,
            ),
        )?;
    }
    if v.levels.len() >= 2 && v.time_steps.is_empty() && outs.len() == specs.len() {
        let e: Vec<&ErrorNorms> = outs
            .iter()
            .filter_map(|s| s.solve.errors.as_ref())
            .collect();
        rep.orders = e
            .windows(2)
            .map(|w| (w[0].l2 / w[1].l2).ln() / (w[0].h / w[1].h).ln())
            .collect();
        lg.line(&format!("observed orders in space {:.3?}", rep.orders));
        if let (Some([lo, hi]), Some(&o)) = (v.expected_order, rep.orders.last())
            && !(lo..=hi).contains(&o)
        {
            rep.failures.push(format!(
                "observed order {o:.3} on the finest pair, outside [{lo}, {hi}]"
            ));
        }
        let pts: Vec<(f64, f64)> = e.iter().map(|e| (e.h, e.l2)).collect();
        let slope: Vec<(f64, f64)> = match (pts.first(), pts.last()) {
            (Some(&a), Some(&b)) => vec![(a.0, b.1 * (a.0 / b.0).powi(2)), b],
            _ => Vec::new(),
        };
        std::fs::write(
            out.join("convergence.svg"),
            plot::chart(
                &format!("{name}: velocity error over mesh size"),
                "h (m)",
                "L2 error / U",
                &[
                    Series {
                        label: "Foehncraft",
                        points: &pts,
                        markers: true,
                    },
                    Series {
                        label: "slope 2",
                        points: &slope,
                        markers: false,
                    },
                ],
                true,
                true,
            ),
        )?;
    }
    for o in &outs {
        let s = &o.solve;
        if let (Some(max), Some(e)) = (v.max_error, &s.errors)
            && e.l2 > max
        {
            rep.failures.push(format!(
                "{}: relative L2 velocity error {:.3e} above {max:e}",
                s.mesh, e.l2
            ));
        }
        if s.mass_imbalance.abs() > 1e-6 {
            rep.failures.push(format!(
                "{}: net mass flow out of the domain {:.2e} of the largest boundary flow",
                s.mesh, s.mass_imbalance
            ));
        }
        if let (Some(tol), Some(r)) = (v.reference_tolerance, &s.reference)
            && r.u_max_dev.max(r.v_max_dev) > tol
        {
            rep.failures.push(format!(
                "{}: deviation from {} of {:.4} above {tol}",
                s.mesh,
                r.reference,
                r.u_max_dev.max(r.v_max_dev)
            ));
        }
        if let (Some(tol), Some(b)) = (v.benchmark_tolerance, &s.benchmark)
            && b.max_deviation() > tol
        {
            rep.failures.push(format!(
                "{}: largest distance from the {} ranges {:.4} above {tol}",
                s.mesh,
                b.benchmark,
                b.max_deviation()
            ));
        }
        if v.benchmark == Some(Benchmark::SchaferTurek2d2) && s.benchmark.is_none() && !o.stopped {
            rep.failures.push(format!(
                "{}: fewer than two full periods of the lift in the evaluation window",
                s.mesh
            ));
        }
    }
    rep.solves = outs.into_iter().map(|o| o.solve).collect();
    rep.seconds = t0.elapsed().as_secs_f64();
    for f in &rep.failures {
        lg.line(&format!("FAILED: {f}"));
    }
    std::fs::write(out.join("flow.log"), &lg.text)?;
    std::fs::write(out.join("flow.json"), serde_json::to_string_pretty(&rep)?)?;
    std::fs::write(out.join("flow.md"), rep.to_markdown())?;
    Ok(rep)
}

fn header<F: FnMut(&str)>(
    lg: &mut Log<F>,
    case: &FlowCase,
    mesh: &foehncraft_mesh::Mesh,
    desc: &str,
    last: bool,
    dt: f64,
) {
    let t = case.time.unwrap_or_default();
    let s = &case.solver;
    lg.line(&format!(
        "Foehncraft {} flow solver (PISO, transient, laminar)",
        env!("CARGO_PKG_VERSION")
    ));
    lg.line(&format!(
        "fluid  ρ = {} kg/m³, μ = {} Pa s; convection {}",
        case.fluid.density,
        case.fluid.viscosity,
        if case.flow.convection {
            "on"
        } else {
            "off (Stokes)"
        }
    ));
    let steps = t.steps.unwrap_or_else(|| {
        t.end_time
            .map_or(0, |e| ((e - t.start_time) / dt).ceil() as usize)
    });
    lg.line(&format!(
        "time   {}, Δt = {dt:e} s{}, from t = {} s {}; {} corrector(s), {} outer iteration(s)",
        t.scheme.label(),
        if t.adaptive {
            format!(" (adaptive, Courant number {})", t.max_cfl)
        } else {
            String::new()
        },
        t.start_time,
        match (t.end_time, t.steps) {
            (Some(e), _) => format!("to t = {e} s"),
            (None, Some(n)) => format!("for {n} steps"),
            _ => String::new(),
        },
        s.correctors,
        s.outer_iterations
    ));
    if last {
        let nc = mesh.n_cells;
        let o = &case.output;
        lg.line(&format!("output for {desc} (about {steps} steps):"));
        lg.line(&format!(
            "  history.csv: {} rows, about {} kB",
            steps / o.history_every + 1,
            steps / o.history_every * 120 / 1000 + 1
        ));
        let walls = case
            .boundary
            .values()
            .filter(|b| matches!(b, BoundaryCondition::Wall { .. }))
            .count();
        if walls > 0 || o.coefficients.is_some() {
            lg.line(&format!(
                "  forces.csv: {} rows (every {} steps), forces on {walls} wall(s){}",
                steps / o.monitor_every + 1,
                o.monitor_every,
                if o.coefficients.is_some() {
                    " and the force coefficients"
                } else {
                    ""
                }
            ));
        }
        if !o.probes.is_empty() {
            lg.line(&format!(
                "  probes.csv: {} probe(s), {} rows",
                o.probes.len(),
                steps / o.monitor_every + 1
            ));
        }
        let mb =
            (nc * (8 * 6) + mesh.face_nodes.len() * 8 + mesh.points.len() * 24) / 1_000_000 + 1;
        if o.vtkhdf {
            let series = steps.checked_div(o.fields_every).unwrap_or(0) + o.fields_at.len();
            lg.line(&format!(
                "  fields.vtkhdf (final){}: about {mb} MB each",
                if series > 0 {
                    format!(" and {series} fields_<step>.vtkhdf")
                } else {
                    String::new()
                }
            ));
        }
        let r = &case.run;
        if r.restart_final || r.restart_every > 0 || !r.restart_at.is_empty() {
            lg.line(&format!(
                "  restart.h5{}: about {} MB",
                if r.restart_at.is_empty() {
                    String::new()
                } else {
                    format!(" and restart_<step>.h5 at {:?}", r.restart_at)
                },
                (nc * 8 * 22 + mesh.n_faces() * 8 * 4) / 1_000_000 + 1
            ));
        }
        lg.line("  patches.csv, flow.json, flow.md, case.resolved.toml");
        lg.line(&format!(
            "  control files in this directory: .stop, .result_N, .read_toml_N{}",
            if r.reread_every > 0 {
                format!("; the case file is re-read every {} steps", r.reread_every)
            } else {
                String::new()
            }
        ));
    }
}

/// The restart content of the current state.
fn restart_of(
    fl: &mut Flow,
    tr: &Transient,
    mesh: &foehncraft_mesh::Mesh,
    name: &str,
    dt_next: f64,
) -> Result<Restart> {
    Ok(Restart {
        case: name.into(),
        cell_id: mesh.cell_id[..mesh.n_cells].to_vec(),
        face_id: mesh.face_id.clone(),
        state: fl.state()?,
        iteration: fl.iteration(),
        amg_setups: fl.amg_setups,
        amg_matrix: fl.amg_matrix.clone(),
        transient: Some((tr.old_levels(fl)?, tr.clock, dt_next, tr.scheme)),
    })
}

#[allow(clippy::too_many_arguments)]
fn run_one<F: FnMut(&str)>(
    case0: &FlowCase,
    name: &str,
    out: &Path,
    lg: &mut Log<F>,
    control: &Control,
    steer: &mut Steer,
    warnings: &mut Vec<crate::diag::Diagnostic>,
    mesh: &foehncraft_mesh::Mesh,
    geo: &foehncraft_mesh::Geometry,
    desc: &str,
    dt0: f64,
    last: bool,
    single: bool,
) -> Result<RunOut> {
    let mut case = case0.clone();
    let t_setup = Instant::now();
    let phys = Physics {
        rho: case.fluid.density,
        mu: case.fluid.viscosity,
        body_force: case.flow.body_force,
        convection: case.flow.convection,
    };
    let time0 = case.time.context("[time]")?;
    let analytic0 = case.analytic_at(time0.start_time);
    let arrays = bc::arrays(mesh, geo, &case.boundary, phys.rho, analytic0.as_ref())?;
    let ub0 = arrays.ub.clone();
    let disc = simplec::discretisation(mesh, geo, &arrays)?;
    let loc = Locator::new(mesh, geo);
    let reference_cell = if arrays.has_pressure() {
        None
    } else {
        let x = case.pressure.reference_point.unwrap_or_else(|| {
            let (lo, hi) = mesh.bounding_box();
            std::array::from_fn(|i| lo[i] + 1e-6 * (hi[i] - lo[i]))
        });
        Some(locate(&loc, geo, x))
    };
    let target = Target {
        device: case.run.backend.clone(),
        threads: case.run.threads,
        ptx: case.run.ptx.clone(),
    };
    let mut fl = Flow::on(
        disc,
        arrays,
        phys,
        case.solver.clone(),
        true,
        reference_cell,
        &target,
    )?;
    let nc = mesh.n_cells;
    let mut tr = Transient::new(&mut fl, time0.scheme, time0.start_time)?;
    let mut af = AnalyticFaces::new(&case, mesh, &ub0);
    let mut dt = dt0;
    let mut restart_from = None;
    let restart_file = if single {
        case.run.restart.clone()
    } else {
        None
    };
    if let Some(path) = &restart_file {
        let r = restart::read(path)?;
        let (r, same) = r.to_local(&mesh.cell_id[..nc], &mesh.face_id, mesh.n_interior_faces())?;
        let Some((old, clock, dt_next, scheme)) = r.transient.clone() else {
            bail!("{} is the restart file of a steady run", path.display());
        };
        if scheme != time0.scheme {
            lg.line(&format!(
                "note: the restart file was written with {}, the case runs {}",
                scheme.label(),
                time0.scheme.label()
            ));
        }
        fl.load(&r.state)?;
        tr.restore(&mut fl, &old, clock)?;
        fl.iteration = r.iteration;
        if let Some(m) = &r.amg_matrix {
            fl.restore_amg(m, r.amg_setups)?;
        }
        if !same {
            let m = "note: the restart file has another cell or face order; the fields were mapped \
                 by identity and the multigrid rebuilt (the continuation is not bit-identical)";
            lg.line(m);
            warnings.push(crate::diag::Diagnostic::new(
                crate::diag::RUN_RESTART_MAPPED,
                m,
            ));
        }
        if !time0.adaptive {
            dt = dt0;
        } else {
            dt = dt_next;
        }
        lg.line(&format!(
            "restart from {} at step {}, t = {:.6e} s",
            path.display(),
            clock.step,
            clock.time
        ));
        restart_from = Some(path.display().to_string());
    } else {
        let (u, p) = match analytic0.filter(|a| a.is_transient()) {
            Some(a) => {
                let mut u = Soa3::zeros(nc);
                let mut p = vec![0.0; nc];
                for (c, pc) in p.iter_mut().enumerate() {
                    let x = geo.cell_centroid[c];
                    u.set(c, a.velocity(x));
                    *pc = a.pressure(x);
                }
                (u, p)
            }
            None => {
                let u0 = case.initial.velocity;
                (
                    Soa3 {
                        x: vec![u0[0]; nc],
                        y: vec![u0[1]; nc],
                        z: vec![u0[2]; nc],
                    },
                    vec![
                        if reference_cell.is_some() {
                            case.pressure.reference_value
                        } else {
                            case.initial.pressure
                        };
                        nc
                    ],
                )
            }
        };
        af.update(&case, mesh, geo, &mut fl, time0.start_time)?;
        fl.initialise(&u, &p)?;
        tr.start(&mut fl)?;
    }
    let setup = t_setup.elapsed().as_secs_f64();
    let first_step = tr.clock.step;

    // monitors
    let bench = case.verification.benchmark;
    let extra = benchmark_points(&case, mesh);
    let probes = Probes::new(&case, &extra, &loc, geo);
    let cont = restart_from.as_ref().map(|_| first_step);
    let hist_dir = out;
    let mut history = History::new(
        hist_dir,
        &suffix_name("history", desc, last),
        &crate::case::csv_preamble(
            "history",
            name,
            None,
            "step -, time s, dt s, cfl - (largest Courant number), momentum_x momentum_y \
             momentum_z continuity - (scaled residuals), outer cg - (iterations)",
        ),
        "step,time,dt,cfl,momentum_x,momentum_y,momentum_z,continuity,outer,cg",
        cont,
    )?;
    let walls: Vec<String> = mesh
        .patches
        .iter()
        .filter(|p| {
            matches!(
                case.boundary.get(&p.name),
                Some(BoundaryCondition::Wall { .. })
            )
        })
        .map(|p| p.name.clone())
        .collect();
    let mut fh = String::from("step,time");
    for w in &walls {
        let _ = write!(fh, ",{w}_fx,{w}_fy,{w}_fz");
    }
    if case.output.coefficients.is_some() {
        fh.push_str(",cd,cl");
    }
    let monitors = last && (!walls.is_empty() || !probes.names.is_empty());
    let mut forces = if monitors {
        Some(History::new(
            out,
            "forces.csv",
            &crate::case::csv_preamble(
                "forces",
                name,
                None,
                "step -, time s, <wall>_fx <wall>_fy <wall>_fz N (force of the fluid on the wall), cd cl -",
            ),
            &fh,
            cont,
        )?)
    } else {
        None
    };
    let mut ph = String::from("step,time");
    for n in &probes.names {
        let _ = write!(ph, ",{n}_u,{n}_v,{n}_w,{n}_p");
    }
    let mut probe_csv = if monitors && !probes.names.is_empty() {
        Some(History::new(
            out,
            "probes.csv",
            &crate::case::csv_preamble(
                "probes",
                name,
                None,
                "step -, time s, <probe>_u <probe>_v <probe>_w m/s, <probe>_p Pa",
            ),
            &ph,
            cont,
        )?)
    } else {
        None
    };
    // benchmark histories
    let (mut bt, mut bcd, mut bcl, mut bdp) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());

    let ts = Instant::now();
    let mut stopped = false;
    let (mut dt_min, mut dt_max, mut cfl_max) = (f64::INFINITY, 0.0f64, 0.0f64);
    let (mut outer_sum, mut cg_sum, mut nsteps) = (0usize, 0usize, 0usize);
    let mut last_rep = StepReport::default();
    loop {
        let t = case.time.context("[time]")?;
        let step = tr.clock.step;
        if t.steps.is_some_and(|n| step >= n) {
            break;
        }
        let mut h = dt;
        if let Some(end) = t.end_time {
            let rest = end - tr.clock.time;
            if rest <= 1e-9 * dt {
                break;
            }
            if h > rest || rest - h < 1e-6 * h {
                h = rest;
            }
        }
        af.update(&case, mesh, geo, &mut fl, tr.clock.time + h)?;
        let rp = tr.step(&mut fl, h)?;
        last_rep = rp;
        nsteps += 1;
        outer_sum += rp.outer;
        cg_sum += rp.residuals.pressure_iterations;
        dt_min = dt_min.min(h);
        dt_max = dt_max.max(h);
        cfl_max = cfl_max.max(rp.cfl);
        let n = rp.step;
        let r = &rp.residuals;
        if n.is_multiple_of(case.output.history_every) {
            history.row(&format!(
                "{n},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{},{}",
                rp.time,
                h,
                rp.cfl,
                r.momentum[0],
                r.momentum[1],
                r.momentum[2],
                r.continuity,
                rp.outer,
                r.pressure_iterations
            ));
            if steer.record.on() {
                steer.record.row(&serde_json::json!({
                    "type": "step",
                    "step": n,
                    "time": rp.time,
                    "dt": h,
                    "cfl": rp.cfl,
                    "momentum": r.momentum,
                    "continuity": r.continuity,
                    "outer": rp.outer,
                    "cg": r.pressure_iterations,
                }));
            }
        }
        if n == first_step + 1 || n.is_multiple_of(case.solver.log_every) {
            steer.record.flush();
            lg.line(&format!(
                "step {n:6}  t {:.6e}  Δt {h:.3e}  Courant {:.3}  momentum {:.2e} {:.2e} {:.2e}  \
                 continuity {:.2e}  outer {}  cg {}",
                rp.time,
                rp.cfl,
                r.momentum[0],
                r.momentum[1],
                r.momentum[2],
                r.continuity,
                rp.outer,
                r.pressure_iterations
            ));
            // the histories reach the disk with every log line, so a run
            // that is killed keeps them up to there
            history.flush()?;
            if let Some(f) = forces.as_mut() {
                f.flush()?;
            }
            if let Some(f) = probe_csv.as_mut() {
                f.flush()?;
            }
        }
        let req = control.check(n, case.run.reread_every);
        let at_monitor = monitors && n.is_multiple_of(case.output.monitor_every);
        let o = &case.output;
        let fields_now = last
            && ((o.fields_every > 0 && n.is_multiple_of(o.fields_every))
                || o.fields_at.contains(&n)
                || req.result);
        let r = &case.run;
        let asked = steer.restart_requested();
        let restart_now = last
            && ((r.restart_every > 0 && n.is_multiple_of(r.restart_every))
                || r.restart_at.contains(&n)
                || asked);
        if at_monitor || fields_now {
            let sol = Solution {
                state: fl.state()?,
                disc: &fl.disc,
                phys: fl.phys,
            };
            if at_monitor {
                let patches = patch_summaries(&sol, mesh, geo, &case);
                let mut row = format!("{n},{:e}", rp.time);
                for w in &walls {
                    if let Some(p) = patches.iter().find(|p| &p.name == w) {
                        let _ = write!(row, ",{:e},{:e},{:e}", p.force[0], p.force[1], p.force[2]);
                    }
                }
                let cc = monitor::coefficients(&case, &patches);
                if let Some([cd, cl]) = cc {
                    let _ = write!(row, ",{cd:e},{cl:e}");
                }
                if let Some(fh) = forces.as_mut() {
                    fh.row(&row);
                }
                let vals = probes.values(&sol, geo);
                if let Some(pc) = probe_csv.as_mut() {
                    let mut row = format!("{n},{:e}", rp.time);
                    for (u, p) in &vals {
                        let _ = write!(row, ",{:e},{:e},{:e},{p:e}", u[0], u[1], u[2]);
                    }
                    pc.row(&row);
                }
                if let (Some(_), Some([cd, cl])) =
                    (bench, benchmark_coefficients(&case, &patches, mesh))
                {
                    let (i, j) = (
                        probes.index("benchmark_front").unwrap_or(0),
                        probes.index("benchmark_back").unwrap_or(0),
                    );
                    bt.push(rp.time);
                    bcd.push(cd);
                    bcl.push(cl);
                    bdp.push(vals[i].1 - vals[j].1);
                }
            }
            if fields_now {
                write_fields(
                    &case,
                    &out.join(format!("fields_{n:06}.vtkhdf")),
                    mesh,
                    geo,
                    &sol,
                    case.analytic_at(rp.time).as_ref(),
                    Some((n, Some(rp.time))),
                )?;
                lg.line(&format!(
                    "      fields of step {n} written (fields_{n:06}.vtkhdf)"
                ));
            }
        }
        // next time step, from the step just taken (h, which is shorter
        // than dt when it was cut to end at end_time): the Courant number
        // scales with the step, and the growth limit holds for the ratio
        // of two successive steps
        let tnow = case.time.context("[time]")?;
        if tnow.adaptive {
            dt = h;
            let grow = if rp.cfl > 0.0 {
                (tnow.max_cfl / rp.cfl).min(tnow.max_growth)
            } else {
                tnow.max_growth
            };
            dt = (dt * grow).max(tnow.min_time_step);
            if let Some(m) = tnow.max_time_step {
                dt = dt.min(m);
            }
        }
        if restart_now && restart::AVAILABLE {
            let rs = restart_of(&mut fl, &tr, mesh, name, dt)?;
            restart::write(&out.join("restart.h5"), &rs)?;
            if case.run.restart_at.contains(&n) {
                restart::write(&out.join(format!("restart_{n:06}.h5")), &rs)?;
            }
            if asked {
                lg.line(&format!("      step {n}: restart.h5 written (SIGUSR1)"));
            }
        }
        if req.reread
            && let Some(old) =
                crate::case::reread_case(&mut case, &mut fl, lg, steer, warnings, n, out)?
            && let (Some(a), Some(b)) = (old.time, case.time)
            && a.time_step != b.time_step
            && !b.adaptive
        {
            dt = b.time_step;
        }
        if let Some(s) = steer.poll(&case.run.limits, nsteps, true, req.stop) {
            lg.line(&format!("      step {n}: {}", s.message));
            stopped = true;
            break;
        }
    }
    history.flush()?;
    if let Some(f) = forces.as_mut() {
        f.flush()?;
    }
    if let Some(f) = probe_csv.as_mut() {
        f.flush()?;
    }
    let secs = ts.elapsed().as_secs_f64();
    let steps_total = tr.clock.step;
    lg.line(&format!(
        "      {} {nsteps} steps to t = {:.6e} s in {secs:.2} s (setup {setup:.2} s), {:.0} ns/cell/step; \
         Δt {:.3e}–{:.3e} s, largest Courant number {:.3}",
        if stopped { "stopped after" } else { "done:" },
        tr.clock.time,
        1e9 * secs / (nc as f64 * nsteps.max(1) as f64),
        if nsteps > 0 { dt_min } else { 0.0 },
        dt_max,
        cfl_max
    ));
    if last && (case.run.restart_final || stopped) && restart::AVAILABLE {
        let rs = restart_of(&mut fl, &tr, mesh, name, dt)?;
        restart::write(&out.join("restart.h5"), &rs)?;
    }
    let sol = Solution {
        state: fl.state()?,
        disc: &fl.disc,
        phys: fl.phys,
    };
    let analytic_end = case.analytic_at(tr.clock.time);
    let errors = analytic_end.map(|a| error_norms(&sol, mesh, geo, &a, dimension(&case)));
    if let Some(e) = &errors {
        lg.line(&format!(
            "      at t = {:.6e} s: velocity error L2 {:.3e}, Linf {:.3e} (relative), pressure Linf {:.3e}",
            tr.clock.time, e.l2, e.linf, e.pressure_linf
        ));
    }
    let reference = case
        .verification
        .reference
        .map(|r| crate::case::compare_reference(&sol, mesh, geo, &loc, r, &case));
    if let Some(r) = &reference {
        lg.line(&format!(
            "      {}: largest deviation u {:.4}, v {:.4} (lid speeds)",
            r.reference, r.u_max_dev, r.v_max_dev
        ));
        if last {
            crate::case::write_reference_svg(out, name, r)?;
        }
    }
    let patches = patch_summaries(&sol, mesh, geo, &case);
    let total: f64 = patches.iter().map(|p| p.mass_flow).sum();
    let umax = (0..nc).map(|c| norm(sol.u.get(c))).fold(0.0f64, f64::max);
    let big = patches
        .iter()
        .map(|p| p.mass_flow.abs().max(case.fluid.density * umax * p.area))
        .fold(0.0f64, f64::max);
    let mass_imbalance = if big > 0.0 { total / big } else { 0.0 };
    let coefficients = monitor::coefficients(&case, &patches);
    if let Some([cd, cl]) = coefficients {
        lg.line(&format!(
            "      force coefficients C_D {cd:.5}, C_L {cl:.5}"
        ));
    }
    // the chart of the coefficients shows the evaluation window of a periodic
    // benchmark (the impulsive start is orders of magnitude larger)
    let mut chart_start = f64::NEG_INFINITY;
    let benchmark = match bench {
        Some(b @ Benchmark::SchaferTurek2d1) => {
            benchmark_coefficients(&case, &patches, mesh).map(|[cd, cl]| {
                let vals = probes.values(&sol, geo);
                let (i, j) = (
                    probes.index("benchmark_front").unwrap_or(0),
                    probes.index("benchmark_back").unwrap_or(0),
                );
                monitor::steady_benchmark(b, cd, cl, vals[i].1 - vals[j].1)
            })
        }
        Some(b @ Benchmark::SchaferTurek2d2) => {
            let start = case
                .verification
                .evaluation_start
                .unwrap_or(0.5 * (time0.start_time + tr.clock.time));
            chart_start = start;
            monitor::periodic_benchmark(b, &bt, &bcd, &bcl, &bdp, start)
        }
        None => None,
    };
    if let Some(b) = &benchmark {
        log_benchmark(b, |l| lg.line(l));
    }
    if last || stopped {
        write_outputs(
            &case,
            out,
            mesh,
            geo,
            &sol,
            analytic_end.as_ref(),
            &patches,
            &loc,
            Some((tr.clock.step, Some(tr.clock.time))),
        )?;
        if !bt.is_empty() {
            let window = |c: &[f64]| -> Vec<(f64, f64)> {
                bt.iter()
                    .copied()
                    .zip(c.iter().copied())
                    .filter(|&(t, _)| t >= chart_start)
                    .collect()
            };
            let (cd, cl) = (window(&bcd), window(&bcl));
            let title = if chart_start.is_finite() {
                format!("{name}: force coefficients from t = {chart_start} s")
            } else {
                format!("{name}: force coefficients")
            };
            std::fs::write(
                out.join("coefficients.svg"),
                plot::chart(
                    &title,
                    "t (s)",
                    "C",
                    &[
                        Series {
                            label: "C_D",
                            points: &cd,
                            markers: false,
                        },
                        Series {
                            label: "C_L",
                            points: &cl,
                            markers: false,
                        },
                    ],
                    false,
                    false,
                ),
            )?;
        }
    }
    let solve = Solve {
        mesh: desc.into(),
        mesh_hash: crate::engine::mesh_hash(mesh),
        device: fl.device(),
        cells: nc,
        iterations: nsteps,
        converged: !stopped,
        residuals: last_rep.residuals,
        ns_per_cell_iteration: 1e9 * secs / (nc as f64 * nsteps.max(1) as f64),
        timings: fl.timings.clone(),
        amg_setups: fl.amg_setups,
        errors,
        reference,
        patches,
        mass_imbalance,
        coefficients,
        transient: Some(TransientSummary {
            scheme: time0.scheme.label().into(),
            steps: nsteps,
            last_step: steps_total,
            time: tr.clock.time,
            dt_min: if nsteps > 0 { dt_min } else { 0.0 },
            dt_max,
            cfl_max,
            outer_mean: outer_sum as f64 / nsteps.max(1) as f64,
            cg_mean: cg_sum as f64 / nsteps.max(1) as f64,
            restart_from,
            stopped,
            ns_per_cell_step: 1e9 * secs / (nc as f64 * nsteps.max(1) as f64),
        }),
        benchmark,
    };
    let volume = geo.cell_volume.clone();
    Ok(RunOut {
        u: sol.state.u.clone(),
        solve,
        volume,
        stopped,
    })
}

/// C_D and C_L of the benchmark cylinder: 2 F / (ρ Ū² D H) with H the depth
/// of the slab.
/// The front and back points of the benchmark (Δp), as extra probes.
pub(crate) fn benchmark_points(
    case: &FlowCase,
    mesh: &foehncraft_mesh::Mesh,
) -> Vec<(String, [f64; 3])> {
    let (lo, hi) = mesh.bounding_box();
    let zc = 0.5 * (lo[2] + hi[2]);
    match case.verification.benchmark {
        Some(_) => vec![
            (
                "benchmark_front".into(),
                [Benchmark::FRONT[0], Benchmark::FRONT[1], zc],
            ),
            (
                "benchmark_back".into(),
                [Benchmark::BACK[0], Benchmark::BACK[1], zc],
            ),
        ],
        None => Vec::new(),
    }
}

/// The steady benchmark (2D-1) of a converged solution.
pub(crate) fn steady_benchmark_of(
    case: &FlowCase,
    sol: &Solution,
    mesh: &foehncraft_mesh::Mesh,
    geo: &foehncraft_mesh::Geometry,
    loc: &Locator,
    patches: &[crate::case::PatchSummary],
) -> Option<monitor::BenchmarkResult> {
    let b = case.verification.benchmark?;
    let [cd, cl] = benchmark_coefficients(case, patches, mesh)?;
    let probes = Probes::new(case, &benchmark_points(case, mesh), loc, geo);
    let vals = probes.values(sol, geo);
    let (i, j) = (
        probes.index("benchmark_front")?,
        probes.index("benchmark_back")?,
    );
    Some(monitor::steady_benchmark(b, cd, cl, vals[i].1 - vals[j].1))
}

/// Log lines of a benchmark result.
pub(crate) fn log_benchmark(b: &monitor::BenchmarkResult, mut line: impl FnMut(&str)) {
    line(&format!("      {}:", b.benchmark));
    for v in &b.values {
        line(&format!(
            "        {:10} {:.5}  reference [{}, {}]  distance {:.4}",
            v.range.name, v.value, v.range.lower, v.range.upper, v.deviation
        ));
    }
}

fn benchmark_coefficients(
    case: &FlowCase,
    patches: &[crate::case::PatchSummary],
    mesh: &foehncraft_mesh::Mesh,
) -> Option<[f64; 2]> {
    let b = case.verification.benchmark?;
    let name = case
        .verification
        .benchmark_patch
        .clone()
        .unwrap_or_else(|| "cylinder".into());
    let p = patches.iter().find(|p| p.name == name)?;
    let (lo, hi) = mesh.bounding_box();
    let depth = hi[2] - lo[2];
    let u = b.mean_velocity();
    let q = 0.5 * case.fluid.density * u * u * Benchmark::DIAMETER * depth;
    Some([p.force[0] / q, p.force[1] / q])
}

fn suffix_name(base: &str, _desc: &str, last: bool) -> String {
    if last {
        format!("{base}.csv")
    } else {
        format!("{base}_run.csv")
    }
}

/// The time scheme label (for reports).
pub fn scheme_label(s: TimeScheme) -> &'static str {
    s.label()
}

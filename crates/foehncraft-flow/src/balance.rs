//! Conservation ledger (D88), balance reports (D89) and the kinetic-energy
//! form of the second-law check (D90 part 4).
//!
//! The ledger books, for each conserved quantity, the storage (the time
//! derivative of its integral over the domain), the boundary fluxes per
//! patch and the volume sources, and closes it:
//!
//! imbalance = storage − Σ boundary fluxes into the domain − Σ sources.
//!
//! **Mass** (kg/s). The density is constant, so the storage is zero and
//! the imbalance is the net mass flow out of the domain. The interior face
//! fluxes cancel exactly, so it is the sum of the continuity residuals of
//! all cells: a check of the pressure-correction tolerance.
//!
//! **Kinetic energy** (W). The momentum equation dotted with u gives, for
//! constant ρ and μ and the Laplacian form of the viscous term that the
//! solver discretises,
//!
//! d/dt ∫ ½ρ|u|² dV = Σ_patches [ −∮ ½ρ|u|² u·dS − ∮ p u·dS + ∮ (μ∇u·dS)·u ]
//!                    + ∫ f·u dV − ∫ μ ∇u:∇u dV,
//!
//! boundary power (convected kinetic energy, pressure work, viscous work),
//! the power of the body force and the viscous dissipation Φ = μ ∇u:∇u,
//! which is never negative. The discrete solution satisfies this only up
//! to the dissipation of the schemes (upwind-biased convection, the
//! interpolations, the pressure–velocity coupling) and, in transient runs,
//! of the time integration. The ledger computes every physical term from
//! the solution and books the remainder as the **numerical dissipation**;
//! the second law asks that it is not negative (D90 part 4), and it falls
//! as the mesh is refined.
//!
//! The ledger splits it in two. The momentum rows hold the pressure as
//! the volume source −∇p_P V_P with the least-squares cell gradient, so
//! their energy budget has −Σ u_P·∇p_P V_P where the exact one has the
//! boundary pressure work; the difference, the discrete −∫ p ∇·u dV of the
//! cell velocities, is the **pressure part**. The momentum interpolation
//! of the face fluxes (Rhie–Chow) makes it a dissipation on regular cells;
//! on distorted cells it holds a consistency error of either sign that
//! falls with the mesh. The rest is the **scheme part**, of the convection
//! and time schemes: positive for upwind convection, of either sign for
//! the linear-upwind correction near singular points.
//!
//! **Mass in transient runs.** The discrete continuity equation holds at
//! the new time level, so the mass flows of a step are those at its end.
//!
//! In a transient run the storage of kinetic energy over a step is
//! (Kⁿ⁺¹ − Kⁿ)/Δt and its rates are the mean of their values at tⁿ and tⁿ⁺¹ (trapezoidal rule,
//! second order), so that the cumulative ledger of the run is
//! ΔK = ∫ (boundary power + body-force power − Φ − numerical dissipation) dt.
//!
//! All sums are pairwise sums over the cells or faces in a fixed order
//! (D31): the ledger is bitwise reproducible and does not depend on the
//! thread count or the device. It is computed from the fields on the host
//! and changes nothing of the solution.
//!
//! The volumes are passed per evaluation (D91): with a deforming mesh
//! (L22) each time level brings its own; today they are the volumes of the
//! fixed mesh.

use crate::case::{FlowCase, wall_face_values};
use crate::diag::{self, Diagnostic};
use crate::simplec::Solution;
use foehncraft_core::repro::pairwise_sum;
use foehncraft_kernels::flow::{MOM_FIXED, MOM_SYMMETRY};
use foehncraft_mesh::geometry::{dot, norm};
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};

/// What an exceeded tolerance does.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum BalanceAction {
    /// A warning with a code (D78) in the log and the result.
    #[default]
    Warn,
    /// Stop the run as at a limit: restart file, final output, exit code 3.
    Stop,
}

/// `[balance]` (D89).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct BalanceSection {
    /// The full balance table in the log every N outer iterations (steady)
    /// or time steps (transient); 0: at the end only. A transient run also
    /// gives one compact balance line with every progress line.
    pub every: usize,
    /// Largest relative mass imbalance of one time step, or of the final
    /// state of a steady run (-).
    pub tolerance_step: f64,
    /// Largest relative mass imbalance over the run (-).
    pub tolerance_run: f64,
    /// Largest negative numerical dissipation of kinetic energy, relative
    /// to the largest term of the ledger, of the final steady state or over
    /// the run (-). Below it the second-law check warns; a transient step
    /// below it counts as negative. The default leaves room for the
    /// consistency errors of distorted cells, which fall with the mesh.
    pub tolerance_energy: f64,
    /// What an exceeded tolerance does: `warn`, or `stop` (a transient run
    /// stops at the step as at a limit, exit code 3; a steady run checks
    /// its final state only, so `stop` acts as `warn` there).
    pub action: BalanceAction,
}

impl Default for BalanceSection {
    fn default() -> Self {
        Self {
            every: 0,
            tolerance_step: 1e-8,
            tolerance_run: 1e-4,
            tolerance_energy: 1e-2,
            action: BalanceAction::Warn,
        }
    }
}

impl BalanceSection {
    /// Checks the ranges.
    pub fn validate(&self) -> anyhow::Result<()> {
        crate::need!(
            self.tolerance_step > 0.0,
            CASE_RANGE,
            "[balance] tolerance_step must be > 0"
        );
        crate::need!(
            self.tolerance_run > 0.0,
            CASE_RANGE,
            "[balance] tolerance_run must be > 0"
        );
        crate::need!(
            self.tolerance_energy > 0.0,
            CASE_RANGE,
            "[balance] tolerance_energy must be > 0"
        );
        Ok(())
    }
}

/// The boundary terms of one patch.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct PatchTerms {
    /// Patch name.
    pub patch: String,
    /// Boundary condition.
    pub kind: String,
    /// Mass flow out of the domain (kg/s).
    pub mass_flow: f64,
    /// Kinetic energy carried out of the domain, Σ ṁ ½|u_b|² (W).
    pub convective_power: f64,
    /// Work of the pressure on the outflow, Σ p_b u_b·S (W).
    pub pressure_power: f64,
    /// Work of the viscous stress on the fluid, Σ (μ ∇u·S)·u_b (W).
    pub viscous_power: f64,
}

impl PatchTerms {
    /// Power into the fluid: viscous work − convected kinetic energy −
    /// pressure work (W).
    pub fn power(&self) -> f64 {
        self.viscous_power - self.convective_power - self.pressure_power
    }

    fn scaled(&self, a: f64) -> Self {
        Self {
            patch: self.patch.clone(),
            kind: self.kind.clone(),
            mass_flow: a * self.mass_flow,
            convective_power: a * self.convective_power,
            pressure_power: a * self.pressure_power,
            viscous_power: a * self.viscous_power,
        }
    }

    fn add(&mut self, o: &Self, a: f64) {
        self.mass_flow += a * o.mass_flow;
        self.convective_power += a * o.convective_power;
        self.pressure_power += a * o.pressure_power;
        self.viscous_power += a * o.viscous_power;
    }
}

/// The instantaneous terms of one state.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Rates {
    /// Kinetic energy in the domain, ∫ ½ρ|u|² dV (J).
    pub kinetic_energy: f64,
    /// Boundary terms per patch.
    pub patches: Vec<PatchTerms>,
    /// Power of the body force, ∫ f·u dV (W).
    pub body_force_power: f64,
    /// Power of the pressure source of the momentum rows on the cell
    /// velocities, −Σ u_P·∇p_P V_P (W). With the boundary pressure work it
    /// gives the pressure part of the numerical dissipation.
    pub pressure_source_power: f64,
    /// Viscous dissipation: of the discrete viscous operator, which
    /// converges to ∫ μ ∇u:∇u dV (W).
    pub viscous_dissipation: f64,
    /// Scale of the mass flows: the largest patch mass flow, or in a closed
    /// domain ρ |u|max times the largest patch area (kg/s).
    pub mass_scale: f64,
}

/// The terms of a state on the host. `volume` holds the cell volumes of
/// the time level of the state (D91); `nonorth` is the non-orthogonal
/// correction of the solver, which the viscous boundary flux follows.
pub fn rates(
    sol: &Solution,
    mesh: &Mesh,
    geo: &Geometry,
    case: &FlowCase,
    mom_kind: &[u32],
    volume: &[f64],
    nonorth: bool,
) -> Rates {
    let nc = mesh.n_cells;
    let ni = mesh.n_interior_faces();
    let (rho, mu) = (sol.phys.rho, sol.phys.mu);
    let f = sol.phys.body_force;
    let mut ke = vec![0.0; nc];
    let mut body = vec![0.0; nc];
    let mut press = vec![0.0; nc];
    let mut umax = 0.0f64;
    for c in 0..nc {
        let u = sol.u.get(c);
        let v = volume[c];
        let uu = dot(u, u);
        umax = umax.max(uu.sqrt());
        ke[c] = 0.5 * rho * uu * v;
        body[c] = dot(f, u) * v;
        // the pressure source of the momentum rows, −∇p_P V_P
        // (foehncraft_kernels::flow::momentum_row)
        press[c] = -dot(sol.grad_p.get(c), u) * v;
    }
    let d = sol.disc;
    let on = if nonorth { 1.0 } else { 0.0 };
    // the dissipation of the discrete viscous operator, face by face:
    // F_f·(u_N − u_P) with the diffusive flux F_f of the momentum rows
    // (foehncraft_kernels::flow::momentum_face, non-orthogonal and
    // skewness corrections included); it converges to ∫ μ ∇u:∇u dV and
    // is what the discrete momentum equation dotted with u removes
    let mut phi = vec![0.0; d.n_faces];
    for (fc, ph) in phi.iter_mut().enumerate().take(ni) {
        let (p, n) = (d.owner[fc] as usize, d.neighbour[fc] as usize);
        let (up, un) = (sol.u.get(p), sol.u.get(n));
        let w = d.weight[fc];
        let xf = d.face_centroid.get(fc);
        let (xp, xn) = (d.cell_centroid.get(p), d.cell_centroid.get(n));
        let rp = [xf[0] - xp[0], xf[1] - xp[1], xf[2] - xp[2]];
        let rm = [
            xf[0] - 0.5 * (xn[0] + xf[0] - rp[0]),
            xf[1] - 0.5 * (xn[1] + xf[1] - rp[1]),
            xf[2] - 0.5 * (xn[2] + xf[2] - rp[2]),
        ];
        let k = d.k.get(fc);
        let dd = mu * d.dgeo[fc];
        let mut e = 0.0;
        for i in 0..3 {
            let (gp, gn) = (sol.grad_u[i].get(p), sol.grad_u[i].get(n));
            let gf = [
                w * gp[0] + (1.0 - w) * gn[0],
                w * gp[1] + (1.0 - w) * gn[1],
                w * gp[2] + (1.0 - w) * gn[2],
            ];
            let skew = dot([gn[0] - gp[0], gn[1] - gp[1], gn[2] - gp[2]], rm);
            let flux = dd * (un[i] - up[i]) + on * (mu * dot(k, gf) + dd * skew);
            e += flux * (un[i] - up[i]);
        }
        *ph = e;
    }
    let mut patches = Vec::with_capacity(mesh.patches.len());
    let mut scale = 0.0f64;
    for p in &mesh.patches {
        let bc = &case.boundary[&p.name];
        let wall = matches!(bc, crate::bc::BoundaryCondition::Wall { .. });
        let n = p.len;
        let (mut m, mut kc, mut pw, mut vw, mut area) = (
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
        );
        for fc in p.range() {
            let fb = fc - ni;
            let s = geo.face_area[fc];
            let ub = sol.ub.get(fb);
            let mf = sol.mass_flux[fc];
            m.push(mf);
            area.push(norm(s));
            kc.push(mf * 0.5 * dot(ub, ub));
            let pf = if wall {
                wall_face_values(sol, geo, fc, ni).0
            } else {
                sol.pb[fb]
            };
            pw.push(pf * dot(ub, s));
            let kind = mom_kind[fb];
            let c = d.owner[fc] as usize;
            let up = sol.u.get(c);
            let dd = mu * d.dgeo[fc];
            if kind == MOM_FIXED {
                // the diffusive flux of the momentum rows
                // (foehncraft_kernels::flow::momentum_face): work on the
                // fluid F·u_b, dissipation F·(u_b − u_P)
                let k = d.k.get(fc);
                let (mut w, mut e) = (0.0, 0.0);
                for i in 0..3 {
                    let flux = dd * (ub[i] - up[i]) + on * mu * dot(k, sol.grad_u[i].get(c));
                    w += flux * ub[i];
                    e += flux * (ub[i] - up[i]);
                }
                vw.push(w);
                phi[fc] = e;
            } else if kind == MOM_SYMMETRY {
                // the normal velocity goes to zero at the face: no work,
                // dissipation μ d (u_P·n̂)²
                let un = dot(up, s) / norm(s);
                vw.push(0.0);
                phi[fc] = dd * un * un;
            } else {
                // outflow: zero normal gradient, no viscous flux
                vw.push(0.0);
            }
        }
        let t = PatchTerms {
            patch: p.name.clone(),
            kind: bc.type_name().into(),
            mass_flow: pairwise_sum(&m),
            convective_power: pairwise_sum(&kc),
            pressure_power: pairwise_sum(&pw),
            viscous_power: pairwise_sum(&vw),
        };
        scale = scale
            .max(t.mass_flow.abs())
            .max(rho * umax * pairwise_sum(&area));
        patches.push(t);
    }
    Rates {
        kinetic_energy: pairwise_sum(&ke),
        patches,
        body_force_power: pairwise_sum(&body),
        pressure_source_power: pairwise_sum(&press),
        viscous_dissipation: pairwise_sum(&phi),
        mass_scale: scale,
    }
}

/// Mass ledger (kg/s, or kg over a run).
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct MassBalance {
    /// Storage (zero for constant density).
    pub storage: f64,
    /// Net mass flow out of the domain.
    pub outflow: f64,
    /// storage − inflow = storage + outflow.
    pub imbalance: f64,
    /// Imbalance over the scale of the mass flows (-).
    pub relative: f64,
}

/// Kinetic-energy ledger (W, or J over a run).
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct EnergyBalance {
    /// Kinetic energy in the domain at the end (J).
    pub energy: f64,
    /// Storage: its change per time (W) or over the run (J).
    pub storage: f64,
    /// Power of all patches into the fluid.
    pub boundary: f64,
    /// Power of the body force.
    pub body_force: f64,
    /// Viscous dissipation (physical, ≥ 0).
    pub viscous_dissipation: f64,
    /// Numerical dissipation: boundary + body force − viscous dissipation −
    /// storage; ≥ 0 by the second law (D90 part 4). It is the sum of the
    /// two parts below.
    pub numerical_dissipation: f64,
    /// Its pressure part: the boundary pressure work minus the power of the
    /// pressure source −∇p_P V_P of the momentum rows on the cell
    /// velocities, the discrete −∫ p ∇·u dV of the cell velocities. The
    /// momentum interpolation of the face fluxes makes it a dissipation;
    /// it is not sign-definite on distorted cells.
    pub pressure_dissipation: f64,
    /// Its scheme part: of the convection and time schemes (the rest).
    pub scheme_dissipation: f64,
    /// Numerical dissipation over the largest term (-).
    pub relative: f64,
}

/// One closed ledger: of a state (steady), a time step or a whole run.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct Balance {
    /// Outer iteration (steady).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub iteration: Option<usize>,
    /// Time step (transient; the last one for a run).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub step: Option<usize>,
    /// Time at the end (s).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub time: Option<f64>,
    /// Length of the step or of the run (s).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub duration: Option<f64>,
    /// Mass.
    pub mass: MassBalance,
    /// Kinetic energy.
    pub kinetic_energy: EnergyBalance,
    /// Boundary terms per patch (rates; integrals over a run).
    pub patches: Vec<PatchTerms>,
}

#[allow(clippy::too_many_arguments)]
fn close(
    storage_ke: f64,
    energy: f64,
    patches: Vec<PatchTerms>,
    body: f64,
    source: f64,
    phi: f64,
    scale: f64,
) -> (MassBalance, EnergyBalance, Vec<PatchTerms>) {
    let out = pairwise_sum(&patches.iter().map(|p| p.mass_flow).collect::<Vec<_>>());
    let boundary = pairwise_sum(&patches.iter().map(PatchTerms::power).collect::<Vec<_>>());
    let wall_pressure = pairwise_sum(&patches.iter().map(|p| p.pressure_power).collect::<Vec<_>>());
    let pressure = -(source + wall_pressure);
    let num = boundary + body - phi - storage_ke;
    let mut big = storage_ke
        .abs()
        .max(body.abs())
        .max(phi.abs())
        .max(source.abs());
    for p in &patches {
        big = big
            .max(p.viscous_power.abs())
            .max(p.convective_power.abs())
            .max(p.pressure_power.abs());
    }
    (
        MassBalance {
            storage: 0.0,
            outflow: out,
            imbalance: out,
            relative: if scale > 0.0 { out / scale } else { 0.0 },
        },
        EnergyBalance {
            energy,
            storage: storage_ke,
            boundary,
            body_force: body,
            viscous_dissipation: phi,
            numerical_dissipation: num,
            pressure_dissipation: pressure,
            scheme_dissipation: num - pressure,
            relative: if big > 0.0 { num / big } else { 0.0 },
        },
        patches,
    )
}

/// The ledger of a steady state (no storage).
pub fn steady(r: &Rates, iteration: usize) -> Balance {
    let (mass, kinetic_energy, patches) = close(
        0.0,
        r.kinetic_energy,
        r.patches.clone(),
        r.body_force_power,
        r.pressure_source_power,
        r.viscous_dissipation,
        r.mass_scale,
    );
    Balance {
        iteration: Some(iteration),
        mass,
        kinetic_energy,
        patches,
        ..Default::default()
    }
}

/// The ledger of a transient run: one balance per step and the totals.
#[derive(Debug, Clone, Default)]
pub struct Ledger {
    prev: Option<Rates>,
    start_energy: f64,
    start_time: f64,
    time: f64,
    step: usize,
    patches: Vec<PatchTerms>,
    body: f64,
    source: f64,
    phi: f64,
    scale: f64,
    energy: f64,
    /// A step counts as negative below −tolerance (relative).
    pub tolerance: f64,
    /// Steps booked.
    pub steps: usize,
    /// Steps with a negative numerical dissipation (below −tolerance).
    pub negative_steps: usize,
    /// Largest relative mass imbalance of a step.
    pub largest_mass_relative: f64,
}

impl Ledger {
    /// Starts at the state of time `t` (start fields or a restart). A step
    /// counts as negative when its relative numerical dissipation is below
    /// −`tolerance` (round-off is not a finding).
    pub fn start(r: Rates, t: f64, step: usize, tolerance: f64) -> Self {
        Self {
            tolerance,
            start_energy: r.kinetic_energy,
            energy: r.kinetic_energy,
            start_time: t,
            time: t,
            step,
            patches: r.patches.iter().map(|p| p.scaled(0.0)).collect(),
            prev: Some(r),
            ..Default::default()
        }
    }

    /// Books the step to the state `r` at time `t`.
    pub fn step(&mut self, r: Rates, step: usize, t: f64) -> Balance {
        let prev = self.prev.take().unwrap_or_else(|| r.clone());
        let dt = t - self.time;
        let half = 0.5;
        let patches: Vec<PatchTerms> = prev
            .patches
            .iter()
            .zip(&r.patches)
            .map(|(a, b)| {
                let mut m = a.scaled(half);
                m.add(b, half);
                // the discrete continuity equation holds at the new level:
                // the mass flows of the step are those at its end
                m.mass_flow = b.mass_flow;
                m
            })
            .collect();
        let body = half * (prev.body_force_power + r.body_force_power);
        let source = half * (prev.pressure_source_power + r.pressure_source_power);
        let phi = half * (prev.viscous_dissipation + r.viscous_dissipation);
        let scale = r.mass_scale;
        let storage = if dt > 0.0 {
            (r.kinetic_energy - prev.kinetic_energy) / dt
        } else {
            0.0
        };
        let (mass, kinetic_energy, patches) =
            close(storage, r.kinetic_energy, patches, body, source, phi, scale);
        for (acc, p) in self.patches.iter_mut().zip(&patches) {
            acc.add(p, dt);
        }
        self.body += dt * body;
        self.source += dt * source;
        self.phi += dt * phi;
        self.scale += dt * scale;
        self.energy = r.kinetic_energy;
        self.time = t;
        self.step = step;
        self.steps += 1;
        if kinetic_energy.relative < -self.tolerance {
            self.negative_steps += 1;
        }
        self.largest_mass_relative = self.largest_mass_relative.max(mass.relative.abs());
        self.prev = Some(r);
        Balance {
            step: Some(step),
            time: Some(t),
            duration: Some(dt),
            mass,
            kinetic_energy,
            patches,
            ..Default::default()
        }
    }

    /// The ledger of the run so far: integrals over time (kg, J).
    pub fn run(&self) -> Balance {
        let (mass, kinetic_energy, patches) = close(
            self.energy - self.start_energy,
            self.energy,
            self.patches.clone(),
            self.body,
            self.source,
            self.phi,
            self.scale,
        );
        Balance {
            step: Some(self.step),
            time: Some(self.time),
            duration: Some(self.time - self.start_time),
            mass,
            kinetic_energy,
            patches,
            ..Default::default()
        }
    }
}

/// The ledger of a solve in the report and in `result.json`.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct BalanceSummary {
    /// The last balance: of the final state (steady) or of the last step
    /// (transient), in kg/s and W.
    pub last: Balance,
    /// Transient: the whole run, in kg and J.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run: Option<Balance>,
    /// Transient: steps booked.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub steps: Option<usize>,
    /// Transient: steps whose numerical dissipation is negative (a health
    /// number: anti-diffusive steps, D90).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub negative_steps: Option<usize>,
    /// Largest relative mass imbalance of a step or of the final state.
    pub largest_mass_relative: f64,
}

/// The tolerance checks of D89 on a finished ledger: the warnings with
/// their codes.
pub fn check(s: &BalanceSummary, set: &BalanceSection) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    if s.largest_mass_relative > set.tolerance_step {
        out.push(
            Diagnostic::new(
                diag::BAL_MASS_STEP,
                format!(
                    "the relative mass imbalance reached {:.2e}, above [balance] tolerance_step {:.1e} \
                     (tighten [solver] pressure_tolerance)",
                    s.largest_mass_relative, set.tolerance_step
                ),
            )
            .key("balance.tolerance_step"),
        );
    }
    let whole = s.run.as_ref().unwrap_or(&s.last);
    if s.run.is_some() && whole.mass.relative.abs() > set.tolerance_run {
        out.push(
            Diagnostic::new(
                diag::BAL_MASS_RUN,
                format!(
                    "the mass imbalance of the run is {:.2e} of the mass flows, above [balance] \
                     tolerance_run {:.1e}",
                    whole.mass.relative, set.tolerance_run
                ),
            )
            .key("balance.tolerance_run"),
        );
    }
    if whole.kinetic_energy.relative < -set.tolerance_energy {
        out.push(
            Diagnostic::new(
                diag::BAL_NEGATIVE_DISSIPATION,
                format!(
                    "the numerical dissipation of kinetic energy is negative: {:.3e} {u} \
                     ({:.2e} of the largest term, below [balance] tolerance_energy {:.1e}; \
                     pressure part {:+.3e} {u}, schemes {:+.3e} {u}); the discretisation adds \
                     energy (D90)",
                    whole.kinetic_energy.numerical_dissipation,
                    whole.kinetic_energy.relative,
                    set.tolerance_energy,
                    whole.kinetic_energy.pressure_dissipation,
                    whole.kinetic_energy.scheme_dissipation,
                    u = if s.run.is_some() { "J" } else { "W" },
                ),
            )
            .key("balance.tolerance_energy"),
        );
    }
    out
}

/// One compact line for the progress log (transient).
pub fn compact(b: &Balance) -> String {
    format!(
        "      balance  mass {:+.2e}  kinetic energy: dK/dt {:+.3e} W, boundary {:+.3e} W, \
         dissipation {:.3e} W, numerical {:+.3e} W ({:+.1e})",
        b.mass.relative,
        b.kinetic_energy.storage,
        b.kinetic_energy.boundary + b.kinetic_energy.body_force,
        b.kinetic_energy.viscous_dissipation,
        b.kinetic_energy.numerical_dissipation,
        b.kinetic_energy.relative
    )
}

/// The full table (log, report): per quantity, patch and term.
pub fn table(title: &str, b: &Balance, integrated: bool) -> String {
    let (m, e) = if integrated {
        ("kg", "J")
    } else {
        ("kg/s", "W")
    };
    // the rows are indented as the title is
    let ind = &title[..title.len() - title.trim_start().len()];
    let mut s = format!("{title}\n");
    let row = |s: &mut String, name: &str, v: f64, unit: &str| {
        s.push_str(&format!("{ind}  {name:<44} {v:>+14.6e} {unit}\n"));
    };
    s.push_str(&format!("{ind}  mass\n"));
    for p in &b.patches {
        row(
            &mut s,
            &format!("  {} ({}) out", p.patch, p.kind),
            p.mass_flow,
            m,
        );
    }
    row(&mut s, "  storage", b.mass.storage, m);
    row(
        &mut s,
        &format!("  imbalance ({:+.2e})", b.mass.relative),
        b.mass.imbalance,
        m,
    );
    s.push_str(&format!("{ind}  kinetic energy\n"));
    for p in &b.patches {
        if p.viscous_power != 0.0 {
            row(
                &mut s,
                &format!("  {} viscous work in", p.patch),
                p.viscous_power,
                e,
            );
        }
        if p.convective_power != 0.0 {
            row(
                &mut s,
                &format!("  {} convected out", p.patch),
                p.convective_power,
                e,
            );
        }
        if p.pressure_power != 0.0 {
            row(
                &mut s,
                &format!("  {} pressure work out", p.patch),
                p.pressure_power,
                e,
            );
        }
    }
    let k = &b.kinetic_energy;
    if k.body_force != 0.0 {
        row(&mut s, "  body force", k.body_force, e);
    }
    row(&mut s, "  viscous dissipation", k.viscous_dissipation, e);
    row(&mut s, "  storage", k.storage, e);
    row(
        &mut s,
        &format!("  numerical dissipation ({:+.2e})", k.relative),
        k.numerical_dissipation,
        e,
    );
    row(&mut s, "    of it pressure", k.pressure_dissipation, e);
    row(&mut s, "    of it schemes", k.scheme_dissipation, e);
    s
}

/// The ledger in the report (`flow.md`): one table per solve.
pub fn markdown(s: &BalanceSummary, cells: usize) -> String {
    let (b, unit_m, unit_e, what) = match &s.run {
        Some(r) => (r, "kg", "J", "over the run"),
        None => (&s.last, "kg/s", "W", "of the final state"),
    };
    let mut o = format!(
        "\n{cells} cells, conservation ledger {what}{}:\n\n| Term | Value |\n|---|---|\n",
        match (s.steps, s.negative_steps) {
            (Some(n), Some(k)) => format!(" ({n} steps, {k} with negative numerical dissipation)"),
            _ => String::new(),
        }
    );
    let mut row =
        |name: String, v: f64, u: &str| o.push_str(&format!("| {name} | {v:+.6e} {u} |\n"));
    row("mass imbalance".into(), b.mass.imbalance, unit_m);
    row(
        "largest relative mass imbalance (-)".into(),
        s.largest_mass_relative,
        "",
    );
    for p in &b.patches {
        if p.power() != 0.0 {
            row(
                format!("power into the fluid through {}", p.patch),
                p.power(),
                unit_e,
            );
        }
    }
    let k = &b.kinetic_energy;
    if k.body_force != 0.0 {
        row("body force".into(), k.body_force, unit_e);
    }
    row("viscous dissipation".into(), k.viscous_dissipation, unit_e);
    row(
        "storage (change of kinetic energy)".into(),
        k.storage,
        unit_e,
    );
    row(
        "numerical dissipation".into(),
        k.numerical_dissipation,
        unit_e,
    );
    row("numerical dissipation, relative (-)".into(), k.relative, "");
    row(
        "of it pressure (momentum interpolation)".into(),
        k.pressure_dissipation,
        unit_e,
    );
    row(
        "of it convection and time schemes".into(),
        k.scheme_dissipation,
        unit_e,
    );
    o
}

const CSV_TERMS: &str = "mass_outflow,mass_relative,kinetic_energy,storage,boundary_power,\
body_force_power,viscous_dissipation,numerical_dissipation,numerical_relative,pressure_dissipation,\
scheme_dissipation";

/// The units line (of the D75 preamble) and the header of
/// `balance.csv`: per time step with the time (transient), or per outer
/// iteration (steady).
pub fn csv_head(transient: bool) -> (String, String) {
    let units = "mass_outflow kg/s, mass_relative -, kinetic_energy J, storage boundary_power \
                 body_force_power viscous_dissipation numerical_dissipation W, numerical_relative -, \
                 pressure_dissipation scheme_dissipation W";
    if transient {
        (
            format!("step -, time s, {units}"),
            format!("step,time,{CSV_TERMS}"),
        )
    } else {
        (
            format!("iteration -, {units}"),
            format!("iteration,{CSV_TERMS}"),
        )
    }
}

/// A row of `balance.csv` (with the time when the balance has one).
pub fn csv_row(b: &Balance) -> String {
    let k = &b.kinetic_energy;
    let mut s = format!("{}", b.step.or(b.iteration).unwrap_or(0));
    if let Some(t) = b.time {
        s.push_str(&format!(",{t:e}"));
    }
    for v in [
        b.mass.outflow,
        b.mass.relative,
        k.energy,
        k.storage,
        k.boundary,
        k.body_force,
        k.viscous_dissipation,
        k.numerical_dissipation,
        k.relative,
        k.pressure_dissipation,
        k.scheme_dissipation,
    ] {
        s.push_str(&format!(",{v:e}"));
    }
    s
}

/// The `balance` row of the run record (D69).
pub fn record_row(b: &Balance) -> serde_json::Value {
    let k = &b.kinetic_energy;
    let mut v = serde_json::json!({
        "type": "balance",
        "mass_outflow": b.mass.outflow,
        "mass_relative": b.mass.relative,
        "kinetic_energy": k.energy,
        "storage": k.storage,
        "boundary_power": k.boundary,
        "body_force_power": k.body_force,
        "viscous_dissipation": k.viscous_dissipation,
        "numerical_dissipation": k.numerical_dissipation,
        "numerical_relative": k.relative,
        "pressure_dissipation": k.pressure_dissipation,
        "scheme_dissipation": k.scheme_dissipation,
    });
    if let Some(m) = v.as_object_mut() {
        if let Some(s) = b.step {
            m.insert("step".into(), s.into());
        }
        if let Some(i) = b.iteration {
            m.insert("iteration".into(), i.into());
        }
        if let Some(t) = b.time {
            m.insert("time".into(), t.into());
        }
    }
    v
}

/// The end-of-run chart: bars of the kinetic-energy terms (W, or J over a
/// transient run) with the numerical dissipation.
pub fn chart(name: &str, b: &Balance, integrated: bool) -> String {
    let e = if integrated { "J" } else { "W" };
    let mut bars: Vec<(String, f64)> = Vec::new();
    for p in &b.patches {
        let w = p.power();
        if w != 0.0 {
            bars.push((format!("{} (power in)", p.patch), w));
        }
    }
    let k = &b.kinetic_energy;
    if k.body_force != 0.0 {
        bars.push(("body force".into(), k.body_force));
    }
    bars.push(("viscous dissipation".into(), -k.viscous_dissipation));
    bars.push(("numerical: pressure".into(), -k.pressure_dissipation));
    bars.push(("numerical: schemes".into(), -k.scheme_dissipation));
    bars.push(("storage".into(), -k.storage));
    crate::plot::bars(
        &format!(
            "{name}: kinetic energy{}",
            if integrated { " over the run" } else { "" }
        ),
        &format!("{e} (sources +, sinks −)"),
        &bars,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rates(ke: f64, lid: f64, phi: f64) -> Rates {
        Rates {
            kinetic_energy: ke,
            patches: vec![PatchTerms {
                patch: "lid".into(),
                kind: "wall".into(),
                viscous_power: lid,
                ..Default::default()
            }],
            body_force_power: 0.0,
            pressure_source_power: 0.0,
            viscous_dissipation: phi,
            mass_scale: 1.0,
        }
    }

    #[test]
    fn steady_ledger_closes() {
        let b = steady(&rates(2.0, 1.0, 0.9), 7);
        assert_eq!(b.kinetic_energy.storage, 0.0);
        assert!((b.kinetic_energy.numerical_dissipation - 0.1).abs() < 1e-15);
        assert_eq!(b.mass.imbalance, 0.0);
    }

    #[test]
    fn transient_ledger_integrates_with_the_trapezoidal_rule() {
        // K(t) = t², power in 2t + 1, dissipation 1: exact, no numerical
        // dissipation, for the linear rates
        let mut l = Ledger::start(rates(0.0, 1.0, 1.0), 0.0, 0, 1e-4);
        let mut t: f64 = 0.0;
        for n in 1..=10 {
            t = 0.1 * n as f64;
            let b = l.step(rates(t * t, 2.0 * t + 1.0, 1.0), n, t);
            assert!(b.kinetic_energy.numerical_dissipation.abs() < 1e-12);
        }
        let r = l.run();
        assert!((r.kinetic_energy.storage - t * t).abs() < 1e-12);
        assert!((r.patches[0].viscous_power - (t * t + t)).abs() < 1e-12);
        assert!(r.kinetic_energy.numerical_dissipation.abs() < 1e-12);
        assert_eq!(l.negative_steps, 0);
    }

    #[test]
    fn negative_dissipation_and_mass_are_flagged() {
        let s = BalanceSummary {
            last: steady(&rates(1.0, 1.0, 1.1), 3),
            largest_mass_relative: 1e-6,
            ..Default::default()
        };
        let d = check(&s, &BalanceSection::default());
        let codes: Vec<&str> = d.iter().map(|d| d.code).collect();
        assert_eq!(codes, [diag::BAL_MASS_STEP, diag::BAL_NEGATIVE_DISSIPATION]);
    }
}

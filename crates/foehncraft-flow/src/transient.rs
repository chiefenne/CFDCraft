//! Transient incompressible flow: the PISO algorithm (Issa 1986) with
//! backward differences of first and second order (BDF1, BDF2 with variable
//! time steps) or Crank–Nicolson (D28), on the state of [`Flow`].
//!
//! One time step from tⁿ to tⁿ⁺¹ = tⁿ + Δt:
//!
//! 1. the old levels move on (uⁿ⁻¹ ← uⁿ, ṁⁿ⁻¹ ← ṁⁿ); for Crank–Nicolson the
//!    spatial operator R(uⁿ) of the old level is evaluated with the old flux;
//! 2. per outer iteration (one for plain PISO):
//!    - boundary values (the prescribed ones of tⁿ⁺¹), gradients, limiters;
//!    - the momentum predictor: the rows of [`crate::simplec`] with the time
//!      term ρV a₀/Δt on the diagonal and the old levels on the right-hand
//!      side, convected by the flux of the latest iterate;
//!    - in the first outer iteration of a second-order step the explicit
//!      quantities are extrapolated to tⁿ⁺¹ from the levels n and n−1
//!      (φ* = (1 + ω)φⁿ − ωφⁿ⁻¹): the convecting flux and the velocity
//!      gradients and limiters of the deferred corrections, so that plain
//!      PISO keeps second order (`extrapolation`; a steady state is not
//!      changed by it);
//!    - the predicted Rhie–Chow flux with the time-step term γ_f Y_f of
//!      `foehncraft_kernels::flow::face_mass_flux` (D27: a steady state of the
//!      time stepping is exactly the steady solution, for any Δt);
//!    - `correctors` pressure corrections: Σ c_f (p'_P − p'_N) = −Σ ±ṁ*
//!      with c_f from D = V/(a_P/α) (Issa's PISO), the fluxes corrected to
//!      continuity, u' = δ − D∇p', p += p'; between two correctors the
//!      neighbour correction δ = H(u')/a_P (the part the first corrector
//!      neglects) is added to u and its flux to ṁ;
//! 3. boundary values and gradients of the new level, the old-level
//!    difference Xⁿ⁺¹ = ṁ − ρ ū·S of the Rhie–Chow time-step term and the
//!    largest Courant number.
//!
//! Backward differences: du/dt ≈ (a₀ uⁿ⁺¹ + a₁ uⁿ + a₂ uⁿ⁻¹)/Δt with, for
//! ω = Δtₙ/Δtₙ₋₁, a₀ = (1 + 2ω)/(1 + ω), a₁ = −(1 + ω), a₂ = ω²/(1 + ω)
//! (BDF2; the first step of a run is BDF1: a₀ = 1, a₁ = −1, a₂ = 0).
//! Crank–Nicolson: a₀ = 1, a₁ = −1, θ = ½ on the spatial terms; the pressure
//! is then the one of tⁿ⁺¹ᐟ². Its first step is a BDF1 step as well:
//! Crank–Nicolson does not damp the part of the initial fields that is not
//! consistent with the discrete equations (the face fluxes, the pressure),
//! which otherwise leaves an error of first order in Δt (the start of
//! Rannacher 1984; measured with the Taylor–Green vortex).

use crate::device::{Buf, Buf3, TimeTerms};
use crate::simplec::{Flow, Residuals};
use anyhow::{Result, ensure};
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// Time integration (D28).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "kebab-case")]
pub enum TimeScheme {
    /// First-order backward differences (backward Euler).
    Bdf1,
    /// Second-order backward differences with variable steps (default);
    /// the first step is BDF1.
    #[default]
    Bdf2,
    /// Crank–Nicolson, second order; the first step is BDF1.
    CrankNicolson,
}

impl TimeScheme {
    /// Name for the log.
    pub fn label(self) -> &'static str {
        match self {
            Self::Bdf1 => "BDF1 (backward Euler)",
            Self::Bdf2 => "BDF2",
            Self::CrankNicolson => "Crank–Nicolson",
        }
    }
}

/// The clock of a transient run (part of the restart).
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq)]
pub struct Clock {
    /// Time of the current level (s).
    pub time: f64,
    /// Time steps done.
    pub step: usize,
    /// Size of the last step (s), 0 before the first.
    pub dt_prev: f64,
    /// Steps behind the current level, up to 2: 0 at the start, 1 after
    /// the first step (an old level exists: BDF2, Crank–Nicolson), 2 later
    /// (the old level is itself a solution of the discrete equations: the
    /// extrapolation, which would otherwise carry an initial field that is
    /// not divergence-free into the convecting flux).
    pub levels: usize,
}

/// Report of one time step.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq)]
pub struct StepReport {
    /// Step number (from 1).
    pub step: usize,
    /// Time at the end of the step (s).
    pub time: f64,
    /// Time step (s).
    pub dt: f64,
    /// Residuals of the first outer iteration (momentum before the solve,
    /// continuity of the predicted fluxes); pressure iterations of all
    /// correctors and outer iterations.
    pub residuals: Residuals,
    /// Residuals at the start of the last outer iteration.
    pub last: Residuals,
    /// Outer iterations done.
    pub outer: usize,
    /// Largest Courant number of the new level.
    pub cfl: f64,
}

/// Old levels and work buffers of the time integration.
#[derive(Debug, Clone, Copy)]
struct Levels {
    u_old: Buf3,
    u_old2: Buf3,
    m_old: Buf,
    m_old2: Buf,
    x_cur: Buf,
    x_old: Buf,
    y: Buf,
    conv: Buf,
    r_old: Buf3,
    ucorr: Buf3,
    delta: Buf3,
    dflux: Buf,
    cfl: Buf,
    /// Velocity gradients and limiters of level n−1.
    g_old: [Buf; N_EXPLICIT],
}

/// Number of cell fields in [`explicit_buffers`].
const N_EXPLICIT: usize = 12;

/// The velocity gradients and limiters of `f` (the explicit parts of the
/// deferred corrections).
fn explicit_buffers(f: &crate::simplec::Fields) -> [Buf; N_EXPLICIT] {
    let g = &f.grad_u;
    [
        g[0].x, g[0].y, g[0].z, g[1].x, g[1].y, g[1].z, g[2].x, g[2].y, g[2].z, f.psi.x, f.psi.y,
        f.psi.z,
    ]
}

fn set_explicit_buffers(f: &mut crate::simplec::Fields, b: [Buf; N_EXPLICIT]) {
    let b3 = |i: usize| Buf3 {
        x: b[i],
        y: b[i + 1],
        z: b[i + 2],
    };
    f.grad_u = [b3(0), b3(3), b3(6)];
    f.psi = b3(9);
}

/// The transient driver on a [`Flow`].
pub struct Transient {
    /// Time integration.
    pub scheme: TimeScheme,
    /// The clock.
    pub clock: Clock,
    l: Levels,
}

/// Host copy of the old levels (for a restart).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OldLevels {
    /// uⁿ⁻¹ per component.
    pub u_old: [Vec<f64>; 3],
    /// ṁⁿ⁻¹.
    pub m_old: Vec<f64>,
    /// Xⁿ = ṁⁿ − ρ ūⁿ·S.
    pub x_cur: Vec<f64>,
    /// Xⁿ⁻¹.
    pub x_old: Vec<f64>,
    /// Velocity gradients (9) and limiters (3) of level n−1.
    pub g_old: Vec<Vec<f64>>,
}

impl Transient {
    /// Buffers on the device of `fl`, the clock at `start_time`. Call
    /// [`start`](Self::start) for new start fields, or
    /// [`restore`](Self::restore) after a restart.
    pub fn new(fl: &mut Flow, scheme: TimeScheme, start_time: f64) -> Result<Self> {
        ensure!(
            fl.transient,
            "the flow solver was created for a steady case (no [time] table)"
        );
        let nc = fl.disc.n_cells;
        let nf = fl.disc.n_faces;
        let dev = &mut fl.dev;
        let mut c3 = || -> Result<Buf3> {
            Ok(Buf3 {
                x: dev.alloc(nc)?,
                y: dev.alloc(nc)?,
                z: dev.alloc(nc)?,
            })
        };
        let (u_old, u_old2, r_old, ucorr, delta) = (c3()?, c3()?, c3()?, c3()?, c3()?);
        let l = Levels {
            u_old,
            u_old2,
            r_old,
            ucorr,
            delta,
            m_old: dev.alloc(nf)?,
            m_old2: dev.alloc(nf)?,
            x_cur: dev.alloc(nf)?,
            x_old: dev.alloc(nf)?,
            y: dev.alloc(nf)?,
            conv: dev.alloc(nf)?,
            dflux: dev.alloc(nf)?,
            cfl: dev.alloc(nc)?,
            g_old: {
                let mut g = [Buf(0); N_EXPLICIT];
                for b in &mut g {
                    *b = dev.alloc(nc)?;
                }
                g
            },
        };
        Ok(Self {
            scheme,
            clock: Clock {
                time: start_time,
                ..Default::default()
            },
            l,
        })
    }

    /// Boundary values, gradients and Xⁿ of new start fields.
    pub fn start(&mut self, fl: &mut Flow) -> Result<()> {
        self.finish_level(fl)
    }

    /// Boundary values, gradients and Xⁿ of the current fields.
    fn finish_level(&mut self, fl: &mut Flow) -> Result<()> {
        fl.update_boundaries()?;
        fl.gradients()?;
        let a = fl.flux_args(None, fl.w.zero_face, None);
        fl.dev.mass_flux(&a, self.l.x_cur)?;
        fl.dev
            .scaled_difference(self.l.x_cur, 1.0, fl.f.mass_flux, self.l.x_cur)
    }

    /// The old levels on the host.
    pub fn old_levels(&self, fl: &mut Flow) -> Result<OldLevels> {
        let l = self.l;
        Ok(OldLevels {
            u_old: [
                fl.dev.download(l.u_old.x)?,
                fl.dev.download(l.u_old.y)?,
                fl.dev.download(l.u_old.z)?,
            ],
            m_old: fl.dev.download(l.m_old)?,
            x_cur: fl.dev.download(l.x_cur)?,
            x_old: fl.dev.download(l.x_old)?,
            g_old: l
                .g_old
                .iter()
                .map(|&b| fl.dev.download(b))
                .collect::<Result<_>>()?,
        })
    }

    /// Restores the old levels and the clock (after [`Flow::load`]).
    pub fn restore(&mut self, fl: &mut Flow, o: &OldLevels, clock: Clock) -> Result<()> {
        let l = self.l;
        fl.dev.upload(l.u_old.x, &o.u_old[0])?;
        fl.dev.upload(l.u_old.y, &o.u_old[1])?;
        fl.dev.upload(l.u_old.z, &o.u_old[2])?;
        fl.dev.upload(l.m_old, &o.m_old)?;
        fl.dev.upload(l.x_cur, &o.x_cur)?;
        fl.dev.upload(l.x_old, &o.x_old)?;
        ensure!(
            o.g_old.len() == l.g_old.len(),
            "the restart file holds {} old gradient fields, expected {}",
            o.g_old.len(),
            l.g_old.len()
        );
        for (&b, g) in l.g_old.iter().zip(&o.g_old) {
            fl.dev.upload(b, g)?;
        }
        self.clock = clock;
        Ok(())
    }

    /// Largest Courant number of the current fluxes for a step `dt`.
    pub fn courant(&mut self, fl: &mut Flow, dt: f64) -> Result<f64> {
        fl.dev.courant(fl.f.mass_flux, dt, self.l.cfl)?;
        Ok(fl.dev.min_max(self.l.cfl)?.1)
    }

    /// One time step of size `dt`. The prescribed boundary values of the new
    /// time are set beforehand ([`Flow::set_boundary_velocity`]).
    pub fn step(&mut self, fl: &mut Flow, dt: f64) -> Result<StepReport> {
        ensure!(dt > 0.0 && dt.is_finite(), "time step {dt} must be > 0");
        let s = fl.settings.clone();
        let (l, w, f) = (self.l, fl.w, fl.f);
        let rho = fl.phys.rho;
        let second = self.scheme != TimeScheme::Bdf1 && self.clock.levels >= 1;
        let omega = if second { dt / self.clock.dt_prev } else { 1.0 };
        let (a0, a1, a2, theta) = match self.scheme {
            TimeScheme::CrankNicolson if self.clock.levels >= 1 => (1.0, -1.0, 0.0, 0.5),
            TimeScheme::Bdf2 if second => (
                (1.0 + 2.0 * omega) / (1.0 + omega),
                -(1.0 + omega),
                omega * omega / (1.0 + omega),
                1.0,
            ),
            _ => (1.0, -1.0, 0.0, 1.0),
        };
        let alpha = s.alpha_u(true);
        let alpha_p = s.alpha_p(true);
        let extrapolate = second && s.extrapolation && self.clock.levels >= 2;
        let t = Instant::now();

        // 1. old levels
        if self.clock.levels >= 1 {
            for (a, b) in l.u_old.all().into_iter().zip(l.u_old2.all()) {
                fl.dev.copy(a, b)?;
            }
            fl.dev.copy(l.m_old, l.m_old2)?;
        }
        for (a, b) in f.u.all().into_iter().zip(l.u_old.all()) {
            fl.dev.copy(a, b)?;
        }
        fl.dev.copy(f.mass_flux, l.m_old)?;
        let ew = 1.0 - theta;
        if ew != 0.0 {
            // R(uⁿ) with the old flux, gradients and boundary values
            fl.momentum_assembly_faces(f.mass_flux)?;
            fl.dev.momentum_operator(w.cp, w.cn, w.bf, f.u, l.r_old)?;
        }
        // G* = (1 + ω)Gⁿ − ωGⁿ⁻¹ for the velocity gradients and limiters,
        // into the buffers of Gⁿ⁻¹, which become the current ones; Gⁿ⁻¹ ←
        // Gⁿ. The gradients of level n are those of the end of the last
        // step (with the boundary values of tⁿ), so the extrapolation is
        // consistent also for time-dependent boundary values. The pressure
        // is not extrapolated: p* = 2pⁿ − pⁿ⁻¹ amplifies the pressure modes
        // that the corrections do not remove on skewed meshes (unstable on
        // the cylinder O-grid), and the corrections recover pⁿ⁺¹ anyway.
        let cur = explicit_buffers(&fl.f);
        let old = self.l.g_old;
        if extrapolate {
            for (&a, &b) in cur.iter().zip(&old) {
                fl.dev.lincomb(b, -omega, b, 1.0 + omega, a)?;
            }
            set_explicit_buffers(&mut fl.f, old);
            self.l.g_old = cur;
        } else {
            for (&a, &b) in cur.iter().zip(&old) {
                fl.dev.copy(a, b)?;
            }
        }
        // Y = Σ (−a_k/a₀) X^k
        fl.dev.lincomb(l.y, -a1 / a0, l.x_cur, -a2 / a0, l.x_old)?;
        // the convecting flux of the new level (also for Crank–Nicolson: the
        // implicit half is R(uⁿ⁺¹) with ṁⁿ⁺¹, the explicit one R(uⁿ) with ṁⁿ,
        // as in the later outer iterations)
        let conv = if extrapolate {
            fl.dev
                .lincomb(l.conv, 1.0 + omega, l.m_old, -omega, l.m_old2)?;
            l.conv
        } else {
            f.mass_flux
        };
        let time = TimeTerms {
            tau: rho * a0 / dt,
            c_old: [rho * (-a1) / dt, rho * (-a2) / dt],
            u_old: l.u_old,
            u_old2: l.u_old2,
            explicit: (ew, l.r_old),
        };
        fl.dev.sync()?;
        fl.timings.momentum_assembly += t.elapsed().as_secs_f64();

        let mut first = Residuals::default();
        let mut last = Residuals::default();
        let mut outer = 0;
        let rebuild_step = self.clock.step.is_multiple_of(s.amg_rebuild);
        for k in 1..=s.outer_iterations {
            outer = k;
            fl.iteration += 1;
            let mut res = Residuals {
                iteration: fl.iteration,
                ..Default::default()
            };
            let t = Instant::now();
            fl.update_boundaries()?;
            // first outer iteration of an extrapolated step: the velocity
            // gradients are G*, the pressure gradient is that of pⁿ
            if !(k == 1 && extrapolate) {
                fl.gradients()?;
            }
            fl.dev.sync()?;
            fl.timings.gradients += t.elapsed().as_secs_f64();
            let t = Instant::now();
            let relax = if alpha < 1.0 {
                fl.relaxation_term(alpha)?;
                w.relax
            } else {
                w.zero_face
            };
            let c = if k == 1 { conv } else { fl.f.mass_flux };
            let scale = fl.momentum_assembly(c, alpha, theta, Some(time))?;
            fl.dev.sync()?;
            fl.timings.momentum_assembly += t.elapsed().as_secs_f64();
            let t = Instant::now();
            fl.momentum_solve(scale, &mut res)?;
            fl.dev.sync()?;
            fl.timings.momentum_solve += t.elapsed().as_secs_f64();

            // predicted fluxes with the time-step term
            let t = Instant::now();
            let a = fl.flux_args(Some(w.d), relax, Some((w.ds, l.y, alpha)));
            fl.dev.mass_flux(&a, w.m)?;
            res.continuity = fl.continuity(w.m)?;
            fl.dev.sync()?;
            fl.timings.pressure_assembly += t.elapsed().as_secs_f64();

            for j in 1..=s.correctors {
                let rebuild = rebuild_step && k == 1 && j == 1;
                res.pressure_iterations += fl.pressure_correction(w.m, rebuild)?;
                let t = Instant::now();
                fl.corrections(w.m, alpha_p)?;
                if j < s.correctors {
                    // u' = δ − D∇p' of this corrector (δ = 0 in the first)
                    if j == 1 {
                        for b in l.ucorr.all() {
                            fl.dev.fill(b, 0.0)?;
                        }
                    }
                    fl.dev.correct_velocity(w.dc, w.gpc, l.ucorr)?;
                    // δ = H(u')/a_P, added to u and, as a flux, to ṁ
                    for i in 0..3 {
                        let (x, dg, dl) = (l.ucorr.all()[i], w.diag.all()[i], l.delta.all()[i]);
                        fl.dev.neighbour_correction(w.val, x, dg, dl)?;
                    }
                    for (u, d) in f.u.all().into_iter().zip(l.delta.all()) {
                        fl.dev.lincomb(u, 1.0, u, 1.0, d)?;
                    }
                    fl.dev.delta_flux(l.delta, l.dflux)?;
                    fl.dev.lincomb(w.m, 1.0, w.mc, 1.0, l.dflux)?;
                    for (d, u) in l.delta.all().into_iter().zip(l.ucorr.all()) {
                        fl.dev.copy(d, u)?;
                    }
                } else {
                    std::mem::swap(&mut fl.f.mass_flux, &mut fl.w.mc);
                }
                fl.dev.sync()?;
                fl.timings.pressure_assembly += t.elapsed().as_secs_f64();
            }
            ensure!(
                res.max().is_finite(),
                "the solution diverged in time step {} (t = {:.6e} s, Δt = {dt:.3e} s; residuals \
                 {:?}, continuity {:.3e}); lower [time] time_step or max_cfl",
                self.clock.step + 1,
                self.clock.time + dt,
                res.momentum,
                res.continuity
            );
            if k == 1 {
                first = res;
            } else {
                first.pressure_iterations += res.pressure_iterations;
            }
            last = res;
            if k > 1 && res.max() < s.outer_tolerance {
                break;
            }
        }
        // the mass flux buffer may have moved: f is stale from here
        let t = Instant::now();
        fl.dev.copy(self.l.x_cur, self.l.x_old)?;
        self.finish_level(fl)?;
        let cfl = self.courant(fl, dt)?;
        fl.dev.sync()?;
        fl.timings.gradients += t.elapsed().as_secs_f64();
        self.clock.time += dt;
        self.clock.step += 1;
        self.clock.dt_prev = dt;
        self.clock.levels = (self.clock.levels + 1).min(2);
        Ok(StepReport {
            step: self.clock.step,
            time: self.clock.time,
            dt,
            residuals: first,
            last,
            outer,
            cfl,
        })
    }
}

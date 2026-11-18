//! The coupled pressure–velocity algorithm (L11, D19): the momentum
//! equations and the continuity equation of one outer iteration are solved
//! as one linear system for the increments of the velocity and the
//! pressure, on the discretisation of [`crate::simplec`].
//!
//! One outer iteration, from the fields u, p and the face fluxes ṁ:
//!
//! 1. boundary values, gradients, limiters and the momentum rows A (relaxed
//!    with α_u) as in SIMPLEC;
//! 2. the residuals of the current fields: r_u = b − A u − V ∇p per
//!    component and r_c = −Σ_f ±ṁ*_f, with ṁ* the Rhie–Chow fluxes of u
//!    and p ([`foehncraft_kernels::flow::face_mass_flux`]);
//! 3. the increments δu, δp from
//!
//!    ```text
//!    | A  G | |δu|   |r_u|        G δp = V ∇(δp),
//!    | D  L | |δp| = |r_c|        D δu = Σ_f ± ρ (δu)_f·S,
//!    ```
//!
//!    L δp = Σ_f c_f (δp_P − δp_N) with c_f = ρ D̄_f |S|²/(d·S),
//!    D = V/(a_P/α_u): the linearisation of the same flux formula (the
//!    interpolated pressure gradient of the Rhie–Chow term and the
//!    skewness corrections stay at the old fields; they are in r_c);
//! 4. u += δu, p += δp, ṁ = ṁ* + ρ (δu)_f·S − c_f (δp_N − δp_P): the
//!    fluxes satisfy continuity to the tolerance of the linear solve;
//! 5. turbulence and scalars as in SIMPLEC.
//!
//! The residuals are those of the segregated algorithm, so a converged
//! result is the same (the Rhie–Chow flux does not depend on α_u,
//! Majumdar 1988).
//!
//! **Linear solve.** The operator is applied matrix-free with the
//! operations of the device (the momentum matrix, the least-squares or
//! Gauss gradient with the mask of a pressure correction, the flux of a
//! velocity increment, the pressure matrix), so the system needs no block
//! matrix. Flexible GMRES (Saad 1993) with right preconditioning; the
//! preconditioner is one SIMPLE-type sweep on the increments (Patankar &
//! Spalding 1972; as a preconditioner: Klaij & Vuik 2013):
//!
//! ```text
//! z_u = Ã⁻¹ r_u,   Ŝ z_p = r_p − D z_u,   z_u −= D̂ ∇z_p,
//! ```
//!
//! with Ã⁻¹ a few BiCGStab iterations per component, D̂ = D^C (SIMPLEC) or
//! D (SIMPLE), and Ŝ the pressure matrix with the coefficient
//! ρ (D̄ + D̂̄)_f |S|²/(d·S) (the Rhie–Chow part L plus the compact form of
//! −D D̂ G), solved by CG with the multigrid of the pressure.

use crate::device::{Buf, Buf3, Lsq, Mask};
use crate::simplec::{Flow, Residuals};
use anyhow::Result;
use foehncraft_fv::solver::Preconditioner;
use foehncraft_kernels::flow as k;
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// The approximation D̂ of the inverse momentum diagonal in the Schur
/// complement of the preconditioner.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Schur {
    /// D^C = V / (a_P/α − Σ|a_N|) (Van Doormaal & Raithby 1984).
    #[default]
    Simplec,
    /// D = V / (a_P/α).
    Simple,
}

/// `[solver.coupled]`: the linear solve of the coupled algorithm.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct CoupledSettings {
    /// Relative reduction of the residual of the coupled system per outer
    /// iteration.
    pub tolerance: f64,
    /// Largest number of Krylov iterations per outer iteration.
    pub max_iter: usize,
    /// The momentum inverse in the Schur complement of the preconditioner.
    pub schur: Schur,
    /// Relative tolerance of the momentum sweeps of the preconditioner.
    pub momentum_tolerance: f64,
    /// Largest number of BiCGStab iterations of a momentum sweep.
    pub momentum_max_iter: usize,
    /// Preconditioner of the momentum sweeps (multigrid: slot 5, built
    /// from the symmetric part of the momentum matrix).
    pub momentum_preconditioner: Preconditioner,
    /// Relative tolerance of the pressure solve of the preconditioner.
    pub pressure_tolerance: f64,
    /// Largest number of CG iterations of that pressure solve.
    pub pressure_max_iter: usize,
    /// After the Krylov solve one pressure correction makes the face
    /// fluxes satisfy continuity to the tolerance of `[solver.pressure]`.
    pub projection: bool,
    /// Factor on the increments of the velocity (explicit relaxation).
    pub velocity_factor: f64,
    /// Factor on the increment of the pressure (explicit relaxation).
    pub pressure_factor: f64,
}

impl Default for CoupledSettings {
    fn default() -> Self {
        Self {
            tolerance: 0.1,
            max_iter: 30,
            schur: Schur::Simplec,
            momentum_tolerance: 0.1,
            momentum_max_iter: 5,
            momentum_preconditioner: Preconditioner::Jacobi,
            pressure_tolerance: 0.1,
            pressure_max_iter: 20,
            projection: true,
            velocity_factor: 1.0,
            pressure_factor: 1.0,
        }
    }
}

impl CoupledSettings {
    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        for (name, t) in [
            ("tolerance", self.tolerance),
            ("momentum_tolerance", self.momentum_tolerance),
            ("pressure_tolerance", self.pressure_tolerance),
        ] {
            crate::need!(
                t > 0.0 && t < 1.0,
                CASE_RANGE,
                "[solver.coupled] {name} must lie in (0, 1)"
            );
        }
        for (name, n) in [
            ("max_iter", self.max_iter),
            ("momentum_max_iter", self.momentum_max_iter),
            ("pressure_max_iter", self.pressure_max_iter),
        ] {
            crate::need!(n >= 1, CASE_RANGE, "[solver.coupled] {name} must be >= 1");
        }
        crate::need!(
            self.max_iter <= 200,
            CASE_RANGE,
            "[solver.coupled] max_iter must be <= 200 (the Krylov vectors are stored)"
        );
        for (name, f) in [
            ("velocity_factor", self.velocity_factor),
            ("pressure_factor", self.pressure_factor),
        ] {
            crate::need!(
                f > 0.0 && f <= 1.0,
                CASE_RANGE,
                "[solver.coupled] {name} must lie in (0, 1]"
            );
        }
        Ok(())
    }
}

/// A vector of the coupled system: three velocity components and the
/// pressure, one value per cell each.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Vec4 {
    pub u: Buf3,
    pub p: Buf,
}

impl Vec4 {
    fn all(self) -> [Buf; 4] {
        [self.u.x, self.u.y, self.u.z, self.p]
    }
}

/// Multigrid slot of the momentum sweeps (0 pressure, 1 to 4 turbulence).
const MOMENTUM_SLOT: usize = 5;

/// Buffers of the coupled solve.
#[derive(Debug)]
pub(crate) struct CoupledWork {
    /// −V per cell.
    neg_vol: Buf,
    /// 1 per cell.
    one: Buf,
    /// ρ Σ|S| over the interior faces plus half of that of the boundary
    /// faces (the scale of the continuity row).
    area: f64,
    /// The Rhie–Chow pressure matrix L and its face coefficients.
    pval: Buf,
    cf: Buf,
    /// D + D̂.
    dsum: Buf,
    /// D̂.
    dhat: Buf,
    /// Right-hand side, solution, work vector.
    r: Vec4,
    x: Vec4,
    w: Vec4,
    /// Gradient, cell and face work.
    g: Buf3,
    t: Buf,
    mf: Buf,
    /// ṁ − ρ ū·S of the current fields (pseudo-time term of the flux).
    y: Buf,
    /// The multigrid of the momentum sweeps exists.
    momentum_amg: bool,
    /// Krylov basis and preconditioned basis.
    v: Vec<Vec4>,
    z: Vec<Vec4>,
}

/// Counters of one coupled solve.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Counts {
    pub krylov: usize,
    pub momentum: [usize; 3],
    pub pressure: usize,
    /// Reached relative residual.
    pub reduction: f64,
}

impl Flow {
    fn vec4(&mut self) -> Result<Vec4> {
        let nc = self.disc.n_cells;
        Ok(Vec4 {
            u: Buf3 {
                x: self.dev.alloc(nc)?,
                y: self.dev.alloc(nc)?,
                z: self.dev.alloc(nc)?,
            },
            p: self.dev.alloc(nc)?,
        })
    }

    /// The buffers of the coupled solve, allocated at the first coupled
    /// iteration (and again when `max_iter` grows).
    fn coupled_work(&mut self) -> Result<()> {
        let m = self.settings.coupled.max_iter;
        if let Some(cw) = self.coupled.as_ref()
            && cw.z.len() >= m
        {
            return Ok(());
        }
        let (nc, nf) = (self.disc.n_cells, self.disc.n_faces);
        let stored = self.disc.sell.stored();
        let mut cw = match self.coupled.take() {
            Some(cw) => cw,
            None => {
                let neg_vol = self.dev.alloc(nc)?;
                let v: Vec<f64> = self.disc.cell_volume.iter().map(|v| -v).collect();
                self.dev.upload(neg_vol, &v)?;
                let one = self.dev.alloc(nc)?;
                self.dev.fill(one, 1.0)?;
                let d = &self.disc;
                let mut area = 0.0;
                for f in 0..nf {
                    let a = (d.face_area.x[f].powi(2)
                        + d.face_area.y[f].powi(2)
                        + d.face_area.z[f].powi(2))
                    .sqrt();
                    area += if f < d.n_interior { a } else { 0.5 * a };
                }
                CoupledWork {
                    neg_vol,
                    one,
                    area: self.phys.rho * area,
                    pval: self.dev.alloc(stored)?,
                    cf: self.dev.alloc(nf)?,
                    dsum: self.dev.alloc(nc)?,
                    dhat: self.dev.alloc(nc)?,
                    r: self.vec4()?,
                    x: self.vec4()?,
                    w: self.vec4()?,
                    g: Buf3 {
                        x: self.dev.alloc(nc)?,
                        y: self.dev.alloc(nc)?,
                        z: self.dev.alloc(nc)?,
                    },
                    t: self.dev.alloc(nc)?,
                    mf: self.dev.alloc(nf)?,
                    y: self.dev.alloc(nf)?,
                    momentum_amg: false,
                    v: Vec::new(),
                    z: Vec::new(),
                }
            }
        };
        while cw.z.len() < m {
            let z = self.vec4()?;
            cw.z.push(z);
        }
        while cw.v.len() < m + 1 {
            let v = self.vec4()?;
            cw.v.push(v);
        }
        self.coupled = Some(cw);
        Ok(())
    }

    fn dot4(&mut self, a: Vec4, b: Vec4) -> Result<f64> {
        let mut s = 0.0;
        for (x, y) in a.all().into_iter().zip(b.all()) {
            s += self.dev.dot(x, y)?;
        }
        Ok(s)
    }

    fn norm4(&mut self, a: Vec4) -> Result<f64> {
        let mut s = 0.0;
        for x in a.all() {
            s += self.dev.sum_sq(x)?;
        }
        Ok(s.sqrt())
    }

    /// out = a x + b y (out may be x or y).
    fn lincomb4(&mut self, out: Vec4, a: f64, x: Vec4, b: f64, y: Vec4) -> Result<()> {
        for ((o, x), y) in out.all().into_iter().zip(x.all()).zip(y.all()) {
            self.dev.lincomb(o, a, x, b, y)?;
        }
        Ok(())
    }

    /// The momentum matrix of component `i` in `w.val` (the diagonal of the
    /// symmetry planes differs per component).
    fn momentum_matrix(&mut self, i: usize, sym: bool) -> Result<()> {
        if sym {
            self.dev.set_diagonal(self.w.val, self.w.diag.all()[i])?;
        }
        Ok(())
    }

    /// out = K x: the coupled operator, the continuity row scaled with
    /// `sc`.
    fn coupled_operator(&mut self, x: Vec4, out: Vec4, sc: f64, sym: bool) -> Result<()> {
        let (g, t, mf, neg_vol, pval) = {
            let cw = self.coupled.as_ref().expect("coupled work");
            (cw.g, cw.t, cw.mf, cw.neg_vol, cw.pval)
        };
        let w = self.w;
        self.dev.gradient(
            Lsq::Pressure,
            x.p,
            Mask::PressureCorrection,
            w.zero_boundary,
            g,
        )?;
        for i in 0..3 {
            self.momentum_matrix(i, sym)?;
            self.dev
                .matrix_product(w.val, x.u.all()[i], out.u.all()[i])?;
        }
        // + V ∇x_p
        self.dev.correct_velocity(neg_vol, g, out.u)?;
        self.dev.delta_flux(x.u, mf)?;
        self.dev.net_outflow(mf, t)?;
        self.dev.matrix_product(pval, x.p, out.p)?;
        self.dev.lincomb(out.p, sc, out.p, sc, t)
    }

    /// z = M⁻¹ r: one SIMPLE-type sweep on the increments.
    fn coupled_preconditioner(
        &mut self,
        r: Vec4,
        z: Vec4,
        sc: f64,
        sym: bool,
        dhat: Buf,
        n: &mut Counts,
    ) -> Result<()> {
        let (g, t, mf) = {
            let cw = self.coupled.as_ref().expect("coupled work");
            (cw.g, cw.t, cw.mf)
        };
        let w = self.w;
        let s = self.settings.coupled.clone();
        let tm = Instant::now();
        let mpc = s.momentum_preconditioner;
        if mpc == Preconditioner::Amg {
            self.dev.amg_select(MOMENTUM_SLOT)?;
        }
        for i in 0..3 {
            self.momentum_matrix(i, sym)?;
            self.dev.fill(z.u.all()[i], 0.0)?;
            let (it, _) = self.dev.solve(
                w.val,
                r.u.all()[i],
                z.u.all()[i],
                mpc,
                false,
                s.momentum_tolerance,
                s.momentum_max_iter,
            )?;
            n.momentum[i] += it;
        }
        if mpc == Preconditioner::Amg {
            self.dev.amg_select(0)?;
        }
        self.dev.sync()?;
        self.timings.momentum_solve += tm.elapsed().as_secs_f64();
        let tm = Instant::now();
        self.dev.delta_flux(z.u, mf)?;
        self.dev.net_outflow(mf, t)?;
        self.dev.lincomb(t, 1.0 / sc, r.p, -1.0, t)?;
        self.dev.fill(z.p, 0.0)?;
        let (it, _) = self.dev.solve(
            w.pval,
            t,
            z.p,
            self.settings.pressure.preconditioner,
            true,
            s.pressure_tolerance,
            s.pressure_max_iter,
        )?;
        n.pressure += it;
        if let Some(rc) = self.reference_cell {
            let shift = self.dev.value(z.p, rc)?;
            self.dev.subtract(z.p, shift)?;
        }
        self.dev.gradient(
            Lsq::Pressure,
            z.p,
            Mask::PressureCorrection,
            w.zero_boundary,
            g,
        )?;
        self.dev.correct_velocity(dhat, g, z.u)?;
        self.dev.sync()?;
        self.timings.pressure_solve += tm.elapsed().as_secs_f64();
        Ok(())
    }

    /// Flexible GMRES on K x = r from x = 0 (right preconditioning,
    /// modified Gram–Schmidt, Givens rotations); the solution in `cw.x`.
    fn coupled_solve(&mut self, sc: f64, sym: bool, dhat: Buf) -> Result<Counts> {
        let (r, x, wv, v, z) = {
            let cw = self.coupled.as_ref().expect("coupled work");
            (cw.r, cw.x, cw.w, cw.v.clone(), cw.z.clone())
        };
        let s = self.settings.coupled.clone();
        let m = s.max_iter;
        let mut n = Counts::default();
        for b in x.all() {
            self.dev.fill(b, 0.0)?;
        }
        let beta = self.norm4(r)?;
        if beta == 0.0 || !beta.is_finite() {
            n.reduction = if beta == 0.0 { 0.0 } else { f64::NAN };
            return Ok(n);
        }
        self.lincomb4(v[0], 1.0 / beta, r, 0.0, r)?;
        // Hessenberg matrix by columns, rotations, right-hand side
        let mut h = vec![vec![0.0; m + 1]; m];
        let (mut cs, mut sn) = (vec![0.0; m], vec![0.0; m]);
        let mut gk = vec![0.0; m + 1];
        gk[0] = beta;
        let mut kk = 0;
        for j in 0..m {
            self.coupled_preconditioner(v[j], z[j], sc, sym, dhat, &mut n)?;
            self.coupled_operator(z[j], wv, sc, sym)?;
            for i in 0..=j {
                h[j][i] = self.dot4(wv, v[i])?;
                self.lincomb4(wv, 1.0, wv, -h[j][i], v[i])?;
            }
            let hn = self.norm4(wv)?;
            h[j][j + 1] = hn;
            for i in 0..j {
                let (a, b) = (h[j][i], h[j][i + 1]);
                h[j][i] = cs[i] * a + sn[i] * b;
                h[j][i + 1] = -sn[i] * a + cs[i] * b;
            }
            let (a, b) = (h[j][j], h[j][j + 1]);
            let d = a.hypot(b);
            if d == 0.0 || !d.is_finite() {
                break;
            }
            cs[j] = a / d;
            sn[j] = b / d;
            h[j][j] = d;
            h[j][j + 1] = 0.0;
            gk[j + 1] = -sn[j] * gk[j];
            gk[j] *= cs[j];
            kk = j + 1;
            n.krylov = kk;
            n.reduction = gk[j + 1].abs() / beta;
            if n.reduction <= s.tolerance || hn == 0.0 {
                break;
            }
            self.lincomb4(v[j + 1], 1.0 / hn, wv, 0.0, wv)?;
        }
        // back substitution, x = Σ y_j z_j
        let mut y = vec![0.0; kk];
        for i in (0..kk).rev() {
            let mut a = gk[i];
            for j in i + 1..kk {
                a -= h[j][i] * y[j];
            }
            y[i] = a / h[i][i];
        }
        for (j, yj) in y.iter().enumerate() {
            self.lincomb4(x, 1.0, x, *yj, z[j])?;
        }
        Ok(n)
    }

    /// Development aid: writes the coupled matrix of this iteration (by
    /// applying the operator to every unit vector), the right-hand side and
    /// the cell data of the preconditioner as text files into `dir`.
    fn coupled_dump(&mut self, dir: &std::path::Path, sym: bool, dhat: Buf) -> Result<()> {
        use std::io::Write;
        std::fs::create_dir_all(dir)?;
        let nc = self.disc.n_cells;
        let (r, x, wv) = {
            let cw = self.coupled.as_ref().expect("coupled work");
            (cw.r, cw.x, cw.w)
        };
        let mut out = std::io::BufWriter::new(std::fs::File::create(dir.join("K.txt"))?);
        let mut e = vec![0.0; nc];
        for b in x.all() {
            self.dev.fill(b, 0.0)?;
        }
        for (k, b) in x.all().into_iter().enumerate() {
            for c in 0..nc {
                e[c] = 1.0;
                self.dev.upload(b, &e)?;
                e[c] = 0.0;
                self.coupled_operator(x, wv, 1.0, sym)?;
                for (kk, o) in wv.all().into_iter().enumerate() {
                    let v = self.dev.download(o)?;
                    for (i, a) in v.iter().enumerate() {
                        if *a != 0.0 {
                            writeln!(out, "{} {} {:e}", kk * nc + i, k * nc + c, a)?;
                        }
                    }
                }
            }
            self.dev.fill(b, 0.0)?;
        }
        out.flush()?;
        let w = self.w;
        let put = |name: &str, bufs: &[Buf], dev: &mut Box<dyn crate::device::FlowDevice>| -> Result<()> {
            let mut f = std::io::BufWriter::new(std::fs::File::create(dir.join(name))?);
            let cols: Vec<Vec<f64>> = bufs
                .iter()
                .map(|b| dev.download(*b))
                .collect::<Result<_>>()?;
            for i in 0..cols[0].len() {
                let line: Vec<String> = cols.iter().map(|c| format!("{:e}", c[i])).collect();
                writeln!(f, "{}", line.join(" "))?;
            }
            Ok(())
        };
        put("rhs.txt", &r.all(), &mut self.dev)?;
        put("cells.txt", &[w.d, w.dc, dhat, w.diag.x, w.diag.y, w.diag.z], &mut self.dev)?;
        put("state.txt", &[self.f.u.x, self.f.u.y, self.f.u.z, self.f.p], &mut self.dev)?;
        let d = &self.disc;
        let mut f = std::io::BufWriter::new(std::fs::File::create(dir.join("geom.txt"))?);
        for c in 0..nc {
            writeln!(
                f,
                "{:e} {:e} {:e} {:e}",
                d.cell_centroid.x[c], d.cell_centroid.y[c], d.cell_centroid.z[c], d.cell_volume[c]
            )?;
        }
        f.flush()?;
        // the matrix of the preconditioner's pressure solve
        let s = d.matrix_csr(&self.dev.download(w.pval)?);
        let mut f = std::io::BufWriter::new(std::fs::File::create(dir.join("S.txt"))?);
        for i in 0..s.n_rows {
            for (j, a) in s.row(i) {
                writeln!(f, "{i} {j} {a:e}")?;
            }
        }
        f.flush()?;
        std::fs::write(
            dir.join("meta.txt"),
            format!(
                "cells {nc}\nalpha {}\nrho {}\niteration {}\n",
                self.settings.alpha_u(false),
                self.phys.rho,
                self.iteration
            ),
        )?;
        eprintln!("coupled: matrix of iteration {} written to {}", self.iteration, dir.display());
        Ok(())
    }

    /// One steady outer iteration of the coupled algorithm; returns the
    /// residuals of the fields before the solve.
    pub(crate) fn iterate_coupled(&mut self) -> Result<Residuals> {
        self.iteration += 1;
        let mut res = Residuals {
            iteration: self.iteration,
            ..Default::default()
        };
        let s = self.settings.clone();
        let alpha = s.alpha_u(false);
        self.coupled_work()?;
        let (r, x, t, mf, cf, pval, dsum, one, area, dhat, y) = {
            let cw = self.coupled.as_ref().expect("coupled work");
            (
                cw.r, cw.x, cw.t, cw.mf, cw.cf, cw.pval, cw.dsum, cw.one, cw.area, cw.dhat, cw.y,
            )
        };
        // prototype: pseudo-time step Δτ (the same in every cell) from an
        // iteration on
        let env = |k: &str| std::env::var(k).ok().and_then(|v| v.parse::<f64>().ok());
        let pseudo = env("FOEHN_PSEUDO_DT")
            .filter(|_| self.iteration as f64 >= env("FOEHN_PSEUDO_FROM").unwrap_or(1.0));
        let tau = pseudo.map_or(0.0, |dt| self.phys.rho / dt);
        let tm = Instant::now();
        self.update_boundaries()?;
        self.gradients()?;
        self.dev.sync()?;
        self.timings.gradients += tm.elapsed().as_secs_f64();

        // momentum rows and the residual of the current fields
        let tm = Instant::now();
        self.relaxation_term(alpha)?;
        let time = pseudo.map(|_| crate::device::TimeTerms {
            tau,
            c_old: [tau, 0.0],
            u_old: self.f.u,
            u_old2: self.f.u,
            explicit: (0.0, self.f.u),
        });
        if pseudo.is_some() {
            let relax = self.w.relax;
            self.dev
                .lincomb(y, 1.0 / (1.0 - alpha), relax, 0.0, relax)?;
        }
        let scale = self.momentum_assembly(self.f.mass_flux, alpha, 1.0, time)?;
        let (f, w) = (self.f, self.w);
        let sym = self.bc.mom_kind.contains(&k::MOM_SYMMETRY);
        for i in 0..3 {
            self.momentum_matrix(i, sym)?;
            self.dev.matrix_product(w.val, f.u.all()[i], t)?;
            self.dev
                .lincomb(r.u.all()[i], 1.0, w.rhs.all()[i], -1.0, t)?;
            res.momentum[i] = self.dev.sum_sq(r.u.all()[i])?.sqrt() / scale;
        }
        self.dev.sync()?;
        self.timings.momentum_assembly += tm.elapsed().as_secs_f64();

        // Rhie–Chow fluxes of the current fields, the continuity residual,
        // the matrix L and the matrix Ŝ of the preconditioner
        let tm = Instant::now();
        let a = self.flux_args(Some(w.d), w.relax, pseudo.map(|_| (w.ds, y, alpha)));
        self.dev.mass_flux(&a, w.m)?;
        self.close_backflow(w.m)?;
        res.continuity = self.continuity(w.m)?;
        self.dev.pressure_faces(w.d, cf)?;
        self.dev.pressure_rows(cf, w.m, pval, r.p)?;
        // D̂: D^C bounded by its value for fluxes that satisfy continuity,
        // D/(1 − α) (prototype: on the host)
        {
            let d = self.dev.download(w.d)?;
            let dc = self.dev.download(w.dc)?;
            let ds = self.dev.download(w.ds)?;
            // V / ((1 − α) a_P/α + τ V/α)
            let v: Vec<f64> = match s.coupled.schur {
                Schur::Simplec => dc
                    .iter()
                    .zip(&ds)
                    .map(|(c, s)| c.min(1.0 / ((1.0 - alpha) / s + tau / alpha)))
                    .collect(),
                Schur::Simple => d,
            };
            self.dev.upload(dhat, &v)?;
        }
        self.dev.lincomb(dsum, 1.0, w.d, 1.0, dhat)?;
        self.dev.pressure_faces(dsum, w.cf)?;
        self.dev.pressure_rows(w.cf, w.m, w.pval, w.prhs)?;
        if let Some(rc) = self.reference_cell {
            let pos = self.disc.sell.pos(rc, 0);
            self.dev.scale_entry(pval, pos, 2.0)?;
            self.dev.scale_entry(w.pval, pos, 2.0)?;
        }
        self.dev.sync()?;
        self.timings.pressure_assembly += tm.elapsed().as_secs_f64();
        let rebuild = (self.iteration - 1).is_multiple_of(s.amg_rebuild);
        self.pressure_amg(rebuild)?;
        if s.coupled.momentum_preconditioner == Preconditioner::Amg
            && (rebuild || !self.coupled.as_ref().is_some_and(|c| c.momentum_amg))
        {
            let tm = Instant::now();
            self.momentum_matrix(0, sym)?;
            self.dev
                .amg_setup_slot(MOMENTUM_SLOT, w.val, &s.pressure.amg, true)?;
            self.dev.amg_select(0)?;
            if let Some(c) = self.coupled.as_mut() {
                c.momentum_amg = true;
            }
            self.timings.amg_setup += tm.elapsed().as_secs_f64();
        }

        if let Some(dir) = std::env::var_os("FOEHN_COUPLED_DUMP") {
            let at: usize = std::env::var("FOEHN_COUPLED_DUMP_AT")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(20);
            if self.iteration == at {
                self.coupled_dump(std::path::Path::new(&dir), sym, dhat)?;
            }
        }
        // the continuity row in the units of the momentum rows: times the
        // velocity Σ a_P / (ρ Σ|S|), the mean diagonal coefficient over the
        // face areas (the mean normal speed where convection dominates,
        // ν/h in a fluid at rest)
        let sc = self.dev.dot(w.diag.x, one)? / area;
        self.dev.lincomb(r.p, sc, r.p, 0.0, r.p)?;
        let n = self.coupled_solve(sc, sym, dhat)?;
        res.momentum_iterations = n.momentum;
        res.pressure_iterations = n.pressure;
        res.coupled_iterations = n.krylov;
        if std::env::var_os("FOEHN_COUPLED_DEBUG").is_some() {
            let xn = self.norm4(x)?;
            let un = self.dev.sum_sq(f.u.x)?.sqrt();
            eprintln!(
                "coupled {}: krylov {} reduction {:.3e} |x| {:.3e} |u_x| {:.3e} sc {:.3e} scale {:.3e} mom {:?} p {}",
                self.iteration, n.krylov, n.reduction, xn, un, sc, scale, n.momentum, n.pressure
            );
        }

        // the increments
        let tm = Instant::now();
        let (fu, fp) = (s.coupled.velocity_factor, s.coupled.pressure_factor);
        if fu != 1.0 || fp != 1.0 {
            for b in x.u.all() {
                self.dev.lincomb(b, fu, b, 0.0, b)?;
            }
            self.dev.lincomb(x.p, fp, x.p, 0.0, x.p)?;
        }
        if let Some(rc) = self.reference_cell {
            let shift = self.dev.value(x.p, rc)?;
            self.dev.subtract(x.p, shift)?;
        }
        for (u, d) in f.u.all().into_iter().zip(x.u.all()) {
            self.dev.lincomb(u, 1.0, u, 1.0, d)?;
        }
        self.dev.correct_fluxes(cf, w.m, x.p, w.mc)?;
        self.dev.delta_flux(x.u, mf)?;
        self.dev.lincomb(w.mc, 1.0, w.mc, 1.0, mf)?;
        self.dev.lincomb(f.p, 1.0, f.p, 1.0, x.p)?;
        if s.coupled.projection {
            // one pressure correction on the remaining continuity defect
            self.dev.pressure_rows(w.cf, w.mc, w.pval, w.prhs)?;
            if let Some(rc) = self.reference_cell {
                let pos = self.disc.sell.pos(rc, 0);
                self.dev.scale_entry(w.pval, pos, 2.0)?;
            }
            self.dev.fill(w.pc, 0.0)?;
            let (it, _) = self.dev.solve(
                w.pval,
                w.prhs,
                w.pc,
                s.pressure.preconditioner,
                true,
                s.pressure_tol(false),
                s.pressure.max_iter,
            )?;
            res.pressure_iterations += it;
            if let Some(rc) = self.reference_cell {
                let shift = self.dev.value(w.pc, rc)?;
                self.dev.subtract(w.pc, shift)?;
            }
            self.dev.correct_fluxes(w.cf, w.mc, w.pc, w.m)?;
            self.dev.copy(w.m, w.mc)?;
            self.dev.gradient(
                Lsq::Pressure,
                w.pc,
                Mask::PressureCorrection,
                w.zero_boundary,
                w.gpc,
            )?;
            self.dev.correct_velocity(dhat, w.gpc, f.u)?;
            self.dev.lincomb(f.p, 1.0, f.p, 1.0, w.pc)?;
        }
        self.drive_mass_flow()?;
        self.dev.sync()?;
        self.timings.pressure_assembly += tm.elapsed().as_secs_f64();
        if !res.finite() || !n.reduction.is_finite() {
            self.diverged = Some(res.clone());
            return Err(crate::diag::err(crate::diag::Diagnostic::new(
                crate::diag::RUN_DIVERGED,
                format!(
                    "the solution diverged in outer iteration {} (residuals {:?}, continuity {:.3e})",
                    self.iteration, res.momentum, res.continuity
                ),
            )));
        }
        std::mem::swap(&mut self.f.mass_flux, &mut self.w.mc);
        if let Some(t) = self.turbulence.as_mut() {
            t.amg_rebuild = rebuild;
            t.pseudo_dt = pseudo.and(env("FOEHN_PSEUDO_DT_TURB"));
            t.pseudo_ct = env("FOEHN_PSEUDO_CT");
        }
        self.solve_scalars(None, &mut res)?;
        Ok(res)
    }
}

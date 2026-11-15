//! CPU device: the kernel bodies of `cfdcraft_kernels::flow` in
//! multi-threaded loops (rayon), the linear solves on the CPU backend of
//! `cfdcraft-fv`. Each output is written by one index, so the result is the
//! same for any number of threads (D31).

use super::{
    BoundaryArgs, Buf, Buf3, Eps2, FlowDevice, FluxArgs, Lsq, Mask, MomentumFaceArgs,
    MomentumRowArgs, Setup,
};
use crate::bc::BoundaryArrays;
use anyhow::{Result, ensure};
use cfdcraft_fv::amg::AmgSettings;
use cfdcraft_fv::backend::{Backend, V};
use cfdcraft_fv::disc::{Discretisation, Soa3};
use cfdcraft_fv::solver::{self as linear, Preconditioner};
use cfdcraft_kernels::flow::{self as k, FaceGeometry, VelocityGradient};
use cfdcraft_kernels::fv::{self, SellLayout, Vec3};
use rayon::prelude::*;

/// Indices per task of the loops.
const GRAIN: usize = 1 << 12;

/// Meshes with fewer cells run on one thread.
pub const SERIAL_BELOW: usize = 20_000;

#[derive(Clone, Copy)]
struct Ptr(*mut f64);
// SAFETY: every loop below writes disjoint entries, one index each.
unsafe impl Send for Ptr {}
// SAFETY: as above.
unsafe impl Sync for Ptr {}
impl Ptr {
    fn of(v: &mut [f64]) -> Self {
        Self(v.as_mut_ptr())
    }
    // methods (not field access) so that closures capture the whole `Ptr`
    fn get(&self) -> *mut f64 {
        self.0
    }
    unsafe fn at(&self, i: usize) -> *mut f64 {
        // SAFETY: guaranteed by the caller (i inside the slice).
        unsafe { self.0.add(i) }
    }
}

/// Runs `f(i)` for `i < n` on the pool, in chunks.
fn par_for(pool: &rayon::ThreadPool, n: usize, f: impl Fn(usize) + Sync) {
    pool.install(|| {
        (0..n.div_ceil(GRAIN)).into_par_iter().for_each(|ch| {
            for i in ch * GRAIN..((ch + 1) * GRAIN).min(n) {
                f(i);
            }
        })
    });
}

/// `out[i] = f(i, out[i])`.
fn par_vec(pool: &rayon::ThreadPool, out: &mut [f64], f: impl Fn(usize, f64) -> f64 + Sync) {
    pool.install(|| {
        out.par_chunks_mut(GRAIN).enumerate().for_each(|(c, ch)| {
            for (j, o) in ch.iter_mut().enumerate() {
                *o = f(c * GRAIN + j, *o);
            }
        })
    });
}

/// Reproducible sum of squares (fixed blocks, fixed combination).
fn sum_sq(v: &[f64]) -> f64 {
    let blocks: Vec<f64> = v
        .chunks(cfdcraft_core::repro::BLOCK)
        .map(|b| b.iter().map(|x| x * x).sum())
        .collect();
    cfdcraft_core::repro::combine_block_sums(&blocks)
}

/// The CPU device.
pub struct CpuFlow {
    disc: Discretisation,
    lsq_p: Soa3,
    bc: BoundaryArrays,
    rho: f64,
    mu: f64,
    eps2: Vec<f64>,
    mu_face: Vec<f64>,
    /// Masks of the gradients: velocity, pressure, pressure correction.
    masks: [Vec<f64>; 3],
    scratch: Vec<f64>,
    zero_cells: Vec<f64>,
    bufs: Vec<Vec<f64>>,
    pool: rayon::ThreadPool,
    be: Box<dyn Backend>,
}

impl CpuFlow {
    /// A device for the static data of `s`.
    pub fn new(s: &Setup) -> Result<Self> {
        let nc = s.disc.n_cells;
        let nb = s.disc.n_faces - s.disc.n_interior;
        // small meshes run on one thread: the latency of waking the pool
        // would dominate (the result is the same for any number of threads)
        let threads = if nc < SERIAL_BELOW {
            Some(1)
        } else {
            s.threads
        };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.unwrap_or(0))
            .build()?;
        let be = cfdcraft_fv::backend::create("cpu", s.disc, threads, std::path::Path::new(""))?;
        let pca: Vec<f64> = (0..nb)
            .map(|fb| {
                if s.bc.flux_kind[fb] == k::FLUX_PRESSURE {
                    0.0
                } else {
                    1.0
                }
            })
            .collect();
        Ok(Self {
            masks: [
                s.bc.u_zero_gradient.clone(),
                s.bc.p_zero_gradient.clone(),
                pca,
            ],
            disc: s.disc.clone(),
            lsq_p: s.lsq_p.clone(),
            bc: s.bc.clone(),
            rho: s.rho,
            mu: s.mu,
            eps2: s.eps2.to_vec(),
            mu_face: Vec::new(),
            scratch: vec![0.0; nb],
            zero_cells: vec![0.0; nc],
            bufs: Vec::new(),
            pool,
            be,
        })
    }

    fn b(&self, h: Buf) -> &[f64] {
        &self.bufs[h.0 as usize]
    }

    fn v3(&self, h: Buf3) -> Vec3<'_, f64> {
        Vec3 {
            x: self.b(h.x),
            y: self.b(h.y),
            z: self.b(h.z),
        }
    }

    fn grad(&self, g: &[Buf3; 3]) -> VelocityGradient<'_, f64> {
        VelocityGradient {
            x: self.v3(g[0]),
            y: self.v3(g[1]),
            z: self.v3(g[2]),
        }
    }

    fn geometry(&self) -> FaceGeometry<'_, f64> {
        let d = &self.disc;
        FaceGeometry {
            area: d.face_area.view(),
            face_centroid: d.face_centroid.view(),
            cell_centroid: d.cell_centroid.view(),
            dgeo: &d.dgeo,
            k: d.k.view(),
            weight: &d.weight,
        }
    }

    fn layout(&self) -> SellLayout<'_> {
        SellLayout {
            chunk: self.disc.sell.chunk,
            slice_ptr: &self.disc.sell.slice_ptr,
        }
    }

    /// Takes the output buffers out of the table (they must differ from each
    /// other and from `inputs`).
    fn take(&mut self, outs: &[Buf], inputs: &[Buf]) -> Result<Vec<Vec<f64>>> {
        for (i, o) in outs.iter().enumerate() {
            ensure!(
                !inputs.contains(o) && !outs[..i].contains(o),
                "flow device: output buffer {o:?} is also an input"
            );
        }
        Ok(outs
            .iter()
            .map(|o| std::mem::take(&mut self.bufs[o.0 as usize]))
            .collect())
    }

    fn put(&mut self, outs: &[Buf], v: Vec<Vec<f64>>) {
        for (o, x) in outs.iter().zip(v) {
            self.bufs[o.0 as usize] = x;
        }
    }

    fn mask(&self, m: Mask) -> &[f64] {
        &self.masks[match m {
            Mask::Velocity => 0,
            Mask::Pressure => 1,
            Mask::PressureCorrection => 2,
        }]
    }

    /// Boundary constants of `value` with mask `m` into the scratch array.
    fn constants(&mut self, m: Mask, value: Buf) {
        let mut s = std::mem::take(&mut self.scratch);
        {
            let a = self.mask(m);
            let v = self.b(value);
            // SAFETY: fb is a boundary face; both arrays have an entry each.
            par_vec(&self.pool, &mut s, |fb, _| unsafe {
                k::boundary_constant(fb, a, v)
            });
        }
        self.scratch = s;
    }
}

impl FlowDevice for CpuFlow {
    fn name(&self) -> String {
        "cpu".into()
    }

    fn device(&self) -> String {
        format!("CPU, {} threads", self.pool.current_num_threads())
    }

    fn alloc(&mut self, n: usize) -> Result<Buf> {
        self.bufs.push(vec![0.0; n]);
        Ok(Buf(self.bufs.len() as u32 - 1))
    }

    fn upload(&mut self, b: Buf, v: &[f64]) -> Result<()> {
        let dst = &mut self.bufs[b.0 as usize];
        ensure!(
            dst.len() == v.len(),
            "upload of {} values into {}",
            v.len(),
            dst.len()
        );
        dst.copy_from_slice(v);
        Ok(())
    }

    fn download(&mut self, b: Buf) -> Result<Vec<f64>> {
        Ok(self.b(b).to_vec())
    }

    fn copy(&mut self, src: Buf, dst: Buf) -> Result<()> {
        if src != dst {
            let s = std::mem::take(&mut self.bufs[src.0 as usize]);
            self.bufs[dst.0 as usize].copy_from_slice(&s);
            self.bufs[src.0 as usize] = s;
        }
        Ok(())
    }

    fn fill(&mut self, b: Buf, v: f64) -> Result<()> {
        self.bufs[b.0 as usize].fill(v);
        Ok(())
    }

    fn scaled_difference(&mut self, out: Buf, a: f64, x: Buf, y: Buf) -> Result<()> {
        let mut o = std::mem::take(&mut self.bufs[out.0 as usize]);
        {
            let (xs, ys) = (self.b(x), self.b(y));
            par_vec(&self.pool, &mut o, |i, oi| {
                let xi = if x == out { oi } else { xs[i] };
                let yi = if y == out { oi } else { ys[i] };
                a * (xi - yi)
            });
        }
        self.bufs[out.0 as usize] = o;
        Ok(())
    }

    fn lincomb(&mut self, out: Buf, a: f64, x: Buf, b: f64, y: Buf) -> Result<()> {
        let mut o = std::mem::take(&mut self.bufs[out.0 as usize]);
        {
            let (xs, ys) = (self.b(x), self.b(y));
            par_vec(&self.pool, &mut o, |i, oi| {
                let xi = if x == out { oi } else { xs[i] };
                let yi = if y == out { oi } else { ys[i] };
                a * xi + b * yi
            });
        }
        self.bufs[out.0 as usize] = o;
        Ok(())
    }

    fn sum_sq(&mut self, x: Buf) -> Result<f64> {
        Ok(sum_sq(self.b(x)))
    }

    fn min_max(&mut self, x: Buf) -> Result<(f64, f64)> {
        let v = self.b(x);
        Ok((
            v.iter().copied().fold(f64::INFINITY, f64::min),
            v.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        ))
    }

    fn value(&mut self, b: Buf, i: usize) -> Result<f64> {
        Ok(self.b(b)[i])
    }

    fn subtract(&mut self, b: Buf, s: f64) -> Result<()> {
        self.bufs[b.0 as usize].iter_mut().for_each(|v| *v -= s);
        Ok(())
    }

    fn scale_entry(&mut self, b: Buf, i: usize, s: f64) -> Result<()> {
        self.bufs[b.0 as usize][i] *= s;
        Ok(())
    }

    fn boundary_values(&mut self, a: &BoundaryArgs) -> Result<()> {
        let nb = self.disc.n_faces - self.disc.n_interior;
        let outs = [a.ub.x, a.ub.y, a.ub.z, a.pb];
        let mut ins = vec![a.p, a.mass_flux];
        ins.extend(a.u.all());
        ins.extend(a.grad_p.all());
        ins.extend(a.ub_set.all());
        for g in &a.grad_u {
            ins.extend(g.all());
        }
        let mut o = self.take(&outs, &ins)?;
        {
            let inp = k::BoundaryInputs {
                topo: self.disc.topology(),
                geo: self.geometry(),
                rho: self.rho,
                mass_flux: self.b(a.mass_flux),
            };
            let (u, g, ub0) = (self.v3(a.u), self.grad(&a.grad_u), self.v3(a.ub_set));
            let (px, py, pz) = (Ptr::of(&mut o[0]), Ptr::of(&mut o[1]), Ptr::of(&mut o[2]));
            let bc = &self.bc;
            par_for(&self.pool, nb, |fb| {
                // SAFETY: fb < nb; the arrays were built for this mesh.
                unsafe {
                    let v = k::boundary_velocity(fb, &inp, bc.ub_rule[fb], ub0.at(fb), &u, &g);
                    *px.at(fb) = v[0];
                    *py.at(fb) = v[1];
                    *pz.at(fb) = v[2];
                }
            });
            let (p, gp) = (self.b(a.p), self.v3(a.grad_p));
            par_vec(&self.pool, &mut o[3], |fb, _| {
                // SAFETY: as above.
                unsafe { k::boundary_pressure(fb, &inp, bc.pb_rule[fb], bc.pset[fb], p, &gp) }
            });
        }
        self.put(&outs, o);
        Ok(())
    }

    fn gradient(&mut self, w: Lsq, phi: Buf, mask: Mask, value: Buf, out: Buf3) -> Result<()> {
        self.constants(mask, value);
        let outs = out.all();
        let mut o = self.take(&outs, &[phi])?;
        {
            let topo = self.disc.topology();
            let w = match w {
                Lsq::Velocity => self.disc.lsq.view(),
                Lsq::Pressure => self.lsq_p.view(),
            };
            let bc = fv::Boundary {
                a: self.mask(mask),
                b: &self.scratch,
            };
            let phi = self.b(phi);
            let (px, py, pz) = (Ptr::of(&mut o[0]), Ptr::of(&mut o[1]), Ptr::of(&mut o[2]));
            par_for(&self.pool, self.disc.n_cells, |c| {
                // SAFETY: c is a local cell; the weights were built for this mesh.
                unsafe {
                    let g = fv::lsq_gradient(c, &topo, &w, phi, &bc);
                    *px.at(c) = g[0];
                    *py.at(c) = g[1];
                    *pz.at(c) = g[2];
                }
            });
        }
        self.put(&outs, o);
        Ok(())
    }

    fn limiter(&mut self, phi: Buf, value: Buf, grad: Buf3, eps2: Eps2, out: Buf) -> Result<()> {
        self.constants(Mask::Velocity, value);
        let mut ins = vec![phi];
        ins.extend(grad.all());
        let mut o = self.take(&[out], &ins)?;
        {
            let single;
            let e2: &[f64] = match eps2 {
                Eps2::PerCell => &self.eps2,
                Eps2::Single(e) => {
                    single = [e];
                    &single
                }
            };
            let topo = self.disc.topology();
            let bc = fv::Boundary {
                a: self.mask(Mask::Velocity),
                b: &self.scratch,
            };
            let (phi, g) = (self.b(phi), self.v3(grad));
            let (cc, fc) = (
                self.disc.cell_centroid.view(),
                self.disc.face_centroid.view(),
            );
            par_vec(&self.pool, &mut o[0], |c, _| {
                // SAFETY: as in `gradient`.
                unsafe { fv::venkat_limiter(c, &topo, phi, &bc, &g, &cc, &fc, e2) }
            });
        }
        self.put(&[out], o);
        Ok(())
    }

    fn mass_flux(&mut self, a: &FluxArgs, out: Buf) -> Result<()> {
        let mut ins = vec![a.p, a.pb, a.relax];
        ins.extend(a.u.all());
        ins.extend(a.grad_p.all());
        ins.extend(a.ub.all());
        for g in &a.grad_u {
            ins.extend(g.all());
        }
        ins.extend(a.d);
        if let Some((ds, y, _)) = a.time {
            ins.push(ds);
            ins.push(y);
        }
        let mut o = self.take(&[out], &ins)?;
        {
            let (ds, time, alpha) = match a.time {
                Some((ds, y, al)) => (self.b(ds), self.b(y), al),
                None => (&[][..], &[][..], 1.0),
            };
            let inp = k::FluxInputs {
                topo: self.disc.topology(),
                geo: self.geometry(),
                rho: self.rho,
                u: self.v3(a.u),
                grad: self.grad(&a.grad_u),
                p: self.b(a.p),
                grad_p: self.v3(a.grad_p),
                pb: self.b(a.pb),
                ub: self.v3(a.ub),
                flux_kind: &self.bc.flux_kind,
                mom_kind: &self.bc.mom_kind,
                d: match a.d {
                    Some(d) => self.b(d),
                    None => &self.zero_cells,
                },
                relax: self.b(a.relax),
                transient: a.time.is_some(),
                ds,
                time,
                alpha,
            };
            par_vec(&self.pool, &mut o[0], |f, _| {
                // SAFETY: f is a face; all arrays were built for this mesh.
                unsafe { k::face_mass_flux(f, &inp) }
            });
        }
        self.put(&[out], o);
        Ok(())
    }

    fn momentum_faces(&mut self, a: &MomentumFaceArgs) -> Result<()> {
        let nf = self.disc.n_faces;
        if self.mu_face.len() != nf {
            // the viscosity is constant in laminar flow (stage 6: per face)
            self.mu_face = vec![self.mu; nf];
        }
        let outs = [a.cp, a.cn, a.bf.x, a.bf.y, a.bf.z];
        let mut ins = vec![a.mass_flux];
        ins.extend(a.u.all());
        ins.extend(a.psi.all());
        ins.extend(a.ub.all());
        for g in &a.grad_u {
            ins.extend(g.all());
        }
        let mut o = self.take(&outs, &ins)?;
        {
            let inp = k::MomentumInputs {
                topo: self.disc.topology(),
                geo: self.geometry(),
                mass_flux: self.b(a.mass_flux),
                mu: &self.mu_face,
                u: self.v3(a.u),
                grad: self.grad(&a.grad_u),
                psi: self.v3(a.psi),
                ub: self.v3(a.ub),
                mom_kind: &self.bc.mom_kind,
                blend: a.blend,
                nonorth: a.nonorth,
                convection: a.convection,
            };
            let p: Vec<Ptr> = o.iter_mut().map(|v| Ptr::of(v)).collect();
            par_for(&self.pool, nf, |f| {
                // SAFETY: f is a face of the mesh.
                unsafe {
                    let (cp, cn, b) = k::momentum_face(f, &inp);
                    *p[0].at(f) = cp;
                    *p[1].at(f) = cn;
                    *p[2].at(f) = b[0];
                    *p[3].at(f) = b[1];
                    *p[4].at(f) = b[2];
                }
            });
        }
        self.put(&outs, o);
        Ok(())
    }

    fn momentum_rows(&mut self, a: &MomentumRowArgs) -> Result<()> {
        let outs = [
            a.val, a.rhs.x, a.rhs.y, a.rhs.z, a.d, a.ds, a.dc, a.diag.x, a.diag.y, a.diag.z,
        ];
        let mut ins = vec![a.cp, a.cn];
        ins.extend(a.bf.all());
        ins.extend(a.grad_p.all());
        ins.extend(a.u.all());
        if let Some(t) = &a.time {
            ins.extend(t.u_old.all());
            ins.extend(t.u_old2.all());
            ins.extend(t.explicit.1.all());
        }
        let mut o = self.take(&outs, &ins)?;
        {
            let empty = Vec3 {
                x: &[][..],
                y: &[][..],
                z: &[][..],
            };
            let (transient, tau, c_old, u_old, u_old2, ew, ex) = match &a.time {
                Some(t) => (
                    true,
                    t.tau,
                    t.c_old,
                    self.v3(t.u_old),
                    self.v3(t.u_old2),
                    t.explicit.0,
                    self.v3(t.explicit.1),
                ),
                None => (false, 0.0, [0.0; 2], empty, empty, 0.0, empty),
            };
            let inp = k::MomentumRowInputs {
                topo: self.disc.topology(),
                cp: self.b(a.cp),
                cn: self.b(a.cn),
                bf: self.v3(a.bf),
                grad_p: self.v3(a.grad_p),
                body_force: a.body_force,
                volume: &self.disc.cell_volume,
                u: self.v3(a.u),
                alpha: a.alpha,
                layout: self.layout(),
                mom_kind: &self.bc.mom_kind,
                face_area: self.disc.face_area.view(),
                theta: a.theta,
                transient,
                tau,
                c_old,
                u_old,
                u_old2,
                explicit_weight: ew,
                explicit: ex,
            };
            let p: Vec<Ptr> = o.iter_mut().map(|v| Ptr::of(v)).collect();
            par_for(&self.pool, self.disc.n_cells, |c| {
                let out = k::MomentumRowOut {
                    val: p[0].get(),
                    rhs: [p[1].get(), p[2].get(), p[3].get()],
                    d: p[4].get(),
                    ds: p[5].get(),
                    dc: p[6].get(),
                    diag: [p[7].get(), p[8].get(), p[9].get()],
                };
                // SAFETY: c is a local cell; row c is written only here.
                unsafe { k::momentum_row(c, &inp, &out) }
            });
        }
        self.put(&outs, o);
        Ok(())
    }

    fn momentum_operator(&mut self, cp: Buf, cn: Buf, bf: Buf3, u: Buf3, out: Buf3) -> Result<()> {
        let outs = out.all();
        let mut ins = vec![cp, cn];
        ins.extend(bf.all());
        ins.extend(u.all());
        let mut o = self.take(&outs, &ins)?;
        {
            let inp = k::OperatorInputs {
                topo: self.disc.topology(),
                cp: self.b(cp),
                cn: self.b(cn),
                bf: self.v3(bf),
                mom_kind: &self.bc.mom_kind,
                face_area: self.disc.face_area.view(),
                u: self.v3(u),
            };
            let p: Vec<Ptr> = o.iter_mut().map(|v| Ptr::of(v)).collect();
            par_for(&self.pool, self.disc.n_cells, |c| {
                // SAFETY: c is a local cell, written once.
                unsafe {
                    let r = k::momentum_operator(c, &inp);
                    *p[0].at(c) = r[0];
                    *p[1].at(c) = r[1];
                    *p[2].at(c) = r[2];
                }
            });
        }
        self.put(&outs, o);
        Ok(())
    }

    fn residual_scale(&mut self, u: Buf3, rhs: Buf3, diag: Buf3, alpha: f64) -> Result<f64> {
        let nc = self.disc.n_cells;
        let mut v = vec![0.0; 3 * nc];
        let mut w = vec![0.0; 3 * nc];
        {
            let (u, b, dg) = (self.v3(u), self.v3(rhs), self.v3(diag));
            let (pv, pw) = (Ptr::of(&mut v), Ptr::of(&mut w));
            par_for(&self.pool, nc, |c| {
                // SAFETY: c is a local cell; entries 3c..3c+3 written once.
                unsafe {
                    let (a, b) = k::residual_scale(c, &u, &b, &dg, alpha);
                    for i in 0..3 {
                        *pv.at(3 * c + i) = a[i];
                        *pw.at(3 * c + i) = b[i];
                    }
                }
            });
        }
        Ok((sum_sq(&v).sqrt() + sum_sq(&w).sqrt()).max(f64::MIN_POSITIVE))
    }

    fn set_diagonal(&mut self, val: Buf, diag: Buf) -> Result<()> {
        let mut o = self.take(&[val], &[diag])?;
        {
            let sell = &self.disc.sell;
            let dg = self.b(diag);
            for (c, &d) in dg.iter().enumerate() {
                o[0][sell.pos(c, 0)] = d;
            }
        }
        self.put(&[val], o);
        Ok(())
    }

    fn net_outflow(&mut self, m: Buf, out: Buf) -> Result<()> {
        let mut o = self.take(&[out], &[m])?;
        {
            let (topo, m) = (self.disc.topology(), self.b(m));
            // SAFETY: c is a local cell, m has one entry per face.
            par_vec(&self.pool, &mut o[0], |c, _| unsafe {
                k::net_outflow(c, &topo, m)
            });
        }
        self.put(&[out], o);
        Ok(())
    }

    fn throughflow(&mut self, m: Buf, out: Buf) -> Result<()> {
        let mut o = self.take(&[out], &[m])?;
        {
            let (topo, m) = (self.disc.topology(), self.b(m));
            // SAFETY: as above.
            par_vec(&self.pool, &mut o[0], |c, _| unsafe {
                k::throughflow(c, &topo, m)
            });
        }
        self.put(&[out], o);
        Ok(())
    }

    fn courant(&mut self, m: Buf, dt: f64, out: Buf) -> Result<()> {
        let mut o = self.take(&[out], &[m])?;
        {
            let (topo, m, vol) = (self.disc.topology(), self.b(m), &self.disc.cell_volume);
            let factor = dt / (2.0 * self.rho);
            // SAFETY: as above.
            par_vec(&self.pool, &mut o[0], |c, _| unsafe {
                k::courant(c, &topo, m, vol, factor)
            });
        }
        self.put(&[out], o);
        Ok(())
    }

    fn pressure_faces(&mut self, dc: Buf, out: Buf) -> Result<()> {
        let mut o = self.take(&[out], &[dc])?;
        {
            let (topo, geo) = (self.disc.topology(), self.geometry());
            let (fk, dc, rho) = (&self.bc.flux_kind, self.b(dc), self.rho);
            // SAFETY: f is a face.
            par_vec(&self.pool, &mut o[0], |f, _| unsafe {
                k::pressure_face(f, &topo, &geo, fk, rho, dc)
            });
        }
        self.put(&[out], o);
        Ok(())
    }

    fn pressure_rows(&mut self, cf: Buf, m: Buf, val: Buf, rhs: Buf) -> Result<()> {
        let outs = [val, rhs];
        let mut o = self.take(&outs, &[cf, m])?;
        {
            let (topo, layout) = (self.disc.topology(), self.layout());
            let (pv, pr) = (Ptr::of(&mut o[0]), Ptr::of(&mut o[1]));
            let (cf, m) = (self.b(cf), self.b(m));
            par_for(&self.pool, self.disc.n_cells, |c| {
                // SAFETY: c is a local cell; row c is written only here.
                unsafe { k::pressure_row(c, &topo, cf, m, &layout, pv.get(), pr.get()) }
            });
        }
        self.put(&outs, o);
        Ok(())
    }

    fn correct_fluxes(&mut self, cf: Buf, m: Buf, pc: Buf, out: Buf) -> Result<()> {
        let mut o = self.take(&[out], &[cf, m, pc])?;
        {
            let topo = self.disc.topology();
            let (cf, m, pc) = (self.b(cf), self.b(m), self.b(pc));
            // SAFETY: f is a face; cf, m have one entry per face, pc per cell.
            par_vec(&self.pool, &mut o[0], |f, _| unsafe {
                k::correct_flux(f, &topo, cf, m, pc)
            });
        }
        self.put(&[out], o);
        Ok(())
    }

    fn correct_velocity(&mut self, dc: Buf, grad: Buf3, u: Buf3) -> Result<()> {
        let outs = u.all();
        let mut ins = vec![dc];
        ins.extend(grad.all());
        let mut o = self.take(&outs, &ins)?;
        {
            let (dc, g) = (self.b(dc), self.v3(grad));
            let p: Vec<Ptr> = o.iter_mut().map(|v| Ptr::of(v)).collect();
            par_for(&self.pool, self.disc.n_cells, |c| {
                // SAFETY: c is a local cell, written once.
                unsafe {
                    *p[0].at(c) -= dc[c] * g.x[c];
                    *p[1].at(c) -= dc[c] * g.y[c];
                    *p[2].at(c) -= dc[c] * g.z[c];
                }
            });
        }
        self.put(&outs, o);
        Ok(())
    }

    fn neighbour_correction(&mut self, val: Buf, x: Buf, diag: Buf, out: Buf) -> Result<()> {
        let mut o = self.take(&[out], &[val, x, diag])?;
        {
            let (topo, layout) = (self.disc.topology(), self.layout());
            let (val, x, dg) = (self.b(val), self.b(x), self.b(diag));
            // SAFETY: c is a local cell; val is the SELL storage of the mesh.
            par_vec(&self.pool, &mut o[0], |c, _| unsafe {
                -k::neighbour_sum(c, &topo, &layout, val, x) / dg[c]
            });
        }
        self.put(&[out], o);
        Ok(())
    }

    fn delta_flux(&mut self, delta: Buf3, out: Buf) -> Result<()> {
        let mut o = self.take(&[out], &delta.all())?;
        {
            let (topo, geo) = (self.disc.topology(), self.geometry());
            let (fk, rho, d) = (&self.bc.flux_kind, self.rho, self.v3(delta));
            // SAFETY: f is a face, delta has an entry per cell.
            par_vec(&self.pool, &mut o[0], |f, _| unsafe {
                k::delta_flux(f, &topo, &geo, fk, rho, &d)
            });
        }
        self.put(&[out], o);
        Ok(())
    }

    fn solve(
        &mut self,
        val: Buf,
        rhs: Buf,
        x: Buf,
        pc: Preconditioner,
        cg: bool,
        tol: f64,
        max_iter: usize,
    ) -> Result<(usize, f64)> {
        let be = self.be.as_mut();
        be.set_matrix(&self.bufs[val.0 as usize])?;
        be.upload(V::Rhs, &self.bufs[rhs.0 as usize])?;
        be.upload(V::Phi, &self.bufs[x.0 as usize])?;
        be.residual(V::Phi, V::Rhs, V::R)?;
        let r0 = be.dot(V::R, V::R)?.sqrt();
        let it = if r0 == 0.0 {
            0
        } else if cg {
            linear::cg(be, pc, max_iter, tol, r0)?.0
        } else {
            linear::bicgstab(be, pc, max_iter, tol, r0)?
        };
        self.bufs[x.0 as usize].copy_from_slice(&be.download(V::Phi)?);
        Ok((it, r0))
    }

    fn amg_setup(&mut self, val: Buf, s: &AmgSettings) -> Result<()> {
        self.be.set_matrix(&self.bufs[val.0 as usize])?;
        linear::setup_amg(self.be.as_mut(), s).map(|_| ())
    }

    fn sync(&mut self) -> Result<()> {
        Ok(())
    }
}

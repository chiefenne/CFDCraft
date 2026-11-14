//! CPU backend: the kernel bodies in multi-threaded loops (rayon).
//!
//! Each output value is computed by one iteration from its index, so the
//! results are bitwise the same for any number of threads (D31).

use super::{Backend, FaceSettings, Limiter, V};
use crate::disc::Discretisation;
use cfdcraft_core::repro::{BLOCK, combine_block_sums};
use cfdcraft_kernels::fv::{self, FaceInputs, SellLayout, Vec3};
use cfdcraft_kernels::{blas, spmv};
use rayon::prelude::*;

/// Indices per task of the parallel loops.
const GRAIN: usize = 1 << 14;

/// CPU backend.
pub struct CpuBackend {
    d: Discretisation,
    pool: rayon::ThreadPool,
    vecs: Vec<Vec<f64>>,
    scratch: Vec<f64>,
    val: Vec<f64>,
    eps2: Option<(f64, Vec<f64>)>,
}

/// A raw pointer that may be shared by the threads of a gather loop in
/// which every index writes its own, disjoint outputs.
#[derive(Clone, Copy)]
struct Shared(*mut f64);
// SAFETY: used only where each thread writes disjoint entries.
unsafe impl Send for Shared {}
// SAFETY: as above.
unsafe impl Sync for Shared {}

impl CpuBackend {
    /// New backend on `threads` threads (`None`: all cores).
    pub fn new(d: Discretisation, threads: Option<usize>) -> anyhow::Result<Self> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.unwrap_or(0))
            .build()?;
        let mut vecs = vec![Vec::new(); V::COUNT];
        for v in V::CELL {
            vecs[v.index()] = vec![0.0; d.n_cells];
        }
        for v in V::FACE {
            vecs[v.index()] = vec![0.0; d.n_faces];
        }
        vecs[V::Psi.index()].fill(1.0);
        let val = vec![0.0; d.sell.stored()];
        Ok(Self {
            d,
            pool,
            vecs,
            scratch: Vec::new(),
            val,
            eps2: None,
        })
    }

    fn v(&self, v: V) -> &[f64] {
        &self.vecs[v.index()]
    }

    /// Computes `out[i] = body(self, i)` into the scratch buffer and swaps it
    /// in, so `out` may also be an input.
    fn map_into(&mut self, out: V, body: impl Fn(&Self, usize) -> f64 + Sync) {
        let n = self.vecs[out.index()].len();
        let mut s = std::mem::take(&mut self.scratch);
        s.resize(n, 0.0);
        let this = &*self;
        this.pool.install(|| {
            s.par_chunks_mut(GRAIN).enumerate().for_each(|(c, ch)| {
                for (j, o) in ch.iter_mut().enumerate() {
                    *o = body(this, c * GRAIN + j);
                }
            })
        });
        std::mem::swap(&mut s, &mut self.vecs[out.index()]);
        self.scratch = s;
    }

    fn grad_view(&self) -> Vec3<'_, f64> {
        Vec3 {
            x: self.v(V::Gx),
            y: self.v(V::Gy),
            z: self.v(V::Gz),
        }
    }

    fn sell_row(&self, r: usize, x: &[f64]) -> f64 {
        let s = &self.d.sell;
        // SAFETY: the structure was built by `sell_structure` for this
        // mesh; every column is a cell index.
        unsafe { spmv::sell_row(r, s.chunk, &s.slice_ptr, &s.slice_len, &s.col, &self.val, x) }
    }
}

impl Backend for CpuBackend {
    fn name(&self) -> String {
        if self.pool.current_num_threads() == 1 {
            "cpu-serial".into()
        } else {
            "cpu".into()
        }
    }

    fn device(&self) -> String {
        format!("CPU, {} threads", self.pool.current_num_threads())
    }

    fn n_cells(&self) -> usize {
        self.d.n_cells
    }

    fn upload(&mut self, v: V, data: &[f64]) -> anyhow::Result<()> {
        let t = &mut self.vecs[v.index()];
        anyhow::ensure!(t.len() == data.len(), "upload of {v:?}: wrong length");
        t.copy_from_slice(data);
        Ok(())
    }

    fn download(&mut self, v: V) -> anyhow::Result<Vec<f64>> {
        Ok(self.vecs[v.index()].clone())
    }

    fn sync(&mut self) -> anyhow::Result<()> {
        Ok(())
    }

    fn gradient(&mut self) -> anyhow::Result<()> {
        let mut gx = std::mem::take(&mut self.vecs[V::Gx.index()]);
        let mut gy = std::mem::take(&mut self.vecs[V::Gy.index()]);
        let mut gz = std::mem::take(&mut self.vecs[V::Gz.index()]);
        let (d, phi) = (&self.d, &self.vecs[V::Phi.index()]);
        let (topo, bc, w) = (d.topology(), d.boundary(), d.lsq.view());
        self.pool.install(|| {
            gx.par_chunks_mut(GRAIN)
                .zip(gy.par_chunks_mut(GRAIN))
                .zip(gz.par_chunks_mut(GRAIN))
                .enumerate()
                .for_each(|(ci, ((x, y), z))| {
                    for j in 0..x.len() {
                        // SAFETY: the discretisation was built from a
                        // validated mesh; c is a local cell.
                        let g = unsafe { fv::lsq_gradient(ci * GRAIN + j, &topo, &w, phi, &bc) };
                        x[j] = g[0];
                        y[j] = g[1];
                        z[j] = g[2];
                    }
                })
        });
        self.vecs[V::Gx.index()] = gx;
        self.vecs[V::Gy.index()] = gy;
        self.vecs[V::Gz.index()] = gz;
        Ok(())
    }

    fn limiter(&mut self, l: Limiter) -> anyhow::Result<()> {
        let eps2: Vec<f64> = match l {
            Limiter::None => {
                self.vecs[V::Psi.index()].fill(1.0);
                return Ok(());
            }
            Limiter::Venkatakrishnan { k } => match &self.eps2 {
                Some((kk, e)) if *kk == k => e.clone(),
                _ => {
                    let e = self.d.venkat_eps2(k);
                    self.eps2 = Some((k, e.clone()));
                    e
                }
            },
            Limiter::VenkatakrishnanWang { k } => {
                let (mn, mx) = self.min_max(V::Phi)?;
                let e = k * (mx - mn);
                vec![e * e]
            }
        };
        self.map_into(V::Psi, |s, c| {
            let d = &s.d;
            // SAFETY: as in `gradient`; eps2 has one entry or one per cell.
            unsafe {
                fv::venkat_limiter(
                    c,
                    &d.topology(),
                    s.v(V::Phi),
                    &d.boundary(),
                    &s.grad_view(),
                    &d.cell_centroid.view(),
                    &d.face_centroid.view(),
                    &eps2,
                )
            }
        });
        Ok(())
    }

    fn faces(&mut self, st: FaceSettings) -> anyhow::Result<()> {
        let mut cp = std::mem::take(&mut self.vecs[V::Cp.index()]);
        let mut cn = std::mem::take(&mut self.vecs[V::Cn.index()]);
        let mut bf = std::mem::take(&mut self.vecs[V::Bf.index()]);
        let d = &self.d;
        let inp = FaceInputs {
            topo: d.topology(),
            bc: d.boundary(),
            mass_flux: &d.mass_flux,
            gamma: &d.gamma,
            dgeo: &d.dgeo,
            k: d.k.view(),
            weight: &d.weight,
            face_centroid: d.face_centroid.view(),
            cell_centroid: d.cell_centroid.view(),
            phi: self.v(V::Phi),
            grad: self.grad_view(),
            psi: self.v(V::Psi),
            blend: st.blend,
            nonorth: if st.nonorth { 1.0 } else { 0.0 },
            nonorth_limit: st.nonorth_limit,
        };
        self.pool.install(|| {
            cp.par_chunks_mut(GRAIN)
                .zip(cn.par_chunks_mut(GRAIN))
                .zip(bf.par_chunks_mut(GRAIN))
                .enumerate()
                .for_each(|(ci, ((a, b), c))| {
                    for j in 0..a.len() {
                        // SAFETY: f is a face of the discretisation.
                        let (x, y, z) = unsafe { fv::face_coefficients(ci * GRAIN + j, &inp) };
                        a[j] = x;
                        b[j] = y;
                        c[j] = z;
                    }
                })
        });
        self.vecs[V::Cp.index()] = cp;
        self.vecs[V::Cn.index()] = cn;
        self.vecs[V::Bf.index()] = bf;
        Ok(())
    }

    fn assemble(&mut self, alpha: f64) -> anyhow::Result<()> {
        let mut rhs = std::mem::take(&mut self.vecs[V::Rhs.index()]);
        let mut dinv = std::mem::take(&mut self.vecs[V::Dinv.index()]);
        let (pv, pr, pd) = (
            Shared(self.val.as_mut_ptr()),
            Shared(rhs.as_mut_ptr()),
            Shared(dinv.as_mut_ptr()),
        );
        let d = &self.d;
        let topo = d.topology();
        let layout = SellLayout {
            chunk: d.sell.chunk,
            slice_ptr: &d.sell.slice_ptr,
        };
        let (cp, cn, bf, phi) = (self.v(V::Cp), self.v(V::Cn), self.v(V::Bf), self.v(V::Phi));
        let n = d.n_cells;
        self.pool.install(|| {
            (0..n.div_ceil(GRAIN)).into_par_iter().for_each(|ci| {
                let (pv, pr, pd) = (pv, pr, pd);
                for c in ci * GRAIN..((ci + 1) * GRAIN).min(n) {
                    // SAFETY: row c and rhs[c], dinv[c] are written by this
                    // iteration only; the SELL storage has a slot for every
                    // interior face of c.
                    unsafe {
                        fv::assemble_row(
                            c, &topo, cp, cn, bf, &d.source, phi, alpha, &layout, pv.0, pr.0, pd.0,
                        )
                    }
                }
            })
        });
        self.vecs[V::Rhs.index()] = rhs;
        self.vecs[V::Dinv.index()] = dinv;
        Ok(())
    }

    fn spmv(&mut self, x: V, out: V) -> anyhow::Result<()> {
        self.map_into(out, |s, r| s.sell_row(r, s.v(x)));
        Ok(())
    }

    fn residual(&mut self, x: V, b: V, out: V) -> anyhow::Result<()> {
        self.map_into(out, |s, r| s.v(b)[r] - s.sell_row(r, s.v(x)));
        Ok(())
    }

    fn dot(&mut self, a: V, b: V) -> anyhow::Result<f64> {
        let (x, y) = (self.v(a), self.v(b));
        let nb = x.len().div_ceil(BLOCK);
        let sums: Vec<f64> = self.pool.install(|| {
            (0..nb)
                .into_par_iter()
                // SAFETY: blk * BLOCK < len.
                .map(|blk| unsafe { blas::block_dot(blk, x, y) })
                .collect()
        });
        Ok(combine_block_sums(&sums))
    }

    fn min_max(&mut self, a: V) -> anyhow::Result<(f64, f64)> {
        let x = self.v(a);
        let nb = x.len().div_ceil(BLOCK);
        let mm: Vec<(f64, f64)> = self.pool.install(|| {
            (0..nb)
                .into_par_iter()
                // SAFETY: blk * BLOCK < len.
                .map(|blk| unsafe { blas::block_min_max(blk, x) })
                .collect()
        });
        Ok(mm
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &(c, d)| {
                (a.min(c), b.max(d))
            }))
    }

    fn update_p(&mut self, beta: f64, omega: f64) -> anyhow::Result<()> {
        // SAFETY: all cell vectors have the same length.
        self.map_into(V::P, |s, i| unsafe {
            blas::update_p(i, s.v(V::R), s.v(V::P), s.v(V::Vv), beta, omega)
        });
        Ok(())
    }

    fn sub_scaled(&mut self, x: V, z: V, alpha: f64, out: V) -> anyhow::Result<()> {
        // SAFETY: as above.
        self.map_into(out, |s, i| unsafe {
            blas::sub_scaled(i, s.v(x), s.v(z), alpha)
        });
        Ok(())
    }

    fn update_x(&mut self, alpha: f64, omega: f64) -> anyhow::Result<()> {
        // SAFETY: as above.
        self.map_into(V::Phi, |s, i| unsafe {
            blas::update_x(i, s.v(V::Phi), s.v(V::Ph), s.v(V::Sh), alpha, omega)
        });
        Ok(())
    }

    fn precondition(&mut self, x: V, out: V) -> anyhow::Result<()> {
        // SAFETY: as above.
        self.map_into(out, |s, i| unsafe { blas::scale(i, s.v(V::Dinv), s.v(x)) });
        Ok(())
    }

    fn copy(&mut self, src: V, dst: V) -> anyhow::Result<()> {
        let s = self.vecs[src.index()].clone();
        self.vecs[dst.index()] = s;
        Ok(())
    }
}

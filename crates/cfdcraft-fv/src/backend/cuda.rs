//! NVIDIA GPU backend (feature `cuda`): the kernel bodies of
//! `cfdcraft_kernels` compiled by rustc to PTX (`gpu/ptx-kernels`, entry
//! points `fv_*` and `blas_*`), launched through the CUDA driver. All data
//! stay on the GPU; only dot products (one value per block of 256) and the
//! limiter range come back to the host.
//!
//! The PTX has no fused multiply-adds and IEEE-rounded division, and the
//! block sums follow the CPU order, so CPU and GPU results are expected to
//! agree bitwise (checked by `cfdcraft scalar` with `backends = ["cpu",
//! "cuda"]`).

use super::{Backend, FaceSettings, Limiter, V};
use crate::disc::Discretisation;
use anyhow::Context as _;
use cfdcraft_core::repro::{BLOCK, combine_block_sums};
use cudarc::driver::{
    CudaContext, CudaFunction, CudaModule, CudaSlice, CudaStream, LaunchConfig, PushKernelArg,
};
use cudarc::nvrtc::Ptx;
use std::path::Path;
use std::sync::Arc;

/// Threads per block.
const THREADS: u32 = 256;

const ENTRIES: [&str; 12] = [
    "fv_gradient",
    "fv_limiter",
    "fv_faces",
    "fv_assemble",
    "fv_residual",
    "g1_spmv_sell",
    "blas_block_dot",
    "blas_block_min_max",
    "blas_update_p",
    "blas_sub_scaled",
    "blas_update_x",
    "blas_scale",
];

fn cfg(n: usize) -> LaunchConfig {
    let n = n as u32;
    LaunchConfig {
        grid_dim: (n.div_ceil(THREADS).max(1), 1, 1),
        block_dim: (THREADS, 1, 1),
        shared_mem_bytes: 0,
    }
}

struct Soa {
    x: CudaSlice<f64>,
    y: CudaSlice<f64>,
    z: CudaSlice<f64>,
}

/// Device copy of the discretisation.
struct Dev {
    owner: CudaSlice<u32>,
    neighbour: CudaSlice<u32>,
    cell_ptr: CudaSlice<u32>,
    cell_faces: CudaSlice<u32>,
    cc: Soa,
    fc: Soa,
    k: Soa,
    lsq: Soa,
    dgeo: CudaSlice<f64>,
    weight: CudaSlice<f64>,
    bc_a: CudaSlice<f64>,
    bc_b: CudaSlice<f64>,
    mass_flux: CudaSlice<f64>,
    gamma: CudaSlice<f64>,
    source: CudaSlice<f64>,
    slice_ptr: CudaSlice<u32>,
    slice_len: CudaSlice<u32>,
    col: CudaSlice<u32>,
    val: CudaSlice<f64>,
    eps2_cell: CudaSlice<f64>,
    eps2_one: CudaSlice<f64>,
    sums: CudaSlice<f64>,
    sums2: CudaSlice<f64>,
}

/// CUDA backend.
pub struct CudaBackend {
    ctx: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    _module: Arc<CudaModule>,
    f: Vec<CudaFunction>,
    d: Dev,
    vecs: Vec<CudaSlice<f64>>,
    nc: u32,
    nf: u32,
    ni: u32,
    ne: u32,
    chunk: u32,
    n_slices: u32,
    stored: u32,
    venkat_k: Option<f64>,
    host: Discretisation,
}

impl CudaBackend {
    /// Uploads the discretisation to GPU `ordinal` and loads the PTX module.
    pub fn new(disc: &Discretisation, ordinal: usize, ptx: &Path) -> anyhow::Result<Self> {
        anyhow::ensure!(
            ptx.exists(),
            "PTX module `{}` not found; build it with gpu/ptx-kernels/build.sh",
            ptx.display()
        );
        let ctx =
            CudaContext::new(ordinal).context("CUDA device (is the NVIDIA driver installed?)")?;
        let src = std::fs::read_to_string(ptx)?;
        let module = ctx
            .load_module(Ptx::from_src(src))
            .context("loading the PTX module")?;
        let mut f = Vec::new();
        for e in ENTRIES {
            f.push(
                module
                    .load_function(e)
                    .with_context(|| format!("entry point {e} (rebuild the PTX module)"))?,
            );
        }
        let stream = ctx.default_stream();
        let s = &stream;
        let soa = |v: &crate::disc::Soa3| -> anyhow::Result<Soa> {
            Ok(Soa {
                x: s.clone_htod(&v.x)?,
                y: s.clone_htod(&v.y)?,
                z: s.clone_htod(&v.z)?,
            })
        };
        let nc = disc.n_cells;
        let nb = nc.div_ceil(BLOCK);
        let d = Dev {
            owner: s.clone_htod(&disc.owner)?,
            neighbour: s.clone_htod(&nonempty_u32(&disc.neighbour))?,
            cell_ptr: s.clone_htod(&disc.cell_ptr)?,
            cell_faces: s.clone_htod(&disc.cell_faces)?,
            cc: soa(&disc.cell_centroid)?,
            fc: soa(&disc.face_centroid)?,
            k: soa(&disc.k)?,
            lsq: soa(&disc.lsq)?,
            dgeo: s.clone_htod(&disc.dgeo)?,
            weight: s.clone_htod(&disc.weight)?,
            bc_a: s.clone_htod(&nonempty_f64(&disc.bc_a))?,
            bc_b: s.clone_htod(&nonempty_f64(&disc.bc_b))?,
            mass_flux: s.clone_htod(&disc.mass_flux)?,
            gamma: s.clone_htod(&disc.gamma)?,
            source: s.clone_htod(&disc.source)?,
            slice_ptr: s.clone_htod(&disc.sell.slice_ptr)?,
            slice_len: s.clone_htod(&disc.sell.slice_len)?,
            col: s.clone_htod(&disc.sell.col)?,
            val: s.alloc_zeros(disc.sell.stored())?,
            eps2_cell: s.alloc_zeros(nc)?,
            eps2_one: s.alloc_zeros(1)?,
            sums: s.alloc_zeros(nb)?,
            sums2: s.alloc_zeros(nb)?,
        };
        let mut vecs = Vec::with_capacity(V::COUNT);
        for i in 0..V::COUNT {
            let n = if V::FACE.iter().any(|v| v.index() == i) {
                disc.n_faces
            } else {
                nc
            };
            vecs.push(s.alloc_zeros(n)?);
        }
        let mut be = Self {
            ctx,
            stream,
            _module: module,
            f,
            d,
            vecs,
            nc: nc as u32,
            nf: disc.n_faces as u32,
            ni: disc.n_interior as u32,
            ne: disc.cell_faces.len() as u32,
            chunk: disc.sell.chunk as u32,
            n_slices: disc.sell.slice_ptr.len() as u32,
            stored: disc.sell.stored() as u32,
            venkat_k: None,
            host: disc.clone(),
        };
        be.upload(V::Psi, &vec![1.0; nc])?;
        be.stream.synchronize()?;
        Ok(be)
    }

    fn func(&self, name: &str) -> &CudaFunction {
        let i = ENTRIES.iter().position(|e| *e == name).expect("entry");
        &self.f[i]
    }

    /// Two vectors at once, the second mutable (distinct).
    fn pair(&mut self, a: V, out: V) -> (&CudaSlice<f64>, &mut CudaSlice<f64>) {
        assert_ne!(a, out);
        let (i, j) = (a.index(), out.index());
        if i < j {
            let (l, r) = self.vecs.split_at_mut(j);
            (&l[i], &mut r[0])
        } else {
            let (l, r) = self.vecs.split_at_mut(i);
            (&r[0], &mut l[j])
        }
    }

    /// Takes a vector out of the table (put it back with `put`).
    fn take(&mut self, v: V) -> anyhow::Result<CudaSlice<f64>> {
        let empty = self.stream.alloc_zeros(1)?;
        Ok(std::mem::replace(&mut self.vecs[v.index()], empty))
    }

    fn put(&mut self, v: V, s: CudaSlice<f64>) {
        self.vecs[v.index()] = s;
    }

    fn block_dot(&mut self, a: V, b: V) -> anyhow::Result<f64> {
        let nb = (self.nc as usize).div_ceil(BLOCK);
        let nbu = nb as u32;
        let s = self.stream.clone();
        let f = self.func("blas_block_dot").clone();
        let (va, vb) = (&self.vecs[a.index()], &self.vecs[b.index()]);
        let mut l = s.launch_builder(&f);
        l.arg(&self.nc)
            .arg(&nbu)
            .arg(va)
            .arg(vb)
            .arg(&mut self.d.sums);
        // SAFETY: argument list of blas_block_dot; arrays of length n_cells
        // and n_blocks.
        unsafe { l.launch(cfg(nb)) }?;
        let sums = s.clone_dtoh(&self.d.sums)?;
        Ok(combine_block_sums(&sums))
    }
}

fn nonempty_u32(v: &[u32]) -> Vec<u32> {
    if v.is_empty() { vec![0] } else { v.to_vec() }
}

fn nonempty_f64(v: &[f64]) -> Vec<f64> {
    if v.is_empty() { vec![0.0] } else { v.to_vec() }
}

impl Backend for CudaBackend {
    fn name(&self) -> String {
        "cuda".into()
    }

    fn device(&self) -> String {
        let name = self.ctx.name().unwrap_or_else(|_| "NVIDIA GPU".into());
        let (ma, mi) = self.ctx.compute_capability().unwrap_or((0, 0));
        format!("{name} (sm_{ma}{mi})")
    }

    fn n_cells(&self) -> usize {
        self.nc as usize
    }

    fn upload(&mut self, v: V, data: &[f64]) -> anyhow::Result<()> {
        let s = self.stream.clone();
        let t = &mut self.vecs[v.index()];
        anyhow::ensure!(t.len() == data.len(), "upload of {v:?}: wrong length");
        s.memcpy_htod(data, t)?;
        Ok(())
    }

    fn download(&mut self, v: V) -> anyhow::Result<Vec<f64>> {
        Ok(self.stream.clone_dtoh(&self.vecs[v.index()])?)
    }

    fn sync(&mut self) -> anyhow::Result<()> {
        Ok(self.stream.synchronize()?)
    }

    fn gradient(&mut self) -> anyhow::Result<()> {
        let mut gx = self.take(V::Gx)?;
        let mut gy = self.take(V::Gy)?;
        let mut gz = self.take(V::Gz)?;
        {
            let d = &self.d;
            let f = self.func("fv_gradient");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc).arg(&self.nf).arg(&self.ni).arg(&self.ne);
            l.arg(&d.owner)
                .arg(&d.neighbour)
                .arg(&d.cell_ptr)
                .arg(&d.cell_faces);
            l.arg(&d.lsq.x).arg(&d.lsq.y).arg(&d.lsq.z);
            l.arg(&self.vecs[V::Phi.index()]).arg(&d.bc_a).arg(&d.bc_b);
            l.arg(&mut gx).arg(&mut gy).arg(&mut gz);
            // SAFETY: argument list of fv_gradient (gpu/ptx-kernels/src/fv.rs).
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(V::Gx, gx);
        self.put(V::Gy, gy);
        self.put(V::Gz, gz);
        Ok(())
    }

    fn limiter(&mut self, lim: Limiter) -> anyhow::Result<()> {
        let n_eps2: u32 = match lim {
            Limiter::None => {
                let ones = vec![1.0; self.nc as usize];
                return self.upload(V::Psi, &ones);
            }
            Limiter::Venkatakrishnan { k } => {
                if self.venkat_k != Some(k) {
                    let e = self.host.venkat_eps2(k);
                    self.stream.memcpy_htod(&e, &mut self.d.eps2_cell)?;
                    self.venkat_k = Some(k);
                }
                self.nc
            }
            Limiter::VenkatakrishnanWang { k } => {
                let (mn, mx) = self.min_max(V::Phi)?;
                let e = k * (mx - mn);
                self.stream.memcpy_htod(&[e * e], &mut self.d.eps2_one)?;
                1
            }
        };
        let mut psi = self.take(V::Psi)?;
        {
            let d = &self.d;
            let eps2 = if n_eps2 == 1 {
                &d.eps2_one
            } else {
                &d.eps2_cell
            };
            let f = self.func("fv_limiter");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc)
                .arg(&self.nf)
                .arg(&self.ni)
                .arg(&self.ne)
                .arg(&n_eps2);
            l.arg(&d.owner)
                .arg(&d.neighbour)
                .arg(&d.cell_ptr)
                .arg(&d.cell_faces);
            l.arg(&self.vecs[V::Phi.index()]).arg(&d.bc_a).arg(&d.bc_b);
            l.arg(&self.vecs[V::Gx.index()])
                .arg(&self.vecs[V::Gy.index()])
                .arg(&self.vecs[V::Gz.index()]);
            l.arg(&d.cc.x).arg(&d.cc.y).arg(&d.cc.z);
            l.arg(&d.fc.x).arg(&d.fc.y).arg(&d.fc.z);
            l.arg(eps2).arg(&mut psi);
            // SAFETY: argument list of fv_limiter.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(V::Psi, psi);
        Ok(())
    }

    fn faces(&mut self, st: FaceSettings) -> anyhow::Result<()> {
        let mut cp = self.take(V::Cp)?;
        let mut cn = self.take(V::Cn)?;
        let mut bf = self.take(V::Bf)?;
        let nonorth = if st.nonorth { 1.0f64 } else { 0.0 };
        {
            let d = &self.d;
            let f = self.func("fv_faces");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc).arg(&self.nf).arg(&self.ni);
            l.arg(&st.blend).arg(&nonorth).arg(&st.nonorth_limit);
            l.arg(&d.owner).arg(&d.neighbour).arg(&d.bc_a).arg(&d.bc_b);
            l.arg(&d.mass_flux).arg(&d.gamma).arg(&d.dgeo);
            l.arg(&d.k.x).arg(&d.k.y).arg(&d.k.z).arg(&d.weight);
            l.arg(&d.fc.x).arg(&d.fc.y).arg(&d.fc.z);
            l.arg(&d.cc.x).arg(&d.cc.y).arg(&d.cc.z);
            l.arg(&self.vecs[V::Phi.index()]);
            l.arg(&self.vecs[V::Gx.index()])
                .arg(&self.vecs[V::Gy.index()])
                .arg(&self.vecs[V::Gz.index()]);
            l.arg(&self.vecs[V::Psi.index()]);
            l.arg(&mut cp).arg(&mut cn).arg(&mut bf);
            // SAFETY: argument list of fv_faces.
            unsafe { l.launch(cfg(self.nf as usize)) }?;
        }
        self.put(V::Cp, cp);
        self.put(V::Cn, cn);
        self.put(V::Bf, bf);
        Ok(())
    }

    fn assemble(&mut self, alpha: f64) -> anyhow::Result<()> {
        let mut rhs = self.take(V::Rhs)?;
        let mut dinv = self.take(V::Dinv)?;
        {
            let d = &mut self.d;
            let i = ENTRIES.iter().position(|e| *e == "fv_assemble").unwrap();
            let mut l = self.stream.launch_builder(&self.f[i]);
            l.arg(&self.nc).arg(&self.nf).arg(&self.ni).arg(&self.ne);
            l.arg(&self.chunk).arg(&self.n_slices).arg(&alpha);
            l.arg(&d.owner)
                .arg(&d.neighbour)
                .arg(&d.cell_ptr)
                .arg(&d.cell_faces);
            l.arg(&self.vecs[V::Cp.index()])
                .arg(&self.vecs[V::Cn.index()])
                .arg(&self.vecs[V::Bf.index()]);
            l.arg(&d.source).arg(&self.vecs[V::Phi.index()]);
            l.arg(&d.slice_ptr)
                .arg(&mut d.val)
                .arg(&mut rhs)
                .arg(&mut dinv);
            // SAFETY: argument list of fv_assemble.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(V::Rhs, rhs);
        self.put(V::Dinv, dinv);
        Ok(())
    }

    fn spmv(&mut self, x: V, out: V) -> anyhow::Result<()> {
        let mut y = self.take(out)?;
        {
            let d = &self.d;
            let f = self.func("g1_spmv_sell");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc)
                .arg(&self.chunk)
                .arg(&self.n_slices)
                .arg(&self.stored);
            l.arg(&d.slice_ptr)
                .arg(&d.slice_len)
                .arg(&d.col)
                .arg(&d.val);
            l.arg(&self.vecs[x.index()]).arg(&mut y);
            // SAFETY: argument list of g1_spmv_sell; σ = 1, so the sorted
            // rows are the cell rows.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(out, y);
        Ok(())
    }

    fn residual(&mut self, x: V, b: V, out: V) -> anyhow::Result<()> {
        anyhow::ensure!(out != x && out != b, "residual: output aliases an input");
        let mut y = self.take(out)?;
        {
            let d = &self.d;
            let f = self.func("fv_residual");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc)
                .arg(&self.chunk)
                .arg(&self.n_slices)
                .arg(&self.stored);
            l.arg(&d.slice_ptr)
                .arg(&d.slice_len)
                .arg(&d.col)
                .arg(&d.val);
            l.arg(&self.vecs[x.index()])
                .arg(&self.vecs[b.index()])
                .arg(&mut y);
            // SAFETY: argument list of fv_residual.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(out, y);
        Ok(())
    }

    fn dot(&mut self, a: V, b: V) -> anyhow::Result<f64> {
        self.block_dot(a, b)
    }

    fn min_max(&mut self, a: V) -> anyhow::Result<(f64, f64)> {
        let nb = (self.nc as usize).div_ceil(BLOCK);
        let nbu = nb as u32;
        {
            let d = &mut self.d;
            let i = ENTRIES
                .iter()
                .position(|e| *e == "blas_block_min_max")
                .unwrap();
            let mut l = self.stream.launch_builder(&self.f[i]);
            l.arg(&self.nc)
                .arg(&nbu)
                .arg(&self.vecs[a.index()])
                .arg(&mut d.sums)
                .arg(&mut d.sums2);
            // SAFETY: argument list of blas_block_min_max.
            unsafe { l.launch(cfg(nb)) }?;
        }
        let mins = self.stream.clone_dtoh(&self.d.sums)?;
        let maxs = self.stream.clone_dtoh(&self.d.sums2)?;
        Ok((
            mins.iter().copied().fold(f64::INFINITY, f64::min),
            maxs.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        ))
    }

    fn update_p(&mut self, beta: f64, omega: f64) -> anyhow::Result<()> {
        let mut p = self.take(V::P)?;
        {
            let f = self.func("blas_update_p");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc).arg(&beta).arg(&omega);
            l.arg(&self.vecs[V::R.index()])
                .arg(&self.vecs[V::Vv.index()])
                .arg(&mut p);
            // SAFETY: argument list of blas_update_p.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(V::P, p);
        Ok(())
    }

    fn sub_scaled(&mut self, x: V, z: V, alpha: f64, out: V) -> anyhow::Result<()> {
        anyhow::ensure!(out != x && out != z, "sub_scaled: output aliases an input");
        let mut y = self.take(out)?;
        {
            let f = self.func("blas_sub_scaled");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc).arg(&alpha);
            l.arg(&self.vecs[x.index()])
                .arg(&self.vecs[z.index()])
                .arg(&mut y);
            // SAFETY: argument list of blas_sub_scaled.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(out, y);
        Ok(())
    }

    fn update_x(&mut self, alpha: f64, omega: f64) -> anyhow::Result<()> {
        let mut x = self.take(V::Phi)?;
        {
            let f = self.func("blas_update_x");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc).arg(&alpha).arg(&omega);
            l.arg(&self.vecs[V::Ph.index()])
                .arg(&self.vecs[V::Sh.index()])
                .arg(&mut x);
            // SAFETY: argument list of blas_update_x.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(V::Phi, x);
        Ok(())
    }

    fn precondition(&mut self, x: V, out: V) -> anyhow::Result<()> {
        anyhow::ensure!(
            out != x && out != V::Dinv,
            "precondition: output aliases an input"
        );
        let mut y = self.take(out)?;
        {
            let f = self.func("blas_scale");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc);
            l.arg(&self.vecs[V::Dinv.index()])
                .arg(&self.vecs[x.index()])
                .arg(&mut y);
            // SAFETY: argument list of blas_scale.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(out, y);
        Ok(())
    }

    fn copy(&mut self, src: V, dst: V) -> anyhow::Result<()> {
        let s = self.stream.clone();
        let (a, b) = self.pair(src, dst);
        s.memcpy_dtod(a, b)?;
        Ok(())
    }
}

//! NVIDIA GPU backends (feature `cuda`).
//!
//! - `cuda-reference`: the hand-written CUDA C kernels of
//!   `kernels/g1_reference.cu`, compiled at run time by NVRTC for the
//!   installed GPU;
//! - `cuda-rust`: the single-source kernel bodies of `inncraft-kernels`,
//!   compiled by rustc for `nvptx64-nvidia-cuda` (`gpu/ptx-kernels`) and
//!   loaded from the PTX file.
//!
//! Both modules export the same entry points with the same arguments, so one
//! launch path serves both and the comparison is only the kernel code. The
//! CUDA driver and NVRTC are loaded dynamically: a build needs no CUDA
//! toolkit, a run needs the NVIDIA driver (and NVRTC for the reference).

use super::Backend;
use crate::workload::{Kernel, Workload};
use anyhow::Context as _;
use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, LaunchConfig, PushKernelArg,
};
use cudarc::nvrtc::{CompileOptions, Ptx};
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

const REFERENCE_SRC: &str = include_str!("../../kernels/g1_reference.cu");

/// Threads per block. 256 is a good default on every NVIDIA architecture
/// since Maxwell for memory-bound one-thread-per-output kernels.
const BLOCK: u32 = 256;

struct Device {
    owner: CudaSlice<u32>,
    neighbour: CudaSlice<u32>,
    phi: CudaSlice<f64>,
    gamma: CudaSlice<f64>,
    mass_flux: CudaSlice<f64>,
    cell_ptr: CudaSlice<u32>,
    cell_faces: CudaSlice<u32>,
    flux: CudaSlice<f64>,
    row_ptr: CudaSlice<u32>,
    col: CudaSlice<u32>,
    val: CudaSlice<f64>,
    s_ptr: CudaSlice<u32>,
    s_len: CudaSlice<u32>,
    s_col: CudaSlice<u32>,
    s_val: CudaSlice<f64>,
    x: CudaSlice<f64>,
    out: CudaSlice<f64>,
    n_cells: u32,
    n_faces: u32,
    nnz: u32,
    chunk: u32,
    n_slices: u32,
    stored: u32,
}

/// CUDA backend.
pub struct CudaBackend {
    name: &'static str,
    ctx: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    functions: HashMap<Kernel, CudaFunction>,
    dev: Option<Device>,
}

impl CudaBackend {
    /// The hand-written CUDA C reference, compiled by NVRTC.
    pub fn reference(ordinal: usize) -> anyhow::Result<Self> {
        let ctx =
            CudaContext::new(ordinal).context("CUDA device (is the NVIDIA driver installed?)")?;
        let (major, minor) = ctx.compute_capability()?;
        let arch: &'static str = Box::leak(format!("compute_{major}{minor}").into_boxed_str());
        let ptx = cudarc::nvrtc::compile_ptx_with_opts(
            REFERENCE_SRC,
            CompileOptions {
                arch: Some(arch),
                // no a*b+c contraction: bitwise comparable with the CPU and
                // the Rust PTX (rustc does not contract by default)
                fmad: Some(false),
                ..Default::default()
            },
        )
        .context("NVRTC compilation of g1_reference.cu")?;
        Self::with_module("cuda-reference", ctx, ptx)
    }

    /// The single-source Rust kernels, from the PTX file built by
    /// `gpu/ptx-kernels/build.sh`.
    pub fn rust(ordinal: usize, ptx_path: &Path) -> anyhow::Result<Self> {
        anyhow::ensure!(
            ptx_path.exists(),
            "PTX module {} not found; build it with gpu/ptx-kernels/build.sh",
            ptx_path.display()
        );
        let ctx =
            CudaContext::new(ordinal).context("CUDA device (is the NVIDIA driver installed?)")?;
        let src = std::fs::read_to_string(ptx_path)?;
        Self::with_module("cuda-rust", ctx, Ptx::from_src(src))
    }

    fn with_module(name: &'static str, ctx: Arc<CudaContext>, ptx: Ptx) -> anyhow::Result<Self> {
        let module = ctx.load_module(ptx).context("loading the PTX module")?;
        let mut functions = HashMap::new();
        for k in Kernel::ALL {
            functions.insert(
                k,
                module.load_function(k.entry()).with_context(|| k.entry())?,
            );
        }
        let stream = ctx.default_stream();
        Ok(Self {
            name,
            ctx,
            stream,
            functions,
            dev: None,
        })
    }
}

impl Backend for CudaBackend {
    fn name(&self) -> &'static str {
        self.name
    }

    fn device(&self) -> String {
        let name = self.ctx.name().unwrap_or_else(|_| "NVIDIA GPU".into());
        let (ma, mi) = self.ctx.compute_capability().unwrap_or((0, 0));
        let mem = self.ctx.total_mem().unwrap_or(0) as f64 / (1u64 << 30) as f64;
        format!("{name} (sm_{ma}{mi}, {mem:.1} GiB)")
    }

    fn setup(&mut self, w: &Workload) -> anyhow::Result<()> {
        let s = &self.stream;
        let m = &w.mesh;
        let max_out = m
            .n_faces()
            .max(m.n_cells)
            .max(w.sell.slice_ptr.len() * w.sell.chunk);
        self.dev = Some(Device {
            owner: s.clone_htod(&m.owner)?,
            neighbour: s.clone_htod(&m.neighbour)?,
            phi: s.clone_htod(&w.phi)?,
            gamma: s.clone_htod(&w.gamma)?,
            mass_flux: s.clone_htod(&w.mass_flux)?,
            cell_ptr: s.clone_htod(&m.cell_ptr)?,
            cell_faces: s.clone_htod(&m.cell_faces)?,
            flux: s.clone_htod(&w.flux)?,
            row_ptr: s.clone_htod(&w.csr.row_ptr)?,
            col: s.clone_htod(&w.csr.col)?,
            val: s.clone_htod(&w.csr.val)?,
            s_ptr: s.clone_htod(&w.sell.slice_ptr)?,
            s_len: s.clone_htod(&w.sell.slice_len)?,
            s_col: s.clone_htod(&w.sell.col)?,
            s_val: s.clone_htod(&w.sell.val)?,
            x: s.clone_htod(&w.x)?,
            out: s.alloc_zeros(max_out)?,
            n_cells: m.n_cells as u32,
            n_faces: m.n_faces() as u32,
            nnz: w.csr.nnz() as u32,
            chunk: w.sell.chunk as u32,
            n_slices: w.sell.slice_ptr.len() as u32,
            stored: w.sell.stored() as u32,
        });
        s.synchronize()?;
        Ok(())
    }

    fn run(&mut self, _w: &Workload, k: Kernel) -> anyhow::Result<f64> {
        let d = self.dev.as_mut().context("setup before run")?;
        let f = &self.functions[&k];
        let n = match k {
            Kernel::FaceFlux => d.n_faces,
            _ => d.n_cells,
        };
        let cfg = LaunchConfig {
            grid_dim: (n.div_ceil(BLOCK), 1, 1),
            block_dim: (BLOCK, 1, 1),
            shared_mem_bytes: 0,
        };
        let start = self
            .ctx
            .new_event(Some(cudarc::driver::sys::CUevent_flags::CU_EVENT_DEFAULT))?;
        let stop = self
            .ctx
            .new_event(Some(cudarc::driver::sys::CUevent_flags::CU_EVENT_DEFAULT))?;
        start.record(&self.stream)?;
        let mut b = self.stream.launch_builder(f);
        let n_entries = 2 * d.n_faces;
        match k {
            Kernel::FaceFlux => {
                b.arg(&d.n_faces)
                    .arg(&d.n_cells)
                    .arg(&d.owner)
                    .arg(&d.neighbour);
                b.arg(&d.phi)
                    .arg(&d.gamma)
                    .arg(&d.mass_flux)
                    .arg(&mut d.out);
            }
            Kernel::CellGather => {
                b.arg(&d.n_cells).arg(&n_entries).arg(&d.n_faces);
                b.arg(&d.cell_ptr)
                    .arg(&d.cell_faces)
                    .arg(&d.flux)
                    .arg(&mut d.out);
            }
            Kernel::SpmvCsr => {
                b.arg(&d.n_cells).arg(&d.nnz).arg(&d.row_ptr).arg(&d.col);
                b.arg(&d.val).arg(&d.x).arg(&mut d.out);
            }
            Kernel::SpmvSell => {
                b.arg(&d.n_cells)
                    .arg(&d.chunk)
                    .arg(&d.n_slices)
                    .arg(&d.stored);
                b.arg(&d.s_ptr).arg(&d.s_len).arg(&d.s_col).arg(&d.s_val);
                b.arg(&d.x).arg(&mut d.out);
            }
        }
        // SAFETY: argument lists match the kernel signatures in
        // g1_reference.cu and gpu/ptx-kernels; every array holds the sizes
        // passed with it.
        unsafe { b.launch(cfg) }?;
        stop.record(&self.stream)?;
        stop.synchronize()?;
        Ok(start.elapsed_ms(&stop)? as f64 * 1e-3)
    }

    fn output(&mut self, w: &Workload, k: Kernel) -> anyhow::Result<Vec<f64>> {
        let d = self.dev.as_ref().context("setup before run")?;
        let all = self.stream.clone_dtoh(&d.out)?;
        let n = w.output_len(k);
        if k == Kernel::SpmvSell {
            let mut y = vec![0.0; n];
            for r in 0..n {
                y[w.sell.perm[r] as usize] = all[r];
            }
            return Ok(y);
        }
        Ok(all[..n].to_vec())
    }
}

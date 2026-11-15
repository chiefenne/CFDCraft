//! The acceptance checks of stage 2 as functions, used by the tests, the
//! `ennecraft scalar` verification cases and the gallery:
//!
//! - [`gradient_linear`]: the least-squares gradient of a linear field is
//!   exact on every mesh (D24);
//! - [`manufactured`]: convection–diffusion with a manufactured solution,
//!   error norms against the exact cell-centroid values;
//! - [`order_sweep`]: the observed order of accuracy over nested refinements;
//! - [`oblique_step`]: convection of a discontinuity at an angle to the mesh,
//!   for the bounds of the limited schemes (D25).

use crate::backend::{self, Backend, V};
use crate::disc::Discretisation;
use crate::problem::{Boundary, BoundaryKind, Problem, Profile};
use crate::solver::{self, Outcome, Settings, Timings};
use anyhow::Result;
use ennecraft_mesh::generate::{BOX_PATCHES, BoxSpec, Kind};
use ennecraft_mesh::renumber::{Ordering, renumber};
use ennecraft_mesh::{Geometry, Mesh};
use serde::Serialize;

/// A generated test mesh with RCM cell order, and its geometry.
pub fn test_mesh(
    kind: Kind,
    dimension: u8,
    cells: usize,
    refine: u32,
    perturb: f64,
) -> Result<(Mesh, Geometry)> {
    let mut spec = BoxSpec::new(kind, dimension, [cells, cells, cells]);
    spec.refine = refine;
    spec.perturb = perturb;
    let mut mesh = spec.generate().map_err(anyhow::Error::msg)?;
    renumber(&mut mesh, Ordering::Rcm, 0).map_err(anyhow::Error::msg)?;
    let geo = Geometry::compute(&mesh);
    Ok((mesh, geo))
}

/// Largest relative error of the least-squares gradient of the linear field
/// φ = 1 + g·x over all cells. Patches in `gradient_patches` take the
/// fixed normal gradient of the field, the others its value; with
/// `dimension = 2` the z patches are zero-gradient (g_z = 0).
pub fn gradient_linear(
    mesh: &Mesh,
    geo: &Geometry,
    dimension: u8,
    gradient_patches: &[&str],
) -> Result<f64> {
    let g = if dimension == 2 {
        [0.7, -1.3, 0.0]
    } else {
        [0.7, -1.3, 2.1]
    };
    let mut p = Problem::default();
    for name in BOX_PATCHES {
        let kind = if dimension == 2 && (name == "zmin" || name == "zmax") {
            BoundaryKind::ZeroGradient
        } else if gradient_patches.contains(&name) {
            BoundaryKind::FixedGradient
        } else {
            BoundaryKind::FixedValue
        };
        p.boundary.insert(
            name.into(),
            Boundary {
                kind,
                profile: Profile::Linear {
                    value: 1.0,
                    gradient: g,
                },
                exact: None,
            },
        );
    }
    let d = Discretisation::new(mesh, geo, &p)?;
    let phi: Vec<f64> = (0..mesh.n_cells)
        .map(|c| {
            let x = geo.cell_centroid[c];
            1.0 + g[0] * x[0] + g[1] * x[1] + g[2] * x[2]
        })
        .collect();
    let mut be = backend::cpu::CpuBackend::new(d, None)?;
    be.upload(V::Phi, &phi)?;
    be.gradient()?;
    let (gx, gy, gz) = (
        be.download(V::Gx)?,
        be.download(V::Gy)?,
        be.download(V::Gz)?,
    );
    let gn = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    let mut worst = 0.0f64;
    for c in 0..mesh.n_cells {
        let e = ((gx[c] - g[0]).powi(2) + (gy[c] - g[1]).powi(2) + (gz[c] - g[2]).powi(2)).sqrt();
        worst = worst.max(e / gn);
    }
    Ok(worst)
}

/// Result of one manufactured-solution run.
#[derive(Debug, Clone, Serialize)]
pub struct MmsResult {
    /// Mesh type.
    pub kind: String,
    /// Cells.
    pub cells: usize,
    /// Mesh size h = (V / N)^(1/3) (2D: (A / N)^(1/2)).
    pub h: f64,
    /// Volume-weighted L2 norm of the error.
    pub l2: f64,
    /// Largest error.
    pub linf: f64,
    /// Largest non-orthogonality (degrees).
    pub max_nonorth_deg: f64,
    /// Outer iterations.
    pub outer: usize,
    /// Converged.
    pub converged: bool,
    /// Phase timings.
    pub timings: Timings,
    /// ns per cell and outer iteration.
    pub ns_per_cell_iteration: f64,
}

/// The manufactured-solution problem on a box mesh: fixed values on all
/// patches (in 2D zero gradient on the z patches and a planar solution).
///
/// Γ = 0.2 kg/(m s) with |u| ≈ 1.2 m/s gives a Péclet number of about 6 on
/// the unit box, so convection and diffusion both matter. At smaller Γ the
/// convective and diffusive truncation errors have opposite signs and
/// partly cancel on the coarse meshes, which spoils the observed order
/// until much finer meshes (docs/verification/stage2.md).
pub fn mms_problem(dimension: u8) -> Problem {
    let mut p = if dimension == 2 {
        Problem::manufactured(&["xmin", "xmax", "ymin", "ymax"], &["zmin", "zmax"], true)
    } else {
        Problem::manufactured(&BOX_PATCHES, &[], false)
    };
    p.diffusivity = MMS_DIFFUSIVITY;
    p
}

/// Diffusion coefficient of [`mms_problem`].
pub const MMS_DIFFUSIVITY: f64 = 0.2;

/// Solves `problem` (with a manufactured solution) on `mesh` with backend
/// `be_name` and compares with the exact solution. Returns the result and
/// the final cell values.
#[allow(clippy::too_many_arguments)]
pub fn manufactured(
    kind: &str,
    mesh: &Mesh,
    geo: &Geometry,
    dimension: u8,
    problem: &Problem,
    settings: &Settings,
    be_name: &str,
    threads: Option<usize>,
    ptx: &std::path::Path,
) -> Result<(MmsResult, Vec<f64>)> {
    let m = problem
        .manufactured
        .ok_or_else(|| anyhow::anyhow!("the problem has no manufactured solution"))?;
    let d = Discretisation::new(mesh, geo, problem)?;
    let max_nonorth_deg = d.max_nonorth_deg;
    let mut be = backend::create(be_name, &d, threads, ptx)?;
    let phi0 = vec![m.offset; mesh.n_cells];
    be.upload(V::Phi, &phi0)?;
    let out: Outcome = solver::solve(be.as_mut(), settings, |_| {})?;
    let phi = be.download(V::Phi)?;
    let (mut e2, mut einf, mut vol) = (0.0, 0.0f64, 0.0);
    for (c, &v) in phi.iter().enumerate() {
        let e = v - m.phi(geo.cell_centroid[c]);
        e2 += e * e * geo.cell_volume[c];
        einf = einf.max(e.abs());
        vol += geo.cell_volume[c];
    }
    let n = mesh.n_cells as f64;
    let h = if dimension == 2 {
        let (lo, hi) = mesh.bounding_box();
        (vol / (hi[2] - lo[2]) / n).sqrt()
    } else {
        (vol / n).cbrt()
    };
    let r = MmsResult {
        kind: kind.into(),
        cells: mesh.n_cells,
        h,
        l2: (e2 / vol).sqrt(),
        linf: einf,
        max_nonorth_deg,
        outer: out.outer,
        converged: out.converged,
        ns_per_cell_iteration: out.ns_per_cell_iteration(mesh.n_cells),
        timings: out.timings,
    };
    Ok((r, phi))
}

/// Observed order between successive results, from the L2 errors.
pub fn observed_orders(r: &[MmsResult], dimension: u8) -> Vec<f64> {
    let _ = dimension;
    r.windows(2)
        .map(|w| (w[0].l2 / w[1].l2).ln() / (w[0].h / w[1].h).ln())
        .collect()
}

/// Manufactured-solution runs on nested refinements `levels` of a base mesh
/// of `cells` intervals per direction.
#[allow(clippy::too_many_arguments)]
pub fn order_sweep(
    kind: Kind,
    dimension: u8,
    cells: usize,
    levels: &[u32],
    perturb: f64,
    settings: &Settings,
    be_name: &str,
    ptx: &std::path::Path,
) -> Result<Vec<MmsResult>> {
    let p = mms_problem(dimension);
    let mut out = Vec::new();
    for &l in levels {
        let (mesh, geo) = test_mesh(kind, dimension, cells, l, perturb)?;
        let (r, _) = manufactured(
            kind.name(),
            &mesh,
            &geo,
            dimension,
            &p,
            settings,
            be_name,
            None,
            ptx,
        )?;
        out.push(r);
    }
    Ok(out)
}

/// The oblique-step problem (a classic test of convection schemes, e.g.
/// Leonard 1979): pure convection at 30° to the x axis across the unit
/// square, φ = 1 entering through `xmin` above y = 0.25... in fact through
/// the lower part of `xmin` and 0 through `ymin`, zero gradient at the
/// outflow and on the z patches; nearly no diffusion. The exact solution is
/// a discontinuity along the characteristic, so φ must stay in [0, 1].
pub fn oblique_step_problem(diffusivity: f64) -> Problem {
    let a = 30f64.to_radians();
    let mut p = Problem {
        diffusivity,
        velocity: [a.cos(), a.sin(), 0.0],
        ..Default::default()
    };
    p.boundary.insert(
        "xmin".into(),
        Boundary {
            kind: BoundaryKind::FixedValue,
            profile: Profile::Step {
                axis: 1,
                at: 0.2,
                below: 1.0,
                above: 0.0,
            },
            exact: None,
        },
    );
    p.boundary.insert("ymin".into(), Boundary::value(1.0));
    for name in ["xmax", "ymax", "zmin", "zmax"] {
        p.boundary.insert(name.into(), Boundary::zero_gradient());
    }
    p
}

/// Bounds of the oblique-step solution: (min, max, outcome).
pub fn oblique_step(
    mesh: &Mesh,
    geo: &Geometry,
    settings: &Settings,
    be_name: &str,
    ptx: &std::path::Path,
) -> Result<(f64, f64, Outcome, Vec<f64>)> {
    let p = oblique_step_problem(1e-6);
    let d = Discretisation::new(mesh, geo, &p)?;
    let mut be: Box<dyn Backend> = backend::create(be_name, &d, None, ptx)?;
    be.upload(V::Phi, &vec![0.0; mesh.n_cells])?;
    let out = solver::solve(be.as_mut(), settings, |_| {})?;
    let phi = be.download(V::Phi)?;
    let mn = phi.iter().copied().fold(f64::INFINITY, f64::min);
    let mx = phi.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    Ok((mn, mx, out, phi))
}

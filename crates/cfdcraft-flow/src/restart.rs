//! Restart files (D42): the solution and the solver state of a flow run in
//! one HDF5 file, so that a run continued from it gives bitwise the result
//! of the run without the interruption (same mesh, same decomposition and
//! device).
//!
//! Layout (format `cfdcraft-restart`, version 1):
//!
//! - attributes `format`, `version`, `case`, `iteration`, `amg_setups`,
//!   `transient`, and for a transient run `time`, `step`, `dt_prev`,
//!   `levels`, `dt_next`, `scheme`;
//! - per cell (local order, with the global identities in `cell_id`):
//!   `u_x u_y u_z p`, the gradients `grad_u_<i><j>` and `grad_p_<j>`, the
//!   limiters `psi_<i>`, and for a transient run the old levels `u_old_<i>`
//!   and `g_old_<k>` (velocity gradients and limiters of the level
//!   before, for the extrapolation);
//! - per face (with `face_id`): `mass_flux`, and transient `m_old`, `x_cur`,
//!   `x_old`; per boundary face `ub_<i>`, `pb`;
//! - `amg_matrix`: the pressure matrix of the last multigrid setup (the
//!   hierarchy is rebuilt from it, so the preconditioner is the same).
//!
//! A file whose identities are in another order (another renumbering) is
//! mapped by identity; the multigrid matrix is then dropped and the
//! hierarchy rebuilt, so the continuation is no longer bit-identical.

use crate::simplec::State;
use crate::transient::{Clock, OldLevels, TimeScheme};
use anyhow::Result;
#[cfg(feature = "hdf5")]
use anyhow::{Context, bail, ensure};
use cfdcraft_fv::disc::Soa3;
use std::path::Path;

/// The content of a restart file.
#[derive(Debug, Clone)]
pub struct Restart {
    /// Case name.
    pub case: String,
    /// Global cell identities in local order.
    pub cell_id: Vec<u64>,
    /// Global face identities in local order.
    pub face_id: Vec<u64>,
    /// The fields.
    pub state: State,
    /// Outer iterations done.
    pub iteration: usize,
    /// Multigrid setups done.
    pub amg_setups: usize,
    /// Pressure matrix of the last multigrid setup.
    pub amg_matrix: Option<Vec<f64>>,
    /// Transient state: old levels, clock, next time step, scheme.
    pub transient: Option<(OldLevels, Clock, f64, TimeScheme)>,
}

#[cfg(feature = "hdf5")]
mod io {
    use super::*;
    use hdf5_metno as h5;

    fn ds(g: &h5::File, name: &str, data: &[f64]) -> Result<()> {
        g.new_dataset::<f64>()
            .shape([data.len()])
            .create(name)?
            .write_raw(data)?;
        Ok(())
    }

    fn ds_u64(g: &h5::File, name: &str, data: &[u64]) -> Result<()> {
        g.new_dataset::<u64>()
            .shape([data.len()])
            .create(name)?
            .write_raw(data)?;
        Ok(())
    }

    fn attr<T: h5::H5Type>(g: &h5::File, name: &str, v: &T) -> Result<()> {
        g.new_attr::<T>().create(name)?.write_scalar(v)?;
        Ok(())
    }

    fn sattr(g: &h5::File, name: &str, v: &str) -> Result<()> {
        let s: h5::types::VarLenUnicode = v.parse().map_err(|e| anyhow::anyhow!("{e:?}"))?;
        g.new_attr::<h5::types::VarLenUnicode>()
            .create(name)?
            .write_scalar(&s)?;
        Ok(())
    }

    fn rd(g: &h5::File, name: &str) -> Result<Vec<f64>> {
        g.dataset(name)
            .with_context(|| format!("restart file: dataset `{name}`"))?
            .read_raw::<f64>()
            .map_err(Into::into)
    }

    fn rattr<T: h5::H5Type>(g: &h5::File, name: &str) -> Result<T> {
        g.attr(name)
            .with_context(|| format!("restart file: attribute `{name}`"))?
            .read_scalar::<T>()
            .map_err(Into::into)
    }

    fn rsattr(g: &h5::File, name: &str) -> Result<String> {
        Ok(rattr::<h5::types::VarLenUnicode>(g, name)?
            .as_str()
            .to_owned())
    }

    fn put3(g: &h5::File, name: &str, v: &Soa3) -> Result<()> {
        ds(g, &format!("{name}x"), &v.x)?;
        ds(g, &format!("{name}y"), &v.y)?;
        ds(g, &format!("{name}z"), &v.z)
    }

    fn get3(g: &h5::File, name: &str) -> Result<Soa3> {
        Ok(Soa3 {
            x: rd(g, &format!("{name}x"))?,
            y: rd(g, &format!("{name}y"))?,
            z: rd(g, &format!("{name}z"))?,
        })
    }

    /// Writes `r` to `path` (through a temporary file, so a crash never
    /// leaves a broken restart behind).
    pub fn write(path: &Path, r: &Restart) -> Result<()> {
        let tmp = path.with_extension("h5.tmp");
        {
            let g =
                h5::File::create(&tmp).with_context(|| format!("creating {}", tmp.display()))?;
            sattr(&g, "format", "cfdcraft-restart")?;
            attr(&g, "version", &1i64)?;
            sattr(&g, "case", &r.case)?;
            attr(&g, "iteration", &(r.iteration as u64))?;
            attr(&g, "amg_setups", &(r.amg_setups as u64))?;
            attr(&g, "transient", &u8::from(r.transient.is_some()))?;
            ds_u64(&g, "cell_id", &r.cell_id)?;
            ds_u64(&g, "face_id", &r.face_id)?;
            let s = &r.state;
            put3(&g, "u_", &s.u)?;
            ds(&g, "p", &s.p)?;
            ds(&g, "mass_flux", &s.mass_flux)?;
            put3(&g, "ub_", &s.ub)?;
            ds(&g, "pb", &s.pb)?;
            for (i, c) in ["x", "y", "z"].iter().enumerate() {
                put3(&g, &format!("grad_u_{c}"), &s.grad_u[i])?;
            }
            put3(&g, "grad_p_", &s.grad_p)?;
            put3(&g, "psi_", &s.psi)?;
            if let Some(m) = &r.amg_matrix {
                ds(&g, "amg_matrix", m)?;
            }
            if let Some((o, clock, dt_next, scheme)) = &r.transient {
                attr(&g, "time", &clock.time)?;
                attr(&g, "step", &(clock.step as u64))?;
                attr(&g, "dt_prev", &clock.dt_prev)?;
                attr(&g, "levels", &(clock.levels as u64))?;
                attr(&g, "dt_next", dt_next)?;
                sattr(&g, "scheme", scheme_name(*scheme))?;
                ds(&g, "u_old_x", &o.u_old[0])?;
                ds(&g, "u_old_y", &o.u_old[1])?;
                ds(&g, "u_old_z", &o.u_old[2])?;
                ds(&g, "m_old", &o.m_old)?;
                ds(&g, "x_cur", &o.x_cur)?;
                ds(&g, "x_old", &o.x_old)?;
                for (i, v) in o.g_old.iter().enumerate() {
                    ds(&g, &format!("g_old_{i}"), v)?;
                }
            }
        }
        std::fs::rename(&tmp, path).with_context(|| format!("writing {}", path.display()))?;
        Ok(())
    }

    fn scheme_name(s: TimeScheme) -> &'static str {
        match s {
            TimeScheme::Bdf1 => "bdf1",
            TimeScheme::Bdf2 => "bdf2",
            TimeScheme::CrankNicolson => "crank-nicolson",
        }
    }

    /// Reads a restart file.
    pub fn read(path: &Path) -> Result<Restart> {
        let g = h5::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let fmt = rsattr(&g, "format")?;
        ensure!(
            fmt == "cfdcraft-restart",
            "{} is not a CFDCraft restart file (format `{fmt}`)",
            path.display()
        );
        let version: i64 = rattr(&g, "version")?;
        ensure!(
            version == 1,
            "restart file version {version} is not supported"
        );
        let transient = rattr::<u8>(&g, "transient")? == 1;
        let state = State {
            u: get3(&g, "u_")?,
            p: rd(&g, "p")?,
            mass_flux: rd(&g, "mass_flux")?,
            ub: get3(&g, "ub_")?,
            pb: rd(&g, "pb")?,
            grad_u: [
                get3(&g, "grad_u_x")?,
                get3(&g, "grad_u_y")?,
                get3(&g, "grad_u_z")?,
            ],
            grad_p: get3(&g, "grad_p_")?,
            psi: get3(&g, "psi_")?,
        };
        let amg_matrix = if g.link_exists("amg_matrix") {
            Some(rd(&g, "amg_matrix")?)
        } else {
            None
        };
        let tr = if transient {
            let scheme = match rsattr(&g, "scheme")?.as_str() {
                "bdf1" => TimeScheme::Bdf1,
                "bdf2" => TimeScheme::Bdf2,
                "crank-nicolson" => TimeScheme::CrankNicolson,
                other => bail!("unknown time scheme `{other}` in the restart file"),
            };
            Some((
                OldLevels {
                    u_old: [rd(&g, "u_old_x")?, rd(&g, "u_old_y")?, rd(&g, "u_old_z")?],
                    m_old: rd(&g, "m_old")?,
                    x_cur: rd(&g, "x_cur")?,
                    x_old: rd(&g, "x_old")?,
                    g_old: (0..12)
                        .map(|i| rd(&g, &format!("g_old_{i}")))
                        .collect::<Result<_>>()?,
                },
                Clock {
                    time: rattr(&g, "time")?,
                    step: rattr::<u64>(&g, "step")? as usize,
                    dt_prev: rattr(&g, "dt_prev")?,
                    levels: rattr::<u64>(&g, "levels")? as usize,
                },
                rattr(&g, "dt_next")?,
                scheme,
            ))
        } else {
            None
        };
        Ok(Restart {
            case: rsattr(&g, "case")?,
            cell_id: g.dataset("cell_id")?.read_raw::<u64>()?,
            face_id: g.dataset("face_id")?.read_raw::<u64>()?,
            state,
            iteration: rattr::<u64>(&g, "iteration")? as usize,
            amg_setups: rattr::<u64>(&g, "amg_setups")? as usize,
            amg_matrix,
            transient: tr,
        })
    }
}

/// Restart files are available (the `hdf5` feature); without it a run
/// writes none.
pub const AVAILABLE: bool = cfg!(feature = "hdf5");

/// Writes a restart file.
pub fn write(path: &Path, r: &Restart) -> Result<()> {
    #[cfg(feature = "hdf5")]
    {
        io::write(path, r)
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = (path, r);
        anyhow::bail!("restart files need the `hdf5` feature")
    }
}

/// Reads a restart file.
pub fn read(path: &Path) -> Result<Restart> {
    #[cfg(feature = "hdf5")]
    {
        io::read(path)
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = path;
        anyhow::bail!("restart files need the `hdf5` feature")
    }
}

/// Position in the file of each local entry, by global identity; `None`
/// when the order is the same.
fn mapping(file_ids: &[u64], local_ids: &[u64]) -> Result<Option<Vec<usize>>> {
    anyhow::ensure!(
        file_ids.len() == local_ids.len(),
        "the restart file has {} entries, the mesh {} (another mesh?)",
        file_ids.len(),
        local_ids.len()
    );
    if file_ids == local_ids {
        return Ok(None);
    }
    let pos: std::collections::HashMap<u64, usize> =
        file_ids.iter().enumerate().map(|(i, &g)| (g, i)).collect();
    let m: Option<Vec<usize>> = local_ids.iter().map(|g| pos.get(g).copied()).collect();
    m.map(Some)
        .ok_or_else(|| anyhow::anyhow!("the restart file belongs to another mesh"))
}

fn remap(v: &[f64], m: &[usize]) -> Vec<f64> {
    m.iter().map(|&i| v[i]).collect()
}

fn remap3(v: &Soa3, m: &[usize]) -> Soa3 {
    Soa3 {
        x: remap(&v.x, m),
        y: remap(&v.y, m),
        z: remap(&v.z, m),
    }
}

impl Restart {
    /// The content in the local order of a mesh with the identities
    /// `cell_id`, `face_id` (n_interior interior faces); false in the second
    /// value when the order differed (the multigrid matrix is then dropped).
    pub fn to_local(
        mut self,
        cell_id: &[u64],
        face_id: &[u64],
        n_interior: usize,
    ) -> Result<(Self, bool)> {
        let n_cells = cell_id.len();
        let mc = mapping(&self.cell_id, cell_id)?;
        let mf = mapping(&self.face_id, face_id)?;
        let same = mc.is_none() && mf.is_none();
        let ident = |n: usize| (0..n).collect::<Vec<usize>>();
        if !same {
            let mc = mc.unwrap_or_else(|| ident(n_cells));
            let mf = mf.unwrap_or_else(|| ident(face_id.len()));
            // boundary faces: the positions among the boundary faces
            let mb: Vec<usize> = mf[n_interior..].iter().map(|&i| i - n_interior).collect();
            anyhow::ensure!(
                mf[n_interior..].iter().all(|&i| i >= n_interior),
                "the restart file has another split of interior and boundary faces"
            );
            let s = &self.state;
            self.state = State {
                u: remap3(&s.u, &mc),
                p: remap(&s.p, &mc),
                mass_flux: remap(&s.mass_flux, &mf),
                ub: remap3(&s.ub, &mb),
                pb: remap(&s.pb, &mb),
                grad_u: [
                    remap3(&s.grad_u[0], &mc),
                    remap3(&s.grad_u[1], &mc),
                    remap3(&s.grad_u[2], &mc),
                ],
                grad_p: remap3(&s.grad_p, &mc),
                psi: remap3(&s.psi, &mc),
            };
            if let Some((o, ..)) = &mut self.transient {
                *o = OldLevels {
                    u_old: [
                        remap(&o.u_old[0], &mc),
                        remap(&o.u_old[1], &mc),
                        remap(&o.u_old[2], &mc),
                    ],
                    m_old: remap(&o.m_old, &mf),
                    x_cur: remap(&o.x_cur, &mf),
                    x_old: remap(&o.x_old, &mf),
                    g_old: o.g_old.iter().map(|v| remap(v, &mc)).collect(),
                };
            }
            self.amg_matrix = None;
            self.cell_id = cell_id.to_vec();
            self.face_id = face_id.to_vec();
        }
        Ok((self, same))
    }
}

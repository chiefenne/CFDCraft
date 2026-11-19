//! The sliding mesh with conformal stepping (D56, D182, variant A): a zone
//! of the mesh (the rotor) that turns about an axis through the time steps
//! of a transient run, past the cells of the rest of the mesh (the stator).
//! Case file: one `[[sliding]]` table.
//!
//! # The method
//!
//! The rotor mesh turns **rigidly**: it keeps its cells, volumes and faces,
//! so the geometric conservation law holds exactly. When it turns by a whole
//! number k of its sectors per time step (a rotor of N identical sectors,
//! Δt = k · 2π / (N |Ω|)), every rotor cell lands on the place of another
//! rotor cell, and the cells of the interface meet new neighbours without
//! any interpolation: the mesh stays **conformal**. The solver keeps the
//! mesh where it is: the fields are the fields at the places of the mesh
//! (Eulerian), also the old levels, from which the gradients, the
//! Rhie–Chow differences and the extrapolations are made. The only thing
//! that sees the motion is the **time term**: the old level of a cell that
//! moved is what that cell held, which is the Eulerian old level at the
//! place that it came from, k sectors back for the level n and 2k for the
//! level n−1. So the time term of the momentum and scalar rows of a rotor
//! cell reads the old levels with the rotor cells shifted by k and 2k
//! sectors ([`Plan::cell`]); everything else is the solver of a fixed mesh.
//!
//! # Equations
//!
//! The solved velocity is the absolute one **u** in the inertial basis
//! (the vectors of the shifted old levels do not turn with their cells: the
//! pressure gradient turns them, as it must), the equations are those of a
//! moving mesh,
//!
//! d/dt ∫ ρ **u** dV + ∮ ρ (**u** − **w**)·**S** **u** = −∮ p **S** + ∮ τ·**S**,
//!
//! with the grid velocity **w** = **Ω** × (**x** − **x**₀) of the rotor and
//! **no Coriolis source**. The convecting flux is the relative flux ṁ −
//! ρ **S**·**w**, the frame flux of the moving reference frame
//! ([`crate::frame`]), which is the grid flux of a rigid turn: Σ_f ρ
//! **S**·**w** = 0 for every closed rotor cell (a rigid turn keeps volumes),
//! and ρ **S**·**w** = 0 on an interface that is a surface of revolution
//! about the axis. Continuity, the Rhie–Chow flux and the pressure
//! correction work on the absolute flux. Scalars and turbulence variables
//! are convected with the relative flux; a scalar that is carried by the
//! rotor (relative flux 0) is carried exactly by the shift, without
//! numerical diffusion.
//!
//! # Limits
//!
//! - The rotor mesh must be a ring of identical sectors: it maps onto itself
//!   when turned by 2π/N about the axis (cells, faces, volumes). The
//!   sectors are found from the geometry (`SLD004` says what failed).
//! - The step is fixed and Δt = k · 2π / (N |Ω|) (`SLD003`).
//! - The boundary faces of the rotor must be walls that turn with it, or
//!   symmetry planes (`SLD004`); a wall with a prescribed velocity expression
//!   (D53) of the rigid turn does that.
//! - One rotor per case; laminar flow and passive scalars (a turbulence
//!   model is refused so far).
//! - On the CPU, on one process (the shifted old levels are made through the
//!   host).

use crate::case::FlowCase;
use crate::diag::{self, Diagnostic};
use crate::frame::{MrfSection, assign_cells_in};
use crate::simplec::Flow;
use crate::transient::Transient;
use anyhow::Result;
use foehncraft_fv::Discretisation;
use foehncraft_mesh::Mesh;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::f64::consts::TAU;

/// `[[sliding]]`: a zone that turns about an axis through the time steps
/// (D56, D182).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct SlidingSection {
    /// The cell zone of the mesh that turns (a zone name of the mesh;
    /// without it all cells turn).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub zone: Option<String>,
    /// A point on the axis of rotation (m).
    pub origin: [f64; 3],
    /// Direction of the axis (any length; the rotation is right-handed about
    /// it for a positive speed).
    pub axis: [f64; 3],
    /// Angular velocity (rad/s). Give this or `rpm`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub omega: Option<f64>,
    /// Angular velocity in revolutions per minute. Give this or `omega`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rpm: Option<f64>,
    /// Sectors of the rotor that it turns in one time step (default 1); the
    /// time step is `shift` × 2π / (sectors × |Ω|).
    pub shift: usize,
}

impl Default for SlidingSection {
    fn default() -> Self {
        Self {
            zone: None,
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            omega: None,
            rpm: None,
            shift: 1,
        }
    }
}

impl SlidingSection {
    /// The rotating zone as a frame of [`crate::frame`] (the same axis,
    /// speed and zone).
    pub fn frame(&self) -> MrfSection {
        MrfSection {
            zone: self.zone.clone(),
            origin: self.origin,
            axis: self.axis,
            omega: self.omega,
            rpm: self.rpm,
        }
    }

    /// The angular velocity vector (rad/s).
    pub fn omega_vector(&self) -> Result<[f64; 3]> {
        self.frame().omega_vector()
    }

    /// The time step that turns the rotor by `shift` of `sectors` sectors.
    pub fn time_step(&self, sectors: usize) -> Result<f64> {
        let o = self.omega_vector()?;
        Ok(self.shift as f64 * TAU / (sectors as f64 * norm(o)))
    }
}

/// Relative deviation of the time step from `shift` × 2π / (N |Ω|) that is
/// still accepted (a time step given with a few digits).
const STEP_TOLERANCE: f64 = 1e-6;

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn err(code: &'static str, key: &str, msg: String) -> anyhow::Error {
    let mut d = Diagnostic::new(code, msg);
    d.key = Some(key.into());
    diag::err(d)
}

/// The checks of the `[[sliding]]` tables that need no mesh (D182).
pub fn validate(case: &FlowCase) -> Result<()> {
    if case.sliding.is_empty() {
        return Ok(());
    }
    let Some(t) = case.time else {
        return Err(err(
            diag::SLIDING_STEADY,
            "sliding",
            "[[sliding]] moves the mesh in time (D182): the case needs a [time] table; a steady \
             case with a rotating zone uses [[mrf]] (D55)"
                .into(),
        ));
    };
    if case.sliding.len() > 1 {
        return Err(err(
            diag::SLIDING_ENTRY,
            "sliding",
            "[[sliding]]: one rotor per case so far".into(),
        ));
    }
    let s = &case.sliding[0];
    if let Err(e) = s.omega_vector() {
        return Err(err(
            diag::SLIDING_ENTRY,
            "sliding[0]",
            format!("[[sliding]]: {e}"),
        ));
    }
    crate::need!(
        s.omega_vector().is_ok_and(|o| norm(o) > 0.0),
        SLIDING_ENTRY,
        "[[sliding]]: the angular velocity must not be 0 (a rotor at rest needs no sliding mesh)"
    );
    crate::need!(
        s.shift >= 1,
        SLIDING_ENTRY,
        "[[sliding]]: `shift` (sectors per time step) must be at least 1"
    );
    crate::need!(
        !t.adaptive && t.error_tolerance.is_none(),
        SLIDING_STEPPING,
        "[[sliding]]: the time step is fixed (shift × 2π / (sectors × Ω)); `adaptive` and \
         `error_tolerance` are not available"
    );
    let v = &case.verification;
    crate::need!(
        v.levels.len() <= 1 && v.time_steps.is_empty(),
        SLIDING_STEPPING,
        "[[sliding]]: the time step follows from the mesh, so a sweep over levels or time steps \
         is a series of runs, not a [verification] table"
    );
    Ok(())
}

/// A rotor zone of the mesh and the shift of its cells per time step.
#[derive(Debug, Clone, PartialEq)]
pub struct Plan {
    /// Identical sectors of the rotor (N).
    pub sectors: usize,
    /// Sectors per step (k).
    pub shift: usize,
    /// The angle turned in one step (rad), in the sense of the rotation.
    pub angle: f64,
    /// Cells of the rotor.
    pub rotor_cells: usize,
    /// Entry i of the shifted array is entry `cell[i]` of the old one: the
    /// cell that is now at the place of cell i was at the place of cell
    /// `cell[i]` one step ago (identity outside the rotor).
    pub cell: Vec<usize>,
}

/// Points that can be found again within a tolerance.
struct PointIndex<'a> {
    pts: &'a [[f64; 3]],
    bin: f64,
    tol: f64,
    map: HashMap<[i64; 3], Vec<u32>>,
}

impl<'a> PointIndex<'a> {
    fn new(pts: &'a [[f64; 3]], tol: f64) -> Self {
        let bin = 16.0 * tol;
        let mut map: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
        for (i, p) in pts.iter().enumerate() {
            map.entry(Self::key(*p, bin)).or_default().push(i as u32);
        }
        Self { pts, bin, tol, map }
    }

    fn key(p: [f64; 3], bin: f64) -> [i64; 3] {
        [
            (p[0] / bin).floor() as i64,
            (p[1] / bin).floor() as i64,
            (p[2] / bin).floor() as i64,
        ]
    }

    /// The nearest point within the tolerance of `q`.
    fn find(&self, q: [f64; 3]) -> Option<usize> {
        let k = Self::key(q, self.bin);
        let mut best: Option<(f64, usize)> = None;
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(v) = self.map.get(&[k[0] + dx, k[1] + dy, k[2] + dz]) {
                        for &i in v {
                            let d = norm(sub(self.pts[i as usize], q));
                            if d <= self.tol && best.is_none_or(|(b, _)| d < b) {
                                best = Some((d, i as usize));
                            }
                        }
                    }
                }
            }
        }
        best.map(|(_, i)| i)
    }
}

/// Rodrigues' rotation of `v` about the unit axis `a` by `angle`.
fn rotate(a: [f64; 3], angle: f64, v: [f64; 3]) -> [f64; 3] {
    let (s, c) = angle.sin_cos();
    let axv = cross(a, v);
    let adv = dot(a, v);
    std::array::from_fn(|i| v[i] * c + axv[i] * s + a[i] * adv * (1.0 - c))
}

/// The shift of the data of the cells `rotor` of `disc` by the turn of
/// `entry` (see the module documentation).
pub fn plan(entry: &SlidingSection, rotor: &[bool], disc: &Discretisation) -> Result<Plan> {
    let o = entry.origin;
    let omega = entry.omega_vector()?;
    let w = norm(omega);
    let a = omega.map(|x| x / w);
    let nc = disc.n_cells;
    let cells: Vec<usize> = (0..nc).filter(|&c| rotor[c]).collect();
    let fail = |msg: String| err(diag::SLIDING_ZONE, "sliding", msg);
    if cells.is_empty() {
        return Err(fail("[[sliding]]: the rotor has no cells".into()));
    }
    let centre = |c: usize| disc.cell_centroid.get(c);
    // the size of the rotor and the tolerance of the matches
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for &c in &cells {
        let x = centre(c);
        for i in 0..3 {
            lo[i] = lo[i].min(x[i]);
            hi[i] = hi[i].max(x[i]);
        }
    }
    let size = norm(sub(hi, lo)).max(f64::MIN_POSITIVE);
    let tol = 1e-6 * size;
    // the sectors: the cells at the radius and height of the outermost cell
    let cyl = |x: [f64; 3]| {
        let d = sub(x, o);
        let z = dot(d, a);
        let r = norm(sub(d, a.map(|ai| ai * z)));
        (z, r)
    };
    let c0 = *cells
        .iter()
        .max_by(|&&p, &&q| cyl(centre(p)).1.total_cmp(&cyl(centre(q)).1))
        .expect("a cell");
    let (z0, r0) = cyl(centre(c0));
    let sectors = cells
        .iter()
        .filter(|&&c| {
            let (z, r) = cyl(centre(c));
            (z - z0).abs() <= tol && (r - r0).abs() <= tol
        })
        .count();
    if sectors < 2 || r0 <= tol {
        return Err(fail(format!(
            "[[sliding]]: the rotor is not a ring of identical sectors about the axis {:?} through \
             {:?}: no second cell lies at the radius and height of the outermost cell",
            entry.axis, entry.origin
        )));
    }
    if entry.shift.is_multiple_of(sectors) {
        return Err(err(
            diag::SLIDING_ENTRY,
            "sliding[0].shift",
            format!(
                "[[sliding]]: a shift of {} sectors is a whole turn of the {sectors} sectors of the rotor",
                entry.shift
            ),
        ));
    }
    let angle = TAU * entry.shift as f64 / sectors as f64;
    let turn = |x: [f64; 3]| {
        let r = rotate(a, angle, sub(x, o));
        [o[0] + r[0], o[1] + r[1], o[2] + r[2]]
    };
    // cells: the image of each rotor cell
    let pts: Vec<[f64; 3]> = cells.iter().map(|&c| centre(c)).collect();
    let index = PointIndex::new(&pts, tol);
    let mut image = vec![usize::MAX; nc];
    let mut cell: Vec<usize> = (0..nc).collect();
    let mut missing = 0usize;
    let mut first_missing = None;
    let mut taken = vec![false; nc];
    for &c in &cells {
        let target = turn(centre(c));
        match index.find(target) {
            Some(i) => {
                let t = cells[i];
                let dv = (disc.cell_volume[t] - disc.cell_volume[c]).abs();
                if dv > 1e-8 * disc.cell_volume[c] || taken[t] {
                    missing += 1;
                    first_missing.get_or_insert(c);
                } else {
                    taken[t] = true;
                    image[c] = t;
                    cell[t] = c;
                }
            }
            None => {
                missing += 1;
                first_missing.get_or_insert(c);
            }
        }
    }
    if missing > 0 {
        let c = first_missing.expect("a cell");
        return Err(fail(format!(
            "[[sliding]]: {missing} of {} rotor cells have no image when the rotor is turned by \
             {:.4}° (1/{sectors} of a turn, the first one is the cell at {:?}): the rotor mesh is not \
             a ring of {sectors} identical sectors about the axis {:?} through {:?}",
            cells.len(),
            angle.to_degrees(),
            centre(c),
            entry.axis,
            entry.origin
        )));
    }
    Ok(Plan {
        sectors,
        shift: entry.shift,
        angle,
        rotor_cells: cells.len(),
        cell,
    })
}

/// What [`attach`] reports for the log and the diagnostics.
#[derive(Debug, Default)]
pub struct SlidingNotes {
    /// Log lines.
    pub lines: Vec<String>,
}

/// The sliding mesh of a run.
#[derive(Debug, Clone)]
pub struct Sliding {
    /// The shift of the data.
    pub plan: Plan,
    /// The time step of the run (s).
    pub dt: f64,
}

/// The sliding mesh of `case` on `disc`, or `None` without `[[sliding]]`;
/// `dt` is the time step of the run.
pub fn attach(
    case: &FlowCase,
    mesh: &Mesh,
    disc: &Discretisation,
    dt: f64,
) -> Result<Option<(Sliding, SlidingNotes)>> {
    validate(case)?;
    let Some(entry) = case.sliding.first() else {
        return Ok(None);
    };
    if !case.turbulence.is_laminar() {
        return Err(err(
            diag::SLIDING_ZONE,
            "turbulence.model",
            "[[sliding]]: the sliding mesh runs laminar flow and passive scalars so far, not a \
             turbulence model"
                .into(),
        ));
    }
    if let Some(name) = &entry.zone
        && !mesh.zones.contains(name)
    {
        let hint = crate::bc::did_you_mean(name, mesh.zones.iter().map(String::as_str))
            .map(|c| format!(" (did you mean `{c}`?)"))
            .unwrap_or_default();
        return Err(err(
            diag::SLIDING_ENTRY,
            "sliding[0].zone",
            format!(
                "[[sliding]]: the mesh has no cell zone `{name}`{hint}; its zones are {:?}",
                mesh.zones
            ),
        ));
    }
    let cell_zone = assign_cells_in(&[entry.frame()], mesh, "[[sliding]]")?;
    let rotor: Vec<bool> = cell_zone.iter().map(Option::is_some).collect();
    // the boundaries of the rotor: walls that turn with it, or symmetry
    for b in &mesh.boundaries {
        if !b
            .range()
            .any(|f| rotor.get(disc.owner[f] as usize).copied().unwrap_or(false))
        {
            continue;
        }
        let Some(bc) = case.boundary.get(&b.name) else {
            continue;
        };
        use crate::bc::BoundaryCondition as Bc;
        if !matches!(
            bc,
            Bc::Wall { .. } | Bc::SlipWall | Bc::Symmetry | Bc::Wedge
        ) {
            return Err(err(
                diag::SLIDING_ZONE,
                &format!("boundary.{}", b.name),
                format!(
                    "[[sliding]]: the boundary `{}` of the rotor is a `{}`; the boundary faces of \
                     the rotor must be walls that turn with it, slip walls or symmetry planes",
                    b.name,
                    bc.type_name()
                ),
            ));
        }
    }
    let plan = plan(entry, &rotor, disc)?;
    let want = entry.time_step(plan.sectors)?;
    if (dt - want).abs() > STEP_TOLERANCE * want {
        return Err(err(
            diag::SLIDING_STEPPING,
            "time.time_step",
            format!(
                "[[sliding]]: the rotor turns {} of its {} sectors in one step, so the time step \
                 must be {want:.10e} s, not {dt:.10e} s",
                plan.shift, plan.sectors
            ),
        ));
    }
    if let Some(t) = case.time
        && let Some(e) = t.end_time
    {
        let n = (e - t.start_time) / dt;
        if (n - n.round()).abs() > 1e-6 {
            return Err(err(
                diag::SLIDING_STEPPING,
                "time.end_time",
                format!(
                    "[[sliding]]: `end_time` must be a whole number of steps after `start_time` \
                     ({n:.6} steps of {dt:.6e} s); give `steps`, or {:.10e} s",
                    t.start_time + n.round() * dt
                ),
            ));
        }
    }
    let notes = SlidingNotes {
        lines: vec![format!(
            "sliding mesh (D182): {} ({} cells) turns about the axis {:?} through {:?} at {:.6e} \
             rad/s; {} sectors, {} per step = {:.4}° per step, conformal stepping",
            entry
                .zone
                .as_deref()
                .map_or("all cells".to_string(), |z| format!("zone `{z}`")),
            plan.rotor_cells,
            entry.axis,
            entry.origin,
            norm(entry.omega_vector()?),
            plan.sectors,
            plan.shift,
            plan.angle.to_degrees()
        )],
    };
    Ok(Some((Sliding { plan, dt }, notes)))
}

impl Sliding {
    /// Tells the transient run how the old levels of the time terms move
    /// (call once, before the first step).
    pub fn install(&self, fl: &mut Flow, tr: &mut Transient) -> Result<()> {
        tr.set_sliding(fl, self.plan.cell.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn time_step_from_sectors_and_speed() {
        let s = SlidingSection {
            omega: Some(10.0),
            shift: 2,
            ..Default::default()
        };
        let dt = s.time_step(64).unwrap();
        assert!((dt - 2.0 * TAU / (64.0 * 10.0)).abs() < 1e-15);
        let rpm = SlidingSection {
            rpm: Some(60.0),
            ..Default::default()
        };
        assert!((rpm.time_step(36).unwrap() - 1.0 / 36.0).abs() < 1e-15);
        assert!(SlidingSection::default().time_step(8).is_err());
    }

    #[test]
    fn rotation_about_an_axis() {
        let r = rotate(
            [0.0, 0.0, 1.0],
            std::f64::consts::FRAC_PI_2,
            [1.0, 0.0, 0.5],
        );
        assert!(r[0].abs() < 1e-15 && (r[1] - 1.0).abs() < 1e-15 && (r[2] - 0.5).abs() < 1e-15);
    }

    #[test]
    fn points_are_found_within_the_tolerance() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0 - 1e-9]];
        let idx = PointIndex::new(&pts, 1e-6);
        assert_eq!(idx.find([1.0, 1e-8, 0.0]), Some(1));
        assert_eq!(idx.find([1.0, 1.0, 1.0]), Some(2));
        assert_eq!(idx.find([0.5, 0.5, 0.5]), None);
    }
}

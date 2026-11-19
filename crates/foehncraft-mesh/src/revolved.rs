//! The built-in mesh of revolution (`[mesh.revolved]`, D186): a block grid in
//! the meridional (r, z) half-plane turned about the z axis, with a rotor and
//! a stator that are two separate meshes of their own pitch around the axis.
//!
//! The blocks are the rectangles between the breakpoints `r` and `z`; the
//! blocks listed in `rotor` form the cell zone `rotor`, the others the zone
//! `stator`. Where a rotor block meets a stator block, the faces are the two
//! boundaries `rotor_interface` (with `n_theta` intervals around the axis)
//! and `stator_interface` (with `n_theta_stator`): an interface that is a
//! surface of revolution of cylinders, discs and, with `cone`, cones, for the
//! sliding mesh by imprinting (`[[sliding]]`). With `r[0] = 0` the grid
//! reaches the axis; the cells there are wedges about one node per z on the
//! axis, so a disc of the interface can pass through the centre.
//!
//! `cone` (degrees) tilts the inner z breakpoints: the line z = z[b] becomes
//! z = z[b] + (r − r_mid) tan(cone), r_mid the middle of the r range, so
//! interfaces on those lines are cones.
//!
//! Boundaries: `inner` (r = r[0], only when it is not the axis), `outer`,
//! `zmin`, `zmax`, and `wall` around the blocks listed in `solid`; the faces
//! of rotor cells on them get the suffix `_rotor` (`inner_rotor`, a hub that
//! turns with the rotor; `wall_rotor`, an impeller disc); and the two
//! interface boundaries. A boundary appears when it has faces.

use crate::builder::MeshBuilder;
use crate::geometry::{dot, polygon_area_centroid, sub};
use crate::mesh::Mesh;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::f64::consts::TAU;

/// The two boundaries of the interface: the rotor's faces and the stator's.
pub const INTERFACE_BOUNDARIES: [&str; 2] = ["rotor_interface", "stator_interface"];

/// Specification of a mesh of revolution (`[mesh.revolved]`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct RevolvedSpec {
    /// Radial breakpoints (m), increasing; `r[0] = 0` reaches the axis.
    pub r: Vec<f64>,
    /// Axial breakpoints (m), increasing.
    pub z: Vec<f64>,
    /// Radial intervals of every r block (`r.len() - 1` values).
    pub n_r: Vec<usize>,
    /// Axial intervals of every z block (`z.len() - 1` values).
    pub n_z: Vec<usize>,
    /// The blocks of the rotor, as [r block, z block] (counted from 0).
    pub rotor: Vec<[usize; 2]>,
    /// Blocks without cells (a solid body: a disc, a hub, a baffle); their
    /// faces are the boundary `wall` (`wall_rotor` on rotor cells, turning
    /// with the rotor).
    pub solid: Vec<[usize; 2]>,
    /// Intervals around the axis of the rotor.
    pub n_theta: usize,
    /// Intervals around the axis of the stator (0: as `n_theta`; the
    /// interface is split into two boundaries also then).
    pub n_theta_stator: usize,
    /// Tilt of the inner z breakpoints (degrees): interfaces on them are
    /// cones (0: discs).
    pub cone: f64,
    /// Refinement level: every interval is halved `refine` times.
    pub refine: u32,
}

impl Default for RevolvedSpec {
    fn default() -> Self {
        // an axial stage: a rotor row on the hub between two stator parts,
        // the interface a disc, a cylinder and a disc
        Self {
            r: vec![0.05, 0.09, 0.10],
            z: vec![0.0, 0.05, 0.08, 0.13],
            n_r: vec![8, 2],
            n_z: vec![10, 6, 10],
            rotor: vec![[0, 1]],
            solid: Vec::new(),
            n_theta: 48,
            n_theta_stator: 64,
            cone: 0.0,
            refine: 0,
        }
    }
}

impl RevolvedSpec {
    /// Intervals around the axis of the rotor and of the stator after
    /// refinement.
    pub fn theta_intervals(&self) -> (usize, usize) {
        let f = 1usize << self.refine;
        let ns = if self.n_theta_stator == 0 {
            self.n_theta
        } else {
            self.n_theta_stator
        };
        (self.n_theta * f, ns * f)
    }

    /// Checks the specification.
    pub fn validate(&self) -> Result<(), String> {
        let e = |m: String| Err(format!("[mesh.revolved]: {m}"));
        if self.r.len() < 2 || self.z.len() < 2 {
            return e("need at least two breakpoints in r and in z".into());
        }
        if self.r.windows(2).any(|w| !(w[1] > w[0])) || !(self.r[0] >= 0.0) {
            return e(format!("r must increase from r[0] >= 0, got {:?}", self.r));
        }
        if self.z.windows(2).any(|w| !(w[1] > w[0])) {
            return e(format!("z must increase, got {:?}", self.z));
        }
        if self.n_r.len() != self.r.len() - 1 || self.n_z.len() != self.z.len() - 1 {
            return e(format!(
                "need {} values in n_r and {} in n_z (one per block)",
                self.r.len() - 1,
                self.z.len() - 1
            ));
        }
        if self.n_r.iter().chain(&self.n_z).any(|&n| n == 0) {
            return e("every block needs at least one interval".into());
        }
        if self.n_theta < 3 || (self.n_theta_stator > 0 && self.n_theta_stator < 3) {
            return e("need n_theta >= 3 and n_theta_stator >= 3".into());
        }
        let nb = (self.r.len() - 1) * (self.z.len() - 1);
        if self.rotor.is_empty() || self.rotor.len() >= nb {
            return e("`rotor` must name some, but not all, of the blocks".into());
        }
        if let Some(b) = self
            .rotor
            .iter()
            .chain(&self.solid)
            .find(|b| b[0] + 1 >= self.r.len() || b[1] + 1 >= self.z.len())
        {
            return e(format!("the block {b:?} is outside the grid"));
        }
        if let Some(b) = self.solid.iter().find(|b| self.rotor.contains(b)) {
            return e(format!("the block {b:?} is both rotor and solid"));
        }
        if self.rotor.len() + self.solid.len() >= nb {
            return e("the stator needs at least one block".into());
        }
        if !(self.cone.abs() < 80.0) {
            return e(format!(
                "cone must lie within ±80 degrees, got {}",
                self.cone
            ));
        }
        // the tilted breakpoints must stay in order at both ends of r
        for &r in [self.r[0], self.r[self.r.len() - 1]].iter() {
            let z: Vec<f64> = (0..self.z.len()).map(|b| self.z_line(b, r)).collect();
            if z.windows(2).any(|w| !(w[1] > w[0])) {
                return e(format!(
                    "with cone = {} degrees the z breakpoints cross at r = {r} m",
                    self.cone
                ));
            }
        }
        Ok(())
    }

    /// The breakpoint line `b` at radius `r`.
    fn z_line(&self, b: usize, r: f64) -> f64 {
        if b == 0 || b + 1 == self.z.len() || self.cone == 0.0 {
            return self.z[b];
        }
        let r_mid = 0.5 * (self.r[0] + self.r[self.r.len() - 1]);
        self.z[b] + (r - r_mid) * self.cone.to_radians().tan()
    }

    /// Number of cells.
    pub fn estimated_cells(&self) -> usize {
        let f = 1usize << self.refine;
        let (nr, ns) = self.theta_intervals();
        let mut n = 0;
        for (i, &a) in self.n_r.iter().enumerate() {
            for (k, &b) in self.n_z.iter().enumerate() {
                if self.solid.contains(&[i, k]) {
                    continue;
                }
                let t = if self.rotor.contains(&[i, k]) { nr } else { ns };
                n += a * b * f * f * t;
            }
        }
        n
    }

    /// One line for the log.
    pub fn describe(&self) -> String {
        let (nr, ns) = self.theta_intervals();
        format!(
            "built-in mesh of revolution: {} rotor block(s) with {nr} and the stator with {ns} \
             intervals around the axis{}{}, refine {}",
            self.rotor.len(),
            if self.r[0] == 0.0 {
                ", through the axis"
            } else {
                ""
            },
            if self.cone != 0.0 {
                format!(", cones at {} degrees", self.cone)
            } else {
                String::new()
            },
            self.refine
        )
    }

    /// Builds the mesh.
    pub fn generate(&self) -> Result<Mesh, String> {
        self.validate()?;
        let f = 1usize << self.refine;
        // the meridional nodes: r_i, and z of (i, k)
        let mut rs: Vec<f64> = vec![self.r[0]];
        let mut r_block: Vec<usize> = Vec::new();
        for (b, &n) in self.n_r.iter().enumerate() {
            let n = n * f;
            for a in 1..=n {
                rs.push(self.r[b] + (self.r[b + 1] - self.r[b]) * a as f64 / n as f64);
                r_block.push(b);
            }
        }
        let mut z_frac: Vec<(usize, f64)> = vec![(0, 0.0)];
        let mut z_block: Vec<usize> = Vec::new();
        for (b, &n) in self.n_z.iter().enumerate() {
            let n = n * f;
            for a in 1..=n {
                z_frac.push((b, a as f64 / n as f64));
                z_block.push(b);
            }
        }
        let zz = |i: usize, k: usize| -> f64 {
            let (b, t) = z_frac[k];
            let (z0, z1) = (self.z_line(b, rs[i]), self.z_line(b + 1, rs[i]));
            z0 + t * (z1 - z0)
        };
        let (ni, nk) = (rs.len() - 1, z_frac.len() - 1);
        let axis = self.r[0] == 0.0;
        let (ntr, nts) = self.theta_intervals();
        // the zone of a meridional cell: 0 stator, 1 rotor, None solid
        let kind = |i: usize, k: usize| -> Option<u32> {
            let b = [r_block[i], z_block[k]];
            (!self.solid.contains(&b)).then(|| u32::from(self.rotor.contains(&b)))
        };
        let n_theta = |zone: u32| if zone == 1 { ntr } else { nts };
        // nodes, per zone, made when first used
        let mut points: Vec<[f64; 3]> = Vec::new();
        let mut ids: HashMap<(u32, usize, usize, usize), u32> = HashMap::new();
        let mut node = |zone: u32, i: usize, k: usize, j: usize| -> u32 {
            let n = n_theta(zone);
            let j = if axis && i == 0 { 0 } else { j % n };
            *ids.entry((zone, i, k, j)).or_insert_with(|| {
                let th = TAU * j as f64 / n as f64;
                points.push([rs[i] * th.cos(), rs[i] * th.sin(), zz(i, k)]);
                (points.len() - 1) as u32
            })
        };
        // cells: meridional cell (i, k), sector j
        let mut first = vec![0u32; ni * nk];
        let mut cell_zone: Vec<u32> = Vec::new();
        let mut centre: Vec<[f64; 3]> = Vec::new();
        for i in 0..ni {
            for k in 0..nk {
                let Some(z) = kind(i, k) else { continue };
                first[i * nk + k] = cell_zone.len() as u32;
                let n = n_theta(z);
                let (rc, zc) = (
                    0.5 * (rs[i] + rs[i + 1]),
                    0.25 * (zz(i, k) + zz(i + 1, k) + zz(i, k + 1) + zz(i + 1, k + 1)),
                );
                for j in 0..n {
                    let th = TAU * (j as f64 + 0.5) / n as f64;
                    cell_zone.push(z);
                    centre.push([rc * th.cos(), rc * th.sin(), zc]);
                }
            }
        }
        let cell = |i: usize, k: usize, j: usize| -> u32 {
            let n = n_theta(kind(i, k).expect("a fluid cell"));
            first[i * nk + k] + (j % n) as u32
        };
        // faces: (nodes, owner, neighbour or boundary name)
        enum Side {
            Cell(u32),
            Named(String),
        }
        let mut faces: Vec<(Vec<u32>, u32, Side)> = Vec::new();
        let named = |base: &str, zone: u32| {
            if zone == 1 {
                format!("{base}_rotor")
            } else {
                base.to_string()
            }
        };
        let iface = |zone: u32| INTERFACE_BOUNDARIES[usize::from(zone == 0)].to_string();
        for i in 0..ni {
            for k in 0..nk {
                let Some(z) = kind(i, k) else { continue };
                let n = n_theta(z);
                // θ faces between the sectors j − 1 and j
                for j in 0..n {
                    let q = vec![
                        node(z, i, k, j),
                        node(z, i + 1, k, j),
                        node(z, i + 1, k + 1, j),
                        node(z, i, k + 1, j),
                    ];
                    faces.push((q, cell(i, k, j + n - 1), Side::Cell(cell(i, k, j))));
                }
            }
        }
        // r faces at the node line i, between the cells (i − 1, k) and (i, k)
        for i in 0..=ni {
            if axis && i == 0 {
                continue;
            }
            for k in 0..nk {
                // None: outside the grid; Some(None): a solid block
                let left = (i > 0).then(|| kind(i - 1, k));
                let right = (i < ni).then(|| kind(i, k));
                let quad = |node: &mut dyn FnMut(u32, usize, usize, usize) -> u32, z, j| {
                    vec![
                        node(z, i, k, j),
                        node(z, i, k, j + 1),
                        node(z, i, k + 1, j + 1),
                        node(z, i, k + 1, j),
                    ]
                };
                match (left, right) {
                    (Some(Some(a)), Some(Some(b))) if a == b => {
                        for j in 0..n_theta(a) {
                            let q = quad(&mut node, a, j);
                            faces.push((q, cell(i - 1, k, j), Side::Cell(cell(i, k, j))));
                        }
                    }
                    _ => {
                        for (side, zone) in [(i.wrapping_sub(1), left), (i, right)] {
                            let Some(Some(zone)) = zone else { continue };
                            let other = if side == i { left } else { right };
                            let name = match other {
                                Some(Some(_)) => iface(zone),
                                Some(None) => named("wall", zone),
                                None if side == i => named("inner", zone),
                                None => named("outer", zone),
                            };
                            for j in 0..n_theta(zone) {
                                let q = quad(&mut node, zone, j);
                                faces.push((q, cell(side, k, j), Side::Named(name.clone())));
                            }
                        }
                    }
                }
            }
        }
        // z faces at the node line k, between the cells (i, k − 1) and (i, k)
        for k in 0..=nk {
            for i in 0..ni {
                let below = (k > 0).then(|| kind(i, k - 1));
                let above = (k < nk).then(|| kind(i, k));
                let poly = |node: &mut dyn FnMut(u32, usize, usize, usize) -> u32, z, j| {
                    if axis && i == 0 {
                        vec![
                            node(z, i, k, 0),
                            node(z, i + 1, k, j),
                            node(z, i + 1, k, j + 1),
                        ]
                    } else {
                        vec![
                            node(z, i, k, j),
                            node(z, i + 1, k, j),
                            node(z, i + 1, k, j + 1),
                            node(z, i, k, j + 1),
                        ]
                    }
                };
                match (below, above) {
                    (Some(Some(a)), Some(Some(b))) if a == b => {
                        for j in 0..n_theta(a) {
                            let q = poly(&mut node, a, j);
                            faces.push((q, cell(i, k - 1, j), Side::Cell(cell(i, k, j))));
                        }
                    }
                    _ => {
                        for (side, zone) in [(k.wrapping_sub(1), below), (k, above)] {
                            let Some(Some(zone)) = zone else { continue };
                            let other = if side == k { below } else { above };
                            let name = match other {
                                Some(Some(_)) => iface(zone),
                                Some(None) => named("wall", zone),
                                None if side == k => named("zmin", zone),
                                None => named("zmax", zone),
                            };
                            for j in 0..n_theta(zone) {
                                let q = poly(&mut node, zone, j);
                                faces.push((q, cell(i, side, j), Side::Named(name.clone())));
                            }
                        }
                    }
                }
            }
        }
        let n_cells = cell_zone.len();
        // orientation: out of the owner
        let mut b = MeshBuilder::new(points.clone(), n_cells);
        b.set_zones(vec!["stator".into(), "rotor".into()], cell_zone);
        for (mut q, owner, side) in faces {
            let (s, c) = polygon_area_centroid(&points, &q);
            if dot(s, sub(c, centre[owner as usize])) < 0.0 {
                q.reverse();
            }
            match side {
                Side::Cell(nb) => b.add_interior(&q, owner, nb),
                Side::Named(name) => b.add_boundary(&q, owner, &name),
            }
        }
        b.finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Geometry;

    fn closes(m: &Mesh) -> f64 {
        let g = Geometry::compute(m);
        (0..m.n_cells)
            .map(|c| {
                let s = g.closure(m, c);
                let tot: f64 = m
                    .cell_face_iter(c)
                    .map(|(f, _)| crate::geometry::norm(g.face_area[f]))
                    .sum();
                crate::geometry::norm(s) / tot
            })
            .fold(0.0, f64::max)
    }

    #[test]
    fn the_stage_has_its_boundaries_and_closed_cells() {
        let s = RevolvedSpec::default();
        let m = s.generate().unwrap();
        assert_eq!(m.n_cells, s.estimated_cells());
        let names: Vec<&str> = m.boundaries.iter().map(|p| p.name.as_str()).collect();
        for want in [
            "inner",
            "outer",
            "zmin",
            "zmax",
            "inner_rotor",
            "rotor_interface",
            "stator_interface",
        ] {
            assert!(names.contains(&want), "{names:?}");
        }
        assert!(closes(&m) < 1e-13);
        let g = Geometry::compute(&m);
        assert!(g.cell_volume[..m.n_cells].iter().all(|&v| v > 0.0));
    }

    #[test]
    fn a_pot_through_the_axis_has_wedges() {
        let s = RevolvedSpec {
            r: vec![0.0, 0.04, 0.1],
            z: vec![0.0, 0.04, 0.08, 0.15],
            n_r: vec![4, 4],
            n_z: vec![3, 3, 4],
            rotor: vec![[0, 1]],
            n_theta: 12,
            n_theta_stator: 18,
            ..Default::default()
        };
        let m = s.generate().unwrap();
        assert_eq!(m.n_cells, s.estimated_cells());
        let names: Vec<&str> = m.boundaries.iter().map(|p| p.name.as_str()).collect();
        assert!(!names.contains(&"inner"), "{names:?}");
        assert!(closes(&m) < 1e-13);
        let g = Geometry::compute(&m);
        let vol: f64 = g.cell_volume[..m.n_cells].iter().sum();
        // the cylinder of the chords: (n / 2) sin(2 pi / n) r^2 h
        let n = 18.0;
        let want = 0.5 * n * (TAU / n).sin() * 0.1f64.powi(2) * 0.15;
        // the rotor has 12 intervals: its polygon is smaller
        assert!((vol - want).abs() < 0.03 * want, "{vol} {want}");
    }

    #[test]
    fn a_solid_block_is_a_wall() {
        // a disc impeller in the rotor of a pot
        let s = RevolvedSpec {
            r: vec![0.0, 0.03, 0.045, 0.1],
            z: vec![0.0, 0.03, 0.045, 0.055, 0.07, 0.1],
            n_r: vec![3, 2, 4],
            n_z: vec![3, 2, 1, 2, 3],
            rotor: vec![[0, 1], [0, 3], [1, 1], [1, 2], [1, 3]],
            solid: vec![[0, 2]],
            n_theta: 12,
            n_theta_stator: 18,
            ..Default::default()
        };
        let m = s.generate().unwrap();
        assert_eq!(m.n_cells, s.estimated_cells());
        let names: Vec<&str> = m.boundaries.iter().map(|p| p.name.as_str()).collect();
        assert!(
            names.contains(&"wall_rotor") && !names.contains(&"wall"),
            "{names:?}"
        );
        assert!(closes(&m) < 1e-13);
        let mut bad = s.clone();
        bad.solid.push([1, 1]);
        assert!(bad.validate().unwrap_err().contains("both rotor and solid"));
    }

    #[test]
    fn cones_tilt_the_inner_breakpoints() {
        let s = RevolvedSpec {
            cone: 20.0,
            ..Default::default()
        };
        let m = s.generate().unwrap();
        assert!(closes(&m) < 1e-13);
        let bad = RevolvedSpec {
            cone: 79.0,
            ..Default::default()
        };
        assert!(bad.validate().unwrap_err().contains("cross"));
    }
}

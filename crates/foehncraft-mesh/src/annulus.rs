//! Built-in mesh of an annulus (`[mesh.annulus]`): the gap between two
//! coaxial cylinders, the geometry of Taylor–Couette flow and the simplest
//! rotor–stator arrangement. One layer of cells in z by default (a 2D case),
//! `n_z` layers for a 3D annulus.
//!
//! A structured ring of hexahedra: `n_theta` intervals around the axis (the
//! z axis through the origin), `n_r` intervals from the inner to the outer
//! cylinder, optionally graded (`radial_ratio`, `symmetric`). `refine` halves
//! every interval; the gradings follow the continuous mapping of
//! [`crate::generate::Grading`], so the refined meshes are nested in the
//! radial direction.
//!
//! With `rotor_layers > 0` the innermost `rotor_layers` radial intervals form
//! the cell zone `rotor` and the rest the zone `stator`, so a rotating
//! zone (`[[mrf]]`, a sliding mesh) has an interface that is a surface of
//! revolution about the axis; without it all cells are in the zone `fluid`.
//!
//! With `n_theta_stator > 0` the rotor and the stator are two separate rings
//! with their own number of intervals around the axis (a pitch ratio that is
//! not one: the sliding mesh by imprinting, D182 variant B): `n_theta` is the
//! rotor's, `n_theta_stator` the stator's, and the faces of the interface are
//! the two boundaries `rotor_interface` and `stator_interface` instead of
//! interior faces.
//!
//! Boundaries: `inner` (r = `r_inner`), `outer` (r = `r_outer`), `zmin`,
//! `zmax`, and with a split interface `rotor_interface`, `stator_interface`.

use crate::builder::{BoundaryFace, CellSoup};
use crate::generate::{Grading, hex_faces};
use crate::mesh::Mesh;
use serde::{Deserialize, Serialize};

/// Boundary names of the annulus mesh.
pub const ANNULUS_BOUNDARIES: [&str; 4] = ["inner", "outer", "zmin", "zmax"];

/// The two boundaries of a split interface (`n_theta_stator > 0`): the
/// rotor's faces and the stator's.
pub const INTERFACE_BOUNDARIES: [&str; 2] = ["rotor_interface", "stator_interface"];

/// Specification of an annulus (`[mesh.annulus]`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct AnnulusSpec {
    /// Radius of the inner cylinder (m).
    pub r_inner: f64,
    /// Radius of the outer cylinder (m).
    pub r_outer: f64,
    /// Extent in z (m): one cell layer for `n_z = 1`.
    pub depth: f64,
    /// Intervals around the axis (of the rotor when `n_theta_stator` is set).
    pub n_theta: usize,
    /// Intervals around the axis of the stator ring (0: as `n_theta`, one
    /// conformal mesh). With a value (also equal to `n_theta`) the interface
    /// between rotor and stator is split into the boundaries
    /// `rotor_interface` and `stator_interface`, for a sliding mesh with a
    /// pitch ratio that is not one (needs `rotor_layers` > 0).
    pub n_theta_stator: usize,
    /// Intervals from the inner to the outer cylinder.
    pub n_r: usize,
    /// Cell layers in z (1: a 2D case; `refine` halves them for more).
    pub n_z: usize,
    /// Ratio of the last to the first radial interval (1: uniform); with
    /// `symmetric` of the middle to the end intervals.
    pub radial_ratio: f64,
    /// Cluster the radial intervals at both cylinders.
    pub symmetric: bool,
    /// Radial intervals (of the unrefined mesh) next to the inner cylinder
    /// that form the cell zone `rotor` (0: one zone `fluid`).
    pub rotor_layers: usize,
    /// Refinement level: every interval is halved `refine` times.
    pub refine: u32,
}

impl Default for AnnulusSpec {
    fn default() -> Self {
        Self {
            r_inner: 0.05,
            r_outer: 0.10,
            depth: 0.01,
            n_theta: 64,
            n_theta_stator: 0,
            n_r: 16,
            n_z: 1,
            radial_ratio: 1.0,
            symmetric: false,
            rotor_layers: 0,
            refine: 0,
        }
    }
}

impl AnnulusSpec {
    /// Intervals (around, radial, z, rotor layers) after refinement.
    pub fn intervals(&self) -> [usize; 4] {
        let f = 1usize << self.refine;
        // one layer in z (a 2D case) stays one layer
        let nz = if self.n_z == 1 { 1 } else { self.n_z * f };
        [self.n_theta * f, self.n_r * f, nz, self.rotor_layers * f]
    }

    /// Intervals around the axis of the stator ring after refinement.
    pub fn stator_intervals(&self) -> usize {
        let f = 1usize << self.refine;
        if self.n_theta_stator > 0 {
            self.n_theta_stator * f
        } else {
            self.n_theta * f
        }
    }

    /// The interface is split into two boundaries (`n_theta_stator` set).
    pub fn split(&self) -> bool {
        self.n_theta_stator > 0
    }

    /// Number of cells.
    pub fn estimated_cells(&self) -> usize {
        let [nt, nr, nz, nrot] = self.intervals();
        nt * nrot * nz + self.stator_intervals() * (nr - nrot) * nz
    }

    /// 2 for one cell layer (a 2D case), else 3.
    pub fn dimension(&self) -> u8 {
        if self.n_z == 1 { 2 } else { 3 }
    }

    /// Checks the spec.
    pub fn validate(&self) -> Result<(), String> {
        if !(self.r_inner > 0.0 && self.r_outer > self.r_inner && self.depth > 0.0) {
            return Err("[mesh.annulus]: need 0 < r_inner < r_outer and depth > 0".into());
        }
        if self.n_theta < 3 || self.n_r < 1 || self.n_z < 1 {
            return Err("[mesh.annulus]: need n_theta >= 3, n_r >= 1 and n_z >= 1".into());
        }
        if self.radial_ratio <= 0.0 {
            return Err("[mesh.annulus]: radial_ratio must be > 0".into());
        }
        if self.n_theta_stator > 0 && self.n_theta_stator < 3 {
            return Err("[mesh.annulus]: need n_theta_stator >= 3".into());
        }
        if self.n_theta_stator > 0 && self.rotor_layers == 0 {
            return Err(
                "[mesh.annulus]: n_theta_stator needs rotor_layers > 0 (the rotor is the part \
                 with the pitch n_theta)"
                    .into(),
            );
        }
        if self.rotor_layers >= self.n_r {
            return Err(format!(
                "[mesh.annulus]: rotor_layers ({}) must be smaller than n_r ({}): the stator \
                 needs at least one layer",
                self.rotor_layers, self.n_r
            ));
        }
        if self.refine > 6 {
            return Err("[mesh.annulus]: refine must be <= 6".into());
        }
        Ok(())
    }

    /// Radius of the interface between rotor and stator, if there is a
    /// rotor.
    pub fn rotor_radius(&self) -> Option<f64> {
        if self.rotor_layers == 0 {
            return None;
        }
        let [_, nr, _, nrot] = self.intervals();
        let s = self.radial_nodes();
        Some(self.r_inner + s[nrot.min(nr)] * (self.r_outer - self.r_inner))
    }

    fn radial_nodes(&self) -> Vec<f64> {
        let [_, nr, _, _] = self.intervals();
        Grading {
            ratio: self.radial_ratio,
            symmetric: self.symmetric,
        }
        .nodes(nr, self.refine)
    }

    /// Generates the mesh.
    pub fn generate(&self) -> Result<Mesh, String> {
        self.validate()?;
        if self.split() {
            return self.generate_split();
        }
        let [nt, nr, nz, nrot] = self.intervals();
        let s = self.radial_nodes();
        let tau = std::f64::consts::TAU;
        let (ri, ro) = (self.r_inner, self.r_outer);
        // 2D nodes, index j + (nr + 1) a
        let mut pts2: Vec<[f64; 2]> = Vec::with_capacity(nt * (nr + 1));
        for a in 0..nt {
            let th = tau * a as f64 / nt as f64;
            let (c, sn) = (th.cos(), th.sin());
            for &sj in &s {
                let r = ri + sj * (ro - ri);
                pts2.push([r * c, r * sn]);
            }
        }
        let node = |a: usize, j: usize| (j + (nr + 1) * (a % nt)) as u32;
        let np = pts2.len() as u32;
        let mut pts = Vec::with_capacity(pts2.len() * (nz + 1));
        for k in 0..=nz {
            let z = self.depth * k as f64 / nz as f64;
            pts.extend(pts2.iter().map(|p| [p[0], p[1], z]));
        }
        let mut soup = CellSoup::new(pts);
        if nrot > 0 {
            soup.set_zone_names(vec!["stator".to_string(), "rotor".to_string()]);
        }
        for k in 0..nz {
            for a in 0..nt {
                for j in 0..nr {
                    let q = [
                        node(a, j),
                        node(a, j + 1),
                        node(a + 1, j + 1),
                        node(a + 1, j),
                    ];
                    let (lo, hi) = (k as u32 * np, (k as u32 + 1) * np);
                    let hx = [
                        q[0] + lo,
                        q[1] + lo,
                        q[2] + lo,
                        q[3] + lo,
                        q[0] + hi,
                        q[1] + hi,
                        q[2] + hi,
                        q[3] + hi,
                    ];
                    let faces = hex_faces(&hx);
                    let zone = u32::from(j < nrot);
                    soup.add_cell_in_zone(faces.iter().map(|f| &f[..]), zone);
                }
            }
        }
        let depth = self.depth;
        let rmid = 0.5 * (ri + ro);
        soup.finish(&ANNULUS_BOUNDARIES, move |b: BoundaryFace| {
            let c = b.centroid;
            let name = if b.normal[2].abs() > 0.5 {
                if c[2] < 0.5 * depth { "zmin" } else { "zmax" }
            } else if c[0].hypot(c[1]) < rmid {
                "inner"
            } else {
                "outer"
            };
            name.to_string()
        })
    }

    /// The rotor ring and the stator ring as two separate meshes: the
    /// interface faces of both are boundaries (`rotor_interface`,
    /// `stator_interface`), the pitch of the two rings is `n_theta` and
    /// `n_theta_stator`.
    fn generate_split(&self) -> Result<Mesh, String> {
        let [ntr, nr, nz, nrot] = self.intervals();
        let nts = self.stator_intervals();
        let s = self.radial_nodes();
        let tau = std::f64::consts::TAU;
        let (ri, ro) = (self.r_inner, self.r_outer);
        let radius = |j: usize| ri + s[j] * (ro - ri);
        // 2D nodes: the rotor ring (j = 0..=nrot), then the stator ring
        // (j = nrot..=nr)
        let nj_r = nrot + 1;
        let nj_s = nr - nrot + 1;
        let mut pts2: Vec<[f64; 2]> = Vec::with_capacity(ntr * nj_r + nts * nj_s);
        for a in 0..ntr {
            let th = tau * a as f64 / ntr as f64;
            let (c, sn) = (th.cos(), th.sin());
            for j in 0..=nrot {
                pts2.push([radius(j) * c, radius(j) * sn]);
            }
        }
        for a in 0..nts {
            let th = tau * a as f64 / nts as f64;
            let (c, sn) = (th.cos(), th.sin());
            for j in nrot..=nr {
                pts2.push([radius(j) * c, radius(j) * sn]);
            }
        }
        let node_r = |a: usize, j: usize| (j + nj_r * (a % ntr)) as u32;
        let off = (ntr * nj_r) as u32;
        let node_s = |a: usize, j: usize| off + (j - nrot + nj_s * (a % nts)) as u32;
        let np = pts2.len() as u32;
        let mut pts = Vec::with_capacity(pts2.len() * (nz + 1));
        for k in 0..=nz {
            let z = self.depth * k as f64 / nz as f64;
            pts.extend(pts2.iter().map(|p| [p[0], p[1], z]));
        }
        let mut soup = CellSoup::new(pts);
        soup.set_zone_names(vec!["stator".to_string(), "rotor".to_string()]);
        for k in 0..nz {
            let (lo, hi) = (k as u32 * np, (k as u32 + 1) * np);
            let mut add = |q: [u32; 4], zone: u32| {
                let hx = [
                    q[0] + lo,
                    q[1] + lo,
                    q[2] + lo,
                    q[3] + lo,
                    q[0] + hi,
                    q[1] + hi,
                    q[2] + hi,
                    q[3] + hi,
                ];
                let faces = hex_faces(&hx);
                soup.add_cell_in_zone(faces.iter().map(|f| &f[..]), zone);
            };
            for a in 0..ntr {
                for j in 0..nrot {
                    add(
                        [
                            node_r(a, j),
                            node_r(a, j + 1),
                            node_r(a + 1, j + 1),
                            node_r(a + 1, j),
                        ],
                        1,
                    );
                }
            }
            for a in 0..nts {
                for j in nrot..nr {
                    add(
                        [
                            node_s(a, j),
                            node_s(a, j + 1),
                            node_s(a + 1, j + 1),
                            node_s(a + 1, j),
                        ],
                        0,
                    );
                }
            }
        }
        let depth = self.depth;
        let r_if = radius(nrot);
        // the faces are planar chords: their centroids lie at r cos(pi / n)
        let chord = |r: f64, n: usize| r * (std::f64::consts::PI / n as f64).cos();
        let (c_rotor_if, c_outer) = (chord(r_if, ntr), chord(ro, nts));
        let (c_inner, c_stator_if) = (chord(ri, ntr), chord(r_if, nts));
        let mut names: Vec<&str> = ANNULUS_BOUNDARIES.to_vec();
        names.extend(INTERFACE_BOUNDARIES);
        soup.finish(&names, move |b: BoundaryFace| {
            let c = b.centroid;
            let r = c[0].hypot(c[1]);
            let name = if b.normal[2].abs() > 0.5 {
                if c[2] < 0.5 * depth { "zmin" } else { "zmax" }
            } else if (b.normal[0] * c[0] + b.normal[1] * c[1]) > 0.0 {
                // out of the domain, away from the axis: the interface
                // side of the rotor or the outer cylinder
                if (r - c_rotor_if).abs() < (r - c_outer).abs() {
                    "rotor_interface"
                } else {
                    "outer"
                }
            } else if (r - c_inner).abs() < (r - c_stator_if).abs() {
                "inner"
            } else {
                "stator_interface"
            };
            name.to_string()
        })
    }

    /// One line for the log.
    pub fn describe(&self) -> String {
        let [nt, nr, nz, nrot] = self.intervals();
        let rotor = if self.split() {
            format!(
                ", rotor zone: the inner {nrot} radial layers with {nt} intervals around the \
                 axis, stator {} (split interface)",
                self.stator_intervals()
            )
        } else if nrot > 0 {
            format!(", rotor zone: the inner {nrot} radial layers")
        } else {
            String::new()
        };
        format!(
            "generated: annulus {} to {} m, {} x {} x {} (around x radial x z), refine {}{}",
            self.r_inner, self.r_outer, nt, nr, nz, self.refine, rotor
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Geometry;

    #[test]
    fn volume_boundaries_and_zones() {
        let s = AnnulusSpec {
            n_theta: 24,
            n_r: 6,
            rotor_layers: 2,
            radial_ratio: 3.0,
            ..Default::default()
        };
        let m = s.generate().unwrap();
        assert_eq!(m.n_cells, s.estimated_cells());
        let g = Geometry::compute(&m);
        assert!(g.cell_volume.iter().all(|&v| v > 0.0));
        // the polygons inscribed in the circles have n_theta sides
        let n = 24.0_f64;
        let ring = 0.5 * n * (std::f64::consts::TAU / n).sin();
        let exact = ring * (0.10_f64.powi(2) - 0.05_f64.powi(2)) * 0.01;
        let vol: f64 = g.cell_volume.iter().sum();
        assert!((vol - exact).abs() < 1e-12 * exact, "{vol} {exact}");
        let names: Vec<&str> = m.boundaries.iter().map(|p| p.name.as_str()).collect();
        for p in ANNULUS_BOUNDARIES {
            assert!(names.contains(&p), "{p} missing in {names:?}");
        }
        let len = |name: &str| m.boundaries.iter().find(|p| p.name == name).unwrap().len;
        assert_eq!(len("inner"), 24);
        assert_eq!(len("outer"), 24);
        assert_eq!(len("zmin"), 24 * 6);
        assert_eq!(len("zmax"), 24 * 6);
        // zones: the inner two layers form the rotor
        assert_eq!(m.zones, vec!["stator".to_string(), "rotor".to_string()]);
        let rotor = m.cell_zone.iter().filter(|&&z| z == 1).count();
        assert_eq!(rotor, 24 * 2);
        let r_if = s.rotor_radius().unwrap();
        for c in 0..m.n_cells {
            let x = g.cell_centroid[c];
            let r = x[0].hypot(x[1]);
            assert_eq!(m.cell_zone[c] == 1, r < r_if, "cell {c} at r = {r}");
        }
    }

    #[test]
    fn a_split_interface_has_two_boundaries_and_unequal_pitch() {
        let s = AnnulusSpec {
            n_theta: 12,
            n_theta_stator: 18,
            n_r: 4,
            rotor_layers: 2,
            ..Default::default()
        };
        let m = s.generate().unwrap();
        assert_eq!(m.n_cells, s.estimated_cells());
        assert_eq!(m.n_cells, 12 * 2 + 18 * 2);
        let len = |name: &str| m.boundaries.iter().find(|p| p.name == name).unwrap().len;
        assert_eq!(len("rotor_interface"), 12);
        assert_eq!(len("stator_interface"), 18);
        assert_eq!(len("inner"), 12);
        assert_eq!(len("outer"), 18);
        let g = Geometry::compute(&m);
        assert!(g.cell_volume.iter().all(|&v| v > 0.0));
        // the faces of the two interface boundaries lie on one circle and
        // point away from the rotor and into it
        let r_if = s.rotor_radius().unwrap();
        for (name, sign) in [("rotor_interface", 1.0), ("stator_interface", -1.0)] {
            let p = m.boundaries.iter().find(|p| p.name == name).unwrap();
            for f in p.range() {
                for &n in m.face(f) {
                    let x = m.points[n as usize];
                    assert!((x[0].hypot(x[1]) - r_if).abs() < 1e-12);
                }
                let c = g.face_centroid[f];
                let (n, _) = g.normal_area(f);
                assert!(sign * (n[0] * c[0] + n[1] * c[1]) > 0.0, "{name} {f}");
            }
        }
        // every cell closes
        for c in 0..m.n_cells {
            let s = g.closure(&m, c);
            assert!(s.iter().all(|x| x.abs() < 1e-12), "cell {c}: {s:?}");
        }
        assert_eq!(m.zones, vec!["stator".to_string(), "rotor".to_string()]);
        // an equal pitch also splits the interface
        let e = AnnulusSpec {
            n_theta: 12,
            n_theta_stator: 12,
            n_r: 4,
            rotor_layers: 2,
            ..Default::default()
        };
        assert_eq!(
            e.generate()
                .unwrap()
                .boundaries
                .iter()
                .find(|p| p.name == "rotor_interface")
                .unwrap()
                .len,
            12
        );
        // refinement doubles both pitches
        let r = AnnulusSpec {
            refine: 1,
            ..s.clone()
        };
        assert_eq!(r.generate().unwrap().n_cells, (24 * 4 + 36 * 4));
        // bad specs
        let bad = AnnulusSpec {
            n_theta_stator: 18,
            rotor_layers: 0,
            ..Default::default()
        };
        assert!(bad.validate().unwrap_err().contains("rotor_layers"));
    }

    #[test]
    fn one_zone_without_rotor_and_layers_in_z() {
        let s = AnnulusSpec {
            n_theta: 12,
            n_r: 3,
            n_z: 4,
            ..Default::default()
        };
        let m = s.generate().unwrap();
        assert_eq!(m.n_cells, 12 * 3 * 4);
        assert_eq!(m.zones, vec!["fluid".to_string()]);
        assert_eq!(s.dimension(), 3);
        let g = Geometry::compute(&m);
        let n = 12.0_f64;
        let ring = 0.5 * n * (std::f64::consts::TAU / n).sin();
        let exact = ring * (0.10_f64.powi(2) - 0.05_f64.powi(2)) * 0.01;
        let vol: f64 = g.cell_volume.iter().sum();
        assert!((vol - exact).abs() < 1e-12 * exact);
    }

    #[test]
    fn the_rotor_interface_is_the_same_radius_on_every_level() {
        let s0 = AnnulusSpec {
            n_theta: 8,
            n_r: 4,
            rotor_layers: 1,
            radial_ratio: 4.0,
            ..Default::default()
        };
        let s2 = AnnulusSpec {
            refine: 2,
            ..s0.clone()
        };
        let (r0, r2) = (s0.rotor_radius().unwrap(), s2.rotor_radius().unwrap());
        assert!((r0 - r2).abs() < 1e-14, "{r0} {r2}");
        assert_eq!(s2.generate().unwrap().n_cells, 32 * 16);
    }

    #[test]
    fn bad_specs_are_refused() {
        let bad = |s: AnnulusSpec| s.validate().unwrap_err();
        assert!(
            bad(AnnulusSpec {
                r_outer: 0.04,
                ..Default::default()
            })
            .contains("r_inner")
        );
        assert!(
            bad(AnnulusSpec {
                rotor_layers: 16,
                ..Default::default()
            })
            .contains("stator")
        );
        assert!(
            bad(AnnulusSpec {
                n_theta: 2,
                ..Default::default()
            })
            .contains("n_theta")
        );
    }
}

//! Built-in mesh of a channel with a circular cylinder (D34): the geometry
//! of the cylinder benchmark of Schäfer & Turek (1996) by default, one cell
//! layer thick in z (a 2D case).
//!
//! Two blocks of hexahedra:
//!
//! - an O-grid around the cylinder that fills the rectangle
//!   `0 ≤ x ≤ height` × `0 ≤ y ≤ height` (the inlet is its left side, the
//!   channel walls its lower and upper sides): `4 n` intervals around the
//!   cylinder (`n` per side of the rectangle, the circle split at the
//!   diagonals), `radial` intervals from the cylinder to the rectangle,
//!   graded with `radial_ratio` (the last to the first interval). The node
//!   rings run from the circle to the rectangle in polar coordinates about
//!   the cylinder centre with the angle blended quadratically, so the grid
//!   lines leave the cylinder normal to it;
//! - the wake block from `x = height` to the outlet, `wake` intervals in x
//!   graded with `wake_ratio`, the `n` intervals of the O-grid side in y.
//!
//! `refine` halves every interval (the gradings follow the continuous
//! mapping of [`crate::generate::Grading`], so the refined meshes are
//! nested in the radial and wake directions).
//!
//! Boundaries: `inlet` (x = 0), `outlet` (x = length), `walls` (y = 0 and
//! y = height), `cylinder`, `zmin`, `zmax`.

use crate::builder::{BoundaryFace, CellSoup};
use crate::generate::{Grading, hex_faces};
use crate::mesh::Mesh;
use serde::{Deserialize, Serialize};

/// Boundary names of the cylinder mesh.
pub const CYLINDER_BOUNDARIES: [&str; 6] = ["inlet", "outlet", "walls", "cylinder", "zmin", "zmax"];

/// Specification of a channel with a cylinder (`[mesh.cylinder]`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct CylinderSpec {
    /// Channel length (m).
    pub length: f64,
    /// Channel height (m).
    pub height: f64,
    /// Cell layer thickness in z (m).
    pub depth: f64,
    /// Cylinder centre (x, y) (m).
    pub center: [f64; 2],
    /// Cylinder diameter (m).
    pub diameter: f64,
    /// Intervals per side of the O-grid (a quarter of the circle).
    pub n: usize,
    /// Intervals from the cylinder to the O-grid rectangle.
    pub radial: usize,
    /// Ratio of the last to the first radial interval.
    pub radial_ratio: f64,
    /// Intervals in x of the wake block.
    pub wake: usize,
    /// Ratio of the last to the first interval of the wake block.
    pub wake_ratio: f64,
    /// Refinement level: every interval is halved `refine` times.
    pub refine: u32,
}

impl Default for CylinderSpec {
    fn default() -> Self {
        Self {
            length: 2.2,
            height: 0.41,
            depth: 0.01,
            center: [0.2, 0.2],
            diameter: 0.1,
            n: 32,
            radial: 32,
            radial_ratio: 6.0,
            wake: 64,
            wake_ratio: 8.0,
            refine: 0,
        }
    }
}

impl CylinderSpec {
    /// Intervals (n, radial, wake) after refinement.
    pub fn intervals(&self) -> [usize; 3] {
        let f = 1usize << self.refine;
        [self.n * f, self.radial * f, self.wake * f]
    }

    /// Number of cells.
    pub fn estimated_cells(&self) -> usize {
        let [n, r, w] = self.intervals();
        4 * n * r + w * n
    }

    /// Checks the spec.
    pub fn validate(&self) -> Result<(), String> {
        let r = 0.5 * self.diameter;
        let [cx, cy] = self.center;
        let h = self.height;
        if !(self.length > h && h > 0.0 && self.depth > 0.0 && self.diameter > 0.0) {
            return Err(
                "[mesh.cylinder]: need length > height > 0, depth > 0 and diameter > 0".into(),
            );
        }
        if !(cx - r > 0.0 && cx + r < h && cy - r > 0.0 && cy + r < h) {
            return Err(format!(
                "[mesh.cylinder]: the cylinder must lie inside the O-grid square 0 < x, y < {h} \
                 (the height); centre {:?}, radius {r}",
                self.center
            ));
        }
        if self.n < 2 || self.radial < 1 || self.wake < 1 {
            return Err("[mesh.cylinder]: need n >= 2, radial >= 1 and wake >= 1".into());
        }
        if self.radial_ratio <= 0.0 || self.wake_ratio <= 0.0 {
            return Err("[mesh.cylinder]: the ratios must be > 0".into());
        }
        if self.refine > 6 {
            return Err("[mesh.cylinder]: refine must be <= 6".into());
        }
        Ok(())
    }

    /// Generates the mesh.
    pub fn generate(&self) -> Result<Mesh, String> {
        self.validate()?;
        let [n, nr, nw] = self.intervals();
        let r = 0.5 * self.diameter;
        let [cx, cy] = self.center;
        let (h, len) = (self.height, self.length);
        let na = 4 * n; // nodes around (periodic)
        let rad = Grading {
            ratio: self.radial_ratio,
            symmetric: false,
        }
        .nodes(nr, self.refine);
        let wake = Grading {
            ratio: self.wake_ratio,
            symmetric: false,
        }
        .nodes(nw, self.refine);
        // point on the rectangle at position a/n along the sides,
        // counter-clockwise from the lower right corner
        let side = |a: usize| -> [f64; 2] {
            let (s, t) = (a / n, (a % n) as f64 / n as f64);
            match s % 4 {
                0 => [h, t * h],
                1 => [h - t * h, h],
                2 => [0.0, h - t * h],
                _ => [t * h, 0.0],
            }
        };
        let tau = std::f64::consts::TAU;
        let mut pts2: Vec<[f64; 2]> = Vec::with_capacity(na * (nr + 1) + nw * (n + 1));
        // O-grid nodes, index j + (nr + 1) a
        for a in 0..na {
            let th_a = tau * (a as f64 / na as f64 - 0.125);
            let b = side(a);
            let (dx, dy) = (b[0] - cx, b[1] - cy);
            let rho_b = dx.hypot(dy);
            let mut th_b = dy.atan2(dx);
            while th_b - th_a > 0.5 * tau {
                th_b -= tau;
            }
            while th_a - th_b > 0.5 * tau {
                th_b += tau;
            }
            for (j, &s) in rad.iter().enumerate() {
                if j == nr {
                    pts2.push(b);
                    continue;
                }
                let rho = r + s * (rho_b - r);
                let th = th_a + s * s * (th_b - th_a);
                pts2.push([cx + rho * th.cos(), cy + rho * th.sin()]);
            }
        }
        let og = |a: usize, j: usize| (j + (nr + 1) * (a % na)) as u32;
        // wake nodes (columns 1..=nw), the first column is the right side
        // of the O-grid (a = 0..=n)
        let w0 = pts2.len();
        for (i, &s) in wake.iter().enumerate().skip(1) {
            for jy in 0..=n {
                let x = if i == nw { len } else { h + s * (len - h) };
                pts2.push([x, side(jy)[1].max(0.0)]);
            }
        }
        let wk = |i: usize, jy: usize| -> u32 {
            if i == 0 {
                og(jy, nr)
            } else {
                (w0 + (i - 1) * (n + 1) + jy) as u32
            }
        };
        // the top corner of the right side: a = n is on side 1 at t = 0
        // (the same point (h, h)); side(n) gives it exactly
        let np = pts2.len() as u32;
        let mut pts = Vec::with_capacity(2 * pts2.len());
        for z in [0.0, self.depth] {
            pts.extend(pts2.iter().map(|p| [p[0], p[1], z]));
        }
        let mut soup = CellSoup::new(pts);
        let hex = |soup: &mut CellSoup, q: [u32; 4]| {
            let hx = [
                q[0],
                q[1],
                q[2],
                q[3],
                q[0] + np,
                q[1] + np,
                q[2] + np,
                q[3] + np,
            ];
            let faces = hex_faces(&hx);
            soup.add_cell(faces.iter().map(|f| &f[..]));
        };
        // O-grid: (radial, around) is right-handed in the plane
        for a in 0..na {
            for j in 0..nr {
                hex(
                    &mut soup,
                    [og(a, j), og(a, j + 1), og(a + 1, j + 1), og(a + 1, j)],
                );
            }
        }
        for i in 0..nw {
            for jy in 0..n {
                hex(
                    &mut soup,
                    [wk(i, jy), wk(i + 1, jy), wk(i + 1, jy + 1), wk(i, jy + 1)],
                );
            }
        }
        let tol = 1e-9 * h;
        let (depth, diameter) = (self.depth, self.diameter);
        soup.finish(&CYLINDER_BOUNDARIES, move |b: BoundaryFace| {
            let c = b.centroid;
            let name = if b.normal[2].abs() > 0.5 {
                if c[2] < 0.5 * depth { "zmin" } else { "zmax" }
            } else if c[0] < tol {
                "inlet"
            } else if c[0] > len - tol {
                "outlet"
            } else if c[1] < tol || c[1] > h - tol {
                "walls"
            } else {
                debug_assert!(
                    ((c[0] - cx).hypot(c[1] - cy) - 0.5 * diameter).abs() < 0.05 * diameter,
                    "boundary face at {c:?} is not on the cylinder"
                );
                "cylinder"
            };
            name.to_string()
        })
    }

    /// One line for the log.
    pub fn describe(&self) -> String {
        let [n, nr, nw] = self.intervals();
        format!(
            "generated: channel {} x {} m with a cylinder D = {} m at {:?}, O-grid {} x {} \
             (around x radial), wake {} x {}, refine {}",
            self.length,
            self.height,
            self.diameter,
            self.center,
            4 * n,
            nr,
            nw,
            n,
            self.refine
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Geometry;

    #[test]
    fn volume_and_boundaries() {
        let s = CylinderSpec {
            n: 8,
            radial: 6,
            wake: 10,
            ..Default::default()
        };
        let m = s.generate().unwrap();
        assert_eq!(m.n_cells, s.estimated_cells());
        let g = Geometry::compute(&m);
        let vol: f64 = g.cell_volume.iter().sum();
        assert!(g.cell_volume.iter().all(|&v| v > 0.0));
        // the polygon inscribed in the circle has 4n sides
        let na = 32.0;
        let r = 0.05_f64;
        let hole = 0.5 * na * r * r * (std::f64::consts::TAU / na).sin();
        let exact = (2.2 * 0.41 - hole) * 0.01;
        assert!((vol - exact).abs() < 1e-12 * exact, "{vol} {exact}");
        let names: Vec<&str> = m.boundaries.iter().map(|p| p.name.as_str()).collect();
        for p in CYLINDER_BOUNDARIES {
            assert!(names.contains(&p), "{p} missing in {names:?}");
        }
        let cyl = m.boundaries.iter().find(|p| p.name == "cylinder").unwrap();
        assert_eq!(cyl.len, 32);
        let inlet = m.boundaries.iter().find(|p| p.name == "inlet").unwrap();
        assert_eq!(inlet.len, 8);
    }
}

//! Mesh quality report.
//!
//! Metrics (all from the geometry of [`crate::geometry`]):
//!
//! - **closure** of a cell: |Σ S_f| / Σ |S_f|, zero up to round-off for a
//!   closed cell (stage-1 acceptance test);
//! - **non-orthogonality** of a face: the angle between its area vector and
//!   the vector d from the owner centroid to the neighbour centroid (for a
//!   boundary face: to the face centroid), in degrees;
//! - **skewness** of an interior face: the distance between the face
//!   centroid and the point where the line through both cell centroids
//!   crosses the face plane, divided by |d|;
//! - **face warp**: the largest distance of a face node from the plane
//!   through the face centroid normal to its area vector, divided by the
//!   square root of the face area;
//! - **volume ratio** of an interior face: larger over smaller cell volume;
//! - **aspect ratio** of a cell: the largest over the smallest distance from
//!   the cell centroid to its face centroids (1 for a cube or a regular
//!   tetrahedron);
//! - **inverted** faces: S_f · d ≤ 0, so the neighbour centroid lies behind
//!   the face; and non-positive cell volumes. Both are errors.

use crate::geometry::{Geometry, add, dot, norm, scale, sub};
use crate::mesh::{CellShape, Mesh};
use rayon::prelude::*;
use serde::Serialize;
use std::fmt::Write as _;

/// Minimum, mean and maximum of a metric, with the item where the maximum
/// (or, for `min_at`, the minimum) occurs.
#[derive(Debug, Clone, Copy, Serialize, Default)]
pub struct Stat {
    /// Smallest value.
    pub min: f64,
    /// Mean value.
    pub mean: f64,
    /// Largest value.
    pub max: f64,
    /// Index of the cell or face with the largest value.
    pub max_at: usize,
    /// Index of the cell or face with the smallest value.
    pub min_at: usize,
}

impl Stat {
    fn of(v: &[f64]) -> Self {
        if v.is_empty() {
            return Self::default();
        }
        let mut s = Self {
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            ..Default::default()
        };
        for (i, &x) in v.iter().enumerate() {
            if x > s.max {
                s.max = x;
                s.max_at = i;
            }
            if x < s.min {
                s.min = x;
                s.min_at = i;
            }
        }
        s.mean = foehncraft_core::repro::pairwise_sum(v) / v.len() as f64;
        s
    }
}

/// Per-cell quality fields (for output).
#[derive(Debug, Clone, Default)]
pub struct CellQuality {
    /// |Σ S_f| / Σ |S_f| per cell.
    pub closure: Vec<f64>,
    /// Largest non-orthogonality of the faces of each cell, degrees.
    pub non_orthogonality: Vec<f64>,
    /// Largest skewness of the faces of each cell.
    pub skewness: Vec<f64>,
    /// Aspect ratio.
    pub aspect_ratio: Vec<f64>,
    /// Largest volume ratio to a neighbour.
    pub volume_ratio: Vec<f64>,
}

/// The mesh quality report.
#[derive(Debug, Clone, Serialize)]
pub struct QualityReport {
    /// Number of local cells.
    pub cells: usize,
    /// Number of faces.
    pub faces: usize,
    /// Number of interior faces.
    pub interior_faces: usize,
    /// Number of boundary faces.
    pub boundary_faces: usize,
    /// Number of nodes.
    pub points: usize,
    /// Cell count per shape.
    pub shapes: Vec<(String, usize)>,
    /// Faces per cell.
    pub faces_per_cell: Stat,
    /// Boundaries: name, number of faces, area.
    pub boundaries: Vec<(String, usize, f64)>,
    /// Bounding box.
    pub bounding_box: ([f64; 3], [f64; 3]),
    /// Sum of the cell volumes (pairwise sum, D31).
    pub total_volume: f64,
    /// Volume enclosed by the boundary faces, (1/3) Σ x_f · S_f over boundary
    /// faces; equals `total_volume` for a watertight mesh.
    pub boundary_volume: f64,
    /// Cell volumes.
    pub volume: Stat,
    /// Closure.
    pub closure: Stat,
    /// Non-orthogonality, degrees, over all faces.
    pub non_orthogonality: Stat,
    /// Skewness over interior faces.
    pub skewness: Stat,
    /// Face warp.
    pub warp: Stat,
    /// Aspect ratio.
    pub aspect_ratio: Stat,
    /// Volume ratio over interior faces.
    pub volume_ratio: Stat,
    /// Number of cells with a non-positive volume.
    pub negative_volumes: usize,
    /// Number of faces with S_f · d ≤ 0.
    pub inverted_faces: usize,
    /// Faces with non-orthogonality above 70°.
    pub non_orthogonal_above_70: usize,
    /// Interior faces with skewness above 0.9: on these the flow solver
    /// uses its cell-quality remediation (D50, `solver.skewness_limit`).
    pub skewness_above_0_9: usize,
}

impl QualityReport {
    /// Computes the report and the per-cell fields.
    pub fn compute(mesh: &Mesh, geo: &Geometry) -> (Self, CellQuality) {
        let nf = mesh.n_faces();
        let ni = mesh.n_interior_faces();
        let nc = mesh.n_cells;
        let n_owned_face = |f: usize| f < ni && (mesh.neighbour[f] as usize) < nc + mesh.n_halo;
        // periodic faces: the neighbour centroid moved next to the face
        let periodic_shift = |f: usize| {
            mesh.periodic
                .as_ref()
                .filter(|p| f >= p.start && f < ni)
                .map(|p| p.shift[f - p.start])
        };

        // per-face metrics
        let face: Vec<(f64, f64, f64, f64, bool)> = (0..nf)
            .into_par_iter()
            .map(|f| {
                let p = mesh.owner[f] as usize;
                let xp = geo.cell_centroid[p];
                let s = geo.face_area[f];
                let fc = geo.face_centroid[f];
                let d = if let Some(s) = periodic_shift(f) {
                    sub(add(geo.cell_centroid[mesh.neighbour[f] as usize], s), xp)
                } else if n_owned_face(f) {
                    sub(geo.cell_centroid[mesh.neighbour[f] as usize], xp)
                } else {
                    sub(fc, xp)
                };
                let sd = dot(s, d);
                let cosang = (sd / (norm(s) * norm(d))).clamp(-1.0, 1.0);
                let nonorth = cosang.acos().to_degrees();
                let (skew, vr) = if f < ni {
                    // intersection of the centroid line with the face plane
                    let t = dot(sub(fc, xp), s) / sd;
                    let xi = [xp[0] + t * d[0], xp[1] + t * d[1], xp[2] + t * d[2]];
                    let vp = geo.cell_volume[p];
                    let vn = geo.cell_volume[mesh.neighbour[f] as usize];
                    let vr = if vp > 0.0 && vn > 0.0 {
                        vp.max(vn) / vp.min(vn)
                    } else {
                        f64::INFINITY
                    };
                    (norm(sub(fc, xi)) / norm(d), vr)
                } else {
                    (0.0, 1.0)
                };
                let area = norm(s);
                let n = scale(s, 1.0 / area);
                let warp = mesh
                    .face(f)
                    .iter()
                    .map(|&k| dot(sub(mesh.points[k as usize], fc), n).abs())
                    .fold(0.0, f64::max)
                    / area.sqrt();
                (nonorth, skew, warp, vr, sd <= 0.0)
            })
            .collect();

        let mut cq = CellQuality {
            closure: vec![0.0; nc],
            non_orthogonality: vec![0.0; nc],
            skewness: vec![0.0; nc],
            aspect_ratio: vec![0.0; nc],
            volume_ratio: vec![1.0; nc],
        };
        let mut faces_per_cell = vec![0.0; nc];
        (0..nc)
            .into_par_iter()
            .zip(cq.closure.par_iter_mut())
            .zip(cq.non_orthogonality.par_iter_mut())
            .zip(cq.skewness.par_iter_mut())
            .zip(cq.aspect_ratio.par_iter_mut())
            .zip(cq.volume_ratio.par_iter_mut())
            .zip(faces_per_cell.par_iter_mut())
            .for_each(|((((((c, clo), no), sk), ar), vr), fpc)| {
                let s = geo.closure(mesh, c);
                let mut sum_area = 0.0;
                let (mut dmin, mut dmax) = (f64::INFINITY, 0.0f64);
                let x = geo.cell_centroid[c];
                for (f, _) in mesh.cell_face_iter(c) {
                    sum_area += norm(geo.face_area[f]);
                    *no = no.max(face[f].0);
                    *sk = sk.max(face[f].1);
                    *vr = vr.max(face[f].3);
                    let d = norm(sub(geo.face_centroid[f], x));
                    dmin = dmin.min(d);
                    dmax = dmax.max(d);
                    *fpc += 1.0;
                }
                *clo = norm(s) / sum_area;
                *ar = dmax / dmin;
            });

        let mut shapes = vec![0usize; CellShape::ALL.len()];
        for c in 0..nc {
            let s = mesh.cell_shape(c);
            shapes[CellShape::ALL.iter().position(|&x| x == s).unwrap()] += 1;
        }
        let shapes = CellShape::ALL
            .iter()
            .zip(shapes)
            .filter(|(_, n)| *n > 0)
            .map(|(s, n)| (s.name().to_string(), n))
            .collect();

        let boundaries = mesh
            .boundaries
            .iter()
            .map(|p| {
                let a: Vec<f64> = p.range().map(|f| norm(geo.face_area[f])).collect();
                (
                    p.name.clone(),
                    p.len,
                    foehncraft_core::repro::pairwise_sum(&a),
                )
            })
            .collect();
        let bvol: Vec<f64> = (ni..nf)
            .map(|f| dot(geo.face_centroid[f], geo.face_area[f]) / 3.0)
            .collect();
        let vols = &geo.cell_volume[..nc];
        let nonorth: Vec<f64> = face.iter().map(|x| x.0).collect();
        let skew: Vec<f64> = face[..ni].iter().map(|x| x.1).collect();
        let warp: Vec<f64> = face.iter().map(|x| x.2).collect();
        let vr: Vec<f64> = face[..ni].iter().map(|x| x.3).collect();
        let report = Self {
            cells: nc,
            faces: nf,
            interior_faces: ni,
            boundary_faces: nf - ni,
            points: mesh.n_points(),
            shapes,
            faces_per_cell: Stat::of(&faces_per_cell),
            boundaries,
            bounding_box: mesh.bounding_box(),
            total_volume: foehncraft_core::repro::pairwise_sum(vols),
            boundary_volume: foehncraft_core::repro::pairwise_sum(&bvol),
            volume: Stat::of(vols),
            closure: Stat::of(&cq.closure),
            non_orthogonality: Stat::of(&nonorth),
            skewness: Stat::of(&skew),
            warp: Stat::of(&warp),
            aspect_ratio: Stat::of(&cq.aspect_ratio),
            volume_ratio: Stat::of(&vr),
            negative_volumes: vols.iter().filter(|&&v| v <= 0.0).count(),
            inverted_faces: face.iter().filter(|x| x.4).count(),
            non_orthogonal_above_70: nonorth.iter().filter(|&&a| a > 70.0).count(),
            skewness_above_0_9: skew.iter().filter(|&&a| a > 0.9).count(),
        };
        (report, cq)
    }

    /// `true` if the mesh has no errors (no inverted faces, no non-positive
    /// volumes, every cell closed to `1e-10`).
    pub fn ok(&self) -> bool {
        self.negative_volumes == 0 && self.inverted_faces == 0 && self.closure.max < 1e-10
    }

    /// The report as Markdown.
    pub fn to_markdown(&self, title: &str) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "# Mesh quality: {title}\n");
        let _ = writeln!(
            s,
            "Verdict: **{}**\n",
            if self.ok() {
                "valid"
            } else {
                "INVALID (see errors below)"
            }
        );
        let _ = writeln!(s, "| | |\n|---|---|");
        let _ = writeln!(s, "| cells | {} |", self.cells);
        for (name, n) in &self.shapes {
            let _ = writeln!(s, "| &nbsp;&nbsp;{name} | {n} |");
        }
        let _ = writeln!(
            s,
            "| faces | {} ({} interior, {} boundary) |",
            self.faces, self.interior_faces, self.boundary_faces
        );
        let _ = writeln!(s, "| nodes | {} |", self.points);
        let _ = writeln!(
            s,
            "| faces per cell | {:.0} to {:.0}, mean {:.2} |",
            self.faces_per_cell.min, self.faces_per_cell.max, self.faces_per_cell.mean
        );
        let (lo, hi) = self.bounding_box;
        let _ = writeln!(
            s,
            "| bounding box | ({:.6}, {:.6}, {:.6}) to ({:.6}, {:.6}, {:.6}) |",
            lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]
        );
        let _ = writeln!(s, "| total volume | {:.12e} |", self.total_volume);
        let _ = writeln!(
            s,
            "| volume enclosed by the boundary | {:.12e} (relative difference {:.1e}) |",
            self.boundary_volume,
            ((self.total_volume - self.boundary_volume) / self.total_volume).abs()
        );
        let _ = writeln!(
            s,
            "\n## Boundaries\n\n| boundary | faces | area |\n|---|---:|---:|"
        );
        for (n, k, a) in &self.boundaries {
            let _ = writeln!(s, "| {n} | {k} | {a:.6e} |");
        }
        let _ = writeln!(
            s,
            "\n## Metrics\n\n| metric | min | mean | max | at |\n|---|---:|---:|---:|---:|"
        );
        let row = |s: &mut String, name: &str, st: &Stat, what: &str| {
            let _ = writeln!(
                s,
                "| {name} | {:.4e} | {:.4e} | {:.4e} | {what} {} |",
                st.min, st.mean, st.max, st.max_at
            );
        };
        row(&mut s, "cell volume", &self.volume, "cell");
        row(
            &mut s,
            "closure &#124;ΣS&#124;/Σ&#124;S&#124;",
            &self.closure,
            "cell",
        );
        row(
            &mut s,
            "non-orthogonality (deg)",
            &self.non_orthogonality,
            "face",
        );
        row(&mut s, "skewness", &self.skewness, "face");
        row(&mut s, "face warp", &self.warp, "face");
        row(&mut s, "aspect ratio", &self.aspect_ratio, "cell");
        row(&mut s, "volume ratio", &self.volume_ratio, "face");
        let _ = writeln!(s, "\n## Errors and warnings\n");
        let _ = writeln!(
            s,
            "- cells with non-positive volume: {}",
            self.negative_volumes
        );
        let _ = writeln!(s, "- inverted faces (S·d ≤ 0): {}", self.inverted_faces);
        let _ = writeln!(
            s,
            "- faces with non-orthogonality above 70°: {}",
            self.non_orthogonal_above_70
        );
        let _ = writeln!(
            s,
            "- faces with skewness above 0.9: {}{}",
            self.skewness_above_0_9,
            if self.skewness_above_0_9 > 0 {
                " (the flow solver treats them with its cell-quality remediation, D50; \
                 an optimised mesh avoids them, Gmsh: `Mesh.OptimizeNetgen = 1`)"
            } else {
                ""
            }
        );
        s
    }
}

//! Mesh geometry: face area vectors and centroids, cell volumes and
//! centroids.
//!
//! Faces may be non-planar polygons. A face is split into triangles fanning
//! out from the average of its nodes; its area vector is the sum of the
//! triangle area vectors, which depends only on the boundary loop of the
//! face, so the area vectors of a closed cell sum to zero up to round-off
//! for any face shape. The face centroid is the average of the triangle
//! centroids weighted by the triangle areas projected on the face normal.
//!
//! A cell is split into tetrahedra spanned by an interior point (the average
//! of its face centroids), the fan point of a face (its node average) and one
//! edge of the face, so cell volumes use the same face surfaces as the area
//! vectors and neighbouring cells tile the domain exactly. Volume and centroid are the sums over these tetrahedra (Gauss
//! divergence theorem). Geometry is always computed in float64, whatever the
//! floating-point type of the solver (D7).
//!
//! The same decomposition is used by point location ([`crate::locate`]).

use crate::mesh::Mesh;
use rayon::prelude::*;

/// Vector sum.
#[inline(always)]
pub fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Vector difference.
#[inline(always)]
pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scaled vector.
#[inline(always)]
pub fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
#[inline(always)]
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
#[inline(always)]
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm.
#[inline(always)]
pub fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Average of the nodes of a polygon.
#[inline]
pub fn node_average(points: &[[f64; 3]], nodes: &[u32]) -> [f64; 3] {
    let mut c = [0.0; 3];
    for &n in nodes {
        c = add(c, points[n as usize]);
    }
    scale(c, 1.0 / nodes.len() as f64)
}

/// Area vector and centroid of a (possibly non-planar) polygon.
pub fn polygon_area_centroid(points: &[[f64; 3]], nodes: &[u32]) -> ([f64; 3], [f64; 3]) {
    let k = nodes.len();
    if k == 3 {
        let (a, b, c) = (
            points[nodes[0] as usize],
            points[nodes[1] as usize],
            points[nodes[2] as usize],
        );
        let s = scale(cross(sub(b, a), sub(c, a)), 0.5);
        return (s, scale(add(add(a, b), c), 1.0 / 3.0));
    }
    let m = node_average(points, nodes);
    let mut s = [0.0; 3];
    let mut tri_s = [[0.0; 3]; 32];
    let mut big: Vec<[f64; 3]> = Vec::new();
    for i in 0..k {
        let p = points[nodes[i] as usize];
        let q = points[nodes[(i + 1) % k] as usize];
        let t = scale(cross(sub(p, m), sub(q, m)), 0.5);
        s = add(s, t);
        if k <= 32 {
            tri_s[i] = t;
        } else {
            big.push(t);
        }
    }
    let ns = norm(s);
    if ns == 0.0 {
        return (s, m);
    }
    let n = scale(s, 1.0 / ns);
    let mut c = [0.0; 3];
    let mut w = 0.0;
    for i in 0..k {
        let p = points[nodes[i] as usize];
        let q = points[nodes[(i + 1) % k] as usize];
        let t = if k <= 32 { tri_s[i] } else { big[i] };
        let a = dot(t, n);
        c = add(c, scale(add(add(p, q), m), a / 3.0));
        w += a;
    }
    (s, scale(c, 1.0 / w))
}

/// Geometry of a mesh, for the local cells and faces.
#[derive(Debug, Clone, Default)]
pub struct Geometry {
    /// Area vector of each face (normal out of the owner, length = area).
    pub face_area: Vec<[f64; 3]>,
    /// Centroid of each face.
    pub face_centroid: Vec<[f64; 3]>,
    /// Volume of each cell (local cells first; halo cells, if any, after
    /// them, filled by the halo exchange or by the decomposition).
    pub cell_volume: Vec<f64>,
    /// Centroid of each cell (local and halo).
    pub cell_centroid: Vec<[f64; 3]>,
}

impl Geometry {
    /// Computes the geometry of all faces and local cells. Halo entries are
    /// left at zero.
    pub fn compute(mesh: &Mesh) -> Self {
        let (face_area, face_centroid): (Vec<_>, Vec<_>) = (0..mesh.n_faces())
            .into_par_iter()
            .map(|f| polygon_area_centroid(&mesh.points, mesh.face(f)))
            .unzip();
        let mut cell_volume = vec![0.0; mesh.n_cells + mesh.n_halo];
        let mut cell_centroid = vec![[0.0; 3]; mesh.n_cells + mesh.n_halo];
        cell_volume[..mesh.n_cells]
            .par_iter_mut()
            .zip(cell_centroid[..mesh.n_cells].par_iter_mut())
            .enumerate()
            .for_each(|(c, (v, x))| {
                let (vol, cen) = cell_volume_centroid(mesh, &face_centroid, c);
                *v = vol;
                *x = cen;
            });
        Self {
            face_area,
            face_centroid,
            cell_volume,
            cell_centroid,
        }
    }

    /// Sum of the outward area vectors of local cell `c` (zero for a closed
    /// cell).
    pub fn closure(&self, mesh: &Mesh, c: usize) -> [f64; 3] {
        let mut s = [0.0; 3];
        for (f, is_nb) in mesh.cell_face_iter(c) {
            let a = self.face_area[f];
            s = if is_nb { sub(s, a) } else { add(s, a) };
        }
        s
    }

    /// Unit normal and area of face `f`.
    #[inline]
    pub fn normal_area(&self, f: usize) -> ([f64; 3], f64) {
        let a = norm(self.face_area[f]);
        (scale(self.face_area[f], 1.0 / a), a)
    }
}

/// Interior point of a local cell used for its tetrahedral decomposition:
/// the average of its face centroids.
pub fn cell_apex(mesh: &Mesh, face_centroid: &[[f64; 3]], c: usize) -> [f64; 3] {
    let mut x = [0.0; 3];
    let mut k = 0.0;
    for (f, _) in mesh.cell_face_iter(c) {
        x = add(x, face_centroid[f]);
        k += 1.0;
    }
    scale(x, 1.0 / k)
}

/// Volume and centroid of local cell `c`.
pub fn cell_volume_centroid(mesh: &Mesh, face_centroid: &[[f64; 3]], c: usize) -> (f64, [f64; 3]) {
    let apex = cell_apex(mesh, face_centroid, c);
    let mut vol = 0.0;
    let mut cen = [0.0; 3];
    for (f, is_nb) in mesh.cell_face_iter(c) {
        let nodes = mesh.face(f);
        let k = nodes.len();
        // fan about the node average, the same surface as the area vector
        let fc = if k == 3 {
            face_centroid[f]
        } else {
            node_average(&mesh.points, nodes)
        };
        for i in 0..k {
            let p = mesh.points[nodes[i] as usize];
            let q = mesh.points[nodes[(i + 1) % k] as usize];
            // tetrahedron (apex, fc, p, q), outward orientation of (fc, p, q)
            let mut v = dot(sub(fc, apex), cross(sub(p, apex), sub(q, apex))) / 6.0;
            if is_nb {
                v = -v;
            }
            vol += v;
            let tc = scale(add(add(apex, fc), add(p, q)), 0.25);
            cen = add(cen, scale(tc, v));
        }
    }
    (vol, scale(cen, 1.0 / vol))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_square_and_warped_quad() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.2],
        ];
        let (s, c) = polygon_area_centroid(&pts, &[0, 1, 2, 3]);
        assert!((s[2] - 1.0).abs() < 1e-15 && s[0] == 0.0 && s[1] == 0.0);
        assert!((c[0] - 0.5).abs() < 1e-15 && (c[1] - 0.5).abs() < 1e-15);
        // a warped quad has the same projected area as its loop
        let (s, _) = polygon_area_centroid(&pts, &[0, 1, 4, 3]);
        assert!((s[2] - 1.0).abs() < 1e-15);
        // triangle centroid
        let (s, c) = polygon_area_centroid(&pts, &[0, 1, 3]);
        assert!((s[2] - 0.5).abs() < 1e-15);
        assert!((c[0] - 1.0 / 3.0).abs() < 1e-15);
    }
}

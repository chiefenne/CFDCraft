//! Translational periodic boundary pairs (D29, D109 (1)).
//!
//! [`pair`] matches the faces of two boundaries that are translated copies
//! of each other, x_b = x_a + T, and turns them into the periodic faces of
//! [`crate::mesh::Periodic`]: each side keeps its own face, owned by its
//! cell, with the cell of the partner face as neighbour and the shift ∓T
//! that places that cell next to the face. The faces are matched by their
//! centroids (moved by T) and checked by their area vectors (opposite), so
//! the two boundaries must carry conformal copies of the same face set, as
//! a periodic mesh from the built-in generator or from Gmsh's
//! `Periodic Surface` has.
//!
//! The translation is the difference of the area-weighted centroids of the
//! two boundaries unless the case gives it. Pairing is the last step of
//! the mesh set-up, after renumbering, on a single partition.

use crate::geometry::{add, dot, norm, polygon_area_centroid, sub};
use crate::mesh::{Boundary, Mesh, Periodic, PeriodicPair};
use std::collections::HashMap;

/// One requested pair.
#[derive(Debug, Clone, PartialEq)]
pub struct PairSpec {
    /// First boundary.
    pub a: String,
    /// Partner boundary.
    pub b: String,
    /// Translation from `a` to `b` (m); `None`: from the face centroids.
    pub translation: Option<[f64; 3]>,
}

/// Relative tolerance of the face matching (to the face size √|S|).
pub const MATCH_TOLERANCE: f64 = 1e-6;

struct Side {
    faces: Vec<usize>,
    centroid: Vec<[f64; 3]>,
    area: Vec<[f64; 3]>,
}

fn side(mesh: &Mesh, b: &Boundary) -> Side {
    let faces: Vec<usize> = b.range().collect();
    let (mut area, mut centroid) = (Vec::new(), Vec::new());
    for &f in &faces {
        let (s, c) = polygon_area_centroid(&mesh.points, mesh.face(f));
        area.push(s);
        centroid.push(c);
    }
    Side {
        faces,
        centroid,
        area,
    }
}

fn mean_centroid(s: &Side) -> [f64; 3] {
    let mut c = [0.0; 3];
    let mut w = 0.0;
    for (x, a) in s.centroid.iter().zip(&s.area) {
        let m = norm(*a);
        for d in 0..3 {
            c[d] += m * x[d];
        }
        w += m;
    }
    c.map(|v| v / w)
}

/// The mesh with the faces of the requested pairs made periodic.
pub fn pair(mesh: &Mesh, specs: &[PairSpec]) -> Result<Mesh, String> {
    if specs.is_empty() {
        return Ok(mesh.clone());
    }
    if mesh.periodic.is_some() {
        return Err("the mesh is paired already".into());
    }
    if mesh.n_halo != 0 || mesh.n_inner_faces != mesh.n_interior_faces() {
        return Err("periodic pairs need the whole mesh on one partition".into());
    }
    let ni = mesh.n_interior_faces();
    let mut used: Vec<&str> = Vec::new();
    let mut matched: Vec<(Vec<usize>, Vec<usize>, [f64; 3])> = Vec::new();
    for sp in specs {
        for n in [&sp.a, &sp.b] {
            if used.contains(&n.as_str()) {
                return Err(format!("boundary `{n}` is in more than one periodic pair"));
            }
            used.push(n);
        }
        if sp.a == sp.b {
            return Err(format!("boundary `{}` cannot be its own periodic partner", sp.a));
        }
        let ba = mesh
            .boundary(&sp.a)
            .ok_or_else(|| format!("no boundary `{}` for the periodic pair", sp.a))?;
        let bb = mesh
            .boundary(&sp.b)
            .ok_or_else(|| format!("no boundary `{}` for the periodic pair", sp.b))?;
        let (sa, sb) = (side(mesh, ba), side(mesh, bb));
        if sa.faces.len() != sb.faces.len() {
            return Err(format!(
                "periodic pair `{}`–`{}`: {} and {} faces; the two boundaries must carry \
                 conformal copies of the same faces",
                sp.a,
                sp.b,
                sa.faces.len(),
                sb.faces.len()
            ));
        }
        if sa.faces.is_empty() {
            return Err(format!("periodic pair `{}`–`{}` has no faces", sp.a, sp.b));
        }
        let t = sp
            .translation
            .unwrap_or_else(|| sub(mean_centroid(&sb), mean_centroid(&sa)));
        // a hash grid on the centroids of side b, cell size the smallest
        // face size
        let h = sb
            .area
            .iter()
            .map(|a| norm(*a).sqrt())
            .fold(f64::INFINITY, f64::min);
        if h <= 0.0 || !h.is_finite() {
            return Err(format!("periodic pair `{}`–`{}`: degenerate faces", sp.a, sp.b));
        }
        let key = |x: [f64; 3]| x.map(|v| (v / h).floor() as i64);
        let mut grid: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        for (i, &x) in sb.centroid.iter().enumerate() {
            grid.entry(key(x)).or_default().push(i);
        }
        let mut taken = vec![false; sb.faces.len()];
        let mut bmatch = Vec::with_capacity(sa.faces.len());
        for i in 0..sa.faces.len() {
            let x = add(sa.centroid[i], t);
            let k = key(x);
            let size = norm(sa.area[i]).sqrt();
            let mut best: Option<(f64, usize)> = None;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(v) = grid.get(&[k[0] + dx, k[1] + dy, k[2] + dz]) {
                            for &j in v {
                                let e = norm(sub(sb.centroid[j], x));
                                if best.is_none_or(|b| e < b.0) {
                                    best = Some((e, j));
                                }
                            }
                        }
                    }
                }
            }
            let Some((e, j)) = best.filter(|b| b.0 <= MATCH_TOLERANCE * size) else {
                return Err(format!(
                    "periodic pair `{}`–`{}`: face {} at {:?} has no partner at {:?} \
                     (translation {:?}); the faces must be conformal translated copies",
                    sp.a, sp.b, sa.faces[i], sa.centroid[i], x, t
                ));
            };
            let _ = e;
            if taken[j] {
                return Err(format!(
                    "periodic pair `{}`–`{}`: two faces meet the same partner",
                    sp.a, sp.b
                ));
            }
            let (s1, s2) = (sa.area[i], sb.area[j]);
            if norm(add(s1, s2)) > 1e-6 * norm(s1) {
                return Err(format!(
                    "periodic pair `{}`–`{}`: faces {} and {} differ in area or orientation",
                    sp.a, sp.b, sa.faces[i], sb.faces[j]
                ));
            }
            taken[j] = true;
            bmatch.push(sb.faces[j]);
        }
        if dot(t, t) == 0.0 {
            return Err(format!("periodic pair `{}`–`{}`: zero translation", sp.a, sp.b));
        }
        matched.push((sa.faces, bmatch, t));
    }

    // the new face order: inner faces, the periodic faces pair by pair
    // (side a, then side b in the order of their partners), the other
    // boundaries
    let mut order: Vec<usize> = (0..ni).collect();
    let mut pairs = Vec::new();
    let mut shift = Vec::new();
    let mut pair_of = Vec::new();
    let mut partner = Vec::new();
    for (k, ((af, bf, t), sp)) in matched.iter().zip(specs).enumerate() {
        let a0 = order.len();
        order.extend(af);
        let b0 = order.len();
        order.extend(bf);
        let n = af.len();
        for i in 0..n {
            shift.push(t.map(|v| -v));
            pair_of.push(k as u32);
            partner.push((b0 + i) as u32);
        }
        for i in 0..n {
            shift.push(*t);
            pair_of.push(k as u32);
            partner.push((a0 + i) as u32);
        }
        pairs.push(PeriodicPair {
            a: sp.a.clone(),
            b: sp.b.clone(),
            translation: *t,
            a_faces: a0..b0,
            b_faces: b0..b0 + n,
        });
    }
    let n_interior = order.len();
    let mut boundaries = Vec::new();
    for b in &mesh.boundaries {
        if used.contains(&b.name.as_str()) {
            continue;
        }
        let start = order.len();
        order.extend(b.range());
        boundaries.push(Boundary {
            name: b.name.clone(),
            start,
            len: b.len,
        });
    }
    debug_assert_eq!(order.len(), mesh.n_faces());

    let mut m = Mesh {
        points: mesh.points.clone(),
        n_inner_faces: ni,
        n_cells: mesh.n_cells,
        n_inner_cells: mesh.n_inner_cells,
        n_halo: 0,
        boundaries,
        cell_zone: mesh.cell_zone.clone(),
        zones: mesh.zones.clone(),
        cell_id: mesh.cell_id.clone(),
        point_id: mesh.point_id.clone(),
        ..Default::default()
    };
    m.face_ptr.push(0);
    for &f in &order {
        m.face_nodes.extend_from_slice(mesh.face(f));
        m.face_ptr.push(m.face_nodes.len() as u32);
        m.owner.push(mesh.owner[f]);
        m.face_id.push(mesh.face_id[f]);
    }
    m.neighbour = mesh.neighbour[..ni].to_vec();
    for i in 0..n_interior - ni {
        let q = partner[i] as usize;
        let n = m.owner[q];
        if n == m.owner[ni + i] {
            let pr = &pairs[pair_of[i] as usize];
            return Err(format!(
                "periodic pair `{}`–`{}`: a cell touches both sides; the mesh needs at least \
                 two cells across the periodic direction",
                pr.a, pr.b
            ));
        }
        m.neighbour.push(n);
    }
    m.periodic = Some(Periodic {
        pairs,
        start: ni,
        shift,
        partner,
        pair: pair_of,
    });
    m.build_cell_faces();
    m.validate()?;
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Geometry;
    use crate::generate::{BoxSpec, Kind};

    fn spec(a: &str, b: &str) -> PairSpec {
        PairSpec {
            a: a.into(),
            b: b.into(),
            translation: None,
        }
    }

    #[test]
    fn pairs_the_box_in_x_and_z() {
        for kind in [Kind::Hex, Kind::Poly, Kind::Tet] {
            let m = BoxSpec::new(kind, 3, [4, 3, 5]).generate().unwrap();
            let names: Vec<String> = m.boundaries.iter().map(|b| b.name.clone()).collect();
            assert!(names.contains(&"xmin".to_string()), "{names:?}");
            let p = pair(&m, &[spec("xmin", "xmax"), spec("zmin", "zmax")]).unwrap();
            p.validate().unwrap();
            let pd = p.periodic.as_ref().unwrap();
            assert_eq!(pd.pairs.len(), 2);
            assert!((pd.pairs[0].translation[0] - 1.0).abs() < 1e-12);
            assert!(p.boundary("xmin").is_none() && p.boundary("ymin").is_some());
            // geometry unchanged (to round-off: the faces of a cell are
            // summed in another order)
            let (g0, g1) = (Geometry::compute(&m), Geometry::compute(&p));
            for c in 0..m.n_cells {
                assert!((g0.cell_volume[c] - g1.cell_volume[c]).abs() < 1e-15);
                assert!(norm(sub(g0.cell_centroid[c], g1.cell_centroid[c])) < 1e-14);
            }
            // the shifted neighbour lies on the other side of the face
            for f in pd.start..p.n_interior_faces() {
                let i = f - pd.start;
                let xn = add(g1.cell_centroid[p.neighbour[f] as usize], pd.shift[i]);
                let xp = g1.cell_centroid[p.owner[f] as usize];
                assert!(dot(sub(xn, xp), g1.face_area[f]) > 0.0);
                assert_eq!(pd.partner[pd.partner[i] as usize - pd.start] as usize, f);
            }
        }
    }

    #[test]
    fn refuses_bad_pairs() {
        let m = BoxSpec::new(Kind::Hex, 2, [3, 3, 1]).generate().unwrap();
        assert!(pair(&m, &[spec("xmin", "ymax")]).is_err());
        assert!(pair(&m, &[spec("xmin", "xmin")]).is_err());
        let one = BoxSpec::new(Kind::Hex, 2, [1, 3, 1]).generate().unwrap();
        let e = pair(&one, &[spec("xmin", "xmax")]).unwrap_err();
        assert!(e.contains("two cells"), "{e}");
    }
}

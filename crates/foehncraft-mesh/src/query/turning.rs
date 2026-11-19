//! The exact wall distance (D48) of a mesh with a part that turns rigidly
//! about an axis, for every angle of the turn without a rebuild: the sliding
//! mesh by imprinting (D182 variant B), whose rotor walls move past the
//! walls of the stator.
//!
//! **One hierarchy per rigid part.** The walls of the rotor (at the angle 0
//! of the base mesh) and the walls of the stator each get a [`Bvh`], built
//! once. A query point of the other part is turned into the frame of the
//! part (the instancing of rigid objects of the geometry-query libraries of
//! computer graphics; one hierarchy per object with its transform, as in
//! Sawhney's FCPW), so the hierarchies are never rebuilt or refitted.
//!
//! **Culling.** The distance of a cell to the walls of its own part does
//! not change with the angle and is kept. A wall of the other part can only
//! be nearer than that where the boundary of the own part is: a straight
//! line from a cell to any point outside its part leaves the part through
//! one of the part's boundary faces, either a wall of the part (then that
//! wall is nearer) or one of its other, **open** faces (the interface,
//! symmetry planes, inlets, outlets). A cell whose own-wall distance is not
//! larger than its distance to the open faces of its part keeps it, at
//! every angle. Only the other cells search the other part's hierarchy,
//! with a **capped** query (the own-wall distance as the largest radius:
//! boxes farther away are never opened) that starts from the face that was
//! the nearest at the previous angle (temporal coherence of a turn of about
//! a cell per step).
//!
//! The distances are those of [`super::wall_distance`] on the mesh at that
//! angle to round-off (a distance to a turned triangle is the distance of
//! the point turned back), and every cell is an independent query: any
//! number of threads gives the same result (D93).

use super::bvh::{Bvh, SurfaceFaces, closest_on_triangle};
use crate::imprint::rotate;
use crate::{Geometry, Mesh};
use rayon::prelude::*;

/// The wall distance of a mesh with rigidly turning parts (see the module
/// documentation).
#[derive(Debug, Clone)]
pub struct TurningWallDistance {
    /// Axis of every rotor: a point on it and the unit direction.
    turns: Vec<([f64; 3], [f64; 3])>,
    /// The part of every cell: rotor i < n, the stator n.
    part: Vec<u32>,
    /// The walls of every part (the rotors at the angle 0).
    walls: Vec<Option<Bvh>>,
    /// The triangles of every face, per part (start, count).
    face_tris: Vec<Vec<(u32, u32)>>,
    /// The distance of every cell to the walls of its own part.
    own: Vec<f64>,
    /// The cells that can see a wall of another part nearer than their
    /// own, and per part the face that was the nearest at the last angles
    /// (`u32::MAX`: none yet).
    candidates: Vec<(u32, Vec<u32>)>,
}

fn face_tris(s: &SurfaceFaces) -> Vec<(u32, u32)> {
    let mut r = vec![(0u32, 0u32); s.len()];
    for (t, &f) in s.tri_face.iter().enumerate() {
        let e = &mut r[f as usize];
        if e.1 == 0 {
            e.0 = t as u32;
        }
        e.1 += 1;
    }
    r
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

impl TurningWallDistance {
    /// From the parts: per cell its centroid at the angle 0 and whether it
    /// belongs to the rotor; per part (0 rotor, 1 stator) its wall faces and
    /// its open faces (every other face of the boundary of the part); the
    /// axis of the turn through `origin` along the unit vector `axis`.
    pub fn new(
        centroid: &[[f64; 3]],
        rotor: &[bool],
        walls: [SurfaceFaces; 2],
        open: [SurfaceFaces; 2],
        origin: [f64; 3],
        axis: [f64; 3],
    ) -> Self {
        let part: Vec<u32> = rotor.iter().map(|&r| u32::from(!r)).collect();
        Self::of_parts(
            centroid,
            &part,
            walls.into(),
            open.into(),
            vec![(origin, axis)],
        )
    }

    /// From `n` rotors and the stator (D185): per cell its centroid at the
    /// angle 0 and its part (rotor i < n, the stator n); per part its wall
    /// faces and its open faces; per rotor the axis of its turn (a point on
    /// it and the unit direction).
    pub fn of_parts(
        centroid: &[[f64; 3]],
        part: &[u32],
        walls: Vec<SurfaceFaces>,
        open: Vec<SurfaceFaces>,
        turns: Vec<([f64; 3], [f64; 3])>,
    ) -> Self {
        let np = turns.len() + 1;
        assert_eq!(centroid.len(), part.len(), "a part per cell");
        assert!(
            walls.len() == np && open.len() == np,
            "walls and open faces per part"
        );
        let ft: Vec<_> = walls.iter().map(face_tris).collect();
        let walls: Vec<Option<Bvh>> = walls
            .into_iter()
            .map(|w| (!w.is_empty()).then(|| Bvh::build(w)))
            .collect();
        let open: Vec<Option<Bvh>> = open
            .into_iter()
            .map(|o| (!o.is_empty()).then(|| Bvh::build(o)))
            .collect();
        let others_have_walls = |p: usize| (0..np).any(|q| q != p && walls[q].is_some());
        let pairs: Vec<(f64, bool)> = (0..centroid.len())
            .into_par_iter()
            .map(|c| {
                let x = centroid[c];
                let p = part[c] as usize;
                let own = walls[p]
                    .as_ref()
                    .and_then(|b| b.closest(x))
                    .map_or(f64::INFINITY, |q| q.distance);
                // the other walls can be nearer only through an open face
                let reach = others_have_walls(p)
                    && open[p].as_ref().is_some_and(|o| {
                        o.closest_within(x, own).is_some_and(|q| q.distance < own)
                    });
                (own, reach)
            })
            .collect();
        let own = pairs.iter().map(|p| p.0).collect();
        let candidates = pairs
            .iter()
            .enumerate()
            .filter(|(_, p)| p.1)
            .map(|(c, _)| (c as u32, vec![u32::MAX; np]))
            .collect();
        Self {
            turns,
            part: part.to_vec(),
            walls,
            face_tris: ft,
            own,
            candidates,
        }
    }

    /// The parts of `mesh` (the base mesh, at the angle 0): the cells
    /// `rotor` turn about the axis through `origin` along the unit vector
    /// `axis`; the faces of the boundaries `walls` (indices into
    /// `mesh.boundaries`) are walls of the part of their cell, every other
    /// boundary face is an open face of its part, and a face between a
    /// rotor and a stator cell is an open face of both. A planar boundary
    /// face whose plane has the whole domain on one side at every angle
    /// (the two planes of a mesh of one cell layer, a flat outlet) is not
    /// open: a straight line between two points of the domain cannot cross
    /// it.
    pub fn of_mesh(
        mesh: &Mesh,
        geo: &Geometry,
        rotor: &[bool],
        walls: &[usize],
        origin: [f64; 3],
        axis: [f64; 3],
    ) -> Self {
        let part: Vec<u32> = rotor.iter().map(|&r| u32::from(!r)).collect();
        Self::of_mesh_parts(mesh, geo, &part, walls, vec![(origin, axis)])
    }

    /// [`Self::of_mesh`] for `n` rotors (D185): `part` of every cell (rotor
    /// i < n, the stator n), `turns` the axis of every rotor.
    pub fn of_mesh_parts(
        mesh: &Mesh,
        geo: &Geometry,
        part: &[u32],
        walls: &[usize],
        turns: Vec<([f64; 3], [f64; 3])>,
    ) -> Self {
        let np = turns.len() + 1;
        let mut w = vec![SurfaceFaces::default(); np];
        let mut o = vec![SurfaceFaces::default(); np];
        let mut support = Support::new(mesh, part, &turns);
        for f in 0..mesh.n_interior_faces() {
            let (a, b) = (
                part[mesh.owner[f] as usize],
                part[mesh.neighbour[f] as usize],
            );
            if a != b {
                o[a as usize].push_face(mesh, geo, f);
                o[b as usize].push_face(mesh, geo, f);
            }
        }
        for (bi, b) in mesh.boundaries.iter().enumerate() {
            let wall = walls.contains(&bi);
            for f in b.range() {
                let p = part[mesh.owner[f] as usize] as usize;
                if wall {
                    w[p].push_face(mesh, geo, f);
                } else if !support.holds(mesh, geo, f) {
                    o[p].push_face(mesh, geo, f);
                }
            }
        }
        Self::of_parts(&geo.cell_centroid[..mesh.n_cells], part, w, o, turns)
    }

    /// Cells that may see a wall of another part nearer than their own
    /// (the queries of every angle; the others keep their distance).
    pub fn candidates(&self) -> usize {
        self.candidates.len()
    }

    /// The wall distance of every cell with the rotor turned by `angle`
    /// (rad, about the axis) from the base mesh; `centroid` holds the cell
    /// centroids of the mesh at that angle (the rotor cells turned).
    pub fn at(&mut self, angle: f64, centroid: &[[f64; 3]]) -> Vec<f64> {
        self.at_angles(&[angle], centroid)
    }

    /// The wall distance of every cell with rotor i turned by `angles[i]`
    /// from the base mesh; `centroid` holds the cell centroids of the mesh
    /// at these angles.
    pub fn at_angles(&mut self, angles: &[f64], centroid: &[[f64; 3]]) -> Vec<f64> {
        assert_eq!(centroid.len(), self.own.len(), "a centroid per cell");
        assert_eq!(angles.len(), self.turns.len(), "an angle per rotor");
        let mut d = self.own.clone();
        let np = self.turns.len() + 1;
        let (walls, ft, part, own, turns) = (
            &self.walls,
            &self.face_tris,
            &self.part,
            &self.own,
            &self.turns,
        );
        let found: Vec<f64> = self
            .candidates
            .par_iter_mut()
            .map(|(c, hints)| {
                let c = *c as usize;
                let p = part[c] as usize;
                let mut best = f64::INFINITY;
                for q in (0..np).filter(|&q| q != p) {
                    let Some(bvh) = walls[q].as_ref() else {
                        continue;
                    };
                    // the point in the frame of the walls of part q
                    let x = if q < turns.len() {
                        let (o, a) = turns[q];
                        rotate(centroid[c], o, a, -angles[q])
                    } else {
                        centroid[c]
                    };
                    let hint = &mut hints[q];
                    // the face nearest at the last angle bounds the search
                    let mut near = f64::INFINITY;
                    if *hint != u32::MAX {
                        let (s, n) = ft[q][*hint as usize];
                        let tris = &bvh.faces().tris;
                        for t in &tris[s as usize..(s + n) as usize] {
                            near = near.min(dist(x, closest_on_triangle(x, t)));
                        }
                    }
                    let r = own[c].min(near);
                    // (a hair wider than the bound, which is a square root of
                    // a squared distance: the face of the hint is found again)
                    match bvh.closest_within(x, r * (1.0 + 1e-12)) {
                        Some(hit) => {
                            *hint = hit.face as u32;
                            best = best.min(hit.distance);
                        }
                        None => {
                            if near < own[c] {
                                best = best.min(near);
                            }
                        }
                    }
                }
                best
            })
            .collect();
        for (&(c, _), f) in self.candidates.iter().zip(found) {
            let c = c as usize;
            d[c] = d[c].min(f);
        }
        d
    }
}

/// Planes that have the whole domain on one side at every rotor angle.
struct Support {
    /// Per node: the rotor it turns with (`u32::MAX`: it stands).
    turns: Vec<u32>,
    /// Axis of every rotor.
    axes: Vec<([f64; 3], [f64; 3])>,
    tol: f64,
    /// Planes already decided (normal and offset in units of `tol`).
    seen: std::collections::HashMap<[i64; 4], bool>,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl Support {
    fn new(mesh: &Mesh, part: &[u32], axes: &[([f64; 3], [f64; 3])]) -> Self {
        let n = axes.len() as u32;
        let mut turns = vec![u32::MAX; mesh.points.len()];
        for f in 0..mesh.n_faces() {
            let mut mark = |c: u32| {
                let p = part[c as usize];
                if p < n {
                    for &v in mesh.face(f) {
                        turns[v as usize] = p;
                    }
                }
            };
            mark(mesh.owner[f]);
            if f < mesh.n_interior_faces() {
                mark(mesh.neighbour[f]);
            }
        }
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for p in &mesh.points {
            for i in 0..3 {
                lo[i] = lo[i].min(p[i]);
                hi[i] = hi[i].max(p[i]);
            }
        }
        Self {
            turns,
            axes: axes.to_vec(),
            tol: 1e-9 * dist(lo, hi).max(f64::MIN_POSITIVE),
            seen: Default::default(),
        }
    }

    /// The boundary face `f` is planar and its plane has every point of the
    /// domain, the rotor at any angle, on its inner side.
    fn holds(&mut self, mesh: &Mesh, geo: &Geometry, f: usize) -> bool {
        let s = geo.face_area[f];
        let a = dot(s, s).sqrt();
        if a == 0.0 {
            return false;
        }
        let n = s.map(|x| x / a);
        let c = geo.face_centroid[f];
        let tol = self.tol;
        if mesh
            .face(f)
            .iter()
            .any(|&v| dot(n, sub(mesh.points[v as usize], c)).abs() > tol)
        {
            return false;
        }
        let off = dot(n, c);
        let key = [
            (n[0] * 1e9).round() as i64,
            (n[1] * 1e9).round() as i64,
            (n[2] * 1e9).round() as i64,
            (off / tol).round() as i64,
        ];
        if let Some(&h) = self.seen.get(&key) {
            return h;
        }
        // the largest n·p of a point p turned about the axis of its rotor
        let h = mesh.points.iter().zip(&self.turns).all(|(&p, &t)| {
            let reach = if t != u32::MAX {
                let (o, ax) = self.axes[t as usize];
                let na = dot(n, ax);
                let n_perp = (dot(n, n) - na * na).max(0.0).sqrt();
                let v = sub(p, o);
                let va = dot(v, ax);
                let v_perp = (dot(v, v) - va * va).max(0.0).sqrt();
                dot(n, o) + na * va + n_perp * v_perp
            } else {
                dot(n, p)
            };
            reach <= off + tol
        });
        self.seen.insert(key, h);
        h
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(s: &mut u64) -> f64 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*s >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn surface(tris: Vec<[[f64; 3]; 3]>) -> SurfaceFaces {
        let n = tris.len();
        SurfaceFaces {
            tris,
            tri_face: (0..n as u32).collect(),
            global_id: (0..n as u64).collect(),
            mesh_face: (0..n).collect(),
            normal: vec![[0.0, 0.0, 1.0]; n],
        }
    }

    /// A cylinder r = `r` for |z| ≤ 1 as `n` × 2 triangles.
    fn cylinder(r: f64, n: usize) -> Vec<[[f64; 3]; 3]> {
        let p = |i: usize, z: f64| {
            let t = std::f64::consts::TAU * i as f64 / n as f64;
            [r * t.cos(), r * t.sin(), z]
        };
        (0..n)
            .flat_map(|i| {
                [
                    [p(i, -1.0), p(i + 1, -1.0), p(i + 1, 1.0)],
                    [p(i, -1.0), p(i + 1, 1.0), p(i, 1.0)],
                ]
            })
            .collect()
    }

    /// A flat plate (two triangles) from radius r0 to r1 at the angle t,
    /// for |z| ≤ h.
    fn plate(r0: f64, r1: f64, t: f64, h: f64) -> Vec<[[f64; 3]; 3]> {
        let (s, c) = t.sin_cos();
        let p = |r: f64, z: f64| [r * c, r * s, z];
        vec![
            [p(r0, -h), p(r1, -h), p(r1, h)],
            [p(r0, -h), p(r1, h), p(r0, h)],
        ]
    }

    /// A rotor of three blades on a hub inside an interface at r = 1, a
    /// stator ring with baffles and an outer wall at r = 2, between two
    /// planes: at every angle the distances are the
    /// brute-force minimum over the turned rotor walls and the stator
    /// walls, and only some of the cells search.
    #[test]
    fn turned_parts_give_the_brute_force_distance() {
        let tau = std::f64::consts::TAU;
        let mut rw = cylinder(0.15, 24);
        for k in 0..3 {
            rw.extend(plate(0.15, 0.8, k as f64 * tau / 3.0, 0.7));
        }
        let mut sw = cylinder(2.0, 64);
        for k in 0..4 {
            sw.extend(plate(1.5, 2.0, 0.3 + k as f64 * tau / 4.0, 1.0));
        }
        // (the planes z = ±1 that close both parts have the domain on one
        // side: no open faces, as TurningWallDistance::of_mesh finds)
        let ro = cylinder(1.0, 48);
        let so = cylinder(1.0, 48);
        let mut s = 7u64;
        let mut x0 = Vec::new();
        let mut rotor = Vec::new();
        while x0.len() < 3000 {
            let r = 0.16 + 1.83 * lcg(&mut s);
            let t = tau * lcg(&mut s);
            let z = 1.98 * lcg(&mut s) - 0.99;
            x0.push([r * t.cos(), r * t.sin(), z]);
            rotor.push(r < 1.0);
        }
        let (o, a) = ([0.0; 3], [0.0, 0.0, 1.0]);
        let mut twd = TurningWallDistance::new(
            &x0,
            &rotor,
            [surface(rw.clone()), surface(sw.clone())],
            [surface(ro), surface(so)],
            o,
            a,
        );
        assert!(
            twd.candidates() > 0 && twd.candidates() < x0.len() / 2,
            "{} of {}",
            twd.candidates(),
            x0.len()
        );
        let brute = |x: [f64; 3], tris: &[[[f64; 3]; 3]]| {
            tris.iter()
                .map(|t| dist(x, closest_on_triangle(x, t)))
                .fold(f64::INFINITY, f64::min)
        };
        for angle in [0.0, 0.1, 0.25, 1.3, 2.9, -0.7] {
            let turn = |p: [f64; 3]| rotate(p, o, a, angle);
            let x: Vec<[f64; 3]> = x0
                .iter()
                .zip(&rotor)
                .map(|(&p, &r)| if r { turn(p) } else { p })
                .collect();
            let rw_t: Vec<[[f64; 3]; 3]> = rw.iter().map(|t| t.map(turn)).collect();
            let d = twd.at(angle, &x);
            for (c, &xc) in x.iter().enumerate() {
                let want = brute(xc, &rw_t).min(brute(xc, &sw));
                assert!(
                    (d[c] - want).abs() <= 1e-12 * want.max(1.0),
                    "angle {angle}, cell {c}: {} vs {want}",
                    d[c]
                );
            }
        }
    }

    /// The split ring of the sliding mesh by imprinting: the planes of the
    /// cell layer are not open, only cells nearer to the interface than to
    /// their own wall search, and the distances are those of the wall
    /// distance on the imprinted mesh of every angle.
    #[test]
    fn the_imprinted_ring_agrees_with_the_wall_distance_of_its_mesh() {
        use crate::annulus::AnnulusSpec;
        use crate::imprint::{Interface, InterfaceSpec};
        let mesh = AnnulusSpec {
            n_theta: 12,
            n_theta_stator: 18,
            n_r: 8,
            rotor_layers: 4,
            radial_ratio: 3.0,
            ..Default::default()
        }
        .generate()
        .unwrap();
        let base = Geometry::compute(&mesh);
        let spec = InterfaceSpec {
            zone: "rotor".into(),
            rotor_patch: "rotor_interface".into(),
            stator_patch: "stator_interface".into(),
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
        };
        let ifc = Interface::new(&mesh, &spec).unwrap();
        let walls_of = |m: &Mesh| -> Vec<usize> {
            m.boundaries
                .iter()
                .enumerate()
                .filter(|(_, b)| b.name == "inner" || b.name == "outer")
                .map(|(i, _)| i)
                .collect()
        };
        let mut twd = TurningWallDistance::of_mesh(
            &mesh,
            &base,
            ifc.rotor_cells(),
            &walls_of(&mesh),
            ifc.origin(),
            ifc.axis(),
        );
        let n = mesh.n_cells;
        // (the half of each part nearer to the interface; with the planes
        // of the layer open, every cell would search)
        assert!(
            twd.candidates() > 0 && twd.candidates() <= n / 2,
            "{} of {n}",
            twd.candidates()
        );
        for angle in [0.0, 0.05, 0.4, 1.7, -2.2] {
            let a = ifc.assemble(&mesh, angle).unwrap();
            let geo = ifc.geometry(&a, &base);
            let want = super::super::wall_distance(&a.mesh, &geo, &walls_of(&a.mesh)).unwrap();
            let d = twd.at(angle, &geo.cell_centroid[..n]);
            for c in 0..n {
                let w = want[c].distance;
                assert!(
                    (d[c] - w).abs() <= 1e-12 * w.max(1.0),
                    "angle {angle}, cell {c}: {} vs {w}",
                    d[c]
                );
            }
        }
    }

    #[test]
    fn two_rotors_agree_with_the_wall_distance_of_their_mesh() {
        use crate::annulus::AnnulusSpec;
        use crate::imprint::{Imprint, Interface, InterfaceSpec};
        let mesh = AnnulusSpec {
            n_theta: 12,
            n_theta_stator: 18,
            n_theta_outer: 15,
            n_r: 9,
            n_z: 2,
            rotor_layers: 2,
            outer_rotor_layers: 5,
            ..Default::default()
        }
        .generate()
        .unwrap();
        let base = Geometry::compute(&mesh);
        let spec = InterfaceSpec {
            zone: "rotor".into(),
            rotor_patch: "rotor_interface".into(),
            stator_patch: "stator_interface".into(),
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
        };
        let outer = InterfaceSpec {
            zone: "rotor_outer".into(),
            rotor_patch: "rotor_interface_outer".into(),
            stator_patch: "stator_interface_outer".into(),
            ..spec.clone()
        };
        let parts = vec![
            Interface::new(&mesh, &spec).unwrap(),
            Interface::new(&mesh, &outer).unwrap(),
        ];
        let part: Vec<u32> = (0..mesh.n_cells)
            .map(|c| {
                if parts[0].rotor_cells()[c] {
                    0
                } else if parts[1].rotor_cells()[c] {
                    1
                } else {
                    2
                }
            })
            .collect();
        let imp = Imprint::new(&mesh, parts).unwrap();
        let walls_of = |m: &Mesh| -> Vec<usize> {
            m.boundaries
                .iter()
                .enumerate()
                .filter(|(_, b)| b.name == "inner" || b.name == "outer")
                .map(|(i, _)| i)
                .collect()
        };
        let axis = ([0.0; 3], [0.0, 0.0, 1.0]);
        let mut twd = TurningWallDistance::of_mesh_parts(
            &mesh,
            &base,
            &part,
            &walls_of(&mesh),
            vec![axis, axis],
        );
        let n = mesh.n_cells;
        assert!(
            twd.candidates() > 0 && twd.candidates() < n,
            "{}",
            twd.candidates()
        );
        for angles in [[0.0, 0.0], [0.05, -0.3], [0.4, 1.0], [1.7, -2.2]] {
            let a = imp.assemble(&mesh, &angles).unwrap();
            let geo = imp.geometry(&a, &base);
            let want = super::super::wall_distance(&a.mesh, &geo, &walls_of(&a.mesh)).unwrap();
            let d = twd.at_angles(&angles, &geo.cell_centroid[..n]);
            for c in 0..n {
                let w = want[c].distance;
                assert!(
                    (d[c] - w).abs() <= 1e-12 * w.max(1.0),
                    "angles {angles:?}, cell {c}: {} vs {w}",
                    d[c]
                );
            }
        }
    }
}

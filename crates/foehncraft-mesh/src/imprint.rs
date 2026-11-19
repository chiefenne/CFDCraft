//! The sliding interface by imprinting (D58, D182 variant B): a rotor and a
//! stator whose faces on the interface differ in angular pitch.
//!
//! The base mesh holds the interface as two boundaries, `rotor_interface`
//! (faces of rotor cells) and `stator_interface` (faces of stator cells),
//! on one circular cylinder about the axis, and no interior face between the
//! two zones. [`Interface::assemble`] gives the conformal mesh of one rotor
//! angle: the nodes of the rotor are turned, and the two boundaries are
//! replaced by interior faces, one for every non-empty intersection of a rotor
//! face with a stator face (the supermesh of Farrell et al. 2009). Every
//! intersection polygon is an ordinary internal face between a rotor cell and
//! a stator cell, so gradients, Rhie–Chow, the multigrid and the halo plans
//! see a conformal polyhedral mesh and need no interpolation weights.
//!
//! **The polygons.** The surface is unrolled into the coordinates (θ, z),
//! θ the angle about the axis and z the axial position. A face is a convex
//! polygon there (its edges are circular arcs or axial lines); the rotor
//! faces are shifted by the rotor angle, an intersection is the clip of two
//! convex polygons (Sutherland–Hodgman), and its corners are mapped back onto
//! the cylinder. The area vector of a face made of such polygons is the sum of
//! theirs only along the axis of the extrusion (the polyline of sub-chords
//! against the chord of the face differ by a sliver whose area vector points
//! along z, and the same sliver at both ends of an axial edge cancels): the
//! cells of an extruded mesh (every cell of `[mesh.annulus]`, a Gmsh mesh
//! extruded along the axis) close to round-off, and the closure defect of
//! other faces is reported by [`Interface::closure`].
//!
//! **Cell geometry.** A rotor cell moves rigidly, so its volume is the one of
//! the base mesh and its centroid is turned with the nodes
//! ([`Interface::geometry`]): the geometric conservation law holds exactly. The
//! faces of an assembled mesh, those of the intersections included, are
//! computed from their nodes.

use crate::builder::MeshBuilder;
use crate::geometry::{Geometry, add, cross, dot, norm, polygon_area_centroid, scale, sub};
use crate::mesh::Mesh;
use std::f64::consts::{PI, TAU};

/// Smallest intersection kept, relative to the smaller of the two faces (the
/// rest is a sliver of two nearly coincident edges).
const SLIVER: f64 = 1e-10;

/// What the base mesh names.
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceSpec {
    /// The cell zone that turns.
    pub zone: String,
    /// The boundary of the rotor's interface faces.
    pub rotor_patch: String,
    /// The boundary of the stator's interface faces.
    pub stator_patch: String,
    /// A point on the axis.
    pub origin: [f64; 3],
    /// Direction of the axis (any length).
    pub axis: [f64; 3],
}

/// A face of one of the two interface boundaries in the unrolled surface.
#[derive(Debug, Clone)]
struct IfFace {
    /// The face of the base mesh.
    face: u32,
    /// Its owner cell.
    cell: u32,
    /// Corners (θ, z), counter-clockwise, θ unwrapped around the first.
    poly: Vec<[f64; 2]>,
    /// +1 when the order of `nodes` is counter-clockwise in (θ, z), else −1
    /// (the normal of the face points away from the axis for +1).
    sign: f64,
    /// θ range of the corners.
    lo: f64,
    hi: f64,
    /// Area in (θ, z).
    area: f64,
}

/// The imprint of a base mesh: prepared once, assembled for each rotor angle.
#[derive(Debug, Clone)]
pub struct Interface {
    origin: [f64; 3],
    axis: [f64; 3],
    e1: [f64; 3],
    e2: [f64; 3],
    radius: f64,
    rotor_cell: Vec<bool>,
    rotor_node: Vec<bool>,
    rotor: Vec<IfFace>,
    stator: Vec<IfFace>,
    /// Indices of the two patches in `Mesh::boundaries`.
    patches: [usize; 2],
    /// Stator faces by angular bin.
    bins: Vec<Vec<u32>>,
}

/// Where a face of an assembled mesh comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceOrigin {
    /// The face of the base mesh with this index.
    Base(u32),
    /// The intersection polygon with this index in [`Assembled::pieces`].
    Piece(u32),
}

/// One intersection of a rotor face and a stator face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Piece {
    /// The rotor face (of the base mesh).
    pub rotor_face: u32,
    /// The stator face (of the base mesh).
    pub stator_face: u32,
    /// The rotor cell.
    pub rotor_cell: u32,
    /// The stator cell.
    pub stator_cell: u32,
    /// Area of the polygon (m²).
    pub area: f64,
}

/// The mesh of one rotor angle.
#[derive(Debug, Clone)]
pub struct Assembled {
    /// The conformal mesh; the cells are those of the base mesh.
    pub mesh: Mesh,
    /// The origin of each face of `mesh`.
    pub origin: Vec<FaceOrigin>,
    /// The intersections, in order of creation.
    pub pieces: Vec<Piece>,
    /// The rotor angle (rad; of the first rotor).
    pub angle: f64,
    /// The angle of every rotor (rad).
    pub angles: Vec<f64>,
}

/// Rotation of `p` about the axis through `o` along the unit vector `a` by
/// `angle` (Rodrigues).
pub fn rotate(p: [f64; 3], o: [f64; 3], a: [f64; 3], angle: f64) -> [f64; 3] {
    let v = sub(p, o);
    let (s, c) = angle.sin_cos();
    let r = add(
        add(scale(v, c), scale(cross(a, v), s)),
        scale(a, dot(a, v) * (1.0 - c)),
    );
    add(o, r)
}

/// The angle in (−π, π].
fn wrap(x: f64) -> f64 {
    let y = x - TAU * (x / TAU).round();
    if y <= -PI { y + TAU } else { y }
}

/// Twice the signed area of a polygon.
fn area2(p: &[[f64; 2]]) -> f64 {
    let n = p.len();
    (0..n)
        .map(|i| {
            let (a, b) = (p[i], p[(i + 1) % n]);
            a[0] * b[1] - a[1] * b[0]
        })
        .sum()
}

/// The convex polygon `subject` clipped to the convex polygon `clip` (both
/// counter-clockwise; Sutherland–Hodgman).
fn clip(subject: &[[f64; 2]], clip: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut out: Vec<[f64; 2]> = subject.to_vec();
    let m = clip.len();
    for i in 0..m {
        if out.is_empty() {
            break;
        }
        let (a, b) = (clip[i], clip[(i + 1) % m]);
        let side = |p: [f64; 2]| (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
        let input = std::mem::take(&mut out);
        let n = input.len();
        for j in 0..n {
            let (p, q) = (input[j], input[(j + 1) % n]);
            let (sp, sq) = (side(p), side(q));
            if sp >= 0.0 {
                out.push(p);
            }
            if (sp > 0.0 && sq < 0.0) || (sp < 0.0 && sq > 0.0) {
                let t = sp / (sp - sq);
                out.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
            }
        }
    }
    out
}

/// The polygon without repeated corners.
fn dedupe(p: Vec<[f64; 2]>, tol: f64) -> Vec<[f64; 2]> {
    let mut out: Vec<[f64; 2]> = Vec::with_capacity(p.len());
    for q in p {
        if out
            .last()
            .is_none_or(|l| (l[0] - q[0]).abs() > tol || (l[1] - q[1]).abs() > tol)
        {
            out.push(q);
        }
    }
    while out.len() > 1 {
        let (f, l) = (out[0], out[out.len() - 1]);
        if (f[0] - l[0]).abs() <= tol && (f[1] - l[1]).abs() <= tol {
            out.pop();
        } else {
            break;
        }
    }
    out
}

/// The bin of the angle `t` among `nb` bins of the circle.
fn bin_of(t: f64, nb: usize) -> usize {
    let u = (t / TAU - (t / TAU).floor()).clamp(0.0, 1.0 - 1e-15);
    ((u * nb as f64) as usize).min(nb - 1)
}

/// The bins of the circle that the angles from `lo` to `hi` cover.
fn bin_range(lo: f64, hi: f64, nb: usize) -> Vec<usize> {
    let (b0, b1) = (bin_of(lo, nb), bin_of(hi, nb));
    let n = (((hi - lo) / TAU * nb as f64).ceil() as usize + 2).min(nb);
    let mut out = Vec::with_capacity(n);
    let mut b = b0;
    for _ in 0..n {
        out.push(b);
        if b == b1 {
            break;
        }
        b = (b + 1) % nb;
    }
    out
}

impl Interface {
    /// Prepares the imprint of `mesh`: checks that the two boundaries lie on
    /// one circular cylinder about the axis, cover it, and separate a rotor
    /// from a stator that share no node.
    pub fn new(mesh: &Mesh, spec: &InterfaceSpec) -> Result<Self, String> {
        let n = norm(spec.axis);
        if !(n > 0.0 && n.is_finite()) {
            return Err("the axis has no direction".into());
        }
        let axis = scale(spec.axis, 1.0 / n);
        let ex = if axis[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let e1 = {
            let v = sub(ex, scale(axis, dot(ex, axis)));
            scale(v, 1.0 / norm(v))
        };
        let e2 = cross(axis, e1);
        let zi = mesh
            .zones
            .iter()
            .position(|z| *z == spec.zone)
            .ok_or_else(|| {
                format!(
                    "the mesh has no cell zone `{}` (zones: {})",
                    spec.zone,
                    mesh.zones.join(", ")
                )
            })? as u32;
        let rotor_cell: Vec<bool> = mesh.cell_zone.iter().map(|&z| z == zi).collect();
        let n_rotor = rotor_cell.iter().filter(|&&r| r).count();
        if n_rotor == 0 || n_rotor == mesh.n_cells {
            return Err(format!(
                "the zone `{}` holds {n_rotor} of {} cells: a rotor and a stator are needed",
                spec.zone, mesh.n_cells
            ));
        }
        let patch = |name: &str| -> Result<usize, String> {
            mesh.boundaries
                .iter()
                .position(|p| p.name == name)
                .ok_or_else(|| {
                    format!(
                        "the mesh has no boundary `{name}` (boundaries: {})",
                        mesh.boundaries
                            .iter()
                            .map(|p| p.name.as_str())
                            .collect::<Vec<_>>()
                            .join(", ")
                    )
                })
        };
        let patches = [patch(&spec.rotor_patch)?, patch(&spec.stator_patch)?];
        if patches[0] == patches[1] {
            return Err("the rotor and the stator interface are the same boundary".into());
        }
        // the zones are separated by the two boundaries
        for f in 0..mesh.n_interior_faces() {
            let (o, nb) = (mesh.owner[f] as usize, mesh.neighbour[f] as usize);
            if rotor_cell[o] != rotor_cell[nb] {
                return Err(format!(
                    "the cells {o} and {nb} of the rotor and the stator share the face {f}: the \
                     interface must be the two boundaries `{}` and `{}` of two separate meshes \
                     (a conformal interface slides by whole sectors, `[[sliding]]` without \
                     interface boundaries)",
                    spec.rotor_patch, spec.stator_patch
                ));
            }
        }
        // the nodes of the rotor turn, those of the stator stand
        let mut rotor_node = vec![false; mesh.n_points()];
        let mut stator_node = vec![false; mesh.n_points()];
        for (c, &is_rotor) in rotor_cell.iter().enumerate() {
            let mark = if is_rotor {
                &mut rotor_node
            } else {
                &mut stator_node
            };
            for (f, _) in mesh.cell_face_iter(c) {
                for &p in mesh.face(f) {
                    mark[p as usize] = true;
                }
            }
        }
        if let Some(p) = (0..mesh.n_points()).find(|&p| rotor_node[p] && stator_node[p]) {
            return Err(format!(
                "the node {p} belongs to rotor and stator cells: the two meshes must not share \
                 nodes"
            ));
        }
        // the unrolled faces
        let to_theta_z = |p: [f64; 3]| -> (f64, f64, f64) {
            let v = sub(p, spec.origin);
            let (x, y, z) = (dot(v, e1), dot(v, e2), dot(v, axis));
            (y.atan2(x), z, x.hypot(y))
        };
        let mut radii = Vec::new();
        let mut build = |pi: usize, want_rotor: bool| -> Result<Vec<IfFace>, String> {
            let b = &mesh.boundaries[pi];
            let mut out = Vec::with_capacity(b.len);
            for f in b.range() {
                let cell = mesh.owner[f];
                if rotor_cell[cell as usize] != want_rotor {
                    return Err(format!(
                        "the face {f} of the boundary `{}` belongs to a cell of the {}",
                        b.name,
                        if want_rotor { "stator" } else { "rotor" }
                    ));
                }
                let nodes = mesh.face(f);
                let mut poly: Vec<[f64; 2]> = Vec::with_capacity(nodes.len());
                let mut t0 = 0.0;
                for (i, &p) in nodes.iter().enumerate() {
                    let (t, z, r) = to_theta_z(mesh.points[p as usize]);
                    radii.push(r);
                    if i == 0 {
                        t0 = t;
                        poly.push([t, z]);
                    } else {
                        poly.push([t0 + wrap(t - t0), z]);
                    }
                }
                let a2 = area2(&poly);
                let sign = if a2 >= 0.0 { 1.0 } else { -1.0 };
                if sign < 0.0 {
                    poly.reverse();
                }
                let lo = poly.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
                let hi = poly.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
                if hi - lo >= 0.99 * PI {
                    return Err(format!(
                        "the face {f} of `{}` spans {:.1} degrees around the axis (at most 178 \
                         are supported)",
                        b.name,
                        (hi - lo).to_degrees()
                    ));
                }
                if a2.abs() <= 0.0 {
                    return Err(format!("the face {f} of `{}` has no area", b.name));
                }
                out.push(IfFace {
                    face: f as u32,
                    cell,
                    poly,
                    sign,
                    lo,
                    hi,
                    area: 0.5 * a2.abs(),
                });
            }
            Ok(out)
        };
        let rotor = build(patches[0], true)?;
        let stator = build(patches[1], false)?;
        if rotor.is_empty() || stator.is_empty() {
            return Err("an interface boundary has no faces".into());
        }
        let r0 = radii.iter().sum::<f64>() / radii.len() as f64;
        if let Some(r) = radii.iter().find(|&&r| (r - r0).abs() > 1e-6 * r0) {
            return Err(format!(
                "the interface is not a circular cylinder about the axis: nodes at the \
                 distances {r0:.6e} and {r:.6e} m (cones and discs are not supported yet)"
            ));
        }
        let (ar, as_): (f64, f64) = (
            rotor.iter().map(|f| f.area).sum(),
            stator.iter().map(|f| f.area).sum(),
        );
        if (ar - as_).abs() > 1e-8 * ar.max(as_) {
            return Err(format!(
                "the two interface boundaries do not cover the same surface (areas {ar:.8e} \
                 and {as_:.8e} in (θ, z)): a partial ring is not supported"
            ));
        }
        // stator faces by angular bin, for the search of the partners
        let nb = (2 * stator.len()).clamp(8, 1 << 20);
        let mut bins = vec![Vec::new(); nb];
        for (i, f) in stator.iter().enumerate() {
            for b in bin_range(f.lo, f.hi, nb) {
                bins[b].push(i as u32);
            }
        }
        Ok(Self {
            origin: spec.origin,
            axis,
            e1,
            e2,
            radius: r0,
            rotor_cell,
            rotor_node,
            rotor,
            stator,
            patches,
            bins,
        })
    }

    /// Number of faces of the rotor's and of the stator's boundary.
    pub fn faces(&self) -> (usize, usize) {
        (self.rotor.len(), self.stator.len())
    }

    /// The cells of the rotor.
    pub fn rotor_cells(&self) -> &[bool] {
        &self.rotor_cell
    }

    /// The unit axis.
    pub fn axis(&self) -> [f64; 3] {
        self.axis
    }

    /// A point on the axis.
    pub fn origin(&self) -> [f64; 3] {
        self.origin
    }

    /// Radius of the interface (m).
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Indices of the two interface boundaries in the base mesh.
    pub fn patches(&self) -> [usize; 2] {
        self.patches
    }

    /// The point of the cylinder at (θ, z).
    fn point(&self, t: f64, z: f64) -> [f64; 3] {
        let (s, c) = t.sin_cos();
        add(
            add(
                self.origin,
                scale(add(scale(self.e1, c), scale(self.e2, s)), self.radius),
            ),
            scale(self.axis, z),
        )
    }

    /// The intersections of the rotor faces, turned by `angle`, with the
    /// stator faces.
    fn intersections(&self, angle: f64) -> Vec<(u32, u32, Vec<[f64; 2]>, f64)> {
        let nb = self.bins.len();
        let mut out = Vec::new();
        let mut seen: Vec<u32> = Vec::new();
        for (ri, r) in self.rotor.iter().enumerate() {
            let (lo, hi) = (r.lo + angle, r.hi + angle);
            // the stator faces in the bins of [lo, hi]
            seen.clear();
            for b in bin_range(lo, hi, nb) {
                seen.extend_from_slice(&self.bins[b]);
            }
            seen.sort_unstable();
            seen.dedup();
            let r_mid = 0.5 * (lo + hi);
            for &si in &seen {
                let s = &self.stator[si as usize];
                let shift = angle + TAU * ((0.5 * (s.lo + s.hi) - r_mid) / TAU).round();
                if r.hi + shift <= s.lo || r.lo + shift >= s.hi {
                    continue;
                }
                let shifted: Vec<[f64; 2]> = r.poly.iter().map(|p| [p[0] + shift, p[1]]).collect();
                let tol = 1e-13 * (1.0 + s.lo.abs().max(s.hi.abs()));
                let c = dedupe(clip(&shifted, &s.poly), tol);
                if c.len() < 3 {
                    continue;
                }
                let a = 0.5 * area2(&c);
                if a <= SLIVER * r.area.min(s.area) {
                    continue;
                }
                out.push((ri as u32, si, c, a));
            }
        }
        out
    }

    /// The conformal mesh with the rotor turned by `angle` (rad, right-handed
    /// about the axis): the interface boundaries are replaced by the
    /// intersection polygons. The cells and their numbering are those of
    /// `mesh`, the other boundaries keep their order.
    pub fn assemble(&self, mesh: &Mesh, angle: f64) -> Result<Assembled, String> {
        assemble_parts(&[self], mesh, &[angle])
    }

    /// The geometry of an assembled mesh: faces from their nodes, cells as
    /// the rigid bodies of the base mesh (volumes unchanged, the centroids of
    /// the rotor turned by the angle).
    pub fn geometry(&self, a: &Assembled, base: &Geometry) -> Geometry {
        geometry_parts(&[self], a, base)
    }

    /// The largest closure defect of a cell of `a`: |Σ S| relative to Σ |S|.
    pub fn closure(&self, a: &Assembled, geo: &Geometry) -> f64 {
        (0..a.mesh.n_cells)
            .map(|c| {
                let s = geo.closure(&a.mesh, c);
                let tot: f64 = a
                    .mesh
                    .cell_face_iter(c)
                    .map(|(f, _)| norm(geo.face_area[f]))
                    .sum();
                norm(s) / tot
            })
            .fold(0.0, f64::max)
    }
}

/// Several rotors of one base mesh, each with its own interface to the
/// stator, assembled together (D185): every rotor turns by its own angle
/// about its own axis, and the stator stands. A rotor that borders another
/// rotor (a rotor–rotor interface) is refused.
#[derive(Debug, Clone)]
pub struct Imprint {
    parts: Vec<Interface>,
}

impl Imprint {
    /// The rotors `parts` of `mesh` (prepared with [`Interface::new`]): their
    /// cells and their interface boundaries must be distinct, and the stator
    /// side of every interface must be a cell of the stator.
    pub fn new(mesh: &Mesh, parts: Vec<Interface>) -> Result<Self, String> {
        if parts.is_empty() {
            return Err("no rotor".into());
        }
        for (i, p) in parts.iter().enumerate() {
            for (j, q) in parts.iter().enumerate().skip(i + 1) {
                if (0..mesh.n_cells).any(|c| p.rotor_cell[c] && q.rotor_cell[c]) {
                    return Err(format!("the rotors {} and {} share cells", i + 1, j + 1));
                }
                if p.patches.iter().any(|b| q.patches.contains(b)) {
                    return Err(format!(
                        "the rotors {} and {} name the same interface boundary",
                        i + 1,
                        j + 1
                    ));
                }
            }
            for (j, q) in parts.iter().enumerate() {
                if i == j {
                    continue;
                }
                if let Some(f) = p.stator.iter().find(|f| q.rotor_cell[f.cell as usize]) {
                    let b = &mesh.boundaries[p.patches[1]].name;
                    return Err(format!(
                        "the face {} of `{b}`, the stator side of the interface of rotor {}, \
                         belongs to a cell of rotor {}: an interface between two rotors is not \
                         available yet (each rotor needs its interface with the stator)",
                        f.face,
                        i + 1,
                        j + 1
                    ));
                }
            }
        }
        Ok(Self { parts })
    }

    /// The rotors.
    pub fn parts(&self) -> &[Interface] {
        &self.parts
    }

    /// The conformal mesh with rotor i turned by `angles[i]`.
    pub fn assemble(&self, mesh: &Mesh, angles: &[f64]) -> Result<Assembled, String> {
        let p: Vec<&Interface> = self.parts.iter().collect();
        assemble_parts(&p, mesh, angles)
    }

    /// The geometry of an assembled mesh (see [`Interface::geometry`]).
    pub fn geometry(&self, a: &Assembled, base: &Geometry) -> Geometry {
        let p: Vec<&Interface> = self.parts.iter().collect();
        geometry_parts(&p, a, base)
    }
}

/// The conformal mesh with the rotor of `parts[i]` turned by `angles[i]`.
fn assemble_parts(parts: &[&Interface], mesh: &Mesh, angles: &[f64]) -> Result<Assembled, String> {
    assert_eq!(parts.len(), angles.len(), "an angle per rotor");
    let mut points = mesh.points.clone();
    for (ifc, &angle) in parts.iter().zip(angles) {
        if angle != 0.0 {
            for (p, &turn) in points.iter_mut().zip(&ifc.rotor_node) {
                if turn {
                    *p = rotate(*p, ifc.origin, ifc.axis, angle);
                }
            }
        }
    }
    // the corners of the polygons are new nodes
    let mut piece_nodes: Vec<Vec<u32>> = Vec::new();
    let mut pieces = Vec::new();
    for (ifc, &angle) in parts.iter().zip(angles) {
        for (ri, si, poly, _) in &ifc.intersections(angle) {
            let (r, s) = (&ifc.rotor[*ri as usize], &ifc.stator[*si as usize]);
            // θ of the polygon is in the frame of the stator
            let ids: Vec<u32> = poly
                .iter()
                .map(|q| {
                    let id = points.len() as u32;
                    points.push(ifc.point(q[0], q[1]));
                    id
                })
                .collect();
            // the orientation of the rotor face: out of the rotor cell
            let ordered: Vec<u32> = if r.sign > 0.0 {
                ids
            } else {
                ids.into_iter().rev().collect()
            };
            piece_nodes.push(ordered);
            pieces.push(Piece {
                rotor_face: r.face,
                stator_face: s.face,
                rotor_cell: r.cell,
                stator_cell: s.cell,
                area: 0.0,
            });
        }
    }
    let mut b = MeshBuilder::new(points, mesh.n_cells);
    b.set_zones(mesh.zones.clone(), mesh.cell_zone.clone());
    let skipped: Vec<usize> = parts.iter().flat_map(|p| p.patches).collect();
    let skip = |name: &str| skipped.iter().any(|&i| mesh.boundaries[i].name == name);
    for p in &mesh.boundaries {
        if !skip(&p.name) {
            b.declare_boundary(&p.name);
        }
    }
    let mut origin: Vec<FaceOrigin> = Vec::new();
    for f in 0..mesh.n_interior_faces() {
        b.add_interior(mesh.face(f), mesh.owner[f], mesh.neighbour[f]);
        origin.push(FaceOrigin::Base(f as u32));
    }
    for p in &mesh.boundaries {
        if skip(&p.name) {
            continue;
        }
        for f in p.range() {
            b.add_boundary(mesh.face(f), mesh.owner[f], &p.name);
            origin.push(FaceOrigin::Base(f as u32));
        }
    }
    for (k, nodes) in piece_nodes.iter().enumerate() {
        b.add_interior(nodes, pieces[k].rotor_cell, pieces[k].stator_cell);
        origin.push(FaceOrigin::Piece(k as u32));
    }
    let (assembled, order) = b.finish_with_order()?;
    let origin: Vec<FaceOrigin> = order.iter().map(|&i| origin[i as usize]).collect();
    let mut out = Assembled {
        mesh: assembled,
        origin,
        pieces,
        angle: angles[0],
        angles: angles.to_vec(),
    };
    // the true areas
    for (f, o) in out.origin.iter().enumerate() {
        if let FaceOrigin::Piece(k) = o {
            let (s, _) = polygon_area_centroid(&out.mesh.points, out.mesh.face(f));
            out.pieces[*k as usize].area = norm(s);
        }
    }
    Ok(out)
}

/// The geometry of a mesh assembled from `parts` (see
/// [`Interface::geometry`]).
fn geometry_parts(parts: &[&Interface], a: &Assembled, base: &Geometry) -> Geometry {
    let (face_area, face_centroid) = Geometry::faces(&a.mesh);
    let n = a.mesh.n_cells;
    let mut cell_volume = vec![0.0; n + a.mesh.n_halo];
    let mut cell_centroid = vec![[0.0; 3]; n + a.mesh.n_halo];
    cell_volume[..n].copy_from_slice(&base.cell_volume[..n]);
    cell_centroid[..n].copy_from_slice(&base.cell_centroid[..n]);
    for (ifc, &angle) in parts.iter().zip(&a.angles) {
        for (c, x) in cell_centroid[..n].iter_mut().enumerate() {
            if ifc.rotor_cell[c] {
                *x = rotate(base.cell_centroid[c], ifc.origin, ifc.axis, angle);
            }
        }
    }
    Geometry {
        face_area,
        face_centroid,
        cell_volume,
        cell_centroid,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annulus::AnnulusSpec;

    fn spec() -> InterfaceSpec {
        InterfaceSpec {
            zone: "rotor".into(),
            rotor_patch: "rotor_interface".into(),
            stator_patch: "stator_interface".into(),
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
        }
    }

    fn annulus(nr: usize, ns: usize, nz: usize) -> Mesh {
        AnnulusSpec {
            n_theta: nr,
            n_theta_stator: ns,
            n_r: 4,
            n_z: nz,
            rotor_layers: 2,
            ..Default::default()
        }
        .generate()
        .unwrap()
    }

    #[test]
    fn the_pieces_cover_the_faces_and_every_cell_closes() {
        for (nr, ns, nz) in [(12, 18, 1), (12, 12, 1), (7, 11, 3), (16, 5, 2)] {
            let m = annulus(nr, ns, nz);
            let ifc = Interface::new(&m, &spec()).unwrap();
            let base = Geometry::compute(&m);
            for angle in [0.0, 0.1, 0.3333, 1.0, 2.5, -0.7, 7.0] {
                let a = ifc.assemble(&m, angle).unwrap();
                let geo = ifc.geometry(&a, &base);
                let defect = ifc.closure(&a, &geo);
                assert!(
                    defect < 1e-12,
                    "{nr}/{ns}/{nz} angle {angle}: closure {defect:e}"
                );
                // the area vectors of the pieces of a rotor face add up to
                // the area vector of the face
                let mut sum = vec![[0.0; 3]; m.n_faces()];
                let mut sum_s = vec![[0.0; 3]; m.n_faces()];
                for (f, o) in a.origin.iter().enumerate() {
                    if let FaceOrigin::Piece(k) = o {
                        let p = a.pieces[*k as usize];
                        // normal out of the owner of the assembled face
                        let s = geo.face_area[f];
                        let from_rotor = a.mesh.owner[f] == p.rotor_cell;
                        let s_r = if from_rotor { s } else { scale(s, -1.0) };
                        let rotated = |v: [f64; 3], inv: bool| {
                            rotate(
                                v,
                                [0.0; 3],
                                [0.0, 0.0, 1.0],
                                if inv { -angle } else { angle },
                            )
                        };
                        sum[p.rotor_face as usize] =
                            add(sum[p.rotor_face as usize], rotated(s_r, true));
                        sum_s[p.stator_face as usize] =
                            add(sum_s[p.stator_face as usize], scale(s_r, -1.0));
                    }
                }
                let mut covered = 0;
                for r in &ifc.rotor {
                    let f = r.face as usize;
                    let want = base.face_area[f];
                    let got = sum[f];
                    let d = norm(sub(want, got));
                    assert!(d < 1e-12 * norm(want), "rotor face {f}: {want:?} {got:?}");
                    covered += 1;
                }
                for s in &ifc.stator {
                    let f = s.face as usize;
                    let want = base.face_area[f];
                    let got = sum_s[f];
                    let d = norm(sub(want, got));
                    assert!(d < 1e-12 * norm(want), "stator face {f}: {want:?} {got:?}");
                    covered += 1;
                }
                assert_eq!(covered, nr * nz + ns * nz);
            }
        }
    }

    #[test]
    fn the_number_of_pieces_follows_the_vertices() {
        // 12 and 18 sectors: the vertices coincide every 2 pi / 6
        let m = annulus(12, 18, 1);
        let ifc = Interface::new(&m, &spec()).unwrap();
        let count = |angle: f64| ifc.assemble(&m, angle).unwrap().pieces.len();
        // generic angle: every vertex of both sides cuts: 12 + 18
        assert_eq!(count(0.01), 30);
        // aligned: the 6 coincident vertices save one piece each
        assert_eq!(count(0.0), 24);
        assert_eq!(count(TAU / 6.0), 24);
        assert_eq!(count(-TAU / 6.0 + TAU), 24);
        // the same pitch: aligned one piece per face, otherwise two
        let m = annulus(12, 12, 1);
        let ifc = Interface::new(&m, &spec()).unwrap();
        let count = |angle: f64| ifc.assemble(&m, angle).unwrap().pieces.len();
        assert_eq!(count(0.0), 12);
        assert_eq!(count(TAU / 12.0), 12);
        assert_eq!(count(0.1), 24);
    }

    #[test]
    fn the_rotor_is_turned_and_its_volumes_stay() {
        let m = annulus(12, 18, 1);
        let ifc = Interface::new(&m, &spec()).unwrap();
        let base = Geometry::compute(&m);
        let a = ifc.assemble(&m, 0.4).unwrap();
        let geo = ifc.geometry(&a, &base);
        for c in 0..m.n_cells {
            assert_eq!(geo.cell_volume[c], base.cell_volume[c]);
            let (x, y) = (base.cell_centroid[c][0], base.cell_centroid[c][1]);
            let (ca, sa) = (0.4f64.cos(), 0.4f64.sin());
            let want = if ifc.rotor_cell[c] {
                [ca * x - sa * y, sa * x + ca * y]
            } else {
                [x, y]
            };
            assert!((geo.cell_centroid[c][0] - want[0]).abs() < 1e-14);
            assert!((geo.cell_centroid[c][1] - want[1]).abs() < 1e-14);
        }
        // a full turn is the base mesh
        let full = ifc.assemble(&m, TAU).unwrap();
        let a0 = ifc.assemble(&m, 0.0).unwrap();
        assert_eq!(full.pieces.len(), a0.pieces.len());
        // the other boundaries keep their names and order
        let names: Vec<&str> = a.mesh.boundaries.iter().map(|p| p.name.as_str()).collect();
        assert_eq!(names, ["inner", "outer", "zmin", "zmax"]);
        assert_eq!(a.mesh.n_cells, m.n_cells);
    }

    #[test]
    fn two_rotors_turn_at_their_own_angles() {
        let m = AnnulusSpec {
            n_theta: 12,
            n_theta_stator: 18,
            n_theta_outer: 14,
            n_r: 6,
            n_z: 2,
            rotor_layers: 2,
            outer_rotor_layers: 2,
            ..Default::default()
        }
        .generate()
        .unwrap();
        let inner = Interface::new(&m, &spec()).unwrap();
        let outer = Interface::new(
            &m,
            &InterfaceSpec {
                zone: "rotor_outer".into(),
                rotor_patch: "rotor_interface_outer".into(),
                stator_patch: "stator_interface_outer".into(),
                ..spec()
            },
        )
        .unwrap();
        let imp = Imprint::new(&m, vec![inner.clone(), outer.clone()]).unwrap();
        let base = Geometry::compute(&m);
        for angles in [[0.0, 0.0], [0.3, -0.7], [1.1, 2.9], [-2.0, 0.05]] {
            let a = imp.assemble(&m, &angles).unwrap();
            let geo = imp.geometry(&a, &base);
            assert!(inner.closure(&a, &geo) < 1e-12, "{angles:?}");
            // the pieces are those of each interface alone
            let single = inner.assemble(&m, angles[0]).unwrap().pieces.len()
                + outer.assemble(&m, angles[1]).unwrap().pieces.len();
            assert_eq!(a.pieces.len(), single);
            let names: Vec<&str> = a.mesh.boundaries.iter().map(|p| p.name.as_str()).collect();
            assert_eq!(names, ["inner", "outer", "zmin", "zmax"]);
            for c in 0..m.n_cells {
                let turn = if inner.rotor_cell[c] {
                    angles[0]
                } else if outer.rotor_cell[c] {
                    angles[1]
                } else {
                    0.0
                };
                let want = rotate(base.cell_centroid[c], [0.0; 3], [0.0, 0.0, 1.0], turn);
                assert!(norm(sub(geo.cell_centroid[c], want)) < 1e-14);
                assert_eq!(geo.cell_volume[c], base.cell_volume[c]);
            }
        }
        // the same rotor twice
        let e = Imprint::new(&m, vec![inner.clone(), inner]).unwrap_err();
        assert!(e.contains("share cells"), "{e}");
    }

    #[test]
    fn what_does_not_fit_is_refused() {
        // a conformal mesh has no interface boundaries
        let conformal = AnnulusSpec {
            n_theta: 12,
            n_r: 4,
            rotor_layers: 2,
            ..Default::default()
        }
        .generate()
        .unwrap();
        let e = Interface::new(&conformal, &spec()).unwrap_err();
        assert!(e.contains("no boundary `rotor_interface`"), "{e}");
        let m = annulus(12, 18, 1);
        let mut s = spec();
        s.zone = "rotr".into();
        assert!(Interface::new(&m, &s).unwrap_err().contains("no cell zone"));
        let mut s = spec();
        s.stator_patch = "rotor_interface".into();
        assert!(
            Interface::new(&m, &s)
                .unwrap_err()
                .contains("same boundary")
        );
        let mut s = spec();
        s.rotor_patch = "stator_interface".into();
        s.stator_patch = "rotor_interface".into();
        let e = Interface::new(&m, &s).unwrap_err();
        assert!(e.contains("belongs to a cell of the"), "{e}");
        // an axis off the centre: the interface is not a cylinder about it
        let mut s = spec();
        s.origin = [0.01, 0.0, 0.0];
        let e = Interface::new(&m, &s).unwrap_err();
        assert!(e.contains("not a circular cylinder"), "{e}");
        let mut s = spec();
        s.axis = [0.0; 3];
        assert!(Interface::new(&m, &s).unwrap_err().contains("no direction"));
    }
}

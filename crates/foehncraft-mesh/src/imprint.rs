//! The sliding interface by imprinting (D58, D182 variant B, D186): a rotor
//! and a stator whose faces on the interface differ in pitch.
//!
//! The base mesh holds the interface as two boundaries, `rotor_interface`
//! (faces of rotor cells) and `stator_interface` (faces of stator cells),
//! on one surface of revolution about the axis, and no interior face between
//! the two zones. [`Interface::assemble`] gives the conformal mesh of one
//! rotor angle: the nodes of the rotor are turned, and the two boundaries are
//! replaced by interior faces, one for every non-empty intersection of a rotor
//! face with a stator face (the supermesh of Farrell et al. 2009). Every
//! intersection polygon is an ordinary internal face between a rotor cell and
//! a stator cell, so gradients, Rhie–Chow, the multigrid and the halo plans
//! see a conformal polyhedral mesh and need no interpolation weights.
//!
//! **The surface (D186).** One interface per rotor, of any shape that is
//! rotationally symmetric about the rotor's axis and whose meridian (its
//! outline in the (r, z) half-plane) is made of straight pieces: a cylinder
//! (a radial gap), a disc (an axial gap between two blade rows), a cone, or
//! any chain of them (a pot of disc, cylinder and disc around an impeller).
//! The meridian is found from the nodes of the two boundaries: every face lies
//! on one straight piece (the meshes have a ring of nodes at every kink),
//! collinear faces form the pieces, and the pieces are chained end to end.
//! The meridian may have several chains (two discs before and after a rotor
//! row); a chain may end on the axis (a disc through the centre), and a
//! branched or closed meridian is refused.
//!
//! **The polygons.** The surface is unrolled into the coordinates (θ, m), θ
//! the angle about the axis and m the arc length along the meridian. A face
//! is a convex polygon there; the rotor faces are shifted by the rotor angle,
//! an intersection is the clip of two convex polygons (Sutherland–Hodgman),
//! and its corners are mapped back onto the surface. A node on the axis has
//! no angle: in the faces around it, it becomes the edge m = m_axis between
//! the angles of its two neighbours (the polar chart), and the corners of a
//! polygon on that edge map back to the one point on the axis.
//!
//! **Closed cells.** The corners of the pieces lie on the surface, and the
//! edges of the faces of the base mesh are chords. Where a piece corner lies
//! on an edge of a rotor or a stator face, it is inserted into the other faces
//! of the base mesh that share the edge (the faces next to the interface get
//! the corners on their edges), so every cell of an assembled mesh is a closed
//! polyhedron and closes to round-off on every surface ([`Interface::closure`]).
//!
//! **Cell geometry.** A rotor cell moves rigidly, so its volume is the one of
//! the base mesh and its centroid is turned with the nodes
//! ([`Interface::geometry`]): the geometric conservation law holds exactly. The
//! faces of an assembled mesh, those of the intersections included, are
//! computed from their nodes.

use crate::builder::MeshBuilder;
use crate::geometry::{Geometry, add, cross, dot, norm, polygon_area_centroid, scale, sub};
use crate::mesh::Mesh;
use std::collections::{HashMap, HashSet};
use std::f64::consts::{PI, TAU};

/// Smallest intersection kept, relative to the smaller of the two faces (the
/// rest is a sliver of two nearly coincident edges).
const SLIVER: f64 = 1e-10;

/// Distance (relative to the size of the surface) within which nodes lie on
/// the meridian, on the axis, or at a common end of two pieces.
const ON_SURFACE: f64 = 1e-6;

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
    /// Corners (θ, m / size), counter-clockwise, θ unwrapped around the
    /// first.
    poly: Vec<[f64; 2]>,
    /// The node of the base mesh at each corner (a node on the axis gives
    /// two corners).
    corner_node: Vec<u32>,
    /// +1 when the order of `nodes` is counter-clockwise in (θ, m), else −1.
    sign: f64,
    /// θ range of the corners.
    lo: f64,
    hi: f64,
    /// m range of the corners.
    mlo: f64,
    mhi: f64,
    /// Area in (θ, m).
    area: f64,
}

/// A straight piece of the meridian in the (r, z) half-plane.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Segment {
    /// (r, z) at the start, in the direction of increasing m.
    a: [f64; 2],
    /// Unit direction in (r, z).
    d: [f64; 2],
    /// Length (m).
    len: f64,
    /// m at the start.
    m0: f64,
}

impl Segment {
    fn at(&self, m: f64) -> [f64; 2] {
        let s = m - self.m0;
        [self.a[0] + s * self.d[0], self.a[1] + s * self.d[1]]
    }

    fn m_of(&self, p: [f64; 2]) -> f64 {
        self.m0 + (p[0] - self.a[0]) * self.d[0] + (p[1] - self.a[1]) * self.d[1]
    }
}

/// The meridian of the interface: chains of straight pieces, numbered by
/// arc length m (the chains apart by a gap).
#[derive(Debug, Clone)]
struct Meridian {
    /// The pieces in the order of m.
    pieces: Vec<Segment>,
    /// Number of chains.
    chains: usize,
    /// m of the ends of the chains that lie on the axis.
    axis_m: Vec<f64>,
}

impl Meridian {
    /// The piece that holds m.
    fn piece(&self, m: f64) -> &Segment {
        let i = self.pieces.partition_point(|s| s.m0 <= m);
        let i = i.saturating_sub(1);
        // between two chains (only by round-off): the nearer end
        let s = &self.pieces[i];
        if m > s.m0 + s.len
            && let Some(n) = self.pieces.get(i + 1)
            && n.m0 - m < m - (s.m0 + s.len)
        {
            return n;
        }
        s
    }

    /// (r, z) at m.
    fn rz(&self, m: f64) -> [f64; 2] {
        let s = self.piece(m);
        let mut p = s.at(m);
        p[0] = p[0].max(0.0);
        p
    }

    /// A description, e.g. "a disc at z = 0.04 m from r = 0 to 0.06 m".
    fn describe(&self) -> String {
        let f = |x: f64| format!("{x:.6e}");
        if self.pieces.len() == 1 {
            let s = &self.pieces[0];
            let b = s.at(s.m0 + s.len);
            let (r0, r1) = (s.a[0].min(b[0]), s.a[0].max(b[0]));
            let (z0, z1) = (s.a[1].min(b[1]), s.a[1].max(b[1]));
            return if s.d[0].abs() < 1e-9 {
                format!(
                    "a cylinder at r = {} m from z = {} to {} m",
                    f(s.a[0]),
                    f(z0),
                    f(z1)
                )
            } else if s.d[1].abs() < 1e-9 {
                format!(
                    "a disc at z = {} m from r = {} to {} m",
                    f(s.a[1]),
                    f(r0),
                    f(r1)
                )
            } else {
                format!(
                    "a cone from (r, z) = ({}, {}) to ({}, {}) m",
                    f(s.a[0]),
                    f(s.a[1]),
                    f(b[0]),
                    f(b[1])
                )
            };
        }
        format!(
            "a surface of revolution of {} straight pieces in {} part{}{}",
            self.pieces.len(),
            self.chains,
            if self.chains == 1 { "" } else { "s" },
            if self.axis_m.is_empty() {
                String::new()
            } else {
                format!(", {} through the axis", self.axis_m.len())
            }
        )
    }
}

/// The meridian of the faces whose nodes are at `faces` in (r, z), and the
/// piece of every face.
fn meridian(faces: &[Vec<[f64; 2]>], tol: f64) -> Result<(Meridian, Vec<usize>), String> {
    // the line of every face: its two farthest nodes
    struct Line {
        d: [f64; 2],
        c: f64,
        s: [f64; 2],
    }
    let mut lines = Vec::with_capacity(faces.len());
    for (k, pts) in faces.iter().enumerate() {
        let (mut best, mut ends) = (-1.0, (0, 0));
        for i in 0..pts.len() {
            for j in i + 1..pts.len() {
                let l = (pts[i][0] - pts[j][0]).hypot(pts[i][1] - pts[j][1]);
                if l > best {
                    best = l;
                    ends = (i, j);
                }
            }
        }
        if best <= tol {
            return Err(format!(
                "the interface face {k} lies on one circle about the axis (a face of a surface \
                 of revolution must extend along the meridian)"
            ));
        }
        let (p, q) = (pts[ends.0], pts[ends.1]);
        let mut d = [(q[0] - p[0]) / best, (q[1] - p[1]) / best];
        // one direction per line: z up, or r out for a disc
        if d[1] < -1e-9 || (d[1].abs() <= 1e-9 && d[0] < 0.0) {
            d = [-d[0], -d[1]];
        }
        let c = d[0] * p[1] - d[1] * p[0];
        for x in pts {
            if (d[0] * x[1] - d[1] * x[0] - c).abs() > tol {
                return Err(format!(
                    "the nodes of the interface face {k} do not lie on one straight piece of the \
                     meridian: the interface must be a surface of revolution about the axis \
                     made of cylinders, cones and discs, with a ring of nodes at every kink"
                ));
            }
        }
        let (s0, s1) = (d[0] * p[0] + d[1] * p[1], d[0] * q[0] + d[1] * q[1]);
        lines.push(Line {
            d,
            c,
            s: [s0.min(s1), s0.max(s1)],
        });
    }
    // collinear faces: a group per line
    let mut groups: Vec<(usize, Vec<usize>)> = Vec::new();
    for (k, l) in lines.iter().enumerate() {
        let g = groups.iter_mut().find(|(r, _)| {
            let m = &lines[*r];
            let sin = (m.d[0] * l.d[1] - m.d[1] * l.d[0]).abs();
            let mid = 0.5 * (l.s[0] + l.s[1]);
            sin * mid.abs().max(1.0) <= tol.max(1e-9) && (m.c - l.c).abs() <= tol
        });
        match g {
            Some((_, v)) => v.push(k),
            None => groups.push((k, vec![k])),
        }
    }
    // the pieces: overlapping intervals of a line
    let mut raw: Vec<([f64; 2], [f64; 2])> = Vec::new();
    let mut face_piece = vec![usize::MAX; faces.len()];
    for (r, members) in &groups {
        let rep = &lines[*r];
        let n = [-rep.d[1], rep.d[0]];
        let proj = |k: usize| {
            let l = &lines[k];
            // the interval along the representative's direction
            let sgn = if rep.d[0] * l.d[0] + rep.d[1] * l.d[1] >= 0.0 {
                1.0
            } else {
                -1.0
            };
            let (a, b) = (sgn * l.s[0], sgn * l.s[1]);
            [a.min(b), a.max(b)]
        };
        let mut order: Vec<usize> = members.clone();
        order.sort_by(|&a, &b| proj(a)[0].total_cmp(&proj(b)[0]));
        let mut cur: Option<[f64; 2]> = None;
        let mut cur_faces: Vec<usize> = Vec::new();
        let mut flush = |iv: [f64; 2], fs: &mut Vec<usize>, raw: &mut Vec<([f64; 2], [f64; 2])>| {
            let pt = |s: f64| [s * rep.d[0] + rep.c * n[0], s * rep.d[1] + rep.c * n[1]];
            let id = raw.len();
            raw.push((pt(iv[0]), pt(iv[1])));
            for &f in fs.iter() {
                face_piece[f] = id;
            }
            fs.clear();
        };
        for k in order {
            let iv = proj(k);
            match cur {
                Some(ref mut c) if iv[0] <= c[1] + tol => c[1] = c[1].max(iv[1]),
                Some(c) => {
                    flush(c, &mut cur_faces, &mut raw);
                    cur = Some(iv);
                }
                None => cur = Some(iv),
            }
            cur_faces.push(k);
        }
        if let Some(c) = cur {
            flush(c, &mut cur_faces, &mut raw);
        }
    }
    // the ends of the pieces, joined where they meet
    let mut ends: Vec<[f64; 2]> = Vec::new();
    let mut end_of = vec![[0usize; 2]; raw.len()];
    for (i, (a, b)) in raw.iter().enumerate() {
        for (j, p) in [a, b].into_iter().enumerate() {
            let e = match ends
                .iter()
                .position(|q| (q[0] - p[0]).hypot(q[1] - p[1]) <= 10.0 * tol)
            {
                Some(e) => e,
                None => {
                    ends.push(*p);
                    ends.len() - 1
                }
            };
            end_of[i][j] = e;
        }
    }
    let mut at: Vec<Vec<usize>> = vec![Vec::new(); ends.len()];
    for (i, e) in end_of.iter().enumerate() {
        at[e[0]].push(i);
        at[e[1]].push(i);
    }
    if let Some(e) = at.iter().position(|v| v.len() > 2) {
        return Err(format!(
            "the meridian of the interface branches at (r, z) = ({:.6e}, {:.6e}) m: one \
             interface per rotor must be one surface (without branches)",
            ends[e][0], ends[e][1]
        ));
    }
    // chains from their free ends: the end nearer the axis first, then the
    // lower z
    let mut starts: Vec<usize> = (0..ends.len()).filter(|&e| at[e].len() == 1).collect();
    starts.sort_by(|&a, &b| {
        let (p, q) = (ends[a], ends[b]);
        p[0].total_cmp(&q[0]).then(p[1].total_cmp(&q[1]))
    });
    let mut used = vec![false; raw.len()];
    let mut pieces: Vec<Segment> = Vec::new();
    let mut piece_of_raw = vec![usize::MAX; raw.len()];
    let mut axis_m = Vec::new();
    let mut chains = 0;
    let span = raw
        .iter()
        .map(|(a, b)| (a[0] - b[0]).hypot(a[1] - b[1]))
        .sum::<f64>();
    let mut m = 0.0;
    for &s in &starts {
        if at[s].iter().all(|&i| used[i]) {
            continue;
        }
        chains += 1;
        if ends[s][0] <= tol {
            axis_m.push(m);
        }
        let mut e = s;
        while let Some(&i) = at[e].iter().find(|&&i| !used[i]) {
            used[i] = true;
            let (fwd, next) = if end_of[i][0] == e {
                (true, end_of[i][1])
            } else {
                (false, end_of[i][0])
            };
            let (a, b) = if fwd { raw[i] } else { (raw[i].1, raw[i].0) };
            let len = (b[0] - a[0]).hypot(b[1] - a[1]);
            piece_of_raw[i] = pieces.len();
            pieces.push(Segment {
                a: ends[e],
                d: [(b[0] - a[0]) / len, (b[1] - a[1]) / len],
                len,
                m0: m,
            });
            m += len;
            e = next;
        }
        if ends[e][0] <= tol {
            axis_m.push(m);
        }
        // the next chain starts after a gap
        m += span;
    }
    if used.iter().any(|u| !u) {
        return Err(
            "the meridian of the interface is a closed loop (a ring that does not reach a wall \
             or the axis): not supported"
                .into(),
        );
    }
    let face_piece = face_piece.into_iter().map(|r| piece_of_raw[r]).collect();
    Ok((
        Meridian {
            pieces,
            chains,
            axis_m,
        },
        face_piece,
    ))
}

/// The imprint of a base mesh: prepared once, assembled for each rotor angle.
#[derive(Debug, Clone)]
pub struct Interface {
    origin: [f64; 3],
    axis: [f64; 3],
    e1: [f64; 3],
    e2: [f64; 3],
    meridian: Meridian,
    /// Length scale of the surface (m): the largest distance of a node from
    /// the axis.
    size: f64,
    rotor_cell: Vec<bool>,
    rotor_node: Vec<bool>,
    rotor: Vec<IfFace>,
    stator: Vec<IfFace>,
    /// Indices of the two patches in `Mesh::boundaries`.
    patches: [usize; 2],
    /// Stator faces by angular bin.
    bins: Vec<Vec<u32>>,
    /// The other faces of the base mesh that share an edge with an
    /// interface face (they get the corners of the pieces on that edge).
    touching: Vec<u32>,
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
        if a == b {
            continue;
        }
        // a corner within round-off of the line lies on it (no spurious
        // crossing of two nearly coincident edges)
        let tol = 1e-12 * (b[0] - a[0]).hypot(b[1] - a[1]);
        let side = |p: [f64; 2]| {
            let s = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
            if s.abs() <= tol { 0.0 } else { s }
        };
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

/// The polygon without repeated corners (tolerances in θ and in m).
fn dedupe(p: Vec<[f64; 2]>, tol: [f64; 2]) -> Vec<[f64; 2]> {
    let same =
        |a: [f64; 2], b: [f64; 2]| (a[0] - b[0]).abs() <= tol[0] && (a[1] - b[1]).abs() <= tol[1];
    let mut out: Vec<[f64; 2]> = Vec::with_capacity(p.len());
    for q in p {
        if out.last().is_none_or(|&l| !same(l, q)) {
            out.push(q);
        }
    }
    while out.len() > 1 && same(out[0], out[out.len() - 1]) {
        out.pop();
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

/// The key of the edge between two nodes.
fn edge_key(a: u32, b: u32) -> (u32, u32) {
    (a.min(b), a.max(b))
}

/// An intersection polygon in (θ, m) of the stator frame.
struct Cut {
    rotor: u32,
    stator: u32,
    poly: Vec<[f64; 2]>,
    /// The rotor shifted into the stator frame by this angle.
    shift: f64,
}

impl Interface {
    /// Prepares the imprint of `mesh`: checks that the two boundaries lie on
    /// one surface of revolution about the axis, cover it alike, and separate
    /// a rotor from a stator that share no node.
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
        // the faces of both sides, owner checked
        let mut faces: Vec<(usize, u32)> = Vec::new();
        for (side, &pi) in patches.iter().enumerate() {
            let b = &mesh.boundaries[pi];
            for f in b.range() {
                let cell = mesh.owner[f];
                let want_rotor = side == 0;
                if rotor_cell[cell as usize] != want_rotor {
                    return Err(format!(
                        "the face {f} of the boundary `{}` belongs to a cell of the {}",
                        b.name,
                        if want_rotor { "stator" } else { "rotor" }
                    ));
                }
                faces.push((f, cell));
            }
        }
        if mesh.boundaries[patches[0]].len == 0 || mesh.boundaries[patches[1]].len == 0 {
            return Err("an interface boundary has no faces".into());
        }
        // cylindrical coordinates (θ, r, z) of the nodes
        let cyl = |p: [f64; 3]| -> [f64; 3] {
            let v = sub(p, spec.origin);
            let (x, y, z) = (dot(v, e1), dot(v, e2), dot(v, axis));
            [y.atan2(x), x.hypot(y), z]
        };
        let coords: Vec<Vec<[f64; 3]>> = faces
            .iter()
            .map(|&(f, _)| {
                mesh.face(f)
                    .iter()
                    .map(|&p| cyl(mesh.points[p as usize]))
                    .collect()
            })
            .collect();
        let size = coords.iter().flatten().map(|c| c[1]).fold(0.0, f64::max);
        if size <= 0.0 {
            return Err("the interface lies on the axis".into());
        }
        let tol = ON_SURFACE * size;
        let rz: Vec<Vec<[f64; 2]>> = coords
            .iter()
            .map(|v| v.iter().map(|c| [c[1], c[2]]).collect())
            .collect();
        let (meridian, face_piece) = meridian(&rz, tol)?;
        // the unrolled faces
        let n_r = mesh.boundaries[patches[0]].len;
        let mut rotor = Vec::with_capacity(n_r);
        let mut stator = Vec::with_capacity(faces.len() - n_r);
        for (k, &(f, cell)) in faces.iter().enumerate() {
            let name = &mesh.boundaries[patches[usize::from(k >= n_r)]].name;
            let nodes = mesh.face(f);
            let seg = &meridian.pieces[face_piece[k]];
            let c = &coords[k];
            let on_axis: Vec<bool> = c.iter().map(|x| x[1] <= tol).collect();
            let Some(first) = on_axis.iter().position(|&a| !a) else {
                return Err(format!("the face {f} of `{name}` lies on the axis"));
            };
            let nn = nodes.len();
            if (0..nn).any(|i| on_axis[i] && on_axis[(i + 1) % nn]) {
                return Err(format!(
                    "the face {f} of `{name}` has an edge on the axis (a face of the surface \
                     meets the axis in one node)"
                ));
            }
            // θ unwrapped about the first node off the axis
            let t0 = c[first][0];
            let theta: Vec<f64> = c.iter().map(|x| t0 + wrap(x[0] - t0)).collect();
            let mut poly: Vec<[f64; 2]> = Vec::with_capacity(nn + 1);
            let mut corner_node: Vec<u32> = Vec::with_capacity(nn + 1);
            for i in 0..nn {
                // m in units of the size of the surface
                let m = seg.m_of([c[i][1], c[i][2]]) / size;
                if on_axis[i] {
                    // the polar chart: the edge m = m_axis between the angles
                    // of the two neighbours
                    let (prev, next) = ((i + nn - 1) % nn, (i + 1) % nn);
                    let m_axis = meridian
                        .axis_m
                        .iter()
                        .map(|a| a / size)
                        .min_by(|a, b| (a - m).abs().total_cmp(&(b - m).abs()))
                        .filter(|a| (a - m).abs() <= 10.0 * ON_SURFACE)
                        .ok_or_else(|| {
                            format!(
                                "the node {} of the face {f} of `{name}` lies on the axis away \
                                 from an end of the meridian",
                                nodes[i]
                            )
                        })?;
                    poly.push([theta[prev], m_axis]);
                    poly.push([theta[next], m_axis]);
                    corner_node.extend([nodes[i], nodes[i]]);
                } else {
                    poly.push([theta[i], m]);
                    corner_node.push(nodes[i]);
                }
            }
            let a2 = area2(&poly);
            if a2.abs() <= 0.0 {
                return Err(format!("the face {f} of `{name}` has no area"));
            }
            let sign = if a2 >= 0.0 { 1.0 } else { -1.0 };
            if sign < 0.0 {
                poly.reverse();
                corner_node.reverse();
            }
            let lo = poly.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
            let hi = poly.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
            if hi - lo >= 0.99 * PI {
                return Err(format!(
                    "the face {f} of `{name}` spans {:.1} degrees around the axis (at most 178 \
                     are supported)",
                    (hi - lo).to_degrees()
                ));
            }
            let mlo = poly.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min);
            let mhi = poly.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max);
            let face = IfFace {
                face: f as u32,
                cell,
                poly,
                corner_node,
                sign,
                lo,
                hi,
                mlo,
                mhi,
                area: 0.5 * a2.abs(),
            };
            if k < n_r {
                rotor.push(face);
            } else {
                stator.push(face);
            }
        }
        let (ar, as_): (f64, f64) = (
            rotor.iter().map(|f| f.area).sum(),
            stator.iter().map(|f| f.area).sum(),
        );
        if (ar - as_).abs() > 1e-8 * ar.max(as_) {
            return Err(format!(
                "the two interface boundaries do not cover the same surface (areas {ar:.8e} \
                 and {as_:.8e} in (θ, m)): a partial ring is not supported"
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
        // the faces next to the interface: an edge in common with it
        let mut edges: HashSet<(u32, u32)> = HashSet::new();
        for f in rotor.iter().chain(&stator) {
            let k = f.corner_node.len();
            for i in 0..k {
                let (a, b) = (f.corner_node[i], f.corner_node[(i + 1) % k]);
                if a != b {
                    edges.insert(edge_key(a, b));
                }
            }
        }
        let touching: Vec<u32> = (0..mesh.n_faces())
            .filter(|&f| {
                let in_if = mesh.boundaries[patches[0]].range().contains(&f)
                    || mesh.boundaries[patches[1]].range().contains(&f);
                let nodes = mesh.face(f);
                let k = nodes.len();
                !in_if && (0..k).any(|i| edges.contains(&edge_key(nodes[i], nodes[(i + 1) % k])))
            })
            .map(|f| f as u32)
            .collect();
        Ok(Self {
            origin: spec.origin,
            axis,
            e1,
            e2,
            meridian,
            size,
            rotor_cell,
            rotor_node,
            rotor,
            stator,
            patches,
            bins,
            touching,
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

    /// The shape of the surface in words ("a cylinder at r = … m", "a disc
    /// at z = … m", "a cone …", "a surface of revolution of n pieces").
    pub fn surface(&self) -> String {
        self.meridian.describe()
    }

    /// The straight pieces of the meridian: (r, z) of their two ends (m), in
    /// the order of the arc length.
    pub fn meridian(&self) -> Vec<([f64; 2], [f64; 2])> {
        self.meridian
            .pieces
            .iter()
            .map(|s| (s.a, s.at(s.m0 + s.len)))
            .collect()
    }

    /// Indices of the two interface boundaries in the base mesh.
    pub fn patches(&self) -> [usize; 2] {
        self.patches
    }

    /// The point of the surface at (θ, u), u = m / size.
    fn point(&self, t: f64, u: f64) -> [f64; 3] {
        let [r, z] = self.meridian.rz(u * self.size);
        let (s, c) = t.sin_cos();
        add(
            add(
                self.origin,
                scale(add(scale(self.e1, c), scale(self.e2, s)), r),
            ),
            scale(self.axis, z),
        )
    }

    /// Whether u = m / size is an end of the meridian on the axis.
    fn on_axis(&self, u: f64) -> bool {
        let tol = 1e-12 * (1.0 + u.abs());
        self.meridian
            .axis_m
            .iter()
            .any(|a| (a / self.size - u).abs() <= tol)
    }

    /// The intersections of the rotor faces, turned by `angle`, with the
    /// stator faces.
    fn intersections(&self, angle: f64) -> Vec<Cut> {
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
                if r.mhi <= s.mlo || r.mlo >= s.mhi {
                    continue;
                }
                let shift = angle + TAU * ((0.5 * (s.lo + s.hi) - r_mid) / TAU).round();
                if r.hi + shift <= s.lo || r.lo + shift >= s.hi {
                    continue;
                }
                let shifted: Vec<[f64; 2]> = r.poly.iter().map(|p| [p[0] + shift, p[1]]).collect();
                let tol = [
                    1e-13 * (1.0 + s.lo.abs().max(s.hi.abs())),
                    1e-13 * (1.0 + s.mlo.abs().max(s.mhi.abs())),
                ];
                let c = dedupe(clip(&shifted, &s.poly), tol);
                if c.len() < 3 {
                    continue;
                }
                let a = 0.5 * area2(&c);
                if a <= SLIVER * r.area.min(s.area) {
                    continue;
                }
                out.push(Cut {
                    rotor: ri as u32,
                    stator: si,
                    poly: c,
                    shift,
                });
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

/// The corners of `poly` (in θ, u) that lie on an edge of `face` (shifted by
/// `shift` in θ): (edge key, parameter from the smaller node, corner index);
/// the ends of the edge have the parameters 0 and 1.
fn on_edges(face: &IfFace, shift: f64, poly: &[[f64; 2]], out: &mut Vec<((u32, u32), f64, usize)>) {
    let k = face.poly.len();
    for i in 0..k {
        let (na, nb) = (face.corner_node[i], face.corner_node[(i + 1) % k]);
        if na == nb {
            continue;
        }
        let (a, b) = (face.poly[i], face.poly[(i + 1) % k]);
        let (ax, ay) = (a[0] + shift, a[1]);
        let (ex, ey) = (b[0] - a[0], b[1] - a[1]);
        let l2 = ex * ex + ey * ey;
        if l2 <= 0.0 {
            continue;
        }
        for (j, q) in poly.iter().enumerate() {
            let (qx, qy) = (q[0] - ax, q[1] - ay);
            let t = (qx * ex + qy * ey) / l2;
            if !(-1e-9..=1.0 + 1e-9).contains(&t) {
                continue;
            }
            let d = (qx * ey - qy * ex).abs() / l2.sqrt();
            if d > 1e-10 * l2.sqrt() {
                continue;
            }
            let t = t.clamp(0.0, 1.0);
            let (key, t) = if na < nb {
                ((na, nb), t)
            } else {
                ((nb, na), 1.0 - t)
            };
            out.push((key, t, j));
        }
    }
}

/// Whether a parameter on an edge lies inside it (not at a node).
fn inside(t: f64) -> bool {
    t > 1e-9 && t < 1.0 - 1e-9
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
    // the corners of the polygons are new nodes; those inside an edge of a
    // rotor or a stator face also go into the other faces of that edge
    type Hits = Vec<((u32, u32), f64, usize)>;
    let mut piece_nodes: Vec<(Vec<u32>, Hits, bool)> = Vec::new();
    let mut pieces = Vec::new();
    let mut inserted: HashMap<(u32, u32), Vec<(f64, u32)>> = HashMap::new();
    for (ifc, &angle) in parts.iter().zip(angles) {
        for cut in &ifc.intersections(angle) {
            let (r, s) = (
                &ifc.rotor[cut.rotor as usize],
                &ifc.stator[cut.stator as usize],
            );
            // corners on the axis are one node
            let mut poly: Vec<[f64; 2]> = Vec::with_capacity(cut.poly.len());
            for &q in &cut.poly {
                let dup = poly
                    .last()
                    .is_some_and(|l: &[f64; 2]| ifc.on_axis(l[1]) && ifc.on_axis(q[1]));
                if !dup {
                    poly.push(q);
                }
            }
            while poly.len() > 1 && ifc.on_axis(poly[0][1]) && ifc.on_axis(poly[poly.len() - 1][1])
            {
                poly.pop();
            }
            if poly.len() < 3 {
                continue;
            }
            // θ of the polygon is in the frame of the stator
            let ids: Vec<u32> = poly
                .iter()
                .map(|q| {
                    let id = points.len() as u32;
                    points.push(ifc.point(q[0], q[1]));
                    id
                })
                .collect();
            let mut hits = Vec::new();
            on_edges(r, cut.shift, &poly, &mut hits);
            on_edges(s, 0.0, &poly, &mut hits);
            for &(key, t, j) in &hits {
                inserted.entry(key).or_default().push((t, ids[j]));
            }
            piece_nodes.push((ids, hits, r.sign < 0.0));
            pieces.push(Piece {
                rotor_face: r.face,
                stator_face: s.face,
                rotor_cell: r.cell,
                stator_cell: s.cell,
                area: 0.0,
            });
        }
    }
    // one node per point of an edge (the pieces on both sides of it give it
    // twice)
    for v in inserted.values_mut() {
        v.sort_by(|a, b| a.0.total_cmp(&b.0));
        v.dedup_by(|b, a| (b.0 - a.0).abs() <= 1e-9);
        v.retain(|x| inside(x.0));
    }
    // the pieces next to each other on an edge of a face (two faces of one
    // side that share an edge) get each other's corners on it
    let piece_nodes: Vec<Vec<u32>> = piece_nodes
        .into_iter()
        .map(|(ids, hits, reverse)| {
            let k = ids.len();
            let mut out = Vec::with_capacity(k + 2);
            for (j, &id) in ids.iter().enumerate() {
                out.push(id);
                let jn = (j + 1) % k;
                for &(key, ta, _) in hits.iter().filter(|h| h.2 == j) {
                    let Some(&(_, tb, _)) = hits.iter().find(|h| h.2 == jn && h.0 == key) else {
                        continue;
                    };
                    let Some(v) = inserted.get(&key) else {
                        continue;
                    };
                    let (lo, hi) = (ta.min(tb) + 1e-9, ta.max(tb) - 1e-9);
                    let between = v.iter().filter(|x| x.0 > lo && x.0 < hi);
                    if ta < tb {
                        out.extend(between.map(|x| x.1));
                    } else {
                        out.extend(between.rev().map(|x| x.1));
                    }
                }
            }
            // the orientation of the rotor face: out of the rotor cell
            if reverse {
                out.reverse();
            }
            out
        })
        .collect();
    let mut touching = vec![false; mesh.n_faces()];
    for ifc in parts {
        for &f in &ifc.touching {
            touching[f as usize] = true;
        }
    }
    let mut buf: Vec<u32> = Vec::new();
    let mut nodes_of = |f: usize| -> Vec<u32> {
        let nodes = mesh.face(f);
        if !touching[f] {
            return nodes.to_vec();
        }
        buf.clear();
        let k = nodes.len();
        for i in 0..k {
            let (a, b) = (nodes[i], nodes[(i + 1) % k]);
            buf.push(a);
            if let Some(v) = inserted.get(&edge_key(a, b)) {
                if a < b {
                    buf.extend(v.iter().map(|x| x.1));
                } else {
                    buf.extend(v.iter().rev().map(|x| x.1));
                }
            }
        }
        buf.clone()
    };
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
        b.add_interior(&nodes_of(f), mesh.owner[f], mesh.neighbour[f]);
        origin.push(FaceOrigin::Base(f as u32));
    }
    for p in &mesh.boundaries {
        if skip(&p.name) {
            continue;
        }
        for f in p.range() {
            b.add_boundary(&nodes_of(f), mesh.owner[f], &p.name);
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

    /// The meshes of revolution of D186: a stage on a hub (disc, cylinder,
    /// disc), its cone, a pot through the axis, and a rotor row between two
    /// discs (two chains).
    fn revolved() -> Vec<(&'static str, Mesh)> {
        use crate::revolved::RevolvedSpec;
        let small = |s: RevolvedSpec| {
            RevolvedSpec {
                n_r: s.n_r.iter().map(|_| 2).collect(),
                n_z: s.n_z.iter().map(|_| 2).collect(),
                n_theta: 12,
                n_theta_stator: 17,
                ..s
            }
            .generate()
            .unwrap()
        };
        vec![
            ("stage", small(RevolvedSpec::default())),
            (
                "cone",
                small(RevolvedSpec {
                    cone: 25.0,
                    ..Default::default()
                }),
            ),
            (
                "pot",
                small(RevolvedSpec {
                    r: vec![0.0, 0.04, 0.1],
                    z: vec![0.0, 0.04, 0.08, 0.15],
                    n_r: vec![1, 1],
                    n_z: vec![1, 1, 1],
                    rotor: vec![[0, 1]],
                    ..Default::default()
                }),
            ),
            (
                "row",
                small(RevolvedSpec {
                    r: vec![0.05, 0.1],
                    n_r: vec![1],
                    cone: -15.0,
                    ..Default::default()
                }),
            ),
        ]
    }

    #[test]
    fn surfaces_of_revolution_close_at_every_angle() {
        let want = [
            (
                "stage",
                "a surface of revolution of 3 straight pieces in 1 part",
            ),
            (
                "cone",
                "a surface of revolution of 3 straight pieces in 1 part",
            ),
            (
                "pot",
                "a surface of revolution of 3 straight pieces in 1 part, 2 through the axis",
            ),
            (
                "row",
                "a surface of revolution of 2 straight pieces in 2 parts",
            ),
        ];
        for ((name, m), (n2, text)) in revolved().into_iter().zip(want) {
            assert_eq!(name, n2);
            let ifc = Interface::new(&m, &spec()).unwrap_or_else(|e| panic!("{name}: {e}"));
            assert_eq!(ifc.surface(), text, "{name}");
            let base = Geometry::compute(&m);
            let rotor_area: f64 = m.boundaries[ifc.patches[0]]
                .range()
                .map(|f| norm(base.face_area[f]))
                .sum();
            for angle in [0.0, 0.05, 0.3333, 1.0, 2.5, -0.7, 7.0] {
                let a = ifc.assemble(&m, angle).unwrap();
                let geo = ifc.geometry(&a, &base);
                let defect = ifc.closure(&a, &geo);
                assert!(defect < 1e-12, "{name} angle {angle}: closure {defect:e}");
                // the pieces cover the rotor's faces (up to the chords)
                let pieces: f64 = a.pieces.iter().map(|p| p.area).sum();
                assert!(
                    (pieces - rotor_area).abs() < 5e-2 * rotor_area,
                    "{name} angle {angle}: {pieces} {rotor_area}"
                );
                // every rotor face has pieces of about its area
                let mut per = vec![0.0; m.n_faces()];
                for p in &a.pieces {
                    per[p.rotor_face as usize] += p.area;
                }
                for f in m.boundaries[ifc.patches[0]].range() {
                    let s = norm(base.face_area[f]);
                    assert!(
                        (per[f] - s).abs() < 1e-1 * s,
                        "{name} face {f}: {} {s}",
                        per[f]
                    );
                }
                assert!(
                    a.mesh
                        .boundaries
                        .iter()
                        .all(|b| !b.name.contains("interface"))
                );
            }
        }
    }

    #[test]
    fn a_cylinder_keeps_its_description() {
        let m = annulus(12, 18, 2);
        let ifc = Interface::new(&m, &spec()).unwrap();
        assert!(
            ifc.surface().starts_with("a cylinder at r = "),
            "{}",
            ifc.surface()
        );
        let mer = ifc.meridian();
        assert_eq!(mer.len(), 1);
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
        assert!(e.contains("surface of revolution about the axis"), "{e}");
        let mut s = spec();
        s.axis = [0.0; 3];
        assert!(Interface::new(&m, &s).unwrap_err().contains("no direction"));
    }
}

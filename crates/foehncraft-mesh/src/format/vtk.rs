//! Reader and writer for VTK unstructured grids (D192): legacy `.vtk` and
//! XML `.vtu`.
//!
//! Written from the public descriptions of the two file formats (the VTK file
//! formats guide).
//!
//! **Reading.** Legacy files in ASCII or BINARY (big-endian) form, also the
//! 5.1 layout with `OFFSETS` and `CONNECTIVITY`; XML files with the data in
//! ASCII, inline base64 or appended (raw or base64) form, little- or
//! big-endian, 32- or 64-bit headers, uncompressed or zlib-compressed (the
//! default of ParaView; inflated by the decoder of this crate). Cell types
//! tetrahedron (10), hexahedron (12), voxel (11), wedge (13), pyramid (14) and
//! polyhedron (42) are the cells; triangles (5), polygons (7) and
//! quadrilaterals (9) are boundary faces; vertices and lines are ignored; a
//! higher-order or 2D cell is refused (a 2D VTK mesh: extrude it in the
//! mesher).
//!
//! **Names.** A grid has numbers where Foehncraft has names. The region
//! numbers are one scalar cell-data array, the first found of `region`,
//! `gmsh:physical`, `CellEntityIds`, `PatchID`, `patch`, `boundary`, `zone`,
//! `MaterialIds` (also in other capitals): the volume cells give the cell
//! zones, the surface cells the boundaries. Names for the numbers come from
//! the string arrays `zone_names` and `boundary_names` of the field data
//! (entries `id=name`, as [`write_vtu`] writes them); else a zone is `fluid`
//! (one zone) or `zone<id>`, a boundary `boundary<id>`. Surface cells without
//! a region array form the boundary `default`; boundary faces that no surface
//! cell names go to `default` as for every format.
//!
//! **Writing.** [`write_vtu`] writes the mesh as a `.vtu` file: every cell a
//! polyhedron, every boundary face a polygon, the region array and the names
//! as above, ASCII or appended binary. The file reads back as the same mesh
//! (also in ParaView).

use super::{ElementMesh, MeshFile, MeshFormat, ReadOptions, inflate};
use crate::mesh::Mesh;
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

/// The VTK format in the registry of [`crate::format`].
#[derive(Debug, Clone, Copy)]
pub struct Vtk;

impl MeshFormat for Vtk {
    fn name(&self) -> &'static str {
        "vtk"
    }

    fn title(&self) -> &'static str {
        "VTK unstructured grid, legacy .vtk and XML .vtu"
    }

    fn extensions(&self) -> &'static [&'static str] {
        &["vtk", "vtu"]
    }

    fn sniff(&self, head: &[u8]) -> bool {
        let start = head
            .iter()
            .position(|c| !c.is_ascii_whitespace())
            .unwrap_or(0);
        let h = &head[start..];
        h.starts_with(b"# vtk DataFile") || h.starts_with(b"<?xml") || h.starts_with(b"<VTKFile")
    }

    fn parse(&self, bytes: &[u8], _options: &ReadOptions) -> Result<MeshFile, String> {
        let probe = &bytes[..bytes.len().min(4096)];
        let (grid, how) = if probe.starts_with(b"# vtk") {
            parse_legacy(bytes)?
        } else {
            parse_xml(bytes)?
        };
        grid.into_mesh_file(how)
    }
}

/// What a file holds, before it becomes a mesh.
#[derive(Default)]
struct Grid {
    points: Vec<[f64; 3]>,
    /// VTK cell type and node list (for a polyhedron: the nodes of its faces
    /// are in `poly`).
    cells: Vec<(u8, Vec<u32>)>,
    poly: HashMap<usize, Vec<Vec<u32>>>,
    /// Scalar cell-data arrays.
    arrays: Vec<(String, Vec<f64>)>,
    /// String arrays of the field data.
    names: HashMap<String, Vec<String>>,
}

const REGION_NAMES: [&str; 12] = [
    "region",
    "gmsh:physical",
    "cellentityids",
    "regionid",
    "regionids",
    "patchid",
    "patch",
    "boundary",
    "boundaryid",
    "zone",
    "zoneid",
    "materialids",
];

fn id_names(list: Option<&Vec<String>>) -> HashMap<i64, String> {
    list.into_iter()
        .flatten()
        .filter_map(|e| {
            let (i, n) = e.split_once('=')?;
            Some((i.trim().parse().ok()?, n.trim().to_string()))
        })
        .collect()
}

impl Grid {
    fn into_mesh_file(self, how: String) -> Result<MeshFile, String> {
        let region: Option<&Vec<f64>> = REGION_NAMES.iter().find_map(|want| {
            self.arrays
                .iter()
                .find(|(n, _)| n.eq_ignore_ascii_case(want))
                .map(|(_, v)| v)
        });
        let id = |i: usize| region.and_then(|r| r.get(i)).map(|v| v.round() as i64);
        let zone_names = id_names(self.names.get("zone_names"));
        let boundary_names = id_names(self.names.get("boundary_names"));

        // the zones and boundaries, in the order of their numbers
        let mut zone_ids: BTreeMap<i64, ()> = BTreeMap::new();
        let mut bound_ids: BTreeMap<i64, ()> = BTreeMap::new();
        for (i, (t, _)) in self.cells.iter().enumerate() {
            match t {
                10 | 11 | 12 | 13 | 14 | 42 => {
                    zone_ids.insert(id(i).unwrap_or(0), ());
                }
                5 | 7 | 9 => {
                    bound_ids.insert(id(i).unwrap_or(i64::MIN), ());
                }
                _ => {}
            }
        }
        let mut em = ElementMesh {
            points: self.points.clone(),
            face_word: "surface cells",
            group_word: "region",
            ..Default::default()
        };
        let single = zone_ids.len() == 1;
        let zone_of: HashMap<i64, u32> = zone_ids
            .keys()
            .map(|&k| {
                let name = zone_names.get(&k).cloned().unwrap_or_else(|| {
                    if single {
                        "fluid".to_string()
                    } else {
                        format!("zone{k}")
                    }
                });
                (k, em.zone(&name))
            })
            .collect();
        let bound_of: HashMap<i64, u32> = bound_ids
            .keys()
            .map(|&k| {
                let name = boundary_names.get(&k).cloned().unwrap_or_else(|| {
                    if k == i64::MIN {
                        "default".to_string()
                    } else {
                        format!("boundary{k}")
                    }
                });
                (k, em.boundary(&name))
            })
            .collect();

        let n = self.points.len() as u32;
        let mut n_volume = 0usize;
        for (i, (t, nodes)) in self.cells.iter().enumerate() {
            if let Some(&v) = nodes.iter().find(|&&v| v >= n) {
                return Err(format!("cell {i} refers to the point {v}, beyond {n}"));
            }
            let need = |k: usize| -> Result<(), String> {
                if nodes.len() == k {
                    Ok(())
                } else {
                    Err(format!(
                        "cell {i} of type {t} has {} points, expected {k}",
                        nodes.len()
                    ))
                }
            };
            let z = *zone_of.get(&id(i).unwrap_or(0)).unwrap_or(&0);
            match t {
                10 => {
                    need(4)?;
                    em.add_cell(super::CellKind::Tetrahedron, nodes, z);
                }
                12 => {
                    need(8)?;
                    em.add_cell(super::CellKind::Hexahedron, nodes, z);
                }
                11 => {
                    need(8)?;
                    let v = [0, 1, 3, 2, 4, 5, 7, 6].map(|k| nodes[k]);
                    em.add_cell(super::CellKind::Hexahedron, &v, z);
                }
                13 => {
                    need(6)?;
                    em.add_cell(super::CellKind::Prism, nodes, z);
                }
                14 => {
                    need(5)?;
                    em.add_cell(super::CellKind::Pyramid, nodes, z);
                }
                42 => {
                    let faces = self
                        .poly
                        .get(&i)
                        .ok_or_else(|| format!("polyhedron {i} has no faces"))?;
                    em.add_polyhedron(faces.clone(), z);
                }
                5 | 7 | 9 => {
                    if nodes.len() < 3 {
                        return Err(format!("surface cell {i} has {} points", nodes.len()));
                    }
                    let b = *bound_of.get(&id(i).unwrap_or(i64::MIN)).unwrap_or(&0);
                    em.add_face(nodes, b);
                    continue;
                }
                1 | 2 | 3 | 4 | 6 | 8 => continue, // points, lines, strips, pixels
                other => {
                    return Err(format!(
                        "cell {i} has the type {other}: only linear tetrahedra, hexahedra, wedges, \
                         pyramids and polyhedra (and triangles and quadrilaterals as boundary \
                         faces) are read"
                    ));
                }
            }
            n_volume += 1;
        }
        if n_volume == 0 {
            return Err("the grid has no 3D cells (a 2D mesh: extrude it in the mesher)".into());
        }
        let built = em.into_mesh()?;
        Ok(MeshFile {
            mesh: built.mesh,
            reader: "vtk",
            format: how,
            warnings: built.warnings,
        })
    }
}

// ---------------------------------------------------------------------
// numbers in binary form

#[derive(Clone, Copy)]
enum Ty {
    I8,
    U8,
    I16,
    U16,
    I32,
    U32,
    I64,
    U64,
    F32,
    F64,
}

impl Ty {
    fn legacy(name: &str) -> Result<Self, String> {
        Ok(match name.to_ascii_lowercase().as_str() {
            "char" => Ty::I8,
            "unsigned_char" => Ty::U8,
            "short" => Ty::I16,
            "unsigned_short" => Ty::U16,
            "int" => Ty::I32,
            "unsigned_int" => Ty::U32,
            "long" | "vtktypeint64" | "vtkidtype" => Ty::I64,
            "unsigned_long" | "vtktypeuint64" => Ty::U64,
            "float" => Ty::F32,
            "double" => Ty::F64,
            other => return Err(format!("unknown data type '{other}'")),
        })
    }

    fn xml(name: &str) -> Result<Self, String> {
        Ok(match name {
            "Int8" => Ty::I8,
            "UInt8" => Ty::U8,
            "Int16" => Ty::I16,
            "UInt16" => Ty::U16,
            "Int32" => Ty::I32,
            "UInt32" => Ty::U32,
            "Int64" => Ty::I64,
            "UInt64" => Ty::U64,
            "Float32" => Ty::F32,
            "Float64" => Ty::F64,
            other => return Err(format!("unknown data type '{other}'")),
        })
    }

    fn size(self) -> usize {
        match self {
            Ty::I8 | Ty::U8 => 1,
            Ty::I16 | Ty::U16 => 2,
            Ty::I32 | Ty::U32 | Ty::F32 => 4,
            Ty::I64 | Ty::U64 | Ty::F64 => 8,
        }
    }

    fn values(self, bytes: &[u8], big: bool) -> Vec<f64> {
        let k = self.size();
        bytes
            .chunks_exact(k)
            .map(|c| {
                let mut a = [0u8; 8];
                a[..k].copy_from_slice(c);
                if big {
                    a[..k].reverse();
                }
                match self {
                    Ty::I8 => a[0] as i8 as f64,
                    Ty::U8 => a[0] as f64,
                    Ty::I16 => i16::from_le_bytes([a[0], a[1]]) as f64,
                    Ty::U16 => u16::from_le_bytes([a[0], a[1]]) as f64,
                    Ty::I32 => i32::from_le_bytes([a[0], a[1], a[2], a[3]]) as f64,
                    Ty::U32 => u32::from_le_bytes([a[0], a[1], a[2], a[3]]) as f64,
                    Ty::I64 => i64::from_le_bytes(a) as f64,
                    Ty::U64 => u64::from_le_bytes(a) as f64,
                    Ty::F32 => f32::from_le_bytes([a[0], a[1], a[2], a[3]]) as f64,
                    Ty::F64 => f64::from_le_bytes(a),
                }
            })
            .collect()
    }
}

fn to_u32(v: f64, what: &str) -> Result<u32, String> {
    if v >= 0.0 && v <= u32::MAX as f64 && v.fract() == 0.0 {
        Ok(v as u32)
    } else {
        Err(format!("{what}: {v} is not a valid index"))
    }
}

// ---------------------------------------------------------------------
// legacy files

struct Cursor<'a> {
    b: &'a [u8],
    p: usize,
}

impl<'a> Cursor<'a> {
    fn line(&mut self) -> &'a str {
        let s = self.p;
        while self.p < self.b.len() && self.b[self.p] != b'\n' {
            self.p += 1;
        }
        let l = &self.b[s..self.p];
        self.p = (self.p + 1).min(self.b.len());
        std::str::from_utf8(l).unwrap_or("").trim()
    }

    fn word(&mut self) -> Option<&'a str> {
        while self.p < self.b.len() && self.b[self.p].is_ascii_whitespace() {
            self.p += 1;
        }
        let s = self.p;
        while self.p < self.b.len() && !self.b[self.p].is_ascii_whitespace() {
            self.p += 1;
        }
        (s < self.p).then(|| std::str::from_utf8(&self.b[s..self.p]).unwrap_or(""))
    }

    fn need(&mut self, what: &str) -> Result<&'a str, String> {
        self.word()
            .ok_or_else(|| format!("the file ends where {what} is expected"))
    }

    fn number(&mut self, what: &str) -> Result<usize, String> {
        let w = self.need(what)?;
        w.parse()
            .map_err(|_| format!("{what}: '{w}' is not a number"))
    }

    fn values(&mut self, n: usize, ty: Ty, binary: bool, what: &str) -> Result<Vec<f64>, String> {
        if binary {
            // one line break follows the header of binary data
            if self.b.get(self.p) == Some(&b'\r') {
                self.p += 1;
            }
            if self.b.get(self.p) == Some(&b'\n') {
                self.p += 1;
            }
            let len = n * ty.size();
            let data = self
                .b
                .get(self.p..self.p + len)
                .ok_or_else(|| format!("{what}: the binary data ends too early"))?;
            self.p += len;
            Ok(ty.values(data, true))
        } else {
            (0..n)
                .map(|_| {
                    let w = self.need(what)?;
                    w.parse::<f64>()
                        .map_err(|_| format!("{what}: '{w}' is not a number"))
                })
                .collect()
        }
    }
}

fn parse_legacy(bytes: &[u8]) -> Result<(Grid, String), String> {
    let mut c = Cursor { b: bytes, p: 0 };
    let first = c.line();
    let version = first
        .rsplit(' ')
        .next()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(0.0);
    c.line(); // title
    let mode = c.line().to_ascii_uppercase();
    let binary = match mode.as_str() {
        "ASCII" => false,
        "BINARY" => true,
        other => return Err(format!("the file is '{other}', expected ASCII or BINARY")),
    };
    let mut g = Grid::default();
    let mut cell_data = false;
    let mut n_cells = 0usize;
    let mut cell_types: Vec<f64> = Vec::new();
    let mut raw_cells: Vec<Vec<u32>> = Vec::new();
    let mut count_data = 0usize;
    while let Some(w) = c.word() {
        match w.to_ascii_uppercase().as_str() {
            "DATASET" => {
                let k = c.need("the dataset type")?;
                if !k.eq_ignore_ascii_case("UNSTRUCTURED_GRID") {
                    return Err(format!(
                        "the dataset is {k}: only UNSTRUCTURED_GRID is read"
                    ));
                }
            }
            "POINTS" => {
                let n = c.number("the number of points")?;
                let ty = Ty::legacy(c.need("the point type")?)?;
                let v = c.values(3 * n, ty, binary, "POINTS")?;
                g.points = v.as_chunks::<3>().0.to_vec();
            }
            "CELLS" => {
                let n1 = c.number("the number of cells")?;
                let n2 = c.number("the size of the cell list")?;
                raw_cells.clear();
                if version >= 5.0 && c.word_is("OFFSETS") {
                    let ty = Ty::legacy(c.need("the offsets type")?)?;
                    let off = c.values(n1, ty, binary, "OFFSETS")?;
                    c.need("CONNECTIVITY")?;
                    let ty = Ty::legacy(c.need("the connectivity type")?)?;
                    let con = c.values(n2, ty, binary, "CONNECTIVITY")?;
                    for w in off.windows(2) {
                        let (a, b) = (w[0] as usize, w[1] as usize);
                        let s = con
                            .get(a..b)
                            .ok_or("the offsets go beyond the connectivity")?;
                        raw_cells.push(
                            s.iter()
                                .map(|&v| to_u32(v, "CELLS"))
                                .collect::<Result<_, _>>()?,
                        );
                    }
                } else {
                    let v = c.values(n2, Ty::I32, binary, "CELLS")?;
                    let mut i = 0;
                    for _ in 0..n1 {
                        let k = *v.get(i).ok_or("CELLS ends too early")? as usize;
                        let s = v.get(i + 1..i + 1 + k).ok_or("CELLS ends too early")?;
                        raw_cells.push(
                            s.iter()
                                .map(|&v| to_u32(v, "CELLS"))
                                .collect::<Result<_, _>>()?,
                        );
                        i += 1 + k;
                    }
                }
                n_cells = raw_cells.len();
            }
            "CELL_TYPES" => {
                let n = c.number("the number of cell types")?;
                cell_types = c.values(n, Ty::I32, binary, "CELL_TYPES")?;
            }
            "CELL_DATA" => {
                count_data = c.number("the number of cell data")?;
                cell_data = true;
            }
            "POINT_DATA" => {
                count_data = c.number("the number of point data")?;
                cell_data = false;
            }
            "SCALARS" => {
                let name = c.need("the array name")?.to_string();
                let ty = Ty::legacy(c.need("the array type")?)?;
                // an optional number of components, then LOOKUP_TABLE
                let mut comps = 1usize;
                let mut next = c.need("LOOKUP_TABLE")?;
                if let Ok(k) = next.parse::<usize>() {
                    comps = k;
                    next = c.need("LOOKUP_TABLE")?;
                }
                if !next.eq_ignore_ascii_case("LOOKUP_TABLE") {
                    return Err(format!(
                        "SCALARS {name}: expected LOOKUP_TABLE, found {next}"
                    ));
                }
                c.need("the table name")?;
                let v = c.values(count_data * comps, ty, binary, &name)?;
                if cell_data && comps == 1 {
                    g.arrays.push((name, v));
                }
            }
            "VECTORS" | "NORMALS" => {
                c.need("the array name")?;
                let ty = Ty::legacy(c.need("the array type")?)?;
                c.values(3 * count_data, ty, binary, "VECTORS")?;
            }
            "FIELD" => {
                c.need("the field name")?;
                let k = c.number("the number of arrays")?;
                for _ in 0..k {
                    let name = c.need("the array name")?.to_string();
                    let comps = c.number("the number of components")?;
                    let tuples = c.number("the number of tuples")?;
                    let ty = Ty::legacy(c.need("the array type")?)?;
                    let v = c.values(comps * tuples, ty, binary, &name)?;
                    if cell_data && comps == 1 {
                        g.arrays.push((name, v));
                    }
                }
            }
            other => {
                return Err(format!(
                    "the keyword {other} is not supported (read: POINTS, CELLS, CELL_TYPES and \
                     scalar cell data)"
                ));
            }
        }
    }
    if g.points.is_empty() {
        return Err("the file has no POINTS".into());
    }
    if cell_types.len() != n_cells {
        return Err(format!(
            "{n_cells} cells but {} cell types",
            cell_types.len()
        ));
    }
    for (i, (nodes, t)) in raw_cells.into_iter().zip(cell_types).enumerate() {
        let t = t as u8;
        if t == 42 {
            g.poly.insert(i, face_stream(&nodes)?);
        }
        g.cells.push((t, nodes));
    }
    let how = format!("VTK legacy {}", if binary { "binary" } else { "ASCII" });
    Ok((g, how))
}

impl Cursor<'_> {
    fn word_is(&mut self, w: &str) -> bool {
        let at = self.p;
        let found = self.word().is_some_and(|x| x.eq_ignore_ascii_case(w));
        if !found {
            self.p = at;
        }
        found
    }
}

/// The faces of a polyhedron from its face stream `n_faces, n0, ids, n1, ids, ...`.
fn face_stream(s: &[u32]) -> Result<Vec<Vec<u32>>, String> {
    let n = *s.first().ok_or("a polyhedron without faces")? as usize;
    let mut faces = Vec::with_capacity(n);
    let mut i = 1;
    for _ in 0..n {
        let k = *s.get(i).ok_or("a polyhedron face stream ends too early")? as usize;
        faces.push(
            s.get(i + 1..i + 1 + k)
                .ok_or("a polyhedron face stream ends too early")?
                .to_vec(),
        );
        i += 1 + k;
    }
    Ok(faces)
}

// ---------------------------------------------------------------------
// XML files

#[derive(Default)]
struct El {
    name: String,
    attrs: HashMap<String, String>,
    kids: Vec<El>,
    text: String,
}

impl El {
    fn kid(&self, name: &str) -> Option<&El> {
        self.kids.iter().find(|k| k.name == name)
    }

    fn attr(&self, name: &str) -> Option<&str> {
        self.attrs.get(name).map(|s| s.as_str())
    }
}

fn parse_attrs(s: &str) -> Result<HashMap<String, String>, String> {
    let mut m = HashMap::new();
    let mut rest = s.trim();
    while !rest.is_empty() {
        let eq = rest
            .find('=')
            .ok_or_else(|| format!("a bad attribute near '{rest}'"))?;
        let key = rest[..eq].trim().to_string();
        let after = rest[eq + 1..].trim_start();
        let q = after
            .chars()
            .next()
            .filter(|&c| c == '"' || c == '\'')
            .ok_or_else(|| format!("the attribute {key} has no quotes"))?;
        let end = after[1..]
            .find(q)
            .ok_or_else(|| format!("the attribute {key} is not closed"))?;
        m.insert(key, after[1..1 + end].to_string());
        rest = after[end + 2..].trim_start();
    }
    Ok(m)
}

fn parse_tree(s: &str) -> Result<El, String> {
    let mut stack: Vec<El> = vec![El::default()];
    let mut i = 0;
    while let Some(off) = s[i..].find('<') {
        let at = i + off;
        if let Some(top) = stack.last_mut() {
            top.text.push_str(&s[i..at]);
        }
        let rest = &s[at..];
        if rest.starts_with("<!--") {
            let e = rest.find("-->").ok_or("a comment is not closed")?;
            i = at + e + 3;
            continue;
        }
        if rest.starts_with("<?") {
            let e = rest.find("?>").ok_or("a declaration is not closed")?;
            i = at + e + 2;
            continue;
        }
        // the end of the tag, outside quotes
        let mut q = None;
        let mut end = None;
        for (k, ch) in rest.char_indices() {
            match (q, ch) {
                (None, '"' | '\'') => q = Some(ch),
                (Some(c), x) if c == x => q = None,
                (None, '>') => {
                    end = Some(k);
                    break;
                }
                _ => {}
            }
        }
        let end = end.ok_or("a tag is not closed")?;
        let inner = &rest[1..end];
        i = at + end + 1;
        if let Some(name) = inner.strip_prefix('/') {
            let done = stack.pop().ok_or("an end tag without a start")?;
            if done.name != name.trim() || stack.is_empty() {
                return Err(format!("the end tag </{}> does not match", name.trim()));
            }
            stack.last_mut().unwrap().kids.push(done);
        } else {
            let selfclose = inner.ends_with('/');
            let inner = inner.trim_end_matches('/');
            let (name, attrs) = match inner.find(|c: char| c.is_ascii_whitespace()) {
                Some(k) => (&inner[..k], &inner[k..]),
                None => (inner, ""),
            };
            let el = El {
                name: name.to_string(),
                attrs: parse_attrs(attrs)?,
                ..Default::default()
            };
            if selfclose {
                stack.last_mut().unwrap().kids.push(el);
            } else {
                stack.push(el);
            }
        }
    }
    if stack.len() != 1 {
        return Err("an element is not closed".into());
    }
    stack
        .pop()
        .unwrap()
        .kids
        .into_iter()
        .next()
        .ok_or_else(|| "the file has no elements".into())
}

fn b64(text: &[u8]) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(text.len() / 4 * 3);
    let (mut acc, mut n) = (0u32, 0);
    for &c in text {
        let v = match c {
            b'A'..=b'Z' => c - b'A',
            b'a'..=b'z' => c - b'a' + 26,
            b'0'..=b'9' => c - b'0' + 52,
            b'+' => 62,
            b'/' => 63,
            b'=' => break,
            c if c.is_ascii_whitespace() => continue,
            _ => return Err("bad base64 data".into()),
        };
        acc = (acc << 6) | v as u32;
        n += 1;
        if n == 4 {
            out.extend_from_slice(&[(acc >> 16) as u8, (acc >> 8) as u8, acc as u8]);
            acc = 0;
            n = 0;
        }
    }
    match n {
        2 => out.push((acc >> 4) as u8),
        3 => out.extend_from_slice(&[(acc >> 10) as u8, (acc >> 2) as u8]),
        _ => {}
    }
    Ok(out)
}

struct Xml<'a> {
    /// Header integer size in bytes (4 or 8).
    hl: usize,
    compressed: bool,
    big: bool,
    appended: &'a [u8],
    appended_b64: bool,
}

impl Xml<'_> {
    fn int(&self, b: &[u8], i: usize) -> Result<usize, String> {
        let k = self.hl;
        let s = b
            .get(i * k..(i + 1) * k)
            .ok_or("the binary data ends too early")?;
        let mut a = [0u8; 8];
        a[..k].copy_from_slice(s);
        if self.big {
            a[..k].reverse();
        }
        Ok(u64::from_le_bytes(a) as usize)
    }

    /// The bytes of one binary array: `src` starts at its header.
    fn payload(&self, src: &[u8], base64: bool) -> Result<Vec<u8>, String> {
        let hl = self.hl;
        let chars = |bytes: usize| bytes.div_ceil(3) * 4;
        if !self.compressed {
            if !base64 {
                let n = self.int(src, 0)?;
                return Ok(src
                    .get(hl..hl + n)
                    .ok_or("the binary data ends too early")?
                    .to_vec());
            }
            let head = src
                .get(..chars(hl).min(src.len()))
                .ok_or("the binary data ends too early")?;
            // header and data are encoded separately (padding in the header)
            // or together
            if head.contains(&b'=') {
                let h = b64(head)?;
                let n = self.int(&h, 0)?;
                let start = chars(hl);
                let d = src
                    .get(start..(start + chars(n)).min(src.len()))
                    .ok_or("the binary data ends too early")?;
                let mut d = b64(d)?;
                d.truncate(n);
                Ok(d)
            } else {
                let all = b64(src)?;
                let n = self.int(&all, 0)?;
                Ok(all
                    .get(hl..hl + n)
                    .ok_or("the binary data ends too early")?
                    .to_vec())
            }
        } else {
            // blocks: number, block size, size of the last block, the
            // compressed sizes, then the compressed blocks
            let (head, data): (Vec<u8>, Vec<u8>) = if base64 {
                let h3 = b64(src
                    .get(..chars(3 * hl))
                    .ok_or("the binary data ends too early")?)?;
                let nb = self.int(&h3, 0)?;
                let hbytes = (3 + nb) * hl;
                let hend = chars(hbytes);
                let h = b64(src.get(..hend).ok_or("the binary data ends too early")?)?;
                (h, b64(&src[hend..])?)
            } else {
                let nb = self.int(src, 0)?;
                let hbytes = (3 + nb) * hl;
                (
                    src.get(..hbytes)
                        .ok_or("the binary data ends too early")?
                        .to_vec(),
                    src[hbytes..].to_vec(),
                )
            };
            let nb = self.int(&head, 0)?;
            let (bs, last) = (self.int(&head, 1)?, self.int(&head, 2)?);
            let mut out = Vec::new();
            let mut at = 0usize;
            for k in 0..nb {
                let sz = self.int(&head, 3 + k)?;
                let blk = data
                    .get(at..at + sz)
                    .ok_or("the compressed data ends too early")?;
                at += sz;
                let mut d = inflate::zlib(blk)?;
                let want = if k + 1 == nb && last != 0 { last } else { bs };
                if d.len() != want {
                    return Err(format!(
                        "a compressed block has {} bytes, expected {want}",
                        d.len()
                    ));
                }
                out.append(&mut d);
            }
            Ok(out)
        }
    }

    fn array(&self, el: &El) -> Result<Vec<f64>, String> {
        let name = el.attr("Name").unwrap_or("?");
        let ty = Ty::xml(el.attr("type").unwrap_or("Float32"))?;
        let big = self.big;
        match el.attr("format").unwrap_or("ascii") {
            "ascii" => el
                .text
                .split_ascii_whitespace()
                .map(|w| {
                    w.parse::<f64>()
                        .map_err(|_| format!("array {name}: '{w}' is not a number"))
                })
                .collect(),
            "binary" => {
                let p = self.payload(el.text.trim().as_bytes(), true)?;
                Ok(ty.values(&p, big))
            }
            "appended" => {
                let off: usize = el
                    .attr("offset")
                    .and_then(|o| o.parse().ok())
                    .ok_or_else(|| format!("array {name}: no offset"))?;
                let src = self
                    .appended
                    .get(off..)
                    .ok_or_else(|| format!("array {name}: the offset is beyond the data"))?;
                let p = self.payload(src, self.appended_b64)?;
                Ok(ty.values(&p, big))
            }
            other => Err(format!("array {name}: the format '{other}' is not read")),
        }
    }
}

fn parse_xml(bytes: &[u8]) -> Result<(Grid, String), String> {
    // the appended data are not XML: cut them out
    let (xml_bytes, appended, app_b64) = match find(bytes, b"<AppendedData") {
        Some(a) => {
            let gt = a + find(&bytes[a..], b">").ok_or("AppendedData is not closed")? + 1;
            let tag = std::str::from_utf8(&bytes[a..gt]).unwrap_or("");
            let end = rfind(bytes, b"</AppendedData>").ok_or("AppendedData is not closed")?;
            let mut data = &bytes[gt..end];
            while data.first().is_some_and(|c| c.is_ascii_whitespace()) {
                data = &data[1..];
            }
            if data.first() == Some(&b'_') {
                data = &data[1..];
            }
            let mut x = bytes[..gt].to_vec();
            x.extend_from_slice(&bytes[end..]);
            (x, data, tag.contains("base64"))
        }
        None => (bytes.to_vec(), &bytes[..0], false),
    };
    let text = std::str::from_utf8(&xml_bytes).map_err(|_| "the XML is not UTF-8 text")?;
    let root = parse_tree(text)?;
    if root.name != "VTKFile" {
        return Err("not a VTKFile".into());
    }
    let ty = root.attr("type").unwrap_or("");
    if ty != "UnstructuredGrid" {
        return Err(format!(
            "the file is a {ty}: only UnstructuredGrid is read (a parallel .pvtu: read its pieces)"
        ));
    }
    let compressor = root.attr("compressor");
    if let Some(c) = compressor
        && c != "vtkZLibDataCompressor"
    {
        return Err(format!(
            "the compressor {c} is not read (zlib is; save the file without compression or with zlib)"
        ));
    }
    let x = Xml {
        hl: if root.attr("header_type") == Some("UInt64") {
            8
        } else {
            4
        },
        compressed: compressor.is_some(),
        big: root.attr("byte_order") == Some("BigEndian"),
        appended,
        appended_b64: app_b64,
    };
    let grid = root.kid("UnstructuredGrid").ok_or("no UnstructuredGrid")?;
    let pieces: Vec<&El> = grid.kids.iter().filter(|k| k.name == "Piece").collect();
    if pieces.len() != 1 {
        return Err(format!(
            "{} pieces: only a single piece is read",
            pieces.len()
        ));
    }
    let piece = pieces[0];
    let mut g = Grid::default();
    let data_of = |parent: &El, name: &str| -> Result<Option<Vec<f64>>, String> {
        parent
            .kids
            .iter()
            .find(|k| k.name == "DataArray" && k.attr("Name") == Some(name))
            .map(|a| x.array(a))
            .transpose()
    };
    let pts = piece
        .kid("Points")
        .and_then(|p| p.kids.iter().find(|k| k.name == "DataArray"))
        .ok_or("the piece has no Points")?;
    let v = x.array(pts)?;
    g.points = v.as_chunks::<3>().0.to_vec();

    let cells = piece.kid("Cells").ok_or("the piece has no Cells")?;
    let need = |n: &str| -> Result<Vec<f64>, String> {
        data_of(cells, n)?.ok_or_else(|| format!("Cells has no {n}"))
    };
    let (conn, offs, types) = (need("connectivity")?, need("offsets")?, need("types")?);
    let faces = data_of(cells, "faces")?;
    let faceoffs = data_of(cells, "faceoffsets")?;
    if offs.len() != types.len() {
        return Err(format!(
            "{} offsets for {} cell types",
            offs.len(),
            types.len()
        ));
    }
    let mut start = 0usize;
    let mut fstart = 0usize;
    for (i, (&end, &t)) in offs.iter().zip(&types).enumerate() {
        let end = end as usize;
        let nodes = conn
            .get(start..end)
            .ok_or("the offsets go beyond the connectivity")?
            .iter()
            .map(|&v| to_u32(v, "connectivity"))
            .collect::<Result<Vec<_>, _>>()?;
        start = end;
        let t = t as u8;
        if t == 42 {
            let (f, fo) = faces
                .as_ref()
                .zip(faceoffs.as_ref())
                .ok_or("polyhedra need the faces and faceoffsets arrays")?;
            let fend = *fo.get(i).ok_or("faceoffsets ends too early")? as usize;
            let s = f
                .get(fstart..fend)
                .ok_or("faceoffsets go beyond the faces")?
                .iter()
                .map(|&v| to_u32(v, "faces"))
                .collect::<Result<Vec<_>, _>>()?;
            fstart = fend;
            g.poly.insert(i, face_stream(&s)?);
        }
        g.cells.push((t, nodes));
    }
    if let Some(cd) = piece.kid("CellData") {
        for a in cd.kids.iter().filter(|k| k.name == "DataArray") {
            let comps = a.attr("NumberOfComponents").unwrap_or("1");
            if comps == "1" {
                g.arrays
                    .push((a.attr("Name").unwrap_or("").to_string(), x.array(a)?));
            }
        }
    }
    // the field data may sit in the grid or in the piece
    for fd in [grid.kid("FieldData"), piece.kid("FieldData")]
        .into_iter()
        .flatten()
    {
        for a in fd.kids.iter().filter(|k| k.name == "DataArray") {
            if a.attr("type") != Some("String") {
                continue;
            }
            let codes = x.array_ascii_ints(a)?;
            let mut list = Vec::new();
            let mut cur = String::new();
            for c in codes {
                if c == 0 {
                    list.push(std::mem::take(&mut cur));
                } else {
                    cur.push(c as u8 as char);
                }
            }
            g.names
                .insert(a.attr("Name").unwrap_or("").to_string(), list);
        }
    }
    let enc = match (
        compressor.is_some(),
        piece
            .kid("Points")
            .and_then(|p| p.kids.first())
            .and_then(|a| a.attr("format")),
    ) {
        (true, Some(f)) => format!("{f}, zlib"),
        (false, Some(f)) => f.to_string(),
        _ => String::new(),
    };
    Ok((g, format!("VTK XML ({enc})")))
}

impl Xml<'_> {
    /// A string array: character codes, written as numbers of any encoding.
    fn array_ascii_ints(&self, el: &El) -> Result<Vec<u32>, String> {
        // String arrays hold UInt8 codes; read them like a UInt8 array
        let mut e = El {
            name: el.name.clone(),
            attrs: el.attrs.clone(),
            kids: Vec::new(),
            text: el.text.clone(),
        };
        e.attrs.insert("type".into(), "UInt8".into());
        Ok(self.array(&e)?.into_iter().map(|v| v as u32).collect())
    }
}

fn find(h: &[u8], n: &[u8]) -> Option<usize> {
    h.windows(n.len()).position(|w| w == n)
}

fn rfind(h: &[u8], n: &[u8]) -> Option<usize> {
    h.windows(n.len()).rposition(|w| w == n)
}

// ---------------------------------------------------------------------
// writing

/// Writes `mesh` as a VTK XML unstructured grid (`.vtu`): every cell a
/// polyhedron, every boundary face a polygon, the cell-data array `region`
/// (cell zone of a cell, boundary of a face) and the names of the zones and
/// boundaries as field data. `binary` selects appended binary data (little
/// endian, uncompressed) instead of ASCII.
pub fn write_vtu(mesh: &Mesh, path: &Path, binary: bool) -> Result<(), String> {
    let n_cells = mesh.n_cells;
    // faces of each cell, outward
    let mut cell_faces: Vec<Vec<Vec<u32>>> = vec![Vec::new(); n_cells];
    let nf = mesh.n_faces();
    let ni = mesh.n_interior_faces();
    let nodes =
        |f: usize| &mesh.face_nodes[mesh.face_ptr[f] as usize..mesh.face_ptr[f + 1] as usize];
    for f in 0..nf {
        let n = nodes(f).to_vec();
        if f < ni {
            let mut r = n.clone();
            r.reverse();
            cell_faces[mesh.neighbour[f] as usize].push(r);
        }
        cell_faces[mesh.owner[f] as usize].push(n);
    }
    let mut conn: Vec<i64> = Vec::new();
    let mut offs: Vec<i64> = Vec::new();
    let mut types: Vec<u8> = Vec::new();
    let mut faces: Vec<i64> = Vec::new();
    let mut faceoffs: Vec<i64> = Vec::new();
    let mut region: Vec<i32> = Vec::new();
    for (c, fs) in cell_faces.iter().enumerate() {
        let mut uniq: Vec<u32> = fs.iter().flatten().copied().collect();
        uniq.sort_unstable();
        uniq.dedup();
        conn.extend(uniq.iter().map(|&v| v as i64));
        offs.push(conn.len() as i64);
        types.push(42);
        faces.push(fs.len() as i64);
        for f in fs {
            faces.push(f.len() as i64);
            faces.extend(f.iter().map(|&v| v as i64));
        }
        faceoffs.push(faces.len() as i64);
        region.push(mesh.cell_zone[c] as i32);
    }
    for (bi, b) in mesh.boundaries.iter().enumerate() {
        for f in b.range() {
            let n = nodes(f);
            conn.extend(n.iter().map(|&v| v as i64));
            offs.push(conn.len() as i64);
            types.push(match n.len() {
                3 => 5,
                4 => 9,
                _ => 7,
            });
            faceoffs.push(-1);
            region.push(bi as i32);
        }
    }
    let n_out = types.len();
    let mut pts: Vec<f64> = Vec::with_capacity(3 * mesh.points.len());
    for p in &mesh.points {
        pts.extend_from_slice(p);
    }
    let names = |list: Vec<String>| -> Vec<u8> {
        let mut v = Vec::new();
        for (i, n) in list.iter().enumerate() {
            v.extend(format!("{i}={n}").bytes());
            v.push(0);
        }
        v
    };
    let zone_names = names(mesh.zones.clone());
    let bound_names = names(mesh.boundaries.iter().map(|b| b.name.clone()).collect());

    enum Data<'a> {
        F64(&'a [f64]),
        I64(&'a [i64]),
        I32(&'a [i32]),
        U8(&'a [u8]),
    }
    impl Data<'_> {
        fn ty(&self) -> &'static str {
            match self {
                Data::F64(_) => "Float64",
                Data::I64(_) => "Int64",
                Data::I32(_) => "Int32",
                Data::U8(_) => "UInt8",
            }
        }
        fn len(&self) -> usize {
            match self {
                Data::F64(v) => v.len(),
                Data::I64(v) => v.len(),
                Data::I32(v) => v.len(),
                Data::U8(v) => v.len(),
            }
        }
        fn bytes(&self) -> Vec<u8> {
            match self {
                Data::F64(v) => v.iter().flat_map(|x| x.to_le_bytes()).collect(),
                Data::I64(v) => v.iter().flat_map(|x| x.to_le_bytes()).collect(),
                Data::I32(v) => v.iter().flat_map(|x| x.to_le_bytes()).collect(),
                Data::U8(v) => v.to_vec(),
            }
        }
        fn ascii(&self) -> String {
            let mut s = String::new();
            let mut put = |x: String| {
                s.push_str(&x);
                s.push(' ');
            };
            match self {
                Data::F64(v) => v.iter().for_each(|x| put(format!("{x:e}"))),
                Data::I64(v) => v.iter().for_each(|x| put(x.to_string())),
                Data::I32(v) => v.iter().for_each(|x| put(x.to_string())),
                Data::U8(v) => v.iter().for_each(|x| put(x.to_string())),
            }
            s
        }
    }
    let mut appended: Vec<u8> = Vec::new();
    let mut array = |name: &str, comps: usize, d: Data| -> String {
        if binary {
            let off = appended.len();
            let b = d.bytes();
            appended.extend((b.len() as u64).to_le_bytes());
            appended.extend(b);
            format!(
                "<DataArray type=\"{}\" Name=\"{name}\" NumberOfComponents=\"{comps}\" \
                 format=\"appended\" offset=\"{off}\"/>\n",
                d.ty()
            )
        } else {
            format!(
                "<DataArray type=\"{}\" Name=\"{name}\" NumberOfComponents=\"{comps}\" \
                 format=\"ascii\">\n{}\n</DataArray>\n",
                d.ty(),
                d.ascii()
            )
        }
    };
    let string_array = |name: &str, bytes: &[u8]| -> String {
        let codes = Data::U8(bytes).ascii();
        format!(
            "<DataArray type=\"String\" Name=\"{name}\" NumberOfTuples=\"1\" format=\"ascii\">\
             {codes}</DataArray>\n"
        )
    };
    let mut x = String::new();
    x.push_str(
        "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" \
         byte_order=\"LittleEndian\" header_type=\"UInt64\">\n<UnstructuredGrid>\n",
    );
    x.push_str("<FieldData>\n");
    x.push_str(&string_array("zone_names", &zone_names));
    x.push_str(&string_array("boundary_names", &bound_names));
    x.push_str("</FieldData>\n");
    x.push_str(&format!(
        "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{n_out}\">\n<Points>\n",
        mesh.points.len()
    ));
    x.push_str(&array("Points", 3, Data::F64(&pts)));
    x.push_str("</Points>\n<Cells>\n");
    x.push_str(&array("connectivity", 1, Data::I64(&conn)));
    x.push_str(&array("offsets", 1, Data::I64(&offs)));
    x.push_str(&array("types", 1, Data::U8(&types)));
    x.push_str(&array("faces", 1, Data::I64(&faces)));
    x.push_str(&array("faceoffsets", 1, Data::I64(&faceoffs)));
    x.push_str("</Cells>\n<CellData>\n");
    x.push_str(&array("region", 1, Data::I32(&region)));
    x.push_str("</CellData>\n</Piece>\n</UnstructuredGrid>\n");
    let _ = Data::len;
    let mut out = x.into_bytes();
    if binary {
        out.extend_from_slice(b"<AppendedData encoding=\"raw\">\n_");
        out.extend(appended);
        out.extend_from_slice(b"\n</AppendedData>\n");
    }
    out.extend_from_slice(b"</VTKFile>\n");
    std::fs::write(path, out).map_err(|e| format!("{}: {e}", path.display()))
}

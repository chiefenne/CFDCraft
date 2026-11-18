//! Reduction multigrid for the non-symmetric matrices of transport
//! equations with strong convection (D179, L44): the sweep of the momentum
//! equations in the preconditioner of the coupled algorithm (D176).
//!
//! Aggregation multigrid and point smoothers rely on a symmetric positive
//! definite matrix. An upwind-dominated matrix is close to triangular
//! instead: information travels along the flow. A reduction method uses
//! that structure (Ries, Trottenberg, Winter 1983; approximate ideal
//! restriction: Manteuffel, Ruge, Southworth 2018; Manteuffel, Münzenmaier,
//! Ruge, Southworth 2019). The unknowns are split into coarse points C and
//! fine points F such that the block A_FF is strongly diagonally dominant;
//! then one Jacobi sweep on the F points and an exact solve of the Schur
//! complement on the C points would solve the system. The method
//! approximates that Schur complement by a Galerkin operator R A P with
//!
//! - the restriction R = [Z I], Z ≈ −A_CF A_FF⁻¹ (the *ideal* restriction,
//!   which makes the coarse residual equation independent of the F error)
//!   with the two-term Neumann series A_FF⁻¹ ≈ 2 D_F⁻¹ − D_F⁻¹ A_FF D_F⁻¹
//!   on the strong C–F connections;
//! - the interpolation P = [W; I] with the direct interpolation of
//!   classical multigrid (Stüben 2001), which keeps the row sums and so
//!   also serves where diffusion dominates.
//!
//! **Setup** (host, float64, deterministic), per level:
//!
//! 1. strong dependencies, the classical measure: j is strong for i when
//!    |a_ij| ≥ θ max_{k≠i} |a_ik|;
//! 2. the first pass of the Ruge–Stüben splitting: the undecided point on
//!    which most points depend becomes C (ties: the smallest index), the
//!    points that depend on it become F; points on which nobody depends are
//!    F from the start;
//! 3. R and P as above; the coarse operator A_c = R A P, of which the
//!    entries below a fraction of the diagonal of their row are dropped
//!    (not lumped);
//! 4. down to at most `max_coarse` unknowns, solved by a dense inverse
//!    (Gauss–Jordan with partial pivoting).
//!
//! **Cycle**: the V-cycle of [`super::cycle`] on the level operations of
//! the backends, unchanged: the hierarchy stores an inverse diagonal that is
//! zero on the C points, so that the Jacobi step of the cycle relaxes the F
//! points only, and one step before and after the coarse correction. Two F
//! sweeps per level in total fit the two terms of the Neumann series: an
//! odd number of sweeps leaves error modes with a factor near −1.

use super::Precision;
use super::csr::{Csr, RowAccumulator};
use super::setup::{Hierarchy, HierarchyInfo, Level};
use anyhow::{Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::time::Instant;

/// Settings of the reduction multigrid (`[solver.coupled.multigrid]`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct ReductionSettings {
    /// Strength threshold θ of the splitting and the interpolation: j is
    /// strong for i when |a_ij| ≥ θ max |a_ik| over the neighbours k.
    pub strength: f64,
    /// The same threshold for the connections that enter the restriction.
    pub restriction_strength: f64,
    /// Entries of a coarse operator below this fraction of the diagonal of
    /// their row are dropped.
    pub drop_tolerance: f64,
    /// Largest coarsest level (solved by a dense inverse).
    pub max_coarse: usize,
    /// Largest number of levels.
    pub max_levels: usize,
    /// Precision of the hierarchy.
    pub precision: Precision,
}

impl Default for ReductionSettings {
    fn default() -> Self {
        Self {
            strength: 0.3,
            restriction_strength: 0.05,
            drop_tolerance: 1e-4,
            max_coarse: 64,
            max_levels: 30,
            precision: Precision::Mixed,
        }
    }
}

impl ReductionSettings {
    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.strength > 0.0 && self.strength <= 1.0,
            "[multigrid] strength must lie in (0, 1]"
        );
        ensure!(
            self.restriction_strength > 0.0 && self.restriction_strength <= 1.0,
            "[multigrid] restriction_strength must lie in (0, 1]"
        );
        ensure!(
            (0.0..0.5).contains(&self.drop_tolerance),
            "[multigrid] drop_tolerance must lie in [0, 0.5)"
        );
        ensure!(
            (1..=4000).contains(&self.max_coarse),
            "[multigrid] max_coarse must lie in 1..=4000 (dense inverse)"
        );
        ensure!(self.max_levels >= 1, "[multigrid] max_levels must be >= 1");
        Ok(())
    }
}

/// Largest |a_ik| over the neighbours k ≠ i of row `i`.
fn largest_neighbour(a: &Csr, i: usize) -> f64 {
    a.row(i)
        .filter(|&(j, _)| j != i)
        .fold(0.0f64, |m, (_, v)| m.max(v.abs()))
}

/// A graph in compressed rows: the points `idx[ptr[i]..ptr[i + 1]]` of
/// point i, ascending.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Graph {
    /// Start of each row (length: points + 1).
    pub ptr: Vec<usize>,
    /// The points of all rows.
    pub idx: Vec<u32>,
}

impl Graph {
    /// Number of points.
    pub fn len(&self) -> usize {
        self.ptr.len().saturating_sub(1)
    }

    /// Whether the graph has no points.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The points of row `i`.
    pub fn row(&self, i: usize) -> &[u32] {
        &self.idx[self.ptr[i]..self.ptr[i + 1]]
    }
}

/// Strong dependencies of every row, columns ascending: j ≠ i with
/// |a_ij| ≥ θ max_{k≠i} |a_ik|.
pub fn strong(a: &Csr, theta: f64) -> Graph {
    // the pattern of a CSR matrix with the values unused
    let g = Csr::build_rows(
        a.n_rows,
        a.n_cols,
        || (),
        |_, i, col, _| {
            let m = largest_neighbour(a, i);
            if m == 0.0 {
                return;
            }
            col.extend(
                a.row(i)
                    .filter(|&(j, v)| j != i && v.abs() >= theta * m)
                    .map(|(j, _)| j as u32),
            );
        },
    );
    Graph {
        ptr: g.row_ptr,
        idx: g.col,
    }
}

/// First pass of the Ruge–Stüben splitting on the graph of strong
/// dependencies (`strong[i]`: the points i depends on). Returns `true` for
/// the C points.
///
/// The measure λ_i of an undecided point is the number of undecided points
/// that depend on it plus the number of F points that do. The point with
/// the largest measure (ties: the smallest index) becomes C; the undecided
/// points that depend on it become F, which raises the measure of the
/// other points they depend on; the measure of the points the new C point
/// depends on falls by one. Points on which nobody depends are F from the
/// start. The result depends on the graph only, not on threads.
pub fn split(strong: &Graph) -> Vec<bool> {
    const U: u8 = 0;
    const F: u8 = 1;
    const C: u8 = 2;
    let n = strong.len();
    // the transposed graph: the points that depend on i
    let mut ptr = vec![0usize; n + 1];
    for &j in &strong.idx {
        ptr[j as usize + 1] += 1;
    }
    for i in 0..n {
        ptr[i + 1] += ptr[i];
    }
    let mut next = ptr.clone();
    let mut dependants = vec![0u32; ptr[n]];
    for i in 0..n {
        for &j in strong.row(i) {
            dependants[next[j as usize]] = i as u32;
            next[j as usize] += 1;
        }
    }
    let mut lambda: Vec<u32> = (0..n).map(|i| (ptr[i + 1] - ptr[i]) as u32).collect();
    let mut state: Vec<u8> = lambda.iter().map(|&l| if l == 0 { F } else { U }).collect();
    // lazy heap: an entry is valid while its measure is the current one
    let mut heap: BinaryHeap<(u32, Reverse<u32>)> = (0..n)
        .filter(|&i| state[i] == U)
        .map(|i| (lambda[i], Reverse(i as u32)))
        .collect();
    while let Some((l, Reverse(i))) = heap.pop() {
        let i = i as usize;
        if state[i] != U || l != lambda[i] {
            continue;
        }
        state[i] = C;
        for &j in &dependants[ptr[i]..ptr[i + 1]] {
            let j = j as usize;
            if state[j] != U {
                continue;
            }
            state[j] = F;
            for &k in strong.row(j) {
                let k = k as usize;
                if state[k] == U {
                    lambda[k] += 1;
                    heap.push((lambda[k], Reverse(k as u32)));
                }
            }
        }
        for &j in strong.row(i) {
            let j = j as usize;
            if state[j] == U && lambda[j] > 0 {
                lambda[j] -= 1;
                heap.push((lambda[j], Reverse(j as u32)));
            }
        }
    }
    state.iter().map(|&s| s == C).collect()
}

/// The interpolation P = [W; I] (rows: all points, columns: C points in
/// their order): a C point takes its coarse value; an F point i the direct
/// interpolation from its strong C neighbours C_i,
/// w_ij = −α a_ij⁻ / a_ii − β a_ij⁺ / a_ii with
/// α = Σ_{k≠i} a_ik⁻ / Σ_{k∈C_i} a_ik⁻ and β likewise for the positive
/// entries, so that the row sums of A are kept. Couplings of a sign
/// without a strong C neighbour are added to the diagonal when they have
/// its sign and left to the relaxation otherwise. An F point without
/// strong C neighbours is not interpolated (it is relaxed only).
fn interpolation(
    a: &Csr,
    diag: &[f64],
    strong: &Graph,
    coarse: &[bool],
    cidx: &[u32],
    nc: usize,
) -> Csr {
    // the selected entries of a row come in column order (cidx is
    // monotone), each once
    Csr::build_rows(a.n_rows, nc, Vec::<(u32, f64)>::new, |sel, i, col, val| {
        if coarse[i] {
            col.push(cidx[i]);
            val.push(1.0);
            return;
        }
        let (mut neg_all, mut pos_all, mut neg_c, mut pos_c) = (0.0, 0.0, 0.0, 0.0);
        let mut s = strong.row(i).iter().map(|&j| j as usize).peekable();
        sel.clear();
        for (j, v) in a.row(i) {
            if j == i {
                continue;
            }
            if v < 0.0 {
                neg_all += v;
            } else {
                pos_all += v;
            }
            while s.peek().is_some_and(|&sj| sj < j) {
                s.next();
            }
            if s.peek() == Some(&j) && coarse[j] {
                if v < 0.0 {
                    neg_c += v;
                } else {
                    pos_c += v;
                }
                sel.push((cidx[j], v));
            }
        }
        if sel.is_empty() {
            return;
        }
        let mut d = diag[i];
        let alpha = if neg_c != 0.0 {
            neg_all / neg_c
        } else {
            if d < 0.0 {
                d += neg_all;
            }
            0.0
        };
        let beta = if pos_c != 0.0 {
            pos_all / pos_c
        } else {
            if d > 0.0 {
                d += pos_all;
            }
            0.0
        };
        for &(c, v) in sel.iter() {
            col.push(c);
            val.push(-(if v < 0.0 { alpha } else { beta }) * v / d);
        }
    })
}

/// The restriction R = [Z I] (rows: C points in their order, columns: all
/// points): Z = −A_CF^s (2 D_F⁻¹ − D_F⁻¹ A_FF D_F⁻¹), where A_CF^s keeps the
/// C–F entries with |a_cf| ≥ θ max_{k≠c} |a_ck|.
fn restriction(a: &Csr, diag: &[f64], coarse: &[bool], theta: f64) -> Csr {
    let c_points: Vec<usize> = (0..a.n_rows).filter(|&i| coarse[i]).collect();
    // entries of a column summed in the order they arrive, as
    // `Csr::from_rows` sums them
    Csr::build_rows(
        c_points.len(),
        a.n_rows,
        || RowAccumulator::new(a.n_rows),
        |acc, r, col, val| {
            let c = c_points[r];
            let m = largest_neighbour(a, c);
            acc.add(c as u32, 1.0);
            for (f, s) in a.row(c) {
                if f == c || coarse[f] || s == 0.0 || s.abs() < theta * m {
                    continue;
                }
                let df = diag[f];
                acc.add(f as u32, -2.0 * s / df);
                for (g, v) in a.row(f) {
                    if !coarse[g] {
                        acc.add(g as u32, s * v / (df * diag[g]));
                    }
                }
            }
            acc.drain_sorted(col, val);
        },
    )
}

/// Removes the entries with |a_ij| ≤ `tol` |a_ii| (the diagonal stays).
fn drop_small(a: &Csr, tol: f64) -> Csr {
    if tol <= 0.0 {
        return a.clone();
    }
    let diag = a.diagonal();
    let mut out = Csr {
        n_rows: a.n_rows,
        n_cols: a.n_cols,
        row_ptr: vec![0],
        ..Csr::default()
    };
    for (i, d) in diag.iter().enumerate() {
        let lim = tol * d.abs();
        for (j, v) in a.row(i) {
            if j == i || v.abs() > lim {
                out.col.push(j as u32);
                out.val.push(v);
            }
        }
        out.row_ptr.push(out.col.len());
    }
    out
}

/// Dense inverse of a general matrix (row-major) by Gauss–Jordan
/// elimination with partial pivoting; fails when a pivot vanishes.
pub fn dense_inverse(a: &Csr) -> Result<Vec<f64>> {
    let n = a.n_rows;
    let mut m = vec![0.0; n * n];
    for i in 0..n {
        for (j, v) in a.row(i) {
            m[i * n + j] = v;
        }
    }
    let scale = m.iter().fold(0.0f64, |s, v| s.max(v.abs()));
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        inv[i * n + i] = 1.0;
    }
    for c in 0..n {
        let (mut pr, mut pv) = (c, m[c * n + c].abs());
        for r in c + 1..n {
            let v = m[r * n + c].abs();
            if v > pv {
                (pr, pv) = (r, v);
            }
        }
        // (a NaN fails the comparison too)
        let regular = pv > 1e-13 * scale;
        if !regular {
            bail!(
                "the coarsest multigrid level ({n} unknowns) of the transport matrix is \
                 singular or not finite"
            );
        }
        if pr != c {
            for k in 0..n {
                m.swap(c * n + k, pr * n + k);
                inv.swap(c * n + k, pr * n + k);
            }
        }
        let d = 1.0 / m[c * n + c];
        for k in 0..n {
            m[c * n + k] *= d;
            inv[c * n + k] *= d;
        }
        for r in 0..n {
            let f = m[r * n + c];
            if r == c || f == 0.0 {
                continue;
            }
            for k in 0..n {
                m[r * n + k] -= f * m[c * n + k];
                inv[r * n + k] -= f * inv[c * n + k];
            }
        }
    }
    Ok(inv)
}

/// Builds the reduction hierarchy of the matrix `a` (square, with a
/// non-zero diagonal; it need not be symmetric).
pub fn build(a: Csr, s: &ReductionSettings) -> Result<Hierarchy> {
    s.validate()?;
    let t0 = Instant::now();
    ensure!(
        a.n_rows == a.n_cols && a.n_rows > 0,
        "multigrid: square matrix needed"
    );
    let mut levels = Vec::new();
    let mut cur = a;
    while cur.n_rows > s.max_coarse && levels.len() + 1 < s.max_levels {
        let n = cur.n_rows;
        let diag = cur.diagonal();
        ensure!(
            diag.iter().all(|&d| d != 0.0 && d.is_finite()),
            "multigrid of a transport matrix: a diagonal entry is zero or not finite on \
             level {}",
            levels.len()
        );
        let st = strong(&cur, s.strength);
        let mut coarse = split(&st);
        if !coarse.iter().any(|&c| c) {
            // a matrix without couplings: one sweep solves it; the cycle
            // needs a coarse level all the same
            coarse[0] = true;
        }
        let mut cidx = vec![u32::MAX; n];
        let mut nc = 0usize;
        for i in 0..n {
            if coarse[i] {
                cidx[i] = nc as u32;
                nc += 1;
            }
        }
        if nc as f64 > 0.9 * n as f64 {
            break;
        }
        let p = interpolation(&cur, &diag, &st, &coarse, &cidx, nc);
        let r = restriction(&cur, &diag, &coarse, s.restriction_strength);
        let ac = drop_small(&r.matmul(&cur.matmul(&p)), s.drop_tolerance);
        levels.push(Level {
            // relaxation on the F points only
            dinv: (0..n)
                .map(|i| if coarse[i] { 0.0 } else { 1.0 / diag[i] })
                .collect(),
            steps: vec![(0.0, 1.0)],
            lambda_max: 1.0,
            a: cur,
            p,
            r,
        });
        cur = ac;
    }
    ensure!(
        cur.n_rows <= 4000,
        "multigrid of a transport matrix: coarsening stalled at {} unknowns (more than the \
         dense limit of 4000)",
        cur.n_rows
    );
    let coarse_inv = dense_inverse(&cur)?;
    let tb = match s.precision {
        Precision::Mixed => 4,
        Precision::Double => 8,
    };
    let mut info = HierarchyInfo {
        precision: s.precision.name().into(),
        ..HierarchyInfo::default()
    };
    for l in &levels {
        info.rows.push(l.a.n_rows);
        info.nnz.push(l.a.nnz());
        info.bytes += (l.a.nnz() + l.p.nnz() + l.r.nnz()) * (tb + 4) + l.a.n_rows * tb;
    }
    info.rows.push(cur.n_rows);
    info.nnz.push(cur.nnz());
    info.bytes += cur.n_rows * cur.n_rows * tb;
    info.grid_complexity = info.rows.iter().sum::<usize>() as f64 / info.rows[0] as f64;
    info.operator_complexity = info.nnz.iter().sum::<usize>() as f64 / info.nnz[0] as f64;
    info.setup_seconds = t0.elapsed().as_secs_f64();
    Ok(Hierarchy {
        levels,
        coarse_n: cur.n_rows,
        coarse_inv,
        precision: s.precision,
        info,
    })
}

/// One V-cycle z = M⁻¹ r of a hierarchy on the host in float64, the cycle
/// of [`super::cycle::vcycle`] (for tests and studies of the setup; the
/// solvers run the cycle on a backend).
pub fn host_cycle(h: &Hierarchy, r: &[f64]) -> Vec<f64> {
    fn level(h: &Hierarchy, l: usize, b: &[f64]) -> Vec<f64> {
        if l == h.levels.len() {
            let n = h.coarse_n;
            return (0..n)
                .map(|i| (0..n).map(|j| h.coarse_inv[i * n + j] * b[j]).sum())
                .collect();
        }
        let lv = &h.levels[l];
        let n = lv.a.n_rows;
        let mut x = vec![0.0; n];
        let mut d = vec![0.0; n];
        let step = |x: &mut Vec<f64>, d: &mut Vec<f64>, c1: f64, c2: f64| {
            let ax = lv.a.spmv(x);
            for i in 0..n {
                d[i] = c1 * d[i] + c2 * lv.dinv[i] * (b[i] - ax[i]);
                x[i] += d[i];
            }
        };
        for &(c1, c2) in &lv.steps {
            step(&mut x, &mut d, c1, c2);
        }
        let ax = lv.a.spmv(&x);
        let w: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
        let xc = level(h, l + 1, &lv.r.spmv(&w));
        let px = lv.p.spmv(&xc);
        for i in 0..n {
            x[i] += px[i];
        }
        for &(c1, c2) in &lv.steps {
            step(&mut x, &mut d, c1, c2);
        }
        x
    }
    level(h, 0, r)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Upwind convection–diffusion on an m × m grid of unit cells with the
    /// velocity field of a solid-body rotation around the centre, scaled so
    /// that the cell Péclet number is `peclet` at the rim, plus `shift` on
    /// the diagonal (a pseudo-time term). Dirichlet boundaries.
    fn rotation(m: usize, peclet: f64, shift: f64) -> Csr {
        let id = |i: usize, j: usize| (i * m + j) as u32;
        let c = 0.5 * (m as f64 - 1.0);
        let rows = (0..m * m)
            .map(|k| {
                let (i, j) = (k / m, k % m);
                // velocity (u along j, v along i) at the faces of the cell
                let u = |ii: f64| -peclet * (ii - c) / c;
                let v = |jj: f64| peclet * (jj - c) / c;
                let mut diag = shift;
                let mut r = Vec::new();
                // (neighbour exists, neighbour id, outflow through the face)
                let faces = [
                    (j + 1 < m, id(i, (j + 1).min(m - 1)), u(i as f64)),
                    (j > 0, id(i, j.saturating_sub(1)), -u(i as f64)),
                    (i + 1 < m, id((i + 1).min(m - 1), j), v(j as f64)),
                    (i > 0, id(i.saturating_sub(1), j), -v(j as f64)),
                ];
                for (inside, nb, out) in faces {
                    // diffusion 1 per face, upwind convection
                    diag += 1.0 + out.max(0.0);
                    if inside {
                        r.push((nb, -1.0 - (-out).max(0.0)));
                    }
                }
                r.push((k as u32, diag));
                r
            })
            .collect();
        Csr::from_rows(m * m, rows)
    }

    /// ‖b − A x‖ after `cycles` stationary cycles, relative to ‖b‖.
    fn reduction(a: &Csr, h: &Hierarchy, cycles: usize) -> f64 {
        let n = a.n_rows;
        let b: Vec<f64> = (0..n).map(|i| 1.0 + ((i * 7919) % 13) as f64).collect();
        let norm = |v: &[f64]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
        let mut x = vec![0.0; n];
        for _ in 0..cycles {
            let ax = a.spmv(&x);
            let r: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
            let z = host_cycle(h, &r);
            for i in 0..n {
                x[i] += z[i];
            }
        }
        let ax = a.spmv(&x);
        let r: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
        norm(&r) / norm(&b)
    }

    #[test]
    fn splitting_of_a_laplacian() {
        let a = super::super::setup::tests::lap2d(24);
        let c = split(&strong(&a, 0.3));
        let nc = c.iter().filter(|&&c| c).count();
        // red-black-like: about half of the points, never two neighbours
        // that both stay undecided
        assert!(nc > 24 * 24 / 4 && nc < 24 * 24 * 2 / 3, "{nc} C points");
    }

    #[test]
    fn splitting_follows_the_flow() {
        // pure upwind transport along a line: every other point is C
        let n = 40;
        let rows = (0..n)
            .map(|i| {
                let mut r = vec![(i as u32, 1.0)];
                if i > 0 {
                    r.push((i as u32 - 1, -1.0));
                }
                r
            })
            .collect();
        let a = Csr::from_rows(n, rows);
        let c = split(&strong(&a, 0.3));
        for i in 0..n - 1 {
            assert!(c[i] != c[i + 1], "points {i} and {} alike", i + 1);
        }
        // a triangular matrix: the reduction is exact up to the Neumann
        // series, the cycle converges fast
        let h = build(
            a.clone(),
            &ReductionSettings {
                max_coarse: 4,
                ..ReductionSettings::default()
            },
        )
        .unwrap();
        assert!(h.levels.len() >= 3);
        assert!(reduction(&a, &h, 2) < 1e-10);
    }

    #[test]
    fn dense_inverse_of_a_nonsymmetric_matrix() {
        let a = rotation(5, 3.0, 0.0);
        let inv = dense_inverse(&a).unwrap();
        let n = a.n_rows;
        for c in 0..n {
            let col: Vec<f64> = (0..n).map(|i| inv[i * n + c]).collect();
            for (i, v) in a.spmv(&col).iter().enumerate() {
                let want = if i == c { 1.0 } else { 0.0 };
                assert!((v - want).abs() < 1e-12);
            }
        }
        // a singular matrix is reported
        let z = Csr::from_rows(2, vec![vec![(0, 1.0), (1, 1.0)], vec![(0, 1.0), (1, 1.0)]]);
        assert!(dense_inverse(&z).is_err());
    }

    #[test]
    fn cycle_converges_from_diffusion_to_convection() {
        // (cells per side, cell Péclet number at the rim, relative shift of
        // the diagonal, bound of the convergence factor per cycle). With
        // the shift of the relaxation cap 0.98 of the coupled algorithm the
        // factor is 0.07 to 0.22 on every mesh; without any shift (closed
        // streamlines, nothing but cross-stream diffusion damps an error
        // along them) it is 0.25 to 0.53.
        for (m, pe, rel, bound) in [
            (32, 0.0, 0.0, 0.5),
            (64, 30.0, 0.0, 0.5),
            (128, 30.0, 0.0, 0.7),
            (128, 5.0, 0.0, 0.7),
            (64, 1.0, 0.02, 0.35),
            (128, 30.0, 0.02, 0.35),
            (256, 30.0, 0.02, 0.35),
            (128, 300.0, 0.02, 0.35),
        ] {
            let mut a = rotation(m, pe, 0.0);
            for i in 0..a.n_rows {
                for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                    if a.col[k] as usize == i {
                        a.val[k] *= 1.0 + rel;
                    }
                }
            }
            let h = build(a.clone(), &ReductionSettings::default()).unwrap();
            let factor = reduction(&a, &h, 10).powf(0.1);
            println!(
                "m {m} Pe {pe} shift {rel}: factor {factor:.2}, levels {}, complexity {:.2}",
                h.levels.len(),
                h.info.operator_complexity
            );
            assert!(
                factor < bound,
                "m {m} Pe {pe} shift {rel}: factor {factor:.2}"
            );
            assert!(
                h.info.operator_complexity < 4.0,
                "m {m} Pe {pe}: operator complexity {}",
                h.info.operator_complexity
            );
            assert!(h.coarse_n <= 64);
        }
    }

    #[test]
    fn matrix_without_couplings() {
        let n = 100;
        let a = Csr::from_rows(
            n,
            (0..n).map(|i| vec![(i as u32, 2.0 + i as f64)]).collect(),
        );
        let h = build(a.clone(), &ReductionSettings::default()).unwrap();
        assert!(reduction(&a, &h, 1) < 1e-14);
    }

    #[test]
    fn setup_does_not_depend_on_threads() {
        let a = rotation(48, 20.0, 0.5);
        let one = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(|| build(a.clone(), &ReductionSettings::default()).unwrap());
        let four = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(|| build(a.clone(), &ReductionSettings::default()).unwrap());
        assert_eq!(one.levels.len(), four.levels.len());
        for (x, y) in one.levels.iter().zip(&four.levels) {
            assert_eq!(x.a, y.a);
            assert_eq!(x.p, y.p);
            assert_eq!(x.r, y.r);
        }
        assert_eq!(one.coarse_inv, four.coarse_inv);
    }
}

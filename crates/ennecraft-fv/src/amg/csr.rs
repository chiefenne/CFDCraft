//! Host-side sparse matrices of the multigrid setup: compressed sparse rows
//! in float64 for the setup algebra, and the sliced-ELL form in which the
//! backends store and apply the hierarchy.
//!
//! Every product and transpose works row by row with a fixed order of the
//! entries (columns ascending), so the setup is deterministic and does not
//! depend on the number of threads (D31).

use crate::disc::CHUNK;
use ennecraft_core::Real;
use rayon::prelude::*;

/// A sparse matrix in compressed-sparse-row form, columns ascending in each
/// row.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Csr {
    /// Rows.
    pub n_rows: usize,
    /// Columns.
    pub n_cols: usize,
    /// Start of each row (length `n_rows + 1`).
    pub row_ptr: Vec<usize>,
    /// Column of each entry.
    pub col: Vec<u32>,
    /// Value of each entry.
    pub val: Vec<f64>,
}

impl Csr {
    /// Builds a matrix from rows of (column, value) pairs; each row is sorted
    /// by column and duplicate columns are summed in their given order.
    pub fn from_rows(n_cols: usize, mut rows: Vec<Vec<(u32, f64)>>) -> Self {
        let n_rows = rows.len();
        // sort (stable: duplicates keep their order) and merge each row
        rows.par_iter_mut().with_min_len(4096).for_each(|r| {
            r.sort_by_key(|e| e.0);
            let mut w = 0;
            for k in 0..r.len() {
                if w > 0 && r[w - 1].0 == r[k].0 {
                    r[w - 1].1 += r[k].1;
                } else {
                    r[w] = r[k];
                    w += 1;
                }
            }
            r.truncate(w);
        });
        let nnz = rows.iter().map(Vec::len).sum();
        let mut row_ptr = Vec::with_capacity(n_rows + 1);
        row_ptr.push(0);
        let mut col = Vec::with_capacity(nnz);
        let mut val = Vec::with_capacity(nnz);
        for r in rows {
            for (c, v) in r {
                col.push(c);
                val.push(v);
            }
            row_ptr.push(col.len());
        }
        Self {
            n_rows,
            n_cols,
            row_ptr,
            col,
            val,
        }
    }

    /// Number of stored entries.
    pub fn nnz(&self) -> usize {
        self.col.len()
    }

    /// Entries of row `i`.
    pub fn row(&self, i: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let (lo, hi) = (self.row_ptr[i], self.row_ptr[i + 1]);
        self.col[lo..hi]
            .iter()
            .zip(&self.val[lo..hi])
            .map(|(&c, &v)| (c as usize, v))
    }

    /// Diagonal entries (zero where a row has none).
    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.n_rows)
            .map(|i| self.row(i).find(|&(c, _)| c == i).map_or(0.0, |e| e.1))
            .collect()
    }

    /// y = A x (rows in parallel, each in column order).
    pub fn spmv(&self, x: &[f64]) -> Vec<f64> {
        (0..self.n_rows)
            .into_par_iter()
            .with_min_len(4096)
            .map(|i| self.row(i).map(|(c, v)| v * x[c]).sum())
            .collect()
    }

    /// The transpose; the entries of a column keep the row order.
    pub fn transpose(&self) -> Csr {
        let mut count = vec![0usize; self.n_cols + 1];
        for &c in &self.col {
            count[c as usize + 1] += 1;
        }
        for i in 0..self.n_cols {
            count[i + 1] += count[i];
        }
        let row_ptr = count.clone();
        let mut next = count;
        let mut col = vec![0u32; self.nnz()];
        let mut val = vec![0.0; self.nnz()];
        for i in 0..self.n_rows {
            for (c, v) in self.row(i) {
                let k = next[c];
                col[k] = i as u32;
                val[k] = v;
                next[c] += 1;
            }
        }
        Csr {
            n_rows: self.n_cols,
            n_cols: self.n_rows,
            row_ptr,
            col,
            val,
        }
    }

    /// The product self · b (Gustavson's row-by-row algorithm with a dense
    /// accumulator per row; the entries of a row are summed in the order of
    /// the row of `self`, then of the rows of `b`).
    pub fn matmul(&self, b: &Csr) -> Csr {
        assert_eq!(self.n_cols, b.n_rows, "matmul: inner dimensions");
        let rows: Vec<Vec<(u32, f64)>> = (0..self.n_rows)
            .into_par_iter()
            .with_min_len(1024)
            .map_init(
                || (vec![usize::MAX; b.n_cols], Vec::<(u32, f64)>::new()),
                |(pos, acc), i| {
                    acc.clear();
                    for (k, a) in self.row(i) {
                        for (j, v) in b.row(k) {
                            if pos[j] == usize::MAX {
                                pos[j] = acc.len();
                                acc.push((j as u32, a * v));
                            } else {
                                acc[pos[j]].1 += a * v;
                            }
                        }
                    }
                    for e in acc.iter() {
                        pos[e.0 as usize] = usize::MAX;
                    }
                    let mut r = acc.clone();
                    r.sort_by_key(|e| e.0);
                    r
                },
            )
            .collect();
        let mut out = Csr {
            n_rows: self.n_rows,
            n_cols: b.n_cols,
            row_ptr: vec![0],
            ..Csr::default()
        };
        for r in rows {
            for (c, v) in r {
                out.col.push(c);
                out.val.push(v);
            }
            out.row_ptr.push(out.col.len());
        }
        out
    }

    /// Largest |a_ij − a_ji| relative to the largest |a_ij|.
    pub fn asymmetry(&self) -> f64 {
        if self.n_rows != self.n_cols {
            return f64::INFINITY;
        }
        let t = self.transpose();
        let scale = self.val.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        if scale == 0.0 {
            return 0.0;
        }
        let mut d = 0.0f64;
        for i in 0..self.n_rows {
            // both rows have ascending columns: merge
            let (mut a, mut b) = (self.row(i).peekable(), t.row(i).peekable());
            loop {
                match (a.peek().copied(), b.peek().copied()) {
                    (None, None) => break,
                    (Some((ca, va)), Some((cb, vb))) if ca == cb => {
                        d = d.max((va - vb).abs());
                        a.next();
                        b.next();
                    }
                    (Some((ca, va)), Some((cb, _))) if ca < cb => {
                        d = d.max(va.abs());
                        a.next();
                    }
                    (Some((_, va)), None) => {
                        d = d.max(va.abs());
                        a.next();
                    }
                    (_, Some((_, vb))) => {
                        d = d.max(vb.abs());
                        b.next();
                    }
                }
            }
        }
        d / scale
    }
}

/// A sliced-ELL matrix (SELL-C-σ, C = [`CHUNK`], σ = 1) owned by the host,
/// in the precision `T` of the hierarchy.
#[derive(Debug, Clone, Default)]
pub struct SellMatrix<T> {
    /// Rows.
    pub n_rows: usize,
    /// Columns.
    pub n_cols: usize,
    /// Slice height C.
    pub chunk: usize,
    /// Start of each slice.
    pub slice_ptr: Vec<u32>,
    /// Slots per row of each slice.
    pub slice_len: Vec<u32>,
    /// Column per slot (padding: column 0, value 0).
    pub col: Vec<u32>,
    /// Value per slot.
    pub val: Vec<T>,
    /// Real entries.
    pub nnz: usize,
}

impl<T: Real> SellMatrix<T> {
    /// Converts a CSR matrix (values rounded to `T`).
    pub fn from_csr(a: &Csr) -> Self {
        let n = a.n_rows;
        let n_slices = n.div_ceil(CHUNK).max(1);
        let mut slice_ptr = Vec::with_capacity(n_slices);
        let mut slice_len = Vec::with_capacity(n_slices);
        let mut total = 0usize;
        for s in 0..n_slices {
            let len = (s * CHUNK..((s + 1) * CHUNK).min(n))
                .map(|r| a.row_ptr[r + 1] - a.row_ptr[r])
                .max()
                .unwrap_or(0)
                .max(1);
            slice_ptr.push(total as u32);
            slice_len.push(len as u32);
            total += len * CHUNK;
        }
        let mut col = vec![0u32; total];
        let mut val = vec![T::ZERO; total];
        for r in 0..n {
            let s = r / CHUNK;
            for (j, (c, v)) in a.row(r).enumerate() {
                let k = slice_ptr[s] as usize + j * CHUNK + r % CHUNK;
                col[k] = c as u32;
                val[k] = T::from_f64(v);
            }
        }
        Self {
            n_rows: n,
            n_cols: a.n_cols,
            chunk: CHUNK,
            slice_ptr,
            slice_len,
            col,
            val,
            nnz: a.nnz(),
        }
    }

    /// Borrowed form for the kernel bodies.
    pub fn view(&self) -> ennecraft_kernels::amg::Sell<'_, T> {
        ennecraft_kernels::amg::Sell {
            chunk: self.chunk,
            slice_ptr: &self.slice_ptr,
            slice_len: &self.slice_len,
            col: &self.col,
            val: &self.val,
        }
    }

    /// Stored slots including padding.
    pub fn stored(&self) -> usize {
        self.col.len()
    }

    /// Bytes of the stored matrix (values and columns).
    pub fn bytes(&self) -> usize {
        self.stored() * (T::BYTES + 4) + 8 * self.slice_ptr.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lap1d(n: usize) -> Csr {
        let rows = (0..n)
            .map(|i| {
                let mut r = vec![(i as u32, 2.0)];
                if i > 0 {
                    r.push((i as u32 - 1, -1.0));
                }
                if i + 1 < n {
                    r.push((i as u32 + 1, -1.0));
                }
                r
            })
            .collect();
        Csr::from_rows(n, rows)
    }

    #[test]
    fn transpose_matmul_symmetry() {
        let a = lap1d(7);
        assert_eq!(a.transpose(), a);
        assert_eq!(a.asymmetry(), 0.0);
        let a2 = a.matmul(&a);
        assert_eq!(a2.row(3).collect::<Vec<_>>().len(), 5);
        let x: Vec<f64> = (0..7).map(|i| i as f64).collect();
        let y1 = a2.spmv(&x);
        let y2 = a.spmv(&a.spmv(&x));
        assert_eq!(y1, y2);
    }

    #[test]
    fn sell_round_trip() {
        let a = lap1d(70);
        let s = SellMatrix::<f64>::from_csr(&a);
        let x: Vec<f64> = (0..70).map(|i| (i as f64).sin()).collect();
        let y = a.spmv(&x);
        for (i, yi) in y.iter().enumerate() {
            // SAFETY: built by from_csr.
            let v = unsafe { s.view().row(i, &x) };
            assert_eq!(v, *yi);
        }
    }
}

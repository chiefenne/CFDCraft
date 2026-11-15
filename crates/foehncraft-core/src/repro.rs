//! Reductions with a fixed summation order (requirement D31).
//!
//! Floating-point addition is not associative, so a sum depends on its order.
//! The functions here fix the order as a function of the data length only:
//! the result is bitwise the same from run to run and does not depend on the
//! number of threads that compute it (see `foehncraft_parallel::repro`).

use crate::Real;

/// Length of the blocks that are summed sequentially before the pairwise
/// tree starts. Part of the definition of the summation order: changing it
/// changes results in the last bits.
pub const BLOCK: usize = 256;

/// Sum of `x` in a fixed order: sequential sums over blocks of [`BLOCK`]
/// values, combined by a pairwise (binary) tree.
///
/// The rounding error grows like O(log n) instead of O(n) for the plain loop
/// (Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., §4.2).
pub fn pairwise_sum<T: Real>(x: &[T]) -> T {
    if x.len() <= BLOCK {
        let mut s = T::ZERO;
        for &v in x {
            s += v;
        }
        return s;
    }
    let half = split_point(x.len());
    pairwise_sum(&x[..half]) + pairwise_sum(&x[half..])
}

/// Split point of the pairwise tree: the largest multiple of [`BLOCK`] times
/// a power of two below `n`. Splitting on block boundaries lets a caller sum
/// whole blocks independently (in parallel) and combine the block sums with
/// [`combine_block_sums`] to the bitwise same result.
#[inline]
pub fn split_point(n: usize) -> usize {
    let blocks = n.div_ceil(BLOCK);
    debug_assert!(blocks > 1);
    // largest power of two strictly below `blocks`
    let p = 1usize << (usize::BITS - 1 - (blocks - 1).leading_zeros());
    p * BLOCK
}

/// Combines the sums of consecutive blocks of [`BLOCK`] values (the last one
/// possibly shorter) in the tree order of [`pairwise_sum`].
pub fn combine_block_sums<T: Real>(block_sums: &[T]) -> T {
    match block_sums.len() {
        0 => T::ZERO,
        1 => block_sums[0],
        n => {
            let p = 1usize << (usize::BITS - 1 - (n - 1).leading_zeros());
            combine_block_sums(&block_sums[..p]) + combine_block_sums(&block_sums[p..])
        }
    }
}

/// Sequential sum of one block, the leaf of the pairwise tree.
#[inline]
pub fn block_sum<T: Real>(x: &[T]) -> T {
    let mut s = T::ZERO;
    for &v in x {
        s += v;
    }
    s
}

#[cfg(test)]
mod tests {
    extern crate std;
    use super::*;
    use std::vec::Vec;

    fn data(n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| 1.0 / (1.0 + i as f64) + 1e-3 * ((i % 7) as f64))
            .collect()
    }

    #[test]
    fn small_and_empty() {
        assert_eq!(pairwise_sum::<f64>(&[]), 0.0);
        assert_eq!(pairwise_sum(&[1.0, 2.0, 3.0]), 6.0);
    }

    #[test]
    fn block_decomposition_gives_identical_bits() {
        for n in [1, 255, 256, 257, 1000, 4096, 10_001, 70_000, 1_000_003] {
            let v = data(n);
            let direct = pairwise_sum(&v);
            let sums: Vec<f64> = v.chunks(BLOCK).map(block_sum).collect();
            let combined = combine_block_sums(&sums);
            assert_eq!(direct.to_bits(), combined.to_bits(), "n = {n}");
        }
    }

    #[test]
    fn pairwise_is_accurate() {
        let v = data(1_000_000);
        // Kahan-Babuska (Neumaier) compensated sum as the reference.
        let (mut s, mut c) = (0.0f64, 0.0f64);
        for &x in &v {
            let t = s + x;
            c += if s.abs() >= x.abs() {
                (s - t) + x
            } else {
                (x - t) + s
            };
            s = t;
        }
        let exact = s + c;
        let pw = pairwise_sum(&v);
        assert!((pw - exact).abs() <= 4.0 * f64::EPSILON * exact.abs());
    }
}

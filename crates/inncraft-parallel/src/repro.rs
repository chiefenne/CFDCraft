//! Multi-threaded reductions with a result independent of the thread count.
//!
//! The data are cut into blocks of [`inncraft_core::repro::BLOCK`] values,
//! each block is summed sequentially (in parallel over blocks), and the block
//! sums are combined in the fixed tree order of
//! [`inncraft_core::repro::pairwise_sum`]. The result is therefore bitwise
//! equal to the sequential `pairwise_sum` for any number of threads (D31).

use inncraft_core::Real;
use inncraft_core::repro::{BLOCK, block_sum, combine_block_sums};
use rayon::prelude::*;

/// Sum of `x`, bitwise equal to `inncraft_core::repro::pairwise_sum(x)` for
/// any number of threads.
pub fn par_sum<T: Real>(x: &[T]) -> T {
    let sums: Vec<T> = x.par_chunks(BLOCK).map(block_sum).collect();
    combine_block_sums(&sums)
}

/// Dot product of `a` and `b` in the same fixed order as [`par_sum`].
pub fn par_dot<T: Real>(a: &[T], b: &[T]) -> T {
    assert_eq!(a.len(), b.len());
    let sums: Vec<T> = a
        .par_chunks(BLOCK)
        .zip(b.par_chunks(BLOCK))
        .map(|(x, y)| {
            let mut s = T::ZERO;
            for (&u, &v) in x.iter().zip(y) {
                s += u * v;
            }
            s
        })
        .collect();
    combine_block_sums(&sums)
}

#[cfg(test)]
mod tests {
    use super::*;
    use inncraft_core::repro::pairwise_sum;

    fn data(n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| inncraft_core::rng::uniform([i as u64, 1, 2]) - 0.3)
            .collect()
    }

    #[test]
    fn independent_of_thread_count() {
        let v = data(1_234_567);
        let reference = pairwise_sum(&v).to_bits();
        for threads in [1, 2, 3, 4, 7] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let s = pool.install(|| par_sum(&v));
            assert_eq!(s.to_bits(), reference, "{threads} threads");
        }
    }

    #[test]
    fn dot_matches_sum_of_products() {
        let a = data(100_001);
        let b: Vec<f64> = data(100_001).iter().map(|x| x * 0.5 + 1.0).collect();
        let prod: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x * y).collect();
        assert_eq!(par_dot(&a, &b).to_bits(), pairwise_sum(&prod).to_bits());
    }
}

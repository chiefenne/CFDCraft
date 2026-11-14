//! Counter-based random numbers (requirement D47).
//!
//! A random number is a pure function of a key (for example particle
//! identity, time step and event number), so its value does not depend on
//! which thread, GPU or rank draws it, nor on the order of the draws. This is
//! what keeps Lagrangian particle runs reproducible (D31, D47).
//!
//! The mixing function is SplitMix64 (Steele, Lea & Flood, *Fast splittable
//! pseudorandom number generators*, OOPSLA 2014), applied in a chain over the
//! key words. It passes the usual statistical batteries as a hash of
//! consecutive counters; a Philox-type generator (Salmon et al., SC 2011) can
//! replace it later behind the same functions.

/// SplitMix64 finaliser: a bijective mixing of a 64-bit word.
#[inline(always)]
pub const fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// 64 random bits for a key of three words (identity, step, event).
#[inline(always)]
pub const fn bits(key: [u64; 3]) -> u64 {
    splitmix64(key[2] ^ splitmix64(key[1] ^ splitmix64(key[0])))
}

/// Uniform number in [0, 1) with 53 random bits, for a key of three words.
#[inline(always)]
pub fn uniform(key: [u64; 3]) -> f64 {
    (bits(key) >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_key_same_value() {
        assert_eq!(uniform([7, 3, 1]).to_bits(), uniform([7, 3, 1]).to_bits());
        assert_ne!(bits([7, 3, 1]), bits([7, 3, 2]));
        assert_ne!(bits([7, 3, 1]), bits([3, 7, 1]));
    }

    #[test]
    fn mean_and_range() {
        let n = 100_000u64;
        let mut s = 0.0;
        for i in 0..n {
            let u = uniform([i, 0, 0]);
            assert!((0.0..1.0).contains(&u));
            s += u;
        }
        let mean = s / n as f64;
        // standard error of the mean of U(0,1) is 1/sqrt(12 n) ~ 9e-4
        assert!((mean - 0.5).abs() < 5e-3, "mean = {mean}");
    }
}

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

/// 64 random bits for a key of three words (identity, step, event) of
/// the stream [`stream::MESH`]; the same as
/// `keyed([stream::MESH, identity, step, event])`.
#[inline(always)]
pub const fn bits(key: [u64; 3]) -> u64 {
    splitmix64(key[2] ^ splitmix64(key[1] ^ splitmix64(key[0])))
}

/// The streams of the generator (D121): the first word of a key names the
/// module and purpose, so two models that draw for the same identity, step
/// and event (breakup and dispersion of one parcel) never get the same
/// number. Every stream of the code is listed here, once; a number is never
/// reused for another purpose.
pub mod stream {
    /// The mesh: perturbation of the test meshes, random renumbering,
    /// start cells of the point location tests. Stream 0 keeps the numbers
    /// of the three-word keys of [`bits`](super::bits).
    pub const MESH: u64 = 0;
    /// Spray injection: parcel sizes, directions and positions (L12).
    pub const SPRAY_INJECTION: u64 = 1;
    /// Spray breakup (L12).
    pub const SPRAY_BREAKUP: u64 = 2;
    /// Turbulent dispersion of parcels (L12).
    pub const SPRAY_DISPERSION: u64 = 3;
    /// Wall impingement and splash of parcels (D123).
    pub const SPRAY_WALL: u64 = 4;
    /// Parcel collision (L12).
    pub const SPRAY_COLLISION: u64 = 5;
    /// Start vectors of the multigrid eigenvalue estimate (D93 (5)).
    pub const AMG_START: u64 = 6;
    /// Monte Carlo radiation: photon bundles (L15).
    pub const RADIATION: u64 = 7;
    /// Every stream with its name, for the documentation and the test that
    /// keeps the numbers distinct.
    pub const ALL: &[(u64, &str)] = &[
        (MESH, "mesh"),
        (SPRAY_INJECTION, "spray-injection"),
        (SPRAY_BREAKUP, "spray-breakup"),
        (SPRAY_DISPERSION, "spray-dispersion"),
        (SPRAY_WALL, "spray-wall"),
        (SPRAY_COLLISION, "spray-collision"),
        (AMG_START, "amg-start"),
        (RADIATION, "radiation"),
    ];
}

/// 64 random bits for a key of four words `[stream, identity, step, event]`
/// (D121). Stream [`stream::MESH`] gives the numbers of [`bits`].
#[inline(always)]
pub const fn keyed(key: [u64; 4]) -> u64 {
    let s = if key[0] == stream::MESH {
        0
    } else {
        splitmix64(key[0] ^ 0x5354_5245_414D_0000)
    };
    splitmix64(key[3] ^ splitmix64(key[2] ^ splitmix64(key[1] ^ s)))
}

/// Uniform number in [0, 1) with 53 random bits, for a key of four words.
#[inline(always)]
pub fn uniform_keyed(key: [u64; 4]) -> f64 {
    (keyed(key) >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
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
    fn streams() {
        // stream 0 is the three-word generator, bit for bit
        for i in 0..1000u64 {
            assert_eq!(keyed([stream::MESH, i, 3, 7]), bits([i, 3, 7]));
        }
        // the same identity, step and event differ between streams
        for &(a, _) in stream::ALL {
            for &(b, _) in stream::ALL {
                if a != b {
                    assert_ne!(keyed([a, 11, 2, 0]), keyed([b, 11, 2, 0]));
                }
            }
        }
        let mut ids: [u64; 8] = [0; 8];
        for (i, &(s, _)) in stream::ALL.iter().enumerate() {
            ids[i] = s;
        }
        ids.sort_unstable();
        assert!(
            ids.windows(2).all(|w| w[0] != w[1]),
            "stream numbers must be distinct"
        );
        let u = uniform_keyed([stream::SPRAY_BREAKUP, 5, 1, 0]);
        assert!((0.0..1.0).contains(&u));
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

//! The generic floating-point type.

use core::fmt::{Debug, Display};
use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

/// Floating-point type of the solver.
///
/// All numerical code is generic over `Real`, so that a float32 or
/// mixed-precision mode (P2) needs no second code base. The reference
/// solution path uses `f64` (D7). Only operations that `core` provides on
/// every target are part of the trait, so it is usable inside GPU kernels.
pub trait Real:
    Copy
    + Send
    + Sync
    + Default
    + Debug
    + Display
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + 'static
{
    /// Additive identity.
    const ZERO: Self;
    /// Multiplicative identity.
    const ONE: Self;
    /// Size in bytes, used by the bandwidth model of the benchmarks.
    const BYTES: usize;
    /// Short name of the type (`"f64"`, `"f32"`).
    const NAME: &'static str;

    /// Conversion from `f64` (rounding for narrower types).
    fn from_f64(x: f64) -> Self;
    /// Conversion to `f64`.
    fn to_f64(self) -> f64;

    /// Absolute value, without the `std` math library.
    #[inline(always)]
    fn abs(self) -> Self {
        if self < Self::ZERO { -self } else { self }
    }

    /// Larger of two values, branch-free (IEEE 754 maxNum: a NaN operand
    /// is ignored).
    fn max_of(self, other: Self) -> Self;

    /// Smaller of two values, branch-free (IEEE 754 minNum).
    fn min_of(self, other: Self) -> Self;
}

impl Real for f64 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const BYTES: usize = 8;
    const NAME: &'static str = "f64";

    #[inline(always)]
    fn from_f64(x: f64) -> Self {
        x
    }
    #[inline(always)]
    fn to_f64(self) -> f64 {
        self
    }
    #[inline(always)]
    fn max_of(self, other: Self) -> Self {
        f64::max(self, other)
    }
    #[inline(always)]
    fn min_of(self, other: Self) -> Self {
        f64::min(self, other)
    }
}

impl Real for f32 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const BYTES: usize = 4;
    const NAME: &'static str = "f32";

    #[inline(always)]
    fn from_f64(x: f64) -> Self {
        x as f32
    }
    #[inline(always)]
    fn to_f64(self) -> f64 {
        self as f64
    }
    #[inline(always)]
    fn max_of(self, other: Self) -> Self {
        f32::max(self, other)
    }
    #[inline(always)]
    fn min_of(self, other: Self) -> Self {
        f32::min(self, other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn generic_abs<T: Real>(x: T) -> T {
        x.abs()
    }

    #[test]
    fn abs_and_extrema() {
        assert_eq!(generic_abs(-2.5_f64), 2.5);
        assert_eq!(generic_abs(-2.5_f32), 2.5);
        assert_eq!(1.0_f64.max_of(3.0), 3.0);
        assert_eq!(1.0_f64.min_of(3.0), 1.0);
    }

    #[test]
    fn conversions() {
        assert_eq!(f32::from_f64(0.1).to_f64(), 0.1_f32 as f64);
        assert_eq!(<f64 as Real>::BYTES, 8);
    }
}

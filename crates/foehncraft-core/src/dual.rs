//! Forward-mode dual numbers (requirement D59): a value and `N` directional
//! derivatives, `x + Σ_k d_k ε_k` with `ε_j ε_k = 0`.
//!
//! [`Dual`] implements [`Real`], so every kernel body generic over `Real`
//! also runs on dual numbers and returns, next to its result, the exact
//! derivatives of that result in `N` seeded directions: the Jacobian
//! contributions come from the same source as the residual, by
//! instantiation, with no hand-derived or copied code.
//!
//! **Non-differentiable operations.** Comparisons look at the value only.
//! [`Real::max_of`] and [`Real::min_of`] return the operand that the value
//! selects, together with its derivatives; at a tie the first operand wins.
//! So a switch (upwind direction, limiter branch, clipping) is *frozen* at
//! the state of the linearisation, which is the derivative of the branch
//! that the primal computation takes (exact away from the switching
//! surfaces). [`Real::abs`] follows from the comparison in the same way.
//!
//! The derivative parts carry the same IEEE operations as the values, so a
//! dual evaluation is reproducible like the primal one (D31).

use crate::Real;
use core::fmt::{Debug, Display, Formatter};
use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

/// A value and `N` derivatives.
#[derive(Clone, Copy)]
pub struct Dual<const N: usize> {
    /// Value.
    pub v: f64,
    /// Derivatives in the `N` seeded directions.
    pub d: [f64; N],
}

impl<const N: usize> Dual<N> {
    /// A constant (all derivatives zero).
    #[inline(always)]
    pub const fn constant(v: f64) -> Self {
        Self { v, d: [0.0; N] }
    }

    /// A value with the given derivatives.
    #[inline(always)]
    pub const fn new(v: f64, d: [f64; N]) -> Self {
        Self { v, d }
    }

    /// The independent variable of direction `k`: derivative 1 there, 0 in
    /// the other directions.
    #[inline(always)]
    pub fn variable(v: f64, k: usize) -> Self {
        let mut d = [0.0; N];
        d[k] = 1.0;
        Self { v, d }
    }

    /// The value with derivatives scaled by `s`, the derivatives of `self`
    /// multiplied by the chain factor `s` (f(x) with f'(x) = s).
    #[inline(always)]
    fn chain(v: f64, s: f64, x: &Self) -> Self {
        let mut d = [0.0; N];
        let mut k = 0;
        while k < N {
            d[k] = s * x.d[k];
            k += 1;
        }
        Self { v, d }
    }
}

impl<const N: usize> Default for Dual<N> {
    fn default() -> Self {
        Self::constant(0.0)
    }
}

impl<const N: usize> Debug for Dual<N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "{:?} + {:?}ε", self.v, self.d)
    }
}

impl<const N: usize> Display for Dual<N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        Display::fmt(&self.v, f)
    }
}

/// Equality of the values (the derivatives do not take part).
impl<const N: usize> PartialEq for Dual<N> {
    fn eq(&self, o: &Self) -> bool {
        self.v == o.v
    }
}

/// Order of the values (the derivatives do not take part).
impl<const N: usize> PartialOrd for Dual<N> {
    fn partial_cmp(&self, o: &Self) -> Option<core::cmp::Ordering> {
        self.v.partial_cmp(&o.v)
    }
}

macro_rules! zip {
    ($a:expr, $b:expr, |$x:ident, $y:ident| $e:expr) => {{
        let mut d = [0.0; N];
        let mut k = 0;
        while k < N {
            let ($x, $y) = ($a.d[k], $b.d[k]);
            d[k] = $e;
            k += 1;
        }
        d
    }};
}

impl<const N: usize> Add for Dual<N> {
    type Output = Self;
    #[inline(always)]
    fn add(self, o: Self) -> Self {
        Self::new(self.v + o.v, zip!(self, o, |a, b| a + b))
    }
}

impl<const N: usize> Sub for Dual<N> {
    type Output = Self;
    #[inline(always)]
    fn sub(self, o: Self) -> Self {
        Self::new(self.v - o.v, zip!(self, o, |a, b| a - b))
    }
}

impl<const N: usize> Mul for Dual<N> {
    type Output = Self;
    #[inline(always)]
    fn mul(self, o: Self) -> Self {
        let (u, v) = (self.v, o.v);
        Self::new(u * v, zip!(self, o, |a, b| a * v + u * b))
    }
}

impl<const N: usize> Div for Dual<N> {
    type Output = Self;
    #[inline(always)]
    fn div(self, o: Self) -> Self {
        let q = self.v / o.v;
        let inv = 1.0 / o.v;
        Self::new(q, zip!(self, o, |a, b| (a - q * b) * inv))
    }
}

impl<const N: usize> Neg for Dual<N> {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self::chain(-self.v, -1.0, &self)
    }
}

impl<const N: usize> AddAssign for Dual<N> {
    #[inline(always)]
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}

impl<const N: usize> SubAssign for Dual<N> {
    #[inline(always)]
    fn sub_assign(&mut self, o: Self) {
        *self = *self - o;
    }
}

impl<const N: usize> MulAssign for Dual<N> {
    #[inline(always)]
    fn mul_assign(&mut self, o: Self) {
        *self = *self * o;
    }
}

impl<const N: usize> Real for Dual<N> {
    const ZERO: Self = Self::constant(0.0);
    const ONE: Self = Self::constant(1.0);
    const BYTES: usize = 8 * (N + 1);
    const NAME: &'static str = "dual";

    #[inline(always)]
    fn from_f64(x: f64) -> Self {
        Self::constant(x)
    }

    #[inline(always)]
    fn to_f64(self) -> f64 {
        self.v
    }

    /// The operand the value selects (IEEE maxNum: a NaN value is
    /// ignored), with its derivatives; `self` at a tie.
    #[inline(always)]
    fn max_of(self, o: Self) -> Self {
        if o.v > self.v || self.v.is_nan() {
            o
        } else {
            self
        }
    }

    /// As [`Self::max_of`] for the smaller value.
    #[inline(always)]
    fn min_of(self, o: Self) -> Self {
        if o.v < self.v || self.v.is_nan() {
            o
        } else {
            self
        }
    }
    #[inline(always)]
    fn sqrt_of(self) -> Self {
        let r = crate::math::sqrt(self.v);
        Self::chain(r, 0.5 / r, &self)
    }

    #[inline(always)]
    fn exp_of(self) -> Self {
        let e = crate::math::exp(self.v);
        Self::chain(e, e, &self)
    }

    #[inline(always)]
    fn ln_of(self) -> Self {
        Self::chain(crate::math::ln(self.v), 1.0 / self.v, &self)
    }

    #[inline(always)]
    fn tanh_of(self) -> Self {
        let t = crate::math::tanh(self.v);
        Self::chain(t, 1.0 - t * t, &self)
    }

    #[inline(always)]
    fn cos_acos3_of(self) -> Self {
        Self::chain(
            crate::math::cos_acos3(self.v),
            crate::math::cos_acos3_derivative(self.v),
            &self,
        )
    }

    #[inline(always)]
    fn sin_of(self) -> Self {
        Self::chain(crate::math::sin(self.v), crate::math::cos(self.v), &self)
    }

    #[inline(always)]
    fn atan_of(self) -> Self {
        Self::chain(
            crate::math::atan(self.v),
            1.0 / (1.0 + self.v * self.v),
            &self,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type D = Dual<2>;

    #[test]
    fn arithmetic_derivatives() {
        let x = D::variable(3.0, 0);
        let y = D::variable(2.0, 1);
        // f = x y / (x − y) + x²
        let f = x * y / (x - y) + x * x;
        assert_eq!(f.v, 6.0 + 9.0);
        // ∂f/∂x = (y (x − y) − x y) / (x − y)² + 2x = −y²/(x−y)² + 2x = 2
        // ∂f/∂y = (x (x − y) + x y) / (x − y)² = x²/(x−y)² = 9
        assert!((f.d[0] - 2.0).abs() < 1e-14);
        assert!((f.d[1] - 9.0).abs() < 1e-14);
        let g = -x + D::from_f64(1.0);
        assert_eq!((g.v, g.d), (-2.0, [-1.0, 0.0]));
    }

    #[test]
    fn elementary_functions_carry_derivatives() {
        let x = D::variable(0.7, 0);
        let h = 1e-6;
        for (f, g) in [
            (x.sqrt_of(), crate::math::sqrt as fn(f64) -> f64),
            (x.exp_of(), crate::math::exp),
            (x.ln_of(), crate::math::ln),
            (x.tanh_of(), crate::math::tanh),
        ] {
            let fd = (g(0.7 + h) - g(0.7 - h)) / (2.0 * h);
            assert_eq!(f.v, g(0.7));
            assert!((f.d[0] - fd).abs() < 1e-8, "{} vs {}", f.d[0], fd);
            assert_eq!(f.d[1], 0.0);
        }
    }

    #[test]
    fn switches_are_frozen() {
        let x = D::variable(-1.0, 0);
        let y = D::variable(2.0, 1);
        assert_eq!(x.max_of(y).d, [0.0, 1.0]);
        assert_eq!(x.min_of(y).d, [1.0, 0.0]);
        assert_eq!(x.abs().v, 1.0);
        assert_eq!(x.abs().d, [-1.0, 0.0]);
        // tie: the first operand
        let z = D::variable(2.0, 0);
        assert_eq!(y.max_of(z).d, [0.0, 1.0]);
        // NaN is ignored like f64::max
        let n = D::constant(f64::NAN);
        assert_eq!(n.max_of(y).v, 2.0);
        assert_eq!(y.max_of(n).v, 2.0);
    }
}

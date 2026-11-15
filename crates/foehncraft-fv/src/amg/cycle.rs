//! The V-cycle, written once over the level operations of a backend.
//!
//! The backend holds, for every level l, the vectors x_l (correction),
//! b_l (right-hand side), w_l (residual) and d_l (smoother update), in the
//! precision of the hierarchy. One cycle z = M⁻¹ r:
//!
//! 1. b_0 = r (rounded to the hierarchy precision);
//! 2. on levels 0 .. L−1: pre-smooth from x_l = 0, w_l = b_l − A_l x_l,
//!    b_{l+1} = R_l w_l;
//! 3. x_L = A_L⁻¹ b_L (dense inverse);
//! 4. on levels L−1 .. 0: x_l += P_l x_{l+1}, post-smooth (the same steps);
//! 5. z = x_0 (back to float64).
//!
//! With the same polynomial smoother before and after the coarse correction
//! the cycle is a symmetric operator (exactly so in float64; in float32 up
//! to rounding, which the flexible CG absorbs).

use crate::backend::V;

/// Level operations a backend provides for the V-cycle.
pub trait LevelOps {
    /// Number of levels with a matrix (the coarsest dense level not
    /// counted).
    fn amg_levels(&self) -> usize;
    /// Smoother steps (c₁, c₂) of level `l`.
    fn amg_steps(&self, l: usize) -> Vec<(f64, f64)>;
    /// b_0 = r.
    fn amg_load_rhs(&mut self, r: V) -> anyhow::Result<()>;
    /// z = x_0.
    fn amg_store(&mut self, z: V) -> anyhow::Result<()>;
    /// First smoother step from zero: d_l = x_l = c₂ D⁻¹ b_l.
    fn amg_smooth_first(&mut self, l: usize, c2: f64) -> anyhow::Result<()>;
    /// Smoother step: d_l = c₁ d_l + c₂ D⁻¹(b_l − A_l x_l); x_l += d_l.
    fn amg_smooth_step(&mut self, l: usize, c1: f64, c2: f64) -> anyhow::Result<()>;
    /// w_l = b_l − A_l x_l; b_{l+1} = R_l w_l.
    fn amg_restrict(&mut self, l: usize) -> anyhow::Result<()>;
    /// x_L = A_L⁻¹ b_L.
    fn amg_coarse(&mut self) -> anyhow::Result<()>;
    /// x_l += P_l x_{l+1}.
    fn amg_prolong(&mut self, l: usize) -> anyhow::Result<()>;
}

/// One V-cycle z = M⁻¹ r.
pub fn vcycle<B: LevelOps + ?Sized>(be: &mut B, r: V, z: V) -> anyhow::Result<()> {
    let nl = be.amg_levels();
    be.amg_load_rhs(r)?;
    for l in 0..nl {
        let steps = be.amg_steps(l);
        be.amg_smooth_first(l, steps[0].1)?;
        for &(c1, c2) in &steps[1..] {
            be.amg_smooth_step(l, c1, c2)?;
        }
        be.amg_restrict(l)?;
    }
    be.amg_coarse()?;
    for l in (0..nl).rev() {
        be.amg_prolong(l)?;
        for &(c1, c2) in &be.amg_steps(l) {
            be.amg_smooth_step(l, c1, c2)?;
        }
    }
    be.amg_store(z)
}

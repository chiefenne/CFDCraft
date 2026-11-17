//! Strain-rate and rotation tensors of the velocity gradient and their
//! invariants (stage 6 specification 2.1), and the dimensionless,
//! Galilean- and rotation-invariant features of D166.
//!
//! The velocity gradient is `g[i][j] = ∂u_i/∂x_j`.
//!
//! - S_ij = ½(g_ij + g_ji), W_ij = ½(g_ij − g_ji);
//! - S = (2 S_ij S_ij)^½ (the strain-rate invariant of SST 2003),
//!   Ω = (2 W_ij W_ij)^½ (the vorticity magnitude);
//! - the five invariants of Pope (1975) of the normalised tensors
//!   Ŝ = τ S, Ŵ = τ W with the turbulent time scale τ:
//!   λ₁ = tr Ŝ², λ₂ = tr Ŵ², λ₃ = tr Ŝ³, λ₄ = tr Ŵ²Ŝ, λ₅ = tr Ŵ²Ŝ²
//!   (S. B. Pope, *A more general effective-viscosity hypothesis*, J. Fluid
//!   Mech. 72, 331–340 (1975); the minimal integrity basis of a symmetric
//!   and an antisymmetric tensor in 3D);
//! - the Q criterion Q = ½(W_ij W_ij − S_ij S_ij) (D72).

use foehncraft_core::Real;

/// A 3 × 3 tensor, row i, column j.
pub type Tensor<T> = [[T; 3]; 3];

/// S_ij and W_ij of the velocity gradient `g` (g_ij = ∂u_i/∂x_j).
#[inline(always)]
pub fn strain_rotation<T: Real>(g: &Tensor<T>) -> (Tensor<T>, Tensor<T>) {
    let half = T::from_f64(0.5);
    let mut s = [[T::ZERO; 3]; 3];
    let mut w = [[T::ZERO; 3]; 3];
    let mut i = 0;
    while i < 3 {
        let mut j = 0;
        while j < 3 {
            s[i][j] = half * (g[i][j] + g[j][i]);
            w[i][j] = half * (g[i][j] - g[j][i]);
            j += 1;
        }
        i += 1;
    }
    (s, w)
}

/// A_ij B_ij.
#[inline(always)]
pub fn contract<T: Real>(a: &Tensor<T>, b: &Tensor<T>) -> T {
    let mut r = T::ZERO;
    let mut i = 0;
    while i < 3 {
        let mut j = 0;
        while j < 3 {
            r += a[i][j] * b[i][j];
            j += 1;
        }
        i += 1;
    }
    r
}

/// The matrix product A B.
#[inline(always)]
pub fn product<T: Real>(a: &Tensor<T>, b: &Tensor<T>) -> Tensor<T> {
    let mut c = [[T::ZERO; 3]; 3];
    let mut i = 0;
    while i < 3 {
        let mut j = 0;
        while j < 3 {
            let mut k = 0;
            while k < 3 {
                c[i][j] += a[i][k] * b[k][j];
                k += 1;
            }
            j += 1;
        }
        i += 1;
    }
    c
}

/// tr A.
#[inline(always)]
pub fn trace<T: Real>(a: &Tensor<T>) -> T {
    a[0][0] + a[1][1] + a[2][2]
}

/// The invariants S² = 2 S_ij S_ij and Ω² = 2 W_ij W_ij of a velocity
/// gradient (their squares: no square root where only the square is
/// needed).
#[inline(always)]
pub fn s2_omega2<T: Real>(g: &Tensor<T>) -> (T, T) {
    let (s, w) = strain_rotation(g);
    let two = T::from_f64(2.0);
    (two * contract(&s, &s), two * contract(&w, &w))
}

/// The five invariants λ₁ … λ₅ of Pope (1975) of Ŝ = τ S and Ŵ = τ W.
#[inline(always)]
pub fn pope<T: Real>(g: &Tensor<T>, tau: T) -> [T; 5] {
    let (mut s, mut w) = strain_rotation(g);
    let mut i = 0;
    while i < 3 {
        let mut j = 0;
        while j < 3 {
            s[i][j] *= tau;
            w[i][j] *= tau;
            j += 1;
        }
        i += 1;
    }
    let s2 = product(&s, &s);
    let w2 = product(&w, &w);
    [
        trace(&s2),
        trace(&w2),
        trace(&product(&s2, &s)),
        trace(&product(&w2, &s)),
        trace(&product(&w2, &s2)),
    ]
}

/// The Q criterion ½(W_ij W_ij − S_ij S_ij) (1/s²).
#[inline(always)]
pub fn q_criterion<T: Real>(g: &Tensor<T>) -> T {
    let (s, w) = strain_rotation(g);
    T::from_f64(0.5) * (contract(&w, &w) - contract(&s, &s))
}

/// The shielding function of delayed DES (Spalart, Deck, Shur, Squires,
/// Strelets & Travin, Theor. Comput. Fluid Dyn. 20, 181–195 (2006)):
/// f_d = 1 − tanh((8 r_d)³), r_d = (ν_t + ν)/(√(g_ij g_ij) κ² d²),
/// κ = 0.41. 0 inside attached boundary layers, 1 away from walls.
#[inline(always)]
pub fn shielding<T: Real>(g: &Tensor<T>, nu_t: T, nu: T, d: T) -> T {
    let kappa = T::from_f64(0.41);
    let gg = contract(g, g).sqrt_of().max_of(T::from_f64(1e-300));
    let rd = (nu_t + nu) / (gg * kappa * kappa * d * d);
    let a = T::from_f64(8.0) * rd;
    T::ONE - (a * a * a).tanh_of()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn shear(a: f64) -> Tensor<f64> {
        // u = a y: ∂u/∂y = a
        let mut g = [[0.0; 3]; 3];
        g[0][1] = a;
        g
    }

    #[test]
    fn pure_shear_has_equal_strain_and_vorticity() {
        let (s2, w2) = s2_omega2(&shear(3.0));
        assert!((s2 - 9.0).abs() < 1e-14);
        assert!((w2 - 9.0).abs() < 1e-14);
        assert!(q_criterion(&shear(3.0)).abs() < 1e-14);
    }

    #[test]
    fn pope_invariants_of_shear_and_rotation() {
        // pure shear: Ŝ = τa/2 (e_xy + e_yx), Ŵ = τa/2 (e_xy − e_yx)
        let l = pope(&shear(2.0), 0.5);
        // τa/2 = 0.5: tr Ŝ² = 2·0.25, tr Ŵ² = −2·0.25, tr Ŝ³ = 0,
        // Ŵ² = −0.25 diag(1,1,0): tr Ŵ²Ŝ = 0, tr Ŵ²Ŝ² = −0.25·0.25·2
        assert!((l[0] - 0.5).abs() < 1e-15);
        assert!((l[1] + 0.5).abs() < 1e-15);
        assert!(l[2].abs() < 1e-15);
        assert!(l[3].abs() < 1e-15);
        assert!((l[4] + 0.125).abs() < 1e-15);
        // solid-body rotation: no strain
        let mut g = [[0.0; 3]; 3];
        g[0][1] = -1.0;
        g[1][0] = 1.0;
        let (s2, w2) = s2_omega2(&g);
        assert_eq!(s2, 0.0);
        assert!((w2 - 4.0).abs() < 1e-15);
        assert!((q_criterion(&g) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn invariants_do_not_depend_on_the_frame() {
        // rotate a general gradient by 30° about z: G' = R G Rᵀ
        let g = [[0.3, 1.2, -0.4], [0.1, -0.5, 0.7], [0.9, 0.2, 0.2]];
        let (c, s) = (0.5_f64.sqrt() * (1.5_f64).sqrt(), 0.5);
        let r = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        let mut rt = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                rt[i][j] = r[j][i];
            }
        }
        let gp = product(&product(&r, &g), &rt);
        let (a, b) = (pope(&g, 0.7), pope(&gp, 0.7));
        for i in 0..5 {
            assert!((a[i] - b[i]).abs() < 1e-14, "λ{} {} {}", i + 1, a[i], b[i]);
        }
        let (s2a, w2a) = s2_omega2(&g);
        let (s2b, w2b) = s2_omega2(&gp);
        assert!((s2a - s2b).abs() < 1e-14 && (w2a - w2b).abs() < 1e-14);
    }

    #[test]
    fn shielding_is_one_far_from_walls_and_zero_near_them() {
        let g = shear(100.0);
        assert!(shielding(&g, 1e-3, 1e-5, 10.0) > 0.999);
        assert!(shielding(&g, 1e-3, 1e-5, 1e-4) < 1e-6);
    }
}

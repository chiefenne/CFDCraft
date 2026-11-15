//! Finite-volume operators on a scalar φ (stage 2): least-squares gradient
//! (D24), Venkatakrishnan limiters (D25), face flux coefficients with linear
//! upwind convection and over-relaxed non-orthogonal diffusion (D25, D26),
//! and the gather-based row assembly (D31).
//!
//! Data are structure-of-arrays slices (one slice per vector component), so
//! the same bodies run in a CPU loop and as one GPU thread per index.
//!
//! **Boundary faces.** Every boundary face carries two coefficients `a`, `b`
//! of its face value, φ_f = a φ_P + b:
//!
//! - fixed value φ_b: a = 0, b = φ_b; the gradient stencil point is the face
//!   centroid;
//! - fixed normal gradient g (zero for symmetry and zero-gradient faces):
//!   a = 1, b = g d_n with d_n = (x_f − x_P)·n̂; the stencil point is the
//!   projection x_P + d_n n̂, so the gradient of a linear field stays exact.
//!
//! With the diffusion factor of the face, D = Γ |S|²/(d·S), the same formula
//! −D (φ_f − φ_P) gives the diffusive flux of both kinds, because
//! D g d_n = Γ g |S| for a fixed gradient.
//!
//! **Face flux.** The flux of φ out of the owner P through face f is linear
//! in the cell values, F_f = c_P φ_P + c_N φ_N + b_f. The implicit part is
//! first-order upwind convection and the orthogonal part of the diffusion;
//! `b_f` carries the deferred corrections: the linear-upwind
//! reconstruction ṁ (φ_f − φ_U), blended by `blend` (0 upwind, 1 linear
//! upwind), and the non-orthogonal diffusion −Γ k·∇φ_f (Jasak 1996, the
//! over-relaxed decomposition S = Δ + k with Δ = d |S|²/(d·S)), switched by
//! `nonorth` and limited by `nonorth_limit`.

use crate::{ld, unpack_face};
use inncraft_core::Real;

/// Three components of a vector field, one slice each.
#[derive(Clone, Copy)]
pub struct Vec3<'a, T> {
    /// x components.
    pub x: &'a [T],
    /// y components.
    pub y: &'a [T],
    /// z components.
    pub z: &'a [T],
}

impl<T: Real> Vec3<'_, T> {
    /// Vector `i`.
    ///
    /// # Safety
    /// `i` is inside all three slices.
    #[inline(always)]
    pub unsafe fn at(&self, i: usize) -> [T; 3] {
        // SAFETY: guaranteed by the caller.
        unsafe { [ld(self.x, i), ld(self.y, i), ld(self.z, i)] }
    }
}

#[inline(always)]
fn dot<T: Real>(a: [T; 3], b: [T; 3]) -> T {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline(always)]
fn sub<T: Real>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Face-based connectivity of a partition (see `inncraft_mesh::Mesh`).
#[derive(Clone, Copy)]
pub struct Topology<'a> {
    /// Owner cell of each face.
    pub owner: &'a [u32],
    /// Neighbour cell of each interior face.
    pub neighbour: &'a [u32],
    /// Start of the face list of each local cell.
    pub cell_ptr: &'a [u32],
    /// Faces of each local cell, packed with the orientation bit.
    pub cell_faces: &'a [u32],
    /// Number of interior (inner and processor) faces; boundary faces
    /// follow.
    pub n_interior: usize,
}

/// Boundary coefficients of the face values, φ_f = a φ_P + b, one entry per
/// boundary face (index `f − n_interior`).
#[derive(Clone, Copy)]
pub struct Boundary<'a, T> {
    /// Coefficient of the owner value.
    pub a: &'a [T],
    /// Constant part.
    pub b: &'a [T],
}

impl<T: Real> Boundary<'_, T> {
    /// Face value of boundary face `fb` for owner value `phi_p`.
    ///
    /// # Safety
    /// `fb` is inside both slices.
    #[inline(always)]
    pub unsafe fn value(&self, fb: usize, phi_p: T) -> T {
        // SAFETY: guaranteed by the caller.
        unsafe { ld(self.a, fb) * phi_p + ld(self.b, fb) }
    }
}

/// Weighted least-squares gradient of cell `c` (D24):
/// ∇φ_c = Σ_k w_k (φ_k − φ_c) over the cell's faces in stored order, with
/// the weight vectors `w` precomputed on the host per cell-face entry,
/// w_k = G⁻¹ ω_k d_k, G = Σ ω_k d_k d_kᵀ.
///
/// # Safety
/// `c` is a local cell of `topo`, every index stored in `topo` is inside the
/// arrays it addresses, `w` has one entry per cell-face entry, `bc` one per
/// boundary face.
#[inline(always)]
pub unsafe fn lsq_gradient<T: Real>(
    c: usize,
    topo: &Topology,
    w: &Vec3<T>,
    phi: &[T],
    bc: &Boundary<T>,
) -> [T; 3] {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let pc = ld(phi, c);
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut g = [T::ZERO; 3];
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            let d = if f < topo.n_interior {
                let other = if is_nb {
                    ld(topo.owner, f)
                } else {
                    ld(topo.neighbour, f)
                };
                ld(phi, other as usize) - pc
            } else {
                bc.value(f - topo.n_interior, pc) - pc
            };
            let wk = w.at(k);
            g[0] += wk[0] * d;
            g[1] += wk[1] * d;
            g[2] += wk[2] * d;
            k += 1;
        }
        g
    }
}

/// Venkatakrishnan's limiter function of one face (Venkatakrishnan 1995):
/// for the unlimited increment Δ₋ = ∇φ·(x_f − x_P) and the allowed
/// increment Δ₊ (to the neighbourhood maximum or minimum),
///
/// ψ = (Δ₊² + ε² + 2 Δ₋ Δ₊) / (Δ₊² + 2 Δ₋² + Δ₋ Δ₊ + ε²),
///
/// clipped to 1. With ε = 0 the limited increment ψ Δ₋ never exceeds Δ₊, so
/// no new extremum is created; ε > 0 lets smooth regions keep the unlimited
/// gradient.
#[inline(always)]
pub fn venkat_psi<T: Real>(dplus: T, dminus: T, eps2: T) -> T {
    let num = dplus * dplus + eps2 + T::from_f64(2.0) * dminus * dplus;
    let den = dplus * dplus + T::from_f64(2.0) * dminus * dminus + dminus * dplus + eps2;
    (num / den).min_of(T::ONE)
}

/// Limiter value ψ_c of cell `c`: the minimum of [`venkat_psi`] over its
/// faces, with the bounds φ_min, φ_max from the cell and its face
/// neighbours (boundary faces count with their face value).
///
/// `eps2` is ε² per cell (Venkatakrishnan: (K h_c)³ with h_c the cube root
/// of the cell volume), or a single value when `eps2.len() == 1`
/// (Venkatakrishnan–Wang: ε² = (K' (φ_max − φ_min)_domain)², Wang 2000).
///
/// # Safety
/// As [`lsq_gradient`]; `grad` and `cell_centroid` have an entry for `c`,
/// `face_centroid` one per face, `eps2` one per local cell or exactly one.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub unsafe fn venkat_limiter<T: Real>(
    c: usize,
    topo: &Topology,
    phi: &[T],
    bc: &Boundary<T>,
    grad: &Vec3<T>,
    cell_centroid: &Vec3<T>,
    face_centroid: &Vec3<T>,
    eps2: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let pc = ld(phi, c);
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let (mut mn, mut mx) = (pc, pc);
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            let v = if f < topo.n_interior {
                let other = if is_nb {
                    ld(topo.owner, f)
                } else {
                    ld(topo.neighbour, f)
                };
                ld(phi, other as usize)
            } else {
                bc.value(f - topo.n_interior, pc)
            };
            mn = mn.min_of(v);
            mx = mx.max_of(v);
            k += 1;
        }
        let e2 = if eps2.len() == 1 {
            ld(eps2, 0)
        } else {
            ld(eps2, c)
        };
        let g = grad.at(c);
        let xc = cell_centroid.at(c);
        let mut psi = T::ONE;
        k = lo;
        while k < hi {
            let (f, _) = unpack_face(ld(topo.cell_faces, k));
            let dm = dot(g, sub(face_centroid.at(f), xc));
            let p = if dm > T::ZERO {
                venkat_psi(mx - pc, dm, e2)
            } else if dm < T::ZERO {
                venkat_psi(mn - pc, dm, e2)
            } else {
                T::ONE
            };
            psi = psi.min_of(p);
            k += 1;
        }
        psi
    }
}

/// Per-face inputs of [`face_coefficients`].
#[derive(Clone, Copy)]
pub struct FaceInputs<'a, T> {
    /// Connectivity.
    pub topo: Topology<'a>,
    /// Boundary face values.
    pub bc: Boundary<'a, T>,
    /// Mass flux ṁ_f (kg/s), positive from owner to neighbour (out of the
    /// domain on a boundary face).
    pub mass_flux: &'a [T],
    /// Diffusion coefficient Γ_f on the face (kg/(m s)).
    pub gamma: &'a [T],
    /// Geometric diffusion factor |S|²/(d·S) per face (m), d from owner to
    /// neighbour or to the boundary face centroid.
    pub dgeo: &'a [T],
    /// Non-orthogonal part k = S − d |S|²/(d·S) per face (m²).
    pub k: Vec3<'a, T>,
    /// Weight of the owner in the face interpolation of the gradient.
    pub weight: &'a [T],
    /// Face centroids.
    pub face_centroid: Vec3<'a, T>,
    /// Cell centroids (local and halo cells).
    pub cell_centroid: Vec3<'a, T>,
    /// Cell values.
    pub phi: &'a [T],
    /// Cell gradients.
    pub grad: Vec3<'a, T>,
    /// Limiter value per cell (1 for no limiter).
    pub psi: &'a [T],
    /// Blending of the linear-upwind correction: 0 first-order upwind, 1
    /// linear upwind.
    pub blend: T,
    /// 1 to apply the non-orthogonal diffusion correction, 0 to omit it.
    pub nonorth: T,
    /// Limit of the non-orthogonal correction as a multiple of the
    /// orthogonal part (Jasak 1996); a negative value means no limit.
    pub nonorth_limit: T,
}

/// Coefficients `(c_P, c_N, b_f)` of the flux out of the owner through face
/// `f`, F_f = c_P φ_P + c_N φ_N + b_f (`c_N = 0` on boundary faces).
///
/// # Safety
/// `f` is a face of the topology and every array of `inp` is long enough
/// for the indices it is read at.
#[inline(always)]
pub unsafe fn face_coefficients<T: Real>(f: usize, inp: &FaceInputs<T>) -> (T, T, T) {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let t = &inp.topo;
        let p = ld(t.owner, f) as usize;
        let m = ld(inp.mass_flux, f);
        let gam = ld(inp.gamma, f);
        let dd = gam * ld(inp.dgeo, f);
        let mp = m.max_of(T::ZERO);
        let mm = m.min_of(T::ZERO);
        let xf = inp.face_centroid.at(f);
        let phi_p = ld(inp.phi, p);
        let gp = inp.grad.at(p);
        let kf = inp.k.at(f);
        if f < t.n_interior {
            let n = ld(t.neighbour, f) as usize;
            let gn = inp.grad.at(n);
            // linear-upwind deferred correction ṁ (φ_f − φ_U)
            let ho = if m >= T::ZERO {
                ld(inp.psi, p) * dot(gp, sub(xf, inp.cell_centroid.at(p)))
            } else {
                ld(inp.psi, n) * dot(gn, sub(xf, inp.cell_centroid.at(n)))
            };
            let mut b = inp.blend * m * ho;
            // non-orthogonal correction −Γ k·∇φ_f, ∇φ_f interpolated
            let w = ld(inp.weight, f);
            let wn = T::ONE - w;
            let gf = [
                w * gp[0] + wn * gn[0],
                w * gp[1] + wn * gn[1],
                w * gp[2] + wn * gn[2],
            ];
            let corr = gam * dot(kf, gf);
            let orth = dd * (ld(inp.phi, n) - phi_p);
            b -= inp.nonorth * limit_corr(corr, orth, inp.nonorth_limit);
            (mp + dd, mm - dd, b)
        } else {
            let fb = f - t.n_interior;
            let a = ld(inp.bc.a, fb);
            let bb = ld(inp.bc.b, fb);
            let one_a = T::ONE - a;
            let ho = ld(inp.psi, p) * dot(gp, sub(xf, inp.cell_centroid.at(p)));
            let corr = gam * dot(kf, gp);
            let orth = dd * (inp.bc.value(fb, phi_p) - phi_p);
            let b = mm * bb - dd * bb + inp.blend * mp * ho
                - inp.nonorth * one_a * limit_corr(corr, orth, inp.nonorth_limit);
            (mp + a * mm + one_a * dd, T::ZERO, b)
        }
    }
}

/// Limits the non-orthogonal correction `corr` to `lim` times the size of
/// the orthogonal part `orth`; `lim < 0` leaves it unchanged.
#[inline(always)]
fn limit_corr<T: Real>(corr: T, orth: T, lim: T) -> T {
    if lim < T::ZERO {
        corr
    } else {
        let bound = lim * orth.abs();
        corr.max_of(-bound).min_of(bound)
    }
}

/// Matrix storage of the assembly: sliced ELL (SELL-C-σ with σ = 1, rows in
/// cell order), slot 0 of each row is the diagonal, slots 1.. the interior
/// faces of the cell in stored order.
#[derive(Clone, Copy)]
pub struct SellLayout<'a> {
    /// Slice height C.
    pub chunk: usize,
    /// Start of each slice in `col`/`val`.
    pub slice_ptr: &'a [u32],
}

impl SellLayout<'_> {
    /// Storage position of slot `j` of row `r`.
    ///
    /// # Safety
    /// `r / chunk < slice_ptr.len()`.
    #[inline(always)]
    pub unsafe fn pos(&self, r: usize, j: usize) -> usize {
        let s = r / self.chunk;
        // SAFETY: guaranteed by the caller.
        unsafe { ld(self.slice_ptr, s) as usize + j * self.chunk + (r - s * self.chunk) }
    }
}

/// Assembly of row `c` by a gather over the cell's faces (D31): the row of
/// the conservation equation Σ_f ±F_f = S_c (out of the cell positive),
/// with implicit under-relaxation α (Patankar 1980):
/// a_P/α φ_P + Σ a_N φ_N = rhs + (1 − α)/α a_P φ_P.
///
/// Writes the diagonal and the off-diagonals of the row into `val` (SELL
/// slots), the right-hand side into `rhs[c]` and 1/a_P(relaxed) into
/// `dinv[c]`. Faces are summed in the stored order, so the result is the
/// same for any loop order over the rows.
///
/// # Safety
/// As [`lsq_gradient`]; `cp`, `cn`, `bf` have one entry per face; `val`
/// holds the SELL storage of `layout` with room for every interior face of
/// the cell; `rhs`, `dinv`, `source`, `phi` have an entry for `c`. No other
/// thread writes row `c`.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub unsafe fn assemble_row<T: Real>(
    c: usize,
    topo: &Topology,
    cp: &[T],
    cn: &[T],
    bf: &[T],
    source: &[T],
    phi: &[T],
    alpha: T,
    layout: &SellLayout,
    val: *mut T,
    rhs: *mut T,
    dinv: *mut T,
) {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut diag = T::ZERO;
        let mut b = ld(source, c);
        let mut j = 1;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            if is_nb {
                diag -= ld(cn, f);
                b += ld(bf, f);
                *val.add(layout.pos(c, j)) = -ld(cp, f);
                j += 1;
            } else {
                diag += ld(cp, f);
                b -= ld(bf, f);
                if f < topo.n_interior {
                    *val.add(layout.pos(c, j)) = ld(cn, f);
                    j += 1;
                }
            }
            k += 1;
        }
        let dr = diag / alpha;
        b += (dr - diag) * ld(phi, c);
        *val.add(layout.pos(c, 0)) = dr;
        *rhs.add(c) = b;
        *dinv.add(c) = T::ONE / dr;
    }
}

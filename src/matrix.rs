// Ported from PaulStoffregen/MotionCal matrix.c.
// SPDX-License-Identifier: BSD-3-Clause
//
// Copyright (c) 2014, Freescale Semiconductor, Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Freescale Semiconductor, Inc. nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL FREESCALE SEMICONDUCTOR, INC. BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//! Direct ports of the relevant routines in [`matrix.c`][upstream] from
//! PaulStoffregen/MotionCal.
//!
//! [upstream]: https://github.com/PaulStoffregen/MotionCal/blob/master/matrix.c

use libm::{fabsf, sqrtf};

/// `f3x3matrixAeqI` in `matrix.c`: set A to the 3x3 identity.
pub(crate) fn f3x3_matrix_a_eq_i(a: &mut [[f32; 3]; 3]) {
    for r in 0..3 {
        for c in 0..3 {
            a[r][c] = 0.0;
        }
        a[r][r] = 1.0;
    }
}

/// `f3x3matrixAeqAxScalar` in `matrix.c`: scale every entry of A by k.
pub(crate) fn f3x3_matrix_a_eq_a_x_scalar(a: &mut [[f32; 3]; 3], k: f32) {
    for r in 0..3 {
        for c in 0..3 {
            a[r][c] *= k;
        }
    }
}

/// `f3x3matrixAeqMinusA` in `matrix.c`: negate A in place.
pub(crate) fn f3x3_matrix_a_eq_minus_a(a: &mut [[f32; 3]; 3]) {
    for r in 0..3 {
        for c in 0..3 {
            a[r][c] = -a[r][c];
        }
    }
}

/// `f3x3matrixDetA` in `matrix.c`: 3x3 determinant.
pub(crate) fn f3x3_matrix_det_a(a: &[[f32; 3]; 3]) -> f32 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        + a[0][1] * (a[1][2] * a[2][0] - a[1][0] * a[2][2])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// `f3x3matrixAeqInvSymB` in `matrix.c`: A = inverse(B) for symmetric
/// B (only on/above-diagonal entries of B are read). On singular B, A
/// is set to the identity (matches the C).
pub(crate) fn f3x3_matrix_a_eq_inv_sym_b(a: &mut [[f32; 3]; 3], b: &[[f32; 3]; 3]) {
    let f_b11_b22_m_b12_b12 = b[1][1] * b[2][2] - b[1][2] * b[1][2];
    let f_b12_b02_m_b01_b22 = b[1][2] * b[0][2] - b[0][1] * b[2][2];
    let f_b01_b12_m_b11_b02 = b[0][1] * b[1][2] - b[1][1] * b[0][2];

    let mut ftmp = b[0][0] * f_b11_b22_m_b12_b12
        + b[0][1] * f_b12_b02_m_b01_b22
        + b[0][2] * f_b01_b12_m_b11_b02;

    if ftmp != 0.0 {
        ftmp = 1.0 / ftmp;
        a[0][0] = f_b11_b22_m_b12_b12 * ftmp;
        a[1][0] = f_b12_b02_m_b01_b22 * ftmp;
        a[0][1] = a[1][0];
        a[2][0] = f_b01_b12_m_b11_b02 * ftmp;
        a[0][2] = a[2][0];
        a[1][1] = (b[0][0] * b[2][2] - b[0][2] * b[0][2]) * ftmp;
        a[2][1] = (b[0][2] * b[0][1] - b[0][0] * b[1][2]) * ftmp;
        a[1][2] = a[2][1];
        a[2][2] = (b[0][0] * b[1][1] - b[0][1] * b[0][1]) * ftmp;
    } else {
        f3x3_matrix_a_eq_i(a);
    }
}

/// `fmatrixAeqInvA` in `matrix.c`, specialised to 4x4: in-place
/// Gauss-Jordan inverse with full pivoting. On singular A, A is set
/// to the identity (matches the C).
pub(crate) fn fmatrix_inverse_4x4(a: &mut [[f32; 4]; 4]) {
    const N: usize = 4;
    let mut pivot_row = 0usize;
    let mut pivot_col = 0usize;
    let mut col_ind = [0usize; N];
    let mut row_ind = [0usize; N];
    let mut pivots = [0u8; N];

    for i in 0..N {
        let mut largest: f32 = 0.0;
        for j in 0..N {
            if pivots[j] == 1 {
                continue;
            }
            for k in 0..N {
                if pivots[k] == 0 {
                    let abs = fabsf(a[j][k]);
                    if abs >= largest {
                        pivot_row = j;
                        pivot_col = k;
                        largest = abs;
                    }
                } else if pivots[k] > 1 {
                    // singular
                    for r in 0..N {
                        for c in 0..N {
                            a[r][c] = if r == c { 1.0 } else { 0.0 };
                        }
                    }
                    return;
                }
            }
        }
        pivots[pivot_col] += 1;

        if pivot_row != pivot_col {
            for l in 0..N {
                let t = a[pivot_row][l];
                a[pivot_row][l] = a[pivot_col][l];
                a[pivot_col][l] = t;
            }
        }
        row_ind[i] = pivot_row;
        col_ind[i] = pivot_col;

        if a[pivot_col][pivot_col] == 0.0 {
            for r in 0..N {
                for c in 0..N {
                    a[r][c] = if r == c { 1.0 } else { 0.0 };
                }
            }
            return;
        }

        let recip = 1.0 / a[pivot_col][pivot_col];
        a[pivot_col][pivot_col] = 1.0;
        for l in 0..N {
            a[pivot_col][l] *= recip;
        }
        for m in 0..N {
            if m == pivot_col {
                continue;
            }
            let scale = a[m][pivot_col];
            a[m][pivot_col] = 0.0;
            for l in 0..N {
                a[m][l] -= a[pivot_col][l] * scale;
            }
        }
    }

    for l in (0..N).rev() {
        let i = row_ind[l];
        let j = col_ind[l];
        if i != j {
            for k in 0..N {
                let t = a[k][i];
                a[k][i] = a[k][j];
                a[k][j] = t;
            }
        }
    }
}

/// `eigencompute` in `matrix.c`: in-place Jacobi-rotation eigenvalue
/// solver for a symmetric real matrix stored in the top-left `n x n`
/// block of a 10x10 buffer. On exit `eigval[0..n]` holds the
/// eigenvalues and `eigvec[i][j]` (for `i, j < n`) holds the j-th
/// normalised eigenvector down rows. `mat_a` is clobbered.
pub(crate) fn eigencompute(
    mat_a: &mut [[f32; 10]; 10],
    eigval: &mut [f32; 10],
    eigvec: &mut [[f32; 10]; 10],
    n: usize,
) {
    const N_ITERATIONS: usize = 15;

    for ir in 0..n {
        for ic in 0..n {
            eigvec[ir][ic] = 0.0;
        }
        eigvec[ir][ir] = 1.0;
        eigval[ir] = mat_a[ir][ir];
    }

    let mut ctr = 0;
    loop {
        let mut residue: f32 = 0.0;
        for ir in 0..n.saturating_sub(1) {
            for ic in (ir + 1)..n {
                residue += fabsf(mat_a[ir][ic]);
            }
        }
        if residue <= 0.0 {
            break;
        }

        for ir in 0..n.saturating_sub(1) {
            for ic in (ir + 1)..n {
                if fabsf(mat_a[ir][ic]) <= 0.0 {
                    continue;
                }
                let cot2phi = 0.5 * (eigval[ic] - eigval[ir]) / mat_a[ir][ic];
                let mut tanphi = 1.0 / (fabsf(cot2phi) + sqrtf(1.0 + cot2phi * cot2phi));
                if cot2phi < 0.0 {
                    tanphi = -tanphi;
                }
                let cosphi = 1.0 / sqrtf(1.0 + tanphi * tanphi);
                let sinphi = tanphi * cosphi;
                let tanhalfphi = sinphi / (1.0 + cosphi);

                let ftmp = tanphi * mat_a[ir][ic];
                eigval[ir] -= ftmp;
                eigval[ic] += ftmp;
                mat_a[ir][ic] = 0.0;

                for j in 0..n {
                    let t = eigvec[j][ir];
                    eigvec[j][ir] = t - sinphi * (eigvec[j][ic] + tanhalfphi * t);
                    eigvec[j][ic] = eigvec[j][ic] + sinphi * (t - tanhalfphi * eigvec[j][ic]);
                }
                if ir > 0 {
                    for j in 0..ir {
                        let t = mat_a[j][ir];
                        mat_a[j][ir] = t - sinphi * (mat_a[j][ic] + tanhalfphi * t);
                        mat_a[j][ic] = mat_a[j][ic] + sinphi * (t - tanhalfphi * mat_a[j][ic]);
                    }
                }
                if ic > ir + 1 {
                    for j in (ir + 1)..ic {
                        let t = mat_a[ir][j];
                        mat_a[ir][j] = t - sinphi * (mat_a[j][ic] + tanhalfphi * t);
                        mat_a[j][ic] = mat_a[j][ic] + sinphi * (t - tanhalfphi * mat_a[j][ic]);
                    }
                }
                for j in (ic + 1)..n {
                    let t = mat_a[ir][j];
                    mat_a[ir][j] = t - sinphi * (mat_a[ic][j] + tanhalfphi * t);
                    mat_a[ic][j] = mat_a[ic][j] + sinphi * (t - tanhalfphi * mat_a[ic][j]);
                }
            }
        }

        ctr += 1;
        if ctr >= N_ITERATIONS {
            break;
        }
    }
}

// Ported from PaulStoffregen/MotionCal magcal.c.
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

#![allow(clippy::needless_range_loop)]

use libm::{fabsf, powf, sqrtf};

use crate::matrix::{
    eigencompute, f3x3_matrix_a_eq_a_x_scalar, f3x3_matrix_a_eq_i, f3x3_matrix_a_eq_inv_sym_b,
    f3x3_matrix_a_eq_minus_a, f3x3_matrix_det_a, fmatrix_inverse_4x4,
};
use crate::{FitError, MagCal, SolverTier};

const ONE_THIRD: f32 = 0.333_333_33;
const ONE_SIXTH: f32 = 0.166_666_67;

/// Numerical-conditioning scale. Applied to inputs and reversed on outputs;
/// the choice doesn't change the result, only how `f32` rounding accumulates.
pub const DEFAULT_SCALE: f32 = 500.0;

pub const MIN_SAMPLES_SPHERE: usize = 40;
pub const MIN_SAMPLES_AXIS_ALIGNED: usize = 100;
pub const MIN_SAMPLES_ELLIPSOID: usize = 150;

const SPHERE_PARAMS: usize = 4;
const AXIS_ALIGNED_PARAMS: usize = 7;
const ELLIPSOID_PARAMS: usize = 10;

/// Streaming magnetometer-calibration solver. Accumulates samples via
/// [`push_sample`](Self::push_sample) and, on demand, runs one of
/// the three fitters ported from [`magcal.c`][upstream]:
/// [`solve_sphere`](Self::solve_sphere) (`fUpdateCalibration4INV`),
/// [`solve_axis_aligned`](Self::solve_axis_aligned)
/// (`fUpdateCalibration7EIG`), or
/// [`solve_ellipsoid`](Self::solve_ellipsoid)
/// (`fUpdateCalibration10EIG`). [`solve`](Self::solve) picks the
/// highest tier whose recommended sample count is met.
///
/// [upstream]: https://github.com/PaulStoffregen/MotionCal/blob/master/magcal.c
#[derive(Clone)]
pub struct Solver {
    scale: f32,
    sample_count: u32,
    i_offset: [i16; 3],

    sph_mat_a: [[f32; 4]; 4],
    sph_vec_b: [f32; 4],
    sph_sum_bp4: f32,

    axis_mat_a: [[f32; 7]; 7],

    ell_mat_a: [[f32; 10]; 10],

    tr_v: [f32; 3],
    tr_inv_w: [[f32; 3]; 3],
    tr_b: f32,
    tr_fit_error_percent: f32,
    last_tier: Option<SolverTier>,
}

impl Solver {
    /// New solver with the default scale ([`DEFAULT_SCALE`]).
    pub const fn new() -> Self {
        Self::with_scale(DEFAULT_SCALE)
    }

    /// New solver with a caller-chosen numerical conditioning scale.
    pub const fn with_scale(scale: f32) -> Self {
        Self {
            scale,
            sample_count: 0,
            i_offset: [0; 3],
            sph_mat_a: [[0.0; 4]; 4],
            sph_vec_b: [0.0; 4],
            sph_sum_bp4: 0.0,
            axis_mat_a: [[0.0; 7]; 7],
            ell_mat_a: [[0.0; 10]; 10],
            tr_v: [0.0; 3],
            tr_inv_w: [[0.0; 3]; 3],
            tr_b: 0.0,
            tr_fit_error_percent: 0.0,
            last_tier: None,
        }
    }

    /// Zero the accumulators and clear the last fit, so a long-lived
    /// (e.g. `static`) [`Solver`] can be reused across cal cycles.
    pub fn reset(&mut self) {
        self.sample_count = 0;
        self.i_offset = [0; 3];
        self.sph_mat_a = [[0.0; 4]; 4];
        self.sph_vec_b = [0.0; 4];
        self.sph_sum_bp4 = 0.0;
        self.axis_mat_a = [[0.0; 7]; 7];
        self.ell_mat_a = [[0.0; 10]; 10];
        self.last_tier = None;
    }

    /// Fold one raw sample into the running per-tier accumulators.
    /// `raw` is in whatever input unit the caller uses; see the
    /// crate-level [Units](crate#units) note.
    pub fn push_sample(&mut self, raw: [i16; 3]) {
        if self.sample_count == 0 {
            self.i_offset = raw;
        }
        let fscaling = 1.0 / self.scale;
        let sx = ((raw[0] as i32 - self.i_offset[0] as i32) as f32) * fscaling;
        let sy = ((raw[1] as i32 - self.i_offset[1] as i32) as f32) * fscaling;
        let sz = ((raw[2] as i32 - self.i_offset[2] as i32) as f32) * fscaling;

        let bp2 = sx * sx + sy * sy + sz * sz;
        self.sph_sum_bp4 += bp2 * bp2;
        self.sph_vec_b[0] += sx * bp2;
        self.sph_vec_b[1] += sy * bp2;
        self.sph_vec_b[2] += sz * bp2;
        self.sph_vec_b[3] += bp2;
        self.sph_mat_a[0][0] += sx * sx;
        self.sph_mat_a[0][1] += sx * sy;
        self.sph_mat_a[0][2] += sx * sz;
        self.sph_mat_a[0][3] += sx;
        self.sph_mat_a[1][1] += sy * sy;
        self.sph_mat_a[1][2] += sy * sz;
        self.sph_mat_a[1][3] += sy;
        self.sph_mat_a[2][2] += sz * sz;
        self.sph_mat_a[2][3] += sz;

        let av = [sx * sx, sy * sy, sz * sz, sx, sy, sz];
        for m in 0..6 {
            self.axis_mat_a[m][6] += av[m];
            for n in m..6 {
                self.axis_mat_a[m][n] += av[m] * av[n];
            }
        }

        let ev = [
            sx * sx,
            2.0 * sx * sy,
            2.0 * sx * sz,
            sy * sy,
            2.0 * sy * sz,
            sz * sz,
            sx,
            sy,
            sz,
        ];
        for m in 0..9 {
            self.ell_mat_a[m][9] += ev[m];
            for n in m..9 {
                self.ell_mat_a[m][n] += ev[m] * ev[n];
            }
        }

        self.sample_count += 1;
    }

    /// Number of samples pushed since construction / last reset.
    pub fn sample_count(&self) -> usize {
        self.sample_count as usize
    }

    /// The most recent fit's result, or `None` if no fit has been run
    /// since the last [`reset`](Self::reset).
    pub fn last_fit(&self) -> Option<MagCal> {
        self.last_tier.map(|tier| MagCal {
            hard_iron: self.tr_v,
            soft_iron: self.tr_inv_w,
            field_strength: self.tr_b,
            fit_error_percent: self.tr_fit_error_percent,
            tier,
        })
    }

    /// Run a fit using the highest-tier solver whose recommended
    /// sample count (see [`SolverTier::min_samples`]) is met by the
    /// current sample count, falling back to a sphere fit otherwise.
    ///
    /// Returns [`FitError::NotEnoughSamples`] only when even the
    /// sphere fit cannot run (fewer than 4 samples).
    ///
    /// For more control, call the direct `solve_*` methods or use
    /// [`MagCal::fit_with`] instead.
    pub fn solve(&mut self) -> Result<MagCal, FitError> {
        let n = self.sample_count();
        if n >= MIN_SAMPLES_ELLIPSOID {
            self.solve_ellipsoid()
        } else if n >= MIN_SAMPLES_AXIS_ALIGNED {
            self.solve_axis_aligned()
        } else {
            self.solve_sphere()
        }
    }

    /// Run a full-ellipsoid fit (hard iron + 3x3 soft iron). Needs at least 10 samples.
    ///
    /// This corresponds with `fUpdateCalibration10EIG` in `magcal.c`.
    pub fn solve_ellipsoid(&mut self) -> Result<MagCal, FitError> {
        let count = self.sample_count();
        if count < ELLIPSOID_PARAMS {
            return Err(FitError::NotEnoughSamples {
                provided: count,
                required: ELLIPSOID_PARAMS,
            });
        }

        let mut mat_a: [[f32; 10]; 10] = self.ell_mat_a;
        mat_a[9][9] = count as f32;
        for m in 1..10 {
            for n in 0..m {
                mat_a[m][n] = mat_a[n][m];
            }
        }

        let mut mat_b = [[0.0_f32; 10]; 10];
        let mut vec_a = [0.0_f32; 10];
        eigencompute(&mut mat_a, &mut vec_a, &mut mat_b, 10);

        let mut j = 0usize;
        for i in 1..10 {
            if vec_a[i] < vec_a[j] {
                j = i;
            }
        }

        let mut a = [[0.0_f32; 3]; 3];
        a[0][0] = mat_b[0][j];
        a[0][1] = mat_b[1][j];
        a[1][0] = a[0][1];
        a[0][2] = mat_b[2][j];
        a[2][0] = a[0][2];
        a[1][1] = mat_b[3][j];
        a[1][2] = mat_b[4][j];
        a[2][1] = a[1][2];
        a[2][2] = mat_b[5][j];

        let mut det = f3x3_matrix_det_a(&a);
        if det < 0.0 {
            f3x3_matrix_a_eq_minus_a(&mut a);
            mat_b[6][j] = -mat_b[6][j];
            mat_b[7][j] = -mat_b[7][j];
            mat_b[8][j] = -mat_b[8][j];
            mat_b[9][j] = -mat_b[9][j];
            det = -det;
        }

        let mut inv_a = [[0.0_f32; 3]; 3];
        f3x3_matrix_a_eq_inv_sym_b(&mut inv_a, &a);

        let mut tr_v = [0.0_f32; 3];
        for k in 0..3 {
            let mut acc = 0.0_f32;
            for m in 0..3 {
                acc += inv_a[k][m] * mat_b[m + 6][j];
            }
            tr_v[k] = acc * -0.5;
        }

        let mut tr_b = sqrtf(fabsf(
            a[0][0] * tr_v[0] * tr_v[0]
                + 2.0 * a[0][1] * tr_v[0] * tr_v[1]
                + 2.0 * a[0][2] * tr_v[0] * tr_v[2]
                + a[1][1] * tr_v[1] * tr_v[1]
                + 2.0 * a[1][2] * tr_v[1] * tr_v[2]
                + a[2][2] * tr_v[2] * tr_v[2]
                - mat_b[9][j],
        ));

        let fit_err = 50.0 * sqrtf(fabsf(vec_a[j]) / count as f32) / (tr_b * tr_b);

        for k in 0..3 {
            tr_v[k] = tr_v[k] * self.scale + self.i_offset[k] as f32;
        }
        tr_b *= self.scale;

        f3x3_matrix_a_eq_a_x_scalar(&mut a, powf(det, -ONE_THIRD));
        tr_b *= powf(det, -ONE_SIXTH);

        let mut a10 = [[0.0_f32; 10]; 10];
        for i in 0..3 {
            for jj in 0..3 {
                a10[i][jj] = a[i][jj];
            }
        }
        eigencompute(&mut a10, &mut vec_a, &mut mat_b, 3);

        for jj in 0..3 {
            let ftmp = sqrtf(sqrtf(fabsf(vec_a[jj])));
            for i in 0..3 {
                mat_b[i][jj] *= ftmp;
            }
        }
        let mut tr_inv_w = [[0.0_f32; 3]; 3];
        for i in 0..3 {
            for jj in i..3 {
                let mut acc = 0.0_f32;
                for k in 0..3 {
                    acc += mat_b[i][k] * mat_b[jj][k];
                }
                tr_inv_w[i][jj] = acc;
                tr_inv_w[jj][i] = acc;
            }
        }

        self.tr_v = tr_v;
        self.tr_inv_w = tr_inv_w;
        self.tr_b = tr_b;
        self.tr_fit_error_percent = fit_err;
        self.last_tier = Some(SolverTier::Ellipsoid);
        Ok(self.last_fit().unwrap())
    }

    /// Run an axis-aligned ellipsoid fit (hard iron + diagonal soft iron). Needs at least 7 samples.
    ///
    /// This corresponds with `fUpdateCalibration7EIG` in `magcal.c`.
    pub fn solve_axis_aligned(&mut self) -> Result<MagCal, FitError> {
        let count = self.sample_count();
        if count < AXIS_ALIGNED_PARAMS {
            return Err(FitError::NotEnoughSamples {
                provided: count,
                required: AXIS_ALIGNED_PARAMS,
            });
        }

        let mut mat_a = [[0.0_f32; 10]; 10];
        for m in 0..7 {
            for n in m..7 {
                mat_a[m][n] = self.axis_mat_a[m][n];
            }
        }
        mat_a[6][6] = count as f32;
        for m in 1..7 {
            for n in 0..m {
                mat_a[m][n] = mat_a[n][m];
            }
        }

        let mut mat_b = [[0.0_f32; 10]; 10];
        let mut vec_a = [0.0_f32; 10];
        eigencompute(&mut mat_a, &mut vec_a, &mut mat_b, 7);

        let mut j = 0usize;
        for i in 1..7 {
            if vec_a[i] < vec_a[j] {
                j = i;
            }
        }

        let mut a = [[0.0_f32; 3]; 3];
        let mut det: f32 = 1.0;
        let mut tr_v = [0.0_f32; 3];
        for k in 0..3 {
            a[k][k] = mat_b[k][j];
            det *= a[k][k];
            tr_v[k] = -0.5 * mat_b[k + 3][j] / a[k][k];
        }
        if det < 0.0 {
            f3x3_matrix_a_eq_minus_a(&mut a);
            mat_b[6][j] = -mat_b[6][j];
            det = -det;
        }

        let mut ftmp = -mat_b[6][j];
        for k in 0..3 {
            ftmp += a[k][k] * tr_v[k] * tr_v[k];
        }

        let fit_err = 50.0 * sqrtf(fabsf(vec_a[j]) / count as f32) / fabsf(ftmp);

        f3x3_matrix_a_eq_a_x_scalar(&mut a, powf(det, -ONE_THIRD));
        let tr_b = sqrtf(fabsf(ftmp)) * self.scale * powf(det, -ONE_SIXTH);

        let mut tr_inv_w = [[0.0_f32; 3]; 3];
        f3x3_matrix_a_eq_i(&mut tr_inv_w);
        for k in 0..3 {
            tr_inv_w[k][k] = sqrtf(fabsf(a[k][k]));
            tr_v[k] = tr_v[k] * self.scale + self.i_offset[k] as f32;
        }

        self.tr_v = tr_v;
        self.tr_inv_w = tr_inv_w;
        self.tr_b = tr_b;
        self.tr_fit_error_percent = fit_err;
        self.last_tier = Some(SolverTier::AxisAligned);
        Ok(self.last_fit().unwrap())
    }

    /// Run a sphere fit (hard iron only; soft iron returned as the identity). Needs at least 4 samples.
    ///
    /// This corresponds with `fUpdateCalibration4INV` in `magcal.c`.
    pub fn solve_sphere(&mut self) -> Result<MagCal, FitError> {
        let count = self.sample_count();
        if count < SPHERE_PARAMS {
            return Err(FitError::NotEnoughSamples {
                provided: count,
                required: SPHERE_PARAMS,
            });
        }

        let mut mat_a4: [[f32; 4]; 4] = self.sph_mat_a;
        mat_a4[3][3] = count as f32;
        for i in 0..4 {
            for j in i + 1..4 {
                mat_a4[j][i] = mat_a4[i][j];
            }
        }
        let mat_a4_full = mat_a4;

        let mut mat_b4 = mat_a4;
        fmatrix_inverse_4x4(&mut mat_b4);

        let mut beta = [0.0_f32; 4];
        for i in 0..4 {
            let mut acc = 0.0_f32;
            for k in 0..4 {
                acc += mat_b4[i][k] * self.sph_vec_b[k];
            }
            beta[i] = acc;
        }

        let mut e = 0.0_f32;
        for i in 0..4 {
            e += beta[i] * self.sph_vec_b[i];
        }
        e = self.sph_sum_bp4 - 2.0 * e;

        let mut xtxb = [0.0_f32; 4];
        for i in 0..4 {
            let mut acc = 0.0_f32;
            for k in 0..4 {
                acc += mat_a4_full[i][k] * beta[k];
            }
            xtxb[i] = acc;
        }
        for i in 0..4 {
            e += xtxb[i] * beta[i];
        }

        let mut tr_v = [0.0_f32; 3];
        for k in 0..3 {
            tr_v[k] = 0.5 * beta[k];
        }

        let mut tr_b = sqrtf(beta[3] + tr_v[0] * tr_v[0] + tr_v[1] * tr_v[1] + tr_v[2] * tr_v[2]);
        let fit_err = sqrtf(e / count as f32) * 100.0 / (2.0 * tr_b * tr_b);

        for k in 0..3 {
            tr_v[k] = tr_v[k] * self.scale + self.i_offset[k] as f32;
        }
        tr_b *= self.scale;

        let mut tr_inv_w = [[0.0_f32; 3]; 3];
        f3x3_matrix_a_eq_i(&mut tr_inv_w);

        self.tr_v = tr_v;
        self.tr_inv_w = tr_inv_w;
        self.tr_b = tr_b;
        self.tr_fit_error_percent = fit_err;
        self.last_tier = Some(SolverTier::Sphere);
        Ok(self.last_fit().unwrap())
    }
}

impl Default for Solver {
    fn default() -> Self {
        Self::new()
    }
}

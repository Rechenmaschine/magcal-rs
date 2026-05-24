//! Magnetometer hard-iron / soft-iron calibration: a Rust port of
//! [`magcal.c`][upstream-magcal] from [PaulStoffregen/MotionCal][upstream]
//! (originally Freescale; algorithm described in NXP AN4246).
//!
//! [upstream]: https://github.com/PaulStoffregen/MotionCal
//! [upstream-magcal]: https://github.com/PaulStoffregen/MotionCal/blob/master/magcal.c
//!
//! Given raw magnetometer samples collected while the sensor was
//! rotated through varied orientations, the crate solves for the
//! offset and transform that map the corrupted point cloud back to a
//! sphere of the local field strength:
//!
//! ```text
//! corrected = soft_iron * (raw - hard_iron)
//! ```
//!
//! Three solver tiers are available, with recommended sample counts
//! (and the absolute minimum the math allows in parentheses):
//!
//! - [`SolverTier::Sphere`]: hard iron only; recommend at least 40 (absolute minimum: 4).
//! - [`SolverTier::AxisAligned`]: plus a diagonal soft iron; recommend at least 100 (absolute minimum: 7).
//! - [`SolverTier::Ellipsoid`]: plus a full 3x3 soft iron; recommend at least 150 (absolute minimum: 10).
//!
//! Recommended counts are not enforced: [`MagCal::fit`] still
//! produces a (lower-quality) sphere fit when the count falls
//! between the absolute floor and the recommended threshold.
//!
//! Sample diversity matters more than raw count. It is recommended to slowly rotate
//! through as many distinct orientations as practical.
//!
//! # Quick start
//!
//! ```no_run
//! # fn main() -> Result<(), magcal::FitError> {
//! use magcal::MagCal;
//!
//! let samples: &[[i16; 3]] = /* raw mag readings */
//! # &[];
//! let cal = MagCal::fit(samples)?;
//! let corrected = cal.apply(samples[0]);
//! # Ok(()) }
//! ```
//!
//! # Units
//!
//! Samples are `i16` in whatever input unit the caller uses (raw chip
//! counts, milligauss, deciuT). [`MagCal::hard_iron`] and
//! [`MagCal::field_strength`] come back in those same units;
//! [`MagCal::soft_iron`] is dimensionless. Multiply by your chip's
//! per-count scale at the boundary to convert to uT.
//!
//! # Embedded use
//!
//! This crate is `#![no_std]` and depends only on `libm`. The same [`Solver`] can
//! be reused across multiple calibrations to save on memory.

#![no_std]
#![deny(rustdoc::broken_intra_doc_links)]

mod matrix;
mod solver;

#[inline]
pub(crate) fn c_fabs(x: f32) -> f64 {
    // Upstream calls C `fabs(float_expr)`, which promotes through double.
    (x as f64).abs()
}

pub use solver::{
    Solver, DEFAULT_SCALE, MIN_SAMPLES_AXIS_ALIGNED, MIN_SAMPLES_ELLIPSOID, MIN_SAMPLES_SPHERE,
};

/// Which solver tier produced a [`MagCal`]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SolverTier {
    /// Hard iron only; soft iron returned as the identity.
    Sphere,
    /// Hard iron + diagonal soft iron.
    AxisAligned,
    /// Hard iron + full 3x3 soft iron.
    Ellipsoid,
}

impl SolverTier {
    /// Recommended minimum sample count for this tier.
    pub const fn min_samples(self) -> usize {
        match self {
            SolverTier::Sphere => MIN_SAMPLES_SPHERE,
            SolverTier::AxisAligned => MIN_SAMPLES_AXIS_ALIGNED,
            SolverTier::Ellipsoid => MIN_SAMPLES_ELLIPSOID,
        }
    }
}

/// Result of a successful fit. Apply via [`apply`](Self::apply) or
/// by hand as `corrected = soft_iron * (raw - hard_iron)`. Units
/// follow the input samples; see the crate-level [Units](crate#units)
/// note.
#[derive(Debug, Clone, Copy)]
pub struct MagCal {
    /// Hard-iron offset, in the same units as the input samples.
    pub hard_iron: [f32; 3],
    /// Soft-iron transform; dimensionless.
    pub soft_iron: [[f32; 3]; 3],
    /// Estimated local field strength, in the same units as the input samples.
    pub field_strength: f32,
    /// Fit residual as a percentage of [`field_strength`](Self::field_strength).
    pub fit_error_percent: f32,
    /// Tier that produced this fit.
    pub tier: SolverTier,
}

impl MagCal {
    /// Solve for a calibration, auto-picking the tier from the
    /// sample count (see [`SolverTier::min_samples`]).
    ///
    /// Coverage of the orientation sphere matters more than raw
    /// count. Below the recommended Sphere threshold, a sphere fit
    /// is still produced but quality suffers.
    ///
    /// # Errors
    ///
    /// Returns [`FitError::NotEnoughSamples`] if there are fewer than
    /// 4 samples (the sphere fit's absolute minimum requirement).
    pub fn fit(samples: &[[i16; 3]]) -> Result<Self, FitError> {
        let mut cal = Solver::new();
        for s in samples {
            cal.push_sample(*s);
        }
        cal.solve()
    }

    /// Solve for a calibration using the given solver tier.
    ///
    /// Like [`fit`](Self::fit), but uses the supplied `tier` instead
    /// of choosing from the sample count, and requires only as many
    /// samples as the tier has parameters (4 / 7 / 10). Useful for
    /// running a higher-tier solver on a thinner dataset than
    /// [`fit`](Self::fit) would accept, at the cost of fit quality.
    ///
    /// # Errors
    ///
    /// Returns [`FitError::NotEnoughSamples`] if there are fewer samples than the
    /// tier's absolute minimum requirement (4 for Sphere, 7 for AxisAligned, 10 for Ellipsoid).
    pub fn fit_with(samples: &[[i16; 3]], tier: SolverTier) -> Result<Self, FitError> {
        let mut cal = Solver::new();
        for s in samples {
            cal.push_sample(*s);
        }
        match tier {
            SolverTier::Sphere => cal.solve_sphere(),
            SolverTier::AxisAligned => cal.solve_axis_aligned(),
            SolverTier::Ellipsoid => cal.solve_ellipsoid(),
        }
    }

    /// Apply this calibration to a raw sample.
    ///
    /// `raw` must be in the same units the calibration was fitted in!
    /// After calibration, `|result|` should be approximately [`field_strength`](Self::field_strength).
    pub fn apply(&self, raw: [i16; 3]) -> [f32; 3] {
        let r = [
            raw[0] as f32 - self.hard_iron[0],
            raw[1] as f32 - self.hard_iron[1],
            raw[2] as f32 - self.hard_iron[2],
        ];
        [
            self.soft_iron[0][0] * r[0] + self.soft_iron[0][1] * r[1] + self.soft_iron[0][2] * r[2],
            self.soft_iron[1][0] * r[0] + self.soft_iron[1][1] * r[1] + self.soft_iron[1][2] * r[2],
            self.soft_iron[2][0] * r[0] + self.soft_iron[2][1] * r[1] + self.soft_iron[2][2] * r[2],
        ]
    }
}

#[derive(Debug, Clone, Copy)]
pub enum FitError {
    /// Sample count is below the tier's absolute minimum (4 / 7 / 10).
    NotEnoughSamples { provided: usize, required: usize },
}

impl core::fmt::Display for FitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            FitError::NotEnoughSamples { provided, required } => write!(
                f,
                "not enough samples: provided {provided}, required at least {required}",
            ),
        }
    }
}

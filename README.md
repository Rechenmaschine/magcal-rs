<h1 align="center">magcal-rs</h1>

<p align="center">
  <em>Magnetometer calibration in pure Rust, optimized for embedded use.</em>
  <br>
  <strong>Minimal dependencies.</strong>
  <strong>Dead-simple API surface.</strong>
  <strong>Very much <code>no_std</code>.</strong>
</p>

This is a Rust port of Freescale/NXP [`magcal.c`](https://github.com/PaulStoffregen/MotionCal/blob/master/magcal.c)
by Paul Stoffregen.

Fits hard-iron and soft-iron magnetometer calibration from a batch of raw samples:

```text
corrected = soft_iron * (raw - hard_iron)
```

There are three solvers available:

* `sphere`: fits hard iron only.
* `axis_aligned`: fits hard iron and diagonal soft iron.
* `ellipsoid`: fits hard iron and full soft iron.

Each tier requires more samples, compute, and memory, but gives increasingly better results.

## Usage

```rust
use magcal::{MagCal, Solver, SolverTier};

let samples: &[[i16; 3]] = todo!();

// Pick the best tier supported by the sample count.
let cal = MagCal::fit(samples).expect("fit failed");

// Or pick a specific tier yourself.
let cal = MagCal::fit_with(samples, SolverTier::Ellipsoid).expect("fit failed");

// If low on memory, reuse a long-lived Solver, for example a static one.
let mut s = Solver::new();

s.push_sample([100, 200, -150]);
// ...

let cal = s.solve().expect("fit failed");
// Alternatively: solve_sphere / solve_axis_aligned / solve_ellipsoid

s.reset(); // Clear samples and reuse.

// Apply to a fresh raw sample.
let corrected = cal.apply([100, 200, -150]);
```

## Contributing

PRs welcome.

## Compatibility

MSRV: **Rust 1.63**. `#![no_std]`, no `alloc`, only depends on `libm`. Builds on stable.

## License

BSD-3-Clause; see [`LICENSE`](LICENSE). The upstream
[`magcal.c`](https://github.com/PaulStoffregen/MotionCal/blob/master/magcal.c)
and [`matrix.c`](https://github.com/PaulStoffregen/MotionCal/blob/master/matrix.c)
that `src/solver.rs` and `src/matrix.rs` port from are BSD-3-Clause by
Freescale Semiconductor, Inc. (2014); the Freescale copyright is preserved
verbatim at the top of each ported file. Their original notice is preserved
verbatim in [`NOTICE-UPSTREAM`](NOTICE-UPSTREAM).

//! Solver-direct demo: 200 synthetic samples drawn from a known
//! ellipsoid (hard iron + a soft-iron deformation), pushed through
//! a [`Solver`] one at a time, then resolved with the full-ellipsoid
//! fitter.
//!
//! ```text
//! cargo run --example ellipsoid_fit
//! ```

use magcal::Solver;

fn main() {
    const HARD_IRON: [f32; 3] = [120.0, -80.0, 45.0];
    const W: [[f32; 3]; 3] = [
        [1.2, 0.05, 0.0],
        [0.05, 0.9, -0.03],
        [0.0, -0.03, 1.1],
    ];
    const RADIUS: f32 = 500.0;
    const N: usize = 200;

    let mut solver = Solver::new();
    for i in 0..N {
        let theta = (i as f32) * 0.3;
        let phi = (i as f32) * 0.7;
        let (st, ct) = (theta.sin(), theta.cos());
        let (sp, cp) = (phi.sin(), phi.cos());
        let p = [RADIUS * sp * ct, RADIUS * sp * st, RADIUS * cp];
        let raw = [
            (W[0][0] * p[0] + W[0][1] * p[1] + W[0][2] * p[2] + HARD_IRON[0]) as i16,
            (W[1][0] * p[0] + W[1][1] * p[1] + W[1][2] * p[2] + HARD_IRON[1]) as i16,
            (W[2][0] * p[0] + W[2][1] * p[1] + W[2][2] * p[2] + HARD_IRON[2]) as i16,
        ];
        solver.push_sample(raw);
    }

    let cal = solver.solve_ellipsoid().expect("ellipsoid fit failed");

    println!("recovered hard iron: {:?}", cal.hard_iron);
    println!("ground truth:        {:?}", HARD_IRON);
    println!("field strength:      {} (truth ~{})", cal.field_strength, RADIUS);
    println!("fit error:           {:.4} %", cal.fit_error_percent);
    println!("recovered soft iron (det-1 normalised):");
    for row in &cal.soft_iron {
        println!("  [{:7.4}, {:7.4}, {:7.4}]", row[0], row[1], row[2]);
    }
}

//! Minimal demo: 60 synthetic samples on a sphere displaced by a
//! known hard-iron offset, fit, and print the recovered offset.
//!
//! ```text
//! cargo run --example sphere_fit
//! ```

use magcal::MagCal;

fn main() {
    const HARD_IRON: [f32; 3] = [120.0, -80.0, 45.0];
    const RADIUS: f32 = 500.0;

    let mut samples = [[0_i16; 3]; 60];
    for (i, slot) in samples.iter_mut().enumerate() {
        let theta = (i as f32) * 0.3;
        let phi = (i as f32) * 0.7;
        let (st, ct) = (theta.sin(), theta.cos());
        let (sp, cp) = (phi.sin(), phi.cos());
        *slot = [
            (RADIUS * sp * ct + HARD_IRON[0]) as i16,
            (RADIUS * sp * st + HARD_IRON[1]) as i16,
            (RADIUS * cp + HARD_IRON[2]) as i16,
        ];
    }

    let cal = MagCal::fit(&samples).expect("fit failed");

    println!("recovered hard iron: {:?}", cal.hard_iron);
    println!("ground truth:        {:?}", HARD_IRON);
    println!("field strength:      {} (truth {})", cal.field_strength, RADIUS);
    println!("tier used:           {:?}", cal.tier);
}

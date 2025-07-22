mod core;
use core::{STDATM, asunoyozora, Meter, Kilogram, Second};

fn main() {
    let t_min = Second(0.0);
    let t_max = Second(65.0);
    let dt = Second(0.01);

    let tol= Meter(0.000001);

    let h_t = Meter(0.1);
    let h_h = Meter(1.6);
    let m = Kilogram(60.0);
    let h_width_rate: f64 = 0.25;

    let stdatm = STDATM::new(h_t, h_h, m, h_width_rate);

    // asunoyozora(t_min, t_max, dt, tol, &stdatm);

    println!("{}", stdatm.analytical_Z(t_max, Meter(20288.15135307835)).0);
}


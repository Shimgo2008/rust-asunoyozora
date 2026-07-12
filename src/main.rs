mod core;
use core::{STDATM, asunoyozora, Meter, Kilogram, Second};

fn main() {
    let t_min = Second(0.0);
    let t_max = Second(240.0);
    let dt = Second(0.01);

    let tol= Meter(0.001);

    let h_t = Meter(0.45);
    let h_h = Meter(1.612);
    let m = Kilogram(60.0);
    let h_width_rate: f64 = 0.23;

    let stdatm = STDATM::new(h_t, h_h, m, h_width_rate);

    asunoyozora(t_min, t_max, dt, tol, &stdatm);
}


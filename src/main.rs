mod core;
use core::{asunoyozora, Kilogram, Meter, Second, STDATM};

fn main() {
    let t_min = Second(0.0);
    let t_max = Second(177.0);
    let dt = Second(0.01);

    let tol = Meter(0.001);

    let h_t = Meter(0.30);
    let h_h = Meter(1.45);
    let m = Kilogram(43.0);
    let h_width_rate: f64 = 0.23;

    let stdatm = STDATM::new(h_t, h_h, m, h_width_rate);

    // print summary of the parameters
    print!("=== Parameters ===\n");
    println!("t_min: {} s", t_min.0);
    println!("t_max: {} s", t_max.0);
    println!("dt: {} s", dt.0);
    println!("tol: {} m", tol.0);

    print!("=== Body Parameters ===\n");
    println!("human thickness: {} m", h_t.0);
    println!("human height: {} m", h_h.0);
    println!("human mass: {} kg", m.0);
    println!("height width rate: {}", h_width_rate);

    asunoyozora(t_min, t_max, dt, tol, &stdatm);
}

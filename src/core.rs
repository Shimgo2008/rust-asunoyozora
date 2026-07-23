use std::ops::{Add, Div, Mul, Sub};

#[derive(Clone, Copy)]
pub struct Meter(pub f64);

#[derive(Clone, Copy)]
pub struct Kelvin(pub f64);

#[derive(Clone, Copy)]
pub struct GPMeter(pub f64);

#[derive(Clone, Copy)]
pub struct Pascal(pub f64);
#[derive(Clone, Copy)]
pub struct Kilogram(pub f64);
#[derive(Clone, Copy)]
pub struct Second(pub f64);
#[derive(Clone, Copy)]
struct SquareMeter(pub f64);

#[derive(Clone, Copy)]
pub struct Vec2(pub [f64; 2]);

impl Meter {
    pub const fn new(v: f64) -> Self {
        Self(v)
    }
}
impl Kelvin {
    pub const fn new(v: f64) -> Self {
        Self(v)
    }
}
impl GPMeter {
    pub const fn new(v: f64) -> Self {
        Self(v)
    }
}
impl Pascal {
    pub const fn new(v: f64) -> Self {
        Self(v)
    }
}
impl Kilogram {
    pub const fn new(v: f64) -> Self {
        Self(v)
    }
}
impl Second {
    pub const fn new(v: f64) -> Self {
        Self(v)
    }
}
impl Vec2 {
    pub const fn new(v: [f64; 2]) -> Self {
        Self(v)
    }
}

impl From<f64> for Second {
    fn from(value: f64) -> Self {
        Self(value)
    }
}

impl Add for Vec2 {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self::Output {
        return Vec2([self.0[0] + rhs.0[0], self.0[1] + rhs.0[1]]);
    }
}
impl Sub for Vec2 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self::Output {
        return Vec2([self.0[0] - rhs.0[0], self.0[1] - rhs.0[1]]);
    }
}
impl Mul for Vec2 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self::Output {
        return Vec2([self.0[0] * rhs.0[0], self.0[1] * rhs.0[1]]);
    }
}
impl Mul<f64> for Vec2 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: f64) -> Self::Output {
        return Vec2([self.0[0] * rhs, self.0[1] * rhs]);
    }
}
impl Div for Vec2 {
    type Output = Self;
    #[inline(always)]
    fn div(self, rhs: Self) -> Self::Output {
        return Vec2([self.0[0] / rhs.0[0], self.0[1] / rhs.0[1]]);
    }
}

const SA_LENGTH: usize = 8;
const R: f64 = 8.314462618e+3; // CODATA
const G_0: f64 = 9.80665; // CODATA
const M_0: f64 = 28.9644; // CODATA
const C_D: f64 = 1.0; // TODO: 難しいから定数 1である妥当性はない
const R_0: Meter = Meter(6356.766e3);
const H_I: [GPMeter; SA_LENGTH] = [
    GPMeter(0.0 * 1000.0),
    GPMeter(11.0 * 1000.0),
    GPMeter(20.0 * 1000.0),
    GPMeter(32.0 * 1000.0),
    GPMeter(47.0 * 1000.0),
    GPMeter(51.0 * 1000.0),
    GPMeter(71.0 * 1000.0),
    GPMeter(84.852 * 1000.0),
];
const P_B: [Pascal; SA_LENGTH] = [
    Pascal(101325.0),
    Pascal(22632.645891130254),
    Pascal(5475.162723334922),
    Pascal(868.0895580343575),
    Pascal(110.91927557768723),
    Pascal(66.94728114099739),
    Pascal(3.9571093634086405),
    Pascal(0.37346372866865546),
];
const T_M_B: [Kelvin; SA_LENGTH] = [
    Kelvin(288.15),
    Kelvin(216.65),
    Kelvin(216.65),
    Kelvin(228.65),
    Kelvin(270.65),
    Kelvin(270.65),
    Kelvin(214.65),
    Kelvin(186.95),
];
const L_M_B: [f64; SA_LENGTH] = [
    -6.5 / 1000.0,
    0.0 / 1000.0,
    1.0 / 1000.0,
    2.8 / 1000.0,
    0.0 / 1000.0,
    -2.8 / 1000.0,
    -2.0 / 1000.0,
    0.0 / 1000.0,
];

pub struct STDATM {
    h_t: Meter,
    h_h: Meter,
    m: Kilogram,
    h_width_rate: f64,
}

#[allow(non_snake_case)]
impl STDATM {
    pub const fn new(h_t: Meter, h_h: Meter, m: Kilogram, h_width_rate: f64) -> Self {
        Self {
            h_t,
            h_h,
            m,
            h_width_rate,
        }
    }

    pub fn analytical_Z(&self, t: Second, base_Z: Meter) -> Meter {
        let local_k = self.k(base_Z);
        let analytical_z_val = f64::sqrt((local_k * G_0) / self.m.0) * t.0;

        if !analytical_z_val.is_finite() || analytical_z_val.is_nan() || analytical_z_val < 0.0 {
            println!("[WARN] Invalid analytical_z_val: {}", analytical_z_val);
            return Meter(0.0);
        }

        // 上限で切る（Wasmでは700超えたらcoshは信用できない）
        if analytical_z_val > 700.0 {
            return Meter((self.m.0 / local_k) * (analytical_z_val - f64::ln(2.0)));
        }

        let cosh_val = f64::cosh(analytical_z_val);
        if !cosh_val.is_finite() || cosh_val <= 0.0 {
            println!("[WARN] Invalid cosh: {}", cosh_val);
            return Meter(0.0);
        }

        Meter((self.m.0 / local_k) * f64::ln(cosh_val))
    }

    fn a(&self, Z: Meter, v: f64) -> f64 {
        -self.g(Z) - (self.k(Z) * v.abs() * v) / self.m.0
    }

    fn g(&self, Z: Meter) -> f64 {
        let z = Z.0.max(0.0);
        G_0 * (R_0.0 / (R_0.0 + z)).powi(2)
    }

    fn k(&self, Z: Meter) -> f64 {
        (0.5) * self.A().0 * C_D * self.rho(Z)
    }

    fn A(&self) -> SquareMeter {
        SquareMeter((self.h_width_rate * self.h_h.0) * self.h_t.0)
    }

    /// 最適化後：Z → H, b を1回だけ求める
    fn rho(&self, Z: Meter) -> f64 {
        let H = self.Z2H(Meter(Z.0.max(0.0)));
        let b = self.get_b(H);
        let T = self.T_from_H(H, b);
        let P = self.P_from_H(H, b);
        let rho = self.rho_from_TP(T, P);

        assert!(
            rho.is_finite() && rho > 0.0,
            "invalid atmosphere: H={}, b={}, T={}, P={}, rho={}",
            H.0,
            b,
            T.0,
            P.0,
            rho
        );

        rho
    }

    /// Hとbが既知のときの温度
    fn T_from_H(&self, H: GPMeter, b: usize) -> Kelvin {
        let dh = H.0 - H_I[b].0;
        Kelvin(T_M_B[b].0 + L_M_B[b] * dh)
    }

    /// Hとbが既知のときの気圧
    fn P_from_H(&self, H: GPMeter, b: usize) -> Pascal {
        let p_b = P_B[b].0;
        let t_b = T_M_B[b].0;
        let lapse = L_M_B[b];
        let dh = H.0 - H_I[b].0;

        if lapse.abs() < f64::EPSILON {
            Pascal(p_b * (-(G_0 * M_0 * dh) / (R * t_b)).exp())
        } else {
            let t = t_b + lapse * dh;
            Pascal(p_b * (t_b / t).powf((G_0 * M_0) / (R * lapse)))
        }
    }

    fn rho_from_TP(&self, T: Kelvin, P: Pascal) -> f64 {
        (P.0 * M_0) / (R * T.0)
    }

    /// b値の取得（Z2Hを省略するため引数はGPMeter）
    fn get_b(&self, H: GPMeter) -> usize {
        assert!(
            H.0 <= H_I[SA_LENGTH - 1].0 + 1.0e-9,
            "geopotential height {} m exceeds model limit {} m",
            H.0,
            H_I[SA_LENGTH - 1].0
        );

        if H.0 <= H_I[0].0 {
            return 0;
        }

        for b in 0..(SA_LENGTH - 1) {
            if H.0 < H_I[b + 1].0 {
                return b;
            }
        }

        SA_LENGTH - 1
    }

    fn Z2H(&self, Z: Meter) -> GPMeter {
        GPMeter((R_0.0 * Z.0) / (R_0.0 + Z.0))
    }

    fn H2Z(&self, H: GPMeter) -> Meter {
        Meter((R_0.0 * H.0) / (R_0.0 - H.0))
    }

    fn max_model_Z(&self) -> Meter {
        self.H2Z(H_I[SA_LENGTH - 1])
    }
}

fn rk4<F>(
    f: F,
    y0: Vec2,
    t0: Second,
    t1: Second,
    dt: Second,
    stdatm: &STDATM,
    is_history: bool,
) -> Vec<Vec2>
where
    F: Fn(&STDATM, Vec2) -> Vec2,
{
    assert!(dt.0 > 0.0);
    assert!(t1.0 >= t0.0);

    let mut results = Vec::new();
    let mut y: Vec2 = y0;
    let mut t = t0.0;
    let mut step = 0usize;

    if is_history {
        results.push(y0);
    }

    while t < t1.0 {
        let h = dt.0.min(t1.0 - t);

        let k1: Vec2 = f(stdatm, y);
        let k2: Vec2 = f(stdatm, y + k1 * (h / 2.0));
        let k3: Vec2 = f(stdatm, y + k2 * (h / 2.0));
        let k4: Vec2 = f(stdatm, y + k3 * h);

        y = y + ((k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0));
        t += h;
        step += 1;

        if is_history && (step % 2 == 0 || t >= t1.0) {
            results.push(y);
        }
    }

    if !is_history {
        results.push(y);
    }

    return results;
}

#[allow(non_snake_case)]
fn vector_ODE(std: &STDATM, y: Vec2) -> Vec2 {
    let Z = Meter(y.0[0]);
    let v = y.0[1];
    let dZ_dt = v;
    let dv_dt = std.a(Z, v);
    return Vec2([dZ_dt, dv_dt]);
}

fn final_altitude(z0: f64, t_min: Second, t_max: Second, dt: Second, stdatm: &STDATM) -> f64 {
    assert!(z0 >= 0.0, "initial altitude must be non-negative: {}", z0);
    assert!(
        z0 <= stdatm.max_model_Z().0,
        "initial altitude {} m exceeds atmosphere model limit {} m",
        z0,
        stdatm.max_model_Z().0
    );

    let y0 = Vec2::new([z0, 0.0]);
    let z = rk4(vector_ODE, y0, t_min, t_max, dt, stdatm, false)
        .last()
        .expect("RK4 returned no result")
        .0[0];

    assert!(z.is_finite(), "final altitude is not finite: {}", z);
    z
}

fn bracket_initial_altitude(
    t_min: Second,
    t_max: Second,
    dt: Second,
    stdatm: &STDATM,
) -> (f64, f64) {
    let model_limit = stdatm.max_model_Z().0;
    let low = 0.0;
    let f_low = final_altitude(low, t_min, t_max, dt, stdatm);

    let mut high = stdatm.analytical_Z(t_max, Meter(0.0)).0;
    if !high.is_finite() || high <= 0.0 {
        high = 1.0;
    }
    high = high.min(model_limit);

    let mut f_high = final_altitude(high, t_min, t_max, dt, stdatm);

    const MAX_BRACKET_ITERATIONS: usize = 64;

    for _ in 0..MAX_BRACKET_ITERATIONS {
        if f_low <= 0.0 && f_high >= 0.0 {
            return (low, high);
        }

        if f_low > 0.0 {
            break;
        }

        if high >= model_limit {
            break;
        }

        high = (high * 2.0).min(model_limit);
        f_high = final_altitude(high, t_min, t_max, dt, stdatm);
    }

    assert!(
        f_low <= 0.0 && f_high >= 0.0,
        "failed to bracket root within atmosphere model: low={}, f(low)={}, high={}, f(high)={}, model_limit={}",
        low,
        f_low,
        high,
        f_high,
        model_limit
    );

    (low, high)
}

#[allow(non_snake_case)]
pub fn asunoyozora(
    t_min: Second,
    t_max: Second,
    dt: Second,
    tol: Meter,
    stdatm: &STDATM,
) -> Vec<Vec2> {
    let (mut low, mut high) = bracket_initial_altitude(t_min, t_max, dt, stdatm);

    const MAX_ITERATIONS: i32 = 100;

    for i in 0..MAX_ITERATIONS {
        println!(
            "[Iter {}] low: {:.4}, high: {:.4}, diff: {:.8}",
            i,
            low,
            high,
            high - low
        );

        if high - low < tol.0 {
            println!("Tolerance reached. Exiting loop.");
            break;
        }

        let mid = low + (high - low) / 2.0;

        if !mid.is_finite() || mid < 0.0 {
            println!("Error: mid is not finite or is negative. Breaking.");

            high = low;
            break;
        }

        let f_mid = final_altitude(mid, t_min, t_max, dt, stdatm);

        if f_mid > 0.0 {
            high = mid;
        } else {
            low = mid;
        }
    }

    let optimal_initial_altitude = (low + high) / 2.0;
    println!(
        "Final estimated altitude: {:.6} m",
        optimal_initial_altitude
    );

    let final_y0 = Vec2::new([optimal_initial_altitude, 0.0]);
    let final_history = rk4(vector_ODE, final_y0, t_min, t_max, dt, &stdatm, true);

    if let Some(last_vec) = final_history.last() {
        if !last_vec.0[1].is_finite() {
            println!("FATAL: Final simulation also resulted in non-finite velocity!");
        }
    }

    println!("{}", final_history[0].0[0]);
    return final_history;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stdatm() -> STDATM {
        STDATM::new(Meter(0.45), Meter(1.612), Kilogram(60.0), 0.23)
    }

    #[test]
    fn get_b_selects_standard_atmosphere_layer() {
        let atmosphere = stdatm();

        assert_eq!(atmosphere.get_b(GPMeter(-1.0)), 0);
        assert_eq!(atmosphere.get_b(GPMeter(0.0)), 0);
        assert_eq!(atmosphere.get_b(GPMeter(10_999.0)), 0);
        assert_eq!(atmosphere.get_b(GPMeter(11_000.0)), 1);
        assert_eq!(atmosphere.get_b(GPMeter(20_000.0)), 2);
        assert_eq!(atmosphere.get_b(H_I[SA_LENGTH - 1]), 7);
    }

    #[test]
    #[should_panic(expected = "exceeds model limit")]
    fn get_b_rejects_height_above_model_limit() {
        stdatm().get_b(GPMeter(90_000.0));
    }

    #[test]
    fn rho_uses_sea_level_standard_atmosphere() {
        let atmosphere = stdatm();
        let rho = atmosphere.rho(Meter(0.0));

        assert!((rho - 1.2249781434738449).abs() < 1.0e-12);
    }

    #[test]
    fn pressure_matches_layer_base_values() {
        let atmosphere = stdatm();

        for b in 0..SA_LENGTH {
            let p = atmosphere.P_from_H(H_I[b], b);
            assert!((p.0 - P_B[b].0).abs() < 1.0e-9);
        }
    }

    #[test]
    fn bracket_initial_altitude_checks_final_altitude_signs() {
        let atmosphere = stdatm();
        let t_min = Second(0.0);
        let t_max = Second(177.0);
        let dt = Second(0.05);

        let (low, high) = bracket_initial_altitude(t_min, t_max, dt, &atmosphere);

        assert!(final_altitude(low, t_min, t_max, dt, &atmosphere) <= 0.0);
        assert!(final_altitude(high, t_min, t_max, dt, &atmosphere) >= 0.0);
        assert!(high <= atmosphere.max_model_Z().0);
    }
}

// for debug
#[allow(dead_code)]
fn main() {
    let t_min = Second(0.0);
    let t_max = Second(62.49);
    let dt = Second(0.01);

    let tol = Meter(0.000001);

    let h_t = Meter(0.1);
    let h_h = Meter(1.6);
    let m = Kilogram(60.0);
    let h_width_rate: f64 = 0.25;

    let stdatm = STDATM::new(h_t, h_h, m, h_width_rate);

    asunoyozora(t_min, t_max, dt, tol, &stdatm);

    // println!("{}", stdatm.analytical_Z(t_max, Meter(20288.15135307835)).0);
}

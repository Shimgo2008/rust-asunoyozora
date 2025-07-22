use std::{ops::{Add, Div, Mul, Sub}};

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
    pub const fn new(v: [f64; 2]) -> Self{
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
const R: f64 = 8.314462618e+3;      // CODATA
const G_0: f64 = 9.80665;           // CODATA
const M_0: f64 = 28.9644;           // CODATA
const C_D: f64 = 1.0;               // TODO: 難しいから定数 1である妥当性はない
const R_0: Meter = Meter(6356.766e3);
const H_I: [GPMeter; SA_LENGTH] = [GPMeter(0.0 * 1000.0), GPMeter(11.0 * 1000.0), GPMeter(20.0 * 1000.0), GPMeter(32.0 * 1000.0), GPMeter(47.0 * 1000.0), GPMeter(51.0 * 1000.0), GPMeter(71.0 * 1000.0), GPMeter(84.852 * 1000.0)];
const P_B: [Pascal; SA_LENGTH] = [Pascal(101325.0), Pascal(17881.924167776844), Pascal(4451.782721266618), Pascal(835.676386546647), Pascal(125.823085025137), Pascal(75.131157658166), Pascal(2.218089806735), Pascal(0.010142033211)];
const T_M_B: [Kelvin; SA_LENGTH] = [Kelvin(288.15), Kelvin(216.65), Kelvin(216.65), Kelvin(228.65), Kelvin(270.65), Kelvin(270.65), Kelvin(214.65), Kelvin(186.95)];
const L_M_B: [f64; SA_LENGTH] = [-6.5 / 1000.0, 0.0 / 1000.0, 1.0 / 1000.0, 2.8 / 1000.0, 0.0 / 1000.0, -2.8 / 1000.0, -2.0 / 1000.0, 0.0 / 1000.0];

pub struct STDATM {
    h_t: Meter,
    h_h: Meter,
    m: Kilogram,
    h_width_rate: f64,
}

#[allow(non_snake_case)]
impl STDATM {
    pub const fn new(h_t: Meter, h_h: Meter, m: Kilogram, h_width_rate: f64) -> Self {
        Self { h_t, h_h, m, h_width_rate }
    }


    pub fn analytical_Z(&self, t: Second, base_Z: Meter) -> Meter {
        let local_k = self.k(base_Z);
        println!("k: {}", local_k);
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
        -G_0 - (self.k(Z) * v * v) / self.m.0
    }

    fn k(&self, Z: Meter) -> f64 {
        (0.5) * self.A().0 * C_D * self.rho(Z)
    }

    fn A(&self) -> SquareMeter {
        SquareMeter((self.h_width_rate * self.h_h.0) * self.h_t.0)
    }

    /// 最適化後：Z → H, b を1回だけ求める
    fn rho(&self, Z: Meter) -> f64 {
        let H = self.Z2H(Z);
        let b = self.get_b(H);
        let T = self.T_from_H(H, b);
        let P = self.P_from_H(H, b);
        println!("H: {}", H.0);
        println!("b: {}", b);
        println!("T: {}", T.0);
        println!("P: {}", P.0);
        self.rho_from_TP(T, P)
    }

    /// Hとbが既知のときの温度
    fn T_from_H(&self, H: GPMeter, b: usize) -> Kelvin {
        if b == 0 {
            return T_M_B[0];  // もしくは定数返すなど
        }
        let rate = L_M_B[b];
        let remainder_kilometer = (H.0 - H_I[b - 1].0) / 1000.0;
        Kelvin(T_M_B[b - 1].0 + remainder_kilometer * rate)
    }

    /// Hとbが既知のときの気圧
    fn P_from_H(&self, H: GPMeter, b: usize) -> Pascal {
        let P_b = P_B[b];
        if L_M_B[b] == 0.0 {
            self.P_33b(H, b, P_b)
        } else {
            self.P_33a(H, b, P_b)
        }
    }

    fn rho_from_TP(&self, T: Kelvin, P: Pascal) -> f64 {
        (P.0 * M_0) / (R * T.0)
    }

    /// 最適化後：bを引数で受け取る（P_33a用）
    fn P_33a(&self, H: GPMeter, b: usize, P_b: Pascal) -> Pascal {
        if b == 0 {
            return P_b; // 変化なし or 別の処理にする
        }
        let P33a_value = (G_0 * M_0) / (R * L_M_B[b]);
        println!("P33a: {}", P33a_value);
        let base = ((T_M_B[b].0) / (T_M_B[b].0 + L_M_B[b] * (H.0 - H_I[b - 1].0)));
        println!("base: {}", base);
        Pascal(P_b.0 * f64::powf(
            base,
            P33a_value,
        ))
    }

    /// 最適化後：bを引数で受け取る（P_33b用）
    fn P_33b(&self, H: GPMeter, b: usize, P_b: Pascal) -> Pascal {
        if b == 0 {
            return P_b; // 変化なし or 別の処理にする
        }
        let P33b_value = (G_0 * M_0) / (R * T_M_B[b].0);
        println!("P33b: {}", P33b_value);
        Pascal(P_b.0 * f64::exp(-1.0 * P33b_value * (H.0 - H_I[b - 1].0)))
    }

    /// b値の取得（Z2Hを省略するため引数はGPMeter）
    fn get_b(&self, H: GPMeter) -> usize {
        if H.0 == 0.0 {
            return 0;
        }
        for i in 1..(SA_LENGTH - 1) {
            if H_I[i].0 <= H.0 && H.0 < H_I[i + 1].0 {
                return i;
            }
        }
        SA_LENGTH - 1
    }

    fn Z2H(&self, Z: Meter) -> GPMeter {
        GPMeter((R_0.0 * Z.0) / (R_0.0 + Z.0))
    }
}

fn rk4<F>(f: F,
    y0: Vec2,
    t0: Second,
    t1: Second,
    dt: Second,
    stdatm: &STDATM,
) -> Vec<Vec2>
where
    F: Fn(&STDATM, Vec2) -> Vec2
{
    let steps = ((t1.0 - t0.0) / dt.0) as usize;
    let mut results = Vec::with_capacity((steps + 1) / 2);
    let mut y: Vec2 = y0;
    results.push(y0);

    for step in 1..steps{
        let k1: Vec2 = f(stdatm, y);
        let k2: Vec2 = f(stdatm, y + k1 * (dt.0 / 2.0));
        let k3: Vec2 = f(stdatm, y + k2 * (dt.0 / 2.0));
        let k4: Vec2 = f(stdatm, y + k3 * dt.0);
        y = y + ((k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt.0 / 6.0));
        if step % 2 == 0 {
            results.push(y);
        }
    }

    return results;
}

#[allow(non_snake_case)]
fn vector_ODE(std: &STDATM, y: Vec2) -> Vec2{
    let Z = Meter(y.0[0]);
    let v = y.0[1];
    let dZ_dt = v;
    let dv_dt = std.a(Z, v);
    return Vec2([dZ_dt, dv_dt]);
}

#[allow(non_snake_case)]
pub fn asunoyozora(t_min: Second, t_max: Second, dt: Second, tol: Meter, stdatm: &STDATM) -> Vec<Vec2> {
    let mut analytical_Zlist = vec![Meter(0.0)];
    for _ in 0..10 {
        let last = analytical_Zlist.last().unwrap();
        let next = stdatm.analytical_Z(t_max, *last);
        println!("next: {}, last: {}", next.0, last.0);
        analytical_Zlist.push(next);
    }

    let mut low = analytical_Zlist.last().unwrap().0 / 2.0;
    let mut high = analytical_Zlist.last().unwrap().0 * 2.0;

    const MAX_ITERATIONS: i32 = 100;

    for i in 0..MAX_ITERATIONS {
        println!("[Iter {}] low: {:.4}, high: {:.4}, diff: {:.8}", i, low, high, high - low);

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

        let y0 = Vec2::new([mid, 0.0]);

        let history = rk4(vector_ODE, y0, t_min, t_max, dt, &stdatm);

        let final_altitude = match history.last() {
            Some(last_vec) if last_vec.0[0].is_finite() => last_vec.0[0],
            _ => {
                println!("Warning: Simulation diverged at mid={:.4}. Assuming height is too high.", mid);
                high = mid;
                continue
            }
        };

        if final_altitude > 0.0 {
            high = mid;
        } else {
            low = mid;
        }
    }

    let optimal_initial_altitude = (low + high) / 2.0;
    println!("Final estimated altitude: {:.6} m", optimal_initial_altitude);

    let final_y0 = Vec2::new([optimal_initial_altitude, 0.0]);
    let final_history = rk4(vector_ODE, final_y0, t_min, t_max, dt, &stdatm);

    if let Some(last_vec) = final_history.last() {
        if !last_vec.0[1].is_finite() {
             println!("FATAL: Final simulation also resulted in non-finite velocity!");
        }
    }

    return final_history;
}

fn main() {
    let t_min = Second(0.0);
    let t_max = Second(50.0);
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


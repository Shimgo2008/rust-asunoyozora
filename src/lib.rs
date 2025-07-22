mod core;

use core::{STDATM, asunoyozora, Meter, Kilogram, Second};
use wasm_bindgen::prelude::*;
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize)]
pub struct AsunoyozoraResult {
    pub altitude: f64,
    pub velocity: f64,
}

#[wasm_bindgen]
extern "C" {
    pub fn alert(s: &str);
}

#[wasm_bindgen]
pub fn asunoyozora_wasm(
    t_min: f64,
    t_max: f64,
    dt: f64,
    tol: f64,
    h_t: f64,
    h_h: f64,
    m: f64,
    h_width_rate: f64,
) -> Result<JsValue, JsValue> {
    let stdatm = STDATM::new(
        Meter(h_t),
        Meter(h_h),
        Kilogram(m),
        h_width_rate,
    );

    let results = asunoyozora(
        Second(t_min),
        Second(t_max),
        Second(dt),
        Meter(tol),
        &stdatm,
    );

    let js_results: Vec<AsunoyozoraResult> = results.into_iter().map(|v| AsunoyozoraResult { altitude: v.0[0], velocity: v.0[1] }).collect();

    Ok(serde_wasm_bindgen::to_value(&js_results)?)
}


// #[pymodule]
// fn rust_asunoyozora(m: &Bound<'_, PyModule>) -> PyResult<()> {
//     m.add_function(wrap_pyfunction!(py_asunoyozora, m)?)
// }

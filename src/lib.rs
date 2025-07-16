mod core;

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::wrap_pyfunction;
#[cfg(feature = "python")]
use pyo3::types::PyModule;

use core::{STDATM, asunoyozora, Meter, Kilogram, Second};

#[cfg(feature = "python")]
#[pyfunction]
fn py_asunoyozora(
    t_min: f64,
    t_max: f64,
    dt: f64,
    tol: f64,
    h_t: f64,
    h_h: f64,
    m: f64,
    h_width_rate: f64,
) -> Vec<(f64, f64)> {
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

    results.into_iter().map(|v| (v.0[0], v.0[1])).collect()
}

#[cfg(feature = "python")]
#[pymodule]
fn pyo3_example(m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_asunoyozora, m)?)?;
    Ok(())
}

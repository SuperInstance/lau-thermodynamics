//! Ideal gas law and real gas corrections (van der Waals).

use crate::constants::R;
use serde::{Serialize, Deserialize};

/// Ideal gas law: PV = nRT → returns pressure.
pub fn ideal_gas_pressure(n: f64, t: f64, v: f64) -> f64 {
    n * R * t / v
}

/// Ideal gas law → returns volume.
pub fn ideal_gas_volume(n: f64, t: f64, p: f64) -> f64 {
    n * R * t / p
}

/// Ideal gas law → returns temperature.
pub fn ideal_gas_temperature(p: f64, v: f64, n: f64) -> f64 {
    p * v / (n * R)
}

/// Ideal gas law → returns number of moles.
pub fn ideal_gas_moles(p: f64, v: f64, t: f64) -> f64 {
    p * v / (R * t)
}

/// Van der Waals equation parameters for a specific gas.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VanDerWaalsParams {
    /// Attraction parameter a (Pa·m⁶/mol²)
    pub a: f64,
    /// Volume parameter b (m³/mol)
    pub b: f64,
}

/// Common van der Waals parameters.
impl VanDerWaalsParams {
    pub fn nitrogen() -> Self {
        Self { a: 0.1370, b: 3.87e-5 }
    }

    pub fn oxygen() -> Self {
        Self { a: 0.1398, b: 3.18e-5 }
    }

    pub fn carbon_dioxide() -> Self {
        Self { a: 0.3643, b: 4.27e-5 }
    }

    pub fn water() -> Self {
        Self { a: 0.5537, b: 3.05e-5 }
    }

    pub fn helium() -> Self {
        Self { a: 0.00346, b: 2.38e-5 }
    }
}

/// Van der Waals pressure: P = nRT/(V - nb) - a(n/V)²
pub fn vdw_pressure(n: f64, t: f64, v: f64, params: &VanDerWaalsParams) -> f64 {
    let nb = n * params.b;
    (n * R * t) / (v - nb) - params.a * (n / v).powi(2)
}

/// Van der Waals volume solved iteratively via Newton's method.
pub fn vdw_volume(n: f64, t: f64, p: f64, params: &VanDerWaalsParams, iterations: usize) -> f64 {
    let a = params.a;
    let b = params.b;
    // Start from ideal gas estimate
    let mut v = n * R * t / p;
    for _ in 0..iterations {
        let nb = n * b;
        let f = (n * R * t) / (v - nb) - a * (n / v).powi(2) - p;
        let df = -(n * R * t) / (v - nb).powi(2) + 2.0 * a * n.powi(2) / v.powi(3);
        if df.abs() < 1e-30 { break; }
        v -= f / df;
        if v <= nb { v = nb + 1e-10; }
    }
    v
}

/// Compression factor Z = PV/(nRT). For ideal gas Z = 1.
pub fn compressibility_factor(p: f64, v: f64, n: f64, t: f64) -> f64 {
    p * v / (n * R * t)
}

/// Boyle temperature for a van der Waals gas: T_B = a/(Rb).
pub fn boyle_temperature(params: &VanDerWaalsParams) -> f64 {
    params.a / (R * params.b)
}

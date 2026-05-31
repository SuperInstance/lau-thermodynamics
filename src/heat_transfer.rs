//! Heat transfer — conduction and convection basics.

use crate::constants::STEFAN_BOLTZMANN;
use serde::{Serialize, Deserialize};

/// Fourier's law of heat conduction: Q = -k·A·dT/dx
/// Steady-state: Q = k·A·(T_hot - T_cold) / L
pub fn conduction_heat(k: f64, area: f64, t_hot: f64, t_cold: f64, thickness: f64) -> f64 {
    if thickness <= 0.0 { return 0.0; }
    k * area * (t_hot - t_cold) / thickness
}

/// Thermal resistance for conduction: R_th = L / (k·A)
pub fn thermal_resistance_conduction(thickness: f64, k: f64, area: f64) -> f64 {
    if k <= 0.0 || area <= 0.0 { return f64::INFINITY; }
    thickness / (k * area)
}

/// Total thermal resistance for series: R_total = Σ R_i
pub fn thermal_resistance_series(resistances: &[f64]) -> f64 {
    resistances.iter().sum()
}

/// Total thermal resistance for parallel: 1/R_total = Σ 1/R_i
pub fn thermal_resistance_parallel(resistances: &[f64]) -> f64 {
    if resistances.is_empty() { return f64::INFINITY; }
    1.0 / resistances.iter().map(|r| 1.0 / r).sum::<f64>()
}

/// Newton's law of cooling (convection): Q = h·A·(T_surface - T_fluid)
pub fn convection_heat(h: f64, area: f64, t_surface: f64, t_fluid: f64) -> f64 {
    h * area * (t_surface - t_fluid)
}

/// Thermal resistance for convection: R_th = 1/(h·A)
pub fn thermal_resistance_convection(h: f64, area: f64) -> f64 {
    if h <= 0.0 || area <= 0.0 { return f64::INFINITY; }
    1.0 / (h * area)
}

/// Radiative heat transfer: Q = ε·σ·A·(T₁⁴ - T₂⁴)
pub fn radiation_heat(emissivity: f64, area: f64, t1: f64, t2: f64) -> f64 {
    emissivity * STEFAN_BOLTZMANN * area * (t1.powi(4) - t2.powi(4))
}

/// Common thermal conductivities (W/(m·K)).
pub fn thermal_conductivity_copper() -> f64 { 401.0 }
pub fn thermal_conductivity_aluminum() -> f64 { 237.0 }
pub fn thermal_conductivity_iron() -> f64 { 80.0 }
pub fn thermal_conductivity_glass() -> f64 { 1.05 }
pub fn thermal_conductivity_water() -> f64 { 0.60 }
pub fn thermal_conductivity_air() -> f64 { 0.026 }

/// Heat transfer result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HeatTransferResult {
    pub conduction: f64,
    pub convection: f64,
    pub radiation: f64,
    pub total: f64,
}

impl HeatTransferResult {
    pub fn new(cond: f64, conv: f64, rad: f64) -> Self {
        Self { conduction: cond, convection: conv, radiation: rad, total: cond + conv + rad }
    }
}

/// Transient conduction: Biot number Bi = hL/k
pub fn biot_number(h: f64, characteristic_length: f64, k: f64) -> f64 {
    if k <= 0.0 { return f64::INFINITY; }
    h * characteristic_length / k
}

/// Fourier number Fo = α·t/L² (dimensionless time for transient conduction)
pub fn fourier_number(thermal_diffusivity: f64, time: f64, characteristic_length: f64) -> f64 {
    if characteristic_length <= 0.0 { return 0.0; }
    thermal_diffusivity * time / (characteristic_length * characteristic_length)
}

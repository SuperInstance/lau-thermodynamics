//! Physical constants used in thermodynamics.

use serde::{Serialize, Deserialize};

/// Universal gas constant R = 8.314462618 J/(mol·K)
pub const R: f64 = 8.314462618;

/// Boltzmann constant k_B = 1.380649e-23 J/K
pub const BOLTZMANN: f64 = 1.380649e-23;

/// Avogadro's number N_A = 6.02214076e23 /mol
pub const AVOGADRO: f64 = 6.02214076e23;

/// Stefan-Boltzmann constant σ = 5.670374419e-8 W/(m²·K⁴)
pub const STEFAN_BOLTZMANN: f64 = 5.670374419e-8;

/// Standard atmospheric pressure (Pa)
pub const ATM: f64 = 101_325.0;

/// Celsius-to-Kelvin offset
pub const C_TO_K: f64 = 273.15;

/// Convert Celsius to Kelvin.
pub fn celsius_to_kelvin(c: f64) -> f64 {
    c + C_TO_K
}

/// Convert Kelvin to Celsius.
pub fn kelvin_to_celsius(k: f64) -> f64 {
    k - C_TO_K
}

/// Molar gas constant helper — returns R.
pub fn gas_constant() -> f64 {
    R
}

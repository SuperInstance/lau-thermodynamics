//! # lau-thermodynamics
//!
//! Classical and statistical thermodynamics — energy, entropy, and equilibrium.
//!
//! Covers the four laws of thermodynamics, ideal and real gas behavior,
//! Carnot cycles, entropy, Maxwell relations, phase transitions,
//! statistical mechanics, heat transfer, and agent energy budgets.

pub mod constants;
pub mod laws;
pub mod gas;
pub mod carnot;
pub mod entropy;
pub mod maxwell;
pub mod phase;
pub mod statistical;
pub mod heat_transfer;
pub mod agent_budget;

pub use constants::*;
pub use laws::*;
pub use gas::*;
pub use carnot::*;
pub use entropy::*;
pub use maxwell::*;
pub use phase::*;
pub use statistical::*;
pub use heat_transfer::*;
pub use agent_budget::*;

#[cfg(test)]
mod tests;

use serde::{Serialize, Deserialize};

/// Core thermodynamic state point.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermoState {
    pub temperature: f64, // Kelvin
    pub pressure: f64,    // Pascals
    pub volume: f64,      // m³
    pub internal_energy: f64, // Joules
    pub entropy: f64,     // J/K
    pub enthalpy: f64,    // Joules
    pub gibbs: f64,       // Joules
}

impl ThermoState {
    /// Compute Gibbs free energy from enthalpy, temperature, and entropy.
    pub fn gibbs_free_energy(h: f64, t: f64, s: f64) -> f64 {
        h - t * s
    }

    /// Compute Helmholtz free energy from internal energy, temperature, and entropy.
    pub fn helmholtz_free_energy(u: f64, t: f64, s: f64) -> f64 {
        u - t * s
    }
}

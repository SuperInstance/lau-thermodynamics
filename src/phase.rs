//! Phase transitions — Clausius-Clapeyron equation.

use crate::constants::R;
use serde::{Serialize, Deserialize};

/// Clausius-Clapeyron equation: dP/dT = ΔH/(T·ΔV)
pub fn clausius_clapeyron_slope(delta_h: f64, t: f64, delta_v: f64) -> f64 {
    if t <= 0.0 || delta_v == 0.0 { return f64::NAN; }
    delta_h / (t * delta_v)
}

/// Integrated Clausius-Clapeyron for liquid-vapor (assuming vapor is ideal gas):
/// ln(P₂/P₁) = -ΔH_vap/R · (1/T₂ - 1/T₁)
pub fn clausius_clapeyron_integrated(
    p1: f64,
    t1: f64,
    t2: f64,
    delta_h_vap: f64,
) -> f64 {
    p1 * ((-delta_h_vap / R) * (1.0 / t2 - 1.0 / t1)).exp()
}

/// Boiling point at a given pressure from reference boiling point.
pub fn boiling_point_at_pressure(
    p_ref: f64,
    t_boil_ref: f64,
    p_target: f64,
    delta_h_vap: f64,
) -> f64 {
    // ln(P_target/P_ref) = -ΔH/R · (1/T_target - 1/T_ref)
    // 1/T_target = 1/T_ref - R/ΔH · ln(P_target/P_ref)
    let inv_t_target = 1.0 / t_boil_ref - (R / delta_h_vap) * (p_target / p_ref).ln();
    1.0 / inv_t_target
}

/// Common phase transition data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhaseTransition {
    pub name: String,
    pub t_transition: f64, // K
    pub delta_h: f64,      // J/mol
    pub delta_s: f64,      // J/(mol·K)
}

impl PhaseTransition {
    /// Water boiling at 1 atm.
    pub fn water_boiling() -> Self {
        Self {
            name: "Water → Steam (1 atm)".into(),
            t_transition: 373.15,
            delta_h: 40_700.0,
            delta_s: 40_700.0 / 373.15,
        }
    }

    /// Water freezing at 1 atm.
    pub fn water_freezing() -> Self {
        Self {
            name: "Water → Ice (1 atm)".into(),
            t_transition: 273.15,
            delta_h: -6_010.0,
            delta_s: -6_010.0 / 273.15,
        }
    }

    /// Verify ΔG = ΔH - TΔS ≈ 0 at equilibrium.
    pub fn verify_equilibrium(&self) -> f64 {
        self.delta_h - self.t_transition * self.delta_s
    }
}

/// Triple point description.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TriplePoint {
    pub temperature: f64,
    pub pressure: f64,
}

impl TriplePoint {
    pub fn water() -> Self {
        Self {
            temperature: 273.16,
            pressure: 611.657,
        }
    }
}

/// Critical point description.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CriticalPoint {
    pub temperature: f64,
    pub pressure: f64,
    pub volume: f64, // Molar volume
}

impl CriticalPoint {
    pub fn water() -> Self {
        Self {
            temperature: 647.096,
            pressure: 22_064_000.0,
            volume: 5.59e-5,
        }
    }

    /// Convert critical point to van der Waals parameters.
    /// At critical point: P_c = a/(27b²), T_c = 8a/(27Rb), V_c = 3b
    pub fn to_vdw_params(&self) -> (f64, f64) {
        let b = self.volume / 3.0;
        let a = 27.0 * b * b * self.pressure;
        (a, b)
    }
}

/// Gibbs phase rule: F = C - P + 2 (degrees of freedom)
pub fn gibbs_phase_rule(components: i32, phases: i32) -> i32 {
    let f = components - phases + 2;
    f.max(0)
}

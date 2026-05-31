//! Carnot cycle and heat engine efficiency.

use serde::{Serialize, Deserialize};

/// Compute Carnot efficiency: η = 1 - T_cold / T_hot
pub fn carnot_efficiency(t_hot: f64, t_cold: f64) -> f64 {
    if t_hot <= 0.0 || t_cold < 0.0 || t_cold >= t_hot {
        return 0.0;
    }
    1.0 - t_cold / t_hot
}

/// Work output of a Carnot engine: W = Q_hot · η
pub fn carnot_work(q_hot: f64, t_hot: f64, t_cold: f64) -> f64 {
    q_hot * carnot_efficiency(t_hot, t_cold)
}

/// Heat rejected to cold reservoir: Q_cold = Q_hot - W
pub fn carnot_heat_rejected(q_hot: f64, t_hot: f64, t_cold: f64) -> f64 {
    q_hot - carnot_work(q_hot, t_hot, t_cold)
}

/// Coefficient of performance for a Carnot refrigerator: COP_R = T_cold / (T_hot - T_cold)
pub fn cop_refrigerator(t_hot: f64, t_cold: f64) -> f64 {
    if t_hot <= t_cold { return 0.0; }
    t_cold / (t_hot - t_cold)
}

/// Coefficient of performance for a Carnot heat pump: COP_HP = T_hot / (T_hot - T_cold)
pub fn cop_heat_pump(t_hot: f64, t_cold: f64) -> f64 {
    if t_hot <= t_cold { return 0.0; }
    t_hot / (t_hot - t_cold)
}

/// A complete Carnot cycle result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CarnotCycle {
    pub t_hot: f64,
    pub t_cold: f64,
    pub efficiency: f64,
    pub q_hot: f64,
    pub q_cold: f64,
    pub work: f64,
    pub entropy_change: f64, // Net entropy change = 0 for reversible
}

impl CarnotCycle {
    /// Build a full Carnot cycle from reservoir temperatures and heat input.
    pub fn new(t_hot: f64, t_cold: f64, q_hot: f64) -> Self {
        let efficiency = carnot_efficiency(t_hot, t_cold);
        let work = q_hot * efficiency;
        let q_cold = q_hot - work;
        Self {
            t_hot,
            t_cold,
            efficiency,
            q_hot,
            q_cold,
            work,
            entropy_change: 0.0, // Reversible cycle
        }
    }
}

/// Actual thermal efficiency given heat input and work output.
pub fn thermal_efficiency(q_in: f64, work_out: f64) -> f64 {
    if q_in <= 0.0 { return 0.0; }
    work_out / q_in
}

/// Check if an efficiency claim violates the Carnot limit.
pub fn violates_carnot(claimed_efficiency: f64, t_hot: f64, t_cold: f64) -> bool {
    claimed_efficiency > carnot_efficiency(t_hot, t_cold)
}

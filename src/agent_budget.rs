//! Agent energy budgets — thermodynamic analysis of agent resource consumption.
//!
//! Models an AI agent as a thermodynamic system with energy inputs (compute),
//! entropy production (information disorder), and efficiency metrics.

use serde::{Serialize, Deserialize};
use crate::entropy::information_entropy;
use crate::carnot::carnot_efficiency;

/// Resource consumption record for an agent step.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentStep {
    pub label: String,
    pub compute_joules: f64,     // Energy consumed
    pub tokens_processed: usize, // Throughput
    pub useful_output_bits: f64, // Useful information produced
    pub wall_time_seconds: f64,  // Time taken
}

/// Agent energy budget analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnergyBudget {
    pub steps: Vec<AgentStep>,
    pub total_energy_j: f64,
    pub total_time_s: f64,
    pub entropy_produced: f64,
    pub efficiency: f64,
    pub useful_work_fraction: f64,
}

impl EnergyBudget {
    /// Build an energy budget from a list of steps.
    pub fn from_steps(steps: Vec<AgentStep>) -> Self {
        let total_energy: f64 = steps.iter().map(|s| s.compute_joules).sum();
        let total_time: f64 = steps.iter().map(|s| s.wall_time_seconds).sum();
        let total_bits: f64 = steps.iter().map(|s| s.useful_output_bits).sum();

        // Compute entropy of the step distribution (energy fractions as probabilities)
        let probs: Vec<f64> = if total_energy > 0.0 {
            steps.iter().map(|s| s.compute_joules / total_energy).collect()
        } else {
            vec![]
        };
        let entropy = information_entropy(&probs);

        // Useful work fraction: information output per unit energy
        // Normalized: bits per joule, relative to a reference
        let useful_fraction = if total_energy > 0.0 {
            (total_bits / total_energy).min(1.0)
        } else {
            0.0
        };

        Self {
            steps,
            total_energy_j: total_energy,
            total_time_s: total_time,
            entropy_produced: entropy,
            efficiency: useful_fraction,
            useful_work_fraction: useful_fraction,
        }
    }

    /// Average power consumption in watts.
    pub fn average_power(&self) -> f64 {
        if self.total_time_s > 0.0 {
            self.total_energy_j / self.total_time_s
        } else {
            0.0
        }
    }

    /// Energy per token.
    pub fn energy_per_token(&self) -> f64 {
        let total_tokens: usize = self.steps.iter().map(|s| s.tokens_processed).sum();
        if total_tokens > 0 {
            self.total_energy_j / total_tokens as f64
        } else {
            0.0
        }
    }

    /// Compare agent efficiency to a theoretical Carnot-like bound.
    /// Uses an analogy where T_hot = peak compute capacity, T_cold = idle baseline.
    pub fn carnet_bound_comparison(&self, t_hot_analog: f64, t_cold_analog: f64) -> f64 {
        carnot_efficiency(t_hot_analog, t_cold_analog)
    }

    /// Most expensive step.
    pub fn most_expensive_step(&self) -> Option<&AgentStep> {
        self.steps.iter().max_by(|a, b| a.compute_joules.partial_cmp(&b.compute_joules).unwrap())
    }

    /// Thermodynamic cost of information: Landauer's principle.
    /// Minimum energy to erase one bit: E = k_B · T · ln(2)
    pub fn landauer_cost_per_bit(temperature: f64) -> f64 {
        let k_b = crate::constants::BOLTZMANN;
        k_b * temperature * 2.0_f64.ln()
    }

    /// Total Landauer cost for all information processed.
    pub fn total_landauer_cost(&self, temperature: f64) -> f64 {
        let total_bits: f64 = self.steps.iter().map(|s| s.useful_output_bits).sum();
        total_bits * Self::landauer_cost_per_bit(temperature)
    }
}

/// Simulate an agent as a heat engine: compute is the "hot reservoir",
/// waste (discarded tokens, errors) is the "cold reservoir".
pub fn agent_heat_engine_analysis(
    compute_energy: f64,
    useful_output_energy: f64,
) -> (f64, f64) {
    let waste = compute_energy - useful_output_energy;
    let efficiency = if compute_energy > 0.0 {
        useful_output_energy / compute_energy
    } else {
        0.0
    };
    (efficiency, waste)
}

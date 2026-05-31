//! Entropy computation — Clausius and statistical formulations.

use crate::constants::{R, BOLTZMANN};
use serde::{Serialize, Deserialize};

/// Clausius entropy change for reversible heat transfer: ΔS = Q/T
pub fn clausius_entropy_change(q: f64, t: f64) -> f64 {
    if t <= 0.0 { return f64::NAN; }
    q / t
}

/// Entropy change for isothermal process: ΔS = nRT·ln(V_f/V_i) / T = nR·ln(V_f/V_i)
pub fn entropy_isothermal(n: f64, v_final: f64, v_initial: f64) -> f64 {
    n * R * (v_final / v_initial).ln()
}

/// Entropy change for isobaric heating: ΔS = n·Cp·ln(T_f/T_i)
pub fn entropy_isobaric(n: f64, cp: f64, t_final: f64, t_initial: f64) -> f64 {
    n * cp * (t_final / t_initial).ln()
}

/// Entropy change for isochoric (constant volume) process: ΔS = n·Cv·ln(T_f/T_i)
pub fn entropy_isochoric(n: f64, cv: f64, t_final: f64, t_initial: f64) -> f64 {
    n * cv * (t_final / t_initial).ln()
}

/// Entropy change for adiabatic process (reversible): ΔS = 0
pub fn entropy_adiabatic_reversible() -> f64 {
    0.0
}

/// Entropy of mixing for ideal gases: ΔS_mix = -nR · Σ(x_i · ln(x_i))
pub fn entropy_of_mixing(mole_fractions: &[f64], n_total: f64) -> f64 {
    -n_total * R * mole_fractions.iter()
        .filter(|&&x| x > 0.0)
        .map(|&x| x * x.ln())
        .sum::<f64>()
}

/// Statistical (Boltzmann) entropy: S = k_B · ln(Ω)
pub fn boltzmann_entropy(omega: f64) -> f64 {
    if omega <= 0.0 { return f64::NAN; }
    BOLTZMANN * omega.ln()
}

/// Gibbs entropy formula: S = -k_B · Σ(p_i · ln(p_i))
pub fn gibbs_entropy(probabilities: &[f64]) -> f64 {
    -BOLTZMANN * probabilities.iter()
        .filter(|&&p| p > 0.0)
        .map(|p| p * p.ln())
        .sum::<f64>()
}

/// Entropy change for phase transition: ΔS = ΔH / T_transition
pub fn entropy_of_phase_transition(delta_h: f64, t_transition: f64) -> f64 {
    if t_transition <= 0.0 { return f64::NAN; }
    delta_h / t_transition
}

/// Von Neumann entropy analog (classical): S = -Σ(p_i · ln(p_i)) without k_B.
pub fn information_entropy(probabilities: &[f64]) -> f64 {
    -probabilities.iter()
        .filter(|&&p| p > 0.0)
        .map(|p| p * p.ln())
        .sum::<f64>()
}

/// Entropy as a function of microstate count and temperature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EntropyResult {
    pub clausius: f64,
    pub statistical: f64,
    pub temperature: f64,
}

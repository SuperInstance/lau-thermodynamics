//! The four laws of thermodynamics with mathematical formulations.

use crate::constants::R;

/// Zeroth law: thermal equilibrium is transitive.
/// If A is in thermal equilibrium with B, and B with C, then A is with C.
/// Returns true if all temperatures are equal within tolerance.
pub fn zeroth_law_equilibrium(temps: &[f64], tolerance: f64) -> bool {
    if temps.is_empty() { return true; }
    let ref_t = temps[0];
    temps.iter().all(|&t| (t - ref_t).abs() < tolerance)
}

/// First law: ΔU = Q - W (energy conservation).
/// Returns the change in internal energy given heat added (Q) and work done by system (W).
pub fn first_law(q: f64, w: f64) -> f64 {
    q - w
}

/// Work done by an ideal gas during isobaric (constant pressure) process: W = P·ΔV
pub fn work_isobaric(p: f64, delta_v: f64) -> f64 {
    p * delta_v
}

/// Work done during isothermal expansion of ideal gas: W = nRT·ln(V_f/V_i)
pub fn work_isothermal(n: f64, t: f64, v_final: f64, v_initial: f64) -> f64 {
    n * R * t * (v_final / v_initial).ln()
}

/// Second law: entropy change of the universe must be ≥ 0.
/// Returns total entropy change (system + surroundings).
pub fn second_law_entropy(delta_s_system: f64, delta_s_surroundings: f64) -> f64 {
    delta_s_system + delta_s_surroundings
}

/// Check if a process is spontaneous: ΔS_universe > 0.
pub fn is_spontaneous(delta_s_universe: f64) -> bool {
    delta_s_universe > 0.0
}

/// Third law: entropy of a perfect crystal at absolute zero is zero.
/// Returns 0.0 for T=0 perfect crystal (by convention).
pub fn third_law_entropy(temperature: f64, is_perfect_crystal: bool) -> f64 {
    if temperature == 0.0 && is_perfect_crystal {
        0.0
    } else if temperature == 0.0 {
        f64::NAN // Residual entropy for non-perfect crystals
    } else {
        f64::NAN // Not at absolute zero
    }
}

/// Specific heat capacity at constant volume for an ideal monatomic gas: Cv = (3/2)R per mole.
pub fn cv_monatomic() -> f64 {
    1.5 * R
}

/// Specific heat capacity at constant pressure for an ideal monatomic gas: Cp = (5/2)R per mole.
pub fn cp_monatomic() -> f64 {
    2.5 * R
}

/// Specific heat ratio γ = Cp/Cv for monatomic ideal gas.
pub fn gamma_monatomic() -> f64 {
    5.0 / 3.0
}

/// Cv for diatomic ideal gas: (5/2)R per mole.
pub fn cv_diatomic() -> f64 {
    2.5 * R
}

/// Cp for diatomic ideal gas: (7/2)R per mole.
pub fn cp_diatomic() -> f64 {
    3.5 * R
}

/// γ for diatomic ideal gas: 7/5.
pub fn gamma_diatomic() -> f64 {
    7.0 / 5.0
}

/// Internal energy change for ideal gas: ΔU = n·Cv·ΔT.
pub fn internal_energy_change_ideal(n: f64, cv: f64, delta_t: f64) -> f64 {
    n * cv * delta_t
}

/// Enthalpy change for ideal gas: ΔH = n·Cp·ΔT.
pub fn enthalpy_change_ideal(n: f64, cp: f64, delta_t: f64) -> f64 {
    n * cp * delta_t
}

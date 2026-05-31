//! Statistical mechanics — partition functions, Boltzmann distribution.

use crate::constants::{R, BOLTZMANN, AVOGADRO};
use serde::{Serialize, Deserialize};

/// Boltzmann factor: probability ratio of two states.
/// P(E_i)/P(E_j) = exp(-(E_i - E_j)/(k_B·T))
pub fn boltzmann_factor(energy_i: f64, energy_j: f64, temperature: f64) -> f64 {
    (-(energy_i - energy_j) / (BOLTZMANN * temperature)).exp()
}

/// Boltzmann probability for a single state: P(E_i) = exp(-E_i/(k_B·T)) / Z
pub fn boltzmann_probability(energy: f64, temperature: f64, partition_function: f64) -> f64 {
    (-energy / (BOLTZMANN * temperature)).exp() / partition_function
}

/// Canonical partition function for a set of energy levels: Z = Σ exp(-E_i/(k_B·T))
pub fn canonical_partition_function(energies: &[f64], temperature: f64) -> f64 {
    energies.iter()
        .map(|&e| (-e / (BOLTZMANN * temperature)).exp())
        .sum()
}

/// Translational partition function for a single particle in 3D box:
/// q_trans = (2πmk_BT/h²)^(3/2) · V
/// Using h = 6.626e-34 J·s
pub fn translational_partition_function(mass: f64, temperature: f64, volume: f64) -> f64 {
    let h = 6.62607015e-34;
    ((2.0 * std::f64::consts::PI * mass * BOLTZMANN * temperature) / (h * h)).powf(1.5) * volume
}

/// Rotational partition function for a linear molecule: q_rot = 8π²Ik_BT/(σh²)
pub fn rotational_partition_function_linear(moment_of_inertia: f64, temperature: f64, sigma: f64) -> f64 {
    let h = 6.62607015e-34;
    8.0 * std::f64::consts::PI.powi(2) * moment_of_inertia * BOLTZMANN * temperature
        / (sigma * h * h)
}

/// Vibrational partition function: q_vib = exp(-θ_v/(2T)) / (1 - exp(-θ_v/T))
/// where θ_v = hν/k_B is the characteristic vibrational temperature.
pub fn vibrational_partition_function(theta_v: f64, temperature: f64) -> f64 {
    if temperature <= 0.0 { return 0.0; }
    let x = theta_v / temperature;
    (-x / 2.0).exp() / (1.0 - (-x).exp())
}

/// Internal energy from partition function: U = k_B·T² · d(ln Z)/dT
/// Numerical derivative version.
pub fn internal_energy_from_partition(
    log_z: fn(f64) -> f64,
    temperature: f64,
    dt: f64,
) -> f64 {
    let dlnz_dt = (log_z(temperature + dt / 2.0) - log_z(temperature - dt / 2.0)) / dt;
    BOLTZMANN * temperature * temperature * dlnz_dt
}

/// Helmholtz free energy from partition function: A = -k_B·T·ln(Z)
pub fn helmholtz_from_partition(partition_function: f64, temperature: f64) -> f64 {
    if partition_function <= 0.0 { return f64::NAN; }
    -BOLTZMANN * temperature * partition_function.ln()
}

/// Entropy from partition function: S = k_B·ln(Z) + U/T
pub fn entropy_from_partition(partition_function: f64, internal_energy: f64, temperature: f64) -> f64 {
    if partition_function <= 0.0 || temperature <= 0.0 { return f64::NAN; }
    BOLTZMANN * partition_function.ln() + internal_energy / temperature
}

/// Molar quantities: multiply per-molecule by Avogadro's number.
pub fn per_molecule_to_per_mole(value: f64) -> f64 {
    value * AVOGADRO
}

/// Mean energy from Boltzmann distribution: <E> = Σ(E_i · exp(-E_i/(k_B·T))) / Z
pub fn mean_energy(energies: &[f64], temperature: f64) -> f64 {
    let z = canonical_partition_function(energies, temperature);
    if z == 0.0 { return 0.0; }
    energies.iter()
        .map(|&e| e * (-e / (BOLTZMANN * temperature)).exp())
        .sum::<f64>() / z
}

/// Heat capacity from energy levels: C_v = d<E>/dT (numerical).
pub fn heat_capacity_from_energies(energies: &[f64], temperature: f64, dt: f64) -> f64 {
    let e1 = mean_energy(energies, temperature - dt / 2.0);
    let e2 = mean_energy(energies, temperature + dt / 2.0);
    (e2 - e1) / dt
}

/// Two-level system partition function and properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TwoLevelSystem {
    pub energy_ground: f64,
    pub energy_excited: f64,
}

impl TwoLevelSystem {
    pub fn new(ground: f64, excited: f64) -> Self {
        Self { energy_ground: ground, energy_excited: excited }
    }

    pub fn partition_function(&self, temperature: f64) -> f64 {
        canonical_partition_function(
            &[self.energy_ground, self.energy_excited],
            temperature,
        )
    }

    pub fn fraction_excited(&self, temperature: f64) -> f64 {
        let z = self.partition_function(temperature);
        (-self.energy_excited / (BOLTZMANN * temperature)).exp() / z
    }

    pub fn fraction_ground(&self, temperature: f64) -> f64 {
        1.0 - self.fraction_excited(temperature)
    }

    /// At high T, both fractions approach 0.5.
    pub fn high_temp_limit_fraction(&self) -> f64 {
        0.5
    }
}

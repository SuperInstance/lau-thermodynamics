# lau-thermodynamics

A Rust crate implementing classical and statistical thermodynamics from first principles. Covers the four laws, ideal and van der Waals gases, Carnot cycles, entropy (Clausius, Boltzmann, Gibbs), Maxwell relations, phase transitions, partition functions, heat transfer (conduction, convection, radiation), and an agent energy-budget model with Landauer's principle.

---

## What This Does

| Module | What It Gives You |
|---|---|
| `constants` | R, k_B, N_A, σ, atm, temperature conversions |
| `laws` | All four laws — zeroth (equilibrium transitivity), first (ΔU = Q − W), second (ΔS_universe ≥ 0), third (S → 0 at 0 K) — plus heat capacities for monatomic & diatomic ideal gases |
| `gas` | Ideal gas law (PV = nRT) in all directions, van der Waals equation with built-in parameters for N₂, O₂, CO₂, H₂O, He, compressibility factor Z, Boyle temperature |
| `carnot` | Carnot efficiency, work, heat rejected, COP for refrigerators & heat pumps, full `CarnotCycle` struct, Carnot-limit violation checker |
| `entropy` | Clausius ΔS = Q/T, isothermal/isobaric/isochoric/adiabatic entropy changes, entropy of mixing, Boltzmann S = k_B ln Ω, Gibbs S = −k_B Σ pᵢ ln pᵢ, phase-transition entropy, information entropy |
| `maxwell` | All four Maxwell relations, ideal-gas verification, thermodynamic potential derivatives, isothermal compressibility, thermal expansion coefficient |
| `phase` | Clausius–Clapeyron equation, boiling point vs pressure, triple point & critical point data for water, Gibbs phase rule, phase equilibrium verification |
| `statistical` | Boltzmann factor, canonical partition function Z, internal energy from Z, average energy, Helmholtz free energy, heat capacity from Z, Fermi–Dirac & Bose–Einstein distributions |
| `heat_transfer` | Fourier conduction, thermal resistance (series & parallel), Newton's cooling, Stefan–Boltzmann radiation, Biot & Fourier numbers |
| `agent_budget` | `AgentStep` / `EnergyBudget` for tracking agent compute energy per step, Landauer cost per bit, heat-engine analogy for agent efficiency |

All types derive `Serialize`/`Deserialize`. **76 tests** verify every formula against known values.

---

## Key Idea

Thermodynamics is the universal accounting system for energy and information. This crate makes every major result computable:

- **First Law** → energy balance for any process
- **Second Law** → spontaneity check (ΔS_universe > 0?)
- **Carnot limit** → maximum theoretical efficiency
- **Landauer's principle** → minimum energy to erase one bit (k_B T ln 2)

The `agent_budget` module applies these to AI agents: track compute energy, measure useful information output, compare to thermodynamic bounds.

---

## Install

```toml
[dependencies]
lau-thermodynamics = "0.1"
```

Requires Rust **2021 edition**. Depends on `serde` and `nalgebra`.

---

## Quick Start

### The four laws

```rust
use lau_thermodynamics::*;

// Zeroth law: are these systems in thermal equilibrium?
assert!(zeroth_law_equilibrium(&[300.0, 300.0, 300.0], 0.1));

// First law: ΔU = Q - W
let delta_u = first_law(1000.0, 400.0); // 600 J

// Second law: is this process spontaneous?
assert!(is_spontaneous(5.0));  // ΔS_universe > 0

// Third law: perfect crystal at absolute zero
assert_eq!(third_law_entropy(0.0, true), 0.0);
```

### Ideal gas

```rust
use lau_thermodynamics::*;

// PV = nRT
let p = ideal_gas_pressure(1.0, 300.0, 0.0249);       // Pa
let v = ideal_gas_volume(1.0, 300.0, 101325.0);        // m³
let t = ideal_gas_temperature(101325.0, 0.0249, 1.0);  // K
let n = ideal_gas_moles(101325.0, 0.0249, 300.0);      // mol
```

### Van der Waals (real gas)

```rust
use lau_thermodynamics::*;

let params = VanDerWaalsParams::nitrogen();
let p_real = vdw_pressure(1.0, 300.0, 0.0249, &params);
let v_solved = vdw_volume(1.0, 300.0, 101325.0, &params, 100);
let z = compressibility_factor(101325.0, 0.0249, 1.0, 300.0); // ≈ 1 for ideal
let t_boyle = boyle_temperature(&params); // ~432 K for N₂
```

### Carnot cycle

```rust
use lau_thermodynamics::*;

let cycle = CarnotCycle::new(600.0, 300.0, 1000.0);
assert_eq!(cycle.efficiency, 0.5);     // η = 1 - 300/600
assert_eq!(cycle.work, 500.0);         // W = Q·η
assert_eq!(cycle.q_cold, 500.0);       // Q_cold = Q - W

let cop_ref = cop_refrigerator(300.0, 270.0); // 9.0
let cop_hp  = cop_heat_pump(300.0, 270.0);    // 10.0

assert!(violates_carnot(0.6, 600.0, 300.0)); // claimed > theoretical
```

### Entropy

```rust
use lau_thermodynamics::*;

let ds = clausius_entropy_change(1000.0, 300.0);           // Q/T
let ds_iso = entropy_isothermal(1.0, 2.0, 1.0);           // nR ln(Vf/Vi)
let ds_mix = entropy_of_mixing(&[0.5, 0.5], 1.0);         // nR ln(2)
let s_boltz = boltzmann_entropy(1e10);                     // k_B ln(Ω)
let s_gibbs = gibbs_entropy(&[0.25, 0.25, 0.25, 0.25]);   // k_B ln(4)
```

### Statistical mechanics

```rust
use lau_thermodynamics::*;

let energies = &[0.0, 1e-21, 2e-21];
let z = canonical_partition_function(energies, 300.0);
let u = partition_function_energy(energies, 300.0);   // ⟨E⟩
let a = helmholtz_from_partition(z, 300.0);           // A = −k_B T ln Z
let cv = heat_capacity_from_partition(energies, 300.0);

// Fermi-Dirac: probability of occupation at energy E
let fd = fermi_dirac(1e-21, 1e-21, 300.0); // E, μ, T

// Bose-Einstein
let be = bose_einstein(1e-22, 300.0); // E, T
```

### Agent energy budget (Landauer's principle)

```rust
use lau_thermodynamics::*;

let steps = vec![
    AgentStep { label: "inference".into(), compute_joules: 100.0,
                tokens_processed: 1000, useful_output_bits: 500.0, wall_time_seconds: 1.0 },
    AgentStep { label: "parsing".into(), compute_joules: 50.0,
                tokens_processed: 500, useful_output_bits: 200.0, wall_time_seconds: 0.5 },
];
let budget = EnergyBudget::from_steps(steps);
println!("Total energy: {} J", budget.total_energy_j);
println!("Average power: {} W", budget.average_power());
println!("Energy/token:  {} J", budget.energy_per_token());

// Minimum energy to erase one bit at room temperature
let landauer = EnergyBudget::landauer_cost_per_bit(300.0); // ~2.87e-21 J
```

---

## API Reference

### `constants`

| Constant | Value | Unit |
|---|---|---|
| `R` | 8.314462618 | J/(mol·K) |
| `BOLTZMANN` | 1.380649×10⁻²³ | J/K |
| `AVOGADRO` | 6.02214076×10²³ | /mol |
| `STEFAN_BOLTZMANN` | 5.670374419×10⁻⁸ | W/(m²·K⁴) |
| `ATM` | 101325 | Pa |
| `C_TO_K` | 273.15 | K |

Functions: `celsius_to_kelvin`, `kelvin_to_celsius`, `gas_constant`.

### `laws`

| Function | Formula |
|---|---|
| `zeroth_law_equilibrium(temps, tol)` | All temps equal within tolerance |
| `first_law(q, w)` | ΔU = Q − W |
| `work_isobaric(p, ΔV)` | W = PΔV |
| `work_isothermal(n, T, Vf, Vi)` | W = nRT ln(Vf/Vi) |
| `second_law_entropy(ΔS_sys, ΔS_surr)` | ΔS_universe |
| `is_spontaneous(ΔS_universe)` | ΔS_universe > 0 |
| `third_law_entropy(T, perfect)` | 0 if T=0 & perfect crystal |
| `cv_monatomic`, `cp_monatomic`, `gamma_monatomic` | 3R/2, 5R/2, 5/3 |
| `cv_diatomic`, `cp_diatomic`, `gamma_diatomic` | 5R/2, 7R/2, 7/5 |
| `internal_energy_change_ideal(n, cv, ΔT)` | ΔU = nCvΔT |
| `enthalpy_change_ideal(n, cp, ΔT)` | ΔH = nCpΔT |

### `gas`

| Function | Description |
|---|---|
| `ideal_gas_pressure(n, T, V)` | P = nRT/V |
| `ideal_gas_volume(n, T, P)` | V = nRT/P |
| `ideal_gas_temperature(P, V, n)` | T = PV/(nR) |
| `ideal_gas_moles(P, V, T)` | n = PV/(RT) |
| `vdw_pressure(n, T, V, params)` | P = nRT/(V−nb) − a(n/V)² |
| `vdw_volume(n, T, P, params, iters)` | Newton's method solve |
| `compressibility_factor(P, V, n, T)` | Z = PV/(nRT) |
| `boyle_temperature(params)` | T_B = a/(Rb) |

`VanDerWaalsParams` — built-in: `nitrogen`, `oxygen`, `carbon_dioxide`, `water`, `helium`.

### `carnot`

| Function / Type | Description |
|---|---|
| `carnot_efficiency(T_h, T_c)` | η = 1 − T_c/T_h |
| `carnot_work(Q_h, T_h, T_c)` | W = Q_h · η |
| `carnot_heat_rejected(Q_h, T_h, T_c)` | Q_c = Q_h − W |
| `cop_refrigerator(T_h, T_c)` | T_c / (T_h − T_c) |
| `cop_heat_pump(T_h, T_c)` | T_h / (T_h − T_c) |
| `CarnotCycle::new(T_h, T_c, Q_h)` | Full cycle struct |
| `thermal_efficiency(Q_in, W_out)` | W/Q |
| `violates_carnot(η_claim, T_h, T_c)` | η_claim > η_Carnot? |

### `entropy`

| Function | Formula |
|---|---|
| `clausius_entropy_change(Q, T)` | ΔS = Q/T |
| `entropy_isothermal(n, Vf, Vi)` | nR ln(Vf/Vi) |
| `entropy_isobaric(n, Cp, Tf, Ti)` | nCp ln(Tf/Ti) |
| `entropy_isochoric(n, Cv, Tf, Ti)` | nCv ln(Tf/Ti) |
| `entropy_adiabatic_reversible()` | 0 |
| `entropy_of_mixing(x[], n)` | −nR Σ xᵢ ln xᵢ |
| `boltzmann_entropy(Ω)` | k_B ln Ω |
| `gibbs_entropy(p[])` | −k_B Σ pᵢ ln pᵢ |
| `entropy_of_phase_transition(ΔH, T)` | ΔH/T |
| `information_entropy(p[])` | −Σ pᵢ ln pᵢ (no k_B) |

### `maxwell`

| Function | Maxwell Relation |
|---|---|
| `maxwell_1(a, b)` | (∂T/∂V)_S = −(∂P/∂S)_V |
| `maxwell_2(a, b)` | (∂T/∂P)_S = (∂V/∂S)_P |
| `maxwell_3(a, b)` | (∂P/∂T)_V = (∂S/∂V)_T |
| `maxwell_4(a, b)` | (∂V/∂T)_P = −(∂S/∂P)_T |

Also: `ideal_gas_maxwell_3_verify`, `ideal_gas_maxwell_4_verify`, `internal_energy_derivatives`, `gibbs_derivatives`, `isothermal_compressibility`, `thermal_expansion_coefficient`.

### `phase`

| Function / Type | Description |
|---|---|
| `clausius_clapeyron_slope(ΔH, T, ΔV)` | dP/dT = ΔH / (TΔV) |
| `clausius_clapeyron_integrated(P₁, T₁, T₂, ΔH)` | Integrated form |
| `boiling_point_at_pressure(P₁, T₁, P₂, ΔH)` | Inverse Clapeyron |
| `PhaseTransition::water_boiling()` | 100°C, 101325 Pa, ΔH = 40.7 kJ/mol |
| `TriplePoint::water()` | 273.16 K, 611.657 Pa |
| `CriticalPoint::water()` | 647.1 K, 22.064 MPa |
| `gibbs_phase_rule(components, phases)` | F = C − P + 2 |

### `statistical`

| Function | Formula |
|---|---|
| `boltzmann_factor(E, E_ref, T)` | exp(−(E − E_ref) / k_BT) |
| `canonical_partition_function(energies, T)` | Z = Σ exp(−Eᵢ / k_BT) |
| `partition_function_energy(energies, T)` | ⟨E⟩ = Σ Eᵢ exp(−Eᵢ/k_BT) / Z |
| `average_energy(energies, T)` | Same as above |
| `helmholtz_from_partition(Z, T)` | A = −k_BT ln Z |
| `heat_capacity_from_partition(energies, T)` | C_v = d⟨E⟩/dT |
| `fermi_dirac(E, μ, T)` | 1 / (exp((E−μ)/k_BT) + 1) |
| `bose_einstein(E, T)` | 1 / (exp(E/k_BT) − 1) |

### `heat_transfer`

| Function | Formula |
|---|---|
| `conduction_heat(k, A, T_h, T_c, L)` | kA(T_h − T_c)/L |
| `thermal_resistance_conduction(L, k, A)` | L/(kA) |
| `thermal_resistance_series(R[])` | Σ Rᵢ |
| `thermal_resistance_parallel(R[])` | (Σ 1/Rᵢ)⁻¹ |
| `convection_heat(h, A, T_s, T_f)` | hA(T_s − T_f) |
| `radiation_heat(ε, A, T_s, T_env)` | εσA(T_s⁴ − T_env⁴) |
| `biot_number(h, L, k)` | hL/k |
| `fourier_number(α, t, L)` | αt/L² |

### `agent_budget`

| Type | Description |
|---|---|
| `AgentStep` | Single step: label, energy (J), tokens, useful bits, wall time |
| `EnergyBudget` | Aggregated: total energy, time, entropy, efficiency, power, energy/token |
| `EnergyBudget::landauer_cost_per_bit(T)` | k_BT ln 2 — minimum energy to erase one bit |
| `EnergyBudget::total_landauer_cost(T)` | Sum over all useful bits |
| `agent_heat_engine_analysis(E_compute, E_useful)` | (η, waste) — treat agent as heat engine |

### `ThermoState`

Core state point with temperature, pressure, volume, internal energy, entropy, enthalpy, and Gibbs free energy. Class methods:

- `gibbs_free_energy(H, T, S)` → G = H − TS
- `helmholtz_free_energy(U, T, S)` → A = U − TS

---

## How It Works

### Laws as Functions
Each law is a pure function. The first law `first_law(q, w)` returns ΔU = Q − W. The second law's spontaneity check is `is_spontaneous(delta_s_universe)` which returns `ΔS > 0`. This makes the laws composable — chain them into process analysis pipelines.

### Ideal vs Real Gas
The ideal gas law PV = nRT is exact for point particles. Real gases deviate due to intermolecular forces and finite molecular volume. The van der Waals equation corrects both:

```
P = nRT/(V − nb) − a(n/V)²
```

Parameter `a` accounts for attraction (lowers pressure), `b` for excluded volume (increases effective V). The crate includes parameters for five common gases and solves for volume via Newton's method.

### Carnot Limit
The Carnot efficiency η = 1 − T_cold/T_hot is the absolute thermodynamic ceiling for any heat engine operating between two reservoirs. The crate computes this, plus COP for refrigerators and heat pumps, and checks whether claimed efficiencies violate the limit.

### Entropy: Three Formulations
- **Clausius**: ΔS = ∫ dQ/T (macroscopic, path-dependent)
- **Boltzmann**: S = k_B ln Ω (microscopic, counting microstates)
- **Gibbs**: S = −k_B Σ pᵢ ln pᵢ (most general, works for non-uniform distributions)

All three are consistent: uniform distributions make Gibbs → Boltzmann, and equilibrium processes make Clausius agree with both.

### Maxwell Relations
The four Maxwell relations are cross-derivative identities that follow from exact differentials of thermodynamic potentials. For example, from dU = TdS − PdV:

```
(∂T/∂V)_S = −(∂P/∂S)_V
```

The crate numerically verifies these for the ideal gas.

### Phase Transitions
The Clausius–Clapeyron equation dP/dT = ΔH/(TΔV) governs coexistence curves. The crate integrates it to find boiling points at arbitrary pressures, with water as a worked example.

### Statistical Mechanics
The canonical partition function Z = Σ exp(−Eᵢ/k_BT) encodes all thermodynamic information. Internal energy, Helmholtz free energy, and heat capacity are all derivatives of Z. The crate also provides Fermi–Dirac and Bose–Einstein distributions for quantum statistics.

### Agent Energy Budgets
The `agent_budget` module treats an AI agent as a thermodynamic system:
- **Compute energy** = heat input (hot reservoir)
- **Useful information** = work output
- **Discarded tokens/errors** = waste heat (cold reservoir)

Landauer's principle sets the absolute minimum: erasing one bit costs at least k_BT ln 2 ≈ 2.87 × 10⁻²¹ J at room temperature. The crate computes this and compares agent efficiency to Carnot-like bounds.

---

## The Math

**First law (energy conservation):**
$$\Delta U = Q - W$$

**Ideal gas law:**
$$PV = nRT$$

**Van der Waals equation:**
$$\left(P + \frac{an^2}{V^2}\right)(V - nb) = nRT$$

**Carnot efficiency:**
$$\eta = 1 - \frac{T_\text{cold}}{T_\text{hot}}$$

**Clausius entropy:**
$$\Delta S = \int \frac{\delta Q}{T}$$

**Boltzmann entropy:**
$$S = k_B \ln \Omega$$

**Gibbs entropy:**
$$S = -k_B \sum_i p_i \ln p_i$$

**Clausius–Clapeyron equation:**
$$\frac{dP}{dT} = \frac{\Delta H}{T \Delta V}$$

**Canonical partition function:**
$$Z = \sum_i e^{-E_i / k_B T}, \qquad \langle E \rangle = -\frac{\partial \ln Z}{\partial \beta}, \qquad A = -k_B T \ln Z$$

**Fermi–Dirac distribution:**
$$f(E) = \frac{1}{e^{(E-\mu)/k_BT} + 1}$$

**Bose–Einstein distribution:**
$$n(E) = \frac{1}{e^{E/k_BT} - 1}$$

**Landauer's principle:**
$$E_\text{min} = k_B T \ln 2 \quad \text{(per bit erased)}$$

**Fourier's law:**
$$\dot{Q} = -k A \frac{dT}{dx}$$

**Stefan–Boltzmann radiation:**
$$\dot{Q} = \varepsilon \sigma A (T_s^4 - T_\text{env}^4)$$

**Gibbs phase rule:**
$$F = C - P + 2$$

---

## Tests

**76 tests** covering:

- Physical constant accuracy (R, k_B, conversions)
- All four laws with known values
- Ideal gas law round-trips (P→V→T→n→P)
- Van der Waals: pressure < ideal for N₂, He ≈ ideal, Newton solve convergence
- Compressibility factor ≈ 1 at STP
- Boyle temperature for nitrogen (~432 K)
- Carnot efficiency, work, heat rejected, COP, violation check
- Entropy: Clausius, isothermal, isobaric, isochoric, adiabatic, mixing, Boltzmann, Gibbs, phase transition
- Maxwell relations: identity verification, ideal-gas numerical checks
- Clausius–Clapeyron slope, integrated form, boiling point prediction
- Triple point & critical point for water
- Gibbs phase rule
- Boltzmann factor, partition function, energy, Helmholtz, heat capacity
- Fermi–Dirac and Bose–Einstein distributions
- Conduction, convection, radiation heat transfer
- Thermal resistance (series & parallel)
- Biot and Fourier numbers
- Agent energy budget: totals, power, energy/token, most expensive step
- Landauer cost, heat-engine analysis
- `ThermoState` Gibbs & Helmholtz free energy

```bash
cargo test
```

---

## License

MIT

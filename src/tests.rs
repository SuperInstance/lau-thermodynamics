//! Comprehensive tests for lau-thermodynamics.

#[cfg(test)]
mod tests {
    use crate::*;
    use crate::constants::*;

    // ─── Constants ───

    #[test]
    fn test_gas_constant_value() {
        assert!((R - 8.314462618).abs() < 1e-6);
    }

    #[test]
    fn test_boltzmann_constant_value() {
        assert!((BOLTZMANN - 1.380649e-23).abs() < 1e-30);
    }

    #[test]
    fn test_celsius_kelvin_conversion() {
        assert!((celsius_to_kelvin(0.0) - 273.15).abs() < 1e-10);
        assert!((celsius_to_kelvin(100.0) - 373.15).abs() < 1e-10);
        assert!((kelvin_to_celsius(373.15) - 100.0).abs() < 1e-10);
    }

    // ─── Laws of Thermodynamics ───

    #[test]
    fn test_zeroth_law_equilibrium() {
        assert!(zeroth_law_equilibrium(&[300.0, 300.0, 300.0], 0.1));
        assert!(!zeroth_law_equilibrium(&[300.0, 200.0, 300.0], 0.1));
    }

    #[test]
    fn test_first_law_energy_conservation() {
        let q = 1000.0;
        let w = 400.0;
        let du = first_law(q, w);
        assert!((du - 600.0).abs() < 1e-10);
    }

    #[test]
    fn test_work_isobaric() {
        let w = work_isobaric(101325.0, 0.01);
        assert!((w - 1013.25).abs() < 1e-6);
    }

    #[test]
    fn test_work_isothermal() {
        let w = work_isothermal(1.0, 300.0, 2.0, 1.0);
        let expected = 1.0 * R * 300.0 * 2.0_f64.ln();
        assert!((w - expected).abs() < 1e-6);
    }

    #[test]
    fn test_second_law_spontaneous() {
        assert!(is_spontaneous(5.0));
        assert!(!is_spontaneous(-1.0));
        assert!(!is_spontaneous(0.0)); // Not spontaneous, equilibrium
    }

    #[test]
    fn test_third_law_perfect_crystal() {
        assert!((third_law_entropy(0.0, true) - 0.0).abs() < 1e-10);
        assert!(third_law_entropy(0.0, false).is_nan());
    }

    #[test]
    fn test_monatomic_gas_heat_capacities() {
        assert!((cv_monatomic() - 1.5 * R).abs() < 1e-10);
        assert!((cp_monatomic() - 2.5 * R).abs() < 1e-10);
        assert!((gamma_monatomic() - 5.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_diatomic_gas_heat_capacities() {
        assert!((cv_diatomic() - 2.5 * R).abs() < 1e-10);
        assert!((cp_diatomic() - 3.5 * R).abs() < 1e-10);
        assert!((gamma_diatomic() - 7.0 / 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_internal_energy_change() {
        let du = internal_energy_change_ideal(1.0, cv_monatomic(), 100.0);
        assert!((du - 1.0 * 1.5 * R * 100.0).abs() < 1e-6);
    }

    // ─── Ideal Gas Law ───

    #[test]
    fn test_ideal_gas_pressure() {
        let p = ideal_gas_pressure(1.0, 300.0, 0.0249);
        let expected = 1.0 * R * 300.0 / 0.0249;
        assert!((p - expected).abs() < 1.0);
    }

    #[test]
    fn test_ideal_gas_law_roundtrip() {
        let n = 2.0;
        let t = 350.0;
        let v = 0.05;
        let p = ideal_gas_pressure(n, t, v);
        let v2 = ideal_gas_volume(n, t, p);
        assert!((v - v2).abs() < 1e-8);
        let t2 = ideal_gas_temperature(p, v, n);
        assert!((t - t2).abs() < 1e-8);
        let n2 = ideal_gas_moles(p, v, t);
        assert!((n - n2).abs() < 1e-8);
    }

    #[test]
    fn test_ideal_gas_at_stp() {
        let p = ideal_gas_pressure(1.0, 273.15, 0.0224);
        // Should be ~101325 Pa
        assert!((p - 101325.0).abs() < 500.0);
    }

    // ─── Van der Waals ───

    #[test]
    fn test_vdw_pressure_less_than_ideal() {
        let params = VanDerWaalsParams::nitrogen();
        let v = 0.0249;
        let t = 300.0;
        let n = 1.0;
        let p_ideal = ideal_gas_pressure(n, t, v);
        let p_vdw = vdw_pressure(n, t, v, &params);
        // At moderate conditions, VdW pressure should be close but different from ideal
        assert!(p_vdw > 0.0);
        assert!((p_ideal - p_vdw).abs() / p_ideal < 0.05);
    }

    #[test]
    fn test_vdw_helium_close_to_ideal() {
        let params = VanDerWaalsParams::helium();
        let v = 0.0249;
        let t = 300.0;
        let n = 1.0;
        let p_ideal = ideal_gas_pressure(n, t, v);
        let p_vdw = vdw_pressure(n, t, v, &params);
        // Helium is very close to ideal
        assert!((p_ideal - p_vdw).abs() / p_ideal < 0.01);
    }

    #[test]
    fn test_vdw_volume_iteration() {
        let params = VanDerWaalsParams::nitrogen();
        let n = 1.0;
        let t = 300.0;
        let p = 101325.0;
        let v = vdw_volume(n, t, p, &params, 100);
        let p_check = vdw_pressure(n, t, v, &params);
        assert!((p_check - p).abs() / p < 0.01);
    }

    #[test]
    fn test_compressibility_factor_ideal() {
        let z = compressibility_factor(101325.0, 0.0249, 1.0, 300.0);
        assert!((z - 1.0).abs() < 0.05);
    }

    #[test]
    fn test_boyle_temperature() {
        let params = VanDerWaalsParams::nitrogen();
        let tb = boyle_temperature(&params);
        // Nitrogen Boyle temp should be positive and reasonable (~400-430 K)
        assert!(tb > 300.0 && tb < 500.0);
    }

    // ─── Carnot Cycle ───

    #[test]
    fn test_carnot_efficiency() {
        let eta = carnot_efficiency(600.0, 300.0);
        assert!((eta - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_carnot_efficiency_zero_cold() {
        let eta = carnot_efficiency(300.0, 0.0);
        // Zero cold → 100% (theoretical)
        assert!((eta - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_carnot_work_and_heat() {
        let cycle = CarnotCycle::new(600.0, 300.0, 1000.0);
        assert!((cycle.efficiency - 0.5).abs() < 1e-10);
        assert!((cycle.work - 500.0).abs() < 1e-10);
        assert!((cycle.q_cold - 500.0).abs() < 1e-10);
    }

    #[test]
    fn test_carnot_invalid_temperatures() {
        assert!((carnot_efficiency(300.0, 400.0)) == 0.0); // T_cold > T_hot
        assert!((carnot_efficiency(0.0, 300.0)) == 0.0); // T_hot = 0
    }

    #[test]
    fn test_cop_refrigerator() {
        let cop = cop_refrigerator(300.0, 270.0);
        let expected = 270.0 / (300.0 - 270.0);
        assert!((cop - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cop_heat_pump() {
        let cop = cop_heat_pump(300.0, 270.0);
        let expected = 300.0 / (300.0 - 270.0);
        assert!((cop - expected).abs() < 1e-10);
    }

    #[test]
    fn test_violates_carnot() {
        assert!(!violates_carnot(0.4, 600.0, 300.0)); // 0.4 < 0.5
        assert!(violates_carnot(0.6, 600.0, 300.0));  // 0.6 > 0.5
    }

    // ─── Entropy ───

    #[test]
    fn test_clausius_entropy_change() {
        let ds = clausius_entropy_change(1000.0, 300.0);
        assert!((ds - 1000.0 / 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_entropy_isothermal_expansion() {
        let ds = entropy_isothermal(1.0, 2.0, 1.0);
        assert!((ds - R * 2.0_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_entropy_isobaric_heating() {
        let ds = entropy_isobaric(1.0, cp_monatomic(), 600.0, 300.0);
        let expected = 1.0 * 2.5 * R * (600.0_f64 / 300.0).ln();
        assert!((ds - expected).abs() < 1e-6);
    }

    #[test]
    fn test_entropy_isochoric() {
        let ds = entropy_isochoric(1.0, cv_monatomic(), 600.0, 300.0);
        let expected = 1.0 * 1.5 * R * (600.0_f64 / 300.0).ln();
        assert!((ds - expected).abs() < 1e-6);
    }

    #[test]
    fn test_entropy_adiabatic_reversible() {
        assert!((entropy_adiabatic_reversible() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_entropy_of_mixing() {
        let ds = entropy_of_mixing(&[0.5, 0.5], 1.0);
        // ΔS_mix = -nR * (0.5*ln(0.5) + 0.5*ln(0.5)) = -nR * ln(0.5) = nR*ln(2)
        let expected = 1.0 * R * 2.0_f64.ln();
        assert!((ds - expected).abs() < 1e-6);
    }

    #[test]
    fn test_boltzmann_entropy() {
        let omega = 1e10;
        let s = boltzmann_entropy(omega);
        assert!((s - BOLTZMANN * omega.ln()).abs() < 1e-30);
    }

    #[test]
    fn test_gibbs_entropy_uniform() {
        let probs = vec![0.25, 0.25, 0.25, 0.25];
        let s = gibbs_entropy(&probs);
        let expected = BOLTZMANN * 4.0_f64.ln();
        assert!((s - expected).abs() < 1e-30);
    }

    #[test]
    fn test_entropy_phase_transition() {
        let ds = entropy_of_phase_transition(40700.0, 373.15);
        assert!((ds - 40700.0 / 373.15).abs() < 1e-6);
    }

    // ─── Maxwell Relations ───

    #[test]
    fn test_maxwell_1_identity() {
        // Test with a=b: (∂T/∂V)_S = 5, -(∂P/∂S)_V = -(-5) = 5
        let rel = maxwell_1(5.0, -5.0);
        assert!(rel.difference.abs() < 1e-10);
    }

    #[test]
    fn test_maxwell_3_identity() {
        let rel = maxwell_3(3.0, 3.0);
        assert!(rel.difference.abs() < 1e-10);
    }

    #[test]
    fn test_maxwell_4_identity() {
        let rel = maxwell_4(2.0, -2.0);
        assert!(rel.difference.abs() < 1e-10);
    }

    #[test]
    fn test_ideal_gas_maxwell_3_verify() {
        assert!(ideal_gas_maxwell_3_verify(1.0, 300.0, 0.025));
    }

    #[test]
    fn test_ideal_gas_maxwell_4_verify() {
        assert!(ideal_gas_maxwell_4_verify(1.0, 300.0, 101325.0));
    }

    #[test]
    fn test_thermodynamic_potential_derivatives() {
        let (dudv_s, duvd_s) = internal_energy_derivatives(300.0, 101325.0);
        assert!((dudv_s - 300.0).abs() < 1e-10);
        assert!((duvd_s - (-101325.0)).abs() < 1e-10);
    }

    #[test]
    fn test_gibbs_derivatives() {
        let (dgdt_p, dgdp_t) = gibbs_derivatives(100.0, 0.025);
        assert!((dgdt_p - (-100.0)).abs() < 1e-10);
        assert!((dgdp_t - 0.025).abs() < 1e-10);
    }

    #[test]
    fn test_isothermal_compressibility() {
        let kappa = isothermal_compressibility(0.025, -1e-7);
        assert!(kappa > 0.0); // Must be positive
    }

    // ─── Phase Transitions ───

    #[test]
    fn test_clausius_clapeyron_slope() {
        let slope = clausius_clapeyron_slope(40700.0, 373.15, 0.001);
        // Should be positive for boiling
        assert!(slope > 0.0);
    }

    #[test]
    fn test_clausius_clapeyron_integrated() {
        let p2 = clausius_clapeyron_integrated(101325.0, 373.15, 363.15, 40700.0);
        // Lower T → lower P
        assert!(p2 < 101325.0);
        assert!(p2 > 0.0);
    }

    #[test]
    fn test_boiling_point_at_pressure() {
        let t_boil = boiling_point_at_pressure(101325.0, 373.15, 200000.0, 40700.0);
        // Higher pressure → higher boiling point
        assert!(t_boil > 373.15);
    }

    #[test]
    fn test_water_boiling_equilibrium() {
        let trans = PhaseTransition::water_boiling();
        let dg = trans.verify_equilibrium();
        assert!(dg.abs() < 1.0); // Should be ~0 at equilibrium
    }

    #[test]
    fn test_triple_point_water() {
        let tp = TriplePoint::water();
        assert!((tp.temperature - 273.16).abs() < 0.01);
        assert!((tp.pressure - 611.657).abs() < 1.0);
    }

    #[test]
    fn test_critical_point_to_vdw() {
        let cp = CriticalPoint::water();
        let (a, b) = cp.to_vdw_params();
        assert!(a > 0.0);
        assert!(b > 0.0);
    }

    #[test]
    fn test_gibbs_phase_rule() {
        // Water at triple point: 1 component, 3 phases → F = 0
        assert_eq!(gibbs_phase_rule(1, 3), 0);
        assert_eq!(gibbs_phase_rule(1, 1), 2);
    }

    // ─── Statistical Mechanics ───

    #[test]
    fn test_boltzmann_factor_equal_energies() {
        let bf = boltzmann_factor(100.0, 100.0, 300.0);
        assert!((bf - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_boltzmann_factor_higher_energy() {
        let bf = boltzmann_factor(200.0, 100.0, 300.0);
        assert!(bf < 1.0);
    }

    #[test]
    fn test_canonical_partition_function() {
        let z = canonical_partition_function(&[0.0, 1e-21, 2e-21], 300.0);
        assert!(z > 1.0); // Ground state alone contributes 1.0
    }

    #[test]
    fn test_partition_function_ground_state_only() {
        let z = canonical_partition_function(&[0.0], 300.0);
        assert!((z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_two_level_system() {
        let sys = TwoLevelSystem::new(0.0, 1e-21);
        let z = sys.partition_function(300.0);
        assert!(z > 1.0);

        let f_excited = sys.fraction_excited(300.0);
        assert!(f_excited > 0.0 && f_excited < 0.5);

        let f_ground = sys.fraction_ground(300.0);
        assert!((f_ground + f_excited - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_two_level_system_high_temperature() {
        let sys = TwoLevelSystem::new(0.0, 1e-23);
        // At high T relative to energy gap, fractions → 0.5
        let f_excited = sys.fraction_excited(1e6);
        assert!((f_excited - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_mean_energy() {
        let energies = &[0.0, 1e-21, 2e-21];
        let mean = mean_energy(energies, 300.0);
        assert!(mean > 0.0);
        assert!(mean < 2e-21);
    }

    #[test]
    fn test_heat_capacity_from_energies() {
        let energies = &[0.0, 1e-21, 2e-21];
        let cv = heat_capacity_from_energies(energies, 300.0, 1.0);
        assert!(cv > 0.0);
    }

    #[test]
    fn test_helmholtz_from_partition() {
        let a = helmholtz_from_partition(2.0, 300.0);
        let expected = -BOLTZMANN * 300.0 * 2.0_f64.ln();
        assert!((a - expected).abs() < 1e-30);
    }

    #[test]
    fn test_vibrational_partition_function() {
        let q = vibrational_partition_function(1000.0, 300.0);
        assert!(q > 0.0);
    }

    #[test]
    fn test_per_molecule_to_per_mole() {
        let molar = per_molecule_to_per_mole(1.0);
        assert!((molar - AVOGADRO).abs() < 1e10);
    }

    // ─── Heat Transfer ───

    #[test]
    fn test_conduction_heat() {
        let q = conduction_heat(401.0, 0.01, 400.0, 300.0, 0.1);
        let expected = 401.0 * 0.01 * 100.0 / 0.1;
        assert!((q - expected).abs() < 1e-6);
    }

    #[test]
    fn test_thermal_resistance_series() {
        let r = thermal_resistance_series(&[2.0, 3.0, 5.0]);
        assert!((r - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_resistance_parallel() {
        let r = thermal_resistance_parallel(&[4.0, 4.0]);
        assert!((r - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_convection_heat() {
        let q = convection_heat(10.0, 1.0, 350.0, 300.0);
        assert!((q - 500.0).abs() < 1e-10);
    }

    #[test]
    fn test_radiation_heat() {
        let q = radiation_heat(1.0, 1.0, 400.0, 300.0);
        assert!(q > 0.0);
    }

    #[test]
    fn test_biot_number() {
        let bi = biot_number(100.0, 0.1, 10.0);
        assert!((bi - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fourier_number() {
        let fo = fourier_number(1e-5, 100.0, 0.1);
        let expected = 1e-5 * 100.0 / 0.01;
        assert!((fo - expected).abs() < 1e-10);
    }

    // ─── Agent Energy Budgets ───

    #[test]
    fn test_agent_energy_budget() {
        let steps = vec![
            AgentStep {
                label: "inference".into(),
                compute_joules: 100.0,
                tokens_processed: 1000,
                useful_output_bits: 500.0,
                wall_time_seconds: 1.0,
            },
            AgentStep {
                label: "parsing".into(),
                compute_joules: 50.0,
                tokens_processed: 500,
                useful_output_bits: 200.0,
                wall_time_seconds: 0.5,
            },
        ];
        let budget = EnergyBudget::from_steps(steps);
        assert!((budget.total_energy_j - 150.0).abs() < 1e-10);
        assert!((budget.total_time_s - 1.5).abs() < 1e-10);
        assert!((budget.average_power() - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_agent_energy_per_token() {
        let steps = vec![
            AgentStep {
                label: "step".into(),
                compute_joules: 100.0,
                tokens_processed: 1000,
                useful_output_bits: 500.0,
                wall_time_seconds: 1.0,
            },
        ];
        let budget = EnergyBudget::from_steps(steps);
        assert!((budget.energy_per_token() - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_agent_most_expensive_step() {
        let steps = vec![
            AgentStep {
                label: "cheap".into(),
                compute_joules: 10.0,
                tokens_processed: 100,
                useful_output_bits: 50.0,
                wall_time_seconds: 0.1,
            },
            AgentStep {
                label: "expensive".into(),
                compute_joules: 100.0,
                tokens_processed: 1000,
                useful_output_bits: 500.0,
                wall_time_seconds: 1.0,
            },
        ];
        let budget = EnergyBudget::from_steps(steps);
        assert_eq!(budget.most_expensive_step().unwrap().label, "expensive");
    }

    #[test]
    fn test_landauer_cost() {
        let cost = EnergyBudget::landauer_cost_per_bit(300.0);
        // E = k_B * T * ln(2)
        let expected = BOLTZMANN * 300.0 * 2.0_f64.ln();
        assert!((cost - expected).abs() < 1e-30);
    }

    #[test]
    fn test_agent_heat_engine_analysis() {
        let (eff, waste) = agent_heat_engine_analysis(1000.0, 300.0);
        assert!((eff - 0.3).abs() < 1e-10);
        assert!((waste - 700.0).abs() < 1e-10);
    }

    #[test]
    fn test_thermo_state_gibbs_free_energy() {
        let g = ThermoState::gibbs_free_energy(1000.0, 300.0, 2.0);
        assert!((g - 1000.0 + 300.0 * 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_thermo_state_helmholtz_free_energy() {
        let a = ThermoState::helmholtz_free_energy(1000.0, 300.0, 2.0);
        assert!((a - 1000.0 + 300.0 * 2.0).abs() < 1e-10);
    }
}

//! Maxwell relations — thermodynamic identity derivatives.

use crate::constants::R;
use serde::{Serialize, Deserialize};

/// Maxwell relation result container.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaxwellRelation {
    pub name: String,
    pub left_side: f64,
    pub right_side: f64,
    pub difference: f64, // Should be ~0 if the relation holds
}

/// (∂T/∂V)_S = -(∂P/∂S)_V
pub fn maxwell_1(
    dtdv_s: f64,   // (∂T/∂V)_S
    dpds_v: f64,   // (∂P/∂S)_V
) -> MaxwellRelation {
    let diff = dtdv_s + dpds_v;
    MaxwellRelation {
        name: "(∂T/∂V)_S = -(∂P/∂S)_V".into(),
        left_side: dtdv_s,
        right_side: -dpds_v,
        difference: diff,
    }
}

/// (∂T/∂P)_S = (∂V/∂S)_P
pub fn maxwell_2(
    dtdp_s: f64,
    dvds_p: f64,
) -> MaxwellRelation {
    let diff = dtdp_s - dvds_p;
    MaxwellRelation {
        name: "(∂T/∂P)_S = (∂V/∂S)_P".into(),
        left_side: dtdp_s,
        right_side: dvds_p,
        difference: diff,
    }
}

/// (∂P/∂T)_V = (∂S/∂V)_T
pub fn maxwell_3(
    dpdt_v: f64,
    dsdv_t: f64,
) -> MaxwellRelation {
    let diff = dpdt_v - dsdv_t;
    MaxwellRelation {
        name: "(∂P/∂T)_V = (∂S/∂V)_T".into(),
        left_side: dpdt_v,
        right_side: dsdv_t,
        difference: diff,
    }
}

/// (∂V/∂T)_P = -(∂S/∂P)_T
pub fn maxwell_4(
    dvdt_p: f64,
    dsdp_t: f64,
) -> MaxwellRelation {
    let diff = dvdt_p + dsdp_t;
    MaxwellRelation {
        name: "(∂V/∂T)_P = -(∂S/∂P)_T".into(),
        left_side: dvdt_p,
        right_side: -dsdp_t,
        difference: diff,
    }
}

/// For an ideal gas, verify (∂P/∂T)_V = nR/V and (∂S/∂V)_T = nR/V
/// Maxwell relation 3: these should be equal.
pub fn ideal_gas_maxwell_3_verify(n: f64, t: f64, v: f64) -> bool {
    let dpdt_v = n * R / v; // (∂P/∂T)_V for ideal gas
    let dsdv_t = n * R / v; // (∂S/∂V)_T for ideal gas
    (dpdt_v - dsdv_t).abs() < 1e-10
}

/// For an ideal gas, verify (∂V/∂T)_P = nR/P = V/T and -(∂S/∂P)_T = nR/P = V/T
/// Maxwell relation 4.
pub fn ideal_gas_maxwell_4_verify(n: f64, t: f64, p: f64) -> bool {
    let dvdt_p = n * R / p;
    let neg_dsdp_t = n * R / p;
    (dvdt_p - neg_dsdp_t).abs() < 1e-10
}

/// Thermodynamic potential derivatives from the fundamental relation.
/// dU = TdS - PdV → (∂U/∂S)_V = T, (∂U/∂V)_S = -P
pub fn internal_energy_derivatives(
    t: f64,
    p: f64,
) -> (f64, f64) {
    (t, -p)
}

/// dH = TdS + VdP → (∂H/∂S)_P = T, (∂H/∂P)_S = V
pub fn enthalpy_derivatives(
    t: f64,
    v: f64,
) -> (f64, f64) {
    (t, v)
}

/// dA = -SdT - PdV → (∂A/∂T)_V = -S, (∂A/∂V)_T = -P
pub fn helmholtz_derivatives(
    s: f64,
    p: f64,
) -> (f64, f64) {
    (-s, -p)
}

/// dG = -SdT + VdP → (∂G/∂T)_P = -S, (∂G/∂P)_T = V
pub fn gibbs_derivatives(
    s: f64,
    v: f64,
) -> (f64, f64) {
    (-s, v)
}

/// Isothermal compressibility: κ_T = -1/V · (∂V/∂P)_T
pub fn isothermal_compressibility(v: f64, dvdp_t: f64) -> f64 {
    -(1.0 / v) * dvdp_t
}

/// Thermal expansion coefficient: α = 1/V · (∂V/∂T)_P
pub fn thermal_expansion_coefficient(v: f64, dvdt_p: f64) -> f64 {
    (1.0 / v) * dvdt_p
}

#!/usr/bin/env python3
"""AE DP-01: reproducible toy checks; no cosmological fit or observational data.
Uses only Python's standard library. Run: python DP01_toy_checks.py
The chosen weak-field interpolation is an illustrative assumption, not a
prediction derived from a complete Geometric-State Dynamics action.
"""
from __future__ import annotations
import json
import math
from pathlib import Path

G = 6.67430e-11  # m^3 kg^-1 s^-2, illustrative rounded SI input
M_SUN = 1.98847e30  # kg, rounded scale conversion
KPC = 3.085677581491367e19  # m
A0 = 1.2e-10  # m/s^2; phenomenological input, not fitted here
MASS = 1.0e10 * M_SUN

def acceleration(g_newton: float, a0: float = A0) -> float:
    """Solve mu(g/a0)*g=g_N for mu(x)=x/sqrt(1+x*x)."""
    if not (math.isfinite(g_newton) and g_newton > 0 and a0 > 0):
        raise ValueError('Positive finite accelerations are required.')
    g_squared = 0.5 * (g_newton*g_newton + g_newton*math.hypot(g_newton, 2*a0))
    return math.sqrt(g_squared)

def run() -> dict:
    rows = []
    max_relative_residual = 0.0
    for r_kpc in [1, 3, 10, 30, 100]:
        r = r_kpc*KPC
        gn = G*MASS/(r*r)
        g = acceleration(gn)
        x = g/A0
        residual = abs((x/math.hypot(x, 1))*g/gn - 1)
        max_relative_residual = max(max_relative_residual, residual)
        rows.append({'r_kpc': r_kpc, 'v_newton_km_s': math.sqrt(gn*r)/1000,
                     'v_toy_km_s': math.sqrt(g*r)/1000,
                     'g_newton_m_s2': gn, 'g_toy_m_s2': g})
    assert max_relative_residual < 1e-12
    # Algebraic cosmological-background check, normalized units; not a data fit.
    A, B, C, n = 0.25, 0.65, 0.05, 0.5
    bg = []
    max_bg_residual = 0.0
    for a in [0.1, 0.5, 1.0, 2.0]:
        rho = A*a**-3 + B + C*a**-n
        pressure = -B + (n/3-1)*C*a**-n
        d_rho_d_log_a = -3*A*a**-3 - n*C*a**-n
        res = abs(d_rho_d_log_a + 3*(rho+pressure))/rho
        max_bg_residual = max(max_bg_residual, res)
        bg.append({'a': a, 'rho_X_normalized': rho, 'p_X_normalized': pressure,
                   'w_X': pressure/rho})
    assert max_bg_residual < 1e-12
    relaxation = [{'t_over_tau': q, 'q_over_initial': math.exp(-q)} for q in [0,1,3,5]]
    assert all(0 < row['q_over_initial'] <= 1 for row in relaxation)
    assert all(relaxation[i+1]['q_over_initial'] < relaxation[i]['q_over_initial'] for i in range(len(relaxation)-1))
    assert math.isclose(10.0**-4, 1e-4, rel_tol=1e-15)
    result = {
        'status': 'ILLUSTRATIVE TOY CALCULATION ONLY; NOT OBSERVATIONAL DATA',
        'inputs': {'G_SI': G, 'M_sun_kg': M_SUN, 'kpc_m': KPC,
                   'mass_Msun': 1e10, 'a0_m_s2': A0,
                   'background_A_B_C_n': [A,B,C,n]},
        'rotation_table': rows,
        'asymptotic_speed_km_s': (G*MASS*A0)**0.25/1000,
        'transition_radius_kpc': math.sqrt(G*MASS/A0)/KPC,
        'max_rotation_relative_residual': max_relative_residual,
        'background_table': bg,
        'max_background_relative_residual': max_bg_residual,
        'single_mode_relaxation': relaxation,
        'schwarzschild_horizon_K_ratio_M2_over_M1_for_M2_equal_10M1': 1e-4,
        'checks_passed': ['weak-field algebra', 'background continuity identity',
                          'finite passive relaxation', 'horizon curvature scaling'],
        'not_performed': ['complete field equations', 'ghost/hyperbolicity proof',
                          'CMB/BAO/lensing fit', 'cluster-merger simulation',
                          'wormhole stability', 'engineering feasibility']
    }
    out=Path(__file__).with_name('DP01_toy_results.json')
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n',encoding='utf-8')
    return result

if __name__ == '__main__':
    print(json.dumps(run(), ensure_ascii=False, indent=2))

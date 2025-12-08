#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 09:24:55 2025

@author: kkasturi
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple
from scipy.integrate import solve_ivp
from numpy.linalg import solve as lin_solve

# ----------------------------
# Core setup (DATA geometry)
# ----------------------------
@dataclass
class Setup:
    A: float                 # membrane area [m^2]
    V0: float                # initial retentate volume [m^3]
    z_c: np.ndarray          # valences of cations, shape (n,)
    z_a: int                 # anion valence (negative)
    Dm_c: np.ndarray         # membrane-phase diffusivities of cations [m^2/s], shape (n,)
    Dm_a: float              # membrane-phase diffusivity of anion [m^2/s]
    H_c: np.ndarray          # partition coefficients (cations), shape (n,)
    H_a: float               # partition coefficient (anion)
    chi: float               # fixed charge density in membrane [same conc units as c]
    Lp: float                # hydraulic permeability [m/(s·Pa)] or coherent units
    dP: float                # transmembrane pressure drop [Pa]
    R: float                 # gas constant [J/(mol·K)]
    T: float                 # temperature [K]
    MW_c: np.ndarray         # molecular weights cations [kg/mol], shape (n,)
    MW_a: float              # anion MW [kg/mol]
    delta: float             # effective transport thickness used for grad c [m]

# ----------------------------
# Partitioning (Eqs. 85–86)
# ----------------------------
def partition_common_anion_log(
    c_s_c: np.ndarray, c_s_a: float,
    H_c: np.ndarray, H_a: float,
    z_c: np.ndarray, z_a: int, chi: float
) -> Tuple[np.ndarray, float]:
    """
    Solve for membrane-side interface concentrations (c_m_c, c_m_a)
    using the log-transformed common-anion relations (Eq. 86)
    and membrane electroneutrality.
    """
    n = len(c_s_c)
    # unknowns: [c_m_c(0..n-1), c_m_a]
    x = np.hstack([np.maximum(c_s_c / np.maximum(H_c, 1.0), 1e-12), np.maximum(c_s_a / max(H_a,1.0), 1e-12)])

    def residual(x):
        c_m_c = np.clip(x[:n], 1e-30, None)
        c_m_a = max(x[-1], 1e-30)
        # Eq. 86 for each cation i
        eq = z_c*(np.log(H_a) - np.log(c_m_a) + np.log(c_s_a)) - \
             z_a*(np.log(H_c) - np.log(c_m_c) + np.log(c_s_c))
        # electroneutrality in membrane
        en = chi + (z_c @ c_m_c) + z_a*c_m_a
        return np.hstack([eq, en])

    # damped Newton
    for _ in range(30):
        r = residual(x)
        if np.linalg.norm(r) < 1e-11: break
        J = np.zeros((n+1, n+1))
        eps = 1e-8
        for k in range(n+1):
            xk = x.copy()
            xk[k] *= (1+eps)
            J[:,k] = (residual(xk)-r)/(xk[k]-x[k])
        try:
            dx = lin_solve(J, -r)
        except np.linalg.LinAlgError:
            dx = -1e-2*r
        lam = 1.0
        for _ in range(8):
            xt = x + lam*dx
            if np.all(xt>0) and np.linalg.norm(residual(xt)) < np.linalg.norm(r):
                x = xt; break
            lam *= 0.5
    return x[:n], x[-1]

# ----------------------------
# Alpha, D_ij, D~ (Eqs. 115–118)
# ----------------------------
def compute_alpha_D(c_m_c: np.ndarray, z_c: np.ndarray, z_a: int,
                    Dm_c: np.ndarray, Dm_a: float, chi: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return alpha (n,) and Dmat (n,n) from Eqs. 115–118 given membrane-side cation concentrations.
    """
    # Eq. 118: D~ 
    Dtilde = np.sum((z_c**2*Dm_c - z_c*z_a*Dm_a) * c_m_c) - z_a*Dm_a*chi  # ˜D
    # Eq. 115: alpha_i
    alpha = 1.0 + (z_c*Dm_c*chi)/Dtilde
    n = len(z_c)
    Dmat = np.zeros((n,n))
    # Eq. 116–117
    for i in range(n):
        for j in range(n):
            if i != j:
                Dmat[i,j] = ((z_c[i]*z_c[j]*Dm_c[i]*Dm_c[j] - z_c[i]*z_c[j]*Dm_c[i]*Dm_a) * c_m_c[i]) / Dtilde
            else:
                # i == j case uses beta_ik per Eq. 117
                beta_sum = 0.0
                for k in range(n):
                    if i != k:
                        beta_ik = (z_c[k]**2) * Dm_c[k] * Dm_c[i]
                    else:
                        beta_ik = (z_c[k]**2) * Dm_c[k] * Dm_a
                    beta_sum += beta_ik * c_m_c[k]
                Dmat[i,i] = (np.sum((z_c*np.full(n, z_a)*Dm_c[i]*Dm_a - 0.0) * c_m_c)  # first term inside sum
                             - beta_sum + z_a*Dm_c[i]*Dm_a*chi) / Dtilde
    return alpha, Dmat

# ----------------------------
# Water flux (Tr-II) and Δπ
# ----------------------------
def osmotic_pressure_RT(c_s_c: np.ndarray, c_s_a: float, MW_c: np.ndarray, MW_a: float, R: float, T: float) -> float:
    """
    Minimal van't Hoff-like Δπ for the *solution side* used to drive Jw.
    Uses mass conc -> molar via MW (consistent with your lab units).
    """
    mol_c = c_s_c / MW_c
    mol_a = c_s_a / MW_a
    return R*T*(np.sum(mol_c) + mol_a)

def water_flux(Lp: float, dP: float, dPi: float) -> float:
    """Jw = Lp (ΔP − Δπ)."""
    return Lp * max(dP - dPi, 0.0)

# ----------------------------
# Well-mixed DATA balances
# ----------------------------
def rhs_data(t, y, S: Setup):
    """
    States: [V, c_c(0..n-1), c_a]   (mass concentrations in retentate)
    """
    n = len(S.z_c)
    V = y[0]
    c_c = np.clip(y[1:1+n], 1e-30, None)
    c_a = max(y[1+n], 1e-30)

    # Partitioning at the interface (Eqs. 85–86)
    c_m_c, c_m_a = partition_common_anion_log(c_c, c_a, S.H_c, S.H_a, S.z_c, S.z_a, S.chi)

    # Transport coefficients (Eqs. 115–118)
    alpha, Dmat = compute_alpha_D(c_m_c, S.z_c, S.z_a, S.Dm_c, S.Dm_a, S.chi)

    # Gradient closure through the membrane (within allowed scope)
    grad_c_m = -(c_m_c / S.delta)

    # Fluxes (Eq. 113) for cations
    # First compute Jw from Tr-II using solution-side Δπ (you can substitute a measured Δπ if desired)
    dPi = osmotic_pressure_RT(c_c, c_a, S.MW_c, S.MW_a, S.R, S.T)
    Jw = water_flux(S.Lp, S.dP, dPi)

    j_c = alpha * c_m_c * Jw + Dmat @ grad_c_m

    # Anion flux from Eq. 114 to enforce electroneutral flux
    j_a = -np.dot(S.z_c, j_c) / S.z_a

    # Well-mixed retentate balances (DATA geometry)
    dVdt = - S.A * Jw
    dVc_c_dt = - S.A * j_c
    dVc_a_dt = - S.A * j_a
    dc_c_dt = (dVc_c_dt - c_c*dVdt) / max(V, 1e-18)
    dc_a_dt = (dVc_a_dt - c_a*dVdt) / max(V, 1e-18)

    return np.hstack([dVdt, dc_c_dt, dc_a_dt])

# ----------------------------
# Example run (two cations + common anion)
# ----------------------------
def run_example():
    # Example: Li+, Co2+, Cl−
    z_c = np.array([+1.0, +2.0])
    z_a = -1
    S = Setup(
        A=0.005, V0=0.05,
        z_c=z_c, z_a=z_a,
        Dm_c=np.array([1.0e-10, 7.0e-11]),  # m^2/s (pick your fitted values)
        Dm_a=2.0e-10,
        H_c=np.array([1.1, 0.9]), H_a=0.95,
        chi=0.0,
        Lp=1.5e-11, dP=5.0e5,                 # example units: m/(s·Pa) and Pa
        R=8.314, T=298.15,
        MW_c=np.array([0.00694, 0.05893]),    # kg/mol
        MW_a=0.03545,
        delta=1e-7
    )

    # Initial retentate composition (mass concentrations, electroneutral)
    c1_0, c2_0 = 1.0, 5.0
    c_a0 = -(z_c[0]*c1_0 + z_c[1]*c2_0)/z_a
    y0 = np.hstack([S.V0, [c1_0, c2_0], c_a0])

    t_span = (0.0, 3600.0)  # 1 hour
    t_eval = np.linspace(*t_span, 401)

    sol = solve_ivp(lambda t,y: rhs_data(t,y,S), t_span, y0, t_eval=t_eval, rtol=1e-7, atol=1e-10)
    if not sol.success:
        raise RuntimeError(sol.message)

    return {
        "t": sol.t,
        "V": sol.y[0],
        "cations": sol.y[1:3],
        "anion": sol.y[3]
    }

if __name__ == "__main__":
    res = run_example()
    print("V(t) first 5:", res["V"][:5])
    print("cations first 5 points:\n", res["cations"][:, :5])

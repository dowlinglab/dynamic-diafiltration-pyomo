#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 13:23:38 2025

@author: kkasturi
"""

import numpy as np

from utility_v3_KKR import (
    list_sheets,
    load_excel_experiment,
    build_data_stru_from_loader,
    solve_model,
    calc_FIM,
    # summarize_FIM,  # uncomment if you added this helper
)

# 1) Path to THIS experiment's Excel file
excel_path = "/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/UnifiedFramework/ExperimentalDataFiles/NF270_MC3.xlsx"

# Show available sheets in the Excel file
print(list_sheets(excel_path))

# 2) Load the experimental data from the desired sheet
ldr = load_excel_experiment(excel_path, measure_sheet="07.09.24_NaCl")
# ldr.load() is already called inside load_excel_experiment

# 3) Plot raw experimental data (optional)
ldr.plot_mass(save_prefix="NF270_MC3_expdata_mass_cumulative", cumulative=True)
ldr.plot_conductivity(save_prefix="NF270_MC3_expdata_conductivity")
# ldr.plot_quicklook(save_prefix="NF270_MC3_expdata_mass_conductivity", cumulative=True)

# 4) Define THIS experiment's metadata / configuration
data_config = {
    'delP'   : 5.0,
    'Temp'   : 298.15,
    'Am'     : 100.0,
    'rho'    : 1.0,
    'ni'     : 1,
    'M_F0'   : 50.0,
    'M_O'    : 5.0,
    'C_D'    : 100.0,
    'C_F0'   : 100.0,
    'Lp0'    : 10.0,
    'B0'     : 1.0,
    'sigma0' : 0.5,
    'namec'  : "KCl",
    'n_v0'   : 1,
    'n_extra': 0,
    'n_h'    : 0,
    'n_A'    : 0,
    # 'n' will be set automatically from len(ldr.times)
}

# 5) Build a data_stru compatible with the existing .mat workflow
data_stru = build_data_stru_from_loader(ldr, data_config, dataset_id=270511)

# 6) Solve the model & get parameter estimates
fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode="Lag", B_form="single")

# 7) Run FIM / sensitivity analysis around the fitted parameters
theta = fit_stru["parameters"]

doe_stru = calc_FIM(
    data_stru,
    mode="Lag",
    theta=theta,        # expansion point
    step=1e-6,          # a bit larger than 1e-8 for numerical stability
    formula="backward", # or 'central'
    B_form="single",
)

# 8) Extract core FIM quantities
Jac = np.array(doe_stru["Jac"])       # sensitivity (Jacobian) matrix
FIM = np.array(doe_stru["FIM"])       # Fisher Information Matrix
eig_vals = np.array(doe_stru["eig_val"])
eig_vecs = np.array(doe_stru["eig_vec"])

print("Jacobian shape:", Jac.shape)
print("FIM shape:", FIM.shape)
print("Eigenvalues(FIM):", eig_vals)

# If you added summarize_FIM in utility_v3_KKR.py:
# summarize_FIM(doe_stru)

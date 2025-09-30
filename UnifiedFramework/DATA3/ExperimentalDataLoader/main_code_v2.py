#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 13:23:38 2025

@author: kkasturi
"""

from experiment_load_v1 import load_experiment

load_experiment("/Users/kkasturi/Desktop/ND/Projects/Diafilteration/Codes/ExperimentalDataFiles/NF270_MC3.xlsx", list_sheets=True)

data = load_experiment(
    "/Users/kkasturi/Desktop/ND/Projects/Diafilteration/Codes/ExperimentalDataFiles/NF270_MC3.xlsx",
    sheet_name="07.09.24_NaCl",
    plot=True,
    save_prefix="NF270_MC3_vialsplit"
)

times = data["times"]           # List[np.ndarray], one per vial
masses = data["masses"]         # Δmass per vial
retent = data["retentate_cond"]
permea  = data["permeate_cond"]
windows = data["windows"]       # [(start_idx, end_idx)]

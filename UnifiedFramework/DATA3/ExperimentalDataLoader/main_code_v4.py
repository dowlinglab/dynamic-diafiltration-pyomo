#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 13:23:38 2025

@author: kkasturi
"""

from experiment_dataload_OOP_v4 import DataLoader, list_sheets

path = "/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/UnifiedFramework/ExperimentalDataFiles/NF270_MC3.xlsx"

# list_sheets=True equivalent:
print(list_sheets(path))

# sheet_name="07.09.24_NaCl" equivalent:
ldr = DataLoader(path, measure_sheet="07.09.24_NaCl")
ldr.load()

times  = ldr.times
masses = ldr.masses
retent = ldr.retentate_cond
perm   = ldr.permeate_cond
wins   = ldr.windows

# ldr.plot_quicklook(save_prefix="NF270_MC3_vialsplit")

# Only mass, per-vial Δmass (default)
# ldr.plot_mass(save_prefix="NF270_MC3_expdata_mass")

# Only mass, cumulative across vial changes
ldr.plot_mass(save_prefix="NF270_MC3_expdata_mass_cumulative", cumulative=True)

# Only conductivity
ldr.plot_conductivity(save_prefix="NF270_MC3_expdata_conductivity")

# Both (quicklook), cumulative mass toggle optional
# ldr.plot_quicklook(save_prefix="NF270_MC3_expdata_mass_conductivity", cumulative=True)

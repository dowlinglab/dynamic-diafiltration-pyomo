#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 13:23:38 2025

@author: kkasturi
"""
import utility_v2_KKR as util
# from utility_v2_KKR import DataLoader, list_sheets, loadmat
from utility_v2_KKR import *
# path = "/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/UnifiedFramework/ExperimentalDataFiles/NF270_MC3.xlsx"

path = loadmat('/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset270511.123.mat')['data_stru']

# list_sheets=True equivalent:
print(list_sheets(path))

# sheet_name="07.09.24_NaCl" equivalent:
ldr = util.DataLoader(path, measure_sheet="07.09.24_NaCl")
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

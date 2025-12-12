#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for: experiment_dataload_xlsx_v2.py

One Excel sheet = one experiment.
Mirrors the MAT test script behavior.
"""

import numpy as np
import matplotlib.pyplot as plt

from experiment_dataload_xlsx_v2 import (
    extract_excel_variables,
    list_sheets,
    plot_mass_and_conductivity,
)

# -----------------------------------------------------------------------------
# 1) Point to an .xlsx file
# -----------------------------------------------------------------------------
path = "/Users/kkasturi/Downloads/NF270_MC2.xlsx"

print("Available experiments (sheets):")
print(list_sheets(path))

# Choose one experiment (sheet)
sheet = list_sheets(path)[0]

# -----------------------------------------------------------------------------
# 2) Load
# -----------------------------------------------------------------------------
vars_ = extract_excel_variables(path, measure_sheet=sheet)

# -----------------------------------------------------------------------------
# 3) Inventory
# -----------------------------------------------------------------------------
print("\n==============================")
print("[info] XLSX variables extracted")
print("==============================")
print("path:", vars_["path"])
print("format:", vars_["format"])
print("dataset:", vars_["dataset"])
print("filename:", vars_["filename"])
print("mode:", vars_["mode"])
print("continuous_cF:", vars_["continuous_cF"])
print("n_obj:", vars_.get("n_obj", None))

def _shape(x):
    return None if x is None else tuple(np.asarray(x).shape)

print("\n[shapes] full-run arrays")
print("time:", _shape(vars_["time"]))
print("mass:", _shape(vars_["mass"]))
print("cF_exp:", _shape(vars_["cF_exp"]))
print("cV_avg:", _shape(vars_["cV_avg"]))

print("\n[vials]")
print("n_vials:", len(vars_["vials"]))
print("vial_ranges:", vars_["vial_ranges"])

print("\n[per-vial shapes]")
print("per_vial_time:", [a.shape for a in vars_["per_vial_time"]])
print("per_vial_mass:", [a.shape for a in vars_["per_vial_mass"]])
print("per_vial_cF_exp:", [a.shape for a in vars_["per_vial_cF_exp"]])
print("per_vial_cV_avg:", [a.shape for a in vars_["per_vial_cV_avg"]])

# -----------------------------------------------------------------------------
# 4) Sanity checks
# -----------------------------------------------------------------------------
time = np.asarray(vars_["time"]).squeeze()
mass = np.asarray(vars_["mass"]).squeeze()
cF = np.asarray(vars_["cF_exp"]).squeeze()
cV = np.asarray(vars_["cV_avg"]).squeeze()
vial_ranges = vars_["vial_ranges"]

print("\n[sanity]")
print("equal lengths (time/mass/cF_exp/cV_avg):",
      len(time) == len(mass) == len(cF) == len(cV))

n = len(time)
ok_ranges = True
for (s, e) in vial_ranges:
    if not (0 <= s <= e < n):
        ok_ranges = False
print("vial_ranges in-bounds:", ok_ranges)

# -----------------------------------------------------------------------------
# 5) Plot (identical style to MAT test)
# -----------------------------------------------------------------------------
plot_mass_and_conductivity(vars_)

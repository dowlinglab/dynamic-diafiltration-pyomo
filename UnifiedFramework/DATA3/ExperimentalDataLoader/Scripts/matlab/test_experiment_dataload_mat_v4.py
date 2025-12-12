#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for: experiment_dataload_mat_v2.py (functional MAT parser)

What it does:
- loads a .mat file using extract_mat_variables()
- prints a compact inventory of returned variables
- sanity-checks array lengths + vial_ranges indices
- (optional) makes quick plots to confirm segmentation
"""

import numpy as np
import matplotlib.pyplot as plt

# Import module (not just the function) so we can toggle loader debugging.
import experiment_dataload_mat_v4 as edm2

# Set True only when debugging schema mismatches
edm2.DEBUG = False

extract_mat_variables = edm2.extract_mat_variables

# -----------------------------------------------------------------------------
# 1) Point to a .mat file
# -----------------------------------------------------------------------------
path = "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset270511.121.mat"

# -----------------------------------------------------------------------------
# 2) Load
# -----------------------------------------------------------------------------
vars_ = extract_mat_variables(path)

# -----------------------------------------------------------------------------
# 3) Inventory
# -----------------------------------------------------------------------------
print("\n==============================")
print("[info] MAT variables extracted")
print("==============================")
print("path:", vars_["path"])
print("format:", vars_["format"])
print("dataset:", vars_["dataset"])
print("filename:", vars_["filename"])
print("mode:", vars_["mode"])
print("continuous_cF:", vars_["continuous_cF"])
if vars_["continuous_cF"] is None:
    print("note: continuous_cF not stored in this .mat (OK under suggestion A)")
print("n_obj:", vars_["n_obj"])

def _shape(x):
    return None if x is None else tuple(np.asarray(x).shape)

print("\n[shapes] full-run arrays")
print("time:", _shape(vars_["time"]))
print("mass:", _shape(vars_["mass"]))
print("cF_exp:", _shape(vars_["cF_exp"]))
print("cV_avg:", _shape(vars_["cV_avg"]))
print("vial_marker:", _shape(vars_.get("vial_marker", None)))

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
# 5) Optional: quick visual checks (MATLAB-style plotting)
# -----------------------------------------------------------------------------
DO_PLOTS = True

def _add_vial_change_lines(ax, time, vial_ranges):
    """Add dashed vertical lines at vial boundaries."""
    for (s, _) in vial_ranges[1:]:
        if 0 <= s < len(time):
            ax.axvline(time[s], linestyle="--", linewidth=1)

if DO_PLOTS and n > 0:
    # --- Mass: single continuous line ---
    fig, ax = plt.subplots()
    ax.plot(time, mass, linestyle="-", label="mass")
    _add_vial_change_lines(ax, time, vial_ranges)
    ax.set_xlabel("time")
    ax.set_ylabel("mass")
    ax.set_title("Mass (raw) with vial changes")
    ax.legend()

    # --- Conductivity / concentration ---
    fig, ax = plt.subplots()
    ax.plot(time, cF, linestyle="-", label="cF_exp (retentate)")
    ax.plot(time, cV, linestyle="--", label="cV_avg (permeate)")
    _add_vial_change_lines(ax, time, vial_ranges)
    ax.set_xlabel("time")
    ax.set_ylabel("signal")
    ax.set_title("cF_exp / cV_avg with vial changes")
    ax.legend()

    plt.show()

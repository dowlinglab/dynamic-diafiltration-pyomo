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
import experiment_dataload_mat_v2 as edm2

# Set True to print detected MATLAB field names for data_stru (helpful when a
# flag like continuous_cF is stored under a different name, e.g. conductivity_cF).
edm2.DEBUG = False

extract_mat_variables = edm2.extract_mat_variables

# -----------------------------------------------------------------------------
# 1) Point to a .mat file (update this)
# -----------------------------------------------------------------------------
path = "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset270511.121.mat"

# -----------------------------------------------------------------------------
# 2) Load
# -----------------------------------------------------------------------------
vars_ = extract_mat_variables(path)

# -----------------------------------------------------------------------------
# 3) Quick inventory (like list_sheets + attribute access checks)
# -----------------------------------------------------------------------------
print("\n==============================")
print("[info] MAT variables extracted")
print("==============================")
print("path:", vars_["path"])
print("format:", vars_["format"])
print("dataset:", vars_["dataset"])
print("filename:", vars_["filename"])
print("mode:", vars_["mode"])
print("continuous_cF:", vars_.get("continuous_cF", None))
print("n_obj:", vars_["n_obj"])

def _shape(x):
    return None if x is None else tuple(np.asarray(x).shape)

print("\n[shapes] full-run arrays")
print("time:", _shape(vars_["time"]))
print("mass:", _shape(vars_["mass"]))
print("cF_exp:", _shape(vars_.get("cF_exp", None)))
print("cV_avg:", _shape(vars_.get("cV_avg", None)))
print("vial_marker:", _shape(vars_.get("vial_marker", None)))

print("\n[vials]")
print("n_vials:", len(vars_["vials"]))
print("vial_ranges:", vars_.get("vial_ranges", None))

print("\n[per-vial shapes]")
print("per_vial_time:", [a.shape for a in vars_["per_vial_time"]])
print("per_vial_mass:", [a.shape for a in vars_["per_vial_mass"]])
print("per_vial_cF_exp:", [a.shape for a in vars_.get("per_vial_cF_exp", [])])
print("per_vial_cV_avg:", [a.shape for a in vars_.get("per_vial_cV_avg", [])])

# -----------------------------------------------------------------------------
# 4) Sanity checks
# -----------------------------------------------------------------------------
time = np.asarray(vars_["time"]).squeeze()
mass = np.asarray(vars_["mass"]).squeeze()

# Conductivity/concentration signals may be missing for some datasets
cF = vars_.get("cF_exp", None)
cV = vars_.get("cV_avg", None)
cF = np.asarray(cF).squeeze() if cF is not None else None
cV = np.asarray(cV).squeeze() if cV is not None else None

vial_ranges = list(vars_.get("vial_ranges", []))

print("\n[sanity]")
if cF is not None and cV is not None:
    print("equal lengths (time/mass/cF_exp/cV_avg):",
          len(time) == len(mass) == len(cF) == len(cV))
else:
    print("equal lengths (time/mass):", len(time) == len(mass))
    if cF is None:
        print("note: cF_exp missing (OK for some files)")
    if cV is None:
        print("note: cV_avg missing (OK for some files)")

n = len(time)
ok_ranges = True
for (s, e) in vial_ranges:
    if not (0 <= s <= e < n):
        ok_ranges = False
        print(f"  BAD vial_range: ({s},{e}) with n={n}")
print("vial_ranges in-bounds:", ok_ranges)

# -----------------------------------------------------------------------------
# 5) Optional: quick visual checks (like your plot_* methods)
# -----------------------------------------------------------------------------
DO_PLOTS = True
if DO_PLOTS and n > 0:
    # Full-run mass (raw) + per-vial mass (re-baselined + spike filtered)
    plt.figure()
    plt.plot(time, mass, label="mass (full-run raw)")
    plt.xlabel("time")
    plt.ylabel("mass")
    plt.title("Full-run mass (raw)")
    plt.legend()

    plt.figure()
    for i, (t_i, m_i) in enumerate(zip(vars_["per_vial_time"], vars_["per_vial_mass"]), start=1):
        plt.plot(t_i, m_i, label=f"vial {i} mass (rebased/filtered)")
    plt.xlabel("time")
    plt.ylabel("Δmass (rebased)")
    plt.title("Per-vial mass (rebased + spike filtered)")
    plt.legend()

    # cF_exp / cV_avg quicklook (if present)
    if cF is not None or cV is not None:
        plt.figure()
        if cF is not None:
            plt.plot(time, cF, label="cF_exp")
        if cV is not None:
            plt.plot(time, cV, label="cV_avg")
        plt.xlabel("time")
        plt.ylabel("signal")
        plt.title("Full-run cF_exp / cV_avg")
        plt.legend()

    # Vial boundaries overlay (if present)
    if vial_ranges:
        plt.figure()
        plt.plot(time, mass, label="mass (raw)")
        for j, (s, e) in enumerate(vial_ranges, start=1):
            plt.axvline(time[s], linestyle="--", label=f"vial{j} start" if j == 1 else None)
            plt.axvline(time[e], linestyle=":",  label=f"vial{j} end" if j == 1 else None)
        plt.xlabel("time")
        plt.ylabel("mass")
        plt.title("Mass with inferred vial boundaries")
        plt.legend()

    plt.show()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for: experiment_dataload_mat_v1.py (functional MAT parser)

What it does:
- loads a .mat file using extract_mat_variables()
- prints a compact inventory of returned variables
- sanity-checks array lengths + window indices
- (optional) makes quick plots to confirm segmentation
"""

import numpy as np
import matplotlib.pyplot as plt

from experiment_dataload_mat_v1 import extract_mat_variables

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
print("conductivity_cF:", vars_["conductivity_cF"])
print("n_obj:", vars_["n_obj"])

def _shape(x):
    return None if x is None else tuple(np.asarray(x).shape)

print("\n[shapes] full-run arrays")
print("time:", _shape(vars_["time"]))
print("mass:", _shape(vars_["mass"]))
print("retentate_cond:", _shape(vars_["retentate_cond"]))
print("permeate_cond:", _shape(vars_["permeate_cond"]))
print("vial_marker:", _shape(vars_["vial_marker"]))

print("\n[vials]")
print("n_vials:", len(vars_["vials"]))
print("windows:", vars_["windows"])

print("\n[per-vial shapes]")
print("per_vial_time:", [a.shape for a in vars_["per_vial_time"]])
print("per_vial_mass:", [a.shape for a in vars_["per_vial_mass"]])
print("per_vial_retentate_cond:", [a.shape for a in vars_["per_vial_retentate_cond"]])
print("per_vial_permeate_cond:", [a.shape for a in vars_["per_vial_permeate_cond"]])

# -----------------------------------------------------------------------------
# 4) Sanity checks
# -----------------------------------------------------------------------------
time = np.asarray(vars_["time"]).squeeze()
mass = np.asarray(vars_["mass"]).squeeze()
rcon = np.asarray(vars_["retentate_cond"]).squeeze()
pcon = np.asarray(vars_["permeate_cond"]).squeeze()
wins = list(vars_["windows"])

print("\n[sanity]")
print("equal lengths (time/mass/rcon/pcon):",
      len(time) == len(mass) == len(rcon) == len(pcon))

n = len(time)
ok_windows = True
for (s, e) in wins:
    if not (0 <= s <= e < n):
        ok_windows = False
        print(f"  BAD window: ({s},{e}) with n={n}")
print("windows in-bounds:", ok_windows)

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

    # Conductivity quicklook
    plt.figure()
    plt.plot(time, rcon, label="retentate_cond")
    plt.plot(time, pcon, label="permeate_cond")
    plt.xlabel("time")
    plt.ylabel("conductivity")
    plt.title("Full-run conductivity")
    plt.legend()

    # Window boundaries overlay (if present)
    if wins:
        plt.figure()
        plt.plot(time, mass, label="mass (raw)")
        for j, (s, e) in enumerate(wins, start=1):
            plt.axvline(time[s], linestyle="--", label=f"win{j} start" if j == 1 else None)
            plt.axvline(time[e], linestyle=":",  label=f"win{j} end" if j == 1 else None)
        plt.xlabel("time")
        plt.ylabel("mass")
        plt.title("Mass with inferred window boundaries")
        plt.legend()

    plt.show()

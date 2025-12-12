#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 13:23:38 2025

Author: kkasturi

Test script for BOTH:
    - Excel experimental files (.xlsx)
    - MATLAB experimental files (.mat, with data_stru)

Uses:
    DataLoader  (existing plotting + vial slicing)
    PyDataStru  (new unified constructors)
"""

from experiment_dataload_OOP_v7 import (
    DataLoader,
    list_sheets,
    PyDataStru,
)

# ---------------------------------------------------------------------------
# CHOOSE YOUR TEST FILE HERE
# ---------------------------------------------------------------------------

# Example Excel experiment file
excel_path = "/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/UnifiedFramework/ExperimentalDataFiles/NF270_MC3.xlsx"

# Example MAT experiment file (data_stru)
mat_path = "/Users/kkasturi/GitHub/dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset270511.121.mat"

# ---------------------------------------------------------------------------
# 1. TEST PyDataStru.from_excel()
# ---------------------------------------------------------------------------

print("\n========== TESTING PyDataStru.from_excel ==========\n")

# list available sheets in the Excel
print("Available sheets in Excel:")
print(list_sheets(excel_path))

sheet = "07.09.24_NaCl"   # choose a sheet from list_sheets()

ds_excel = PyDataStru.from_excel(excel_path, sheet)

print("\n--- Excel PyDataStru ---")
print("filename:", ds_excel.filename)
print("dataset:", ds_excel.dataset)
print("num vials:", len(ds_excel.data_raw))
print("first vial time length:", len(ds_excel.data_raw[0].time))
print("first vial mass length:", len(ds_excel.data_raw[0].mass))

# ---------------------------------------------------------------------------
# 2. TEST PyDataStru.from_mat()
# ---------------------------------------------------------------------------

print("\n========== TESTING PyDataStru.from_mat ==========\n")

ds_mat = PyDataStru.from_mat(mat_path)

print("\n--- MATLAB PyDataStru ---")
print("dataset:", ds_mat.dataset)
print("filename:", ds_mat.filename)
print("mode:", ds_mat.mode)
print("num vials:", len(ds_mat.data_raw))

v0 = ds_mat.data_raw[0]
print("first vial time length:", len(v0.time))
print("first vial mass length:", len(v0.mass))
print("first vial cF_exp exists:", v0.cF_exp is not None)
print("first vial cV_avg exists:", v0.cV_avg is not None)

# ---------------------------------------------------------------------------
# 3. TEST DataLoader on Excel
# ---------------------------------------------------------------------------

print("\n========== TESTING DataLoader (Excel) ==========\n")

ldr_excel = DataLoader(excel_path, measure_sheet=sheet)
ldr_excel.load()

print("Excel windows:", ldr_excel.windows)
print("Excel number of vials:", len(ldr_excel.times))

# Uncomment if you want plots:
# ldr_excel.plot_mass(save_prefix="excel_mass", cumulative=False)
# ldr_excel.plot_conductivity(save_prefix="excel_cond")

# ---------------------------------------------------------------------------
# 4. TEST DataLoader on MATLAB
# ---------------------------------------------------------------------------

print("\n========== TESTING DataLoader (MATLAB data_stru) ==========\n")

ldr_mat = DataLoader(mat_path)
ldr_mat.load()

print("MAT windows:", ldr_mat.windows)
print("MAT number of vials:", len(ldr_mat.times))

# Uncomment if you want plots:
ldr_mat.plot_mass(save_prefix="mat_mass", cumulative=False)
ldr_mat.plot_conductivity(save_prefix="mat_cond")

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

print("\n========== ALL TESTS COMPLETE ==========\n")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 16:19:48 2025

@author: kkasturi
"""

from simple_excel_loader_v3 import SimpleExcelLoader

loader = SimpleExcelLoader('/Users/kkasturi/Desktop/ND/Projects/Diafilteration/Codes/ExperimentalDataFiles/NF270_MC3.xlsx')
loader.load()

# --- ensure .windows exists even if loader is an older version ---
import numpy as np, pandas as pd

def ensure_windows(loader):
    if hasattr(loader, "windows") and isinstance(getattr(loader, "windows"), list):
        return loader.windows

    df = getattr(loader, "_meas_df", None)
    # Fallback: single full-range window if we only have times
    if df is None or getattr(df, "empty", True):
        if loader.times:
            w = [(0, len(loader.times[0]) - 1)]
            setattr(loader, "windows", w)
            return w
        setattr(loader, "windows", [])
        return []

    # Try to infer from a 'Vial' / 'Vial Swap' column
    vial_col = next((c for c in df.columns if "vial" in str(c).lower()), None)
    if vial_col is not None:
        v = pd.to_numeric(df[vial_col], errors="coerce").fillna(-1).to_numpy()
        boundaries = [0] + [i for i in range(1, len(v)) if v[i] != v[i-1]] + [len(v)]
        w = [(s, e - 1) for s, e in zip(boundaries[:-1], boundaries[1:]) if (e - s) >= 3]
        if not w and loader.times:
            w = [(0, len(loader.times[0]) - 1)]
        setattr(loader, "windows", w)
        return w

    # No vial column: fallback
    w = [(0, len(loader.times[0]) - 1)] if loader.times else []
    setattr(loader, "windows", w)
    return w

windows = ensure_windows(loader)
print("windows =", windows[:5], "(total:", len(windows), ")")

# Variables
times = loader.times                  # List[np.ndarray]
masses = loader.masses                # Δmass per vial
retentate = loader.retentate_cond
permeate = loader.permeate_cond
windows = loader.windows              # [(start_idx, end_idx), ...]

# Plots (also saved if you pass a prefix)
loader.plot_quicklook(save_prefix="NF270_MC3_vialsplit")

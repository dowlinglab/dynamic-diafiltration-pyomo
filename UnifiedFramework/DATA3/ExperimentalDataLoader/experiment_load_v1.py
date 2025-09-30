#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 13:20:10 2025

@author: kkasturi
"""

# experiment_loader.py
# --------------------------------------------------------------------------------
# Unified loader for diafiltration-like data that reads Excel (.xlsx/.xls/.xlsm/.xlsb)
# or MATLAB (.mat) and returns simple Python/numpy structures ready for analysis/plotting.
#
# This version:
#   - Excel: list sheets, choose ONE sheet (return single-sheet dict), or MANY sheets
#            (return {"sheets": [...], "by_sheet": {...}}).
#   - Robust header-row detection (handles offset headers).
#   - Vial segmentation from a 'Vial'/'Vial Swap' column (contiguous segments).
#   - Optional quicklook plots in a consistent style (same color, vial boundaries, top axis).
# --------------------------------------------------------------------------------

from __future__ import annotations

from typing import List, Tuple, Optional, Dict, Union
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_experiment(
    path: str,
    plot: bool = False,
    save_prefix: Optional[str] = None,
    sheet_name: Optional[Union[str, List[str]]] = None,
    list_sheets: bool = False,
) -> Dict[str, object]:
    """
    Load experiment data from Excel or MATLAB and (optionally) plot.

    Parameters
    ----------
    path : str
        Path to the input file (.xlsx/.xls/.xlsm/.xlsb or .mat).
    plot : bool, default False
        If True, generate quicklook plots (mass & conductivity).
    save_prefix : str, optional
        If provided, save PNGs (e.g., f"{save_prefix}_mass.png"). When loading multiple
        Excel sheets, the sheet name is appended automatically (e.g., f"{save_prefix}_<sheet>_mass.png").
    sheet_name : str | list[str] | None, default None
        Excel only:
          - None  -> auto-detect a single measurement sheet (time + mass/cond)
          - str   -> process exactly that one sheet; return a single-sheet dict
          - list  -> process multiple sheets; return {"sheets":[...], "by_sheet":{...}}
    list_sheets : bool, default False
        Excel only. If True, print workbook sheet names and return {} immediately.

    Returns
    -------
    dict
        If Excel with multiple sheets requested:
            {
              "sheets": ["SheetA","SheetB",...],
              "by_sheet": {
                  "SheetA": {"times": [...], "masses":[...], "retentate_cond":[...],
                             "permeate_cond":[...], "windows":[(s,e),...]},
                  "SheetB": {...}
              }
            }
        Otherwise (single sheet or MAT):
            {
              "times": [...], "masses":[...], "retentate_cond":[...],
              "permeate_cond":[...], "windows":[(s,e),...]
            }
    """
    # ---------- Basic checks ----------
    ext = os.path.splitext(path)[1].lower()               # e.g., ".xlsx" or ".mat"
    if not os.path.exists(path):                          # fail early if path is wrong
        raise FileNotFoundError(f"File not found: {path}")

    # ---------- Sheet listing mode (Excel only) ----------
    if list_sheets:
        if ext not in {".xlsx", ".xls", ".xlsm", ".xlsb"}:
            print("list_sheets=True is only for Excel workbooks.")
            return {}
        xl = pd.ExcelFile(path)                           # open workbook
        print("Available sheets in workbook:")
        for sh in xl.sheet_names:
            print("  -", sh)
        return {}                                         # nothing else to return

    # ---------- Excel branch ----------
    if ext in {".xlsx", ".xls", ".xlsm", ".xlsb"}:
        xl = pd.ExcelFile(path)                           # open workbook once

        # 1) sheet_name is None → autodetect ONE measurement sheet and return single-sheet dict
        if sheet_name is None:
            sel = _autodetect_measure_sheet(xl)
            if sel is None:
                raise ValueError("Could not auto-detect a measurement-like sheet.")
            return _load_excel_one_sheet(path, sel, plot=plot, save_prefix=save_prefix)

        # 2) sheet_name is a string → process exactly that sheet and return single-sheet dict
        if isinstance(sheet_name, str):
            if sheet_name not in xl.sheet_names:
                raise ValueError(f"Sheet '{sheet_name}' not found. Available: {xl.sheet_names}")
            return _load_excel_one_sheet(path, sheet_name, plot=plot, save_prefix=save_prefix)

        # 3) sheet_name is a list → process multiple sheets and return multi-sheet dict
        names: List[str] = list(sheet_name)
        missing = [nm for nm in names if nm not in xl.sheet_names]
        if missing:
            raise ValueError(f"Sheet(s) not found: {missing}. Available: {xl.sheet_names}")

        by_sheet: Dict[str, dict] = {}
        for nm in names:
            # If saving, append sheet name to the prefix for uniqueness
            sp = f"{save_prefix}_{nm}" if save_prefix else None
            by_sheet[nm] = _load_excel_one_sheet(path, nm, plot=plot, save_prefix=sp)
        return {"sheets": names, "by_sheet": by_sheet}

    # ---------- MATLAB branch ----------
    if ext == ".mat":
        return _load_mat_single(path, plot=plot, save_prefix=save_prefix)

    # ---------- Unsupported ----------
    raise ValueError(f"Unsupported file type: {ext}")


# ============================== Excel helpers ==============================

def _autodetect_measure_sheet(xl: pd.ExcelFile) -> Optional[str]:
    """Pick the first sheet that looks like measurement data (has 'time' and ('mass' or 'cond'))."""
    for sh in xl.sheet_names:
        try:
            sniff = xl.parse(sh, nrows=10)                # read a few rows to inspect columns
        except Exception:
            continue
        cols = [str(c).lower() for c in sniff.columns]
        if any("time" in c for c in cols) and (any("mass" in c for c in cols) or any("cond" in c for c in cols)):
            return sh
    return None                                            # none matched


def _detect_header_row(raw_df: pd.DataFrame, max_rows: int = 10) -> int:
    """
    Find which row in the top 'max_rows' looks like the header:
    it should contain 'time' and ( 'mass' or 'cond' or 'pressure' ).
    """
    for i in range(min(max_rows, len(raw_df))):
        row = raw_df.iloc[i].astype(str).str.lower().tolist()
        has_time = any("time" in s for s in row)
        has_mass_or_cond = any(("mass" in s) or ("cond" in s) or ("pressure" in s) for s in row)
        if has_time and has_mass_or_cond:
            return i
    return 0                                               # fallback to the first row if none matched


def _first_col_like(df: pd.DataFrame, needle: str) -> Optional[str]:
    """Return the first column name that contains `needle` (case-insensitive), or None."""
    n = needle.lower()
    for c in df.columns:
        if n in str(c).lower():
            return c
    return None


def _load_excel_one_sheet(
    path: str,
    sheet: str,
    plot: bool = False,
    save_prefix: Optional[str] = None,
) -> Dict[str, object]:
    """
    Load exactly one Excel sheet and return a single-sheet dict with:
      times, masses, retentate_cond, permeate_cond, windows
    Optionally plot the result.
    """
    # --- Read the sheet without header first so we can detect the real header row ---
    xl = pd.ExcelFile(path)
    raw = xl.parse(sheet, header=None)                    # raw 2D values (header may be offset)
    hdr = _detect_header_row(raw)                         # which row is the true header?
    df = raw.iloc[hdr + 1 :].copy()                       # actual data are below the header row
    df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()  # set header row as column names
    df = df.reset_index(drop=True)                        # clean index to 0..N-1

    # --- Identify relevant columns (robust substring matching) ---
    tcol = _first_col_like(df, "time")                    # time column
    mcol = _first_col_like(df, "mass")                    # mass column
    rcol = next((c for c in df.columns if "retent" in str(c).lower() and "cond" in str(c).lower()), None)
    pcol = next((c for c in df.columns if "permeate" in str(c).lower() and "cond" in str(c).lower()), None)

    # --- Compute "vial windows" from a Vial/Vial Swap column if available ---
    windows: List[Tuple[int, int]] = []
    vcol = _first_col_like(df, "vial")
    if vcol is not None:
        v = pd.to_numeric(df[vcol], errors="coerce").fillna(-1).to_numpy()
        if np.isfinite(v).any():
            # Boundaries where vial value changes; add start=0 and end=len(v)
            boundaries = [0] + [i for i in range(1, len(v)) if v[i] != v[i - 1]] + [len(v)]
            # Inclusive windows (s..e-1) of minimum length 3 to skip single-row swaps
            windows = [(int(s), int(e - 1)) for s, e in zip(boundaries[:-1], boundaries[1:]) if (e - s) >= 3]
    if not windows and len(df) > 0:
        windows = [(0, len(df) - 1)]                      # fallback: one full-range window

    # --- Convert chosen columns into numeric arrays (sensible defaults if missing) ---
    to_num = lambda s: pd.to_numeric(s, errors="coerce").to_numpy()
    tvals = to_num(df[tcol]) if tcol else np.arange(len(df), dtype=float)
    mvals = to_num(df[mcol]) if mcol else np.full(len(df), np.nan)
    rvals = to_num(df[rcol]) if rcol else np.full(len(df), np.nan)
    pvals = to_num(df[pcol]) if pcol else np.full(len(df), np.nan)

    # --- Slice arrays into per-vial segments; compute Δmass; filter spikes > 1.2 g ---
    times:   List[np.ndarray] = []
    masses:  List[np.ndarray] = []
    retent:  List[np.ndarray] = []
    permea:  List[np.ndarray] = []
    for (s, e) in windows:
        sl = slice(s, e + 1)                              # inclusive slice [s, e]
        t = tvals[sl]                                     # time segment
        mraw = mvals[sl]                                  # absolute mass segment
        m = mraw - (mraw[0] if len(mraw) and not np.isnan(mraw[0]) else 0.0)  # Δmass per vial
        m = np.where(m > 1.2, np.nan, m)                  # spike filter per your spec
        r = rvals[sl]                                     # retentate conductivity segment
        p = pvals[sl]                                     # permeate conductivity segment
        times.append(t); masses.append(m); retent.append(r); permea.append(p)

    # --- Optional plotting for this sheet ---
    if plot:
        _plot_quicklook(
            times, masses, retent, permea, windows,
            save_prefix=(save_prefix if save_prefix else None)
        )

    # --- Return a single-sheet dictionary ---
    return {
        "times": times,
        "masses": masses,
        "retentate_cond": retent,
        "permeate_cond": permea,
        "windows": windows,
    }


# ============================== MATLAB helper ==============================

def _load_mat_single(
    path: str,
    plot: bool = False,
    save_prefix: Optional[str] = None,
) -> Dict[str, object]:
    """
    Load a single MATLAB .mat file into the same structure as one Excel sheet.
    (Assumes v7 MAT via scipy.io.loadmat. If your files are v7.3/HDF5, we can add h5py support.)
    """
    from scipy.io import loadmat                                    # v7 MAT support
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)     # read into dict-like

    # Helper: pull by any of several likely names
    def _grab(names: List[str]):
        for nm in names:
            if nm in mat:
                return np.array(mat[nm]).squeeze()
        return None

    # Try common keys; adapt if your MAT files differ
    time = _grab(["Time", "time", "t"])
    mass = _grab(["Mass", "mass", "m"])
    rcon = _grab(["RetentateCond", "Retentate_Cond", "retentate_cond", "ret_cond"])
    pcon = _grab(["PermeateCond", "Permeate_Cond", "permeate_cond", "perm_cond"])
    vial = _grab(["Vial", "VialSwap", "vial", "vial_swap"])

    if time is None or mass is None:
        raise ValueError("MAT file missing required arrays: Time and Mass")

    # Ensure consistent lengths based on time; pad/trim aux signals as needed
    n = len(np.array(time).squeeze())
    def _fit_len(x):
        if x is None: return np.full(n, np.nan)
        arr = np.array(x).squeeze()
        return arr[:n] if arr.size >= n else np.pad(arr, (0, n - arr.size), constant_values=np.nan)

    time = _fit_len(time).astype(float)
    mass = _fit_len(mass).astype(float)
    rcon = _fit_len(rcon).astype(float)
    pcon = _fit_len(pcon).astype(float)

    # Windows from vial if available (min length 3), else full-length
    windows: List[Tuple[int, int]] = [(0, n - 1)]
    if vial is not None:
        v = np.array(vial).squeeze()
        vnum = pd.to_numeric(pd.Series(v), errors="coerce").fillna(-1).to_numpy()
        boundaries = [0] + [i for i in range(1, len(vnum)) if vnum[i] != vnum[i - 1]] + [len(vnum)]
        cand = [(s, e - 1) for s, e in zip(boundaries[:-1], boundaries[1:]) if (e - s) >= 3]
        if cand:
            windows = cand

    # Build outputs
    times:   List[np.ndarray] = []
    masses:  List[np.ndarray] = []
    retent:  List[np.ndarray] = []
    permea:  List[np.ndarray] = []
    for (s, e) in windows:
        sl = slice(s, e + 1)
        t = time[sl]
        mraw = mass[sl]
        m = mraw - (mraw[0] if len(mraw) and not np.isnan(mraw[0]) else 0.0)  # Δmass
        m = np.where(m > 1.2, np.nan, m)                                      # spike filter
        r = rcon[sl]
        p = pcon[sl]
        times.append(t); masses.append(m); retent.append(r); permea.append(p)

    # Optional plotting
    if plot:
        _plot_quicklook(times, masses, retent, permea, windows, save_prefix=save_prefix)

    return {
        "times": times,
        "masses": masses,
        "retentate_cond": retent,
        "permeate_cond": permea,
        "windows": windows,
    }


# ============================== Plot helper ==============================

def _plot_quicklook(
    times: List[np.ndarray],
    masses: List[np.ndarray],
    retentate: List[np.ndarray],
    permeate: List[np.ndarray],
    windows: List[Tuple[int, int]],
    save_prefix: Optional[str] = None,
) -> None:
    """
    Make overlay plots using a consistent style:
    - All vials the same line color (black).
    - Vertical dashed lines at vial boundaries.
    - Top x-axis with vial labels placed at the center time of each vial.
    """
    # Boundary x-positions: last time of each vial except the last
    boundary_times = [t[-1] for t in times[:-1] if len(t)]
    # Tick positions for the top axis: mean time within each vial
    centers = [t.mean() for t in times if len(t)]

    # ---------- Plot Δmass vs time ----------
    if any(len(m) > 0 for m in masses):
        fig, ax = plt.subplots()
        for t, m in zip(times, masses):
            mask = ~np.isnan(m)                          # ignore NaNs from spike filter
            if mask.any():
                ax.plot(t[mask], m[mask], color="black") # SAME color for all vials
        for bt in boundary_times:                        # mark vial boundaries
            ax.axvline(bt, color="gray", linestyle=":", linewidth=1)
        ax2 = ax.twiny()                                 # top x-axis (vial labels)
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(centers)
        ax2.set_xticklabels([f"Vial {i+1}" for i in range(len(centers))])
        ax2.set_xlabel("Vials")
        ax.set_xlabel("Time [s]"); ax.set_ylabel("Mass change [g]"); ax.set_title("Mass vs Time")
        fig.tight_layout()
        if save_prefix:
            fig.savefig(f"{save_prefix}_mass.png", dpi=160)

    # ---------- Plot Conductivity vs time ----------
    if any(np.isfinite(arr).any() for arr in retentate) or any(np.isfinite(arr).any() for arr in permeate):
        fig, ax = plt.subplots()
        for t, r in zip(times, retentate):              # retentate = solid black
            if r.size and len(r) == len(t):
                ax.plot(t, r, color="black")
        for t, p in zip(times, permeate):               # permeate = dashed black
            if p.size and len(p) == len(t):
                ax.plot(t, p, linestyle="--", color="black")
        for bt in boundary_times:
            ax.axvline(bt, color="gray", linestyle=":", linewidth=1)
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(centers)
        ax2.set_xticklabels([f"Vial {i+1}" for i in range(len(centers))])
        ax2.set_xlabel("Vials")
        ax.set_xlabel("Time [s]"); ax.set_ylabel("Conductivity [arb]"); ax.set_title("Conductivity vs Time")
        fig.tight_layout()
        if save_prefix:
            fig.savefig(f"{save_prefix}_cond.png", dpi=160)

    plt.show()

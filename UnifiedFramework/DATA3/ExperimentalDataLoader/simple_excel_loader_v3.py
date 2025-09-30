#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 20:29:54 2025

@author: kkasturi
"""

# simple_excel_loader.py
# -------------------------------------------------------------------------------------------------
# Purpose:
#   Unified data loader for Excel (.xlsx/.xls/.xlsm/.xlsb) and MATLAB (.mat) files used in
#   diafiltration experiments. The loader extracts time-series arrays and splits the run into
#   vial-sized windows using either an Index sheet (Start/End or jstart/jend) or a "Vial"/
#   "Vial Swap" marker column embedded in the measurement sheet (or MAT variables).
#
# Public API (after calling .load()):
#   - loader.times            : List[np.ndarray]   # one array per vial/window (time in seconds)
#   - loader.masses           : List[np.ndarray]   # Δmass per vial (mass normalized to window start)
#   - loader.retentate_cond   : List[np.ndarray]   # retentate conductivity per vial (may be NaN)
#   - loader.permeate_cond    : List[np.ndarray]   # permeate conductivity per vial (may be NaN)
#   - loader.windows          : List[(start_idx, end_idx)] for the underlying sheet/vector indices
#
# Plot helper:
#   - loader.plot_quicklook(save_prefix=None)  # creates mass/cond overlay plots; saves if prefix given
#
# Notes:
#   - Matplotlib only (no seaborn). One figure per chart. No explicit colors.
#   - Robust to "offset headers" in Excel (detects the real header row).
#   - For .mat: supports both v7 (scipy.io.loadmat) and v7.3 (HDF5 via h5py) where possible.
#   - Mass spike filter: values of Δmass > 1.2 g are set to NaN (same idea as your draft).
#   - If no windows can be parsed, falls back to a single full-length window.
#
# Usage example:
#   from simple_excel_loader import SimpleExcelLoader
#   loader = SimpleExcelLoader("/path/to/file.xlsx")  # or .mat
#   loader.load()
#   times, masses = loader.times, loader.masses
#   windows = loader.windows
#   loader.plot_quicklook(save_prefix="run1")
# -------------------------------------------------------------------------------------------------

from __future__ import annotations  # for forward type hints on older Python
from typing import List, Tuple, Optional
import os  # filesystem path checks and joins
import numpy as np  # numerical arrays and vector ops
import pandas as pd  # tabular parsing (Excel)
import matplotlib.pyplot as plt  # plotting (one chart per figure)

class SimpleExcelLoader:
    """
    A single class that handles both Excel and .mat files while exposing the same API.

    Constructor:
        SimpleExcelLoader(path: str)

    Important public attributes (populated after .load()):
        times, masses, retentate_cond, permeate_cond, windows
    """

    def __init__(self, path: str):
        # Store the path to the input file (.xlsx/.xls/.xlsm/.xlsb or .mat)
        self.path: str = path

        # Outputs (public):
        self.times: List[np.ndarray] = []            # time arrays per vial/window
        self.masses: List[np.ndarray] = []           # Δmass arrays per vial/window
        self.retentate_cond: List[np.ndarray] = []   # retentate conductivity arrays
        self.permeate_cond: List[np.ndarray] = []    # permeate conductivity arrays
        self.windows: List[Tuple[int, int]] = []     # index windows (start_idx, end_idx)

        # Internals (Excel only; useful for debugging/inspection):
        self._meas_df: Optional[pd.DataFrame] = None  # parsed measurement DataFrame
        self.measure_sheet: Optional[str] = None      # selected measurement sheet name
        self.index_sheet: Optional[str] = None        # selected index sheet name

    # ------------------------------- PUBLIC ENTRY POINT -----------------------------------------

    def load(self) -> None:
        """Detect file type by extension and route to the appropriate loader."""
        # Normalize extension to lowercase for robust comparison
        ext = os.path.splitext(self.path)[1].lower()

        # Validate the file exists early (clear error if not found)
        if not os.path.exists(self.path):
            raise FileNotFoundError(f"File not found: {self.path}")

        # Route by extension
        if ext in {".xlsx", ".xls", ".xlsm", ".xlsb"}:
            self._load_excel()  # Excel pathway
        elif ext == ".mat":
            self._load_mat()    # MATLAB pathway
        else:
            # If the extension is unexpected, provide a helpful error
            raise ValueError(f"Unsupported file type: {ext} (expected Excel or .mat)")

    # ------------------------------------ EXCEL PATHWAY ----------------------------------------

    def _load_excel(self) -> None:
        """Load from an Excel workbook, robust to offset headers, and compute vial windows."""

        # Helper: detect which sheet looks like the measurement sheet.
        def _detect_measurement_sheet(xl: pd.ExcelFile) -> Optional[str]:
            # Iterate over all sheets
            for sh in xl.sheet_names:
                try:
                    # Read only a few rows to sniff columns; faster than full parse
                    df = xl.parse(sh, nrows=10)
                except Exception:
                    # If a sheet fails to parse, skip it
                    continue
                # Lowercased column names (as strings) for robust matching
                cols = [str(c).lower() for c in df.columns]
                # Heuristic: measurement sheet has "time" and either "mass" or "cond" somewhere
                if any("time" in c for c in cols) and (any("mass" in c for c in cols) or any("cond" in c for c in cols)):
                    return sh
            # Fallback: use the first sheet if nothing matches
            return xl.sheet_names[0] if xl.sheet_names else None

        # Helper: detect a likely index sheet (contains Start/End or jstart/jend)
        def _detect_index_sheet(xl: pd.ExcelFile) -> Optional[str]:
            # Check each sheet with a quick sniff
            for sh in xl.sheet_names:
                try:
                    df = xl.parse(sh, nrows=10)
                except Exception:
                    continue
                cols = [str(c).lower() for c in df.columns]
                # Signals that an index sheet might be present
                if any("start" in c for c in cols) and any("end" in c for c in cols):
                    return sh
                if any("jstart" in c for c in cols) and any("jend" in c for c in cols):
                    return sh
            return None

        # Helper: find which row is actually the header row (handles offset headers)
        def _detect_header_row(df: pd.DataFrame, max_rows: int = 10) -> int:
            # Look only in the top 'max_rows' rows
            for i in range(min(max_rows, len(df))):
                # Coerce row values to strings, lowercased, for robust substring search
                row = df.iloc[i].astype(str).str.lower().tolist()
                # A good header row includes 'time' and also 'mass' or 'cond' or 'pressure'
                has_time = any("time" in s for s in row)
                has_mass_or_cond = any(("mass" in s) or ("cond" in s) or ("pressure" in s) for s in row)
                if has_time and has_mass_or_cond:
                    # Return the index of the detected header row
                    return i
            # If not found, assume the very first row is the header
            return 0

        # Helper: find the first column whose name contains a needle (case-insensitive)
        def _first_col_like(df: pd.DataFrame, needle: str):
            n = needle.lower()
            for c in df.columns:
                if n in str(c).lower():
                    return c
            return None

        # Parse the Excel file using pandas
        xl = pd.ExcelFile(self.path)  # engine is auto-selected by pandas based on extension

        # Figure out which sheets to use
        self.measure_sheet = _detect_measurement_sheet(xl)  # measurement sheet
        self.index_sheet = _detect_index_sheet(xl)          # optional index sheet

        # Safety check: we must have a measurement sheet to proceed
        if self.measure_sheet is None:
            raise ValueError("No measurement-like sheet (Time + Mass/Cond) found in Excel file.")

        # Read the entire measurement sheet with no header first (so we can detect it)
        raw = xl.parse(self.measure_sheet, header=None)

        # Detect which row is the actual header row
        hdr = _detect_header_row(raw, max_rows=10)

        # Build a DataFrame with the detected header (row 'hdr'), and data below it
        df = raw.iloc[hdr + 1:].copy()                             # take all rows after the header
        df.columns = raw.iloc[hdr].astype(str).str.strip().tolist() # set column names from header row
        df = df.reset_index(drop=True)                              # clean index to 0..N-1
        self._meas_df = df                                          # store for inspection/debugging

        # Identify key columns by substring matching (robust to naming variations)
        tcol = _first_col_like(df, "time")      # time column
        mcol = _first_col_like(df, "mass")      # mass column
        # Conductivity columns: look for retentate/permeate + 'cond'
        rcol = next((c for c in df.columns if ("retent" in str(c).lower() and "cond" in str(c).lower())), None)
        pcol = next((c for c in df.columns if ("permeate" in str(c).lower() and "cond" in str(c).lower())), None)

        # Compute vial windows:
        # 1) Prefer index sheet if present
        windows: List[Tuple[int, int]] = []
        if self.index_sheet is not None:
            idx_df = xl.parse(self.index_sheet)  # read full index sheet
            c_start = _first_col_like(idx_df, "start") or _first_col_like(idx_df, "jstart")  # start index column
            c_end   = _first_col_like(idx_df, "end")   or _first_col_like(idx_df, "jend")    # end index column
            if c_start is None or c_end is None:
                # Fallback: if there are at least 3 columns, try the 2nd and 3rd as start/end
                if idx_df.shape[1] >= 3:
                    c_start = idx_df.columns[1]
                    c_end = idx_df.columns[2]
            if c_start is not None and c_end is not None:
                # Convert to numeric arrays; drop non-numeric values
                starts = pd.to_numeric(idx_df[c_start], errors="coerce").dropna().astype(int).to_numpy()
                ends   = pd.to_numeric(idx_df[c_end], errors="coerce").dropna().astype(int).to_numpy()
                # Build (start, end) tuples, ensuring shape match and end >= start
                n = min(len(starts), len(ends))
                windows = [(int(starts[i]), int(ends[i])) for i in range(n) if ends[i] >= starts[i]]

        # 2) If no index sheet windows, try to build windows from a 'Vial' / 'Vial Swap' column
        if not windows:
            vcol = _first_col_like(df, "vial")  # find a column with 'vial' in its name
            if vcol is not None:
                v = pd.to_numeric(df[vcol], errors="coerce").fillna(-1).to_numpy(dtype=float)  # coerce to numeric
                # If there are finite values and the value changes, we can segment on changes
                if np.isfinite(v).any() and (np.diff(v[np.isfinite(v)]) != 0).any():
                    # Build boundaries where the vial marker changes value
                    boundaries = [0]
                    for i in range(1, len(v)):
                        if v[i] != v[i - 1]:
                            boundaries.append(i)
                    boundaries.append(len(v))
                    # Turn boundaries into [start, end] windows; ignore tiny 1- or 2-row fragments
                    for s, e in zip(boundaries[:-1], boundaries[1:]):
                        if e - s >= 3:  # require at least 3 rows to avoid swap spikes
                            windows.append((int(s), int(e - 1)))

        # 3) If we still have no windows, use the entire data as one window
        if not windows and len(df) > 0:
            windows = [(0, len(df) - 1)]

        # Persist windows to the public attribute
        self.windows = windows

        # Helper to convert a series to a numeric NumPy array
        def to_num(series):
            return pd.to_numeric(series, errors="coerce").to_numpy()

        # Convert columns to numeric arrays, with sensible fallbacks (NaNs) if missing
        tvals = to_num(df[tcol]) if tcol else np.arange(len(df), dtype=float)  # time defaults to 0..N-1 if absent
        mvals = to_num(df[mcol]) if mcol else np.full(len(df), np.nan)         # mass defaults to NaN if absent
        rvals = to_num(df[rcol]) if rcol else np.full(len(df), np.nan)         # retentate cond defaults to NaN
        pvals = to_num(df[pcol]) if pcol else np.full(len(df), np.nan)         # permeate cond defaults to NaN

        # Build per-window arrays (time, Δmass, retentate, permeate)
        self.times, self.masses, self.retentate_cond, self.permeate_cond = [], [], [], []
        for (s, e) in self.windows:
            sl = slice(s, e + 1)          # slice for this window
            t = tvals[sl]                  # time segment
            mraw = mvals[sl]               # raw mass segment (absolute)
            # Normalize mass to Δmass per window: subtract the first mass value if it is valid
            m = mraw - (mraw[0] if len(mraw) and not np.isnan(mraw[0]) else 0.0)
            # Apply spike filter: Δmass > 1.2 g are considered invalid and set to NaN
            m = np.where(m > 1.2, np.nan, m)
            r = rvals[sl]                  # retentate cond segment
            p = pvals[sl]                  # permeate cond segment
            # Append segments to the outputs
            self.times.append(t)
            self.masses.append(m)
            self.retentate_cond.append(r)
            self.permeate_cond.append(p)

    # ------------------------------------- MATLAB PATHWAY ---------------------------------------

    def _load_mat(self) -> None:
        """Load from a MATLAB .mat file (supports v7 via scipy.io.loadmat and v7.3 via h5py)."""

        # Try: scipy.io.loadmat (works for v7 MAT files)
        def _try_scipy_loadmat(path: str):
            try:
                from scipy.io import loadmat  # import locally to keep dependency optional
                return loadmat(path, squeeze_me=True, struct_as_record=False)
            except Exception:
                return None  # return None to indicate failure

        # Try: h5py (works for v7.3 MAT files which are HDF5-based)
        def _try_h5py(path: str):
            try:
                import h5py  # import locally to keep dependency optional
                return h5py.File(path, "r")
            except Exception:
                return None  # return None to indicate failure

        # Attempt scipy first (v7)
        mat = _try_scipy_loadmat(self.path)

        if mat is None:
            # If scipy failed, attempt h5py (v7.3 HDF5)
            h5 = _try_h5py(self.path)
            if h5 is None:
                # Neither loader could open the file
                raise ValueError("Could not read .mat: neither scipy.io.loadmat (v7) nor h5py (v7.3) succeeded.")

            # Helper to read a dataset if it exists in HDF5
            def _h5_read(name: str):
                # Return a squeezed NumPy array if path exists, else None
                return np.array(h5[name]).squeeze() if name in h5 else None

            # Try common key spellings (adjust these if your MAT keys differ)
            time = _h5_read("Time") or _h5_read("time")
            mass = _h5_read("Mass") or _h5_read("mass")
            rcon = _h5_read("RetentateCond") or _h5_read("Retentate_Cond") or _h5_read("retentate_cond")
            pcon = _h5_read("PermeateCond") or _h5_read("Permeate_Cond") or _h5_read("permeate_cond")
            vial = _h5_read("Vial") or _h5_read("VialSwap") or _h5_read("vial") or _h5_read("vial_swap")

            # Close the HDF5 file to free resources
            h5.close()

        else:
            # scipy path (v7): mat is a dict-like structure of arrays and possibly structs
            # Helper: find a key in a dict regardless of capitalization
            def _find_key(d, *candidates):
                keys = {k.lower(): k for k in d.keys()}
                for cand in candidates:
                    if cand.lower() in keys:
                        return d[keys[cand.lower()]]
                return None

            # Some MAT files store everything inside a top-level struct (e.g., data_stru)
            top = _find_key(mat, "data_stru", "data", "data_struct", "S", "D")

            # Choose the source dict-like object (struct fields vs. flat dict)
            src = top.__dict__ if hasattr(top, "__dict__") else mat

            # Helper: try several candidate names for a variable and return the first found
            def _grab(*names):
                for nm in names:
                    v = _find_key(src, nm)
                    if v is not None:
                        return np.array(v).squeeze()
                return None

            # Pull common signals with flexible naming
            time = _grab("Time", "time", "t")
            mass = _grab("Mass", "mass", "m")
            rcon = _grab("RetentateCond", "Retentate_Cond", "retentate_cond", "ret_cond")
            pcon = _grab("PermeateCond", "Permeate_Cond", "permeate_cond", "perm_cond")
            vial = _grab("Vial", "VialSwap", "vial", "vial_swap")

        # At minimum, we require time and mass to exist
        if time is None or mass is None:
            raise ValueError("MAT file missing required arrays: Time and Mass.")

        # Ensure arrays are 1D and have consistent length; pad/truncate aux arrays to match time length
        n = len(np.array(time).squeeze())  # canonical length from time
        def _fit_len(x):
            # Helper: return a 1D array of length n; if x is None, return NaNs
            if x is None:
                return np.full(n, np.nan)
            arr = np.array(x).squeeze()
            # If longer than n: truncate; if shorter: pad with NaNs
            return arr[:n] if arr.size >= n else np.pad(arr, (0, n - arr.size), constant_values=np.nan)

        time = _fit_len(time).astype(float)  # time array (float)
        mass = _fit_len(mass).astype(float)  # mass array (float)
        rcon = _fit_len(rcon).astype(float)  # retentate cond array (float, may be NaN)
        pcon = _fit_len(pcon).astype(float)  # permeate cond array (float, may be NaN)

        # Build windows from a vial marker if available; else use a single full window
        windows: List[Tuple[int, int]] = []
        if vial is not None:
            v = np.array(vial).squeeze()
            try:
                # Try numeric coercion first (works for 0/1, 1/2, etc.)
                import pandas as pd  # safe local import, pandas already available
                vnum = pd.to_numeric(pd.Series(v), errors="coerce").fillna(-1).to_numpy()
            except Exception:
                # As a fallback, if labels are non-numeric, split on changes in string labels
                v_str = np.array(list(map(str, v)))
                vnum = (v_str != np.roll(v_str, 1)).astype(int).cumsum()

            # Identify change points in the vial stream
            boundaries = [0] + [i for i in range(1, len(vnum)) if vnum[i] != vnum[i - 1]] + [len(vnum)]
            # Convert boundaries into windows; discard tiny fragments (<3 rows)
            for s, e in zip(boundaries[:-1], boundaries[1:]):
                if e - s >= 3:
                    windows.append((int(s), int(e - 1)))

        # If still no windows, use a single full-length window
        if not windows:
            windows = [(0, n - 1)]

        # Persist windows
        self.windows = windows

        # Build outputs per window (same semantics as Excel)
        self.times, self.masses, self.retentate_cond, self.permeate_cond = [], [], [], []
        for s, e in self.windows:
            sl = slice(s, e + 1)             # window slice
            t = time[sl]                      # time segment
            mraw = mass[sl]                   # raw mass segment
            m = mraw - (mraw[0] if len(mraw) and not np.isnan(mraw[0]) else 0.0)  # Δmass
            m = np.where(m > 1.2, np.nan, m)  # spike filter like your draft
            r = rcon[sl]                      # retentate cond segment
            p = pcon[sl]                      # permeate cond segment
            self.times.append(t)
            self.masses.append(m)
            self.retentate_cond.append(r)
            self.permeate_cond.append(p)

    # ------------------------------------------ PLOTS -------------------------------------------

    # def plot_quicklook(self, save_prefix: Optional[str] = None) -> None:
    #     """
    #     Make overlay plots of Δmass vs time and conductivity vs time (one figure per chart).
    #     If save_prefix is provided, saves PNGs to f"{save_prefix}_mass.png" and "_cond.png".
    #     """
    #     # Plot 1: Δmass vs time (overlay all vials/windows)
    #     if any(len(m) > 0 for m in self.masses):
    #         plt.figure()  # create a new figure (no subplots per requirements)
    #         for i, (t, m) in enumerate(zip(self.times, self.masses)):
    #             mask = ~np.isnan(m)  # only plot finite Δmass points
    #             if mask.any():
    #                 plt.plot(t[mask], m[mask], label=f"Vial {i}")
    #         plt.xlabel("Time [s]")           # x-axis label
    #         plt.ylabel("Mass change [g]")    # y-axis label
    #         plt.title("Mass vs Time (by vial)")  # plot title
    #         plt.legend(loc="best")           # legend placement
    #         plt.tight_layout()               # clean layout to avoid clipping
    #         if save_prefix:
    #             # Save figure if requested
    #             plt.savefig(f"{save_prefix}_mass.png", dpi=160)

    #     # Plot 2: conductivity vs time (retentate solid, permeate dashed)
    #     if any(np.isfinite(arr).any() for arr in self.retentate_cond) or any(np.isfinite(arr).any() for arr in self.permeate_cond):
    #         plt.figure()  # new figure
    #         # Plot retentate in solid lines
    #         for i, (t, r) in enumerate(zip(self.times, self.retentate_cond)):
    #             if r.size and len(r) == len(t):
    #                 plt.plot(t, r, label=f"Retentate v{i}")
    #         # Plot permeate in dashed lines
    #         for i, (t, p) in enumerate(zip(self.times, self.permeate_cond)):
    #             if p.size and len(p) == len(t):
    #                 plt.plot(t, p, linestyle="--", label=f"Permeate v{i}")
    #         plt.xlabel("Time [s]")
    #         plt.ylabel("Conductivity [arb]")
    #         plt.title("Conductivity vs Time")
    #         plt.legend(loc="best")
    #         plt.tight_layout()
    #         if save_prefix:
    #             plt.savefig(f"{save_prefix}_cond.png", dpi=160)

    #     # Show all generated figures
    #     plt.show()

    def plot_quicklook(self, save_prefix: Optional[str] = None) -> None:
        """
        Make overlay plots of Δmass vs time and conductivity vs time (one figure per chart).
        Adds a top x-axis with vial numbers and vertical lines at vial boundaries.
        All vials use the same color for clarity.
        """

        # ---------------- Collect vial boundary times ----------------
        # We'll mark where each vial ends (last time point in the window),
        # except the very last vial (no boundary after it).
        boundary_times = []
        for t in self.times[:-1]:        # iterate over all but last vial
            if len(t):                   # make sure the vial has data
                boundary_times.append(t[-1])  # take the last time of that vial

        # ---------------- Plot 1: Δmass vs time ----------------
        if any(len(m) > 0 for m in self.masses):   # only plot if at least one vial has mass data
            fig, ax = plt.subplots()               # create a figure + main axis

            # Plot each vial’s Δmass vs time
            # All vials get the same color (black) for consistency.
            for t, m in zip(self.times, self.masses):
                mask = ~np.isnan(m)                # only plot valid (non-NaN) values
                if mask.any():
                    ax.plot(t[mask], m[mask], color="red")  # plot line in black

            # Set axis labels and title
            ax.set_xlabel("Time [s]")              # bottom x-axis label
            ax.set_ylabel("Mass change [g]")       # y-axis label
            ax.set_title("Mass vs Time")           # title

            # Draw vertical lines at vial boundaries
            for bt in boundary_times:
                ax.axvline(bt, color="gray", linestyle=":", linewidth=1)

            # Add a second x-axis on top showing vial numbers
            ax2 = ax.twiny()                       # create a new top x-axis sharing the same data scale
            ax2.set_xlim(ax.get_xlim())            # align limits with bottom axis
            centers = [t.mean() for t in self.times if len(t)]   # compute average time of each vial
            ax2.set_xticks(centers)                # place ticks at vial centers
            ax2.set_xticklabels([f"Vial {i}" for i in range(len(centers))])  # label them
            ax2.set_xlabel("Vials")                # top axis label

            # Adjust layout and save figure if requested
            fig.tight_layout()
            if save_prefix:
                fig.savefig(f"{save_prefix}_mass.png", dpi=160)

        # ---------------- Plot 2: Conductivity vs time ----------------
        # Only plot if either retentate or permeate has any finite values
        if any(np.isfinite(arr).any() for arr in self.retentate_cond) or any(np.isfinite(arr).any() for arr in self.permeate_cond):
            fig, ax = plt.subplots()               # create figure + main axis

            # Plot retentate curves (all solid black)
            for t, r in zip(self.times, self.retentate_cond):
                if r.size and len(r) == len(t):    # check length match
                    ax.plot(t, r, color="black")

            # Plot permeate curves (all dashed black)
            for t, p in zip(self.times, self.permeate_cond):
                if p.size and len(p) == len(t):    # check length match
                    ax.plot(t, p, linestyle="--", color="black")

            # Axis labels and title
            ax.set_xlabel("Time [s]")
            ax.set_ylabel("Conductivity [arb]")
            ax.set_title("Conductivity vs Time")

            # Add vertical boundary lines
            for bt in boundary_times:
                ax.axvline(bt, color="gray", linestyle=":", linewidth=1)

            # Add top x-axis with vial numbers
            ax2 = ax.twiny()
            ax2.set_xlim(ax.get_xlim())
            centers = [t.mean() for t in self.times if len(t)]
            ax2.set_xticks(centers)
            ax2.set_xticklabels([f"Vial {i}" for i in range(len(centers))])
            ax2.set_xlabel("Vials")

            # Layout + optional save
            fig.tight_layout()
            if save_prefix:
                fig.savefig(f"{save_prefix}_cond.png", dpi=160)

        # Finally, show all generated figures on screen
        plt.show()


"""
simple_excel_loader.py (v2: vial splitting)
-------------------------------------------
Excel → variables → plots for NF270-style workbooks, robust to "offset headers"
and supports vial segmentation from either:
  - an Index sheet with Start/End (or jstart/jend), OR
  - a 'Vial' / 'Vial Swap' marker column in the measurement sheet (contiguous segments).

Usage:
    from simple_excel_loader import SimpleExcelLoader
    loader = SimpleExcelLoader("NF270_MC3.xlsx")
    loader.load()
    # Access variables:
    times = loader.times              # List[np.ndarray]
    masses = loader.masses            # List[np.ndarray]  (Δmass per vial)
    retentate = loader.retentate_cond # List[np.ndarray]
    permeate = loader.permeate_cond   # List[np.ndarray]
    windows = loader.windows          # List[(start_idx, end_idx)] index windows in the sheet

    # Quicklook plots:
    loader.plot_quicklook(save_prefix="NF270_MC3_vialsplit")

Notes:
- Uses matplotlib (no seaborn).
- Each chart is a single plot; no explicit colors set.
"""

from typing import Any, Dict, List, Optional, Tuple
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _first_col_like(df: pd.DataFrame, needle: str):
    n = needle.lower()
    for c in df.columns:
        if n in str(c).lower():
            return c
    return None

def _detect_measurement_sheet(xl: pd.ExcelFile) -> Optional[str]:
    for sh in xl.sheet_names:
        try:
            df = xl.parse(sh, nrows=10)
        except Exception:
            continue
        cols = [str(c).lower() for c in df.columns]
        if any("time" in c for c in cols) and (any("mass" in c for c in cols) or any("cond" in c for c in cols)):
            return sh
    return xl.sheet_names[0] if xl.sheet_names else None

def _detect_index_sheet(xl: pd.ExcelFile) -> Optional[str]:
    for sh in xl.sheet_names:
        try:
            df = xl.parse(sh, nrows=10)
        except Exception:
            continue
        cols = [str(c).lower() for c in df.columns]
        if any("start" in c for c in cols) and any("end" in c for c in cols):
            return sh
        if any("jstart" in c for c in cols) and any("jend" in c for c in cols):
            return sh
    return None

def _detect_header_row(df: pd.DataFrame, max_rows: int = 10) -> int:
    for i in range(min(max_rows, len(df))):
        row = df.iloc[i].astype(str).str.lower().tolist()
        has_time = any("time" in s for s in row)
        has_mass_or_cond = any(("mass" in s) or ("cond" in s) or ("pressure" in s) for s in row)
        if has_time and has_mass_or_cond:
            return i
    return 0


class SimpleExcelLoader:
    def __init__(self, path: str):
        self.path = path
        self.xl: Optional[pd.ExcelFile] = None
        self.measure_sheet: Optional[str] = None
        self.index_sheet: Optional[str] = None

        # Saved variables
        self.times: List[np.ndarray] = []
        self.masses: List[np.ndarray] = []
        self.retentate_cond: List[np.ndarray] = []
        self.permeate_cond: List[np.ndarray] = []
        self.windows: List[Tuple[int,int]] = []

        # Internals
        self._meas_df: Optional[pd.DataFrame] = None

    def load(self) -> None:
        if not os.path.exists(self.path):
            raise FileNotFoundError(f"Excel file not found: {self.path}")
        self.xl = pd.ExcelFile(self.path)

        # Detect sheets
        self.measure_sheet = _detect_measurement_sheet(self.xl)
        self.index_sheet = _detect_index_sheet(self.xl)

        # Read measurement with header-row detection
        if self.measure_sheet is None:
            raise ValueError("Could not detect a measurement-like sheet (with Time + Mass/Conductivity).")
        raw = self.xl.parse(self.measure_sheet, header=None)
        hdr = _detect_header_row(raw, max_rows=10)
        df = raw.iloc[hdr+1:].copy()
        df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()
        df = df.reset_index(drop=True)
        self._meas_df = df

        # Identify relevant columns
        tcol = _first_col_like(df, "time")
        mcol = _first_col_like(df, "mass")
        # Retentate & permeate conductivity columns
        rcol = None
        for c in df.columns:
            cl = str(c).lower()
            if "retent" in cl and "cond" in cl:
                rcol = c; break
        pcol = None
        for c in df.columns:
            cl = str(c).lower()
            if "permeate" in cl and "cond" in cl:
                pcol = c; break

        # Decide windows (vials)
        windows: List[Tuple[int,int]] = []
        if self.index_sheet is not None:
            idx_df = self.xl.parse(self.index_sheet)
            c_start = _first_col_like(idx_df, "start") or _first_col_like(idx_df, "jstart")
            c_end   = _first_col_like(idx_df, "end")   or _first_col_like(idx_df, "jend")
            if c_start is None or c_end is None:
                if idx_df.shape[1] >= 3:
                    c_start = idx_df.columns[1]
                    c_end = idx_df.columns[2]
            if c_start is not None and c_end is not None:
                starts = pd.to_numeric(idx_df[c_start], errors="coerce").dropna().astype(int).to_numpy()
                ends   = pd.to_numeric(idx_df[c_end], errors="coerce").dropna().astype(int).to_numpy()
                n = min(len(starts), len(ends))
                windows = [(int(starts[i]), int(ends[i])) for i in range(n) if ends[i] >= starts[i]]

        if not windows:
            # Try a 'Vial' marker column in the measurement sheet
            vcol = _first_col_like(df, "vial")
            if vcol is not None:
                v = pd.to_numeric(df[vcol], errors="coerce").fillna(-1).to_numpy(dtype=float)
                # If it's all NaN or constant, fallback to full range
                if np.isfinite(v).any() and (np.diff(v[np.isfinite(v)]) != 0).any():
                    # Build boundaries when value changes; skip 1-row spikes (swap rows)
                    boundaries = [0]
                    for i in range(1, len(v)):
                        if v[i] != v[i-1]:
                            boundaries.append(i)
                    boundaries.append(len(v))
                    for s, e in zip(boundaries[:-1], boundaries[1:]):
                        # Exclude single-row swap markers (e-s==1) and require length>=3
                        if e - s >= 3:
                            windows.append((int(s), int(e-1)))

        if not windows:
            # Full range as single window
            if len(df) > 0:
                windows = [(0, len(df)-1)]

        self.windows = windows

        # Extract arrays per window
        def to_num(series):
            return pd.to_numeric(series, errors="coerce").to_numpy()

        tvals = to_num(df[tcol]) if tcol else np.arange(len(df), dtype=float)
        mvals = to_num(df[mcol]) if mcol else np.full(len(df), np.nan)
        rvals = to_num(df[rcol]) if rcol else np.full(len(df), np.nan)
        pvals = to_num(df[pcol]) if pcol else np.full(len(df), np.nan)

        self.times, self.masses, self.retentate_cond, self.permeate_cond = [], [], [], []
        for (s, e) in self.windows:
            sl = slice(s, e+1)
            t = tvals[sl]
            mraw = mvals[sl]
            m = mraw - (mraw[0] if len(mraw) and not np.isnan(mraw[0]) else 0.0)
            m = np.where(m > 1.2, np.nan, m)  # spike filter per your draft
            r = rvals[sl]
            p = pvals[sl]
            self.times.append(t); self.masses.append(m); self.retentate_cond.append(r); self.permeate_cond.append(p)

    def quick_summary(self) -> str:
        lines = []
        lines.append(f"Workbook: {os.path.basename(self.path)}")
        lines.append(f"Measurement sheet: {self.measure_sheet}")
        lines.append(f"Index sheet: {self.index_sheet}")
        lines.append(f"# vials/windows: {len(self.windows)}")
        if self.times:
            lines.append("Lengths: " + ", ".join(str(len(t)) for t in self.times[:10]))
        return "\n".join(lines)

    def plot_quicklook(self, save_prefix: Optional[str] = None) -> None:
        if any(len(m)>0 for m in self.masses):
            plt.figure()
            for i, (t, m) in enumerate(zip(self.times, self.masses)):
                mask = ~np.isnan(m)
                if mask.any():
                    plt.plot(t[mask], m[mask], label=f"Vial {i}")
            plt.xlabel("Time [s]"); plt.ylabel("Mass change [g]")
            plt.title("Mass vs Time (by vial)"); plt.legend(loc="best"); plt.tight_layout()
            if save_prefix:
                plt.savefig(f"{save_prefix}_mass.png", dpi=160)
        if any(np.isfinite(arr).any() for arr in self.retentate_cond) or any(np.isfinite(arr).any() for arr in self.permeate_cond):
            plt.figure()
            for i, (t, r) in enumerate(zip(self.times, self.retentate_cond)):
                if r.size and len(r) == len(t):
                    plt.plot(t, r, label=f"Retentate v{i}")
            for i, (t, p) in enumerate(zip(self.times, self.permeate_cond)):
                if p.size and len(p) == len(t):
                    plt.plot(t, p, linestyle="--", label=f"Permeate v{i}")
            plt.xlabel("Time [s]"); plt.ylabel("Conductivity [arb]")
            plt.title("Conductivity vs Time"); plt.legend(loc="best"); plt.tight_layout()
            if save_prefix:
                plt.savefig(f"{save_prefix}_cond.png", dpi=160)
        plt.show()

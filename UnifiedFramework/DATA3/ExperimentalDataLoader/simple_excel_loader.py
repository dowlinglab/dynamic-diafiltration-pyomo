
"""
simple_excel_loader.py
----------------------
Excel → variables → plots for NF270-style workbooks, robust to "offset headers" and
sheet naming differences.

Usage:
    from simple_excel_loader import SimpleExcelLoader
    loader = SimpleExcelLoader("NF270_MC3.xlsx")
    loader.load()
    # Saved variables
    times = loader.times            # List[np.ndarray]
    masses = loader.masses          # List[np.ndarray]
    retentate = loader.retentate_cond
    permeate = loader.permeate_cond
    # Quick plots
    loader.plot_quicklook(save_prefix="NF270_MC3_quicklook")

Constraints:
- Uses matplotlib (no seaborn).
- One figure per chart; no custom colors.
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
    """
    Find a row within the first max_rows that looks like a header row:
    containing tokens like 'time' and ('mass' or 'cond' or 'pressure').
    Returns row index (int). If not found, returns 0.
    """
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
        self.config: Dict[str, Any] = {}

        # Internals
        self._meas_df: Optional[pd.DataFrame] = None

    def load(self) -> None:
        if not os.path.exists(self.path):
            raise FileNotFoundError(f"Excel file not found: {self.path}")
        self.xl = pd.ExcelFile(self.path)

        # Sheets
        self.measure_sheet = _detect_measurement_sheet(self.xl)
        self.index_sheet = _detect_index_sheet(self.xl)

        # Read measurement sheet with robust header handling
        if self.measure_sheet is None:
            raise ValueError("Could not detect a measurement-like sheet (with Time + Mass/Conductivity).")
        raw = self.xl.parse(self.measure_sheet, header=None)  # inspect to find header row
        hdr = _detect_header_row(raw, max_rows=10)
        df = raw.iloc[hdr+1:].copy()
        df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()
        df = df.reset_index(drop=True)
        self._meas_df = df

        # Identify columns
        tcol = _first_col_like(df, "time")
        mcol = _first_col_like(df, "mass")
        # Conductivity columns
        ret_cand, per_cand, generic_cond = None, None, None
        for c in df.columns:
            cl = str(c).lower()
            if "cond" in cl and "retent" in cl:
                ret_cand = c
            elif "cond" in cl and ("perm" in cl or "permea" in cl):
                per_cand = c
            elif "cond" in cl and generic_cond is None:
                generic_cond = c
        if ret_cand is None and generic_cond is not None:
            ret_cand = generic_cond

        # Windows (index sheet) or single window
        windows: List[Tuple[int, int]] = []
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
            nrows = len(df)
            if nrows > 0:
                windows = [(0, nrows - 1)]

        # Extract arrays
        tvals = pd.to_numeric(df[tcol], errors="coerce").to_numpy() if tcol else None
        mvals = pd.to_numeric(df[mcol], errors="coerce").to_numpy() if mcol else None
        rvals = pd.to_numeric(df[ret_cand], errors="coerce").to_numpy() if ret_cand else None
        pvals = pd.to_numeric(df[per_cand], errors="coerce").to_numpy() if per_cand else None

        self.times, self.masses, self.retentate_cond, self.permeate_cond = [], [], [], []
        for (jstart, jend) in windows:
            sl = slice(jstart, jend + 1)
            t = tvals[sl] if tvals is not None else np.arange(jend - jstart + 1, dtype=float)
            if mvals is not None:
                mraw = mvals[sl].astype(float)
                if len(mraw) and not np.isnan(mraw[0]):
                    m = mraw - mraw[0]
                else:
                    m = mraw
                m = np.where(m > 1.2, np.nan, m)  # spike filter as in draft
            else:
                m = np.array([], dtype=float)
            r = rvals[sl].astype(float) if rvals is not None else np.array([], dtype=float)
            p = pvals[sl].astype(float) if pvals is not None else np.array([], dtype=float)
            self.times.append(t); self.masses.append(m); self.retentate_cond.append(r); self.permeate_cond.append(p)

    def quick_summary(self) -> str:
        lines = []
        lines.append(f"Workbook: {os.path.basename(self.path)}")
        lines.append(f"Measurement sheet: {self.measure_sheet}")
        lines.append(f"Index sheet: {self.index_sheet}")
        lines.append(f"# vials/windows: {len(self.times)}")
        if self.times:
            lines.append("Per vial lengths (first 5): " + ', '.join(str(len(t)) for t in self.times[:5]))
        lines.append(f"Mass available: {any(len(m)>0 for m in self.masses)}")
        lines.append(f"Retentate conductivity available: {any(len(r)>0 for r in self.retentate_cond)}")
        lines.append(f"Permeate conductivity available: {any(len(p)>0 for p in self.permeate_cond)}")
        return '\n'.join(lines)

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
        if any(len(r)>0 for r in self.retentate_cond) or any(len(p)>0 for p in self.permeate_cond):
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

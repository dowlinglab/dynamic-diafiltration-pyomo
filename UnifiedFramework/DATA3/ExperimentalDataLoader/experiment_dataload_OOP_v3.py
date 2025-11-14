#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OOP-structured unified data loader for Excel (.xlsx/.xls/.xlsm/.xlsb) and MATLAB (.mat)
used in diafiltration experiments, with unified plotting that:
- Mass vs time: single solid color for all vials + legend.
- Conductivity vs time: retentate solid (one color), permeate dashed (one color) + legend.
- Top x-axis vial labels and vertical lines at vial boundaries.

Provides:
- DataLoader(path, *, measure_sheet=None, index_sheet=None)
    .load()
    .times, .masses, .retentate_cond, .permeate_cond, .windows
    .plot_quicklook(save_prefix=None)
- list_sheets(path) -> List[str]
- CLI: python ooph_structured_data_loader.py <file> [--quicklook] [--save-prefix PREFIX]
       [--list-sheets] [--dump-windows] [--export-csv OUTDIR]
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple, Optional, Protocol, runtime_checkable
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------- Data container -----------------------------------------------

@dataclass
class DataRun:
    time: np.ndarray
    mass: np.ndarray
    retentate_cond: np.ndarray
    permeate_cond: np.ndarray
    vial_marker: Optional[np.ndarray] = None

    def __post_init__(self):
        n = int(np.asarray(self.time).size)
        def _fit(x, fill=np.nan):
            if x is None:
                return np.full(n, fill)
            arr = np.asarray(x).squeeze()
            arr = arr[:n] if arr.size >= n else np.pad(arr, (0, n - arr.size), constant_values=fill)
            return arr
        self.time = _fit(self.time).astype(float)
        self.mass = _fit(self.mass).astype(float)
        self.retentate_cond = _fit(self.retentate_cond).astype(float)
        self.permeate_cond = _fit(self.permeate_cond).astype(float)
        if self.vial_marker is not None:
            vm = np.asarray(self.vial_marker).squeeze()
            self.vial_marker = vm[:n] if vm.size >= n else np.pad(vm, (0, n - vm.size), constant_values=np.nan)

# -------------------------------- Data sources --------------------------------------------------

@runtime_checkable
class DataSource(Protocol):
    def sniff(self) -> bool: ...
    def parse(self) -> DataRun: ...

class ExcelSource:
    """Parses Excel workbooks with offset headers and optional index sheet."""
    def __init__(self, path: str, measure_sheet_override: Optional[str] = None, index_sheet_override: Optional[str] = None):
        self.path = path
        self._measure_sheet_override = measure_sheet_override
        self._index_sheet_override = index_sheet_override
        self.measure_sheet: Optional[str] = None
        self.index_sheet: Optional[str] = None
        self._meas_df: Optional[pd.DataFrame] = None

    def sniff(self) -> bool:
        return os.path.splitext(self.path)[1].lower() in {".xlsx", ".xls", ".xlsm", ".xlsb"}

    @staticmethod
    def _first_col_like(df: pd.DataFrame, needle: str) -> Optional[str]:
        n = needle.lower()
        for c in df.columns:
            if n in str(c).lower():
                return c
        return None

    @staticmethod
    def _detect_header_row(df: pd.DataFrame, max_rows: int = 10) -> int:
        for i in range(min(max_rows, len(df))):
            row = df.iloc[i].astype(str).str.lower().tolist()
            has_time = any("time" in s for s in row)
            has_mass_or_cond = any(("mass" in s) or ("cond" in s) or ("pressure" in s) for s in row)
            if has_time and has_mass_or_cond:
                return i
        return 0

    def _detect_measurement_sheet(self, xl: pd.ExcelFile) -> Optional[str]:
        for sh in xl.sheet_names:
            try:
                df = xl.parse(sh, nrows=10)
            except Exception:
                continue
            cols = [str(c).lower() for c in df.columns]
            if any("time" in c for c in cols) and (any("mass" in c for c in cols) or any("cond" in c for c in cols)):
                return sh
        return xl.sheet_names[0] if xl.sheet_names else None

    def _detect_index_sheet(self, xl: pd.ExcelFile) -> Optional[str]:
        for sh in xl.sheet_names:
            try:
                df = xl.parse(sh, nrows=10)
            except Exception:
                continue
            cols = [str(c).lower() for c in df.columns]
            if (any("start" in c for c in cols) and any("end" in c for c in cols)) or \
               (any("jstart" in c for c in cols) and any("jend" in c for c in cols)):
                return sh
        return None

    def parse(self) -> DataRun:
        xl = pd.ExcelFile(self.path)
        # Resolve measurement sheet
        if self._measure_sheet_override is not None:
            if self._measure_sheet_override not in xl.sheet_names:
                raise ValueError(f"Measurement sheet '{self._measure_sheet_override}' not found.")
            self.measure_sheet = self._measure_sheet_override
        else:
            self.measure_sheet = self._detect_measurement_sheet(xl)
        if self.measure_sheet is None:
            raise ValueError("No measurement-like sheet (Time + Mass/Cond) found.")
        # Resolve index sheet
        if self._index_sheet_override is not None:
            if self._index_sheet_override not in xl.sheet_names:
                raise ValueError(f"Index sheet '{self._index_sheet_override}' not found.")
            self.index_sheet = self._index_sheet_override
        else:
            self.index_sheet = self._detect_index_sheet(xl)

        raw = xl.parse(self.measure_sheet, header=None)
        hdr = self._detect_header_row(raw, max_rows=10)
        df = raw.iloc[hdr + 1:].copy()
        df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()
        df = df.reset_index(drop=True)
        self._meas_df = df

        tcol = self._first_col_like(df, "time")
        mcol = self._first_col_like(df, "mass")
        rcol = next((c for c in df.columns if ("retent" in str(c).lower() and "cond" in str(c).lower())), None)
        pcol = next((c for c in df.columns if ("permeate" in str(c).lower() and "cond" in str(c).lower())), None)
        vcol = self._first_col_like(df, "vial") or self._first_col_like(df, "vial swap")

        def to_num(series):
            return pd.to_numeric(series, errors="coerce").to_numpy()

        time = to_num(df[tcol]) if tcol else np.arange(len(df), dtype=float)
        mass = to_num(df[mcol]) if mcol else np.full(len(df), np.nan)
        rcon = to_num(df[rcol]) if rcol else np.full(len(df), np.nan)
        pcon = to_num(df[pcol]) if pcol else np.full(len(df), np.nan)
        vial = to_num(df[vcol]) if vcol else None

        return DataRun(time=time, mass=mass, retentate_cond=rcon, permeate_cond=pcon, vial_marker=vial)

class MatSource:
    """Parses MATLAB .mat files (v7 via scipy.io.loadmat; v7.3 via h5py)."""
    def __init__(self, path: str):
        self.path = path

    def sniff(self) -> bool:
        return os.path.splitext(self.path)[1].lower() == ".mat"

    def _try_scipy(self):
        try:
            from scipy.io import loadmat
            return loadmat(self.path, squeeze_me=True, struct_as_record=False)
        except Exception:
            return None

    def _try_h5py(self):
        try:
            import h5py
            return h5py.File(self.path, "r")
        except Exception:
            return None

    def parse(self) -> DataRun:
        mat = self._try_scipy()
        if mat is None:
            h5 = self._try_h5py()
            if h5 is None:
                raise ValueError("Could not read .mat with scipy or h5py.")
            def _h5(name: str):
                return np.array(h5[name]).squeeze() if name in h5 else None
            time = _h5("Time") or _h5("time")
            mass = _h5("Mass") or _h5("mass")
            rcon = _h5("RetentateCond") or _h5("Retentate_Cond") or _h5("retentate_cond")
            pcon = _h5("PermeateCond") or _h5("Permeate_Cond") or _h5("permeate_cond")
            vial = _h5("Vial") or _h5("VialSwap") or _h5("vial") or _h5("vial_swap")
            h5.close()
        else:
            def _find_key(d, *cands):
                keys = {k.lower(): k for k in d.keys()}
                for c in cands:
                    if c.lower() in keys:
                        return d[keys[c.lower()]]
                return None
            top = _find_key(mat, "data_stru", "data", "data_struct", "S", "D")
            src = top.__dict__ if hasattr(top, "__dict__") else mat
            def _grab(*names):
                for nm in names:
                    v = _find_key(src, nm)
                    if v is not None:
                        return np.array(v).squeeze()
                return None
            time = _grab("Time", "time", "t")
            mass = _grab("Mass", "mass", "m")
            rcon = _grab("RetentateCond", "Retentate_Cond", "retentate_cond", "ret_cond")
            pcon = _grab("PermeateCond", "Permeate_Cond", "permeate_cond", "perm_cond")
            vial = _grab("Vial", "VialSwap", "vial", "vial_swap")

        if time is None or mass is None:
            raise ValueError("MAT file missing required arrays: Time and Mass.")
        return DataRun(time=np.asarray(time), mass=np.asarray(mass),
                       retentate_cond=np.asarray(rcon) if rcon is not None else None,
                       permeate_cond=np.asarray(pcon) if pcon is not None else None,
                       vial_marker=np.asarray(vial) if vial is not None else None)

# ------------------------------- Windowing strategies -----------------------------------------

class WindowStrategy(Protocol):
    def infer(self, source: ExcelSource, meas_df: Optional[pd.DataFrame], run: DataRun) -> List[Tuple[int, int]]: ...

class IndexSheetStrategy:
    def infer(self, source: ExcelSource, meas_df: Optional[pd.DataFrame], run: DataRun) -> List[Tuple[int, int]]:
        windows: List[Tuple[int, int]] = []
        if not isinstance(source, ExcelSource) or source.index_sheet is None:
            return windows
        xl = pd.ExcelFile(source.path)
        try:
            idx_df = xl.parse(source.index_sheet)
        except Exception:
            return windows
        def first_col_like(df: pd.DataFrame, needle: str):
            n = needle.lower()
            for c in df.columns:
                if n in str(c).lower():
                    return c
            return None
        c_start = first_col_like(idx_df, "start") or first_col_like(idx_df, "jstart")
        c_end   = first_col_like(idx_df, "end")   or first_col_like(idx_df, "jend")
        if c_start is None or c_end is None:
            if idx_df.shape[1] >= 3:
                c_start = idx_df.columns[1]
                c_end = idx_df.columns[2]
        if c_start is None or c_end is None:
            return windows
        starts = pd.to_numeric(idx_df[c_start], errors="coerce").dropna().astype(int).to_numpy()
        ends   = pd.to_numeric(idx_df[c_end], errors="coerce").dropna().astype(int).to_numpy()
        n = min(len(starts), len(ends))
        windows = [(int(starts[i]), int(ends[i])) for i in range(n) if ends[i] >= starts[i]]
        return windows

class VialMarkerStrategy:
    def infer(self, source: ExcelSource, meas_df: Optional[pd.DataFrame], run: DataRun) -> List[Tuple[int, int]]:
        v = run.vial_marker
        if v is None or np.asarray(v).size == 0:
            return []
        v_arr = np.asarray(v).squeeze()
        try:
            vnum = pd.to_numeric(pd.Series(v_arr), errors="coerce").fillna(-1).to_numpy()
        except Exception:
            v_str = np.array(list(map(str, v_arr)))
            vnum = (v_str != np.roll(v_str, 1)).astype(int).cumsum()
        boundaries = [0] + [i for i in range(1, len(vnum)) if vnum[i] != vnum[i-1]] + [len(vnum)]
        windows: List[Tuple[int, int]] = []
        for s, e in zip(boundaries[:-1], boundaries[1:]):
            if e - s >= 3:
                windows.append((int(s), int(e - 1)))
        return windows

class SingleWindowStrategy:
    def infer(self, source: ExcelSource, meas_df: Optional[pd.DataFrame], run: DataRun) -> List[Tuple[int, int]]:
        n = run.time.size
        return [(0, int(n - 1))] if n > 0 else []

# ------------------------------- Orchestrator + plotting --------------------------------------

class DataLoader:
    """Selects a DataSource, infers windows, exposes per-window arrays, and plots quicklooks."""
    def __init__(self, path: str, *, measure_sheet: Optional[str] = None, index_sheet: Optional[str] = None):
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")
        self.path = path
        self._measure_sheet_override = measure_sheet
        self._index_sheet_override = index_sheet
        # outputs
        self.times: List[np.ndarray] = []
        self.masses: List[np.ndarray] = []
        self.retentate_cond: List[np.ndarray] = []
        self.permeate_cond: List[np.ndarray] = []
        self.windows: List[Tuple[int, int]] = []
        # internals
        self._source: Optional[DataSource] = None
        self._run: Optional[DataRun] = None
        self._meas_df: Optional[pd.DataFrame] = None
        # strategies
        self._strategies: List[WindowStrategy] = [IndexSheetStrategy(), VialMarkerStrategy(), SingleWindowStrategy()]

    def _select_source(self) -> None:
        excel = ExcelSource(self.path, self._measure_sheet_override, self._index_sheet_override)
        mat = MatSource(self.path)
        if excel.sniff():
            self._source = excel
        elif mat.sniff():
            self._source = mat
        else:
            raise ValueError(f"Unsupported file type: {os.path.splitext(self.path)[1].lower()} (expected Excel or .mat)")

    def load(self) -> None:
        self._select_source()
        run = self._source.parse()  # type: ignore[union-attr]
        self._run = run
        if isinstance(self._source, ExcelSource):
            xl = pd.ExcelFile(self.path)
            measure_sheet = self._source.measure_sheet
            raw = xl.parse(measure_sheet, header=None)
            hdr = ExcelSource._detect_header_row(raw)
            df = raw.iloc[hdr + 1:].copy()
            df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()
            df = df.reset_index(drop=True)
            self._meas_df = df
        windows: List[Tuple[int, int]] = []
        for strat in self._strategies:
            windows = strat.infer(self._source, self._meas_df, run)  # type: ignore[arg-type]
            if windows:
                break
        if not windows and run.time.size > 0:
            windows = [(0, int(run.time.size - 1))]
        self.windows = windows
        # slice per-window arrays and normalize mass + spike filter
        self.times, self.masses, self.retentate_cond, self.permeate_cond = [], [], [], []
        for (s, e) in self.windows:
            sl = slice(s, e + 1)
            t = run.time[sl]
            mraw = run.mass[sl]
            m = mraw - (mraw[0] if mraw.size and not np.isnan(mraw[0]) else 0.0)
            m = np.where(m > 1.2, np.nan, m)
            r = run.retentate_cond[sl]
            p = run.permeate_cond[sl]
            self.times.append(t)
            self.masses.append(m)
            self.retentate_cond.append(r)
            self.permeate_cond.append(p)

    # ---------- plotting (with vial labels and unified colors) ----------
    def _vial_boundaries_and_centers(self):
        boundary_times: List[float] = []
        centers: List[float] = []
        labels: List[str] = []
        for i, t in enumerate(self.times):
            if len(t):
                centers.append(float(np.nanmean(t)))
                labels.append(f"Vial {i}")
        for t in self.times[:-1]:
            if len(t):
                boundary_times.append(float(t[-1]))
        return boundary_times, centers, labels

    def _decorate_with_vials(self, ax, boundary_times, centers, labels):
        for bt in boundary_times:
            ax.axvline(bt, color="gray", linestyle=":", linewidth=1)
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        if centers:
            ax2.set_xticks(centers)
            ax2.set_xticklabels(labels)
        ax2.set_xlabel("Vials")

    def plot_quicklook(self, save_prefix: Optional[str] = None) -> None:
        boundary_times, centers, labels = self._vial_boundaries_and_centers()

        # Mass vs Time (additive across vials)
        if any(len(m) > 0 for m in self.masses):
            fig, ax = plt.subplots()
            color_mass = "blue"  # single color, solid

            cumulative_offset = 0.0
            first = True
            for t, m in zip(self.times, self.masses):
                if len(t) == 0 or len(m) == 0:
                    continue
                mask = ~np.isnan(m)
                if mask.any():
                    y = m[mask] + cumulative_offset
                    ax.plot(
                        t[mask], y,
                        color=color_mass, linewidth=1.5,
                        label="Cumulative mass change" if first else None,
                    )
                    # advance offset by the last finite delta in this vial
                    cumulative_offset += float(m[mask][-1])
                    first = False

            ax.set_xlabel("Time [s]")
            ax.set_ylabel("Mass change [g]")
            ax.set_title("Mass vs Time (additive)")
            self._decorate_with_vials(ax, *self._vial_boundaries_and_centers())
            if not first:
                ax.legend(loc="best")
            fig.tight_layout()
            if save_prefix:
                fig.savefig(f"{save_prefix}_mass.png", dpi=160)

        # Conductivity vs Time
        has_ret = any(np.isfinite(arr).any() for arr in self.retentate_cond)
        has_perm = any(np.isfinite(arr).any() for arr in self.permeate_cond)
        if has_ret or has_perm:
            fig, ax = plt.subplots()
            color_ret = "black"
            color_perm = "red"
            first_ret = True
            for t, r in zip(self.times, self.retentate_cond):
                if r.size and len(r) == len(t):
                    mask = np.isfinite(r)
                    if mask.any():
                        ax.plot(t[mask], r[mask], color=color_ret, linewidth=1.5,
                                label="Retentate (solid)" if first_ret else None)
                        first_ret = False
            first_perm = True
            for t, p in zip(self.times, self.permeate_cond):
                if p.size and len(p) == len(t):
                    mask = np.isfinite(p)
                    if mask.any():
                        ax.plot(t[mask], p[mask], linestyle="--", color=color_perm, linewidth=1.5,
                                label="Permeate (dashed)" if first_perm else None)
                        first_perm = False
            ax.set_xlabel("Time [s]")
            ax.set_ylabel("Conductivity [arb]")
            ax.set_title("Conductivity vs Time")
            self._decorate_with_vials(ax, boundary_times, centers, labels)
            if (not first_ret) or (not first_perm):
                ax.legend(loc="best")
            fig.tight_layout()
            if save_prefix:
                fig.savefig(f"{save_prefix}_cond.png", dpi=160)

        plt.show()

    # Optional getters
    @property
    def measure_sheet_name(self) -> Optional[str]:
        return getattr(self._source, "measure_sheet", None) if isinstance(self._source, ExcelSource) else None

    @property
    def index_sheet_name(self) -> Optional[str]:
        return getattr(self._source, "index_sheet", None) if isinstance(self._source, ExcelSource) else None

    @property
    def measurement_dataframe(self) -> Optional[pd.DataFrame]:
        return self._meas_df

# --------------------------------- CLI & helpers -----------------------------------------------

def list_sheets(path: str) -> List[str]:
    ext = os.path.splitext(path)[1].lower()
    if ext not in {".xlsx", ".xls", ".xlsm", ".xlsb"}:
        return []
    xl = pd.ExcelFile(path)
    return list(xl.sheet_names)


def _cli_list_sheets(path: str) -> None:
    ext = os.path.splitext(path)[1].lower()
    if ext not in {".xlsx", ".xls", ".xlsm", ".xlsb"}:
        print("[info] Not an Excel file; no sheets to list.")
        return
    xl = pd.ExcelFile(path)
    print("Sheets:")
    for sh in xl.sheet_names:
        print(f"  - {sh}")


def _cli_export_csv(loader: 'DataLoader', outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    for i, (t, m, r, p) in enumerate(zip(loader.times, loader.masses, loader.retentate_cond, loader.permeate_cond)):
        df = pd.DataFrame({
            'time_s': t,
            'dmass_g': m,
            'retentate_cond': r,
            'permeate_cond': p,
        })
        df.to_csv(os.path.join(outdir, f"vial_{i:02d}.csv"), index=False)
    print(f"[ok] Exported {len(loader.times)} CSV files to {outdir}")


def _main():
    import argparse
    parser = argparse.ArgumentParser(description="Diafiltration data loader (Excel/MAT) with quicklook plots.")
    parser.add_argument('path', help='Path to .xlsx/.xls/.xlsm/.xlsb or .mat file')
    parser.add_argument('--quicklook', action='store_true', help='Generate quicklook plots')
    parser.add_argument('--save-prefix', default=None, help='If set, save plots with this prefix')
    parser.add_argument('--list-sheets', action='store_true', help='List Excel sheet names and exit')
    parser.add_argument('--dump-windows', action='store_true', help='Print inferred (start,end) windows')
    parser.add_argument('--export-csv', default=None, metavar='OUTDIR', help='Export each vial window to CSV files')
    parser.add_argument('--measure-sheet', default=None, help='Override measurement sheet name')
    parser.add_argument('--index-sheet', default=None, help='Override index sheet name')

    args = parser.parse_args()
    if args.list_sheets:
        _cli_list_sheets(args.path)
        return

    loader = DataLoader(args.path, measure_sheet=args.measure_sheet, index_sheet=args.index_sheet)
    loader.load()

    if args.dump_windows:
        print("Windows (start_idx, end_idx):")
        for w in loader.windows:
            print(w)

    if args.quicklook or args.save_prefix:
        loader.plot_quicklook(save_prefix=args.save_prefix)

    if args.export_csv:
        _cli_export_csv(loader, args.export_csv)


if __name__ == '__main__':
    _main()

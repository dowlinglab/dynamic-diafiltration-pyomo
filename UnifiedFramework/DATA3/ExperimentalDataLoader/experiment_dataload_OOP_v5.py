#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 08:45:29 2025

@author: kkasturi

Diafiltration Data Loader
Supports Excel (.xlsx, .xls, .xlsm) and MATLAB (.mat) input.
This version adds PyDataStru (MATLAB-like data_stru in Python) but does not
modify the existing DataLoader workflow.
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
import scipy.io as sio
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Any, Protocol, runtime_checkable
import matplotlib.pyplot as plt

# =============================================================================
# SECTION 1 — PyDataStru (NEW)
# Python analog of MATLAB data_stru (data_config + data_raw)
# =============================================================================

@dataclass
class PyDataConfig:
    """Python mirror of MATLAB data_stru.data_config."""
    M_F0: float | None = None
    C_F0: np.ndarray | None = None
    M_O: float | None = None
    C_D: np.ndarray | None = None
    nc: int | None = None
    namec: Any | None = None
    ni: int | None = None
    delP: float | None = None
    Temp: float | None = None
    Am: float | None = None
    rho: float | None = None
    Lp0: float | None = None
    B0: np.ndarray | None = None
    sigma0: float | None = None
    theta0: np.ndarray | None = None
    n: int | None = None
    nr: int | None = None
    n_extra: int | None = None
    extras: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PyVial:
    """Python mirror of MATLAB data_stru.data_raw(i)."""
    number: int
    time: np.ndarray
    mass: np.ndarray
    cV_avg: np.ndarray | None = None
    cF_exp: np.ndarray | None = None
    nr: int | None = None
    extras: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PyDataStru:
    """
    Python representation of a full diafiltration experiment, mirroring MATLAB data_stru.
    """
    dataset: float | None
    filename: str | None
    mode: str | None
    conductivity_cF: bool | None
    n_obj: int | None
    data_config: PyDataConfig
    data_raw: List[PyVial]

    # -------------------------------------------------------------------------
    # Constructor from MATLAB .mat (load_data.m output)
    # -------------------------------------------------------------------------
    @staticmethod
    def from_mat(path: str):
        mat = sio.loadmat(path, squeeze_me=True, struct_as_record=False)

        if "data_stru" not in mat:
            raise ValueError(f"{path} does not contain 'data_stru'")

        ds = mat["data_stru"]

        # Top-level
        dataset = float(getattr(ds, "dataset", np.nan))
        filename = str(getattr(ds, "filename", None))
        mode = str(getattr(ds, "mode", None))
        conductivity_cF = bool(getattr(ds, "conductivity_cF", False))
        n_obj = int(getattr(ds, "n_obj", 0))

        # data_config
        dc = getattr(ds, "data_config", None)
        cfg_kwargs = {}
        extras = {}

        if dc is not None and hasattr(dc, "_fieldnames"):
            for fn in dc._fieldnames:
                val = getattr(dc, fn)
                # convert arrays from MATLAB
                if isinstance(val, (list, np.ndarray)):
                    val = np.asarray(val).squeeze()
                if fn in PyDataConfig.__dataclass_fields__:
                    cfg_kwargs[fn] = val
                else:
                    extras[fn] = val

        data_config = PyDataConfig(**cfg_kwargs, extras=extras)

        # data_raw (vials)
        dr = getattr(ds, "data_raw", [])
        if isinstance(dr, np.ndarray):
            dr = dr.ravel()

        vials = []
        for i, v in enumerate(dr):
            t = np.asarray(getattr(v, "time", []), dtype=float).squeeze()
            m = np.asarray(getattr(v, "mass", []), dtype=float).squeeze()

            cV = getattr(v, "cV_avg", None)
            if cV is not None:
                cV = np.asarray(cV).squeeze()

            cF = getattr(v, "cF_exp", None)
            if cF is not None:
                cF = np.asarray(cF).squeeze()

            nr = getattr(v, "nr", None)
            num = getattr(v, "number", i + 1)

            vial_extras = {}
            if hasattr(v, "_fieldnames"):
                for fn in v._fieldnames:
                    if fn not in ["time", "mass", "cV_avg", "cF_exp", "nr", "number"]:
                        vial_extras[fn] = getattr(v, fn)

            vials.append(PyVial(num, t, m, cV, cF, nr, vial_extras))

        return PyDataStru(
            dataset=dataset,
            filename=filename,
            mode=mode,
            conductivity_cF=conductivity_cF,
            n_obj=n_obj,
            data_config=data_config,
            data_raw=vials,
        )

    # -------------------------------------------------------------------------
    # Constructor from Excel Sheet
    # NOTE: This version does NOT handle vial segmentation yet.
    # It creates ONE PyVial per sheet.
    # -------------------------------------------------------------------------
    @staticmethod
    def from_excel(path: str, sheet: str):

        df = pd.read_excel(path, sheet_name=sheet)

        # basic required columns
        # tcol = [c for c in df.columns if "time" in str(c).lower()]
        # mcol = [c for c in df.columns if "mass" in str(c).lower()]

        # if not tcol or not mcol:
        #     raise ValueError("Excel sheet must have columns for time and mass.")
        
        # --- Robust time column detection ---
        time_candidates = [
            "time", "t", "record", "sec", "seconds"
        ]
        
        def has_any(name: str, keywords: list[str]):
            s = name.lower()
            return any(k in s for k in keywords)
        
        tcols = [c for c in df.columns if has_any(str(c), time_candidates)]
        if not tcols:
            raise ValueError(f"No time column found in Excel sheet. Columns = {list(df.columns)}")
        
        tcol = tcols[0]
        
        # --- Robust mass column detection ---
        mass_candidates = [
            "mass", "m(g)", "weight", "collection", "permeate mass", "permeate"
        ]
        
        mcols = [c for c in df.columns if has_any(str(c), mass_candidates)]
        if not mcols:
            raise ValueError(f"No mass column found in Excel sheet. Columns = {list(df.columns)}")
        
        mcol = mcols[0]


        time = pd.to_numeric(df[tcol[0]], errors="coerce").to_numpy()
        mass = pd.to_numeric(df[mcol[0]], errors="coerce").to_numpy()

        # optional conductivity columns
        cF = df.filter(regex="retent.*cond", case=False).apply(pd.to_numeric, errors="coerce")
        cV = df.filter(regex="permeat.*cond", case=False).apply(pd.to_numeric, errors="coerce")

        vial = PyVial(
            number=1,
            time=time,
            mass=mass,
            cF_exp=cF.to_numpy().squeeze() if not cF.empty else None,
            cV_avg=cV.to_numpy().squeeze() if not cV.empty else None,
            nr=np.isfinite(mass).sum(),
            extras={}
        )

        cfg = PyDataConfig()

        return PyDataStru(
            dataset=None,
            filename=os.path.basename(path),
            mode=None,
            conductivity_cF=False,
            n_obj=1,
            data_config=cfg,
            data_raw=[vial]
        )

# =============================================================================
# SECTION 2 — EXISTING DataRun, DataLoader, DataSource, plotting, CLI
# (UNCHANGED — your code continues here exactly as before)
# =============================================================================

# ------------- KEEP YOUR ORIGINAL CODE BELOW THIS LINE ---------------

# -------------------------------- Data container -----------------------------------------------

@dataclass
class DataRun:
    """
    Canonical 1D run representation, independent of source (Excel vs .mat).

    - time, mass, retentate_cond, permeate_cond, vial_marker: core series.
    - extra_signals: all other numeric columns (pressure, temps, ICP, etc.).
    - metadata: run-level information (sheet name, header rows, data_stru scalars).
    """
    time: np.ndarray
    mass: np.ndarray
    retentate_cond: np.ndarray
    permeate_cond: np.ndarray
    vial_marker: Optional[np.ndarray] = None

    extra_signals: Dict[str, np.ndarray] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        n = int(np.asarray(self.time).size)

        def _fit(x, fill=np.nan):
            if x is None:
                return np.full(n, fill)
            arr = np.asarray(x).squeeze()
            arr = arr[:n] if arr.size >= n else np.pad(arr, (0, n - arr.size), constant_values=fill)
            return arr

        # Core series
        self.time = _fit(self.time).astype(float)
        self.mass = _fit(self.mass).astype(float)
        self.retentate_cond = _fit(self.retentate_cond).astype(float)
        self.permeate_cond = _fit(self.permeate_cond).astype(float)

        # Vial marker (if present)
        if self.vial_marker is not None:
            vm = np.asarray(self.vial_marker).squeeze()
            self.vial_marker = vm[:n] if vm.size >= n else np.pad(vm, (0, n - vm.size), constant_values=np.nan)

        # Normalize extra signals: all 1D arrays of length n
        norm_extra: Dict[str, np.ndarray] = {}
        for name, arr in self.extra_signals.items():
            a = np.asarray(arr).squeeze()
            a = a[:n] if a.size >= n else np.pad(a, (0, n - a.size), constant_values=np.nan)
            norm_extra[name] = a.astype(float)
        self.extra_signals = norm_extra

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
        """
        Try to find the row containing column labels by looking for 'time' and either 'mass' or 'cond' or 'pressure'.
        """
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
        """
        Try to detect a sheet with start/end indices for vials (like exp_index in MATLAB).
        """
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
        """
        Parse an Excel sheet into a DataRun, extracting:

        - time, mass, retentate_cond, permeate_cond, vial_marker (canonical)
        - extra_signals: all other numeric columns
        - metadata: sheet name + header rows (for notes, feed, datapoints, etc.)
        """
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

        # Read full sheet with no header; then detect header row
        raw = xl.parse(self.measure_sheet, header=None)
        hdr = self._detect_header_row(raw, max_rows=10)

        # Everything above hdr is "header region" and treated as metadata
        header_region = raw.iloc[:hdr].copy()

        # Measurement table: header row + data rows below
        df = raw.iloc[hdr + 1:].copy()
        df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()
        df = df.reset_index(drop=True)
        self._meas_df = df

        # Convert all columns to numeric arrays (NaN where non-numeric)
        def to_num(series: pd.Series) -> np.ndarray:
            return pd.to_numeric(series, errors="coerce").to_numpy()

        numeric_cols: Dict[str, np.ndarray] = {
            str(col): to_num(df[col]) for col in df.columns
        }

        # Identify canonical columns
        tcol = self._first_col_like(df, "time")
        mcol = self._first_col_like(df, "mass")
        rcol = next(
            (c for c in df.columns if ("retent" in str(c).lower() and "cond" in str(c).lower())),
            None,
        )
        pcol = next(
            (c for c in df.columns if ("permeate" in str(c).lower() and "cond" in str(c).lower())),
            None,
        )
        vcol = self._first_col_like(df, "vial") or self._first_col_like(df, "vial swap")

        n_rows = len(df)

        if tcol:
            time = numeric_cols.pop(str(tcol))
        else:
            time = np.arange(n_rows, dtype=float)

        mass = numeric_cols.pop(str(mcol)) if mcol else np.full(n_rows, np.nan, dtype=float)
        rcon = numeric_cols.pop(str(rcol)) if rcol else np.full(n_rows, np.nan, dtype=float)
        pcon = numeric_cols.pop(str(pcol)) if pcol else np.full(n_rows, np.nan, dtype=float)
        vial = numeric_cols.pop(str(vcol)) if vcol else None

        # Everything else stays as extra signals (pressure, temps, ICP, etc.)
        extra_signals = dict(numeric_cols)

        # Metadata: sheet name + header rows (for feed, diafiltrate, notes, datapoints, etc.)
        metadata: Dict[str, Any] = {
            "sheet_name": self.measure_sheet,
            "header_rows": header_region,
        }

        return DataRun(
            time=time,
            mass=mass,
            retentate_cond=rcon,
            permeate_cond=pcon,
            vial_marker=vial,
            extra_signals=extra_signals,
            metadata=metadata,
        )

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

    @staticmethod
    def _data_stru_to_datarun(ds) -> DataRun:
        """
        Convert a MATLAB data_stru struct (output of load_data.m) into a canonical DataRun.

        - Concatenates all vials in data_stru.data_raw into one run.
        - Uses:
            time   <- data_raw(i).time
            mass   <- data_raw(i).mass
            rcon   <- data_raw(i).cF_exp  (retentate, possibly time-series)
            pcon   <- data_raw(i).cV_avg  (permeate, per-vial or time-series)
            vial_marker <- data_raw(i).number
        - Stores top-level scalars (dataset, mode, filename, etc.) in metadata.
        """
        def _to_arr(val):
            if val is None:
                return None
            return np.asarray(val).squeeze()

        raw = getattr(ds, "data_raw", None)
        if raw is None:
            raise ValueError("data_stru has no data_raw field")

        if isinstance(raw, np.ndarray):
            vial_structs = list(raw.ravel())
        else:
            vial_structs = [raw]

        times_list: List[np.ndarray] = []
        mass_list: List[np.ndarray] = []
        ret_list: List[np.ndarray] = []
        perm_list: List[np.ndarray] = []
        vial_markers: List[np.ndarray] = []

        for idx, vial in enumerate(vial_structs):
            t = _to_arr(getattr(vial, "time", None))
            if t is None or t.size == 0:
                continue
            t = t.astype(float)
            n = int(t.size)

            m = _to_arr(getattr(vial, "mass", None))
            if m is None:
                m = np.full(n, np.nan, dtype=float)
            else:
                m = m.astype(float)
                m = m[:n] if m.size >= n else np.pad(m, (0, n - m.size), constant_values=np.nan)

            cF = _to_arr(getattr(vial, "cF_exp", None))
            cV = _to_arr(getattr(vial, "cV_avg", None))

            def _fit_len(arr):
                if arr is None:
                    return np.full(n, np.nan, dtype=float)
                a = np.asarray(arr).squeeze()
                if a.ndim > 1:
                    a = a[:, 0]
                a = a[:n] if a.size >= n else np.pad(a, (0, n - a.size), constant_values=np.nan)
                return a.astype(float)

            r = _fit_len(cF)
            p = _fit_len(cV)

            vnum = getattr(vial, "number", idx + 1)
            vm = np.full(n, float(vnum), dtype=float)

            times_list.append(t)
            mass_list.append(m)
            ret_list.append(r)
            perm_list.append(p)
            vial_markers.append(vm)

        if not times_list:
            raise ValueError("No non-empty vials found in data_stru.data_raw")

        time = np.concatenate(times_list)
        mass = np.concatenate(mass_list)
        rcon = np.concatenate(ret_list)
        pcon = np.concatenate(perm_list)
        vial_marker = np.concatenate(vial_markers)

        metadata: Dict[str, Any] = {}
        for attr in ("dataset", "filename", "mode", "conductivity_cF", "n_obj"):
            if hasattr(ds, attr):
                metadata[attr] = getattr(ds, attr)

        # Optionally stash data_config summary if present
        if hasattr(ds, "data_config"):
            dc = getattr(ds, "data_config")
            if hasattr(dc, "_fieldnames"):
                dc_meta = {}
                for name in dc._fieldnames:
                    val = getattr(dc, name)
                    # light conversion of numpy scalars
                    if isinstance(val, np.generic):
                        val = val.item()
                    metadata[f"data_config.{name}"] = val

        return DataRun(
            time=time,
            mass=mass,
            retentate_cond=rcon,
            permeate_cond=pcon,
            vial_marker=vial_marker,
            extra_signals={},
            metadata=metadata,
        )

    def parse(self) -> DataRun:
        mat = self._try_scipy()
        if mat is None:
            # h5py fallback for v7.3 MAT files with simple arrays
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

            if time is None or mass is None:
                raise ValueError("MAT file missing required arrays: Time and Mass.")

            return DataRun(
                time=np.asarray(time),
                mass=np.asarray(mass),
                retentate_cond=np.asarray(rcon) if rcon is not None else None,
                permeate_cond=np.asarray(pcon) if pcon is not None else None,
                vial_marker=np.asarray(vial) if vial is not None else None,
            )

        # scipy path
        keys_no_meta = {k for k in mat.keys() if not k.startswith("__")}

        # Prefer explicit data_stru struct if present (output of load_data.m)
        if "data_stru" in keys_no_meta:
            ds = mat["data_stru"]
            return self._data_stru_to_datarun(ds)

        # Generic MAT arrays
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

        return DataRun(
            time=np.asarray(time),
            mass=np.asarray(mass),
            retentate_cond=np.asarray(rcon) if rcon is not None else None,
            permeate_cond=np.asarray(pcon) if pcon is not None else None,
            vial_marker=np.asarray(vial) if vial is not None else None,
        )

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
        c_end = first_col_like(idx_df, "end") or first_col_like(idx_df, "jend")

        if c_start is None or c_end is None:
            if idx_df.shape[1] >= 3:
                c_start = idx_df.columns[1]
                c_end = idx_df.columns[2]
        if c_start is None or c_end is None:
            return windows

        starts = pd.to_numeric(idx_df[c_start], errors="coerce").dropna().astype(int).to_numpy()
        ends = pd.to_numeric(idx_df[c_end], errors="coerce").dropna().astype(int).to_numpy()

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
        boundaries = [0] + [i for i in range(1, len(vnum)) if vnum[i] != vnum[i - 1]] + [len(vnum)]
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

        # For Excel, keep a copy of the measurement table for strategies / inspection
        if isinstance(self._source, ExcelSource):
            xl = pd.ExcelFile(self.path)
            measure_sheet = self._source.measure_sheet
            raw = xl.parse(measure_sheet, header=None)
            hdr = ExcelSource._detect_header_row(raw)
            df = raw.iloc[hdr + 1:].copy()
            df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()
            df = df.reset_index(drop=True)
            self._meas_df = df

        # Infer windows
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
            # zero mass per window and spike-filter big jumps
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
                labels.append(f"Vial {i+1}")
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

    def plot_mass(self, save_prefix: Optional[str] = None, *, cumulative: bool = False) -> None:
        """Plot Mass vs Time.
        Args:
            save_prefix: if provided, saves to f"{save_prefix}_mass.png" (per-vial) or
                         f"{save_prefix}_mass_cum.png" (cumulative)
            cumulative: if True, plot a single additive line across vials; otherwise
                        plot per-vial Δmass with a single color and legend.
        """
        boundary_times, centers, labels = self._vial_boundaries_and_centers()
        if not any(len(m) > 0 for m in self.masses):
            return
        fig, ax = plt.subplots()
        color_mass = "blue"
        if cumulative:
            cumulative_offset = 0.0
            first = True
            for t, m in zip(self.times, self.masses):
                if len(t) == 0 or len(m) == 0:
                    continue
                mask = ~np.isnan(m)
                if mask.any():
                    y = m[mask] + cumulative_offset
                    ax.plot(t[mask], y, color=color_mass, linewidth=1.5,
                            label="Cumulative mass change" if first else None)
                    cumulative_offset += float(m[mask][-1])
                    first = False
            ax.set_title("Cumulative Permeate Mass vs Time")
            ax.set_ylabel("Cumulative permeate mass [g]")
        else:
            first = True
            for t, m in zip(self.times, self.masses):
                if len(t) == 0 or len(m) == 0:
                    continue
                mask = ~np.isnan(m)
                if mask.any():
                    ax.plot(t[mask], m[mask], color=color_mass, linewidth=1.5,
                            label="Mass change" if first else None)
                    first = False
            ax.set_title("Mass vs Time")
            ax.set_ylabel("Mass change [g]")
        ax.set_xlabel("Time [s]")
        self._decorate_with_vials(ax, boundary_times, centers, labels)
        if not first:
            ax.legend(loc="best")
        fig.tight_layout()
        if save_prefix:
            suffix = "mass_cum" if cumulative else "mass"
            fig.savefig(f"{save_prefix}_{suffix}.png", dpi=160)

    def plot_conductivity(self, save_prefix: Optional[str] = None) -> None:
        """Plot Conductivity vs Time with unified colors (retentate solid, permeate dashed)."""
        boundary_times, centers, labels = self._vial_boundaries_and_centers()
        has_ret = any(np.isfinite(arr).any() for arr in self.retentate_cond)
        has_perm = any(np.isfinite(arr).any() for arr in self.permeate_cond)
        if not (has_ret or has_perm):
            return
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

    def plot_quicklook(self, save_prefix: Optional[str] = None, *, cumulative: bool = False) -> None:
        """Convenience wrapper: plot mass (configurable cumulative) and conductivity."""
        self.plot_mass(save_prefix=save_prefix, cumulative=cumulative)
        self.plot_conductivity(save_prefix=save_prefix)
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

def _cli_export_csv(loader: "DataLoader", outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    for i, (t, m, r, p) in enumerate(zip(loader.times, loader.masses, loader.retentate_cond, loader.permeate_cond)):
        df = pd.DataFrame({
            "time_s": t,
            "dmass_g": m,
            "retentate_cond": r,
            "permeate_cond": p,
        })
        df.to_csv(os.path.join(outdir, f"vial_{i:02d}.csv"), index=False)
    print(f"[ok] Exported {len(loader.times)} CSV files to {outdir}")

def _main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Diafiltration data loader (Excel/MAT) with quicklook plots and flexible mass/conductivity options."
    )
    parser.add_argument("path", help="Path to .xlsx/.xls/.xlsm/.xlsb or .mat file")
    parser.add_argument("--quicklook", action="store_true", help="Generate plots (mass and conductivity unless filtered)")
    parser.add_argument("--save-prefix", default=None, help="If set, save plots with this prefix")
    parser.add_argument("--list-sheets", action="store_true", help="List Excel sheet names and exit")
    parser.add_argument("--dump-windows", action="store_true", help="Print inferred (start,end) windows")
    parser.add_argument("--export-csv", default=None, metavar="OUTDIR", help="Export each vial window to CSV files")
    parser.add_argument("--measure-sheet", default=None, help="Override measurement sheet name")
    parser.add_argument("--index-sheet", default=None, help="Override index sheet name")
    parser.add_argument("--cumulative", action="store_true", help="Mass plot: additive across vials (cumulative mode)")
    parser.add_argument("--mass-only", action="store_true", help="Only plot mass vs time")
    parser.add_argument("--cond-only", action="store_true", help="Only plot conductivity vs time")

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
        if args.mass_only and not args.cond_only:
            loader.plot_mass(save_prefix=args.save_prefix, cumulative=args.cumulative)
        elif args.cond_only and not args.mass_only:
            loader.plot_conductivity(save_prefix=args.save_prefix)
        else:
            loader.plot_quicklook(save_prefix=args.save_prefix, cumulative=args.cumulative)

    if args.export_csv:
        _cli_export_csv(loader, args.export_csv)

if __name__ == "__main__":
    _main()

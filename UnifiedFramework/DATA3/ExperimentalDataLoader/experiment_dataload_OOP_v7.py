#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unified experimental data loader for diafiltration.

authors: Keshav Kasturi Rangan, Alex Dowling

Supports:
    - MATLAB .mat files containing data_stru (output of load_data.m)
    - Excel .xlsx files where each sheet is an experimental run.

Provides:
    - PyDataStru / PyVial / PyDataConfig:
        PyDataStru.from_mat(path)
        PyDataStru.from_excel(path, sheet)

    - DataLoader:
        DataLoader(path, measure_sheet=None)
        .load()
        .times, .masses, .retentate_cond, .permeate_cond, .windows

    - list_sheets(path): convenience for Excel files
"""

import os # Standard library: path handling, file existence checks, directory joins, etc.
          # Used throughout the loader to detect files, handle Excel/MAT paths, etc.import re

import re # Python's regular expression engine.
          # We use it for case-insensitive pattern matching of Excel column names
          # (e.g., re.compile("retent.*cond", flags=re.I)).

from dataclasses import dataclass, field # 'dataclass' turns a Python class into a light, structured container
                                         # (MATLAB struct equivalent). We use it for PyVial, PyDataStru, DataRun, etc.
                                         # 'field' is used to declare mutable default values (e.g., dicts, lists).

from typing import Any, Dict, List, Optional, Protocol, Tuple, runtime_checkable # Type hints to clarify interfaces and structure.
                                                                                 # - Optional[X] = X or None
                                                                                 # - Protocol = an interface (like an abstract class). We use it for DataSource.
                                                                                 # - runtime_checkable allows checking if an object implements a Protocol.
                                                                                 # - Dict, List, Tuple are just typed containers.

import numpy as np # NumPy: the numerical backbone of the loader.
                   # Used for arrays (time, mass, conductivity), padding, concatenation, NaN handling, etc.

import pandas as pd # Pandas: Excel parsing, column detection, row slicing, measurement-table reconstruction.
                    # We rely on pandas to read Excel data and manipulate tabular structures.

import scipy.io as sio # SciPy I/O: needed for reading MATLAB .mat files.
                       # 'sio.loadmat' reads MATLAB structs (data_stru) into nested Python objects.

import matplotlib.pyplot as plt # Matplotlib: optional plotting for quick visual inspection of time/mass/cond data.
                                # Used in DataLoader.plot_mass() and plot_conductivity().

# =============================================================================
# Shared helper: read Excel measurement table
# =============================================================================

def _excel_get_meas_df(path: str, sheet: str, max_header_scan: int = 25) -> Tuple[pd.DataFrame, int]:
    """
    Read an Excel sheet and return:
        df   = measurement table (data rows with correct column names)
        hdr  = index of the header row in the raw sheet

    We scan the first `max_header_scan` rows looking for a row that:
        - contains 'time'
        - contains either 'mass' or 'cond' (conductivity)
    """
    raw = pd.read_excel(path, sheet_name=sheet, header=None)

    def detect_header_row(df: pd.DataFrame, max_rows: int) -> int:
        for i in range(min(max_rows, len(df))):
            row = df.iloc[i].astype(str).str.lower().tolist()
            has_time = any("time" in s for s in row)
            has_mass_or_cond = any(("mass" in s) or ("cond" in s) for s in row)
            if has_time and has_mass_or_cond:
                return i
        raise ValueError("Could not find a header row containing both time and mass/cond.")

    hdr = detect_header_row(raw, max_header_scan)

    df = raw.iloc[hdr + 1:].copy()
    df.columns = raw.iloc[hdr].astype(str).str.strip().tolist()
    df = df.reset_index(drop=True)

    return df, hdr


# =============================================================================
# SECTION 1 — PyDataStru / PyVial / PyDataConfig
# =============================================================================

@dataclass
class PyDataConfig:
    """Python mirror of MATLAB data_stru.data_config."""
    M_F0: Optional[float] = None
    C_F0: Optional[np.ndarray] = None
    M_O: Optional[float] = None
    C_D: Optional[np.ndarray] = None
    nc: Optional[int] = None
    namec: Any = None
    ni: Optional[int] = None
    delP: Optional[float] = None
    Temp: Optional[float] = None
    Am: Optional[float] = None
    rho: Optional[float] = None
    Lp0: Optional[float] = None
    B0: Optional[np.ndarray] = None
    sigma0: Optional[float] = None
    theta0: Optional[np.ndarray] = None
    n: Optional[int] = None
    nr: Optional[int] = None
    n_extra: Optional[int] = None
    extras: Dict[str, Any] = field(default_factory=dict)

    '''
    meta = _extract_excel_metadata(...)
    cfg_kwargs = {}
    rpm = meta.pop("stirring_rpm", None)
    cfg_kwargs["stirring_rpm"] = rpm
    data_config = PyDataConfig(**cfg_kwargs, extras=meta)
    '''

@dataclass
class PyVial:
    """Python mirror of MATLAB data_stru.data_raw(i)."""
    number: int
    time: np.ndarray
    mass: np.ndarray
    cV_avg: Optional[np.ndarray] = None
    cF_exp: Optional[np.ndarray] = None
    nr: Optional[int] = None
    extras: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PyDataStru:
    """
    Unified experiment structure, mirroring MATLAB data_stru.

    Top-level:
        dataset, filename, mode, conductivity_cF, n_obj
        data_config : PyDataConfig
        data_raw    : List[PyVial]
    """
    dataset: Optional[float]
    filename: Optional[str]
    mode: Optional[str]
    conductivity_cF: Optional[bool]
    n_obj: Optional[int]
    data_config: PyDataConfig
    data_raw: List[PyVial]

    # -------------------------------------------------------------------------
    # Constructor from MATLAB data_stru .mat
    # -------------------------------------------------------------------------
    @staticmethod
    def from_mat(path: str) -> "PyDataStru":
        mat = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
        if "data_stru" not in mat:
            raise ValueError(f"{path} does not contain 'data_stru'")

        ds = mat["data_stru"]

        dataset = float(getattr(ds, "dataset", np.nan))
        filename = str(getattr(ds, "filename", None))
        mode = str(getattr(ds, "mode", None))
        conductivity_cF = bool(getattr(ds, "conductivity_cF", False))
        n_obj = int(getattr(ds, "n_obj", 0))

        # data_config
        dc = getattr(ds, "data_config", None)
        cfg_kwargs: Dict[str, Any] = {}
        extras: Dict[str, Any] = {}

        if dc is not None and hasattr(dc, "_fieldnames"):
            for fn in dc._fieldnames:
                val = getattr(dc, fn)
                if isinstance(val, (list, np.ndarray)):
                    val = np.asarray(val).squeeze()
                if fn in PyDataConfig.__dataclass_fields__:
                    cfg_kwargs[fn] = val
                else:
                    extras[fn] = val

        data_config = PyDataConfig(**cfg_kwargs, extras=extras)

        # data_raw vials
        dr = getattr(ds, "data_raw", [])
        if isinstance(dr, np.ndarray):
            dr = dr.ravel()

        vials: List[PyVial] = []
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

            vial_extras: Dict[str, Any] = {}
            if hasattr(v, "_fieldnames"):
                for fn in v._fieldnames:
                    if fn not in ["time", "mass", "cV_avg", "cF_exp", "nr", "number"]:
                        vial_extras[fn] = getattr(v, fn)

            vials.append(PyVial(num, t, m, cV, cF, nr, vial_extras))

        return PyDataStru(dataset, filename, mode, conductivity_cF, n_obj, data_config, vials)

    # -------------------------------------------------------------------------
    # Constructor from Excel sheet
    # -------------------------------------------------------------------------
    @staticmethod
    def from_excel(path: str, sheet: str) -> "PyDataStru":
        """
        Robust Excel → PyDataStru converter.

        Uses the same header-row detection logic as ExcelSource.
        Currently builds ONE PyVial for the sheet (no vial segmentation yet).
        """
        df, hdr = _excel_get_meas_df(path, sheet, max_header_scan=25)

        # Helper to find column by keywords
        def find_col(keywords: List[str]) -> Optional[str]:
            for c in df.columns:
                lc = str(c).lower()
                if any(k in lc for k in keywords):
                    return c
            return None

        tcol = find_col(["time", "sec", "seconds"])
        mcol = find_col(["mass", "m (g)", "permeate mass", "collection mass", "weight"])

        if tcol is None:
            raise ValueError(f"Time column not found. Columns = {list(df.columns)}")
        if mcol is None:
            raise ValueError(f"Mass column not found. Columns = {list(df.columns)}")

        time = pd.to_numeric(df[tcol], errors="coerce").to_numpy()
        mass = pd.to_numeric(df[mcol], errors="coerce").to_numpy()

        # Optional conductivity signals (case-insensitive)
        cF_df = df.filter(regex=re.compile("retent.*cond", flags=re.I))
        cV_df = df.filter(regex=re.compile("permeat.*cond", flags=re.I))


        cF = cF_df.to_numpy().squeeze() if not cF_df.empty else None
        cV = cV_df.to_numpy().squeeze() if not cV_df.empty else None

        vial = PyVial(
            number=1,
            time=time,
            mass=mass,
            cF_exp=cF,
            cV_avg=cV,
            nr=int(np.isfinite(mass).sum()),
            extras={},
        )

        return PyDataStru(
            dataset=None,
            filename=os.path.basename(path),
            mode=None,
            conductivity_cF=False,
            n_obj=1,
            data_config=PyDataConfig(),
            data_raw=[vial],
        )


# =============================================================================
# SECTION 2 — DataRun (core representation for DataLoader)
# =============================================================================

@dataclass
class DataRun:
    time: np.ndarray
    mass: np.ndarray
    retentate_cond: np.ndarray
    permeate_cond: np.ndarray
    vial_marker: Optional[np.ndarray] = None
    extra_signals: Dict[str, np.ndarray] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        n = len(np.asarray(self.time).squeeze())

        def pad(x: Any) -> np.ndarray:
            if x is None:
                return np.full(n, np.nan)
            arr = np.asarray(x).squeeze()
            if arr.size < n:
                arr = np.pad(arr, (0, n - arr.size), constant_values=np.nan)
            return arr[:n]

        self.time = pad(self.time)
        self.mass = pad(self.mass)
        self.retentate_cond = pad(self.retentate_cond)
        self.permeate_cond = pad(self.permeate_cond)

        if self.vial_marker is not None:
            self.vial_marker = pad(self.vial_marker)

        for k, v in list(self.extra_signals.items()):
            self.extra_signals[k] = pad(v)


# =============================================================================
# SECTION 3 — DataSource protocol, ExcelSource, MatSource
# =============================================================================

@runtime_checkable
class DataSource(Protocol):
    def sniff(self) -> bool: ...
    def parse(self) -> DataRun: ...


class ExcelSource:
    """
    Excel parser used by DataLoader.

    Uses the same header detection logic as PyDataStru.from_excel,
    but returns a flat DataRun for the whole sheet (one window).
    """
    def __init__(self, path: str, measure_sheet: Optional[str] = None):
        self.path = path
        self.measure_sheet = measure_sheet

    def sniff(self) -> bool:
        ext = os.path.splitext(self.path)[1].lower()
        return ext in {".xlsx", ".xls", ".xlsm"}

    def parse(self) -> DataRun:
        xl = pd.ExcelFile(self.path)
        sheet = self.measure_sheet or xl.sheet_names[0]

        df, hdr = _excel_get_meas_df(self.path, sheet, max_header_scan=25)

        def find_col(keywords: List[str]) -> Optional[str]:
            for c in df.columns:
                lc = str(c).lower()
                if any(k in lc for k in keywords):
                    return c
            return None

        tcol = find_col(["time", "sec", "seconds"])
        mcol = find_col(["mass", "m (g)", "permeate mass", "collection mass", "weight"])

        if tcol is None:
            raise ValueError(f"Time column not found in ExcelSource. Columns = {list(df.columns)}")
        if mcol is None:
            raise ValueError(f"Mass column not found in ExcelSource. Columns = {list(df.columns)}")

        time = pd.to_numeric(df[tcol], errors="coerce").to_numpy()
        mass = pd.to_numeric(df[mcol], errors="coerce").to_numpy()

        cF_df = df.filter(regex=re.compile("retent.*cond", flags=re.I))
        cV_df = df.filter(regex=re.compile("permeat.*cond", flags=re.I))


        r = cF_df.to_numpy().squeeze() if not cF_df.empty else None
        p = cV_df.to_numpy().squeeze() if not cV_df.empty else None

        metadata = {"sheet_name": sheet, "header_row_index": hdr}

        return DataRun(
            time=time,
            mass=mass,
            retentate_cond=r,
            permeate_cond=p,
            vial_marker=None,
            extra_signals={},
            metadata=metadata,
        )


class MatSource:
    """
    MATLAB parser.

    Uses PyDataStru.from_mat to fully read data_stru, then flattens
    all vials into a single DataRun with a vial_marker field.
    """
    def __init__(self, path: str):
        self.path = path

    def sniff(self) -> bool:
        return self.path.lower().endswith(".mat")

    def parse(self) -> DataRun:
        ds = PyDataStru.from_mat(self.path)

        times: List[np.ndarray] = []
        masses: List[np.ndarray] = []
        rcons: List[np.ndarray] = []
        pcons: List[np.ndarray] = []
        markers: List[np.ndarray] = []

        for vial in ds.data_raw:
            t = np.asarray(vial.time, dtype=float).squeeze()
            m = np.asarray(vial.mass, dtype=float).squeeze()
            r = np.asarray(vial.cF_exp, dtype=float).squeeze() if vial.cF_exp is not None else np.full_like(t, np.nan)
            p = np.asarray(vial.cV_avg, dtype=float).squeeze() if vial.cV_avg is not None else np.full_like(t, np.nan)
            vnum = float(vial.number)

            n = t.size
            if m.size < n:
                m = np.pad(m, (0, n - m.size), constant_values=np.nan)
            if r.size < n:
                r = np.pad(r, (0, n - r.size), constant_values=np.nan)
            if p.size < n:
                p = np.pad(p, (0, n - p.size), constant_values=np.nan)

            times.append(t)
            masses.append(m)
            rcons.append(r)
            pcons.append(p)
            markers.append(np.full(n, vnum))

        time_all = np.concatenate(times)
        mass_all = np.concatenate(masses)
        r_all = np.concatenate(rcons)
        p_all = np.concatenate(pcons)
        vm_all = np.concatenate(markers)

        metadata = {
            "dataset": ds.dataset,
            "filename": ds.filename,
            "mode": ds.mode,
            "conductivity_cF": ds.conductivity_cF,
            "n_obj": ds.n_obj,
        }

        return DataRun(
            time=time_all,
            mass=mass_all,
            retentate_cond=r_all,
            permeate_cond=p_all,
            vial_marker=vm_all,
            extra_signals={},
            metadata=metadata,
        )


# =============================================================================
# SECTION 4 — DataLoader
# =============================================================================

class DataLoader:
    """
    High-level loader:

        ldr = DataLoader(path, measure_sheet=optional_sheet_name)
        ldr.load()

        ldr.times          : List[np.ndarray] (per window)
        ldr.masses         : List[np.ndarray]
        ldr.retentate_cond : List[np.ndarray]
        ldr.permeate_cond  : List[np.ndarray]
        ldr.windows        : List[(start_idx, end_idx)]

    For Excel:
        - currently one window (whole sheet)

    For MAT (data_stru):
        - windows segmented by vial_marker changes
    """
    def __init__(self, path: str, measure_sheet: Optional[str] = None):
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        self.path = path
        self.measure_sheet = measure_sheet

        self.source: Optional[DataSource] = None
        self.windows: List[Tuple[int, int]] = []
        self.times: List[np.ndarray] = []
        self.masses: List[np.ndarray] = []
        self.retentate_cond: List[np.ndarray] = []
        self.permeate_cond: List[np.ndarray] = []

    def _select_source(self) -> None:
        excel = ExcelSource(self.path, self.measure_sheet)
        mat = MatSource(self.path)
        if excel.sniff():
            self.source = excel
        elif mat.sniff():
            self.source = mat
        else:
            raise ValueError(f"Unsupported file type: {self.path}")

    def load(self) -> None:
        self._select_source()
        assert self.source is not None

        run = self.source.parse()

        # Windowing logic
        if run.vial_marker is not None:
            vm = np.asarray(run.vial_marker).astype(float).squeeze()
            boundaries = [0]
            for i in range(1, len(vm)):
                if vm[i] != vm[i - 1]:
                    boundaries.append(i)
            boundaries.append(len(vm))

            self.windows = [(boundaries[i], boundaries[i + 1] - 1) for i in range(len(boundaries) - 1)]
        else:
            n = len(run.time)
            self.windows = [(0, n - 1)] if n > 0 else []

        self.times = []
        self.masses = []
        self.retentate_cond = []
        self.permeate_cond = []

        for (s, e) in self.windows:
            sl = slice(s, e + 1)
            self.times.append(run.time[sl])
            self.masses.append(run.mass[sl])
            self.retentate_cond.append(run.retentate_cond[sl])
            self.permeate_cond.append(run.permeate_cond[sl])

    def plot_mass(self) -> None:
        for t, m in zip(self.times, self.masses):
            plt.plot(t, m)
        plt.xlabel("Time [s]")
        plt.ylabel("Mass [g]")
        plt.title("Mass vs Time")
        plt.show()

    def plot_conductivity(self) -> None:
        for t, r in zip(self.times, self.retentate_cond):
            plt.plot(t, r, label="Retentate")
        for t, p in zip(self.times, self.permeate_cond):
            plt.plot(t, p, "--", label="Permeate")
        plt.xlabel("Time [s]")
        plt.ylabel("Conductivity [arb]")
        plt.title("Conductivity vs Time")
        plt.legend()
        plt.show()


# =============================================================================
# SECTION 5 — Utility: list_sheets + CLI
# =============================================================================

def list_sheets(path: str) -> List[str]:
    """Return list of sheet names if path is an Excel file; else empty list."""
    ext = os.path.splitext(path)[1].lower()
    if ext not in {".xlsx", ".xls", ".xlsm"}:
        return []
    xl = pd.ExcelFile(path)
    return list(xl.sheet_names)


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Experimental data loader test harness")
    parser.add_argument("path", help="Path to .xlsx or .mat")
    parser.add_argument("--sheet", default=None, help="Excel sheet name (optional)")
    parser.add_argument("--quicklook", action="store_true", help="Plot mass and cond quicklook")
    args = parser.parse_args()

    ldr = DataLoader(args.path, measure_sheet=args.sheet)
    ldr.load()
    print("Windows:", ldr.windows)
    print("Number of windows:", len(ldr.windows))

    if args.quicklook:
        ldr.plot_mass()
        ldr.plot_conductivity()


if __name__ == "__main__":
    _main()
